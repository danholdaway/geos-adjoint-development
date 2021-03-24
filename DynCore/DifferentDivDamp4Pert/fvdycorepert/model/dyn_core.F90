module dyn_core_mod

  use fv_arrays_mod,      only: g_precision, p_precision, pg_precision
  use mpp_domains_mod,    only: CGRID_NE, DGRID_NE, mpp_get_boundary,   &
                                mpp_update_domains
  use fv_mp_mod,          only: domain, isd, ied, jsd, jed, is, ie, js, je, gid, mp_barrier, mp_reduce_max
  use fv_control_mod,     only: hord_mt, hord_vt, hord_tm, hord_dp, hord_ze, hord_tr, n_sponge, m_riem, m_split, &
                                dddmp, fv_d2_bg=>d2_bg, fv_d4_bg=>d4_bg, d_ext, vtdm4, beta1, beta, init_wind_m, m_grad_p, &
                                a2b_ord, ppm_limiter, master, fv_debug, d_con, fv_nord=>nord, max_courant_no, &
                                no_cgrid, fill_dp, nwat, inline_q, breed_vortex_inline, shallow_water, spin_up_hours
  use sw_core_mod,        only: c_sw, d_sw, divergence_corner, d2a2c
  use a2b_edge_mod,       only: a2b_ord2, a2b_ord4
  use nh_core_mod,        only: Riem_Solver_C, Riem_Solver, update_dz_c, update_dz_d
  use fv_grid_tools_mod,  only: rdx, rdy, rdxc, dxc, dyc, rdyc, dx, dy, area, rarea, grid_type
  use fv_grid_utils_mod,  only: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n,  &
                                ec1, ec2, en1, en2, da_min_c
  use fv_timing_mod,      only: timing_on, timing_off
  use fv_diagnostics_mod, only: prt_maxmin
  use mpp_parameter_mod,  only: CORNER
  use mpp_mod,            only: FATAL, mpp_error

  use test_cases_mod,     only: test_case, case9_forcing1, case9_forcing2
#ifndef MAPL_MODE
  use nwp_nudge_mod,      only: breed_slp_inline
#endif

implicit none

private
public :: dyn_core

    real  , allocatable::  delzc(:,:,:)
    real  , allocatable::     ut(:,:,:)
    real  , allocatable::     vt(:,:,:)
    real  , allocatable::    crx(:,:,:), xfx(:,:,:)
    real  , allocatable::    cry(:,:,:), yfx(:,:,:)
    real  , allocatable:: divg_d(:,:,:)
    real  , allocatable::     zh(:,:,:)
    real(pg_precision), allocatable::     du(:,:,:), dv(:,:,:)
    real(pg_precision), allocatable::    pkc(:,:,:)
    real              , allocatable::    pkd(:,:,:)
    real              , allocatable::  delpc(:,:,:)
    real(pg_precision), allocatable::    pk3(:,:,:)
    real              , allocatable::    ptc(:,:,:)
    real(pg_precision), allocatable::     gz(:,:,:)

    real   , save :: myStep = 0
    real   , save :: myTime = 0
    logical, save :: first_step = .true.

contains

!-----------------------------------------------------------------------
!     dyn_core :: FV Lagrangian dynamics driver
!-----------------------------------------------------------------------
 
 subroutine dyn_core(npx, npy, npz, ng, sphum, nq, bdt, n_split, zvir, cp, akap, rdgas, grav, hydrostatic,  &
                     u,  v,  um, vm, w, delz, pt, q, delp, pe, pk, phis, omga, ptop, pfull, ua, va, & 
                     uc, vc, mfx, mfy, cx, cy, pem, pkz, peln, ak, bk, init_step, end_step, time_total)
    integer, intent(IN) :: npx
    integer, intent(IN) :: npy
    integer, intent(IN) :: npz
    integer, intent(IN) :: ng, nq, sphum
    integer, intent(IN) :: n_split
    real   , intent(IN) :: bdt
    real   , intent(IN) :: zvir, cp, akap, rdgas, grav
    real(p_precision) , intent(IN) :: ptop
    real, intent(in), dimension(npz+1):: ak, bk
    logical, intent(IN) :: hydrostatic
    logical, intent(IN) :: init_step, end_step
    real, intent(in) :: pfull(npz)
    real, intent(inout), dimension(isd:ied,jsd:jed+1,npz):: um ! D grid zonal wind (m/s)
    real, intent(inout), dimension(isd:ied+1,jsd:jed,npz):: vm ! D grid meridional wind (m/s)
    real, intent(inout), dimension(isd:ied,jsd:jed+1,npz):: u! D grid zonal wind (m/s)
    real, intent(inout), dimension(isd:ied+1,jsd:jed,npz):: v! D grid meridional wind (m/s)
    real, intent(inout) :: w(   isd:ied  ,jsd:jed  ,npz)  ! vertical vel. (m/s)
    real, intent(inout) :: delz(is :ie   ,js :je   ,npz)  ! delta-height (m)
    real, intent(inout) :: pt(  isd:ied  ,jsd:jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(isd:ied  ,jsd:jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) ::  q(  isd:ied  ,jsd:jed  ,npz, nq)  !
    real, intent(IN), optional:: time_total  ! total time (seconds) since start

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout):: phis(isd:ied,jsd:jed)      ! Surface geopotential (g*Z_surf)
    real(p_precision), intent(inout):: pe(is-1:ie+1, npz+1,js-1:je+1)  ! edge pressure (pascal)
    real, intent(out) :: pem(is-1:ie+1, npz+1,js-1:je+1)
    real(p_precision), intent(out):: peln(is:ie,npz+1,js:je)           ! ln(pe)
    real(p_precision), intent(inout):: pk(is:ie,js:je, npz+1)        ! pe**kappa

!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, intent(inout):: omga(isd:ied,jsd:jed,npz)    ! Vertical pressure velocity (pa/s)
    real, intent(inout):: uc(isd:ied+1,jsd:jed  ,npz)  ! (uc, vc) are mostly used as the C grid winds
    real, intent(inout):: vc(isd:ied  ,jsd:jed+1,npz)
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: ua, va

! The Flux capacitors: accumulated Mass flux arrays
    real, intent(inout)::  mfx(is:ie+1, js:je,   npz)
    real, intent(inout)::  mfy(is:ie  , js:je+1, npz)
! Accumulated Courant number arrays
    real, intent(inout)::  cx(is:ie+1, jsd:jed, npz)
    real, intent(inout)::  cy(isd:ied ,js:je+1, npz)
! Work:
    real(p_precision), intent(inout):: pkz(is:ie,js:je,npz)  ! 

! Auto 1D & 2D arrays:
   !real   wbuffer(npy+2,npz)
   !real   ebuffer(npy+2,npz)
   !real   nbuffer(npx+2,npz)
   !real   sbuffer(npx+2,npz)
    real   wbuffer(js:je,npz)
    real   sbuffer(is:ie,npz)
    real   ebuffer(js:je,npz)
    real   nbuffer(is:ie,npz)
! ----   For external mode:
    real(pg_precision) divg2(is:ie+1,js:je+1)
    real   wk(isd:ied,jsd:jed)
! --- For no_cgrid option ---
    real   u2(isd:ied,jsd:jed+1)
    real   v2(isd:ied+1,jsd:jed)
!-------------------------------------

    real :: d4_bg, d2_bg
    integer :: hord_m, hord_v, hord_t, hord_p, nord, nord_k
    integer :: i,j,k, it, iq
    integer :: ism1, iep1, jsm1, jep1
    integer :: ieb1, jeb1
    real    :: alpha, damp_k
    real    :: dt2, dt, rdt, rgrav
    real    :: d2_divg, dd_divg
    logical :: do_omega
    real(p_precision)  :: ptk

   !--------------------------------------------
   ! sub-cycle again to recover from instability
   !--------------------------------------------
    real :: cmax
    real :: ndt, elapsed_dt
    integer :: subt
    integer :: subcycle

      d2_bg = fv_d2_bg
      d4_bg = fv_d4_bg
      nord = fv_nord
      if (myTime < 3600.0*spin_up_hours) then
         nord = 0
         d2_bg = 0.075
         d4_bg = 0.0
         myStep = myStep+1.0
         myTime = myStep*bdt
         if (gid==0) print*, 'Spinning Up with NORD:', nord, 'myTime:', myTime
      endif

      ptk = ptop** akap
      ndt = bdt / real(n_split)
      dt  = ndt
      dt2 = 0.5*dt
      rdt = 1./dt
      rgrav = 1./grav

! Indexes:
      ism1 = is - 1;  iep1 = ie + 1
      jsm1 = js - 1;  jep1 = je + 1
                                                 call timing_on('COMM_TOTAL')
      if ( npz>1 )   &
      call mpp_update_domains( pt,   domain, complete=.false. )
      call mpp_update_domains( delp, domain, complete=.true. )
      call mpp_update_domains( u, v, domain, gridtype=DGRID_NE, complete=.true. )

                                                call timing_off('COMM_TOTAL')

      if ( init_step ) then

           allocate(    gz(isd:ied, jsd:jed ,npz+1) )
           allocate(   ptc(isd:ied, jsd:jed ,npz ) )
           allocate( delzc(is:ie, js:je ,npz ) )
           allocate( crx(is :ie+1, jsd:jed,  npz) )
           allocate( xfx(is :ie+1, jsd:jed,  npz) )
           allocate( cry(isd:ied,  js :je+1, npz) )
           allocate( yfx(isd:ied,  js :je+1, npz) )
           allocate( divg_d(isd:ied+1,jsd:jed+1,npz) )
           allocate(   pkc(isd:ied, jsd:jed  ,npz+1) )
           allocate( delpc(isd:ied, jsd:jed  ,npz  ) )
           allocate(   pkd(isd:ied, jsd:jed  ,npz+1) )

          if ( .not. no_cgrid ) then
               allocate( ut(isd:ied, jsd:jed, npz) )
               allocate( vt(isd:ied, jsd:jed, npz) )
               ut(:,:,:) = 0.
               vt(:,:,:) = 0.
          endif
          if ( .not. hydrostatic ) then 
               allocate( zh(isd:ied, jsd:jed, npz) )
               if ( m_grad_p==0 ) allocate ( pk3(isd:ied,jsd:jed,npz+1) )
          endif

          if ( beta > 1.E-4 .or. no_cgrid ) then
               allocate( du(isd:ied,  jsd:jed+1,npz) )
               allocate( dv(isd:ied+1,jsd:jed,  npz) )
          endif
      endif

     if ( beta > 1.E-4 .or. (no_cgrid .and. init_wind_m) ) then
          call geopk(ptop, pe, peln, delp, pkc, gz, phis, pt, pkz, npz, akap, .false.)
          call grad1_p(du, dv, pkc, gz, delp, dt, ng, npx, npy, npz, ptop, ptk, hydrostatic)

          if( init_wind_m )   &
          call mpp_update_domains(du, dv, domain, gridtype=DGRID_NE)
     endif

! Empty the "flux capacitors"
     mfx(:,:,:) = 0.;  mfy(:,:,:) = 0.
      cx(:,:,:) = 0.;   cy(:,:,:) = 0.

  elapsed_dt = 0.0
  subcycle=1
  if (.not. hydrostatic) then
     m_split = max(1, CEILING(0.5 + abs(dt)/3.0))
    !m_riem = 0
    !if(m_split >= 2) m_riem = 2
    !if(m_split >= 4) m_riem = 4
    !if(gid==0) print*, 'm_split: ', m_split, 'm_riem: ', m_riem
  endif
!-----------------------------------------------------
  do it=1,n_split
!-----------------------------------------------------

   !--------------------------------------------
   ! Dynamic sub-cycling to prevent instability
   !--------------------------------------------
   if (max_courant_no >= 0.0) then
    cmax = 0.0
    do k=1,npz
      do j=js,je
        do i=is,ie+1
          crx(i,j,k) = v(i,j,k)*ndt*rdy(i,j)  ! D-Grid meridonal courant number stored in crx
          cmax = MAX(cmax, ABS(crx(i,j,k)))
        enddo
      enddo
      do j=js,je+1
        do i=is,ie
          cry(i,j,k) = u(i,j,k)*ndt*rdx(i,j)  ! D-Grid zonal courant number stored in cry
          cmax = MAX(cmax, ABS(cry(i,j,k)))
        enddo
      enddo
    enddo
    call mp_reduce_max(cmax)
    call mp_barrier()
    subt = CEILING(cmax/max_courant_no)
    if ( subcycle /= subt) then
       subcycle = subt
       dt  = (bdt - elapsed_dt) / real((n_split-(it-1))*subcycle)
       dt2 = 0.5*dt
       rdt = 1./dt
       if (gid==0) write(*,101) 'Sub-Cycling for stability = ', it, subcycle, dt, cmax
101    format(A,i3,2x,i3,2x,f10.6,2x,f10.6)
       if (.not. hydrostatic) then
         m_split = max(1, CEILING(0.5 + abs(dt)/3.0))
        !m_riem = 0
        !if(m_split >= 2) m_riem = 2
        !if(m_split >= 4) m_riem = 4
        !if(gid==0) print*, 'm_split: ', m_split, 'm_riem: ', m_riem
       endif
    endif
   endif

    do subt=1,subcycle
   !--------------------------------------------

     if ( .not. hydrostatic ) then
        do j=js,je
           do i=is,ie
              zh(i,j,npz) = phis(i,j)*rgrav - delz(i,j,npz)
           enddo
           do k=npz-1,1,-1
              do i=is,ie
                 zh(i,j,k) = zh(i,j,k+1) - delz(i,j,k)
              enddo
           enddo
        enddo
                                 call timing_on('COMM_TOTAL')
        call mpp_update_domains(zh, domain, complete=.false.)
        call mpp_update_domains(w,  domain, complete=.true.)
                                call timing_off('COMM_TOTAL')
     endif

     if (shallow_water) then
        do_omega  = .false.
        if (test_case==9) call case9_forcing1(phis, time_total)
     else
        if ( ((it==n_split) .and. (subt==subcycle)) ) then
!$omp parallel do default(shared) private(i, j, k)
             do j=jsm1,jep1
                do i=ism1,iep1
                   pem(i,1,j) = ptop
                enddo
                do k=1,npz
                   do i=ism1,iep1
                      pem(i,k+1,j) = pem(i,k,j) + delp(i,j,k)
                   enddo
               enddo
             enddo
             do_omega  = .true.
        else
             do_omega  = .false.
        endif
     endif

     if ( no_cgrid ) then
!---------------------------------------------------------------
! Using time extrapolated wind in place of computed C Grid wind
!---------------------------------------------------------------
                                               call timing_on('no_cgrid')
        do k=1,npz
           if ( init_wind_m ) then
              do j=jsd,jed+1
                 do i=isd,ied
                    u2(i,j) = u(i,j,k) + 0.5*du(i,j,k)
                 enddo
              enddo
              do j=jsd,jed
                 do i=isd,ied+1
                    v2(i,j) = v(i,j,k) + 0.5*dv(i,j,k)
                 enddo
              enddo
           else

           do j=jsd,jed+1
              do i=isd,ied
                 u2(i,j) = 1.5*u(i,j,k) - 0.5*um(i,j,k)
              enddo
           enddo
           do j=jsd,jed
              do i=isd,ied+1
                 v2(i,j) = 1.5*v(i,j,k) - 0.5*vm(i,j,k)
              enddo
           enddo
           endif
           call d2a2c( u(isd,jsd,k),  v(isd,jsd,k), u2, v2, ua(isd,jsd,k),   &
                      va(isd,jsd,k), uc(isd,jsd,k), vc(isd,jsd,k), nord>0 )
        enddo

        if ( .not. hydrostatic ) delpc(:,:,:) = delp(:,:,:)
        um(:,:,:) = u(:,:,:)
        vm(:,:,:) = v(:,:,:)
                                               call timing_off('no_cgrid')
    else
                                                     call timing_on('c_sw')
!$omp parallel do default(shared) private(i,j,k)
      do k=1,npz
         call c_sw(delpc(isd,jsd,k), delp(isd,jsd,k),  ptc(isd,jsd,k),  &
                      pt(isd,jsd,k),    u(isd,jsd,k),    v(isd,jsd,k),  &
                       w(isd,jsd,k),   uc(isd,jsd,k),   vc(isd,jsd,k),  &
                      ua(isd,jsd,k),   va(isd,jsd,k), omga(isd,jsd,k),  &
                      ut(isd,jsd,k),   vt(isd,jsd,k), dt2, hydrostatic, nord>0 )
! on output omga is updated w
      enddo
                                                     call timing_off('c_sw')

      if( fill_dp ) call mix_dp(hydrostatic, omga, delpc, ptc, npz, ak, bk, .true.)

      if ( hydrostatic ) then
           if ( beta1 > 0.001 ) then
                      alpha = 1. - beta1
               delpc(:,:,:) = beta1*delp(:,:,:) + alpha*delpc(:,:,:)  
                 ptc(:,:,:) = beta1*  pt(:,:,:) + alpha*  ptc(:,:,:)  
           endif
           call geopk(ptop, pe, peln, delpc, pkc, gz, phis, ptc, pkz, npz, akap, .true.)
      else
           call update_dz_c(is,   ie, js, je,  npz,    ng,    &
                            area, zh, ut, vt, delz, delzc, gz)
                                               call timing_on('Riem_C')
           call Riem_Solver_C( dt2,   is,  ie,   js,   je,   npz,   ng,   &
                               akap, rdgas, grav,  cp,  ptop, phis, omga, delzc, ptc,  &
                               delpc, gz,  pkc,  1 )
                                               call timing_off('Riem_C')
! pkc is full non-hydro pressure
                                               call timing_on('COMM_TOTAL')
!          call mpp_update_domains(pkc,domain,whalo=1,ehalo=1,shalo=1, nhalo=1, complete=.false.)
!          call mpp_update_domains(gz ,domain,whalo=1,ehalo=1,shalo=1, nhalo=1, complete=.true.)
           call mpp_update_domains(pkc, domain, complete=.false.)
           call mpp_update_domains(gz , domain, complete=.true.)
                                               call timing_off('COMM_TOTAL')
      endif

!-----------------------------------------
! Update time-centered winds on the C-Grid
!-----------------------------------------
      ieb1 = ie+1;   jeb1 = je+1

!$omp parallel do default(shared) private(i, j, k, wk)
      do k=1,npz
         if ( hydrostatic ) then
              do j=jsm1,jeb1
                 do i=ism1,ieb1
                    wk(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
                 enddo
              enddo
         else
              do j=jsm1,jeb1
                 do i=ism1,ieb1
                       wk(i,j) = delpc(i,j,k)
                 enddo
              enddo
              delpc(:,:,k) = delp(:,:,k) ! Save delp for update_dz_d
         endif

         do j=js,je
            do i=is,ieb1
               uc(i,j,k) = uc(i,j,k) + dt2*rdxc(i,j) / (wk(i-1,j)+wk(i,j)) *   &
                      ( (gz(i-1,j,k+1)-gz(i,j,k  ))*(pkc(i,j,k+1)-pkc(i-1,j,k))  &
                      + (gz(i-1,j,k) - gz(i,j,k+1))*(pkc(i-1,j,k+1)-pkc(i,j,k)) )
            enddo
         enddo
         do j=js,jeb1
            do i=is,ie
               vc(i,j,k) = vc(i,j,k) + dt2*rdyc(i,j) / (wk(i,j-1)+wk(i,j)) *   &
                      ( (gz(i,j-1,k+1)-gz(i,j,k  ))*(pkc(i,j,k+1)-pkc(i,j-1,k))  &
                      + (gz(i,j-1,k) - gz(i,j,k+1))*(pkc(i,j-1,k+1)-pkc(i,j,k)) )
            enddo
         enddo
      enddo

    endif      ! end no_cgrid section
                                                     call timing_on('COMM_TOTAL')
      call mpp_update_domains(uc, vc, domain, gridtype=CGRID_NE, complete=.true.)
                                                     call timing_off('COMM_TOTAL')

    if (shallow_water) then
       if (test_case==9) call case9_forcing2(phis)
    endif
    if ( inline_q ) then
                                   call timing_on('COMM_TOTAL')
         call mpp_update_domains( q,  domain, complete=.true. )
                                   call timing_off('COMM_TOTAL')
    endif

    if ( nord>0 ) then
         call divergence_corner(u, v, ua, va, divg_d, npz)
                                            call timing_on('COMM_TOTAL')
         call mpp_update_domains(divg_d, domain, position=CORNER)
                                            call timing_off('COMM_TOTAL')
    endif

                                                     call timing_on('d_sw')
!$omp parallel do default(shared) private(i, j, k, nord_k, damp_k, d2_divg, dd_divg, hord_m, hord_v, hord_t, hord_p, wk)
    do k=1,npz
       hord_m = hord_mt
       hord_t = hord_tm
       hord_v = hord_vt
       hord_p = hord_dp
       nord_k = nord
       damp_k = dddmp
       d2_divg = min(0.20, d2_bg*(1.-3.*tanh(0.1*log(pfull(k)/pfull(npz)))))
#ifdef NEW
       if ( n_sponge==-1 .or. npz==1 ) then
! Constant divg damping coefficient:
           d2_divg = d2_bg
           dd_divg = d4_bg
       else
         dd_divg = d4_bg
         if ( n_sponge/=0 .and. k<=2*n_sponge) then
            dd_divg = MIN(0.20,d4_bg*(1.-3.*tanh(0.1*log(pfull(k)/pfull(2*n_sponge)))))
            hord_m = 2
            hord_v = 2
            hord_t = 2
            hord_p = 2
         endif
         if ( n_sponge/=0 .and. k<=n_sponge+1) then
            if (k<=n_sponge+1) then
               hord_m = 1
               hord_v = 1
               hord_t = 1
               hord_p = 1
               nord_k = 0
               damp_k = 0.2
               d2_divg = 0.20
            endif
         endif
         damp_k = max(damp_k, dddmp)
       endif
#else
       d2_divg = min(0.20, d2_bg*(1.-3.*tanh(0.1*log(pfull(k)/pfull(npz)))))
       if ( n_sponge==-1 .or. npz==1 ) then
! Constant divg damping coefficient:
           d2_divg = d2_bg
       else
         if ( n_sponge==0 .and. k==1 ) then
               hord_v = 2
               hord_t = 2
               hord_p = 2
               nord_k = max(0, nord-1)
               damp_k = 0.2
               d2_divg = min(0.20, 2.*d2_bg)
               d2_divg = max(0.06 ,  d2_divg)
         else
           if( k <= n_sponge .and. npz>16 ) then
! Apply first order scheme for damping the sponge layer
               hord_m = 1
               hord_v = 1
               hord_t = 1
               hord_p = 1
               nord_k = 0
               damp_k = 0.2
               d2_divg = min(0.20, 4.*d2_bg)   ! 0.25 is the stability limit
               d2_divg = max(0.15,  d2_divg)
           elseif( k == n_sponge+1 .and. npz>24 ) then
               hord_v = 2
               hord_t = 2
               hord_p = 2
               nord_k = max(0, nord-1)
               d2_divg = min(0.20, 2.*d2_bg)
               d2_divg = max(0.08, d2_divg)
               if ( nord > 1 ) then
                    damp_k = 0.
               else
                    damp_k = 0.12
               endif
           endif
         endif
       endif
       dd_divg = d4_bg
        damp_k = max(damp_k, dddmp)
#endif
     ! if (gid==0 .and. first_step) then
     !    print*, k, nord_k, damp_k, d2_divg, dd_divg 
     ! endif

!--- external mode divergence damping ---
       if ( d_ext > 0. ) then 
            call a2b_ord2(delp(:,:,k), wk, npx, npy, is, ie, js, je, ng, .false.)
       endif

       call d_sw(pkd(isd,jsd,k), delp(isd,jsd,k), ptc(isd,jsd,k),  pt(isd,jsd,k),    &
                   u(isd,jsd,k),    v(isd,jsd,k),   w(isd,jsd,k),  uc(isd,jsd,k),    &
                  vc(isd,jsd,k),   ua(isd,jsd,k),  va(isd,jsd,k), divg_d(isd,jsd,k), &
                 mfx(is, js, k),  mfy(is, js, k),  cx(is, jsd,k),  cy(isd,js, k),    &
                 crx(is, jsd,k),  cry(isd,js, k), xfx(is, jsd,k), yfx(isd,js, k),    &
                 zvir, sphum, nq, q, k, npz,   inline_q,        &
                 pkz(is,js,k),    dt, hord_tr, hord_m, hord_v,  &
                 hord_t, hord_p, nord_k, damp_k, d2_divg, dd_divg, vtdm4, &
                 d_con, hydrostatic, ppm_limiter )

       if ( d_ext > 0. ) then
            do j=js,jep1
               do i=is,iep1
                  ptc(i,j,k) = wk(i,j)   ! delp at cell corners
               enddo
            enddo
       endif
    enddo
                                                     call timing_off('d_sw')

    if( fill_dp ) call mix_dp(hydrostatic, w, delp, pt, npz, ak, bk, .false.)

    if ( fv_debug ) then
         call prt_maxmin('DELP', delp, is, ie, js, je, ng, npz, 1.E-2, master)
    endif

    if ( d_ext > 0. ) then
          d2_divg = d_ext * da_min_c
! pkd() is 3D field of horizontal divergence
! ptc is "delp" at cell corners
!$omp parallel do default(shared) private(i,j,k)
          do j=js,jep1
              do i=is,iep1
                    wk(i,j) = ptc(i,j,1)
                 divg2(i,j) = wk(i,j)*pkd(i,j,1)
              enddo
              do k=2,npz
                 do i=is,iep1
                       wk(i,j) =    wk(i,j) + ptc(i,j,k)
                    divg2(i,j) = divg2(i,j) + ptc(i,j,k)*pkd(i,j,k)
                 enddo
              enddo
              do i=is,iep1
                 divg2(i,j) = d2_divg*divg2(i,j)/wk(i,j)
              enddo
          enddo
    else
        divg2 = 0.
    endif
                             call timing_on('COMM_TOTAL')
    call mpp_update_domains(  pt, domain, complete=.false.)
    call mpp_update_domains(delp, domain, complete=.true.)
                             call timing_off('COMM_TOTAL')

     if ( hydrostatic ) then
          call geopk(ptop, pe, peln, delp, pkc, gz, phis, pt, pkz, npz, akap, .false.)
     else
                                            call timing_on('UPDATE_DZ')
          call update_dz_d(hord_tm, is, ie, js, je, npz, ng, npx, npy, area,  &
                           zh, crx, cry, xfx, yfx, delz, delzc, delpc, n_sponge)
                                            call timing_off('UPDATE_DZ')
                                                          call timing_on('Riem_D')
!-----------------------------------------------------------
! mgrad_p = 1: pkc is full pressure
! mgrad_p = 0: pkc is non-hydrostatic perturbation pressure
!-----------------------------------------------------------
          call Riem_Solver(dt,   is,   ie,   js,   je, npz,  ng,  &
                           akap, rdgas, grav, cp,   ptop, phis, peln, w,  delz,      &
                           pt,   delp, gz,   pkc,  pk, pe, ((it==n_split) .and. (subt==subcycle)), m_grad_p)
                                                 call timing_off('Riem_D')

                                       call timing_on('COMM_TOTAL')
          if ( m_grad_p==0 ) then
             do k=1,npz+1
                do j=js,je
                   do i=is,ie
                      pk3(i,j,k) = pk(i,j,k)
                   enddo
                enddo
             enddo
             call mpp_update_domains(pk3, domain, complete=.false.)
          endif

          call mpp_update_domains(pkc, domain, complete=.true.)
          call mpp_update_domains(gz , domain, complete=.true.)
                                       call timing_off('COMM_TOTAL')
     endif    ! end hydro case


     if (.not. shallow_water) then
      if ( breed_vortex_inline .or. (((it==n_split) .and. (subt==subcycle)) .and. hydrostatic) ) then
!$omp parallel do default(shared) private(i, j, k)
           do k=1,npz+1
              do j=js,je
                 do i=is,ie
                    pk(i,j,k) = pkc(i,j,k)
                 enddo
              enddo
           enddo
     endif
    
      if ( do_omega ) then
!------------------------------
! Compute time tendency
!------------------------------
         do k=1,npz
            do j=js,je
               do i=is,ie
                  omga(i,j,k) = (pe(i,k+1,j) - pem(i,k+1,j)) * rdt 
               enddo
            enddo
         enddo
!------------------------------
! Compute the "advective term"
!------------------------------
         call adv_pe(ua, va, pem, omga, npx, npy,  npz, ng)
      endif
    endif

      if ( .not.hydrostatic .and. m_grad_p == 0 ) then
           call two_grad_p(u, v, pkc, gz, delp, pk3, divg2, dt, ng, npx, npy,   &
                           npz, ptk)  
      else
       if ( beta > 1.E-4 ) then
          do k=1,npz
             do j=js,je+1
                do i=is,ie
                   u(i,j,k) = (u(i,j,k)+divg2(i,j)-divg2(i+1,j))*rdx(i,j) + beta*du(i,j,k)
                enddo
             enddo
          enddo
          do k=1,npz
             do j=js,je
                do i=is,ie+1
                   v(i,j,k) = (v(i,j,k)+divg2(i,j)-divg2(i,j+1))*rdy(i,j) + beta*dv(i,j,k)
                enddo
             enddo
          enddo
          call grad1_p(du, dv, pkc, gz, delp, dt, ng, npx, npy, npz, ptop, ptk, hydrostatic)
          alpha = 1. - beta
          do k=1,npz
             do j=js,je+1
                do i=is,ie
                   u(i,j,k) = u(i,j,k) + alpha*du(i,j,k)
                enddo
             enddo
          enddo
          do k=1,npz
             do j=js,je
                do i=is,ie+1
                   v(i,j,k) = v(i,j,k) + alpha*dv(i,j,k)
                enddo
             enddo
          enddo
       else
          call one_grad_p(u, v, pkc, gz, divg2, delp, dt, ng, npx, npy, npz,   &
                          ptop, ptk, hydrostatic)  
       endif
      endif

!-------------------------------------------------------------------------------------------------------
#ifndef MAPL_MODE
      if ( breed_vortex_inline )     &
      call breed_slp_inline(it, dt, npz, ak, bk, phis, pe, pk, peln, delp, u, v, pt, q, nwat, zvir)
#endif
!-------------------------------------------------------------------------------------------------------

                                                                call timing_on('COMM_TOTAL')
      if( ((it==n_split) .and. (subt==subcycle)) .and. grid_type<4 ) then
! Prevent accumulation of rounding errors at overlapped domain edges:
          call mpp_get_boundary(u, v, domain, wbuffery=wbuffer, ebuffery=ebuffer,  &
                               sbufferx=sbuffer, nbufferx=nbuffer, gridtype=DGRID_NE )
          do k=1,npz
             do i=is,ie
                u(i,je+1,k) = nbuffer(i,k)
             enddo
          enddo
          do k=1,npz
             do j=js,je
                v(ie+1,j,k) = ebuffer(j,k)
             enddo
          enddo
      else
          call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
      endif
                                                                call timing_off('COMM_TOTAL')
      init_wind_m = .false.
      first_step = .false.

   elapsed_dt = elapsed_dt + dt
   enddo  ! sub-cycle

!-----------------------------------------------------
  enddo   ! time split loop
!-----------------------------------------------------


   if ( end_step ) then
        deallocate(    gz )
        deallocate(   ptc )
        deallocate( delzc )
        deallocate(   crx )
        deallocate(   xfx )
        deallocate(   cry )
        deallocate(   yfx )
        deallocate ( divg_d )
        deallocate(   pkc )
        deallocate( delpc )
        deallocate(   pkd )

        if ( .not. no_cgrid ) then
              deallocate( ut )
              deallocate( vt )
        endif
        if ( .not. hydrostatic ) then
              deallocate( zh )
              if ( m_grad_p==0 ) deallocate ( pk3 )
        endif
        if ( beta > 1.E-4 .or. no_cgrid ) then
             deallocate( du )
             deallocate( dv )
        endif
   endif

 end subroutine dyn_core



 subroutine two_grad_p(u, v, pk, gh, delp, pkt, divg2, dt, ng, npx, npy, npz, ptk)  

    integer, intent(IN) :: ng, npx, npy, npz
    real  ,    intent(IN) :: dt
    real(p_precision), intent(IN) :: ptk
    real(pg_precision),    intent(in) :: divg2(is:ie+1, js:je+1)
    real  , intent(inout) ::  delp(isd:ied, jsd:jed, npz)
    real(pg_precision), intent(inout) ::    pk(isd:ied, jsd:jed, npz+1)  ! perturbation pressure
    real(pg_precision), intent(inout) ::   pkt(isd:ied, jsd:jed, npz+1)  ! p**kappa
    real(pg_precision), intent(inout) ::    gh(isd:ied, jsd:jed, npz+1)  ! g * zh
    real  , intent(inout) ::     u(isd:ied,  jsd:jed+1,npz) 
    real  , intent(inout) ::     v(isd:ied+1,jsd:jed,  npz)
! Local:
    real(pg_precision)   dp(isd:ied, jsd:jed)
    real                 wk1(isd:ied, jsd:jed)
    real                  wk(is: ie+1,js: je+1)

    integer iep1, jep1
    integer i,j,k

    iep1 = ie + 1
    jep1 = je + 1

    do j=js,jep1
       do i=is,iep1
           pk(i,j,1) = 0.
          pkt(i,j,1) = ptk
       enddo
    enddo

    do k=1,npz+1
       if ( k/=1 ) then
         if ( a2b_ord==4 ) then
           call a2b_ord4( pk(isd:,jsd:,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
           call a2b_ord4(pkt(isd:,jsd:,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
         else
           call a2b_ord2( pk(isd:,jsd:,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
           call a2b_ord2(pkt(isd:,jsd:,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
         endif
       endif

       if ( a2b_ord==4 ) then
           call a2b_ord4( gh(isd:,jsd:,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
       else
           call a2b_ord2( gh(isd:,jsd:,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
       endif
    enddo

    do k=1,npz
       if ( a2b_ord==4 ) then
            call a2b_ord4(delp(:,:,k), wk1, npx, npy, is, ie, js, je, ng)
       else
            call a2b_ord2(delp(:,:,k), wk1, npx, npy, is, ie, js, je, ng)
       endif

       do j=js,jep1
          do i=is,iep1
             wk(i,j) = pkt(i,j,k+1) - pkt(i,j,k)
          enddo
       enddo

       do j=js,jep1
          do i=is,ie
!------------------
! Perturbation term:
!------------------
             u(i,j,k) = u(i,j,k) + dt/(wk1(i,j)+wk1(i+1,j)) *   &
                   ((gh(i,j,k+1)-gh(i+1,j,k  ))*(pk(i+1,j,k+1)-pk(i,j,k  )) &
                  + (gh(i,j,k  )-gh(i+1,j,k+1))*(pk(i,j,k+1  )-pk(i+1,j,k)))
!-----------------
! Hydrostatic term
!-----------------
             u(i,j,k) = rdx(i,j)*(divg2(i,j)-divg2(i+1,j)+u(i,j,k) + dt/(wk(i,j)+wk(i+1,j)) *  &
                   ((gh(i,j,k+1)-gh(i+1,j,k  ))*(pkt(i+1,j,k+1)-pkt(i  ,j,k)) &
                  + (gh(i,j,k  )-gh(i+1,j,k+1))*(pkt(i  ,j,k+1)-pkt(i+1,j,k))))
          enddo
       enddo

       do j=js,je
          do i=is,iep1
!------------------
! Perturbation term:
!------------------
             v(i,j,k) = v(i,j,k) + dt/(wk1(i,j)+wk1(i,j+1)) *   &
                   ((gh(i,j,k+1)-gh(i,j+1,k  ))*(pk(i,j+1,k+1)-pk(i,j  ,k)) &
                  + (gh(i,j,k  )-gh(i,j+1,k+1))*(pk(i,j  ,k+1)-pk(i,j+1,k)))
!-----------------
! Hydrostatic term
!-----------------
             v(i,j,k) = rdy(i,j)*(divg2(i,j)-divg2(i,j+1)+v(i,j,k) + dt/(wk(i,j)+wk(i,j+1)) * &
                   ((gh(i,j,k+1)-gh(i,j+1,k  ))*(pkt(i,j+1,k+1)-pkt(i,j  ,k)) &
                  + (gh(i,j,k  )-gh(i,j+1,k+1))*(pkt(i,j  ,k+1)-pkt(i,j+1,k))))
          enddo
       enddo
    enddo    ! end k-loop

 end subroutine two_grad_p

 subroutine one_grad_p(u, v, pk, gh, divg2, delp, dt, ng, npx, npy, npz,  &
                       ptop, ptk,  hydrostatic)  

    integer, intent(IN) :: ng, npx, npy, npz
    real  ,    intent(IN) :: dt
    real(p_precision), intent(IN) :: ptop, ptk
    logical, intent(in) :: hydrostatic
    real(pg_precision),    intent(in) :: divg2(is:ie+1,js:je+1)
    real(pg_precision), intent(inout) ::    pk(isd:ied,  jsd:jed  ,npz+1)
    real(pg_precision), intent(inout) ::    gh(isd:ied,  jsd:jed  ,npz+1)
    real  , intent(inout) ::  delp(isd:ied,  jsd:jed  ,npz)
    real  , intent(inout) ::     u(isd:ied  ,jsd:jed+1,npz) 
    real  , intent(inout) ::     v(isd:ied+1,jsd:jed  ,npz)
! Local:
    real, dimension(isd:ied,jsd:jed):: wk
    real:: wk1(is:ie+1,js:je)
    real:: wk2(is:ie,js:je+1)
    real top_value
    integer :: iep1, jep1
    integer i,j,k

    real(pg_precision)  :: dt_pg

    dt_pg = dt

    iep1 = ie + 1
    jep1 = je + 1

    if ( hydrostatic ) then
! pk is pe**kappa if hydrostatic
         top_value = ptk
    else
! pk is full pressure if non-hydrostatic
         top_value = ptop
    endif

    do j=js,jep1
       do i=is,iep1
          pk(i,j,1) = top_value
       enddo
    enddo

    do k=2,npz+1
       if ( a2b_ord==4 ) then
         call a2b_ord4(pk(isd:,jsd:,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       else
         call a2b_ord2(pk(isd:,jsd:,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       endif
    enddo

    do k=1,npz+1
       if ( a2b_ord==4 ) then
         call a2b_ord4( gh(isd:,jsd:,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       else
         call a2b_ord2( gh(isd:,jsd:,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       endif
    enddo

    do j=js,jep1
       do i=is,ie
          wk2(i,j) = divg2(i,j)-divg2(i+1,j)
       enddo
    enddo
    do j=js,je
       do i=is,iep1
          wk1(i,j) = divg2(i,j)-divg2(i,j+1)
       enddo
    enddo

    do k=1,npz

       if ( hydrostatic ) then
            do j=js,jep1
               do i=is,iep1
                  wk(i,j) = pk(i,j,k+1) - pk(i,j,k)
               enddo
            enddo
       else
         if ( a2b_ord==4 ) then
            call a2b_ord4(delp(:,:,k), wk, npx, npy, is, ie, js, je, ng)
         else
            call a2b_ord2(delp(:,:,k), wk, npx, npy, is, ie, js, je, ng)
         endif
       endif

       do j=js,jep1
          do i=is,ie
             u(i,j,k) = rdx(i,j)*(wk2(i,j)+u(i,j,k) + dt_pg/(wk(i,j)+wk(i+1,j)) *  &
                         ((gh(i,j,k+1)-gh(i+1,j,k  ))*(pk(i+1,j,k+1)-pk(i  ,j,k)) &
                        + (gh(i,j,k  )-gh(i+1,j,k+1))*(pk(i  ,j,k+1)-pk(i+1,j,k))))
          enddo
       enddo
       do j=js,je
          do i=is,iep1
             v(i,j,k) = rdy(i,j)*(wk1(i,j)+v(i,j,k) + dt_pg/(wk(i,j)+wk(i,j+1)) *  &
                         ((gh(i,j,k+1)-gh(i,j+1,k  ))*(pk(i,j+1,k+1)-pk(i,j  ,k)) &
                        + (gh(i,j,k  )-gh(i,j+1,k+1))*(pk(i,j  ,k+1)-pk(i,j+1,k))))
          enddo
       enddo
    enddo    ! end k-loop

 end subroutine one_grad_p

 subroutine grad1_p(delu, delv, pk, gh,  delp, dt, ng, npx, npy, npz,  &
                    ptop, ptk,  hydrostatic)  

    integer, intent(in) :: ng, npx, npy, npz
    real  ,    intent(in) :: dt
    real(p_precision), intent(IN) :: ptk
    real(p_precision),    intent(in) :: ptop
    logical, intent(in) :: hydrostatic
    real(pg_precision), intent(inout) ::    pk(isd:ied,  jsd:jed  ,npz+1)
    real(pg_precision), intent(inout) ::    gh(isd:ied,  jsd:jed  ,npz+1)
    real  , intent(inout) ::  delp(isd:ied,  jsd:jed  ,npz)

    real(pg_precision), intent(out) ::    delu(isd:ied  ,jsd:jed+1,npz) 
    real(pg_precision), intent(out) ::    delv(isd:ied+1,jsd:jed  ,npz)
! Local:
    real:: wk(isd:ied,jsd:jed)
    real top_value
    integer i,j,k
    real(pg_precision) :: dt_pg
  
    dt_pg = dt

    if ( hydrostatic ) then
! pk is pe**kappa if hydrostatic
         top_value = ptk
    else
! pk is full pressure if non-hydrostatic
         top_value = ptop
    endif

    do j=js,je+1
       do i=is,ie+1
          pk(i,j,1) = top_value
       enddo
    enddo

    do k=2,npz+1
       if ( a2b_ord==4 ) then
         call a2b_ord4(pk(isd:,jsd:,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       else
         call a2b_ord2(pk(isd:,jsd:,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       endif
    enddo

    do k=1,npz+1
       if ( a2b_ord==4 ) then
         call a2b_ord4( gh(isd:,jsd:,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       else
         call a2b_ord2( gh(isd:,jsd:,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       endif
    enddo


    do k=1,npz

       if ( hydrostatic ) then
            do j=js,je+1
               do i=is,ie+1
                  wk(i,j) = pk(i,j,k+1) - pk(i,j,k)
               enddo
            enddo
       else
         if ( a2b_ord==4 ) then
            call a2b_ord4(delp(:,:,k), wk, npx, npy, is, ie, js, je, ng)
         else
            call a2b_ord2(delp(:,:,k), wk, npx, npy, is, ie, js, je, ng)
         endif
       endif

       do j=js,je+1
          do i=is,ie
             delu(i,j,k) = rdx(i,j) * dt_pg/(wk(i,j)+wk(i+1,j)) *  &
                         ((gh(i,j,k+1)-gh(i+1,j,k  ))*(pk(i+1,j,k+1)-pk(i  ,j,k)) &
                        + (gh(i,j,k  )-gh(i+1,j,k+1))*(pk(i  ,j,k+1)-pk(i+1,j,k)))
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             delv(i,j,k) = rdy(i,j) * dt_pg/(wk(i,j)+wk(i,j+1)) *  &
                         ((gh(i,j,k+1)-gh(i,j+1,k  ))*(pk(i,j+1,k+1)-pk(i,j  ,k)) &
                        + (gh(i,j,k  )-gh(i,j+1,k+1))*(pk(i,j  ,k+1)-pk(i,j+1,k)))
          enddo
       enddo
    enddo    ! end k-loop

 end subroutine grad1_p

 subroutine mix_dp(hydrostatic, w, delp, pt, km, ak, bk, CG)
  integer, intent(IN) :: km
  real, intent(IN) :: ak(km+1), bk(km+1)
  real  , intent(INOUT), dimension(isd:ied,jsd:jed,km):: pt, delp, w
  logical, intent(IN) :: hydrostatic, CG
! Local:
     real dp, dpmin
     integer i, j, k, ip
     integer ifirst, ilast
     integer jfirst, jlast


     if ( CG ) then
          ifirst = is-1; ilast = ie+1
          jfirst = js-1; jlast = je+1
     else
          ifirst = is; ilast = ie
          jfirst = js; jlast = je
     endif


!$omp parallel do default(shared) private(i, j, k, ip, dpmin, dp)
     do 1000 j=jfirst,jlast

          ip = 0

          do k=1, km-1
             dpmin = 0.005 * ( ak(k+1)-ak(k) + (bk(k+1)-bk(k))*1.E5 )
             do i=ifirst, ilast
              if(delp(i,j,k) < dpmin) then
! Remap from below and mix pt
                dp = dpmin - delp(i,j,k)
                pt(i,j,k) = (pt(i,j,k)*delp(i,j,k) + pt(i,j,k+1)*dp) / dpmin
                if ( .not.hydrostatic ) w(i,j,k) = (w(i,j,k)*delp(i,j,k) + w(i,j,k+1)*dp) / dpmin
                delp(i,j,k) = dpmin
                delp(i,j,k+1) = delp(i,j,k+1) - dp
                ip = ip + 1
              endif
            enddo
          enddo

! Bottom (k=km):
          dpmin = 0.005 * ( ak(km+1)-ak(km) + (bk(km+1)-bk(km))*1.E5 )
          do i=ifirst, ilast
            if(delp(i,j,km) < dpmin) then
! Remap from above and mix pt
              dp = dpmin - delp(i,j,km)
              pt(i,j,km) = (pt(i,j,km)*delp(i,j,km) + pt(i,j,km-1)*dp)/dpmin
              if ( .not.hydrostatic ) w(i,j,km) = (w(i,j,km)*delp(i,j,km) + w(i,j,km-1)*dp) / dpmin
              delp(i,j,km) = dpmin
              delp(i,j,km-1) = delp(i,j,km-1) - dp
              ip = ip + 1
            endif
          enddo
!      if ( fv_debug .and. ip/=0 ) write(*,*) 'Warning: Mix_dp', gid, j, ip 
       if ( ip/=0 ) write(*,*) 'Warning: Mix_dp', gid, j, ip 
1000   continue

 end subroutine  mix_dp



 subroutine geopk(ptop, pe, peln, delp, pk, gh, hs, pt, pkz, km, akap, CG)

     integer, intent(IN) :: km
     real             , intent(IN) :: akap
     real(p_precision), intent(IN) :: ptop
     real             , intent(IN) :: hs(isd:ied,jsd:jed)
     real             , intent(IN), dimension(isd:ied,jsd:jed,km):: pt, delp
     logical, intent(IN) :: CG
! !OUTPUT PARAMETERS
     real(pg_precision), intent(OUT), dimension(isd:ied,jsd:jed,km+1):: gh
     real(pg_precision), intent(OUT), dimension(isd:ied,jsd:jed,km+1):: pk
     real(p_precision) , intent(OUT)::   pe(is-1:ie+1,km+1,js-1:je+1)
     real(p_precision) , intent(out):: peln(is:ie,km+1,js:je)           ! ln(pe)
     real(p_precision) , intent(out):: pkz(is:ie,js:je,km)
! !DESCRIPTION:
!    Calculates geopotential and pressure to the kappa.
! Local:
     real(pg_precision)  p1d(is-2:ie+2)
     real(pg_precision) logp(is-2:ie+2)
     real(pg_precision) ptk
     integer i, j, k
     integer ifirst, ilast
     integer jfirst, jlast

     ptk  = ptop ** akap

     if ( .not. CG .and. a2b_ord==4 ) then   ! D-Grid
          ifirst = is-2; ilast = ie+2
          jfirst = js-2; jlast = je+2
     else
          ifirst = is-1; ilast = ie+1
          jfirst = js-1; jlast = je+1
     endif

!$omp parallel do default(shared) private(i, j, k, p1d, dp, logp)
     do 2000 j=jfirst,jlast

        do i=ifirst, ilast
           p1d(i) = ptop
           pk(i,j,1) = ptk
           gh(i,j,km+1) = hs(i,j)
        enddo

        if( j>(js-2) .and. j<(je+2) ) then
           do i=max(ifirst,is-1), min(ilast,ie+1) 
              pe(i,1,j) = ptop
           enddo
        endif

! Top down
        do k=2,km+1
          do i=ifirst, ilast
               p1d(i)  = p1d(i) + delp(i,j,k-1)
!            pk(i,j,k) = p1d(i) ** akap
! Optimized form:
              logp(i)  = log(p1d(i))
             pk(i,j,k) = exp( akap*logp(i) )
          enddo

          if( j>(js-2) .and. j<(je+2) ) then
             do i=max(ifirst,is-1), min(ilast,ie+1) 
                pe(i,k,j) = p1d(i)
             enddo
             if ( j>=js .and. j<=je ) then
                  do i=is,ie
                     peln(i,k,j) = logp(i)
                  enddo
             endif
          endif

        enddo

! Bottom up
        do k=km,1,-1
           do i=ifirst, ilast
              gh(i,j,k) = gh(i,j,k+1) + pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
           enddo
        enddo

2000  continue

      if ( .not. CG ) then
! This is for hydrostatic only
         do k=1,km
            do j=js,je
               do i=is,ie
                  pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
               enddo
            enddo
         enddo
      endif

 end subroutine geopk

 
 subroutine adv_pe(ua, va, pem, om, npx, npy, npz, ng)

 integer, intent(in) :: npx, npy, npz, ng
! Contra-variant wind components:
 real, intent(in), dimension(isd:ied,jsd:jed,npz):: ua, va
! Pressure at edges:
 real, intent(in) :: pem(is-1:ie+1,1:npz+1,js-1:je+1)
 real, intent(inout) :: om(isd:ied,jsd:jed,npz)

! Local:
 real, dimension(is:ie,js:je):: up, vp
 real v3(3,is:ie,js:je)

 real pin(isd:ied,jsd:jed)
 real  pb(isd:ied,jsd:jed)

 real grad(3,is:ie,js:je)
 real pdx(3,is:ie,js:je+1)
 real pdy(3,is:ie+1,js:je)
 integer :: i,j,k, n

!$omp parallel do default(shared) private(i, j, k, n, pdx, pdy, pin, pb, up, vp, grad, v3)
 do k=1,npz
    if ( k==npz ) then
       do j=js,je
          do i=is,ie
             up(i,j) = ua(i,j,npz)
             vp(i,j) = va(i,j,npz)
          enddo
       enddo
    else
       do j=js,je
          do i=is,ie
             up(i,j) = 0.5*(ua(i,j,k)+ua(i,j,k+1))
             vp(i,j) = 0.5*(va(i,j,k)+va(i,j,k+1))
          enddo
       enddo
    endif

! Compute Vect wind:
    do j=js,je
       do i=is,ie
          do n=1,3
             v3(n,i,j) = up(i,j)*ec1(n,i,j) + vp(i,j)*ec2(n,i,j) 
          enddo
       enddo
    enddo

    do j=js-1,je+1
       do i=is-1,ie+1
          pin(i,j) = pem(i,k+1,j)
       enddo
    enddo

! Compute pe at 4 cell corners:
    call a2b_ord2(pin, pb, npx, npy, is, ie, js, je, ng)


    do j=js,je+1
       do i=is,ie
          do n=1,3
             pdx(n,i,j) = (pb(i,j)+pb(i+1,j))*dx(i,j)*en1(n,i,j)
          enddo
       enddo
    enddo
    do j=js,je
       do i=is,ie+1
          do n=1,3
             pdy(n,i,j) = (pb(i,j)+pb(i,j+1))*dy(i,j)*en2(n,i,j)
          enddo
       enddo
    enddo

! Compute grad (pe) by Green's theorem
    do j=js,je
       do i=is,ie
          do n=1,3
             grad(n,i,j) = pdx(n,i,j+1) - pdx(n,i,j) - pdy(n,i,j) + pdy(n,i+1,j)
          enddo
       enddo
    enddo

! Compute inner product: V3 * grad (pe)
       do j=js,je
          do i=is,ie
             om(i,j,k) = om(i,j,k) + 0.5*rarea(i,j)*(v3(1,i,j)*grad(1,i,j) +   &
                         v3(2,i,j)*grad(2,i,j) + v3(3,i,j)*grad(3,i,j))
          enddo
       enddo
 enddo

 end subroutine adv_pe 

end module dyn_core_mod
