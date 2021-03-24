! SJL: Apr 12, 2012
! This revision may actually produce rounding level differences due to the elimination of KS to compute
! pressure level for remapping.
module fv_mapz_mod

!#ifndef MAPL_MODE
!  use constants_mod,     only: radius, pi, rvgas, rdgas, grav, hlv, hlf, cp_air
!#else
  use MAPL_Mod, only: MAPL_KAPPA, MAPL_RADIUS, MAPL_PI_R8, &
                      MAPL_RVAP, MAPL_RGAS, MAPL_GRAV
  use constants_mod,     only: hlv, hlf, cp_air
  use fv_arrays_mod,     only: REAL4, REAL8, FVPRC
!#endif
  use fv_arrays_mod,  only: REAL4, REAL8, FVPRC
  use tracer_manager_mod,only: get_tracer_index
  use field_manager_mod, only: MODEL_ATMOS
  use fv_grid_utils_mod, only: g_sum, ptop_min
!  use fv_fill_mod,       only: fillz
!  use fv_eta_mod ,       only: compute_dz_L32, hybrid_z_dz
  use mpp_domains_mod,   only: mpp_update_domains, domain2d
!  use mpp_mod,           only: FATAL, mpp_error, get_unit, mpp_root_pe, mpp_pe
  use fv_arrays_mod,     only: fv_grid_type
  use fv_timing_mod,     only: timing_on, timing_off
!#ifndef MAPL_MODE
!  use lin_cld_microphys_mod, only: sat_adj2
!#endif

  implicit none
!#ifdef MAPL_MODE
  real(REAL8), parameter :: RADIUS       = MAPL_RADIUS
  real(REAL8), parameter :: PI           = MAPL_PI_R8
  real(REAL8), parameter :: RVGAS        = MAPL_RVAP
  real(REAL8), parameter :: RDGAS        = MAPL_RGAS
  real(REAL8), parameter :: GRAV         = MAPL_GRAV
!#endif
  real(FVPRC), parameter:: t_min= 184.   ! below which applies stricter constraint
  real(FVPRC), parameter:: r3 = 1./3., r23 = 2./3., r12 = 1./12.
  real(REAL8), parameter :: r2=1./2., r0=0.0
  real(FVPRC), parameter:: cv_vap = 3.*rvgas  ! 1384.5
  real(FVPRC), parameter:: cv_air =  cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68
  real(FVPRC), parameter:: c_ice = 2106.           ! heat capacity of ice at 0.C
  real(FVPRC), parameter:: c_liq = 4190.
  real(kind=4) :: E_Flux, E_Flux_Nest
  private

  public compute_total_energy, Lagrangian_to_Eulerian,    &
         rst_remap, mappm, E_Flux, E_Flux_nest, map1_cubic, &
         map1_q2

!---- version number -----
!  character(len=128) :: version = '$Id: fv_mapz.F90,v 1.1.2.2.2.2.30.1.22.2.6.2.34.5.4.4.2.2.2.1.4.2.2.2 2017/02/21 20:43:38 bmauer Exp $'
!  character(len=128) :: tagname = '$Name: Heracles-UNSTABLE_ncepdyn_Feb222017 $'

contains

 subroutine Lagrangian_to_Eulerian(last_step, consv, ps, pe, delp, pkz, pk,   &
                      mdt, pdt, km, is,ie,js,je, isd,ied,jsd,jed,       &
                      nq, nwat, sphum, q_con, u, v, w, delz, pt, q, hs, r_vir, cp,  &
                      akap, cappa, kord_mt, kord_wz, kord_tr, kord_tm, kord_mt_pert, &
                      kord_wz_pert, kord_tr_pert, kord_tm_pert,  peln, te0_2d,        &
                      ng, ua, va, omga, te, ws, fill, reproduce_sum, out_dt, dtdt,      &
                      ptop, ak, bk, gridstruct, domain, ze0, gmao_cubic,remap_t, do_sat_adj, &
                      hydrostatic, hybrid_z, do_omega, do_adiabatic_init)
  logical, intent(in):: last_step
  real(FVPRC),    intent(in):: mdt                   ! remap time step
  real(FVPRC),    intent(in):: pdt                   ! phys time step
  integer, intent(in):: km
  integer, intent(in):: nq                    ! number of tracers (including h2o)
  integer, intent(in):: nwat
  integer, intent(in):: sphum                 ! index for water vapor (specific humidity)
  integer, intent(in):: ng
  integer, intent(in):: is,ie,isd,ied         ! starting & ending X-Dir index
  integer, intent(in):: js,je,jsd,jed         ! starting & ending Y-Dir index
  integer, intent(in):: kord_mt               ! Mapping order for the vector winds
  integer, intent(in):: kord_wz               ! Mapping order/option for w
  integer, intent(in):: kord_tr(nq)           ! Mapping order for tracers
  integer, intent(in):: kord_tm               ! Mapping order for thermodynamics
  integer, intent(in):: kord_mt_pert          ! Mapping order for the vector winds
  integer, intent(in):: kord_wz_pert          ! Mapping order/option for w
  integer, intent(in):: kord_tr_pert(nq)      ! Mapping order for tracers
  integer, intent(in):: kord_tm_pert          ! Mapping order for thermodynamics

  real(FVPRC), intent(in):: consv                 ! factor for TE conservation
  real(FVPRC), intent(in):: r_vir
  real(FVPRC), intent(in):: cp
  real(FVPRC), intent(in):: akap
  real(FVPRC), intent(in):: hs(isd:ied,jsd:jed)  ! surface geopotential
  real(FVPRC), intent(inout):: te0_2d(is:ie,js:je)
  real(FVPRC), intent(in):: ws(is:ie,js:je)

  logical, intent(in):: do_sat_adj
  logical, intent(in):: fill                  ! fill negative tracers
  logical, intent(in):: reproduce_sum
  logical, intent(in):: do_omega, do_adiabatic_init
  real(FVPRC), intent(in) :: ptop
  real(FVPRC), intent(in) :: ak(km+1)
  real(FVPRC), intent(in) :: bk(km+1)
  type(fv_grid_type), intent(IN), target :: gridstruct
  type(domain2d), intent(INOUT) :: domain

! !INPUT/OUTPUT
  real(FVPRC), intent(inout):: pk(is:ie,js:je,km+1) ! pe to the kappa
  real(FVPRC), intent(inout):: q(isd:ied,jsd:jed,km,nq)
  real(FVPRC), intent(inout):: delp(isd:ied,jsd:jed,km) ! pressure thickness
  real(FVPRC), intent(inout)::  pe(is-1:ie+1,km+1,js-1:je+1) ! pressure at layer edges
  real(FVPRC), intent(inout):: ps(isd:ied,jsd:jed)      ! surface pressure
  real(FVPRC), intent(inout):: ze0(is:ie,js:je,km+1)    ! Specified height at edges (m)

! u-wind will be ghosted one latitude to the north upon exit
  real(FVPRC), intent(inout)::  u(isd:ied  ,jsd:jed+1,km)   ! u-wind (m/s)
  real(FVPRC), intent(inout)::  v(isd:ied+1,jsd:jed  ,km)   ! v-wind (m/s)
  real(FVPRC), intent(inout)::  w(isd:ied,jsd:jed,km)   ! vertical velocity (m/s)
  real(FVPRC), intent(inout):: pt(isd:ied  ,jsd:jed  ,km)   ! cp*virtual potential temperature 
                                                     ! as input; output: temperature
  real(FVPRC), intent(inout), dimension(isd:ied,jsd:jed,km)::delz
  real(FVPRC), intent(inout), dimension(isd:isd,jsd:jsd,1)::cappa, q_con
  logical, intent(inout):: remap_t
  logical, intent(in):: gmao_cubic
  logical, intent(in):: hydrostatic
  logical, intent(in):: hybrid_z
  logical, intent(in):: out_dt

  real(FVPRC), intent(inout)::   ua(isd:ied,jsd:jed,km)   ! u-wind (m/s) on physics grid
  real(FVPRC), intent(inout)::   va(isd:ied,jsd:jed,km)   ! v-wind (m/s) on physics grid
  real(FVPRC), intent(inout):: omga(isd:ied,jsd:jed,km)   ! vertical press. velocity (pascal/sec)

  real(FVPRC), intent(inout)::   dtdt(is:ie,js:je,km)
  real(FVPRC), intent(out)::    pkz(is:ie,js:je,km)       ! layer-mean pk for converting t to pt
  real(FVPRC), intent(out)::     te(is:ie,js:je,km)
  real(FVPRC), intent(inout)::   peln(is:ie,km+1,js:je)     ! log(pe)

! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!-----------------------------------------------------------------------
  integer :: i,j,k 
! real(FVPRC) q_source(is:ie,js:je,nq)    ! numerical tracer source from surface
                                   ! in case fillz is not sufficient
  real(FVPRC), dimension(is:ie,js:je):: te_2d, zsum0, zsum1, dpeln
  real(FVPRC), dimension(is:ie,km)  :: q2, dp2, deng
  real(FVPRC), dimension(is:ie,km+1):: ze1, ze2, pe1, pe2, pk1, pk2, pn2, phis
     real(FVPRC)  pe0(is:ie+1,km+1)
     real(FVPRC)  pe3(is:ie+1,km+1)
     real(FVPRC)   gz(is:ie)
     real(FVPRC) dz1(km)
     real(FVPRC) rcp, rg, tmp, tpe, cvm, rgama, rrg, bkh, dtmp, dlnp, ztop, z_rat
     real(FVPRC) k1k, kapag, q_liq, q_sol
     integer liq_wat, ice_wat, rainwat, snowwat, cld_amt, graupel, iq, n, kp, k_next
     logical te_map
!     real(FVPRC), pointer, dimension(:,:) :: cosa_s, rsin2

  integer :: isp1, iep1, jsp1, jep1
  integer :: isdp1, iedp1, jsdp1, jedp1
  integer :: abs_kord_tm, abs_kord_tm_pert

  real(FVPRC) ::    u_tj(isd:ied  ,jsd:jed+1,km)
  real(FVPRC) ::    v_tj(isd:ied+1,jsd:jed  ,km)
  real(FVPRC) ::    w_tj(isd:ied,jsd:jed,km)
  real(FVPRC) ::   pt_tj(isd:ied  ,jsd:jed  ,km)
  real(FVPRC) ::   q2_tj(is:ie,km)
  real(FVPRC) :: deng_tj(is:ie,km)
  real(FVPRC) ::   te_tj(is:ie,js:je,km)
  real(FVPRC) ::    q_tj(isd:ied,jsd:jed,km,nq)
  real(FVPRC) :: delz_tj(isd:ied,jsd:jed,km)

  isp1 = is+1; iep1 = ie+1; jsp1 = js+1; jep1 = je+1
  isdp1 = isd+1; iedp1 = ied+1; jsdp1 = jsd+1; jedp1 = jed+1

!      cosa_s => gridstruct%cosa_s
!      rsin2 => gridstruct%rsin2

        k1k = akap / (1.-akap)    ! rg/Cv=0.4
      kapag = -akap / grav
         rg = akap * cp
      rgama = 1. - akap           ! cv/cp
        rcp = 1./ cp
        rrg = - rdgas/grav

      if ( kord_tm < 0 ) then
           te_map = .false.
      else
           te_map = .true.
      endif

      abs_kord_tm = abs(kord_tm)
      abs_kord_tm_pert = abs(kord_tm_pert)

      if ( nwat>=3 ) then
           liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
           ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
           cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
      endif
      if ( nwat.eq.6 ) then
           rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
           snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
           graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
      endif

!$OMP parallel do default(none) shared(is,ie,js,je,km,pe,ptop,kord_tm,remap_t,hydrostatic, &
!$OMP                                  pt,pk,rg,peln,q,liq_wat,rainwat,ice_wat,snowwat,    &
!$OMP                                  graupel,q_con,sphum,cappa,r_vir,rcp,k1k,kapag,delp, &
!$OMP                                  delz,akap,pkz,te,rsin2,u,v,cosa_s,hybrid_z,ztop,ps, &
!$OMP                                  ze0,ak,bk,nq,isd,ied,jsd,jed,kord_tr,fill,te_map,   &
!$OMP                                  hs,gz,w,ws,kord_wz,do_omega,omga,rrg,kord_mt,ua)    &
!$OMP                          private(cvm,dz1,q_liq,q_sol,z_rat,kp,k_next,bkh,deng,dp2,   &
!$OMP                                  pe0,pe1,pe2,pe3,pk1,pk2,pn2,phis,q2,ze1,ze2)
  do 1000 j=js,je+1

!$AD II-LOOP
      do k=1,km+1
        do i=is,ie
           pe1(i,k) = pe(i,k,j)
        enddo
     enddo

     do i=is,ie
        pe2(i,   1) = ptop
        pe2(i,km+1) = pe(i,km+1,j)
     enddo

  if ( j /= (je+1) ) then
       if ( kord_tm < 0 ) then
          if ( remap_t ) then
! Note: pt at this stage is cp*Theta_v
! Transform virtual pt to virtual Temp
             if ( hydrostatic ) then
             do k=1,km
                   do i=is,ie
                      pt(i,j,k) = pt(i,j,k) * (pk(i,j,k+1)-pk(i,j,k)) /  &
                                 (rg*(peln(i,k+1,j)-peln(i,k,j)) )
                   enddo
             enddo
             else
               do k=1,km
                  do i=is,ie
!#ifdef MOIST_CAPPA
!#ifdef USE_NWAT3
!                     q_liq = q(i,j,k,liq_wat)
!                     q_sol = q(i,j,k,ice_wat)
!#else
!                     q_liq = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
!                     q_sol = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
!#endif
!!                    cappa(i,j,k) = rdgas/(rdgas + ((1.-(q(i,j,k,sphum)+q_liq+q_sol))*cv_air + q(i,j,k,sphum)*cv_vap +  &
!!                                                        q_liq*c_liq + q_sol*c_ice)/(1.+r_vir*q(i,j,k,sphum)))
!                     q_con(i,j,k) = q_liq + q_sol
!                     cvm = (1.-(q(i,j,k,sphum)+q_Con(i,j,k)))*cv_air+q(i,j,k,sphum)*cv_vap+q_liq*c_liq+q_sol*c_ice
!                     cappa(i,j,k) = rdgas/(rdgas + cvm*(1.-q_con(i,j,k))/(1.+r_vir*q(i,j,k,sphum)-q_con(i,j,k)))
!                     pt(i,j,k) = rcp*pt(i,j,k)*exp(cappa(i,j,k)/(1.-cappa(i,j,k))*log(kapag*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
!#else
                     pt(i,j,k) = rcp*pt(i,j,k)*exp(k1k*log(kapag*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
!#endif
                  enddo
               enddo
             endif         ! hydro test
          endif            ! remap_t test
       else
           remap_t = .false.
           call pkez(km, is, ie, js, je, j, pe, pk, akap, peln, pkz, ptop)
! Compute cp*T + KE
           do k=1,km
                 do i=is,ie
                    te(i,j,k) = 0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                                 v(i,j,k)**2+v(i+1,j,k)**2 -  &
                               (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))  &
                              +  pt(i,j,k)*pkz(i,j,k)
                 enddo
           enddo
       endif

     if ( .not. hydrostatic ) then
        if ( hybrid_z ) then
           if ( km==32 ) then
                !call compute_dz_L32(km, ztop, dz1)
!          else
!               ztop = 18.E3
!               call hybrid_z_dz(km, dz1, ztop, 1.)
!               call mpp_error(FATAL,'==> Error from fv_mapz: hybrid_z works only with L32')
           endif
        else
           do k=1,km
                 do i=is,ie
                    delz(i,j,k) = -delz(i,j,k) / delp(i,j,k) ! ="specific volume"/grav
                 enddo
           enddo
        endif
      endif

! update ps
      do i=is,ie
         ps(i,j) = pe1(i,km+1)
      enddo

   if ( hybrid_z ) then
!--------------------------
! hybrid z_p coordinate
!--------------------------

        do i=is,ie
           ze1(i,km+1) = ze0(i,j,km+1)
        enddo

!!$AD II-LOOP
         do k=km,1,-1
           do i=is,ie
              ze1(i,k) = ze1(i,k+1) - delz(i,j,k)   ! current height
           enddo
        enddo
!
! Copy ztop; the top layer must be thick enough to prevent numerical problems.
!
        do i=is,ie
           ze2(i,  1) = ze1(i,1)
           ze0(i,j,1) = ze1(i,1)      ! Note: ze0 (top) updated
        enddo

!$AD II-LOOP
         do k=2,km+1
           do i=is,ie
              ze2(i,k) = ze0(i,j,k)   ! specified height
           enddo
        enddo

! Check/fix monotonicity of top 3 layers thickness
!-----------------------------------------------
!       if ( km==32 ) then
!        do i=is,ie
!           do k=1,3
!              z_rat = (ze2(i,1)-ze2(i,4)) / (dz1(1)+dz1(2)+dz1(3))
!              ze2(i,3) = ze2(i,4) + dz1(3)*z_rat
!              ze2(i,2) = ze2(i,3) + dz1(2)*z_rat
!           enddo
!        enddo
!       endif
!-----------------------------------------------

!$AD II-LOOP
         do k=1,km
           do i=is,ie
              deng(i,k) = -delp(i,j,k)/delz(i,j,k)  ! density * grav
           enddo
        enddo

        if (abs_kord_tm == abs_kord_tm_pert) then
           call remap_z_fb(km, ze1,       km, ze2, deng, is, ie, 1, abs_kord_tm)
        else
           deng_tj = deng
           !UNCOM      call remap_z(km, ze1,       km, ze2, deng, is, ie, 1, abs_kord_tm)
           call remap_z(km, ze1,       km, ze2, deng, is, ie, 1, abs_kord_tm_pert)
        endif

!-------------
! Update delz
!-------------
        do k=1,km
           do i=is,ie
              delz(i,j,k) = ze2(i,k+1) - ze2(i,k)
           enddo
        enddo

!------------
! update delp
!------------
        do k=1,km-1
           do i=is,ie
               dp2(i,k  ) = -delz(i,j,k) * deng(i,k)
               pe2(i,k+1) =     pe2(i,k) +  dp2(i,k)
           enddo
        enddo

        do i=is,ie
           dp2(i,km) = pe2(i,km+1) - pe2(i,km)  ! to reduce rounding error
        enddo
   else
!
! Hybrid sigma-P coordinate:
!
        do k=2,km
           do i=is,ie
              pe2(i,k) = ak(k) + bk(k)*pe(i,km+1,j)
           enddo
        enddo
!$AD II-LOOP
         do k=1,km
           do i=is,ie
              dp2(i,k) = pe2(i,k+1) - pe2(i,k)
           enddo
        enddo
   endif

!------------
! update delp
!------------
      do k=1,km
         do i=is,ie
            delp(i,j,k) = dp2(i,k)
         enddo
      enddo

!----------------
! Map constituents
!----------------
      if( nq > 5 ) then
         if (kord_tr(1) == kord_tr_pert(1)) then
           call mapn_tracer_fb(nq, km, km, pe1, pe2, q, dp2, kord_tr, j,     &
                            is, ie, isd, ied, jsd, jed, 0.0_FVPRC, fill)
         else
           q_tj(:,j,:,:) = q(:,j,:,:)
           !UNCOM      call mapn_tracer(nq, km, km, pe1, pe2, q, dp2, kord_tr, j,     &
           !UNCOM                       is, ie, isd, ied, jsd, jed, 0.0_FVPRC, fill)
           call mapn_tracer(nq, km, km, pe1, pe2, q, dp2, kord_tr_pert, j,     &
                            is, ie, isd, ied, jsd, jed, 0.0_FVPRC, fill)
         endif
      elseif ( nq > 0 ) then
! Remap one tracer at a time
          do iq=1,nq
             if (kord_tr(iq) == kord_tr_pert(iq)) then
                call map1_q2_fb(km, pe1, q(isd:ied,jsd:jed,1:km,iq),     &
                             km, pe2, q2, dp2,             &
                             is, ie, 0, kord_tr(iq), j, isd, ied, jsd, jed, 0.0_FVPRC)
             else
                q2_tj = q2
                !UNCOM        call map1_q2(km, pe1, q(isd:ied,jsd:jed,1:km,iq),     &
                !UNCOM                     km, pe2, q2, dp2,             &
                !UNCOM                     is, ie, 0, kord_tr(iq), j, isd, ied, jsd, jed, 0.0_FVPRC)
                call map1_q2(km, pe1, q(isd:ied,jsd:jed,1:km,iq),     &
                             km, pe2, q2, dp2,             &
                             is, ie, 0, kord_tr_pert(iq), j, isd, ied, jsd, jed, 0.0_FVPRC)
             endif
            !if (fill) call fillz(ie-is+1, km, 1, q2, dp2)
            do k=1,km
               do i=is,ie
                  q(i,j,k,iq) = q2(i,k)
               enddo
            enddo
          enddo
      endif

!------------------
! Compute p**Kappa
!------------------
!$AD II-LOOP
    do k=1,km+1
      do i=is,ie
         pk1(i,k) = pk(i,j,k)
      enddo
   enddo

   do i=is,ie
      pn2(i,   1) = peln(i,   1,j)
      pn2(i,km+1) = peln(i,km+1,j)
      pk2(i,   1) = pk1(i,   1)
      pk2(i,km+1) = pk1(i,km+1)
   enddo

   do k=2,km
      do i=is,ie
         pn2(i,k) = log(pe2(i,k))
         pk2(i,k) = exp(akap*pn2(i,k))
      enddo
   enddo

   if ( te_map ) then
!---------------------
! Compute Total Energy
!---------------------
        do i=is,ie
           phis(i,km+1) = hs(i,j)
        enddo
!!$AD II-LOOP
         do k=km,1,-1
           do i=is,ie
              phis(i,k) = phis(i,k+1) + pt(i,j,k)*(pk1(i,k+1)-pk1(i,k))
           enddo
        enddo
!$AD II-LOOP
         do k=1,km+1
           do i=is,ie
              phis(i,k) = phis(i,k) * pe1(i,k)
           enddo
        enddo
        do k=1,km
           do i=is,ie
              te(i,j,k) = te(i,j,k)+(phis(i,k+1)-phis(i,k))/(pe1(i,k+1)-pe1(i,k))
           enddo
        enddo
!----------------
! Map Total Energy
!----------------
        if (gmao_cubic) then
        call map1_cubic (km,   pe1,            &
                         km,   pe2,  te,       &
                         is, ie, j, is, ie, js, je, akap, T_VAR=1, conserv=.true.)
        else 
        if (kord_tm == kord_tm_pert) then
           call map1_ppm_fb (km,   pe1,      gz,   &
                          km,   pe2,  te,       &
                          is, ie, j, is, ie, js, je, 1, kord_tm)
        else
           te_tj(:,j,:) = te(:,j,:)
           !UNCOM      call map1_ppm (km,   pe1,      gz,   &
           !UNCOM                     km,   pe2,  te,       &
           !UNCOM                     is, ie, j, is, ie, js, je, 1, kord_tm)
           call map1_ppm (km,   pe1,      gz,   &
                          km,   pe2,  te,       &
                          is, ie, j, is, ie, js, je, 1, kord_tm_pert)
        endif
        endif
   else
     if ( remap_t ) then
!----------------------------------
! Map t using ze1 (hybrid_z) or logp 
!Note: the '1' in remap_z and map1_ppm was a '2' in siena_prerelease. It has been changed to be its value in siena -- 19mar12 lmh
!----------------------------------
       if ( hybrid_z ) then
!$AD II-LOOP
            do k=1,km
              do i=is,ie
                 deng(i,k) = pt(i,j,k)
              enddo
           enddo
           if (abs_kord_tm == abs_kord_tm_pert) then
              call remap_z_fb(km, ze1,       km, ze2, deng, is, ie, 1, abs_kord_tm)
           else
              deng_tj = deng
              !UNCOM        call remap_z(km, ze1,       km, ze2, deng, is, ie, 1, abs_kord_tm)
              call remap_z(km, ze1,       km, ze2, deng, is, ie, 1, abs_kord_tm_pert)
           endif
           do k=1,km
              do i=is,ie
                 pt(i,j,k) = deng(i,k)
              enddo
           enddo
       else
         if (abs_kord_tm == abs_kord_tm_pert) then
            call map_scalar_fb(km,  peln(is:ie,1:km+1,j),      gz,   &
                            km,  pn2,           pt,              &
                            is, ie, j, isd, ied, jsd, jed, 1, abs_kord_tm, t_min)
         else
            pt_tj(:,j,:) = pt(:,j,:)
            !UNCOM      call map_scalar(km,  peln(is:ie,1:km+1,j),      gz,   &
            !UNCOM                      km,  pn2,           pt,              &
            !UNCOM                      is, ie, j, isd, ied, jsd, jed, 1, abs_kord_tm, t_min)
            call map_scalar(km,  peln(is:ie,1:km+1,j),      gz,   &
                            km,  pn2,           pt,              &
                            is, ie, j, isd, ied, jsd, jed, 1, abs_kord_tm_pert, t_min)
         endif
       endif
     else
!----------------
! Map pt using pk
!----------------
       if (gmao_cubic) then
       call map1_cubic (km,   pe1,            &
                        km,   pe2,  pt,       &
                        is, ie, j, isd, ied, jsd, jed, akap, T_VAR=3, conserv=.false.)
       else
       if (abs_kord_tm == abs_kord_tm_pert) then
          call map1_ppm_fb (km,  pk1,       gz,       &
                         km,  pk2,  pt,                  &
                         is, ie, j, isd, ied, jsd, jed, 1, abs_kord_tm)
       else
          pt_tj(:,j,:) = pt(:,j,:)          
          !UNCOM    call map1_ppm (km,  pk1,       gz,       &
          !UNCOM                   km,  pk2,  pt,                  &
          !UNCOM                   is, ie, j, isd, ied, jsd, jed, 1, abs_kord_tm)
          call map1_ppm (km,  pk1,       gz,       &
                         km,  pk2,  pt,                  &
                         is, ie, j, isd, ied, jsd, jed, 1, abs_kord_tm_pert)
       endif
       end if
     endif
   endif

   if ( .not. hydrostatic ) then
! Remap vertical wind:
        if (kord_wz == kord_wz_pert) then
           call map1_ppm_fb (km,   pe1,      ws(is:ie,j),   &
                          km,   pe2,  w,              &
                          is, ie, j, isd, ied, jsd, jed, -2, kord_wz)
        else
           w_tj(:,j,:) = w(:,j,:)
           !UNCOM      call map1_ppm (km,   pe1,      ws(is:ie,j),   &
           !UNCOM                     km,   pe2,  w,              &
           !UNCOM                     is, ie, j, isd, ied, jsd, jed, -2, kord_wz)
           call map1_ppm (km,   pe1,      ws(is:ie,j),   &
                          km,   pe2,  w,              &
                          is, ie, j, isd, ied, jsd, jed, -2, kord_wz_pert)
        endif
     if ( .not. hybrid_z ) then
! Remap delz for hybrid sigma-p coordinate
        if (abs_kord_tm == abs_kord_tm_pert) then
           call map1_ppm_fb (km,   pe1,        gz,   &
                          km,   pe2, delz,              &
                          is, ie, j, isd,  ied,  jsd,  jed,  1, abs_kord_tm)
        else
           delz_tj(:,j,:) = delz(:,j,:)
           !UNCOM        call map1_ppm (km,   pe1,        gz,   &
           !UNCOM                       km,   pe2, delz,              &
           !UNCOM                       is, ie, j, isd,  ied,  jsd,  jed,  1, abs_kord_tm)
           call map1_ppm (km,   pe1,        gz,   &
                          km,   pe2, delz,              &
                          is, ie, j, isd,  ied,  jsd,  jed,  1, abs_kord_tm_pert)
        endif
        do k=1,km
           do i=is,ie
              delz(i,j,k) = -delz(i,j,k)*dp2(i,k)
           enddo
        enddo
     endif
   endif

!----------
! Update pk
!----------
   do k=2,km
      do i=is,ie
         pk(i,j,k) = pk2(i,k)
      enddo
   enddo

!----------------
   if ( do_omega ) then
      do i=is,ie
         pe3(i,1) = 0.
      enddo
!$AD II-LOOP
       do k=2,km+1
         do i=is,ie
            pe3(i,k) = omga(i,j,k-1)
         enddo
      enddo
   endif

   do k=1,km+1
      do i=is,ie
          pe0(i,k)   = peln(i,k,j)
         peln(i,k,j) =  pn2(i,k)
      enddo
   enddo

!------------
! Compute pkz
!------------
   if ( hydrostatic ) then
      do k=1,km
         do i=is,ie
            pkz(i,j,k) = (pk2(i,k+1)-pk2(i,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
         enddo
      enddo
   else
      if ( remap_t ) then
! Note: pt at this stage is T_v or T_m
         do k=1,km
         do i=is,ie
!#ifdef MOIST_CAPPA
!#ifdef USE_NWAT3
!            q_liq = q(i,j,k,liq_wat)
!            q_sol = q(i,j,k,ice_wat)
!#else
!            q_liq = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
!            q_sol = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
!#endif
!            q_con(i,j,k) = q_liq + q_sol
!!           cappa(i,j,k) = rdgas/(rdgas+((1.-q(i,j,k,sphum)-q_con(i,j,k))*cv_air + q(i,j,k,sphum)*cv_vap +  &
!!                                            q_liq*c_liq + q_sol*c_ice)/(1.+r_vir*q(i,j,k,sphum)))
!            cvm = (1.-(q(i,j,k,sphum)+q_con(i,j,k)))*cv_air+q(i,j,k,sphum)*cv_vap+q_liq*c_liq+q_sol*c_ice
!            cappa(i,j,k) = rdgas/(rdgas + cvm*(1.-q_con(i,j,k))/(1.+r_vir*q(i,j,k,sphum)-q_con(i,j,k)))
!            pkz(i,j,k) = exp(cappa(i,j,k)*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
!#else
            pkz(i,j,k) = exp(akap*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
!#endif
         enddo
         enddo
      else
! Note: pt at this stage is cp*Theta_v
         do k=1,km
         do i=is,ie
!#ifdef MOIST_CAPPA
!            pkz(i,j,k) = exp(cappa(i,j,k)/(1.-cappa(i,j,k))*log(kapag*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)) )
!#else
            pkz(i,j,k) = exp( k1k*log(kapag*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)) )
!#endif
         enddo
         enddo
      endif
   endif

! Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
   if ( do_omega ) then
!$AD II-LOOP
    do k=1,km
      do i=is,ie
         dp2(i,k) = 0.5*(peln(i,k,j) + peln(i,k+1,j))
      enddo
   enddo
   do i=is,ie
       k_next = 1
       do n=1,km
          kp = k_next
          do k=kp,km
             if( dp2(i,n) <= pe0(i,k+1) .and. dp2(i,n) >= pe0(i,k) ) then
                 omga(i,j,n) = pe3(i,k)  +  (pe3(i,k+1) - pe3(i,k)) *    &
                       (dp2(i,n)-pe0(i,k)) / (pe0(i,k+1)-pe0(i,k) )
                 k_next = k
                 exit
             endif
          enddo
       enddo
   enddo
   endif     ! end do_omega

  endif !(j < je+1)

 if ( .not.hybrid_z ) then
      do i=is,ie+1
         pe0(i,1) = pe(i,1,j)
      enddo
!------
! map u
!------
!$AD II-LOOP
       do k=2,km+1
         do i=is,ie
            pe0(i,k) = 0.5*(pe(i,k,j-1)+pe1(i,k))
         enddo
      enddo

      do k=1,km+1
         bkh = 0.5*bk(k)
         do i=is,ie
            pe3(i,k) = ak(k) + bkh*(pe(i,km+1,j-1)+pe1(i,km+1))
         enddo
      enddo

      if (kord_mt == kord_mt_pert) then
         call map1_ppm_fb( km, pe0(is:ie,1:km+1),        gz,   &
                        km, pe3(is:ie,1:km+1),   u,               &
                        is, ie, j, isd, ied, jsd, jedp1, -1, kord_mt)
      else
         u_tj(:,j,:) = u(:,j,:)
         !UNCOM    call map1_ppm( km, pe0(is:ie,1:km+1),        gz,   &
         !UNCOM                   km, pe3(is:ie,1:km+1),   u,               &
         !UNCOM                   is, ie, j, isd, ied, jsd, jedp1, -1, kord_mt)
         call map1_ppm( km, pe0(is:ie,1:km+1),        gz,   &
                        km, pe3(is:ie,1:km+1),   u,               &
                        is, ie, j, isd, ied, jsd, jedp1, -1, kord_mt_pert)
      endif

   if (j < je+1) then
!------
! map v
!------
       do i=is,ie+1
          pe3(i,1) = ak(1)
       enddo

       do k=2,km+1
          bkh = 0.5*bk(k)
          do i=is,ie+1
             pe0(i,k) =         0.5*(pe(i-1,k,   j)+pe(i,k,   j))
             pe3(i,k) = ak(k) + bkh*(pe(i-1,km+1,j)+pe(i,km+1,j))
          enddo
       enddo

       if (kord_mt == kord_mt_pert) then
          call map1_ppm_fb (km, pe0,     gz,    &
                         km, pe3,  v, is, iep1,    &
                         j, isd, iedp1, jsd, jed, -1, kord_mt)
       else
          v_tj(:,j,:) = v(:,j,:)
          !UNCOM      call map1_ppm (km, pe0,     gz,    &
          !UNCOM                     km, pe3,  v, is, iep1,    &
          !UNCOM                     j, isd, iedp1, jsd, jed, -1, kord_mt)
          call map1_ppm (km, pe0,     gz,    &
                         km, pe3,  v, is, iep1,    &
                         j, isd, iedp1, jsd, jed, -1, kord_mt_pert)
       endif
   endif ! (j < je+1)
 endif    ! end hybrid_z check

     do k=1,km
        do i=is,ie
           ua(i,j,k) = pe2(i,k+1)
        enddo
     enddo

1000  continue

if ( hybrid_z ) then   
!------- Hybrid_z section ---------------
   if ( gridstruct%square_domain ) then
     call mpp_update_domains(ua , domain,  whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
   else
     call mpp_update_domains(ua , domain, complete=.true.)
   endif

! u-wind


!$OMP parallel do default(none) shared(is,ie,js,je,km,ptop,u,gz,pe,ua,isd,ied,jsd,jed,kord_mt) &
!$OMP                          private(pe1, pe2)
   do j=js,je+1
      do i=is,ie
         pe1(i,1) = ptop
         pe2(i,1) = ptop
      enddo
      do k=2,km+1
         do i=is,ie
            pe1(i,k) = 0.5*(pe(i,k,  j-1) + pe(i,k,j  ))
            pe2(i,k) = 0.5*(ua(i,j-1,k-1) + ua(i,j,k-1))
         enddo
      enddo

      if (kord_mt == kord_mt_pert) then
         call map1_ppm_fb( km, pe1,       gz, &
                        km, pe2,   u,            &
                        is, ie, j, isd, ied, jsd, jedp1, -1, kord_mt)
      else
         u_tj(:,j,:) = u(:,j,:)
         !UNCOM    call map1_ppm( km, pe1,       gz, &
         !UNCOM                   km, pe2,   u,            &
         !UNCOM                   is, ie, j, isd, ied, jsd, jedp1, -1, kord_mt)
         call map1_ppm( km, pe1,       gz, &
                        km, pe2,   u,            &
                        is, ie, j, isd, ied, jsd, jedp1, -1, kord_mt_pert)
      endif
   enddo

! v-wind
!$OMP parallel do default(none) shared(is,ie,js,je,ptop,pe,ua,km,v,gz,isd,ied,jsd,jed,kord_mt) &
!$OMP                          private(pe0, pe3)
   do j=js,je
      do i=is,ie+1
         pe0(i,1) = ptop
         pe3(i,1) = ptop
      enddo

      do k=2,km+1
         do i=is,ie+1
            pe0(i,k) = 0.5*(pe(i-1,k,j  ) + pe(i,k,j  ))
            pe3(i,k) = 0.5*(ua(i-1,j,k-1) + ua(i,j,k-1))
         enddo
      enddo

      if (kord_mt == kord_mt_pert) then
         call map1_ppm_fb (km, pe0,      gz,   &
                        km, pe3,  v, is, iep1,    &
                        j, isd, iedp1, jsd, jed, -1, kord_mt)
      else
         v_tj(:,j,:) = v(:,j,:)
         !UNCOM    call map1_ppm (km, pe0,      gz,   &
         !UNCOM                   km, pe3,  v, is, iep1,    &
         !UNCOM                   j, isd, iedp1, jsd, jed, -1, kord_mt)
         call map1_ppm (km, pe0,      gz,   &
                        km, pe3,  v, is, iep1,    &
                        j, isd, iedp1, jsd, jed, -1, kord_mt_pert)
      endif
   enddo
endif 
!------------- Hybrid_z section ----------------------

!$OMP parallel do default(none) shared(is,ie,js,je,km,pe,ua)
!$AD II-LOOP
   do k=2,km
     do j=js,je
        do i=is,ie
           pe(i,k,j) = ua(i,j,k-1)
        enddo
     enddo
  enddo

!  call cubed_to_latlon(u,  v, ua, va, dx, dy, rdxa, rdya, km, 1, flagstruct%c2l_ord)

  if( last_step .and. consv > 0.  .and. (.not.do_adiabatic_init)  ) then

    if ( te_map ) then
!$OMP parallel do default(none) shared(is,ie,js,je,km,te_2d,te,delp)
      do j=js,je
          do i=is,ie
             te_2d(i,j) = te(i,j,1)*delp(i,j,1)
          enddo
          do k=2,km
             do i=is,ie
                te_2d(i,j) = te_2d(i,j) + te(i,j,k)*delp(i,j,k)
             enddo
          enddo
      enddo
    else
!$OMP parallel do default(none) shared(is,ie,js,je,km,remap_t,hydrostatic,hs,rg,pt,peln, &
!$OMP                                  te_2d,pe,delp,cp,rsin2,u,v,cosa_s,delz,rainwat,   &
!$OMP                                  liq_wat,ice_wat,snowwat,graupel,q_con,r_vir,      &
!$OMP                                  sphum,w,pk,pkz,rgama)                             &
!$OMP                          private(q_liq, q_sol, cvm, gz, phis)
      do j=js,je
      if ( remap_t ) then
         if ( hydrostatic ) then
         do i=is,ie
            gz(i) = hs(i,j)
            do k=1,km
               gz(i) = gz(i) + rg*pt(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))
            enddo
         enddo
         do i=is,ie
            te_2d(i,j) = pe(i,km+1,j)*hs(i,j) - pe(i,1,j)*gz(i)
         enddo

         do k=1,km
            do i=is,ie
               te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cp*pt(i,j,k) +   &
                            0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                             v(i,j,k)**2+v(i+1,j,k)**2 -  &
                           (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j)))
            enddo
         enddo
         else
! non-hydrostatic & remap_t
           do i=is,ie
              phis(i,km+1) = hs(i,j)
              do k=km,1,-1
                 phis(i,k) = phis(i,k+1) - grav*delz(i,j,k)
              enddo
           enddo
           do i=is,ie
              te_2d(i,j) = 0.
           enddo
           do k=1,km
              do i=is,ie
! KE using 3D winds:
!#ifdef USE_COND
!#ifdef USE_NWAT3
!           q_liq = q(i,j,k,liq_wat)
!           q_sol = q(i,j,k,ice_wat)
!!#else
!           q_liq = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
!           q_sol = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
!#endif
!           q_con(i,j,k) = q_liq + q_sol
!#ifdef MOIST_CAPPA
!           cvm = (1.-(q(i,j,k,sphum)+q_con(i,j,k)))*cv_air+q(i,j,k,sphum)*cv_vap+q_liq*c_liq+q_sol*c_ice
!#else
!           cvm = cv_air
!#endif
!                 te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cvm*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum)-q_con(i,j,k)) + &
!#else
                 te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cv_air*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum)) + &
!#endif
                             0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*gridstruct%rsin2(i,j)*( &
                              u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                             (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))))
              enddo
           enddo
         endif
      else
         if ( hydrostatic ) then
            do i=is,ie
               gz(i) = hs(i,j)
               do k=1,km
                  gz(i) = gz(i) + pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
               enddo
            enddo

            do i=is,ie
               te_2d(i,j) = pe(i,km+1,j)*hs(i,j) - pe(i,1,j)*gz(i)
            enddo
            do k=1,km
               do i=is,ie
                  te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(pt(i,j,k)*pkz(i,j,k) +   &
                               0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                                v(i,j,k)**2+v(i+1,j,k)**2 -  &
                            (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j)))
               enddo
            enddo
         else
!-----------------
! Non-hydrostatic:
!-----------------
           do i=is,ie
              phis(i,km+1) = hs(i,j)
              do k=km,1,-1
                 phis(i,k) = phis(i,k+1) - grav*delz(i,j,k)
              enddo
           enddo
           do i=is,ie
              te_2d(i,j) = 0.
           enddo
           do k=1,km
              do i=is,ie
! KE using 3D winds:
!#ifdef USE_COND
!#ifdef USE_NWAT3
!                 q_con(i,j,k) = q(i,j,k,liq_wat) + q(i,j,k,ice_wat)
!#else
!                 q_con(i,j,k) = q(i,j,k,liq_wat) + q(i,j,k,rainwat) + q(i,j,k,ice_wat)  &
!                              + q(i,j,k,snowwat) + q(i,j,k,graupel)
!#endif
!                 te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(rgama*pt(i,j,k)*pkz(i,j,k)/(1.+r_vir*q(i,j,k,sphum)-q_con(i,j,k)) + &
!#else
                 te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(rgama*pt(i,j,k)*pkz(i,j,k)/(1.+r_vir*q(i,j,k,sphum)) + &
!#endif
                              0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*gridstruct%rsin2(i,j)*( &
                              u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                             (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))))
              enddo
           enddo
         endif
      endif
      enddo   ! openMP j-loop
    endif

!$OMP parallel do default(none) shared(is,ie,js,je,zsum1,pkz,delp,km,zsum0,ptop,pk,te_2d,te0_2d)
      do j=js,je
         do i=is,ie
            zsum1(i,j) = pkz(i,j,1)*delp(i,j,1)
         enddo
         do k=2,km
            do i=is,ie
               zsum1(i,j) = zsum1(i,j) + pkz(i,j,k)*delp(i,j,k)
            enddo
         enddo

         do i=is,ie
            zsum0(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1)) + zsum1(i,j)
            te_2d(i,j) = te0_2d(i,j) - te_2d(i,j)
         enddo
      enddo

         tpe = consv*g_sum(domain, te_2d, is, ie, js, je, ng, gridstruct%area_64, 0)
      E_Flux = tpe / (grav*pdt*4.*pi*radius**2)    ! unit: W/m**2
                                                   ! Note pdt is "phys" time step

      if ( hydrostatic ) then
           dtmp = tpe / (cp*g_sum(domain, zsum0,  is, ie, js, je, ng, gridstruct%area_64, 0))
      else
           dtmp = tpe / (cv_air*g_sum(domain, zsum1,  is, ie, js, je, ng, gridstruct%area_64, 0))
      endif
!-------------------------------------------------------------------------------
! One may use this quick fix to ensure reproducibility at the expense of a lower
! floating precision; this is fine for the TE correction
!-------------------------------------------------------------------------------
      if ( reproduce_sum ) dtmp = real(dtmp, 4) ! convert to 4-byte real
  else
      dtmp   = 0.
      E_Flux = 0.
!#ifdef USE_COND
!!$OMP parallel do default(none) shared(is,ie,js,je,km,q_con,q,liq_wat,ice_wat, &
!!$OMP                                  rainwat,snowwat,graupel)
!      do k=1,km
!         do j=js,je
!         do i=is,ie
!#ifdef USE_NWAT3
!            q_con(i,j,k) = q(i,j,k,liq_wat)+q(i,j,k,ice_wat)
!#else
!            q_con(i,j,k) = q(i,j,k,liq_wat)+q(i,j,k,rainwat)+q(i,j,k,ice_wat)+q(i,j,k,snowwat)+q(i,j,k,graupel)
!#endif
!         enddo
!         enddo
!      enddo
!#endif
  endif        ! end consv check

  if ( te_map ) then !kord_mt >= 0
!$OMP parallel do default(none) shared(is,ie,js,je,km,hs,te,rsin2,u,v,cosa_s,rg,peln, &
!$OMP                                  cp,pe,delp,pt,pkz,dtmp) &
!$OMP                          private(gz, tpe, tmp, dlnp)
      do j=js,je
         do i=is,ie
            gz(i) = hs(i,j)
         enddo
         do k=km,1,-1
            do i=is,ie
               tpe = te(i,j,k) - gz(i) - 0.25*gridstruct%rsin2(i,j)*(    &
                     u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                    (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j) )
               dlnp = rg*(peln(i,k+1,j) - peln(i,k,j))
               tmp = tpe / (cp - pe(i,k,j)*dlnp/delp(i,j,k))
               pt(i,j,k) = cp*(tmp/pkz(i,j,k) + dtmp)
               gz(i) = gz(i) + dlnp*tmp
            enddo
         enddo           ! end k-loop
      enddo
  else
    if ( remap_t ) then
!#ifndef MAPL_MODE
!        ! Note: pt at this stage is T_v
!        if ( .not. do_adiabatic_init ) then
!#ifdef FAST_LIN_MP
!! Notes: must set consv=0
!                                           call timing_on('FAST_LIN_MP')
!      call lin_fast_mp(u, v, w, pt, q, sphum, liq_wat, rainwat, ice_wat, snowwat, graupel, cld_amt,  & 
!                       delp, delz, gridstruct, mdt, lprec, fprec, f_land,      &
!                       hydrostatic, phys_hydrostatic, time)
!                                           call timing_off('FAST_LIN_MP')
!#else
!          if ( do_sat_adj ) then
!                                           call timing_on('sat_adj2')
!!$OMP parallel do default(none) shared(is,ie,js,je,km,peln,ng,mdt,hydrostatic,consv,te, &
!!$OMP                                  q,isd,jsd,sphum,liq_wat,ice_wat,rainwat,snowwat, &
!!$OMP                                  graupel,cld_amt,gridstruct,delz,pt,delp,         &
!!$OMP                                  q_con,cappa,dtdt,out_dt,last_step,pkz,rrg,akap ) &
!!$OMP                          private(dpeln)
!           do k=1,km
!              do j=js,je
!                 do i=is,ie
!                    dpeln(i,j) = peln(i,k+1,j) - peln(i,k,j)
!                 enddo
!              enddo
!              call sat_adj2(mdt, is, ie, js, je, ng, hydrostatic, consv>0., &
!                            te(is,js,k), q(isd,jsd,k,sphum), q(isd,jsd,k,liq_wat),   &
!                            q(isd,jsd,k,ice_wat), q(isd,jsd,k,rainwat),    &
!                            q(isd,jsd,k,snowwat), q(isd,jsd,k,graupel), q(isd,jsd,k,cld_amt), gridstruct%area(isd,jsd), &
!                            dpeln, delz(isd:,jsd:,k), pt(isd,jsd,k), delp(isd,jsd,k),           &
!                            q_con(isd:,jsd:,k), cappa(isd:,jsd:,k), dtdt(is,js,k), out_dt, last_step)
!              if ( .not. hydrostatic  ) then
!                 do j=js,je
!                    do i=is,ie
!#ifdef MOIST_CAPPA
!                       pkz(i,j,k) = exp(cappa(i,j,k)*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
!#else
!                       pkz(i,j,k) = exp(akap*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
!#endif
!                    enddo
!                 enddo
!              endif
!           enddo    ! OpenMP k-loop
!           if ( consv > 0. ) then
!!$OMP parallel do default(none) shared(is,ie,js,je,km,te_2d,te,te0_2d)
!                do j=js,je
!                   do k=1,km
!                      do i=is,ie
!                         te0_2d(i,j) = te0_2d(i,j) + te(i,j,k)
!                      enddo
!                   enddo
!                enddo
!           endif
!                                           call timing_off('sat_adj2')
!         endif
!#endif
!        endif
!#endif
        if ( last_step ) then
         if ( hydrostatic ) then
!$OMP parallel do default(none) shared(is,ie,js,je,km,pt,dtmp,pkz,r_vir,q,sphum)
!$AD II-LOOP
            do k=1,km
              do j=js,je
                 do i=is,ie
                    pt(i,j,k) = (pt(i,j,k)+dtmp*pkz(i,j,k))/(1.+r_vir*q(i,j,k,sphum))
                 enddo
              enddo
           enddo
         else
!$OMP parallel do default(none) shared(is,ie,js,je,km,q,liq_wat,ice_wat,rainwat,snowwat, &
!$OMP                                  graupel,q_con,pt,dtmp,pkz,r_vir,sphum )           &
!$OMP                          private(q_liq, q_sol, cvm)
!$AD II-LOOP
            do k=1,km
              do j=js,je
                 do i=is,ie
                    ! Output temperature if last_step
!#ifdef USE_COND
!#ifdef MOIST_CAPPA
!#ifdef USE_NWAT3
!           q_liq = q(i,j,k,liq_wat)
!           q_sol = q(i,j,k,ice_wat)
!#else
!           q_liq = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
!           q_sol = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
!#endif
!           cvm = (1.-(q(i,j,k,sphum)+q_con(i,j,k)))*cv_air+q(i,j,k,sphum)*cv_vap+q_liq*c_liq+q_sol*c_ice
!
!                    pt(i,j,k) = (pt(i,j,k)+dtmp*cv_air/cvm*pkz(i,j,k))/(1.+r_vir*q(i,j,k,sphum)-q_con(i,j,k))
!#else
!                    pt(i,j,k) = (pt(i,j,k)+dtmp*pkz(i,j,k))/(1.+r_vir*q(i,j,k,sphum)-q_con(i,j,k))
!#endif
!#else
                    pt(i,j,k) = (pt(i,j,k)+dtmp*pkz(i,j,k))/(1.+r_vir*q(i,j,k,sphum))
!#endif
                 enddo
              enddo
           enddo
         endif
        else
!$OMP parallel do default(none) shared(is,ie,js,je,km,pt,cp,pkz,dtmp,q,liq_wat,ice_wat, &
!$OMP                                  rainwat,snowwat,graupel,sphum,q_con) &
!$OMP                          private(q_liq, q_sol, cvm)
!$AD II-LOOP
            do k=1,km
              do j=js,je
                 do i=is,ie
!#ifdef MOIST_CAPPA
!#ifdef USE_NWAT3
!           q_liq = q(i,j,k,liq_wat)
!           q_sol = q(i,j,k,ice_wat)
!#else
!           q_liq = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
!           q_sol = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
!#endif
!           cvm = (1.-(q(i,j,k,sphum)+q_con(i,j,k)))*cv_air+q(i,j,k,sphum)*cv_vap+q_liq*c_liq+q_sol*c_ice
!                    pt(i,j,k) = cp*(pt(i,j,k)/pkz(i,j,k) + dtmp*cv_air/cvm)
!#else
                    pt(i,j,k) = cp*(pt(i,j,k)/pkz(i,j,k) + dtmp)
!#endif
                 enddo
              enddo
           enddo
        endif
    else
!$OMP parallel do default(none) shared(is,ie,js,je,km,pt,cp,dtmp)
!$AD II-LOOP
         do k=1,km
           do j=js,je
              do i=is,ie
                 pt(i,j,k) = pt(i,j,k) + cp*dtmp
              enddo
           enddo
        enddo
     endif
  endif

! call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km, 1, flagstruct%c2l_ord)

!  nullify(cosa_s)
!  nullify(rsin2)

 end subroutine Lagrangian_to_Eulerian


 subroutine compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, km,       &
                                 u, v, w, delz, pt, delp, q, qc, pe, peln, hs, &
                                 rsin2_l, cosa_s_l, &
                                 r_vir,  cp, rg, hlv, te_2d, ua, va, teq, &
                                 moist_phys, sphum, liq_wat, rainwat, ice_wat, snowwat, graupel, hydrostatic, id_te)
!------------------------------------------------------
! Compute vertically integrated total energy per column
!------------------------------------------------------
! !INPUT PARAMETERS:
   integer,  intent(in):: km, is, ie, js, je, isd, ied, jsd, jed, id_te
   integer,  intent(in):: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel
   real(FVPRC), intent(inout), dimension(isd:ied,jsd:jed,km):: ua, va
   real(FVPRC), intent(in), dimension(isd:ied,jsd:jed,km):: pt, delp
   real(FVPRC), intent(in), dimension(isd:ied,jsd:jed,km,*):: q
   real(FVPRC), intent(in), dimension(is:ie,js:je,km):: qc
   real(FVPRC), intent(inout)::  u(isd:ied,  jsd:jed+1,km)
   real(FVPRC), intent(inout)::  v(isd:ied+1,jsd:jed,  km)
   real(FVPRC), intent(in)::  w(isd:ied,jsd:jed,km)   ! vertical velocity (m/s)
   real(FVPRC), intent(in):: delz(isd:ied,jsd:jed,km)
   real(FVPRC), intent(in):: hs(isd:ied,jsd:jed)  ! surface geopotential
   real(FVPRC), intent(in)::   pe(is-1:ie+1,km+1,js-1:je+1) ! pressure at layer edges
   real(FVPRC), intent(in):: peln(is:ie,km+1,js:je)  ! log(pe)
   real(FVPRC), intent(in):: cp, rg, r_vir, hlv
   real(FVPRC), intent(in) :: rsin2_l(isd:ied, jsd:jed)
   real(FVPRC), intent(in) :: cosa_s_l(isd:ied, jsd:jed)
   logical, intent(in):: moist_phys, hydrostatic
! Output:
   real(FVPRC), intent(out):: te_2d(is:ie,js:je)   ! vertically integrated TE
   real(FVPRC), intent(out)::   teq(is:ie,js:je)   ! Moist TE
! Local
   real(FVPRC), dimension(is:ie,km):: tv
   real(FVPRC)  phiz(is:ie,km+1)
   real(FVPRC) cvm, q_liq, q_sol
   integer i, j, k

!----------------------
! Output lat-lon winds:
!----------------------
!  call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km, flagstruct%c2l_ord)

!$OMP parallel do default(none) shared(is,ie,js,je,km,hydrostatic,hs,pt,qc,rg,peln,te_2d, &
!$OMP                                  q,pe,delp,cp,rsin2_l,u,v,cosa_s_l,delz,moist_phys,w, &
!$OMP                                  liq_wat,rainwat,ice_wat,snowwat,graupel,sphum)   &
!$OMP                          private(phiz, tv, cvm, q_sol, q_liq)
  do j=js,je

     if ( hydrostatic ) then

        do i=is,ie
           phiz(i,km+1) = hs(i,j)
        enddo
        do k=km,1,-1
           do i=is,ie
                tv(i,k) = pt(i,j,k)*(1.+qc(i,j,k))
              phiz(i,k) = phiz(i,k+1) + rg*tv(i,k)*(peln(i,k+1,j)-peln(i,k,j))
           enddo
        enddo

        do i=is,ie
           te_2d(i,j) = pe(i,km+1,j)*phiz(i,km+1) - pe(i,1,j)*phiz(i,1)
        enddo

        do k=1,km
           do i=is,ie
              te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cp*tv(i,k) +            &
                           0.25*rsin2_l(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +      &
                                            v(i,j,k)**2+v(i+1,j,k)**2 -      &
                       (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s_l(i,j)))
           enddo
        enddo

     else
!-----------------
! Non-hydrostatic:
!-----------------
     do i=is,ie
        phiz(i,km+1) = hs(i,j)
        do k=km,1,-1
           phiz(i,k) = phiz(i,k+1) - grav*delz(i,j,k)
        enddo
     enddo
     do i=is,ie
        te_2d(i,j) = 0.
     enddo
     if ( moist_phys ) then
     do k=1,km
        do i=is,ie
!#if defined(USE_COND) && defined(MOIST_CAPPA)
!#ifdef USE_NWAT3
!           q_liq = q(i,j,k,liq_wat)
!           q_sol = q(i,j,k,ice_wat)
!#else
!           q_liq = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
!           q_sol = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
!#endif
!           cvm = (1.-(q(i,j,k,sphum)+q_liq+q_sol))*cv_air+q(i,j,k,sphum)*cv_vap+q_liq*c_liq+q_sol*c_ice
!           te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cvm*pt(i,j,k) +  &
!#else
           te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv_air*pt(i,j,k) +  &
!#endif
                        0.5*(phiz(i,k)+phiz(i,k+1)+w(i,j,k)**2+0.5*rsin2_l(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                        v(i,j,k)**2+v(i+1,j,k)**2-(u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s_l(i,j))))
        enddo
     enddo
     else
       do k=1,km
          do i=is,ie
             te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv_air*pt(i,j,k) +  &
                          0.5*(phiz(i,k)+phiz(i,k+1)+w(i,j,k)**2+0.5*rsin2_l(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                          v(i,j,k)**2+v(i+1,j,k)**2-(u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s_l(i,j))))
          enddo
       enddo
     endif
     endif
  enddo

!-------------------------------------
! Doganostics computation for moist TE
!-------------------------------------
  if( id_te>0 ) then
!$OMP parallel do default(none) shared(is,ie,js,je,teq,te_2d,moist_phys,km,hlv,sphum,q,delp)
      do j=js,je
         do i=is,ie
            teq(i,j) = te_2d(i,j)
         enddo
         if ( moist_phys ) then
           do k=1,km
              do i=is,ie
                 teq(i,j) = teq(i,j) + hlv*q(i,j,k,sphum)*delp(i,j,k)
              enddo
           enddo
         endif
      enddo
   endif

  end subroutine compute_total_energy


  subroutine pkez(km, ifirst, ilast, jfirst, jlast, j, &
                  pe, pk, akap, peln, pkz, ptop)

! !INPUT PARAMETERS:
   integer, intent(in):: km, j
   integer, intent(in):: ifirst, ilast        ! Latitude strip
   integer, intent(in):: jfirst, jlast        ! Latitude strip
   real(FVPRC), intent(in):: akap
   real(FVPRC), intent(in):: pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1)
   real(FVPRC), intent(in):: pk(ifirst:ilast,jfirst:jlast,km+1)
   real(FVPRC), intent(IN):: ptop
! !OUTPUT
   real(FVPRC), intent(out):: pkz(ifirst:ilast,jfirst:jlast,km)
   real(FVPRC), intent(inout):: peln(ifirst:ilast, km+1, jfirst:jlast)   ! log (pe)
! Local
   real(FVPRC) pk2(ifirst:ilast, km+1)
   real(FVPRC) pek
   real(FVPRC) lnp
   real(FVPRC) ak1
   integer i, k

   ak1 = (akap + 1.) / akap

        pek = pk(ifirst,j,1)
        do i=ifirst, ilast
           pk2(i,1) = pek
        enddo

        do k=2,km+1
           do i=ifirst, ilast
!             peln(i,k,j) =  log(pe(i,k,j))
              pk2(i,k) =  pk(i,j,k)
           enddo
        enddo

!---- GFDL modification
       if( ptop < ptop_min ) then
           do i=ifirst, ilast
               peln(i,1,j) = peln(i,2,j) - ak1
           enddo
       else
           lnp = log( ptop )
           do i=ifirst, ilast
              peln(i,1,j) = lnp
           enddo
       endif
!---- GFDL modification

       do k=1,km
          do i=ifirst, ilast
             pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k) )  /  &
                          (akap*(peln(i,k+1,j) - peln(i,k,j)) )
          enddo
       enddo

 end subroutine pkez



 subroutine remap_z(km, pe1,     kn, pe2, q2, i1, i2, iv, kord)

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension
      integer, intent(in) :: iv

      real(FVPRC), intent(in) ::  pe1(i1:i2,km+1)     ! height at layer edges 
                                               ! (from model top to bottom surface)
      real(FVPRC), intent(in) ::  pe2(i1:i2,kn+1)     ! hieght at layer edges 
                                               ! (from model top to bottom surface)
!      real(FVPRC), intent(in) ::  q1(i1:i2,km)        ! Field input

! !INPUT/OUTPUT PARAMETERS:
      real(FVPRC), intent(inout)::  q2(i1:i2,kn)      ! Field output

! !LOCAL VARIABLES:
      real(FVPRC)   qs(i1:i2)
      real(FVPRC)  dp1(  i1:i2,km)
      real(FVPRC)   q4(4,i1:i2,km)
      real(FVPRC)   pl, pr, qsum, delp, esl
      integer i, k, l, m, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)      ! negative
            q4(1,i,k) = q2(i,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        if (kord == 111) then
           call  cs_profile_linear( qs, q4, dp1, km, i1, i2, iv )
        else
           call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
        endif
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) <= pe1(i,l) .and. pe2(i,k) >= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) >= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) < pe1(i,m+1) ) then
! Whole layer..
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                    delp = pe2(i,k+1)-pe1(i,m)
                    esl = delp / dp1(i,m)
                    qsum = qsum + delp*(q4(2,i,m)+0.5*esl*               &
                         (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                    k0 = m
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue

 end subroutine remap_z

 subroutine map_scalar( km,   pe1,          qs,           &
                      kn,   pe2,    q2,   i1, i2,       &
                      j,    ibeg, iend, jbeg, jend, iv,  kord, q_min)
! iv=1
 integer, intent(in) :: i1                ! Starting longitude
 integer, intent(in) :: i2                ! Finishing longitude
 integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == temp
                                          !       2 == remap temp with cs scheme
 integer, intent(in) :: kord              ! Method order
 integer, intent(in) :: j                 ! Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                ! Original vertical dimension
 integer, intent(in) :: kn                ! Target vertical dimension
 real(FVPRC), intent(in) ::   qs(i1:i2)       ! bottom BC
 real(FVPRC), intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the original vertical coordinate
 real(FVPRC), intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the new vertical coordinate
! real(FVPRC), intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
 real(FVPRC), intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) ! Field output
 real(FVPRC), intent(in):: q_min

! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
   real(FVPRC)    dp1(i1:i2,km)
   real(FVPRC)   q4(4,i1:i2,km)
   real(FVPRC)    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q2(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        if (kord == 111) then
           call scalar_profile_linear( qs, q4, dp1, km, i1, i2, iv, q_min )
        else
           call scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
        endif
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map_scalar


 subroutine map1_ppm( km,   pe1,          qs,           &
                      kn,   pe2,    q2,   i1, i2,       &
                      j,    ibeg, iend, jbeg, jend, iv,  kord)
 integer, intent(in) :: i1                ! Starting longitude
 integer, intent(in) :: i2                ! Finishing longitude
 integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == ???
                                          !       2 == remap temp with cs scheme
 integer, intent(in) :: kord              ! Method order
 integer, intent(in) :: j                 ! Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                ! Original vertical dimension
 integer, intent(in) :: kn                ! Target vertical dimension
 real(FVPRC), intent(in) ::   qs(i1:i2)       ! bottom BC
 real(FVPRC), intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the original vertical coordinate
 real(FVPRC), intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the new vertical coordinate
! real(FVPRC), intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
 real(FVPRC), intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) ! Field output

! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
   real(FVPRC)    dp1(i1:i2,km)
   real(FVPRC)   q4(4,i1:i2,km)
   real(FVPRC)    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q2(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        if (kord == 111) then
           call  cs_profile_linear( qs, q4, dp1, km, i1, i2, iv )
        else
           call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
        endif
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map1_ppm


 subroutine mapn_tracer(nq, km, kn, pe1, pe2, q1, dp2, kord, j,     &
                        i1, i2, isd, ied, jsd, jed, q_min, fill)
! !INPUT PARAMETERS:
      integer, intent(in):: km                ! Original vertical dimension
      integer, intent(in):: kn                ! Target vertical dimension
      integer, intent(in):: j, nq, i1, i2
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: kord(nq)
      real(FVPRC), intent(in)::  pe1(i1:i2,km+1)     ! pressure at layer edges 
                                              ! (from model top to bottom surface)
                                              ! in the original vertical coordinate
      real(FVPRC), intent(in)::  pe2(i1:i2,kn+1)     ! pressure at layer edges 
                                              ! (from model top to bottom surface)
                                              ! in the new vertical coordinate
      real(FVPRC), intent(in)::  dp2(i1:i2,kn)
      real(FVPRC), intent(in)::  q_min
      logical, intent(in):: fill
      real(FVPRC), intent(inout):: q1(isd:ied,jsd:jed,km,nq) ! Field input
! !LOCAL VARIABLES:
      real(FVPRC):: q4(4,i1:i2,km,nq)
      real(FVPRC):: q2(i1:i2,kn,nq) ! Field output
      real(FVPRC):: qsum(nq)
      real(FVPRC):: dp1(i1:i2,km)
      real(FVPRC):: qs(i1:i2)
      real(FVPRC):: pl, pr, dp, esl, fac1, fac2
      integer:: i, k, l, m, k0, iq

      do k=1,km
         do i=i1,i2
            dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         enddo
      enddo

      do iq=1,nq
         do k=1,km
            do i=i1,i2
               q4(1,i,k,iq) = q1(i,j,k,iq)
            enddo
         enddo
	 if (kord(iq) == 111) then
            call scalar_profile_linear( qs, q4(1,i1:i2,1,iq), dp1, km, i1, i2, 0, q_min )
         else
            call scalar_profile( qs, q4(1,i1:i2,1,iq), dp1, km, i1, i2, 0, kord(iq), q_min )
         endif
      enddo

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            fac1 = pr + pl
            fac2 = r3*(pr*fac1 + pl*pl) 
            fac1 = 0.5*fac1
            do iq=1,nq
               q2(i,k,iq) = q4(2,i,l,iq) + (q4(4,i,l,iq)+q4(3,i,l,iq)-q4(2,i,l,iq))*fac1  &
                                         -  q4(4,i,l,iq)*fac2
            enddo
            k0 = l
            goto 555
          else
! Fractional area...
            dp = pe1(i,l+1) - pe2(i,k)
            fac1 = 1. + pl
            fac2 = r3*(1.+pl*fac1)
            fac1 = 0.5*fac1
            do iq=1,nq
               qsum(iq) = dp*(q4(2,i,l,iq) + (q4(4,i,l,iq)+   &
                              q4(3,i,l,iq) - q4(2,i,l,iq))*fac1 - q4(4,i,l,iq)*fac2)
            enddo
            do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
               if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                  do iq=1,nq
                     qsum(iq) = qsum(iq) + dp1(i,m)*q4(1,i,m,iq)
                  enddo
               else
                  dp = pe2(i,k+1)-pe1(i,m)
                  esl = dp / dp1(i,m)
                  fac1 = 0.5*esl
                  fac2 = 1.-r23*esl
                  do iq=1,nq
                     qsum(iq) = qsum(iq) + dp*( q4(2,i,m,iq) + fac1*(         &
                                q4(3,i,m,iq)-q4(2,i,m,iq)+q4(4,i,m,iq)*fac2 ) )
                  enddo
                  k0 = m
                  goto 123
               endif
            enddo
            goto 123
          endif
      endif
100   continue
123   continue
      do iq=1,nq
         q2(i,k,iq) = qsum(iq) / dp2(i,k)
      enddo
555   continue
1000  continue

  !if (fill) call fillz(i2-i1+1, kn, nq, q2, dp2)

  do iq=1,nq
!    if (fill) call fillz(i2-i1+1, kn, 1, q2(i1,1,iq), dp2)
     do k=1,kn
        do i=i1,i2
           q1(i,j,k,iq) = q2(i,k,iq)
        enddo
     enddo
  enddo

 end subroutine mapn_tracer


 subroutine map1_q2(km,   pe1,   q1,            &
                    kn,   pe2,   q2,   dp2,     &
                    i1,   i2,    iv,   kord, j, &
                    ibeg, iend, jbeg, jend, q_min )


! !INPUT PARAMETERS:
      integer, intent(in) :: j
      integer, intent(in) :: i1, i2
      integer, intent(in) :: ibeg, iend, jbeg, jend
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents 1 == ???
      integer, intent(in) :: kord
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real(FVPRC), intent(in) ::  pe1(i1:i2,km+1)     ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real(FVPRC), intent(in) ::  pe2(i1:i2,kn+1)     ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real(FVPRC), intent(in) ::  q1(ibeg:iend,jbeg:jend,km) ! Field input
      real(FVPRC), intent(in) ::  dp2(i1:i2,kn)
      real(FVPRC), intent(in) ::  q_min
! !INPUT/OUTPUT PARAMETERS:
      real(FVPRC), intent(inout):: q2(i1:i2,kn) ! Field output
! !LOCAL VARIABLES:
      real(FVPRC)   qs(i1:i2)
      real(FVPRC)   dp1(i1:i2,km)
      real(FVPRC)   q4(4,i1:i2,km)
      real(FVPRC)   pl, pr, qsum, dp, esl

      integer i, k, l, m, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        if (kord == 111) then
           call  scalar_profile_linear( qs, q4, dp1, km, i1, i2, iv, q_min )
        else
           call  scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
        endif
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                    esl = dp / dp1(i,m)
                   qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                       (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                   k0 = m
                   goto 123
                 endif
              enddo
              goto 123
          endif
      endif
100   continue
123   q2(i,k) = qsum / dp2(i,k)
555   continue
1000  continue

 end subroutine map1_q2

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  map1_cubic --- Cubic Interpolation for vertical re-mapping
!
! !INTERFACE:
  subroutine map1_cubic( km,   pe1,                        &
                         kn,   pe2,    q2,   i1, i2,       &
                         j,    ibeg, iend, jbeg, jend, akap, T_VAR, conserv)
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      real(FVPRC), intent(in) :: akap
      integer, intent(in) :: T_VAR             ! Thermodynamic variable to remap
                                               !     1:TE  2:T  3:PT 
      logical, intent(in) :: conserv
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: ibeg, iend, jbeg, jend
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real(FVPRC), intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real(FVPRC), intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate

!      real(FVPRC), intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
      real(FVPRC), intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) ! Field output

! !DESCRIPTION:
!
!     Perform Cubic Interpolation a given latitude
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY:
!    2005.11.14   Takacs    Initial Code
!    2016.07.20   Putman    Modified to make genaric for any thermodynamic variable
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real(FVPRC)       qx(i1:i2,km)
      real(FVPRC)   logpl1(i1:i2,km)
      real(FVPRC)   logpl2(i1:i2,kn)
      real(FVPRC)   dlogp1(i1:i2,km)
      real(FVPRC)    vsum1(i1:i2)
      real(FVPRC)    vsum2(i1:i2)
      real(FVPRC)   am2,am1,ap0,ap1,P,PLP1,PLP0,PLM1,PLM2,DLP0,DLM1,DLM2

      integer i, k, LM2,LM1,LP0,LP1

! Initialization
! --------------

      select case (T_VAR)
      case(1)
       ! Total Energy Remapping in Log(P)
        do k=1,km
            qx(:,k) = q2(i1:i2,j,k)
        logpl1(:,k) = log( r2*(pe1(:,k)+pe1(:,k+1)) )
        enddo
        do k=1,kn
        logpl2(:,k) = log( r2*(pe2(:,k)+pe2(:,k+1)) )
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

      case(2)
       ! Temperature Remapping in Log(P)
        do k=1,km
            qx(:,k) = q2(i1:i2,j,k)
        logpl1(:,k) = log( r2*(pe1(:,k)+pe1(:,k+1)) )
        enddo
        do k=1,kn
        logpl2(:,k) = log( r2*(pe2(:,k)+pe2(:,k+1)) )
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

      case(3)
       ! Potential Temperature Remapping in P^KAPPA
        do k=1,km
            qx(:,k) = q2(i1:i2,j,k)
        logpl1(:,k) = exp( akap*log( r2*(pe1(:,k)+pe1(:,k+1))) )
        enddo
        do k=1,kn
        logpl2(:,k) = exp( akap*log( r2*(pe2(:,k)+pe2(:,k+1))) )
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

      end select


      if (conserv) then
! Compute vertical integral of Input TE
! -------------------------------------
        vsum1(:) = r0
        do i=i1,i2
        do k=1,km
        vsum1(i) = vsum1(i) + qx(i,k)*( pe1(i,k+1)-pe1(i,k) )
        enddo
        vsum1(i) = vsum1(i) / ( pe1(i,km+1)-pe1(i,1) )
        enddo

      endif

! Interpolate TE onto target Pressures
! ------------------------------------
      do i=i1,i2
      do k=1,kn
         LM1 = 1
         LP0 = 1
         do while( LP0.le.km )
            if (logpl1(i,LP0).lt.logpl2(i,k)) then
               LP0 = LP0+1
            else
               exit
            endif
         enddo
         LM1 = max(LP0-1,1)
         LP0 = min(LP0, km)

! Extrapolate Linearly in LogP above first model level
! ----------------------------------------------------
         if( LM1.eq.1 .and. LP0.eq.1 ) then
             q2(i,j,k) = qx(i,1) + ( qx(i,2)-qx(i,1) )*( logpl2(i,k)-logpl1(i,1) ) &
                                                      /( logpl1(i,2)-logpl1(i,1) )

! Extrapolate Linearly in LogP below last model level
! ---------------------------------------------------
         else if( LM1.eq.km .and. LP0.eq.km ) then
             q2(i,j,k) = qx(i,km) + ( qx(i,km)-qx(i,km-1) )*( logpl2(i,k )-logpl1(i,km  ) ) &
                                                           /( logpl1(i,km)-logpl1(i,km-1) )

! Interpolate Linearly in LogP between levels 1 => 2 and km-1 => km
! -----------------------------------------------------------------
         else if( LM1.eq.1 .or. LP0.eq.km ) then
             q2(i,j,k) = qx(i,LP0) + ( qx(i,LM1)-qx(i,LP0) )*( logpl2(i,k  )-logpl1(i,LP0) ) &
                                                            /( logpl1(i,LM1)-logpl1(i,LP0) )
! Interpolate Cubicly in LogP between other model levels
! ------------------------------------------------------
         else
              LP1 = LP0+1
              LM2 = LM1-1
             P    = logpl2(i,k)
             PLP1 = logpl1(i,LP1)
             PLP0 = logpl1(i,LP0)
             PLM1 = logpl1(i,LM1)
             PLM2 = logpl1(i,LM2)
             DLP0 = dlogp1(i,LP0)
             DLM1 = dlogp1(i,LM1)
             DLM2 = dlogp1(i,LM2)

              ap1 = (P-PLP0)*(P-PLM1)*(P-PLM2)/( DLP0*(DLP0+DLM1)*(DLP0+DLM1+DLM2) )
              ap0 = (PLP1-P)*(P-PLM1)*(P-PLM2)/( DLP0*      DLM1 *(     DLM1+DLM2) )
              am1 = (PLP1-P)*(PLP0-P)*(P-PLM2)/( DLM1*      DLM2 *(DLP0+DLM1     ) )
              am2 = (PLP1-P)*(PLP0-P)*(PLM1-P)/( DLM2*(DLM1+DLM2)*(DLP0+DLM1+DLM2) )

             q2(i,j,k) = ap1*qx(i,LP1) + ap0*qx(i,LP0) + am1*qx(i,LM1) + am2*qx(i,LM2)

         endif

      enddo
      enddo
      if (conserv) then

! Compute vertical integral of Output TE
! --------------------------------------
        vsum2(:) = r0
        do i=i1,i2
        do k=1,kn
        vsum2(i) = vsum2(i) + q2(i,j,k)*( pe2(i,k+1)-pe2(i,k) )
        enddo
        vsum2(i) = vsum2(i) / ( pe2(i,kn+1)-pe2(i,1) )
        enddo

! Adjust Final TE to conserve
! ---------------------------
        do i=i1,i2
        do k=1,kn
           q2(i,j,k) = q2(i,j,k) + vsum1(i)-vsum2(i)
!          q2(i,j,k) = q2(i,j,k) * vsum1(i)/vsum2(i)
        enddo
        enddo

      endif

      return
!EOC
 end subroutine map1_cubic
!-----------------------------------------------------------------------



 subroutine remap_2d(km,   pe1,   q1,        &
                     kn,   pe2,   q2,        &
                     i1,   i2,    iv,   kord)
   integer, intent(in):: i1, i2
   integer, intent(in):: iv               ! Mode: 0 ==  constituents 1 ==others
   integer, intent(in):: kord
   integer, intent(in):: km               ! Original vertical dimension
   integer, intent(in):: kn               ! Target vertical dimension
   real(FVPRC), intent(in):: pe1(i1:i2,km+1)     ! pressure at layer edges 
                                          ! (from model top to bottom surface)
                                          ! in the original vertical coordinate
   real(FVPRC), intent(in):: pe2(i1:i2,kn+1)     ! pressure at layer edges 
                                          ! (from model top to bottom surface)
                                          ! in the new vertical coordinate
   real(FVPRC), intent(in) :: q1(i1:i2,km) ! Field input
   real(FVPRC), intent(out):: q2(i1:i2,kn) ! Field output
! !LOCAL VARIABLES:
   real(FVPRC)   qs(i1:i2)
   real(FVPRC)   dp1(i1:i2,km)
   real(FVPRC)   q4(4,i1:i2,km)
   real(FVPRC)   pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
          dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        if (kord == 111) then
           call  cs_profile_linear( qs, q4, dp1, km, i1, i2, iv )
        else
           call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
        endif
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

   do i=i1,i2
      k0 = 1
      do 555 k=1,kn
!#ifdef OLD_TOP_EDGE
!         if( pe2(i,k+1) <= pe1(i,1) ) then
!! Entire grid above old ptop
!             q2(i,k) = q4(2,i,1)
!         elseif( pe2(i,k) < pe1(i,1) .and. pe2(i,k+1)>pe1(i,1) ) then
!! Partially above old ptop:
!             q2(i,k) = q1(i,1)
!#else
         if( pe2(i,k) <= pe1(i,1) ) then
! above old ptop:
             q2(i,k) = q1(i,1)
!#endif
         else
           do l=k0,km
! locate the top edge: pe2(i,k)
           if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
               pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
               if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
                  pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
                  q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                          *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
                  k0 = l
                  goto 555
               else
! Fractional area...
                 qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                         q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                        (r3*(1.+pl*(1.+pl))))
                 do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                    if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                       qsum = qsum + dp1(i,m)*q4(1,i,m)
                    else
                       dp = pe2(i,k+1)-pe1(i,m)
                      esl = dp / dp1(i,m)
                      qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                            (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                      k0 = m
                      goto 123
                    endif
                 enddo
                 goto 123
               endif
           endif
           enddo
123        q2(i,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
         endif
555   continue
   enddo

 end subroutine remap_2d

 subroutine scalar_profile_linear(qs, a4, delp, km, i1, i2, iv, qmin)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 real(FVPRC), intent(in)   ::   qs(i1:i2)
 real(FVPRC), intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real(FVPRC), intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
 real(FVPRC), intent(in):: qmin
!-----------------------------------------------------------------------
 logical:: extm(i1:i2,km) 
 real(FVPRC)  gam(i1:i2,km)
 real(FVPRC)    q(i1:i2,km+1)
 real(FVPRC)   d4(i1:i2)
 real(FVPRC)   bet, a_bot, grat 
 real(FVPRC)   pmp_1, lac_1, pmp_2, lac_2
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km) 
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo
 
  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif

   do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo

 end subroutine scalar_profile_linear

 subroutine scalar_profile(qs, a4, delp, km, i1, i2, iv, kord, qmin)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: kord
 real(FVPRC), intent(in)   ::   qs(i1:i2)
 real(FVPRC), intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real(FVPRC), intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
 real(FVPRC), intent(in):: qmin
!-----------------------------------------------------------------------
 logical:: extm(i1:i2,km) 
 real(FVPRC)  gam(i1:i2,km)
 real(FVPRC)    q(i1:i2,km+1)
 real(FVPRC)   d4(i1:i2)
 real(FVPRC)   bet, a_bot, grat 
 real(FVPRC)   pmp_1, lac_1, pmp_2, lac_2
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km) 
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo
 
  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif

!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints 
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=2,km-1
     do i=i1,i2
        if ( gam(i,k)*gam(i,k+1) > 0.0 ) then
             extm(i,k) = .false. 
        else
             extm(i,k) = .true.
        endif
     enddo
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then 
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( iv==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( iv/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. (extm(i,k-1).or.extm(i,k+1).or.a4(1,i,k)<qmin) ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( extm(i,k) ) then
              if( a4(1,i,k)<qmin .or. extm(i,k-1) .or. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected; or q is too small -> ehance vertical mixing
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
                   a4(4,i,k) = 0.
              else
! True local extremum
                a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
              endif
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==12 ) then
       do i=i1,i2
          if( extm(i,k) ) then
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==14 ) then
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     else      ! kord = 11, 12, 13
       do i=i1,i2
         if ( extm(i,k) .and. (extm(i,k-1) .or. extm(i,k+1)) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv .eq. -1 ) then 
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine scalar_profile

 subroutine cs_profile_linear(qs, a4, delp, km, i1, i2, iv)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 real(FVPRC), intent(in)   ::   qs(i1:i2)
 real(FVPRC), intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real(FVPRC), intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
!-----------------------------------------------------------------------
 logical:: extm(i1:i2,km) 
 real(FVPRC)  gam(i1:i2,km)
 real(FVPRC)    q(i1:i2,km+1)
 real(FVPRC)   d4(i1:i2)
 real(FVPRC)   bet, a_bot, grat 
 real(FVPRC)   pmp_1, lac_1, pmp_2, lac_2
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km) 
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo
 
  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo

 end subroutine cs_profile_linear

 subroutine cs_profile(qs, a4, delp, km, i1, i2, iv, kord)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: kord
 real(FVPRC), intent(in)   ::   qs(i1:i2)
 real(FVPRC), intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real(FVPRC), intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
!-----------------------------------------------------------------------
 logical:: extm(i1:i2,km) 
 real(FVPRC)  gam(i1:i2,km)
 real(FVPRC)    q(i1:i2,km+1)
 real(FVPRC)   d4(i1:i2)
 real(FVPRC)   bet, a_bot, grat 
 real(FVPRC)   pmp_1, lac_1, pmp_2, lac_2
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km) 
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo
 
  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif

!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints 
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=2,km-1
     do i=i1,i2
        if ( gam(i,k)*gam(i,k+1) > 0.0 ) then
             extm(i,k) = .false. 
        else
             extm(i,k) = .true.
        endif
     enddo
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then 
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( iv==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( iv/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. (extm(i,k-1).or.extm(i,k+1)) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( extm(i,k) ) then
              if( extm(i,k-1) .or. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
                   a4(4,i,k) = 0.
              else
! True local extremum
                a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
              endif
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==12 ) then
       do i=i1,i2
          if( extm(i,k) ) then
! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==14 ) then
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     else      ! kord = 11, 12, 13
       do i=i1,i2
         if ( extm(i,k) .and. (extm(i,k-1) .or. extm(i,k+1)) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv .eq. -1 ) then 
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine cs_profile


 subroutine cs_limiters(im, extm, a4, iv)
 integer, intent(in) :: im
 integer, intent(in) :: iv
 logical, intent(in) :: extm(im)
 real(FVPRC) , intent(inout) :: a4(4,im)   ! PPM array
! !LOCAL VARIABLES:
 real(FVPRC)  da1, da2, a6da
 integer i

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
    if( a4(1,i)<=0.) then
        a4(2,i) = a4(1,i)
        a4(3,i) = a4(1,i)
        a4(4,i) = 0.
    else
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
         if( (a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12) < 0. ) then
! local minimum is negative
             if( a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i) ) then
                 a4(3,i) = a4(1,i)
                 a4(2,i) = a4(1,i)
                 a4(4,i) = 0.
             elseif( a4(3,i) > a4(2,i) ) then
                 a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                 a4(3,i) = a4(2,i) - a4(4,i)
             else
                 a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                 a4(2,i) = a4(3,i) - a4(4,i)
             endif
         endif
      endif
    endif
    enddo
 elseif ( iv==1 ) then
    do i=1,im
      if( (a4(1,i)-a4(2,i))*(a4(1,i)-a4(3,i))>=0. ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
      if( extm(i) ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 endif
 end subroutine cs_limiters



 subroutine ppm_profile(a4, delp, km, i1, i2, iv, kord)

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
                               ! iv = 2: temp (if remap_t) and w (iv=-2)
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real(FVPRC) , intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real(FVPRC) , intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
! 
! !REVISION HISTORY: 
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
      real(FVPRC)    dc(i1:i2,km)
      real(FVPRC)    h2(i1:i2,km)
      real(FVPRC)  delq(i1:i2,km)
      real(FVPRC)   df2(i1:i2,km)
      real(FVPRC)    d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt, it
      real(FVPRC)  fac
      real(FVPRC)  a1, a2, c1, c2, c3, d1, d2
      real(FVPRC)  qm, dq, lac, qmp, pmp

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
                 c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
                 c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            df2(i,k) = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            dc(i,k) = sign( min(abs(df2(i,k)),              &
                            max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))-a4(1,i,k),  &
                  a4(1,i,k)-min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))), df2(i,k) )
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=3,km1
         do i=i1,i2
            c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
            a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
            a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
            a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                      ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                        delp(i,k-1)*a1*dc(i,k  ) )
         enddo
      enddo

      if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
!-------------------------------------------------------
!        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
!-------------------------------------------------------
! No over- and undershoot condition
         a4(2,i,2) = max( a4(2,i,2), min(a4(1,i,1), a4(1,i,2)) )
         a4(2,i,2) = min( a4(2,i,2), max(a4(1,i,1), a4(1,i,2)) )
         dc(i,1) =  0.5*(a4(2,i,2) - a4(1,i,1))
      enddo

! Enforce monotonicity  within the top layer

      if( iv==0 ) then
         do i=i1,i2
            a4(2,i,1) = max(0., a4(2,i,1))
            a4(2,i,2) = max(0., a4(2,i,2))
         enddo 
      elseif( iv==-1 ) then
         do i=i1,i2
            if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
         enddo
      elseif( abs(iv)==2 ) then
         do i=i1,i2
            a4(2,i,1) = a4(1,i,1)
            a4(3,i,1) = a4(1,i,1)
         enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
!        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
!-----------------------------------------------------
!        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
! No over- and under-shoot condition
         a4(2,i,km) = max( a4(2,i,km), min(a4(1,i,km), a4(1,i,km1)) )
         a4(2,i,km) = min( a4(2,i,km), max(a4(1,i,km), a4(1,i,km1)) )
         dc(i,km) = 0.5*(a4(1,i,km) - a4(2,i,km))
      enddo


! Enforce constraint on the "slope" at the surface

!#ifdef BOT_MONO
!      do i=i1,i2
!         a4(4,i,km) = 0
!         if( a4(3,i,km) * a4(1,i,km) <= 0. ) a4(3,i,km) = 0.
!         d1 = a4(1,i,km) - a4(2,i,km)
!         d2 = a4(3,i,km) - a4(1,i,km)
!         if ( d1*d2 < 0. ) then
!              a4(2,i,km) = a4(1,i,km)
!              a4(3,i,km) = a4(1,i,km)
!         else
!              dq = sign(min(abs(d1),abs(d2),0.5*abs(delq(i,km-1))), d1)
!              a4(2,i,km) = a4(1,i,km) - dq
!              a4(3,i,km) = a4(1,i,km) + dq
!         endif
!      enddo
!#else
      if( iv==0 ) then
          do i=i1,i2
             a4(2,i,km) = max(0.,a4(2,i,km))
             a4(3,i,km) = max(0.,a4(3,i,km))
          enddo
      elseif( iv<0 ) then
          do i=i1,i2
             if( a4(1,i,km)*a4(3,i,km) <= 0. )  a4(3,i,km) = 0.
          enddo
      endif
!#endif

   do k=1,km1
      do i=i1,i2
         a4(3,i,k) = a4(2,i,k+1)
      enddo
   enddo

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=2,km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
            h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))  &
                     / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )        &
                     * delp(i,k)**2 
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      fac = 1.5           ! original quasi-monotone

      do k=3,km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), qmp, lac)),    &
                                        max(a4(1,i,k), qmp, lac) )
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         a4(2,i,k) = min(max(a4(2,i,k),  min(a4(1,i,k), qmp, lac)),   &
                     max(a4(1,i,k), qmp, lac))
!-------------
! Recompute A6
!-------------
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord >= 6 )                      &
             call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 2)
      enddo

      else

         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

         do k=3,km-2
            if( kord /= 4) then
              do i=i1,i2
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              enddo
            endif
            if(kord/=6) call ppm_limiters(dc(i1,k), a4(1,i1,k), it, lmt)
         enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

 end subroutine ppm_profile


 subroutine ppm_limiters(dm, a4, itot, lmt)

! !INPUT PARAMETERS:
      real(FVPRC) , intent(in):: dm(*)     ! the linear slope
      integer, intent(in) :: itot      ! Total Longitudes
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real(FVPRC) , intent(inout) :: a4(4,*)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)
! !LOCAL VARIABLES:
      real(FVPRC)  qmp
      real(FVPRC)  da1, da2, a6da
      real(FVPRC)  fmin
      integer i

! Developer: S.-J. Lin

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) < 10.e-10) then !dh fix to prevent tap complaint about real equality test
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

 end subroutine ppm_limiters



 subroutine steepz(i1, i2, km, a4, df2, dm, dq, dp, d4)
 integer, intent(in) :: km, i1, i2
   real(FVPRC) , intent(in) ::  dp(i1:i2,km)       ! grid size
   real(FVPRC) , intent(in) ::  dq(i1:i2,km)       ! backward diff of q
   real(FVPRC) , intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
   real(FVPRC) , intent(in) :: df2(i1:i2,km)       ! first guess mismatch
   real(FVPRC) , intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch
! !INPUT/OUTPUT PARAMETERS:
      real(FVPRC) , intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened
! !LOCAL VARIABLES:
      integer i, k
      real(FVPRC)  alfa(i1:i2,km)
      real(FVPRC)     f(i1:i2,km)
      real(FVPRC)   rat(i1:i2,km)
      real(FVPRC)   dg2

! Compute ratio of dq/dp
      do k=2,km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=2,km-1
         do i=i1,i2
            f(i,k) =   (rat(i,k+1) - rat(i,k))                          &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=3,km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<0. .and. abs(df2(i,k))>10.e-10) then !dh fix to prevenet tap complaint on real equality test
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(0., min(0.5, -0.1875*dg2/df2(i,k)))
         else
            alfa(i,k) = 0.
         endif
         enddo
      enddo

      do k=4,km-2
         do i=i1,i2
            a4(2,i,k) = (1.-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

 end subroutine steepz



 subroutine rst_remap(km, kn, is,ie,js,je, isd,ied,jsd,jed, nq,  &
                      delp_r, u_r, v_r, w_r, delz_r, pt_r, q_r,      &
                      delp,   u,   v,   w,   delz,   pt,   q,        &
                      ak_r, bk_r, ptop, ak, bk, hydrostatic, make_nh, &
                      domain, square_domain)
!------------------------------------
! Assuming hybrid sigma-P coordinate:
!------------------------------------
! !INPUT PARAMETERS:
  integer, intent(in):: km                    ! Restart z-dimension
  integer, intent(in):: kn                    ! Run time dimension
  integer, intent(in):: nq                    ! number of tracers (including h2o)
  integer, intent(in):: is,ie,isd,ied         ! starting & ending X-Dir index
  integer, intent(in):: js,je,jsd,jed         ! starting & ending Y-Dir index
  logical, intent(in):: hydrostatic, make_nh, square_domain
  real(FVPRC), intent(IN) :: ptop
  real(FVPRC), intent(in) :: ak_r(km+1)
  real(FVPRC), intent(in) :: bk_r(km+1)
  real(FVPRC), intent(in) :: ak(kn+1)
  real(FVPRC), intent(in) :: bk(kn+1)
  real(FVPRC), intent(in):: delp_r(is:ie,js:je,km) ! pressure thickness
  real(FVPRC), intent(in)::   u_r(is:ie,  js:je+1,km)   ! u-wind (m/s)
  real(FVPRC), intent(in)::   v_r(is:ie+1,js:je  ,km)   ! v-wind (m/s)
  real(FVPRC), intent(inout)::  pt_r(is:ie,js:je,km)
  real(FVPRC), intent(in)::   w_r(is:ie,js:je,km)
  real(FVPRC), intent(in)::   q_r(is:ie,js:je,km,*)
  real(FVPRC), intent(inout)::delz_r(is:ie,js:je,km)
  type(domain2d), intent(INOUT) :: domain
! Output:
  real(FVPRC), intent(out):: delp(isd:ied,jsd:jed,kn) ! pressure thickness
  real(FVPRC), intent(out)::  u(isd:ied  ,jsd:jed+1,kn)   ! u-wind (m/s)
  real(FVPRC), intent(out)::  v(isd:ied+1,jsd:jed  ,kn)   ! v-wind (m/s)
  real(FVPRC), intent(out)::  w(isd:,jsd:,1:)   ! vertical velocity (m/s)
  real(FVPRC), intent(out):: pt(isd:ied  ,jsd:jed  ,kn)   ! temperature
  real(FVPRC), intent(out):: q(isd:ied,jsd:jed,kn,*)
  real(FVPRC), intent(out):: delz(isd:,jsd:,1:)   ! delta-height (m)
!-----------------------------------------------------------------------
  real(FVPRC) r_vir
  real(FVPRC) ps(isd:ied,jsd:jed)  ! surface pressure
  real(FVPRC)  pe1(is:ie,km+1)
  real(FVPRC)  pe2(is:ie,kn+1)
  real(FVPRC)  pv1(is:ie+1,km+1)
  real(FVPRC)  pv2(is:ie+1,kn+1)

  integer i,j,k , iq
  integer, parameter:: kord=4

  r_vir = rvgas/rdgas - 1.

!$OMP parallel do default(none) shared(is,ie,js,je,ps,ak_r)
  do j=js,je
     do i=is,ie
        ps(i,j) = ak_r(1)
     enddo
  enddo

! this OpenMP do-loop setup cannot work in it's current form....
!$OMP parallel do default(none) shared(is,ie,js,je,km,ps,delp_r)
  do j=js,je
     do k=1,km
        do i=is,ie
           ps(i,j) = ps(i,j) + delp_r(i,j,k)
        enddo
     enddo
  enddo

! only one cell is needed
  if ( square_domain ) then
      call mpp_update_domains(ps, domain,  whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
  else
      call mpp_update_domains(ps, domain, complete=.true.)
  endif

! Compute virtual Temp
!$OMP parallel do default(none) shared(is,ie,js,je,km,pt_r,r_vir,q_r)
  do k=1,km
     do j=js,je
        do i=is,ie
           pt_r(i,j,k) = pt_r(i,j,k) * (1.+r_vir*q_r(i,j,k,1))
        enddo
     enddo
  enddo

!$OMP parallel do default(none) shared(is,ie,js,je,km,ak_r,bk_r,ps,kn,ak,bk,u_r,u,delp, &
!$OMP                                  nq,hydrostatic,make_nh,w_r,w,delz_r,delp_r,delz, &
!$OMP                                  pt_r,pt,v_r,v) &
!$OMP                          private(pe1,  pe2, pv1, pv2)
  do 1000 j=js,je+1
!------
! map u
!------
     do k=1,km+1
        do i=is,ie
           pe1(i,k) = ak_r(k) + 0.5*bk_r(k)*(ps(i,j-1)+ps(i,j))
        enddo
     enddo

     do k=1,kn+1
        do i=is,ie
           pe2(i,k) = ak(k) + 0.5*bk(k)*(ps(i,j-1)+ps(i,j))
        enddo
     enddo

     call remap_2d(km, pe1, u_r(is:ie,j:j,1:km),       &
                   kn, pe2,   u(is:ie,j:j,1:kn),       &
                   is, ie, -1, kord)

  if ( j /= (je+1) )  then 

!---------------
! Hybrid sigma-p
!---------------
     do k=1,km+1
        do i=is,ie
           pe1(i,k) = ak_r(k) + bk_r(k)*ps(i,j)
        enddo
     enddo

     do k=1,kn+1
        do i=is,ie
           pe2(i,k) =   ak(k) + bk(k)*ps(i,j)
        enddo
     enddo

!-------------
! Compute delp
!-------------
      do k=1,kn
         do i=is,ie
            delp(i,j,k) = pe2(i,k+1) - pe2(i,k)
         enddo
      enddo

!----------------
! Map constituents
!----------------
      if( nq /= 0 ) then
          do iq=1,nq
             call remap_2d(km, pe1, q_r(is:ie,j:j,1:km,iq:iq),  &
                           kn, pe2,   q(is:ie,j:j,1:kn,iq:iq),  &
                           is, ie, 0, kord)
          enddo
      endif

      if ( .not. hydrostatic .and. .not. make_nh) then
! Remap vertical wind:
         call remap_2d(km, pe1, w_r(is:ie,j:j,1:km),       &
                       kn, pe2,   w(is:ie,j:j,1:kn),       &
                       is, ie, -1, kord)
! Remap delz for hybrid sigma-p coordinate
         do k=1,km
            do i=is,ie
               delz_r(i,j,k) = -delz_r(i,j,k)/delp_r(i,j,k) ! ="specific volume"/grav
            enddo
         enddo
         call remap_2d(km, pe1, delz_r(is:ie,j:j,1:km),       &
                       kn, pe2,   delz(is:ie,j:j,1:kn),       &
                       is, ie, 1, kord)
         do k=1,kn
            do i=is,ie
               delz(i,j,k) = -delz(i,j,k)*delp(i,j,k)
            enddo
         enddo
      endif

! Geopotential conserving remap of virtual temperature:
       do k=1,km+1
          do i=is,ie
             pe1(i,k) = log(pe1(i,k))
          enddo
       enddo
       do k=1,kn+1
          do i=is,ie
             pe2(i,k) = log(pe2(i,k))
          enddo
       enddo

       call remap_2d(km, pe1, pt_r(is:ie,j:j,1:km),       &
                     kn, pe2,   pt(is:ie,j:j,1:kn),       &
                     is, ie, 1, kord)
!------
! map v
!------
       do k=1,km+1
          do i=is,ie+1
             pv1(i,k) = ak_r(k) + 0.5*bk_r(k)*(ps(i-1,j)+ps(i,j))
          enddo
       enddo
       do k=1,kn+1
          do i=is,ie+1
             pv2(i,k) = ak(k) + 0.5*bk(k)*(ps(i-1,j)+ps(i,j))
          enddo
       enddo

       call remap_2d(km, pv1, v_r(is:ie+1,j:j,1:km),       &
                     kn, pv2,   v(is:ie+1,j:j,1:kn),       &
                     is, ie+1, -1, kord)

  endif !(j < je+1)
1000  continue

!$OMP parallel do default(none) shared(is,ie,js,je,kn,pt,r_vir,q)
  do k=1,kn
     do j=js,je
        do i=is,ie
           pt(i,j,k) = pt(i,j,k) / (1.+r_vir*q(i,j,k,1))
        enddo
     enddo   
  enddo

 end subroutine rst_remap

!#ifdef ANA_REMAP
! subroutine ana_remap_(km, kn, is,ie,js,je, isd,ied,jsd,jed, nq,  &
!                       delp_r, u_r, v_r, ptv_r, q_r,      &
!                       delp,   u  , v  , ptv  , q  ,      &
!                       ak, bk, istatus, phis_r, phis)
!  implicit none
!
!! Todling: Stripped from original S.J.Lin's routine;
!!          remove non-hydrostatic support
!!------------------------------------
!! Assuming hybrid sigma-P coordinate:
!!------------------------------------
!! !INPUT PARAMETERS:
!  integer, intent(in):: km                    ! Restart z-dimension
!  integer, intent(in):: kn                    ! Run time dimension
!  integer, intent(in):: nq                    ! number of tracers (including h2o)
!  integer, intent(in):: is,ie,isd,ied         ! starting & ending X-Dir index
!  integer, intent(in):: js,je,jsd,jed         ! starting & ending Y-Dir index
!  real(REAL8), intent(in) :: ak(kn+1)
!  real(REAL8), intent(in) :: bk(kn+1)
!  real(REAL8), intent(in), optional:: phis_r(is:ie,js:je)   ! topography
!  real(REAL8), intent(in):: delp_r(is:ie  ,js:je  ,km)      ! pressure thickness
!  real(REAL8), intent(in)::    u_r(is:ie  ,js:je+1,km)      ! u-wind (m/s)
!  real(REAL8), intent(in)::    v_r(is:ie+1,js:je  ,km)      ! v-wind (m/s)
!  real(REAL8), intent(in)::  ptv_r(is:ie  ,js:je  ,km)      ! virtual-potential temperature  
!  real(REAL8), intent(in)::    q_r(is:ie  ,js:je  ,km,nq)
!! Input related to desired remapped output:
!  real(REAL8), intent(in), optional:: phis(is:ie,js:je)       ! topography (unghosted)
!! Output:
!  integer, intent(out) :: istatus
!  real(REAL8), intent(inout):: delp(isd:ied  ,jsd:jed  ,kn)   ! pressure thickness
!  real(REAL8), intent(out)::      u(isd:ied  ,jsd:jed+1,kn)   ! u-wind (m/s)
!  real(REAL8), intent(out)::      v(isd:ied+1,jsd:jed  ,kn)   ! v-wind (m/s)
!  real(REAL8), intent(out)::    ptv(isd:ied  ,jsd:jed  ,kn)   ! virtual-potential temperature
!  real(REAL8), intent(out)::      q(isd:ied  ,jsd:jed  ,kn,nq)
!!-----------------------------------------------------------------------
!  real(REAL8) r_vir,cp,kappa
!  real(REAL8) ps_r(isd:ied,jsd:jed)  ! surface pressure
!  real(REAL8) ps  (isd:ied,jsd:jed)  ! surface pressure
!  real(REAL8)  pe1(is:ie,km+1)
!  real(REAL8)  pe2(is:ie,kn+1)
!  real(REAL8)  pv1(is:ie+1,km+1)
!  real(REAL8)  pv2(is:ie+1,kn+1)
!  real(REAL8),allocatable:: ple(:,:,:)
!  real(REAL8),allocatable:: pke(:,:,:)
!  real(REAL8),allocatable:: phi(:,:,:)
!
!  integer i,j,k , iq
!  integer, parameter:: kord=4
!  kappa = MAPL_KAPPA
!  r_vir = rvgas/rdgas - 1.
!  cp    = rvgas/kappa
!
!  istatus=0
!! For now, take only same eta levels between input and output
!  if(km/=kn) then
!     istatus=1
!     return
!  endif
!
!!$omp parallel do default(shared)
!  do j=js,je
!     do i=is,ie
!        ps_r(i,j) = ak(1)
!     enddo
!  enddo
!
!!$omp parallel do default(shared)
!  do k=1,km
!     do j=js,je
!        do i=is,ie
!           ps_r(i,j) = ps_r(i,j) + delp_r(i,j,k)
!        enddo
!     enddo
!  enddo
!  if ( present(phis) .and. present(phis_r) ) then
!
!!    Construct Input Heights
!!    -----------------------
!     allocate(ple(is:ie,js:je,km+1))
!     allocate(pke(is:ie,js:je,km+1))
!     allocate(phi(is:ie,js:je,km+1))
!     ple(:,:,1) = ak(1)
!     do k=2,km+1
!        ple(:,:,k) = ple(:,:,k-1) + delp_r(:,:,k-1)
!     enddo
!     pke = ple**kappa
!
!     phi(:,:,km+1) = phis_r
!     do k=km,1,-1
!        phi(:,:,k) = phi(:,:,k+1) + cp*ptv_r(:,:,k)*( pke(:,:,k+1)-pke(:,:,k) )
!     enddo
!
!!    Compute new surface pressure consistent with output topography
!!    --------------------------------------------------------------
!     do j=js,je
!        do i=is,ie
!           k = km
!           do while ( phi(i,j,k).lt.phis(i,j) )
!              k = k-1
!           enddo
!           ps(i,j) = ple(i,j,k+1)*( 1+(phi(i,j,k+1)-phis(i,j))/(cp*ptv_r(i,j,k)*pke(i,j,k+1)) )**(1.0/kappa)
!        enddo
!     enddo
!     deallocate(phi)
!     deallocate(pke)
!     deallocate(ple)
!
!  else ! when no phis available ...
!!$omp parallel do default(shared)
!     do j=js,je
!        do i=is,ie
!           ps(i,j) = ak(1)
!        enddo
!     enddo
!
!!$omp parallel do default(shared)
!     do k=1,kn
!        do j=js,je
!           do i=is,ie
!              ps(i,j) = ps(i,j) + delp(i,j,k)
!           enddo
!        enddo
!     enddo
!
!  endif ! <phis check>
!
!! only one cell is needed
!  call mpp_update_domains(ps  , domain, complete=.true.)
!  call mpp_update_domains(ps_r, domain, complete=.true.)
!
!!$omp parallel do default(shared) private(pe1,  pe2, pv1, pv2)
!  do 1000 j=js,je+1
!!------
!! map u
!!------
!     do k=1,km+1
!        do i=is,ie
!           pe1(i,k) = ak(k) + 0.5*bk(k)*(ps_r(i,j-1)+ps_r(i,j))
!        enddo
!     enddo
!     do k=1,kn+1
!        do i=is,ie
!           pe2(i,k) = ak(k) + 0.5*bk(k)*(ps(i,j-1)+ps(i,j))
!        enddo
!     enddo
!
!     call remap_2d(km, pe1, u_r(is:ie,j:j,1:km),       &
!                   kn, pe2,   u(is:ie,j:j,1:kn),       &
!                   is, ie, -1, kord)
!
!  if ( j /= (je+1) )  then
!!---------------
!! Hybrid sigma-p
!!---------------
!     do k=1,km+1
!        do i=is,ie
!           pe1(i,k) = ak(k) + bk(k)*ps_r(i,j)
!        enddo
!     enddo
!
!     do k=1,kn+1
!        do i=is,ie
!           pe2(i,k) =   ak(k) + bk(k)*ps(i,j)
!        enddo
!     enddo
!
!!-------------
!! Compute delp
!!-------------
!      do k=1,kn
!         do i=is,ie
!            delp(i,j,k) = pe2(i,k+1) - pe2(i,k)
!         enddo
!      enddo
!
!!----------------
!! Map constituents
!!----------------
!      if( nq /= 0 ) then
!          do iq=1,nq
!             call remap_2d(km, pe1, q_r(is:ie,j:j,1:km,iq:iq),  &
!                           kn, pe2,   q(is:ie,j:j,1:kn,iq:iq),  &
!                           is, ie, 0, kord)
!          enddo
!      endif
!
!! Geopotential conserving remap of virtual temperature:
!       do k=1,km+1
!          do i=is,ie
!             pe1(i,k) = log(pe1(i,k))
!          enddo
!       enddo
!       do k=1,kn+1
!          do i=is,ie
!             pe2(i,k) = log(pe2(i,k))
!          enddo
!       enddo
!
!       call remap_2d(km, pe1,ptv_r(is:ie,j:j,1:km),       &
!                     kn, pe2,  ptv(is:ie,j:j,1:kn),       &
!                     is, ie, 1, kord)
!!------
!! map v
!!------
!       do k=1,km+1
!          do i=is,ie+1
!             pv1(i,k) = ak(k) + 0.5*bk(k)*(ps_r(i-1,j)+ps_r(i,j))
!          enddo
!       enddo
!       do k=1,kn+1
!          do i=is,ie+1
!             pv2(i,k) = ak(k) + 0.5*bk(k)*(ps(i-1,j)+ps(i,j))
!          enddo
!       enddo
!
!       call remap_2d(km, pv1, v_r(is:ie+1,j:j,1:km),       &
!                     kn, pv2,   v(is:ie+1,j:j,1:kn),       &
!                     is, ie+1, -1, kord)
!
!  endif !(j < je+1)
!1000  continue
!
! end subroutine ana_remap_
!#endif



 subroutine mappm(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord, ptop)

! IV = 0: constituents
! IV = 1: potential temp
! IV =-1: winds
 
! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)
 
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

 integer, intent(in):: i1, i2, km, kn, kord, iv
 real(FVPRC), intent(in ):: pe1(i1:i2,km+1), pe2(i1:i2,kn+1)
 real(FVPRC), intent(in )::  q1(i1:i2,km)
 real(FVPRC), intent(out)::  q2(i1:i2,kn)
 real(FVPRC), intent(IN) :: ptop
! local
      real(FVPRC)  qs(i1:i2)
      real(FVPRC) dp1(i1:i2,km)
      real(FVPRC) a4(4,i1:i2,km)
      integer i, k, l
      integer k0, k1
      real(FVPRC) pl, pr, tt, delp, qsum, dpsum, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            a4(1,i,k) = q1(i,k)
         enddo
      enddo

      if ( kord >7 ) then
           if (kord == 111) then
              call  cs_profile_linear( qs, a4, dp1, km, i1, i2, iv )
           else
              call  cs_profile( qs, a4, dp1, km, i1, i2, iv, kord )
           endif
      else
           call ppm_profile( a4, dp1, km, i1, i2, iv, kord )
      endif

!------------------------------------
! Lowest layer: constant distribution
!------------------------------------
      do i=i1,i2
         a4(2,i,km) = q1(i,km)
         a4(3,i,km) = q1(i,km)
         a4(4,i,km) = 0.
      enddo

      do 5555 i=i1,i2
         k0 = 1
      do 555 k=1,kn

         if(pe2(i,k+1) .le. pe1(i,1)) then
! Entire grid above old ptop
            q2(i,k) = a4(2,i,1)
         elseif(pe2(i,k) .ge. pe1(i,km+1)) then
! Entire grid below old ps
            q2(i,k) = a4(3,i,km)
         elseif(pe2(i,k  ) .lt. pe1(i,1) .and.   &
                pe2(i,k+1) .gt. pe1(i,1))  then
! Part of the grid above ptop
            q2(i,k) = a4(1,i,1)
         else

         do 45 L=k0,km
! locate the top edge at pe2(i,k)
         if( pe2(i,k) .ge. pe1(i,L) .and.        &
             pe2(i,k) .le. pe1(i,L+1)    ) then
             k0 = L
             PL = (pe2(i,k)-pe1(i,L)) / dp1(i,L)
             if(pe2(i,k+1) .le. pe1(i,L+1)) then

! entire new grid is within the original grid
               PR = (pe2(i,k+1)-pe1(i,L)) / dp1(i,L)
               TT = r3*(PR*(PR+PL)+PL**2)
               q2(i,k) = a4(2,i,L) + 0.5*(a4(4,i,L)+a4(3,i,L)  &
                       - a4(2,i,L))*(PR+PL) - a4(4,i,L)*TT
              goto 555
             else
! Fractional area...
              delp = pe1(i,L+1) - pe2(i,k)
              TT   = r3*(1.+PL*(1.+PL))
              qsum = delp*(a4(2,i,L)+0.5*(a4(4,i,L)+            &
                     a4(3,i,L)-a4(2,i,L))*(1.+PL)-a4(4,i,L)*TT)
              dpsum = delp
              k1 = L + 1
             goto 111
             endif
         endif
45       continue

111      continue
         do 55 L=k1,km
         if( pe2(i,k+1) .gt. pe1(i,L+1) ) then

! Whole layer..

            qsum  =  qsum + dp1(i,L)*q1(i,L)
            dpsum = dpsum + dp1(i,L)
         else
           delp = pe2(i,k+1)-pe1(i,L)
           esl  = delp / dp1(i,L)
           qsum = qsum + delp * (a4(2,i,L)+0.5*esl*            &
                 (a4(3,i,L)-a4(2,i,L)+a4(4,i,L)*(1.-r23*esl)) )
          dpsum = dpsum + delp
           k0 = L
           goto 123
         endif
55       continue
        delp = pe2(i,k+1) - pe1(i,km+1)
        if(delp > 0.) then
! Extended below old ps
           qsum = qsum + delp * a4(3,i,km)
          dpsum = dpsum + delp
        endif
123     q2(i,k) = qsum / dpsum
      endif
555   continue
5555  continue

 end subroutine mappm


 subroutine remap_z_fb(km, pe1,     kn, pe2, q2, i1, i2, iv, kord)

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension
      integer, intent(in) :: iv

      real(FVPRC), intent(in) ::  pe1(i1:i2,km+1)     ! height at layer edges 
                                               ! (from model top to bottom surface)
      real(FVPRC), intent(in) ::  pe2(i1:i2,kn+1)     ! hieght at layer edges 
                                               ! (from model top to bottom surface)
!      real(FVPRC), intent(in) ::  q1(i1:i2,km)        ! Field input

! !INPUT/OUTPUT PARAMETERS:
      real(FVPRC), intent(inout)::  q2(i1:i2,kn)      ! Field output

! !LOCAL VARIABLES:
      real(FVPRC)   qs(i1:i2)
      real(FVPRC)  dp1(  i1:i2,km)
      real(FVPRC)   q4(4,i1:i2,km)
      real(FVPRC)   pl, pr, qsum, delp, esl
      integer i, k, l, m, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)      ! negative
            q4(1,i,k) = q2(i,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        if (kord == 111) then
           call  cs_profile_linear_fb( qs, q4, dp1, km, i1, i2, iv )
        else
           !call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
        endif
   else
        !call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) <= pe1(i,l) .and. pe2(i,k) >= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) >= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) < pe1(i,m+1) ) then
! Whole layer..
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                    delp = pe2(i,k+1)-pe1(i,m)
                    esl = delp / dp1(i,m)
                    qsum = qsum + delp*(q4(2,i,m)+0.5*esl*               &
                         (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                    k0 = m
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue

 end subroutine remap_z_fb

 subroutine map_scalar_fb( km,   pe1,          qs,           &
                      kn,   pe2,    q2,   i1, i2,       &
                      j,    ibeg, iend, jbeg, jend, iv,  kord, q_min)
! iv=1
 integer, intent(in) :: i1                ! Starting longitude
 integer, intent(in) :: i2                ! Finishing longitude
 integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == temp
                                          !       2 == remap temp with cs scheme
 integer, intent(in) :: kord              ! Method order
 integer, intent(in) :: j                 ! Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                ! Original vertical dimension
 integer, intent(in) :: kn                ! Target vertical dimension
 real(FVPRC), intent(in) ::   qs(i1:i2)       ! bottom BC
 real(FVPRC), intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the original vertical coordinate
 real(FVPRC), intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the new vertical coordinate
! real(FVPRC), intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
 real(FVPRC), intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) ! Field output
 real(FVPRC), intent(in):: q_min

! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
   real(FVPRC)    dp1(i1:i2,km)
   real(FVPRC)   q4(4,i1:i2,km)
   real(FVPRC)    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q2(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        if (kord == 111) then
           call scalar_profile_linear_fb( qs, q4, dp1, km, i1, i2, iv, q_min )
        else
           !call scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
        endif
   else
        !call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map_scalar_fb


 subroutine map1_ppm_fb( km,   pe1,          qs,           &
                      kn,   pe2,    q2,   i1, i2,       &
                      j,    ibeg, iend, jbeg, jend, iv,  kord)
 integer, intent(in) :: i1                ! Starting longitude
 integer, intent(in) :: i2                ! Finishing longitude
 integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == ???
                                          !       2 == remap temp with cs scheme
 integer, intent(in) :: kord              ! Method order
 integer, intent(in) :: j                 ! Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                ! Original vertical dimension
 integer, intent(in) :: kn                ! Target vertical dimension
 real(FVPRC), intent(in) ::   qs(i1:i2)       ! bottom BC
 real(FVPRC), intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the original vertical coordinate
 real(FVPRC), intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the new vertical coordinate
! real(FVPRC), intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
 real(FVPRC), intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) ! Field output

! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
   real(FVPRC)    dp1(i1:i2,km)
   real(FVPRC)   q4(4,i1:i2,km)
   real(FVPRC)    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q2(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        if (kord == 111) then
           call  cs_profile_linear_fb( qs, q4, dp1, km, i1, i2, iv )
        else
           !call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
        endif
   else
        !call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map1_ppm_fb


 subroutine mapn_tracer_fb(nq, km, kn, pe1, pe2, q1, dp2, kord, j,     &
                        i1, i2, isd, ied, jsd, jed, q_min, fill)
! !INPUT PARAMETERS:
      integer, intent(in):: km                ! Original vertical dimension
      integer, intent(in):: kn                ! Target vertical dimension
      integer, intent(in):: j, nq, i1, i2
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: kord(nq)
      real(FVPRC), intent(in)::  pe1(i1:i2,km+1)     ! pressure at layer edges 
                                              ! (from model top to bottom surface)
                                              ! in the original vertical coordinate
      real(FVPRC), intent(in)::  pe2(i1:i2,kn+1)     ! pressure at layer edges 
                                              ! (from model top to bottom surface)
                                              ! in the new vertical coordinate
      real(FVPRC), intent(in)::  dp2(i1:i2,kn)
      real(FVPRC), intent(in)::  q_min
      logical, intent(in):: fill
      real(FVPRC), intent(inout):: q1(isd:ied,jsd:jed,km,nq) ! Field input
! !LOCAL VARIABLES:
      real(FVPRC):: q4(4,i1:i2,km,nq)
      real(FVPRC):: q2(i1:i2,kn,nq) ! Field output
      real(FVPRC):: qsum(nq)
      real(FVPRC):: dp1(i1:i2,km)
      real(FVPRC):: qs(i1:i2)
      real(FVPRC):: pl, pr, dp, esl, fac1, fac2
      integer:: i, k, l, m, k0, iq

      do k=1,km
         do i=i1,i2
            dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         enddo
      enddo

      do iq=1,nq
         do k=1,km
            do i=i1,i2
               q4(1,i,k,iq) = q1(i,j,k,iq)
            enddo
         enddo
	 if (kord(iq) == 111) then
            call scalar_profile_linear_fb( qs, q4(1,i1:i2,1,iq), dp1, km, i1, i2, 0, q_min )
         else
            !call scalar_profile( qs, q4(1,i1:i2,1,iq), dp1, km, i1, i2, 0, kord(iq), q_min )
         endif
      enddo

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            fac1 = pr + pl
            fac2 = r3*(pr*fac1 + pl*pl) 
            fac1 = 0.5*fac1
            do iq=1,nq
               q2(i,k,iq) = q4(2,i,l,iq) + (q4(4,i,l,iq)+q4(3,i,l,iq)-q4(2,i,l,iq))*fac1  &
                                         -  q4(4,i,l,iq)*fac2
            enddo
            k0 = l
            goto 555
          else
! Fractional area...
            dp = pe1(i,l+1) - pe2(i,k)
            fac1 = 1. + pl
            fac2 = r3*(1.+pl*fac1)
            fac1 = 0.5*fac1
            do iq=1,nq
               qsum(iq) = dp*(q4(2,i,l,iq) + (q4(4,i,l,iq)+   &
                              q4(3,i,l,iq) - q4(2,i,l,iq))*fac1 - q4(4,i,l,iq)*fac2)
            enddo
            do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
               if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                  do iq=1,nq
                     qsum(iq) = qsum(iq) + dp1(i,m)*q4(1,i,m,iq)
                  enddo
               else
                  dp = pe2(i,k+1)-pe1(i,m)
                  esl = dp / dp1(i,m)
                  fac1 = 0.5*esl
                  fac2 = 1.-r23*esl
                  do iq=1,nq
                     qsum(iq) = qsum(iq) + dp*( q4(2,i,m,iq) + fac1*(         &
                                q4(3,i,m,iq)-q4(2,i,m,iq)+q4(4,i,m,iq)*fac2 ) )
                  enddo
                  k0 = m
                  goto 123
               endif
            enddo
            goto 123
          endif
      endif
100   continue
123   continue
      do iq=1,nq
         q2(i,k,iq) = qsum(iq) / dp2(i,k)
      enddo
555   continue
1000  continue

  !if (fill) call fillz(i2-i1+1, kn, nq, q2, dp2)

  do iq=1,nq
!    if (fill) call fillz(i2-i1+1, kn, 1, q2(i1,1,iq), dp2)
     do k=1,kn
        do i=i1,i2
           q1(i,j,k,iq) = q2(i,k,iq)
        enddo
     enddo
  enddo

 end subroutine mapn_tracer_fb


 subroutine map1_q2_fb(km,   pe1,   q1,            &
                    kn,   pe2,   q2,   dp2,     &
                    i1,   i2,    iv,   kord, j, &
                    ibeg, iend, jbeg, jend, q_min )


! !INPUT PARAMETERS:
      integer, intent(in) :: j
      integer, intent(in) :: i1, i2
      integer, intent(in) :: ibeg, iend, jbeg, jend
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents 1 == ???
      integer, intent(in) :: kord
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real(FVPRC), intent(in) ::  pe1(i1:i2,km+1)     ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real(FVPRC), intent(in) ::  pe2(i1:i2,kn+1)     ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real(FVPRC), intent(in) ::  q1(ibeg:iend,jbeg:jend,km) ! Field input
      real(FVPRC), intent(in) ::  dp2(i1:i2,kn)
      real(FVPRC), intent(in) ::  q_min
! !INPUT/OUTPUT PARAMETERS:
      real(FVPRC), intent(inout):: q2(i1:i2,kn) ! Field output
! !LOCAL VARIABLES:
      real(FVPRC)   qs(i1:i2)
      real(FVPRC)   dp1(i1:i2,km)
      real(FVPRC)   q4(4,i1:i2,km)
      real(FVPRC)   pl, pr, qsum, dp, esl

      integer i, k, l, m, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        if (kord == 111) then
           call  scalar_profile_linear_fb( qs, q4, dp1, km, i1, i2, iv, q_min )
        else
           !call  scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
        endif
   else
        !call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                    esl = dp / dp1(i,m)
                   qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                       (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                   k0 = m
                   goto 123
                 endif
              enddo
              goto 123
          endif
      endif
100   continue
123   q2(i,k) = qsum / dp2(i,k)
555   continue
1000  continue

 end subroutine map1_q2_fb

 subroutine cs_profile_linear_fb(qs, a4, delp, km, i1, i2, iv)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 real(FVPRC), intent(in)   ::   qs(i1:i2)
 real(FVPRC), intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real(FVPRC), intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
!-----------------------------------------------------------------------
 logical:: extm(i1:i2,km) 
 real(FVPRC)  gam(i1:i2,km)
 real(FVPRC)    q(i1:i2,km+1)
 real(FVPRC)   d4(i1:i2)
 real(FVPRC)   bet, a_bot, grat 
 real(FVPRC)   pmp_1, lac_1, pmp_2, lac_2
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km) 
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo
 
  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo

 end subroutine cs_profile_linear_fb

 subroutine scalar_profile_linear_fb(qs, a4, delp, km, i1, i2, iv, qmin)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 real(FVPRC), intent(in)   ::   qs(i1:i2)
 real(FVPRC), intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real(FVPRC), intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
 real(FVPRC), intent(in):: qmin
!-----------------------------------------------------------------------
 logical:: extm(i1:i2,km) 
 real(FVPRC)  gam(i1:i2,km)
 real(FVPRC)    q(i1:i2,km+1)
 real(FVPRC)   d4(i1:i2)
 real(FVPRC)   bet, a_bot, grat 
 real(FVPRC)   pmp_1, lac_1, pmp_2, lac_2
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km) 
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo
 
  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif

   do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo

 end subroutine scalar_profile_linear_fb

end module fv_mapz_mod
