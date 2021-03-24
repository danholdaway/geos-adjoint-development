module dyn_core_tlm_mod

  use fv_arrays_mod,      only: g_precision, p_precision, pg_precision
  use mpp_domains_mod,    only: CGRID_NE, DGRID_NE, mpp_get_boundary,   &
                                mpp_update_domains
  use fv_mp_mod,          only: domain, isd, ied, jsd, jed, is, ie, js, je, gid, mp_barrier, mp_reduce_max
  use fv_control_mod,     only: hord_mt, hord_vt, hord_tm, hord_dp, hord_ze, hord_tr, n_sponge, m_riem, m_split, &
                                hord_mt_pert, hord_vt_pert, hord_tm_pert, hord_dp_pert, hord_tr_pert, &
                                dddmp, fv_d2_bg=>d2_bg, fv_d4_bg=>d4_bg, d_ext, vtdm4, beta1, beta, init_wind_m, m_grad_p, &
                                a2b_ord, ppm_limiter, master, fv_debug, d_con, fv_nord=>nord, max_courant_no, &
                                no_cgrid, fill_dp, nwat, inline_q, breed_vortex_inline, shallow_water, spin_up_hours
  use sw_core_mod,        only: c_sw, d_sw, divergence_corner, d2a2c
  use sw_core_tlm_mod,    only: c_sw_tlm, d_sw_tlm, divergence_corner_tlm, d2a2c_tlm
  use a2b_edge_mod,       only: a2b_ord2, a2b_ord4
  use a2b_edge_tlm_mod,   only: a2b_ord2_tlm, a2b_ord4_tlm
  use nh_core_mod,        only: Riem_Solver_C, Riem_Solver, update_dz_c, update_dz_d
  use fv_grid_tools_mod,  only: rdx, rdy, rdxc, dxc, dyc, rdyc, dx, dy, area, rarea, grid_type
  use fv_grid_utils_mod,  only: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n,  &
                                ec1, ec2, en1, en2, da_min_c
  use fv_timing_mod,      only: timing_on, timing_off
  use fv_diagnostics_mod, only: prt_maxmin
  use mpp_parameter_mod,  only: CORNER
  use mpp_mod,            only: FATAL, mpp_error

  use test_cases_mod,     only: test_case, case9_forcing1, case9_forcing2
!#ifndef MAPL_MODE
!  use nwp_nudge_mod,      only: breed_slp_inline
!#endif

implicit none

private
public :: dyn_core_tlm

    real   , save :: myStep = 0
    real   , save :: myTime = 0
    logical, save :: first_step = .true.

CONTAINS
!  Differentiation of dyn_core in forward (tangent) mode (with options r8):
!   variations   of useful results: peln q u v w delp ua mfx mfy
!                pkz pe pk pt cx cy
!   with respect to varying inputs: peln q u v w delp ua va pkz
!                pe pk pt
!-----------------------------------------------------------------------
!     dyn_core :: FV Lagrangian dynamics driver
!-----------------------------------------------------------------------
!, init_step, end_step, time_total)
  SUBROUTINE DYN_CORE_TLM(npx, npy, npz, ng, sphum, nq, bdt, n_split, &
&   zvir, cp, akap, rdgas, grav, hydrostatic, u, u_tl, v, v_tl, um, &
&   um_tl, vm, vm_tl, w, w_tl, delz, pt, pt_tl, q, q_tl, delp, delp_tl, &
&   pe, pe_tl, pk, pk_tl, phis, omga, omga_tl, ptop, pfull, ua, ua_tl, &
&   va, va_tl, uc, uc_tl, vc, vc_tl, mfx, mfx_tl, mfy, mfy_tl, cx, cx_tl&
&   , cy, cy_tl, pem, pkz, pkz_tl, peln, peln_tl, ak, bk)
    IMPLICIT NONE
!call timing_off('COMM_TOTAL')
! sub-cycle
!-----------------------------------------------------
! time split loop
!-----------------------------------------------------
!   if ( end_step ) then
!        deallocate(    gz )
!        deallocate(   ptc )
!        deallocate( delzc )
!        deallocate(   crx )
!        deallocate(   xfx )
!        deallocate(   cry )
!        deallocate(   yfx )
!        deallocate ( divg_d )
!        deallocate(   pkc )
!        deallocate( delpc )
!        deallocate(   pkd )
!
!        if ( .not. no_cgrid ) then
!              deallocate( ut )
!              deallocate( vt )
!        endif
!        if ( .not. hydrostatic ) then
!              deallocate( zh )
!              if ( m_grad_p==0 ) deallocate ( pk3 )
!        endif
!        if ( beta > 1.E-4 .or. no_cgrid ) then
!             deallocate( du )
!             deallocate( dv )
!        endif
!   endif
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: ng, nq, sphum
    INTEGER, INTENT(IN) :: n_split
    REAL, INTENT(IN) :: bdt
    REAL, INTENT(IN) :: zvir, cp, akap, rdgas, grav
    REAL(p_precision), INTENT(IN) :: ptop
    REAL, DIMENSION(npz + 1), INTENT(IN) :: ak, bk
    LOGICAL, INTENT(IN) :: hydrostatic
!    logical, intent(IN) :: init_step, end_step
    REAL, INTENT(IN) :: pfull(npz)
! D grid zonal wind (m/s)
    REAL, DIMENSION(isd:ied, jsd:jed + 1, npz), INTENT(INOUT) :: um
    REAL, DIMENSION(isd:ied, jsd:jed+1, npz), INTENT(INOUT) :: um_tl
! D grid meridional wind (m/s)
    REAL, DIMENSION(isd:ied + 1, jsd:jed, npz), INTENT(INOUT) :: vm
    REAL, DIMENSION(isd:ied+1, jsd:jed, npz), INTENT(INOUT) :: vm_tl
! D grid zonal wind (m/s)
    REAL, DIMENSION(isd:ied, jsd:jed + 1, npz), INTENT(INOUT) :: u
    REAL, DIMENSION(isd:ied, jsd:jed+1, npz), INTENT(INOUT) :: u_tl
! D grid meridional wind (m/s)
    REAL, DIMENSION(isd:ied + 1, jsd:jed, npz), INTENT(INOUT) :: v
    REAL, DIMENSION(isd:ied+1, jsd:jed, npz), INTENT(INOUT) :: v_tl
! vertical vel. (m/s)
    REAL, INTENT(INOUT) :: w(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: w_tl(isd:ied, jsd:jed, npz)
! delta-height (m)
    REAL, INTENT(INOUT) :: delz(is:ie, js:je, npz)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: pt_tl(isd:ied, jsd:jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp_tl(isd:ied, jsd:jed, npz)
!
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, npz, nq)
    REAL, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed, npz, nq)
!    real, intent(IN), optional:: time_total  ! total time (seconds) since start
!THESE ARE ACTUALLY IN THE MODULE
    REAL :: delzc(is:ie, js:je, npz)
    REAL :: ut(isd:ied, jsd:jed, npz)
    REAL :: ut_tl(isd:ied, jsd:jed, npz)
    REAL :: vt(isd:ied, jsd:jed, npz)
    REAL :: vt_tl(isd:ied, jsd:jed, npz)
    REAL :: crx(is:ie+1, jsd:jed, npz), xfx(is:ie+1, jsd:jed, npz)
    REAL :: crx_tl(is:ie+1, jsd:jed, npz), xfx_tl(is:ie+1, jsd:jed, npz)
    REAL :: cry(isd:ied, js:je+1, npz), yfx(isd:ied, js:je+1, npz)
    REAL :: cry_tl(isd:ied, js:je+1, npz), yfx_tl(isd:ied, js:je+1, npz)
    REAL :: divg_d(isd:ied+1, jsd:jed+1, npz)
    REAL :: divg_d_tl(isd:ied+1, jsd:jed+1, npz)
    REAL :: zh(isd:ied, jsd:jed, npz)
    REAL(pg_precision) :: du(isd:ied, jsd:jed+1, npz), dv(isd:ied+1, jsd&
&   :jed, npz)
    REAL(pg_precision) :: du_tl(isd:ied, jsd:jed+1, npz), dv_tl(isd:ied+&
&   1, jsd:jed, npz)
    REAL(pg_precision) :: pkc(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision) :: pkc_tl(isd:ied, jsd:jed, npz+1)
    REAL :: pkd(isd:ied, jsd:jed, npz+1)
    REAL :: pkd_tl(isd:ied, jsd:jed, npz+1)
    REAL :: delpc(isd:ied, jsd:jed, npz)
    REAL :: delpc_tl(isd:ied, jsd:jed, npz)
    REAL(pg_precision) :: pk3(isd:ied, jsd:jed, npz+1)
    REAL :: ptc(isd:ied, jsd:jed, npz)
    REAL :: ptc_tl(isd:ied, jsd:jed, npz)
    REAL(pg_precision) :: gz(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision) :: gz_tl(isd:ied, jsd:jed, npz+1)
!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
! Surface geopotential (g*Z_surf)
    REAL, INTENT(INOUT) :: phis(isd:ied, jsd:jed)
! edge pressure (pascal)
    REAL(p_precision), INTENT(INOUT) :: pe(is-1:ie+1, npz+1, js-1:je+1)
    REAL(p_precision), INTENT(INOUT) :: pe_tl(is-1:ie+1, npz+1, js-1:je+&
&   1)
    REAL, INTENT(OUT) :: pem(is-1:ie+1, npz+1, js-1:je+1)
! ln(pe)
    REAL(p_precision), INTENT(OUT) :: peln(is:ie, npz+1, js:je)
    REAL(p_precision), INTENT(OUT) :: peln_tl(is:ie, npz+1, js:je)
! pe**kappa
    REAL(p_precision), INTENT(INOUT) :: pk(is:ie, js:je, npz+1)
    REAL(p_precision), INTENT(INOUT) :: pk_tl(is:ie, js:je, npz+1)
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
! Vertical pressure velocity (pa/s)
    REAL, INTENT(INOUT) :: omga(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: omga_tl(isd:ied, jsd:jed, npz)
! (uc, vc) are mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: uc_tl(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: vc(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: vc_tl(isd:ied, jsd:jed+1, npz)
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua_tl, &
&   va_tl
! The Flux capacitors: accumulated Mass flux arrays
    REAL, INTENT(INOUT) :: mfx(is:ie+1, js:je, npz)
    REAL, INTENT(INOUT) :: mfx_tl(is:ie+1, js:je, npz)
    REAL, INTENT(INOUT) :: mfy(is:ie, js:je+1, npz)
    REAL, INTENT(INOUT) :: mfy_tl(is:ie, js:je+1, npz)
! Accumulated Courant number arrays
    REAL, INTENT(INOUT) :: cx(is:ie+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: cx_tl(is:ie+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: cy(isd:ied, js:je+1, npz)
    REAL, INTENT(INOUT) :: cy_tl(isd:ied, js:je+1, npz)
! Work:
! 
    REAL(p_precision), INTENT(INOUT) :: pkz(is:ie, js:je, npz)
    REAL(p_precision), INTENT(INOUT) :: pkz_tl(is:ie, js:je, npz)
! Auto 1D & 2D arrays:
!real   wbuffer(npy+2,npz)
!real   ebuffer(npy+2,npz)
!real   nbuffer(npx+2,npz)
!real   sbuffer(npx+2,npz)
    REAL :: wbuffer(js:je, npz)
    REAL :: sbuffer(is:ie, npz)
    REAL :: ebuffer(js:je, npz)
    REAL :: nbuffer(is:ie, npz)
! ----   For external mode:
    REAL(pg_precision) :: divg2(is:ie+1, js:je+1)
    REAL(pg_precision) :: divg2_tl(is:ie+1, js:je+1)
    REAL :: wk(isd:ied, jsd:jed)
    REAL :: wk_tl(isd:ied, jsd:jed)
! --- For no_cgrid option ---
    REAL :: u2(isd:ied, jsd:jed+1)
    REAL :: u2_tl(isd:ied, jsd:jed+1)
    REAL :: v2(isd:ied+1, jsd:jed)
    REAL :: v2_tl(isd:ied+1, jsd:jed)
!-------------------------------------
    REAL :: d4_bg, d2_bg
    INTEGER :: hord_m, hord_v, hord_t, hord_p, nord, nord_k
    INTEGER :: hord_m_pert, hord_v_pert, hord_t_pert, hord_p_pert
    INTEGER :: i, j, k, it, iq
    INTEGER :: ism1, iep1, jsm1, jep1
    INTEGER :: ieb1, jeb1
    REAL :: alpha, damp_k
    REAL :: dt2, dt, rdt, rgrav
    REAL :: d2_divg, dd_divg
    LOGICAL :: do_omega
    REAL(p_precision) :: ptk
!--------------------------------------------
! sub-cycle again to recover from instability
!--------------------------------------------
    REAL :: cmax
    REAL :: ndt, elapsed_dt
    INTEGER :: subt
    INTEGER :: subcycle
    INTRINSIC REAL
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC NINT
    INTRINSIC LOG
    INTRINSIC TANH
    INTRINSIC MIN
    REAL :: arg1
    REAL :: arg2
    REAL :: y4
    REAL :: y3
    REAL :: y2
    REAL :: y1
    d2_bg = fv_d2_bg
    d4_bg = fv_d4_bg
    nord = fv_nord
!      if (myTime < 3600.0*spin_up_hours) then
!         nord = 0
!         d2_bg = 0.075
!         d4_bg = 0.0
!         myStep = myStep+1.0
!         myTime = myStep*bdt
!         if (gid==0) print*, 'Spinning Up with NORD:', nord, 'myTime:', myTime
!      endif
    ptk = ptop**akap
    ndt = bdt/REAL(n_split)
    dt = ndt
    dt2 = 0.5*dt
    rdt = 1./dt
    rgrav = 1./grav
! Indexes:
    ism1 = is - 1
    iep1 = ie + 1
    jsm1 = js - 1
    jep1 = je + 1
!call timing_on('COMM_TOTAL')
    IF (npz .GT. 1) then
       call mpp_update_domains( pt,   domain, complete=.false. )
       call mpp_update_domains( pt_tl,   domain, complete=.false. )
    ENDIF
    call mpp_update_domains( delp, domain, complete=.true. )
    call mpp_update_domains( delp_tl, domain, complete=.true. )
    call mpp_update_domains( u, v, domain, gridtype=DGRID_NE, complete=.true. )
    call mpp_update_domains( u_tl, v_tl, domain, gridtype=DGRID_NE, complete=.true. )
!call timing_off('COMM_TOTAL')
!Initiliaze vars
    gz = 0.0
    ptc = 0.0
    delzc = 0.0
    crx = 0.0
    xfx = 0.0
    cry = 0.0
    yfx = 0.0
    divg_d = 0.0
    pkc = 0.0
    delpc = 0.0
    pkd = 0.0
    IF (.NOT.no_cgrid) THEN
      ut(:, :, :) = 0.
      vt(:, :, :) = 0.
    END IF
!         if ( .not. hydrostatic ) then 
!            zh = 0.0
!            if ( m_grad_p==0 ) pk3 = 0.0
!         endif
    IF (beta .GT. 1.e-4 .OR. no_cgrid) THEN
      du = 0.0
      dv = 0.0
    END IF
!   if ( init_step ) then
!
!           allocate(    gz(isd:ied, jsd:jed ,npz+1) )
!           allocate(   ptc(isd:ied, jsd:jed ,npz ) )
!           allocate( delzc(is:ie, js:je ,npz ) )
!           allocate( crx(is :ie+1, jsd:jed,  npz) )
!           allocate( xfx(is :ie+1, jsd:jed,  npz) )
!           allocate( cry(isd:ied,  js :je+1, npz) )
!           allocate( yfx(isd:ied,  js :je+1, npz) )
!           allocate( divg_d(isd:ied+1,jsd:jed+1,npz) )
!           allocate(   pkc(isd:ied, jsd:jed  ,npz+1) )
!           allocate( delpc(isd:ied, jsd:jed  ,npz  ) )
!           allocate(   pkd(isd:ied, jsd:jed  ,npz+1) )
!
!          if ( .not. no_cgrid ) then
!               allocate( ut(isd:ied, jsd:jed, npz) )
!               allocate( vt(isd:ied, jsd:jed, npz) )
!               ut(:,:,:) = 0.
!               vt(:,:,:) = 0.
!          endif
!          if ( .not. hydrostatic ) then 
!               allocate( zh(isd:ied, jsd:jed, npz) )
!               if ( m_grad_p==0 ) allocate ( pk3(isd:ied,jsd:jed,npz+1) )
!          endif
!
!          if ( beta > 1.E-4 .or. no_cgrid ) then
!               allocate( du(isd:ied,  jsd:jed+1,npz) )
!               allocate( dv(isd:ied+1,jsd:jed,  npz) )
!          endif
!      endif
    IF (beta .GT. 1.e-4 .OR. (no_cgrid .AND. init_wind_m)) THEN
      gz_tl = 0.0_8
      pkc_tl = 0.0_8
      CALL GEOPK_TLM(ptop, pe, pe_tl, peln, peln_tl, delp, delp_tl, pkc&
&              , pkc_tl, gz, gz_tl, phis, pt, pt_tl, pkz, pkz_tl, npz, &
&              akap, .false.)
      dv_tl = 0.0_8
      du_tl = 0.0_8
      CALL GRAD1_P_TLM(du, du_tl, dv, dv_tl, pkc, pkc_tl, gz, gz_tl, &
&                delp, delp_tl, dt, ng, npx, npy, npz, ptop, ptk, &
&                hydrostatic)
      IF (init_wind_m) call mpp_update_domains(du, dv, domain, gridtype=DGRID_NE)
      IF (init_wind_m) call mpp_update_domains(du_tl, dv_tl, domain, gridtype=DGRID_NE)
!      IF (init_wind_m) CALL MPP_UPDATE_DOMAINS_DUMMY3_TLM(du, du_tl, isd&
!&                                                   , ied, jsd, jed + 1&
!&                                                   , npz)
!      CALL MPP_UPDATE_DOMAINS_DUMMY3_TLM(dv, dv_tl, isd, ied + 1, jsd, &
!&                                  jed, npz)
    ELSE
      gz_tl = 0.0_8
      du_tl = 0.0_8
      dv_tl = 0.0_8
      pkc_tl = 0.0_8
    END IF
! Empty the "flux capacitors"
    mfx(:, :, :) = 0.
    mfy(:, :, :) = 0.
    cx(:, :, :) = 0.
    cy(:, :, :) = 0.
    elapsed_dt = 0.0
    subcycle = 1
    uc_tl = 0.0_8
    mfx_tl = 0.0_8
    mfy_tl = 0.0_8
    um_tl = 0.0_8
    vc_tl = 0.0_8
    vm_tl = 0.0_8
    cx_tl = 0.0_8
    cy_tl = 0.0_8
    xfx_tl = 0.0_8
    v2_tl = 0.0_8
    ptc_tl = 0.0_8
    ut_tl = 0.0_8
    pkd_tl = 0.0_8
    delpc_tl = 0.0_8
    yfx_tl = 0.0_8
    vt_tl = 0.0_8
    divg_d_tl = 0.0_8
    u2_tl = 0.0_8
    crx_tl = 0.0_8
    cry_tl = 0.0_8
    wk_tl = 0.0_8
    divg2_tl = 0.0_8
!  if (.not. hydrostatic) then
!    ! m_split = max(1, NINT(0.5 + abs(dt)/3.0)) !NINT2CEILING
!    !m_riem = 0
!    !if(m_split >= 2) m_riem = 2
!    !if(m_split >= 4) m_riem = 4
!    !if(gid==0) print*, 'm_split: ', m_split, 'm_riem: ', m_riem
!  endif
!-----------------------------------------------------
    DO it=1,n_split
!-----------------------------------------------------
!--------------------------------------------
! Dynamic sub-cycling to prevent instability
!--------------------------------------------
      IF (max_courant_no .GE. 0.0) THEN
        cmax = 0.0
        DO k=1,npz
          DO j=js,je
            DO i=is,ie+1
! D-Grid meridonal courant number stored in crx
              crx_tl(i, j, k) = ndt*rdy(i, j)*v_tl(i, j, k)
              crx(i, j, k) = v(i, j, k)*ndt*rdy(i, j)
              IF (crx(i, j, k) .GE. 0.) THEN
                y1 = crx(i, j, k)
              ELSE
                y1 = -crx(i, j, k)
              END IF
              IF (cmax .LT. y1) THEN
                cmax = y1
              ELSE
                cmax = cmax
              END IF
            END DO
          END DO
          DO j=js,je+1
            DO i=is,ie
! D-Grid zonal courant number stored in cry
              cry_tl(i, j, k) = ndt*rdx(i, j)*u_tl(i, j, k)
              cry(i, j, k) = u(i, j, k)*ndt*rdx(i, j)
              IF (cry(i, j, k) .GE. 0.) THEN
                y2 = cry(i, j, k)
              ELSE
                y2 = -cry(i, j, k)
              END IF
              IF (cmax .LT. y2) THEN
                cmax = y2
              ELSE
                cmax = cmax
              END IF
            END DO
          END DO
        END DO
        call mp_reduce_max(cmax)
        call mp_barrier()
!NINT2CEILING
        subt = CEILING(cmax/max_courant_no)
        IF (subcycle .NE. subt) THEN
          subcycle = subt
          dt = (bdt-elapsed_dt)/REAL((n_split-(it-1))*subcycle)
          dt2 = 0.5*dt
          rdt = 1./dt
!if (gid==0) write(*,101) 'Sub-Cycling for stability = ', it, subcycle, dt, cmax
!101    format(A,i3,2x,i3,2x,f10.6,2x,f10.6)
!       if (.not. hydrostatic) then
!        ! m_split = max(1, NINT(0.5 + abs(dt)/3.0)) !NINT2CEILING
!        !m_riem = 0
!        !if(m_split >= 2) m_riem = 2
!        !if(m_split >= 4) m_riem = 4
!        !if(gid==0) print*, 'm_split: ', m_split, 'm_riem: ', m_riem
!       endif
        END IF
      END IF
      DO subt=1,subcycle
!--------------------------------------------
!     if ( .not. hydrostatic ) then
!        do j=js,je
!           do i=is,ie
!              zh(i,j,npz) = phis(i,j)*rgrav - delz(i,j,npz)
!           enddo
!           do k=npz-1,1,-1
!              do i=is,ie
!                 zh(i,j,k) = zh(i,j,k+1) - delz(i,j,k)
!              enddo
!           enddo
!        enddo
!                                 !call timing_on('COMM_TOTAL')
!        !oncall mpp_update_domains(zh, domain, complete=.false.)
!        call mpp_update_domains_dummy3(zh,isd,ied, jsd,jed, npz)
!        !oncall mpp_update_domains(w,  domain, complete=.true.)
!        call mpp_update_domains_dummy3(w,isd,ied  ,jsd,jed  ,npz)
!                                !call timing_off('COMM_TOTAL')
!     endif
        IF (shallow_water) THEN
          do_omega = .false.
          !if (test_case==9) call case9_forcing1(phis, time_total)
        ELSE IF (it .EQ. n_split .AND. subt .EQ. subcycle) THEN
!$omp parallel do default(shared) private(i, j, k)
          DO j=jsm1,jep1
            DO i=ism1,iep1
              pem(i, 1, j) = ptop
            END DO
            DO k=1,npz
              DO i=ism1,iep1
                pem(i, k+1, j) = pem(i, k, j) + delp(i, j, k)
              END DO
            END DO
          END DO
          do_omega = .true.
        ELSE
          do_omega = .false.
        END IF
! end no_cgrid section
        IF (no_cgrid) THEN
!---------------------------------------------------------------
! Using time extrapolated wind in place of computed C Grid wind
!---------------------------------------------------------------
!call timing_on('no_cgrid')
          DO k=1,npz
            IF (init_wind_m) THEN
              DO j=jsd,jed+1
                DO i=isd,ied
                  u2_tl(i, j) = u_tl(i, j, k) + 0.5*du_tl(i, j, k)
                  u2(i, j) = u(i, j, k) + 0.5*du(i, j, k)
                END DO
              END DO
              DO j=jsd,jed
                DO i=isd,ied+1
                  v2_tl(i, j) = v_tl(i, j, k) + 0.5*dv_tl(i, j, k)
                  v2(i, j) = v(i, j, k) + 0.5*dv(i, j, k)
                END DO
              END DO
            ELSE
              DO j=jsd,jed+1
                DO i=isd,ied
                  u2_tl(i, j) = 1.5*u_tl(i, j, k) - 0.5*um_tl(i, j, k)
                  u2(i, j) = 1.5*u(i, j, k) - 0.5*um(i, j, k)
                END DO
              END DO
              DO j=jsd,jed
                DO i=isd,ied+1
                  v2_tl(i, j) = 1.5*v_tl(i, j, k) - 0.5*vm_tl(i, j, k)
                  v2(i, j) = 1.5*v(i, j, k) - 0.5*vm(i, j, k)
                END DO
              END DO
            END IF
            CALL D2A2C_TLM(u(isd, jsd, k), u_tl(isd, jsd, k), v(isd, jsd&
&                    , k), v_tl(isd, jsd, k), u2, u2_tl, v2, v2_tl, ua(&
&                    isd, jsd, k), ua_tl(isd, jsd, k), va(isd, jsd, k), &
&                    va_tl(isd, jsd, k), uc(isd, jsd, k), uc_tl(isd, jsd&
&                    , k), vc(isd, jsd, k), vc_tl(isd, jsd, k), nord &
&                    .GT. 0)
          END DO
!if ( .not. hydrostatic ) delpc(:,:,:) = delp(:,:,:)
          um_tl(:, :, :) = u_tl(:, :, :)
          um(:, :, :) = u(:, :, :)
          vm_tl(:, :, :) = v_tl(:, :, :)
          vm(:, :, :) = v(:, :, :)
!call timing_off('no_cgrid')
        ELSE
!call timing_on('c_sw')
!$omp parallel do default(shared) private(i,j,k)
          DO k=1,npz
            CALL C_SW_TLM(delpc(isd, jsd, k), delpc_tl(isd, jsd, k), &
&                   delp(isd, jsd, k), delp_tl(isd, jsd, k), ptc(isd, &
&                   jsd, k), ptc_tl(isd, jsd, k), pt(isd, jsd, k), pt_tl&
&                   (isd, jsd, k), u(isd, jsd, k), u_tl(isd, jsd, k), v(&
&                   isd, jsd, k), v_tl(isd, jsd, k), w(isd, jsd, k), &
&                   w_tl(isd, jsd, k), uc(isd, jsd, k), uc_tl(isd, jsd, &
&                   k), vc(isd, jsd, k), vc_tl(isd, jsd, k), ua(isd, jsd&
&                   , k), ua_tl(isd, jsd, k), va(isd, jsd, k), va_tl(isd&
&                   , jsd, k), omga(isd, jsd, k), ut(isd, jsd, k), ut_tl&
&                   (isd, jsd, k), vt(isd, jsd, k), vt_tl(isd, jsd, k), &
&                   dt2, hydrostatic, nord .GT. 0)
          END DO
! on output omga is updated w
!call timing_off('c_sw')
          IF (fill_dp) THEN
            omga_tl = 0.0_8
            CALL MIX_DP_TLM(hydrostatic, omga, omga_tl, delpc, delpc_tl&
&                     , ptc, ptc_tl, npz, ak, bk, .true.)
          END IF
!      if ( hydrostatic ) then
          IF (beta1 .GT. 0.001) THEN
            alpha = 1. - beta1
            delpc_tl(:, :, :) = beta1*delp_tl(:, :, :) + alpha*delpc_tl(&
&             :, :, :)
            delpc(:, :, :) = beta1*delp(:, :, :) + alpha*delpc(:, :, :)
            ptc_tl(:, :, :) = beta1*pt_tl(:, :, :) + alpha*ptc_tl(:, :, &
&             :)
            ptc(:, :, :) = beta1*pt(:, :, :) + alpha*ptc(:, :, :)
          END IF
          CALL GEOPK_TLM(ptop, pe, pe_tl, peln, peln_tl, delpc, delpc_tl&
&                  , pkc, pkc_tl, gz, gz_tl, phis, ptc, ptc_tl, pkz, &
&                  pkz_tl, npz, akap, .true.)
!      else
!           call update_dz_c(is,   ie, js, je,  npz,    ng,    &
!                            area, zh, ut, vt, delz, delzc, gz)
!                                               !call timing_on('Riem_C')
!           call Riem_Solver_C( dt2,   is,  ie,   js,   je,   npz,   ng,   &
!                               akap, rdgas, grav,  cp,  ptop, phis, omga, delzc, ptc,  &
!                               delpc, gz,  pkc,  1 )
!                                               !call timing_off('Riem_C')
!! pkc is full non-hydro pressure
!                                               !call timing_on('COMM_TOTAL')
!!          call mpp_update_domains(pkc,domain,whalo=1,ehalo=1,shalo=1, nhalo=1, complete=.false.)
!!          call mpp_update_domains(gz ,domain,whalo=1,ehalo=1,shalo=1, nhalo=1, complete=.true.)
!           call mpp_update_domains(pkc, domain, complete=.false.)
!           call mpp_update_domains(gz , domain, complete=.true.)
!                                               !call timing_off('COMM_TOTAL')
!      endif
!-----------------------------------------
! Update time-centered winds on the C-Grid
!-----------------------------------------
          ieb1 = ie + 1
          jeb1 = je + 1
!$omp parallel do default(shared) private(i, j, k, wk)
          DO k=1,npz
!         if ( hydrostatic ) then
            DO j=jsm1,jeb1
              DO i=ism1,ieb1
                wk_tl(i, j) = pkc_tl(i, j, k+1) - pkc_tl(i, j, k)
                wk(i, j) = pkc(i, j, k+1) - pkc(i, j, k)
              END DO
            END DO
!         else
!              do j=jsm1,jeb1
!                 do i=ism1,ieb1
!                       wk(i,j) = delpc(i,j,k)
!                 enddo
!              enddo
!              delpc(:,:,k) = delp(:,:,k) ! Save delp for update_dz_d
!         endif
            DO j=js,je
              DO i=is,ieb1
                uc_tl(i, j, k) = uc_tl(i, j, k) + dt2*rdxc(i, j)*((gz_tl&
&                 (i-1, j, k+1)-gz_tl(i, j, k))*(pkc(i, j, k+1)-pkc(i-1&
&                 , j, k))+(gz(i-1, j, k+1)-gz(i, j, k))*(pkc_tl(i, j, k&
&                 +1)-pkc_tl(i-1, j, k))+(gz_tl(i-1, j, k)-gz_tl(i, j, k&
&                 +1))*(pkc(i-1, j, k+1)-pkc(i, j, k))+(gz(i-1, j, k)-gz&
&                 (i, j, k+1))*(pkc_tl(i-1, j, k+1)-pkc_tl(i, j, k)))/(&
&                 wk(i-1, j)+wk(i, j)) - dt2*rdxc(i, j)*(wk_tl(i-1, j)+&
&                 wk_tl(i, j))*((gz(i-1, j, k+1)-gz(i, j, k))*(pkc(i, j&
&                 , k+1)-pkc(i-1, j, k))+(gz(i-1, j, k)-gz(i, j, k+1))*(&
&                 pkc(i-1, j, k+1)-pkc(i, j, k)))/(wk(i-1, j)+wk(i, j))&
&                 **2
                uc(i, j, k) = uc(i, j, k) + dt2*rdxc(i, j)/(wk(i-1, j)+&
&                 wk(i, j))*((gz(i-1, j, k+1)-gz(i, j, k))*(pkc(i, j, k+&
&                 1)-pkc(i-1, j, k))+(gz(i-1, j, k)-gz(i, j, k+1))*(pkc(&
&                 i-1, j, k+1)-pkc(i, j, k)))
              END DO
            END DO
            DO j=js,jeb1
              DO i=is,ie
                vc_tl(i, j, k) = vc_tl(i, j, k) + dt2*rdyc(i, j)*((gz_tl&
&                 (i, j-1, k+1)-gz_tl(i, j, k))*(pkc(i, j, k+1)-pkc(i, j&
&                 -1, k))+(gz(i, j-1, k+1)-gz(i, j, k))*(pkc_tl(i, j, k+&
&                 1)-pkc_tl(i, j-1, k))+(gz_tl(i, j-1, k)-gz_tl(i, j, k+&
&                 1))*(pkc(i, j-1, k+1)-pkc(i, j, k))+(gz(i, j-1, k)-gz(&
&                 i, j, k+1))*(pkc_tl(i, j-1, k+1)-pkc_tl(i, j, k)))/(wk&
&                 (i, j-1)+wk(i, j)) - dt2*rdyc(i, j)*(wk_tl(i, j-1)+&
&                 wk_tl(i, j))*((gz(i, j-1, k+1)-gz(i, j, k))*(pkc(i, j&
&                 , k+1)-pkc(i, j-1, k))+(gz(i, j-1, k)-gz(i, j, k+1))*(&
&                 pkc(i, j-1, k+1)-pkc(i, j, k)))/(wk(i, j-1)+wk(i, j))&
&                 **2
                vc(i, j, k) = vc(i, j, k) + dt2*rdyc(i, j)/(wk(i, j-1)+&
&                 wk(i, j))*((gz(i, j-1, k+1)-gz(i, j, k))*(pkc(i, j, k+&
&                 1)-pkc(i, j-1, k))+(gz(i, j-1, k)-gz(i, j, k+1))*(pkc(&
&                 i, j-1, k+1)-pkc(i, j, k)))
              END DO
            END DO
          END DO
        END IF
!call timing_on('COMM_TOTAL')
        call mpp_update_domains(uc, vc, domain, gridtype=CGRID_NE, complete=.true.)
        call mpp_update_domains(uc_tl, vc_tl, domain, gridtype=CGRID_NE, complete=.true.)
        !if (test_case==9) call case9_forcing2(phis)
!call timing_on('COMM_TOTAL')
        IF (inline_q) call mpp_update_domains( q,  domain, complete=.true. )
        IF (inline_q) call mpp_update_domains( q_tl,  domain, complete=.true. )
!call timing_off('COMM_TOTAL')
        IF (nord .GT. 0) THEN
          CALL DIVERGENCE_CORNER_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, &
&                              va_tl, divg_d, divg_d_tl, npz)
!call timing_on('COMM_TOTAL')
          call mpp_update_domains(divg_d, domain, position=CORNER)
          call mpp_update_domains(divg_d_tl, domain, position=CORNER)
!call timing_off('COMM_TOTAL')
        END IF
!call timing_on('d_sw')
!$omp parallel do default(shared) private(i, j, k, nord_k, damp_k, d2_divg, dd_divg, hord_m, hord_v, hord_t, hord_p, wk)
        DO k=1,npz
          hord_m = hord_mt
          hord_t = hord_tm
          hord_v = hord_vt
          hord_p = hord_dp
          hord_m_pert = hord_mt_pert
          hord_t_pert = hord_tm_pert
          hord_v_pert = hord_vt_pert
          hord_p_pert = hord_dp_pert
          nord_k = nord
          damp_k = dddmp
          arg1 = pfull(k)/pfull(npz)
          arg2 = 0.1*LOG(arg1)
          y3 = d2_bg*(1.-3.*TANH(arg2))
          IF (0.20 .GT. y3) THEN
            d2_divg = y3
          ELSE
            d2_divg = 0.20
          END IF
          arg1 = pfull(k)/pfull(npz)
          arg2 = 0.1*LOG(arg1)
          y4 = d2_bg*(1.-3.*TANH(arg2))
          IF (0.20 .GT. y4) THEN
            d2_divg = y4
          ELSE
            d2_divg = 0.20
          END IF
          IF (n_sponge .EQ. -1 .OR. npz .EQ. 1) THEN
! Constant divg damping coefficient:
            d2_divg = d2_bg
          ELSE IF (n_sponge .EQ. 0 .AND. k .EQ. 1) THEN
            hord_v = 2
            hord_t = 2
            hord_p = 2
            hord_v_pert = 1
            hord_t_pert = 1
            hord_p_pert = 1
            IF (0 .LT. nord - 1) THEN
              nord_k = nord - 1
            ELSE
              nord_k = 0
            END IF
            damp_k = 0.2
            IF (0.20 .GT. 2.*d2_bg) THEN
              d2_divg = 2.*d2_bg
            ELSE
              d2_divg = 0.20
            END IF
            IF (0.06 .LT. d2_divg) THEN
              d2_divg = d2_divg
            ELSE
              d2_divg = 0.06
            END IF
          ELSE IF (k .LE. n_sponge .AND. npz .GT. 16) THEN
! Apply first order scheme for damping the sponge layer
            hord_m = 1
            hord_v = 1
            hord_t = 1
            hord_p = 1
            hord_m_pert = 1
            hord_v_pert = 1
            hord_t_pert = 1
            hord_p_pert = 1
            nord_k = 0
            damp_k = 0.2
            IF (0.20 .GT. 4.*d2_bg) THEN
              d2_divg = 4.*d2_bg
            ELSE
              d2_divg = 0.20
            END IF
            IF (0.15 .LT. d2_divg) THEN
              d2_divg = d2_divg
            ELSE
              d2_divg = 0.15
            END IF
          ELSE IF (k .EQ. n_sponge + 1 .AND. npz .GT. 24) THEN
            hord_v = 2
            hord_t = 2
            hord_p = 2
            hord_v_pert = 1
            hord_t_pert = 1
            hord_p_pert = 1
            IF (0 .LT. nord - 1) THEN
              nord_k = nord - 1
            ELSE
              nord_k = 0
            END IF
            IF (0.20 .GT. 2.*d2_bg) THEN
              d2_divg = 2.*d2_bg
            ELSE
              d2_divg = 0.20
            END IF
            IF (0.08 .LT. d2_divg) THEN
              d2_divg = d2_divg
            ELSE
              d2_divg = 0.08
            END IF
            IF (nord .GT. 1) THEN
              damp_k = 0.
            ELSE
              damp_k = 0.12
            END IF
          END IF
          dd_divg = d4_bg
          IF (damp_k .LT. dddmp) THEN
            damp_k = dddmp
          ELSE
            damp_k = damp_k
          END IF
!#endif
! if (gid==0 .and. first_step) then
!    print*, k, nord_k, damp_k, d2_divg, dd_divg 
! endif
!--- external mode divergence damping ---
          IF (d_ext .GT. 0.) CALL A2B_ORD2_TLM(delp(:, :, k), delp_tl(:&
&                                        , :, k), wk, wk_tl, npx, npy, &
&                                        is, ie, js, je, ng, .false.)
          CALL D_SW_TLM(pkd(isd, jsd, k), pkd_tl(isd, jsd, k), delp(isd&
&                 , jsd, k), delp_tl(isd, jsd, k), ptc(isd, jsd, k), &
&                 ptc_tl(isd, jsd, k), pt(isd, jsd, k), pt_tl(isd, jsd, &
&                 k), u(isd, jsd, k), u_tl(isd, jsd, k), v(isd, jsd, k)&
&                 , v_tl(isd, jsd, k), w(isd, jsd, k), w_tl(isd, jsd, k)&
&                 , uc(isd, jsd, k), uc_tl(isd, jsd, k), vc(isd, jsd, k)&
&                 , vc_tl(isd, jsd, k), ua(isd, jsd, k), ua_tl(isd, jsd&
&                 , k), va(isd, jsd, k), va_tl(isd, jsd, k), divg_d(isd&
&                 , jsd, k), divg_d_tl(isd, jsd, k), mfx(is, js, k), &
&                 mfx_tl(is, js, k), mfy(is, js, k), mfy_tl(is, js, k), &
&                 cx(is, jsd, k), cx_tl(is, jsd, k), cy(isd, js, k), &
&                 cy_tl(isd, js, k), crx(is, jsd, k), crx_tl(is, jsd, k)&
&                 , cry(isd, js, k), cry_tl(isd, js, k), xfx(is, jsd, k)&
&                 , xfx_tl(is, jsd, k), yfx(isd, js, k), yfx_tl(isd, js&
&                 , k), zvir, sphum, nq, q, q_tl, k, npz, inline_q, pkz(&
&                 is, js, k), pkz_tl(is, js, k), dt, hord_tr, hord_m, &
&                 hord_v, hord_t, hord_p, hord_tr_pert, hord_m_pert, &
&                 hord_v_pert, hord_t_pert, hord_p_pert, nord_k, damp_k, d2_divg, &
&                 dd_divg, vtdm4, d_con, hydrostatic, ppm_limiter)
          IF (d_ext .GT. 0.) THEN
            DO j=js,jep1
              DO i=is,iep1
! delp at cell corners
                ptc_tl(i, j, k) = wk_tl(i, j)
                ptc(i, j, k) = wk(i, j)
              END DO
            END DO
          END IF
        END DO
!call timing_off('d_sw')
        IF (fill_dp) CALL MIX_DP_TLM(hydrostatic, w, w_tl, delp, delp_tl&
&                              , pt, pt_tl, npz, ak, bk, .false.)
!    if ( fv_debug ) then
!         call prt_maxmin('DELP', delp, is, ie, js, je, ng, npz, 1.E-2, master)
!    endif
        IF (d_ext .GT. 0.) THEN
          d2_divg = d_ext*da_min_c
! pkd() is 3D field of horizontal divergence
! ptc is "delp" at cell corners
!$omp parallel do default(shared) private(i,j,k)
          DO j=js,jep1
            DO i=is,iep1
              wk_tl(i, j) = ptc_tl(i, j, 1)
              wk(i, j) = ptc(i, j, 1)
              divg2_tl(i, j) = wk_tl(i, j)*pkd(i, j, 1) + wk(i, j)*&
&               pkd_tl(i, j, 1)
              divg2(i, j) = wk(i, j)*pkd(i, j, 1)
            END DO
            DO k=2,npz
              DO i=is,iep1
                wk_tl(i, j) = wk_tl(i, j) + ptc_tl(i, j, k)
                wk(i, j) = wk(i, j) + ptc(i, j, k)
                divg2_tl(i, j) = divg2_tl(i, j) + ptc_tl(i, j, k)*pkd(i&
&                 , j, k) + ptc(i, j, k)*pkd_tl(i, j, k)
                divg2(i, j) = divg2(i, j) + ptc(i, j, k)*pkd(i, j, k)
              END DO
            END DO
            DO i=is,iep1
              divg2_tl(i, j) = (d2_divg*divg2_tl(i, j)*wk(i, j)-d2_divg*&
&               divg2(i, j)*wk_tl(i, j))/wk(i, j)**2
              divg2(i, j) = d2_divg*divg2(i, j)/wk(i, j)
            END DO
          END DO
        ELSE
          divg2 = 0.
          divg2_tl = 0.0_8
        END IF
!call timing_on('COMM_TOTAL')
        call mpp_update_domains(  pt, domain, complete=.false.)
        call mpp_update_domains(  pt_tl, domain, complete=.false.)
        call mpp_update_domains(delp, domain, complete=.true.)
        call mpp_update_domains(delp_tl, domain, complete=.true.)
! end hydro case
!call timing_off('COMM_TOTAL')
!     if ( hydrostatic ) then
        CALL GEOPK_TLM(ptop, pe, pe_tl, peln, peln_tl, delp, delp_tl, &
&                pkc, pkc_tl, gz, gz_tl, phis, pt, pt_tl, pkz, pkz_tl, &
&                npz, akap, .false.)
!     else
!                                            !call timing_on('UPDATE_DZ')
!          call update_dz_d(hord_tm, is, ie, js, je, npz, ng, npx, npy, area,  &
!                           zh, crx, cry, xfx, yfx, delz, delzc, delpc, n_sponge)
!                                            !call timing_off('UPDATE_DZ')
!                                                          !call timing_on('Riem_D')
!!-----------------------------------------------------------
!! mgrad_p = 1: pkc is full pressure
!! mgrad_p = 0: pkc is non-hydrostatic perturbation pressure
!!-----------------------------------------------------------
!          call Riem_Solver(dt,   is,   ie,   js,   je, npz,  ng,  &
!                           akap, rdgas, grav, cp,   ptop, phis, peln, w,  delz,      &
!                           pt,   delp, gz,   pkc,  pk, pe, ((it==n_split) .and. (subt==subcycle)), m_grad_p)
!                                                 !call timing_off('Riem_D')
!
!                                       !call timing_on('COMM_TOTAL')
!          if ( m_grad_p==0 ) then
!             do k=1,npz+1
!                do j=js,je
!                   do i=is,ie
!                      pk3(i,j,k) = pk(i,j,k)
!                   enddo
!                enddo
!             enddo
!             call mpp_update_domains(pk3, domain, complete=.false.)
!          endif
!
!          call mpp_update_domains(pkc, domain, complete=.true.)
!          call mpp_update_domains(gz , domain, complete=.true.)
!                                       !call timing_off('COMM_TOTAL')
!     endif    ! end hydro case
        IF (.NOT.shallow_water) THEN
          IF (breed_vortex_inline .OR. (it .EQ. n_split .AND. subt .EQ. &
&             subcycle .AND. hydrostatic)) THEN
!$omp parallel do default(shared) private(i, j, k)
            DO k=1,npz+1
              DO j=js,je
                DO i=is,ie
                  pk_tl(i, j, k) = pkc_tl(i, j, k)
                  pk(i, j, k) = pkc(i, j, k)
                END DO
              END DO
            END DO
          END IF
          IF (do_omega) THEN
!------------------------------
! Compute time tendency
!------------------------------
            DO k=1,npz
              DO j=js,je
                DO i=is,ie
                  omga(i, j, k) = (pe(i, k+1, j)-pem(i, k+1, j))*rdt
                END DO
              END DO
            END DO
!------------------------------
! Compute the "advective term"
!------------------------------
            CALL ADV_PE(ua, va, pem, omga, npx, npy, npz, ng)
          END IF
        END IF
!      if ( .not.hydrostatic .and. m_grad_p == 0 ) then
!           call two_grad_p(u, v, pkc, gz, delp, pk3, divg2, dt, ng, npx, npy,   &
!                           npz, ptk)  
!      else
        IF (beta .GT. 1.e-4) THEN
          DO k=1,npz
            DO j=js,je+1
              DO i=is,ie
                u_tl(i, j, k) = rdx(i, j)*(u_tl(i, j, k)+divg2_tl(i, j)-&
&                 divg2_tl(i+1, j)) + beta*du_tl(i, j, k)
                u(i, j, k) = (u(i, j, k)+divg2(i, j)-divg2(i+1, j))*rdx(&
&                 i, j) + beta*du(i, j, k)
              END DO
            END DO
          END DO
          DO k=1,npz
            DO j=js,je
              DO i=is,ie+1
                v_tl(i, j, k) = rdy(i, j)*(v_tl(i, j, k)+divg2_tl(i, j)-&
&                 divg2_tl(i, j+1)) + beta*dv_tl(i, j, k)
                v(i, j, k) = (v(i, j, k)+divg2(i, j)-divg2(i, j+1))*rdy(&
&                 i, j) + beta*dv(i, j, k)
              END DO
            END DO
          END DO
          CALL GRAD1_P_TLM(du, du_tl, dv, dv_tl, pkc, pkc_tl, gz, gz_tl&
&                    , delp, delp_tl, dt, ng, npx, npy, npz, ptop, ptk, &
&                    hydrostatic)
          alpha = 1. - beta
          DO k=1,npz
            DO j=js,je+1
              DO i=is,ie
                u_tl(i, j, k) = u_tl(i, j, k) + alpha*du_tl(i, j, k)
                u(i, j, k) = u(i, j, k) + alpha*du(i, j, k)
              END DO
            END DO
          END DO
          DO k=1,npz
            DO j=js,je
              DO i=is,ie+1
                v_tl(i, j, k) = v_tl(i, j, k) + alpha*dv_tl(i, j, k)
                v(i, j, k) = v(i, j, k) + alpha*dv(i, j, k)
              END DO
            END DO
          END DO
        ELSE
          CALL ONE_GRAD_P_TLM(u, u_tl, v, v_tl, pkc, pkc_tl, gz, gz_tl, &
&                       divg2, divg2_tl, delp, delp_tl, dt, ng, npx, npy&
&                       , npz, ptop, ptk, hydrostatic)
        END IF
!      endif
!-------------------------------------------------------------------------------------------------------
!#ifndef MAPL_MODE
!      if ( breed_vortex_inline )     &
!      call breed_slp_inline(it, dt, npz, ak, bk, phis, pe, pk, peln, delp, u, v, pt, q, nwat, zvir)
!#endif
!-------------------------------------------------------------------------------------------------------
!call timing_on('COMM_TOTAL')
        IF (it .EQ. n_split .AND. subt .EQ. subcycle .AND. grid_type &
&           .LT. 4) THEN
! Prevent accumulation of rounding errors at overlapped domain edges:
          call mpp_get_boundary(u, v, domain, wbuffery=wbuffer, ebuffery=ebuffer,  &
                                sbufferx=sbuffer, nbufferx=nbuffer, gridtype=DGRID_NE )
          DO k=1,npz
            DO i=is,ie
              u(i, je+1, k) = nbuffer(i, k)
            END DO
          END DO
          DO k=1,npz
            DO j=js,je
              v(ie+1, j, k) = ebuffer(j, k)
            END DO
          END DO

          call mpp_get_boundary(u_tl, v_tl, domain, wbuffery=wbuffer, ebuffery=ebuffer,  &
                                sbufferx=sbuffer, nbufferx=nbuffer, gridtype=DGRID_NE )
          DO k=1,npz
            DO i=is,ie
              u_tl(i, je+1, k) = nbuffer(i, k)
            END DO
          END DO
          DO k=1,npz
            DO j=js,je
              v_tl(ie+1, j, k) = ebuffer(j, k)
            END DO
          END DO

        ELSE
          call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
          call mpp_update_domains(u_tl, v_tl, domain, gridtype=DGRID_NE)
        END IF
!call timing_off('COMM_TOTAL')
!init_wind_m = .false. !dh always false
!first_step = .false.
        elapsed_dt = elapsed_dt + dt
      END DO
    END DO
  END SUBROUTINE DYN_CORE_TLM

  SUBROUTINE TWO_GRAD_P_TLM(u, u_tl, v, v_tl, pk, pk_tl, gh, gh_tl, delp&
&   , delp_tl, pkt, pkt_tl, divg2, divg2_tl, dt, ng, npx, npy, npz, ptk)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    REAL(p_precision), INTENT(IN) :: ptk
    REAL(pg_precision), INTENT(IN) :: divg2(is:ie+1, js:je+1)
    REAL(pg_precision), INTENT(IN) :: divg2_tl(is:ie+1, js:je+1)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp_tl(isd:ied, jsd:jed, npz)
! perturbation pressure
    REAL(pg_precision), INTENT(INOUT) :: pk(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: pk_tl(isd:ied, jsd:jed, npz+1)
! p**kappa
    REAL(pg_precision), INTENT(INOUT) :: pkt(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: pkt_tl(isd:ied, jsd:jed, npz+1)
! g * zh
    REAL(pg_precision), INTENT(INOUT) :: gh(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: gh_tl(isd:ied, jsd:jed, npz+1)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: u_tl(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: v_tl(isd:ied+1, jsd:jed, npz)
! Local:
    REAL(pg_precision) :: dp(isd:ied, jsd:jed)
    REAL :: wk1(isd:ied, jsd:jed)
    REAL :: wk1_tl(isd:ied, jsd:jed)
    REAL :: wk(is:ie+1, js:je+1)
    REAL :: wk_tl(is:ie+1, js:je+1)
    INTEGER :: iep1, jep1
    INTEGER :: i, j, k
    iep1 = ie + 1
    jep1 = je + 1
    DO j=js,jep1
      DO i=is,iep1
        pk_tl(i, j, 1) = 0.0
        pk(i, j, 1) = 0.
        pkt_tl(i, j, 1) = 0.0
        pkt(i, j, 1) = ptk
      END DO
    END DO
    wk1_tl = 0.0
    DO k=1,npz+1
      IF (k .NE. 1) THEN
        IF (a2b_ord .EQ. 4) THEN
          CALL A2B_ORD4_TLM(pk(isd:, jsd:, k), pk_tl(isd:, jsd:, k), wk1&
&                     , wk1_tl, npx, npy, is, ie, js, je, ng, .true.)
          CALL A2B_ORD4_TLM(pkt(isd:, jsd:, k), pkt_tl(isd:, jsd:, k), &
&                     wk1, wk1_tl, npx, npy, is, ie, js, je, ng, .true.)
        ELSE
          CALL A2B_ORD2_TLM(pk(isd:, jsd:, k), pk_tl(isd:, jsd:, k), wk1&
&                     , wk1_tl, npx, npy, is, ie, js, je, ng, .true.)
          CALL A2B_ORD2_TLM(pkt(isd:, jsd:, k), pkt_tl(isd:, jsd:, k), &
&                     wk1, wk1_tl, npx, npy, is, ie, js, je, ng, .true.)
        END IF
      END IF
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(gh(isd:, jsd:, k), gh_tl(isd:, jsd:, k), wk1, &
&                   wk1_tl, npx, npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2_TLM(gh(isd:, jsd:, k), gh_tl(isd:, jsd:, k), wk1, &
&                   wk1_tl, npx, npy, is, ie, js, je, ng, .true.)
      END IF
    END DO
    wk_tl = 0.0
    DO k=1,npz
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(delp(:, :, k), delp_tl(:, :, k), wk1, wk1_tl, &
&                   npx, npy, is, ie, js, je, ng)
      ELSE
        CALL A2B_ORD2_TLM(delp(:, :, k), delp_tl(:, :, k), wk1, wk1_tl, &
&                   npx, npy, is, ie, js, je, ng)
      END IF
      DO j=js,jep1
        DO i=is,iep1
          wk_tl(i, j) = pkt_tl(i, j, k+1) - pkt_tl(i, j, k)
          wk(i, j) = pkt(i, j, k+1) - pkt(i, j, k)
        END DO
      END DO
      DO j=js,jep1
        DO i=is,ie
!------------------
! Perturbation term:
!------------------
          u_tl(i, j, k) = u_tl(i, j, k) + dt*((gh_tl(i, j, k+1)-gh_tl(i+&
&           1, j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gh(i, j, k+1)-gh(i+&
&           1, j, k))*(pk_tl(i+1, j, k+1)-pk_tl(i, j, k))+(gh_tl(i, j, k&
&           )-gh_tl(i+1, j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k))+(gh(i, j&
&           , k)-gh(i+1, j, k+1))*(pk_tl(i, j, k+1)-pk_tl(i+1, j, k)))/(&
&           wk1(i, j)+wk1(i+1, j)) - dt*(wk1_tl(i, j)+wk1_tl(i+1, j))*((&
&           gh(i, j, k+1)-gh(i+1, j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(&
&           gh(i, j, k)-gh(i+1, j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))/&
&           (wk1(i, j)+wk1(i+1, j))**2
          u(i, j, k) = u(i, j, k) + dt/(wk1(i, j)+wk1(i+1, j))*((gh(i, j&
&           , k+1)-gh(i+1, j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gh(i, j&
&           , k)-gh(i+1, j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))
!-----------------
! Hydrostatic term
!-----------------
          u_tl(i, j, k) = rdx(i, j)*(divg2_tl(i, j)-divg2_tl(i+1, j)+&
&           u_tl(i, j, k)+dt*((gh_tl(i, j, k+1)-gh_tl(i+1, j, k))*(pkt(i&
&           +1, j, k+1)-pkt(i, j, k))+(gh(i, j, k+1)-gh(i+1, j, k))*(&
&           pkt_tl(i+1, j, k+1)-pkt_tl(i, j, k))+(gh_tl(i, j, k)-gh_tl(i&
&           +1, j, k+1))*(pkt(i, j, k+1)-pkt(i+1, j, k))+(gh(i, j, k)-gh&
&           (i+1, j, k+1))*(pkt_tl(i, j, k+1)-pkt_tl(i+1, j, k)))/(wk(i&
&           , j)+wk(i+1, j))-dt*(wk_tl(i, j)+wk_tl(i+1, j))*((gh(i, j, k&
&           +1)-gh(i+1, j, k))*(pkt(i+1, j, k+1)-pkt(i, j, k))+(gh(i, j&
&           , k)-gh(i+1, j, k+1))*(pkt(i, j, k+1)-pkt(i+1, j, k)))/(wk(i&
&           , j)+wk(i+1, j))**2)
          u(i, j, k) = rdx(i, j)*(divg2(i, j)-divg2(i+1, j)+u(i, j, k)+&
&           dt/(wk(i, j)+wk(i+1, j))*((gh(i, j, k+1)-gh(i+1, j, k))*(pkt&
&           (i+1, j, k+1)-pkt(i, j, k))+(gh(i, j, k)-gh(i+1, j, k+1))*(&
&           pkt(i, j, k+1)-pkt(i+1, j, k))))
        END DO
      END DO
      DO j=js,je
        DO i=is,iep1
!------------------
! Perturbation term:
!------------------
          v_tl(i, j, k) = v_tl(i, j, k) + dt*((gh_tl(i, j, k+1)-gh_tl(i&
&           , j+1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gh(i, j, k+1)-gh(i&
&           , j+1, k))*(pk_tl(i, j+1, k+1)-pk_tl(i, j, k))+(gh_tl(i, j, &
&           k)-gh_tl(i, j+1, k+1))*(pk(i, j, k+1)-pk(i, j+1, k))+(gh(i, &
&           j, k)-gh(i, j+1, k+1))*(pk_tl(i, j, k+1)-pk_tl(i, j+1, k)))/&
&           (wk1(i, j)+wk1(i, j+1)) - dt*(wk1_tl(i, j)+wk1_tl(i, j+1))*(&
&           (gh(i, j, k+1)-gh(i, j+1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+&
&           (gh(i, j, k)-gh(i, j+1, k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))&
&           /(wk1(i, j)+wk1(i, j+1))**2
          v(i, j, k) = v(i, j, k) + dt/(wk1(i, j)+wk1(i, j+1))*((gh(i, j&
&           , k+1)-gh(i, j+1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gh(i, j&
&           , k)-gh(i, j+1, k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))
!-----------------
! Hydrostatic term
!-----------------
          v_tl(i, j, k) = rdy(i, j)*(divg2_tl(i, j)-divg2_tl(i, j+1)+&
&           v_tl(i, j, k)+dt*((gh_tl(i, j, k+1)-gh_tl(i, j+1, k))*(pkt(i&
&           , j+1, k+1)-pkt(i, j, k))+(gh(i, j, k+1)-gh(i, j+1, k))*(&
&           pkt_tl(i, j+1, k+1)-pkt_tl(i, j, k))+(gh_tl(i, j, k)-gh_tl(i&
&           , j+1, k+1))*(pkt(i, j, k+1)-pkt(i, j+1, k))+(gh(i, j, k)-gh&
&           (i, j+1, k+1))*(pkt_tl(i, j, k+1)-pkt_tl(i, j+1, k)))/(wk(i&
&           , j)+wk(i, j+1))-dt*(wk_tl(i, j)+wk_tl(i, j+1))*((gh(i, j, k&
&           +1)-gh(i, j+1, k))*(pkt(i, j+1, k+1)-pkt(i, j, k))+(gh(i, j&
&           , k)-gh(i, j+1, k+1))*(pkt(i, j, k+1)-pkt(i, j+1, k)))/(wk(i&
&           , j)+wk(i, j+1))**2)
          v(i, j, k) = rdy(i, j)*(divg2(i, j)-divg2(i, j+1)+v(i, j, k)+&
&           dt/(wk(i, j)+wk(i, j+1))*((gh(i, j, k+1)-gh(i, j+1, k))*(pkt&
&           (i, j+1, k+1)-pkt(i, j, k))+(gh(i, j, k)-gh(i, j+1, k+1))*(&
&           pkt(i, j, k+1)-pkt(i, j+1, k))))
        END DO
      END DO
    END DO
  END SUBROUTINE TWO_GRAD_P_TLM

  SUBROUTINE ONE_GRAD_P_TLM(u, u_tl, v, v_tl, pk, pk_tl, gh, gh_tl, &
&   divg2, divg2_tl, delp, delp_tl, dt, ng, npx, npy, npz, ptop, ptk, &
&   hydrostatic)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    REAL(p_precision), INTENT(IN) :: ptop, ptk
    LOGICAL, INTENT(IN) :: hydrostatic
    REAL(pg_precision), INTENT(IN) :: divg2(is:ie+1, js:je+1)
    REAL(pg_precision), INTENT(IN) :: divg2_tl(is:ie+1, js:je+1)
    REAL(pg_precision), INTENT(INOUT) :: pk(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: pk_tl(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: gh(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: gh_tl(isd:ied, jsd:jed, npz+1)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp_tl(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: u_tl(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: v_tl(isd:ied+1, jsd:jed, npz)
! Local:
    REAL, DIMENSION(isd:ied, jsd:jed) :: wk
    REAL, DIMENSION(isd:ied, jsd:jed) :: wk_tl
    REAL :: wk1(is:ie+1, js:je)
    REAL :: wk1_tl(is:ie+1, js:je)
    REAL :: wk2(is:ie, js:je+1)
    REAL :: wk2_tl(is:ie, js:je+1)
    REAL :: top_value
    INTEGER :: iep1, jep1
    INTEGER :: i, j, k
    REAL(pg_precision) :: dt_pg
    dt_pg = dt
    iep1 = ie + 1
    jep1 = je + 1
    IF (hydrostatic) THEN
! pk is pe**kappa if hydrostatic
      top_value = ptk
    ELSE
! pk is full pressure if non-hydrostatic
      top_value = ptop
    END IF
    DO j=js,jep1
      DO i=is,iep1
        pk_tl(i, j, 1) = 0.0_8
        pk(i, j, 1) = top_value
      END DO
    END DO
    wk_tl = 0.0_8
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(pk(isd:, jsd:, k), pk_tl(isd:, jsd:, k), wk, &
&                   wk_tl, npx, npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2_TLM(pk(isd:, jsd:, k), pk_tl(isd:, jsd:, k), wk, &
&                   wk_tl, npx, npy, is, ie, js, je, ng, .true.)
      END IF
    END DO
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(gh(isd:, jsd:, k), gh_tl(isd:, jsd:, k), wk, &
&                   wk_tl, npx, npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2_TLM(gh(isd:, jsd:, k), gh_tl(isd:, jsd:, k), wk, &
&                   wk_tl, npx, npy, is, ie, js, je, ng, .true.)
      END IF
    END DO
    wk2_tl = 0.0_8
    DO j=js,jep1
      DO i=is,ie
        wk2_tl(i, j) = divg2_tl(i, j) - divg2_tl(i+1, j)
        wk2(i, j) = divg2(i, j) - divg2(i+1, j)
      END DO
    END DO
    wk1_tl = 0.0_8
    DO j=js,je
      DO i=is,iep1
        wk1_tl(i, j) = divg2_tl(i, j) - divg2_tl(i, j+1)
        wk1(i, j) = divg2(i, j) - divg2(i, j+1)
      END DO
    END DO
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js,jep1
          DO i=is,iep1
            wk_tl(i, j) = pk_tl(i, j, k+1) - pk_tl(i, j, k)
            wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
          END DO
        END DO
      ELSE IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(delp(:, :, k), delp_tl(:, :, k), wk, wk_tl, &
&                   npx, npy, is, ie, js, je, ng)
      ELSE
        CALL A2B_ORD2_TLM(delp(:, :, k), delp_tl(:, :, k), wk, wk_tl, &
&                   npx, npy, is, ie, js, je, ng)
      END IF
      DO j=js,jep1
        DO i=is,ie
          u_tl(i, j, k) = rdx(i, j)*(wk2_tl(i, j)+u_tl(i, j, k)+dt_pg*((&
&           gh_tl(i, j, k+1)-gh_tl(i+1, j, k))*(pk(i+1, j, k+1)-pk(i, j&
&           , k))+(gh(i, j, k+1)-gh(i+1, j, k))*(pk_tl(i+1, j, k+1)-&
&           pk_tl(i, j, k))+(gh_tl(i, j, k)-gh_tl(i+1, j, k+1))*(pk(i, j&
&           , k+1)-pk(i+1, j, k))+(gh(i, j, k)-gh(i+1, j, k+1))*(pk_tl(i&
&           , j, k+1)-pk_tl(i+1, j, k)))/(wk(i, j)+wk(i+1, j))-dt_pg*(&
&           wk_tl(i, j)+wk_tl(i+1, j))*((gh(i, j, k+1)-gh(i+1, j, k))*(&
&           pk(i+1, j, k+1)-pk(i, j, k))+(gh(i, j, k)-gh(i+1, j, k+1))*(&
&           pk(i, j, k+1)-pk(i+1, j, k)))/(wk(i, j)+wk(i+1, j))**2)
          u(i, j, k) = rdx(i, j)*(wk2(i, j)+u(i, j, k)+dt_pg/(wk(i, j)+&
&           wk(i+1, j))*((gh(i, j, k+1)-gh(i+1, j, k))*(pk(i+1, j, k+1)-&
&           pk(i, j, k))+(gh(i, j, k)-gh(i+1, j, k+1))*(pk(i, j, k+1)-pk&
&           (i+1, j, k))))
        END DO
      END DO
      DO j=js,je
        DO i=is,iep1
          v_tl(i, j, k) = rdy(i, j)*(wk1_tl(i, j)+v_tl(i, j, k)+dt_pg*((&
&           gh_tl(i, j, k+1)-gh_tl(i, j+1, k))*(pk(i, j+1, k+1)-pk(i, j&
&           , k))+(gh(i, j, k+1)-gh(i, j+1, k))*(pk_tl(i, j+1, k+1)-&
&           pk_tl(i, j, k))+(gh_tl(i, j, k)-gh_tl(i, j+1, k+1))*(pk(i, j&
&           , k+1)-pk(i, j+1, k))+(gh(i, j, k)-gh(i, j+1, k+1))*(pk_tl(i&
&           , j, k+1)-pk_tl(i, j+1, k)))/(wk(i, j)+wk(i, j+1))-dt_pg*(&
&           wk_tl(i, j)+wk_tl(i, j+1))*((gh(i, j, k+1)-gh(i, j+1, k))*(&
&           pk(i, j+1, k+1)-pk(i, j, k))+(gh(i, j, k)-gh(i, j+1, k+1))*(&
&           pk(i, j, k+1)-pk(i, j+1, k)))/(wk(i, j)+wk(i, j+1))**2)
          v(i, j, k) = rdy(i, j)*(wk1(i, j)+v(i, j, k)+dt_pg/(wk(i, j)+&
&           wk(i, j+1))*((gh(i, j, k+1)-gh(i, j+1, k))*(pk(i, j+1, k+1)-&
&           pk(i, j, k))+(gh(i, j, k)-gh(i, j+1, k+1))*(pk(i, j, k+1)-pk&
&           (i, j+1, k))))
        END DO
      END DO
    END DO
  END SUBROUTINE ONE_GRAD_P_TLM
  SUBROUTINE GRAD1_P_TLM(delu, delu_tl, delv, delv_tl, pk, pk_tl, gh, &
&   gh_tl, delp, delp_tl, dt, ng, npx, npy, npz, ptop, ptk, hydrostatic)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    REAL(p_precision), INTENT(IN) :: ptk
    REAL(p_precision), INTENT(IN) :: ptop
    LOGICAL, INTENT(IN) :: hydrostatic
    REAL(pg_precision), INTENT(INOUT) :: pk(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: pk_tl(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: gh(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: gh_tl(isd:ied, jsd:jed, npz+1)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp_tl(isd:ied, jsd:jed, npz)
    REAL(pg_precision), INTENT(OUT) :: delu(isd:ied, jsd:jed+1, npz)
    REAL(pg_precision), INTENT(OUT) :: delu_tl(isd:ied, jsd:jed+1, npz)
    REAL(pg_precision), INTENT(OUT) :: delv(isd:ied+1, jsd:jed, npz)
    REAL(pg_precision), INTENT(OUT) :: delv_tl(isd:ied+1, jsd:jed, npz)
! Local:
    REAL :: wk(isd:ied, jsd:jed)
    REAL :: wk_tl(isd:ied, jsd:jed)
    REAL :: top_value
    INTEGER :: i, j, k
    REAL(pg_precision) :: dt_pg
    dt_pg = dt
    IF (hydrostatic) THEN
! pk is pe**kappa if hydrostatic
      top_value = ptk
    ELSE
! pk is full pressure if non-hydrostatic
      top_value = ptop
    END IF
    DO j=js,je+1
      DO i=is,ie+1
        pk_tl(i, j, 1) = 0.0_8
        pk(i, j, 1) = top_value
      END DO
    END DO
    wk_tl = 0.0_8
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(pk(isd:, jsd:, k), pk_tl(isd:, jsd:, k), wk, &
&                   wk_tl, npx, npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2_TLM(pk(isd:, jsd:, k), pk_tl(isd:, jsd:, k), wk, &
&                   wk_tl, npx, npy, is, ie, js, je, ng, .true.)
      END IF
    END DO
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(gh(isd:, jsd:, k), gh_tl(isd:, jsd:, k), wk, &
&                   wk_tl, npx, npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD2_TLM(gh(isd:, jsd:, k), gh_tl(isd:, jsd:, k), wk, &
&                   wk_tl, npx, npy, is, ie, js, je, ng, .true.)
      END IF
    END DO
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js,je+1
          DO i=is,ie+1
            wk_tl(i, j) = pk_tl(i, j, k+1) - pk_tl(i, j, k)
            wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
          END DO
        END DO
      ELSE IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4_TLM(delp(:, :, k), delp_tl(:, :, k), wk, wk_tl, &
&                   npx, npy, is, ie, js, je, ng)
      ELSE
        CALL A2B_ORD2_TLM(delp(:, :, k), delp_tl(:, :, k), wk, wk_tl, &
&                   npx, npy, is, ie, js, je, ng)
      END IF
      DO j=js,je+1
        DO i=is,ie
          delu_tl(i, j, k) = rdx(i, j)*dt_pg*((gh_tl(i, j, k+1)-gh_tl(i+&
&           1, j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gh(i, j, k+1)-gh(i+&
&           1, j, k))*(pk_tl(i+1, j, k+1)-pk_tl(i, j, k))+(gh_tl(i, j, k&
&           )-gh_tl(i+1, j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k))+(gh(i, j&
&           , k)-gh(i+1, j, k+1))*(pk_tl(i, j, k+1)-pk_tl(i+1, j, k)))/(&
&           wk(i, j)+wk(i+1, j)) - rdx(i, j)*dt_pg*(wk_tl(i, j)+wk_tl(i+&
&           1, j))*((gh(i, j, k+1)-gh(i+1, j, k))*(pk(i+1, j, k+1)-pk(i&
&           , j, k))+(gh(i, j, k)-gh(i+1, j, k+1))*(pk(i, j, k+1)-pk(i+1&
&           , j, k)))/(wk(i, j)+wk(i+1, j))**2
          delu(i, j, k) = rdx(i, j)*dt_pg/(wk(i, j)+wk(i+1, j))*((gh(i, &
&           j, k+1)-gh(i+1, j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gh(i, &
&           j, k)-gh(i+1, j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          delv_tl(i, j, k) = rdy(i, j)*dt_pg*((gh_tl(i, j, k+1)-gh_tl(i&
&           , j+1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gh(i, j, k+1)-gh(i&
&           , j+1, k))*(pk_tl(i, j+1, k+1)-pk_tl(i, j, k))+(gh_tl(i, j, &
&           k)-gh_tl(i, j+1, k+1))*(pk(i, j, k+1)-pk(i, j+1, k))+(gh(i, &
&           j, k)-gh(i, j+1, k+1))*(pk_tl(i, j, k+1)-pk_tl(i, j+1, k)))/&
&           (wk(i, j)+wk(i, j+1)) - rdy(i, j)*dt_pg*(wk_tl(i, j)+wk_tl(i&
&           , j+1))*((gh(i, j, k+1)-gh(i, j+1, k))*(pk(i, j+1, k+1)-pk(i&
&           , j, k))+(gh(i, j, k)-gh(i, j+1, k+1))*(pk(i, j, k+1)-pk(i, &
&           j+1, k)))/(wk(i, j)+wk(i, j+1))**2
          delv(i, j, k) = rdy(i, j)*dt_pg/(wk(i, j)+wk(i, j+1))*((gh(i, &
&           j, k+1)-gh(i, j+1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gh(i, &
&           j, k)-gh(i, j+1, k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))
        END DO
      END DO
    END DO
  END SUBROUTINE GRAD1_P_TLM
  SUBROUTINE MIX_DP_TLM(hydrostatic, w, w_tl, delp, delp_tl, pt, pt_tl, &
&   km, ak, bk, cg)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL, INTENT(IN) :: ak(km+1), bk(km+1)
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: pt, delp, w
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: pt_tl, &
&   delp_tl, w_tl
    LOGICAL, INTENT(IN) :: hydrostatic, cg
! Local:
    REAL :: dp, dpmin
    REAL :: dp_tl
    INTEGER :: i, j, k, ip
    INTEGER :: ifirst, ilast
    INTEGER :: jfirst, jlast
    IF (cg) THEN
      ifirst = is - 1
      ilast = ie + 1
      jfirst = js - 1
      jlast = je + 1
    ELSE
      ifirst = is
      ilast = ie
      jfirst = js
      jlast = je
    END IF
!$omp parallel do default(shared) private(i, j, k, ip, dpmin, dp)
    DO j=jfirst,jlast
      ip = 0
      DO k=1,km-1
        dpmin = 0.005*(ak(k+1)-ak(k)+(bk(k+1)-bk(k))*1.e5)
        DO i=ifirst,ilast
          IF (delp(i, j, k) .LT. dpmin) THEN
! Remap from below and mix pt
            dp_tl = -delp_tl(i, j, k)
            dp = dpmin - delp(i, j, k)
            pt_tl(i, j, k) = (pt_tl(i, j, k)*delp(i, j, k)+pt(i, j, k)*&
&             delp_tl(i, j, k)+pt_tl(i, j, k+1)*dp+pt(i, j, k+1)*dp_tl)/&
&             dpmin
            pt(i, j, k) = (pt(i, j, k)*delp(i, j, k)+pt(i, j, k+1)*dp)/&
&             dpmin
            IF (.NOT.hydrostatic) THEN
              w_tl(i, j, k) = (w_tl(i, j, k)*delp(i, j, k)+w(i, j, k)*&
&               delp_tl(i, j, k)+w_tl(i, j, k+1)*dp+w(i, j, k+1)*dp_tl)/&
&               dpmin
              w(i, j, k) = (w(i, j, k)*delp(i, j, k)+w(i, j, k+1)*dp)/&
&               dpmin
            END IF
            delp_tl(i, j, k) = 0.0_8
            delp(i, j, k) = dpmin
            delp_tl(i, j, k+1) = delp_tl(i, j, k+1) - dp_tl
            delp(i, j, k+1) = delp(i, j, k+1) - dp
            ip = ip + 1
          END IF
        END DO
      END DO
! Bottom (k=km):
      dpmin = 0.005*(ak(km+1)-ak(km)+(bk(km+1)-bk(km))*1.e5)
      DO i=ifirst,ilast
        IF (delp(i, j, km) .LT. dpmin) THEN
! Remap from above and mix pt
          dp_tl = -delp_tl(i, j, km)
          dp = dpmin - delp(i, j, km)
          pt_tl(i, j, km) = (pt_tl(i, j, km)*delp(i, j, km)+pt(i, j, km)&
&           *delp_tl(i, j, km)+pt_tl(i, j, km-1)*dp+pt(i, j, km-1)*dp_tl&
&           )/dpmin
          pt(i, j, km) = (pt(i, j, km)*delp(i, j, km)+pt(i, j, km-1)*dp)&
&           /dpmin
          IF (.NOT.hydrostatic) THEN
            w_tl(i, j, km) = (w_tl(i, j, km)*delp(i, j, km)+w(i, j, km)*&
&             delp_tl(i, j, km)+w_tl(i, j, km-1)*dp+w(i, j, km-1)*dp_tl)&
&             /dpmin
            w(i, j, km) = (w(i, j, km)*delp(i, j, km)+w(i, j, km-1)*dp)/&
&             dpmin
          END IF
          delp_tl(i, j, km) = 0.0_8
          delp(i, j, km) = dpmin
          delp_tl(i, j, km-1) = delp_tl(i, j, km-1) - dp_tl
          delp(i, j, km-1) = delp(i, j, km-1) - dp
          ip = ip + 1
        END IF
      END DO
!      if ( fv_debug .and. ip/=0 ) write(*,*) 'Warning: Mix_dp', gid, j, ip 
      IF (ip .NE. 0) WRITE(*, *) 'Warning: Mix_dp', gid, j, ip
    END DO
  END SUBROUTINE MIX_DP_TLM
  SUBROUTINE GEOPK_TLM(ptop, pe, pe_tl, peln, peln_tl, delp, delp_tl, pk&
&   , pk_tl, gh, gh_tl, hs, pt, pt_tl, pkz, pkz_tl, km, akap, cg)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL, INTENT(IN) :: akap
    REAL(p_precision), INTENT(IN) :: ptop
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: pt, delp
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: pt_tl, delp_tl
    LOGICAL, INTENT(IN) :: cg
! !OUTPUT PARAMETERS
    REAL(pg_precision), DIMENSION(isd:ied, jsd:jed, km + 1), INTENT(OUT)&
&   :: gh
    REAL(pg_precision), DIMENSION(isd:ied, jsd:jed, km+1), INTENT(OUT) &
&   :: gh_tl
    REAL(pg_precision), DIMENSION(isd:ied, jsd:jed, km + 1), INTENT(OUT)&
&   :: pk
    REAL(pg_precision), DIMENSION(isd:ied, jsd:jed, km+1), INTENT(OUT) &
&   :: pk_tl
    REAL(p_precision), INTENT(OUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL(p_precision), INTENT(OUT) :: pe_tl(is-1:ie+1, km+1, js-1:je+1)
! ln(pe)
    REAL(p_precision), INTENT(OUT) :: peln(is:ie, km+1, js:je)
    REAL(p_precision), INTENT(OUT) :: peln_tl(is:ie, km+1, js:je)
    REAL(p_precision), INTENT(OUT) :: pkz(is:ie, js:je, km)
    REAL(p_precision), INTENT(OUT) :: pkz_tl(is:ie, js:je, km)
! !DESCRIPTION:
!    Calculates geopotential and pressure to the kappa.
! Local:
    REAL(pg_precision) :: p1d(is-2:ie+2)
    REAL(pg_precision) :: p1d_tl(is-2:ie+2)
    REAL(pg_precision) :: logp(is-2:ie+2)
    REAL(pg_precision) :: logp_tl(is-2:ie+2)
    REAL(pg_precision) :: ptk
    INTEGER :: i, j, k
    INTEGER :: ifirst, ilast
    INTEGER :: jfirst, jlast
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC LOG
    INTRINSIC EXP
    INTEGER :: min2
    INTEGER :: min1
    INTEGER :: max2
    INTEGER :: max1
    ptk = ptop**akap
    IF (.NOT.cg .AND. a2b_ord .EQ. 4) THEN
! D-Grid
      ifirst = is - 2
      ilast = ie + 2
      jfirst = js - 2
      jlast = je + 2
      logp_tl = 0.0_8
      p1d_tl = 0.0_8
    ELSE
      ifirst = is - 1
      ilast = ie + 1
      jfirst = js - 1
      jlast = je + 1
      logp_tl = 0.0_8
      p1d_tl = 0.0_8
    END IF
!$omp parallel do default(shared) private(i, j, k, p1d, dp, logp)
    DO j=jfirst,jlast
      DO i=ifirst,ilast
        p1d_tl(i) = 0.0_8
        p1d(i) = ptop
        pk_tl(i, j, 1) = 0.0_8
        pk(i, j, 1) = ptk
        gh_tl(i, j, km+1) = 0.0_8
        gh(i, j, km+1) = hs(i, j)
      END DO
      IF (j .GT. js - 2 .AND. j .LT. je + 2) THEN
        IF (ifirst .LT. is - 1) THEN
          max1 = is - 1
        ELSE
          max1 = ifirst
        END IF
        IF (ilast .GT. ie + 1) THEN
          min1 = ie + 1
        ELSE
          min1 = ilast
        END IF
        DO i=max1,min1
          pe_tl(i, 1, j) = 0.0_8
          pe(i, 1, j) = ptop
        END DO
      END IF
! Top down
      DO k=2,km+1
        DO i=ifirst,ilast
          p1d_tl(i) = p1d_tl(i) + delp_tl(i, j, k-1)
          p1d(i) = p1d(i) + delp(i, j, k-1)
!            pk(i,j,k) = p1d(i) ** akap
! Optimized form:
          logp_tl(i) = p1d_tl(i)/p1d(i)
          logp(i) = LOG(p1d(i))
          pk_tl(i, j, k) = akap*logp_tl(i)*EXP(akap*logp(i))
          pk(i, j, k) = EXP(akap*logp(i))
        END DO
        IF (j .GT. js - 2 .AND. j .LT. je + 2) THEN
          IF (ifirst .LT. is - 1) THEN
            max2 = is - 1
          ELSE
            max2 = ifirst
          END IF
          IF (ilast .GT. ie + 1) THEN
            min2 = ie + 1
          ELSE
            min2 = ilast
          END IF
          DO i=max2,min2
            pe_tl(i, k, j) = p1d_tl(i)
            pe(i, k, j) = p1d(i)
          END DO
          IF (j .GE. js .AND. j .LE. je) THEN
            DO i=is,ie
              peln_tl(i, k, j) = logp_tl(i)
              peln(i, k, j) = logp(i)
            END DO
          END IF
        END IF
      END DO
! Bottom up
      DO k=km,1,-1
        DO i=ifirst,ilast
          gh_tl(i, j, k) = gh_tl(i, j, k+1) + pt_tl(i, j, k)*(pk(i, j, k&
&           +1)-pk(i, j, k)) + pt(i, j, k)*(pk_tl(i, j, k+1)-pk_tl(i, j&
&           , k))
          gh(i, j, k) = gh(i, j, k+1) + pt(i, j, k)*(pk(i, j, k+1)-pk(i&
&           , j, k))
        END DO
      END DO
    END DO
    IF (.NOT.cg) THEN
! This is for hydrostatic only
      DO k=1,km
        DO j=js,je
          DO i=is,ie
            pkz_tl(i, j, k) = ((pk_tl(i, j, k+1)-pk_tl(i, j, k))*akap*(&
&             peln(i, k+1, j)-peln(i, k, j))-(pk(i, j, k+1)-pk(i, j, k))&
&             *akap*(peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(akap*(peln(i&
&             , k+1, j)-peln(i, k, j)))**2
            pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(akap*(peln(i, k+&
&             1, j)-peln(i, k, j)))
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE GEOPK_TLM
  SUBROUTINE ADV_PE(ua, va, pem, om, npx, npy, npz, ng)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, npz, ng
! Contra-variant wind components:
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: ua, va
! Pressure at edges:
    REAL, INTENT(IN) :: pem(is-1:ie+1, npz+1, js-1:je+1)
    REAL, INTENT(INOUT) :: om(isd:ied, jsd:jed, npz)
! Local:
    REAL, DIMENSION(is:ie, js:je) :: up, vp
    REAL :: v3(3, is:ie, js:je)
    REAL :: pin(isd:ied, jsd:jed)
    REAL :: pb(isd:ied, jsd:jed)
    REAL :: grad(3, is:ie, js:je)
    REAL :: pdx(3, is:ie, js:je+1)
    REAL :: pdy(3, is:ie+1, js:je)
    INTEGER :: i, j, k, n
!$omp parallel do default(shared) private(i, j, k, n, pdx, pdy, pin, pb, up, vp, grad, v3)
    DO k=1,npz
      IF (k .EQ. npz) THEN
        DO j=js,je
          DO i=is,ie
            up(i, j) = ua(i, j, npz)
            vp(i, j) = va(i, j, npz)
          END DO
        END DO
      ELSE
        DO j=js,je
          DO i=is,ie
            up(i, j) = 0.5*(ua(i, j, k)+ua(i, j, k+1))
            vp(i, j) = 0.5*(va(i, j, k)+va(i, j, k+1))
          END DO
        END DO
      END IF
! Compute Vect wind:
      DO j=js,je
        DO i=is,ie
          DO n=1,3
            v3(n, i, j) = up(i, j)*ec1(n, i, j) + vp(i, j)*ec2(n, i, j)
          END DO
        END DO
      END DO
      DO j=js-1,je+1
        DO i=is-1,ie+1
          pin(i, j) = pem(i, k+1, j)
        END DO
      END DO
! Compute pe at 4 cell corners:
      CALL A2B_ORD2(pin, pb, npx, npy, is, ie, js, je, ng)
      DO j=js,je+1
        DO i=is,ie
          DO n=1,3
            pdx(n, i, j) = (pb(i, j)+pb(i+1, j))*dx(i, j)*en1(n, i, j)
          END DO
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          DO n=1,3
            pdy(n, i, j) = (pb(i, j)+pb(i, j+1))*dy(i, j)*en2(n, i, j)
          END DO
        END DO
      END DO
! Compute grad (pe) by Green's theorem
      DO j=js,je
        DO i=is,ie
          DO n=1,3
            grad(n, i, j) = pdx(n, i, j+1) - pdx(n, i, j) - pdy(n, i, j)&
&             + pdy(n, i+1, j)
          END DO
        END DO
      END DO
! Compute inner product: V3 * grad (pe)
      DO j=js,je
        DO i=is,ie
          om(i, j, k) = om(i, j, k) + 0.5*rarea(i, j)*(v3(1, i, j)*grad(&
&           1, i, j)+v3(2, i, j)*grad(2, i, j)+v3(3, i, j)*grad(3, i, j)&
&           )
        END DO
      END DO
    END DO
  END SUBROUTINE ADV_PE

END MODULE DYN_CORE_TLM_MOD
