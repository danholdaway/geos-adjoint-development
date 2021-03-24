module dyn_core_adm_mod

  use fv_arrays_mod,      only: g_precision, p_precision, pg_precision
  use mpp_domains_mod,    only: CGRID_NE, DGRID_NE, mpp_get_boundary,   &
                                mpp_update_domains
  use fv_mp_mod,          only: domain, isd, ied, jsd, jed, is, ie, js, je, gid, mp_barrier, mp_reduce_max
  use fv_control_mod,     only: hord_mt, hord_vt, hord_tm, hord_dp, hord_ze, hord_tr, n_sponge, m_riem, m_split, &
                                hord_mt_pert, hord_vt_pert, hord_tm_pert, hord_dp_pert, hord_tr_pert, &
                                dddmp, fv_d2_bg=>d2_bg, fv_d4_bg=>d4_bg, d_ext, vtdm4, beta1, beta, init_wind_m, m_grad_p, &
                                a2b_ord, ppm_limiter, master, fv_debug, d_con, fv_nord=>nord, max_courant_no, &
                                no_cgrid, fill_dp, nwat, inline_q, breed_vortex_inline, shallow_water, spin_up_hours
  use sw_core_adm_mod,    only: c_sw, d_sw, divergence_corner, d2a2c, c_sw_adm, d_sw_adm, divergence_corner_adm, d2a2c_adm
  use a2b_edge_adm_mod,   only: a2b_ord2, a2b_ord4, a2b_ord2_adm, a2b_ord4_adm
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
public :: dyn_core, dyn_core_adm

    real   , save :: myStep = 0
    real   , save :: myTime = 0
    logical, save :: first_step = .true.

contains

  SUBROUTINE DYN_CORE_ADM(npx, npy, npz, ng, sphum, nq, bdt, n_split, &
&   zvir, cp, akap, rdgas, grav, hydrostatic, u, u_ad, v, v_ad, um, &
&   um_ad, vm, vm_ad, w, w_ad, delz, pt, pt_ad, q, q_ad, delp, delp_ad, &
&   pe, pe_ad, pk, pk_ad, phis, omga, omga_ad, ptop, pfull, ua, ua_ad, &
&   va, va_ad, uc, uc_ad, vc, vc_ad, mfx, mfx_ad, mfy, mfy_ad, cx, cx_ad&
&   , cy, cy_ad, pem, pkz, pkz_ad, peln, peln_ad, ak, bk)
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
    REAL, DIMENSION(isd:ied, jsd:jed+1, npz), INTENT(INOUT) :: um_ad
! D grid meridional wind (m/s)
    REAL, DIMENSION(isd:ied + 1, jsd:jed, npz), INTENT(INOUT) :: vm
    REAL, DIMENSION(isd:ied+1, jsd:jed, npz), INTENT(INOUT) :: vm_ad
! D grid zonal wind (m/s)
    REAL, DIMENSION(isd:ied, jsd:jed + 1, npz), INTENT(INOUT) :: u
    REAL, DIMENSION(isd:ied, jsd:jed+1, npz), INTENT(INOUT) :: u_ad
! D grid meridional wind (m/s)
    REAL, DIMENSION(isd:ied + 1, jsd:jed, npz), INTENT(INOUT) :: v
    REAL, DIMENSION(isd:ied+1, jsd:jed, npz), INTENT(INOUT) :: v_ad
! vertical vel. (m/s)
    REAL, INTENT(INOUT) :: w(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: w_ad(isd:ied, jsd:jed, npz)
! delta-height (m)
    REAL, INTENT(INOUT) :: delz(is:ie, js:je, npz)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: pt_ad(isd:ied, jsd:jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp_ad(isd:ied, jsd:jed, npz)
!
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, npz, nq)
    REAL, INTENT(INOUT) :: q_ad(isd:ied, jsd:jed, npz, nq)
!    real, intent(IN), optional:: time_total  ! total time (seconds) since start
!THESE ARE ACTUALLY IN THE MODULE
    REAL :: delzc(is:ie, js:je, npz)
    REAL :: ut(isd:ied, jsd:jed, npz)
    REAL :: ut_ad(isd:ied, jsd:jed, npz)
    REAL :: vt(isd:ied, jsd:jed, npz)
    REAL :: vt_ad(isd:ied, jsd:jed, npz)
    REAL :: crx(is:ie+1, jsd:jed, npz), xfx(is:ie+1, jsd:jed, npz)
    REAL :: crx_ad(is:ie+1, jsd:jed, npz), xfx_ad(is:ie+1, jsd:jed, npz)
    REAL :: cry(isd:ied, js:je+1, npz), yfx(isd:ied, js:je+1, npz)
    REAL :: cry_ad(isd:ied, js:je+1, npz), yfx_ad(isd:ied, js:je+1, npz)
    REAL :: divg_d(isd:ied+1, jsd:jed+1, npz)
    REAL :: divg_d_ad(isd:ied+1, jsd:jed+1, npz)
    REAL :: zh(isd:ied, jsd:jed, npz)
    REAL(pg_precision) :: du(isd:ied, jsd:jed+1, npz), dv(isd:ied+1, jsd&
&   :jed, npz)
    REAL(pg_precision) :: du_ad(isd:ied, jsd:jed+1, npz), dv_ad(isd:ied+&
&   1, jsd:jed, npz)
    REAL(pg_precision) :: pkc(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision) :: pkc_ad(isd:ied, jsd:jed, npz+1)
    REAL :: pkd(isd:ied, jsd:jed, npz+1)
    REAL :: pkd_ad(isd:ied, jsd:jed, npz+1)
    REAL :: delpc(isd:ied, jsd:jed, npz)
    REAL :: delpc_ad(isd:ied, jsd:jed, npz)
    REAL(pg_precision) :: pk3(isd:ied, jsd:jed, npz+1)
    REAL :: ptc(isd:ied, jsd:jed, npz)
    REAL :: ptc_ad(isd:ied, jsd:jed, npz)
    REAL(pg_precision) :: gz(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision) :: gz_ad(isd:ied, jsd:jed, npz+1)
!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
! Surface geopotential (g*Z_surf)
    REAL, INTENT(INOUT) :: phis(isd:ied, jsd:jed)
! edge pressure (pascal)
    REAL(p_precision), INTENT(INOUT) :: pe(is-1:ie+1, npz+1, js-1:je+1)
    REAL(p_precision) :: pe_ad(is-1:ie+1, npz+1, js-1:je+1)
    REAL, INTENT(OUT) :: pem(is-1:ie+1, npz+1, js-1:je+1)
! ln(pe)
    REAL(p_precision) :: peln(is:ie, npz+1, js:je)
    REAL(p_precision) :: peln_ad(is:ie, npz+1, js:je)
! pe**kappa
    REAL(p_precision), INTENT(INOUT) :: pk(is:ie, js:je, npz+1)
    REAL(p_precision) :: pk_ad(is:ie, js:je, npz+1)
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
! Vertical pressure velocity (pa/s)
    REAL, INTENT(INOUT) :: omga(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: omga_ad(isd:ied, jsd:jed, npz)
! (uc, vc) are mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: uc_ad(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: vc(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: vc_ad(isd:ied, jsd:jed+1, npz)
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua_ad
! The Flux capacitors: accumulated Mass flux arrays
    REAL, INTENT(INOUT) :: mfx(is:ie+1, js:je, npz)
    REAL, INTENT(INOUT) :: mfx_ad(is:ie+1, js:je, npz)
    REAL, INTENT(INOUT) :: mfy(is:ie, js:je+1, npz)
    REAL, INTENT(INOUT) :: mfy_ad(is:ie, js:je+1, npz)
! Accumulated Courant number arrays
    REAL, INTENT(INOUT) :: cx(is:ie+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: cx_ad(is:ie+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: cy(isd:ied, js:je+1, npz)
    REAL, INTENT(INOUT) :: cy_ad(isd:ied, js:je+1, npz)
! Work:
! 
    REAL(p_precision), INTENT(INOUT) :: pkz(is:ie, js:je, npz)
    REAL(p_precision) :: pkz_ad(is:ie, js:je, npz)
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
    REAL(pg_precision) :: divg2_ad(is:ie+1, js:je+1)
    REAL :: wk(isd:ied, jsd:jed)
    REAL :: wk_ad(isd:ied, jsd:jed)
! --- For no_cgrid option ---
    REAL :: u2(isd:ied, jsd:jed+1)
    REAL :: u2_ad(isd:ied, jsd:jed+1)
    REAL :: v2(isd:ied+1, jsd:jed)
    REAL :: v2_ad(isd:ied+1, jsd:jed)
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
    INTEGER :: arg1
    LOGICAL :: arg10
    INTEGER :: arg2
    INTEGER :: branch
    INTEGER :: ad_to
    INTEGER :: ad_to0
    INTEGER :: ad_to1
    INTEGER :: ad_to2
    INTEGER :: ad_to3
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: va_ad
    REAL(pg_precision) :: temp3
    REAL(pg_precision) :: temp2
    REAL(pg_precision) :: temp1
    REAL(pg_precision) :: temp0
    REAL(pg_precision) :: temp_ad5
    REAL(pg_precision) :: temp_ad4
    REAL(pg_precision) :: temp_ad3
    REAL :: temp_ad2
    REAL(pg_precision) :: temp_ad1
    REAL :: temp_ad0
    REAL(pg_precision) :: temp_ad
    REAL :: temp
    REAL(pg_precision) :: temp8
    REAL :: y4
    REAL(pg_precision) :: temp7
    REAL :: y3
    REAL(pg_precision) :: temp6
    REAL :: y2
    REAL(pg_precision) :: temp5
    REAL :: y1
    REAL :: temp4
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
! Indexes:
    ism1 = is - 1
    iep1 = ie + 1
    jsm1 = js - 1
    jep1 = je + 1
!call timing_on('COMM_TOTAL')
    IF (npz .GT. 1) THEN
!oncall mpp_update_domains( pt,   domain, complete=.false. )
      CALL PUSHREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1)*npz)
      CALL mpp_update_domains( pt,   domain, complete=.false. )
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
!oncall mpp_update_domains( delp, domain, complete=.true. )
    CALL PUSHREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1)*npz)
    CALL mpp_update_domains( delp, domain, complete=.true. )
!oncall mpp_update_domains( u, v, domain, gridtype=DGRID_NE, complete=.true. )
    arg1 = jed + 1
    CALL PUSHREAL8ARRAY(u, (ied-isd+1)*(jed-jsd+2)*npz)
    CALL PUSHREAL8ARRAY(v, (ied-isd+2)*(jed-jsd+1)*npz)
    call mpp_update_domains( u, v, domain, gridtype=DGRID_NE, complete=.true. )
!call timing_off('COMM_TOTAL')
!Initiliaze vars
    gz = 0.0
    ptc = 0.0
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
      CALL PUSHREAL8ARRAY(pkz, (ie-is+1)*(je-js+1)*npz)
      CALL PUSHREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*(npz&
&                   +1)/8)
      CALL PUSHREAL8ARRAY(peln, (ie-is+1)*(npz+1)*(je-js+1))
      CALL PUSHREAL8ARRAY(pe, (ie-is+3)*(npz+1)*(je-js+3))
      CALL GEOPK(ptop, pe, peln, delp, pkc, gz, phis, pt, pkz, npz, akap&
&          , .false.)
      CALL PUSHREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1)*npz)
      CALL PUSHREAL8ARRAY(gz, pg_precision*(ied-isd+1)*(jed-jsd+1)*(npz+&
&                   1)/8)
      CALL PUSHREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*(npz&
&                   +1)/8)
      CALL GRAD1_P(du, dv, pkc, gz, delp, dt, ng, npx, npy, npz, ptop, &
&            ptk, hydrostatic)
      IF (init_wind_m) THEN
!oncall mpp_update_domains(du, dv, domain, gridtype=DGRID_NE)
        call mpp_update_domains(du, dv, domain, gridtype=DGRID_NE)
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
! Empty the "flux capacitors"
    mfx(:, :, :) = 0.
    mfy(:, :, :) = 0.
    cx(:, :, :) = 0.
    cy(:, :, :) = 0.
    elapsed_dt = 0.0
    subcycle = 1
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
              crx(i, j, k) = v(i, j, k)*ndt*rdy(i, j)
              IF (crx(i, j, k) .GE. 0.) THEN
                y1 = crx(i, j, k)
              ELSE
                y1 = -crx(i, j, k)
              END IF
              IF (cmax .LT. y1) THEN
                CALL PUSHCONTROL1B(0)
                cmax = y1
              ELSE
                CALL PUSHCONTROL1B(1)
                cmax = cmax
              END IF
            END DO
          END DO
          DO j=js,je+1
            DO i=is,ie
! D-Grid zonal courant number stored in cry
              cry(i, j, k) = u(i, j, k)*ndt*rdx(i, j)
              IF (cry(i, j, k) .GE. 0.) THEN
                y2 = cry(i, j, k)
              ELSE
                y2 = -cry(i, j, k)
              END IF
              IF (cmax .LT. y2) THEN
                CALL PUSHCONTROL1B(0)
                cmax = y2
              ELSE
                CALL PUSHCONTROL1B(1)
                cmax = cmax
              END IF
            END DO
          END DO
        END DO
!oncall mp_reduce_max(cmax)
!oncall mp_barrier()
        call mp_reduce_max(cmax)
        call mp_barrier()
!NINT2CEILING
        subt = CEILING(cmax/max_courant_no)
        IF (subcycle .NE. subt) THEN
          subcycle = subt
          CALL PUSHREAL8(dt)
          dt = (bdt-elapsed_dt)/REAL((n_split-(it-1))*subcycle)
          CALL PUSHREAL8(dt2)
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
          CALL PUSHCONTROL2B(2)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(0)
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
          CALL PUSHCONTROL1B(0)
          do_omega = .false.
!onif (test_case==9) call case9_forcing1(phis, time_total)
        ELSE IF (it .EQ. n_split .AND. subt .EQ. subcycle) THEN
          CALL PUSHCONTROL1B(1)
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
          CALL PUSHCONTROL1B(1)
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
                  u2(i, j) = u(i, j, k) + 0.5*du(i, j, k)
                END DO
              END DO
              DO j=jsd,jed
                DO i=isd,ied+1
                  v2(i, j) = v(i, j, k) + 0.5*dv(i, j, k)
                END DO
              END DO
              CALL PUSHCONTROL1B(0)
            ELSE
              DO j=jsd,jed+1
                DO i=isd,ied
                  u2(i, j) = 1.5*u(i, j, k) - 0.5*um(i, j, k)
                END DO
              END DO
              DO j=jsd,jed
                DO i=isd,ied+1
                  v2(i, j) = 1.5*v(i, j, k) - 0.5*vm(i, j, k)
                END DO
              END DO
              CALL PUSHCONTROL1B(1)
            END IF
            arg10 = nord .GT. 0
            CALL PUSHREAL8ARRAY(va(isd, jsd, k), ied - isd + 1)
            CALL PUSHREAL8ARRAY(ua(isd, jsd, k), ied - isd + 1)
            CALL D2A2C(u(isd, jsd, k), v(isd, jsd, k), u2, v2, ua(isd, &
&                jsd, k), va(isd, jsd, k), uc(isd, jsd, k), vc(isd, jsd&
&                , k), arg10)
          END DO
!if ( .not. hydrostatic ) delpc(:,:,:) = delp(:,:,:)
          um(:, :, :) = u(:, :, :)
          vm(:, :, :) = v(:, :, :)
!call timing_off('no_cgrid')
          CALL PUSHCONTROL1B(0)
        ELSE
!call timing_on('c_sw')
!$omp parallel do default(shared) private(i,j,k)
          DO k=1,npz
            arg10 = nord .GT. 0
            CALL PUSHREAL8ARRAY(vt(isd, jsd, k), ied - isd + 1)
            CALL PUSHREAL8ARRAY(ut(isd, jsd, k), ied - isd + 1)
            CALL PUSHREAL8ARRAY(va(isd, jsd, k), ied - isd + 1)
            CALL PUSHREAL8ARRAY(ua(isd, jsd, k), ied - isd + 1)
            CALL PUSHREAL8ARRAY(vc(isd, jsd, k), ied - isd + 1)
            CALL PUSHREAL8ARRAY(uc(isd, jsd, k), ied - isd + 2)
            CALL PUSHREAL8ARRAY(pt(isd, jsd, k), ied - isd + 1)
            CALL PUSHREAL8ARRAY(ptc(isd, jsd, k), ied - isd + 1)
            CALL PUSHREAL8ARRAY(delp(isd, jsd, k), ied - isd + 1)
            CALL PUSHREAL8ARRAY(delpc(isd, jsd, k), ied - isd + 1)
            CALL C_SW(delpc(isd, jsd, k), delp(isd, jsd, k), ptc(isd, &
&               jsd, k), pt(isd, jsd, k), u(isd, jsd, k), v(isd, jsd, k)&
&               , w(isd, jsd, k), uc(isd, jsd, k), vc(isd, jsd, k), ua(&
&               isd, jsd, k), va(isd, jsd, k), omga(isd, jsd, k), ut(isd&
&               , jsd, k), vt(isd, jsd, k), dt2, hydrostatic, arg10)
          END DO
! on output omga is updated w
!call timing_off('c_sw')
          IF (fill_dp) THEN
            CALL PUSHREAL8ARRAY(ptc, (ied-isd+1)*(jed-jsd+1)*npz)
            CALL PUSHREAL8ARRAY(delpc, (ied-isd+1)*(jed-jsd+1)*npz)
            CALL PUSHREAL8ARRAY(omga, (ied-isd+1)*(jed-jsd+1)*npz)
            CALL MIX_DP(hydrostatic, omga, delpc, ptc, npz, ak, bk, &
&                 .true.)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!      if ( hydrostatic ) then
          IF (beta1 .GT. 0.001) THEN
            CALL PUSHREAL8(alpha)
            alpha = 1. - beta1
            CALL PUSHREAL8ARRAY(delpc, (ied-isd+1)*(jed-jsd+1)*npz)
            delpc(:, :, :) = beta1*delp(:, :, :) + alpha*delpc(:, :, :)
            CALL PUSHREAL8ARRAY(ptc, (ied-isd+1)*(jed-jsd+1)*npz)
            ptc(:, :, :) = beta1*pt(:, :, :) + alpha*ptc(:, :, :)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8ARRAY(pkz, (ie-is+1)*(je-js+1)*npz)
          CALL PUSHREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*&
&                       (npz+1)/8)
          CALL PUSHREAL8ARRAY(peln, (ie-is+1)*(npz+1)*(je-js+1))
          CALL PUSHREAL8ARRAY(pe, (ie-is+3)*(npz+1)*(je-js+3))
          CALL GEOPK(ptop, pe, peln, delpc, pkc, gz, phis, ptc, pkz, npz&
&              , akap, .true.)
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
                CALL PUSHREAL8(wk(i, j))
                wk(i, j) = pkc(i, j, k+1) - pkc(i, j, k)
              END DO
              CALL PUSHINTEGER4(i - 1)
            END DO
            CALL PUSHINTEGER4(j - 1)
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
                uc(i, j, k) = uc(i, j, k) + dt2*rdxc(i, j)/(wk(i-1, j)+&
&                 wk(i, j))*((gz(i-1, j, k+1)-gz(i, j, k))*(pkc(i, j, k+&
&                 1)-pkc(i-1, j, k))+(gz(i-1, j, k)-gz(i, j, k+1))*(pkc(&
&                 i-1, j, k+1)-pkc(i, j, k)))
              END DO
              CALL PUSHINTEGER4(i - 1)
            END DO
            DO j=js,jeb1
              DO i=is,ie
                vc(i, j, k) = vc(i, j, k) + dt2*rdyc(i, j)/(wk(i, j-1)+&
&                 wk(i, j))*((gz(i, j-1, k+1)-gz(i, j, k))*(pkc(i, j, k+&
&                 1)-pkc(i, j-1, k))+(gz(i, j-1, k)-gz(i, j, k+1))*(pkc(&
&                 i, j-1, k+1)-pkc(i, j, k)))
              END DO
            END DO
            CALL PUSHINTEGER4(j - 1)
          END DO
          CALL PUSHCONTROL1B(1)
        END IF
!call timing_on('COMM_TOTAL')
!oncall mpp_update_domains(uc, vc, domain, gridtype=CGRID_NE, complete=.true.)
        call mpp_update_domains(uc, vc, domain, gridtype=CGRID_NE, complete=.true.)
!onif (test_case==9) call case9_forcing2(phis)
        IF (inline_q) THEN
!call timing_on('COMM_TOTAL')
!oncall mpp_update_domains( q,  domain, complete=.true. )
          CALL PUSHREAL8ARRAY(q, (ied-isd+1)*(jed-jsd+1)*npz*nq)
          call mpp_update_domains( q,  domain, complete=.true. )
!call timing_off('COMM_TOTAL')
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (nord .GT. 0) THEN
          CALL DIVERGENCE_CORNER(u, v, ua, va, divg_d, npz)
!call timing_on('COMM_TOTAL')
!oncall mpp_update_domains(divg_d, domain, position=CORNER)
          call mpp_update_domains(divg_d, domain, position=CORNER)
!call timing_off('COMM_TOTAL')
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
!call timing_on('d_sw')
!$omp parallel do default(shared) private(i, j, k, nord_k, damp_k, d2_divg, dd_divg, hord_m, hord_v, hord_t, hord_p, wk)
        DO k=1,npz
          CALL PUSHINTEGER4(hord_m)
          hord_m = hord_mt
          CALL PUSHINTEGER4(hord_t)
          hord_t = hord_tm
          CALL PUSHINTEGER4(hord_v)
          hord_v = hord_vt
          CALL PUSHINTEGER4(hord_p)
          hord_p = hord_dp
          CALL PUSHINTEGER4(hord_m_pert)
          hord_m_pert = hord_mt_pert
          CALL PUSHINTEGER4(hord_t_pert)
          hord_t_pert = hord_tm_pert
          CALL PUSHINTEGER4(hord_v_pert)
          hord_v_pert = hord_vt_pert
          CALL PUSHINTEGER4(hord_p_pert)
          hord_p_pert = hord_dp_pert
          CALL PUSHINTEGER4(nord_k)
          nord_k = nord
          CALL PUSHREAL8(damp_k)
          damp_k = dddmp
          y4 = d2_bg*(1.-3.*TANH(0.1*LOG(pfull(k)/pfull(npz))))
          IF (0.20 .GT. y4) THEN
            CALL PUSHREAL8(d2_divg)
            d2_divg = y4
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL8(d2_divg)
            d2_divg = 0.20
            CALL PUSHCONTROL1B(1)
          END IF
          IF (n_sponge .EQ. -1 .OR. npz .EQ. 1) THEN
            CALL PUSHCONTROL3B(2)
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
              CALL PUSHCONTROL3B(3)
              d2_divg = d2_divg
            ELSE
              CALL PUSHCONTROL3B(3)
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
              CALL PUSHCONTROL3B(4)
              d2_divg = d2_divg
            ELSE
              CALL PUSHCONTROL3B(4)
              d2_divg = 0.15
            END IF
          ELSE IF (k .EQ. n_sponge + 1 .AND. npz .GT. 24) THEN
            hord_v = 2
            hord_t = 2
            hord_p = 2
            hord_v_pert = 2
            hord_t_pert = 2
            hord_p_pert = 2
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
              CALL PUSHCONTROL3B(1)
              damp_k = 0.
            ELSE
              CALL PUSHCONTROL3B(1)
              damp_k = 0.12
            END IF
          ELSE
            CALL PUSHCONTROL3B(0)
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
          IF (d_ext .GT. 0.) THEN
            CALL PUSHREAL8ARRAY(wk, (ied-isd+1)*(jed-jsd+1))
            CALL A2B_ORD2(delp(:, :, k), wk, npx, npy, is, ie, js, je, &
&                   ng, .false.)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8ARRAY(q, (ied-isd+1)*(jed-jsd+1)*npz*nq)
          CALL PUSHREAL8ARRAY(yfx(isd, js, k), ied - isd + 1)
          CALL PUSHREAL8ARRAY(xfx(is, jsd, k), ie - is + 2)
          CALL PUSHREAL8ARRAY(cry(isd, js, k), ied - isd + 1)
          CALL PUSHREAL8ARRAY(crx(is, jsd, k), ie - is + 2)
          CALL PUSHREAL8ARRAY(divg_d(isd, jsd, k), ied - isd + 2)
          CALL PUSHREAL8ARRAY(vc(isd, jsd, k), ied - isd + 1)
          CALL PUSHREAL8ARRAY(uc(isd, jsd, k), ied - isd + 2)
          CALL PUSHREAL8ARRAY(w(isd, jsd, k), ied - isd + 1)
          CALL PUSHREAL8ARRAY(v(isd, jsd, k), ied - isd + 2)
          CALL PUSHREAL8ARRAY(u(isd, jsd, k), ied - isd + 1)
          CALL PUSHREAL8ARRAY(pt(isd, jsd, k), ied - isd + 1)
          CALL PUSHREAL8ARRAY(ptc(isd, jsd, k), ied - isd + 1)
          CALL PUSHREAL8ARRAY(delp(isd, jsd, k), ied - isd + 1)
          CALL PUSHREAL8ARRAY(pkd(isd, jsd, k), ied - isd + 1)
          CALL D_SW(pkd(isd, jsd, k), delp(isd, jsd, k), ptc(isd, jsd, k&
&             ), pt(isd, jsd, k), u(isd, jsd, k), v(isd, jsd, k), w(isd&
&             , jsd, k), uc(isd, jsd, k), vc(isd, jsd, k), ua(isd, jsd, &
&             k), va(isd, jsd, k), divg_d(isd, jsd, k), mfx(is, js, k), &
&             mfy(is, js, k), cx(is, jsd, k), cy(isd, js, k), crx(is, &
&             jsd, k), cry(isd, js, k), xfx(is, jsd, k), yfx(isd, js, k)&
&             , zvir, sphum, nq, q, k, npz, inline_q, pkz(is, js, k), dt&
&             , hord_tr, hord_m, hord_v, hord_t, hord_p, nord_k, damp_k&
&             , d2_divg, dd_divg, vtdm4, d_con, hydrostatic, ppm_limiter&
&            )
          IF (d_ext .GT. 0.) THEN
            DO j=js,jep1
              DO i=is,iep1
! delp at cell corners
                CALL PUSHREAL8(ptc(i, j, k))
                ptc(i, j, k) = wk(i, j)
              END DO
            END DO
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
!call timing_off('d_sw')
        IF (fill_dp) THEN
          CALL PUSHREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1)*npz)
          CALL PUSHREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1)*npz)
          CALL PUSHREAL8ARRAY(w, (ied-isd+1)*(jed-jsd+1)*npz)
          CALL MIX_DP(hydrostatic, w, delp, pt, npz, ak, bk, .false.)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!    if ( fv_debug ) then
!         call prt_maxmin('DELP', delp, is, ie, js, je, ng, npz, 1.E-2, master)
!    endif
        IF (d_ext .GT. 0.) THEN
          CALL PUSHREAL8(d2_divg)
          d2_divg = d_ext*da_min_c
! pkd() is 3D field of horizontal divergence
! ptc is "delp" at cell corners
!$omp parallel do default(shared) private(i,j,k)
          DO j=js,jep1
            DO i=is,iep1
              CALL PUSHREAL8(wk(i, j))
              wk(i, j) = ptc(i, j, 1)
              divg2(i, j) = wk(i, j)*pkd(i, j, 1)
            END DO
            DO k=2,npz
              DO i=is,iep1
                CALL PUSHREAL8(wk(i, j))
                wk(i, j) = wk(i, j) + ptc(i, j, k)
                divg2(i, j) = divg2(i, j) + ptc(i, j, k)*pkd(i, j, k)
              END DO
            END DO
            DO i=is,iep1
              CALL PUSHREAL8ARRAY(divg2(i, j), pg_precision/8)
              divg2(i, j) = d2_divg*divg2(i, j)/wk(i, j)
            END DO
          END DO
          CALL PUSHCONTROL1B(0)
        ELSE
          divg2 = 0.
          CALL PUSHCONTROL1B(1)
        END IF
!call timing_on('COMM_TOTAL')
!oncall mpp_update_domains(  pt, domain, complete=.false.)
        CALL PUSHREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1)*npz)
        call mpp_update_domains(  pt, domain, complete=.false.)
!oncall mpp_update_domains(delp, domain, complete=.true.)
        call mpp_update_domains(delp, domain, complete=.true.)
!call timing_off('COMM_TOTAL')
!     if ( hydrostatic ) then
        CALL PUSHREAL8ARRAY(pkz, (ie-is+1)*(je-js+1)*npz)
        CALL PUSHREAL8ARRAY(gz, pg_precision*(ied-isd+1)*(jed-jsd+1)*(&
&                     npz+1)/8)
        CALL PUSHREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*(&
&                     npz+1)/8)
        CALL PUSHREAL8ARRAY(peln, (ie-is+1)*(npz+1)*(je-js+1))
        CALL PUSHREAL8ARRAY(pe, (ie-is+3)*(npz+1)*(je-js+3))
        CALL GEOPK(ptop, pe, peln, delp, pkc, gz, phis, pt, pkz, npz, &
&            akap, .false.)
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
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (do_omega) THEN
            CALL PUSHCONTROL1B(0)
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
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!      if ( .not.hydrostatic .and. m_grad_p == 0 ) then
!           call two_grad_p(u, v, pkc, gz, delp, pk3, divg2, dt, ng, npx, npy,   &
!                           npz, ptk)  
!      else
        IF (beta .GT. 1.e-4) THEN
          DO k=1,npz
            DO j=js,je+1
              DO i=is,ie
                CALL PUSHREAL8(u(i, j, k))
                u(i, j, k) = (u(i, j, k)+divg2(i, j)-divg2(i+1, j))*rdx(&
&                 i, j) + beta*du(i, j, k)
              END DO
            END DO
          END DO
          DO k=1,npz
            DO j=js,je
              DO i=is,ie+1
                CALL PUSHREAL8(v(i, j, k))
                v(i, j, k) = (v(i, j, k)+divg2(i, j)-divg2(i, j+1))*rdy(&
&                 i, j) + beta*dv(i, j, k)
              END DO
            END DO
          END DO
          CALL PUSHREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1)*npz)
          CALL PUSHREAL8ARRAY(gz, pg_precision*(ied-isd+1)*(jed-jsd+1)*(&
&                       npz+1)/8)
          CALL PUSHREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*&
&                       (npz+1)/8)
          CALL GRAD1_P(du, dv, pkc, gz, delp, dt, ng, npx, npy, npz, &
&                ptop, ptk, hydrostatic)
          CALL PUSHREAL8(alpha)
          alpha = 1. - beta
          DO k=1,npz
            DO j=js,je+1
              DO i=is,ie
                CALL PUSHREAL8(u(i, j, k))
                u(i, j, k) = u(i, j, k) + alpha*du(i, j, k)
              END DO
            END DO
          END DO
          DO k=1,npz
            DO j=js,je
              DO i=is,ie+1
                CALL PUSHREAL8(v(i, j, k))
                v(i, j, k) = v(i, j, k) + alpha*dv(i, j, k)
              END DO
            END DO
          END DO
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1)*npz)
          CALL PUSHREAL8ARRAY(gz, pg_precision*(ied-isd+1)*(jed-jsd+1)*(&
&                       npz+1)/8)
          CALL PUSHREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*&
&                       (npz+1)/8)
          CALL PUSHREAL8ARRAY(v, (ied-isd+2)*(jed-jsd+1)*npz)
          CALL PUSHREAL8ARRAY(u, (ied-isd+1)*(jed-jsd+2)*npz)
          CALL ONE_GRAD_P(u, v, pkc, gz, divg2, delp, dt, ng, npx, npy, &
&                   npz, ptop, ptk, hydrostatic)
          CALL PUSHCONTROL1B(1)
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
!oncall mpp_get_boundary(u, v, domain, wbuffery=wbuffer, ebuffery=ebuffer,  &
!                     sbufferx=sbuffer, nbufferx=nbuffer, gridtype=DGRID_NE )
          !arg1 = jed + 1
          CALL PUSHREAL8ARRAY(u, (ied-isd+1)*(jed-jsd+2)*npz)
          !CALL MPP_GET_BOUNDARY_DUMMY3(u, isd, ied, jsd, arg1, npz)
          !nbuffer = 1.0
          !arg1 = ied + 1
          CALL PUSHREAL8ARRAY(v, (ied-isd+2)*(jed-jsd+1)*npz)
          !CALL MPP_GET_BOUNDARY_DUMMY3(v, isd, arg1, jsd, jed, npz)
          call mpp_get_boundary(u, v, domain, wbuffery=wbuffer, ebuffery=ebuffer, &
                                sbufferx=sbuffer, nbufferx=nbuffer, gridtype=DGRID_NE )
          !ebuffer = 1.0
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
          CALL PUSHCONTROL1B(0)
        ELSE
!oncall mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
          !arg1 = jed + 1
          CALL PUSHREAL8ARRAY(u, (ied-isd+1)*(jed-jsd+2)*npz)
          !CALL MPP_UPDATE_DOMAINS_DUMMY3(u, isd, ied, jsd, arg1, npz)
          !arg1 = ied + 1
          CALL PUSHREAL8ARRAY(v, (ied-isd+2)*(jed-jsd+1)*npz)
          call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
          CALL PUSHCONTROL1B(1)
        END IF
!call timing_off('COMM_TOTAL')
!init_wind_m = .false. !dh always false
!first_step = .false.
        elapsed_dt = elapsed_dt + dt
      END DO
      CALL PUSHINTEGER4(subt - 1)
    END DO
    uc_ad = 0.0_8
    um_ad = 0.0_8
    va_ad = 0.0_8
    vc_ad = 0.0_8
    vm_ad = 0.0_8
    xfx_ad = 0.0_8
    v2_ad = 0.0_8
    gz_ad = 0.0_8
    du_ad = 0.0_8
    dv_ad = 0.0_8
    ptc_ad = 0.0_8
    ut_ad = 0.0_8
    pkc_ad = 0.0_8
    pkd_ad = 0.0_8
    delpc_ad = 0.0_8
    yfx_ad = 0.0_8
    vt_ad = 0.0_8
    divg_d_ad = 0.0_8
    u2_ad = 0.0_8
    crx_ad = 0.0_8
    cry_ad = 0.0_8
    wk_ad = 0.0_8
    divg2_ad = 0.0_8
    DO it=n_split,1,-1
      CALL POPINTEGER4(ad_to3)
      DO subt=ad_to3,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          call mpp_get_boundary(u_ad, v_ad, domain, wbuffery=wbuffer, ebuffery=ebuffer,  &
                                sbufferx=sbuffer, nbufferx=nbuffer, gridtype=DGRID_NE )
          DO k=1,npz
            DO i=is,ie
              u_ad(i, je+1, k) = nbuffer(i, k)
            END DO
          END DO
          DO k=1,npz
            DO j=js,je
              v_ad(ie+1, j, k) = ebuffer(j, k)
            END DO
          END DO
          CALL POPREAL8ARRAY(v, (ied-isd+2)*(jed-jsd+1)*npz)
          CALL POPREAL8ARRAY(u, (ied-isd+1)*(jed-jsd+2)*npz)
        ELSE
          CALL POPREAL8ARRAY(v, (ied-isd+2)*(jed-jsd+1)*npz)
          !CALL MPP_UPDATE_DOMAINS_DUMMY3_ADM(v, v_ad, isd, arg1, jsd,jed, npz)
          CALL POPREAL8ARRAY(u, (ied-isd+1)*(jed-jsd+2)*npz)
          call mpp_update_domains(u_ad, v_ad, domain, gridtype=DGRID_NE)
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO k=npz,1,-1
            DO j=je,js,-1
              DO i=ie+1,is,-1
                CALL POPREAL8(v(i, j, k))
                dv_ad(i, j, k) = dv_ad(i, j, k) + alpha*v_ad(i, j, k)
              END DO
            END DO
          END DO
          DO k=npz,1,-1
            DO j=je+1,js,-1
              DO i=ie,is,-1
                CALL POPREAL8(u(i, j, k))
                du_ad(i, j, k) = du_ad(i, j, k) + alpha*u_ad(i, j, k)
              END DO
            END DO
          END DO
          CALL POPREAL8(alpha)
          CALL POPREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*(&
&                      npz+1)/8)
          CALL POPREAL8ARRAY(gz, pg_precision*(ied-isd+1)*(jed-jsd+1)*(&
&                      npz+1)/8)
          CALL POPREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1)*npz)
          CALL GRAD1_P_ADM(du, du_ad, dv, dv_ad, pkc, pkc_ad, gz, gz_ad&
&                    , delp, delp_ad, dt, ng, npx, npy, npz, ptop, ptk, &
&                    hydrostatic)
          DO k=npz,1,-1
            DO j=je,js,-1
              DO i=ie+1,is,-1
                CALL POPREAL8(v(i, j, k))
                temp_ad5 = rdy(i, j)*v_ad(i, j, k)
                divg2_ad(i, j) = divg2_ad(i, j) + temp_ad5
                divg2_ad(i, j+1) = divg2_ad(i, j+1) - temp_ad5
                dv_ad(i, j, k) = dv_ad(i, j, k) + beta*v_ad(i, j, k)
                v_ad(i, j, k) = temp_ad5
              END DO
            END DO
          END DO
          DO k=npz,1,-1
            DO j=je+1,js,-1
              DO i=ie,is,-1
                CALL POPREAL8(u(i, j, k))
                temp_ad4 = rdx(i, j)*u_ad(i, j, k)
                divg2_ad(i, j) = divg2_ad(i, j) + temp_ad4
                divg2_ad(i+1, j) = divg2_ad(i+1, j) - temp_ad4
                du_ad(i, j, k) = du_ad(i, j, k) + beta*u_ad(i, j, k)
                u_ad(i, j, k) = temp_ad4
              END DO
            END DO
          END DO
        ELSE
          CALL POPREAL8ARRAY(u, (ied-isd+1)*(jed-jsd+2)*npz)
          CALL POPREAL8ARRAY(v, (ied-isd+2)*(jed-jsd+1)*npz)
          CALL POPREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*(&
&                      npz+1)/8)
          CALL POPREAL8ARRAY(gz, pg_precision*(ied-isd+1)*(jed-jsd+1)*(&
&                      npz+1)/8)
          CALL POPREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1)*npz)
          CALL ONE_GRAD_P_ADM(u, u_ad, v, v_ad, pkc, pkc_ad, gz, gz_ad, &
&                       divg2, divg2_ad, delp, delp_ad, dt, ng, npx, npy&
&                       , npz, ptop, ptk, hydrostatic)
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            DO k=npz+1,1,-1
              DO j=je,js,-1
                DO i=ie,is,-1
                  pkc_ad(i, j, k) = pkc_ad(i, j, k) + pk_ad(i, j, k)
                  pk_ad(i, j, k) = 0.0_8
                END DO
              END DO
            END DO
          END IF
        END IF
        CALL POPREAL8ARRAY(pe, (ie-is+3)*(npz+1)*(je-js+3))
        CALL POPREAL8ARRAY(peln, (ie-is+1)*(npz+1)*(je-js+1))
        CALL POPREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*(&
&                    npz+1)/8)
        CALL POPREAL8ARRAY(gz, pg_precision*(ied-isd+1)*(jed-jsd+1)*(npz&
&                    +1)/8)
        CALL POPREAL8ARRAY(pkz, (ie-is+1)*(je-js+1)*npz)
        CALL GEOPK_ADM(ptop, pe, pe_ad, peln, peln_ad, delp, delp_ad, &
&                pkc, pkc_ad, gz, gz_ad, phis, pt, pt_ad, pkz, pkz_ad, &
&                npz, akap, .false.)
        call mpp_update_domains(delp_ad, domain, complete=.true.)
        CALL POPREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1)*npz)
        call mpp_update_domains(  pt_ad, domain, complete=.false.)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO j=jep1,js,-1
            DO i=iep1,is,-1
              CALL POPREAL8ARRAY(divg2(i, j), pg_precision/8)
              temp_ad3 = d2_divg*divg2_ad(i, j)/wk(i, j)
              wk_ad(i, j) = wk_ad(i, j) - divg2(i, j)*temp_ad3/wk(i, j)
              divg2_ad(i, j) = temp_ad3
            END DO
            DO k=npz,2,-1
              DO i=iep1,is,-1
                ptc_ad(i, j, k) = ptc_ad(i, j, k) + wk_ad(i, j) + pkd(i&
&                 , j, k)*divg2_ad(i, j)
                pkd_ad(i, j, k) = pkd_ad(i, j, k) + ptc(i, j, k)*&
&                 divg2_ad(i, j)
                CALL POPREAL8(wk(i, j))
              END DO
            END DO
            DO i=iep1,is,-1
              wk_ad(i, j) = wk_ad(i, j) + pkd(i, j, 1)*divg2_ad(i, j)
              pkd_ad(i, j, 1) = pkd_ad(i, j, 1) + wk(i, j)*divg2_ad(i, j&
&               )
              divg2_ad(i, j) = 0.0_8
              CALL POPREAL8(wk(i, j))
              ptc_ad(i, j, 1) = ptc_ad(i, j, 1) + wk_ad(i, j)
              wk_ad(i, j) = 0.0_8
            END DO
          END DO
          CALL POPREAL8(d2_divg)
        ELSE
          divg2_ad = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8ARRAY(w, (ied-isd+1)*(jed-jsd+1)*npz)
          CALL POPREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1)*npz)
          CALL POPREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1)*npz)
          CALL MIX_DP_ADM(hydrostatic, w, w_ad, delp, delp_ad, pt, pt_ad&
&                   , npz, ak, bk, .false.)
        END IF
        DO k=npz,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            DO j=jep1,js,-1
              DO i=iep1,is,-1
                CALL POPREAL8(ptc(i, j, k))
                wk_ad(i, j) = wk_ad(i, j) + ptc_ad(i, j, k)
                ptc_ad(i, j, k) = 0.0_8
              END DO
            END DO
          END IF
          dd_divg = d4_bg
          CALL POPREAL8ARRAY(pkd(isd, jsd, k), ied - isd + 1)
          CALL POPREAL8ARRAY(delp(isd, jsd, k), ied - isd + 1)
          CALL POPREAL8ARRAY(ptc(isd, jsd, k), ied - isd + 1)
          CALL POPREAL8ARRAY(pt(isd, jsd, k), ied - isd + 1)
          CALL POPREAL8ARRAY(u(isd, jsd, k), ied - isd + 1)
          CALL POPREAL8ARRAY(v(isd, jsd, k), ied - isd + 2)
          CALL POPREAL8ARRAY(w(isd, jsd, k), ied - isd + 1)
          CALL POPREAL8ARRAY(uc(isd, jsd, k), ied - isd + 2)
          CALL POPREAL8ARRAY(vc(isd, jsd, k), ied - isd + 1)
          CALL POPREAL8ARRAY(divg_d(isd, jsd, k), ied - isd + 2)
          CALL POPREAL8ARRAY(crx(is, jsd, k), ie - is + 2)
          CALL POPREAL8ARRAY(cry(isd, js, k), ied - isd + 1)
          CALL POPREAL8ARRAY(xfx(is, jsd, k), ie - is + 2)
          CALL POPREAL8ARRAY(yfx(isd, js, k), ied - isd + 1)
          CALL POPREAL8ARRAY(q, (ied-isd+1)*(jed-jsd+1)*npz*nq)
          CALL D_SW_ADM(pkd(isd:, jsd, k), pkd_ad(isd:, jsd, k), delp(&
&                 isd:, jsd, k), delp_ad(isd:, jsd, k), ptc(isd:, jsd, k&
&                 ), ptc_ad(isd:, jsd, k), pt(isd:, jsd, k), pt_ad(isd:&
&                 , jsd, k), u(isd:, jsd, k), u_ad(isd:, jsd, k), v(isd:&
&                 , jsd, k), v_ad(isd:, jsd, k), w(isd:, jsd, k), w_ad(&
&                 isd:, jsd, k), uc(isd:, jsd, k), uc_ad(isd:, jsd, k), &
&                 vc(isd:, jsd, k), vc_ad(isd:, jsd, k), ua(isd, jsd, k)&
&                 , ua_ad(isd, jsd, k), va(isd, jsd, k), va_ad(isd, jsd&
&                 , k), divg_d(isd:, jsd, k), divg_d_ad(isd:, jsd, k), &
&                 mfx(is, js, k), mfx_ad(is, js, k), mfy(is, js, k), &
&                 mfy_ad(is, js, k), cx(is, jsd, k), cx_ad(is, jsd, k), &
&                 cy(isd, js, k), cy_ad(isd, js, k), crx(is:, jsd, k), &
&                 crx_ad(is:, jsd, k), cry(isd:, js, k), cry_ad(isd:, js&
&                 , k), xfx(is:, jsd, k), xfx_ad(is:, jsd, k), yfx(isd:&
&                 , js, k), yfx_ad(isd:, js, k), zvir, sphum, nq, q, &
&                 q_ad, k, npz, inline_q, pkz(is, js, k), pkz_ad(is, js&
&                 , k), dt, hord_tr, hord_m, hord_v, hord_t, hord_p, &
&                 hord_tr_pert, hord_m_pert, hord_v_pert, hord_t_pert, hord_p_pert, &
&                 nord_k, damp_k, d2_divg, dd_divg, vtdm4, d_con, &
&                 hydrostatic, ppm_limiter)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8ARRAY(wk, (ied-isd+1)*(jed-jsd+1))
            CALL A2B_ORD2_ADM(delp(:, :, k), delp_ad(:, :, k), wk, wk_ad&
&                       , npx, npy, is, ie, js, je, ng, .false.)
          END IF
          CALL POPCONTROL3B(branch)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(d2_divg)
          ELSE
            CALL POPREAL8(d2_divg)
          END IF
          CALL POPREAL8(damp_k)
          CALL POPINTEGER4(nord_k)
          CALL POPINTEGER4(hord_p_pert)
          CALL POPINTEGER4(hord_v_pert)
          CALL POPINTEGER4(hord_t_pert)
          CALL POPINTEGER4(hord_m_pert)
          CALL POPINTEGER4(hord_p)
          CALL POPINTEGER4(hord_v)
          CALL POPINTEGER4(hord_t)
          CALL POPINTEGER4(hord_m)
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          call mpp_update_domains(divg_d_ad, domain, position=CORNER)
          CALL DIVERGENCE_CORNER_ADM(u, u_ad, v, v_ad, ua, ua_ad, va, &
&                              va_ad, divg_d, divg_d_ad, npz)
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8ARRAY(q, (ied-isd+1)*(jed-jsd+1)*npz*nq)
          call mpp_update_domains( q_ad,  domain, complete=.true. )
        END IF
        call mpp_update_domains(uc_ad, vc_ad, domain, gridtype=CGRID_NE, complete=.true.)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          v_ad = v_ad + vm_ad
          u_ad = u_ad + um_ad
          um_ad = 0.0_8
          vm_ad = 0.0_8
          DO k=npz,1,-1
            arg10 = nord .GT. 0
            CALL POPREAL8ARRAY(ua(isd, jsd, k), ied - isd + 1)
            CALL POPREAL8ARRAY(va(isd, jsd, k), ied - isd + 1)
            CALL D2A2C_ADM(u(isd, jsd, k), u_ad(isd, jsd, k), v(isd, jsd&
&                    , k), v_ad(isd, jsd, k), u2, u2_ad, v2, v2_ad, ua(&
&                    isd:, jsd, k), ua_ad(isd:, jsd, k), va(isd:, jsd, k&
&                    ), va_ad(isd:, jsd, k), uc(isd, jsd, k), uc_ad(isd&
&                    , jsd, k), vc(isd, jsd, k), vc_ad(isd, jsd, k), &
&                    arg10)
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              DO j=jed,jsd,-1
                DO i=ied+1,isd,-1
                  v_ad(i, j, k) = v_ad(i, j, k) + v2_ad(i, j)
                  dv_ad(i, j, k) = dv_ad(i, j, k) + 0.5*v2_ad(i, j)
                  v2_ad(i, j) = 0.0_8
                END DO
              END DO
              DO j=jed+1,jsd,-1
                DO i=ied,isd,-1
                  u_ad(i, j, k) = u_ad(i, j, k) + u2_ad(i, j)
                  du_ad(i, j, k) = du_ad(i, j, k) + 0.5*u2_ad(i, j)
                  u2_ad(i, j) = 0.0_8
                END DO
              END DO
            ELSE
              DO j=jed,jsd,-1
                DO i=ied+1,isd,-1
                  v_ad(i, j, k) = v_ad(i, j, k) + 1.5*v2_ad(i, j)
                  vm_ad(i, j, k) = vm_ad(i, j, k) - 0.5*v2_ad(i, j)
                  v2_ad(i, j) = 0.0_8
                END DO
              END DO
              DO j=jed+1,jsd,-1
                DO i=ied,isd,-1
                  u_ad(i, j, k) = u_ad(i, j, k) + 1.5*u2_ad(i, j)
                  um_ad(i, j, k) = um_ad(i, j, k) - 0.5*u2_ad(i, j)
                  u2_ad(i, j) = 0.0_8
                END DO
              END DO
            END IF
          END DO
        ELSE
          DO k=npz,1,-1
            CALL POPINTEGER4(ad_to2)
            DO j=ad_to2,js,-1
              DO i=ie,is,-1
                temp4 = wk(i, j-1) + wk(i, j)
                temp8 = pkc(i, j-1, k+1) - pkc(i, j, k)
                temp7 = gz(i, j-1, k) - gz(i, j, k+1)
                temp6 = pkc(i, j, k+1) - pkc(i, j-1, k)
                temp5 = gz(i, j-1, k+1) - gz(i, j, k)
                temp_ad1 = dt2*rdyc(i, j)*vc_ad(i, j, k)/temp4
                temp_ad2 = -((temp5*temp6+temp7*temp8)*temp_ad1/temp4)
                gz_ad(i, j-1, k+1) = gz_ad(i, j-1, k+1) + temp6*temp_ad1
                gz_ad(i, j, k) = gz_ad(i, j, k) - temp6*temp_ad1
                pkc_ad(i, j, k+1) = pkc_ad(i, j, k+1) + temp5*temp_ad1
                pkc_ad(i, j-1, k) = pkc_ad(i, j-1, k) - temp5*temp_ad1
                gz_ad(i, j-1, k) = gz_ad(i, j-1, k) + temp8*temp_ad1
                gz_ad(i, j, k+1) = gz_ad(i, j, k+1) - temp8*temp_ad1
                pkc_ad(i, j-1, k+1) = pkc_ad(i, j-1, k+1) + temp7*&
&                 temp_ad1
                pkc_ad(i, j, k) = pkc_ad(i, j, k) - temp7*temp_ad1
                wk_ad(i, j-1) = wk_ad(i, j-1) + temp_ad2
                wk_ad(i, j) = wk_ad(i, j) + temp_ad2
              END DO
            END DO
            DO j=je,js,-1
              CALL POPINTEGER4(ad_to1)
              DO i=ad_to1,is,-1
                temp = wk(i-1, j) + wk(i, j)
                temp3 = pkc(i-1, j, k+1) - pkc(i, j, k)
                temp2 = gz(i-1, j, k) - gz(i, j, k+1)
                temp1 = pkc(i, j, k+1) - pkc(i-1, j, k)
                temp0 = gz(i-1, j, k+1) - gz(i, j, k)
                temp_ad = dt2*rdxc(i, j)*uc_ad(i, j, k)/temp
                temp_ad0 = -((temp0*temp1+temp2*temp3)*temp_ad/temp)
                gz_ad(i-1, j, k+1) = gz_ad(i-1, j, k+1) + temp1*temp_ad
                gz_ad(i, j, k) = gz_ad(i, j, k) - temp1*temp_ad
                pkc_ad(i, j, k+1) = pkc_ad(i, j, k+1) + temp0*temp_ad
                pkc_ad(i-1, j, k) = pkc_ad(i-1, j, k) - temp0*temp_ad
                gz_ad(i-1, j, k) = gz_ad(i-1, j, k) + temp3*temp_ad
                gz_ad(i, j, k+1) = gz_ad(i, j, k+1) - temp3*temp_ad
                pkc_ad(i-1, j, k+1) = pkc_ad(i-1, j, k+1) + temp2*&
&                 temp_ad
                pkc_ad(i, j, k) = pkc_ad(i, j, k) - temp2*temp_ad
                wk_ad(i-1, j) = wk_ad(i-1, j) + temp_ad0
                wk_ad(i, j) = wk_ad(i, j) + temp_ad0
              END DO
            END DO
            CALL POPINTEGER4(ad_to0)
            DO j=ad_to0,jsm1,-1
              CALL POPINTEGER4(ad_to)
              DO i=ad_to,ism1,-1
                CALL POPREAL8(wk(i, j))
                pkc_ad(i, j, k+1) = pkc_ad(i, j, k+1) + wk_ad(i, j)
                pkc_ad(i, j, k) = pkc_ad(i, j, k) - wk_ad(i, j)
                wk_ad(i, j) = 0.0_8
              END DO
            END DO
          END DO
          CALL POPREAL8ARRAY(pe, (ie-is+3)*(npz+1)*(je-js+3))
          CALL POPREAL8ARRAY(peln, (ie-is+1)*(npz+1)*(je-js+1))
          CALL POPREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*(&
&                      npz+1)/8)
          CALL POPREAL8ARRAY(pkz, (ie-is+1)*(je-js+1)*npz)
          CALL GEOPK_ADM(ptop, pe, pe_ad, peln, peln_ad, delpc, delpc_ad&
&                  , pkc, pkc_ad, gz, gz_ad, phis, ptc, ptc_ad, pkz, &
&                  pkz_ad, npz, akap, .true.)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            alpha = 1. - beta1
            CALL POPREAL8ARRAY(ptc, (ied-isd+1)*(jed-jsd+1)*npz)
            pt_ad = pt_ad + beta1*ptc_ad
            ptc_ad(:, :, :) = alpha*ptc_ad(:, :, :)
            CALL POPREAL8ARRAY(delpc, (ied-isd+1)*(jed-jsd+1)*npz)
            delp_ad = delp_ad + beta1*delpc_ad
            delpc_ad(:, :, :) = alpha*delpc_ad(:, :, :)
            CALL POPREAL8(alpha)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8ARRAY(omga, (ied-isd+1)*(jed-jsd+1)*npz)
            CALL POPREAL8ARRAY(delpc, (ied-isd+1)*(jed-jsd+1)*npz)
            CALL POPREAL8ARRAY(ptc, (ied-isd+1)*(jed-jsd+1)*npz)
            omga_ad = 0.0_8
            CALL MIX_DP_ADM(hydrostatic, omga, omga_ad, delpc, delpc_ad&
&                     , ptc, ptc_ad, npz, ak, bk, .true.)
          END IF
          DO k=npz,1,-1
            arg10 = nord .GT. 0
            CALL POPREAL8ARRAY(delpc(isd, jsd, k), ied - isd + 1)
            CALL POPREAL8ARRAY(delp(isd, jsd, k), ied - isd + 1)
            CALL POPREAL8ARRAY(ptc(isd, jsd, k), ied - isd + 1)
            CALL POPREAL8ARRAY(pt(isd, jsd, k), ied - isd + 1)
            CALL POPREAL8ARRAY(uc(isd, jsd, k), ied - isd + 2)
            CALL POPREAL8ARRAY(vc(isd, jsd, k), ied - isd + 1)
            CALL POPREAL8ARRAY(ua(isd, jsd, k), ied - isd + 1)
            CALL POPREAL8ARRAY(va(isd, jsd, k), ied - isd + 1)
            CALL POPREAL8ARRAY(ut(isd, jsd, k), ied - isd + 1)
            CALL POPREAL8ARRAY(vt(isd, jsd, k), ied - isd + 1)
            CALL C_SW_ADM(delpc(isd:, jsd, k), delpc_ad(isd:, jsd, k), &
&                   delp(isd:, jsd, k), delp_ad(isd:, jsd, k), ptc(isd:&
&                   , jsd, k), ptc_ad(isd:, jsd, k), pt(isd:, jsd, k), &
&                   pt_ad(isd:, jsd, k), u(isd, jsd, k), u_ad(isd, jsd, &
&                   k), v(isd, jsd, k), v_ad(isd, jsd, k), w(isd, jsd, k&
&                   ), w_ad(isd, jsd, k), uc(isd:, jsd, k), uc_ad(isd:, &
&                   jsd, k), vc(isd:, jsd, k), vc_ad(isd:, jsd, k), ua(&
&                   isd:, jsd, k), ua_ad(isd:, jsd, k), va(isd:, jsd, k)&
&                   , va_ad(isd:, jsd, k), omga(isd, jsd, k), ut(isd:, &
&                   jsd, k), ut_ad(isd:, jsd, k), vt(isd:, jsd, k), &
&                   vt_ad(isd:, jsd, k), dt2, hydrostatic, arg10)
          END DO
        END IF
        CALL POPCONTROL1B(branch)
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .NE. 0) THEN
        IF (branch .NE. 1) THEN
          CALL POPREAL8(dt2)
          CALL POPREAL8(dt)
        END IF
        DO k=npz,1,-1
          DO j=je+1,js,-1
            DO i=ie,is,-1
              CALL POPCONTROL1B(branch)
              u_ad(i, j, k) = u_ad(i, j, k) + ndt*rdx(i, j)*cry_ad(i, j&
&               , k)
              cry_ad(i, j, k) = 0.0_8
            END DO
          END DO
          DO j=je,js,-1
            DO i=ie+1,is,-1
              CALL POPCONTROL1B(branch)
              v_ad(i, j, k) = v_ad(i, j, k) + ndt*rdy(i, j)*crx_ad(i, j&
&               , k)
              crx_ad(i, j, k) = 0.0_8
            END DO
          END DO
        END DO
      END IF
    END DO
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      call mpp_update_domains(du_ad, dv_ad, domain, gridtype=DGRID_NE)
    ELSE IF (branch .NE. 1) THEN
      GOTO 100
    END IF
    CALL POPREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*(npz+1)&
&                /8)
    CALL POPREAL8ARRAY(gz, pg_precision*(ied-isd+1)*(jed-jsd+1)*(npz+1)/&
&                8)
    CALL POPREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1)*npz)
    CALL GRAD1_P_ADM(du, du_ad, dv, dv_ad, pkc, pkc_ad, gz, gz_ad, delp&
&              , delp_ad, dt, ng, npx, npy, npz, ptop, ptk, hydrostatic)
    CALL POPREAL8ARRAY(pe, (ie-is+3)*(npz+1)*(je-js+3))
    CALL POPREAL8ARRAY(peln, (ie-is+1)*(npz+1)*(je-js+1))
    CALL POPREAL8ARRAY(pkc, pg_precision*(ied-isd+1)*(jed-jsd+1)*(npz+1)&
&                /8)
    CALL POPREAL8ARRAY(pkz, (ie-is+1)*(je-js+1)*npz)
    CALL GEOPK_ADM(ptop, pe, pe_ad, peln, peln_ad, delp, delp_ad, pkc, &
&            pkc_ad, gz, gz_ad, phis, pt, pt_ad, pkz, pkz_ad, npz, akap&
&            , .false.)
 100 CALL POPREAL8ARRAY(v, (ied-isd+2)*(jed-jsd+1)*npz)
    CALL POPREAL8ARRAY(u, (ied-isd+1)*(jed-jsd+2)*npz)
    call mpp_update_domains( u_ad, v_ad, domain, gridtype=DGRID_NE, complete=.true. )
    CALL POPREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1)*npz)
    call mpp_update_domains( delp_ad, domain, complete=.true. )    
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1)*npz)
      IF (npz .GT. 1) call mpp_update_domains( pt,   domain, complete=.false. )      
    END IF
  END SUBROUTINE DYN_CORE_ADM
!-----------------------------------------------------------------------
!     dyn_core :: FV Lagrangian dynamics driver
!-----------------------------------------------------------------------
!, init_step, end_step, time_total)
  SUBROUTINE DYN_CORE(npx, npy, npz, ng, sphum, nq, bdt, n_split, zvir, &
&   cp, akap, rdgas, grav, hydrostatic, u, v, um, vm, w, delz, pt, q, &
&   delp, pe, pk, phis, omga, ptop, pfull, ua, va, uc, vc, mfx, mfy, cx&
&   , cy, pem, pkz, peln, ak, bk)
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
! D grid meridional wind (m/s)
    REAL, DIMENSION(isd:ied + 1, jsd:jed, npz), INTENT(INOUT) :: vm
! D grid zonal wind (m/s)
    REAL, DIMENSION(isd:ied, jsd:jed + 1, npz), INTENT(INOUT) :: u
! D grid meridional wind (m/s)
    REAL, DIMENSION(isd:ied + 1, jsd:jed, npz), INTENT(INOUT) :: v
! vertical vel. (m/s)
    REAL, INTENT(INOUT) :: w(isd:ied, jsd:jed, npz)
! delta-height (m)
    REAL, INTENT(INOUT) :: delz(is:ie, js:je, npz)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
!
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, npz, nq)
!    real, intent(IN), optional:: time_total  ! total time (seconds) since start
!THESE ARE ACTUALLY IN THE MODULE
    REAL :: delzc(is:ie, js:je, npz)
    REAL :: ut(isd:ied, jsd:jed, npz)
    REAL :: vt(isd:ied, jsd:jed, npz)
    REAL :: crx(is:ie+1, jsd:jed, npz), xfx(is:ie+1, jsd:jed, npz)
    REAL :: cry(isd:ied, js:je+1, npz), yfx(isd:ied, js:je+1, npz)
    REAL :: divg_d(isd:ied+1, jsd:jed+1, npz)
    REAL :: zh(isd:ied, jsd:jed, npz)
    REAL(pg_precision) :: du(isd:ied, jsd:jed+1, npz), dv(isd:ied+1, jsd&
&   :jed, npz)
    REAL(pg_precision) :: pkc(isd:ied, jsd:jed, npz+1)
    REAL :: pkd(isd:ied, jsd:jed, npz+1)
    REAL :: delpc(isd:ied, jsd:jed, npz)
    REAL(pg_precision) :: pk3(isd:ied, jsd:jed, npz+1)
    REAL :: ptc(isd:ied, jsd:jed, npz)
    REAL(pg_precision) :: gz(isd:ied, jsd:jed, npz+1)
!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
! Surface geopotential (g*Z_surf)
    REAL, INTENT(INOUT) :: phis(isd:ied, jsd:jed)
! edge pressure (pascal)
    REAL(p_precision), INTENT(INOUT) :: pe(is-1:ie+1, npz+1, js-1:je+1)
    REAL, INTENT(OUT) :: pem(is-1:ie+1, npz+1, js-1:je+1)
! ln(pe)
    REAL(p_precision), INTENT(OUT) :: peln(is:ie, npz+1, js:je)
! pe**kappa
    REAL(p_precision), INTENT(INOUT) :: pk(is:ie, js:je, npz+1)
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
! Vertical pressure velocity (pa/s)
    REAL, INTENT(INOUT) :: omga(isd:ied, jsd:jed, npz)
! (uc, vc) are mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: vc(isd:ied, jsd:jed+1, npz)
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua, va
! The Flux capacitors: accumulated Mass flux arrays
    REAL, INTENT(INOUT) :: mfx(is:ie+1, js:je, npz)
    REAL, INTENT(INOUT) :: mfy(is:ie, js:je+1, npz)
! Accumulated Courant number arrays
    REAL, INTENT(INOUT) :: cx(is:ie+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: cy(isd:ied, js:je+1, npz)
! Work:
! 
    REAL(p_precision), INTENT(INOUT) :: pkz(is:ie, js:je, npz)
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
    REAL :: wk(isd:ied, jsd:jed)
! --- For no_cgrid option ---
    REAL :: u2(isd:ied, jsd:jed+1)
    REAL :: v2(isd:ied+1, jsd:jed)
!-------------------------------------
    REAL :: d4_bg, d2_bg
    INTEGER :: hord_m, hord_v, hord_t, hord_p, nord, nord_k
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
    INTEGER :: arg1
    LOGICAL :: arg10
    INTEGER :: arg2
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
    IF (npz .GT. 1) call mpp_update_domains( pt,   domain, complete=.false. )
!oncall mpp_update_domains( pt,   domain, complete=.false. )
!oncall mpp_update_domains( delp, domain, complete=.true. )
    call mpp_update_domains( delp, domain, complete=.true. )
!oncall mpp_update_domains( u, v, domain, gridtype=DGRID_NE, complete=.true. )
    call mpp_update_domains( u, v, domain, gridtype=DGRID_NE, complete=.true. )
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
      CALL GEOPK(ptop, pe, peln, delp, pkc, gz, phis, pt, pkz, npz, akap&
&          , .false.)
      CALL GRAD1_P(du, dv, pkc, gz, delp, dt, ng, npx, npy, npz, ptop, &
&            ptk, hydrostatic)
      IF (init_wind_m) THEN
!oncall mpp_update_domains(du, dv, domain, gridtype=DGRID_NE)
         call mpp_update_domains(du, dv, domain, gridtype=DGRID_NE)
      END IF
    END IF
! Empty the "flux capacitors"
    mfx(:, :, :) = 0.
    mfy(:, :, :) = 0.
    cx(:, :, :) = 0.
    cy(:, :, :) = 0.
    elapsed_dt = 0.0
    subcycle = 1
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
!oncall mp_reduce_max(cmax)
!oncall mp_barrier()
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
!onif (test_case==9) call case9_forcing1(phis, time_total)
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
                  u2(i, j) = u(i, j, k) + 0.5*du(i, j, k)
                END DO
              END DO
              DO j=jsd,jed
                DO i=isd,ied+1
                  v2(i, j) = v(i, j, k) + 0.5*dv(i, j, k)
                END DO
              END DO
            ELSE
              DO j=jsd,jed+1
                DO i=isd,ied
                  u2(i, j) = 1.5*u(i, j, k) - 0.5*um(i, j, k)
                END DO
              END DO
              DO j=jsd,jed
                DO i=isd,ied+1
                  v2(i, j) = 1.5*v(i, j, k) - 0.5*vm(i, j, k)
                END DO
              END DO
            END IF
            arg10 = nord .GT. 0
            CALL D2A2C(u(isd, jsd, k), v(isd, jsd, k), u2, v2, ua(isd, &
&                jsd, k), va(isd, jsd, k), uc(isd, jsd, k), vc(isd, jsd&
&                , k), arg10)
          END DO
!if ( .not. hydrostatic ) delpc(:,:,:) = delp(:,:,:)
          um(:, :, :) = u(:, :, :)
          vm(:, :, :) = v(:, :, :)
!call timing_off('no_cgrid')
        ELSE
!call timing_on('c_sw')
!$omp parallel do default(shared) private(i,j,k)
          DO k=1,npz
            arg10 = nord .GT. 0
            CALL C_SW(delpc(isd, jsd, k), delp(isd, jsd, k), ptc(isd, &
&               jsd, k), pt(isd, jsd, k), u(isd, jsd, k), v(isd, jsd, k)&
&               , w(isd, jsd, k), uc(isd, jsd, k), vc(isd, jsd, k), ua(&
&               isd, jsd, k), va(isd, jsd, k), omga(isd, jsd, k), ut(isd&
&               , jsd, k), vt(isd, jsd, k), dt2, hydrostatic, arg10)
          END DO
! on output omga is updated w
!call timing_off('c_sw')
          IF (fill_dp) CALL MIX_DP(hydrostatic, omga, delpc, ptc, npz, &
&                            ak, bk, .true.)
!      if ( hydrostatic ) then
          IF (beta1 .GT. 0.001) THEN
            alpha = 1. - beta1
            delpc(:, :, :) = beta1*delp(:, :, :) + alpha*delpc(:, :, :)
            ptc(:, :, :) = beta1*pt(:, :, :) + alpha*ptc(:, :, :)
          END IF
          CALL GEOPK(ptop, pe, peln, delpc, pkc, gz, phis, ptc, pkz, npz&
&              , akap, .true.)
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
                uc(i, j, k) = uc(i, j, k) + dt2*rdxc(i, j)/(wk(i-1, j)+&
&                 wk(i, j))*((gz(i-1, j, k+1)-gz(i, j, k))*(pkc(i, j, k+&
&                 1)-pkc(i-1, j, k))+(gz(i-1, j, k)-gz(i, j, k+1))*(pkc(&
&                 i-1, j, k+1)-pkc(i, j, k)))
              END DO
            END DO
            DO j=js,jeb1
              DO i=is,ie
                vc(i, j, k) = vc(i, j, k) + dt2*rdyc(i, j)/(wk(i, j-1)+&
&                 wk(i, j))*((gz(i, j-1, k+1)-gz(i, j, k))*(pkc(i, j, k+&
&                 1)-pkc(i, j-1, k))+(gz(i, j-1, k)-gz(i, j, k+1))*(pkc(&
&                 i, j-1, k+1)-pkc(i, j, k)))
              END DO
            END DO
          END DO
        END IF
!call timing_on('COMM_TOTAL')
!oncall mpp_update_domains(uc, vc, domain, gridtype=CGRID_NE, complete=.true.)
        call mpp_update_domains(uc, vc, domain, gridtype=CGRID_NE, complete=.true.)
!onif (test_case==9) call case9_forcing2(phis)
        call mpp_update_domains( q,  domain, complete=.true. )
!call timing_on('COMM_TOTAL')
!oncall mpp_update_domains( q,  domain, complete=.true. )
!call timing_off('COMM_TOTAL')
        IF (nord .GT. 0) THEN
          CALL DIVERGENCE_CORNER(u, v, ua, va, divg_d, npz)
!call timing_on('COMM_TOTAL')
!oncall mpp_update_domains(divg_d, domain, position=CORNER)
          call mpp_update_domains(divg_d, domain, position=CORNER)
!call timing_off('COMM_TOTAL')
        END IF
!call timing_on('d_sw')
!$omp parallel do default(shared) private(i, j, k, nord_k, damp_k, d2_divg, dd_divg, hord_m, hord_v, hord_t, hord_p, wk)
        DO k=1,npz
          hord_m = hord_mt
          hord_t = hord_tm
          hord_v = hord_vt
          hord_p = hord_dp
          nord_k = nord
          damp_k = dddmp
          y3 = d2_bg*(1.-3.*TANH(0.1*LOG(pfull(k)/pfull(npz))))
          IF (0.20 .GT. y3) THEN
            d2_divg = y3
          ELSE
            d2_divg = 0.20
          END IF
          y4 = d2_bg*(1.-3.*TANH(0.1*LOG(pfull(k)/pfull(npz))))
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
          IF (d_ext .GT. 0.) CALL A2B_ORD2(delp(:, :, k), wk, npx, npy, &
&                                    is, ie, js, je, ng, .false.)
          CALL D_SW(pkd(isd, jsd, k), delp(isd, jsd, k), ptc(isd, jsd, k&
&             ), pt(isd, jsd, k), u(isd, jsd, k), v(isd, jsd, k), w(isd&
&             , jsd, k), uc(isd, jsd, k), vc(isd, jsd, k), ua(isd, jsd, &
&             k), va(isd, jsd, k), divg_d(isd, jsd, k), mfx(is, js, k), &
&             mfy(is, js, k), cx(is, jsd, k), cy(isd, js, k), crx(is, &
&             jsd, k), cry(isd, js, k), xfx(is, jsd, k), yfx(isd, js, k)&
&             , zvir, sphum, nq, q, k, npz, inline_q, pkz(is, js, k), dt&
&             , hord_tr, hord_m, hord_v, hord_t, hord_p, nord_k, damp_k&
&             , d2_divg, dd_divg, vtdm4, d_con, hydrostatic, ppm_limiter&
&            )
          IF (d_ext .GT. 0.) THEN
            DO j=js,jep1
              DO i=is,iep1
! delp at cell corners
                ptc(i, j, k) = wk(i, j)
              END DO
            END DO
          END IF
        END DO
!call timing_off('d_sw')
        IF (fill_dp) CALL MIX_DP(hydrostatic, w, delp, pt, npz, ak, bk, &
&                          .false.)
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
              wk(i, j) = ptc(i, j, 1)
              divg2(i, j) = wk(i, j)*pkd(i, j, 1)
            END DO
            DO k=2,npz
              DO i=is,iep1
                wk(i, j) = wk(i, j) + ptc(i, j, k)
                divg2(i, j) = divg2(i, j) + ptc(i, j, k)*pkd(i, j, k)
              END DO
            END DO
            DO i=is,iep1
              divg2(i, j) = d2_divg*divg2(i, j)/wk(i, j)
            END DO
          END DO
        ELSE
          divg2 = 0.
        END IF
!call timing_on('COMM_TOTAL')
!oncall mpp_update_domains(  pt, domain, complete=.false.)
        call mpp_update_domains(  pt, domain, complete=.false.)
!oncall mpp_update_domains(delp, domain, complete=.true.)
        call mpp_update_domains(delp, domain, complete=.true.)
!call timing_off('COMM_TOTAL')
!     if ( hydrostatic ) then
        CALL GEOPK(ptop, pe, peln, delp, pkc, gz, phis, pt, pkz, npz, &
&            akap, .false.)
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
                u(i, j, k) = (u(i, j, k)+divg2(i, j)-divg2(i+1, j))*rdx(&
&                 i, j) + beta*du(i, j, k)
              END DO
            END DO
          END DO
          DO k=1,npz
            DO j=js,je
              DO i=is,ie+1
                v(i, j, k) = (v(i, j, k)+divg2(i, j)-divg2(i, j+1))*rdy(&
&                 i, j) + beta*dv(i, j, k)
              END DO
            END DO
          END DO
          CALL GRAD1_P(du, dv, pkc, gz, delp, dt, ng, npx, npy, npz, &
&                ptop, ptk, hydrostatic)
          alpha = 1. - beta
          DO k=1,npz
            DO j=js,je+1
              DO i=is,ie
                u(i, j, k) = u(i, j, k) + alpha*du(i, j, k)
              END DO
            END DO
          END DO
          DO k=1,npz
            DO j=js,je
              DO i=is,ie+1
                v(i, j, k) = v(i, j, k) + alpha*dv(i, j, k)
              END DO
            END DO
          END DO
        ELSE
          CALL ONE_GRAD_P(u, v, pkc, gz, divg2, delp, dt, ng, npx, npy, &
&                   npz, ptop, ptk, hydrostatic)
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
!oncall mpp_get_boundary(u, v, domain, wbuffery=wbuffer, ebuffery=ebuffer,  &
!                     sbufferx=sbuffer, nbufferx=nbuffer, gridtype=DGRID_NE )
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
        ELSE
!oncall mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
          call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
        END IF
!call timing_off('COMM_TOTAL')
!init_wind_m = .false. !dh always false
!first_step = .false.
        elapsed_dt = elapsed_dt + dt
      END DO
    END DO
  END SUBROUTINE DYN_CORE
  SUBROUTINE TWO_GRAD_P(u, v, pk, gh, delp, pkt, divg2, dt, ng, npx, npy&
&   , npz, ptk)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    REAL(p_precision), INTENT(IN) :: ptk
    REAL(pg_precision), INTENT(IN) :: divg2(is:ie+1, js:je+1)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
! perturbation pressure
    REAL(pg_precision), INTENT(INOUT) :: pk(isd:ied, jsd:jed, npz+1)
! p**kappa
    REAL(pg_precision), INTENT(INOUT) :: pkt(isd:ied, jsd:jed, npz+1)
! g * zh
    REAL(pg_precision), INTENT(INOUT) :: gh(isd:ied, jsd:jed, npz+1)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, npz)
! Local:
    REAL(pg_precision) :: dp(isd:ied, jsd:jed)
    REAL :: wk1(isd:ied, jsd:jed)
    REAL :: wk(is:ie+1, js:je+1)
    INTEGER :: iep1, jep1
    INTEGER :: i, j, k
    iep1 = ie + 1
    jep1 = je + 1
    DO j=js,jep1
      DO i=is,iep1
        pk(i, j, 1) = 0.
        pkt(i, j, 1) = ptk
      END DO
    END DO
    DO k=1,npz+1
      IF (k .NE. 1) THEN
        IF (a2b_ord .EQ. 4) THEN
          CALL A2B_ORD4(pk(isd:, jsd:, k), wk1, npx, npy, is, ie, js, je&
&                 , ng, .true.)
          CALL A2B_ORD4(pkt(isd:, jsd:, k), wk1, npx, npy, is, ie, js, &
&                 je, ng, .true.)
        ELSE
          CALL A2B_ORD2(pk(isd:, jsd:, k), wk1, npx, npy, is, ie, js, je&
&                 , ng, .true.)
          CALL A2B_ORD2(pkt(isd:, jsd:, k), wk1, npx, npy, is, ie, js, &
&                 je, ng, .true.)
        END IF
      END IF
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(gh(isd:, jsd:, k), wk1, npx, npy, is, ie, js, je, &
&               ng, .true.)
      ELSE
        CALL A2B_ORD2(gh(isd:, jsd:, k), wk1, npx, npy, is, ie, js, je, &
&               ng, .true.)
      END IF
    END DO
    DO k=1,npz
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(delp(:, :, k), wk1, npx, npy, is, ie, js, je, ng)
      ELSE
        CALL A2B_ORD2(delp(:, :, k), wk1, npx, npy, is, ie, js, je, ng)
      END IF
      DO j=js,jep1
        DO i=is,iep1
          wk(i, j) = pkt(i, j, k+1) - pkt(i, j, k)
        END DO
      END DO
      DO j=js,jep1
        DO i=is,ie
!------------------
! Perturbation term:
!------------------
          u(i, j, k) = u(i, j, k) + dt/(wk1(i, j)+wk1(i+1, j))*((gh(i, j&
&           , k+1)-gh(i+1, j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gh(i, j&
&           , k)-gh(i+1, j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))
!-----------------
! Hydrostatic term
!-----------------
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
          v(i, j, k) = v(i, j, k) + dt/(wk1(i, j)+wk1(i, j+1))*((gh(i, j&
&           , k+1)-gh(i, j+1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gh(i, j&
&           , k)-gh(i, j+1, k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))
!-----------------
! Hydrostatic term
!-----------------
          v(i, j, k) = rdy(i, j)*(divg2(i, j)-divg2(i, j+1)+v(i, j, k)+&
&           dt/(wk(i, j)+wk(i, j+1))*((gh(i, j, k+1)-gh(i, j+1, k))*(pkt&
&           (i, j+1, k+1)-pkt(i, j, k))+(gh(i, j, k)-gh(i, j+1, k+1))*(&
&           pkt(i, j, k+1)-pkt(i, j+1, k))))
        END DO
      END DO
    END DO
  END SUBROUTINE TWO_GRAD_P
!  Differentiation of one_grad_p in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: u v delp pk divg2 gh
!   with respect to varying inputs: u v delp pk divg2 gh
  SUBROUTINE ONE_GRAD_P_ADM(u, u_ad, v, v_ad, pk, pk_ad, gh, gh_ad, &
&   divg2, divg2_ad, delp, delp_ad, dt, ng, npx, npy, npz, ptop, ptk, &
&   hydrostatic)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    REAL(p_precision), INTENT(IN) :: ptop, ptk
    LOGICAL, INTENT(IN) :: hydrostatic
    REAL(pg_precision), INTENT(IN) :: divg2(is:ie+1, js:je+1)
    REAL(pg_precision) :: divg2_ad(is:ie+1, js:je+1)
    REAL(pg_precision), INTENT(INOUT) :: pk(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision) :: pk_ad(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: gh(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision) :: gh_ad(isd:ied, jsd:jed, npz+1)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp_ad(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: u_ad(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: v_ad(isd:ied+1, jsd:jed, npz)
! Local:
    REAL, DIMENSION(isd:ied, jsd:jed) :: wk
    REAL, DIMENSION(isd:ied, jsd:jed) :: wk_ad
    REAL :: wk1(is:ie+1, js:je)
    REAL :: wk1_ad(is:ie+1, js:je)
    REAL :: wk2(is:ie, js:je+1)
    REAL :: wk2_ad(is:ie, js:je+1)
    REAL :: top_value
    INTEGER :: iep1, jep1
    INTEGER :: i, j, k
    REAL(pg_precision) :: dt_pg
    INTEGER :: branch
    REAL(pg_precision) :: temp3
    REAL(pg_precision) :: temp2
    REAL(pg_precision) :: temp1
    REAL(pg_precision) :: temp0
    REAL :: temp_ad4
    REAL(pg_precision) :: temp_ad3
    REAL(pg_precision) :: temp_ad2
    REAL :: temp_ad1
    REAL(pg_precision) :: temp_ad0
    REAL(pg_precision) :: temp_ad
    REAL :: temp
    REAL(pg_precision) :: temp8
    REAL(pg_precision) :: temp7
    REAL(pg_precision) :: temp6
    REAL(pg_precision) :: temp5
    REAL :: temp4
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
        pk(i, j, 1) = top_value
      END DO
    END DO
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(pk(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL A2B_ORD2(pk(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(gh(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL A2B_ORD2(gh(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js,jep1
          DO i=is,iep1
            CALL PUSHREAL8(wk(i, j))
            wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
          END DO
        END DO
        CALL PUSHCONTROL2B(2)
      ELSE IF (a2b_ord .EQ. 4) THEN
        CALL PUSHREAL8ARRAY(wk, (ied-isd+1)*(jed-jsd+1))
        CALL PUSHREAL8ARRAY(delp(:, :, k), (ied-isd+1)*(jed-jsd+1))
        CALL A2B_ORD4(delp(:, :, k), wk, npx, npy, is, ie, js, je, ng)
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHREAL8ARRAY(wk, (ied-isd+1)*(jed-jsd+1))
        CALL PUSHREAL8ARRAY(delp(:, :, k), (ied-isd+1)*(jed-jsd+1))
        CALL A2B_ORD2(delp(:, :, k), wk, npx, npy, is, ie, js, je, ng)
        CALL PUSHCONTROL2B(0)
      END IF
    END DO
    wk1_ad = 0.0_8
    wk2_ad = 0.0_8
    wk_ad = 0.0_8
    DO k=npz,1,-1
      DO j=je,js,-1
        DO i=iep1,is,-1
          temp4 = wk(i, j) + wk(i, j+1)
          temp8 = pk(i, j, k+1) - pk(i, j+1, k)
          temp7 = gh(i, j, k) - gh(i, j+1, k+1)
          temp6 = pk(i, j+1, k+1) - pk(i, j, k)
          temp5 = gh(i, j, k+1) - gh(i, j+1, k)
          temp_ad2 = rdy(i, j)*v_ad(i, j, k)
          temp_ad3 = dt_pg*temp_ad2/temp4
          temp_ad4 = -((temp5*temp6+temp7*temp8)*temp_ad3/temp4)
          wk1_ad(i, j) = wk1_ad(i, j) + temp_ad2
          gh_ad(i, j, k+1) = gh_ad(i, j, k+1) + temp6*temp_ad3
          gh_ad(i, j+1, k) = gh_ad(i, j+1, k) - temp6*temp_ad3
          pk_ad(i, j+1, k+1) = pk_ad(i, j+1, k+1) + temp5*temp_ad3
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp5*temp_ad3
          gh_ad(i, j, k) = gh_ad(i, j, k) + temp8*temp_ad3
          gh_ad(i, j+1, k+1) = gh_ad(i, j+1, k+1) - temp8*temp_ad3
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp7*temp_ad3
          pk_ad(i, j+1, k) = pk_ad(i, j+1, k) - temp7*temp_ad3
          wk_ad(i, j) = wk_ad(i, j) + temp_ad4
          wk_ad(i, j+1) = wk_ad(i, j+1) + temp_ad4
          v_ad(i, j, k) = temp_ad2
        END DO
      END DO
      DO j=jep1,js,-1
        DO i=ie,is,-1
          temp = wk(i, j) + wk(i+1, j)
          temp3 = pk(i, j, k+1) - pk(i+1, j, k)
          temp2 = gh(i, j, k) - gh(i+1, j, k+1)
          temp1 = pk(i+1, j, k+1) - pk(i, j, k)
          temp0 = gh(i, j, k+1) - gh(i+1, j, k)
          temp_ad = rdx(i, j)*u_ad(i, j, k)
          temp_ad0 = dt_pg*temp_ad/temp
          temp_ad1 = -((temp0*temp1+temp2*temp3)*temp_ad0/temp)
          wk2_ad(i, j) = wk2_ad(i, j) + temp_ad
          gh_ad(i, j, k+1) = gh_ad(i, j, k+1) + temp1*temp_ad0
          gh_ad(i+1, j, k) = gh_ad(i+1, j, k) - temp1*temp_ad0
          pk_ad(i+1, j, k+1) = pk_ad(i+1, j, k+1) + temp0*temp_ad0
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp0*temp_ad0
          gh_ad(i, j, k) = gh_ad(i, j, k) + temp3*temp_ad0
          gh_ad(i+1, j, k+1) = gh_ad(i+1, j, k+1) - temp3*temp_ad0
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp2*temp_ad0
          pk_ad(i+1, j, k) = pk_ad(i+1, j, k) - temp2*temp_ad0
          wk_ad(i, j) = wk_ad(i, j) + temp_ad1
          wk_ad(i+1, j) = wk_ad(i+1, j) + temp_ad1
          u_ad(i, j, k) = temp_ad
        END DO
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8ARRAY(delp(:, :, k), (ied-isd+1)*(jed-jsd+1))
        CALL POPREAL8ARRAY(wk, (ied-isd+1)*(jed-jsd+1))
        CALL A2B_ORD2_ADM(delp(:, :, k), delp_ad(:, :, k), wk, wk_ad, &
&                   npx, npy, is, ie, js, je, ng)
      ELSE IF (branch .EQ. 1) THEN
        CALL POPREAL8ARRAY(delp(:, :, k), (ied-isd+1)*(jed-jsd+1))
        CALL POPREAL8ARRAY(wk, (ied-isd+1)*(jed-jsd+1))
        CALL A2B_ORD4_ADM(delp(:, :, k), delp_ad(:, :, k), wk, wk_ad, &
&                   npx, npy, is, ie, js, je, ng)
      ELSE
        DO j=jep1,js,-1
          DO i=iep1,is,-1
            CALL POPREAL8(wk(i, j))
            pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + wk_ad(i, j)
            pk_ad(i, j, k) = pk_ad(i, j, k) - wk_ad(i, j)
            wk_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
    END DO
    DO j=je,js,-1
      DO i=iep1,is,-1
        divg2_ad(i, j) = divg2_ad(i, j) + wk1_ad(i, j)
        divg2_ad(i, j+1) = divg2_ad(i, j+1) - wk1_ad(i, j)
        wk1_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=jep1,js,-1
      DO i=ie,is,-1
        divg2_ad(i, j) = divg2_ad(i, j) + wk2_ad(i, j)
        divg2_ad(i+1, j) = divg2_ad(i+1, j) - wk2_ad(i, j)
        wk2_ad(i, j) = 0.0_8
      END DO
    END DO
    DO k=npz+1,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL A2B_ORD2_ADM(gh(isd:, jsd:, k), gh_ad(isd:, jsd:, k), wk, &
&                   wk_ad, npx, npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD4_ADM(gh(isd:, jsd:, k), gh_ad(isd:, jsd:, k), wk, &
&                   wk_ad, npx, npy, is, ie, js, je, ng, .true.)
      END IF
    END DO
    DO k=npz+1,2,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL A2B_ORD2_ADM(pk(isd:, jsd:, k), pk_ad(isd:, jsd:, k), wk, &
&                   wk_ad, npx, npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD4_ADM(pk(isd:, jsd:, k), pk_ad(isd:, jsd:, k), wk, &
&                   wk_ad, npx, npy, is, ie, js, je, ng, .true.)
      END IF
    END DO
    DO j=jep1,js,-1
      DO i=iep1,is,-1
        pk_ad(i, j, 1) = 0.0_8
      END DO
    END DO
  END SUBROUTINE ONE_GRAD_P_ADM
  SUBROUTINE ONE_GRAD_P(u, v, pk, gh, divg2, delp, dt, ng, npx, npy, npz&
&   , ptop, ptk, hydrostatic)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    REAL(p_precision), INTENT(IN) :: ptop, ptk
    LOGICAL, INTENT(IN) :: hydrostatic
    REAL(pg_precision), INTENT(IN) :: divg2(is:ie+1, js:je+1)
    REAL(pg_precision), INTENT(INOUT) :: pk(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: gh(isd:ied, jsd:jed, npz+1)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, npz)
! Local:
    REAL, DIMENSION(isd:ied, jsd:jed) :: wk
    REAL :: wk1(is:ie+1, js:je)
    REAL :: wk2(is:ie, js:je+1)
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
        pk(i, j, 1) = top_value
      END DO
    END DO
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(pk(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
      ELSE
        CALL A2B_ORD2(pk(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
      END IF
    END DO
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(gh(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
      ELSE
        CALL A2B_ORD2(gh(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
      END IF
    END DO
    DO j=js,jep1
      DO i=is,ie
        wk2(i, j) = divg2(i, j) - divg2(i+1, j)
      END DO
    END DO
    DO j=js,je
      DO i=is,iep1
        wk1(i, j) = divg2(i, j) - divg2(i, j+1)
      END DO
    END DO
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js,jep1
          DO i=is,iep1
            wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
          END DO
        END DO
      ELSE IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(delp(:, :, k), wk, npx, npy, is, ie, js, je, ng)
      ELSE
        CALL A2B_ORD2(delp(:, :, k), wk, npx, npy, is, ie, js, je, ng)
      END IF
      DO j=js,jep1
        DO i=is,ie
          u(i, j, k) = rdx(i, j)*(wk2(i, j)+u(i, j, k)+dt_pg/(wk(i, j)+&
&           wk(i+1, j))*((gh(i, j, k+1)-gh(i+1, j, k))*(pk(i+1, j, k+1)-&
&           pk(i, j, k))+(gh(i, j, k)-gh(i+1, j, k+1))*(pk(i, j, k+1)-pk&
&           (i+1, j, k))))
        END DO
      END DO
      DO j=js,je
        DO i=is,iep1
          v(i, j, k) = rdy(i, j)*(wk1(i, j)+v(i, j, k)+dt_pg/(wk(i, j)+&
&           wk(i, j+1))*((gh(i, j, k+1)-gh(i, j+1, k))*(pk(i, j+1, k+1)-&
&           pk(i, j, k))+(gh(i, j, k)-gh(i, j+1, k+1))*(pk(i, j, k+1)-pk&
&           (i, j+1, k))))
        END DO
      END DO
    END DO
  END SUBROUTINE ONE_GRAD_P
!  Differentiation of grad1_p in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: delp delu delv pk gh
!   with respect to varying inputs: delp delu delv pk gh
  SUBROUTINE GRAD1_P_ADM(delu, delu_ad, delv, delv_ad, pk, pk_ad, gh, &
&   gh_ad, delp, delp_ad, dt, ng, npx, npy, npz, ptop, ptk, hydrostatic)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    REAL(p_precision), INTENT(IN) :: ptk
    REAL(p_precision), INTENT(IN) :: ptop
    LOGICAL, INTENT(IN) :: hydrostatic
    REAL(pg_precision), INTENT(INOUT) :: pk(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision) :: pk_ad(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: gh(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision) :: gh_ad(isd:ied, jsd:jed, npz+1)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp_ad(isd:ied, jsd:jed, npz)
    REAL(pg_precision) :: delu(isd:ied, jsd:jed+1, npz)
    REAL(pg_precision) :: delu_ad(isd:ied, jsd:jed+1, npz)
    REAL(pg_precision) :: delv(isd:ied+1, jsd:jed, npz)
    REAL(pg_precision) :: delv_ad(isd:ied+1, jsd:jed, npz)
! Local:
    REAL :: wk(isd:ied, jsd:jed)
    REAL :: wk_ad(isd:ied, jsd:jed)
    REAL :: top_value
    INTEGER :: i, j, k
    REAL(pg_precision) :: dt_pg
    INTEGER :: branch
    REAL(pg_precision) :: temp3
    REAL(pg_precision) :: temp2
    REAL(pg_precision) :: temp1
    REAL(pg_precision) :: temp0
    REAL :: temp_ad2
    REAL(pg_precision) :: temp_ad1
    REAL :: temp_ad0
    REAL(pg_precision) :: temp_ad
    REAL :: temp
    REAL(pg_precision) :: temp8
    REAL(pg_precision) :: temp7
    REAL(pg_precision) :: temp6
    REAL(pg_precision) :: temp5
    REAL :: temp4
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
        pk(i, j, 1) = top_value
      END DO
    END DO
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(pk(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL A2B_ORD2(pk(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(gh(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL A2B_ORD2(gh(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js,je+1
          DO i=is,ie+1
            CALL PUSHREAL8(wk(i, j))
            wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
          END DO
        END DO
        CALL PUSHCONTROL2B(2)
      ELSE IF (a2b_ord .EQ. 4) THEN
        CALL PUSHREAL8ARRAY(wk, (ied-isd+1)*(jed-jsd+1))
        CALL PUSHREAL8ARRAY(delp(:, :, k), (ied-isd+1)*(jed-jsd+1))
        CALL A2B_ORD4(delp(:, :, k), wk, npx, npy, is, ie, js, je, ng)
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHREAL8ARRAY(wk, (ied-isd+1)*(jed-jsd+1))
        CALL PUSHREAL8ARRAY(delp(:, :, k), (ied-isd+1)*(jed-jsd+1))
        CALL A2B_ORD2(delp(:, :, k), wk, npx, npy, is, ie, js, je, ng)
        CALL PUSHCONTROL2B(0)
      END IF
    END DO
    wk_ad = 0.0_8
    DO k=npz,1,-1
      DO j=je,js,-1
        DO i=ie+1,is,-1
          temp4 = wk(i, j) + wk(i, j+1)
          temp8 = pk(i, j, k+1) - pk(i, j+1, k)
          temp7 = gh(i, j, k) - gh(i, j+1, k+1)
          temp6 = pk(i, j+1, k+1) - pk(i, j, k)
          temp5 = gh(i, j, k+1) - gh(i, j+1, k)
          temp_ad1 = rdy(i, j)*dt_pg*delv_ad(i, j, k)/temp4
          temp_ad2 = -((temp5*temp6+temp7*temp8)*temp_ad1/temp4)
          gh_ad(i, j, k+1) = gh_ad(i, j, k+1) + temp6*temp_ad1
          gh_ad(i, j+1, k) = gh_ad(i, j+1, k) - temp6*temp_ad1
          pk_ad(i, j+1, k+1) = pk_ad(i, j+1, k+1) + temp5*temp_ad1
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp5*temp_ad1
          gh_ad(i, j, k) = gh_ad(i, j, k) + temp8*temp_ad1
          gh_ad(i, j+1, k+1) = gh_ad(i, j+1, k+1) - temp8*temp_ad1
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp7*temp_ad1
          pk_ad(i, j+1, k) = pk_ad(i, j+1, k) - temp7*temp_ad1
          wk_ad(i, j) = wk_ad(i, j) + temp_ad2
          wk_ad(i, j+1) = wk_ad(i, j+1) + temp_ad2
          delv_ad(i, j, k) = 0.0_8
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie,is,-1
          temp = wk(i, j) + wk(i+1, j)
          temp3 = pk(i, j, k+1) - pk(i+1, j, k)
          temp2 = gh(i, j, k) - gh(i+1, j, k+1)
          temp1 = pk(i+1, j, k+1) - pk(i, j, k)
          temp0 = gh(i, j, k+1) - gh(i+1, j, k)
          temp_ad = rdx(i, j)*dt_pg*delu_ad(i, j, k)/temp
          temp_ad0 = -((temp0*temp1+temp2*temp3)*temp_ad/temp)
          gh_ad(i, j, k+1) = gh_ad(i, j, k+1) + temp1*temp_ad
          gh_ad(i+1, j, k) = gh_ad(i+1, j, k) - temp1*temp_ad
          pk_ad(i+1, j, k+1) = pk_ad(i+1, j, k+1) + temp0*temp_ad
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp0*temp_ad
          gh_ad(i, j, k) = gh_ad(i, j, k) + temp3*temp_ad
          gh_ad(i+1, j, k+1) = gh_ad(i+1, j, k+1) - temp3*temp_ad
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp2*temp_ad
          pk_ad(i+1, j, k) = pk_ad(i+1, j, k) - temp2*temp_ad
          wk_ad(i, j) = wk_ad(i, j) + temp_ad0
          wk_ad(i+1, j) = wk_ad(i+1, j) + temp_ad0
          delu_ad(i, j, k) = 0.0_8
        END DO
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8ARRAY(delp(:, :, k), (ied-isd+1)*(jed-jsd+1))
        CALL POPREAL8ARRAY(wk, (ied-isd+1)*(jed-jsd+1))
        CALL A2B_ORD2_ADM(delp(:, :, k), delp_ad(:, :, k), wk, wk_ad, &
&                   npx, npy, is, ie, js, je, ng)
      ELSE IF (branch .EQ. 1) THEN
        CALL POPREAL8ARRAY(delp(:, :, k), (ied-isd+1)*(jed-jsd+1))
        CALL POPREAL8ARRAY(wk, (ied-isd+1)*(jed-jsd+1))
        CALL A2B_ORD4_ADM(delp(:, :, k), delp_ad(:, :, k), wk, wk_ad, &
&                   npx, npy, is, ie, js, je, ng)
      ELSE
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPREAL8(wk(i, j))
            pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + wk_ad(i, j)
            pk_ad(i, j, k) = pk_ad(i, j, k) - wk_ad(i, j)
            wk_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
    END DO
    DO k=npz+1,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL A2B_ORD2_ADM(gh(isd:, jsd:, k), gh_ad(isd:, jsd:, k), wk, &
&                   wk_ad, npx, npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD4_ADM(gh(isd:, jsd:, k), gh_ad(isd:, jsd:, k), wk, &
&                   wk_ad, npx, npy, is, ie, js, je, ng, .true.)
      END IF
    END DO
    DO k=npz+1,2,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL A2B_ORD2_ADM(pk(isd:, jsd:, k), pk_ad(isd:, jsd:, k), wk, &
&                   wk_ad, npx, npy, is, ie, js, je, ng, .true.)
      ELSE
        CALL A2B_ORD4_ADM(pk(isd:, jsd:, k), pk_ad(isd:, jsd:, k), wk, &
&                   wk_ad, npx, npy, is, ie, js, je, ng, .true.)
      END IF
    END DO
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        pk_ad(i, j, 1) = 0.0_8
      END DO
    END DO
  END SUBROUTINE GRAD1_P_ADM
  SUBROUTINE GRAD1_P(delu, delv, pk, gh, delp, dt, ng, npx, npy, npz, &
&   ptop, ptk, hydrostatic)
    IMPLICIT NONE
! end k-loop
    INTEGER, INTENT(IN) :: ng, npx, npy, npz
    REAL, INTENT(IN) :: dt
    REAL(p_precision), INTENT(IN) :: ptk
    REAL(p_precision), INTENT(IN) :: ptop
    LOGICAL, INTENT(IN) :: hydrostatic
    REAL(pg_precision), INTENT(INOUT) :: pk(isd:ied, jsd:jed, npz+1)
    REAL(pg_precision), INTENT(INOUT) :: gh(isd:ied, jsd:jed, npz+1)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL(pg_precision), INTENT(OUT) :: delu(isd:ied, jsd:jed+1, npz)
    REAL(pg_precision), INTENT(OUT) :: delv(isd:ied+1, jsd:jed, npz)
! Local:
    REAL :: wk(isd:ied, jsd:jed)
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
        pk(i, j, 1) = top_value
      END DO
    END DO
    DO k=2,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(pk(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
      ELSE
        CALL A2B_ORD2(pk(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
      END IF
    END DO
    DO k=1,npz+1
      IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(gh(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
      ELSE
        CALL A2B_ORD2(gh(isd:, jsd:, k), wk, npx, npy, is, ie, js, je, &
&               ng, .true.)
      END IF
    END DO
    DO k=1,npz
      IF (hydrostatic) THEN
        DO j=js,je+1
          DO i=is,ie+1
            wk(i, j) = pk(i, j, k+1) - pk(i, j, k)
          END DO
        END DO
      ELSE IF (a2b_ord .EQ. 4) THEN
        CALL A2B_ORD4(delp(:, :, k), wk, npx, npy, is, ie, js, je, ng)
      ELSE
        CALL A2B_ORD2(delp(:, :, k), wk, npx, npy, is, ie, js, je, ng)
      END IF
      DO j=js,je+1
        DO i=is,ie
          delu(i, j, k) = rdx(i, j)*dt_pg/(wk(i, j)+wk(i+1, j))*((gh(i, &
&           j, k+1)-gh(i+1, j, k))*(pk(i+1, j, k+1)-pk(i, j, k))+(gh(i, &
&           j, k)-gh(i+1, j, k+1))*(pk(i, j, k+1)-pk(i+1, j, k)))
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          delv(i, j, k) = rdy(i, j)*dt_pg/(wk(i, j)+wk(i, j+1))*((gh(i, &
&           j, k+1)-gh(i, j+1, k))*(pk(i, j+1, k+1)-pk(i, j, k))+(gh(i, &
&           j, k)-gh(i, j+1, k+1))*(pk(i, j, k+1)-pk(i, j+1, k)))
        END DO
      END DO
    END DO
  END SUBROUTINE GRAD1_P
!  Differentiation of mix_dp in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: w delp pt
!   with respect to varying inputs: w delp pt
  SUBROUTINE MIX_DP_ADM(hydrostatic, w, w_ad, delp, delp_ad, pt, pt_ad, &
&   km, ak, bk, cg)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL, INTENT(IN) :: ak(km+1), bk(km+1)
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: pt, delp, w
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: pt_ad
    LOGICAL, INTENT(IN) :: hydrostatic, cg
! Local:
    REAL :: dp, dpmin
    REAL :: dp_ad
    INTEGER :: i, j, k, ip
    INTEGER :: ifirst, ilast
    INTEGER :: jfirst, jlast
    INTEGER :: branch
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: w_ad
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: delp_ad
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
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
      DO k=1,km-1
        CALL PUSHREAL8(dpmin)
        dpmin = 0.005*(ak(k+1)-ak(k)+(bk(k+1)-bk(k))*1.e5)
        DO i=ifirst,ilast
          IF (delp(i, j, k) .LT. dpmin) THEN
! Remap from below and mix pt
            CALL PUSHREAL8(dp)
            dp = dpmin - delp(i, j, k)
            CALL PUSHREAL8(pt(i, j, k))
            pt(i, j, k) = (pt(i, j, k)*delp(i, j, k)+pt(i, j, k+1)*dp)/&
&             dpmin
            IF (.NOT.hydrostatic) THEN
              CALL PUSHREAL8(w(i, j, k))
              w(i, j, k) = (w(i, j, k)*delp(i, j, k)+w(i, j, k+1)*dp)/&
&               dpmin
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREAL8(delp(i, j, k))
            delp(i, j, k) = dpmin
            CALL PUSHREAL8(delp(i, j, k+1))
            delp(i, j, k+1) = delp(i, j, k+1) - dp
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
! Bottom (k=km):
      CALL PUSHREAL8(dpmin)
      dpmin = 0.005*(ak(km+1)-ak(km)+(bk(km+1)-bk(km))*1.e5)
      DO i=ifirst,ilast
        IF (delp(i, j, km) .LT. dpmin) THEN
! Remap from above and mix pt
          CALL PUSHREAL8(dp)
          dp = dpmin - delp(i, j, km)
          IF (.NOT.hydrostatic) THEN
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=jlast,jfirst,-1
      dpmin = 0.005*(ak(km+1)-ak(km)+(bk(km+1)-bk(km))*1.e5)
      DO i=ilast,ifirst,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          dp_ad = -delp_ad(i, j, km-1)
          delp_ad(i, j, km) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            dp = dpmin - delp(i, j, km)
            temp_ad2 = w_ad(i, j, km)/dpmin
            delp_ad(i, j, km) = delp_ad(i, j, km) + w(i, j, km)*temp_ad2
            w_ad(i, j, km-1) = w_ad(i, j, km-1) + dp*temp_ad2
            dp_ad = dp_ad + w(i, j, km-1)*temp_ad2
            w_ad(i, j, km) = delp(i, j, km)*temp_ad2
          END IF
          temp_ad1 = pt_ad(i, j, km)/dpmin
          pt_ad(i, j, km-1) = pt_ad(i, j, km-1) + dp*temp_ad1
          dp_ad = dp_ad + pt(i, j, km-1)*temp_ad1
          delp_ad(i, j, km) = delp_ad(i, j, km) + pt(i, j, km)*temp_ad1 &
&           - dp_ad
          pt_ad(i, j, km) = delp(i, j, km)*temp_ad1
          CALL POPREAL8(dp)
        END IF
      END DO
      CALL POPREAL8(dpmin)
      DO k=km-1,1,-1
        DO i=ilast,ifirst,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            CALL POPREAL8(delp(i, j, k+1))
            dp_ad = -delp_ad(i, j, k+1)
            CALL POPREAL8(delp(i, j, k))
            delp_ad(i, j, k) = 0.0_8
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL8(w(i, j, k))
              temp_ad0 = w_ad(i, j, k)/dpmin
              delp_ad(i, j, k) = delp_ad(i, j, k) + w(i, j, k)*temp_ad0
              w_ad(i, j, k+1) = w_ad(i, j, k+1) + dp*temp_ad0
              dp_ad = dp_ad + w(i, j, k+1)*temp_ad0
              w_ad(i, j, k) = delp(i, j, k)*temp_ad0
            END IF
            CALL POPREAL8(pt(i, j, k))
            temp_ad = pt_ad(i, j, k)/dpmin
            pt_ad(i, j, k+1) = pt_ad(i, j, k+1) + dp*temp_ad
            dp_ad = dp_ad + pt(i, j, k+1)*temp_ad
            delp_ad(i, j, k) = delp_ad(i, j, k) + pt(i, j, k)*temp_ad - &
&             dp_ad
            pt_ad(i, j, k) = delp(i, j, k)*temp_ad
            CALL POPREAL8(dp)
          END IF
        END DO
        CALL POPREAL8(dpmin)
      END DO
    END DO
  END SUBROUTINE MIX_DP_ADM
  SUBROUTINE MIX_DP(hydrostatic, w, delp, pt, km, ak, bk, cg)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL, INTENT(IN) :: ak(km+1), bk(km+1)
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: pt, delp, w
    LOGICAL, INTENT(IN) :: hydrostatic, cg
! Local:
    REAL :: dp, dpmin
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
            dp = dpmin - delp(i, j, k)
            pt(i, j, k) = (pt(i, j, k)*delp(i, j, k)+pt(i, j, k+1)*dp)/&
&             dpmin
            IF (.NOT.hydrostatic) w(i, j, k) = (w(i, j, k)*delp(i, j, k)&
&               +w(i, j, k+1)*dp)/dpmin
            delp(i, j, k) = dpmin
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
          dp = dpmin - delp(i, j, km)
          pt(i, j, km) = (pt(i, j, km)*delp(i, j, km)+pt(i, j, km-1)*dp)&
&           /dpmin
          IF (.NOT.hydrostatic) w(i, j, km) = (w(i, j, km)*delp(i, j, km&
&             )+w(i, j, km-1)*dp)/dpmin
          delp(i, j, km) = dpmin
          delp(i, j, km-1) = delp(i, j, km-1) - dp
          ip = ip + 1
        END IF
      END DO
!      if ( fv_debug .and. ip/=0 ) write(*,*) 'Warning: Mix_dp', gid, j, ip 
      IF (ip .NE. 0) WRITE(*, *) 'Warning: Mix_dp', gid, j, ip
    END DO
  END SUBROUTINE MIX_DP
!  Differentiation of geopk in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: peln delp pkz pe pk pt gh
!   with respect to varying inputs: peln delp pkz pe pk pt gh
  SUBROUTINE GEOPK_ADM(ptop, pe, pe_ad, peln, peln_ad, delp, delp_ad, pk&
&   , pk_ad, gh, gh_ad, hs, pt, pt_ad, pkz, pkz_ad, km, akap, cg)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL, INTENT(IN) :: akap
    REAL(p_precision), INTENT(IN) :: ptop
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: pt, delp
    REAL, DIMENSION(isd:ied, jsd:jed, km) :: pt_ad, delp_ad
    LOGICAL, INTENT(IN) :: cg
! !OUTPUT PARAMETERS
    REAL(pg_precision), DIMENSION(isd:ied, jsd:jed, km + 1) :: gh
    REAL(pg_precision), DIMENSION(isd:ied, jsd:jed, km+1) :: gh_ad
    REAL(pg_precision), DIMENSION(isd:ied, jsd:jed, km + 1) :: pk
    REAL(pg_precision), DIMENSION(isd:ied, jsd:jed, km+1) :: pk_ad
    REAL(p_precision) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL(p_precision) :: pe_ad(is-1:ie+1, km+1, js-1:je+1)
! ln(pe)
    REAL(p_precision) :: peln(is:ie, km+1, js:je)
    REAL(p_precision) :: peln_ad(is:ie, km+1, js:je)
    REAL(p_precision) :: pkz(is:ie, js:je, km)
    REAL(p_precision) :: pkz_ad(is:ie, js:je, km)
! !DESCRIPTION:
!    Calculates geopotential and pressure to the kappa.
! Local:
    REAL(pg_precision) :: p1d(is-2:ie+2)
    REAL(pg_precision) :: p1d_ad(is-2:ie+2)
    REAL(pg_precision) :: logp(is-2:ie+2)
    REAL(pg_precision) :: logp_ad(is-2:ie+2)
    REAL(pg_precision) :: ptk
    INTEGER :: i, j, k
    INTEGER :: ifirst, ilast
    INTEGER :: jfirst, jlast
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC LOG
    INTRINSIC EXP
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_from0
    INTEGER :: ad_to0
    INTEGER :: branch
    INTEGER :: min2
    INTEGER :: min1
    REAL(p_precision) :: temp_ad1
    REAL*8 :: temp_ad0
    REAL(pg_precision) :: temp_ad
    REAL*8 :: temp
    INTEGER :: max2
    INTEGER :: max1
    ptk = ptop**akap
    IF (.NOT.cg .AND. a2b_ord .EQ. 4) THEN
! D-Grid
      ifirst = is - 2
      ilast = ie + 2
      jfirst = js - 2
      jlast = je + 2
    ELSE
      ifirst = is - 1
      ilast = ie + 1
      jfirst = js - 1
      jlast = je + 1
    END IF
!$omp parallel do default(shared) private(i, j, k, p1d, dp, logp)
    DO j=jfirst,jlast
      DO i=ifirst,ilast
        CALL PUSHREAL8ARRAY(p1d(i), pg_precision/8)
        p1d(i) = ptop
        CALL PUSHREAL8ARRAY(pk(i, j, 1), pg_precision/8)
        pk(i, j, 1) = ptk
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
        ad_from = max1
        i = min1 + 1
        CALL PUSHINTEGER4(i - 1)
        CALL PUSHINTEGER4(ad_from)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
! Top down
      DO k=2,km+1
        DO i=ifirst,ilast
          CALL PUSHREAL8ARRAY(p1d(i), pg_precision/8)
          p1d(i) = p1d(i) + delp(i, j, k-1)
!            pk(i,j,k) = p1d(i) ** akap
! Optimized form:
          CALL PUSHREAL8ARRAY(logp(i), pg_precision/8)
          logp(i) = LOG(p1d(i))
          CALL PUSHREAL8ARRAY(pk(i, j, k), pg_precision/8)
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
          ad_from0 = max2
          i = min2 + 1
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from0)
          IF (j .GE. js .AND. j .LE. je) THEN
            DO i=is,ie
              CALL PUSHREAL8(peln(i, k, j))
              peln(i, k, j) = logp(i)
            END DO
            CALL PUSHCONTROL2B(2)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(0)
        END IF
      END DO
    END DO
    IF (.NOT.cg) THEN
      DO k=km,1,-1
        DO j=je,js,-1
          DO i=ie,is,-1
            temp = akap*(peln(i, k+1, j)-peln(i, k, j))
            temp_ad0 = pkz_ad(i, j, k)/temp
            temp_ad1 = -((pk(i, j, k+1)-pk(i, j, k))*akap*temp_ad0/temp)
            pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp_ad0
            pk_ad(i, j, k) = pk_ad(i, j, k) - temp_ad0
            peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad1
            peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad1
            pkz_ad(i, j, k) = 0.0_8
          END DO
        END DO
      END DO
    END IF
    logp_ad = 0.0_8
    p1d_ad = 0.0_8
    DO j=jlast,jfirst,-1
      DO k=1,km,1
        DO i=ilast,ifirst,-1
          temp_ad = pt(i, j, k)*gh_ad(i, j, k)
          gh_ad(i, j, k+1) = gh_ad(i, j, k+1) + gh_ad(i, j, k)
          pt_ad(i, j, k) = pt_ad(i, j, k) + (pk(i, j, k+1)-pk(i, j, k))*&
&           gh_ad(i, j, k)
          pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp_ad
          pk_ad(i, j, k) = pk_ad(i, j, k) - temp_ad
          gh_ad(i, j, k) = 0.0_8
        END DO
      END DO
      DO k=km+1,2,-1
        CALL POPCONTROL2B(branch)
        IF (branch .NE. 0) THEN
          IF (branch .NE. 1) THEN
            DO i=ie,is,-1
              CALL POPREAL8(peln(i, k, j))
              logp_ad(i) = logp_ad(i) + peln_ad(i, k, j)
              peln_ad(i, k, j) = 0.0_8
            END DO
          END IF
          CALL POPINTEGER4(ad_from0)
          CALL POPINTEGER4(ad_to0)
          DO i=ad_to0,ad_from0,-1
            p1d_ad(i) = p1d_ad(i) + pe_ad(i, k, j)
            pe_ad(i, k, j) = 0.0_8
          END DO
        END IF
        DO i=ilast,ifirst,-1
          CALL POPREAL8ARRAY(pk(i, j, k), pg_precision/8)
          logp_ad(i) = logp_ad(i) + EXP(akap*logp(i))*akap*pk_ad(i, j, k&
&           )
          pk_ad(i, j, k) = 0.0_8
          CALL POPREAL8ARRAY(logp(i), pg_precision/8)
          p1d_ad(i) = p1d_ad(i) + logp_ad(i)/p1d(i)
          logp_ad(i) = 0.0_8
          CALL POPREAL8ARRAY(p1d(i), pg_precision/8)
          delp_ad(i, j, k-1) = delp_ad(i, j, k-1) + p1d_ad(i)
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        CALL POPINTEGER4(ad_from)
        CALL POPINTEGER4(ad_to)
        DO i=ad_to,ad_from,-1
          pe_ad(i, 1, j) = 0.0_8
        END DO
      END IF
      DO i=ilast,ifirst,-1
        gh_ad(i, j, km+1) = 0.0_8
        CALL POPREAL8ARRAY(pk(i, j, 1), pg_precision/8)
        pk_ad(i, j, 1) = 0.0_8
        CALL POPREAL8ARRAY(p1d(i), pg_precision/8)
        p1d_ad(i) = 0.0_8
      END DO
    END DO
  END SUBROUTINE GEOPK_ADM
  SUBROUTINE GEOPK(ptop, pe, peln, delp, pk, gh, hs, pt, pkz, km, akap, &
&   cg)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL, INTENT(IN) :: akap
    REAL(p_precision), INTENT(IN) :: ptop
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: pt, delp
    LOGICAL, INTENT(IN) :: cg
! !OUTPUT PARAMETERS
    REAL(pg_precision), DIMENSION(isd:ied, jsd:jed, km + 1), INTENT(OUT)&
&   :: gh
    REAL(pg_precision), DIMENSION(isd:ied, jsd:jed, km + 1), INTENT(OUT)&
&   :: pk
    REAL(p_precision), INTENT(OUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
! ln(pe)
    REAL(p_precision), INTENT(OUT) :: peln(is:ie, km+1, js:je)
    REAL(p_precision), INTENT(OUT) :: pkz(is:ie, js:je, km)
! !DESCRIPTION:
!    Calculates geopotential and pressure to the kappa.
! Local:
    REAL(pg_precision) :: p1d(is-2:ie+2)
    REAL(pg_precision) :: logp(is-2:ie+2)
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
    ELSE
      ifirst = is - 1
      ilast = ie + 1
      jfirst = js - 1
      jlast = je + 1
    END IF
!$omp parallel do default(shared) private(i, j, k, p1d, dp, logp)
    DO j=jfirst,jlast
      DO i=ifirst,ilast
        p1d(i) = ptop
        pk(i, j, 1) = ptk
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
          pe(i, 1, j) = ptop
        END DO
      END IF
! Top down
      DO k=2,km+1
        DO i=ifirst,ilast
          p1d(i) = p1d(i) + delp(i, j, k-1)
!            pk(i,j,k) = p1d(i) ** akap
! Optimized form:
          logp(i) = LOG(p1d(i))
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
            pe(i, k, j) = p1d(i)
          END DO
          IF (j .GE. js .AND. j .LE. je) THEN
            DO i=is,ie
              peln(i, k, j) = logp(i)
            END DO
          END IF
        END IF
      END DO
! Bottom up
      DO k=km,1,-1
        DO i=ifirst,ilast
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
            pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(akap*(peln(i, k+&
&             1, j)-peln(i, k, j)))
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE GEOPK
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

end module dyn_core_adm_mod
