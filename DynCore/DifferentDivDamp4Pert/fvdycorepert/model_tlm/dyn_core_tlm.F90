module dyn_core_tlm_mod

  use fv_arrays_mod,      only: g_precision, p_precision, pg_precision  
  use mpp_domains_mod,    only: CGRID_NE, DGRID_NE, mpp_get_boundary,   &
                                mpp_update_domains
  use fv_mp_mod,          only: domain, isd, ied, jsd, jed, is, ie, js, je, gid, mp_barrier, mp_reduce_max
  use fv_control_mod,     only: hord_mt, hord_vt, hord_tm, hord_dp, hord_ze, hord_tr, &
                                hord_mt_pert, hord_vt_pert, hord_tm_pert, hord_dp_pert, hord_tr_pert, &
                                n_sponge, m_riem, m_split, &
                                dddmp, dddmp_pert, fv_d2_bg=>d2_bg, fv_d4_bg=>d4_bg, d_ext, vtdm4, beta1, beta, init_wind_m, m_grad_p, &
                                fv_d2_bg_pert=>d2_bg_pert, fv_d4_bg_pert=>d4_bg_pert, &
                                a2b_ord, ppm_limiter, master, fv_debug, d_con, fv_nord=>nord, fv_nord_pert=>nord_pert, max_courant_no, &
                                no_cgrid, fill_dp, nwat, inline_q, breed_vortex_inline, shallow_water, spin_up_hours
  use sw_core_tlm_mod,    only: c_sw_tlm, d_sw_tlm, divergence_corner_tlm, d2a2c_tlm, divergence_corner_onlytlm
  use sw_core_mod,        only: divergence_corner
  use a2b_edge_tlm_mod,   only: a2b_ord2_tlm, a2b_ord4_tlm
  use fv_grid_tools_mod,  only: rdx, rdy, rdxc, dxc, dyc, rdyc, dx, dy, area, rarea, grid_type
  use fv_grid_utils_mod,  only: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n,  &
                                ec1, ec2, en1, en2, da_min_c
  use mpp_parameter_mod,  only: CORNER

implicit none

private
public :: dyn_core_tlm

contains

  SUBROUTINE DYN_CORE_TLM(npx, npy, npz, ng, sphum, nq, bdt, n_split, &
&   zvir, cp, akap, rdgas, grav, hydrostatic, u, u_tl, v, v_tl, um, &
&   um_tl, vm, vm_tl, w, w_tl, delz, pt, pt_tl, q, q_tl, delp, delp_tl, &
&   pe, pe_tl, pk, pk_tl, phis, omga, omga_tl, ptop, pfull, ua, ua_tl, &
&   va, va_tl, uc, uc_tl, vc, vc_tl, mfx, mfx_tl, mfy, mfy_tl, cx, cx_tl&
&   , cy, cy_tl, pem, pem_tl, pkz, pkz_tl, peln, peln_tl, ak, bk, &
&   init_step, end_step, proc_id, time_total)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, proc_id
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
    INTEGER, INTENT(IN) :: ng, nq, sphum
    INTEGER, INTENT(IN) :: n_split
    REAL, INTENT(IN) :: bdt
    REAL, INTENT(IN) :: zvir, cp, akap, rdgas, grav
    REAL(p_precision), INTENT(IN) :: ptop
    REAL, DIMENSION(npz + 1), INTENT(IN) :: ak, bk
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: init_step, end_step
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
! total time (seconds) since start
    REAL, INTENT(IN), OPTIONAL :: time_total
!THESE ARE ACTUALLY IN THE MODULE
    REAL, ALLOCATABLE :: delzc(:, :, :)
    REAL, ALLOCATABLE :: ut(:, :, :)
    REAL, ALLOCATABLE :: ut_tl(:, :, :)
    REAL, ALLOCATABLE :: vt(:, :, :)
    REAL, ALLOCATABLE :: vt_tl(:, :, :)
    REAL, ALLOCATABLE :: crx(:, :, :), xfx(:, :, :)
    REAL, ALLOCATABLE :: crx_tl(:, :, :), xfx_tl(:, :, :)
    REAL, ALLOCATABLE :: cry(:, :, :), yfx(:, :, :)
    REAL, ALLOCATABLE :: cry_tl(:, :, :), yfx_tl(:, :, :)
    REAL, ALLOCATABLE :: divg_d(:, :, :)
    REAL, ALLOCATABLE :: divg_d_tl(:, :, :)
    REAL, ALLOCATABLE :: zh(:, :, :)
    REAL, ALLOCATABLE :: zh_tl(:, :, :)
    REAL(pg_precision), ALLOCATABLE :: du(:, :, :), dv(:, :, :)
    REAL(pg_precision), ALLOCATABLE :: du_tl(:, :, :), dv_tl(:, :, :)
    REAL(pg_precision), ALLOCATABLE :: pkc(:, :, :)
    REAL(pg_precision), ALLOCATABLE :: pkc_tl(:, :, :)
    REAL, ALLOCATABLE :: pkd(:, :, :)
    REAL, ALLOCATABLE :: pkd_tj(:, :, :)
    REAL, ALLOCATABLE :: pkd_tl(:, :, :)
    REAL, ALLOCATABLE :: delpc(:, :, :)
    REAL, ALLOCATABLE :: delpc_tl(:, :, :)
    REAL(pg_precision), ALLOCATABLE :: pk3(:, :, :)
    REAL(pg_precision), ALLOCATABLE :: pk3_tl(:, :, :)
    REAL, ALLOCATABLE :: ptc(:, :, :)
    REAL, ALLOCATABLE :: ptc_tj(:, :, :)
    REAL, ALLOCATABLE :: ptc_tl(:, :, :)
    REAL(pg_precision), ALLOCATABLE :: gz(:, :, :)
    REAL(pg_precision), ALLOCATABLE :: gz_tl(:, :, :)
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
    REAL, INTENT(OUT) :: pem_tl(is-1:ie+1, npz+1, js-1:je+1)
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
    REAL(pg_precision) :: divg2_tj(is:ie+1, js:je+1)
    REAL(pg_precision) :: divg2_tl(is:ie+1, js:je+1)
    REAL :: wk(isd:ied, jsd:jed)
    REAL :: wk_tj(isd:ied, jsd:jed)
    REAL :: wk_tl(isd:ied, jsd:jed)
! --- For no_cgrid option ---
    REAL :: u2(isd:ied, jsd:jed+1)
    REAL :: u2_tl(isd:ied, jsd:jed+1)
    REAL :: v2(isd:ied+1, jsd:jed)
    REAL :: v2_tl(isd:ied+1, jsd:jed)
!-------------------------------------
    REAL :: d4_bg, d2_bg
    REAL :: d4_bg_pert, d2_bg_pert
    INTEGER :: hord_m, hord_v, hord_t, hord_p, nord, nord_k
    INTEGER :: hord_m_pert, hord_v_pert, hord_t_pert, hord_p_pert, nord_pert, nord_k_pert
    INTEGER :: i, j, k, it, iq
    INTEGER :: ism1, iep1, jsm1, jep1
    INTEGER :: ieb1, jeb1
    REAL :: alpha, damp_k, damp_k_pert
    REAL :: dt2, dt, rdt, rgrav
    REAL :: d2_divg, d2_divg_pert, dd_divg, dd_divg_pert
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
    INTRINSIC LOG
    INTRINSIC TANH
    INTRINSIC MIN
    REAL :: arg1
    REAL :: arg2
    REAL :: abs1
    REAL :: abs0
    REAL :: y5, y5_pert
    INTEGER :: y4
    REAL :: y3
    REAL :: y2
    INTEGER :: y1

    integer :: ioffset, joffset, ii, jj, kk, idh, jdh, kdh


    d2_bg = fv_d2_bg
    d4_bg = fv_d4_bg
    d4_bg_pert = fv_d4_bg_pert
    d2_bg_pert = fv_d2_bg_pert
    nord = fv_nord
    nord_pert = fv_nord_pert
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
if ( npz>1 ) then
call mpp_update_domains( pt_tl,   domain, complete=.false. )
call mpp_update_domains( pt,      domain, complete=.false. )
endif
call mpp_update_domains( delp_tl, domain, complete=.true. )
call mpp_update_domains( delp, domain, complete=.true. )
call mpp_update_domains( u_tl, v_tl, domain, gridtype=DGRID_NE, complete=.true. )
call mpp_update_domains( u, v, domain, gridtype=DGRID_NE, complete=.true. )



    IF (init_step) THEN
      ALLOCATE(gz_tl(isd:ied, jsd:jed, npz+1))
      ALLOCATE(gz(isd:ied, jsd:jed, npz+1))
      ALLOCATE(ptc_tl(isd:ied, jsd:jed, npz))
      ALLOCATE(ptc(isd:ied, jsd:jed, npz))
      ALLOCATE(ptc_tj(isd:ied, jsd:jed, npz))
      ALLOCATE(delzc(is:ie, js:je, npz))
      ALLOCATE(crx_tl(is:ie+1, jsd:jed, npz))
      ALLOCATE(crx(is:ie+1, jsd:jed, npz))
      ALLOCATE(xfx_tl(is:ie+1, jsd:jed, npz))
      ALLOCATE(xfx(is:ie+1, jsd:jed, npz))
      ALLOCATE(cry_tl(isd:ied, js:je+1, npz))
      ALLOCATE(cry(isd:ied, js:je+1, npz))
      ALLOCATE(yfx_tl(isd:ied, js:je+1, npz))
      ALLOCATE(yfx(isd:ied, js:je+1, npz))
      ALLOCATE(divg_d_tl(isd:ied+1, jsd:jed+1, npz))
      ALLOCATE(divg_d(isd:ied+1, jsd:jed+1, npz))
      ALLOCATE(pkc_tl(isd:ied, jsd:jed, npz+1))
      ALLOCATE(pkc(isd:ied, jsd:jed, npz+1))
      ALLOCATE(delpc_tl(isd:ied, jsd:jed, npz))
      ALLOCATE(delpc(isd:ied, jsd:jed, npz))
      ALLOCATE(pkd_tl(isd:ied, jsd:jed, npz+1))
      ALLOCATE(pkd(isd:ied, jsd:jed, npz+1))
      ALLOCATE(pkd_tj(isd:ied, jsd:jed, npz+1))
      IF (.NOT.no_cgrid) THEN
        ALLOCATE(ut_tl(isd:ied, jsd:jed, npz))
        ALLOCATE(ut(isd:ied, jsd:jed, npz))
        ALLOCATE(vt_tl(isd:ied, jsd:jed, npz))
        ALLOCATE(vt(isd:ied, jsd:jed, npz))
        ut_tl(:, :, :) = 0.0
        ut(:, :, :) = 0.
        vt_tl(:, :, :) = 0.0
        vt(:, :, :) = 0.
      END IF
      IF (.NOT.hydrostatic) THEN
        ALLOCATE(zh_tl(isd:ied, jsd:jed, npz))
        ALLOCATE(zh(isd:ied, jsd:jed, npz))
        IF (m_grad_p .EQ. 0) THEN
          ALLOCATE(pk3_tl(isd:ied, jsd:jed, npz+1))
          ALLOCATE(pk3(isd:ied, jsd:jed, npz+1))
        END IF
      END IF
      IF (beta .GT. 1.e-4 .OR. no_cgrid) THEN
        ALLOCATE(du_tl(isd:ied, jsd:jed+1, npz))
        ALLOCATE(du(isd:ied, jsd:jed+1, npz))
        ALLOCATE(dv_tl(isd:ied+1, jsd:jed, npz))
        ALLOCATE(dv(isd:ied+1, jsd:jed, npz))
      END IF
    END IF
    IF (beta .GT. 1.e-4 .OR. (no_cgrid .AND. init_wind_m)) THEN
      CALL GEOPK_TLM(ptop, pe, pe_tl, peln, peln_tl, delp, delp_tl, pkc&
&              , pkc_tl, gz, gz_tl, phis, pt, pt_tl, pkz, pkz_tl, npz, &
&              akap, .false.)
      CALL GRAD1_P_TLM(du, du_tl, dv, dv_tl, pkc, pkc_tl, gz, gz_tl, &
&                delp, delp_tl, dt, ng, npx, npy, npz, ptop, ptk, &
&                hydrostatic)
if( init_wind_m ) then
call mpp_update_domains(du_tl, dv_tl, domain, gridtype=DGRID_NE)
call mpp_update_domains(du,    dv,    domain, gridtype=DGRID_NE)
endif
    END IF
! Empty the "flux capacitors"
    mfx(:, :, :) = 0.
    mfy(:, :, :) = 0.
    cx(:, :, :) = 0.
    cy(:, :, :) = 0.
    elapsed_dt = 0.0
    subcycle = 1
    IF (.NOT.hydrostatic) THEN
      IF (dt .GE. 0.) THEN
        abs0 = dt
      ELSE
        abs0 = -dt
      END IF
      y1 = CEILING(0.5 + abs0/3.0)
      IF (1 .LT. y1) THEN
        m_split = y1
      ELSE
        m_split = 1
      END IF
      mfx_tl = 0.0
      mfy_tl = 0.0
      cx_tl = 0.0
      cy_tl = 0.0
      v2_tl = 0.0
      u2_tl = 0.0
      wk_tl = 0.0
      divg2_tl = 0.0
    ELSE
      mfx_tl = 0.0
      mfy_tl = 0.0
      cx_tl = 0.0
      cy_tl = 0.0
      v2_tl = 0.0
      u2_tl = 0.0
      wk_tl = 0.0
      divg2_tl = 0.0
    END IF
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
                y2 = crx(i, j, k)
              ELSE
                y2 = -crx(i, j, k)
              END IF
              IF (cmax .LT. y2) THEN
                cmax = y2
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
                y3 = cry(i, j, k)
              ELSE
                y3 = -cry(i, j, k)
              END IF
              IF (cmax .LT. y3) THEN
                cmax = y3
              ELSE
                cmax = cmax
              END IF
            END DO
          END DO
        END DO
        call mp_reduce_max(cmax)
        call mp_barrier()
!NINTTOCEILING
        subt = CEILING(cmax/max_courant_no)
        IF (subcycle .NE. subt) THEN
          subcycle = subt
          dt = (bdt-elapsed_dt)/REAL((n_split-(it-1))*subcycle)
          dt2 = 0.5*dt
          rdt = 1./dt
!       if (gid==0) write(*,101) 'Sub-Cycling for stability = ', it, subcycle, dt, cmax
!101    format(A,i3,2x,i3,2x,f10.6,2x,f10.6)
          IF (.NOT.hydrostatic) THEN
            IF (dt .GE. 0.) THEN
              abs1 = dt
            ELSE
              abs1 = -dt
            END IF
            y4 = CEILING(0.5 + abs1/3.0)
            IF (1 .LT. y4) THEN
              m_split = y4
            ELSE
              m_split = 1
            END IF
          END IF
        END IF
      END IF
      DO subt=1,subcycle
!--------------------------------------------
        IF (.NOT.hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              zh(i, j, npz) = phis(i, j)*rgrav - delz(i, j, npz)
            END DO
            DO k=npz-1,1,-1
              DO i=is,ie
                zh(i, j, k) = zh(i, j, k+1) - delz(i, j, k)
              END DO
            END DO
          END DO
          call mpp_update_domains(zh_tl, domain, complete=.false.)
          call mpp_update_domains(zh, domain, complete=.false.)
          call mpp_update_domains(w_tl,  domain, complete=.true.)
          call mpp_update_domains(w,  domain, complete=.true.)
        END IF
        IF (shallow_water) THEN
          do_omega = .false.
        ELSE IF (it .EQ. n_split .AND. subt .EQ. subcycle) THEN
          DO j=jsm1,jep1
            DO i=ism1,iep1
              pem_tl(i, 1, j) = 0.0
              pem(i, 1, j) = ptop
            END DO
            DO k=1,npz
              DO i=ism1,iep1
                pem_tl(i, k+1, j) = pem_tl(i, k, j) + delp_tl(i, j, k)
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
&                    .GT. 0, nord_pert .GT. 0 )
          END DO
          IF (.NOT.hydrostatic) THEN
            delpc_tl(:, :, :) = delp_tl(:, :, :)
            delpc(:, :, :) = delp(:, :, :)
          END IF
          um_tl(:, :, :) = u_tl(:, :, :)
          um(:, :, :) = u(:, :, :)
          vm_tl(:, :, :) = v_tl(:, :, :)
          vm(:, :, :) = v(:, :, :)
        ELSE
          DO k=1,npz
            CALL C_SW_TLM(delpc(isd, jsd, k), delpc_tl(isd, jsd, k), &
&                   delp(isd, jsd, k), delp_tl(isd, jsd, k), ptc(isd, &
&                   jsd, k), ptc_tl(isd, jsd, k), pt(isd, jsd, k), pt_tl&
&                   (isd, jsd, k), u(isd, jsd, k), u_tl(isd, jsd, k), v(&
&                   isd, jsd, k), v_tl(isd, jsd, k), w(isd, jsd, k), &
&                   w_tl(isd, jsd, k), uc(isd, jsd, k), uc_tl(isd, jsd, &
&                   k), vc(isd, jsd, k), vc_tl(isd, jsd, k), ua(isd, jsd&
&                   , k), ua_tl(isd, jsd, k), va(isd, jsd, k), va_tl(isd&
&                   , jsd, k), omga(isd, jsd, k), omga_tl(isd, jsd, k), &
&                   ut(isd, jsd, k), ut_tl(isd, jsd, k), vt(isd, jsd, k)&
&                   , vt_tl(isd, jsd, k), dt2, hydrostatic, nord .GT. 0, &
                    nord_pert .GT. 0 )
          END DO


! on output omga is updated w
          IF (fill_dp) CALL MIX_DP_TLM(hydrostatic, omga, omga_tl, delpc&
&                                , delpc_tl, ptc, ptc_tl, npz, ak, bk, &
&                                .true.)
          IF (hydrostatic) THEN
            IF (beta1 .GT. 0.001) THEN
              alpha = 1. - beta1
              delpc_tl(:, :, :) = beta1*delp_tl(:, :, :) + alpha*&
&               delpc_tl(:, :, :)
              delpc(:, :, :) = beta1*delp(:, :, :) + alpha*delpc(:, :, :&
&               )
              ptc_tl(:, :, :) = beta1*pt_tl(:, :, :) + alpha*ptc_tl(:, :&
&               , :)
              ptc(:, :, :) = beta1*pt(:, :, :) + alpha*ptc(:, :, :)
            END IF
            CALL GEOPK_TLM(ptop, pe, pe_tl, peln, peln_tl, delpc, &
&                    delpc_tl, pkc, pkc_tl, gz, gz_tl, phis, ptc, ptc_tl&
&                    , pkz, pkz_tl, npz, akap, .true.)
          END IF
!DELETED NH_CORE CALLS -dh
!-----------------------------------------
! Update time-centered winds on the C-Grid
!-----------------------------------------
          ieb1 = ie + 1
          jeb1 = je + 1
          DO k=1,npz
            IF (hydrostatic) THEN
              DO j=jsm1,jeb1
                DO i=ism1,ieb1
                  wk_tl(i, j) = pkc_tl(i, j, k+1) - pkc_tl(i, j, k)
                  wk(i, j) = pkc(i, j, k+1) - pkc(i, j, k)
                END DO
              END DO
            ELSE
              DO j=jsm1,jeb1
                DO i=ism1,ieb1
                  wk_tl(i, j) = delpc_tl(i, j, k)
                  wk(i, j) = delpc(i, j, k)
                END DO
              END DO
! Save delp for update_dz_d
              delpc_tl(:, :, k) = delp_tl(:, :, k)
              delpc(:, :, k) = delp(:, :, k)
            END IF
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
        call mpp_update_domains(uc_tl, vc_tl, domain, gridtype=CGRID_NE, complete=.true.)
        call mpp_update_domains(uc, vc, domain, gridtype=CGRID_NE, complete=.true.)

        IF (inline_q) then
          call mpp_update_domains( q_tl,  domain, complete=.true. )
          call mpp_update_domains( q,  domain, complete=.true. )
        endif
        IF (nord .GT. 0) THEN
          CALL DIVERGENCE_CORNER(u, v, ua, va, divg_d, npz)
          call mpp_update_domains(divg_d, domain, position=CORNER)
        END IF
        IF (nord_pert .GT. 0) THEN
          CALL DIVERGENCE_CORNER_ONLYTLM(u_tl, v_tl, ua_tl, va_tl, divg_d_tl, npz)
          call mpp_update_domains(divg_d_tl, domain, position=CORNER)
        END IF
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
          nord_k_pert = nord_pert
          damp_k = dddmp
          damp_k_pert = dddmp_pert

          arg1 = pfull(k)/pfull(npz)
          arg2 = 0.1*LOG(arg1)
          y5 = d2_bg*(1.-3.*TANH(arg2))
          y5_pert = d2_bg_pert*(1.-3.*TANH(arg2))
          IF (0.20 .GT. y5) THEN
            d2_divg = y5
          ELSE
            d2_divg = 0.20
          END IF
           IF (0.20 .GT. y5_pert) THEN
             d2_divg_pert = y5_pert
           ELSE
             d2_divg_pert = 0.20
           END IF
          IF (n_sponge .EQ. -1 .OR. npz .EQ. 1) THEN
! Constant divg damping coefficient:
            d2_divg = d2_bg
            d2_divg_pert = d2_bg_pert
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
             IF (0 .LT. nord_pert - 1) THEN
               nord_k_pert = nord_pert - 1
             ELSE
               nord_k_pert = 0
             END IF
            damp_k = 0.2
            damp_k_pert = 0.2
            IF (0.20 .GT. 2.*d2_bg) THEN
              d2_divg = 2.*d2_bg
            ELSE
              d2_divg = 0.20
            END IF
             IF (0.20 .GT. 2.*d2_bg_pert) THEN
               d2_divg_pert = 2.*d2_bg_pert
             ELSE
               d2_divg_pert = 0.20
             END IF
            IF (0.06 .LT. d2_divg) THEN
              d2_divg = d2_divg
            ELSE
              d2_divg = 0.06
            END IF
             IF (0.06 .LT. d2_divg_pert) THEN
               d2_divg_pert = d2_divg_pert
             ELSE
               d2_divg_pert = 0.06
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
            nord_k_pert = 0
            damp_k = 0.2
            damp_k_pert = 0.2

            IF (0.20 .GT. 4.*d2_bg) THEN
              d2_divg = 4.*d2_bg
            ELSE
              d2_divg = 0.20
            END IF
             IF (0.20 .GT. 4.*d2_bg_pert) THEN
               d2_divg_pert = 4.*d2_bg_pert
             ELSE
               d2_divg_pert = 0.20
             END IF
            IF (0.15 .LT. d2_divg) THEN
              d2_divg = d2_divg
            ELSE
              d2_divg = 0.15
            END IF
             IF (0.15 .LT. d2_divg_pert) THEN
               d2_divg_pert = d2_divg_pert
             ELSE
               d2_divg_pert = 0.15
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
             IF (0 .LT. nord_pert - 1) THEN
               nord_k_pert = nord_pert - 1
             ELSE
               nord_k_pert = 0
             END IF
            IF (0.20 .GT. 2.*d2_bg) THEN
              d2_divg = 2.*d2_bg
            ELSE
              d2_divg = 0.20
            END IF
             IF (0.20 .GT. 2.*d2_bg_pert) THEN
               d2_divg_pert = 2.*d2_bg_pert
             ELSE
               d2_divg_pert = 0.20
             END IF
            IF (0.08 .LT. d2_divg) THEN
              d2_divg = d2_divg
            ELSE
              d2_divg = 0.08
            END IF
             IF (0.08 .LT. d2_divg_pert) THEN
               d2_divg_pert = d2_divg_pert
             ELSE
               d2_divg_pert = 0.08
             END IF
            IF (nord .GT. 1) THEN
              damp_k = 0.
            ELSE
              damp_k = 0.12
            END IF
             IF (nord_pert .GT. 1) THEN
               damp_k_pert = 0.
             ELSE
               damp_k_pert = 0.12
             END IF
          END IF
          dd_divg = d4_bg
          dd_divg_pert = d4_bg_pert
          IF (damp_k .LT. dddmp) THEN
            damp_k = dddmp
          ELSE
            damp_k = damp_k
          END IF
           IF (damp_k_pert .LT. dddmp_pert) THEN
             damp_k_pert = dddmp_pert
           ELSE
             damp_k_pert = damp_k_pert
           END IF
! if (gid==0 .and. first_step) then
!    print*, k, nord_k, damp_k, d2_divg, dd_divg 
! endif
!--- external mode divergence damping ---
          IF (d_ext .GT. 0.) CALL A2B_ORD2_TLM(delp(:, :, k), delp_tl(:&
&                                        , :, k), wk, wk_tl, npx, npy, &
&                                        is, ie, js, je, ng, .false.)
          CALL D_SW_TLM(pkd(isd, jsd, k), pkd_tj(isd, jsd, k), &
&                 pkd_tl(isd, jsd, k), delp(isd&
&                 , jsd, k), delp_tl(isd, jsd, k), ptc(isd, jsd, k), ptc_tj(isd, jsd, k), &
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
&                 hord_v_pert, hord_t_pert, hord_p_pert, nord_k, nord_k_pert, &
&                 damp_k, damp_k_pert, d2_divg, d2_divg_pert, &
&                 dd_divg, dd_divg_pert, vtdm4, d_con, hydrostatic, ppm_limiter,proc_id,it)

          IF (d_ext .GT. 0.) THEN
            DO j=js,jep1
              DO i=is,iep1
! delp at cell corners
                ptc_tl(i, j, k) = wk_tl(i, j)
                ptc_tj(i, j, k) = wk(i, j)
                ptc(i, j, k) = wk(i, j)
              END DO
            END DO
          END IF
        END DO
        IF (fill_dp) CALL MIX_DP_TLM(hydrostatic, w, w_tl, delp, delp_tl&
&                              , pt, pt_tl, npz, ak, bk, .false.)
        IF (d_ext .GT. 0.) THEN
          d2_divg = d_ext*da_min_c
! pkd() is 3D field of horizontal divergence
! ptc is "delp" at cell corners
          DO j=js,jep1
            DO i=is,iep1
              wk_tl(i, j) = ptc_tl(i, j, 1)
              wk_tj(i, j) = ptc_tj(i, j, 1)
              wk(i, j) = ptc(i, j, 1)
              divg2_tl(i, j) = wk_tl(i, j)*pkd_tj(i, j, 1) + wk_tj(i, j)*&
&               pkd_tl(i, j, 1)
              divg2_tj(i, j) = wk_tj(i, j)*pkd_tj(i, j, 1)
              divg2(i, j) = wk(i, j)*pkd(i, j, 1)
            END DO
            DO k=2,npz
              DO i=is,iep1
                wk_tl(i, j) = wk_tl(i, j) + ptc_tl(i, j, k)
                wk_tj(i, j) = wk_tj(i, j) + ptc_tj(i, j, k)
                wk(i, j) = wk(i, j) + ptc(i, j, k)
                divg2_tl(i, j) = divg2_tl(i, j) + ptc_tl(i, j, k)*pkd_tj(i&
&                 , j, k) + ptc_tj(i, j, k)*pkd_tl(i, j, k)
                divg2_tj(i, j) = divg2_tj(i, j) + ptc_tj(i, j, k)*pkd_tj(i, j, k)
                divg2(i, j) = divg2(i, j) + ptc(i, j, k)*pkd(i, j, k)
              END DO
            END DO
            DO i=is,iep1
              divg2_tl(i, j) = (d2_divg*divg2_tl(i, j)*wk_tj(i, j)-d2_divg*&
&               divg2_tj(i, j)*wk_tl(i, j))/wk_tj(i, j)**2
              divg2_tj(i, j) = d2_divg*divg2_tj(i, j)/wk_tj(i, j)
              divg2(i, j) = d2_divg*divg2(i, j)/wk(i, j)
            END DO
          END DO
        ELSE
          divg2 = 0.
          divg2_tj = 0.
          divg2_tl = 0.0
        END IF
        call mpp_update_domains(  pt_tl, domain, complete=.false.)
        call mpp_update_domains(  pt, domain, complete=.false.)
        call mpp_update_domains(delp_tl, domain, complete=.true.)
        call mpp_update_domains(delp, domain, complete=.true.)
! end hydro case
        IF (hydrostatic) CALL GEOPK_TLM(ptop, pe, pe_tl, peln, peln_tl, &
&                                 delp, delp_tl, pkc, pkc_tl, gz, gz_tl&
&                                 , phis, pt, pt_tl, pkz, pkz_tl, npz, &
&                                 akap, .false.)
!DELETED NH_CORE CALLS - dh
        IF (.NOT.shallow_water) THEN
          IF (breed_vortex_inline .OR. (it .EQ. n_split .AND. subt .EQ. &
&             subcycle .AND. hydrostatic)) THEN
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
                  omga_tl(i, j, k) = rdt*(pe_tl(i, k+1, j)-pem_tl(i, k+1&
&                   , j))
                  omga(i, j, k) = (pe(i, k+1, j)-pem(i, k+1, j))*rdt
                END DO
              END DO
            END DO
!------------------------------
! Compute the "advective term"
!------------------------------
            CALL ADV_PE_TLM(ua, ua_tl, va, va_tl, pem, pem_tl, omga, &
&                     omga_tl, npx, npy, npz, ng)
          END IF
        END IF
        IF (.NOT.hydrostatic .AND. m_grad_p .EQ. 0) THEN
          CALL TWO_GRAD_P_TLM(u, u_tl, v, v_tl, pkc, pkc_tl, gz, gz_tl, &
&                       delp, delp_tl, pk3, pk3_tl, divg2, divg2_tl, dt&
&                       , ng, npx, npy, npz, ptk)

        ELSE IF (beta .GT. 1.e-4) THEN
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
          call mpp_update_domains(u_tl, v_tl, domain, gridtype=DGRID_NE)
          call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
        END IF
        init_wind_m = .false.
!first_step = .false.
        elapsed_dt = elapsed_dt + dt
      END DO
    END DO !n_split
! sub-cycle
!-----------------------------------------------------
! time split loop
!-----------------------------------------------------
    IF (end_step) THEN
      DEALLOCATE(gz_tl)
      DEALLOCATE(gz)
      DEALLOCATE(ptc_tl)
      DEALLOCATE(ptc)
      DEALLOCATE(ptc_tj)
      DEALLOCATE(delzc)
      DEALLOCATE(crx_tl)
      DEALLOCATE(crx)
      DEALLOCATE(xfx_tl)
      DEALLOCATE(xfx)
      DEALLOCATE(cry_tl)
      DEALLOCATE(cry)
      DEALLOCATE(yfx_tl)
      DEALLOCATE(yfx)
      DEALLOCATE(divg_d_tl)
      DEALLOCATE(divg_d)
      DEALLOCATE(pkc_tl)
      DEALLOCATE(pkc)
      DEALLOCATE(delpc_tl)
      DEALLOCATE(delpc)
      DEALLOCATE(pkd_tl)
      DEALLOCATE(pkd)
      DEALLOCATE(pkd_tj)
      IF (.NOT.no_cgrid) THEN
        DEALLOCATE(ut_tl)
        DEALLOCATE(ut)
        DEALLOCATE(vt_tl)
        DEALLOCATE(vt)
      END IF
      IF (.NOT.hydrostatic) THEN
        DEALLOCATE(zh_tl)
        DEALLOCATE(zh)
        IF (m_grad_p .EQ. 0) THEN
          DEALLOCATE(pk3_tl)
          DEALLOCATE(pk3)
        END IF
      END IF
      IF (beta .GT. 1.e-4 .OR. no_cgrid) THEN
        DEALLOCATE(du_tl)
        DEALLOCATE(du)
        DEALLOCATE(dv_tl)
        DEALLOCATE(dv)
      END IF
    END IF
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
        pk_tl(i, j, 1) = 0.0
        pk(i, j, 1) = top_value
      END DO
    END DO
    wk_tl = 0.0
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
    wk2_tl = 0.0
    DO j=js,jep1
      DO i=is,ie
        wk2_tl(i, j) = divg2_tl(i, j) - divg2_tl(i+1, j)
        wk2(i, j) = divg2(i, j) - divg2(i+1, j)
      END DO
    END DO
    wk1_tl = 0.0
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
        pk_tl(i, j, 1) = 0.0
        pk(i, j, 1) = top_value
      END DO
    END DO
    wk_tl = 0.0
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
            delp_tl(i, j, k) = 0.0
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
          delp_tl(i, j, km) = 0.0
          delp(i, j, km) = dpmin
          delp_tl(i, j, km-1) = delp_tl(i, j, km-1) - dp_tl
          delp(i, j, km-1) = delp(i, j, km-1) - dp
          ip = ip + 1
        END IF
      END DO
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
      logp_tl = 0.0
      p1d_tl = 0.0
    ELSE
      ifirst = is - 1
      ilast = ie + 1
      jfirst = js - 1
      jlast = je + 1
      logp_tl = 0.0
      p1d_tl = 0.0
    END IF
    DO j=jfirst,jlast
      DO i=ifirst,ilast
        p1d_tl(i) = 0.0
        p1d(i) = ptop
        pk_tl(i, j, 1) = 0.0
        pk(i, j, 1) = ptk
        gh_tl(i, j, km+1) = 0.0
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

  SUBROUTINE ADV_PE_TLM(ua, ua_tl, va, va_tl, pem, pem_tl, om, om_tl, &
&   npx, npy, npz, ng)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, npz, ng
! Contra-variant wind components:
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(IN) :: ua_tl, va_tl
! Pressure at edges:
    REAL, INTENT(IN) :: pem(is-1:ie+1, npz+1, js-1:je+1)
    REAL, INTENT(IN) :: pem_tl(is-1:ie+1, npz+1, js-1:je+1)
    REAL, INTENT(INOUT) :: om(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: om_tl(isd:ied, jsd:jed, npz)
! Local:
    REAL, DIMENSION(is:ie, js:je) :: up, vp
    REAL, DIMENSION(is:ie, js:je) :: up_tl, vp_tl
    REAL :: v3(3, is:ie, js:je)
    REAL :: v3_tl(3, is:ie, js:je)
    REAL :: pin(isd:ied, jsd:jed)
    REAL :: pin_tl(isd:ied, jsd:jed)
    REAL :: pb(isd:ied, jsd:jed)
    REAL :: pb_tl(isd:ied, jsd:jed)
    REAL :: grad(3, is:ie, js:je)
    REAL :: grad_tl(3, is:ie, js:je)
    REAL :: pdx(3, is:ie, js:je+1)
    REAL :: pdx_tl(3, is:ie, js:je+1)
    REAL :: pdy(3, is:ie+1, js:je)
    REAL :: pdy_tl(3, is:ie+1, js:je)
    INTEGER :: i, j, k, n
    v3_tl = 0.0
    grad_tl = 0.0
    up_tl = 0.0
    pdx_tl = 0.0
    pdy_tl = 0.0
    pb_tl = 0.0
    vp_tl = 0.0
    pin_tl = 0.0
!$omp parallel do default(shared) private(i, j, k, n, pdx, pdy, pin, pb, up, vp, grad, v3)
    DO k=1,npz
      IF (k .EQ. npz) THEN
        DO j=js,je
          DO i=is,ie
            up_tl(i, j) = ua_tl(i, j, npz)
            up(i, j) = ua(i, j, npz)
            vp_tl(i, j) = va_tl(i, j, npz)
            vp(i, j) = va(i, j, npz)
          END DO
        END DO
      ELSE
        DO j=js,je
          DO i=is,ie
            up_tl(i, j) = 0.5*(ua_tl(i, j, k)+ua_tl(i, j, k+1))
            up(i, j) = 0.5*(ua(i, j, k)+ua(i, j, k+1))
            vp_tl(i, j) = 0.5*(va_tl(i, j, k)+va_tl(i, j, k+1))
            vp(i, j) = 0.5*(va(i, j, k)+va(i, j, k+1))
          END DO
        END DO
      END IF
! Compute Vect wind:
      DO j=js,je
        DO i=is,ie
          DO n=1,3
            v3_tl(n, i, j) = ec1(n, i, j)*up_tl(i, j) + ec2(n, i, j)*&
&             vp_tl(i, j)
            v3(n, i, j) = up(i, j)*ec1(n, i, j) + vp(i, j)*ec2(n, i, j)
          END DO
        END DO
      END DO
      DO j=js-1,je+1
        DO i=is-1,ie+1
          pin_tl(i, j) = pem_tl(i, k+1, j)
          pin(i, j) = pem(i, k+1, j)
        END DO
      END DO
! Compute pe at 4 cell corners:
      CALL A2B_ORD2_TLM(pin, pin_tl, pb, pb_tl, npx, npy, is, ie, js, je&
&                 , ng)
      DO j=js,je+1
        DO i=is,ie
          DO n=1,3
            pdx_tl(n, i, j) = dx(i, j)*en1(n, i, j)*(pb_tl(i, j)+pb_tl(i&
&             +1, j))
            pdx(n, i, j) = (pb(i, j)+pb(i+1, j))*dx(i, j)*en1(n, i, j)
          END DO
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          DO n=1,3
            pdy_tl(n, i, j) = dy(i, j)*en2(n, i, j)*(pb_tl(i, j)+pb_tl(i&
&             , j+1))
            pdy(n, i, j) = (pb(i, j)+pb(i, j+1))*dy(i, j)*en2(n, i, j)
          END DO
        END DO
      END DO
! Compute grad (pe) by Green's theorem
      DO j=js,je
        DO i=is,ie
          DO n=1,3
            grad_tl(n, i, j) = pdx_tl(n, i, j+1) - pdx_tl(n, i, j) - &
&             pdy_tl(n, i, j) + pdy_tl(n, i+1, j)
            grad(n, i, j) = pdx(n, i, j+1) - pdx(n, i, j) - pdy(n, i, j)&
&             + pdy(n, i+1, j)
          END DO
        END DO
      END DO
! Compute inner product: V3 * grad (pe)
      DO j=js,je
        DO i=is,ie
          om_tl(i, j, k) = om_tl(i, j, k) + 0.5*rarea(i, j)*(v3_tl(1, i&
&           , j)*grad(1, i, j)+v3(1, i, j)*grad_tl(1, i, j)+v3_tl(2, i, &
&           j)*grad(2, i, j)+v3(2, i, j)*grad_tl(2, i, j)+v3_tl(3, i, j)&
&           *grad(3, i, j)+v3(3, i, j)*grad_tl(3, i, j))
          om(i, j, k) = om(i, j, k) + 0.5*rarea(i, j)*(v3(1, i, j)*grad(&
&           1, i, j)+v3(2, i, j)*grad(2, i, j)+v3(3, i, j)*grad(3, i, j)&
&           )
        END DO
      END DO
    END DO
  END SUBROUTINE ADV_PE_TLM

end module dyn_core_tlm_mod
