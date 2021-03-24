module fv_dynamics_tlm_mod

   use fv_arrays_mod,         only: p_precision
   use dyn_core_tlm_mod,      only: dyn_core_tlm
   use fv_mapz_tlm_mod,       only: compute_total_energy_tlm, Lagrangian_to_Eulerian_tlm
   use fv_tracer2d_tlm_mod,   only: tracer_2d_tlm, tracer_2d_1L_tlm
   use fv_control_mod,        only: hord_mt, hord_vt, hord_tm, hord_tr, &
                                    hord_mt_pert, hord_vt_pert, hord_tm_pert, hord_tr_pert, &
                                    kord_mt, kord_tm, kord_tr, moist_phys, &
                                    kord_mt_pert, kord_tm_pert, kord_tr_pert, &
                                    range_warn, inline_q, z_tracer, tau, rf_center, nf_omega,   &
                                    te_method, remap_t,  k_top, p_ref, nwat, k_split, &
                                    check_surface_pressure, shallow_water

   use fv_grid_utils_mod,     only: g_sum
   use fv_grid_utils_mod,     only: sina_u, sina_v, sw_corner, se_corner, &
                                    ne_corner, nw_corner, da_min, ptop
   use fv_grid_utils_tlm_mod, only: c2l_ord2_tlm
   use fv_grid_utils_mod,     only: cubed_to_latlon
   use fv_grid_tools_mod,     only: dx, dy, rdxa, rdya, rdxc, rdyc, area, rarea
   use fv_mp_mod,             only: is,js,ie,je, isd,jsd,ied,jed, gid, domain
   use fv_timing_mod,         only: timing_on, timing_off

   use diag_manager_mod,      only: send_data
   use fv_diagnostics_mod,    only: id_divg, id_te, fv_time, prt_maxmin, range_check
   use mpp_domains_mod,       only: DGRID_NE, mpp_update_domains
   use field_manager_mod,     only: MODEL_ATMOS
   use tracer_manager_mod,    only: get_tracer_index
   use fv_sg_mod,             only: neg_adj3
   use tp_core_mod,           only: copy_corners
   use fv_grid_tools_mod,     only: agrid


implicit none
   logical :: RF_initialized = .false.
   logical :: bad_range
   real, allocatable ::  rf(:), rf_tl(:)
   integer:: kmax=1
private
public :: fv_dynamics_tlm

contains
 
  SUBROUTINE FV_DYNAMICS_TLM(npx, npy, npz, nq, ng, bdt, consv_te, fill&
&   , reproduce_sum, kappa, cp_air, zvir, ks, ncnst, n_split, q_split, u&
&   , u_tl, v, v_tl, um, um_tl, vm, vm_tl, w, w_tl, delz, delz_tl, &
&   hydrostatic, pt, pt_tl, delp, delp_tl, q, q_tl, ps, pe, pe_tl, pk, &
&   pk_tl, peln, peln_tl, pkz, pkz_tl, phis, omga, omga_tl, ua, ua_tl, &
&   va, va_tl, uc, uc_tl, vc, vc_tl, ak, bk, mfx, mfx_tl, mfy, mfy_tl, &
&   cx, cx_tl, cy, cy_tl, ze0, ze0_tl, hybrid_z, rdgas, grav, pi, radius&
&   , hlv, proc_id, time_total, elapsed_time, advection_test_case3)
    IMPLICIT NONE
! Large time-step
    REAL, INTENT(IN) :: bdt
    INTEGER, INTENT(IN) :: proc_id
    REAL, INTENT(IN) :: consv_te
    REAL, INTENT(IN) :: kappa, cp_air
    REAL, INTENT(IN) :: zvir
    REAL, INTENT(IN), OPTIONAL :: time_total
    REAL, INTENT(IN), OPTIONAL :: elapsed_time
    LOGICAL, INTENT(IN), OPTIONAL :: advection_test_case3
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
! transported tracers
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: ng
    INTEGER, INTENT(IN) :: ks
    INTEGER, INTENT(IN) :: ncnst
! small-step horizontal dynamics
    INTEGER, INTENT(IN) :: n_split
! tracer
    INTEGER, INTENT(IN) :: q_split
    LOGICAL, INTENT(IN) :: fill
    LOGICAL, INTENT(IN) :: reproduce_sum
    LOGICAL, INTENT(IN) :: hydrostatic
! Using hybrid_z for remapping
    LOGICAL, INTENT(IN) :: hybrid_z
! D grid zonal wind (m/s)
    REAL, DIMENSION(isd:ied, jsd:jed + 1, npz), INTENT(INOUT) :: u, um
    REAL, DIMENSION(isd:ied, jsd:jed+1, npz), INTENT(INOUT) :: u_tl, &
&   um_tl
! D grid meridional wind (m/s)
    REAL, DIMENSION(isd:ied + 1, jsd:jed, npz), INTENT(INOUT) :: v, vm
    REAL, DIMENSION(isd:ied+1, jsd:jed, npz), INTENT(INOUT) :: v_tl, &
&   vm_tl
    REAL, INTENT(IN) :: rdgas, grav, pi, radius, hlv
!  W (m/s)
    REAL, INTENT(INOUT) :: w(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: w_tl(isd:ied, jsd:jed, npz)
! temperature (K)
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: pt_tl(isd:ied, jsd:jed, npz)
! pressure thickness (pascal)
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: delp_tl(isd:ied, jsd:jed, npz)
! specific humidity and constituents
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, npz, ncnst)
    REAL, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed, npz, ncnst)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(is:ie, js:je, npz)
    REAL, INTENT(INOUT) :: delz_tl(is:ie, js:je, npz)
! height at edges (m); non-hydrostatic
    REAL, INTENT(INOUT) :: ze0(is:ie, js:je, npz+1)
    REAL, INTENT(INOUT) :: ze0_tl(is:ie, js:je, npz+1)
!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
! Surface pressure (pascal)
    REAL(p_precision), INTENT(INOUT) :: ps(isd:ied, jsd:jed)
! edge pressure (pascal)
    REAL(p_precision), INTENT(INOUT) :: pe(is-1:ie+1, npz+1, js-1:je+1)
    REAL(p_precision), INTENT(INOUT) :: pe_tl(is-1:ie+1, npz+1, js-1:je+&
&   1)
! pe**cappa
    REAL(p_precision), INTENT(INOUT) :: pk(is:ie, js:je, npz+1)
    REAL(p_precision), INTENT(INOUT) :: pk_tl(is:ie, js:je, npz+1)
! ln(pe)
    REAL(p_precision), INTENT(INOUT) :: peln(is:ie, npz+1, js:je)
    REAL(p_precision), INTENT(INOUT) :: peln_tl(is:ie, npz+1, js:je)
! finite-volume mean pk
    REAL(p_precision), INTENT(INOUT) :: pkz(is:ie, js:je, npz)
    REAL(p_precision), INTENT(INOUT) :: pkz_tl(is:ie, js:je, npz)
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
! Surface geopotential (g*Z_surf)
    REAL, INTENT(INOUT) :: phis(isd:ied, jsd:jed)
! Vertical pressure velocity (pa/s)
    REAL, INTENT(INOUT) :: omga(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: omga_tl(isd:ied, jsd:jed, npz)
! (uc,vc) mostly used as the C grid winds
    REAL, INTENT(INOUT) :: uc(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: uc_tl(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: vc(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: vc_tl(isd:ied, jsd:jed+1, npz)
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, npz), INTENT(INOUT) :: ua_tl, &
&   va_tl
    REAL, DIMENSION(npz + 1), INTENT(IN) :: ak, bk
! Accumulated Mass flux arrays: the "Flux Capacitor"
    REAL, INTENT(INOUT) :: mfx(is:ie+1, js:je, npz)
    REAL, INTENT(INOUT) :: mfx_tl(is:ie+1, js:je, npz)
    REAL, INTENT(INOUT) :: mfy(is:ie, js:je+1, npz)
    REAL, INTENT(INOUT) :: mfy_tl(is:ie, js:je+1, npz)
! Accumulated Courant number arrays
    REAL, INTENT(INOUT) :: cx(is:ie+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: cx_tl(is:ie+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: cy(isd:ied, js:je+1, npz)
    REAL, INTENT(INOUT) :: cy_tl(isd:ied, js:je+1, npz)
! Local Arrays
    REAL :: q2(isd:ied, jsd:jed, nq)
    REAL :: q2_tl(isd:ied, jsd:jed, nq)
    REAL :: te_2d(is:ie, js:je)
    REAL :: te_2d_tl(is:ie, js:je)
    REAL :: teq(is:ie, js:je)
    REAL :: pfull(npz)
    REAL :: gz(is:ie)
    REAL, ALLOCATABLE :: dp1(:, :, :)
    REAL, ALLOCATABLE :: dp1_tl(:, :, :)
    REAL, ALLOCATABLE :: pem(:, :, :)
    REAL, ALLOCATABLE :: pem_tl(:, :, :)
    REAL :: akap, rg, ph1, ph2, mdt
    INTEGER :: i, j, k, iq, n_map
    INTEGER :: sphum
    LOGICAL :: used, last_step
    REAL :: ps_sum
    INTRINSIC LOG
    INTRINSIC REAL
    INTRINSIC ALL

    ALLOCATE(dp1_tl(is:ie, js:je, 1:npz))
    ALLOCATE(dp1(is:ie, js:je, 1:npz))
    ALLOCATE(pem_tl(is-1:ie+1, 1:npz+1, js-1:je+1))
    ALLOCATE(pem(is-1:ie+1, 1:npz+1, js-1:je+1))
! shallow_water
    IF (shallow_water) THEN
      akap = 1.
      pfull(1) = 0.5*p_ref
      ua_tl = 0.0
      va_tl = 0.0
      te_2d_tl = 0.0
    ELSE
      sphum = 1
      akap = kappa
      rg = kappa*cp_air
      DO k=1,npz
        ph1 = ak(k) + bk(k)*p_ref
        ph2 = ak(k+1) + bk(k+1)*p_ref
        pfull(k) = (ph2-ph1)/LOG(ph2/ph1)
      END DO
!Compute Total Energy
!--------------------
      IF (consv_te .GT. 0.) THEN
        CALL COMPUTE_TOTAL_ENERGY_TLM(is, ie, js, je, isd, ied, jsd, jed&
&                               , npz, u, u_tl, v, v_tl, w, delz, &
&                               delz_tl, pt, pt_tl, delp, delp_tl, q, &
&                               q_tl, pe, pe_tl, peln, peln_tl, phis, &
&                               grav, zvir, cp_air, rg, hlv, te_2d, &
&                               te_2d_tl, ua, va, teq, moist_phys, sphum&
&                               , hydrostatic, id_te)
      ELSE
        te_2d_tl = 0.0
      END IF
      IF (tau .GT. 0.) THEN
        CALL RAYLEIGH_FRICTION_TLM(bdt, npx, npy, npz, ks, pfull, tau, &
&                            rf_center, u, u_tl, v, v_tl, w, w_tl, pt, &
&                            pt_tl, ua, ua_tl, va, va_tl, delz, delz_tl&
&                            , cp_air, rg, hydrostatic, .true.)
      ELSE
        ua_tl = 0.0
        va_tl = 0.0
      END IF
!Convert pt to virtual potential temperature * CP
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            pt_tl(i, j, k) = (cp_air*pt_tl(i, j, k)*pkz(i, j, k)-cp_air*&
&             pt(i, j, k)*pkz_tl(i, j, k))*(1.+zvir*q(i, j, k, sphum))/&
&             pkz(i, j, k)**2 + cp_air*pt(i, j, k)*zvir*q_tl(i, j, k, &
&             sphum)/pkz(i, j, k)
            pt(i, j, k) = cp_air*pt(i, j, k)/pkz(i, j, k)*(1.+zvir*q(i, &
&             j, k, sphum))
          END DO
        END DO
      END DO
    END IF
    mdt = bdt/REAL(k_split)
    uc_tl = 0.0
    um_tl = 0.0
    omga_tl = 0.0
    vc_tl = 0.0
    vm_tl = 0.0
    pk_tl = 0.0_8
    q2_tl = 0.0
! first level of time-split
    DO n_map=1,k_split
      IF (n_map .EQ. k_split) THEN
        last_step = .true.
      ELSE
        last_step = .false.
      END IF
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            dp1_tl(i, j, k) = delp_tl(i, j, k)
            dp1(i, j, k) = delp(i, j, k)
          END DO
        END DO
      END DO
      CALL DYN_CORE_TLM(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir&
&                 , cp_air, akap, rdgas, grav, hydrostatic, u, u_tl, v, &
&                 v_tl, um, um_tl, vm, vm_tl, w, w_tl, delz, pt, pt_tl, &
&                 q, q_tl, delp, delp_tl, pe, pe_tl, pk, pk_tl, phis, &
&                 omga, omga_tl, ptop, pfull, ua, ua_tl, va, va_tl, uc, &
&                 uc_tl, vc, vc_tl, mfx, mfx_tl, mfy, mfy_tl, cx, cx_tl&
&                 , cy, cy_tl, pem, pem_tl, pkz, pkz_tl, peln, peln_tl, &
&                 ak, bk, n_map .EQ. 1, last_step, proc_id, time_total)
! shallow_water
      IF (shallow_water) THEN
        DO j=js,je
          DO i=is,ie
            ps(i, j) = delp(i, j, 1)/grav
          END DO
        END DO
      ELSE
        IF (.NOT.inline_q) THEN
!Diagnose divergence:
          IF (nq .NE. 0) THEN
!Perform large-time-step scalar transport using the accumulated CFL and mass fluxes
            IF (.NOT.ALL(q(is:ie, js:je, 1:npz, 1:nq) .EQ. 0.0)) THEN
              IF (z_tracer) THEN
                DO k=1,npz
                  DO iq=1,nq
                    DO j=js,je
                      DO i=is,ie
                        q2_tl(i, j, iq) = q_tl(i, j, k, iq)
                        q2(i, j, iq) = q(i, j, k, iq)
                      END DO
                    END DO
                  END DO
                  CALL TRACER_2D_1L_TLM(q2, q2_tl, dp1(is, js, k), &
&                                 dp1_tl(is, js, k), mfx(is, js, k), &
&                                 mfx_tl(is, js, k), mfy(is, js, k), &
&                                 mfy_tl(is, js, k), cx(is, jsd, k), &
&                                 cx_tl(is, jsd, k), cy(isd, js, k), &
&                                 cy_tl(isd, js, k), npx, npy, npz, nq, &
&                                 hord_tr_pert, q_split, k, q, q_tl, mdt, &
&                                 id_divg)
                END DO
              ELSE
                CALL TRACER_2D_TLM(q, q_tl, dp1, dp1_tl, mfx, mfx_tl, &
&                            mfy, mfy_tl, cx, cx_tl, cy, cy_tl, npx, npy&
&                            , npz, nq, hord_tr_pert, q_split, mdt, id_divg)
              END IF
            END IF
          END IF
        END IF
        IF (npz .GT. 4) CALL LAGRANGIAN_TO_EULERIAN_TLM(last_step, &
&                                                 consv_te, ps, pe, &
&                                                 pe_tl, delp, delp_tl, &
&                                                 pkz, pkz_tl, pk, pk_tl&
&                                                 , bdt, npz, is, ie, js&
&                                                 , je, isd, ied, jsd, &
&                                                 jed, nq, sphum, u, &
&                                                 u_tl, v, v_tl, w, w_tl&
&                                                 , delz, delz_tl, pt, &
&                                                 pt_tl, q, q_tl, phis, &
&                                                 zvir, cp_air, akap, pi&
&                                                 , radius, grav, &
&                                                 kord_mt, kord_tr, &
&                                                 kord_tm, & 
&                                                 kord_mt_pert, kord_tr_pert, &
&                                                 kord_tm_pert, &
&                                                 peln, peln_tl&
&                                                 , te_2d, te_2d_tl, ng&
&                                                 , ua, ua_tl, va, omga&
&                                                 , omga_tl, dp1, dp1_tl&
&                                                 , pem, fill, &
&                                                 reproduce_sum, ak, bk&
&                                                 , ks, ze0, ze0_tl, &
&                                                 te_method, remap_t, &
&                                                 hydrostatic, hybrid_z&
&                                                 , last_step, k_top, &
&                                                 ncnst)
!Perform vertical remapping
      END IF
    END DO
! n_map loop
    IF (.NOT.shallow_water) THEN
!Convert back to temperature
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            pt_tl(i, j, k) = ((pt_tl(i, j, k)*pkz(i, j, k)+pt(i, j, k)*&
&             pkz_tl(i, j, k))*cp_air*(1.+zvir*q(i, j, k, sphum))-pt(i, &
&             j, k)*pkz(i, j, k)*cp_air*zvir*q_tl(i, j, k, sphum))/(&
&             cp_air*(1.+zvir*q(i, j, k, sphum)))**2
            pt(i, j, k) = pt(i, j, k)*pkz(i, j, k)/(cp_air*(1.+zvir*q(i&
&             , j, k, sphum)))
          END DO
        END DO
      END DO
    END IF
    DEALLOCATE(dp1_tl)
    DEALLOCATE(dp1)
    DEALLOCATE(pem_tl)
    DEALLOCATE(pem)
  END SUBROUTINE FV_DYNAMICS_TLM


  SUBROUTINE RAYLEIGH_FRICTION_TLM(dt, npx, npy, npz, ks, pm, tau, p_c, &
&   u, u_tl, v, v_tl, w, w_tl, pt, pt_tl, ua, ua_tl, va, va_tl, delz, &
&   delz_tl, cp, rg, hydrostatic, conserve)
    IMPLICIT NONE
    REAL, INTENT(IN) :: dt
! time scale (days)
    REAL, INTENT(IN) :: tau
    REAL, INTENT(IN) :: p_c
    REAL, INTENT(IN) :: cp, rg
    INTEGER, INTENT(IN) :: npx, npy, npz, ks
    REAL, DIMENSION(npz), INTENT(IN) :: pm
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: conserve
! D grid zonal wind (m/s)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, npz)
    REAL, INTENT(INOUT) :: u_tl(isd:ied, jsd:jed+1, npz)
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: v_tl(isd:ied+1, jsd:jed, npz)
! cell center vertical wind (m/s)
    REAL, INTENT(INOUT) :: w(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: w_tl(isd:ied, jsd:jed, npz)
! temp
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: pt_tl(isd:ied, jsd:jed, npz)
! 
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: ua_tl(isd:ied, jsd:jed, npz)
! 
    REAL, INTENT(INOUT) :: va(isd:ied, jsd:jed, npz)
    REAL, INTENT(INOUT) :: va_tl(isd:ied, jsd:jed, npz)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(is:ie, js:je, npz)
    REAL, INTENT(INOUT) :: delz_tl(is:ie, js:je, npz)
! local:
    REAL, ALLOCATABLE :: u2f(:, :, :)
    REAL, ALLOCATABLE :: u2f_tl(:, :, :)
    REAL, PARAMETER :: sday=86400.
! scaling velocity  **2
    REAL, PARAMETER :: u000=4900.
    REAL :: c1, pc, fac
    INTEGER :: i, j, k
    INTRINSIC LOG10
    INTRINSIC TANH
    INTRINSIC SQRT
    REAL :: arg1
    REAL :: arg1_tl
    REAL :: result1
    REAL :: result1_tl
    IF (.NOT.rf_initialized) THEN
      ALLOCATE(rf_tl(npz))
      ALLOCATE(rf(npz))
      !ALLOCATE(rw(npz))
      IF (p_c .LE. 0.) THEN
        pc = pm(1)
      ELSE
        pc = p_c
      END IF
      IF (gid .EQ. 0) WRITE(6, *) &
&                     'Rayleigh friction E-folding time [days]:'
      c1 = 1./(tau*sday)
      kmax = 1
      DO k=1,npz
        IF (pm(k) .LT. 40.e2) THEN
          arg1 = LOG10(pc/pm(k))
          rf(k) = c1*(1.+TANH(arg1))
          kmax = k
          IF (gid .EQ. 0) WRITE(6, *) k, 0.01*pm(k), 1./(rf(k)*sday)
        ELSE
          GOTO 100
        END IF
      END DO
 100  IF (gid .EQ. 0) WRITE(6, *) 'Rayleigh Friction kmax=', kmax
      rf_initialized = .true.
    END IF
    ALLOCATE(u2f_tl(isd:ied, jsd:jed, kmax))
    ALLOCATE(u2f(isd:ied, jsd:jed, kmax))
    CALL C2L_ORD2_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, dx, dy, &
&               rdxa, rdya, npz)
    u2f_tl = 0.0
    u2f = 0.
    DO k=1,kmax
      IF (hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            u2f_tl(i, j, k) = 2*ua(i, j, k)*ua_tl(i, j, k) + 2*va(i, j, &
&             k)*va_tl(i, j, k)
            u2f(i, j, k) = ua(i, j, k)**2 + va(i, j, k)**2
          END DO
        END DO
      ELSE
        DO j=js,je
          DO i=is,ie
            u2f_tl(i, j, k) = 2*ua(i, j, k)*ua_tl(i, j, k) + 2*va(i, j, &
&             k)*va_tl(i, j, k) + 2*w(i, j, k)*w_tl(i, j, k)
            u2f(i, j, k) = ua(i, j, k)**2 + va(i, j, k)**2 + w(i, j, k)&
&             **2
          END DO
        END DO
      END IF
    END DO

    call mpp_update_domains(u2f   , domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
    call mpp_update_domains(u2f_tl, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)

    DO k=1,kmax
      IF (conserve) THEN
        IF (hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              arg1_tl = u2f_tl(i, j, k)/u000
              arg1 = u2f(i, j, k)/u000
              IF (arg1 .EQ. 0.0) THEN
                result1_tl = 0.0
              ELSE
                result1_tl = arg1_tl/(2.0*SQRT(arg1))
              END IF
              result1 = SQRT(arg1)
              pt_tl(i, j, k) = pt_tl(i, j, k) + 0.5*u2f_tl(i, j, k)*(1.-&
&               1./(1.+dt*rf(k)*result1)**2)/(cp-rg*ptop/pm(k)) + 0.5*&
&               u2f(i, j, k)*2*dt*rf(k)*result1_tl/((cp-rg*ptop/pm(k))*(&
&               1.+dt*rf(k)*result1)**3)
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)/(cp-rg*ptop/&
&               pm(k))*(1.-1./(1.+dt*rf(k)*result1)**2)
            END DO
          END DO
        ELSE
          DO j=js,je
            DO i=is,ie
              delz_tl(i, j, k) = (delz_tl(i, j, k)*pt(i, j, k)-delz(i, j&
&               , k)*pt_tl(i, j, k))/pt(i, j, k)**2
              delz(i, j, k) = delz(i, j, k)/pt(i, j, k)
              arg1_tl = u2f_tl(i, j, k)/u000
              arg1 = u2f(i, j, k)/u000
              IF (arg1 .EQ. 0.0) THEN
                result1_tl = 0.0
              ELSE
                result1_tl = arg1_tl/(2.0*SQRT(arg1))
              END IF
              result1 = SQRT(arg1)
              pt_tl(i, j, k) = pt_tl(i, j, k) + 0.5*u2f_tl(i, j, k)*(1.-&
&               1./(1.+dt*rf(k)*result1)**2)/(cp-rg*ptop/pm(k)) + 0.5*&
&               u2f(i, j, k)*2*dt*rf(k)*result1_tl/((cp-rg*ptop/pm(k))*(&
&               1.+dt*rf(k)*result1)**3)
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)/(cp-rg*ptop/&
&               pm(k))*(1.-1./(1.+dt*rf(k)*result1)**2)
              delz_tl(i, j, k) = delz_tl(i, j, k)*pt(i, j, k) + delz(i, &
&               j, k)*pt_tl(i, j, k)
              delz(i, j, k) = delz(i, j, k)*pt(i, j, k)
            END DO
          END DO
        END IF
      END IF
      DO j=js-1,je+1
        DO i=is-1,ie+1
          arg1_tl = u2f_tl(i, j, k)/u000
          arg1 = u2f(i, j, k)/u000
          IF (arg1 .EQ. 0.0) THEN
            result1_tl = 0.0
          ELSE
            result1_tl = arg1_tl/(2.0*SQRT(arg1))
          END IF
          result1 = SQRT(arg1)
          u2f_tl(i, j, k) = dt*rf(k)*result1_tl
          u2f(i, j, k) = dt*rf(k)*result1
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          u_tl(i, j, k) = (u_tl(i, j, k)*(1.+0.5*(u2f(i, j-1, k)+u2f(i, &
&           j, k)))-u(i, j, k)*0.5*(u2f_tl(i, j-1, k)+u2f_tl(i, j, k)))/&
&           (1.+0.5*(u2f(i, j-1, k)+u2f(i, j, k)))**2
          u(i, j, k) = u(i, j, k)/(1.+0.5*(u2f(i, j-1, k)+u2f(i, j, k)))
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v_tl(i, j, k) = (v_tl(i, j, k)*(1.+0.5*(u2f(i-1, j, k)+u2f(i, &
&           j, k)))-v(i, j, k)*0.5*(u2f_tl(i-1, j, k)+u2f_tl(i, j, k)))/&
&           (1.+0.5*(u2f(i-1, j, k)+u2f(i, j, k)))**2
          v(i, j, k) = v(i, j, k)/(1.+0.5*(u2f(i-1, j, k)+u2f(i, j, k)))
        END DO
      END DO
      IF (.NOT.hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            w_tl(i, j, k) = (w_tl(i, j, k)*(1.+u2f(i, j, k))-w(i, j, k)*&
&             u2f_tl(i, j, k))/(1.+u2f(i, j, k))**2
            w(i, j, k) = w(i, j, k)/(1.+u2f(i, j, k))
          END DO
        END DO
      END IF
    END DO
    DEALLOCATE(u2f_tl)
    DEALLOCATE(u2f)
  END SUBROUTINE RAYLEIGH_FRICTION_TLM

end module fv_dynamics_tlm_mod
