module fv_dynamics_tlm_mod

   use fv_arrays_mod,      only: p_precision
!#ifndef MAPL_MODE
!   use constants_mod,   only: grav, pi, radius, hlv, rdgas, kappa
!#endif
!   use dyn_core_mod,    only: dyn_core
!   use fv_mapz_mod,     only: compute_total_energy, Lagrangian_to_Eulerian
!   use fv_tracer2d_mod, only: tracer_2d, tracer_2d_1L
   use dyn_core_tlm_mod,    only: dyn_core_tlm
   use fv_mapz_tlm_mod,     only: compute_total_energy_tlm, Lagrangian_to_Eulerian_tlm
   use dyn_core_adm_mod,    only: dyn_core
   use fv_mapz_adm_mod,     only: compute_total_energy
   use fv_tracer2d_adm_mod, only: tracer_2d, tracer_2d_1L
   use fv_tracer2d_tlm_mod, only: tracer_2d_tlm, tracer_2d_1L_tlm
   use fv_control_mod,  only: hord_mt, hord_vt, hord_tm, hord_tr, hord_tr_pert, &
                              kord_mt, kord_tm, kord_tr, moist_phys, &
                              kord_mt_pert, kord_tm_pert, kord_tr_pert, &
                              range_warn, inline_q, z_tracer, tau, rf_center, nf_omega,   &
                              te_method, remap_t,  k_top, p_ref, nwat, fv_debug, k_split, &
                              check_surface_pressure, shallow_water

   use fv_grid_utils_mod, only: g_sum
   use fv_grid_utils_mod, only: sina_u, sina_v, sw_corner, se_corner, &
                                ne_corner, nw_corner, da_min, ptop, c2l_ord2
   use fv_grid_utils_tlm_mod, only: c2l_ord2_tlm
   use fv_grid_utils_mod, only: cubed_to_latlon
   use fv_grid_tools_mod, only: dx, dy, rdxa, rdya, rdxc, rdyc, area, rarea
   use fv_mp_mod,         only: is,js,ie,je, isd,jsd,ied,jed, gid, domain
   use fv_timing_mod,     only: timing_on, timing_off

   use diag_manager_mod,   only: send_data
   use fv_diagnostics_mod, only: id_divg, id_te, fv_time, prt_maxmin, range_check
   use mpp_domains_mod,    only: DGRID_NE, mpp_update_domains
   use field_manager_mod,  only: MODEL_ATMOS
   use tracer_manager_mod, only: get_tracer_index
   use fv_sg_mod,          only: neg_adj3
   use tp_core_mod,        only: copy_corners
   use fv_grid_tools_mod,  only: agrid



!#ifdef WAVE_MAKER
!   use time_manager_mod,  only: get_time
!#endif


implicit none
   logical :: RF_initialized = .false.
   logical :: bad_range
   real, allocatable ::  rf(:), rw(:)
   real, allocatable ::  rf_tl(:), rw_tl(:)
   integer:: kmax=1
private
public :: fv_dynamics_tlm, FVState2PertState, PertState2FVState

contains

  SUBROUTINE FV_DYNAMICS_TLM(npx, npy, npz, nq, ng, bdt, consv_te, fill&
&   , reproduce_sum, kappa, cp_air, zvir, ks, ncnst, n_split, q_split, u&
&   , u_tl, v, v_tl, um, um_tl, vm, vm_tl, w, w_tl, delz, delz_tl, &
&   hydrostatic, pt, pt_tl, delp, delp_tl, q, q_tl, ps, pe, pe_tl, pk, &
&   pk_tl, peln, peln_tl, pkz, pkz_tl, phis, omga, omga_tl, ua, ua_tl, &
&   va, va_tl, uc, uc_tl, vc, vc_tl, ak, bk, mfx, mfx_tl, mfy, mfy_tl, &
&   cx, cx_tl, cy, cy_tl, ze0, hybrid_z, rdgas, grav, pi, radius, hlv, &
&   time_total, elapsed_time, advection_test_case3)
    IMPLICIT NONE
! Large time-step
    REAL, INTENT(IN) :: bdt
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
!#ifdef MAPL_MODE        
    REAL, INTENT(IN) :: rdgas, grav, pi, radius, hlv
!#endif
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
    REAL :: dp1(is:ie, js:je, npz)
    REAL :: dp1_tl(is:ie, js:je, npz)
    REAL :: pem(is-1:ie+1, npz+1, js-1:je+1)
    REAL :: akap, rg, ph1, ph2, mdt
    INTEGER :: i, j, k, iq, n_map
! GFDL physics
    INTEGER :: sphum, liq_wat, ice_wat
    INTEGER :: rainwat, snowwat, graupel, cld_amt
    LOGICAL :: used, last_step
    REAL :: ps_sum
    INTRINSIC LOG
    INTRINSIC REAL
!     real te_den
!#ifdef WAVE_MAKER
!      integer seconds, days
!      real  r0, stime
!
!         call get_time (fv_time, seconds,  days)
!         r0 = pi/30.
!         stime = real(seconds)/86400.*2.*pi
!         do j=jsd,jed
!            do i=isd,ied
!               phis(i,j) = grav*250.*sin(agrid(i,j,1))*sin(stime) / exp( (agrid(i,j,2)/r0)**2 )
!            enddo
!         enddo
!#endif
!      allocate ( dp1(is:ie, js:je, 1:npz) )
!      allocate ( pem(is-1:ie+1, 1:npz+1, js-1:je+1) )
!call timing_on('FV_DYNAMICS')
!Compute the FV variables internally, for checkpointing purposes
    CALL PERTSTATE2FVSTATE_TLM(delp, delp_tl, pe, pe_tl, pk, pk_tl, pkz&
&                        , pkz_tl, peln, peln_tl, ptop, pt, pt_tl, isd, &
&                        ied, jsd, jed, is, ie, js, je, npz, kappa)
! shallow_water
    IF (shallow_water) THEN
      akap = 1.
      pfull(1) = 0.5*p_ref
      ua_tl = 0.0
      va_tl = 0.0
      te_2d_tl = 0.0
    ELSE
!      if ( nwat==6 ) then
!             sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
!           liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
!           ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
!           rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
!           snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
!           graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
!           cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
!      else
      sphum = 1
!      endif
      akap = kappa
      rg = kappa*cp_air
      DO k=1,npz
        ph1 = ak(k) + bk(k)*p_ref
        ph2 = ak(k+1) + bk(k+1)*p_ref
        pfull(k) = (ph2-ph1)/LOG(ph2/ph1)
      END DO
!      if ( fv_debug ) then
!         call prt_maxmin('T_dyn_b',   pt, is, ie, js, je, ng, npz, 1., gid==0)
!         call prt_maxmin('delp_b ', delp, is, ie, js, je, ng, npz, 0.01, gid==0)
!         call prt_maxmin('pkz_b  ',  REAL(pkz), is, ie, js, je,  0, npz, 1., gid==0)
!      endif
!---------------------
! Compute Total Energy
!---------------------
      IF (consv_te .GT. 0.) THEN
!#ifdef MAPL_MODE
!#endif
!        CALL COMPUTE_TOTAL_ENERGY_TLM(is, ie, js, je, isd, ied, jsd, jed&
!&                               , npz, u, u_tl, v, v_tl, w, delz, &
!&                               delz_tl, pt, pt_tl, delp, delp_tl, q, &
!&                               q_tl, pe, pe_tl, peln, peln_tl, phis, &
!&                               grav, zvir, cp_air, rg, hlv, te_2d, &
!&                               te_2d_tl, ua, teq, moist_phys, sphum, &
!&                               hydrostatic, id_te)
        CALL COMPUTE_TOTAL_ENERGY(is, ie, js, je, isd, ied, jsd, jed, &
&                           npz, u, v, w, delz, pt, delp, q, pe, peln, &
&                           phis, grav, zvir, cp_air, rg, hlv, te_2d, ua&
&                           , teq, moist_phys, sphum, hydrostatic, id_te&
&                          )
!#ifndef MAPL_MODE
!           if( id_te>0 ) then
!               used = send_data(id_te, teq, fv_time)
!!              te_den=1.E-9*g_sum(teq, is, ie, js, je, ng, area, 0)/(grav*4.*pi*radius**2)
!!              if(gid==0)  write(*,*) 'Total Energy Density (Giga J/m**2)=',te_den
!           endif
!#endif
      ELSE
        te_2d_tl = 0.0
      END IF
      IF (tau .GT. 0.) THEN
!        CALL RAYLEIGH_FRICTION_TLM(bdt, npx, npy, npz, ks, pfull, tau, &
!&                            rf_center, u, u_tl, v, v_tl, w, w_tl, pt, &
!&                            pt_tl, ua, ua_tl, va, va_tl, delz, delz_tl&
!&                            , cp_air, rg, hydrostatic, .true.)
        CALL RAYLEIGH_FRICTION(bdt, npx, npy, npz, ks, pfull, tau, &
&                        rf_center, u, v, w, pt, ua, va, delz, cp_air, &
&                        rg, hydrostatic, .true.)
      ELSE
        ua_tl = 0.0
        va_tl = 0.0
      END IF
! Convert pt to virtual potential temperature * CP
!$omp parallel do default(shared) private(i, j, k)
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
!      if ( fv_debug ) then
!         call prt_maxmin('PT_dyn_b',   pt, is, ie, js, je, ng, npz, 1., gid==0)
!      endif
    mdt = bdt/REAL(k_split)
!  do n_map=1, k_split   ! first level of time-split
    n_map = 1
!if ( n_map==k_split )  then
    last_step = .true.
    dp1_tl = 0.0
!else
!     last_step = .false.
!endif
!$omp parallel do default(shared) private(i, j, k)
    DO k=1,npz
      DO j=js,je
        DO i=is,ie
          dp1_tl(i, j, k) = delp_tl(i, j, k)
          dp1(i, j, k) = delp(i, j, k)
        END DO
      END DO
    END DO
!    CALL DYN_CORE_TLM(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir, &
!&               cp_air, akap, rdgas, grav, hydrostatic, u, u_tl, v, v_tl&
!&               , um, um_tl, vm, vm_tl, w, w_tl, delz, pt, pt_tl, q, &
!&               q_tl, delp, delp_tl, pe, pe_tl, pk, pk_tl, phis, omga, &
!&               omga_tl, ptop, pfull, ua, ua_tl, va, va_tl, uc, uc_tl, &
!&               vc, vc_tl, mfx, mfx_tl, mfy, mfy_tl, cx, cx_tl, cy, &
!&               cy_tl, pem, pkz, pkz_tl, peln, peln_tl, ak, bk)
    CALL DYN_CORE(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir, &
&           cp_air, akap, rdgas, grav, hydrostatic, u, v, um, vm, w, &
&           delz, pt, q, delp, pe, pk, phis, omga, ptop, pfull, ua, va, &
&           uc, vc, mfx, mfy, cx, cy, pem, pkz, peln, ak, bk)
! shallow_water
!, n_map==1, last_step,  time_total)
    IF (.NOT.shallow_water) THEN
!do j=js,je
!   do i=is,ie
!      ps(i,j) = delp(i,j,1) / grav
!   enddo
!enddo
      IF (.NOT.inline_q) THEN
! diagnose divergence:
        IF (nq .NE. 0) THEN
!--------------------------------------------------------
! Perform large-time-step scalar transport using the accumulated CFL and
! mass fluxes
!call timing_on('tracer_2d')
!     if (.not. all(q(is:ie,js:je,1:npz,1:nq) == 0.0)) then
          IF (z_tracer) THEN
            q2_tl = 0.0
            DO k=1,npz
              DO iq=1,nq
                DO j=js,je
! To_do list:
                  DO i=is,ie
! The data copying can be avoided if q is
                    q2_tl(i, j, iq) = q_tl(i, j, k, iq)
                    q2(i, j, iq) = q(i, j, k, iq)
                  END DO
                END DO
              END DO
! re-dimensioned as q(i,j,nq,k)
!              CALL TRACER_2D_1L_TLM(q2, q2_tl, dp1(is, js, k), dp1_tl(is&
!&                             , js, k), mfx(is, js, k), mfx_tl(is, js, k&
!&                             ), mfy(is, js, k), mfy_tl(is, js, k), cx(&
!&                             is, jsd, k), cx_tl(is, jsd, k), cy(isd, js&
!&                             , k), cy_tl(isd, js, k), npx, npy, npz, nq&
!&                             , hord_tr_pert, q_split, k, q, q_tl, mdt, &
!&                             id_divg)
            END DO
          ELSE
! if ( fv_debug ) then
!   do iq=1,nq
!    call prt_maxmin('Q ',  q(:,:,k,iq), is, ie, js, je,  ng, 1., gid==0)
!   enddo
! endif
!            CALL TRACER_2D_TLM(q, q_tl, dp1, dp1_tl, mfx, mfx_tl, mfy, &
!&                        mfy_tl, cx, cx_tl, cy, cy_tl, npx, npy, npz, nq&
!&                        , hord_tr_pert, q_split, mdt, id_divg)
! if ( fv_debug ) then
!   do iq=1,nq
!    call prt_maxmin('Q ',  q(:,:,:,iq), is, ie, js, je,  ng, npz, 1., gid==0)
!   enddo
! endif
          END IF
        END IF
      END IF
!      endif
!call timing_off('tracer_2d')
!#ifndef MAPL_MODE
!         if( last_step .and. id_divg>0 ) used = send_data(id_divg, dp1, fv_time) 
!#endif
!      IF (npz .GT. 4) CALL LAGRANGIAN_TO_EULERIAN_TLM(last_step, &
!&                                               consv_te, pe, pe_tl, &
!&                                               delp, delp_tl, pkz, &
!&                                               pkz_tl, pk, pk_tl, bdt, &
!&                                               npz, is, ie, js, je, isd&
!&                                               , ied, jsd, jed, nq, &
!&                                               sphum, u, u_tl, v, v_tl&
!&                                               , pt, pt_tl, q, q_tl, &
!&                                               phis, zvir, cp_air, akap&
!&                                               , pi, radius, grav, &
!&                                               kord_mt, kord_tr, &
!&                                               kord_tm, peln, peln_tl, &
!&                                               te_2d, te_2d_tl, ng, ua&
!&                                               , ua_tl, dp1, dp1_tl, &
!&                                               pem, fill, reproduce_sum&
!&                                               , ak, bk, ks, te_method&
!&                                               , remap_t, ncnst)
!------------------------------------------------------------------------
! Peroform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.
!------------------------------------------------------------------------
!call timing_on('Remapping')
!      if ( fv_debug ) then
!         call prt_maxmin('pkz_a1 ',  REAL(pkz), is, ie, js, je,  0, npz, 1., gid==0)
!      endif
!#ifdef MAPL_MODE
!#endif
!      if ( fv_debug ) then
!         call prt_maxmin('pkz_a2 ',  REAL(pkz), is, ie, js, je,  0, npz, 1., gid==0)
!      endif 
!call timing_off('Remapping')
!--------------------------
! Filter omega for physics:
!--------------------------
!         if( last_step .and. nf_omega>0 )  then
!             call del2_cubed(omga, 0.20*REAL(da_min), npx, npy, npz, nf_omega)
!         endif
    END IF
!  enddo         ! n_map loop
    IF (.NOT.shallow_water) THEN
! Convert back to temperature
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
!  deallocate ( dp1 )
!  deallocate ( pem )
!  if ( fv_debug ) then
!         call prt_maxmin('delp_a',  delp, is, ie, js, je, ng, npz, 0.01, gid==0)
!         call prt_maxmin('T_dyn_a',  pt, is, ie, js, je, ng, npz, 1., gid==0)
!         call prt_maxmin('pk_a',   REAL(pk), is, ie, js, je, 0, npz+1, 1., gid==0)
!         call prt_maxmin('pkz_a',  REAL(pkz), is, ie, js, je, 0, npz, 1., gid==0)
!  endif
!  if( nwat==6 ) then
!      call neg_adj3(is, ie, js, je, ng, npz,        &
!                    pt, delp, q(isd,jsd,1,sphum),   &
!                              q(isd,jsd,1,liq_wat), &
!                              q(isd,jsd,1,rainwat), &
!                              q(isd,jsd,1,ice_wat), &
!                              q(isd,jsd,1,snowwat), &
!                              q(isd,jsd,1,graupel), &
!                              q(isd,jsd,1,cld_amt)  )
!     if ( fv_debug ) then
!       call prt_maxmin('SPHUM_dyn',   q(isd:,jsd:,1:,sphum  ), is, ie, js, je, ng, npz, 1., gid==0)
!       call prt_maxmin('liq_wat_dyn', q(isd:,jsd:,1:,liq_wat), is, ie, js, je, ng, npz, 1., gid==0)
!       call prt_maxmin('ice_wat_dyn', q(isd:,jsd:,1:,ice_wat), is, ie, js, je, ng, npz, 1., gid==0)
!       call prt_maxmin('snowwat_dyn', q(isd:,jsd:,1:,snowwat), is, ie, js, je, ng, npz, 1., gid==0)
!       call prt_maxmin('graupel_dyn', q(isd:,jsd:,1:,graupel), is, ie, js, je, ng, npz, 1., gid==0)
!!      call prt_maxmin('cld_amt_dyn', q(isd:,jsd:,1:,cld_amt), is, ie, js, je, ng, npz, 1., gid==0)
!     endif
!  endif
!
!  call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, npz, 1)
!
!  if ( range_warn ) then
!       call range_check('UA_dyn', ua, is, ie, js, je, ng, npz, agrid,   &
!                         gid==0, -220., 260., bad_range)
!       call range_check('VA_dyn', ua, is, ie, js, je, ng, npz, agrid,   &
!                         gid==0, -220., 220., bad_range)
!       if (.not. shallow_water) then
!          call range_check('TA_dyn', pt, is, ie, js, je, ng, npz, agrid,   &
!                            gid==0, 150., 350., bad_range)
!       endif
!  endif
!
!  if ( check_surface_pressure ) then
!    ps_sum=g_sum(ps(is:ie,js:je), is, ie, js, je, ng, area, 1, .true.)
!    if(gid==0)  write(*,*) 'Surface Pressure Global Sum (PA)=',ps_sum
!  endif
!call timing_off('FV_DYNAMICS')
!Convert back to potential temperature
    CALL FVSTATE2PERTSTATE_TLM(pkz, pkz_tl, pt, pt_tl, isd, ied, jsd, &
&                        jed, is, ie, js, je, npz)
  END SUBROUTINE FV_DYNAMICS_TLM
!-----------------------------------------------------------------------
!     fv_dynamics :: FV dynamical core driver
!-----------------------------------------------------------------------
!#ifdef MAPL_MODE
!#endif

!  Differentiation of rayleigh_friction in forward (tangent) mode:
!   variations   of useful results: u v w ua delz va pt
!   with respect to varying inputs: u v w delz pt
! subroutine del2_cubed(q, cd, npx, npy, km, nmax)
!!---------------------------------------------------------------
!! This routine is for filtering the omega field for the physics
!!---------------------------------------------------------------
!   integer, intent(in):: npx, npy, km, nmax
!   real,    intent(in):: cd            ! cd = K * da_min;   0 < K < 0.25
!   real, intent(inout):: q(isd:ied,jsd:jed,km)
!   real, parameter:: r3  = 1./3.
!   real :: fx(isd:ied+1,jsd:jed), fy(isd:ied,jsd:jed+1)
!   real :: q2(isd:ied,jsd:jed)
!   integer i,j,k, n, nt, ntimes
!
!   ntimes = min(3, nmax)
!
!                     !call timing_on('COMM_TOTAL')
!   call mpp_update_domains_dummy(q, domain, complete=.true.)
!                     !call timing_off('COMM_TOTAL')
!
!
!   do n=1,ntimes
!      nt = ntimes - n
!
!   do k=1,km
!
!      if ( sw_corner ) then
!           q(1,1,k) = (q(1,1,k)+q(0,1,k)+q(1,0,k)) * r3
!           q(0,1,k) =  q(1,1,k)
!           q(1,0,k) =  q(1,1,k)
!      endif
!      if ( se_corner ) then
!           q(ie, 1,k) = (q(ie,1,k)+q(npx,1,k)+q(ie,0,k)) * r3
!           q(npx,1,k) =  q(ie,1,k)
!           q(ie, 0,k) =  q(ie,1,k)
!      endif
!      if ( ne_corner ) then
!           q(ie, je,k) = (q(ie,je,k)+q(npx,je,k)+q(ie,npy,k)) * r3
!           q(npx,je,k) =  q(ie,je,k)
!           q(ie,npy,k) =  q(ie,je,k)
!      endif
!      if ( nw_corner ) then
!           q(1, je,k) = (q(1,je,k)+q(0,je,k)+q(1,npy,k)) * r3
!           q(0, je,k) =  q(1,je,k)
!           q(1,npy,k) =  q(1,je,k)
!      endif
!
!      if(nt>0) call copy_corners(q(isd:,jsd:,k), npx, npy, 1)
!      do j=js-nt,je+nt
!         do i=is-nt,ie+1+nt
!            fx(i,j) = dy(i,j)*sina_u(i,j)*(q(i-1,j,k)-q(i,j,k))*rdxc(i,j)
!         enddo
!      enddo
!
!      if(nt>0) call copy_corners(q(isd:,jsd:,k), npx, npy, 2)
!      do j=js-nt,je+1+nt
!         do i=is-nt,ie+nt
!            fy(i,j) = dx(i,j)*sina_v(i,j)*(q(i,j-1,k)-q(i,j,k))*rdyc(i,j)
!         enddo
!      enddo
!
!      do j=js-nt,je+nt
!         do i=is-nt,ie+nt
!            q(i,j,k) = q(i,j,k) + cd*rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))
!         enddo
!      enddo
!   enddo
!   enddo
!
! end subroutine del2_cubed
!
! subroutine del2_cubed_old(q, cd, npx, npy, km, ntimes)
!!---------------------------------------------------------------
!! This routine is for filtering the omega field for the physics
!!---------------------------------------------------------------
!   integer, intent(in):: npx, npy, km, ntimes
!   real,    intent(in):: cd            ! cd = K * da_min;   0 < K < 0.25
!   real, intent(inout):: q(isd:ied,jsd:jed,km)
!   real, parameter:: r3  = 1./3.
!   real :: fx(is:ie+1,js:je), fy(is:ie,js:je+1)
!   integer i,j,k, n
!
!   do n=1,ntimes
!                     !call timing_on('COMM_TOTAL')
!!  call mpp_update_domains_dummy(q, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
!   call mpp_update_domains_dummy(q, domain, complete=.true.)
!                     !call timing_off('COMM_TOTAL')
!   do k=1,km
!      if ( sw_corner ) then
!           q(1,1,k) = (q(1,1,k)+q(0,1,k)+q(1,0,k)) * r3
!           q(0,1,k) =  q(1,1,k)
!           q(1,0,k) =  q(1,1,k)
!      endif
!      if ( se_corner ) then
!           q(ie, 1,k) = (q(ie,1,k)+q(npx,1,k)+q(ie,0,k)) * r3
!           q(npx,1,k) =  q(ie,1,k)
!           q(ie, 0,k) =  q(ie,1,k)
!      endif
!      if ( ne_corner ) then
!           q(ie, je,k) = (q(ie,je,k)+q(npx,je,k)+q(ie,npy,k)) * r3
!           q(npx,je,k) =  q(ie,je,k)
!           q(ie,npy,k) =  q(ie,je,k)
!      endif
!      if ( nw_corner ) then
!           q(1, je,k) = (q(1,je,k)+q(0,je,k)+q(1,npy,k)) * r3
!           q(0, je,k) =  q(1,je,k)
!           q(1,npy,k) =  q(1,je,k)
!      endif
!
!      do j=js,je
!         do i=is,ie+1
!            fx(i,j) = cd*dy(i,j)*sina_u(i,j)*(q(i-1,j,k)-q(i,j,k))*rdxc(i,j)
!         enddo
!      enddo
!
!      do j=js,je+1
!         do i=is,ie
!            fy(i,j) = cd*dx(i,j)*sina_v(i,j)*(q(i,j-1,k)-q(i,j,k))*rdyc(i,j)
!         enddo
!      enddo
!
!      do j=js,je
!         do i=is,ie
!            q(i,j,k) = q(i,j,k) + rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))
!         enddo
!      enddo
!   enddo
!   enddo
!
! end subroutine del2_cubed_old
!
!
!
!#ifdef OLD_RAYF
!
! subroutine Rayleigh_Friction(dt, npx, npy, npz, ks, pm, tau, p_c, u, v, w, pt,  &
!                              ua, va, delz, cp, rg, hydrostatic, conserve)
!    real, intent(in):: dt
!    real, intent(in):: tau              ! time scale (days)
!    real, intent(in):: p_c
!    real, intent(in):: cp, rg
!    real, intent(in),  dimension(npz):: pm
!    integer, intent(in):: npx, npy, npz, ks
!    logical, intent(in):: hydrostatic
!    logical, intent(in):: conserve
!    real, intent(inout):: u(isd:ied  ,jsd:jed+1,npz) ! D grid zonal wind (m/s)
!    real, intent(inout):: v(isd:ied+1,jsd:jed,npz) ! D grid meridional wind (m/s)
!    real, intent(inout)::  w(isd:ied,jsd:jed,npz) ! cell center vertical wind (m/s)
!    real, intent(inout):: pt(isd:ied,jsd:jed,npz) ! temp
!    real, intent(inout):: ua(isd:ied,jsd:jed,npz) ! 
!    real, intent(inout):: va(isd:ied,jsd:jed,npz) ! 
!    real, intent(inout):: delz(is:ie,js:je,npz)   ! delta-height (m); non-hydrostatic only
!    real, parameter:: sday = 86400.
!    real, parameter:: wfac = 10.     ! factor to amplify the drag on w
!    real c1, pc, fac
!    integer i, j, k
!
!    kmax = max(npz/3+1, ks)
!
!    if ( .not. RF_initialized ) then
!          allocate( rf(npz) )
!          allocate( rw(npz) )
!
!          if ( p_c <= 0. ) then
!               pc = pm(1)
!          else
!               pc = p_c
!          endif
!
!          if( gid==0 ) write(6,*) 'Rayleigh friction E-folding time [days]:'
!          c1 = 1. / (tau*sday)
!          do k=1,kmax
!             if ( pm(k) < 30.E2 ) then
!                  rf(k) = c1*(1.+tanh(log10(pc/pm(k))))
!                  if( gid==0 ) write(6,*) k, 0.01*pm(k), 1./(rf(k)*sday)
!                  rf(k) = 1./(1.+dt*rf(k))
!                  rw(k) = 1./(1.+dt*rf(k)*wfac)
!             endif
!          enddo
!          RF_initialized = .true.
!    endif
!
!     if(conserve) call c2l_ord2(u, v, ua, va, dx, dy, rdxa, rdya, npz)
!
!     do k=1,kmax
!        if ( pm(k) < 30.E2 ) then
!! Add heat so as to conserve TE
!          if ( conserve ) then
!               fac = 0.5*(1.-rf(k)**2) / (cp-rg*ptop/pm(k))
!               do j=js,je
!                  do i=is,ie
!                     pt(i,j,k) = pt(i,j,k) + fac*(ua(i,j,k)**2 + va(i,j,k)**2)
!                  enddo
!               enddo
!          endif
!             do j=js,je+1
!                do i=is,ie
!                   u(i,j,k) = u(i,j,k)*rf(k)
!                enddo
!             enddo
!             do j=js,je
!                do i=is,ie+1
!                   v(i,j,k) = v(i,j,k)*rf(k)
!                enddo
!             enddo
!          if ( .not. hydrostatic ) then
!             do j=js,je
!                do i=is,ie
!                   w(i,j,k) = w(i,j,k)*rw(k)
!                enddo
!             enddo
!          endif
!        endif
!     enddo
!
! end subroutine Rayleigh_Friction
!
!#else
  SUBROUTINE RAYLEIGH_FRICTION_TLM(dt, npx, npy, npz, ks, pm, tau, p_c, &
&   u, u_tl, v, v_tl, w, w_tl, pt, pt_tl, ua, ua_tl, va, va_tl, delz, &
&   delz_tl, cp, rg, hydrostatic, conserve)
    IMPLICIT NONE
!     deallocate ( u2f )
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
    REAL :: u2f(isd:ied, jsd:jed, kmax)
    REAL :: u2f_tl(isd:ied, jsd:jed, kmax)
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
      ALLOCATE(rw(npz))
      IF (p_c .LE. 0.) THEN
        pc = pm(1)
      ELSE
        pc = p_c
      END IF
!if( gid==0 ) write(6,*) 'Rayleigh friction E-folding time [days]:'
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
!    allocate( u2f(isd:ied,jsd:jed,kmax) )
    CALL C2L_ORD2_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, dx, dy, &
&               rdxa, rdya, npz)
    u2f = 0.
    u2f_tl = 0.0
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
!call timing_on('COMM_TOTAL')
    call mpp_update_domains(u2f   , domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
    call mpp_update_domains(u2f_tl, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
    !CALL MPP_UPDATE_DOMAINS_DUMMY3_TLM(u2f, u2f_tl, isd, ied, jsd, jed, &
!&                                kmax)
!call timing_off('COMM_TOTAL')
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
  END SUBROUTINE RAYLEIGH_FRICTION_TLM
! subroutine del2_cubed(q, cd, npx, npy, km, nmax)
!!---------------------------------------------------------------
!! This routine is for filtering the omega field for the physics
!!---------------------------------------------------------------
!   integer, intent(in):: npx, npy, km, nmax
!   real,    intent(in):: cd            ! cd = K * da_min;   0 < K < 0.25
!   real, intent(inout):: q(isd:ied,jsd:jed,km)
!   real, parameter:: r3  = 1./3.
!   real :: fx(isd:ied+1,jsd:jed), fy(isd:ied,jsd:jed+1)
!   real :: q2(isd:ied,jsd:jed)
!   integer i,j,k, n, nt, ntimes
!
!   ntimes = min(3, nmax)
!
!                     !call timing_on('COMM_TOTAL')
!   call mpp_update_domains_dummy(q, domain, complete=.true.)
!                     !call timing_off('COMM_TOTAL')
!
!
!   do n=1,ntimes
!      nt = ntimes - n
!
!   do k=1,km
!
!      if ( sw_corner ) then
!           q(1,1,k) = (q(1,1,k)+q(0,1,k)+q(1,0,k)) * r3
!           q(0,1,k) =  q(1,1,k)
!           q(1,0,k) =  q(1,1,k)
!      endif
!      if ( se_corner ) then
!           q(ie, 1,k) = (q(ie,1,k)+q(npx,1,k)+q(ie,0,k)) * r3
!           q(npx,1,k) =  q(ie,1,k)
!           q(ie, 0,k) =  q(ie,1,k)
!      endif
!      if ( ne_corner ) then
!           q(ie, je,k) = (q(ie,je,k)+q(npx,je,k)+q(ie,npy,k)) * r3
!           q(npx,je,k) =  q(ie,je,k)
!           q(ie,npy,k) =  q(ie,je,k)
!      endif
!      if ( nw_corner ) then
!           q(1, je,k) = (q(1,je,k)+q(0,je,k)+q(1,npy,k)) * r3
!           q(0, je,k) =  q(1,je,k)
!           q(1,npy,k) =  q(1,je,k)
!      endif
!
!      if(nt>0) call copy_corners(q(isd:,jsd:,k), npx, npy, 1)
!      do j=js-nt,je+nt
!         do i=is-nt,ie+1+nt
!            fx(i,j) = dy(i,j)*sina_u(i,j)*(q(i-1,j,k)-q(i,j,k))*rdxc(i,j)
!         enddo
!      enddo
!
!      if(nt>0) call copy_corners(q(isd:,jsd:,k), npx, npy, 2)
!      do j=js-nt,je+1+nt
!         do i=is-nt,ie+nt
!            fy(i,j) = dx(i,j)*sina_v(i,j)*(q(i,j-1,k)-q(i,j,k))*rdyc(i,j)
!         enddo
!      enddo
!
!      do j=js-nt,je+nt
!         do i=is-nt,ie+nt
!            q(i,j,k) = q(i,j,k) + cd*rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))
!         enddo
!      enddo
!   enddo
!   enddo
!
! end subroutine del2_cubed
!
! subroutine del2_cubed_old(q, cd, npx, npy, km, ntimes)
!!---------------------------------------------------------------
!! This routine is for filtering the omega field for the physics
!!---------------------------------------------------------------
!   integer, intent(in):: npx, npy, km, ntimes
!   real,    intent(in):: cd            ! cd = K * da_min;   0 < K < 0.25
!   real, intent(inout):: q(isd:ied,jsd:jed,km)
!   real, parameter:: r3  = 1./3.
!   real :: fx(is:ie+1,js:je), fy(is:ie,js:je+1)
!   integer i,j,k, n
!
!   do n=1,ntimes
!                     !call timing_on('COMM_TOTAL')
!!  call mpp_update_domains_dummy(q, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
!   call mpp_update_domains_dummy(q, domain, complete=.true.)
!                     !call timing_off('COMM_TOTAL')
!   do k=1,km
!      if ( sw_corner ) then
!           q(1,1,k) = (q(1,1,k)+q(0,1,k)+q(1,0,k)) * r3
!           q(0,1,k) =  q(1,1,k)
!           q(1,0,k) =  q(1,1,k)
!      endif
!      if ( se_corner ) then
!           q(ie, 1,k) = (q(ie,1,k)+q(npx,1,k)+q(ie,0,k)) * r3
!           q(npx,1,k) =  q(ie,1,k)
!           q(ie, 0,k) =  q(ie,1,k)
!      endif
!      if ( ne_corner ) then
!           q(ie, je,k) = (q(ie,je,k)+q(npx,je,k)+q(ie,npy,k)) * r3
!           q(npx,je,k) =  q(ie,je,k)
!           q(ie,npy,k) =  q(ie,je,k)
!      endif
!      if ( nw_corner ) then
!           q(1, je,k) = (q(1,je,k)+q(0,je,k)+q(1,npy,k)) * r3
!           q(0, je,k) =  q(1,je,k)
!           q(1,npy,k) =  q(1,je,k)
!      endif
!
!      do j=js,je
!         do i=is,ie+1
!            fx(i,j) = cd*dy(i,j)*sina_u(i,j)*(q(i-1,j,k)-q(i,j,k))*rdxc(i,j)
!         enddo
!      enddo
!
!      do j=js,je+1
!         do i=is,ie
!            fy(i,j) = cd*dx(i,j)*sina_v(i,j)*(q(i,j-1,k)-q(i,j,k))*rdyc(i,j)
!         enddo
!      enddo
!
!      do j=js,je
!         do i=is,ie
!            q(i,j,k) = q(i,j,k) + rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))
!         enddo
!      enddo
!   enddo
!   enddo
!
! end subroutine del2_cubed_old
!
!
!
!#ifdef OLD_RAYF
!
! subroutine Rayleigh_Friction(dt, npx, npy, npz, ks, pm, tau, p_c, u, v, w, pt,  &
!                              ua, va, delz, cp, rg, hydrostatic, conserve)
!    real, intent(in):: dt
!    real, intent(in):: tau              ! time scale (days)
!    real, intent(in):: p_c
!    real, intent(in):: cp, rg
!    real, intent(in),  dimension(npz):: pm
!    integer, intent(in):: npx, npy, npz, ks
!    logical, intent(in):: hydrostatic
!    logical, intent(in):: conserve
!    real, intent(inout):: u(isd:ied  ,jsd:jed+1,npz) ! D grid zonal wind (m/s)
!    real, intent(inout):: v(isd:ied+1,jsd:jed,npz) ! D grid meridional wind (m/s)
!    real, intent(inout)::  w(isd:ied,jsd:jed,npz) ! cell center vertical wind (m/s)
!    real, intent(inout):: pt(isd:ied,jsd:jed,npz) ! temp
!    real, intent(inout):: ua(isd:ied,jsd:jed,npz) ! 
!    real, intent(inout):: va(isd:ied,jsd:jed,npz) ! 
!    real, intent(inout):: delz(is:ie,js:je,npz)   ! delta-height (m); non-hydrostatic only
!    real, parameter:: sday = 86400.
!    real, parameter:: wfac = 10.     ! factor to amplify the drag on w
!    real c1, pc, fac
!    integer i, j, k
!
!    kmax = max(npz/3+1, ks)
!
!    if ( .not. RF_initialized ) then
!          allocate( rf(npz) )
!          allocate( rw(npz) )
!
!          if ( p_c <= 0. ) then
!               pc = pm(1)
!          else
!               pc = p_c
!          endif
!
!          if( gid==0 ) write(6,*) 'Rayleigh friction E-folding time [days]:'
!          c1 = 1. / (tau*sday)
!          do k=1,kmax
!             if ( pm(k) < 30.E2 ) then
!                  rf(k) = c1*(1.+tanh(log10(pc/pm(k))))
!                  if( gid==0 ) write(6,*) k, 0.01*pm(k), 1./(rf(k)*sday)
!                  rf(k) = 1./(1.+dt*rf(k))
!                  rw(k) = 1./(1.+dt*rf(k)*wfac)
!             endif
!          enddo
!          RF_initialized = .true.
!    endif
!
!     if(conserve) call c2l_ord2(u, v, ua, va, dx, dy, rdxa, rdya, npz)
!
!     do k=1,kmax
!        if ( pm(k) < 30.E2 ) then
!! Add heat so as to conserve TE
!          if ( conserve ) then
!               fac = 0.5*(1.-rf(k)**2) / (cp-rg*ptop/pm(k))
!               do j=js,je
!                  do i=is,ie
!                     pt(i,j,k) = pt(i,j,k) + fac*(ua(i,j,k)**2 + va(i,j,k)**2)
!                  enddo
!               enddo
!          endif
!             do j=js,je+1
!                do i=is,ie
!                   u(i,j,k) = u(i,j,k)*rf(k)
!                enddo
!             enddo
!             do j=js,je
!                do i=is,ie+1
!                   v(i,j,k) = v(i,j,k)*rf(k)
!                enddo
!             enddo
!          if ( .not. hydrostatic ) then
!             do j=js,je
!                do i=is,ie
!                   w(i,j,k) = w(i,j,k)*rw(k)
!                enddo
!             enddo
!          endif
!        endif
!     enddo
!
! end subroutine Rayleigh_Friction
!
!#else
  SUBROUTINE RAYLEIGH_FRICTION(dt, npx, npy, npz, ks, pm, tau, p_c, u, v&
&   , w, pt, ua, va, delz, cp, rg, hydrostatic, conserve)
    IMPLICIT NONE
!     deallocate ( u2f )
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
! D grid meridional wind (m/s)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, npz)
! cell center vertical wind (m/s)
    REAL, INTENT(INOUT) :: w(isd:ied, jsd:jed, npz)
! temp
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, npz)
! 
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, npz)
! 
    REAL, INTENT(INOUT) :: va(isd:ied, jsd:jed, npz)
! delta-height (m); non-hydrostatic only
    REAL, INTENT(INOUT) :: delz(is:ie, js:je, npz)
! local:
    REAL :: u2f(isd:ied, jsd:jed, kmax)
    REAL, PARAMETER :: sday=86400.
! scaling velocity  **2
    REAL, PARAMETER :: u000=4900.
    REAL :: c1, pc, fac
    INTEGER :: i, j, k
    INTRINSIC LOG10
    INTRINSIC TANH
    INTRINSIC SQRT
    REAL :: arg1
    REAL :: result1
    IF (.NOT.rf_initialized) THEN
      ALLOCATE(rf(npz))
      ALLOCATE(rw(npz))
      IF (p_c .LE. 0.) THEN
        pc = pm(1)
      ELSE
        pc = p_c
      END IF
!if( gid==0 ) write(6,*) 'Rayleigh friction E-folding time [days]:'
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
!    allocate( u2f(isd:ied,jsd:jed,kmax) )
    CALL C2L_ORD2(u, v, ua, va, dx, dy, rdxa, rdya, npz)
    u2f = 0.
    DO k=1,kmax
      IF (hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            u2f(i, j, k) = ua(i, j, k)**2 + va(i, j, k)**2
          END DO
        END DO
      ELSE
        DO j=js,je
          DO i=is,ie
            u2f(i, j, k) = ua(i, j, k)**2 + va(i, j, k)**2 + w(i, j, k)&
&             **2
          END DO
        END DO
      END IF
    END DO
!call timing_on('COMM_TOTAL')
    call mpp_update_domains(u2f, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
    !CALL MPP_UPDATE_DOMAINS_DUMMY3(u2f, isd, ied, jsd, jed, kmax)
!call timing_off('COMM_TOTAL')
    DO k=1,kmax
      IF (conserve) THEN
        IF (hydrostatic) THEN
          DO j=js,je
            DO i=is,ie
              arg1 = u2f(i, j, k)/u000
              result1 = SQRT(arg1)
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)/(cp-rg*ptop/&
&               pm(k))*(1.-1./(1.+dt*rf(k)*result1)**2)
            END DO
          END DO
        ELSE
          DO j=js,je
            DO i=is,ie
              delz(i, j, k) = delz(i, j, k)/pt(i, j, k)
              arg1 = u2f(i, j, k)/u000
              result1 = SQRT(arg1)
              pt(i, j, k) = pt(i, j, k) + 0.5*u2f(i, j, k)/(cp-rg*ptop/&
&               pm(k))*(1.-1./(1.+dt*rf(k)*result1)**2)
              delz(i, j, k) = delz(i, j, k)*pt(i, j, k)
            END DO
          END DO
        END IF
      END IF
      DO j=js-1,je+1
        DO i=is-1,ie+1
          arg1 = u2f(i, j, k)/u000
          result1 = SQRT(arg1)
          u2f(i, j, k) = dt*rf(k)*result1
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          u(i, j, k) = u(i, j, k)/(1.+0.5*(u2f(i, j-1, k)+u2f(i, j, k)))
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v(i, j, k) = v(i, j, k)/(1.+0.5*(u2f(i-1, j, k)+u2f(i, j, k)))
        END DO
      END DO
      IF (.NOT.hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            w(i, j, k) = w(i, j, k)/(1.+u2f(i, j, k))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE RAYLEIGH_FRICTION
!  Differentiation of pertstate2fvstate in forward (tangent) mode:
!   variations   of useful results: peln pkz pe pk pt
!   with respect to varying inputs: peln delp pkz pe pt
!#endif
  SUBROUTINE PERTSTATE2FVSTATE_TLM(delp, delp_tl, pe, pe_tl, pk, pk_tl, &
&   pkz, pkz_tl, peln, peln_tl, ptop, pt, pt_tl, isd, ied, jsd, jed, is&
&   , ie, js, je, km, kappa)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed, is, ie, js, je, km
    REAL, INTENT(IN) :: kappa
    REAL, INTENT(IN) :: ptop
    REAL, INTENT(IN) :: delp(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: delp_tl(isd:ied, jsd:jed, km)
    REAL(p_precision), INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL(p_precision), INTENT(INOUT) :: pe_tl(is-1:ie+1, km+1, js-1:je+1&
&   )
    REAL(p_precision), INTENT(INOUT) :: pk(is:ie, js:je, km+1)
    REAL(p_precision), INTENT(INOUT) :: pk_tl(is:ie, js:je, km+1)
    REAL(p_precision), INTENT(INOUT) :: pkz(is:ie, js:je, km)
    REAL(p_precision), INTENT(INOUT) :: pkz_tl(is:ie, js:je, km)
    REAL(p_precision), INTENT(INOUT) :: peln(is:ie, km+1, js:je)
    REAL(p_precision), INTENT(INOUT) :: peln_tl(is:ie, km+1, js:je)
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: pt_tl(isd:ied, jsd:jed, km)
    INTEGER :: i, j, k, l
    REAL :: pke_local(isd:ied, jsd:jed, km+1)
    REAL :: pke_local_tl(isd:ied, jsd:jed, km+1)
    REAL :: bx_local(isd:ied, jsd:jed, km)
    REAL :: bx_local_tl(isd:ied, jsd:jed, km)
    REAL :: cx_local(isd:ied, jsd:jed, km)
    REAL :: cx_local_tl(isd:ied, jsd:jed, km)
    INTRINSIC LOG
    INTRINSIC EXP
    pe_tl(:, 1, :) = 0.0_8
    pe(:, 1, :) = ptop
    DO l=2,km+1
      DO j=js-1,je+1
        DO i=is-1,ie+1
          pe_tl(i, l, j) = pe_tl(i, l-1, j) + delp_tl(i, j, l-1)
          pe(i, l, j) = pe(i, l-1, j) + delp(i, j, l-1)
        END DO
      END DO
    END DO
    pke_local = 0.
    pke_local_tl = 0.0
    DO l=1,km+1
      DO j=js-1,je+1
        DO i=is-1,ie+1
          IF (pe(i, l, j) .NE. 0.) THEN
            pke_local_tl(i, j, l) = pe_tl(i, l, j)/pe(i, l, j)
            pke_local(i, j, l) = LOG(pe(i, l, j))
          END IF
        END DO
      END DO
    END DO
    DO l=1,km+1
      DO j=js,je
        DO i=is,ie
          peln_tl(i, l, j) = pke_local_tl(i, j, l)
          peln(i, l, j) = pke_local(i, j, l)
        END DO
      END DO
    END DO
    pk_tl = 0.0_8
    DO k=1,km+1
      DO j=js,je
        DO i=is,ie
          pk_tl(i, j, k) = kappa*peln_tl(i, k, j)*EXP(kappa*peln(i, k, j&
&           ))
          pk(i, j, k) = EXP(kappa*peln(i, k, j))
        END DO
      END DO
    END DO
    bx_local_tl = 0.0
    DO j=jsd,jed
      DO i=isd,ied
        bx_local_tl(i, j, :) = pke_local_tl(i, j, 2:km+1) - pke_local_tl&
&         (i, j, 1:km)
        bx_local(i, j, :) = pke_local(i, j, 2:km+1) - pke_local(i, j, 1:&
&         km)
      END DO
    END DO
    pke_local = 0.
    pke_local_tl = 0.0
    DO l=1,km+1
      DO j=js-1,je+1
        DO i=is-1,ie+1
          IF (pe(i, l, j) .GT. 0.0 .OR. (pe(i, l, j) .LT. 0.0 .AND. &
&             kappa .EQ. INT(kappa))) THEN
            pke_local_tl(i, j, l) = kappa*pe(i, l, j)**(kappa-1)*pe_tl(i&
&             , l, j)
          ELSE IF (pe(i, l, j) .EQ. 0.0 .AND. kappa .EQ. 1.0) THEN
            pke_local_tl(i, j, l) = pe_tl(i, l, j)
          ELSE
            pke_local_tl(i, j, l) = 0.0
          END IF
          pke_local(i, j, l) = pe(i, l, j)**kappa
        END DO
      END DO
    END DO
    cx_local_tl = 0.0
    DO j=jsd,jed
      DO i=isd,ied
        cx_local_tl(i, j, :) = pke_local_tl(i, j, 2:km+1) - pke_local_tl&
&         (i, j, 1:km)
        cx_local(i, j, :) = pke_local(i, j, 2:km+1) - pke_local(i, j, 1:&
&         km)
      END DO
    END DO
    DO l=1,km
      DO j=jsd,jed
        DO i=isd,ied
          IF (bx_local(i, j, l) .NE. 0.) THEN
            bx_local_tl(i, j, l) = -(kappa*bx_local_tl(i, j, l)/(kappa*&
&             bx_local(i, j, l))**2)
            bx_local(i, j, l) = 1./(kappa*bx_local(i, j, l))
          END IF
        END DO
      END DO
    END DO
    DO l=1,km
      DO j=js,je
        DO i=is,ie
          pkz_tl(i, j, l) = bx_local_tl(i, j, l)*cx_local(i, j, l) + &
&           bx_local(i, j, l)*cx_local_tl(i, j, l)
          pkz(i, j, l) = bx_local(i, j, l)*cx_local(i, j, l)
        END DO
      END DO
    END DO
!Convert potential temperature to temperature
    DO l=1,km
      DO j=js,je
        DO i=is,ie
          pt_tl(i, j, l) = pt_tl(i, j, l)*pkz(i, j, l) + pt(i, j, l)*&
&           pkz_tl(i, j, l)
          pt(i, j, l) = pt(i, j, l)*pkz(i, j, l)
        END DO
      END DO
    END DO
  END SUBROUTINE PERTSTATE2FVSTATE_TLM
!#endif
  SUBROUTINE PERTSTATE2FVSTATE(delp, pe, pk, pkz, peln, ptop, pt, isd, &
&   ied, jsd, jed, is, ie, js, je, km, kappa)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed, is, ie, js, je, km
    REAL, INTENT(IN) :: kappa
    REAL, INTENT(IN) :: ptop
    REAL, INTENT(IN) :: delp(isd:ied, jsd:jed, km)
    REAL(p_precision), INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL(p_precision), INTENT(INOUT) :: pk(is:ie, js:je, km+1)
    REAL(p_precision), INTENT(INOUT) :: pkz(is:ie, js:je, km)
    REAL(p_precision), INTENT(INOUT) :: peln(is:ie, km+1, js:je)
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, km)
    INTEGER :: i, j, k, l
    REAL :: pke_local(isd:ied, jsd:jed, km+1)
    REAL :: bx_local(isd:ied, jsd:jed, km)
    REAL :: cx_local(isd:ied, jsd:jed, km)
    INTRINSIC LOG
    INTRINSIC EXP
    pe(:, 1, :) = ptop
    DO l=2,km+1
      DO j=js-1,je+1
        DO i=is-1,ie+1
          pe(i, l, j) = pe(i, l-1, j) + delp(i, j, l-1)
        END DO
      END DO
    END DO
    pke_local = 0.
    DO l=1,km+1
      DO j=js-1,je+1
        DO i=is-1,ie+1
          IF (pe(i, l, j) .NE. 0.) pke_local(i, j, l) = LOG(pe(i, l, j))
        END DO
      END DO
    END DO
    DO l=1,km+1
      DO j=js,je
        DO i=is,ie
          peln(i, l, j) = pke_local(i, j, l)
        END DO
      END DO
    END DO
    DO k=1,km+1
      DO j=js,je
        DO i=is,ie
          pk(i, j, k) = EXP(kappa*peln(i, k, j))
        END DO
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied
        bx_local(i, j, :) = pke_local(i, j, 2:km+1) - pke_local(i, j, 1:&
&         km)
      END DO
    END DO
    pke_local = 0.
    DO l=1,km+1
      DO j=js-1,je+1
        DO i=is-1,ie+1
          pke_local(i, j, l) = pe(i, l, j)**kappa
        END DO
      END DO
    END DO
    DO j=jsd,jed
      DO i=isd,ied
        cx_local(i, j, :) = pke_local(i, j, 2:km+1) - pke_local(i, j, 1:&
&         km)
      END DO
    END DO
    DO l=1,km
      DO j=jsd,jed
        DO i=isd,ied
          IF (bx_local(i, j, l) .NE. 0.) bx_local(i, j, l) = 1./(kappa*&
&             bx_local(i, j, l))
        END DO
      END DO
    END DO
    DO l=1,km
      DO j=js,je
        DO i=is,ie
          pkz(i, j, l) = bx_local(i, j, l)*cx_local(i, j, l)
        END DO
      END DO
    END DO
!Convert potential temperature to temperature
    DO l=1,km
      DO j=js,je
        DO i=is,ie
          pt(i, j, l) = pt(i, j, l)*pkz(i, j, l)
        END DO
      END DO
    END DO
  END SUBROUTINE PERTSTATE2FVSTATE
!  Differentiation of fvstate2pertstate in forward (tangent) mode:
!   variations   of useful results: pt
!   with respect to varying inputs: pkz pt
  SUBROUTINE FVSTATE2PERTSTATE_TLM(pkz, pkz_tl, pt, pt_tl, isd, ied, jsd&
&   , jed, is, ie, js, je, km)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed, is, ie, js, je, km
    REAL(p_precision), INTENT(INOUT) :: pkz(is:ie, js:je, km)
    REAL(p_precision), INTENT(INOUT) :: pkz_tl(is:ie, js:je, km)
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: pt_tl(isd:ied, jsd:jed, km)
    INTEGER :: i, j, k, l
    DO l=1,km
      DO j=js,je
        DO i=is,ie
          IF (pkz(i, j, l) .NE. 0.) THEN
            pt_tl(i, j, l) = (pt_tl(i, j, l)*pkz(i, j, l)-pt(i, j, l)*&
&             pkz_tl(i, j, l))/pkz(i, j, l)**2
            pt(i, j, l) = pt(i, j, l)/pkz(i, j, l)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE FVSTATE2PERTSTATE_TLM
  SUBROUTINE FVSTATE2PERTSTATE(pkz, pt, isd, ied, jsd, jed, is, ie, js, &
&   je, km)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: isd, ied, jsd, jed, is, ie, js, je, km
    REAL(p_precision), INTENT(INOUT) :: pkz(is:ie, js:je, km)
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, km)
    INTEGER :: i, j, k, l
    DO l=1,km
      DO j=js,je
        DO i=is,ie
          IF (pkz(i, j, l) .NE. 0.) pt(i, j, l) = pt(i, j, l)/pkz(i, j, &
&             l)
        END DO
      END DO
    END DO
  END SUBROUTINE FVSTATE2PERTSTATE

end module fv_dynamics_tlm_mod
