module fv_dynamics_mod
!#ifndef MAPL_MODE
!   use constants_mod,       only: grav, pi, radius, hlv, rdgas, omega, rvgas, cp_vapor
!#else
   use MAPL_MOD
!#endif
   use fv_arrays_mod,  only: REAL4, REAL8, FVPRC, R_GRID, tiny_number, fv_timing_onoff
   use dyn_core_mod,        only: dyn_core, del2_cubed
   use fv_mapz_mod,         only: compute_total_energy, Lagrangian_to_Eulerian
   use fv_tracer2d_mod,     only: tracer_2d, tracer_2d_1L, tracer_2d_nested
   use fv_grid_utils_mod,   only: cubed_to_latlon, c2l_ord2, g_sum
   use fv_mp_mod,           only: is_master
   use fv_mp_mod,           only: group_halo_update_type
   use fv_mp_mod,           only: start_group_halo_update, complete_group_halo_update
   use fv_timing_mod,       only: timing_on, timing_off
!   use diag_manager_mod,    only: send_data
   use fv_diagnostics_mod,  only: fv_time, prt_mxm, range_check, prt_minmax
   use mpp_domains_mod,     only: DGRID_NE, CGRID_NE, mpp_update_domains, domain2D
   use field_manager_mod,   only: MODEL_ATMOS
   use tracer_manager_mod,  only: get_tracer_index
!   use fv_sg_mod,           only: neg_adj3
!   use fv_nesting_mod,      only: setup_nested_grid_BCs
   use fv_arrays_mod,       only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_grid_bounds_type
   use fv_nwp_nudge_mod,    only: do_adiabatic_init
!#ifdef MAPL_MODE
!   use fv_control_mod,      only: dyn_timer, comm_timer
!#endif

implicit none

!#ifdef MAPL_MODE
  ! Include the MPI library definitons:
!  include 'mpif.h'
!#endif

   logical :: RF_initialized = .false.
   logical :: bad_range
!   real(FVPRC), allocatable ::  rf(:)
   integer :: kmax=1
   real(FVPRC) :: agrav
!#ifdef HIWPP
!   real(FVPRC), allocatable:: u00(:,:,:), v00(:,:,:)
!#endif
private
public :: fv_dynamics

!---- version number -----
!   character(len=128) :: version = '$Id: fv_dynamics.F90,v 1.2.2.2.2.2.30.1.4.1.2.1.20.2.46.4.4.2.2.3.2.1.6.1.2.1 2017/02/16 03:47:47 aoloso Exp $'
!   character(len=128) :: tagname = '$Name: Heracles-UNSTABLE_ncepdyn_Feb222017 $'

!#ifdef MAPL_MODE
  real(FVPRC), parameter :: RADIUS       = MAPL_RADIUS
  real(FVPRC), parameter :: PI           = MAPL_PI_R8
  real(FVPRC), parameter :: RDGAS        = MAPL_RGAS
  real(FVPRC), parameter :: GRAV         = MAPL_GRAV
  real(FVPRC), parameter :: HLV          = MAPL_ALHL
  real(FVPRC), parameter :: CP_VAPOR     = MAPL_CP
  real(FVPRC), parameter :: RVGAS        = MAPL_RVAP
  real(FVPRC), parameter :: OMEGA        = MAPL_OMEGA
!#endif

  logical, save :: IdealTest = .false.

contains

!-----------------------------------------------------------------------
!     fv_dynamics :: FV dynamical core driver
!-----------------------------------------------------------------------
   subroutine fv_dynamics_dummy(npx, npy, npz, nq_tot,  ng, bdt, consv_te, fill,               &
                        reproduce_sum, kappa, cp_air, zvir, ptop, ks, ncnst, n_split,     &
                        q_split, u, v, w, delz, hydrostatic, pt, delp, q,   &
                        ps, pe, pk, peln, pkz, phis, q_con, omga, ua, va, uc, vc,          &
                        ak, bk, mfx, mfy, cx, cy, ze0, hybrid_z, &
                        gridstruct, flagstruct, neststruct, idiag, bd, &
                        parent_grid, domain, time_total)

       real(FVPRC), intent(IN) :: bdt  ! Large time-step
    real(FVPRC), intent(IN) :: consv_te
    real(FVPRC), intent(IN) :: kappa, cp_air
    real(FVPRC), intent(IN) :: zvir, ptop
    real(FVPRC), intent(IN), optional :: time_total

    integer, intent(IN) :: npx
    integer, intent(IN) :: npy
    integer, intent(IN) :: npz
    integer, intent(IN) :: nq_tot             ! transported tracers
    integer, intent(IN) :: ng
    integer, intent(IN) :: ks
    integer, intent(IN) :: ncnst
    integer, intent(IN) :: n_split        ! small-step horizontal dynamics
    integer, intent(IN) :: q_split        ! tracer
    logical, intent(IN) :: fill
    logical, intent(IN) :: reproduce_sum
    logical, intent(IN) :: hydrostatic
    logical, intent(IN) :: hybrid_z       ! Using hybrid_z for remapping

    type(fv_grid_bounds_type), intent(IN) :: bd
    real(FVPRC), intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u ! D grid zonal wind (m/s)
    real(FVPRC), intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v ! D grid meridional wind (m/s)
    !real(FVPRC), intent(inout) :: w( bd%isd:,bd%jsd:,1:)  !  W (m/s)
    real(FVPRC), intent(inout) :: w( bd%isd:bd%ied,bd%jsd:bd%jed,npz)  !  W (m/s) ! replaced previous line bma
    real(FVPRC), intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real(FVPRC), intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real(FVPRC), intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) ! specific humidity and constituents
    !real(FVPRC), intent(inout) :: delz(bd%isd:,bd%jsd:,1:)   ! delta-height (m); non-hydrostatic only
    real(FVPRC), intent(inout) :: delz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! replace previous line bma
    real(FVPRC), intent(inout) ::  ze0(bd%is:bd%ie,bd%js:bd%je,npz+1) ! height at edges (m); non-hydrostatic

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real(FVPRC), intent(inout) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)           ! Surface pressure (pascal)
    real(FVPRC), intent(inout) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)  ! edge pressure (pascal)
    real(FVPRC), intent(inout) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)          ! pe**kappa
    real(FVPRC), intent(inout) :: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)           ! ln(pe)
    real(FVPRC), intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             ! finite-volume mean pk
    !real(FVPRC), intent(inout):: q_con(bd%isd:, bd%jsd:, 1:)
    real(FVPRC), intent(inout):: q_con(bd%isd:bd%isd, bd%jsd:bd%jsd, 1) !bma replaced previous line
    
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real(FVPRC), intent(inout) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)       ! Surface geopotential (g*Z_surf)
    real(FVPRC), intent(inout) :: omga(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! Vertical pressure velocity (pa/s)
    real(FVPRC), intent(inout) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) ! (uc,vc) mostly used as the C grid winds
    real(FVPRC), intent(inout) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)

    real(FVPRC), intent(inout), dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz):: ua, va
    real(FVPRC), intent(in),    dimension(npz+1):: ak, bk

! Accumulated Mass flux arrays: the "Flux Capacitor"
    real(FVPRC), intent(inout) ::  mfx(bd%is:bd%ie+1, bd%js:bd%je,   npz)
    real(FVPRC), intent(inout) ::  mfy(bd%is:bd%ie  , bd%js:bd%je+1, npz)
! Accumulated Courant number arrays
    real(FVPRC), intent(inout) ::  cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    real(FVPRC), intent(inout) ::  cy(bd%isd:bd%ied ,bd%js:bd%je+1, npz)

    type(fv_grid_type),  intent(inout), target :: gridstruct
    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_nest_type),  intent(INOUT) :: neststruct
    type(domain2d), intent(INOUT) :: domain
    type(fv_atmos_type), intent(INOUT) :: parent_grid
    type(fv_diag_type), intent(IN) :: idiag



       call fv_dynamics(npx, npy, npz, nq_tot,  ng, bdt, consv_te, fill,               &
                        reproduce_sum, kappa, cp_air, zvir, ptop, ks, ncnst, n_split,     &
                        q_split, u, v, w, delz, hydrostatic, pt, delp, q,   &
                        ps, pe, pk, peln, pkz, phis, q_con, omga, ua, va, uc, vc,          &
                        ak, bk, mfx, mfy, cx, cy, ze0, hybrid_z, &
                        gridstruct, flagstruct, neststruct, idiag, bd, &
                        parent_grid, domain, time_total)


  end subroutine fv_dynamics_dummy

  subroutine fv_dynamics(npx, npy, npz, nq_tot,  ng, bdt, consv_te, fill,               &
                        reproduce_sum, kappa, cp_air, zvir, ptop, ks, ncnst, n_split,     &
                        q_split, u, v, w, delz, hydrostatic, pt, delp, q,   &
                        ps, pe, pk, peln, pkz, phis, q_con, omga, ua, va, uc, vc,          &
                        ak, bk, mfx, mfy, cx, cy, ze0, hybrid_z, &
                        gridstruct, flagstruct, neststruct, idiag, bd, &
                        parent_grid, domain, time_total)

    real(FVPRC), intent(IN) :: bdt  ! Large time-step
    real(FVPRC), intent(IN) :: consv_te
    real(FVPRC), intent(IN) :: kappa, cp_air
    real(FVPRC), intent(IN) :: zvir, ptop
    real(FVPRC), intent(IN), optional :: time_total

    integer, intent(IN) :: npx
    integer, intent(IN) :: npy
    integer, intent(IN) :: npz
    integer, intent(IN) :: nq_tot             ! transported tracers
    integer, intent(IN) :: ng
    integer, intent(IN) :: ks
    integer, intent(IN) :: ncnst
    integer, intent(IN) :: n_split        ! small-step horizontal dynamics
    integer, intent(IN) :: q_split        ! tracer
    logical, intent(IN) :: fill
    logical, intent(IN) :: reproduce_sum
    logical, intent(IN) :: hydrostatic
    logical, intent(IN) :: hybrid_z       ! Using hybrid_z for remapping

    type(fv_grid_bounds_type), intent(IN) :: bd
    real(FVPRC), intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u ! D grid zonal wind (m/s)
    real(FVPRC), intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v ! D grid meridional wind (m/s)
    !real(FVPRC), intent(inout) :: w( bd%isd:,bd%jsd:,1:)  !  W (m/s)
    real(FVPRC), intent(inout) :: w( bd%isd:bd%ied,bd%jsd:bd%jed,npz)  !  W (m/s) ! replaced previous line bma
    real(FVPRC), intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real(FVPRC), intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real(FVPRC), intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) ! specific humidity and constituents
    !real(FVPRC), intent(inout) :: delz(bd%isd:,bd%jsd:,1:)   ! delta-height (m); non-hydrostatic only
    real(FVPRC), intent(inout) :: delz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! replace previous line bma
    real(FVPRC), intent(inout) ::  ze0(bd%is:bd%ie,bd%js:bd%je,npz+1) ! height at edges (m); non-hydrostatic

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real(FVPRC), intent(inout) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)           ! Surface pressure (pascal)
    real(FVPRC), intent(inout) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)  ! edge pressure (pascal)
    real(FVPRC), intent(inout) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)          ! pe**kappa
    real(FVPRC), intent(inout) :: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)           ! ln(pe)
    real(FVPRC), intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             ! finite-volume mean pk
    !real(FVPRC), intent(inout):: q_con(bd%isd:, bd%jsd:, 1:)
    real(FVPRC), intent(inout):: q_con(bd%isd:bd%isd, bd%jsd:bd%jsd, 1) !bma replaced previous line
    
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real(FVPRC), intent(inout) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)       ! Surface geopotential (g*Z_surf)
    real(FVPRC), intent(inout) :: omga(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! Vertical pressure velocity (pa/s)
    real(FVPRC), intent(inout) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) ! (uc,vc) mostly used as the C grid winds
    real(FVPRC), intent(inout) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)

    real(FVPRC), intent(inout), dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz):: ua, va
    real(FVPRC), intent(in),    dimension(npz+1):: ak, bk

! Accumulated Mass flux arrays: the "Flux Capacitor"
    real(FVPRC), intent(inout) ::  mfx(bd%is:bd%ie+1, bd%js:bd%je,   npz)
    real(FVPRC), intent(inout) ::  mfy(bd%is:bd%ie  , bd%js:bd%je+1, npz)
! Accumulated Courant number arrays
    real(FVPRC), intent(inout) ::  cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    real(FVPRC), intent(inout) ::  cy(bd%isd:bd%ied ,bd%js:bd%je+1, npz)

    type(fv_grid_type),  intent(inout), target :: gridstruct
    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_nest_type),  intent(INOUT) :: neststruct
    type(domain2d), intent(INOUT) :: domain
    type(fv_atmos_type), intent(INOUT) :: parent_grid
    type(fv_diag_type), intent(IN) :: idiag

    real(FVPRC), parameter:: c_liq = 4190.       ! heat capacity of water at 0C
    real(FVPRC), parameter:: c_ice = 2106.       ! heat capacity of ice at 0C: c=c_ice+7.3*(T-Tice) 
    real(FVPRC), parameter:: cv_vap = cp_vapor - rvgas  ! 1384.5
! Local Arrays
      real(FVPRC):: ws(bd%is:bd%ie,bd%js:bd%je)
      real(FVPRC):: te_2d(bd%is:bd%ie,bd%js:bd%je)
      real(FVPRC)::   teq(bd%is:bd%ie,bd%js:bd%je)
      real(FVPRC):: ps2(bd%isd:bd%ied,bd%jsd:bd%jed)
      real(FVPRC):: m_fac(bd%is:bd%ie,bd%js:bd%je)
      real(FVPRC):: pfull(npz)
      !real(FVPRC):: gz(bd%is:bd%ie)
      real(FVPRC) :: dp1(bd%is:bd%ie,bd%js:bd%je,1:npz),dtdt_m(bd%is:bd%ie,bd%js:bd%je,npz),cappa(bd%isd:bd%isd,bd%jsd:bd%jsd,1)
      real(FVPRC):: akap, rg, rdg, ph1, ph2, mdt, gam, amdt, u0
      integer :: i,j,k, n, iq, n_map, nq, nwat, k_split
      integer :: sphum, liq_wat, ice_wat      ! GFDL physics
      integer :: rainwat, snowwat, graupel, cld_amt
      logical used, last_step, consv_am, do_omega
      integer, parameter :: max_packs=24
      type(group_halo_update_type), save :: i_pack(max_packs)
      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed
      real(FVPRC) :: rcv, dt2, consv_fac, q_liq, q_sol, cvm
      real(FVPRC):: cv_air
      real(kind=8) :: t1, t2
      integer :: status

      real(FVPRC) ::     gz(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz+1)
      real(FVPRC) ::    pkc(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz+1)
      real(FVPRC) ::    ptc(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz  )
      real(FVPRC) ::    crx(bd%is :bd%ie +1,bd%jsd:bd%jed  ,npz  )
      real(FVPRC) ::    xfx(bd%is :bd%ie +1,bd%jsd:bd%jed  ,npz  )
      real(FVPRC) ::    cry(bd%isd:bd%ied  ,bd%js :bd%je +1,npz  )
      real(FVPRC) ::    yfx(bd%isd:bd%ied  ,bd%js :bd%je +1,npz  )
      real(FVPRC) ::  divgd(bd%isd:bd%ied+1,bd%jsd:bd%jed+1,npz  )
      real(FVPRC) ::  delpc(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz  )
      real(FVPRC) ::     ut(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz  )
      real(FVPRC) ::     vt(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz  )
      real(FVPRC) ::     zh(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz+1)
      real(FVPRC) ::    pk3(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz+1)
      real(FVPRC) ::     du(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz  )
      real(FVPRC) ::     dv(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz  )
      integer :: hord_tr, hord_tr_pert
      integer :: kord_mt, kord_wz, kord_tr(ncnst), kord_tm
      integer :: kord_mt_pert, kord_wz_pert, kord_tr_pert(ncnst), kord_tm_pert

         gz = 0.0
        pkc = 0.0
        ptc = 0.0
        crx = 0.0
        xfx = 0.0
        cry = 0.0
        yfx = 0.0
      divgd = 0.0
      delpc = 0.0
         ut = 0.0
         vt = 0.0
         zh = 0.0
        pk3 = 0.0
         du = 0.0
         dv = 0.0

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

      !Compute the FV variables internally, for checkpointing purposes
      if ( hydrostatic .and. .not.(IdealTest) ) then
         CALL GEOS_to_FV3(bd, npz, kappa, ptop, delp, pe, pk, pkz, peln, pt)
      endif

!      dyn_timer = 0
!      comm_timer = 0

      cv_air =  cp_air - rdgas
      agrav = 1. / grav
        dt2 = 0.5*bdt
      consv_am = flagstruct%consv_am
      k_split = flagstruct%k_split
      nwat = flagstruct%nwat
      nq = nq_tot - flagstruct%dnats

      if ( flagstruct%no_dycore ) goto 911

!#ifdef MAPL_MODE
!! Begin Dynamics timer for GEOS history processing
!      t1 = MPI_Wtime(status)
!#endif

      !allocate ( dp1(is:ie, js:je, 1:npz) )

!#ifdef MOIST_CAPPA
!      allocate ( cappa(isd:ied,jsd:jed,npz) )
!#else
!      allocate ( cappa(isd:isd,jsd:jsd,1) )
!#endif

!#ifdef SW_DYNAMICS
!      akap  = 1.
!      pfull(1) = 0.5*flagstruct%p_ref
!#ifdef TEST_TRACER
!      sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
!#endif
!#else

      if (nwat>=3 ) then
             sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
           liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
           ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
           cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
      endif
      if ( nwat==6 ) then
           rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
           snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
           graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
      else
           sphum = 1
           cld_amt = -1   ! to cause trouble if (mis)used
      endif

      akap  = kappa
      rg = kappa*cp_air
      rcv = 1./(cp_air-rg)

!$OMP parallel do default(none) shared(npz,ak,bk,flagstruct,pfull) &
!$OMP                          private(ph1, ph2)
      do k=1,npz
         ph1 = ak(k  ) + bk(k  )*flagstruct%p_ref
         ph2 = ak(k+1) + bk(k+1)*flagstruct%p_ref
         pfull(k) = (ph2 - ph1) / log(ph2/ph1)
      enddo

    if ( hydrostatic ) then
!$OMP parallel do default(none) shared(is,ie,js,je,npz,dp1,zvir,q,q_con,sphum,liq_wat, &
!$OMP                                  rainwat,ice_wat,snowwat,graupel )
      do k=1,npz
         do j=js,je
            do i=is,ie
               dp1(i,j,k) = zvir*q(i,j,k,sphum)
!#ifdef USE_COND
!!#ifdef USE_NWAT3
!!               q_con(i,j,k) = q(i,j,k,liq_wat) + q(i,j,k,ice_wat)
!!#else
!               q_con(i,j,k) = q(i,j,k,liq_wat) + q(i,j,k,rainwat) + q(i,j,k,ice_wat)  &
!                            + q(i,j,k,snowwat) + q(i,j,k,graupel)
!!#endif
!#endif
            enddo
         enddo
      enddo
    else
       rdg = -rdgas * agrav
!$OMP parallel do default(none) shared(is,ie,js,je,npz,dp1,zvir,q,q_con,sphum,liq_wat, &
!$OMP                                  rainwat,ice_wat,snowwat,graupel,pkz,cv_air,     & 
!$OMP                                  cappa,kappa,rdg,delp,pt,delz)                         &
!$OMP                          private(cvm, q_liq, q_sol)
!$AD II-LOOP
       do k=1,npz
          do j=js,je
             do i=is,ie
                dp1(i,j,k) = zvir*q(i,j,k,sphum)
!#ifdef USE_COND
!!#ifdef USE_NWAT3
!!               q_liq = q(i,j,k,liq_wat)
!!               q_sol = q(i,j,k,ice_wat)
!!#else
!               q_liq = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
!               q_sol = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
!!#endif
!               q_con(i,j,k) = q_liq + q_Sol
!!#ifdef MOIST_CAPPA
!!!-------Simplified form ---------------------------------------------------------------------------------------
!!!              cappa(i,j,k) = rdgas/(rdgas+((1.-q(i,j,k,sphum)-q_con(i,j,k))*cv_air + q(i,j,k,sphum)*cv_vap + &
!!!                                               q_liq*c_liq + q_sol*c_ice )/(1.+dp1(i,j,k)))
!!!-------Simplified form ---------------------------------------------------------------------------------------
!!               cvm = (1.-(q(i,j,k,sphum)+q_Con(i,j,k)))*cv_air+q(i,j,k,sphum)*cv_vap+q_liq*c_liq+q_sol*c_ice
!!               cappa(i,j,k) = rdgas/(rdgas + cvm*(1.-q_con(i,j,k))/(1.+zvir*q(i,j,k,sphum)-q_con(i,j,k)))
!!               pkz(i,j,k) = exp(cappa(i,j,k)*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
!!                            (1.+dp1(i,j,k)-q_con(i,j,k))/delz(i,j,k)) )
!!#else
!               pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
!                            (1.+dp1(i,j,k)-q_con(i,j,k))/delz(i,j,k)) )
!!#endif
!
!#else
               pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
                            (1.+dp1(i,j,k))/delz(i,j,k)) )
!#endif
             enddo
          enddo
       enddo
    endif

!      if ( flagstruct%fv_debug ) then
!!#ifdef MOIST_CAPPA
!!         call prt_mxm('cappa', cappa, is, ie, js, je, ng, npz, 1._FVPRC, gridstruct%area_64, domain)
!!#endif
!         call prt_mxm('PS',        ps, is, ie, js, je, ng,   1, 0.01_FVPRC, gridstruct%area_64, domain)
!         call prt_mxm('T_dyn_b',   pt, is, ie, js, je, ng, npz, 1._FVPRC,   gridstruct%area_64, domain)
!         call prt_mxm('delz',    delz, is, ie, js, je, ng, npz, 1._FVPRC, gridstruct%area_64, domain)
!         call prt_mxm('delp_b ', delp, is, ie, js, je, ng, npz, 0.01_FVPRC, gridstruct%area_64, domain)
!         call prt_mxm('pk_b',    pk, is, ie, js, je, 0, npz+1, 1._FVPRC,gridstruct%area_64, domain)
!         call prt_mxm('pkz_b',   pkz,is, ie, js, je, 0, npz,   1._FVPRC,gridstruct%area_64, domain)
!      endif

!---------------------
! Compute Total Energy
!---------------------
      if ( consv_te > 0.  .and. (.not.do_adiabatic_init) ) then
           call compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, npz,        &
                                     u, v, w, delz, pt, delp, q, dp1, pe, peln, phis, &
                                     gridstruct%rsin2, gridstruct%cosa_s, &
                                     zvir, cp_air, rg, hlv, te_2d, ua, va, teq,        &
                                     flagstruct%moist_phys, sphum, liq_wat, rainwat,   &
                                     ice_wat, snowwat, graupel, hydrostatic, idiag%id_te)
           if( idiag%id_te>0 ) then
               !used = send_data(idiag%id_te, teq, fv_time)
!              te_den=1.E-9*g_sum(teq, is, ie, js, je, ng, area, 0)/(grav*4.*pi*radius**2)
!              if(is_master())  write(*,*) 'Total Energy Density (Giga J/m**2)=',te_den
           endif
      endif

      if( (consv_am.or.idiag%id_amdt>0) .and. (.not.do_adiabatic_init) ) then
          call compute_aam(npz, is, ie, js, je, isd, ied, jsd, jed, gridstruct, bd,   &
                           ptop, ua, va, u, v, delp, teq, ps2, m_fac)
      endif

      if( flagstruct%tau > 0. ) then
        if ( gridstruct%grid_type<4 ) then
             call Rayleigh_Super(abs(bdt), npx, npy, npz, ks, pfull, phis, flagstruct%tau, u, v, w, pt,  &
                  ua, va, delz, gridstruct%agrid, cp_air, rg, ptop, hydrostatic, .true., &
                  flagstruct%rf_cutoff, gridstruct, domain, bd)
        else
             call Rayleigh_Friction(abs(bdt), npx, npy, npz, ks, pfull, flagstruct%tau, u, v, w, pt,  &
                  ua, va, delz, cp_air, rg, ptop, hydrostatic, .true., flagstruct%rf_cutoff, gridstruct, domain, bd)
        endif
      endif

!#endif

      !We call this BEFORE converting pt to virtual potential temperature, 
      !since we interpolate on (regular) temperature rather than theta.
!      if (gridstruct%nested .or. ANY(neststruct%child_grids)) then
!                                           if (fv_timing_onoff) call timing_on('NEST_BCs')
!         call setup_nested_grid_BCs(npx, npy, npz, cp_air, zvir, ncnst, sphum,     &
!              u, v, w, pt, delp, delz, q, uc, vc, pkz, &
!              neststruct%nested, flagstruct%inline_q, flagstruct%make_nh, ng, &
!              gridstruct, flagstruct, neststruct, &
!              neststruct%nest_timestep, neststruct%tracer_nest_timestep, domain, bd)
!                                           if (fv_timing_onoff) call timing_off('NEST_BCs')
!      endif

!#ifndef SW_DYNAMICS
! Convert pt to virtual potential temperature * CP
  if ( hydrostatic ) then
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,cp_air,dp1,pkz)
!$AD II-LOOP
   do k=1,npz
     do j=js,je
        do i=is,ie
           pt(i,j,k) = cp_air*pt(i,j,k)*(1.+dp1(i,j,k))/pkz(i,j,k)
        enddo
     enddo
  enddo
  else
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,cp_air,dp1,pkz,q_con)
!$AD II-LOOP
   do k=1,npz
     do j=js,je
        do i=is,ie
!#ifdef USE_COND
!           pt(i,j,k) = cp_air*pt(i,j,k)*(1.+dp1(i,j,k)-q_con(i,j,k))/pkz(i,j,k)
!#else
           pt(i,j,k) = cp_air*pt(i,j,k)*(1.+dp1(i,j,k))/pkz(i,j,k)
!#endif
        enddo
     enddo
  enddo
!#endif
  endif

  last_step = .false.
  mdt = bdt / real(k_split)

  if ( idiag%id_mdt > 0 ) then
       !allocate ( dtdt_m(is:ie,js:je,npz) )
!$OMP parallel do default(none) shared(is,ie,js,je,npz,dtdt_m)
       do k=1,npz
          do j=js,je
             do i=is,ie
                dtdt_m(i,j,k) = 0.
             enddo
          enddo
       enddo
  endif


                                                  !if (fv_timing_onoff) call timing_on('  FV_DYN_LOOP')
  do n_map=1, k_split   ! first level of time-split
                                           !if (fv_timing_onoff) call timing_on('COMM_TOTAL')
!#ifdef USE_COND
!      call start_group_halo_update(i_pack(11), i_pack(11+12), q_con, domain)
!!#ifdef MOIST_CAPPA
!!      call start_group_halo_update(i_pack(12), i_pack(12+12), cappa, domain)
!!#endif
!#endif
      call start_group_halo_update(i_pack(1), i_pack(1+12), delp, domain)
      call start_group_halo_update(i_pack(2), i_pack(2+12), pt,   domain)
!#ifndef ROT3
      call start_group_halo_update(i_pack(8), i_pack(8+12), u, v, domain, gridtype=DGRID_NE)
!#endif
                                           !if (fv_timing_onoff) call timing_off('COMM_TOTAL')
!$OMP parallel do default(none) shared(is,ie,js,je,npz,dp1,delp)
!$AD II-LOOP
       do k=1,npz
         do j=js,je
            do i=is,ie
               dp1(i,j,k) = delp(i,j,k)
            enddo
         enddo
      enddo

      if ( n_map==k_split ) last_step = .true.

!#ifdef USE_COND
!                                           if (fv_timing_onoff) call timing_on('COMM_TOTAL')
!     call complete_group_halo_update(i_pack(11), i_pack(11+12), domain)
!!#ifdef MOIST_CAPPA
!!     call complete_group_halo_update(i_pack(12), i_pack(12+12), domain)
!!#endif
!                                           if (fv_timing_onoff) call timing_off('COMM_TOTAL')
!#endif

                                           if (fv_timing_onoff) call timing_on('  DYN_CORE')
      call dyn_core(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir, cp_air, akap, cappa, grav, hydrostatic, &
                    u, v, w, delz, pt, q, delp, pe, pk, phis, ws, omga, ptop, pfull, ua, va,           & 
                    uc, vc, mfx, mfy, cx, cy, pkz, peln, q_con, ak, bk, ks, &
                    gridstruct, flagstruct, neststruct, idiag, bd, &
                    domain, n_map==1, i_pack, last_step, &
                    gz,pkc,ptc,crx,xfx,cry,yfx,divgd,delpc,ut,vt,zh,pk3,du,dv,time_total)
                                           if (fv_timing_onoff) call timing_off('  DYN_CORE')

!  if ( flagstruct%fv_debug ) then
!       call prt_mxm('delp_a1',  delp, is, ie, js, je, ng, npz, 0.01_FVPRC, gridstruct%area_64, domain)
!       call prt_mxm('PT_dyn_a1',  pt, is, ie, js, je, ng, npz, 1._FVPRC, gridstruct%area_64, domain)
!       call prt_mxm('pk_a1',   pk, is, ie, js, je, 0, npz+1, 1._FVPRC, gridstruct%area_64, domain)
!  endif

!#ifdef SW_DYNAMICS
!!$OMP parallel do default(none) shared(is,ie,js,je,delp,agrav)
!      do j=js,je
!         do i=is,ie
!            ps(i,j) = delp(i,j,1) * agrav
!         enddo
!      enddo
!#else
      if( .not. flagstruct%inline_q .and. nq /= 0 ) then    
!--------------------------------------------------------
! Perform large-time-step scalar transport using the accumulated CFL and
! mass fluxes
                                              if (fv_timing_onoff) call timing_on('  tracer_2d')
       hord_tr = flagstruct%hord_tr
       hord_tr_pert = flagstruct%hord_tr_pert
       if (last_step) hord_tr = hord_tr_pert
       !!! CLEANUP: merge these two calls?
       if (gridstruct%nested .or. ANY(neststruct%child_grids)) then
         call tracer_2d_nested(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz, nq,    &
                        hord_tr, hord_tr_pert, q_split, mdt, idiag%id_divg, i_pack(10), i_pack(10+12), &
                        flagstruct%z_tracer, k_split, neststruct, parent_grid)          
       else
         call tracer_2d(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz, 1, nq,    &
                        hord_tr, hord_tr_pert, q_split, mdt, idiag%id_divg, i_pack(10), i_pack(10+12), &
                        flagstruct%z_tracer, k_split)
       endif
                                             if (fv_timing_onoff) call timing_off('  tracer_2d')

         if( last_step .and. idiag%id_divg>0 ) then
             !used = send_data(idiag%id_divg, dp1, fv_time) 
!             if(flagstruct%fv_debug) call prt_mxm('divg',  dp1, is, ie, js, je, 0, npz, 1._FVPRC,gridstruct%area_64, domain)
         endif
      endif

      if ( npz > 4 ) then
!------------------------------------------------------------------------
! Peroform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.
!------------------------------------------------------------------------

         do iq=1,nq
                                kord_tr     (iq) = flagstruct%kord_tr
                                kord_tr_pert(iq) = flagstruct%kord_tr_pert
            if ( iq==cld_amt )  kord_tr     (iq) = 9      ! monotonic
            if ( iq==cld_amt )  kord_tr_pert(iq) = 111
         enddo

         kord_mt = flagstruct%kord_mt
         kord_wz = flagstruct%kord_wz
         kord_tm = flagstruct%kord_tm
         kord_mt_pert = flagstruct%kord_mt_pert
         kord_wz_pert = flagstruct%kord_wz_pert
         kord_tm_pert = flagstruct%kord_tm_pert

         if ((kord_tm < 0 .and. kord_tm_pert >= 0) .or. (kord_tm >= 0 .and. kord_tm_pert < 0)) then
            !This would result in different anticipated paths through the code
            kord_tm = kord_tm_pert
         endif

         if (last_step) then
            !Does not matter how trajectory is remapped as it is about to be overwritten
            kord_mt = kord_mt_pert
            kord_wz = kord_wz_pert
            kord_tr = kord_tr_pert
            kord_tm = kord_tm_pert
         endif

         do_omega = hydrostatic .and. last_step
                                                  if (fv_timing_onoff) call timing_on('  Remapping')

         call Lagrangian_to_Eulerian(last_step, consv_te, ps, pe, delp,          &
                     pkz, pk, mdt, bdt, npz, is,ie,js,je, isd,ied,jsd,jed,       &
                     nq, nwat, sphum, q_con, u,  v, w, delz, pt, q, phis,    &
                     zvir, cp_air, akap, cappa, kord_mt, kord_wz, &
                     kord_tr, kord_tm, kord_mt_pert, kord_wz_pert, &
                     kord_tr_pert, kord_tm_pert, peln, te_2d,               &
                     ng, ua, va, omga, dp1, ws, fill, reproduce_sum,             &
                     idiag%id_mdt>0, dtdt_m, &
                     ptop, ak, bk, gridstruct, domain, ze0, flagstruct%gmao_cubic, flagstruct%remap_t,  &
                     flagstruct%do_sat_adj, hydrostatic, hybrid_z, do_omega, do_adiabatic_init)

                                                  if (fv_timing_onoff) call timing_off('  Remapping')
         if( last_step )  then
            if( hydrostatic ) then
!--------------------------
! Filter omega for physics:
!--------------------------
                if(flagstruct%nf_omega>0)    &
                call del2_cubed(omga, 0.20*gridstruct%da_min, gridstruct, domain, npx, npy, npz, flagstruct%nf_omega, bd)
            else
!$OMP parallel do default(none) shared(is,ie,js,je,npz,omga,delp,delz,w)
!$AD II-LOOP
                do k=1,npz
                  do j=js,je
                     do i=is,ie
                        omga(i,j,k) = delp(i,j,k)/delz(i,j,k)*w(i,j,k)
                     enddo
                  enddo
               enddo
            endif

! Convert back to temperature
            if ( .not. flagstruct%remap_t ) then
              if ( hydrostatic ) then
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,pkz,cp_air,zvir,q,sphum)
!$AD II-LOOP
                do k=1,npz
                  do j=js,je
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)/(cp_air*(1.+zvir*q(i,j,k,sphum)))
                  enddo
                  enddo
               enddo
              else
!$OMP parallel do default(none)  shared(is,ie,js,je,npz,pt,pkz,cp_air,zvir,q,sphum,q_con)
!$AD II-LOOP
                do k=1,npz
                  do j=js,je
                  do i=is,ie
!#ifdef USE_COND
!                     pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)/(cp_air*(1.+zvir*q(i,j,k,sphum)-q_con(i,j,k)))
!#else
                     pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)/(cp_air*(1.+zvir*q(i,j,k,sphum)))
!#endif
                  enddo
                  enddo
               enddo
              endif
            endif

         endif
      end if
!#endif
  enddo    ! n_map loop
                                                  !if (fv_timing_onoff) call timing_off('  FV_DYN_LOOP')
  if ( idiag%id_mdt > 0 .and. (.not.do_adiabatic_init) ) then
! Output temperature tendency due to inline moist physics:
!$OMP parallel do default(none) shared(is,ie,js,je,npz,dtdt_m,bdt)
       do k=1,npz
          do j=js,je
             do i=is,ie
                dtdt_m(i,j,k) = dtdt_m(i,j,k) / bdt
             enddo
          enddo
       enddo
       !used = send_data(idiag%id_mdt, dtdt_m, fv_time)
       !deallocate ( dtdt_m )
  endif

!  if( nwat==6 ) then
!      call neg_adj3(is, ie, js, je, ng, npz,        &
!                    flagstruct%hydrostatic,         &
!                    peln, delz,                     &
!                    pt, delp, q(isd,jsd,1,sphum),   &
!                              q(isd,jsd,1,liq_wat), &
!                              q(isd,jsd,1,rainwat), &
!                              q(isd,jsd,1,ice_wat), &
!                              q(isd,jsd,1,snowwat), &
!                              q(isd,jsd,1,graupel), &
!                              q(isd,jsd,1,cld_amt), flagstruct%check_negative)
!!     if ( flagstruct%fv_debug ) then
!!       call prt_mxm('T_dyn_a3',    pt, is, ie, js, je, ng, npz, 1._FVPRC, gridstruct%area_64, domain)
!!       call prt_mxm('SPHUM_dyn',   q(isd,jsd,1,sphum  ), is, ie, js, je, ng, npz, 1._FVPRC,gridstruct%area_64, domain)
!!       call prt_mxm('liq_wat_dyn', q(isd,jsd,1,liq_wat), is, ie, js, je, ng, npz, 1._FVPRC,gridstruct%area_64, domain)
!!       call prt_mxm('rainwat_dyn', q(isd,jsd,1,rainwat), is, ie, js, je, ng, npz, 1._FVPRC,gridstruct%area_64, domain)
!!       call prt_mxm('ice_wat_dyn', q(isd,jsd,1,ice_wat), is, ie, js, je, ng, npz, 1._FVPRC,gridstruct%area_64, domain)
!!       call prt_mxm('snowwat_dyn', q(isd,jsd,1,snowwat), is, ie, js, je, ng, npz, 1._FVPRC,gridstruct%area_64, domain)
!!       call prt_mxm('graupel_dyn', q(isd,jsd,1,graupel), is, ie, js, je, ng, npz, 1._FVPRC,gridstruct%area_64, domain)
!!     endif
!  endif

  if( consv_am .or. idiag%id_amdt>0 .or. idiag%id_aam>0 .and. (.not.do_adiabatic_init)  ) then
      call compute_aam(npz, is, ie, js, je, isd, ied, jsd, jed, gridstruct, bd,   &
                       ptop, ua, va, u, v, delp, te_2d, ps, m_fac)
      if( idiag%id_aam>0 ) then
          !used = send_data(idiag%id_aam, te_2d, fv_time)
          if ( prt_minmax ) then
             gam = g_sum( domain, te_2d, is, ie, js, je, ng, gridstruct%area_64, 0) 
             !if( is_master() ) write(6,*) 'Total AAM =', gam
          endif
      endif
  endif

  if( consv_am .or. idiag%id_amdt>0 .and. (.not.do_adiabatic_init)  ) then
!$OMP parallel do default(none) shared(is,ie,js,je,te_2d,teq,dt2,ps2,ps,idiag) 
      do j=js,je
         do i=is,ie
! Note: the mountain torque computation contains also numerical error
! The numerical error is mostly from the zonal gradient of the terrain (zxg)
            te_2d(i,j) = te_2d(i,j)-teq(i,j) + dt2*(ps2(i,j)+ps(i,j))*idiag%zxg(i,j)
         enddo
      enddo
      !if( idiag%id_amdt>0 ) used = send_data(idiag%id_amdt, te_2d/bdt, fv_time)

      if ( consv_am .or. prt_minmax ) then
         amdt = g_sum( domain, te_2d, is, ie, js, je, ng, gridstruct%area_64, 0) 
         u0 = -radius*amdt/g_sum( domain, m_fac, is, ie, js, je, ng, gridstruct%area_64, 0)
         u0 = real(u0, 4)    ! truncate to enforce reproducibility
         !if(is_master() .and. prt_minmax)         &
!        write(6,*) 'Dynamic AM tendency =', amdt/(bdt*1.e18), 'del-u (per yr)=', u0*365.*86400./bdt
!         write(6,*) 'Dynamic Angular Momentum tendency (Hadleys)=', amdt/(bdt*1.e18)
      endif

      if( consv_am ) then
!$OMP parallel do default(none) shared(is,ie,js,je,m_fac,u0,gridstruct)
      do j=js,je
         do i=is,ie
            m_fac(i,j) = u0*cos(gridstruct%agrid(i,j,2))
         enddo
      enddo
!$OMP parallel do default(none) shared(is,ie,js,je,npz,hydrostatic,pt,m_fac,ua,cp_air, &
!$OMP                                  rcv,u,u0,gridstruct,v )
      if ( hydrostatic ) then
!$AD II-LOOP
       do k=1,npz
         do j=js,je
         do i=is,ie
            pt(i,j,k) = pt(i,j,k) - m_fac(i,j)*(0.5*m_fac(i,j)+ua(i,j,k))/cp_air
         enddo
         enddo
      enddo
      else
!$AD II-LOOP
       do k=1,npz
         do j=js,je
         do i=is,ie
            pt(i,j,k) = pt(i,j,k) - m_fac(i,j)*(0.5*m_fac(i,j)+ua(i,j,k))*rcv
         enddo
         enddo
      enddo
      endif
      do k=1,npz
      do j=js,je+1
         do i=is,ie
            u(i,j,k) = u(i,j,k) + u0*gridstruct%l2c_u(i,j)
         enddo
      enddo
      do j=js,je
         do i=is,ie+1
            v(i,j,k) = v(i,j,k) + u0*gridstruct%l2c_v(i,j)
         enddo
      enddo
      enddo
      endif   !  consv_am
  endif

911  call cubed_to_latlon(u, v, ua, va, gridstruct, &
          npx, npy, npz, 1, gridstruct%grid_type, domain, gridstruct%nested, flagstruct%c2l_ord, bd)

!     if ( flagstruct%fv_debug ) then
!       call prt_mxm('UA', ua, is, ie, js, je, ng, npz, 1._FVPRC, gridstruct%area_64, domain)
!       call prt_mxm('VA', va, is, ie, js, je, ng, npz, 1._FVPRC, gridstruct%area_64, domain)
!       call prt_mxm('TA', pt, is, ie, js, je, ng, npz, 1._FVPRC, gridstruct%area_64, domain)
!     endif

!  if ( flagstruct%range_warn ) then
!       call range_check('UA_dyn', ua, is, ie, js, je, ng, npz, gridstruct%agrid,   &
!                         -220._FVPRC, 260._FVPRC, bad_range)
!       call range_check('VA_dyn', ua, is, ie, js, je, ng, npz, gridstruct%agrid,   &
!                         -220._FVPRC, 220._FVPRC, bad_range)
!!#ifndef SW_DYNAMICS
!       call range_check('TA_dyn', pt, is, ie, js, je, ng, npz, gridstruct%agrid,   &
!                         160._FVPRC, 330._FVPRC, bad_range)
!       if ( .not. hydrostatic ) &
!       call range_check('W_dyn', w, is, ie, js, je, ng, npz, gridstruct%agrid,   &
!                         -20._FVPRC, 20._FVPRC, bad_range)
!!#endif
!
!  endif

  !deallocate ( dp1 )
  !deallocate ( cappa )

  !Convert back to potential temperature
  if ( hydrostatic .and. .not.(IdealTest)) then
     CALL FV3_to_GEOS(bd, npz, pkz, pt)
  endif

!#ifdef MAPL_MODE
!  t2 = MPI_Wtime(status)
!  dyn_timer = dyn_timer + (t2-t1)
!#endif

  end subroutine fv_dynamics


 subroutine Rayleigh_Super(dt, npx, npy, npz, ks, pm, phis, tau, u, v, w, pt,  &
                           ua, va, delz, agrid, cp, rg, ptop, hydrostatic, conserve, rf_cutoff, gridstruct, domain, bd)
    real(FVPRC), intent(in):: dt
    real(FVPRC), intent(in):: tau              ! time scale (days)
    real(FVPRC), intent(in):: cp, rg, ptop, rf_cutoff
    real(FVPRC), intent(in),  dimension(npz):: pm
    integer, intent(in):: npx, npy, npz, ks
    logical, intent(in):: hydrostatic
    logical, intent(in):: conserve
    type(fv_grid_bounds_type), intent(IN) :: bd
    real(FVPRC), intent(inout):: u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) ! D grid zonal wind (m/s)
    real(FVPRC), intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz) ! D grid meridional wind (m/s)
    real(FVPRC), intent(inout)::  w(bd%isd:bd%isd,bd%jsd:bd%jsd,npz) ! cell center vertical wind (m/s)
    real(FVPRC), intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! temp
    real(FVPRC), intent(inout):: ua(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real(FVPRC), intent(inout):: va(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real(FVPRC), intent(inout):: delz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! delta-height (m); non-hydrostatic only
    real(FVPRC),   intent(in) :: agrid(bd%isd:bd%ied,  bd%jsd:bd%jed,2)
    real(FVPRC), intent(in) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)       ! Surface geopotential (g*Z_surf)
    type(fv_grid_type), intent(IN) :: gridstruct
    type(domain2d), intent(INOUT) :: domain
!
    real(FVPRC)  ::  u2f(bd%isd:bd%ied,bd%jsd:bd%jed,kmax)
    real(FVPRC), parameter:: u0   = 60.   ! scaling velocity
    real(FVPRC), parameter:: sday = 86400.
    real(FVPRC) rcv, tau0
    integer i, j, k

    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

    real(FVPRC) :: rf(npz)

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

    rcv = 1. / (cp - rg)

     if ( .not. RF_initialized ) then
!#ifdef HIWPP
!          allocate ( u00(is:ie,  js:je+1,npz) )
!          allocate ( v00(is:ie+1,js:je  ,npz) )
!!$OMP parallel do default(none) shared(is,ie,js,je,npz,u00,u,v00,v)
!          do k=1,npz
!             do j=js,je+1
!                do i=is,ie
!                   u00(i,j,k) = u(i,j,k)
!                enddo
!             enddo
!             do j=js,je
!                do i=is,ie+1
!                   v00(i,j,k) = v(i,j,k)
!                enddo
!             enddo
!          enddo
!#endif
!#ifdef SMALL_EARTH
!          tau0 = tau
!#else
          tau0 = tau * sday
!#endif
          !allocate( rf(npz) )
          rf(:) = 0.

          !if( is_master() ) write(6,*) 'Rayleigh friction E-folding time (days):'
          do k=1, npz
             if ( pm(k) < rf_cutoff ) then
                  rf(k) = dt/tau0*sin(0.5*pi*log(rf_cutoff/pm(k))/log(rf_cutoff/ptop))**2
                  !if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
                  kmax = k
             else
                  exit
             endif
          enddo
          RF_initialized = .true.
     endif

    call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd)

    !allocate( u2f(isd:ied,jsd:jed,kmax) )

!$OMP parallel do default(none) shared(is,ie,js,je,kmax,pm,rf_cutoff,hydrostatic,ua,va,agrid, &
!$OMP                                  u2f,rf,w)
    do k=1,kmax
       if ( pm(k) < rf_cutoff ) then
       do j=js,je
        if ( hydrostatic ) then
          do i=is,ie
             if ( abs(ua(i,j,k)) > 35.*cos(agrid(i,j,2)) )  then
                  u2f(i,j,k) = 1./(1.+rf(k)*sqrt(ua(i,j,k)**2+va(i,j,k)**2)/u0)
             else
                  u2f(i,j,k) = 1.
             endif
          enddo
        else
          do i=is,ie
             if ( abs(ua(i,j,k)) > 35.*cos(agrid(i,j,2)) .or. abs(w(i,j,k))>7.5 )  then
                  u2f(i,j,k) = 1./(1.+rf(k)*sqrt(ua(i,j,k)**2+va(i,j,k)**2+w(i,j,k)**2)/u0)
             else
                  u2f(i,j,k) = 1.
             endif
          enddo
        endif
       enddo
       endif ! p check
    enddo

                                        !if (fv_timing_onoff) call timing_on('COMM_TOTAL')
    call mpp_update_domains(u2f, domain)
                                        !if (fv_timing_onoff) call timing_off('COMM_TOTAL')


!$OMP parallel do default(none) shared(is,ie,js,je,kmax,pm,rf_cutoff,w,rf,u,v,u00,v00, &
!$OMP                                  conserve,hydrostatic,pt,ua,va,u2f,cp,rg,ptop,rcv)
     do k=1,kmax
        if ( pm(k) < rf_cutoff ) then
!#ifdef HIWPP
!             do j=js,je
!                do i=is,ie
!                   w(i,j,k) = w(i,j,k)/(1.+rf(k))
!                enddo
!             enddo
!             do j=js,je+1
!                do i=is,ie
!                   u(i,j,k) = (u(i,j,k)+rf(k)*u00(i,j,k))/(1.+rf(k))
!                enddo
!             enddo
!             do j=js,je
!                do i=is,ie+1
!                   v(i,j,k) = (v(i,j,k)+rf(k)*v00(i,j,k))/(1.+rf(k))
!                enddo
!             enddo
!#else
! Add heat so as to conserve TE
          if ( conserve ) then
             if ( hydrostatic ) then
               do j=js,je
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k) + 0.5*(ua(i,j,k)**2+va(i,j,k)**2)*(1.-u2f(i,j,k)**2)/(cp-rg*ptop/pm(k))
                  enddo
               enddo
             else
               do j=js,je
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k) + 0.5*(ua(i,j,k)**2+va(i,j,k)**2+w(i,j,k)**2)*(1.-u2f(i,j,k)**2)*rcv
                  enddo
               enddo
             endif
          endif
             do j=js,je+1
                do i=is,ie
                   u(i,j,k) = 0.5*(u2f(i,j-1,k)+u2f(i,j,k))*u(i,j,k)
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   v(i,j,k) = 0.5*(u2f(i-1,j,k)+u2f(i,j,k))*v(i,j,k)
                enddo
             enddo
          if ( .not. hydrostatic ) then
             do j=js,je
                do i=is,ie
                   w(i,j,k) = u2f(i,j,k)*w(i,j,k)
                enddo
             enddo
          endif
!#endif
        endif
     enddo

     !deallocate ( u2f )

 end subroutine Rayleigh_Super


 subroutine Rayleigh_Friction(dt, npx, npy, npz, ks, pm, tau, u, v, w, pt,  &
                              ua, va, delz, cp, rg, ptop, hydrostatic, conserve, &
                              rf_cutoff, gridstruct, domain, bd)
    real(FVPRC), intent(in):: dt
    real(FVPRC), intent(in):: tau              ! time scale (days)
    real(FVPRC), intent(in):: cp, rg, ptop, rf_cutoff
    real(FVPRC), intent(in),  dimension(npz):: pm
    integer, intent(in):: npx, npy, npz, ks
    logical, intent(in):: hydrostatic
    logical, intent(in):: conserve
    type(fv_grid_bounds_type), intent(IN) :: bd
    real(FVPRC), intent(inout):: u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) ! D grid zonal wind (m/s)
    real(FVPRC), intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz) ! D grid meridional wind (m/s)
    real(FVPRC), intent(inout):: w(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! cell center vertical wind (m/s)
    real(FVPRC), intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! temp
    real(FVPRC), intent(inout):: ua(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real(FVPRC), intent(inout):: va(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real(FVPRC), intent(inout):: delz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! delta-height (m); non-hydrostatic only
    type(fv_grid_type), intent(IN) :: gridstruct
    type(domain2d), intent(INOUT) :: domain
! local:
    real(FVPRC)  ::  u2f(bd%isd:bd%ied,bd%jsd:bd%jed,kmax)
    real(FVPRC), parameter:: sday = 86400.
    real(FVPRC), parameter:: u000 = 4900.   ! scaling velocity  **2
    real(FVPRC)  rcv
    integer i, j, k

    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

    real(FVPRC) :: rf(npz)
    
    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed


    rcv = 1. / (cp - rg)

    if ( .not. RF_initialized ) then
          !allocate( rf(npz) )
          !if( is_master() ) write(6,*) 'Rayleigh friction E-folding time (days):'
          do k=1, npz
             if ( pm(k) < rf_cutoff ) then
                  rf(k) = dt/(tau*sday)*sin(0.5*pi*log(rf_cutoff/pm(k))/log(rf_cutoff/ptop))**2
                  !if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
                  kmax = k
             else
                  exit
             endif
          enddo
          RF_initialized = .true.
    endif

    !allocate( u2f(isd:ied,jsd:jed,kmax) )

    call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd)
    u2f = 0.
!$OMP parallel do default(none) shared(is,ie,js,je,kmax,u2f,hydrostatic,ua,va,w)
    do k=1,kmax
        if ( hydrostatic ) then
           do j=js,je
              do i=is,ie
                 u2f(i,j,k) = ua(i,j,k)**2 + va(i,j,k)**2
              enddo
           enddo
        else
           do j=js,je
              do i=is,ie
                 u2f(i,j,k) = ua(i,j,k)**2 + va(i,j,k)**2 + w(i,j,k)**2
              enddo
           enddo
        endif
    enddo
                                        !if (fv_timing_onoff) call timing_on('COMM_TOTAL')
    call mpp_update_domains(u2f, domain)
                                        !if (fv_timing_onoff) call timing_off('COMM_TOTAL')

!$OMP parallel do default(none) shared(is,ie,js,je,kmax,conserve,hydrostatic,pt,u2f,cp,rg, &
!$OMP                                  ptop,pm,rf,delz,rcv,u,v,w)
     do k=1,kmax

        if ( conserve ) then
           if ( hydrostatic ) then
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = pt(i,j,k) + 0.5*u2f(i,j,k)/(cp-rg*ptop/pm(k))      &
                             * ( 1. - 1./(1.+rf(k)*sqrt(u2f(i,j,k)/u000))**2 )
                enddo
             enddo
           else
             do j=js,je
                do i=is,ie
                   delz(i,j,k) = delz(i,j,k) / pt(i,j,k)
                   pt(i,j,k) = pt(i,j,k) + 0.5*u2f(i,j,k) * rcv      &
                             * ( 1. - 1./(1.+rf(k)*sqrt(u2f(i,j,k)/u000))**2 )
                   delz(i,j,k) = delz(i,j,k) * pt(i,j,k)
                enddo
             enddo
           endif
        endif

        do j=js-1,je+1
           do i=is-1,ie+1
              u2f(i,j,k) = rf(k)*sqrt(u2f(i,j,k)/u000)
           enddo
        enddo

        do j=js,je+1
           do i=is,ie
              u(i,j,k) = u(i,j,k) / (1.+0.5*(u2f(i,j-1,k)+u2f(i,j,k)))
           enddo
        enddo
        do j=js,je
           do i=is,ie+1
              v(i,j,k) = v(i,j,k) / (1.+0.5*(u2f(i-1,j,k)+u2f(i,j,k)))
           enddo
        enddo

        if ( .not. hydrostatic ) then
              do j=js,je
                 do i=is,ie
                    w(i,j,k) = w(i,j,k) / (1.+u2f(i,j,k))
                 enddo
              enddo
        endif

     enddo

     !deallocate ( u2f )

 end subroutine Rayleigh_Friction

 subroutine compute_aam(npz, is, ie, js, je, isd, ied, jsd, jed, gridstruct, bd,   &
                        ptop, ua, va, u, v, delp, aam, ps, m_fac)
! Compute vertically (mass) integrated Atmospheric Angular Momentum
    integer, intent(in):: npz
    integer, intent(in):: is,  ie,  js,  je
    integer, intent(in):: isd, ied, jsd, jed
    real(FVPRC), intent(in):: ptop
    real(FVPRC), intent(inout):: u(isd:ied  ,jsd:jed+1,npz) ! D grid zonal wind (m/s)
    real(FVPRC), intent(inout):: v(isd:ied+1,jsd:jed,npz) ! D grid meridional wind (m/s)
    real(FVPRC), intent(inout):: delp(isd:ied,jsd:jed,npz)
    real(FVPRC), intent(inout), dimension(isd:ied,jsd:jed, npz):: ua, va
    real(FVPRC), intent(out):: aam(is:ie,js:je)
    real(FVPRC), intent(out):: m_fac(is:ie,js:je)
    real(FVPRC), intent(out):: ps(isd:ied,jsd:jed)
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(fv_grid_type), intent(IN) :: gridstruct
! local:
    real(FVPRC), dimension(is:ie):: r1, r2, dm
    integer i, j, k

  call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd)
    
!$OMP parallel do default(none) shared(is,ie,js,je,npz,gridstruct,aam,m_fac,ps,ptop,delp,agrav,ua) &
!$OMP                          private(r1, r2, dm)
  do j=js,je
     do i=is,ie
        r1(i) = radius*cos(gridstruct%agrid(i,j,2))
        r2(i) = r1(i)*r1(i)
        aam(i,j) = 0.
        m_fac(i,j) = 0.
        ps(i,j) = ptop
     enddo
     do k=1,npz
        do i=is,ie
           dm(i) = delp(i,j,k)
           ps(i,j) = ps(i,j) + dm(i)
           dm(i) = dm(i)*agrav
           aam(i,j) = aam(i,j) + (r2(i)*omega + r1(i)*ua(i,j,k)) * dm(i)
           m_fac(i,j) = m_fac(i,j) + dm(i)*r2(i)
        enddo
     enddo
  enddo

 end subroutine compute_aam

  subroutine geos_to_fv3(bd, npz, kappa, ptop, delp, pe, pk, pkz, peln, pt)

    implicit none

    !Arguments
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: npz
    real(FVPRC), intent(in) :: kappa, ptop

    real(FVPRC), intent(inout) :: pt  (bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real(FVPRC), intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real(FVPRC), intent(inout) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)
    real(FVPRC), intent(inout) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)
    real(FVPRC), intent(inout) :: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)
    real(FVPRC), intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)

    !Locals
    integer :: i, j, k, is, ie, js, je

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je

    pe(:,:,:) = tiny_number
    pe(:,1,:) = ptop
    do k=2,npz+1
      do j=js,je
        do i=is,ie
          pe(i,k,j)   = pe(i, k-1, j) + delp(i,j,k-1)
        enddo
      enddo
    enddo

    do k=1,npz+1
      do j=js,je
        do i=is,ie
          peln(i,k,j) = log(pe(i,k,j))
        enddo
      enddo
    enddo

    do k=1,npz+1
      do j=js,je
        do i=is,ie
          pk(i,j,k)   = exp( kappa*peln(i,k,j) )
        enddo
      enddo
    enddo

    do k=1,npz
      do j=js,je
        do i=is,ie
          pkz(i,j,k) = (pk(i,j,k+1) - pk(i,j,k) ) / (kappa*(peln(i,k+1,j) - peln(i,k,j)) )
        end do
      end do
    end do

    pt(is:ie,js:je,:) = pt(is:ie,js:je,:)*pkz(is:ie,js:je,:)

  end subroutine geos_to_fv3

  subroutine fv3_to_geos(bd, npz, pkz, pt)

    implicit none

    !Arguments
    type(fv_grid_bounds_type), intent(in) :: bd
    integer, intent(in) :: npz
    real(fvprc), intent(inout) :: pt  (bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real(fvprc), intent(inout) :: pkz (bd%is :bd%ie ,bd%js :bd%je ,npz)

    !Locals
    integer :: i, j, k, is, ie, js, je

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je

    do k=1,npz
      do j=js,je
        do i=is,ie
           pt(i,j,k) = pt(i,j,k)/pkz(i,j,k)
        end do
      end do
    end do

  end subroutine fv3_to_geos

end module fv_dynamics_mod
