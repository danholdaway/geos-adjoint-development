module fv_arrays_mod
#include <fms_platform.h>
 use mpp_domains_mod,  only: domain2d
public

! integer, parameter :: REAL4  = selected_real_kind(P= 6,R=37)
! integer, parameter :: REAL8  = selected_real_kind(P=13,R=300)
  integer, parameter :: REAL4  = kind(1.00)
  integer, parameter :: REAL8  = kind(1.d0)

#ifdef SINGLE_FV
  integer, public, parameter ::  m_precision = REAL4  ! Full Build Precision for the model
  integer, public, parameter ::  i_precision = REAL4  ! Precision of the Vector Interpolations 
  integer, public, parameter ::  p_precision = REAL8  ! Precision of the Pressure/Remapping terms
  integer, public, parameter :: pg_precision = REAL8  ! Precision of the Pressure Gradient terms
  integer, public, parameter :: ke_precision = REAL4  ! Precision of the D-Grid Kinetic Energy Calculation
  integer, public, parameter ::  g_precision = REAL8  ! Precision of the Grid Arrays
#else
  integer, public, parameter ::  m_precision = REAL8  ! Full Build Precision for the model
  integer, public, parameter ::  i_precision = REAL8  ! Precision of the Vector Interpolations 
  integer, public, parameter ::  p_precision = REAL8  ! Precision of the Pressure/Remapping terms
  integer, public, parameter :: pg_precision = REAL8  ! Precision of the Pressure Gradient terms
  integer, public, parameter :: ke_precision = REAL8  ! Precision of the D-Grid Kinetic Energy Calculation
  integer, public, parameter ::  g_precision = REAL8  ! Precision of the Grid Arrays
#endif

  type fv_atmos_type
     type(domain2d), pointer :: domain =>NULL()
!-----------------------------------------------------------------------
! Five prognostic state variables for the f-v dynamics
!-----------------------------------------------------------------------
! dyn_state:
! D-grid prognostatic variables: u, v, and delp (and other scalars)
!
!     o--------u(i,j+1)----------o
!     |           |              |
!     |           |              |
!  v(i,j)------scalar(i,j)----v(i+1,j)
!     |           |              |
!     |           |              |
!     o--------u(i,j)------------o
!
! The C grid component is "diagnostic" in that it is predicted every time step
! from the D grid variables.
    real, _ALLOCATABLE :: u(:,:,:)    _NULL  ! D grid zonal wind (m/s)
    real, _ALLOCATABLE :: v(:,:,:)    _NULL  ! D grid meridional wind (m/s)
    real, _ALLOCATABLE :: um(:,:,:)   _NULL  ! D grid zonal wind (m/s) at n-1
    real, _ALLOCATABLE :: vm(:,:,:)   _NULL  ! D .... meridional.............
    real, _ALLOCATABLE :: pt(:,:,:)   _NULL  ! temperature (K)
    real, _ALLOCATABLE :: delp(:,:,:) _NULL  ! pressure thickness (pascal)
    real, _ALLOCATABLE :: q(:,:,:,:)  _NULL  ! specific humidity and constituents

!----------------------
! non-hydrostatic state:
!----------------------------------------------------------------------
    real, _ALLOCATABLE ::     w(:,:,:)  _NULL  ! cell center vertical wind (m/s)
    real, _ALLOCATABLE ::  delz(:,:,:)  _NULL  ! layer thickness (meters)
    real, _ALLOCATABLE ::   ze0(:,:,:)  _NULL  ! height at layer edges for remapping

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real(p_precision), _ALLOCATABLE :: ps (:,:)      _NULL  ! Surface pressure (pascal)
    real(p_precision), _ALLOCATABLE :: pe (:,:,: )   _NULL  ! edge pressure (pascal)
    real(p_precision), _ALLOCATABLE :: pk  (:,:,:)   _NULL  ! pe**cappa
    real(p_precision), _ALLOCATABLE :: peln(:,:,:)   _NULL  ! ln(pe)
    real(p_precision), _ALLOCATABLE :: pkz (:,:,:)   _NULL  ! finite-volume mean pk

! For phys coupling:
    real, _ALLOCATABLE :: u_srf(:,:)    _NULL  ! Surface u-wind
    real, _ALLOCATABLE :: v_srf(:,:)    _NULL  ! Surface v-wind
    real, _ALLOCATABLE :: sgh(:,:)      _NULL  ! Terrain standard deviation
    real, _ALLOCATABLE :: oro(:,:)      _NULL  ! land fraction (1: all land; 0: all water)
    real, _ALLOCATABLE ::  ts(:,:)      _NULL  ! skin (sst) temperature from NCEP/GFS (K) -- tile
 
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, _ALLOCATABLE :: phis(:,:)     _NULL  ! Surface geopotential (g*Z_surf)
    real, _ALLOCATABLE :: omga(:,:,:)   _NULL  ! Vertical pressure velocity (pa/s)
    real, _ALLOCATABLE :: ua(:,:,:)     _NULL  ! (ua, va) are mostly used as the A grid winds
    real, _ALLOCATABLE :: va(:,:,:)     _NULL
    real, _ALLOCATABLE :: uc(:,:,:)     _NULL  ! (uc, vc) are mostly used as the C grid winds
    real, _ALLOCATABLE :: vc(:,:,:)     _NULL

    real, _ALLOCATABLE :: ak(:)  _NULL
    real, _ALLOCATABLE :: bk(:)  _NULL

! Accumulated Mass flux arrays
    real, _ALLOCATABLE ::  mfx(:,:,:)  _NULL
    real, _ALLOCATABLE ::  mfy(:,:,:)  _NULL
! Accumulated Courant number arrays
    real, _ALLOCATABLE ::  cx(:,:,:)  _NULL
    real, _ALLOCATABLE ::  cy(:,:,:)  _NULL

! Horizontal Grid descriptors
    real(g_precision), pointer :: grid(:,:,:)  _NULL  ! Leave as a pointer for now
    real(g_precision), pointer :: agrid(:,:,:)  _NULL  ! Leave as a pointer for now
    real(g_precision), pointer :: grid_g(:,:,:) _NULL  ! "global" grid (one face of a cube)

    real   :: consv_te

    integer :: isc, iec, jsc, jec
    integer :: isd, ied, jsd, jed
    integer :: ks, npx, npy, npz, npz_rst, ng, ntiles
    integer :: n_sponge    ! Number of sponge layers at the top of the atmosphere
    integer :: k_top       ! Starting layer for non-hydrostatic dynamics
    integer :: ncnst, pnats, ndims, k_split, n_split, m_split, q_split, print_freq
    integer :: nwat        ! water substance
    integer :: fv_sg_adj
    integer :: i_sst, j_sst

! Namelist control values
    logical :: fill
    logical :: check_surface_pressure
    logical :: range_warn
    logical :: z_tracer
    logical :: do_Held_Suarez
    logical :: reproduce_sum
    logical :: moist_phys
    logical :: srf_init
    logical :: mountain
    logical :: non_ortho
    logical :: adjust_dry_mass
    logical :: shallow_water, hydrostatic, phys_hydrostatic
    logical :: hybrid_z, Make_NH, make_hybrid_z
    logical :: external_ic
    logical :: ncep_ic
    logical :: fv_diag_ic
    logical :: fv_land
    logical :: init_wind_m, shift_west, no_cgrid
    logical :: nudge
    logical :: tq_filter
    logical :: warm_start

    character(len=128) :: res_latlon_dynamics  ! restart file from the latlon FV core
    character(len=128) :: res_latlon_tracers   ! tracer restart file from the latlon core

    real    :: dry_mass

  end type fv_atmos_type
end module fv_arrays_mod
