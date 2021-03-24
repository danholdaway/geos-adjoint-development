
!> Fortran module handling geometry for the FV3 model

module fv3jedi_geom_mod

use mpp_domains_mod, only: domain2d, kind_real

implicit none
private
public :: fv3jedi_geom

! ------------------------------------------------------------------------------

!> Skinny version of fv_grid_bounds_type
type fv_grid_bounds_type
    integer :: isd, ied, jsd, jed ! data domain
    integer :: isc, iec, jsc, jec ! compute domain
end type fv_grid_bounds_type

!> Fortran derived type to hold geometry data for the FV3JEDI model
type :: fv3jedi_geom
  !From user, maybe via input.nml
  integer :: npx                                          !x-dir grid edge points per tile
  integer :: npy                                          !y-dir grid edge points per tile
  integer :: npz                                          !z-dir grid points global
  logical :: hydrostatic                                  !Are fields on this geometry hydrostatic
  integer :: layout(2)                                    !Processor layout for computation
  integer :: io_layout(2)                                 !Processor layout for read/write
  integer :: halo                                         !Number of halo points, normally 3
  character(len=255) :: nml_file                          !FV3 nml file associated with this geom
  character(len=255) :: trc_file                          !FV3 field_table associated with this geom
  character(len=255) :: wind_type                         !A-grid or D-grid in the state vector
  !Hardwired or determined
  logical :: am_i_root_pe = .false.                       !Is this the root process 
  integer :: size_cubic_grid                              !Size of cubed sphere grid (cell center)
  integer :: ntracers                                     !Number of tracers
  type(domain2D) :: domain                                !MPP domain
  type(fv_grid_bounds_type) :: bd                         !FV grid bounds
  integer :: ntile                                        !Tile ID
  integer :: ntiles = 6                                   !Number of tiles, always 6
  integer :: stackmax                                     !Stackmax
  real(kind=kind_real), allocatable :: grid_lon(:,:)      !Longitude at cell center
  real(kind=kind_real), allocatable :: grid_lat(:,:)      !Latitude at cell center
  real(kind=kind_real), allocatable :: egrid_lon(:,:)      !Longitude at cell center
  real(kind=kind_real), allocatable :: egrid_lat(:,:)      !Latitude at cell center
  real(kind=kind_real), allocatable :: area(:,:)          !Grid area
  real(kind=kind_real), allocatable :: ak(:),bk(:)        !Model level coefficients
  real(kind=kind_real) :: ptop                            !Pressure at top of domain
  real(kind=kind_real), allocatable :: sin_sg(:,:,:)
  real(kind=kind_real), allocatable :: cos_sg(:,:,:)
  real(kind=kind_real), allocatable :: cosa_u(:,:)
  real(kind=kind_real), allocatable :: cosa_v(:,:)
  real(kind=kind_real), allocatable :: cosa_s(:,:)
  real(kind=kind_real), allocatable :: rsin_u(:,:)
  real(kind=kind_real), allocatable :: rsin_v(:,:)
  real(kind=kind_real), allocatable :: rsin2(:,:)
  real(kind=kind_real), allocatable :: dxa(:,:)
  real(kind=kind_real), allocatable :: dya(:,:)
  real(kind=kind_real), allocatable :: dx(:,:)
  real(kind=kind_real), allocatable :: dy(:,:)
  real(kind=kind_real), allocatable :: dxc(:,:)
  real(kind=kind_real), allocatable :: dyc(:,:)
  real(kind=kind_real), allocatable :: rarea(:,:)
  real(kind=kind_real), allocatable :: rarea_c(:,:)
  real(kind=kind_real), allocatable :: edge_w(:)
  real(kind=kind_real), allocatable :: edge_e(:)
  real(kind=kind_real), allocatable :: edge_s(:)
  real(kind=kind_real), allocatable :: edge_n(:)
  real(kind=kind_real), allocatable :: grid(:,:,:)
  real(kind=kind_real), allocatable :: agrid(:,:,:)

  logical :: sw_corner, se_corner, ne_corner, nw_corner
end type fv3jedi_geom

end module fv3jedi_geom_mod
