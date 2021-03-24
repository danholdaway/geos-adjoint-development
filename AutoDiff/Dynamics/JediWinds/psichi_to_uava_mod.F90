module psichi_to_uava_mod

implicit none

integer, parameter :: kind_real = 8

private
public psichi_to_uava

type domain2D
  integer :: dummy
end type

type :: fv3jedi_geom

  integer :: isd, ied, jsd, jed                           !data domain
  integer :: isc, iec, jsc, jec                           !compute domain
  integer :: npx,npy,npz                                  !x/y/z-dir grid edge points per tile
  integer :: layout(2)                                    !Processor layout for computation
  integer :: io_layout(2)                                 !Processor layout for read/write
  integer :: halo                                         !Number of halo points, normally 3
  character(len=255) :: nml_file                          !FV3 nml file associated with this geom
  integer :: size_cubic_grid                              !Size of cubed sphere grid (cell center)
  type(domain2D) :: domain                                !MPP domain
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
  
  real(kind=kind_real), allocatable :: vlon(:,:,:)
  real(kind=kind_real), allocatable :: vlat(:,:,:)
  real(kind=kind_real), allocatable :: edge_vect_n(:)
  real(kind=kind_real), allocatable :: edge_vect_e(:)
  real(kind=kind_real), allocatable :: edge_vect_s(:)
  real(kind=kind_real), allocatable :: edge_vect_w(:)
  real(kind=kind_real), allocatable :: es(:,:,:,:)
  real(kind=kind_real), allocatable :: ew(:,:,:,:)

end type fv3jedi_geom

contains

subroutine psichi_to_uava(geom,psi,chi,ua,va)

 implicit none
 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(inout) :: psi(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz) !Stream function
 real(kind=kind_real), intent(inout) :: chi(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(out)   ::  ua(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Agrid winds (u)
 real(kind=kind_real), intent(out)   ::  va(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Agrid winds (v)

 integer :: i,j,k

 call mpp_update_domains(psi, geom%domain, complete=.true.)
 call mpp_update_domains(chi, geom%domain, complete=.true.)
 
 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec

        ua(i,j,k) =  (psi(i,j+1,k) - psi(i,j-1,k))/(geom%dyc(i,j) + geom%dyc(i,j+1)) + &
                     (chi(i+1,j,k) - chi(i-1,j,k))/(geom%dxc(i,j) + geom%dxc(i+1,j))
        va(i,j,k) = -(psi(i+1,j,k) - psi(i-1,j,k))/(geom%dxc(i,j) + geom%dxc(i+1,j)) + &
                     (chi(i,j+1,k) - chi(i,j-1,k))/(geom%dyc(i,j) + geom%dyc(i,j+1))

     enddo
   enddo
 enddo 

endsubroutine psichi_to_uava

!-----------------

subroutine mpp_update_domains(field, domain, complete)

real(kind=kind_real), intent(inout) :: field(:,:,:)
type(domain2D), intent(in) :: domain
logical, intent(in) :: complete

if (complete) then
field = field * field * domain%dummy
endif

end subroutine mpp_update_domains

!-----------------

end module psichi_to_uava_mod
