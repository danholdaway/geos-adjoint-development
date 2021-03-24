module d2a_mod

implicit none

integer, parameter :: kind_real = 8

integer, parameter :: DGRID_NE = 1

private
public d2a

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
  real(kind=kind_real), allocatable :: a11(:,:)
  real(kind=kind_real), allocatable :: a12(:,:)
  real(kind=kind_real), allocatable :: a21(:,:)
  real(kind=kind_real), allocatable :: a22(:,:)

end type fv3jedi_geom

contains

 subroutine d2a(geom, u, v, ua, va)

 !c2l_ord4

 !use mpp_domains_mod, only: mpp_update_domains, DGRID_NE

 implicit none

 type(fv3jedi_geom), intent(inout) :: geom

  real(kind=kind_real), intent(inout):: u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,geom%npz)
  real(kind=kind_real), intent(inout):: v(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,geom%npz)
  real(kind=kind_real), intent(inout)::  ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,geom%npz)
  real(kind=kind_real), intent(inout)::  va(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,geom%npz)

! Local 
! 4-pt Lagrange interpolation
  real(kind=kind_real) :: c1 =  1.125
  real(kind=kind_real) :: c2 = -0.125
  real(kind=kind_real) :: utmp(geom%isc:geom%iec,  geom%jsc:geom%jec+1)
  real(kind=kind_real) :: vtmp(geom%isc:geom%iec+1,geom%jsc:geom%jec)
  real(kind=kind_real) :: wu(geom%isc:geom%iec,  geom%jsc:geom%jec+1)
  real(kind=kind_real) :: wv(geom%isc:geom%iec+1,geom%jsc:geom%jec)

  integer i, j, k
  integer :: is,  ie,  js,  je, npx, npy, npz

  is  = geom%isc
  ie  = geom%iec
  js  = geom%jsc
  je  = geom%jec
  npx = geom%npx
  npy = geom%npy
  npz = geom%npz

       call mpp_update_domains(u, v, geom%domain, gridtype=DGRID_NE)

!$OMP parallel do default(none) shared(is,ie,js,je,npz,npx,npy,c2,c1, &
!$OMP                                  u,v,ua,va)         &
!$OMP                          private(utmp, vtmp, wu, wv)
 do k=1,npz

   do j=max(2,js),min(npy-2,je)
     do i=max(2,is),min(npx-2,ie)
       utmp(i,j) = c2*(u(i,j-1,k)+u(i,j+2,k)) + c1*(u(i,j,k)+u(i,j+1,k))
       vtmp(i,j) = c2*(v(i-1,j,k)+v(i+2,j,k)) + c1*(v(i,j,k)+v(i+1,j,k))
     enddo
   enddo

   if ( js==1  ) then
     do i=is,ie+1
       wv(i,1) = v(i,1,k)*geom%dy(i,1)
     enddo
     do i=is,ie
       vtmp(i,1) = 2.*(wv(i,1) + wv(i+1,1)) / (geom%dy(i,1)+geom%dy(i+1,1))
       utmp(i,1) = 2.*(u(i,1,k)*geom%dx(i,1) + u(i,2,k)*geom%dx(i,2))   &
                    / (         geom%dx(i,1) +          geom%dx(i,2))
     enddo
   endif

   if ( (je+1)==npy   ) then
     j = npy-1
     do i=is,ie+1
       wv(i,j) = v(i,j,k)*geom%dy(i,j)
     enddo
     do i=is,ie
       vtmp(i,j) = 2.*(wv(i,j) + wv(i+1,j)) / (geom%dy(i,j)+geom%dy(i+1,j))
       utmp(i,j) = 2.*(u(i,j,k)*geom%dx(i,j) + u(i,j+1,k)*geom%dx(i,j+1))   &
                   / (         geom%dx(i,j) +            geom%dx(i,j+1))
     enddo
   endif

   if ( is==1 ) then
     i = 1
     do j=js,je
       wv(1,j) = v(1,j,k)*geom%dy(1,j)
       wv(2,j) = v(2,j,k)*geom%dy(2,j)
     enddo
     do j=js,je+1
       wu(i,j) = u(i,j,k)*geom%dx(i,j)
     enddo
     do j=js,je
       utmp(i,j) = 2.*(wu(i,j) + wu(i,j+1))/(geom%dx(i,j)+geom%dx(i,j+1))
       vtmp(i,j) = 2.*(wv(1,j) + wv(2,j  ))/(geom%dy(1,j)+geom%dy(2,j))
     enddo
   endif

   if ( (ie+1)==npx) then
     i = npx-1
     do j=js,je
       wv(i,  j) = v(i,  j,k)*geom%dy(i,  j)
       wv(i+1,j) = v(i+1,j,k)*geom%dy(i+1,j)
     enddo
     do j=js,je+1
       wu(i,j) = u(i,j,k)*geom%dx(i,j)
     enddo
     do j=js,je
       utmp(i,j) = 2.*(wu(i,j) + wu(i,  j+1))/(geom%dx(i,j)+geom%dx(i,j+1))
       vtmp(i,j) = 2.*(wv(i,j) + wv(i+1,j  ))/(geom%dy(i,j)+geom%dy(i+1,j))
     enddo
   endif

   !Transform local a-grid winds into latitude-longitude coordinates
   do j=js,je
     do i=is,ie
       ua(i,j,k) = geom%a11(i,j)*utmp(i,j) + geom%a12(i,j)*vtmp(i,j)
       va(i,j,k) = geom%a21(i,j)*utmp(i,j) + geom%a22(i,j)*vtmp(i,j)
     enddo
   enddo

 enddo

end subroutine d2a

!-----------------

subroutine mpp_update_domains(fieldu, fieldv, domain, gridtype)

real(kind=kind_real), intent(inout) :: fieldu(:,:,:)
real(kind=kind_real), intent(inout) :: fieldv(:,:,:)
type(domain2D), intent(in) :: domain
integer, intent(in) :: gridtype

if (gridtype == 1) then
fieldu = fieldu * fieldu * domain%dummy
fieldv = fieldv * fieldv * domain%dummy
endif

end subroutine mpp_update_domains

!-----------------

end module d2a_mod
