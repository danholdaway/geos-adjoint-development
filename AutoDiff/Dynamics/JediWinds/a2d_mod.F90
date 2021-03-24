module a2d_mod

implicit none

integer, parameter :: kind_real = 8

private
public a2d

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
subroutine a2d(geom, ua, va, ud, vd)

!use mpp_domains_mod, only: mpp_update_domains

 implicit none

 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(in)    :: ua(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,geom%npz)
 real(kind=kind_real), intent(in)    :: va(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,geom%npz)
 real(kind=kind_real), intent(inout) :: ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,geom%npz)
 real(kind=kind_real), intent(inout) :: vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,geom%npz)

 integer :: is ,ie , js ,je 
 integer :: npx, npy, npz
 integer :: i,j,k, im2,jm2

 real(kind=kind_real) :: uatemp(geom%isd:geom%ied,geom%jsd:geom%jed,geom%npz)
 real(kind=kind_real) :: vatemp(geom%isd:geom%ied,geom%jsd:geom%jed,geom%npz)

 real(kind=kind_real) :: v3(geom%isc-1:geom%iec+1,geom%jsc-1:geom%jec+1,3)
 real(kind=kind_real) :: ue(geom%isc-1:geom%iec+1,geom%jsc  :geom%jec+1,3)    ! 3D winds at edges
 real(kind=kind_real) :: ve(geom%isc  :geom%iec+1,geom%jsc-1:geom%jec+1,3)    ! 3D winds at edges
 real(kind=kind_real), dimension(geom%isc:geom%iec):: ut1, ut2, ut3
 real(kind=kind_real), dimension(geom%jsc:geom%jec):: vt1, vt2, vt3

 npx = geom%npx
 npy = geom%npy
 npz = geom%npz
 is  = geom%isc
 ie  = geom%iec
 js  = geom%jsc
 je  = geom%jec

 im2 = (npx-1)/2
 jm2 = (npy-1)/2

 uatemp(:,:,:) = 0.0
 vatemp(:,:,:) = 0.0

 uatemp(is:ie,js:je,:) = ua
 vatemp(is:ie,js:je,:) = va

 call mpp_update_domains(uatemp, geom%domain, complete=.true.)
 call mpp_update_domains(vatemp, geom%domain, complete=.true.)

 do k=1, npz

   do j=js-1,je+1
     do i=is-1,ie+1
       v3(i,j,1) = uatemp(i,j,k)*geom%vlon(i,j,1) + vatemp(i,j,k)*geom%vlat(i,j,1)
       v3(i,j,2) = uatemp(i,j,k)*geom%vlon(i,j,2) + vatemp(i,j,k)*geom%vlat(i,j,2)
       v3(i,j,3) = uatemp(i,j,k)*geom%vlon(i,j,3) + vatemp(i,j,k)*geom%vlat(i,j,3)
     enddo
   enddo

   do j=js,je+1
     do i=is-1,ie+1
       ue(i,j,1) = v3(i,j-1,1) + v3(i,j,1)
       ue(i,j,2) = v3(i,j-1,2) + v3(i,j,2)
       ue(i,j,3) = v3(i,j-1,3) + v3(i,j,3)
     enddo
   enddo

   do j=js-1,je+1
     do i=is,ie+1
       ve(i,j,1) = v3(i-1,j,1) + v3(i,j,1)
       ve(i,j,2) = v3(i-1,j,2) + v3(i,j,2)
       ve(i,j,3) = v3(i-1,j,3) + v3(i,j,3)
     enddo
   enddo

   if ( is==1 ) then
     !i = 1
     do j=js,je
       if ( j>jm2 ) then
         vt1(j) = geom%edge_vect_w(j)*ve(1,j-1,1)+(1.-geom%edge_vect_w(j))*ve(1,j,1)
         vt2(j) = geom%edge_vect_w(j)*ve(1,j-1,2)+(1.-geom%edge_vect_w(j))*ve(1,j,2)
         vt3(j) = geom%edge_vect_w(j)*ve(1,j-1,3)+(1.-geom%edge_vect_w(j))*ve(1,j,3)
       else
         vt1(j) = geom%edge_vect_w(j)*ve(1,j+1,1)+(1.-geom%edge_vect_w(j))*ve(1,j,1)
         vt2(j) = geom%edge_vect_w(j)*ve(1,j+1,2)+(1.-geom%edge_vect_w(j))*ve(1,j,2)
         vt3(j) = geom%edge_vect_w(j)*ve(1,j+1,3)+(1.-geom%edge_vect_w(j))*ve(1,j,3)
       endif
     enddo
     do j=js,je
       ve(1,j,1) = vt1(j)
       ve(1,j,2) = vt2(j)
       ve(1,j,3) = vt3(j)
     enddo
   endif

   if ( (ie+1)==npx ) then
     !i = npx
     do j=js,je
       if ( j>jm2 ) then
         vt1(j) = geom%edge_vect_e(j)*ve(npx,j-1,1)+(1.-geom%edge_vect_e(j))*ve(npx,j,1)
         vt2(j) = geom%edge_vect_e(j)*ve(npx,j-1,2)+(1.-geom%edge_vect_e(j))*ve(npx,j,2)
         vt3(j) = geom%edge_vect_e(j)*ve(npx,j-1,3)+(1.-geom%edge_vect_e(j))*ve(npx,j,3)
       else
         vt1(j) = geom%edge_vect_e(j)*ve(npx,j+1,1)+(1.-geom%edge_vect_e(j))*ve(npx,j,1)
         vt2(j) = geom%edge_vect_e(j)*ve(npx,j+1,2)+(1.-geom%edge_vect_e(j))*ve(npx,j,2)
         vt3(j) = geom%edge_vect_e(j)*ve(npx,j+1,3)+(1.-geom%edge_vect_e(j))*ve(npx,j,3)
       endif
     enddo
     do j=js,je
       ve(npx,j,1) = vt1(j)
       ve(npx,j,2) = vt2(j)
       ve(npx,j,3) = vt3(j)
     enddo
   endif

   if ( js==1 ) then
     !j = 1
     do i=is,ie
       if ( i>im2 ) then
         ut1(i) = geom%edge_vect_s(i)*ue(i-1,1,1)+(1.-geom%edge_vect_s(i))*ue(i,1,1)
         ut2(i) = geom%edge_vect_s(i)*ue(i-1,1,2)+(1.-geom%edge_vect_s(i))*ue(i,1,2)
         ut3(i) = geom%edge_vect_s(i)*ue(i-1,1,3)+(1.-geom%edge_vect_s(i))*ue(i,1,3)
       else
         ut1(i) = geom%edge_vect_s(i)*ue(i+1,1,1)+(1.-geom%edge_vect_s(i))*ue(i,1,1)
         ut2(i) = geom%edge_vect_s(i)*ue(i+1,1,2)+(1.-geom%edge_vect_s(i))*ue(i,1,2)
         ut3(i) = geom%edge_vect_s(i)*ue(i+1,1,3)+(1.-geom%edge_vect_s(i))*ue(i,1,3)
       endif
     enddo
     do i=is,ie
       ue(i,1,1) = ut1(i)
       ue(i,1,2) = ut2(i)
       ue(i,1,3) = ut3(i)
     enddo
   endif

   if ( (je+1)==npy ) then
     !j = npy
     do i=is,ie
       if ( i>im2 ) then
         ut1(i) = geom%edge_vect_n(i)*ue(i-1,npy,1)+(1.-geom%edge_vect_n(i))*ue(i,npy,1)
         ut2(i) = geom%edge_vect_n(i)*ue(i-1,npy,2)+(1.-geom%edge_vect_n(i))*ue(i,npy,2)
         ut3(i) = geom%edge_vect_n(i)*ue(i-1,npy,3)+(1.-geom%edge_vect_n(i))*ue(i,npy,3)
       else
         ut1(i) = geom%edge_vect_n(i)*ue(i+1,npy,1)+(1.-geom%edge_vect_n(i))*ue(i,npy,1)
         ut2(i) = geom%edge_vect_n(i)*ue(i+1,npy,2)+(1.-geom%edge_vect_n(i))*ue(i,npy,2)
         ut3(i) = geom%edge_vect_n(i)*ue(i+1,npy,3)+(1.-geom%edge_vect_n(i))*ue(i,npy,3)
       endif
     enddo
     do i=is,ie
       ue(i,npy,1) = ut1(i)
       ue(i,npy,2) = ut2(i)
       ue(i,npy,3) = ut3(i)
     enddo
   endif

   do j=js,je+1
     do i=is,ie
       ud(i,j,k) = 0.5*( ue(i,j,1)*geom%es(1,i,j,1) +  &
                         ue(i,j,2)*geom%es(2,i,j,1) +  &
                         ue(i,j,3)*geom%es(3,i,j,1) )
     enddo
   enddo

   do j=js,je
     do i=is,ie+1
       vd(i,j,k) = 0.5*( ve(i,j,1)*geom%ew(1,i,j,2) +  &
                         ve(i,j,2)*geom%ew(2,i,j,2) +  &
                         ve(i,j,3)*geom%ew(3,i,j,2) )
     enddo
   enddo

 enddo

end subroutine a2d



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

end module a2d_mod
