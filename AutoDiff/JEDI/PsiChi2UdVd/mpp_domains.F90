module mpp_domains_mod

implicit none
private

public mpp_update_domains, domain2d, kind_real

  integer, parameter :: kind_real = 8

  type domain2D
   integer :: dummy
  end type domain2D   

  interface mpp_update_domains
     module procedure mpp_update_domain2D_2d
     module procedure mpp_update_domain2D_2d_r8
     module procedure mpp_update_domain2D_3d
     module procedure mpp_update_domain2D_3d_r8
     module procedure mpp_update_domain2D_4d
     module procedure mpp_update_domain2D_5d
     module procedure mpp_update_domain2D_2dv
     module procedure mpp_update_domain2D_3dv
     module procedure mpp_update_domain2D_4dv
     module procedure mpp_update_domain2D_5dv
  end interface

 contains




! mpp_update_domains
! ------------------

 subroutine mpp_update_domain2D_2d(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real, intent(inout) :: field1(:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_2d
 subroutine mpp_update_domain2D_2d_r8(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(kind=8), intent(inout) :: field1(:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_2d_r8
 subroutine mpp_update_domain2D_3d(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real, intent(inout) :: field1(:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_3d
 subroutine mpp_update_domain2D_3d_r8(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(kind=8), intent(inout) :: field1(:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_3d_r8
 subroutine mpp_update_domain2D_4d(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real, intent(inout) :: field1(:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_4d
 subroutine mpp_update_domain2D_5d(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real, intent(inout) :: field1(:,:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_5d

 subroutine mpp_update_domain2D_2dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real, intent(inout) :: field1(:,:), field2(:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2

 endsubroutine mpp_update_domain2D_2dv
 subroutine mpp_update_domain2D_3dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real, intent(inout) :: field1(:,:,:), field2(:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2


 endsubroutine mpp_update_domain2D_3dv
 subroutine mpp_update_domain2D_4dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real, intent(inout) :: field1(:,:,:,:), field2(:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2

 endsubroutine mpp_update_domain2D_4dv
 subroutine mpp_update_domain2D_5dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real, intent(inout) :: field1(:,:,:,:,:), field2(:,:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2

 endsubroutine mpp_update_domain2D_5dv

end module mpp_domains_mod
