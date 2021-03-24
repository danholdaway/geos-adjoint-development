module mpp_domains_mod

use fv_arrays_mod, only: REAL8

implicit none
private

public :: domain2d, CGRID_NE, mpp_update_domains, mpp_get_boundary

integer :: CGRID_NE

  type domain2D
     private
     logical :: dummy
  end type domain2D

  interface mpp_update_domains
     module procedure mpp_update_domain2D_2d
     module procedure mpp_update_domain2D_3d
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

  real(REAL8), intent(inout) :: field1(:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_2d
 subroutine mpp_update_domain2D_3d(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(REAL8), intent(inout) :: field1(:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_3d
 subroutine mpp_update_domain2D_4d(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(REAL8), intent(inout) :: field1(:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_4d
 subroutine mpp_update_domain2D_5d(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(REAL8), intent(inout) :: field1(:,:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_5d

 subroutine mpp_update_domain2D_2dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(REAL8), intent(inout) :: field1(:,:), field2(:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2

 endsubroutine mpp_update_domain2D_2dv
 subroutine mpp_update_domain2D_3dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(REAL8), intent(inout) :: field1(:,:,:), field2(:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2


 endsubroutine mpp_update_domain2D_3dv
 subroutine mpp_update_domain2D_4dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(REAL8), intent(inout) :: field1(:,:,:,:), field2(:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2

 endsubroutine mpp_update_domain2D_4dv
 subroutine mpp_update_domain2D_5dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(REAL8), intent(inout) :: field1(:,:,:,:,:), field2(:,:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2

 endsubroutine mpp_update_domain2D_5dv


! mpp_get_boundary
! ----------------

 subroutine mpp_get_boundary(field1,field2,domain,wbuffery,ebuffery,sbufferx,nbufferx,wbufferx,ebufferx,sbuffery,nbuffery,gridtype)

  real(REAL8), intent(inout) :: field1(:,:,:),field2(:,:,:)
  type(domain2d), intent(in) :: domain
  real(REAL8), optional, intent(inout) :: wbuffery(:,:),ebuffery(:,:),sbufferx(:,:),nbufferx(:,:)
  real(REAL8), optional, intent(inout) :: wbufferx(:,:),ebufferx(:,:),sbuffery(:,:),nbuffery(:,:)
  integer, optional, intent(in) :: gridtype

  field1 = 2*field1
  field2 = 2*field2

  wbuffery(:,:) = field1(1,:,:)
  ebuffery(:,:) = field1(2,:,:)
  sbufferx(:,:) = field2(1,:,:)
  nbufferx(:,:) = field2(2,:,:)
  
 endsubroutine mpp_get_boundary

endmodule mpp_domains_mod
