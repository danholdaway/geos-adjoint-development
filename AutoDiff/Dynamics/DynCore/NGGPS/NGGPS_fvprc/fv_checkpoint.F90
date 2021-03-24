module fv_checkpoint_mod

! 20170413 | D.Holdaway | Tool for interfacing with Tapenade checkpointing.

 use fv_arrays_mod, only: FVPRC

 implicit none
 private

 public pushrealarray, poprealarray

 interface pushrealarray
   module procedure pushrealarray_sc
   module procedure pushrealarray_r1
   module procedure pushrealarray_r2
   module procedure pushrealarray_r3
   module procedure pushrealarray_r4
 end interface

 interface poprealarray
   module procedure poprealarray_sc
   module procedure poprealarray_r1
   module procedure poprealarray_r2
   module procedure poprealarray_r3
   module procedure poprealarray_r4
 end interface

 contains

! pushrealarray
! -------------

 subroutine pushrealarray_sc(field,dimen)

  real(FVPRC), intent(inout) :: field
  integer    , optional, intent(in   ) :: dimen

#ifdef SINGLE_FV
  CALL PUSHREAL4ARRAY(field,1)
#else
  CALL PUSHREAL8ARRAY(field,1)
#endif

 endsubroutine pushrealarray_sc

 subroutine pushrealarray_r1(field,dimen)

  real(FVPRC), intent(inout) :: field(:)
  integer    , intent(in   ) :: dimen

#ifdef SINGLE_FV
  CALL PUSHREAL4ARRAY(field(:),dimen)
#else
  CALL PUSHREAL8ARRAY(field(:),dimen)
#endif

 endsubroutine pushrealarray_r1


 subroutine pushrealarray_r2(field,dimen)

  real(FVPRC), intent(inout) :: field(:,:)
  integer    , intent(in   ) :: dimen

#ifdef SINGLE_FV
  CALL PUSHREAL4ARRAY(field(:,:),dimen)
#else
  CALL PUSHREAL8ARRAY(field(:,:),dimen)
#endif

 endsubroutine pushrealarray_r2


 subroutine pushrealarray_r3(field,dimen)

  real(FVPRC), intent(inout) :: field(:,:,:)
  integer    , intent(in   ) :: dimen

#ifdef SINGLE_FV
  CALL PUSHREAL4ARRAY(field(:,:,:),dimen)
#else
  CALL PUSHREAL8ARRAY(field(:,:,:),dimen)
#endif

 endsubroutine pushrealarray_r3


 subroutine pushrealarray_r4(field,dimen)

  real(FVPRC), intent(inout) :: field(:,:,:,:)
  integer    , intent(in   ) :: dimen

#ifdef SINGLE_FV
  CALL PUSHREAL4ARRAY(field(:,:,:,:),dimen)
#else
  CALL PUSHREAL8ARRAY(field(:,:,:,:),dimen)
#endif

 endsubroutine pushrealarray_r4



! poprealarray
! -------------

 subroutine poprealarray_sc(field,dimen)

  real(FVPRC), intent(inout) :: field
  integer    , optional, intent(in   ) :: dimen

#ifdef SINGLE_FV
  CALL POPREAL4ARRAY(field,1)
#else
  CALL POPREAL8ARRAY(field,1)
#endif

 endsubroutine poprealarray_sc

 subroutine poprealarray_r1(field,dimen)

  real(FVPRC), intent(inout) :: field(:)
  integer    , intent(in   ) :: dimen

#ifdef SINGLE_FV
  CALL POPREAL4ARRAY(field(:),dimen)
#else
  CALL POPREAL8ARRAY(field(:),dimen)
#endif

 endsubroutine poprealarray_r1


 subroutine poprealarray_r2(field,dimen)

  real(FVPRC), intent(inout) :: field(:,:)
  integer    , intent(in   ) :: dimen

#ifdef SINGLE_FV
  CALL POPREAL4ARRAY(field(:,:),dimen)
#else
  CALL POPREAL8ARRAY(field(:,:),dimen)
#endif

 endsubroutine poprealarray_r2


 subroutine poprealarray_r3(field,dimen)

  real(FVPRC), intent(inout) :: field(:,:,:)
  integer    , intent(in   ) :: dimen

#ifdef SINGLE_FV
  CALL POPREAL4ARRAY(field(:,:,:),dimen)
#else
  CALL POPREAL8ARRAY(field(:,:,:),dimen)
#endif

 endsubroutine poprealarray_r3


 subroutine poprealarray_r4(field,dimen)

  real(FVPRC), intent(inout) :: field(:,:,:,:)
  integer    , intent(in   ) :: dimen

#ifdef SINGLE_FV
  CALL POPREAL4ARRAY(field(:,:,:,:),dimen)
#else
  CALL POPREAL8ARRAY(field(:,:,:,:),dimen)
#endif

 endsubroutine poprealarray_r4



endmodule fv_checkpoint_mod
