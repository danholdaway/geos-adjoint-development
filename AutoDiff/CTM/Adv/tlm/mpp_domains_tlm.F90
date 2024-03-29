!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:30
!
MODULE MPP_DOMAINS_MOD_D
  USE FV_ARRAYS_MOD, ONLY : real8
  IMPLICIT NONE
  PRIVATE 
  PUBLIC domain2d, cgrid_ne, mpp_update_domains, mpp_get_boundary
  PUBLIC mpp_update_domains_tlm
  INTEGER :: cgrid_ne
  TYPE DOMAIN2D
      PRIVATE 
      LOGICAL :: dummy
  END TYPE DOMAIN2D
  INTERFACE MPP_UPDATE_DOMAINS
      MODULE PROCEDURE MPP_UPDATE_DOMAIN2D_2D
      MODULE PROCEDURE MPP_UPDATE_DOMAIN2D_3D
      MODULE PROCEDURE MPP_UPDATE_DOMAIN2D_4D
      MODULE PROCEDURE MPP_UPDATE_DOMAIN2D_5D
      MODULE PROCEDURE MPP_UPDATE_DOMAIN2D_2DV
      MODULE PROCEDURE MPP_UPDATE_DOMAIN2D_3DV
      MODULE PROCEDURE MPP_UPDATE_DOMAIN2D_4DV
      MODULE PROCEDURE MPP_UPDATE_DOMAIN2D_5DV
  END INTERFACE

  INTERFACE MPP_UPDATE_DOMAINS_TLM
      MODULE PROCEDURE MPP_UPDATE_DOMAIN2D_3D_TLM
      MODULE PROCEDURE MPP_UPDATE_DOMAIN2D_4D_TLM
  END INTERFACE


CONTAINS
! mpp_update_domains
! ------------------
  SUBROUTINE MPP_UPDATE_DOMAIN2D_2D(field1, domain, gridtype, complete, &
&   whalo, ehalo, shalo, nhalo)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: field1(:, :)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER, OPTIONAL, INTENT(IN) :: gridtype
    LOGICAL, OPTIONAL, INTENT(IN) :: complete
    INTEGER, OPTIONAL, INTENT(IN) :: whalo, ehalo, shalo, nhalo
    field1 = 2*field1
  END SUBROUTINE MPP_UPDATE_DOMAIN2D_2D
!  Differentiation of mpp_update_domain2d_3d in forward (tangent) mode:
!   variations   of useful results: field1
!   with respect to varying inputs: field1
  SUBROUTINE MPP_UPDATE_DOMAIN2D_3D_TLM(field1, field1_tl, domain, &
&   gridtype, complete, whalo, ehalo, shalo, nhalo)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: field1(:, :, :)
    REAL(real8), INTENT(INOUT) :: field1_tl(:, :, :)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER, OPTIONAL, INTENT(IN) :: gridtype
    LOGICAL, OPTIONAL, INTENT(IN) :: complete
    INTEGER, OPTIONAL, INTENT(IN) :: whalo, ehalo, shalo, nhalo
    field1_tl = 2*field1_tl
    field1 = 2*field1
  END SUBROUTINE MPP_UPDATE_DOMAIN2D_3D_TLM
  SUBROUTINE MPP_UPDATE_DOMAIN2D_3D(field1, domain, gridtype, complete, &
&   whalo, ehalo, shalo, nhalo)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: field1(:, :, :)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER, OPTIONAL, INTENT(IN) :: gridtype
    LOGICAL, OPTIONAL, INTENT(IN) :: complete
    INTEGER, OPTIONAL, INTENT(IN) :: whalo, ehalo, shalo, nhalo
    field1 = 2*field1
  END SUBROUTINE MPP_UPDATE_DOMAIN2D_3D
!  Differentiation of mpp_update_domain2d_4d in forward (tangent) mode:
!   variations   of useful results: field1
!   with respect to varying inputs: field1
  SUBROUTINE MPP_UPDATE_DOMAIN2D_4D_TLM(field1, field1_tl, domain, &
&   gridtype, complete, whalo, ehalo, shalo, nhalo)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: field1(:, :, :, :)
    REAL(real8), INTENT(INOUT) :: field1_tl(:, :, :, :)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER, OPTIONAL, INTENT(IN) :: gridtype
    LOGICAL, OPTIONAL, INTENT(IN) :: complete
    INTEGER, OPTIONAL, INTENT(IN) :: whalo, ehalo, shalo, nhalo
    field1_tl = 2*field1_tl
    field1 = 2*field1
  END SUBROUTINE MPP_UPDATE_DOMAIN2D_4D_TLM
  SUBROUTINE MPP_UPDATE_DOMAIN2D_4D(field1, domain, gridtype, complete, &
&   whalo, ehalo, shalo, nhalo)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: field1(:, :, :, :)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER, OPTIONAL, INTENT(IN) :: gridtype
    LOGICAL, OPTIONAL, INTENT(IN) :: complete
    INTEGER, OPTIONAL, INTENT(IN) :: whalo, ehalo, shalo, nhalo
    field1 = 2*field1
  END SUBROUTINE MPP_UPDATE_DOMAIN2D_4D
  SUBROUTINE MPP_UPDATE_DOMAIN2D_5D(field1, domain, gridtype, complete, &
&   whalo, ehalo, shalo, nhalo)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: field1(:, :, :, :, :)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER, OPTIONAL, INTENT(IN) :: gridtype
    LOGICAL, OPTIONAL, INTENT(IN) :: complete
    INTEGER, OPTIONAL, INTENT(IN) :: whalo, ehalo, shalo, nhalo
    field1 = 2*field1
  END SUBROUTINE MPP_UPDATE_DOMAIN2D_5D
  SUBROUTINE MPP_UPDATE_DOMAIN2D_2DV(field1, field2, domain, gridtype, &
&   complete, whalo, ehalo, shalo, nhalo)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: field1(:, :), field2(:, :)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER, OPTIONAL, INTENT(IN) :: gridtype
    LOGICAL, OPTIONAL, INTENT(IN) :: complete
    INTEGER, OPTIONAL, INTENT(IN) :: whalo, ehalo, shalo, nhalo
    field1 = 2*field1
    field2 = 2*field2
  END SUBROUTINE MPP_UPDATE_DOMAIN2D_2DV
  SUBROUTINE MPP_UPDATE_DOMAIN2D_3DV(field1, field2, domain, gridtype, &
&   complete, whalo, ehalo, shalo, nhalo)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: field1(:, :, :), field2(:, :, :)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER, OPTIONAL, INTENT(IN) :: gridtype
    LOGICAL, OPTIONAL, INTENT(IN) :: complete
    INTEGER, OPTIONAL, INTENT(IN) :: whalo, ehalo, shalo, nhalo
    field1 = 2*field1
    field2 = 2*field2
  END SUBROUTINE MPP_UPDATE_DOMAIN2D_3DV
  SUBROUTINE MPP_UPDATE_DOMAIN2D_4DV(field1, field2, domain, gridtype, &
&   complete, whalo, ehalo, shalo, nhalo)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: field1(:, :, :, :), field2(:, :, :, :)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER, OPTIONAL, INTENT(IN) :: gridtype
    LOGICAL, OPTIONAL, INTENT(IN) :: complete
    INTEGER, OPTIONAL, INTENT(IN) :: whalo, ehalo, shalo, nhalo
    field1 = 2*field1
    field2 = 2*field2
  END SUBROUTINE MPP_UPDATE_DOMAIN2D_4DV
  SUBROUTINE MPP_UPDATE_DOMAIN2D_5DV(field1, field2, domain, gridtype, &
&   complete, whalo, ehalo, shalo, nhalo)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: field1(:, :, :, :, :), field2(:, :, :&
&   , :, :)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER, OPTIONAL, INTENT(IN) :: gridtype
    LOGICAL, OPTIONAL, INTENT(IN) :: complete
    INTEGER, OPTIONAL, INTENT(IN) :: whalo, ehalo, shalo, nhalo
    field1 = 2*field1
    field2 = 2*field2
  END SUBROUTINE MPP_UPDATE_DOMAIN2D_5DV
! mpp_get_boundary
! ----------------
  SUBROUTINE MPP_GET_BOUNDARY(field1, field2, domain, wbuffery, ebuffery&
&   , sbufferx, nbufferx, wbufferx, ebufferx, sbuffery, nbuffery, &
&   gridtype)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: field1(:, :, :), field2(:, :, :)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    REAL(real8), OPTIONAL, INTENT(INOUT) :: wbuffery(:, :), ebuffery(:, &
&   :), sbufferx(:, :), nbufferx(:, :)
    REAL(real8), OPTIONAL, INTENT(INOUT) :: wbufferx(:, :), ebufferx(:, &
&   :), sbuffery(:, :), nbuffery(:, :)
    INTEGER, OPTIONAL, INTENT(IN) :: gridtype
    field1 = 2*field1
    field2 = 2*field2
    wbuffery(:, :) = field1(1, :, :)
    ebuffery(:, :) = field1(2, :, :)
    sbufferx(:, :) = field2(1, :, :)
    nbufferx(:, :) = field2(2, :, :)
  END SUBROUTINE MPP_GET_BOUNDARY
END MODULE MPP_DOMAINS_MOD_D
