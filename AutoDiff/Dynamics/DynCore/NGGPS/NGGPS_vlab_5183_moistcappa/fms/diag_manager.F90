MODULE diag_manager_mod

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: send_data

  INTERFACE send_data
     MODULE PROCEDURE send_data_0d
     MODULE PROCEDURE send_data_1d
     MODULE PROCEDURE send_data_2d
     MODULE PROCEDURE send_data_3d
  END INTERFACE

CONTAINS

  LOGICAL FUNCTION send_data_0d(diag_field_id, field, time, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, INTENT(in) :: field
    INTEGER, INTENT(in), OPTIONAL :: time
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL :: field_out(1, 1, 1)

    send_data_0d = .false.

  END FUNCTION send_data_0d

  LOGICAL FUNCTION send_data_1d(diag_field_id, field, time, is_in, mask, rmask, ie_in, weight, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, DIMENSION(:), INTENT(in) :: field
    REAL, INTENT(in), OPTIONAL :: weight
    REAL, INTENT(in), DIMENSION(:), OPTIONAL :: rmask
    INTEGER, INTENT(in), OPTIONAL :: time
    INTEGER, INTENT(in), OPTIONAL :: is_in, ie_in
    LOGICAL, INTENT(in), DIMENSION(:), OPTIONAL :: mask
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL, DIMENSION(SIZE(field(:)), 1, 1) :: field_out
    LOGICAL, DIMENSION(SIZE(field(:)), 1, 1) ::  mask_out

    send_data_1d = .false.

  END FUNCTION send_data_1d

  LOGICAL FUNCTION send_data_2d(diag_field_id, field, time, is_in, js_in, &
       & mask, rmask, ie_in, je_in, weight, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, INTENT(in), DIMENSION(:,:) :: field
    REAL, INTENT(in), OPTIONAL :: weight
    INTEGER, INTENT(in), OPTIONAL :: time
    INTEGER, INTENT(in), OPTIONAL :: is_in, js_in, ie_in, je_in
    LOGICAL, INTENT(in), DIMENSION(:,:), OPTIONAL :: mask
    REAL, INTENT(in), DIMENSION(:,:),OPTIONAL :: rmask
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL, DIMENSION(SIZE(field,1),SIZE(field,2),1) :: field_out
    LOGICAL, DIMENSION(SIZE(field,1),SIZE(field,2),1) ::  mask_out

    send_data_2d = .false.

  END FUNCTION send_data_2d

  LOGICAL FUNCTION send_data_3d(diag_field_id, field, time, is_in, js_in, ks_in, &
             & mask, rmask, ie_in, je_in, ke_in, weight, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, DIMENSION(:,:,:), INTENT(in) :: field
    REAL, INTENT(in), OPTIONAL :: weight
    INTEGER, INTENT(in), OPTIONAL :: time
    INTEGER, INTENT(in), OPTIONAL :: is_in, js_in, ks_in,ie_in,je_in, ke_in
    LOGICAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask
    REAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: rmask
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    send_data_3d = .false.

  END FUNCTION send_data_3d

END MODULE diag_manager_mod
