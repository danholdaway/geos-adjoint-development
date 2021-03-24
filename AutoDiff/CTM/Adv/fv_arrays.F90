module fv_arrays_mod

 private

  integer, parameter :: REAL4 = 4 !kind(1.00)
  integer, parameter :: REAL8 = 8 !kind(1.d0)
  integer, parameter :: FVPRC = REAL8  ! Full Build Precision for the model

  public:: FVPRC, REAL4, REAL8, CNVT

  INTERFACE CNVT 
    module procedure cnvt0
    module procedure cnvt1
    module procedure cnvt2
    module procedure cnvt3
  END INTERFACE

  contains

  real(FVPRC) function cnvt0(dbl_var)
     real(REAL8), intent(IN) :: dbl_var
     cnvt0 = SIGN(MIN(huge_number,MAX(tiny_number,ABS(dbl_var))),dbl_var)
  end function

  function cnvt1(dbl_var)
     real(REAL8), intent(IN) :: dbl_var(:)
     real(FVPRC)  :: cnvt1(LBOUND(dbl_var,1):UBOUND(dbl_var,1))
     integer :: i
     do i=LBOUND(dbl_var,1),UBOUND(dbl_var,1)
        cnvt1(i) = SIGN(MIN(huge_number,MAX(tiny_number,ABS(dbl_var(i)))),dbl_var(i))
     enddo
  end function

  function cnvt2(dbl_var)
     real(REAL8), intent(IN) :: dbl_var(:,:)
     real(FVPRC)  :: cnvt2(LBOUND(dbl_var,1):UBOUND(dbl_var,1),&
                           LBOUND(dbl_var,2):UBOUND(dbl_var,2))
     integer :: i, j
     do j=LBOUND(dbl_var,2),UBOUND(dbl_var,2)
        do i=LBOUND(dbl_var,1),UBOUND(dbl_var,1)
           cnvt2(i,j) = SIGN(MIN(huge_number,MAX(tiny_number,ABS(dbl_var(i,j)))),dbl_var(i,j))
        enddo
     enddo
  end function

  function cnvt3(dbl_var)
     real(REAL8), intent(IN) :: dbl_var(:,:,:)
     real(FVPRC)  :: cnvt3(LBOUND(dbl_var,1):UBOUND(dbl_var,1),&
                           LBOUND(dbl_var,2):UBOUND(dbl_var,2),&
                           LBOUND(dbl_var,3):UBOUND(dbl_var,3))
     integer :: i, j, k
     do k=LBOUND(dbl_var,3),UBOUND(dbl_var,3)
        do j=LBOUND(dbl_var,2),UBOUND(dbl_var,2)
           do i=LBOUND(dbl_var,1),UBOUND(dbl_var,1)
              cnvt3(i,j,k) = SIGN(MIN(huge_number,MAX(tiny_number,ABS(dbl_var(i,j,k)))),dbl_var(i,j,k))
           enddo
        enddo
     enddo
  end function

end module fv_arrays_mod

