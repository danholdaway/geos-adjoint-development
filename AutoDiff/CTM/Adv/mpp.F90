module mpp_mod

use fv_arrays_mod, only: REAL8

implicit none
private

public mpp_sum

interface mpp_sum
 module procedure mpp_sum_r1
 module procedure mpp_sum_r2
end interface

contains


 subroutine mpp_sum_r1(psums,dims)

  real(REAL8), intent(inout) :: psums(:)
  integer, intent(in) :: dims
  
  psums = sum(psums)

 endsubroutine mpp_sum_r1

 subroutine mpp_sum_r2(psums,dims)

  real(REAL8), intent(inout) :: psums(:,:)
  integer, intent(in) :: dims
  
  psums = sum(psums)

 endsubroutine mpp_sum_r2

end module mpp_mod
