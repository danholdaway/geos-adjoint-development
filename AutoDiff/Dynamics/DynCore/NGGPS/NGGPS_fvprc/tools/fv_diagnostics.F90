module fv_diagnostics_mod

implicit none
private

public fv_time, prt_mxm, range_check, prt_minmax

logical, parameter :: prt_minmax = .false.
integer :: fv_time !dummy

contains


 subroutine range_check(stub)

  integer, intent(in) :: stub

 endsubroutine range_check


 subroutine prt_mxm(stub)

  integer, intent(in) :: stub

 endsubroutine prt_mxm


end module fv_diagnostics_mod

