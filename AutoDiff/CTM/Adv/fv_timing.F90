module fv_timing_mod

implicit none
private

integer :: NAME_LENGTH = 64

public timing_on, timing_off


contains



 subroutine timing_on(timechar)

  character(len=NAME_LENGTH), intent(in) :: timechar

  !print*, timechar

 end subroutine timing_on


 subroutine timing_off(timechar)

  character(len=NAME_LENGTH), intent(in) :: timechar

  !print*, timechar

 end subroutine timing_off



end module fv_timing_mod

