module ESMF

implicit none
private


integer, parameter :: esmf_maxstr = 64


public :: esmf_maxstr, spread


contains

real*8 function spread(array,dimen,length)

real*8, intent(in) :: array(:,:)
integer, intent(in) :: dimen, length

integer :: ik, i1, i2

i1 = size(array,1)
i2 = size(array,2)


   spread = sum(array(:,:))


end function

end module ESMF
