program single_vs_double

implicit none

real(8) :: myvar

real(8), parameter :: cnst_0p20=0.20

real, parameter::  big_number=1.d8
real, parameter:: tiny_number=1.d-8

real, parameter:: ptop_min=1.d-8

myvar = 1.0*5.0


print*, myvar
print*, cnst_0p20

print*, big_number
print*, tiny_number
print*, ptop_min

end program
