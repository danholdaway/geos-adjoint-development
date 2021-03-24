program single_vs_double

implicit none

integer, parameter :: im = 2000 !Grid points
integer, parameter :: nh = 2 !Halo points

real(8) :: x(im+nh,1)
real(8) :: T(im,im+nh)
real(8) :: y(im,1)
real(8) :: xhat(im+nh,1)

real(8) :: dotd(2)
real(16) :: dotq(2)

integer :: i

call random_number(x)
call random_number(T)

y = matmul(T,x)

dotd = 0.0_8
dotq = 0.0_16

do i = 1,im
   dotd(1) = dotd(1) + y(i,1)*y(i,1)
   dotq(1) = dotq(1) + y(i,1)*y(i,1)
enddo

!Introduce adjoint bug
!T(1,1) = T(1,1)+1e-5

xhat = matmul(transpose(T),y)

do i = 1,im+nh
   dotd(2) = dotd(2) + xhat(i,1)*x(i,1)
   dotq(2) = dotq(2) + xhat(i,1)*x(i,1)
enddo


print*, dotd(1)
print*, dotd(2)
print*, (dotd(1) - dotd(2))/dotd(1)

print*, dotq(1)
print*, dotq(2)
print*, (dotq(1) - dotq(2))/dotq(1)


end program
