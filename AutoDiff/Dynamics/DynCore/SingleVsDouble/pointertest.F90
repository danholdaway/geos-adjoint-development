program pointertest

implicit none

integer, parameter :: dim = 10000

real(8), target :: x(dim,dim)
real(8) :: xp1(dim,dim)
real(8), pointer :: xp2(:,:)

integer(8) :: i,j,k,l

real(8) :: calc

real :: t1,t2

call random_number(x)

call cpu_time(t1)

  do i = 1,10
    do j = 1,10

       xp1 => x

    enddo
  enddo

call cpu_time(t2)
print*, t2-t1

call cpu_time(t1)

  do i = 1,10
    do j = 1,10

       xp2 => x

    enddo
  enddo

call cpu_time(t2)
print*, t2-t1


end program
