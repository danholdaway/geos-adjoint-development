module dyn_core_mod

use array_mod

implicit none
private

contains

subroutine dyn_core(bd,npz,u,v,us,vs)

 type(grid), intent(in) :: bd
 integer, intent(in) :: npz
 real, intent(inout) :: u(bd%is:bd%ie,bd%js:bd%ie,npz+1)
 real, intent(inout) :: us(bd%is:bd%ie,bd%js:bd%ie)
 real, intent(inout) :: v(bd%is:bd%ie,bd%js:bd%ie,npz+1)
 real, intent(inout) :: vs(bd%is:bd%ie,bd%js:bd%ie)

 integer :: i,j,k
 integer :: is,ie,js,je

 is = bd%is
 ie = bd%ie
 js = bd%js
 je = bd%je

 do j=js,je
    do i=is,ie
       u(i,j,npz+1) = us(i,j)
    enddo
    do k=npz,1,-1
       do i=is,ie
          u(i,j,k) = u(i,j,k+1) - u(i,j,k)*v(i,j,k)
       enddo
    enddo
 enddo

 do j=js,je
    do i=is,ie
       v(i,j,npz+1) = vs(i,j)
    enddo
    do k=npz,1,-1
       do i=is,ie
          v(i,j,k) = v(i,j,k+1) - v(i,j,k)*u(i,j,k)
       enddo
    enddo
 enddo

 call again(is,ie,js,je,npz,u,v,us,vs)

endsubroutine dyn_core

subroutine again(is,ie,js,je,npz,u,v,us,vs)

 integer, intent(in) :: is,ie,js,je,npz
 real, intent(inout) :: u(is:ie,js:ie,npz+1)
 real, intent(inout) :: us(is:ie,js:ie)
 real, intent(inout) :: v(is:ie,js:ie,npz+1)
 real, intent(inout) :: vs(is:ie,js:ie)

 integer :: i,j,k

 do j=js,je
    do i=is,ie
       u(i,j,npz+1) = us(i,j)
    enddo
    do k=npz,1,-1
       do i=is,ie
          u(i,j,k) = u(i,j,k+1) - u(i,j,k)*v(i,j,k)
       enddo
    enddo
 enddo

 do j=js,je
    do i=is,ie
       v(i,j,npz+1) = vs(i,j)
    enddo
    do k=npz,1,-1
       do i=is,ie
          v(i,j,k) = v(i,j,k+1) - v(i,j,k)*u(i,j,k)
       enddo
    enddo
 enddo 

endsubroutine again

endmodule dyn_core_mod

