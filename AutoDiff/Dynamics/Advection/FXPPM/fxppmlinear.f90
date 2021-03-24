 subroutine fxppmlinear(c, q, flux, iord, ifirst, ilast, jfirst, jlast, npx, npy, ppm_limiter, dxa)

! !INPUT PARAMETERS:
 integer, INTENT(IN) :: ifirst, ilast               !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast               !  Y-Dir strip
 integer, INTENT(IN) :: iord
 integer, INTENT(IN) :: npx, npy
 real   , INTENT(IN) :: q(ifirst-ng:ilast+ng,jfirst:jlast)
 real   , INTENT(IN) :: c(ifirst   :ilast+1 ,jfirst:jlast) ! Courant   N (like FLUX)
 real   , INTENT(IN) :: ppm_limiter
 real   , INTENT(IN) :: dxa(:,:)
! !OUTPUT PARAMETERS:
 real   , INTENT(OUT) :: flux(ifirst:ilast+1,jfirst:jlast) !  Flux
! Local
 real dm1(ifirst-2:ilast+2)
 real  al(ifirst-1:ilast+2)
 real  bl(ifirst-1:ilast+1)
 real  br(ifirst-1:ilast+1)
 real  dq(ifirst-3:ilast+2)
 real dl, dr, pmp, lac, ct, qe
 real xt, x1, x0
 integer i, j, is3, ie3, ie2, it

!dummy
 integer, parameter :: ie=10, is=10
 integer, parameter:: ng    = 3
 integer, parameter:: grid_type = 0 

 is3 = max(3,is-1);   ie3 = min(npx-3,ie+1)
 ie2 = min(npx-2,ie+2)

       do j=jfirst,jlast

          do i=is-2,ie+2
             xt = 0.25*(q(i+1,j) - q(i-1,j))
             dm1(i) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
                               q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
          enddo

          if (grid_type < 3) then

             do i=is3,ie2
                !al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1)-dm1(i))
                al(i) = (7.0/12.0)*(q(i-1,j)+q(i,j)) - (1.0/12.0)*(q(i+1,j)+q(i-2,j))
             enddo

                do i=is3,ie3
                   xt = 2.*dm1(i)
                   bl(i) =-sign(min(abs(xt), abs(al(i)  -q(i,j))), xt)
                   br(i) = sign(min(abs(xt), abs(al(i+1)-q(i,j))), xt)
                enddo
!--------------
! fix the edges
!--------------
             if ( is==1 ) then
                br(2) = al(3) - q(2,j)
!               xt = t11*(q(0,j)+q(1,j)) + t12*(q(-1,j)+q(2,j)) + t13*(dm1(2)-dm1(-1))
!!!             xt = 0.75*(q(0,j)+q(1,j)) - 0.25*(q(-1,j)+q(2,j))
                xt = 0.5*( (2.*dxa(1,j)+dxa(2,j))*(q(0,j)+q(1,j))  &
                   - dxa(1,j)*(q(-1,j)+q(2,j)) ) / ( dxa(1,j)+dxa(2,j) )
                xt = max(0., xt)
                bl(1) = xt - q(1,j)
                br(0) = xt - q(0,j)
                xt = 4./7.*dm1(-1) + 11./14.*q(-1,j) + 3./14.*q(0,j)
                xt = max(0., xt)
                bl(0) =  xt - q(0,j)
                xt = 3./14.*q(1,j) + 11./14.*q(2,j) - 4./7.*dm1(2)
                xt = max(0., xt)
                br(1) = xt - q(1,j)
                bl(2) = xt - q(2,j)
                call pert_ppm_dh(3, q(0:,j), bl(0:), br(0:), 1)
             endif

             if ( (ie+1)==npx ) then
                bl(npx-2) = al(npx-2) - q(npx-2,j)
!               xt = t11*(q(npx-1,j)+q(npx,j)) + t12*(q(npx-2,j)+q(npx+1,j))   &
!                  + t13*(dm1(npx+1)-dm1(npx-2))
!!!             xt = 0.75*(q(npx-1,j)+q(npx,j)) - 0.25*(q(npx-2,j)+q(npx+1,j))
                xt = 0.5*((2.*dxa(npx-1,j)+dxa(npx-2,j))*(q(npx-1,j)+q(npx,j)) -   &
                     dxa(npx-1,j)*(q(npx-2,j)+q(npx+1,j)) )  &
                 / ( dxa(npx-1,j)+dxa(npx-2,j) )
                xt = max(0., xt)
                br(npx-1) = xt - q(npx-1,j)
                bl(npx  ) = xt - q(npx  ,j)
!               br(npx) = 11./14.*q(npx+1,j) + 3./14.*q(npx,j) - 4./7.*dm1(npx+1)
                xt = 11./14.*q(npx+1,j) + 3./14.*q(npx,j) - 4./7.*dm1(npx+1)
                xt = max(0., xt)
                br(npx) = xt - q(npx,j)
                xt = 3./14.*q(npx-1,j) + 11./14.*q(npx-2,j) + 4./7.*dm1(npx-2)
                xt = max(0., xt)
                br(npx-2) = xt - q(npx-2,j)
                bl(npx-1) = xt - q(npx-1,j)
                call pert_ppm_dh(3, q(npx-2:,j), bl(npx-2:), br(npx-2:), 1)
             endif
          else
!--------------
! grid_type >=4
!--------------
             do i=ifirst-1,ilast+2
                !al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1)-dm1(i))
                al(i) = (7.0/12.0)*(q(i-1,j)+q(i,j)) - (1.0/12.0)*(q(i+1,j)+q(i-2,j))
             enddo
             
                do i=ifirst-1,ilast+1
                   xt = 2.*dm1(i)
                   bl(i) =-sign(min(abs(xt), abs(al(i)  -q(i,j))), xt)
                   br(i) = sign(min(abs(xt), abs(al(i+1)-q(i,j))), xt)
                enddo
          endif
          
          do i=ifirst,ilast+1

     
             if( c(i,j)>0. ) then
                flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
             else
                flux(i,j) = q(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
             endif

          enddo

       enddo

 end subroutine fxppmlinear


 subroutine pert_ppm_dh(im, a0, al, ar, iv)

 implicit none

 integer, intent(in):: im
 integer, intent(in):: iv
 real*4, intent(in)   :: a0(im)
 real*4, intent(inout):: al(im), ar(im)
! Local:
 real*4 a4, da1, da2, a6da, fmin
 integer i
 real*4, parameter:: r12 = 1./12.

!-----------------------------------
! Optimized PPM in perturbation form:
!-----------------------------------

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
        a4 = -3.*(ar(i) + al(i))
       da1 =      ar(i) - al(i)
      if( abs(da1) < -a4 ) then
         fmin = a0(i) + 0.25/a4*da1**2 + a4*r12
         if( fmin < 0. ) then
             if( ar(i)>0. .and. al(i)>0. ) then
                 ar(i) = 0.
                 al(i) = 0.
             elseif( da1 > 0. ) then
                 ar(i) = -2.*al(i)
             else
                 al(i) = -2.*ar(i)
             endif
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
       if ( al(i)*ar(i) < 0. ) then
            da1 = al(i) - ar(i)
            da2 = da1**2
            a6da = 3.*(al(i)+ar(i))*da1
            if( a6da < -da2 ) then
                ar(i) = -2.*al(i)
            elseif( a6da > da2 ) then
                al(i) = -2.*ar(i)
            endif
       else
! effect of dm=0 included here
            al(i) = 0.
            ar(i) = 0.
       endif
  enddo
 endif

 end subroutine pert_ppm_dh
