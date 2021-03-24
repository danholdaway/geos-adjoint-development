subroutine fyppmlinear(c, q, flux, jord, ifirst, ilast, jfirst, jlast, npx, npy, dm, ppm_limiter, dya, isd, ied)

! !INPUT PARAMETERS:
 integer, INTENT(IN) :: ifirst, ilast               !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast               !  Y-Dir strip
 integer, INTENT(IN) :: jord
 integer, INTENT(IN) :: npx, npy
 real   , INTENT(IN) :: q(ifirst:ilast,jfirst-ng:jlast+ng)
 real   , intent(in) :: c(isd:ied,js:je+1 )  ! Courant number
 real   , intent(in) :: ppm_limiter
 real   , INTENT(IN) :: dya(:,:)
! !OUTPUT PARAMETERS:
 real   , INTENT(OUT):: flux(ifirst:ilast,jfirst:jlast+1)   !  Flux
 real   , INTENT(OUT)::   dm(ifirst:ilast,jfirst-2:jlast+2)
! Local:
 real al(ifirst:ilast,jfirst-1:jlast+2)
 real bl(ifirst:ilast,jfirst-1:jlast+1)
 real br(ifirst:ilast,jfirst-1:jlast+1)
 real dq(ifirst:ilast,jfirst-3:jlast+2)
 real dl, dr, pmp, lac, ct, qe
 real xt, x0, x1
 integer i, j, js3, je3, je2, jt

!dummy
 integer, parameter :: je=10, js=10
 integer, parameter:: ng    = 3
 integer, parameter:: grid_type = 0

 js3 = max(3,js-1); je3 = min(npy-3,je+1)
 je2 = min(npy-2,je+2)


   do j=js-2,je+2
      do i=ifirst,ilast
         xt = 0.25*(q(i,j+1) - q(i,j-1))
         dm(i,j) = sign(min(abs(xt), max(q(i,j-1), q(i,j), q(i,j+1)) - q(i,j),   &
                            q(i,j) - min(q(i,j-1), q(i,j), q(i,j+1))), xt)
      enddo
   enddo

   if (grid_type < 3) then

      do j=js3,je2
         do i=ifirst,ilast
            !al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
            al(i,j) = (7.0/12.0)*(q(i,j-1)+q(i,j)) -(1.0/12.0)*(q(i,j+1) +q(i,j-2))
         enddo
      enddo

         do j=js3,je3
            do i=ifirst,ilast
               xt = 2.*dm(i,j)
               bl(i,j) = -sign(min(abs(xt), abs(al(i,j  )-q(i,j))), xt)
               br(i,j) =  sign(min(abs(xt), abs(al(i,j+1)-q(i,j))), xt)
            enddo
         enddo
 
!--------------
! Fix the edges:
!--------------
      if( js==1 ) then
         do i=ifirst,ilast
            br(i,2) = al(i,3) - q(i,2)
!           xt = t11*(q(i,0)+q(i,1)) + t12*(q(i,-1)+q(i,2))   &
!              + t13*(dm(i,2)-dm(i,-1))
!!!         xt = 0.75*(q(i,0)+q(i,1)) - 0.25*(q(i,-1)+q(i,2))
            xt = 0.5*((2.*dya(i,1)+dya(i,2))*(q(i,0)+q(i,1))  &
               -dya(i,1)*(q(i,-1)+q(i,2))) / (dya(i,1)+dya(i,2))
            xt = max(0., xt)
            bl(i,1) = xt - q(i,1)
            br(i,0) = xt - q(i,0)
            xt = 4./7.*dm(i,-1) + 11./14.*q(i,-1) + 3./14.*q(i,0)
            xt = max(0., xt)
            bl(i,0) = xt - q(i,0)

            xt = 3./14.*q(i,1) + 11./14.*q(i,2) - 4./7.*dm(i,2)
            xt = max(0., xt)
            br(i,1) = xt - q(i,1)
            bl(i,2) = xt - q(i,2)
         enddo
         do j=0,2
            call pert_ppm_dh(ilast-ifirst+1, q(ifirst:,j), bl(ifirst:,j), br(ifirst:,j), 1)
         enddo
      endif

      if( (je+1)==npy ) then
         do i=ifirst,ilast
            bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)
!           xt = t11*(q(i,npy-1)+q(i,npy)) + t12*(q(i,npy-2)+q(i,npy+1))   &
!               + t13*(dm(i,npy+1)-dm(i,npy-2))
!!!         xt = 0.75*(q(i,npy-1)+q(i,npy)) - 0.25*(q(i,npy-2)+q(i,npy+1))
            xt = 0.5*((2.*dya(i,npy-1)+dya(i,npy-2))*(q(i,npy-1)+q(i,npy)) &
               - dya(i,npy-1)*(q(i,npy-2)+q(i,npy+1)))  &
                / ( dya(i,npy-1)+dya(i,npy-2) )
            xt = max(0., xt)
            br(i,npy-1) = xt - q(i,npy-1)
            bl(i,npy  ) = xt - q(i,npy)
            xt = 3./14.*q(i,npy) + 11./14.*q(i,npy+1) - 4./7.*dm(i,npy+1)
            xt = max(0., xt)
            br(i,npy) = xt - q(i,npy)
            xt = 3./14.*q(i,npy-1) + 11./14.*q(i,npy-2) + 4./7.*dm(i,npy-2)
            xt = max(0., xt)
            br(i,npy-2) = xt - q(i,npy-2)
            bl(i,npy-1) = xt - q(i,npy-1)
         enddo
         do j=npy-2,npy
            call pert_ppm_dh(ilast-ifirst+1, q(ifirst:,j), bl(ifirst:,j), br(ifirst:,j), 1)
         enddo
      endif

   else

      do j=js-1,je+2
         do i=ifirst,ilast
            !al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
            al(i,j) = (7.0/12.0)*(q(i,j-1)+q(i,j)) -(1.0/12.0)*(q(i,j+1) +q(i,j-2))
         enddo
      enddo

         do j=js-1,je+1
            do i=ifirst,ilast
               xt = 2.*dm(i,j)
               bl(i,j) = -sign(min(abs(xt), abs(al(i,j  )-q(i,j))), xt)
               br(i,j) =  sign(min(abs(xt), abs(al(i,j+1)-q(i,j))), xt)
            enddo
         enddo
 
   endif

   do j=js,je+1
      do i=ifirst,ilast

         if( c(i,j)>0. ) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
      
      enddo
   enddo

 end subroutine fyppmlinear


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
