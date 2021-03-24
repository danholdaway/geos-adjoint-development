 subroutine fxppm(c, q, flux, iord, ifirst, ilast, jfirst,   &
                  jlast, npx, npy, ppm_limiter, dxa)

 IMPLICIT NONE

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

!Dummy
 integer, parameter :: ie=10, is=10
 real, parameter:: b1 =   1./30.
 real, parameter:: b2 = -13./60.
 real, parameter:: b3 = -13./60.
 real, parameter:: b4 =  0.45
 real, parameter:: b5 = -0.05
 real, parameter:: t11 = 27./28., t12 = -13./28., t13=3./7.
 real, parameter:: s11 = 11./14., s14 = 4./7.,    s15=3./14.
 real, parameter:: c1 = -2./14.
 real, parameter:: c2 = 11./14.
 real, parameter:: c3 =  5./14.
 real, parameter:: p1 =  7./12.
 real, parameter:: p2 = -1./12.
 real, parameter:: r3 = 1./3.
 integer, parameter:: ng    = 3
 integer, parameter:: grid_type = 0 

! Local
 real dm1(ifirst-2:ilast+2)
 real  al(ifirst-1:ilast+2)
 real  bl(ifirst-1:ilast+1)
 real  br(ifirst-1:ilast+1)
 real  dq(ifirst-3:ilast+2)
 real dl, dr, pmp, lac, ct, qe
 real xt, x1, x0
 integer i, j, is3, ie3, ie2, it

 is3 = max(3,is-1);   ie3 = min(npx-3,ie+1)
 ie2 = min(npx-2,ie+2)

 if (iord<=4) then

     do j=jfirst,jlast

        do i=is-2,ie+2
           xt = 0.25*(q(i+1,j) - q(i-1,j))
           dm1(i) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
                             q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
        enddo

      if (grid_type < 3) then
        do i=max(3,is-1),min(npx-2,ie+2)
           al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1) - dm1(i))
        enddo

! Fix the edges:
        if ( is==1 ) then
             x0 = 0.5*((2.*dxa(1,j)+dxa(2,j))*(q(0,j)+q(1,j))   &
                - dxa(1,j)*(q(-1,j)+q(2,j)))/ ( dxa(1,j)+dxa(2,j))
            al(1) = x0
               x1 = s15*q(0,j) + s11*q(-1,j) + s14*dm1(-1)
           dm1(0) = 0.5*(x0 - x1)
           dm1(0) = sign(min(abs(dm1(0)), max(q(0,j), x0, x1) - q(0,j),    &
                                 q(0,j) - min(q(0,j), x0, x1)), dm1(0) )
            al(0) = 0.5*(q(-1,j)+q(0,j)) + r3*(dm1(-1) - dm1(0))
!
               x1 = s15*q(1,j) + s11*q(2,j) - s14*dm1(2)
           dm1(1) = 0.5*(x1 - x0)
           dm1(1) = sign( min(abs(dm1(1)),  max(q(1,j), x0, x1) - q(1,j),  &
                                   q(1,j) - min(q(1,j), x0, x1) ), dm1(1) )
            al(2) = 0.5*(q(1,j)+q(2,j)) + r3*(dm1(1) - dm1(2))
        endif

        if ( (ie+1)==npx ) then
              x0 = 0.5*( (2.*dxa(npx-1,j)+dxa(npx-2,j))*(q(npx-1,j)+q(npx,j))   &
                - dxa(npx-1,j)*(q(npx-2,j)+q(npx+1,j)))/( dxa(npx-1,j)+dxa(npx-2,j))
           al(npx) = x0
              x1 = s15*q(npx-1,j) + s11*q(npx-2,j) + s14*dm1(npx-2)
           dm1(npx-1) = 0.5*(x0 - x1)
           dm1(npx-1) = sign(min(abs(dm1(npx-1)), max(q(npx-1,j), x0, x1) - q(npx-1,j),   &
                                     q(npx-1,j) - min(q(npx-1,j), x0, x1)), dm1(npx-1) )
           al(npx-1) = 0.5*(q(npx-2,j)+q(npx-1,j)) + r3*(dm1(npx-2) - dm1(npx-1))
!
                 x1 = s15*q(npx,j) + s11*q(npx+1,j) - s14*dm1(npx+1)
           dm1(npx) = 0.5*(x1 - x0)
           dm1(npx) = sign(min(abs(dm1(npx)),  max(q(npx,j), x0, x1) - q(npx,j),   &
                                    q(npx,j) - min(q(npx,j), x0, x1)), dm1(npx))
           al(npx+1) = 0.5*(q(npx,j)+q(npx+1,j)) + r3*(dm1(npx) - dm1(npx+1))
        endif
      else
! For doubly periodic BC
           do i=is-1,ie+2
              al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1) - dm1(i))
           enddo
      endif

      if ( iord==3 ) then
           do i=is-1,ie+1
              bl(i) = al(i  ) - q(i,j)
              br(i) = al(i+1) - q(i,j)
           enddo
           call pert_ppm_dh(ie-is+3, q(is-1:,j), bl(is-1:), br(is-1:), 1)
           do i=is,ie+1
              if(c(i,j)>0.) then
                 flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
              else
                 flux(i,j) = q(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
              endif
        enddo
      else
        do i=is,ie+1
          if( c(i,j)>0. ) then
              xt = ppm_limiter*dm1(i-1)
              dl = sign(min(abs(xt), abs(al(i-1)-q(i-1,j))), xt)
              dr = sign(min(abs(xt), abs(al(i  )-q(i-1,j))), xt)
              flux(i,j) = q(i-1,j) + (1.-c(i,j))*(c(i,j)*(dl-dr) + dr)
          else
              xt = ppm_limiter*dm1(i)
              dl = sign(min(abs(xt), abs(al(i  )-q(i,j))), xt)
              dr = sign(min(abs(xt), abs(al(i+1)-q(i,j))), xt)
              flux(i,j) = q(i,j) - (1.+c(i,j))*(c(i,j)*(dl-dr) + dl)
          endif 
        enddo
      endif
     enddo

 elseif (iord==5) then
! PPM with Hunyh's 2nd constraint
     do j=jfirst,jlast
        do i=ifirst-3,ilast+2
           dq(i) = q(i+1,j) - q(i,j)
        enddo

        do i=ifirst-2,ilast+2
           xt = 0.25*(q(i+1,j) - q(i-1,j))
           dm1(i) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
                             q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
        enddo

        do i=ifirst-1,ilast+2
           al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1) - dm1(i))
        enddo

        do i=ifirst-1,ilast+1
           pmp = -2.*dq(i)
           lac = pmp + 1.5*dq(i+1)
           bl(i) = min(max(0., pmp, lac), max(al(i)-q(i,j), min(0.,pmp, lac)))
           pmp = 2.*dq(i-1)
           lac = pmp - 1.5*dq(i-2)
           br(i) = min(max(0., pmp, lac), max(al(i+1)-q(i,j), min(0.,pmp, lac)))
        enddo

        do i=ifirst,ilast+1
           if(c(i,j)>0.) then
              flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
           else
              flux(i,j) = q(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
           endif
        enddo
     enddo

 elseif ( iord==6 .or. iord==7 ) then

     do j=jfirst,jlast

        if ( iord==6 ) then
! Non-monotonic "5th order" scheme (not really 5th order)
          do i=is3, ie3
             bl(i) = b5*q(i-2,j) + b4*q(i-1,j) + b3*q(i,j) + b2*q(i+1,j) + b1*q(i+2,j)
             br(i) = b1*q(i-2,j) + b2*q(i-1,j) + b3*q(i,j) + b4*q(i+1,j) + b5*q(i+2,j)
          enddo
        else
          do i=is-3,ie+2
             dq(i) = q(i+1,j) - q(i,j)
          enddo
          do i=is3, ie3
!-----------------------------------------------
!- Huynh's 2nd constraint + simple mono limiter
!-----------------------------------------------
             dl = b5*q(i-2,j) + b4*q(i-1,j) + b3*q(i,j) + b2*q(i+1,j) + b1*q(i+2,j)
             dr = b1*q(i-2,j) + b2*q(i-1,j) + b3*q(i,j) + b4*q(i+1,j) + b5*q(i+2,j)
               dl = -sign(min(abs(dl), abs(dq(i-1))), dq(i-1))   ! 1st constraint
              pmp = -2.*dq(i)
              lac = pmp + 1.5*dq(i+1)
            bl(i) = min(max(0., pmp, lac), max(dl, min(0.,pmp, lac)))  ! 2nd constraint
!---
               dr = sign(min(abs(dr), abs(dq(i))), dq(i))   ! 1st constraint
              pmp = 2.*dq(i-1)
              lac = pmp - 1.5*dq(i-2)
            br(i) = min(max(0., pmp, lac), max(dr, min(0.,pmp, lac)))
          enddo
        endif

!--------------
! fix the edges
!--------------
        if ( is==1 ) then
             br(2) = p1*(q(2,j)+q(3,j)) + p2*(q(1,j)+q(4,j)) - q(2,j)
             xt = 0.5*((2.*dxa(1,j)+dxa(2,j))*(q(0,j)+q(1,j))   &
                - dxa(1,j)*(q(-1,j)+q(2,j)))/ ( dxa(1,j)+dxa(2,j))
             bl(1) = xt - q(1,j)
             br(0) = xt - q(0,j)

             xt = c1*q(-2,j) + c2*q(-1,j) + c3*q(0,j)
             xt = max( xt, min(q(-1,j),q(0,j)) )
             xt = min( xt, max(q(-1,j),q(0,j)) )
             bl(0) = xt - q(0,j)

             xt = c3*q(1,j) + c2*q(2,j) +c1*q(3,j)
             xt = max( xt, min(q(1,j),q(2,j)) )
             xt = min( xt, max(q(1,j),q(2,j)) )
             br(1) = xt - q(1,j)
             bl(2) = xt - q(2,j)

             if(iord==7) call pert_ppm_dh(3, q(0:,j), bl(0:), br(0:), 1)
        endif

        if ( (ie+1)==npx ) then
             bl(npx-2) = p1*(q(npx-2,j)+q(npx-3,j)) + p2*(q(npx-4,j)+q(npx-1,j)) - q(npx-2,j)
             xt = 0.5*( (2.*dxa(npx-1,j)+dxa(npx-2,j))*(q(npx-1,j)+q(npx,j))   &
                - dxa(npx-1,j)*(q(npx-2,j)+q(npx+1,j)))/( dxa(npx-1,j)+dxa(npx-2,j))

             br(npx-1) = xt - q(npx-1,j)
             bl(npx  ) = xt - q(npx  ,j)

             xt = c3*q(npx,j) + c2*q(npx+1,j) + c1*q(npx+2,j)
             xt = max( xt, min(q(npx,j),q(npx+1,j)) )
             xt = min( xt, max(q(npx,j),q(npx+1,j)) )
             br(npx) = xt - q(npx,j)

             xt = c1*q(npx-3,j) + c2*q(npx-2,j) + c3*q(npx-1,j)
             xt = max( xt, min(q(npx-2,j),q(npx-1,j)) )
             xt = min( xt, max(q(npx-2,j),q(npx-1,j)) )
             br(npx-2) = xt - q(npx-2,j)
             bl(npx-1) = xt - q(npx-1,j)

             if(iord==7) call pert_ppm_dh(3, q(npx-2:,j), bl(npx-2:), br(npx-2:), 1)
        endif

        do i=ifirst,ilast+1
           if(c(i,j)>0.) then
              flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
           else
              flux(i,j) = q(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
           endif
        enddo
     enddo

 elseif( iord<=10 ) then    ! iord=8, 9, 10

     do j=jfirst,jlast

        if (grid_type < 3) then

        do i=is-3,ie+2
           dq(i) = q(i+1,j) - q(i,j)
        enddo

        do i=is-2,ie+2
               xt = 0.25*(q(i+1,j) - q(i-1,j))
           dm1(i) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
                             q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
        enddo

        do i=is3,min(npx-2,ie+2)
           al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1)-dm1(i))
        enddo

        if ( iord==8 ) then
           do i=is3, ie3
              xt = 2.*dm1(i)
              bl(i) = -sign(min(abs(xt), abs(al(i  )-q(i,j))), xt)
              br(i) =  sign(min(abs(xt), abs(al(i+1)-q(i,j))), xt)
           enddo
        elseif ( iord==9 ) then
           do i=is3, ie3
              pmp = -2.*dq(i)
              lac = pmp + 1.5*dq(i+1)
              bl(i) = min(max(0., pmp, lac), max(al(i  )-q(i,j), min(0.,pmp, lac)))
              pmp = 2.*dq(i-1)
              lac = pmp - 1.5*dq(i-2)
              br(i) = min(max(0., pmp, lac), max(al(i+1)-q(i,j), min(0.,pmp, lac)))
           enddo
        else
           do i=is3, ie3
              bl(i) = al(i  ) - q(i,j)
              br(i) = al(i+1) - q(i,j)
              if ( dq(i-1)*dq(i) <= 0. ) then
                   pmp = -2.*dq(i)
                   lac = pmp + 1.5*dq(i+1)
                   bl(i) = min(max(0., pmp, lac), max(bl(i), min(0.,pmp, lac)))
                   pmp = 2.*dq(i-1)
                   lac = pmp - 1.5*dq(i-2)
                   br(i) = min(max(0., pmp, lac), max(br(i), min(0.,pmp, lac)))
              endif
           enddo
        endif

!--------------
! fix the edges
!--------------
           if ( is==1 ) then
              br(2) = al(3) - q(2,j)
!             xt = t11*(q(0,j)+q(1,j)) + t12*(q(-1,j)+q(2,j)) + t13*(dm1(2)-dm1(-1))
              xt = 0.5*((2.*dxa(1,j)+dxa(2,j))*(q(0,j)+q(1,j))   &
                 - dxa(1,j)*(q(-1,j)+q(2,j)))/ ( dxa(1,j)+dxa(2,j))
              bl(1) = xt - q(1,j)
              br(0) = xt - q(0,j)
              xt = s14*dm1(-1) - s11*dq(-1) + q(0,j)

!             xt = max( xt, min(q(-1,j),q(0,j)) )
!             xt = min( xt, max(q(-1,j),q(0,j)) )

              bl(0) = xt - q(0,j)
              xt = s15*q(1,j) + s11*q( 2,j) - s14*dm1( 2)

!             xt = max( xt, min(q(1,j),q(2,j)) )
!             xt = min( xt, max(q(1,j),q(2,j)) )

              br(1) = xt - q(1,j)
              bl(2) = xt - q(2,j)
              call pert_ppm_dh(3, q(0:,j), bl(0:), br(0:), 1)
           endif

           if ( (ie+1)==npx ) then
              bl(npx-2) = al(npx-2) - q(npx-2,j)
!             xt = t11*(q(npx-1,j)+q(npx,j)) + t12*(q(npx-2,j)+q(npx+1,j))   &
!                                            + t13*(dm1(npx+1)-dm1(npx-2))
              xt = 0.5*( (2.*dxa(npx-1,j)+dxa(npx-2,j))*(q(npx-1,j)+q(npx,j))   &
                 - dxa(npx-1,j)*(q(npx-2,j)+q(npx+1,j)))/( dxa(npx-1,j)+dxa(npx-2,j))

              br(npx-1) = xt - q(npx-1,j)
              bl(npx  ) = xt - q(npx  ,j)
              xt = s11*dq(npx) - s14*dm1(npx+1) + q(npx,j)

!             xt = min( xt, max(q(npx,j), q(npx+1,j)) )
!             xt = max( xt, min(q(npx,j), q(npx+1,j)) )

              br(npx) = xt - q(npx,j)
              xt = s15*q(npx-1,j) + s11*q(npx-2,j) + s14*dm1(npx-2)

!             xt = min( xt, max(q(npx-2,j), q(npx-1,j)) )
!             xt = max( xt, min(q(npx-2,j), q(npx-1,j)) )

              br(npx-2) = xt - q(npx-2,j)
              bl(npx-1) = xt - q(npx-1,j)
              call pert_ppm_dh(3, q(npx-2:,j), bl(npx-2:), br(npx-2:), 1)
           endif
        else
!---------------
! grid_type == 4
!---------------
           do i=ifirst-2,ilast+2
              xt = 0.25*(q(i+1,j) - q(i-1,j))
              dm1(i) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
                                q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
           enddo

           do i=ifirst-1,ilast+2
              al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1)-dm1(i))
           enddo

           do i=ifirst-3,ilast+2
              dq(i) = q(i+1,j) - q(i,j)
           enddo
           
           do i=ifirst-1,ilast+1
              pmp = -2.*dq(i)
              lac = pmp + 1.5*dq(i+1)
              bl(i) = min(max(0., pmp, lac), max(al(i  )-q(i,j), min(0.,pmp, lac)))
              pmp = 2.*dq(i-1)
              lac = pmp - 1.5*dq(i-2)
              br(i) = min(max(0., pmp, lac), max(al(i+1)-q(i,j), min(0.,pmp, lac)))
           enddo

        endif     ! grid_type check

        do i=ifirst,ilast+1

             if( c(i,j)>0. ) then
                flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
             else
                flux(i,j) = q(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
             endif

           enddo


        enddo
    else
!------------------------------
! For positive definite tracers:
!------------------------------
! iord=11: PPM mono constraint (Lin 2004)
! iord=12: Huynh 2nd constraint (Lin 2004) + positive definite (Lin & Rood 1996)
! iord>12: positive definite only (Lin & Rood 1996)

       do j=jfirst,jlast

          do i=is-2,ie+2
             xt = 0.25*(q(i+1,j) - q(i-1,j))
             dm1(i) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
                               q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
          enddo

          if (grid_type < 3) then

             do i=is3,ie2
                al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1)-dm1(i))
             enddo

             if ( iord ==11 ) then
                do i=is3,ie3
                   xt = 2.*dm1(i)
                   bl(i) =-sign(min(abs(xt), abs(al(i)  -q(i,j))), xt)
                   br(i) = sign(min(abs(xt), abs(al(i+1)-q(i,j))), xt)
                enddo
             elseif( iord==12 ) then
                do i=is-3,ie+2
                   dq(i) = q(i+1,j) - q(i,j)
                enddo
                do i=is3,ie3
                   pmp = -2.*dq(i)
                   lac = pmp + 1.5*dq(i+1)
                   bl(i) = min(max(0., pmp, lac), max(al(i  )-q(i,j), min(0.,pmp, lac)))
                   pmp = 2.*dq(i-1)
                   lac = pmp - 1.5*dq(i-2)
                   br(i) = min(max(0., pmp, lac), max(al(i+1)-q(i,j), min(0.,pmp, lac)))
                enddo
             else
                do i=is3,ie3
                   bl(i) = al(i  ) - q(i,j)
                   br(i) = al(i+1) - q(i,j)
                enddo
             endif

! Positive definite constraint:
             if(iord/=11) call pert_ppm_dh(ie3-is3+1, q(is3:,j), bl(is3:), br(is3:), 0)

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
                al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1)-dm1(i))
             enddo
             
             if ( iord ==11 ) then
                do i=ifirst-1,ilast+1
                   xt = 2.*dm1(i)
                   bl(i) =-sign(min(abs(xt), abs(al(i)  -q(i,j))), xt)
                   br(i) = sign(min(abs(xt), abs(al(i+1)-q(i,j))), xt)
                enddo
             elseif( iord==12 ) then
                do i=ifirst-3,ilast+2
                   dq(i) = q(i+1,j) - q(i,j)
                enddo
                do i=ifirst-1,ilast+1
                   pmp = -2.*dq(i)
                   lac = pmp + 1.5*dq(i+1)
                   bl(i) = min(max(0., pmp, lac), max(al(i  )-q(i,j), min(0.,pmp, lac)))
                   pmp = 2.*dq(i-1)
                   lac = pmp - 1.5*dq(i-2)
                   br(i) = min(max(0., pmp, lac), max(al(i+1)-q(i,j), min(0.,pmp, lac)))
                enddo
             else
                do i=is-1,ie+1
                   bl(i) = al(i  ) - q(i,j)
                   br(i) = al(i+1) - q(i,j)
                enddo
             endif

! Positive definite constraint:
             if(iord/=11) call pert_ppm_dh(ie-is+3, q(is-1:,j), bl(is-1:), br(is-1:), 0)

          endif
          
          do i=ifirst,ilast+1

             if( c(i,j)>0. ) then
                flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
             else
                flux(i,j) = q(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
             endif


          enddo

       enddo

    endif

 end subroutine fxppm

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

