module tp_core_mod
!BOP
!
! !MODULE: tp_core --- A collection of routines to support FV transport
!
 use fv_arrays_mod,     only: FVPRC, REAL4, REAL8, CNVT
 use fv_mp_mod,         only: is,js,ie,je, ng, isd,jsd,ied,jed
 use fv_grid_utils_mod, only: sw_corner, se_corner, ne_corner, nw_corner, &
                              sin_sg=>sin_sg, &
                              da_min=>da_min, big_number_R8=>big_number
 use fv_grid_tools_mod, only: grid_type, &
                                dx=>  dx,    dy=>   dy, &
                              rdxc=>rdxc,  rdyc=> rdyc, &
                               dxa=> dxa,   dya=>  dya, &
                              area=>area, rarea=>rarea

 implicit none

 private
 public fv_tp_2d, pert_ppm, copy_corners

 INTERFACE fv_tp_2d
   MODULE PROCEDURE fv_tp_2d_r4
!#ifdef SINGLE_FV
!   MODULE PROCEDURE fv_tp_2d_r8
!#endif
 END INTERFACE

 INTERFACE copy_corners
   MODULE PROCEDURE copy_corners_r4
   MODULE PROCEDURE copy_corners_r8
 END INTERFACE

!#ifdef MAPL_MODE
 logical, parameter :: nested=.false.
!#endif

 real(FVPRC), parameter:: r3 = 1./3.
 real(FVPRC), parameter:: near_zero = 1.E-25
 real(FVPRC), parameter:: ppm_limiter = 2.0
 real(FVPRC), parameter:: big_number = big_number_R8

!#ifdef WAVE_FORM
!! Suresh & Huynh scheme 2.2 (purtabation form)
!! The wave-form is more diffusive than scheme 2.1
! real(FVPRC), parameter:: b1 =   0.0375
! real(FVPRC), parameter:: b2 =  -7./30.
! real(FVPRC), parameter:: b3 =  -23./120.
! real(FVPRC), parameter:: b4 =  13./30.
! real(FVPRC), parameter:: b5 = -11./240.
!#else
! scheme 2.1: perturbation form
 real(FVPRC), parameter:: b1 =   1./30.
 real(FVPRC), parameter:: b2 = -13./60.
 real(FVPRC), parameter:: b3 = -13./60.
 real(FVPRC), parameter:: b4 =  0.45
 real(FVPRC), parameter:: b5 = -0.05
!#endif
 real(FVPRC), parameter:: t11 = 27./28., t12 = -13./28., t13=3./7.
 real(FVPRC), parameter:: s11 = 11./14., s14 = 4./7.,    s15=3./14.
!----------------------------------------------------
! volume-conserving cubic with 2nd drv=0 at end point:
!----------------------------------------------------
! Non-monotonic
  real(FVPRC), parameter:: c1 = -2./14.
  real(FVPRC), parameter:: c2 = 11./14.
  real(FVPRC), parameter:: c3 =  5./14.
!----------------------
! PPM volume mean form:
!----------------------
  real(FVPRC), parameter:: p1 =  7./12.     ! 0.58333333
  real(FVPRC), parameter:: p2 = -1./12.


!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

CONTAINS

!#ifdef SINGLE_FV
! subroutine fv_tp_2d_r8(qIN, crx, cry, npx, npy, hord, fx, fy, xfx, yfx,  &
!                     ra_x, ra_y, mfx, mfy, mass, nord, damp_c)
!   integer, intent(in):: npx, npy
!   integer, intent(in)::hord
!
!   real(FVPRC), intent(in)::  crx(is:ie+1,jsd:jed)  !
!   real(FVPRC), intent(in)::  xfx(is:ie+1,jsd:jed)  !
!   real(FVPRC), intent(in)::  cry(isd:ied,js:je+1 )  !
!   real(FVPRC), intent(in)::  yfx(isd:ied,js:je+1 )  !
!   real(FVPRC), intent(in):: ra_x(is:ie,jsd:jed)
!   real(FVPRC), intent(in):: ra_y(isd:ied,js:je)
!   real(REAL8), intent(in)::  qIN(isd:ied,jsd:jed)  ! transported scalar
!   real(FVPRC), intent(out)::fx(is:ie+1 ,js:je)    ! Flux in x ( E )
!   real(FVPRC), intent(out)::fy(is:ie,   js:je+1 )    ! Flux in y ( N )
!! optional Arguments:
!   real(FVPRC), OPTIONAL, intent(in):: mfx(is:ie+1,js:je  )  ! Mass Flux X-Dir
!   real(FVPRC), OPTIONAL, intent(in):: mfy(is:ie  ,js:je+1)  ! Mass Flux Y-Dir
!   real(FVPRC), OPTIONAL, intent(in):: mass(isd:ied,jsd:jed)
!   real(FVPRC), OPTIONAL, intent(in):: damp_c
!   integer, OPTIONAL, intent(in):: nord
!! Local:
!   real(FVPRC)   q(isd:ied,jsd:jed)
!   q = CNVT(qIN)
!   call fv_tp_2d(q, crx, cry, npx, npy, hord, fx, fy, xfx, yfx,  &
!                     ra_x, ra_y, mfx=mfx, mfy=mfy, mass=mass, nord=nord, damp_c=damp_c)
! end subroutine fv_tp_2d_r8
!#endif

 subroutine fv_tp_2d_r4(q, crx, cry, npx, npy, hord, fx, fy, xfx, yfx,  &
                     ra_x, ra_y, mfx, mfy, mass, nord, damp_c)
   integer, intent(in):: npx, npy
   integer, intent(in)::hord

   real(FVPRC), intent(in)::  crx(is:ie+1,jsd:jed)  !
   real(FVPRC), intent(in)::  xfx(is:ie+1,jsd:jed)  !
   real(FVPRC), intent(in)::  cry(isd:ied,js:je+1 )  !
   real(FVPRC), intent(in)::  yfx(isd:ied,js:je+1 )  !
   real(FVPRC), intent(in):: ra_x(is:ie,jsd:jed)
   real(FVPRC), intent(in):: ra_y(isd:ied,js:je)
   real(FVPRC), intent(inout):: q(isd:ied,jsd:jed)  ! transported scalar
   real(FVPRC), intent(out)::fx(is:ie+1 ,js:je)    ! Flux in x ( E )
   real(FVPRC), intent(out)::fy(is:ie,   js:je+1 )    ! Flux in y ( N )
! optional Arguments:
   real(FVPRC), OPTIONAL, intent(in):: mfx(is:ie+1,js:je  )  ! Mass Flux X-Dir
   real(FVPRC), OPTIONAL, intent(in):: mfy(is:ie  ,js:je+1)  ! Mass Flux Y-Dir
   real(FVPRC), OPTIONAL, intent(in):: mass(isd:ied,jsd:jed)
   real(FVPRC), OPTIONAL, intent(in):: damp_c
   integer, OPTIONAL, intent(in):: nord
! Local:
   integer ord_in, ord_ou
   real(FVPRC) q_i(isd:ied,js:je)
   real(FVPRC) q_j(is:ie,jsd:jed)
   real(FVPRC)   fx1(is :ie+1,js:je  )
   real(FVPRC)   fy1(is :ie  ,js:je+1)
   real(FVPRC)   fx2(is:ie+1,jsd:jed)
   real(FVPRC)   fy2(isd:ied,js:je+1)
   real(FVPRC)   fyy(isd:ied,js:je+1)
   real(FVPRC)   fxx(is:ie+1)
   real(FVPRC)   damp
   integer i, j

   ord_in = hord
   ord_ou = hord

   call copy_corners(q, npx, npy, 2)

   if ( ord_in < 0 ) then
      call yppm0(fy2, q, cry, ord_in, isd, ied, js, je, npx, npy)
   else
      call   ytp(fy2, q, cry, ord_in, isd, ied, js, je, npx, npy)
   endif

   do j=js,je+1
      do i=isd,ied
      fyy(i,j) = yfx(i,j) * fy2(i,j) 
      enddo
   enddo
   do j=js,je
      do i=isd,ied
         q_i(i,j) = (q(i,j)*area(i,j) + fyy(i,j)-fyy(i,j+1))/ra_y(i,j)
      enddo
  enddo
  if( ord_ou < 0 ) then
      call xppm0(fx1, q_i, crx(is,js), ord_ou, is, ie, js, je, npx, npy)
  else
      call   xtp(fx1, q_i, crx(is,js), ord_ou, is, ie, js, je, npx, npy)
  endif

  call copy_corners(q, npx, npy, 1)
  if( ord_in < 0 ) then
      call xppm0(fx2, q, crx, ord_in, is, ie, jsd, jed, npx, npy)
  else
      call   xtp(fx2, q, crx, ord_in, is, ie, jsd, jed, npx, npy)
  endif

  do j=jsd,jed
     do i=is,ie+1
        fxx(i) =  xfx(i,j) * fx2(i,j)
     enddo
     do i=is,ie
        q_j(i,j) = (q(i,j)*area(i,j) + fxx(i)-fxx(i+1))/ra_x(i,j)
     enddo
  enddo
  if ( ord_ou < 0 ) then
      call yppm0(fy1, q_j, cry, ord_ou, is, ie, js, je, npx, npy)
  else
      call   ytp(fy1, q_j, cry, ord_ou, is, ie, js, je, npx, npy)
  endif

!----------------
! Flux averaging:
!----------------

   if ( present(mfx) .and. present(mfy) ) then
!---------------------------------
! For transport of pt and tracers
!---------------------------------
      do j=js,je
         do i=is,ie+1
            fx(i,j) = 0.5*(fx1(i,j) + fx2(i,j)) * mfx(i,j)
         enddo
      enddo
      do j=js,je+1
         do i=is,ie
            fy(i,j) = 0.5*(fy1(i,j) + fy2(i,j)) * mfy(i,j)
         enddo
      enddo
      if ( present(nord) .and. present(damp_c) .and. present(mass) ) then
        if ( damp_c > 1.e-4 ) then
           damp = (damp_c * da_min)**(nord+1)
           call deln_flux( nord, npx, npy, damp, q, fx, fy, mass )
        endif
      endif
   else
!---------------------------------
! For transport of delp, vorticity
!---------------------------------
      do j=js,je
         do i=is,ie+1
            fx(i,j) = 0.5*(fx1(i,j) + fx2(i,j)) * xfx(i,j)
         enddo
      enddo
      do j=js,je+1
         do i=is,ie
            fy(i,j) = 0.5*(fy1(i,j) + fy2(i,j)) * yfx(i,j)
         enddo
      enddo
      if ( present(nord) .and. present(damp_c) ) then
           if ( damp_c > 1.E-4 ) then
                damp = (damp_c * da_min)**(nord+1)
                call deln_flux( nord, npx, npy, damp, q, fx, fy )
           endif
      endif
   endif

  

 end subroutine fv_tp_2d_r4


 subroutine copy_corners_r4(q, npx, npy, dir)
 integer, intent(in):: npx, npy, dir
 real(REAL4), intent(inout):: q(isd:ied,jsd:jed)
 integer  i,j

 if ( dir == 1 ) then
! XDir:
    if ( sw_corner ) then
         do j=1-ng,0
            do i=1-ng,0
               q(i,j) = q(j,1-i)
            enddo
         enddo
    endif
    if ( se_corner ) then
         do j=1-ng,0
            do i=npx,npx+ng-1
               q(i,j) = q(npy-j,i-npx+1)
            enddo
         enddo
    endif
    if ( ne_corner ) then
         do j=npy,npy+ng-1
            do i=npx,npx+ng-1
               q(i,j) = q(j,2*npx-1-i)
            enddo
         enddo
    endif
    if ( nw_corner ) then
         do j=npy,npy+ng-1
            do i=1-ng,0
               q(i,j) = q(npy-j,i-1+npx)
            enddo
         enddo
    endif

 elseif ( dir == 2 ) then
! YDir:

    if ( sw_corner ) then
         do j=1-ng,0
            do i=1-ng,0
               q(i,j) = q(1-j,i)
            enddo
         enddo
    endif
    if ( se_corner ) then
         do j=1-ng,0
            do i=npx,npx+ng-1
               q(i,j) = q(npy+j-1,npx-i)
            enddo
         enddo
    endif
    if ( ne_corner ) then
         do j=npy,npy+ng-1
            do i=npx,npx+ng-1
               q(i,j) = q(2*npy-1-j,i)
            enddo
         enddo
    endif
    if ( nw_corner ) then
         do j=npy,npy+ng-1
            do i=1-ng,0
               q(i,j) = q(j+1-npx,npy-i)
            enddo
         enddo
    endif

 endif
      
 end subroutine copy_corners_r4

 subroutine copy_corners_r8(q, npx, npy, dir)
 integer, intent(in):: npx, npy, dir
 real(REAL8), intent(inout):: q(isd:ied,jsd:jed)
 integer  i,j

 if ( dir == 1 ) then
! XDir:
    if ( sw_corner ) then
         do j=1-ng,0
            do i=1-ng,0
               q(i,j) = q(j,1-i)
            enddo
         enddo
    endif
    if ( se_corner ) then
         do j=1-ng,0
            do i=npx,npx+ng-1
               q(i,j) = q(npy-j,i-npx+1)
            enddo
         enddo
    endif
    if ( ne_corner ) then
         do j=npy,npy+ng-1
            do i=npx,npx+ng-1
               q(i,j) = q(j,2*npx-1-i)
            enddo
         enddo
    endif
    if ( nw_corner ) then
         do j=npy,npy+ng-1
            do i=1-ng,0
               q(i,j) = q(npy-j,i-1+npx)
            enddo
         enddo
    endif

 elseif ( dir == 2 ) then
! YDir:

    if ( sw_corner ) then
         do j=1-ng,0
            do i=1-ng,0
               q(i,j) = q(1-j,i)
            enddo
         enddo
    endif
    if ( se_corner ) then
         do j=1-ng,0
            do i=npx,npx+ng-1
               q(i,j) = q(npy+j-1,npx-i)
            enddo
         enddo
    endif
    if ( ne_corner ) then
         do j=npy,npy+ng-1
            do i=npx,npx+ng-1
               q(i,j) = q(2*npy-1-j,i)
            enddo
         enddo
    endif
    if ( nw_corner ) then
         do j=npy,npy+ng-1
            do i=1-ng,0
               q(i,j) = q(j+1-npx,npy-i)
            enddo
         enddo
    endif

 endif
      
 end subroutine copy_corners_r8

 subroutine xtp(fx,  q,  c, iord, ifirst, ilast, jfirst, jlast, npx, npy)
   integer, intent(IN):: ifirst, ilast   !  X-Dir strip
   integer, intent(IN):: jfirst, jlast   !  Y-Dir strip
   integer, intent(IN):: npx, npy
   integer, intent(IN):: iord
   real(FVPRC)   , intent(in):: c(is :ie+1, jfirst:jlast)      ! Courant numbers
   real(FVPRC)   , intent(in):: q(isd:ied,  jfirst:jlast)
   real(FVPRC)   , intent(out):: fx(ifirst:ilast+1,jfirst:jlast)           
! Local:
   real(FVPRC)   dm(is-2:ie+2)
   real(FVPRC)   x0L, x0R, x1
   integer i, j

   if (iord==1) then

      do j=jfirst,jlast
         do i=ifirst,ilast+1
           if ( c(i,j)>0. ) then
                fx(i,j) = q(i-1,j)
           else
                fx(i,j) = q(i,j)
           endif
         enddo
      enddo

   elseif (iord==333) then

     !Advection using third order scheme
     do j=jfirst,jlast
       do i=ifirst,ilast+1
         if ( c(i,j)>0. ) then
            fx(i,j) = (2.0*q(i,j)+5.0*q(i-1,j)-q(i-2,j) )/6.0 &
                       - 0.5*c(i,j)*(q(i,j)-q(i-1,j)) &
                       + (c(i,j)*c(i,j)/6.0)*(q(i,j)-2.0*q(i-1,j)+q(i-2,j) )
         else
            fx(i,j) = (2.0*q(i-1,j)+5.0*q(i,j)-q(i+1,j) )/6.0 &
                       - 0.5*c(i,j)*(q(i,j)-q(i-1,j)) &
                       + (c(i,j)*c(i,j)/6.0)*(q(i+1,j)-2.0*q(i,j)+q(i-1,j) )
         endif
       enddo
     enddo 

   else
      call fxppm(c, q, fx, iord, ifirst, ilast, jfirst, jlast, npx, npy)
   endif

 end subroutine xtp



 subroutine ytp(fy, q, c, jord, ifirst, ilast, jfirst, jlast, npx, npy)
 integer, intent(in) :: npx, npy
 integer, INTENT(IN) :: ifirst, ilast  !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast  !  Y-Dir strip
 integer, intent(in):: jord
 real(FVPRC), intent(in)::   q(ifirst:ilast,jfirst-ng:jlast+ng) 
 real(FVPRC), intent(in)::   c(isd:ied,js:je+1 )  ! Courant number
 real(FVPRC), intent(out):: fy(ifirst:ilast,jfirst:jlast+1)     !  Flux
! !LOCAL VARIABLES:
 real(FVPRC)   dm(ifirst:ilast,jfirst-2:jlast+2)
 real(FVPRC)   x0L, x0R, x1
 integer i, j

   if(jord==1) then

      do j=jfirst,jlast+1
         do i=ifirst,ilast
            if ( c(i,j)>0. ) then
                 fy(i,j) = q(i,j-1)
            else
                 fy(i,j) = q(i,j)
            endif
         enddo
      enddo

   elseif (jord==333) then

      !Advected using third order scheme
      do j=jfirst,jlast+1
         do i=ifirst,ilast
            if ( c(i,j)>0. ) then
               fy(i,j) = (2.0*q(i,j)+5.0*q(i,j-1)-q(i,j-2) )/6.0 &
                           - 0.5*c(i,j)*(q(i,j)-q(i,j-1)) &
                           + (c(i,j)*c(i,j)/6.0)*(q(i,j)-2.0*q(i,j-1)+q(i,j-2) ) 
            else
               fy(i,j) = (2.0*q(i,j-1)+5.0*q(i,j)-q(i,j+1) )/6.0 &
                           - 0.5*c(i,j)*(q(i,j)-q(i,j-1)) &
                           + (c(i,j)*c(i,j)/6.0)*(q(i,j+1)-2.0*q(i,j)+q(i,j-1) ) 
            endif
         enddo
      enddo

   else
      call fyppm(c, q, fy, jord, ifirst,ilast,jfirst,jlast, npx, npy)
   endif

 end subroutine ytp

!#ifdef XPPM_2D
! subroutine xppm0(flux, q, c, iord, ifirst, ilast, jfirst, jlast, npx, npy, bd, dxa, nested, grid_type)
! type(fv_grid_bounds_type), intent(IN) :: bd
! integer, INTENT(IN) :: ifirst, ilast               !  X-Dir strip
! integer, INTENT(IN) :: jfirst, jlast               !  Y-Dir strip
! integer, INTENT(IN) :: iord
! integer, INTENT(IN) :: npx, npy
! real   , INTENT(IN) :: q(ifirst-ng:ilast+ng,jfirst:jlast)
! real   , INTENT(IN) :: c(ifirst   :ilast+1 ,jfirst:jlast) ! Courant   N (like FLUX)
! real   , intent(IN) :: dxa(bd%isd:bd%ied,bd%jsd:bd%jed)
! logical, intent(IN) :: nested
! integer, intent(IN) :: grid_type
!! !OUTPUT PARAMETERS:
! real  , INTENT(OUT) :: flux(ifirst:ilast+1,jfirst:jlast) !  Flux
!! Local
! real, dimension(3*jfirst:3*jlast+2) :: q_tmp, bl_tmp, br_tmp
! real, dimension(ifirst-1:ilast+1,jfirst:jlast):: bl, br
! real  al(ifirst-1:ilast+2,jfirst:jlast)
! real  dm(ifirst-2:ilast+2,jfirst:jlast)
! real  dq(ifirst-3:ilast+2,jfirst:jlast)
! logical extm(ifirst-2:ilast+2,jfirst:jlast)
! integer i, j, ie3, is1, ie1
! real xt, pmp_1, lac_1, pmp_2, lac_2
!
! if ( .not. nested .and. grid_type<3 ) then
!    is1 = max(3,is-1);  ie3 = min(npx-2,ie+2)
!                        ie1 = min(npx-3,ie+1)
! else
!    is1 = is-1;         ie3 = ie+2
!                        ie1 = ie+1
! end if
!
! if ( abs(iord) < 8 ) then
!
!    ! ord = -5: linear scheme based on PPM 4th order interpolation
!    ! ord = -6: linear PPM with 2-delta filter
!    ! ord = -7: (-6) with additional Positive definite constraint:
!
!    do j=jfirst,jlast
!       do i=is1, ie3
!          al(i,j) = p1*(q(i-1,j)+q(i,j)) + p2*(q(i-2,j)+q(i+1,j))
!       enddo
!    enddo
!    if ( .not. nested .and. grid_type<3 ) then
!       if ( is==1 ) then
!          do j=jfirst,jlast
!             al(0,j) = c1*q(-2,j) + c2*q(-1,j) + c3*q(0,j)
!             al(1,j) = 0.5*(((2.*dxa(0,j)+dxa(-1,j))*q(0,j)-dxa(0,j)*q(-1,j))/(dxa(-1,j)+dxa(0,j)) &
!                  +      ((2.*dxa(1,j)+dxa( 2,j))*q(1,j)-dxa(1,j)*q(2,j))/(dxa(1, j)+dxa(2,j)))
!             al(2,j) = c3*q(1,j) + c2*q(2,j) +c1*q(3,j)
!          enddo
!          if(iord==-7) then
!             do j=jfirst,jlast
!                al(0,j) = max(0., al(0,j))
!                al(1,j) = max(0., al(1,j))
!                al(2,j) = max(0., al(2,j))
!             enddo
!          endif
!       endif
!       if ( (ie+1)==npx ) then
!          do j=jfirst,jlast
!             al(npx-1,j) = c1*q(npx-3,j) + c2*q(npx-2,j) + c3*q(npx-1,j)
!             al(npx,j) = 0.5*(((2.*dxa(npx-1,j)+dxa(npx-2,j))*q(npx-1,j)-dxa(npx-1,j)*q(npx-2,j))/(dxa(npx-2,j)+dxa(npx-1,j)) &
!                  +      ((2.*dxa(npx,  j)+dxa(npx+1,j))*q(npx,j)-dxa(npx,  j)*q(npx+1,j))/(dxa(npx,  j)+dxa(npx+1,j)))
!             al(npx+1,j) = c3*q(npx,j) + c2*q(npx+1,j) + c1*q(npx+2,j)
!          enddo
!          if(iord==-7) then
!             do j=jfirst,jlast
!                al(npx-1,j) = max(0., al(npx-1,j))
!                al(npx,  j) = max(0., al(npx  ,j))
!                al(npx+1,j) = max(0., al(npx+1,j))
!             enddo
!          endif
!       endif
!    endif
!
!    if ( iord==-5 ) then
!       do j=jfirst,jlast
!          do i=ifirst-1,ilast+1
!             bl(i,j) = al(i,j)   - q(i,j)
!             br(i,j) = al(i+1,j) - q(i,j)
!          enddo
!       enddo
!    else
!       do j=jfirst,jlast
!          do i=ifirst-3,ilast+2
!             dq(i,j) = q(i+1,j) - q(i,j)
!          enddo
!       enddo
!       do j=jfirst,jlast
!          do i=ifirst-2, ilast+2
!             if ( dq(i-1,j)*dq(i,j) > 0. ) then
!                extm(i,j) = .false.
!             else
!                extm(i,j) = .true.
!             endif
!          enddo
!       enddo
!       do j=jfirst,jlast
!          do i=ifirst-1,ilast+1
!             if ( extm(i-1,j).and.extm(i,j).and.extm(i+1,j) ) then
!                bl(i,j) = 0.
!                br(i,j) = 0.
!             else
!                bl(i,j) = al(i,j)   - q(i,j)
!                br(i,j) = al(i+1,j) - q(i,j)
!             endif
!          enddo
!       enddo
!       ! Additional positive definite constraint:
!       if(iord==-7) call pert_ppm(ilast-ifirst+3, q(ifirst-1,j), bl(ifirst-1,j), br(ifirst-1,j), 0)
!    endif
! else
!
!! Monotonic constraints:
!! ord = 8: PPM with Lin's PPM fast monotone constraint
!! ord = 10: PPM with Lin's modification of Huynh 2nd constraint
!! ord = 13: 10 plus positive definite constraint
!    do j=jfirst,jlast
!!DEC$ VECTOR ALWAYS
!       do i=is-2,ie+2
!          xt = 0.25*(q(i+1,j) - q(i-1,j))
!          dm(i,j) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
!               q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
!       enddo
!    enddo
!    do j=jfirst,jlast
!!DEC$ VECTOR ALWAYS
!       do i=is1,ie1+1
!          al(i,j) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm(i-1,j)-dm(i,j))
!       enddo
!    enddo
!    if ( iord==-8 ) then
!       do j=jfirst,jlast
!!DEC$ VECTOR ALWAYS
!          do i=is1, ie1
!             xt = 2.*dm(i,j)
!             bl(i,j) = -sign(min(abs(xt), abs(al(i,j  )-q(i,j))), xt)
!             br(i,j) =  sign(min(abs(xt), abs(al(i+1,j)-q(i,j))), xt)
!          enddo
!       enddo
!    else
!       do j=jfirst,jlast
!          do i=is1-2, ie1+1
!             dq(i,j) = q(i+1,j) - q(i,j)
!          enddo
!       enddo
!       do j=jfirst,jlast
!          do i=is1, ie1
!             bl(i,j) = al(i,j  ) - q(i,j)
!             br(i,j) = al(i+1,j) - q(i,j)
!             if ( abs(dm(i-1,j))+abs(dm(i,j))+abs(dm(i+1,j)) < near_zero ) then
!                bl(i,j) = 0.
!                br(i,j) = 0.
!             elseif( abs(3.*(bl(i,j)+br(i,j))) > abs(bl(i,j)-br(i,j)) ) then
!                pmp_1 = -(dq(i,j) + dq(i,j))
!                lac_1 = pmp_1 + 1.5*dq(i+1,j)
!                bl(i,j) = min( max(0., pmp_1, lac_1), max(bl(i,j), min(0., pmp_1, lac_1)) )
!                pmp_2 = dq(i-1,j) + dq(i-1,j)
!                lac_2 = pmp_2 - 1.5*dq(i-2,j)
!                br(i,j) = min( max(0., pmp_2, lac_2), max(br(i,j), min(0., pmp_2, lac_2)) )
!             endif
!          enddo
!       enddo
!    endif
!    ! Positive definite constraint:
!    if(iord==-9 .or. iord==-13) call pert_ppm(ie1-is1+1, q(is1,j), bl(is1,j), br(is1,j), 0)
!
!    if (.not. nested .and. grid_type<3) then
!       if ( is==1 ) then
!          do j=jfirst,jlast
!             bl_tmp(3*j) = s14*dm(-1,j) + s11*(q(-1,j)-q(0,j))
!
!             xt = 0.5*(((2.*dxa(0,j)+dxa(-1,j))*q(0,j)-dxa(0,j)*q(-1,j))/(dxa(-1,j)+dxa(0,j)) &
!                  +      ((2.*dxa(1,j)+dxa( 2,j))*q(1,j)-dxa(1,j)*q(2,j))/(dxa(1, j)+dxa(2,j)))
!             br_tmp(3*j) = xt - q(0,j)
!             bl_tmp(3*j+1) = xt - q(1,j)
!             xt = s15*q(1,j) + s11*q(2,j) - s14*dm(2,j)
!             br_tmp(3*j+1) = xt - q(1,j)
!             bl_tmp(3*j+2) = xt - q(2,j)
!             br_tmp(3*j+2) = al(3,j) - q(2,j)
!             q_tmp(j*3:j*3+2) = q(0:2,j)
!          enddo
!          call pert_ppm(3*(jlast-jfirst+1), q_tmp, bl_tmp, br_tmp, 1)
!          do j=jfirst,jlast
!             bl(0:2,j) = bl_tmp(j*3:j*3+2)
!             br(0:2,j) = br_tmp(j*3:j*3+2)
!          enddo
!      endif
!      if ( (ie+1)==npx ) then
!         do j=jfirst,jlast
!            bl_tmp(3*j) = al(npx-2,j) - q(npx-2,j)
!
!            xt = s15*q(npx-1,j) + s11*q(npx-2,j) + s14*dm(npx-2,j)
!            br_tmp(3*j) = xt - q(npx-2,j)
!            bl_tmp(3*j+1) = xt - q(npx-1,j)
!
!            xt = 0.5*(((2.*dxa(npx-1,j)+dxa(npx-2,j))*q(npx-1,j)-dxa(npx-1,j)*q(npx-2,j))/(dxa(npx-2,j)+dxa(npx-1,j)) &
!               +      ((2.*dxa(npx,  j)+dxa(npx+1,j))*q(npx,j  )-dxa(npx,  j)*q(npx+1,j))/(dxa(npx,  j)+dxa(npx+1,j)))
!            br_tmp(3*j+1) = xt - q(npx-1,j)
!            bl_tmp(3*j+2) = xt - q(npx,j  )
!
!            br_tmp(3*j+2) = s11*(q(npx+1,j)-q(npx,j)) - s14*dm(npx+1,j)
!            q_tmp(j*3:j*3+2) = q(npx-2:npx,j)
!         enddo
!         call pert_ppm(3*(jlast-jfirst+1), q_tmp, bl_tmp, br_tmp, 1)
!
!         do j=jfirst,jlast
!            bl(npx-2:npx,j) = bl_tmp(j*3:j*3+2)
!            br(npx-2:npx,j) = br_tmp(j*3:j*3+2)
!         enddo
!      endif
!    endif
!
!  endif
!
!  do j=jfirst,jlast
!!DEC$ VECTOR ALWAYS
!     do i=ifirst,ilast+1
!        if( c(i,j)>0. ) then
!           flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1,j)-c(i,j)*(bl(i-1,j)+br(i-1,j)))
!        else
!           flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
!        endif
!     enddo
!  enddo
!
! end subroutine xppm0
!#else
 subroutine xppm0(flux, q, c, iord, ifirst, ilast, jfirst, jlast, npx, npy)
 integer, INTENT(IN) :: ifirst, ilast               !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast               !  Y-Dir strip
 integer, INTENT(IN) :: iord
 integer, INTENT(IN) :: npx, npy
 real   , INTENT(IN) :: q(ifirst-ng:ilast+ng,jfirst:jlast)
 real   , INTENT(IN) :: c(ifirst   :ilast+1 ,jfirst:jlast) ! Courant   N (like FLUX)
! !OUTPUT PARAMETERS:
 real  , INTENT(OUT) :: flux(ifirst:ilast+1,jfirst:jlast) !  Flux
! Local
 real, dimension(ifirst-ng:ilast+ng):: q1
 real, dimension(ifirst-1:ilast+1):: bl, br
 real  al(ifirst-1:ilast+2)
 real  dm(ifirst-2:ilast+2)
 real  dq(ifirst-3:ilast+2)
 logical extm(ifirst-2:ilast+2)
 integer i, j, ie3, is1, ie1
 real xt, pmp_1, lac_1, pmp_2, lac_2

 if ( .not. nested .and. grid_type<3 ) then
    is1 = max(3,is-1);  ie3 = min(npx-2,ie+2)
                        ie1 = min(npx-3,ie+1)
 else
    is1 = is-1;         ie3 = ie+2
                        ie1 = ie+1
 end if

 do j=jfirst,jlast

    do i=ifirst-ng, ilast+ng
       q1(i) = q(i,j)
    enddo

  if ( abs(iord) < 8 ) then

! ord = -5: linear scheme based on PPM 4th order interpolation
! ord = -6: linear PPM with 2-delta filter
! ord = -7: (-6) with additional Positive definite constraint:

   do i=is1, ie3
      al(i) = p1*(q1(i-1)+q1(i)) + p2*(q1(i-2)+q1(i+1))
   enddo
   if ( .not. nested .and. grid_type<3 ) then
     if ( is==1 ) then
       al(0) = c1*q1(-2) + c2*q1(-1) + c3*q1(0)
       al(1) = 0.5*(((2.*dxa(0,j)+dxa(-1,j))*q1(0)-dxa(0,j)*q1(-1))/(dxa(-1,j)+dxa(0,j)) &
             +      ((2.*dxa(1,j)+dxa( 2,j))*q1(1)-dxa(1,j)*q1( 2))/(dxa(1, j)+dxa(2,j)))
       al(2) = c3*q1(1) + c2*q1(2) +c1*q1(3)
     endif
     if ( (ie+1)==npx ) then
       al(npx-1) = c1*q1(npx-3) + c2*q1(npx-2) + c3*q1(npx-1)
       al(npx) = 0.5*(((2.*dxa(npx-1,j)+dxa(npx-2,j))*q1(npx-1)-dxa(npx-1,j)*q1(npx-2))/(dxa(npx-2,j)+dxa(npx-1,j)) &
               +      ((2.*dxa(npx,  j)+dxa(npx+1,j))*q1(npx  )-dxa(npx,  j)*q1(npx+1))/(dxa(npx,  j)+dxa(npx+1,j)))
       al(npx+1) = c3*q1(npx) + c2*q1(npx+1) + c1*q1(npx+2)
     endif
   endif

   if ( iord==-5 ) then
      do i=ifirst-1,ilast+1
         bl(i) = al(i)   - q1(i)
         br(i) = al(i+1) - q1(i)
      enddo
   endif

  endif

  do i=ifirst,ilast+1
     if( c(i,j)>0. ) then
         flux(i,j) = q1(i-1) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
     else
         flux(i,j) = q1(i  ) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
     endif
  enddo

 enddo

 end subroutine xppm0

!#endif

 subroutine yppm0(flux, q, c, jord, ifirst, ilast, jfirst, jlast, npx, npy)
 integer, INTENT(IN) :: ifirst, ilast               !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast               !  Y-Dir strip
 integer, INTENT(IN) :: jord
 integer, INTENT(IN) :: npx, npy
 real   , INTENT(IN) :: q(ifirst:ilast,jfirst-ng:jlast+ng)
 real   , intent(in) :: c(isd:ied,js:je+1 )  ! Courant number
 real   , INTENT(OUT):: flux(ifirst:ilast,jfirst:jlast+1)   !  Flux
! Local:
 real:: dm(ifirst:ilast,jfirst-2:jlast+2)
 real:: al(ifirst:ilast,jfirst-1:jlast+2)
 real bl(ifirst:ilast,jfirst-1:jlast+1)
 real br(ifirst:ilast,jfirst-1:jlast+1)
 real dq(ifirst:ilast,jfirst-3:jlast+2)
 logical extm(ifirst:ilast,jfirst-2:jlast+2)
 real xt, pmp_1, lac_1, pmp_2, lac_2
 integer i, j, js1, je3, je1

   if ( .not.nested .and. grid_type < 3 ) then
! Cubed-sphere:
      js1 = max(3,js-1); je3 = min(npy-2,je+2)
                         je1 = min(npy-3,je+1)
   else
! Nested grid OR Doubly periodic domain:
      js1 = js-1;        je3 = je+2
                         je1 = je+1
   endif

if ( abs(jord) < 8 ) then

! ord = -5: linear scheme based on PPM 4th order interpolation
! ord = -6: linear PPM with 2-delta filter
! ord = -7: (-6) with additional Positive definite constraint:

   do j=js1, je3
      do i=ifirst,ilast
         al(i,j) = p1*(q(i,j-1)+q(i,j)) + p2*(q(i,j-2)+q(i,j+1))
      enddo
   enddo
   if ( .not. nested .and. grid_type<3 ) then
      if( js==1 ) then
        do i=ifirst,ilast
           al(i,0) = c1*q(i,-2) + c2*q(i,-1) + c3*q(i,0)
           al(i,1) = 0.5*(((2.*dya(i,0)+dya(i,-1))*q(i,0)-dya(i,0)*q(i,-1))/(dya(i,-1)+dya(i,0))   &
                   +      ((2.*dya(i,1)+dya(i,2))*q(i,1)-dya(i,1)*q(i,2))/(dya(i,1)+dya(i,2)))
           al(i,2) = c3*q(i,1) + c2*q(i,2) + c1*q(i,3)
        enddo
      endif
      if( (je+1)==npy ) then
        do i=ifirst,ilast
         al(i,npy-1) = c1*q(i,npy-3) + c2*q(i,npy-2) + c3*q(i,npy-1)
         al(i,npy) = 0.5*(((2.*dya(i,npy-1)+dya(i,npy-2))*q(i,npy-1)-dya(i,npy-1)*q(i,npy-2))/(dya(i,npy-2)+dya(i,npy-1))  &
                   +      ((2.*dya(i,npy)+dya(i,npy+1))*q(i,npy)-dya(i,npy)*q(i,npy+1))/(dya(i,npy)+dya(i,npy+1)))
         al(i,npy+1) = c3*q(i,npy) + c2*q(i,npy+1) + c1*q(i,npy+2)
        enddo
      endif
   endif

   if ( jord==-5 ) then

      do j=jfirst-1,jlast+1
         do i=ifirst,ilast
            bl(i,j) = al(i,j  ) - q(i,j)
            br(i,j) = al(i,j+1) - q(i,j)
         enddo
      enddo

   endif

endif

  do j=jfirst,jlast+1
     do i=ifirst,ilast
        if( c(i,j)>0. ) then
           flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
        else
           flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
        endif
     enddo
  enddo

 end subroutine yppm0


 subroutine fxppm(c, q, flux, iord, ifirst, ilast, jfirst, jlast, npx, npy)
! !INPUT PARAMETERS:
 integer, INTENT(IN) :: ifirst, ilast               !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast               !  Y-Dir strip
 integer, INTENT(IN) :: iord
 integer, INTENT(IN) :: npx, npy
 real(FVPRC)   , INTENT(IN) :: q(ifirst-ng:ilast+ng,jfirst:jlast)
 real(FVPRC)   , INTENT(IN) :: c(ifirst   :ilast+1 ,jfirst:jlast) ! Courant   N (like FLUX)
! !OUTPUT PARAMETERS:
 real(FVPRC)   , INTENT(OUT) :: flux(ifirst:ilast+1,jfirst:jlast) !  Flux
! Local
 logical extm(ifirst-2:ilast+2)
 real(FVPRC) dm1(ifirst-2:ilast+2)
 real(FVPRC)  al(ifirst-1:ilast+2)
 real(FVPRC)  bl(ifirst-1:ilast+1)
 real(FVPRC)  br(ifirst-1:ilast+1)
 real(FVPRC)  dq(ifirst-3:ilast+2)
 real(FVPRC) dl, dr, pmp, lac, ct, qe
 real(FVPRC) pmp_1, lac_1, pmp_2, lac_2
 real(FVPRC) xt, x1, x0, x0L, x0R
 integer i, j, is3, ie3, it

 x0 = big_number

 is3 = max(3,is-1);   ie3 = min(npx-3,ie+1)

 if (iord == 666) then

     !From NCEP dynamics version

     do j=jfirst,jlast

        do i=is3, ie3
           bl(i) = b5*q(i-2,j) + b4*q(i-1,j) + b3*q(i,j) + b2*q(i+1,j) + b1*q(i+2,j)
           br(i) = b1*q(i-2,j) + b2*q(i-1,j) + b3*q(i,j) + b4*q(i+1,j) + b5*q(i+2,j)
        enddo

!--------------
! fix the edges
!--------------
        if ( is==1 ) then
             x0L = 0.5*((2.*dxa(0,j)+dxa(-1,j))*(q(0,j))   &
                - dxa(0,j)*(q(-1,j)))/ ( dxa(0,j)+dxa(-1,j))
             x0R = 0.5*((2.*dxa(1,j)+dxa(2,j))*(q(1,j))   &
                - dxa(1,j)*(q(2,j)))/ ( dxa(1,j)+dxa(2,j))
             br(2) = p1*(q(2,j)+q(3,j)) + p2*(q(1,j)+q(4,j)) - q(2,j)
             xt = x0L + x0R
             bl(1) = xt - q(1,j)
             br(0) = xt - q(0,j)

             xt = c1*q(-2,j) + c2*q(-1,j) + c3*q(0,j)
             bl(0) = xt - q(0,j)

             xt = c3*q(1,j) + c2*q(2,j) +c1*q(3,j)
             br(1) = xt - q(1,j)
             bl(2) = xt - q(2,j)
        endif

        if ( (ie+1)==npx ) then
           x0L = 0.5*( (2.*dxa(npx-1,j)+dxa(npx-2,j))*(q(npx-1,j))   &
                - dxa(npx-1,j)*(q(npx-2,j)))/( dxa(npx-1,j)+dxa(npx-2,j))
           x0R = 0.5*( (2.*dxa(npx,j)+dxa(npx+1,j))*(q(npx,j))   &
                - dxa(npx,j)*(q(npx+1,j)))/( dxa(npx,j)+dxa(npx+1,j))
             bl(npx-2) = p1*(q(npx-2,j)+q(npx-3,j)) + p2*(q(npx-4,j)+q(npx-1,j)) - q(npx-2,j)
!             xt = x0L*sin_sg(npx-1,j,3) + x0R*sin_sg(npx,j,1)
!             xt = 2.*xt / (sin_sg(npx,j,1) + sin_sg(npx-1,j,3))
             xt = x0L + x0R

             br(npx-1) = xt - q(npx-1,j)
             bl(npx  ) = xt - q(npx  ,j)

             xt = c3*q(npx,j) + c2*q(npx+1,j) + c1*q(npx+2,j)
             br(npx) = xt - q(npx,j)

             xt = c1*q(npx-3,j) + c2*q(npx-2,j) + c3*q(npx-1,j)
             br(npx-2) = xt - q(npx-2,j)
             bl(npx-1) = xt - q(npx-1,j)
        endif

        do i=ifirst,ilast+1
           if(c(i,j)>0.) then
              flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
           else
              flux(i,j) = q(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
           endif
        enddo
     enddo

 endif

 end subroutine fxppm



 subroutine fyppm(c,  q,  flux, jord, ifirst, ilast, jfirst, jlast, npx, npy)
 integer, INTENT(IN) :: ifirst, ilast               !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast               !  Y-Dir strip
 integer, INTENT(IN) :: jord
 integer, INTENT(IN) :: npx, npy
 real(FVPRC)   , INTENT(IN) :: q(ifirst:ilast,jfirst-ng:jlast+ng)
 real(FVPRC)   , intent(in) :: c(isd:ied,js:je+1 )  ! Courant number
 real(FVPRC)   , INTENT(OUT):: flux(ifirst:ilast,jfirst:jlast+1)   !  Flux
! real(FVPRC)   , INTENT(OUT)::   dm(ifirst:ilast,jfirst-2:jlast+2)
! Local:
 logical extm(ifirst:ilast,jfirst-2:jlast+2)
 real(FVPRC) al(ifirst:ilast,jfirst-1:jlast+2)
 real(FVPRC) bl(ifirst:ilast,jfirst-1:jlast+1)
 real(FVPRC) br(ifirst:ilast,jfirst-1:jlast+1)
 real(FVPRC) dq(ifirst:ilast,jfirst-3:jlast+2)
 real(FVPRC) dl, dr, pmp, lac, ct, qe
 real(FVPRC) pmp_1, lac_1, pmp_2, lac_2
 real(FVPRC) xt, x0, x1, x0L, x0R
 integer i, j, js3, je3, jt

 if( jord==666 ) then

   do j=max(3,js-1),min(npy-3,je+1)
      do i=ifirst,ilast
         bl(i,j) = b5*q(i,j-2) + b4*q(i,j-1) + b3*q(i,j) + b2*q(i,j+1) + b1*q(i,j+2)
         br(i,j) = b1*q(i,j-2) + b2*q(i,j-1) + b3*q(i,j) + b4*q(i,j+1) + b5*q(i,j+2)
      enddo
   enddo

   if( js==1) then
         do i=ifirst,ilast
!           br(i,2) = al(i,3) - q(i,2)
            br(i,2) = p1*(q(i,2)+q(i,3)) + p2*(q(i,1)+q(i,4)) - q(i,2)
            x0L = 0.5*((2.*dya(i,0)+dya(i,-1))*(q(i,0))   &
               -dya(i,0)*(q(i,-1))) / ( dya(i,0)+dya(i,-1) )
            x0R = 0.5*((2.*dya(i,1)+dya(i,2))*(q(i,1))   &
               -dya(i,1)*(q(i,2))) / ( dya(i,1)+dya(i,2) )
            xt = x0L + x0R
!            xt = ( x0L*sin_sg(i,0,4) + x0R*sin_sg(i,1,2) ) 
!            xt = 2.*xt / ( sin_sg(i,0,4) + sin_sg(i,1,2) )
            bl(i,1) = xt - q(i,1)
            br(i,0) = xt - q(i,0)

!           xt = s14*0.25*(q(i,0)-q(i,-2)) - s11*(q(i,0)-q(i,-1)) + q(i,0)
            xt = c1*q(i,-2) + c2*q(i,-1) + c3*q(i,0)
            bl(i,0) = xt - q(i,0)

!           xt = s15*q(i,1) + s11*q(i,2) - s14*0.25*(q(i,3)-q(i,1))
            xt = c3*q(i,1) + c2*q(i,2) + c1*q(i,3)
            br(i,1) = xt - q(i,1)
            bl(i,2) = xt - q(i,2)
         enddo
   endif

   if( (je+1)==npy) then
         do i=ifirst,ilast
!           bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)
            bl(i,npy-2) = p1*(q(i,npy-3)+q(i,npy-2)) + p2*(q(i,npy-4)+q(i,npy-1)) - q(i,npy-2)
            x0L = 0.5*((2.*dya(i,npy-1)+dya(i,npy-2))*(q(i,npy-1))  &
               -dya(i,npy-1)*(q(i,npy-2)))/(dya(i,npy-1)+dya(i,npy-2))
            x0R = 0.5*((2.*dya(i,npy)+dya(i,npy+1))*(q(i,npy))  &
               -dya(i,npy)*(q(i,npy+1)))/(dya(i,npy)+dya(i,npy+1))
            xt = x0L + x0R
!            xt = x0L*sin_sg(i,npy-1,4) + x0R*sin_sg(i,npy,2)
!            xt = 2.*xt /( sin_sg(i,npy-1,4) + sin_sg(i,npy,2) )
            br(i,npy-1) = xt - q(i,npy-1)
            bl(i,npy  ) = xt - q(i,npy)

!           xt = s11*(q(i,npy+1)-q(i,npy)) - s14*0.25*(q(i,npy+2)-q(i,npy)) + q(i,npy)
            xt = c3*q(i,npy) + c2*q(i,npy+1) + c1*q(i,npy+2)
            br(i,npy) = xt - q(i,npy)

!           xt = s15*q(i,npy-1) + s11*q(i,npy-2) + s14*0.25*(q(i,npy-1)-q(i,npy-3))
            xt = c1*q(i,npy-3) + c2*q(i,npy-2) + c3*q(i,npy-1)
            br(i,npy-2) = xt - q(i,npy-2)
            bl(i,npy-1) = xt - q(i,npy-1)
         enddo
   endif

   do j=jfirst,jlast+1
      do i=ifirst,ilast
         if(c(i,j)>0.) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo

 endif

 end subroutine fyppm

                       

 subroutine mp_ghost_ew(im, jm, km, nq, ifirst, ilast, jfirst, jlast, &
                              kfirst, klast, ng_w, ng_e, ng_s, ng_n, q_ghst, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: ifirst, ilast
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_e      ! eastern  zones to ghost
      integer, intent(in):: ng_w      ! western  zones to ghost
      integer, intent(in):: ng_s      ! southern zones to ghost
      integer, intent(in):: ng_n      ! northern zones to ghost
      real(FVPRC), intent(inout):: q_ghst(ifirst-ng_w:ilast+ng_e,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
      real(FVPRC), optional, intent(in):: q(ifirst:ilast,jfirst:jlast,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Ghost 4d east/west 
!
! !REVISION HISTORY:
!    2005.08.22   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: i,j,k,n

      if (present(q)) then
         q_ghst(ifirst:ilast,jfirst:jlast,kfirst:klast,1:nq) = &
              q(ifirst:ilast,jfirst:jlast,kfirst:klast,1:nq)
      endif

!      Assume Periodicity in X-dir and not overlapping
      do n=1,nq
         do k=kfirst,klast
            do j=jfirst-ng_s,jlast+ng_n
               do i=1, ng_w
                  q_ghst(ifirst-i,j,k,n) = q_ghst(ilast-i+1,j,k,n)
               enddo
               do i=1, ng_e
                  q_ghst(ilast+i,j,k,n) = q_ghst(ifirst+i-1,j,k,n)
               enddo
            enddo
         enddo
      enddo

 end subroutine mp_ghost_ew



 subroutine pert_ppm(im, a0, al, ar, iv)
 integer, intent(in):: im
 integer, intent(in):: iv
 real(FVPRC), intent(in)   :: a0(im)
 real(FVPRC), intent(inout):: al(im), ar(im)
! Local:
 real(FVPRC) a4, da1, da2, a6da, fmin
 integer i
 real(FVPRC), parameter:: r12 = 1./12.

!-----------------------------------
! Optimized PPM in perturbation form:
!-----------------------------------

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
     if ( a0(i) <= 0. ) then
          al(i) = 0.
          ar(i) = 0.
     else
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

 end subroutine pert_ppm


 subroutine deln_flux( nord, npx, npy, damp, q, fx, fy, mass )
! Del-n damping for the cell-mean values (A grid)
!------------------
! nord = 0:   del-2
! nord = 1:   del-4
! nord = 2:   del-6
! nord = 3:   del-8 --> requires more ghosting than current
!------------------
   integer, intent(in):: nord            ! del-n
   integer, intent(in):: npx, npy
   real(FVPRC), intent(in):: damp
   real(FVPRC), intent(in):: q(is-ng:ie+ng, js-ng:je+ng)  ! q ghosted on input
   real(FVPRC), optional, intent(in):: mass(isd:ied, jsd:jed)  ! q ghosted on input
! diffusive fluxes:
   real(FVPRC), intent(inout):: fx(is:ie+1,js:je), fy(is:ie,js:je+1)
! local:
   real(FVPRC) fx2(isd:ied+1,jsd:jed), fy2(isd:ied,jsd:jed+1)
   real(FVPRC) d2(isd:ied,jsd:jed)
   real(FVPRC) damp2
   integer i,j, n, nt, i1, i2, j1, j2

   i1 = is-1-nord;    i2 = ie+1+nord
   j1 = js-1-nord;    j2 = je+1+nord

   if ( .not. present(mass) ) then
     do j=j1, j2
        do i=i1,i2
           d2(i,j) = damp*q(i,j)
        enddo
     enddo
   else
     do j=j1, j2
        do i=i1,i2
           d2(i,j) = q(i,j)
        enddo
     enddo
   endif

   if( nord>0 ) call copy_corners(d2, npx, npy, 1)

   do j=js-nord,je+nord
      do i=is-nord,ie+nord+1
         fx2(i,j) = 0.5*(sin_sg(i-1,j,3)+sin_sg(i,j,1))*dy(i,j)*(d2(i-1,j)-d2(i,j))*rdxc(i,j)
      enddo
   enddo

   if( nord>0 ) call copy_corners(d2, npx, npy, 2)
   do j=js-nord,je+nord+1
         do i=is-nord,ie+nord
            fy2(i,j) = 0.5*(sin_sg(i,j-1,4)+sin_sg(i,j,2))*dx(i,j)*(d2(i,j-1)-d2(i,j))*rdyc(i,j)
         enddo
   enddo

   if ( nord>0 ) then

!----------
! high-order
!----------

   do n=1, nord

      nt = nord-n

      do j=js-nt-1,je+nt+1
         do i=is-nt-1,ie+nt+1
            d2(i,j) = (fx2(i,j)-fx2(i+1,j)+fy2(i,j)-fy2(i,j+1))*rarea(i,j)
         enddo
      enddo

      call copy_corners(d2, npx, npy, 1)
      do j=js-nt,je+nt
         do i=is-nt,ie+nt+1
            fx2(i,j) = 0.5*(sin_sg(i-1,j,3)+sin_sg(i,j,1))*dy(i,j)*(d2(i,j)-d2(i-1,j))*rdxc(i,j)
         enddo
      enddo

      call copy_corners(d2, npx, npy, 2)
      do j=js-nt,je+nt+1
            do i=is-nt,ie+nt
               fy2(i,j) = dx(i,j)*(d2(i,j)-d2(i,j-1))*rdyc(i,j) &
                         *0.5*(sin_sg(i,j-1,4) + sin_sg(i,j,2) )
            enddo
      enddo
   enddo

   endif

!---------------------------------------------
! Add the diffusive fluxes to the flux arrays:
!---------------------------------------------

   if ( present(mass) ) then
! Apply mass weighting to diffusive fluxes:
        damp2 = 0.5*damp
        do j=js,je
           do i=is,ie+1
              fx(i,j) = fx(i,j) + damp2*(mass(i-1,j)+mass(i,j))*fx2(i,j)
           enddo
        enddo
        do j=js,je+1
           do i=is,ie
              fy(i,j) = fy(i,j) + damp2*(mass(i,j-1)+mass(i,j))*fy2(i,j)
           enddo
        enddo
   else
        do j=js,je
           do i=is,ie+1
              fx(i,j) = fx(i,j) + fx2(i,j)
           enddo
        enddo
        do j=js,je+1
           do i=is,ie
              fy(i,j) = fy(i,j) + fy2(i,j)
           enddo
        enddo
   endif

 end subroutine deln_flux


end module tp_core_mod
