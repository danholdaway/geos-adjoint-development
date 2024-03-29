module a2b_edge_mod

  use fv_grid_utils_mod, only: edge_w, edge_e, edge_s, edge_n, sw_corner, se_corner,  &
                               nw_corner, ne_corner, van2
  use fv_grid_tools_mod, only: dxa, dya, grid_type
  use fv_mp_mod,         only: gid
  implicit none

  real, parameter:: r3 = 1./3.
!----------------------------
! 4-pt Lagrange interpolation
!----------------------------
  real, parameter:: a1 =  0.5625  !  9/16
  real, parameter:: a2 = -0.0625  ! -1/16
!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: b1 =  7./12.     ! 0.58333333
  real, parameter:: b2 = -1./12.

  private
  public :: a2b_ord2, a2b_ord4

      INTERFACE a2b_ord2
        MODULE PROCEDURE a2b_ord2_r4
        MODULE PROCEDURE a2b_ord2_r8
      END INTERFACE

      INTERFACE a2b_ord4
        MODULE PROCEDURE a2b_ord4_r4
        MODULE PROCEDURE a2b_ord4_r8
      END INTERFACE

contains

  subroutine a2b_ord4_r8(qin, qout, npx, npy, is, ie, js, je, ng, replace)
  integer, intent(IN):: npx, npy, is, ie, js, je, ng
  real*8, intent(INOUT)::  qin(is-ng:ie+ng,js-ng:je+ng)   ! A-grid field
  real  , intent(INOUT):: qout(is-ng:ie+ng,js-ng:je+ng)   ! Output  B-grid field
  logical, optional, intent(IN):: replace
! local: compact 4-pt cubic
  real  , parameter:: c1 =  2./3.
  real  , parameter:: c2 = -1./6.
! Parabolic spline
! real  , parameter:: c1 =  0.75
! real  , parameter:: c2 = -0.25
! 6-pt corner interpolation:
  real  , parameter:: d1 =  0.375                   !   0.5
  real  , parameter:: d2 = -1./24.                  !  -1./6.


  real   qx(is:ie+1,js-ng:je+ng)
  real   qy(is-ng:ie+ng,js:je+1)
  real   qxx(is-ng:ie+ng,js-ng:je+ng)
  real   qyy(is-ng:ie+ng,js-ng:je+ng)
  real   fx(is:ie), fy(is-2:ie+2,js:je)
  real   gratio, qt(npy)
  integer :: i, j, is1, js1, is2, js2, ie1, je1
  integer :: im2, jm2

    im2 = (npx-1)/2
    jm2 = (npy-1)/2

  if (grid_type < 3) then

    is1 = max(1,is-1)
    js1 = max(1,js-1)
    is2 = max(2,is)
    js2 = max(2,js)

    ie1 = min(npx-1,ie+1)
    je1 = min(npy-1,je+1)

! Corners:
#ifdef USE_3PT
    if ( sw_corner ) qout(1,    1) = r3*(qin(1,        1)+qin(1,      0)+qin(0,      1))
    if ( se_corner ) qout(npx,  1) = r3*(qin(npx-1,    1)+qin(npx-1,  0)+qin(npx,    1))
    if ( ne_corner ) qout(npx,npy) = r3*(qin(npx-1,npy-1)+qin(npx,npy-1)+qin(npx-1,npy))
    if ( nw_corner ) qout(1,  npy) = r3*(qin(1,    npy-1)+qin(0,  npy-1)+qin(1,    npy))
#else
! 6-point formular:
    if ( sw_corner ) then
        qout(1,1) = d1*(qin(1, 0) + qin( 0,1) + qin(1,1)) +  &
                    d2*(qin(2,-1) + qin(-1,2) + qin(2,2))
    endif
    if ( se_corner ) then
        qout(npx,1) = d1*(qin(npx-1, 0) + qin(npx-1,1) + qin(npx,  1)) +  &
                      d2*(qin(npx-2,-1) + qin(npx-2,2) + qin(npx+1,2))
    endif
    if ( ne_corner ) then
        qout(npx,npy) = d1*(qin(npx-1,npy-1) + qin(npx,  npy-1) + qin(npx-1,npy)) +  &
                        d2*(qin(npx-2,npy-2) + qin(npx+1,npy-2) + qin(npx-2,npy+1))
    endif
    if ( nw_corner ) then
        qout(1,npy) = d1*(qin( 0,npy-1) + qin(1,npy-1) + qin(1,npy)) +   &
                      d2*(qin(-1,npy-2) + qin(2,npy-2) + qin(2,npy+1))
    endif
#endif

!------------
! X-Interior:
!------------
    do j=max(1,js-2),min(npy-1,je+2)
       do i=max(3,is), min(npx-2,ie+1)
          qx(i,j) = b2*(qin(i-2,j)+qin(i+1,j)) + b1*(qin(i-1,j)+qin(i,j))
       enddo
    enddo

! West Edges:
    if ( is==1 ) then

       do j=max(1,js-2),min(npy-1,je+2)
           gratio = dxa(2,j) / dxa(1,j)
          qx(1,j) = 0.5*((2.+gratio)*(qin(0,j)+qin(1,j))    &
                  - (qin(-1,j)+qin(2,j))) / (1.+gratio)
#ifdef TEST2
! Note: Caused noises in test_case-5 for large n_split
          qx(2,j) = (2.*gratio*(gratio+1.)*qin(1,j)+qin(2,j) -     &
                     gratio*(gratio+0.5)*qx(1,j))/(1.+gratio*(gratio+1.5))
#else
          qx(2,j) = (3.*(gratio*qin(1,j)+qin(2,j)) - (gratio*qx(1,j)+qx(3,j)))/(2.+2.*gratio)
#endif
       enddo

       do j=max(3,js),min(npy-2,je+1)
          qout(1,j) = a2*(qx(1,j-2)+qx(1,j+1)) + a1*(qx(1,j-1)+qx(1,j))
       enddo

       if( js==1 )     qout(1,    2) = c1*(qx(1,1)+qx(1,2))         + c2*(qout(1,1)+qout(1,3))
       if((je+1)==npy) qout(1,npy-1) = c1*(qx(1,npy-2)+qx(1,npy-1)) + c2*(qout(1,npy-2)+qout(1,npy))
    endif

! East Edges:
    if ( (ie+1)==npx ) then

       do j=max(1,js-2),min(npy-1,je+2)
               gratio = dxa(npx-2,j) / dxa(npx-1,j)
          qx(npx  ,j) = 0.5*((2.+gratio)*(qin(npx-1,j)+qin(npx,j))   &
                        - (qin(npx-2,j)+qin(npx+1,j))) / (1.+gratio )
#ifdef TEST2
          qx(npx-1,j) = (2.*gratio*(gratio+1.)*qin(npx-1,j)+qin(npx-2,j) -  &
                         gratio*(gratio+0.5)*qx(npx,j))/(1.+gratio*(gratio+1.5))
#else
          qx(npx-1,j) = (3.*(qin(npx-2,j)+gratio*qin(npx-1,j)) - (gratio*qx(npx,j)+qx(npx-2,j)))/(2.+2.*gratio)
#endif
       enddo

       do j=max(3,js),min(npy-2,je+1)
          qout(npx,j) = a2*(qx(npx,j-2)+qx(npx,j+1)) + a1*(qx(npx,j-1)+qx(npx,j))
       enddo

       if(js==1) qout(npx,2) = c1*(qx(npx,1)+qx(npx,2))+c2*(qout(npx,1)+qout(npx,3))
       if((je+1)==npy)   qout(npx,npy-1) =             &
                         c1*(qx(npx,npy-2)+qx(npx,npy-1))+c2*(qout(npx,npy-2)+qout(npx,npy))
    endif

!------------
! Y-Interior:
!------------
    do j=max(3,js),min(npy-2,je+1)
       do i=max(1,is-2), min(npx-1,ie+2)
          qy(i,j) = b2*(qin(i,j-2)+qin(i,j+1)) + b1*(qin(i,j-1)+qin(i,j))
       enddo
    enddo

! South Edges:
    if ( js==1 ) then

       do i=max(1,is-2),min(npx-1,ie+2)
           gratio = dya(i,2) / dya(i,1)
          qy(i,1) = 0.5*((2.+gratio)*(qin(i,0)+qin(i,1))   &
                  - (qin(i,-1)+qin(i,2))) / (1.+gratio )
#ifdef TEST2
          qy(i,2) = (2.*gratio*(gratio+1.)*qin(i,1)+qin(i,2) -     &
                     gratio*(gratio+0.5)*qy(i,1))/(1.+gratio*(gratio+1.5))
#else
          qy(i,2) = (3.*(gratio*qin(i,1)+qin(i,2)) - (gratio*qy(i,1)+qy(i,3)))/(2.+2.*gratio)
#endif
       enddo

       do i=max(3,is),min(npx-2,ie+1)
          qout(i,1) = a2*(qy(i-2,1)+qy(i+1,1)) + a1*(qy(i-1,1)+qy(i,1))
       enddo

       if( is==1 )    qout(2,1) = c1*(qy(1,1)+qy(2,1))+c2*(qout(1,1)+qout(3,1))
       if((ie+1)==npx) qout(npx-1,1) = c1*(qy(npx-2,1)+qy(npx-1,1))+c2*(qout(npx-2,1)+qout(npx,1))
    endif


! North Edges:
    if ( (je+1)==npy ) then
       do i=max(1,is-2),min(npx-1,ie+2)
               gratio = dya(i,npy-2) / dya(i,npy-1)
          qy(i,npy  ) = 0.5*((2.+gratio)*(qin(i,npy-1)+qin(i,npy))  &
                      - (qin(i,npy-2)+qin(i,npy+1))) / (1.+gratio)
#ifdef TEST2
          qy(i,npy-1) = (2.*gratio*(gratio+1.)*qin(i,npy-1)+qin(i,npy-2) - &
                         gratio*(gratio+0.5)*qy(i,npy))/(1.+gratio*(gratio+1.5))
#else
          qy(i,npy-1) = (3.*(qin(i,npy-2)+gratio*qin(i,npy-1)) - (gratio*qy(i,npy)+qy(i,npy-2)))/(2.+2.*gratio)
#endif
       enddo

       do i=max(3,is),min(npx-2,ie+1)
          qout(i,npy) = a2*(qy(i-2,npy)+qy(i+1,npy)) + a1*(qy(i-1,npy)+qy(i,npy))
       enddo

       if( is==1 )  qout(2,npy) = c1*(qy(1,npy)+qy(2,npy))+c2*(qout(1,npy)+qout(3,npy))
       if((ie+1)==npx) qout(npx-1,npy) = c1*(qy(npx-2,npy)+qy(npx-1,npy))+c2*(qout(npx-2,npy)+qout(npx,npy))
    endif
    
    do j=max(3,js),min(npy-2,je+1)
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,j) = a2*(qx(i,j-2)+qx(i,j+1)) + a1*(qx(i,j-1)+qx(i,j))
       enddo
    enddo

    if ( js==1 ) then
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,2) = c1*(qx(i,1)+qx(i,2))+c2*(qout(i,1)+qxx(i,3))
       enddo
    endif
    if ( (je+1)==npy ) then
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,npy-1) = c1*(qx(i,npy-2)+qx(i,npy-1))+c2*(qout(i,npy)+qxx(i,npy-2))
       enddo
    endif

    
    do j=max(2,js),min(npy-1,je+1)
       do i=max(3,is),min(npx-2,ie+1)
          qyy(i,j) = a2*(qy(i-2,j)+qy(i+1,j)) + a1*(qy(i-1,j)+qy(i,j))
       enddo
       if ( is==1 ) qyy(2,j) = c1*(qy(1,j)+qy(2,j))+c2*(qout(1,j)+qyy(3,j))
       if((ie+1)==npx) qyy(npx-1,j) = c1*(qy(npx-2,j)+qy(npx-1,j))+c2*(qout(npx,j)+qyy(npx-2,j))
 
       do i=max(2,is),min(npx-1,ie+1)
          qout(i,j) = 0.5*(qxx(i,j) + qyy(i,j))   ! averaging
       enddo
    enddo

 else  ! grid_type>=3
!------------------------
! Doubly periodic domain:
!------------------------
! X-sweep: PPM
    do j=js-2,je+2
       do i=is,ie+1
          qx(i,j) = b1*(qin(i-1,j)+qin(i,j)) + b2*(qin(i-2,j)+qin(i+1,j))
       enddo
    enddo
! Y-sweep: PPM
    do j=js,je+1
       do i=is-2,ie+2
          qy(i,j) = b1*(qin(i,j-1)+qin(i,j)) + b2*(qin(i,j-2)+qin(i,j+1))
       enddo
    enddo
    
    do j=js,je+1
       do i=is,ie+1
          qout(i,j) = 0.5*( a1*(qx(i,j-1)+qx(i,j  ) + qy(i-1,j)+qy(i,  j)) +  &
                            a2*(qx(i,j-2)+qx(i,j+1) + qy(i-2,j)+qy(i+1,j)) )
       enddo
    enddo
 endif

    if ( present(replace) ) then
       if ( replace ) then
          do j=js,je+1
          do i=is,ie+1
             qin(i,j) = qout(i,j)
          enddo
          enddo
       endif
    endif
    
  end subroutine a2b_ord4_r8


  subroutine a2b_ord2_r8(qin, qout, npx, npy, is, ie, js, je, ng, replace)
    integer, intent(IN   ) :: npx, npy, is, ie, js, je, ng
    real*8   , intent(INOUT) ::  qin(is-ng:ie+ng,js-ng:je+ng)   ! A-grid field
    real     , intent(  OUT) :: qout(is-ng:ie+ng,js-ng:je+ng)   ! Output  B-grid field
    logical, optional, intent(IN) ::  replace
    ! local:
    real   q1(npx), q2(npy)
    integer :: i,j
    integer :: is1, js1, is2, js2, ie1, je1

    if (grid_type < 3) then

    is1 = max(1,is-1)
    js1 = max(1,js-1)
    is2 = max(2,is)
    js2 = max(2,js)

    ie1 = min(npx-1,ie+1)
    je1 = min(npy-1,je+1)

    do j=js2,je1
       do i=is2,ie1
          qout(i,j) = 0.25*(qin(i-1,j-1)+qin(i,j-1)+qin(i-1,j)+qin(i,j))
       enddo
    enddo

! Fix the 4 Corners:
    if ( sw_corner ) qout(1,    1) = r3*(qin(1,        1)+qin(1,      0)+qin(0,      1))
    if ( se_corner ) qout(npx,  1) = r3*(qin(npx-1,    1)+qin(npx-1,  0)+qin(npx,    1))
    if ( ne_corner ) qout(npx,npy) = r3*(qin(npx-1,npy-1)+qin(npx,npy-1)+qin(npx-1,npy))
    if ( nw_corner ) qout(1,  npy) = r3*(qin(1,    npy-1)+qin(0,  npy-1)+qin(1,    npy))

    ! *** West Edges:
    if ( is==1 ) then
       do j=js1, je1
          q2(j) = 0.5*(qin(0,j) + qin(1,j))
       enddo
       do j=js2, je1
          qout(1,j) = edge_w(j)*q2(j-1) + (1.-edge_w(j))*q2(j)
       enddo
    endif

    ! East Edges:
    if ( (ie+1)==npx ) then
       do j=js1, je1
          q2(j) = 0.5*(qin(npx-1,j) + qin(npx,j))
       enddo
       do j=js2, je1
          qout(npx,j) = edge_e(j)*q2(j-1) + (1.-edge_e(j))*q2(j)
       enddo
    endif

    ! South Edges:
    if ( js==1 ) then
       do i=is1, ie1
          q1(i) = 0.5*(qin(i,0) + qin(i,1))
       enddo
       do i=is2, ie1
          qout(i,1) = edge_s(i)*q1(i-1) + (1.-edge_s(i))*q1(i)
       enddo
    endif

    ! North Edges:
    if ( (je+1)==npy ) then
       do i=is1, ie1
          q1(i) = 0.5*(qin(i,npy-1) + qin(i,npy))
       enddo
       do i=is2, ie1
          qout(i,npy) = edge_n(i)*q1(i-1) + (1.-edge_n(i))*q1(i)
       enddo
    endif

 else

    do j=js,je+1
       do i=is,ie+1
          qout(i,j) = 0.25*(qin(i-1,j-1)+qin(i,j-1)+qin(i-1,j)+qin(i,j))
       enddo
    enddo

 endif

    
    if ( present(replace) ) then
       if ( replace ) then
          do j=js,je+1
             do i=is,ie+1
                qin(i,j) = qout(i,j)
             enddo
          enddo
       endif
    endif
    
  end subroutine a2b_ord2_r8

  subroutine a2b_ord4_r4(qin, qout, npx, npy, is, ie, js, je, ng, replace)
  integer, intent(IN):: npx, npy, is, ie, js, je, ng
  real*4, intent(INOUT)::  qin(is-ng:ie+ng,js-ng:je+ng)   ! A-grid field
  real  , intent(INOUT):: qout(is-ng:ie+ng,js-ng:je+ng)   ! Output  B-grid field
  logical, optional, intent(IN):: replace
! local: compact 4-pt cubic
  real  , parameter:: c1 =  2./3.
  real  , parameter:: c2 = -1./6.
! 6-pt corner interpolation:
  real  , parameter:: d1 =  0.375                   !   0.5
  real  , parameter:: d2 = -1./24.                  !  -1./6.

  real   qx(is:ie+1,js-ng:je+ng)
  real   qy(is-ng:ie+ng,js:je+1)
  real   qxx(is-ng:ie+ng,js-ng:je+ng)
  real   qyy(is-ng:ie+ng,js-ng:je+ng)
  real   fx(is:ie), fy(is-2:ie+2,js:je)
  real   gratio, qt(npy)
  integer :: i, j, is1, js1, is2, js2, ie1, je1
  integer :: im2, jm2

    im2 = (npx-1)/2
    jm2 = (npy-1)/2

  if (grid_type < 3) then

    is1 = max(1,is-1)
    js1 = max(1,js-1)
    is2 = max(2,is)
    js2 = max(2,js)

    ie1 = min(npx-1,ie+1)
    je1 = min(npy-1,je+1)

! Corners:
#ifdef USE_3PT
    if ( sw_corner ) qout(1,    1) = r3*(qin(1,        1)+qin(1,      0)+qin(0,      1))
    if ( se_corner ) qout(npx,  1) = r3*(qin(npx-1,    1)+qin(npx-1,  0)+qin(npx,    1))
    if ( ne_corner ) qout(npx,npy) = r3*(qin(npx-1,npy-1)+qin(npx,npy-1)+qin(npx-1,npy))
    if ( nw_corner ) qout(1,  npy) = r3*(qin(1,    npy-1)+qin(0,  npy-1)+qin(1,    npy))
#else
! 6-point formular:
    if ( sw_corner ) then
        qout(1,1) = d1*(qin(1, 0) + qin( 0,1) + qin(1,1)) +  &
                    d2*(qin(2,-1) + qin(-1,2) + qin(2,2))
    endif
    if ( se_corner ) then
        qout(npx,1) = d1*(qin(npx-1, 0) + qin(npx-1,1) + qin(npx,  1)) +  &
                      d2*(qin(npx-2,-1) + qin(npx-2,2) + qin(npx+1,2))
    endif
    if ( ne_corner ) then
        qout(npx,npy) = d1*(qin(npx-1,npy-1) + qin(npx,  npy-1) + qin(npx-1,npy)) +  &
                        d2*(qin(npx-2,npy-2) + qin(npx+1,npy-2) + qin(npx-2,npy+1))
    endif
    if ( nw_corner ) then
        qout(1,npy) = d1*(qin( 0,npy-1) + qin(1,npy-1) + qin(1,npy)) +   &
                      d2*(qin(-1,npy-2) + qin(2,npy-2) + qin(2,npy+1))
    endif
#endif

!------------
! X-Interior:
!------------
    do j=max(1,js-2),min(npy-1,je+2)
       do i=max(3,is), min(npx-2,ie+1)
          qx(i,j) = b2*(qin(i-2,j)+qin(i+1,j)) + b1*(qin(i-1,j)+qin(i,j))
       enddo
    enddo

! West Edges:
    if ( is==1 ) then

       do j=max(1,js-2),min(npy-1,je+2)
           gratio = dxa(2,j) / dxa(1,j)
          qx(1,j) = 0.5*((2.+gratio)*(qin(0,j)+qin(1,j))    &
                  - (qin(-1,j)+qin(2,j))) / (1.+gratio)
#ifdef TEST2
! Note: Caused noises in test_case-5 for large n_split
          qx(2,j) = (2.*gratio*(gratio+1.)*qin(1,j)+qin(2,j) -     &
                     gratio*(gratio+0.5)*qx(1,j))/(1.+gratio*(gratio+1.5))
#else
          qx(2,j) = (3.*(gratio*qin(1,j)+qin(2,j)) - (gratio*qx(1,j)+qx(3,j)))/(2.+2.*gratio)
#endif
       enddo

       do j=max(3,js),min(npy-2,je+1)
          qout(1,j) = a2*(qx(1,j-2)+qx(1,j+1)) + a1*(qx(1,j-1)+qx(1,j))
       enddo

       if( js==1 )     qout(1,    2) = c1*(qx(1,1)+qx(1,2))         + c2*(qout(1,1)+qout(1,3))
       if((je+1)==npy) qout(1,npy-1) = c1*(qx(1,npy-2)+qx(1,npy-1)) + c2*(qout(1,npy-2)+qout(1,npy))
    endif

! East Edges:
    if ( (ie+1)==npx ) then

       do j=max(1,js-2),min(npy-1,je+2)
               gratio = dxa(npx-2,j) / dxa(npx-1,j)
          qx(npx  ,j) = 0.5*((2.+gratio)*(qin(npx-1,j)+qin(npx,j))   &
                        - (qin(npx-2,j)+qin(npx+1,j))) / (1.+gratio )
#ifdef TEST2
          qx(npx-1,j) = (2.*gratio*(gratio+1.)*qin(npx-1,j)+qin(npx-2,j) -  &
                         gratio*(gratio+0.5)*qx(npx,j))/(1.+gratio*(gratio+1.5))
#else
          qx(npx-1,j) = (3.*(qin(npx-2,j)+gratio*qin(npx-1,j)) - (gratio*qx(npx,j)+qx(npx-2,j)))/(2.+2.*gratio)
#endif
       enddo

       do j=max(3,js),min(npy-2,je+1)
          qout(npx,j) = a2*(qx(npx,j-2)+qx(npx,j+1)) + a1*(qx(npx,j-1)+qx(npx,j))
       enddo

       if(js==1) qout(npx,2) = c1*(qx(npx,1)+qx(npx,2))+c2*(qout(npx,1)+qout(npx,3))
       if((je+1)==npy)   qout(npx,npy-1) =             &
                         c1*(qx(npx,npy-2)+qx(npx,npy-1))+c2*(qout(npx,npy-2)+qout(npx,npy))
    endif

!------------
! Y-Interior:
!------------
    do j=max(3,js),min(npy-2,je+1)
       do i=max(1,is-2), min(npx-1,ie+2)
          qy(i,j) = b2*(qin(i,j-2)+qin(i,j+1)) + b1*(qin(i,j-1)+qin(i,j))
       enddo
    enddo

! South Edges:
    if ( js==1 ) then

       do i=max(1,is-2),min(npx-1,ie+2)
           gratio = dya(i,2) / dya(i,1)
          qy(i,1) = 0.5*((2.+gratio)*(qin(i,0)+qin(i,1))   &
                  - (qin(i,-1)+qin(i,2))) / (1.+gratio )
#ifdef TEST2
          qy(i,2) = (2.*gratio*(gratio+1.)*qin(i,1)+qin(i,2) -     &
                     gratio*(gratio+0.5)*qy(i,1))/(1.+gratio*(gratio+1.5))
#else
          qy(i,2) = (3.*(gratio*qin(i,1)+qin(i,2)) - (gratio*qy(i,1)+qy(i,3)))/(2.+2.*gratio)
#endif
       enddo

       do i=max(3,is),min(npx-2,ie+1)
          qout(i,1) = a2*(qy(i-2,1)+qy(i+1,1)) + a1*(qy(i-1,1)+qy(i,1))
       enddo

       if( is==1 )    qout(2,1) = c1*(qy(1,1)+qy(2,1))+c2*(qout(1,1)+qout(3,1))
       if((ie+1)==npx) qout(npx-1,1) = c1*(qy(npx-2,1)+qy(npx-1,1))+c2*(qout(npx-2,1)+qout(npx,1))
    endif

! North Edges:
    if ( (je+1)==npy ) then
       do i=max(1,is-2),min(npx-1,ie+2)
               gratio = dya(i,npy-2) / dya(i,npy-1)
          qy(i,npy  ) = 0.5*((2.+gratio)*(qin(i,npy-1)+qin(i,npy))  &
                      - (qin(i,npy-2)+qin(i,npy+1))) / (1.+gratio)
#ifdef TEST2
          qy(i,npy-1) = (2.*gratio*(gratio+1.)*qin(i,npy-1)+qin(i,npy-2) - &
                         gratio*(gratio+0.5)*qy(i,npy))/(1.+gratio*(gratio+1.5))
#else
          qy(i,npy-1) = (3.*(qin(i,npy-2)+gratio*qin(i,npy-1)) - (gratio*qy(i,npy)+qy(i,npy-2)))/(2.+2.*gratio)
#endif
       enddo

       do i=max(3,is),min(npx-2,ie+1)
          qout(i,npy) = a2*(qy(i-2,npy)+qy(i+1,npy)) + a1*(qy(i-1,npy)+qy(i,npy))
       enddo

       if( is==1 )  qout(2,npy) = c1*(qy(1,npy)+qy(2,npy))+c2*(qout(1,npy)+qout(3,npy))
       if((ie+1)==npx) qout(npx-1,npy) = c1*(qy(npx-2,npy)+qy(npx-1,npy))+c2*(qout(npx-2,npy)+qout(npx,npy))
    endif
    
    do j=max(3,js),min(npy-2,je+1)
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,j) = a2*(qx(i,j-2)+qx(i,j+1)) + a1*(qx(i,j-1)+qx(i,j))
       enddo
    enddo

    if ( js==1 ) then
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,2) = c1*(qx(i,1)+qx(i,2))+c2*(qout(i,1)+qxx(i,3))
       enddo
    endif
    if ( (je+1)==npy ) then
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,npy-1) = c1*(qx(i,npy-2)+qx(i,npy-1))+c2*(qout(i,npy)+qxx(i,npy-2))
       enddo
    endif

    
    do j=max(2,js),min(npy-1,je+1)
       do i=max(3,is),min(npx-2,ie+1)
          qyy(i,j) = a2*(qy(i-2,j)+qy(i+1,j)) + a1*(qy(i-1,j)+qy(i,j))
       enddo
       if ( is==1 ) qyy(2,j) = c1*(qy(1,j)+qy(2,j))+c2*(qout(1,j)+qyy(3,j))
       if((ie+1)==npx) qyy(npx-1,j) = c1*(qy(npx-2,j)+qy(npx-1,j))+c2*(qout(npx,j)+qyy(npx-2,j))
 
       do i=max(2,is),min(npx-1,ie+1)
          qout(i,j) = 0.5*(qxx(i,j) + qyy(i,j))   ! averaging
       enddo
    enddo

 else  ! grid_type>=3
!------------------------
! Doubly periodic domain:
!------------------------
! X-sweep: PPM
    do j=js-2,je+2
       do i=is,ie+1
          qx(i,j) = b1*(qin(i-1,j)+qin(i,j)) + b2*(qin(i-2,j)+qin(i+1,j))
       enddo
    enddo
! Y-sweep: PPM
    do j=js,je+1
       do i=is-2,ie+2
          qy(i,j) = b1*(qin(i,j-1)+qin(i,j)) + b2*(qin(i,j-2)+qin(i,j+1))
       enddo
    enddo
    
    do j=js,je+1
       do i=is,ie+1
          qout(i,j) = 0.5*( a1*(qx(i,j-1)+qx(i,j  ) + qy(i-1,j)+qy(i,  j)) +  &
                            a2*(qx(i,j-2)+qx(i,j+1) + qy(i-2,j)+qy(i+1,j)) )
       enddo
    enddo
 endif

    if ( present(replace) ) then
       if ( replace ) then
          do j=js,je+1
          do i=is,ie+1
             qin(i,j) = qout(i,j)
          enddo
          enddo
       endif
    endif
    
  end subroutine a2b_ord4_r4

  subroutine a2b_ord2_r4(qin, qout, npx, npy, is, ie, js, je, ng, replace)
    integer, intent(IN   ) :: npx, npy, is, ie, js, je, ng
    real*4   , intent(INOUT) ::  qin(is-ng:ie+ng,js-ng:je+ng)   ! A-grid field
    real     , intent(  OUT) :: qout(is-ng:ie+ng,js-ng:je+ng)   ! Output  B-grid field
    logical, optional, intent(IN) ::  replace
    ! local:
    real   q1(npx), q2(npy)
    integer :: i,j
    integer :: is1, js1, is2, js2, ie1, je1

    if (grid_type < 3) then

    is1 = max(1,is-1)
    js1 = max(1,js-1)
    is2 = max(2,is)
    js2 = max(2,js)

    ie1 = min(npx-1,ie+1)
    je1 = min(npy-1,je+1)

    do j=js2,je1
       do i=is2,ie1
          qout(i,j) = 0.25*(qin(i-1,j-1)+qin(i,j-1)+qin(i-1,j)+qin(i,j))
       enddo
    enddo

! Fix the 4 Corners:
    if ( sw_corner ) qout(1,    1) = r3*(qin(1,        1)+qin(1,      0)+qin(0,      1))
    if ( se_corner ) qout(npx,  1) = r3*(qin(npx-1,    1)+qin(npx-1,  0)+qin(npx,    1))
    if ( ne_corner ) qout(npx,npy) = r3*(qin(npx-1,npy-1)+qin(npx,npy-1)+qin(npx-1,npy))
    if ( nw_corner ) qout(1,  npy) = r3*(qin(1,    npy-1)+qin(0,  npy-1)+qin(1,    npy))

    ! *** West Edges:
    if ( is==1 ) then
       do j=js1, je1
          q2(j) = 0.5*(qin(0,j) + qin(1,j))
       enddo
       do j=js2, je1
          qout(1,j) = edge_w(j)*q2(j-1) + (1.-edge_w(j))*q2(j)
       enddo
    endif

    ! East Edges:
    if ( (ie+1)==npx ) then
       do j=js1, je1
          q2(j) = 0.5*(qin(npx-1,j) + qin(npx,j))
       enddo
       do j=js2, je1
          qout(npx,j) = edge_e(j)*q2(j-1) + (1.-edge_e(j))*q2(j)
       enddo
    endif

    ! South Edges:
    if ( js==1 ) then
       do i=is1, ie1
          q1(i) = 0.5*(qin(i,0) + qin(i,1))
       enddo
       do i=is2, ie1
          qout(i,1) = edge_s(i)*q1(i-1) + (1.-edge_s(i))*q1(i)
       enddo
    endif

    ! North Edges:
    if ( (je+1)==npy ) then
       do i=is1, ie1
          q1(i) = 0.5*(qin(i,npy-1) + qin(i,npy))
       enddo
       do i=is2, ie1
          qout(i,npy) = edge_n(i)*q1(i-1) + (1.-edge_n(i))*q1(i)
       enddo
    endif

 else

    do j=js,je+1
       do i=is,ie+1
          qout(i,j) = 0.25*(qin(i-1,j-1)+qin(i,j-1)+qin(i-1,j)+qin(i,j))
       enddo
    enddo

 endif

    
    if ( present(replace) ) then
       if ( replace ) then
          do j=js,je+1
             do i=is,ie+1
                qin(i,j) = qout(i,j)
             enddo
          enddo
       endif
    endif
    
  end subroutine a2b_ord2_r4

  
end module a2b_edge_mod
