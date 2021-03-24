module a2b_edge_tlm_mod

  use fv_grid_utils_mod, only: edge_w, edge_e, edge_s, edge_n, sw_corner, se_corner,  &
                               nw_corner, ne_corner, van2
  use fv_grid_tools_mod, only: dxa, dya, grid_type

  implicit none

  real, parameter:: r3 = 1./3.

  private
  public :: a2b_ord2_tlm, a2b_ord4_tlm

contains

  SUBROUTINE A2B_ORD4_TLM(qin, qin_tl, qout, qout_tl, npx, npy, is, ie, &
&   js, je, ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL*8, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
    REAL*8, INTENT(INOUT) :: qin_tl(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL, INTENT(INOUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: qout_tl(is-ng:ie+ng, js-ng:je+ng)
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
! local: compact 4-pt cubic
    REAL, PARAMETER :: c1=2./3.
    REAL, PARAMETER :: c2=-(1./6.)
! Parabolic spline
! real  , parameter:: c1 =  0.75
! real  , parameter:: c2 = -0.25
! 6-pt corner interpolation:
!   0.5
    REAL, PARAMETER :: d1=0.375
!  -1./6.
    REAL, PARAMETER :: d2=-(1./24.)
    REAL :: qx(is:ie+1, js-ng:je+ng)
    REAL :: qx_tl(is:ie+1, js-ng:je+ng)
    REAL :: qy(is-ng:ie+ng, js:je+1)
    REAL :: qy_tl(is-ng:ie+ng, js:je+1)
    REAL :: qxx(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qxx_tl(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy_tl(is-ng:ie+ng, js-ng:je+ng)
    REAL :: fx(is:ie), fy(is-2:ie+2, js:je)
    REAL :: gratio, qt(npy)
    INTEGER :: i, j, is1, js1, is2, js2, ie1, je1
    INTEGER :: im2, jm2
!  9/16
    REAL, PARAMETER :: a1=0.5625
! -1/16
    REAL, PARAMETER :: a2=-0.0625
! 0.58333333
    REAL, PARAMETER :: b1=7./12.
    REAL, PARAMETER :: b2=-(1./12.)
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC PRESENT
    INTEGER :: min9
    INTEGER :: min8
    INTEGER :: min7
    INTEGER :: min6
    INTEGER :: min5
    INTEGER :: min4
    INTEGER :: min3
    INTEGER :: min2
    INTEGER :: min1
    INTEGER :: min19
    INTEGER :: min18
    INTEGER :: min17
    INTEGER :: min16
    INTEGER :: min15
    INTEGER :: min14
    INTEGER :: min13
    INTEGER :: min12
    INTEGER :: min11
    INTEGER :: min10
    INTEGER :: max9
    INTEGER :: max8
    INTEGER :: max19
    INTEGER :: max7
    INTEGER :: max18
    INTEGER :: max6
    INTEGER :: max17
    INTEGER :: max5
    INTEGER :: max16
    INTEGER :: max4
    INTEGER :: max15
    INTEGER :: max3
    INTEGER :: max14
    INTEGER :: max2
    INTEGER :: max13
    INTEGER :: max1
    INTEGER :: max12
    INTEGER :: max11
    INTEGER :: max10
    im2 = (npx-1)/2
    jm2 = (npy-1)/2
    IF (grid_type .LT. 3) THEN
      IF (1 .LT. is - 1) THEN
        is1 = is - 1
      ELSE
        is1 = 1
      END IF
      IF (1 .LT. js - 1) THEN
        js1 = js - 1
      ELSE
        js1 = 1
      END IF
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (2 .LT. js) THEN
        js2 = js
      ELSE
        js2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 1
      END IF
      IF (npy - 1 .GT. je + 1) THEN
        je1 = je + 1
      ELSE
        je1 = npy - 1
      END IF
! Corners:
! 6-point formular:
      IF (sw_corner) THEN
        qout_tl(1, 1) = d1*(qin_tl(1, 0)+qin_tl(0, 1)+qin_tl(1, 1)) + d2&
&         *(qin_tl(2, -1)+qin_tl(-1, 2)+qin_tl(2, 2))
        qout(1, 1) = d1*(qin(1, 0)+qin(0, 1)+qin(1, 1)) + d2*(qin(2, -1)&
&         +qin(-1, 2)+qin(2, 2))
      END IF
      IF (se_corner) THEN
        qout_tl(npx, 1) = d1*(qin_tl(npx-1, 0)+qin_tl(npx-1, 1)+qin_tl(&
&         npx, 1)) + d2*(qin_tl(npx-2, -1)+qin_tl(npx-2, 2)+qin_tl(npx+1&
&         , 2))
        qout(npx, 1) = d1*(qin(npx-1, 0)+qin(npx-1, 1)+qin(npx, 1)) + d2&
&         *(qin(npx-2, -1)+qin(npx-2, 2)+qin(npx+1, 2))
      END IF
      IF (ne_corner) THEN
        qout_tl(npx, npy) = d1*(qin_tl(npx-1, npy-1)+qin_tl(npx, npy-1)+&
&         qin_tl(npx-1, npy)) + d2*(qin_tl(npx-2, npy-2)+qin_tl(npx+1, &
&         npy-2)+qin_tl(npx-2, npy+1))
        qout(npx, npy) = d1*(qin(npx-1, npy-1)+qin(npx, npy-1)+qin(npx-1&
&         , npy)) + d2*(qin(npx-2, npy-2)+qin(npx+1, npy-2)+qin(npx-2, &
&         npy+1))
      END IF
      IF (nw_corner) THEN
        qout_tl(1, npy) = d1*(qin_tl(0, npy-1)+qin_tl(1, npy-1)+qin_tl(1&
&         , npy)) + d2*(qin_tl(-1, npy-2)+qin_tl(2, npy-2)+qin_tl(2, npy&
&         +1))
        qout(1, npy) = d1*(qin(0, npy-1)+qin(1, npy-1)+qin(1, npy)) + d2&
&         *(qin(-1, npy-2)+qin(2, npy-2)+qin(2, npy+1))
      END IF
      IF (1 .LT. js - 2) THEN
        max1 = js - 2
      ELSE
        max1 = 1
      END IF
      IF (npy - 1 .GT. je + 2) THEN
        min1 = je + 2
        qx_tl = 0.0_8
      ELSE
        min1 = npy - 1
        qx_tl = 0.0_8
      END IF
!------------
! X-Interior:
!------------
      DO j=max1,min1
        IF (3 .LT. is) THEN
          max2 = is
        ELSE
          max2 = 3
        END IF
        IF (npx - 2 .GT. ie + 1) THEN
          min2 = ie + 1
        ELSE
          min2 = npx - 2
        END IF
        DO i=max2,min2
          qx_tl(i, j) = b2*(qin_tl(i-2, j)+qin_tl(i+1, j)) + b1*(qin_tl(&
&           i-1, j)+qin_tl(i, j))
          qx(i, j) = b2*(qin(i-2, j)+qin(i+1, j)) + b1*(qin(i-1, j)+qin(&
&           i, j))
        END DO
      END DO
! West Edges:
      IF (is .EQ. 1) THEN
        IF (1 .LT. js - 2) THEN
          max3 = js - 2
        ELSE
          max3 = 1
        END IF
        IF (npy - 1 .GT. je + 2) THEN
          min3 = je + 2
        ELSE
          min3 = npy - 1
        END IF
        DO j=max3,min3
          gratio = dxa(2, j)/dxa(1, j)
          qx_tl(1, j) = 0.5*((2.+gratio)*(qin_tl(0, j)+qin_tl(1, j))-&
&           qin_tl(-1, j)-qin_tl(2, j))/(1.+gratio)
          qx(1, j) = 0.5*((2.+gratio)*(qin(0, j)+qin(1, j))-(qin(-1, j)+&
&           qin(2, j)))/(1.+gratio)
          qx_tl(2, j) = (3.*(gratio*qin_tl(1, j)+qin_tl(2, j))-gratio*&
&           qx_tl(1, j)-qx_tl(3, j))/(2.+2.*gratio)
          qx(2, j) = (3.*(gratio*qin(1, j)+qin(2, j))-(gratio*qx(1, j)+&
&           qx(3, j)))/(2.+2.*gratio)
        END DO
        IF (3 .LT. js) THEN
          max4 = js
        ELSE
          max4 = 3
        END IF
        IF (npy - 2 .GT. je + 1) THEN
          min4 = je + 1
        ELSE
          min4 = npy - 2
        END IF
        DO j=max4,min4
          qout_tl(1, j) = a2*(qx_tl(1, j-2)+qx_tl(1, j+1)) + a1*(qx_tl(1&
&           , j-1)+qx_tl(1, j))
          qout(1, j) = a2*(qx(1, j-2)+qx(1, j+1)) + a1*(qx(1, j-1)+qx(1&
&           , j))
        END DO
        IF (js .EQ. 1) THEN
          qout_tl(1, 2) = c1*(qx_tl(1, 1)+qx_tl(1, 2)) + c2*(qout_tl(1, &
&           1)+qout_tl(1, 3))
          qout(1, 2) = c1*(qx(1, 1)+qx(1, 2)) + c2*(qout(1, 1)+qout(1, 3&
&           ))
        END IF
        IF (je + 1 .EQ. npy) THEN
          qout_tl(1, npy-1) = c1*(qx_tl(1, npy-2)+qx_tl(1, npy-1)) + c2*&
&           (qout_tl(1, npy-2)+qout_tl(1, npy))
          qout(1, npy-1) = c1*(qx(1, npy-2)+qx(1, npy-1)) + c2*(qout(1, &
&           npy-2)+qout(1, npy))
        END IF
      END IF
! East Edges:
      IF (ie + 1 .EQ. npx) THEN
        IF (1 .LT. js - 2) THEN
          max5 = js - 2
        ELSE
          max5 = 1
        END IF
        IF (npy - 1 .GT. je + 2) THEN
          min5 = je + 2
        ELSE
          min5 = npy - 1
        END IF
        DO j=max5,min5
          gratio = dxa(npx-2, j)/dxa(npx-1, j)
          qx_tl(npx, j) = 0.5*((2.+gratio)*(qin_tl(npx-1, j)+qin_tl(npx&
&           , j))-qin_tl(npx-2, j)-qin_tl(npx+1, j))/(1.+gratio)
          qx(npx, j) = 0.5*((2.+gratio)*(qin(npx-1, j)+qin(npx, j))-(qin&
&           (npx-2, j)+qin(npx+1, j)))/(1.+gratio)
          qx_tl(npx-1, j) = (3.*(qin_tl(npx-2, j)+gratio*qin_tl(npx-1, j&
&           ))-gratio*qx_tl(npx, j)-qx_tl(npx-2, j))/(2.+2.*gratio)
          qx(npx-1, j) = (3.*(qin(npx-2, j)+gratio*qin(npx-1, j))-(&
&           gratio*qx(npx, j)+qx(npx-2, j)))/(2.+2.*gratio)
        END DO
        IF (3 .LT. js) THEN
          max6 = js
        ELSE
          max6 = 3
        END IF
        IF (npy - 2 .GT. je + 1) THEN
          min6 = je + 1
        ELSE
          min6 = npy - 2
        END IF
        DO j=max6,min6
          qout_tl(npx, j) = a2*(qx_tl(npx, j-2)+qx_tl(npx, j+1)) + a1*(&
&           qx_tl(npx, j-1)+qx_tl(npx, j))
          qout(npx, j) = a2*(qx(npx, j-2)+qx(npx, j+1)) + a1*(qx(npx, j-&
&           1)+qx(npx, j))
        END DO
        IF (js .EQ. 1) THEN
          qout_tl(npx, 2) = c1*(qx_tl(npx, 1)+qx_tl(npx, 2)) + c2*(&
&           qout_tl(npx, 1)+qout_tl(npx, 3))
          qout(npx, 2) = c1*(qx(npx, 1)+qx(npx, 2)) + c2*(qout(npx, 1)+&
&           qout(npx, 3))
        END IF
        IF (je + 1 .EQ. npy) THEN
          qout_tl(npx, npy-1) = c1*(qx_tl(npx, npy-2)+qx_tl(npx, npy-1))&
&           + c2*(qout_tl(npx, npy-2)+qout_tl(npx, npy))
          qout(npx, npy-1) = c1*(qx(npx, npy-2)+qx(npx, npy-1)) + c2*(&
&           qout(npx, npy-2)+qout(npx, npy))
        END IF
      END IF
      IF (3 .LT. js) THEN
        max7 = js
      ELSE
        max7 = 3
      END IF
      IF (npy - 2 .GT. je + 1) THEN
        min7 = je + 1
        qy_tl = 0.0_8
      ELSE
        min7 = npy - 2
        qy_tl = 0.0_8
      END IF
!------------
! Y-Interior:
!------------
      DO j=max7,min7
        IF (1 .LT. is - 2) THEN
          max8 = is - 2
        ELSE
          max8 = 1
        END IF
        IF (npx - 1 .GT. ie + 2) THEN
          min8 = ie + 2
        ELSE
          min8 = npx - 1
        END IF
        DO i=max8,min8
          qy_tl(i, j) = b2*(qin_tl(i, j-2)+qin_tl(i, j+1)) + b1*(qin_tl(&
&           i, j-1)+qin_tl(i, j))
          qy(i, j) = b2*(qin(i, j-2)+qin(i, j+1)) + b1*(qin(i, j-1)+qin(&
&           i, j))
        END DO
      END DO
! South Edges:
      IF (js .EQ. 1) THEN
        IF (1 .LT. is - 2) THEN
          max9 = is - 2
        ELSE
          max9 = 1
        END IF
        IF (npx - 1 .GT. ie + 2) THEN
          min9 = ie + 2
        ELSE
          min9 = npx - 1
        END IF
        DO i=max9,min9
          gratio = dya(i, 2)/dya(i, 1)
          qy_tl(i, 1) = 0.5*((2.+gratio)*(qin_tl(i, 0)+qin_tl(i, 1))-&
&           qin_tl(i, -1)-qin_tl(i, 2))/(1.+gratio)
          qy(i, 1) = 0.5*((2.+gratio)*(qin(i, 0)+qin(i, 1))-(qin(i, -1)+&
&           qin(i, 2)))/(1.+gratio)
          qy_tl(i, 2) = (3.*(gratio*qin_tl(i, 1)+qin_tl(i, 2))-gratio*&
&           qy_tl(i, 1)-qy_tl(i, 3))/(2.+2.*gratio)
          qy(i, 2) = (3.*(gratio*qin(i, 1)+qin(i, 2))-(gratio*qy(i, 1)+&
&           qy(i, 3)))/(2.+2.*gratio)
        END DO
        IF (3 .LT. is) THEN
          max10 = is
        ELSE
          max10 = 3
        END IF
        IF (npx - 2 .GT. ie + 1) THEN
          min10 = ie + 1
        ELSE
          min10 = npx - 2
        END IF
        DO i=max10,min10
          qout_tl(i, 1) = a2*(qy_tl(i-2, 1)+qy_tl(i+1, 1)) + a1*(qy_tl(i&
&           -1, 1)+qy_tl(i, 1))
          qout(i, 1) = a2*(qy(i-2, 1)+qy(i+1, 1)) + a1*(qy(i-1, 1)+qy(i&
&           , 1))
        END DO
        IF (is .EQ. 1) THEN
          qout_tl(2, 1) = c1*(qy_tl(1, 1)+qy_tl(2, 1)) + c2*(qout_tl(1, &
&           1)+qout_tl(3, 1))
          qout(2, 1) = c1*(qy(1, 1)+qy(2, 1)) + c2*(qout(1, 1)+qout(3, 1&
&           ))
        END IF
        IF (ie + 1 .EQ. npx) THEN
          qout_tl(npx-1, 1) = c1*(qy_tl(npx-2, 1)+qy_tl(npx-1, 1)) + c2*&
&           (qout_tl(npx-2, 1)+qout_tl(npx, 1))
          qout(npx-1, 1) = c1*(qy(npx-2, 1)+qy(npx-1, 1)) + c2*(qout(npx&
&           -2, 1)+qout(npx, 1))
        END IF
      END IF
! North Edges:
      IF (je + 1 .EQ. npy) THEN
        IF (1 .LT. is - 2) THEN
          max11 = is - 2
        ELSE
          max11 = 1
        END IF
        IF (npx - 1 .GT. ie + 2) THEN
          min11 = ie + 2
        ELSE
          min11 = npx - 1
        END IF
        DO i=max11,min11
          gratio = dya(i, npy-2)/dya(i, npy-1)
          qy_tl(i, npy) = 0.5*((2.+gratio)*(qin_tl(i, npy-1)+qin_tl(i, &
&           npy))-qin_tl(i, npy-2)-qin_tl(i, npy+1))/(1.+gratio)
          qy(i, npy) = 0.5*((2.+gratio)*(qin(i, npy-1)+qin(i, npy))-(qin&
&           (i, npy-2)+qin(i, npy+1)))/(1.+gratio)
          qy_tl(i, npy-1) = (3.*(qin_tl(i, npy-2)+gratio*qin_tl(i, npy-1&
&           ))-gratio*qy_tl(i, npy)-qy_tl(i, npy-2))/(2.+2.*gratio)
          qy(i, npy-1) = (3.*(qin(i, npy-2)+gratio*qin(i, npy-1))-(&
&           gratio*qy(i, npy)+qy(i, npy-2)))/(2.+2.*gratio)
        END DO
        IF (3 .LT. is) THEN
          max12 = is
        ELSE
          max12 = 3
        END IF
        IF (npx - 2 .GT. ie + 1) THEN
          min12 = ie + 1
        ELSE
          min12 = npx - 2
        END IF
        DO i=max12,min12
          qout_tl(i, npy) = a2*(qy_tl(i-2, npy)+qy_tl(i+1, npy)) + a1*(&
&           qy_tl(i-1, npy)+qy_tl(i, npy))
          qout(i, npy) = a2*(qy(i-2, npy)+qy(i+1, npy)) + a1*(qy(i-1, &
&           npy)+qy(i, npy))
        END DO
        IF (is .EQ. 1) THEN
          qout_tl(2, npy) = c1*(qy_tl(1, npy)+qy_tl(2, npy)) + c2*(&
&           qout_tl(1, npy)+qout_tl(3, npy))
          qout(2, npy) = c1*(qy(1, npy)+qy(2, npy)) + c2*(qout(1, npy)+&
&           qout(3, npy))
        END IF
        IF (ie + 1 .EQ. npx) THEN
          qout_tl(npx-1, npy) = c1*(qy_tl(npx-2, npy)+qy_tl(npx-1, npy))&
&           + c2*(qout_tl(npx-2, npy)+qout_tl(npx, npy))
          qout(npx-1, npy) = c1*(qy(npx-2, npy)+qy(npx-1, npy)) + c2*(&
&           qout(npx-2, npy)+qout(npx, npy))
        END IF
      END IF
      IF (3 .LT. js) THEN
        max13 = js
      ELSE
        max13 = 3
      END IF
      IF (npy - 2 .GT. je + 1) THEN
        min13 = je + 1
        qxx_tl = 0.0_8
      ELSE
        min13 = npy - 2
        qxx_tl = 0.0_8
      END IF
      DO j=max13,min13
        IF (2 .LT. is) THEN
          max14 = is
        ELSE
          max14 = 2
        END IF
        IF (npx - 1 .GT. ie + 1) THEN
          min14 = ie + 1
        ELSE
          min14 = npx - 1
        END IF
        DO i=max14,min14
          qxx_tl(i, j) = a2*(qx_tl(i, j-2)+qx_tl(i, j+1)) + a1*(qx_tl(i&
&           , j-1)+qx_tl(i, j))
          qxx(i, j) = a2*(qx(i, j-2)+qx(i, j+1)) + a1*(qx(i, j-1)+qx(i, &
&           j))
        END DO
      END DO
      IF (js .EQ. 1) THEN
        IF (2 .LT. is) THEN
          max15 = is
        ELSE
          max15 = 2
        END IF
        IF (npx - 1 .GT. ie + 1) THEN
          min15 = ie + 1
        ELSE
          min15 = npx - 1
        END IF
        DO i=max15,min15
          qxx_tl(i, 2) = c1*(qx_tl(i, 1)+qx_tl(i, 2)) + c2*(qout_tl(i, 1&
&           )+qxx_tl(i, 3))
          qxx(i, 2) = c1*(qx(i, 1)+qx(i, 2)) + c2*(qout(i, 1)+qxx(i, 3))
        END DO
      END IF
      IF (je + 1 .EQ. npy) THEN
        IF (2 .LT. is) THEN
          max16 = is
        ELSE
          max16 = 2
        END IF
        IF (npx - 1 .GT. ie + 1) THEN
          min16 = ie + 1
        ELSE
          min16 = npx - 1
        END IF
        DO i=max16,min16
          qxx_tl(i, npy-1) = c1*(qx_tl(i, npy-2)+qx_tl(i, npy-1)) + c2*(&
&           qout_tl(i, npy)+qxx_tl(i, npy-2))
          qxx(i, npy-1) = c1*(qx(i, npy-2)+qx(i, npy-1)) + c2*(qout(i, &
&           npy)+qxx(i, npy-2))
        END DO
      END IF
      IF (2 .LT. js) THEN
        max17 = js
      ELSE
        max17 = 2
      END IF
      IF (npy - 1 .GT. je + 1) THEN
        min17 = je + 1
        qyy_tl = 0.0_8
      ELSE
        min17 = npy - 1
        qyy_tl = 0.0_8
      END IF
      DO j=max17,min17
        IF (3 .LT. is) THEN
          max18 = is
        ELSE
          max18 = 3
        END IF
        IF (npx - 2 .GT. ie + 1) THEN
          min18 = ie + 1
        ELSE
          min18 = npx - 2
        END IF
        DO i=max18,min18
          qyy_tl(i, j) = a2*(qy_tl(i-2, j)+qy_tl(i+1, j)) + a1*(qy_tl(i-&
&           1, j)+qy_tl(i, j))
          qyy(i, j) = a2*(qy(i-2, j)+qy(i+1, j)) + a1*(qy(i-1, j)+qy(i, &
&           j))
        END DO
        IF (is .EQ. 1) THEN
          qyy_tl(2, j) = c1*(qy_tl(1, j)+qy_tl(2, j)) + c2*(qout_tl(1, j&
&           )+qyy_tl(3, j))
          qyy(2, j) = c1*(qy(1, j)+qy(2, j)) + c2*(qout(1, j)+qyy(3, j))
        END IF
        IF (ie + 1 .EQ. npx) THEN
          qyy_tl(npx-1, j) = c1*(qy_tl(npx-2, j)+qy_tl(npx-1, j)) + c2*(&
&           qout_tl(npx, j)+qyy_tl(npx-2, j))
          qyy(npx-1, j) = c1*(qy(npx-2, j)+qy(npx-1, j)) + c2*(qout(npx&
&           , j)+qyy(npx-2, j))
        END IF
        IF (2 .LT. is) THEN
          max19 = is
        ELSE
          max19 = 2
        END IF
        IF (npx - 1 .GT. ie + 1) THEN
          min19 = ie + 1
        ELSE
          min19 = npx - 1
        END IF
        DO i=max19,min19
! averaging
          qout_tl(i, j) = 0.5*(qxx_tl(i, j)+qyy_tl(i, j))
          qout(i, j) = 0.5*(qxx(i, j)+qyy(i, j))
        END DO
      END DO
    ELSE
      qx_tl = 0.0_8
! grid_type>=3
!------------------------
! Doubly periodic domain:
!------------------------
! X-sweep: PPM
      DO j=js-2,je+2
        DO i=is,ie+1
          qx_tl(i, j) = b1*(qin_tl(i-1, j)+qin_tl(i, j)) + b2*(qin_tl(i-&
&           2, j)+qin_tl(i+1, j))
          qx(i, j) = b1*(qin(i-1, j)+qin(i, j)) + b2*(qin(i-2, j)+qin(i+&
&           1, j))
        END DO
      END DO
      qy_tl = 0.0_8
! Y-sweep: PPM
      DO j=js,je+1
        DO i=is-2,ie+2
          qy_tl(i, j) = b1*(qin_tl(i, j-1)+qin_tl(i, j)) + b2*(qin_tl(i&
&           , j-2)+qin_tl(i, j+1))
          qy(i, j) = b1*(qin(i, j-1)+qin(i, j)) + b2*(qin(i, j-2)+qin(i&
&           , j+1))
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie+1
          qout_tl(i, j) = 0.5*(a1*(qx_tl(i, j-1)+qx_tl(i, j)+qy_tl(i-1, &
&           j)+qy_tl(i, j))+a2*(qx_tl(i, j-2)+qx_tl(i, j+1)+qy_tl(i-2, j&
&           )+qy_tl(i+1, j)))
          qout(i, j) = 0.5*(a1*(qx(i, j-1)+qx(i, j)+qy(i-1, j)+qy(i, j))&
&           +a2*(qx(i, j-2)+qx(i, j+1)+qy(i-2, j)+qy(i+1, j)))
        END DO
      END DO
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=js,je+1
          DO i=is,ie+1
            qin_tl(i, j) = qout_tl(i, j)
            qin(i, j) = qout(i, j)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE A2B_ORD4_TLM

  SUBROUTINE A2B_ORD2_TLM(qin, qin_tl, qout, qout_tl, npx, npy, is, ie, &
&   js, je, ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL*8, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
    REAL*8, INTENT(INOUT) :: qin_tl(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL, INTENT(OUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(OUT) :: qout_tl(is-ng:ie+ng, js-ng:je+ng)
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
! local:
    REAL :: q1(npx), q2(npy)
    REAL :: q1_tl(npx), q2_tl(npy)
    INTEGER :: i, j
    INTEGER :: is1, js1, is2, js2, ie1, je1
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC PRESENT
    IF (grid_type .LT. 3) THEN
      IF (1 .LT. is - 1) THEN
        is1 = is - 1
      ELSE
        is1 = 1
      END IF
      IF (1 .LT. js - 1) THEN
        js1 = js - 1
      ELSE
        js1 = 1
      END IF
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (2 .LT. js) THEN
        js2 = js
      ELSE
        js2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 1
      END IF
      IF (npy - 1 .GT. je + 1) THEN
        je1 = je + 1
      ELSE
        je1 = npy - 1
      END IF
      DO j=js2,je1
        DO i=is2,ie1
          qout_tl(i, j) = 0.25*(qin_tl(i-1, j-1)+qin_tl(i, j-1)+qin_tl(i&
&           -1, j)+qin_tl(i, j))
          qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin(i&
&           , j))
        END DO
      END DO
! Fix the 4 Corners:
      IF (sw_corner) THEN
        qout_tl(1, 1) = r3*(qin_tl(1, 1)+qin_tl(1, 0)+qin_tl(0, 1))
        qout(1, 1) = r3*(qin(1, 1)+qin(1, 0)+qin(0, 1))
      END IF
      IF (se_corner) THEN
        qout_tl(npx, 1) = r3*(qin_tl(npx-1, 1)+qin_tl(npx-1, 0)+qin_tl(&
&         npx, 1))
        qout(npx, 1) = r3*(qin(npx-1, 1)+qin(npx-1, 0)+qin(npx, 1))
      END IF
      IF (ne_corner) THEN
        qout_tl(npx, npy) = r3*(qin_tl(npx-1, npy-1)+qin_tl(npx, npy-1)+&
&         qin_tl(npx-1, npy))
        qout(npx, npy) = r3*(qin(npx-1, npy-1)+qin(npx, npy-1)+qin(npx-1&
&         , npy))
      END IF
      IF (nw_corner) THEN
        qout_tl(1, npy) = r3*(qin_tl(1, npy-1)+qin_tl(0, npy-1)+qin_tl(1&
&         , npy))
        qout(1, npy) = r3*(qin(1, npy-1)+qin(0, npy-1)+qin(1, npy))
      END IF
! *** West Edges:
      IF (is .EQ. 1) THEN
        q2_tl = 0.0_8
        DO j=js1,je1
          q2_tl(j) = 0.5*(qin_tl(0, j)+qin_tl(1, j))
          q2(j) = 0.5*(qin(0, j)+qin(1, j))
        END DO
        DO j=js2,je1
          qout_tl(1, j) = edge_w(j)*q2_tl(j-1) + (1.-edge_w(j))*q2_tl(j)
          qout(1, j) = edge_w(j)*q2(j-1) + (1.-edge_w(j))*q2(j)
        END DO
      ELSE
        q2_tl = 0.0_8
      END IF
! East Edges:
      IF (ie + 1 .EQ. npx) THEN
        DO j=js1,je1
          q2_tl(j) = 0.5*(qin_tl(npx-1, j)+qin_tl(npx, j))
          q2(j) = 0.5*(qin(npx-1, j)+qin(npx, j))
        END DO
        DO j=js2,je1
          qout_tl(npx, j) = edge_e(j)*q2_tl(j-1) + (1.-edge_e(j))*q2_tl(&
&           j)
          qout(npx, j) = edge_e(j)*q2(j-1) + (1.-edge_e(j))*q2(j)
        END DO
      END IF
! South Edges:
      IF (js .EQ. 1) THEN
        q1_tl = 0.0_8
        DO i=is1,ie1
          q1_tl(i) = 0.5*(qin_tl(i, 0)+qin_tl(i, 1))
          q1(i) = 0.5*(qin(i, 0)+qin(i, 1))
        END DO
        DO i=is2,ie1
          qout_tl(i, 1) = edge_s(i)*q1_tl(i-1) + (1.-edge_s(i))*q1_tl(i)
          qout(i, 1) = edge_s(i)*q1(i-1) + (1.-edge_s(i))*q1(i)
        END DO
      ELSE
        q1_tl = 0.0_8
      END IF
! North Edges:
      IF (je + 1 .EQ. npy) THEN
        DO i=is1,ie1
          q1_tl(i) = 0.5*(qin_tl(i, npy-1)+qin_tl(i, npy))
          q1(i) = 0.5*(qin(i, npy-1)+qin(i, npy))
        END DO
        DO i=is2,ie1
          qout_tl(i, npy) = edge_n(i)*q1_tl(i-1) + (1.-edge_n(i))*q1_tl(&
&           i)
          qout(i, npy) = edge_n(i)*q1(i-1) + (1.-edge_n(i))*q1(i)
        END DO
      END IF
    ELSE
      DO j=js,je+1
        DO i=is,ie+1
          qout_tl(i, j) = 0.25*(qin_tl(i-1, j-1)+qin_tl(i, j-1)+qin_tl(i&
&           -1, j)+qin_tl(i, j))
          qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin(i&
&           , j))
        END DO
      END DO
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=js,je+1
          DO i=is,ie+1
            qin_tl(i, j) = qout_tl(i, j)
            qin(i, j) = qout(i, j)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE A2B_ORD2_TLM
  
end module a2b_edge_tlm_mod
