module a2b_edge_adm_mod

  use fv_grid_utils_mod, only: edge_w, edge_e, edge_s, edge_n, sw_corner, se_corner,  &
                               nw_corner, ne_corner, van2
  use fv_grid_tools_mod, only: dxa, dya, grid_type

  implicit none

  real, parameter:: r3 = 1./3.

  private
  public :: a2b_ord2, a2b_ord4, a2b_ord2_adm, a2b_ord4_adm

contains

  SUBROUTINE A2B_ORD4_ADM(qin, qin_ad, qout, qout_ad, npx, npy, is, ie, &
&   js, je, ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL*8, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
    REAL*8 :: qin_ad(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL, INTENT(INOUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(INOUT) :: qout_ad(is-ng:ie+ng, js-ng:je+ng)
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
    REAL :: qx_ad(is:ie+1, js-ng:je+ng)
    REAL :: qy(is-ng:ie+ng, js:je+1)
    REAL :: qy_ad(is-ng:ie+ng, js:je+1)
    REAL :: qxx(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qxx_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy_ad(is-ng:ie+ng, js-ng:je+ng)
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
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_from0
    INTEGER :: ad_to0
    INTEGER :: ad_from1
    INTEGER :: ad_to1
    INTEGER :: ad_from2
    INTEGER :: ad_to2
    INTEGER :: ad_from3
    INTEGER :: ad_to3
    INTEGER :: branch
    INTEGER :: min9
    INTEGER :: min8
    INTEGER :: min7
    INTEGER :: min6
    INTEGER :: min5
    INTEGER :: min4
    INTEGER :: min3
    INTEGER :: min2
    INTEGER :: min1
    REAL*8 :: temp_ad9
    REAL*8 :: temp_ad8
    REAL*8 :: temp_ad7
    REAL*8 :: temp_ad6
    REAL*8 :: temp_ad5
    REAL*8 :: temp_ad4
    REAL*8 :: temp_ad3
    REAL*8 :: temp_ad2
    REAL*8 :: temp_ad1
    REAL*8 :: temp_ad0
    INTEGER :: min19
    REAL*8 :: temp_ad
    INTEGER :: min18
    INTEGER :: min17
    INTEGER :: min16
    INTEGER :: min15
    INTEGER :: min14
    INTEGER :: min13
    REAL :: temp_ad17
    INTEGER :: min12
    REAL :: temp_ad16
    INTEGER :: min11
    REAL :: temp_ad15
    INTEGER :: min10
    REAL*8 :: temp_ad14
    REAL*8 :: temp_ad13
    REAL*8 :: temp_ad12
    REAL*8 :: temp_ad11
    REAL*8 :: temp_ad10
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
    IF (grid_type .LT. 3) THEN
! Corners:
! 6-point formular:
      IF (sw_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (se_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ne_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (nw_corner) THEN
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      IF (1 .LT. js - 2) THEN
        max1 = js - 2
      ELSE
        max1 = 1
      END IF
      IF (npy - 1 .GT. je + 2) THEN
        min1 = je + 2
      ELSE
        min1 = npy - 1
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
        ad_from = max2
        i = min2 + 1
        CALL PUSHINTEGER4(i - 1)
        CALL PUSHINTEGER4(ad_from)
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
        IF (js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL2B(0)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(2)
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
        IF (js .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (je + 1 .EQ. npy) THEN
          CALL PUSHCONTROL2B(2)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(0)
      END IF
      IF (3 .LT. js) THEN
        max7 = js
      ELSE
        max7 = 3
      END IF
      IF (npy - 2 .GT. je + 1) THEN
        min7 = je + 1
      ELSE
        min7 = npy - 2
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
        ad_from0 = max8
        i = min8 + 1
        CALL PUSHINTEGER4(i - 1)
        CALL PUSHINTEGER4(ad_from0)
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
        IF (is .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (ie + 1 .EQ. npx) THEN
          CALL PUSHCONTROL2B(0)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(2)
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
        IF (is .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (ie + 1 .EQ. npx) THEN
          CALL PUSHCONTROL2B(2)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(0)
      END IF
      IF (3 .LT. js) THEN
        max13 = js
      ELSE
        max13 = 3
      END IF
      IF (npy - 2 .GT. je + 1) THEN
        min13 = je + 1
      ELSE
        min13 = npy - 2
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
        ad_from1 = max14
        i = min14 + 1
        CALL PUSHINTEGER4(i - 1)
        CALL PUSHINTEGER4(ad_from1)
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
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
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
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      IF (2 .LT. js) THEN
        max17 = js
      ELSE
        max17 = 2
      END IF
      IF (npy - 1 .GT. je + 1) THEN
        min17 = je + 1
      ELSE
        min17 = npy - 1
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
        ad_from2 = max18
        i = min18 + 1
        CALL PUSHINTEGER4(i - 1)
        CALL PUSHINTEGER4(ad_from2)
        IF (is .EQ. 1) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (ie + 1 .EQ. npx) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
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
        ad_from3 = max19
        i = min19 + 1
        CALL PUSHINTEGER4(i - 1)
        CALL PUSHINTEGER4(ad_from3)
      END DO
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            qout_ad(i, j) = qout_ad(i, j) + qin_ad(i, j)
            qin_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      qy_ad = 0.0_8
      qxx_ad = 0.0_8
      qyy_ad = 0.0_8
      DO j=min17,max17,-1
        CALL POPINTEGER4(ad_from3)
        CALL POPINTEGER4(ad_to3)
        DO i=ad_to3,ad_from3,-1
          qxx_ad(i, j) = qxx_ad(i, j) + 0.5*qout_ad(i, j)
          qyy_ad(i, j) = qyy_ad(i, j) + 0.5*qout_ad(i, j)
          qout_ad(i, j) = 0.0_8
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          qy_ad(npx-2, j) = qy_ad(npx-2, j) + c1*qyy_ad(npx-1, j)
          qy_ad(npx-1, j) = qy_ad(npx-1, j) + c1*qyy_ad(npx-1, j)
          qout_ad(npx, j) = qout_ad(npx, j) + c2*qyy_ad(npx-1, j)
          qyy_ad(npx-2, j) = qyy_ad(npx-2, j) + c2*qyy_ad(npx-1, j)
          qyy_ad(npx-1, j) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          qy_ad(1, j) = qy_ad(1, j) + c1*qyy_ad(2, j)
          qy_ad(2, j) = qy_ad(2, j) + c1*qyy_ad(2, j)
          qout_ad(1, j) = qout_ad(1, j) + c2*qyy_ad(2, j)
          qyy_ad(3, j) = qyy_ad(3, j) + c2*qyy_ad(2, j)
          qyy_ad(2, j) = 0.0_8
        END IF
        CALL POPINTEGER4(ad_from2)
        CALL POPINTEGER4(ad_to2)
        DO i=ad_to2,ad_from2,-1
          qy_ad(i-2, j) = qy_ad(i-2, j) + a2*qyy_ad(i, j)
          qy_ad(i+1, j) = qy_ad(i+1, j) + a2*qyy_ad(i, j)
          qy_ad(i-1, j) = qy_ad(i-1, j) + a1*qyy_ad(i, j)
          qy_ad(i, j) = qy_ad(i, j) + a1*qyy_ad(i, j)
          qyy_ad(i, j) = 0.0_8
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        qx_ad = 0.0_8
      ELSE
        qx_ad = 0.0_8
        DO i=min16,max16,-1
          qx_ad(i, npy-2) = qx_ad(i, npy-2) + c1*qxx_ad(i, npy-1)
          qx_ad(i, npy-1) = qx_ad(i, npy-1) + c1*qxx_ad(i, npy-1)
          qout_ad(i, npy) = qout_ad(i, npy) + c2*qxx_ad(i, npy-1)
          qxx_ad(i, npy-2) = qxx_ad(i, npy-2) + c2*qxx_ad(i, npy-1)
          qxx_ad(i, npy-1) = 0.0_8
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO i=min15,max15,-1
          qx_ad(i, 1) = qx_ad(i, 1) + c1*qxx_ad(i, 2)
          qx_ad(i, 2) = qx_ad(i, 2) + c1*qxx_ad(i, 2)
          qout_ad(i, 1) = qout_ad(i, 1) + c2*qxx_ad(i, 2)
          qxx_ad(i, 3) = qxx_ad(i, 3) + c2*qxx_ad(i, 2)
          qxx_ad(i, 2) = 0.0_8
        END DO
      END IF
      DO j=min13,max13,-1
        CALL POPINTEGER4(ad_from1)
        CALL POPINTEGER4(ad_to1)
        DO i=ad_to1,ad_from1,-1
          qx_ad(i, j-2) = qx_ad(i, j-2) + a2*qxx_ad(i, j)
          qx_ad(i, j+1) = qx_ad(i, j+1) + a2*qxx_ad(i, j)
          qx_ad(i, j-1) = qx_ad(i, j-1) + a1*qxx_ad(i, j)
          qx_ad(i, j) = qx_ad(i, j) + a1*qxx_ad(i, j)
          qxx_ad(i, j) = 0.0_8
        END DO
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .NE. 0) THEN
        IF (branch .NE. 1) THEN
          qy_ad(npx-2, npy) = qy_ad(npx-2, npy) + c1*qout_ad(npx-1, npy)
          qy_ad(npx-1, npy) = qy_ad(npx-1, npy) + c1*qout_ad(npx-1, npy)
          qout_ad(npx-2, npy) = qout_ad(npx-2, npy) + c2*qout_ad(npx-1, &
&           npy)
          qout_ad(npx, npy) = qout_ad(npx, npy) + c2*qout_ad(npx-1, npy)
          qout_ad(npx-1, npy) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          qy_ad(1, npy) = qy_ad(1, npy) + c1*qout_ad(2, npy)
          qy_ad(2, npy) = qy_ad(2, npy) + c1*qout_ad(2, npy)
          qout_ad(1, npy) = qout_ad(1, npy) + c2*qout_ad(2, npy)
          qout_ad(3, npy) = qout_ad(3, npy) + c2*qout_ad(2, npy)
          qout_ad(2, npy) = 0.0_8
        END IF
        DO i=min12,max12,-1
          qy_ad(i-2, npy) = qy_ad(i-2, npy) + a2*qout_ad(i, npy)
          qy_ad(i+1, npy) = qy_ad(i+1, npy) + a2*qout_ad(i, npy)
          qy_ad(i-1, npy) = qy_ad(i-1, npy) + a1*qout_ad(i, npy)
          qy_ad(i, npy) = qy_ad(i, npy) + a1*qout_ad(i, npy)
          qout_ad(i, npy) = 0.0_8
        END DO
        DO i=min11,max11,-1
          gratio = dya(i, npy-2)/dya(i, npy-1)
          temp_ad13 = qy_ad(i, npy-1)/(gratio*2.+2.)
          qin_ad(i, npy-2) = qin_ad(i, npy-2) + 3.*temp_ad13
          qy_ad(i, npy) = qy_ad(i, npy) - gratio*temp_ad13
          qy_ad(i, npy-2) = qy_ad(i, npy-2) - temp_ad13
          qy_ad(i, npy-1) = 0.0_8
          temp_ad14 = 0.5*qy_ad(i, npy)/(gratio+1.)
          qin_ad(i, npy-1) = qin_ad(i, npy-1) + (gratio+2.)*temp_ad14 + &
&           3.*gratio*temp_ad13
          qin_ad(i, npy) = qin_ad(i, npy) + (gratio+2.)*temp_ad14
          qin_ad(i, npy-2) = qin_ad(i, npy-2) - temp_ad14
          qin_ad(i, npy+1) = qin_ad(i, npy+1) - temp_ad14
          qy_ad(i, npy) = 0.0_8
        END DO
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        qy_ad(npx-2, 1) = qy_ad(npx-2, 1) + c1*qout_ad(npx-1, 1)
        qy_ad(npx-1, 1) = qy_ad(npx-1, 1) + c1*qout_ad(npx-1, 1)
        qout_ad(npx-2, 1) = qout_ad(npx-2, 1) + c2*qout_ad(npx-1, 1)
        qout_ad(npx, 1) = qout_ad(npx, 1) + c2*qout_ad(npx-1, 1)
        qout_ad(npx-1, 1) = 0.0_8
      ELSE IF (branch .NE. 1) THEN
        GOTO 100
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        qy_ad(1, 1) = qy_ad(1, 1) + c1*qout_ad(2, 1)
        qy_ad(2, 1) = qy_ad(2, 1) + c1*qout_ad(2, 1)
        qout_ad(1, 1) = qout_ad(1, 1) + c2*qout_ad(2, 1)
        qout_ad(3, 1) = qout_ad(3, 1) + c2*qout_ad(2, 1)
        qout_ad(2, 1) = 0.0_8
      END IF
      DO i=min10,max10,-1
        qy_ad(i-2, 1) = qy_ad(i-2, 1) + a2*qout_ad(i, 1)
        qy_ad(i+1, 1) = qy_ad(i+1, 1) + a2*qout_ad(i, 1)
        qy_ad(i-1, 1) = qy_ad(i-1, 1) + a1*qout_ad(i, 1)
        qy_ad(i, 1) = qy_ad(i, 1) + a1*qout_ad(i, 1)
        qout_ad(i, 1) = 0.0_8
      END DO
      DO i=min9,max9,-1
        gratio = dya(i, 2)/dya(i, 1)
        temp_ad11 = qy_ad(i, 2)/(gratio*2.+2.)
        qin_ad(i, 1) = qin_ad(i, 1) + 3.*gratio*temp_ad11
        qin_ad(i, 2) = qin_ad(i, 2) + 3.*temp_ad11
        qy_ad(i, 1) = qy_ad(i, 1) - gratio*temp_ad11
        qy_ad(i, 3) = qy_ad(i, 3) - temp_ad11
        qy_ad(i, 2) = 0.0_8
        temp_ad12 = 0.5*qy_ad(i, 1)/(gratio+1.)
        qin_ad(i, 0) = qin_ad(i, 0) + (gratio+2.)*temp_ad12
        qin_ad(i, 1) = qin_ad(i, 1) + (gratio+2.)*temp_ad12
        qin_ad(i, -1) = qin_ad(i, -1) - temp_ad12
        qin_ad(i, 2) = qin_ad(i, 2) - temp_ad12
        qy_ad(i, 1) = 0.0_8
      END DO
 100  DO j=min7,max7,-1
        CALL POPINTEGER4(ad_from0)
        CALL POPINTEGER4(ad_to0)
        DO i=ad_to0,ad_from0,-1
          qin_ad(i, j-2) = qin_ad(i, j-2) + b2*qy_ad(i, j)
          qin_ad(i, j+1) = qin_ad(i, j+1) + b2*qy_ad(i, j)
          qin_ad(i, j-1) = qin_ad(i, j-1) + b1*qy_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qy_ad(i, j)
          qy_ad(i, j) = 0.0_8
        END DO
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .NE. 0) THEN
        IF (branch .NE. 1) THEN
          qx_ad(npx, npy-2) = qx_ad(npx, npy-2) + c1*qout_ad(npx, npy-1)
          qx_ad(npx, npy-1) = qx_ad(npx, npy-1) + c1*qout_ad(npx, npy-1)
          qout_ad(npx, npy-2) = qout_ad(npx, npy-2) + c2*qout_ad(npx, &
&           npy-1)
          qout_ad(npx, npy) = qout_ad(npx, npy) + c2*qout_ad(npx, npy-1)
          qout_ad(npx, npy-1) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          qx_ad(npx, 1) = qx_ad(npx, 1) + c1*qout_ad(npx, 2)
          qx_ad(npx, 2) = qx_ad(npx, 2) + c1*qout_ad(npx, 2)
          qout_ad(npx, 1) = qout_ad(npx, 1) + c2*qout_ad(npx, 2)
          qout_ad(npx, 3) = qout_ad(npx, 3) + c2*qout_ad(npx, 2)
          qout_ad(npx, 2) = 0.0_8
        END IF
        DO j=min6,max6,-1
          qx_ad(npx, j-2) = qx_ad(npx, j-2) + a2*qout_ad(npx, j)
          qx_ad(npx, j+1) = qx_ad(npx, j+1) + a2*qout_ad(npx, j)
          qx_ad(npx, j-1) = qx_ad(npx, j-1) + a1*qout_ad(npx, j)
          qx_ad(npx, j) = qx_ad(npx, j) + a1*qout_ad(npx, j)
          qout_ad(npx, j) = 0.0_8
        END DO
        DO j=min5,max5,-1
          gratio = dxa(npx-2, j)/dxa(npx-1, j)
          temp_ad9 = qx_ad(npx-1, j)/(gratio*2.+2.)
          qin_ad(npx-2, j) = qin_ad(npx-2, j) + 3.*temp_ad9
          qx_ad(npx, j) = qx_ad(npx, j) - gratio*temp_ad9
          qx_ad(npx-2, j) = qx_ad(npx-2, j) - temp_ad9
          qx_ad(npx-1, j) = 0.0_8
          temp_ad10 = 0.5*qx_ad(npx, j)/(gratio+1.)
          qin_ad(npx-1, j) = qin_ad(npx-1, j) + (gratio+2.)*temp_ad10 + &
&           3.*gratio*temp_ad9
          qin_ad(npx, j) = qin_ad(npx, j) + (gratio+2.)*temp_ad10
          qin_ad(npx-2, j) = qin_ad(npx-2, j) - temp_ad10
          qin_ad(npx+1, j) = qin_ad(npx+1, j) - temp_ad10
          qx_ad(npx, j) = 0.0_8
        END DO
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        qx_ad(1, npy-2) = qx_ad(1, npy-2) + c1*qout_ad(1, npy-1)
        qx_ad(1, npy-1) = qx_ad(1, npy-1) + c1*qout_ad(1, npy-1)
        qout_ad(1, npy-2) = qout_ad(1, npy-2) + c2*qout_ad(1, npy-1)
        qout_ad(1, npy) = qout_ad(1, npy) + c2*qout_ad(1, npy-1)
        qout_ad(1, npy-1) = 0.0_8
      ELSE IF (branch .NE. 1) THEN
        GOTO 110
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        qx_ad(1, 1) = qx_ad(1, 1) + c1*qout_ad(1, 2)
        qx_ad(1, 2) = qx_ad(1, 2) + c1*qout_ad(1, 2)
        qout_ad(1, 1) = qout_ad(1, 1) + c2*qout_ad(1, 2)
        qout_ad(1, 3) = qout_ad(1, 3) + c2*qout_ad(1, 2)
        qout_ad(1, 2) = 0.0_8
      END IF
      DO j=min4,max4,-1
        qx_ad(1, j-2) = qx_ad(1, j-2) + a2*qout_ad(1, j)
        qx_ad(1, j+1) = qx_ad(1, j+1) + a2*qout_ad(1, j)
        qx_ad(1, j-1) = qx_ad(1, j-1) + a1*qout_ad(1, j)
        qx_ad(1, j) = qx_ad(1, j) + a1*qout_ad(1, j)
        qout_ad(1, j) = 0.0_8
      END DO
      DO j=min3,max3,-1
        gratio = dxa(2, j)/dxa(1, j)
        temp_ad7 = qx_ad(2, j)/(gratio*2.+2.)
        qin_ad(1, j) = qin_ad(1, j) + 3.*gratio*temp_ad7
        qin_ad(2, j) = qin_ad(2, j) + 3.*temp_ad7
        qx_ad(1, j) = qx_ad(1, j) - gratio*temp_ad7
        qx_ad(3, j) = qx_ad(3, j) - temp_ad7
        qx_ad(2, j) = 0.0_8
        temp_ad8 = 0.5*qx_ad(1, j)/(gratio+1.)
        qin_ad(0, j) = qin_ad(0, j) + (gratio+2.)*temp_ad8
        qin_ad(1, j) = qin_ad(1, j) + (gratio+2.)*temp_ad8
        qin_ad(-1, j) = qin_ad(-1, j) - temp_ad8
        qin_ad(2, j) = qin_ad(2, j) - temp_ad8
        qx_ad(1, j) = 0.0_8
      END DO
 110  DO j=min1,max1,-1
        CALL POPINTEGER4(ad_from)
        CALL POPINTEGER4(ad_to)
        DO i=ad_to,ad_from,-1
          qin_ad(i-2, j) = qin_ad(i-2, j) + b2*qx_ad(i, j)
          qin_ad(i+1, j) = qin_ad(i+1, j) + b2*qx_ad(i, j)
          qin_ad(i-1, j) = qin_ad(i-1, j) + b1*qx_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qx_ad(i, j)
          qx_ad(i, j) = 0.0_8
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        temp_ad5 = d1*qout_ad(1, npy)
        temp_ad6 = d2*qout_ad(1, npy)
        qin_ad(0, npy-1) = qin_ad(0, npy-1) + temp_ad5
        qin_ad(1, npy-1) = qin_ad(1, npy-1) + temp_ad5
        qin_ad(1, npy) = qin_ad(1, npy) + temp_ad5
        qin_ad(-1, npy-2) = qin_ad(-1, npy-2) + temp_ad6
        qin_ad(2, npy-2) = qin_ad(2, npy-2) + temp_ad6
        qin_ad(2, npy+1) = qin_ad(2, npy+1) + temp_ad6
        qout_ad(1, npy) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        temp_ad3 = d1*qout_ad(npx, npy)
        temp_ad4 = d2*qout_ad(npx, npy)
        qin_ad(npx-1, npy-1) = qin_ad(npx-1, npy-1) + temp_ad3
        qin_ad(npx, npy-1) = qin_ad(npx, npy-1) + temp_ad3
        qin_ad(npx-1, npy) = qin_ad(npx-1, npy) + temp_ad3
        qin_ad(npx-2, npy-2) = qin_ad(npx-2, npy-2) + temp_ad4
        qin_ad(npx+1, npy-2) = qin_ad(npx+1, npy-2) + temp_ad4
        qin_ad(npx-2, npy+1) = qin_ad(npx-2, npy+1) + temp_ad4
        qout_ad(npx, npy) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        temp_ad1 = d1*qout_ad(npx, 1)
        temp_ad2 = d2*qout_ad(npx, 1)
        qin_ad(npx-1, 0) = qin_ad(npx-1, 0) + temp_ad1
        qin_ad(npx-1, 1) = qin_ad(npx-1, 1) + temp_ad1
        qin_ad(npx, 1) = qin_ad(npx, 1) + temp_ad1
        qin_ad(npx-2, -1) = qin_ad(npx-2, -1) + temp_ad2
        qin_ad(npx-2, 2) = qin_ad(npx-2, 2) + temp_ad2
        qin_ad(npx+1, 2) = qin_ad(npx+1, 2) + temp_ad2
        qout_ad(npx, 1) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        temp_ad = d1*qout_ad(1, 1)
        temp_ad0 = d2*qout_ad(1, 1)
        qin_ad(1, 0) = qin_ad(1, 0) + temp_ad
        qin_ad(0, 1) = qin_ad(0, 1) + temp_ad
        qin_ad(1, 1) = qin_ad(1, 1) + temp_ad
        qin_ad(2, -1) = qin_ad(2, -1) + temp_ad0
        qin_ad(-1, 2) = qin_ad(-1, 2) + temp_ad0
        qin_ad(2, 2) = qin_ad(2, 2) + temp_ad0
        qout_ad(1, 1) = 0.0_8
      END IF
    ELSE
      qx_ad = 0.0_8
      qy_ad = 0.0_8
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          temp_ad15 = 0.5*qout_ad(i, j)
          temp_ad16 = a1*temp_ad15
          temp_ad17 = a2*temp_ad15
          qx_ad(i, j-1) = qx_ad(i, j-1) + temp_ad16
          qx_ad(i, j) = qx_ad(i, j) + temp_ad16
          qy_ad(i-1, j) = qy_ad(i-1, j) + temp_ad16
          qy_ad(i, j) = qy_ad(i, j) + temp_ad16
          qx_ad(i, j-2) = qx_ad(i, j-2) + temp_ad17
          qx_ad(i, j+1) = qx_ad(i, j+1) + temp_ad17
          qy_ad(i-2, j) = qy_ad(i-2, j) + temp_ad17
          qy_ad(i+1, j) = qy_ad(i+1, j) + temp_ad17
          qout_ad(i, j) = 0.0_8
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+2,is-2,-1
          qin_ad(i, j-1) = qin_ad(i, j-1) + b1*qy_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qy_ad(i, j)
          qin_ad(i, j-2) = qin_ad(i, j-2) + b2*qy_ad(i, j)
          qin_ad(i, j+1) = qin_ad(i, j+1) + b2*qy_ad(i, j)
          qy_ad(i, j) = 0.0_8
        END DO
      END DO
      DO j=je+2,js-2,-1
        DO i=ie+1,is,-1
          qin_ad(i-1, j) = qin_ad(i-1, j) + b1*qx_ad(i, j)
          qin_ad(i, j) = qin_ad(i, j) + b1*qx_ad(i, j)
          qin_ad(i-2, j) = qin_ad(i-2, j) + b2*qx_ad(i, j)
          qin_ad(i+1, j) = qin_ad(i+1, j) + b2*qx_ad(i, j)
          qx_ad(i, j) = 0.0_8
        END DO
      END DO
    END IF
  END SUBROUTINE A2B_ORD4_ADM
  SUBROUTINE A2B_ORD4(qin, qout, npx, npy, is, ie, js, je, ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL*8, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL, INTENT(INOUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
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
    REAL :: qy(is-ng:ie+ng, js:je+1)
    REAL :: qxx(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qyy(is-ng:ie+ng, js-ng:je+ng)
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
      IF (sw_corner) qout(1, 1) = d1*(qin(1, 0)+qin(0, 1)+qin(1, 1)) + &
&         d2*(qin(2, -1)+qin(-1, 2)+qin(2, 2))
      IF (se_corner) qout(npx, 1) = d1*(qin(npx-1, 0)+qin(npx-1, 1)+qin(&
&         npx, 1)) + d2*(qin(npx-2, -1)+qin(npx-2, 2)+qin(npx+1, 2))
      IF (ne_corner) qout(npx, npy) = d1*(qin(npx-1, npy-1)+qin(npx, npy&
&         -1)+qin(npx-1, npy)) + d2*(qin(npx-2, npy-2)+qin(npx+1, npy-2)&
&         +qin(npx-2, npy+1))
      IF (nw_corner) qout(1, npy) = d1*(qin(0, npy-1)+qin(1, npy-1)+qin(&
&         1, npy)) + d2*(qin(-1, npy-2)+qin(2, npy-2)+qin(2, npy+1))
      IF (1 .LT. js - 2) THEN
        max1 = js - 2
      ELSE
        max1 = 1
      END IF
      IF (npy - 1 .GT. je + 2) THEN
        min1 = je + 2
      ELSE
        min1 = npy - 1
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
          qx(1, j) = 0.5*((2.+gratio)*(qin(0, j)+qin(1, j))-(qin(-1, j)+&
&           qin(2, j)))/(1.+gratio)
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
          qout(1, j) = a2*(qx(1, j-2)+qx(1, j+1)) + a1*(qx(1, j-1)+qx(1&
&           , j))
        END DO
        IF (js .EQ. 1) qout(1, 2) = c1*(qx(1, 1)+qx(1, 2)) + c2*(qout(1&
&           , 1)+qout(1, 3))
        IF (je + 1 .EQ. npy) qout(1, npy-1) = c1*(qx(1, npy-2)+qx(1, npy&
&           -1)) + c2*(qout(1, npy-2)+qout(1, npy))
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
          qx(npx, j) = 0.5*((2.+gratio)*(qin(npx-1, j)+qin(npx, j))-(qin&
&           (npx-2, j)+qin(npx+1, j)))/(1.+gratio)
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
          qout(npx, j) = a2*(qx(npx, j-2)+qx(npx, j+1)) + a1*(qx(npx, j-&
&           1)+qx(npx, j))
        END DO
        IF (js .EQ. 1) qout(npx, 2) = c1*(qx(npx, 1)+qx(npx, 2)) + c2*(&
&           qout(npx, 1)+qout(npx, 3))
        IF (je + 1 .EQ. npy) qout(npx, npy-1) = c1*(qx(npx, npy-2)+qx(&
&           npx, npy-1)) + c2*(qout(npx, npy-2)+qout(npx, npy))
      END IF
      IF (3 .LT. js) THEN
        max7 = js
      ELSE
        max7 = 3
      END IF
      IF (npy - 2 .GT. je + 1) THEN
        min7 = je + 1
      ELSE
        min7 = npy - 2
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
          qy(i, 1) = 0.5*((2.+gratio)*(qin(i, 0)+qin(i, 1))-(qin(i, -1)+&
&           qin(i, 2)))/(1.+gratio)
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
          qout(i, 1) = a2*(qy(i-2, 1)+qy(i+1, 1)) + a1*(qy(i-1, 1)+qy(i&
&           , 1))
        END DO
        IF (is .EQ. 1) qout(2, 1) = c1*(qy(1, 1)+qy(2, 1)) + c2*(qout(1&
&           , 1)+qout(3, 1))
        IF (ie + 1 .EQ. npx) qout(npx-1, 1) = c1*(qy(npx-2, 1)+qy(npx-1&
&           , 1)) + c2*(qout(npx-2, 1)+qout(npx, 1))
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
          qy(i, npy) = 0.5*((2.+gratio)*(qin(i, npy-1)+qin(i, npy))-(qin&
&           (i, npy-2)+qin(i, npy+1)))/(1.+gratio)
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
          qout(i, npy) = a2*(qy(i-2, npy)+qy(i+1, npy)) + a1*(qy(i-1, &
&           npy)+qy(i, npy))
        END DO
        IF (is .EQ. 1) qout(2, npy) = c1*(qy(1, npy)+qy(2, npy)) + c2*(&
&           qout(1, npy)+qout(3, npy))
        IF (ie + 1 .EQ. npx) qout(npx-1, npy) = c1*(qy(npx-2, npy)+qy(&
&           npx-1, npy)) + c2*(qout(npx-2, npy)+qout(npx, npy))
      END IF
      IF (3 .LT. js) THEN
        max13 = js
      ELSE
        max13 = 3
      END IF
      IF (npy - 2 .GT. je + 1) THEN
        min13 = je + 1
      ELSE
        min13 = npy - 2
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
      ELSE
        min17 = npy - 1
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
          qyy(i, j) = a2*(qy(i-2, j)+qy(i+1, j)) + a1*(qy(i-1, j)+qy(i, &
&           j))
        END DO
        IF (is .EQ. 1) qyy(2, j) = c1*(qy(1, j)+qy(2, j)) + c2*(qout(1, &
&           j)+qyy(3, j))
        IF (ie + 1 .EQ. npx) qyy(npx-1, j) = c1*(qy(npx-2, j)+qy(npx-1, &
&           j)) + c2*(qout(npx, j)+qyy(npx-2, j))
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
          qout(i, j) = 0.5*(qxx(i, j)+qyy(i, j))
        END DO
      END DO
    ELSE
! grid_type>=3
!------------------------
! Doubly periodic domain:
!------------------------
! X-sweep: PPM
      DO j=js-2,je+2
        DO i=is,ie+1
          qx(i, j) = b1*(qin(i-1, j)+qin(i, j)) + b2*(qin(i-2, j)+qin(i+&
&           1, j))
        END DO
      END DO
! Y-sweep: PPM
      DO j=js,je+1
        DO i=is-2,ie+2
          qy(i, j) = b1*(qin(i, j-1)+qin(i, j)) + b2*(qin(i, j-2)+qin(i&
&           , j+1))
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie+1
          qout(i, j) = 0.5*(a1*(qx(i, j-1)+qx(i, j)+qy(i-1, j)+qy(i, j))&
&           +a2*(qx(i, j-2)+qx(i, j+1)+qy(i-2, j)+qy(i+1, j)))
        END DO
      END DO
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=js,je+1
          DO i=is,ie+1
            qin(i, j) = qout(i, j)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE A2B_ORD4
!  Differentiation of a2b_ord2 in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: qin qout
!   with respect to varying inputs: qin qout
  SUBROUTINE A2B_ORD2_ADM(qin, qin_ad, qout, qout_ad, npx, npy, is, ie, &
&   js, je, ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL*8, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
    REAL*8 :: qin_ad(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL :: qout(is-ng:ie+ng, js-ng:je+ng)
    REAL :: qout_ad(is-ng:ie+ng, js-ng:je+ng)
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
! local:
    REAL :: q1(npx), q2(npy)
    REAL :: q1_ad(npx), q2_ad(npy)
    INTEGER :: i, j
    INTEGER :: is1, js1, is2, js2, ie1, je1
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC PRESENT
    INTEGER :: branch
    REAL*8 :: temp_ad4
    REAL*8 :: temp_ad3
    REAL*8 :: temp_ad2
    REAL*8 :: temp_ad1
    REAL*8 :: temp_ad0
    REAL*8 :: temp_ad
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
! Fix the 4 Corners:
      IF (sw_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (se_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ne_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (nw_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
! *** West Edges:
      IF (is .EQ. 1) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
! East Edges:
      IF (ie + 1 .EQ. npx) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
! South Edges:
      IF (js .EQ. 1) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
! North Edges:
      IF (je + 1 .EQ. npy) THEN
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            qout_ad(i, j) = qout_ad(i, j) + qin_ad(i, j)
            qin_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
    END IF
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      q1_ad = 0.0_8
      DO i=ie1,is2,-1
        q1_ad(i-1) = q1_ad(i-1) + edge_n(i)*qout_ad(i, npy)
        q1_ad(i) = q1_ad(i) + (1.-edge_n(i))*qout_ad(i, npy)
        qout_ad(i, npy) = 0.0_8
      END DO
      DO i=ie1,is1,-1
        qin_ad(i, npy-1) = qin_ad(i, npy-1) + 0.5*q1_ad(i)
        qin_ad(i, npy) = qin_ad(i, npy) + 0.5*q1_ad(i)
        q1_ad(i) = 0.0_8
      END DO
    ELSE IF (branch .EQ. 1) THEN
      q1_ad = 0.0_8
    ELSE
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          temp_ad4 = 0.25*qout_ad(i, j)
          qin_ad(i-1, j-1) = qin_ad(i-1, j-1) + temp_ad4
          qin_ad(i, j-1) = qin_ad(i, j-1) + temp_ad4
          qin_ad(i-1, j) = qin_ad(i-1, j) + temp_ad4
          qin_ad(i, j) = qin_ad(i, j) + temp_ad4
          qout_ad(i, j) = 0.0_8
        END DO
      END DO
      GOTO 100
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO i=ie1,is2,-1
        q1_ad(i-1) = q1_ad(i-1) + edge_s(i)*qout_ad(i, 1)
        q1_ad(i) = q1_ad(i) + (1.-edge_s(i))*qout_ad(i, 1)
        qout_ad(i, 1) = 0.0_8
      END DO
      DO i=ie1,is1,-1
        qin_ad(i, 0) = qin_ad(i, 0) + 0.5*q1_ad(i)
        qin_ad(i, 1) = qin_ad(i, 1) + 0.5*q1_ad(i)
        q1_ad(i) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      q2_ad = 0.0_8
      DO j=je1,js2,-1
        q2_ad(j-1) = q2_ad(j-1) + edge_e(j)*qout_ad(npx, j)
        q2_ad(j) = q2_ad(j) + (1.-edge_e(j))*qout_ad(npx, j)
        qout_ad(npx, j) = 0.0_8
      END DO
      DO j=je1,js1,-1
        qin_ad(npx-1, j) = qin_ad(npx-1, j) + 0.5*q2_ad(j)
        qin_ad(npx, j) = qin_ad(npx, j) + 0.5*q2_ad(j)
        q2_ad(j) = 0.0_8
      END DO
    ELSE
      q2_ad = 0.0_8
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=je1,js2,-1
        q2_ad(j-1) = q2_ad(j-1) + edge_w(j)*qout_ad(1, j)
        q2_ad(j) = q2_ad(j) + (1.-edge_w(j))*qout_ad(1, j)
        qout_ad(1, j) = 0.0_8
      END DO
      DO j=je1,js1,-1
        qin_ad(0, j) = qin_ad(0, j) + 0.5*q2_ad(j)
        qin_ad(1, j) = qin_ad(1, j) + 0.5*q2_ad(j)
        q2_ad(j) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      temp_ad3 = r3*qout_ad(1, npy)
      qin_ad(1, npy-1) = qin_ad(1, npy-1) + temp_ad3
      qin_ad(0, npy-1) = qin_ad(0, npy-1) + temp_ad3
      qin_ad(1, npy) = qin_ad(1, npy) + temp_ad3
      qout_ad(1, npy) = 0.0_8
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      temp_ad2 = r3*qout_ad(npx, npy)
      qin_ad(npx-1, npy-1) = qin_ad(npx-1, npy-1) + temp_ad2
      qin_ad(npx, npy-1) = qin_ad(npx, npy-1) + temp_ad2
      qin_ad(npx-1, npy) = qin_ad(npx-1, npy) + temp_ad2
      qout_ad(npx, npy) = 0.0_8
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      temp_ad1 = r3*qout_ad(npx, 1)
      qin_ad(npx-1, 1) = qin_ad(npx-1, 1) + temp_ad1
      qin_ad(npx-1, 0) = qin_ad(npx-1, 0) + temp_ad1
      qin_ad(npx, 1) = qin_ad(npx, 1) + temp_ad1
      qout_ad(npx, 1) = 0.0_8
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      temp_ad0 = r3*qout_ad(1, 1)
      qin_ad(1, 1) = qin_ad(1, 1) + temp_ad0
      qin_ad(1, 0) = qin_ad(1, 0) + temp_ad0
      qin_ad(0, 1) = qin_ad(0, 1) + temp_ad0
      qout_ad(1, 1) = 0.0_8
    END IF
    DO j=je1,js2,-1
      DO i=ie1,is2,-1
        temp_ad = 0.25*qout_ad(i, j)
        qin_ad(i-1, j-1) = qin_ad(i-1, j-1) + temp_ad
        qin_ad(i, j-1) = qin_ad(i, j-1) + temp_ad
        qin_ad(i-1, j) = qin_ad(i-1, j) + temp_ad
        qin_ad(i, j) = qin_ad(i, j) + temp_ad
        qout_ad(i, j) = 0.0_8
      END DO
    END DO
 100 CONTINUE
  END SUBROUTINE A2B_ORD2_ADM
  SUBROUTINE A2B_ORD2(qin, qout, npx, npy, is, ie, js, je, ng, replace)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL*8, INTENT(INOUT) :: qin(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL, INTENT(OUT) :: qout(is-ng:ie+ng, js-ng:je+ng)
    LOGICAL, OPTIONAL, INTENT(IN) :: replace
! local:
    REAL :: q1(npx), q2(npy)
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
          qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin(i&
&           , j))
        END DO
      END DO
! Fix the 4 Corners:
      IF (sw_corner) qout(1, 1) = r3*(qin(1, 1)+qin(1, 0)+qin(0, 1))
      IF (se_corner) qout(npx, 1) = r3*(qin(npx-1, 1)+qin(npx-1, 0)+qin(&
&         npx, 1))
      IF (ne_corner) qout(npx, npy) = r3*(qin(npx-1, npy-1)+qin(npx, npy&
&         -1)+qin(npx-1, npy))
      IF (nw_corner) qout(1, npy) = r3*(qin(1, npy-1)+qin(0, npy-1)+qin(&
&         1, npy))
! *** West Edges:
      IF (is .EQ. 1) THEN
        DO j=js1,je1
          q2(j) = 0.5*(qin(0, j)+qin(1, j))
        END DO
        DO j=js2,je1
          qout(1, j) = edge_w(j)*q2(j-1) + (1.-edge_w(j))*q2(j)
        END DO
      END IF
! East Edges:
      IF (ie + 1 .EQ. npx) THEN
        DO j=js1,je1
          q2(j) = 0.5*(qin(npx-1, j)+qin(npx, j))
        END DO
        DO j=js2,je1
          qout(npx, j) = edge_e(j)*q2(j-1) + (1.-edge_e(j))*q2(j)
        END DO
      END IF
! South Edges:
      IF (js .EQ. 1) THEN
        DO i=is1,ie1
          q1(i) = 0.5*(qin(i, 0)+qin(i, 1))
        END DO
        DO i=is2,ie1
          qout(i, 1) = edge_s(i)*q1(i-1) + (1.-edge_s(i))*q1(i)
        END DO
      END IF
! North Edges:
      IF (je + 1 .EQ. npy) THEN
        DO i=is1,ie1
          q1(i) = 0.5*(qin(i, npy-1)+qin(i, npy))
        END DO
        DO i=is2,ie1
          qout(i, npy) = edge_n(i)*q1(i-1) + (1.-edge_n(i))*q1(i)
        END DO
      END IF
    ELSE
      DO j=js,je+1
        DO i=is,ie+1
          qout(i, j) = 0.25*(qin(i-1, j-1)+qin(i, j-1)+qin(i-1, j)+qin(i&
&           , j))
        END DO
      END DO
    END IF
    IF (PRESENT(replace)) THEN
      IF (replace) THEN
        DO j=js,je+1
          DO i=is,ie+1
            qin(i, j) = qout(i, j)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE A2B_ORD2
  
end module a2b_edge_adm_mod
