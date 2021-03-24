!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.7 (r4786) - 21 Feb 2013 15:53
!
!  Differentiation of radcouple in reverse (adjoint) mode:
!   gradient     of useful results: rad_cf rad_qi rad_ql rad_ri
!                rad_rl
!   with respect to varying inputs: qcian af rad_cf qclls qcils
!                rad_qi rad_ql cf qclan rad_ri rad_rl te
!   RW status of diff variables: qcian:out af:out rad_cf:in-zero
!                qclls:out qcils:out rad_qi:in-zero rad_ql:in-zero
!                cf:out qclan:out rad_ri:in-zero rad_rl:in-zero
!                te:out
SUBROUTINE RADCOUPLE_B(te, teb, pl, cf, cfb, af, afb, qclls, qcllsb, &
&  qcils, qcilsb, qclan, qclanb, qcian, qcianb, rad_ql, rad_qlb, rad_qi, &
&  rad_qib, rad_cf, rad_cfb, rad_rl, rad_rlb, rad_ri, rad_rib, tempor)
  IMPLICIT NONE
!Inputs
  REAL, INTENT(IN) :: te, pl, tempor
  REAL :: teb
  REAL, INTENT(IN) :: af, cf, qclan, qcian, qclls, qcils
  REAL :: afb, cfb, qclanb, qcianb, qcllsb, qcilsb
! real, intent(in ) :: QRN_ALL, QSN_ALL
!Outputs
  REAL :: rad_ql, rad_qi, rad_cf, rad_rl, rad_ri
  REAL :: rad_qlb, rad_qib, rad_cfb, rad_rlb, rad_rib
! real, intent(out) :: RAD_QR,RAD_QS
!Locals
  REAL :: ss, rad_ri_an, afx, alph
  REAL :: ssb, rad_ri_anb, afxb
  REAL, PARAMETER :: min_ri=20.e-6, max_ri=40.e-6, ri_anv=30.e-6
  INTEGER :: branch
  REAL :: min3
  REAL :: min2
  REAL :: max2b
  REAL :: min1
  REAL :: tempb2
  REAL :: tempb1
  REAL :: tempb0
  INTRINSIC MAX
  REAL :: min1b
  REAL :: x2
  REAL :: x1
  REAL :: x2b
  REAL :: tempb
  REAL :: max3b
  INTRINSIC MIN
  REAL :: temp
  REAL :: max3
  REAL :: max2
  REAL :: max1
  REAL :: min2b
!Initialize outputs
  rad_ri = 0.0
!RAD_QR = 0.0
!RAD_QS = 0.0
! Adjust Anvil fractions for warm clouds
  alph = 0.1
  ss = (280.-te)/20.
  IF (1.0 .GT. ss) THEN
    CALL PUSHCONTROL1B(0)
    ss = ss
  ELSE
    ss = 1.0
    CALL PUSHCONTROL1B(1)
  END IF
  IF (0.0 .LT. ss) THEN
    CALL PUSHCONTROL1B(0)
    ss = ss
  ELSE
    ss = 0.0
    CALL PUSHCONTROL1B(1)
  END IF
  CALL PUSHREAL4(ss)
  ss = alph + ss**3*(1.0-alph)
  afx = af*ss*0.5
  IF (cf + afx .GT. 1.00) THEN
    rad_cf = 1.00
    CALL PUSHCONTROL1B(0)
  ELSE
    rad_cf = cf + afx
    CALL PUSHCONTROL1B(1)
  END IF
!Total In-cloud liquid
  IF (rad_cf .GT. 10.0e-8) THEN
!0 -> 10e-8 FOR LINEARIZATION PROTECTION
    rad_ql = (qclls+qclan)/rad_cf
    CALL PUSHCONTROL1B(1)
  ELSE
    rad_ql = 0.0
    CALL PUSHCONTROL1B(0)
  END IF
  IF (rad_ql .GT. 0.01) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
! Total In-cloud ice
  IF (rad_cf .GT. 10.0e-8) THEN
!0 -> 10e-8 FOR LINEARIZATION PROTECTION
    rad_qi = (qcils+qcian)/rad_cf
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
    rad_qi = 0.0
  END IF
  IF (rad_qi .GT. 0.01) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
! Total In-cloud precipitation
! if (  RAD_CF >0. ) then
!    RAD_QR = ( QRN_ALL ) / RAD_CF
!    RAD_QS = ( QSN_ALL ) / RAD_CF
! else
!    RAD_QR = 0.0
!    RAD_QS = 0.0
! end if
! RAD_QR = MIN( RAD_QR, 0.01 )
! RAD_QS = MIN( RAD_QS, 0.01 )
  IF (pl .LT. 150.) rad_ri = max_ri
  IF (pl .GE. 150.) rad_ri = max_ri*150./pl
! Weigh in a separate R_ice for Anvil Ice according to
  rad_ri_an = rad_ri
  IF (qcils + qcian .GT. 0.0) THEN
    IF (qcils/rad_ri + qcian/ri_anv .GT. 10e-8) THEN
!LINEARIZATION PROTECTION
      rad_ri_an = (qcils+qcian)/(qcils/rad_ri+qcian/ri_anv)
      CALL PUSHCONTROL2B(2)
    ELSE
      CALL PUSHCONTROL2B(1)
    END IF
  ELSE
    CALL PUSHCONTROL2B(0)
  END IF
  IF (rad_ri .GT. rad_ri_an) THEN
    CALL PUSHREAL4(rad_ri)
    rad_ri = rad_ri_an
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHREAL4(rad_ri)
    rad_ri = rad_ri
    CALL PUSHCONTROL1B(1)
  END IF
  IF (rad_ri .LT. min_ri) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 825. .AND. te .LE. 282. .AND. tempor .EQ. 1.) THEN
    IF (0.71*te - 190.25 .LT. 5.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 775. .AND. pl .LT. 825. .AND. te .LE. 282. .AND. te .GT. &
&      275. .AND. tempor .EQ. 1.) THEN
    IF (-(0.1*pl) + 0.71*te - 107.75 .GT. 10.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 825. .AND. te .LE. 275. .AND. tempor .EQ. 1.) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
! Thin low tropical clouds
  IF (pl .GE. 950. .AND. te .GE. 285.) THEN
    IF (2.2*te - 617. .GT. 21.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 925. .AND. te .GE. 290.) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 925. .AND. pl .LT. 950. .AND. te .GT. 285. .AND. te .LT. &
&      290.) THEN
    IF (0.44*pl + 2.2*te - 1035. .GT. 21.) THEN
      x2 = 21.
      CALL PUSHCONTROL1B(0)
    ELSE
      x2 = 0.44*pl + 2.2*te - 1035.
      CALL PUSHCONTROL1B(1)
    END IF
    IF (x2 .LT. 10.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 950. .AND. te .GE. 290.) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (rad_cf .LT. 1.e-5) THEN
    rad_cfb = 0.0
    rad_qib = 0.0
    rad_qlb = 0.0
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rlb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    max3b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      x2b = 0.0
    ELSE
      x2b = max3b
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      teb = 0.0
    ELSE
      teb = 2.2*x2b
    END IF
    rad_rlb = 0.0
  ELSE
    teb = 0.0
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rlb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    min2b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) teb = teb + 2.2*min2b
    rad_rlb = 0.0
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rlb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    min1b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) teb = teb + 0.71*min1b
    rad_rlb = 0.0
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    max2b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) teb = teb + 0.71*max2b
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rib = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    CALL POPREAL4(rad_ri)
    rad_ri_anb = rad_rib
  ELSE
    CALL POPREAL4(rad_ri)
    rad_ri_anb = 0.0
  END IF
  CALL POPCONTROL2B(branch)
  IF (branch .EQ. 0) THEN
    qcianb = 0.0
    qcilsb = 0.0
  ELSE IF (branch .EQ. 1) THEN
    qcianb = 0.0
    qcilsb = 0.0
  ELSE
    temp = qcils/rad_ri + qcian/ri_anv
    tempb1 = rad_ri_anb/temp
    tempb2 = -((qcils+qcian)*tempb1/temp)
    qcilsb = tempb2/rad_ri + tempb1
    qcianb = tempb2/ri_anv + tempb1
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_qib = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) THEN
    tempb0 = rad_qib/rad_cf
    qcilsb = qcilsb + tempb0
    qcianb = qcianb + tempb0
    rad_cfb = rad_cfb - (qcils+qcian)*tempb0/rad_cf
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_qlb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    qcllsb = 0.0
    qclanb = 0.0
  ELSE
    tempb = rad_qlb/rad_cf
    qcllsb = tempb
    qclanb = tempb
    rad_cfb = rad_cfb - (qclls+qclan)*tempb/rad_cf
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    cfb = 0.0
    afxb = 0.0
  ELSE
    cfb = rad_cfb
    afxb = rad_cfb
  END IF
  afb = 0.5*ss*afxb
  ssb = 0.5*af*afxb
  CALL POPREAL4(ss)
  ssb = (1.0-alph)*3*ss**2*ssb
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) ssb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) ssb = 0.0
  teb = teb - ssb/20.
  rad_cfb = 0.0
  rad_qib = 0.0
  rad_qlb = 0.0
  rad_rib = 0.0
  rad_rlb = 0.0
END SUBROUTINE RADCOUPLE_B

