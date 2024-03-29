!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5903) - 14 Dec 2015 10:32
!
MODULE QSAT_UTIL_D
!This module contains the subroutines used to generate the lookup table for computing saturation vapour pressure.
!Variables designated as real*8 should be kept at least as so, even if ESTBLX is required in real*4.
  USE MAPL_CONSTANTSMOD
  IMPLICIT NONE
  PRIVATE 
  INTEGER, PARAMETER :: degsubs=100
  REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
  INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
  PUBLIC esinit, dqsatpert, dqsatscapert, degsubs, tmintbl, tablesize, &
& estblx
  PUBLIC dqsatpert_tlm, dqsatscapert_tlm
  REAL*8, DIMENSION(tablesize) :: estblx

CONTAINS
  SUBROUTINE ESINIT()
    IMPLICIT NONE
!LOCALS
    REAL*8, PARAMETER :: zeroc=273.16, tmix=-20.0
    REAL*8, DIMENSION(tablesize) :: estble, estblw
    INTEGER :: i
    REAL*8 :: t, delta_t
    delta_t = 1.0/degsubs
    DO i=1,tablesize
      t = (i-1)*delta_t + tmintbl
      IF (t .GT. zeroc) THEN
        CALL QSATLQU0(estble(i), t)
      ELSE
        CALL QSATICE0(estble(i), t)
      END IF
      CALL QSATLQU0(estblw(i), t)
      t = t - zeroc
      IF (t .GE. tmix .AND. t .LT. 0.0) THEN
        estblx(i) = t/tmix*(estble(i)-estblw(i)) + estblw(i)
      ELSE
        estblx(i) = estble(i)
      END IF
    END DO
  END SUBROUTINE ESINIT
  SUBROUTINE QSATLQU0(qs, tl)
    IMPLICIT NONE
!INPUTS
!, TMAXTBL
    REAL*8 :: tl
!OUTPUTS
    REAL*8 :: qs
!LOCALS
    REAL*8, PARAMETER :: zeroc=273.16
    REAL*8, PARAMETER :: tminlqu=zeroc-40.0
    REAL*8, PARAMETER :: b6=6.136820929e-11*100.0
    REAL*8, PARAMETER :: b5=2.034080948e-8*100.0
    REAL*8, PARAMETER :: b4=3.031240396e-6*100.0
    REAL*8, PARAMETER :: b3=2.650648471e-4*100.0
    REAL*8, PARAMETER :: b2=1.428945805e-2*100.0
    REAL*8, PARAMETER :: b1=4.436518521e-1*100.0
    REAL*8, PARAMETER :: b0=6.107799961e+0*100.0
    REAL*8 :: tx, ex, ti, tt
    tx = tl
    IF (tx .LT. tminlqu) THEN
      ti = tminlqu
    ELSE IF (tx .GT. tmaxtbl) THEN
      ti = tmaxtbl
    ELSE
      ti = tx
    END IF
!Starr polynomial fit
    tt = ti - zeroc
    ex = tt*(tt*(tt*(tt*(tt*(tt*b6+b5)+b4)+b3)+b2)+b1) + b0
    tl = tx
    qs = ex
    RETURN
  END SUBROUTINE QSATLQU0
  SUBROUTINE QSATICE0(qs, tl)
    IMPLICIT NONE
!INPUTS
    REAL*8 :: tl
!OUTPUTS
    REAL*8 :: qs
!LOCALS
    REAL*8, PARAMETER :: zeroc=273.16, tminstr=-95.0
    REAL*8, PARAMETER :: tminice=zeroc+tminstr
    REAL*8, PARAMETER :: tstarr1=-75.0, tstarr2=-65.0, tstarr3=-50.0, &
&   tstarr4=-40.0
    REAL*8, PARAMETER :: bi6=1.838826904e-10*100.0
    REAL*8, PARAMETER :: bi5=4.838803174e-8*100.0
    REAL*8, PARAMETER :: bi4=5.824720280e-6*100.0
    REAL*8, PARAMETER :: bi3=4.176223716e-4*100.0
    REAL*8, PARAMETER :: bi2=1.886013408e-2*100.0
    REAL*8, PARAMETER :: bi1=5.034698970e-1*100.0
    REAL*8, PARAMETER :: bi0=6.109177956e+0*100.0
    REAL*8, PARAMETER :: s16=0.516000335e-11*100.0
    REAL*8, PARAMETER :: s15=0.276961083e-8*100.0
    REAL*8, PARAMETER :: s14=0.623439266e-6*100.0
    REAL*8, PARAMETER :: s13=0.754129933e-4*100.0
    REAL*8, PARAMETER :: s12=0.517609116e-2*100.0
    REAL*8, PARAMETER :: s11=0.191372282e+0*100.0
    REAL*8, PARAMETER :: s10=0.298152339e+1*100.0
    REAL*8, PARAMETER :: s26=0.314296723e-10*100.0
    REAL*8, PARAMETER :: s25=0.132243858e-7*100.0
    REAL*8, PARAMETER :: s24=0.236279781e-5*100.0
    REAL*8, PARAMETER :: s23=0.230325039e-3*100.0
    REAL*8, PARAMETER :: s22=0.129690326e-1*100.0
    REAL*8, PARAMETER :: s21=0.401390832e+0*100.0
    REAL*8, PARAMETER :: s20=0.535098336e+1*100.0
    REAL*8 :: tx, ti, tt, w, ex
    tx = tl
    IF (tx .LT. tminice) THEN
      ti = tminice
    ELSE IF (tx .GT. zeroc) THEN
      ti = zeroc
    ELSE
      ti = tx
    END IF
    tt = ti - zeroc
    IF (tt .LT. tstarr1) THEN
      ex = tt*(tt*(tt*(tt*(tt*(tt*s16+s15)+s14)+s13)+s12)+s11) + s10
    ELSE IF (tt .GE. tstarr1 .AND. tt .LT. tstarr2) THEN
      w = (tstarr2-tt)/(tstarr2-tstarr1)
      ex = w*(tt*(tt*(tt*(tt*(tt*(tt*s16+s15)+s14)+s13)+s12)+s11)+s10) +&
&       (1.-w)*(tt*(tt*(tt*(tt*(tt*(tt*s26+s25)+s24)+s23)+s22)+s21)+s20)
    ELSE IF (tt .GE. tstarr2 .AND. tt .LT. tstarr3) THEN
      ex = tt*(tt*(tt*(tt*(tt*(tt*s26+s25)+s24)+s23)+s22)+s21) + s20
    ELSE IF (tt .GE. tstarr3 .AND. tt .LT. tstarr4) THEN
      w = (tstarr4-tt)/(tstarr4-tstarr3)
      ex = w*(tt*(tt*(tt*(tt*(tt*(tt*s26+s25)+s24)+s23)+s22)+s21)+s20) +&
&       (1.-w)*(tt*(tt*(tt*(tt*(tt*(tt*bi6+bi5)+bi4)+bi3)+bi2)+bi1)+bi0)
    ELSE
      ex = tt*(tt*(tt*(tt*(tt*(tt*bi6+bi5)+bi4)+bi3)+bi2)+bi1) + bi0
    END IF
    qs = ex
    RETURN
  END SUBROUTINE QSATICE0
!  Differentiation of dqsatpert in forward (tangent) mode (with options r8):
!   variations   of useful results: dqs qss
!   with respect to varying inputs: temp dqs qss
  SUBROUTINE DQSATPERT_TLM(dqs, dqs_tl, qss, qss_tl, temp, temp_tl, plo&
&   , lm)
    IMPLICIT NONE
!Inputs
    INTEGER :: lm
    REAL*8, DIMENSION(lm) :: temp, plo
    REAL*8, DIMENSION(lm) :: temp_tl
!Outputs
    REAL*8, DIMENSION(lm) :: dqs, qss
    REAL*8, DIMENSION(lm) :: dqs_tl, qss_tl
!Locals
    REAL*8, PARAMETER :: max_mixing_ratio=1.0
    REAL*8 :: esfac
    INTEGER :: k
    REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
    REAL*8 :: tl_tl, tt_tl, ti_tl, dqsat_tl, qsat_tl, qq_tl, dd_tl
    INTEGER :: it
    INTEGER, PARAMETER :: degsubs=100
    REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
    INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
    INTRINSIC NINT
    INTRINSIC INT
    esfac = mapl8_h2omw/mapl8_airmw
    DO k=1,lm
      tl_tl = temp_tl(k)
      tl = temp(k)
      pl = plo(k)
      pp = pl*100.0
      IF (tl .LE. tmintbl) THEN
        ti = tmintbl
        ti_tl = 0.0_8
      ELSE IF (tl .GE. tmaxtbl - .001) THEN
        ti_tl = 0.0_8
        ti = tmaxtbl - .001
      ELSE
        ti_tl = tl_tl
        ti = tl
      END IF
      tt_tl = degsubs*ti_tl
      tt = (ti-tmintbl)*degsubs + 1
      it = INT(tt)
      dqq = estblx(it+1) - estblx(it)
      qq_tl = dqq*tt_tl
      qq = (tt-it)*dqq + estblx(it)
      IF (pp .LE. qq) THEN
        qsat = max_mixing_ratio
        dqsat = 0.0
        qsat_tl = 0.0_8
        dqsat_tl = 0.0_8
      ELSE
        dd_tl = -((-((1.0-esfac)*qq_tl))/(pp-(1.0-esfac)*qq)**2)
        dd = 1.0/(pp-(1.0-esfac)*qq)
        qsat_tl = esfac*(qq_tl*dd+qq*dd_tl)
        qsat = esfac*qq*dd
        dqsat_tl = esfac*degsubs*dqq*pp*(dd_tl*dd+dd*dd_tl)
        dqsat = esfac*degsubs*dqq*pp*(dd*dd)
      END IF
      dqs_tl(k) = dqsat_tl
      dqs(k) = dqsat
      qss_tl(k) = qsat_tl
      qss(k) = qsat
    END DO
  END SUBROUTINE DQSATPERT_TLM
  SUBROUTINE DQSATPERT(dqs, qss, temp, plo, lm)
    IMPLICIT NONE
!Inputs
    INTEGER :: lm
    REAL*8, DIMENSION(lm) :: temp, plo
!Outputs
    REAL*8, DIMENSION(lm) :: dqs, qss
!Locals
    REAL*8, PARAMETER :: max_mixing_ratio=1.0
    REAL*8 :: esfac
    INTEGER :: k
    REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
    INTEGER :: it
    INTEGER, PARAMETER :: degsubs=100
    REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
    INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
    INTRINSIC NINT
    INTRINSIC INT
    esfac = mapl8_h2omw/mapl8_airmw
    DO k=1,lm
      tl = temp(k)
      pl = plo(k)
      pp = pl*100.0
      IF (tl .LE. tmintbl) THEN
        ti = tmintbl
      ELSE IF (tl .GE. tmaxtbl - .001) THEN
        ti = tmaxtbl - .001
      ELSE
        ti = tl
      END IF
      tt = (ti-tmintbl)*degsubs + 1
      it = INT(tt)
      dqq = estblx(it+1) - estblx(it)
      qq = (tt-it)*dqq + estblx(it)
      IF (pp .LE. qq) THEN
        qsat = max_mixing_ratio
        dqsat = 0.0
      ELSE
        dd = 1.0/(pp-(1.0-esfac)*qq)
        qsat = esfac*qq*dd
        dqsat = esfac*degsubs*dqq*pp*(dd*dd)
      END IF
      dqs(k) = dqsat
      qss(k) = qsat
    END DO
  END SUBROUTINE DQSATPERT
!  Differentiation of dqsatscapert in forward (tangent) mode (with options r8):
!   variations   of useful results: dqs qss
!   with respect to varying inputs: temp
  SUBROUTINE DQSATSCAPERT_TLM(dqs, dqs_tl, qss, qss_tl, temp, temp_tl, &
&   plo)
    IMPLICIT NONE
!Inputs
    REAL*8 :: temp, plo
    REAL*8 :: temp_tl
!Outputs
    REAL*8 :: dqs, qss
    REAL*8 :: dqs_tl, qss_tl
!Locals
    REAL*8, PARAMETER :: max_mixing_ratio=1.0
    REAL*8 :: esfac
    INTEGER :: k
    REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
    REAL*8 :: tl_tl, tt_tl, ti_tl, dqsat_tl, qsat_tl, qq_tl, dd_tl
    INTEGER :: it
    INTEGER, PARAMETER :: degsubs=100
    REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
    INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
    INTRINSIC NINT
    INTRINSIC INT
    esfac = mapl8_h2omw/mapl8_airmw
    tl_tl = temp_tl
    tl = temp
    pl = plo
    pp = pl*100.0
    IF (tl .LE. tmintbl) THEN
      ti = tmintbl
      ti_tl = 0.0_8
    ELSE IF (tl .GE. tmaxtbl - .001) THEN
      ti = tmaxtbl - .001
      ti_tl = 0.0_8
    ELSE
      ti_tl = tl_tl
      ti = tl
    END IF
    tt_tl = degsubs*ti_tl
    tt = (ti-tmintbl)*degsubs + 1
    it = INT(tt)
    dqq = estblx(it+1) - estblx(it)
    qq_tl = dqq*tt_tl
    qq = (tt-it)*dqq + estblx(it)
    IF (pp .LE. qq) THEN
      qsat = max_mixing_ratio
      dqsat = 0.0
      qsat_tl = 0.0_8
      dqsat_tl = 0.0_8
    ELSE
      dd_tl = -((-((1.0-esfac)*qq_tl))/(pp-(1.0-esfac)*qq)**2)
      dd = 1.0/(pp-(1.0-esfac)*qq)
      qsat_tl = esfac*(qq_tl*dd+qq*dd_tl)
      qsat = esfac*qq*dd
      dqsat_tl = esfac*degsubs*dqq*pp*(dd_tl*dd+dd*dd_tl)
      dqsat = esfac*degsubs*dqq*pp*(dd*dd)
    END IF
    dqs_tl = dqsat_tl
    dqs = dqsat
    qss_tl = qsat_tl
    qss = qsat
  END SUBROUTINE DQSATSCAPERT_TLM
  SUBROUTINE DQSATSCAPERT(dqs, qss, temp, plo)
    IMPLICIT NONE
!Inputs
    REAL*8 :: temp, plo
!Outputs
    REAL*8 :: dqs, qss
!Locals
    REAL*8, PARAMETER :: max_mixing_ratio=1.0
    REAL*8 :: esfac
    INTEGER :: k
    REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
    INTEGER :: it
    INTEGER, PARAMETER :: degsubs=100
    REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
    INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
    INTRINSIC NINT
    INTRINSIC INT
    esfac = mapl8_h2omw/mapl8_airmw
    tl = temp
    pl = plo
    pp = pl*100.0
    IF (tl .LE. tmintbl) THEN
      ti = tmintbl
    ELSE IF (tl .GE. tmaxtbl - .001) THEN
      ti = tmaxtbl - .001
    ELSE
      ti = tl
    END IF
    tt = (ti-tmintbl)*degsubs + 1
    it = INT(tt)
    dqq = estblx(it+1) - estblx(it)
    qq = (tt-it)*dqq + estblx(it)
    IF (pp .LE. qq) THEN
      qsat = max_mixing_ratio
      dqsat = 0.0
    ELSE
      dd = 1.0/(pp-(1.0-esfac)*qq)
      qsat = esfac*qq*dd
      dqsat = esfac*degsubs*dqq*pp*(dd*dd)
    END IF
    dqs = dqsat
    qss = qsat
  END SUBROUTINE DQSATSCAPERT
END MODULE QSAT_UTIL_D
