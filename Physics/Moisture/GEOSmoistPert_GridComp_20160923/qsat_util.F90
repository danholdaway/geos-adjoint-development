module qsat_util

!This module contains the subroutines used to generate the lookup table for computing saturation vapour pressure.
!Variables designated as real*8 should be kept at least as so, even if ESTBLX is required in real*4.

use MAPL_ConstantsMod

IMPLICIT NONE
private

integer, parameter :: DEGSUBS    =  100
real(8), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1
real(8), dimension(TABLESIZE) :: ESTBLX

PUBLIC ESINIT, DQSATpert, DQSATpert_TLM, DQSATpert_ADM, DQSATpert_FWD, DQSATpert_BWD

CONTAINS

subroutine ESINIT

 IMPLICIT NONE

 !LOCALS
 real(8), parameter :: ZEROC = 273.16, TMIX = -20.0

 real(8), dimension(TABLESIZE) :: ESTBLE, ESTBLW

 integer :: I
 real(8)    :: T, DELTA_T

 DELTA_T = 1.0/DEGSUBS

 do I=1,TABLESIZE

    T = (I-1)*DELTA_T + TMINTBL

    if(T>ZEROC) then
       call QSATLQU0(ESTBLE(I),T)
    else
       call QSATICE0(ESTBLE(I),T)
    end if

    call QSATLQU0(ESTBLW(I),T)

    T = T-ZEROC
    if(T>=TMIX .and. T<0.0) then
       ESTBLX(I) = ( T/TMIX )*( ESTBLE(I) - ESTBLW(I) ) + ESTBLW(I)
    else
       ESTBLX(I) = ESTBLE(I)
    end if

 end do

 end subroutine ESINIT


subroutine QSATLQU0(QS,TL)
!SUPERSATURATED AS LIQUID

 IMPLICIT NONE

 !INPUTS
 real(8) :: TL!, TMAXTBL

 !OUTPUTS
 real(8) :: QS

 !LOCALS
 real(8), parameter :: ZEROC   = 273.16
 real(8), parameter :: TMINLQU = ZEROC - 40.0

 real*8,  parameter :: B6 = 6.136820929E-11*100.0
 real*8,  parameter :: B5 = 2.034080948E-8 *100.0
 real*8,  parameter :: B4 = 3.031240396E-6 *100.0
 real*8,  parameter :: B3 = 2.650648471E-4 *100.0
 real*8,  parameter :: B2 = 1.428945805E-2 *100.0
 real*8,  parameter :: B1 = 4.436518521E-1 *100.0
 real*8,  parameter :: B0 = 6.107799961E+0 *100.0

 real(8) :: TX, EX, TI, TT

 TX = TL

 if    (TX<TMINLQU) then
    TI = TMINLQU
 elseif(TX>TMAXTBL) then
    TI = TMAXTBL
 else
    TI = TX
 end if

 TT = TI-ZEROC  !Starr polynomial fit
 EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)

 TL = TX
 QS = EX

 return

end subroutine QSATLQU0


subroutine QSATICE0(QS,TL)
!SUPERSATURATED AS ICE

 IMPLICIT NONE   
      
 !INPUTS
 real(8) :: TL

 !OUTPUTS
 real(8) :: QS

 !LOCALS
 real(8), parameter :: ZEROC = 273.16, TMINSTR = -95.0
 real(8), parameter :: TMINICE = ZEROC + TMINSTR

 real(8), parameter :: TSTARR1 = -75.0, TSTARR2 = -65.0, TSTARR3 = -50.0,  TSTARR4 = -40.0

 real*8,  parameter :: BI6= 1.838826904E-10*100.0
 real*8,  parameter :: BI5= 4.838803174E-8 *100.0
 real*8,  parameter :: BI4= 5.824720280E-6 *100.0
 real*8,  parameter :: BI3= 4.176223716E-4 *100.0
 real*8,  parameter :: BI2= 1.886013408E-2 *100.0
 real*8,  parameter :: BI1= 5.034698970E-1 *100.0
 real*8,  parameter :: BI0= 6.109177956E+0 *100.0
 real*8,  parameter :: S16= 0.516000335E-11*100.0
 real*8,  parameter :: S15= 0.276961083E-8 *100.0
 real*8,  parameter :: S14= 0.623439266E-6 *100.0
 real*8,  parameter :: S13= 0.754129933E-4 *100.0
 real*8,  parameter :: S12= 0.517609116E-2 *100.0
 real*8,  parameter :: S11= 0.191372282E+0 *100.0
 real*8,  parameter :: S10= 0.298152339E+1 *100.0
 real*8,  parameter :: S26= 0.314296723E-10*100.0
 real*8,  parameter :: S25= 0.132243858E-7 *100.0
 real*8,  parameter :: S24= 0.236279781E-5 *100.0
 real*8,  parameter :: S23= 0.230325039E-3 *100.0
 real*8,  parameter :: S22= 0.129690326E-1 *100.0
 real*8,  parameter :: S21= 0.401390832E+0 *100.0
 real*8,  parameter :: S20= 0.535098336E+1 *100.0

 real(8) :: TX, TI, TT, W, EX

 TX = TL

 if (TX<TMINICE) then
    TI = TMINICE
 elseif(TX>ZEROC  ) then
    TI = ZEROC
 else
    TI = TX
 end if

 TT = TI - ZEROC
 if (TT < TSTARR1 ) then
     EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
 elseif(TT >= TSTARR1 .and. TT < TSTARR2) then
     W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
     EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
              + (1.-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
 elseif(TT >= TSTARR2 .and. TT < TSTARR3) then
     EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
 elseif(TT >= TSTARR3 .and. TT < TSTARR4) then
     W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
     EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
              + (1.-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
 else
     EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
 endif

 QS = EX

 return
 
end subroutine QSATICE0

subroutine DQSATpert(DQS,QSS,TEMP,PLO,lm)
!COMPUTES SATURATION VAPOUR PRESSURE QSS AND GRADIENT w.r.t TEMPERATURE DQS.
!INPUTS ARE TEMPERATURE AND PLO (PRESSURE AT T-LEVELS)
!VALES ARE COMPUTED FROM LOOK-UP TALBE (PIECEWISE LINEAR)

 IMPLICIT NONE

 !Inputs
 INTEGER :: lm
 REAL(8), dimension(lm) :: TEMP, PLO

 !Outputs
 REAL(8), dimension(lm) :: DQS, QSS

 !Locals
 REAL(8), parameter :: MAX_MIXING_RATIO = 1.0
 REAL(8) :: ESFAC

 INTEGER :: k

 REAL(8) :: TL, TT, TI, DQSAT, QSAT, DQQ, QQ, PL, PP, DD
 INTEGER :: IT

 ESFAC = MAPL_H2OMW/MAPL_AIRMW

 do K=1,LM

    TL = TEMP(K)
    PL = PLO(K)

    PP = PL*100.0

    if (TL<=TMINTBL) then
       TI = TMINTBL
    elseif(TL>=TMAXTBL-.001) then
       TI = TMAXTBL-.001
    else
       TI = TL
    end if

    TT = (TI - TMINTBL)*DEGSUBS+1
    IT = int(TT)

    DQQ =  ESTBLX(IT+1) - ESTBLX(IT)
    QQ  =  (TT-IT)*DQQ + ESTBLX(IT)

    if (PP <= QQ) then
       QSAT = MAX_MIXING_RATIO
       DQSAT = 0.0
    else
       DD = 1.0/(PP - (1.0-ESFAC)*QQ)
       QSAT = ESFAC*QQ*DD
       DQSAT = (ESFAC*DEGSUBS)*DQQ*PP*(DD*DD)
    end if

    DQS(K) = DQSAT
    QSS(K) = QSAT

 end do

end subroutine DQSATpert

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
    INTRINSIC NINT
    INTRINSIC INT
    esfac = mapl_h2omw/mapl_airmw
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

  SUBROUTINE DQSATPERT_ADM(dqs, dqs_ad, qss, qss_ad, temp, temp_ad, plo&
&   , lm)
    IMPLICIT NONE
!Inputs
    INTEGER :: lm
    REAL*8, DIMENSION(lm) :: temp, plo
    REAL*8, DIMENSION(lm) :: temp_ad
!Outputs
    REAL*8, DIMENSION(lm) :: dqs, qss
    REAL*8, DIMENSION(lm) :: dqs_ad, qss_ad
!Locals
    REAL*8, PARAMETER :: max_mixing_ratio=1.0
    REAL*8 :: esfac
    INTEGER :: k
    REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
    REAL*8 :: tl_ad, tt_ad, ti_ad, dqsat_ad, qsat_ad, qq_ad, dd_ad
    INTEGER :: it
    INTRINSIC NINT
    INTRINSIC INT
    INTEGER :: branch
    REAL*8 :: temp0
    esfac = mapl_h2omw/mapl_airmw
    DO k=1,lm
      tl = temp(k)
      pl = plo(k)
      pp = pl*100.0
      IF (tl .LE. tmintbl) THEN
        ti = tmintbl
        CALL PUSHCONTROL2B(0)
      ELSE IF (tl .GE. tmaxtbl - .001) THEN
        ti = tmaxtbl - .001
        CALL PUSHCONTROL2B(1)
      ELSE
        ti = tl
        CALL PUSHCONTROL2B(2)
      END IF
      tt = (ti-tmintbl)*degsubs + 1
      it = INT(tt)
      CALL PUSHREAL8(dqq)
      dqq = estblx(it+1) - estblx(it)
      CALL PUSHREAL8(qq)
      qq = (tt-it)*dqq + estblx(it)
      IF (pp .LE. qq) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
    END DO
    DO k=lm,1,-1
      qsat_ad = qss_ad(k)
      qss_ad(k) = 0.0_8
      dqsat_ad = dqs_ad(k)
      dqs_ad(k) = 0.0_8
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        qq_ad = 0.0_8
      ELSE
        pl = plo(k)
        pp = pl*100.0
        dd = 1.0/(pp-(1.0-esfac)*qq)
        dd_ad = esfac*qq*qsat_ad + esfac*degsubs*dqq*pp*2*dd*dqsat_ad
        temp0 = pp - (-esfac+1.0)*qq
        qq_ad = (1.0-esfac)*dd_ad/temp0**2 + esfac*dd*qsat_ad
      END IF
      CALL POPREAL8(qq)
      tt_ad = dqq*qq_ad
      CALL POPREAL8(dqq)
      ti_ad = degsubs*tt_ad
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        tl_ad = 0.0_8
      ELSE IF (branch .EQ. 1) THEN
        tl_ad = 0.0_8
      ELSE
        tl_ad = ti_ad
      END IF
      temp_ad(k) = temp_ad(k) + tl_ad
    END DO
  END SUBROUTINE DQSATPERT_ADM

  SUBROUTINE DQSATPERT_FWD(dqs, qss, temp, plo, lm)
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
    esfac = mapl_h2omw/mapl_airmw
    DO k=1,lm
      tl = temp(k)
      pl = plo(k)
      pp = pl*100.0
      IF (tl .LE. tmintbl) THEN
        ti = tmintbl
        CALL PUSHCONTROL2B(0)
      ELSE IF (tl .GE. tmaxtbl - .001) THEN
        ti = tmaxtbl - .001
        CALL PUSHCONTROL2B(1)
      ELSE
        ti = tl
        CALL PUSHCONTROL2B(2)
      END IF
      tt = (ti-tmintbl)*degsubs + 1
      it = INT(tt)
      CALL PUSHREAL8(dqq)
      dqq = estblx(it+1) - estblx(it)
      CALL PUSHREAL8(qq)
      qq = (tt-it)*dqq + estblx(it)
      IF (pp .LE. qq) THEN
        qsat = max_mixing_ratio
        dqsat = 0.0
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHREAL8(dd)
        dd = 1.0/(pp-(1.0-esfac)*qq)
        qsat = esfac*qq*dd
        dqsat = esfac*degsubs*dqq*pp*(dd*dd)
        CALL PUSHCONTROL1B(1)
      END IF
      CALL PUSHREAL8(dqs(k))
      dqs(k) = dqsat
      CALL PUSHREAL8(qss(k))
      qss(k) = qsat
    END DO
    CALL PUSHREAL8(dd)
    CALL PUSHREAL8(esfac)
    CALL PUSHREAL8(qq)
    CALL PUSHREAL8(dqq)
  END SUBROUTINE DQSATPERT_FWD

  SUBROUTINE DQSATPERT_BWD(dqs, dqs_ad, qss, qss_ad, temp, temp_ad, plo&
&   , lm)
    IMPLICIT NONE
    INTEGER :: lm
    REAL*8, DIMENSION(lm) :: temp, plo
    REAL*8, DIMENSION(lm) :: temp_ad
    REAL*8, DIMENSION(lm) :: dqs, qss
    REAL*8, DIMENSION(lm) :: dqs_ad, qss_ad
    REAL*8, PARAMETER :: max_mixing_ratio=1.0
    REAL*8 :: esfac
    INTEGER :: k
    REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
    REAL*8 :: tl_ad, tt_ad, ti_ad, dqsat_ad, qsat_ad, qq_ad, dd_ad
    INTEGER :: it
    INTEGER, PARAMETER :: degsubs=100
    REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
    INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
    INTRINSIC NINT
    INTRINSIC INT
    INTEGER :: branch
    REAL*8 :: temp0
    CALL POPREAL8(dqq)
    CALL POPREAL8(qq)
    CALL POPREAL8(esfac)
    CALL POPREAL8(dd)
    DO k=lm,1,-1
      CALL POPREAL8(qss(k))
      qsat_ad = qss_ad(k)
      qss_ad(k) = 0.0_8
      CALL POPREAL8(dqs(k))
      dqsat_ad = dqs_ad(k)
      dqs_ad(k) = 0.0_8
      pl = plo(k)
      pp = pl*100.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        qq_ad = 0.0_8
      ELSE
        temp0 = pp - (-esfac+1.0)*qq
        dd_ad = esfac*qq*qsat_ad + esfac*degsubs*dqq*pp*2*dd*dqsat_ad
        qq_ad = (1.0-esfac)*dd_ad/temp0**2 + esfac*dd*qsat_ad
        CALL POPREAL8(dd)
      END IF
      CALL POPREAL8(qq)
      tt_ad = dqq*qq_ad
      CALL POPREAL8(dqq)
      ti_ad = degsubs*tt_ad
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        tl_ad = 0.0_8
      ELSE IF (branch .EQ. 1) THEN
        tl_ad = 0.0_8
      ELSE
        tl_ad = ti_ad
      END IF
      temp_ad(k) = temp_ad(k) + tl_ad
    END DO
  END SUBROUTINE DQSATPERT_BWD

end module qsat_util
