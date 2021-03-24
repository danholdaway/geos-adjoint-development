subroutine precipandevap(K, LM, DT, FRLAND, RHCR3, QPl, QPi, QCl, QCi, TE, QV, mass, imass, &
                            PL, dZE, QDDF3, AA, BB, AREA, PFl_above_in, PFl_above_out, PFi_above_in, PFi_above_out, &
                            EVAP_DD_above_in, EVAP_DD_above_out, SUBL_DD_above_in, SUBL_DD_above_out, &
                            FRZ_DIAG, ENVFC, DDRFC,  &
                            CONS_ALHF, CONS_ALHS, CONS_ALHL, CONS_CP, CONS_TICE, REVAP_OFF_P, &
                            C_ACC, C_EV_R, C_EV_S, RHO_W, ESTBLX )

!                            RAIN, SNOW, REVAP_DIAG, RSUBL_DIAG, ACRLL_DIAG, ACRIL_DIAG, &
!                            PFL_DIAG, PFI_DIAG, VFALLSN, VFALLRN )

IMPLICIT NONE

!Inputs
integer, intent(in) :: K,LM
real(8), intent(in) :: DT, mass, imass, PL, AA, BB, RHCR3, dZE, QDDF3, AREA, FRLAND, ENVFC, DDRFC
real(8), intent(in) :: CONS_ALHF, CONS_ALHS, CONS_ALHL, CONS_CP, CONS_TICE, REVAP_OFF_P
real(8), intent(in) :: C_ACC, C_EV_R, C_EV_S, RHO_W
real(8), intent(in) :: ESTBLX(50)

!Prognostics
real(8), intent(inout) :: QV, QPl, QPi, QCl, QCi, TE, FRZ_DIAG
real(8), intent(inout) :: PFl_above_in, Pfl_above_out, PFi_above_in, PFi_above_out
real(8), intent(inout) :: EVAP_DD_above_in, EVAP_DD_above_out, SUBL_DD_above_in, SUBL_DD_above_out

!Outputs - not used for now
!real(8), intent(  out) :: RAIN, SNOW, REVAP_DIAG, RSUBL_DIAG, ACRLL_DIAG, ACRIL_DIAG, PFL_DIAG, PFI_DIAG, VFALLSN, VFALLRN

!Locals
integer :: NS, NSMX, itr,L

real(8) :: PFi, PFl, QS, dQS, ENVFRAC, TKo, QKo, QSTKo, DQSTKo, RH_BOX, T_ED, QPlKo, QPiKo
real(8) :: Ifactor, RAINRAT0, SNOWRAT0, FALLRN, FALLSN, VEsn, VErn, NRAIN, NSNOW, Efactor
real(8) :: TinLAYERrn, DIAMrn, DROPRAD, TinLAYERsn, DIAMsn, FLAKRAD
real(8) :: EVAP, SUBL, ACCR, MLTFRZ, EVAPx, SUBLx, EVAP_DD,SUBL_DD,DDFRACT, LANDSEAF
real(8) :: TAU_FRZ, TAU_MLT

real(8), parameter :: TRMV_L = 1.0   !m/s
logical, parameter :: taneff = .false.

!Fraction of precip falling through "environment" vs through cloud
real(8), parameter :: B_SUB = 1.00

 !Initialize diagnostics
 !rain = 0.0
 !snow = 0.0
 !revap_diag = 0.0
 !rsubl_diag = 0.0
 !acrll_diag = 0.0
 !acril_diag = 0.0
 !pfl_diag = 0.0
 !pfi_diag = 0.0
 !vfallsn = 0.0
 !vfallrn = 0.0

 if (taneff) then

    !envfrac = 1.00

    !if (pl .le. 600.) then
    !   envfrac = 0.25
    !else
    !   envfrac = 0.25 + (1.-0.25)/(19.) * ((atan( (2.*(pl-600.)/(900.-600.)-1.) *       &
    !                                    tan(20.*CONS_PI/21.-0.5*CONS_PI) ) + 0.5*CONS_PI) * 21./CONS_PI - 1.)
    !end if
    !
    !envfrac = min(envfrac,1.)

 else

    ENVFRAC = ENVFC

 endif


 IF ( AREA > 0. ) THEN
    Ifactor = 1./ ( AREA )
 ELSE
    Ifactor = 1.00
 END if

 Ifactor = MAX( Ifactor, 1.) !

 !Start at top of precip column:
 !
 !   a) Accrete                   
 !   b) Evaporate/Sublimate  
 !   c) Rain/Snow-out to next level down 
 !   d) return to (a)

 !Update saturated humidity
 call DQSAT_sub_sca(dQS,QS,TE,PL,ESTBLX)

 DDFRACT = DDRFC 

 IF (K == 1) THEN

    PFl=QPl*MASS
    PFi=QPi*MASS

    EVAP_DD = 0.
    SUBL_DD = 0.

    !VFALLRN = 0.0
    !VFALLSN = 0.0

 ELSE 

    QPl = QPl + PFl_above_in * iMASS
    PFl = 0.00

    QPi = QPi + PFi_above_in * iMASS
    PFi = 0.00

    ACCR = B_SUB * C_ACC * ( QPl*MASS ) *QCl   

    ACCR = MIN(  ACCR , QCl  )

    QPl     = QPl + ACCR
    QCl     = QCl - ACCR

    !ACRLL_DIAG = ACCR / DT

    !Accretion of liquid condensate by falling ice/snow
    ACCR = B_SUB * C_ACC * ( QPi*MASS ) *QCl   

    ACCR = MIN(  ACCR , QCl  )

    QPi = QPi + ACCR
    QCl = QCl - ACCR

    !! Liquid freezes when accreted by snow
    TE  = TE + CONS_ALHF*ACCR/CONS_CP

    !ACRIL_DIAG = ACCR / DT

    RAINRAT0 = Ifactor*QPl*MASS/DT
    SNOWRAT0 = Ifactor*QPi*MASS/DT

    call MARSHPALM(RAINRAT0,PL,DIAMrn,NRAIN,FALLrn,VErn)
    call MARSHPALM(SNOWRAT0,PL,DIAMsn,NSNOW,FALLsn,VEsn)

    !VFALLRN = FALLrn
    !VFALLSN = FALLsn

    TinLAYERrn = dZE / ( FALLrn+0.01 )
    TinLAYERsn = dZE / ( FALLsn+0.01 )

    !Melting of Frozen precipitation      
    TAU_FRZ = 5000.  ! time scale for freezing (s). 

    MLTFRZ = 0.0
    IF ( (TE > CONS_TICE ) .and.(TE <= CONS_TICE+5. ) ) THEN
       MLTFRZ = TinLAYERsn * QPi *( TE - CONS_TICE ) / TAU_FRZ 
       MLTFRZ = MIN( QPi , MLTFRZ )
       TE  = TE  - CONS_ALHF*MLTFRZ/CONS_CP
       QPl = QPl + MLTFRZ
       QPi = QPi - MLTFRZ
    END IF
    FRZ_DIAG = FRZ_DIAG - MLTFRZ / DT

    MLTFRZ = 0.0
    IF ( TE > CONS_TICE+5.  ) THEN  ! Go Ahead and melt any snow/hail left above 5 C 
       MLTFRZ= QPi 
       TE  = TE  - CONS_ALHF*MLTFRZ/CONS_CP
       QPl = QPl + MLTFRZ
       QPi = QPi - MLTFRZ
    END IF
    FRZ_DIAG = FRZ_DIAG - MLTFRZ / DT

    MLTFRZ = 0.0
    if ( K >= LM-1 ) THEN
       IF ( TE > CONS_TICE+0.  ) THEN ! Go Ahead and melt any snow/hail left above 0 C in lowest layers 
          MLTFRZ= QPi 
          TE  = TE  - CONS_ALHF*MLTFRZ/CONS_CP
          QPl = QPl + MLTFRZ
          QPi = QPi - MLTFRZ
       END IF
    endif
    FRZ_DIAG = FRZ_DIAG - MLTFRZ / DT

    !Freezing of liquid precipitation      
    MLTFRZ = 0.0
    IF ( TE <= CONS_TICE ) THEN
       TE  = TE + CONS_ALHF*QPl/CONS_CP
       QPi = QPl + QPi
       MLTFRZ = QPl
       QPl = 0.
    END IF
    FRZ_DIAG = FRZ_DIAG + MLTFRZ / DT

    !In the exp below, evaporation time scale is determined "microphysically" from temp, 
    !press, and drop size. In this context C_EV becomes a dimensionless fudge-fraction. 
    !Also remember that these microphysics are still only for liquid.

    QKo   = QV
    TKo   = TE
    QPlKo = QPl
    QPiKo = QPi

    do itr = 1,3

       DQSTKo = dQS
       QSTKo  = QS  + DQSTKo * ( TKo - TE )
       QSTKo  = MAX( QSTKo , 1.0e-7 )

       RH_BOX = QKo/QSTKo

       QKo   = QV
       TKo   = TE

       IF ( RH_BOX < RHCR3 ) THEN
          Efactor =  RHO_W * ( AA + BB )    / (RHCR3 - RH_BOX )
       else
          Efactor = 9.99e9
       end if

       if ( FRLAND < 0.1 ) then
          LANDSEAF = 0.5  ! Over Ocean
       else
          LANDSEAF = 0.5  ! Over Land
       end if

       LANDSEAF = 1.00

       !Rain falling
       if ( (RH_BOX < RHCR3) .AND. (DIAMrn > 0.00) .AND. (PL > 100.) .AND. (PL < REVAP_OFF_P) ) then
          DROPRAD=0.5*DIAMrn
          T_ED =  Efactor * DROPRAD**2 
          T_ED =  T_ED * ( 1.0 + DQSTKo*CONS_ALHL/CONS_CP )
          EVAP =  QPl*(1.0 - EXP( -C_EV_R * VErn * LANDSEAF * ENVFRAC * TinLAYERrn / T_ED ) )
       ELSE
          EVAP = 0.0
       END if

       !Snow falling
       if ( (RH_BOX < RHCR3) .AND. (DIAMsn > 0.00) .AND. (PL > 100.) .AND. (PL < REVAP_OFF_P) ) then
          FLAKRAD=0.5*DIAMsn
          T_ED =  Efactor * FLAKRAD**2   
          T_ED =  T_ED * ( 1.0 + DQSTKo*CONS_ALHS/CONS_CP )
          SUBL =  QPi*(1.0 - EXP( -C_EV_S * VEsn * LANDSEAF * ENVFRAC * TinLAYERsn / T_ED ) )
       ELSE
          SUBL = 0.0
       END IF

       if (itr == 1) then 
          EVAPx  = EVAP
          SUBLx  = SUBL
       else
          EVAP   = (EVAP+EVAPx) /2.0
          SUBL   = (SUBL+SUBLx) /2.0
       endif

       QKo = QV + EVAP + SUBL
       TKo = TE - EVAP * CONS_ALHL / CONS_CP - SUBL * CONS_ALHS / CONS_CP

    enddo
      
    QPi  = QPi - SUBL
    QPl  = QPl - EVAP

    !Put some re-evap/re-subl precip in to a \quote{downdraft} to be applied later
    EVAP_DD = EVAP_DD_above_in + DDFRACT*EVAP*MASS 
    EVAP    = EVAP          - DDFRACT*EVAP
    SUBL_DD = SUBL_DD_above_in + DDFRACT*SUBL*MASS 
    SUBL    = SUBL          - DDFRACT*SUBL

    QV   = QV  + EVAP + SUBL
    TE   = TE  - EVAP * CONS_ALHL / CONS_CP - SUBL * CONS_ALHS / CONS_CP

    !REVAP_DIAG = EVAP / DT
    !RSUBL_DIAG = SUBL / DT

    PFl  = QPl*MASS
    PFi  = QPi*MASS

    !PFL_DIAG =  PFl/DT
    !PFI_DIAG =  PFi/DT
 end if

 EVAP = QDDF3*EVAP_DD/MASS
 SUBL = QDDF3*SUBL_DD/MASS
 QV   = QV  + EVAP + SUBL
 TE   = TE  - EVAP * CONS_ALHL / CONS_CP - SUBL * CONS_ALHS / CONS_CP
 !REVAP_DIAG = REVAP_DIAG + EVAP / DT
 !RSUBL_DIAG = RSUBL_DIAG + SUBL / DT

 IF (K == LM) THEN
    !RAIN  = PFl/DT
    !SNOW  = PFi/DT
 END IF

 QPi = 0.
 QPl = 0.

 PFl_above_out = PFl
 PFi_above_out = Pfi

 EVAP_DD_above_out = EVAP_DD
 SUBL_DD_above_out = SUBL_DD

end subroutine precipandevap

subroutine DQSAT_sub_sca(DQSi,QSSi,TEMP,PLO,ESTBLX)
!COMPUTES SATURATION VAPOUR PRESSURE QSSi AND GRADIENT w.r.t TEMPERATURE DQSi.
!INPUTS ARE TEMPERATURE AND PLO (PRESSURE AT T-LEVELS)
!VALES ARE COMPUTED FROM LOOK-UP TALBE (PIECEWISE LINEAR)

 IMPLICIT NONE

 !Inputs
 real(8) :: TEMP, PLO
 real(8) :: ESTBLX(:)

 !Outputs
 real(8) :: DQSi, QSSi

 !Locals
 real(8), parameter :: MAX_MIXING_RATIO = 1.0
 real(8), parameter :: CONS_H2OMW = 18.01, CONS_AIRMW = 28.97
 real(8),    parameter :: ESFAC = CONS_H2OMW/CONS_AIRMW

 real(8) :: TL, TT, TI, DQSAT, QSAT, DQQ, QQ, PL, PP, DD
 integer :: IT

 integer, parameter :: DEGSUBS    =  100
 real(8), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
 integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1

 TL = TEMP
 PL = PLO

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

 DQSi = DQSAT
 QSSi = QSAT

end subroutine DQSAT_sub_sca

subroutine marshpalm(RAIN,PR,DIAM3,NTOTAL,W,VE)

Implicit None

!Inputs
real(8), intent(in ) :: RAIN, PR     ! in kg m^-2 s^-1, mbar

!Outputs
real(8), intent(out) :: DIAM3, NTOTAL, W, VE

!Locals
integer :: IQD
real(8), parameter  :: N0 = 0.08  !cm^-3
real(8) :: RAIN_DAY,SLOPR,DIAM1

real(8) :: RX(8) , D3X(8)

 RAIN_DAY = 0.0
 SLOPR = 0.0
 DIAM1 = 0.0

 !Marshall-Palmer sizes at different rain-rates: avg(D^3)
 !RX = (/ 0.   , 5.   , 20.  , 80.  , 320. , 1280., 5120., 20480. /)  ! rain per in mm/day
 RX(1) = 0.
 RX(2) = 5.
 RX(3) = 20.
 RX(4) = 80.
 RX(5) = 320.
 RX(6) = 1280.
 RX(7) = 5120.
 RX(8) = 20480.

 !D3X= (/ 0.019, 0.032, 0.043, 0.057, 0.076, 0.102, 0.137, 0.183  /)
 D3X(1) = 0.019
 D3X(2) = 0.032
 D3X(3) = 0.043
 D3X(4) = 0.057
 D3X(5) = 0.076
 D3X(6) = 0.102
 D3X(7) = 0.137
 D3X(8) = 0.183

 RAIN_DAY = RAIN * 3600. *24.

 IF ( RAIN_DAY <= 0.00 ) THEN
    DIAM1 = 0.00
    DIAM3 = 0.00
    NTOTAL= 0.00
    W     = 0.00
 END IF

 DO IQD = 1,7
    IF ( (RAIN_DAY <= RX(IQD+1)) .AND. (RAIN_DAY > RX(IQD) ) ) THEN
       SLOPR =( D3X(IQD+1)-D3X(IQD) ) / ( RX(IQD+1)-RX(IQD) )
       DIAM3 = D3X(IQD) + (RAIN_DAY-RX(IQD))*SLOPR
    END IF
 END DO

 IF ( RAIN_DAY >= RX(8) ) THEN
    DIAM3=D3X(8)
 END IF

 NTOTAL = 0.019*DIAM3

 DIAM3  = 0.664 * DIAM3  

 W      = (2483.8 * DIAM3 + 80.)*SQRT(1000./PR)

 VE     = MAX( 0.99*W/100. , 1.000 )

 DIAM1  = 3.0*DIAM3

 DIAM1  = DIAM1/100.
 DIAM3  = DIAM3/100.
 W      = W/100.
 NTOTAL = NTOTAL*1.0e6

end subroutine marshpalm

