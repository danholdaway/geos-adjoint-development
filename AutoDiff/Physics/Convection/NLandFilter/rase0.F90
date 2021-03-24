 SUBROUTINE RASE0(IDIM, IRUN, K0, ICMIN, DT,                  &
                  CONS_CP, CONS_ALHL, CONS_GRAV, CONS_RGAS,   &
                  CONS_H2OMW, CONS_AIRMW, CONS_VIREPS,        &
                  SEEDRAS, SIGE,                              &
                  KCBL, WGT0, WGT1, FRLAND, TS,               &
                  THO, QHO,                                   &
                  CO_AUTO, PLE,                               &
                  CLW, FLXD, CNV_PRC3, CNV_UPDFRC,            &
                  RASPARAMS, ESTBLX                           )

 IMPLICIT NONE

 !INPUTS
 INTEGER,                        INTENT(IN   ) ::  IDIM, IRUN, K0, ICMIN
 REAL(8), DIMENSION (IDIM,K0+1), INTENT(IN   ) ::  PLE
 REAL(8), DIMENSION (     K0+1), INTENT(IN   ) ::  SIGE
 REAL(8),                        INTENT(IN   ) ::  DT, CONS_CP, CONS_ALHL, CONS_GRAV, CONS_RGAS
 REAL(8),                        INTENT(IN   ) ::  CONS_H2OMW,CONS_AIRMW,CONS_VIREPS
 INTEGER, DIMENSION (IDIM)  ,    INTENT(IN   ) ::  SEEDRAS
 INTEGER, DIMENSION (IDIM  ),    INTENT(IN   ) ::  KCBL
 REAL(8), DIMENSION (IDIM  ),    INTENT(IN   ) ::  TS, FRLAND
 REAL(8), DIMENSION (IDIM     ), INTENT(IN   ) ::  CO_AUTO
 REAL(8), DIMENSION (IDIM,K0  ), INTENT(IN   ) ::  WGT0,WGT1
 REAL(8), DIMENSION(:),          INTENT(IN   ) ::  RASPARAMS
 REAL(8), DIMENSION(:),          INTENT(IN   ) ::  ESTBLX

 !OUTPUTS
 REAL(8), DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  CLW , FLXD
 REAL(8), DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  CNV_PRC3
 REAL(8), DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  CNV_UPDFRC

 !PROGNOSTIC
 REAL(8), DIMENSION (IDIM,K0  ), INTENT(INOUT) ::  THO, QHO

 !LOCALS
 INTEGER :: I, IC, L, KK, K

 !Parameters
 REAL(8), PARAMETER :: ONEPKAP = 1.+ 2./7., DAYLEN = 86400.0
 REAL(8), PARAMETER :: RHMAX   = 0.9999
 REAL(8), PARAMETER :: CBL_QPERT = 0.0, CBL_TPERT = 1.0
 REAL(8), PARAMETER :: CBL_TPERT_MXOCN = 2.0, CBL_TPERT_MXLND = 4.0

 !Constants
 REAL(8) :: GRAV, CP, ALHL, CPBG, ALHI, CPI, GRAVI, DDT, LBCP 
   
 !Rasparams
 REAL(8) :: FRICFAC, CLI_CRIT, RASAL1, RASAL2 
 REAL(8) :: FRICLAMBDA
 REAL(8) :: SDQV2, SDQV3, SDQVT1
 REAL(8) :: ACRITFAC,  PBLFRAC, AUTORAMPB
 REAL(8) :: MAXDALLOWED, RHMN, RHMX

 REAL(8) :: MXDIAM
 REAL(8) :: TX2, TX3, AKM, ACR, ALM, TTH, QQH, DQX
 REAL(8) :: WFN, TEM, TRG, TRGEXP, EVP, WLQ, QCC
 REAL(8) :: CLI, TE_A, C00_X, CLI_CRIT_X, TOKI
 REAL(8) :: DT_LYR, RATE, CVW_X, CLOSS, F2, F3, F4
 REAL(8) :: WGHT0, PRCBL, RNDU
 REAL(8) :: LAMBDA_MIN, LAMBDA_MAX
 REAL(8) :: TPERT,QPERT

 REAL(8), DIMENSION(K0) :: POI_SV, QOI_SV
 REAL(8), DIMENSION(K0) :: POI, QOI, DQQ, BET, GAM, CLL
 REAL(8), DIMENSION(K0) :: POI_c, QOI_c
 REAL(8), DIMENSION(K0) :: PRH,  PRI,  GHT, DPT, DPB, PKI
 REAL(8), DIMENSION(K0) :: CLN, RNS, POL,DM
 REAL(8), DIMENSION(K0) :: QST, SSL,  RMF, RNN, RN1, RMFC, RMFP
 REAL(8), DIMENSION(K0) :: GMS, ETA, GMH, EHT,  GM1, HCC, RMFD
 REAL(8), DIMENSION(K0) :: HOL, HST, QOL, ZOL, HCLD, CLL0,CLLX,CLLI
 REAL(8), DIMENSION(K0) :: BKE , CVW, UPDFRC
 REAL(8), DIMENSION(K0) :: RASAL, UPDFRP,BK2,DLL0,DLLX
 REAL(8), DIMENSION(K0) :: WGHT, MASSF
 REAL(8), DIMENSION(K0) :: QSS, DQS, Pf, PK, TEMPf, ZLO
 REAL(8), DIMENSION(K0+1) :: PRJ, PRS, QHT, SHT ,ZET, ZLE, PKE

 CNV_PRC3  =0.
 CNV_UPDFRC=0. 
 FLXD = 0.
 CLW = 0.

 FRICFAC      = RASPARAMS(1)     !  ---  1
 CLI_CRIT     = RASPARAMS(4)     !  ---  4
 RASAL1       = RASPARAMS(5)     !  ---  5
 RASAL2       = RASPARAMS(6)     !  ---  6
 FRICLAMBDA   = RASPARAMS(11)    !  --- 11
 SDQV2        = RASPARAMS(14)    !  --- 14
 SDQV3        = RASPARAMS(15)    !  --- 15
 SDQVT1       = RASPARAMS(16)    !  --- 16
 ACRITFAC     = RASPARAMS(17)    !  --- 17
 PBLFRAC      = RASPARAMS(20)    !  --- 20
 AUTORAMPB    = RASPARAMS(21)    !  --- 21
 RHMN         = RASPARAMS(24)    !  --- 24
 MAXDALLOWED  = RASPARAMS(23)    !  --- 24
 RHMX         = RASPARAMS(25)    !  --- 25

 GRAV  = CONS_GRAV
 ALHL  = CONS_ALHL
 CP    = CONS_CP
 CPI   = 1.0/CP      
 ALHI  = 1.0/ALHL
 GRAVI = 1.0/GRAV
 CPBG  = CP*GRAVI
 DDT   = DAYLEN/DT
 LBCP  = ALHL*CPI
 
 DO I = 1,IRUN

    !CALL FINDBASE
    K = KCBL(I)

    IF ( K > 0 ) THEN 

       !Get saturation specific humidity and gradient wrt to T
       PKE  = (PLE(I,:)/1000.)**(CONS_RGAS/CONS_CP)
       Pf   = 0.5*(PLE(I,1:K0) +  PLE(I,2:K0+1  ) )
       PK   = (Pf/1000.)**(CONS_RGAS/CONS_CP)
       TEMPf = THO(I,:)*PK

       ZLE = 0.0
       ZLO = 0.0
       ZLE(K0+1) = 0.
       do L=K0,1,-1
          ZLE(L) = THO(I,L)   * (1.+CONS_VIREPS*QHO(I,L))
          ZLO(L) = ZLE(L+1) + (CONS_CP/CONS_GRAV)*( PKE(L+1)-PK (L  ) ) * ZLE(L)
          ZLE(L) = ZLO(L)   + (CONS_CP/CONS_GRAV)*( PK (L)  -PKE(L)   ) * ZLE(L)
       end do

       TPERT  = CBL_TPERT * ( TS(I) - ( TEMPf(K0)+ CONS_GRAV*ZLO(K0)/CONS_CP )  ) 
       QPERT  = CBL_QPERT !* ( QSSFC - Q(:,:,K0) ) [CBL_QPERT = 0.0]
       TPERT  = MAX( TPERT , 0.0 )
       QPERT  = MAX( QPERT , 0.0 )
      
       if (FRLAND(I) < 0.1) then
          TPERT = MIN( TPERT , CBL_TPERT_MXOCN ) ! ocean
       else
          TPERT = MIN( TPERT , CBL_TPERT_MXLND ) ! land
       endif

       call DQSAT_RAS(DQS,QSS,TEMPf,Pf,K0,ESTBLX,CONS_H2OMW,CONS_AIRMW)

       do kk=icmin,k+1
          PRJ(kk) = PKE(kk)
       enddo

       poi=0.        ! These initialized here in order not to confuse Valgrind debugger
       qoi=0.        ! Do not believe it actually makes any difference.

       PRS(ICMIN:K0+1) = PLE(I,ICMIN:K0+1)
       POI(ICMIN:K)   = THO(I,ICMIN:K)
       QOI(ICMIN:K)   = QHO(I,ICMIN:K)

       QST(ICMIN:K) = QSS(ICMIN:K)
       DQQ(ICMIN:K) = DQS(ICMIN:K)

       !Mass fraction of each layer below cloud base
       MASSF(:) = WGT0(I,:)

       !RESET PRESSURE at bottom edge of CBL 
       PRCBL = PRS(K)
       do l= K,K0
          PRCBL = PRCBL + MASSF(l)*( PRS(l+1)-PRS(l) )
       end do
       PRS(K+1) = PRCBL
       PRJ(K+1) = (PRS(K+1)/1000.)**(CONS_RGAS/CONS_CP)

       DO L=K,ICMIN,-1
          POL(L)  = 0.5*(PRS(L)+PRS(L+1))
          PRH(L)  = (PRS(L+1)*PRJ(L+1)-PRS(L)*PRJ(L)) / (ONEPKAP*(PRS(L+1)-PRS(L)))
          PKI(L)  = 1.0 / PRH(L)
          DPT(L)  = PRH(L  ) - PRJ(L)
          DPB(L)  = PRJ(L+1) - PRH(L)
          PRI(L)  = .01 / (PRS(L+1)-PRS(L))
       ENDDO

       !RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
       if ( K<=K0) then
          POI(K) = 0.
          QOI(K) = 0.

          !SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
          WGHT = 0.
          DO L=K,K0
             WGHT(L)   = MASSF(L) * ( PLE(I,L+1) - PLE(I,L) ) / ( PRS(K+1)   - PRS(K)  )
          END DO      

          DO L=K,K0
             POI(K) = POI(K) + WGHT(L)*THO(I,L)
             QOI(K) = QOI(K) + WGHT(L)*QHO(I,L)
          ENDDO

          call DQSATs_RAS(DQQ(K), QST(K), POI(K)*PRH(K),POL(K),ESTBLX,CONS_H2OMW,CONS_AIRMW)

       endif

       RNDU = max( seedras(I)/1000000., 1e-6 )
       MXDIAM = maxdallowed*( rndu**(-1./2.) )

       DO L=K,ICMIN,-1
          BET(L)  = DQQ(L)*PKI(L)  !*
          GAM(L)  = PKI(L)/(1.0+LBCP*DQQ(L)) !*
          IF (L<K) THEN
             GHT(L+1) = GAM(L)*DPB(L) + GAM(L+1)*DPT(L+1)
             GM1(L+1) = 0.5*LBCP*(DQQ(L  )/(ALHL*(1.0+LBCP*DQQ(L  ))) + DQQ(L+1)/(ALHL*(1.0+LBCP*DQQ(L+1))) )
          ENDIF
       ENDDO

       RNS  = 0.
       CLL  = 0.
       RMF  = 0.
       RMFD = 0.
       RMFC = 0.
       RMFP = 0.
       CLL0 = 0.
       DLL0 = 0.
       CLLX = 0.
       DLLX = 0.
       CLLI = 0.

       POI_SV = POI
       QOI_SV = QOI

       CVW     = 0.0
       UPDFRC  = 0.0
       UPDFRP  = 0.0

       hol=0.            ! HOL initialized here in order not to confuse Valgrind debugger
       ZET(K+1) = 0
       SHT(K+1) = CP*POI(K)*PRJ(K+1)
       DO L=K,ICMIN,-1
          QOL(L)  = MIN(QST(L)*RHMAX,QOI(L))
          QOL(L)  = MAX( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
          SSL(L)  = CP*PRJ(L+1)*POI(L) + GRAV*ZET(L+1)
          HOL(L)  = SSL(L) + QOL(L)*ALHL
          HST(L)  = SSL(L) + QST(L)*ALHL
          TEM     = POI(L)*(PRJ(L+1)-PRJ(L))*CPBG
          ZET(L)  = ZET(L+1) + TEM
          ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI(L)*CPBG
       ENDDO

       DO IC =  K,ICMIN+1,-1
              
          ALM   = 0.
          TRG   = MIN(1.,(QOI(K)/QST(K)-RHMN)/(RHMX-RHMN))

          F4  = MIN(   1.0,  MAX( 0.0 , (AUTORAMPB-SIGE(IC))/0.2 )  )  
     
          IF (TRG <= 1.0E-5) THEN
             CYCLE !================>>
          ENDIF

          !RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL
          POI_c = POI
          QOI_c = QOI
          POI_c(K) =  POI_c(K) + TPERT
          QOI_c(K) =  QOI_c(K) + QPERT

          ZET(K+1) = 0.
          SHT(K+1) = CP*POI_c(K)*PRJ(K+1)
          DO L=K,IC,-1
             QOL(L)  = MIN(QST(L)*RHMAX,QOI_c(L))
             QOL(L)  = MAX( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
             SSL(L)  = CP*PRJ(L+1)*POI_c(L) + GRAV*ZET(L+1)
             HOL(L)  = SSL(L) + QOL(L)*ALHL
             HST(L)  = SSL(L) + QST(L)*ALHL
             TEM     = POI_c(L)*(PRJ(L+1)-PRJ(L))*CPBG
             ZET(L)  = ZET(L+1) + TEM
             ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI_c(L)*CPBG
          ENDDO

          DO L=IC+1,K
             TEM  = (PRJ(L)-PRH(L-1))/(PRH(L)-PRH(L-1))
             SHT(L)  = SSL(L-1) + TEM*(SSL(L)-SSL(L-1)) 
             QHT(L)  = .5*(QOL(L)+QOL(L-1))
          ENDDO

          !CALCULATE LAMBDA, ETA, AND WORKFUNCTION
          LAMBDA_MIN = .2/MXDIAM
          LAMBDA_MAX = .2/  200. 

          IF (HOL(K) <= HST(IC)) THEN
             CYCLE !================>>
          ENDIF

          !LAMBDA CALCULATION: MS-A18

          TEM  =       (HST(IC)-HOL(IC))*(ZOL(IC)-ZET(IC+1)) 
          DO L=IC+1,K-1
             TEM = TEM + (HST(IC)-HOL(L ))*(ZET(L )-ZET(L +1))
          ENDDO

          IF (TEM <= 0.0) THEN
             CYCLE !================>>
          ENDIF

          ALM = (HOL(K)-HST(IC)) / TEM

          IF (ALM > LAMBDA_MAX) THEN
             CYCLE !================>>
          ENDIF

          TOKI=1.0

          IF (ALM < LAMBDA_MIN) THEN
             TOKI = ( ALM/LAMBDA_MIN )**2
          ENDIF

          !ETA CALCULATION: MS-A2
          DO L=IC+1,K
             ETA(L) = 1.0 + ALM * (ZET(L )-ZET(K))
          ENDDO
          ETA(IC) = 1.0 + ALM * (ZOL(IC)-ZET(K))

          !WORKFUNCTION CALCULATION:  MS-A22
          WFN     = 0.0
          HCC(K)  = HOL(K)
          DO L=K-1,IC+1,-1
             HCC(L) = HCC(L+1) + (ETA(L) - ETA(L+1))*HOL(L)
             TEM    = HCC(L+1)*DPB(L) + HCC(L)*DPT(L)
             EHT(L) = ETA(L+1)*DPB(L) + ETA(L)*DPT(L)
             WFN    = WFN + (TEM - EHT(L)*HST(L))*GAM(L)
          ENDDO
          HCC(IC) = HST(IC)*ETA(IC)
          WFN     = WFN + (HCC(IC+1)-HST(IC)*ETA(IC+1))*GAM(IC)*DPB(IC)

          !VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
          BK2(K)   = 0.0
          BKE(K)   = 0.0
          HCLD(K)  = HOL(K)
          DO L=K-1,IC,-1
             HCLD(L) = ( ETA(L+1)*HCLD(L+1) + (ETA(L) - ETA(L+1))*HOL(L) ) / ETA(L)
             TEM     = (HCLD(L)-HST(L) )*(ZET(L)-ZET(L+1))/ (1.0+LBCP*DQQ(L))
             BKE(L)  = BKE(L+1) + GRAV * TEM / ( CP*PRJ(L+1)*POI(L) )
             BK2(L)  = BK2(L+1) + GRAV * MAX(TEM,0.0) / ( CP*PRJ(L+1)*POI(L) )
             CVW(L) = SQRT(  2.0* MAX( BK2(L) , 0.0 )  )    
          ENDDO

          !ALPHA CALCULATION 
          IF ( ZET(IC) <  2000. ) THEN
             RASAL(IC) = RASAL1
          ENDIF
          IF ( ZET(IC)  >= 2000. ) THEN 
             RASAL(IC) = RASAL1 + (RASAL2-RASAL1)*(ZET(IC) - 2000.)/8000.
          ENDIF
          RASAL(IC) = MIN( RASAL(IC) , 1.0e5 )
    
          RASAL(IC) = DT / RASAL(IC)

          CVW(IC:K) = MAX(  CVW(IC:K) , 1.00 )

          CALL ACRITN(POL(IC), PRS(K), ACR, ACRITFAC)

          IF (WFN <= ACR) THEN
             CYCLE !================>>
          ENDIF

          WLQ     = QOL(K)
          RNN(K)  = 0.
          CLL0(K) = 0.

          DO L=K-1,IC,-1
             TEM   = ETA(L) - ETA(L+1)
             WLQ   = WLQ + TEM * QOL(L)

             IF (L>IC) THEN
                TX2   = 0.5*(QST(L)+QST(L-1))*ETA(L)
                TX3   = 0.5*(HST(L)+HST(L-1))*ETA(L)
                QCC   = TX2 + GM1(L)*(HCC(L)-TX3)
                CLL0(L)  = (WLQ-QCC)
             ELSE
                CLL0(L)   = (WLQ-QST(IC)*ETA(IC))
             ENDIF

             CLL0(L)    = MAX( CLL0(L) , 0.00 )

             CLI  = CLL0(L) / ETA(L)
             TE_A = POI(L)*PRH(L)

             CALL SUNDQ3_ICE( TE_A,  SDQV2, SDQV3, SDQVT1, F2 , F3 )

             C00_X  =  CO_AUTO(I) * F2 * F3  * F4 
             CLI_CRIT_X =  CLI_CRIT / ( F2 * F3 )
       
             RATE = C00_X * ( 1.0 - EXP( -(CLI)**2 / CLI_CRIT_X**2 ) )

             CVW_X     = MAX( CVW(L) , 1.00 )  

             DT_LYR  = ( ZET(L)-ZET(L+1) )/CVW_X 

             CLOSS   = CLL0(L) * RATE * DT_LYR
             CLOSS   = MIN( CLOSS , CLL0(L) )

             CLL0(L) = CLL0(L) - CLOSS
             DLL0(L) = CLOSS

             IF (CLOSS>0.) then
                WLQ = WLQ - CLOSS
                RNN(L) = CLOSS 
             else
                RNN(L) = 0.
             ENDIF

          ENDDO

          WLQ = WLQ - QST(IC)*ETA(IC)

          !CALCULATE GAMMAS AND KERNEL
          GMS(K) =          (SHT(K)-SSL(K))*PRI(K)
          GMH(K) = GMS(K) + (QHT(K)-QOL(K))*PRI(K)*ALHL
          AKM    = GMH(K)*GAM(K-1)*DPB(K-1)

          TX2     = GMH(K)
          DO L=K-1,IC+1,-1
             GMS(L) = ( ETA(L  )*(SHT(L)-SSL(L  )) + ETA(L+1)*(SSL(L)-SHT(L+1)) )     *PRI(L)
             GMH(L) = GMS(L) + ( ETA(L  )*(QHT(L)-QOL(L  )) + ETA(L+1)*(QOL(L)-QHT(L+1)) )*ALHL*PRI(L)
             TX2 = TX2 + (ETA(L) - ETA(L+1)) * GMH(L)
             AKM = AKM - GMS(L)*EHT(L)*PKI(L) + TX2*GHT(L)
          ENDDO

          GMS(IC) = ETA(IC+1)*(SSL(IC)-SHT(IC+1))*PRI(IC)
          AKM     = AKM - GMS(IC)*ETA(IC+1)*DPB(IC)*PKI(IC)

          GMH(IC) =   GMS(IC) + ( ETA(IC+1)*(QOL(IC)-QHT(IC+1))*ALHL + ETA(IC)*(HST(IC)-HOL(IC)))*PRI(IC)

          !CLOUD BASE MASS FLUX
          IF (AKM >= 0.0 .OR. WLQ < 0.0)  THEN
             CYCLE !================>>
          ENDIF

          WFN = - (WFN-ACR)/AKM
          WFN = MIN( ( RASAL(IC)*TRG*TOKI )*WFN  ,   (PRS(K+1)-PRS(K) )*(100.*PBLFRAC))
     
          !CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
          TEM      = WFN*GRAVI
          CLL (IC) = CLL (IC) + WLQ*TEM
          RMF (IC) = RMF (IC) +     TEM
          RMFD(IC) = RMFD(IC) +     TEM * ETA(IC)

          DO L=IC+1,K
             RMFP(L) = TEM * ETA(L)
             RMFC(L) = RMFC(L)  +  RMFP(L)
           
             DLLX(L) = DLLX(L)  +  TEM*DLL0(L)        

             IF ( CVW(L) > 0.0 ) THEN
                UPDFRP(L) = rmfp(L) * (DDT/DAYLEN) * 1000. / ( CVW(L) * PRS(L) )
             ELSE 
                UPDFRP(L) = 0.0       
             ENDIF

             CLLI(L) = CLL0(L)/ETA(L)

             UPDFRC(L) =  UPDFRC(L) +  UPDFRP(L)      
 
          ENDDO

          !THETA AND Q CHANGE DUE TO CLOUD TYPE IC
          DO L=IC,K
             RNS(L) = RNS(L) + RNN(L)*TEM
             GMH(L) = GMH(L) * WFN
             GMS(L) = GMS(L) * WFN
             QOI(L) = QOI(L) + (GMH(L) - GMS(L)) * ALHI
             POI(L) = POI(L) + GMS(L)*PKI(L)*CPI
             QST(L) = QST(L) + GMS(L)*BET(L)*CPI
          ENDDO

       ENDDO !CLOUD LOOP

       IF ( SUM( RMF(ICMIN:K) ) > 0.0 ) THEN

          DO L=ICMIN,K
             TEM    = PRI(L)*GRAV
             CNV_PRC3(I,L) = RNS(L)*TEM     
          ENDDO

          THO(I,ICMIN:K-1) = POI(ICMIN:K-1)
          QHO(I,ICMIN:K-1) = QOI(ICMIN:K-1)
          CNV_UPDFRC(I,ICMIN:K-1)   =  UPDFRC(ICMIN:K-1)

          !De-strap tendencies from RAS
          WGHT   = WGT1(I,:)

          !Scale properly by layer masses
          wght0 = 0.
          DO L=K,K0 
             wght0 = wght0 + WGHT(L)* ( PLE(I,L+1) - PLE(I,L) )
          END DO
         
          wght0 = ( PRS(K+1)   - PRS(K)  )/wght0
          WGHT  = wght0 * WGHT

          DO L=K,K0 
             THO(I,L) =  THO(I,L) + WGHT(L)*(POI(K) - POI_SV(K))
             QHO(I,L) =  QHO(I,L) + WGHT(L)*(QOI(K) - QOI_SV(K))
          END DO

          FLXD(I,ICMIN:K) = RMFD(ICMIN:K) * DDT/DAYLEN
          CLW (I,ICMIN:K) = CLL (ICMIN:K) * DDT/DAYLEN

          FLXD(I,1:ICMIN-1) = 0.
          CLW (I,1:ICMIN-1) = 0.

          IF ( K < K0 ) THEN 
             FLXD(I,K:K0) = 0.
             CLW (I,K:K0) = 0.
          END IF

       ELSE

          FLXD(I,:) = 0.
          CLW (I,:) = 0.

       ENDIF 

    ELSE 
     
       FLXD(I,:) = 0.
       CLW (I,:) = 0.

    ENDIF

 ENDDO

END SUBROUTINE RASE0
