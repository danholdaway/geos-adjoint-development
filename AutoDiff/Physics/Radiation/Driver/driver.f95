subroutine radiation_driver_dummy(PTT,QVT,O3T,CFLST,CFCNT,QIT,QLT,DU001T, DU002T, DU003T, DU004T, DU005T, &
                                  IM,JM,LM,RUNalarm,NA, NBCHOU_IR, NBCHOU_SO,DT,PLE,UT,PS,TS, &
                              EMIS, DELT, COSZ, SLR, RGBUV, RGFUV, RGBIR, RGFIR , SC, &
                              TAUA_IR_C, SSAA_IR_C, ASYA_IR_C, TAUA_SO_C, SSAA_SO_C, ASYA_SO_C, &
                              ILSF, ICNF, LLSF, LCNF, &
                              levs925, LCLDMH, LCLDLM, &
                              mapl_p00, mapl_rgas, mapl_cp, mapl_kappa, mapl_grav, &
                              aib_ir, awb_ir, aiw_ir,   &
                              aww_ir, aig_ir, awg_ir,   &
                              xkw, xke, mw, aw, bw, pm, &
                              fkw, gkw, cb, dcb,        &
                              w11, w12, w13, p11, p12,  &
                              p13, dwe, dpe,            &
                              c1,  c2,  c3,             &
                              oo1, oo2, oo3,            &
                              h11, h12, h13,            &
                              h21, h22, h23,            &
                              h81, h82, h83,             &
                              wk_uv, zk_uv, ry_uv,         &
                              xk_ir, ry_ir,                &
                              cah, coa,                    &
                              aig_uv, awg_uv, arg_uv,      &
                              aib_uv, awb_uv, arb_uv,      &
                              aib_nir, awb_nir, arb_nir,   &
                              aia_nir, awa_nir, ara_nir,   &
                              aig_nir, awg_nir, arg_nir,   &
                              caib, caif  ,                 &
                              HK_UV_TEMP, HK_IR_TEMP          )

IMPLICIT NONE

real(8), intent(inout), dimension(IM,JM,LM) :: PTT, QVT, O3T, CFLST, CFCNT, QIT, QLT
real(8), intent(inout), dimension(IM,JM,LM) :: DU001T, DU002T, DU003T, DU004T, DU005T

integer, intent(in) :: IM,JM,LM,RUNalarm
integer, intent(in) :: NA, NBCHOU_IR, NBCHOU_SO
real(8), intent(in) :: SC, DT
real(8), intent(in), dimension(IM,JM,0:LM) :: PLE
real(8), intent(in), dimension(IM,JM,LM) :: UT
real(8), intent(in), dimension(IM,JM) :: PS,TS
real(8), intent(in), dimension(IM,JM) :: EMIS, DELT, COSZ, SLR, RGBUV, RGFUV, RGBIR, RGFIR

real(8), intent(in), dimension(IM,JM,LM) :: ILSF, ICNF, LLSF, LCNF

real(8), intent(in) ::  HK_IR_TEMP(3,10), HK_UV_TEMP(5)

real(8), intent(in), dimension(IM,JM,LM,NBCHOU_IR,NA) :: TAUA_IR_C, SSAA_IR_C, ASYA_IR_C
real(8), intent(in), dimension(IM,JM,LM,NBCHOU_SO,NA) :: TAUA_SO_C, SSAA_SO_C, ASYA_SO_C

integer, intent(in) :: levs925, LCLDMH, LCLDLM

real(8), intent(in) :: mapl_p00, mapl_rgas, mapl_cp, mapl_kappa, mapl_grav

 real(8), intent(IN)    :: aib_ir(3,10), awb_ir(4,10), aiw_ir(4,10)
 real(8), intent(IN)    :: aww_ir(4,10), aig_ir(4,10), awg_ir(4,10)

 integer, intent(IN) :: mw(9)
 real(8), intent(IN)    :: xkw(9), xke(9), aw(9), bw(9), pm(9)
 real(8), intent(IN)    :: fkw(6,9), gkw(6,3), cb(6,10), dcb(5,10)
 real(8), intent(IN)    :: w11, w12, w13, p11, p12
 real(8), intent(IN)    :: p13, dwe, dpe
 real(8), intent(IN)    :: c1(26,30),  c2(26,30),  c3(26,30)
 real(8), intent(IN)    :: oo1(26,21), oo2(26,21), oo3(26,21)
 real(8), intent(IN)    :: h11(26,31), h12(26,31), h13(26,31)
 real(8), intent(IN)    :: h21(26,31), h22(26,31), h23(26,31)
 real(8), intent(IN)    :: h81(26,31), h82(26,31), h83(26,31)

real(8), intent(in) :: wk_uv(5), zk_uv(5), ry_uv(5)
real(8), intent(in) :: xk_ir(10), ry_ir(3)
real(8), intent(in) :: cah(43,37), coa(62,101)

real(8), intent(in) :: aig_uv(3), awg_uv(3), arg_uv(3)
real(8), intent(in) :: aib_uv, awb_uv(2), arb_uv(2)
real(8), intent(in) :: aib_nir, awb_nir(3,2), arb_nir(3,2)
real(8), intent(in) :: aia_nir(3,3), awa_nir(3,3), ara_nir(3,3)
real(8), intent(in) :: aig_nir(3,3), awg_nir(3,3), arg_nir(3,3)
real(8), intent(in) :: caib(11,9,11), caif(9,11)

character(len=50) :: AEROSOLS(NA)

             !Names of aerosols as named in NL model
             AEROSOLS(1) = 'du001'
             AEROSOLS(2) = 'du002'
             AEROSOLS(3) = 'du003'
             AEROSOLS(4) = 'du004'
             AEROSOLS(5) = 'du005'


PTT = PTT*(QVT*10e-4 + O3T*4)
QVT = QVT*(PTT*10e-4 + O3T*4)
O3T = O3T*(QVT*10e-4 + PTT*4)
CFLST = CFLST*(QVT*10e-4 + O3T*4 + 1e-6*(QIT + QLT))
CFCNT = CFCNT*(QVT*10e-4 + O3T*4 + 1e-6*(QIT + QLT))
QIT = QIT*(QVT*10e-4 + O3T*4 + 1e-6*(CFCNT + CFLST))
QLT = QLT*(QVT*10e-4 + O3T*4 + 1e-6*(CFCNT + CFLST))
DU001T = QVT*10e-4 + O3T*4 + 1e-6*(DU001T + DU002T + DU003T + DU004T + DU005T )
DU002T = QVT*10e-4 + O3T*4 + 1e-6*(DU001T + DU002T + DU003T + DU004T + DU005T )
DU003T = QVT*10e-4 + O3T*4 + 1e-6*(DU001T + DU002T + DU003T + DU004T + DU005T )
DU004T = QVT*10e-4 + O3T*4 + 1e-6*(DU001T + DU002T + DU003T + DU004T + DU005T )
DU005T = QVT*10e-4 + O3T*4 + 1e-6*(DU001T + DU002T + DU003T + DU004T + DU005T )



call radiation_driver ( IM,JM,LM,RUNalarm,DT,NA,NBCHOU_IR,NBCHOU_SO,PLE,UT,PTT,QVT,O3T,CFLST,CFCNT,QIT,QLT,PS,TS, &
                              EMIS, DELT, COSZ, SLR, RGBUV, RGFUV, RGBIR, RGFIR , SC, &
                              DU001T, DU002T, DU003T, DU004T, DU005T, AEROSOLS, &
                              TAUA_IR_C, SSAA_IR_C, ASYA_IR_C, TAUA_SO_C, SSAA_SO_C, ASYA_SO_C, &
                              ILSF, ICNF, LLSF, LCNF, &
                              levs925, LCLDMH, LCLDLM, &
                              mapl_p00, mapl_rgas, mapl_cp, mapl_kappa, mapl_grav, &
                              aib_ir, awb_ir, aiw_ir,   &
                              aww_ir, aig_ir, awg_ir,   &
                              xkw, xke, mw, aw, bw, pm, &
                              fkw, gkw, cb, dcb,        &
                              w11, w12, w13, p11, p12,  &
                              p13, dwe, dpe,            &
                              c1,  c2,  c3,             &
                              oo1, oo2, oo3,            &
                              h11, h12, h13,            &
                              h21, h22, h23,            &
                              h81, h82, h83,             &
                              wk_uv, zk_uv, ry_uv,         &
                              xk_ir, ry_ir,                &
                              cah, coa,                    &
                              aig_uv, awg_uv, arg_uv,      &
                              aib_uv, awb_uv, arb_uv,      &
                              aib_nir, awb_nir, arb_nir,   &
                              aia_nir, awa_nir, ara_nir,   &
                              aig_nir, awg_nir, arg_nir,   &
                              caib, caif  ,                 &
                              HK_UV_TEMP, HK_IR_TEMP )


PTT = PTT*(QVT*10e-4 + O3T*4)
QVT = QVT*(PTT*10e-4 + O3T*4)
O3T = O3T*(QVT*10e-4 + PTT*4)
CFLST = CFLST*(QVT*10e-4 + O3T*4 + 1e-6*(QIT + QLT))
CFCNT = CFCNT*(QVT*10e-4 + O3T*4 + 1e-6*(QIT + QLT))
QIT = QIT*(QVT*10e-4 + O3T*4 + 1e-6*(CFCNT + CFLST))
QLT = QLT*(QVT*10e-4 + O3T*4 + 1e-6*(CFCNT + CFLST))
DU001T = QVT*10e-4 + O3T*4 + 1e-6*(DU001T + DU002T + DU003T + DU004T + DU005T )
DU002T = QVT*10e-4 + O3T*4 + 1e-6*(DU001T + DU002T + DU003T + DU004T + DU005T )
DU003T = QVT*10e-4 + O3T*4 + 1e-6*(DU001T + DU002T + DU003T + DU004T + DU005T )
DU004T = QVT*10e-4 + O3T*4 + 1e-6*(DU001T + DU002T + DU003T + DU004T + DU005T )
DU005T = QVT*10e-4 + O3T*4 + 1e-6*(DU001T + DU002T + DU003T + DU004T + DU005T )


end subroutine radiation_driver_dummy

subroutine radiation_driver ( IM,JM,LM,RUNalarm,DT,NA,NBCHOU_IR,NBCHOU_SO,PLE,UT,PTT,QVT,O3T,CFLST,CFCNT,QIT,QLT,PS,TS, &
                              EMIS, DELT, COSZ, SLR, RGBUV, RGFUV, RGBIR, RGFIR , SC, &
                              DU001T, DU002T, DU003T, DU004T, DU005T, AEROSOLS, &
                              TAUA_IR_C, SSAA_IR_C, ASYA_IR_C, TAUA_SO_C, SSAA_SO_C, ASYA_SO_C, &
                              ILSF, ICNF, LLSF, LCNF, &
                              levs925, LCLDMH, LCLDLM, &
                              mapl_p00, mapl_rgas, mapl_cp, mapl_kappa, mapl_grav, &
                              aib_ir, awb_ir, aiw_ir,   &
                              aww_ir, aig_ir, awg_ir,   &
                              xkw, xke, mw, aw, bw, pm, &
                              fkw, gkw, cb, dcb,        &
                              w11, w12, w13, p11, p12,  &
                              p13, dwe, dpe,            &
                              c1,  c2,  c3,             &
                              oo1, oo2, oo3,            &
                              h11, h12, h13,            &
                              h21, h22, h23,            &
                              h81, h82, h83,             &
                              wk_uv, zk_uv, ry_uv,         &
                              xk_ir, ry_ir,                &
                              cah, coa,                    &
                              aig_uv, awg_uv, arg_uv,      &
                              aib_uv, awb_uv, arb_uv,      &
                              aib_nir, awb_nir, arb_nir,   &
                              aia_nir, awa_nir, ara_nir,   &
                              aig_nir, awg_nir, arg_nir,   &
                              caib, caif  ,                 &
                              HK_UV_TEMP, HK_IR_TEMP )


IMPLICIT NONE

!In
integer, intent(in) :: IM,JM,LM,RUNalarm
integer, intent(in) :: NA, NBCHOU_IR, NBCHOU_SO
real(8), intent(in) :: SC, DT
real(8), intent(in), dimension(IM,JM,0:LM) :: PLE
real(8), intent(in), dimension(IM,JM,LM) :: UT, QVT, O3T, CFLST, CFCNT, QIT, QLT
real(8), intent(in), dimension(IM,JM) :: PS,TS
real(8), intent(in), dimension(IM,JM) :: EMIS, DELT, COSZ, SLR, RGBUV, RGFUV, RGBIR, RGFIR
real(8), intent(in), dimension(IM,JM,LM) :: DU001T, DU002T, DU003T, DU004T, DU005T

real(8), intent(in), dimension(IM,JM,LM) :: ILSF, ICNF, LLSF, LCNF

real(8), intent(in) ::  HK_IR_TEMP(3,10), HK_UV_TEMP(5)

real(8), intent(in), dimension(IM,JM,LM,NBCHOU_IR,NA) :: TAUA_IR_C, SSAA_IR_C, ASYA_IR_C
real(8), intent(in), dimension(IM,JM,LM,NBCHOU_SO,NA) :: TAUA_SO_C, SSAA_SO_C, ASYA_SO_C

integer, intent(in) :: levs925, LCLDMH, LCLDLM

real(8), intent(in) :: mapl_p00, mapl_rgas, mapl_cp, mapl_kappa, mapl_grav

 real(8), intent(IN)    :: aib_ir(3,10), awb_ir(4,10), aiw_ir(4,10)
 real(8), intent(IN)    :: aww_ir(4,10), aig_ir(4,10), awg_ir(4,10)

 integer, intent(IN) :: mw(9)
 real(8), intent(IN)    :: xkw(9), xke(9), aw(9), bw(9), pm(9)
 real(8), intent(IN)    :: fkw(6,9), gkw(6,3), cb(6,10), dcb(5,10)
 real(8), intent(IN)    :: w11, w12, w13, p11, p12
 real(8), intent(IN)    :: p13, dwe, dpe
 real(8), intent(IN)    :: c1(26,30),  c2(26,30),  c3(26,30)
 real(8), intent(IN)    :: oo1(26,21), oo2(26,21), oo3(26,21)
 real(8), intent(IN)    :: h11(26,31), h12(26,31), h13(26,31)
 real(8), intent(IN)    :: h21(26,31), h22(26,31), h23(26,31)
 real(8), intent(IN)    :: h81(26,31), h82(26,31), h83(26,31)

real(8), intent(in) :: wk_uv(5), zk_uv(5), ry_uv(5)
real(8), intent(in) :: xk_ir(10), ry_ir(3)
real(8), intent(in) :: cah(43,37), coa(62,101)

real(8), intent(in) :: aig_uv(3), awg_uv(3), arg_uv(3)
real(8), intent(in) :: aib_uv, awb_uv(2), arb_uv(2)
real(8), intent(in) :: aib_nir, awb_nir(3,2), arb_nir(3,2)
real(8), intent(in) :: aia_nir(3,3), awa_nir(3,3), ara_nir(3,3)
real(8), intent(in) :: aig_nir(3,3), awg_nir(3,3), arg_nir(3,3)
real(8), intent(in) :: caib(11,9,11), caif(9,11)

!Inouts
real(8), intent(inout), dimension(IM,JM,LM) :: PTT

!Locals
integer :: i,j,k,l
real(8) :: CO2
logical, parameter :: TRACE = .true., OVERCAST = .false. 
  integer, parameter :: NS = 1

real(8), dimension(IM,JM,LM,10) :: TAUA_IRT, SSAA_IRT, ASYA_IRT
real(8), dimension(IM,JM,LM,8)  :: TAUA_SOT, SSAA_SOT, ASYA_SOT
real(8), dimension(IM,JM,LM,NA)  :: AEROT, SAEROT
character(len=50), intent(in) :: AEROSOLS(NA)
real(8) :: X
integer :: AI

real(8), dimension(IM,JM,LM) :: PTT1, TT, PLOf, PIf, DELTAP, TT_OUT, DTDT
real(8), dimension(IM,JM,LM+1) :: PLESO
real(8), dimension(IM,JM,LM) :: O3mmT
real(8), dimension(IM,JM,LM) :: N2OT,CH4T,CFC11T,CFC12T,CFC22T
real(8), dimension(IM,JM,LM) :: RAD_QLT, RAD_QIT, RAD_CFT, RAD_RLT, RAD_RIT
real(8), dimension(IM,JM,LM,4) :: CWCT, REFFT
real(8), dimension(IM,JM,1)    :: FST, TGT, TVT
real(8), dimension(IM,JM,1,10) :: EGT, EVT, RVT

real(8), dimension(IM,JM,0:LM) :: FLWUT, FLWDT, FLXNT, FLWT_INT, FSWT_INT
real(8), dimension(IM,JM,0:LM) :: flwavet, fswavet, dfdtst_int

real(8), dimension(IM,JM,LM) :: QILST, QLLST, QICNT, QLCNT

real(8), dimension(IM,JM) :: T2MT, TEMPOR, COSZTMP


       DELTAP = ( PLE(:,:,1:LM) - PLE(:,:,0:LM-1) )

       !Pressure (hPa) and Exner pressure at the full levels
       PLOf(:,:,1:LM) = 0.01 * 0.5 * (PLE(:,:,0:LM-1) +  PLE(:,:,1:LM  ) )
       PIf(:,:,1:LM) = (PLOf(:,:,1:LM)/1000.)**(MAPL_RGAS/MAPL_CP)

       !Potential temperature with p0=10000
       PTT1 = PTT*(MAPL_P00**MAPL_KAPPA) !Place in holder so as not to overwrite

       !Temperature
       TT = PTT1*PIf !Save input for later.

IF (RUNALARM == 1) then


          !2m temperature
          T2MT = TT(:,:,LM)*(0.5*(1.0 + PLE(:,:,LM-1)/PLE(:,:,LM)))**(-MAPL_KAPPA)

          !Pressure in mb for SORAD
          PLESO(:,:,1:LM+1) = 0.01 * PLE(:,:,0:LM)

          TEMPOR = 0.
          do l = levs925, lm
             where (UT(:,:,l).gt.4.) TEMPOR(:,:) = 1.
          end do

          !Convert ozone from ppmv to kgkg 
          O3mmT = O3T / 1.e6

             CO2 = 1.0

          !Initialize other gases, not available right now.
          N2OT   = 0.0
          CH4T   = 0.0
          CFC11T = 0.0
          CFC12T = 0.0
          CFC22T = 0.0

          QILST = QIT * ILSF
          QICNT = QIT * ICNF
          QLLST = QLT * LLSF
          QLCNT = QLT * LCNF

          !Initialize RAD/CLOUD variables
          RAD_CFT = 0.0
          RAD_QLT = 0.0
          RAD_QIT = 0.0
          RAD_RLT = 0.0
          RAD_RIT = 0.0


          !Call RADcouple to produce cloud-rad variables
          DO I = 1,IM
             DO J = 1,JM
                DO L = 1,LM
                   call RADCOUPLE( TT(I,J,L), PLOf(I,J,L),                                   & 
                                   CFLST(I,J,L), CFCNT(I,J,L),                               &
                                   QLLST(I,J,L), QILST(I,J,L), QLCNT(I,J,L), QICNT(I,J,L),   & 
                                   RAD_QLT(I,J,L), RAD_QIT(I,J,L),                           &  
                                   RAD_CFT(I,J,L),                                           &
                                   RAD_RLT(I,J,L), RAD_RIT(I,J,L),                           &
                                   TEMPOR(I,J)                                               )
                endDO
             endDO
          endDO

          CWCT(:,:,:,1) = RAD_QIT !QI
          CWCT(:,:,:,2) = RAD_QLT !QL
          CWCT(:,:,:,3) = 0.0     !QR - this is not available right now
          CWCT(:,:,:,4) = 0.0     !QI - this is not available right now

          REFFT(:,:,:,1) = RAD_RIT * 1.0e6 !RI
          REFFT(:,:,:,2) = RAD_RLT * 1.0e6 !RL
          REFFT(:,:,:,3) = 100.e-6 * 1.0e6 !RR
          REFFT(:,:,:,4) = 140.e-6 * 1.0e6 !RS

          !Set surface quantities
          EGT = 0.0
          FST = 0.0
          TGT = 0.0
          TVT = 0.0
          EVT = 0.0
          RVT = 0.0

          do L = 1, 10
             EGT(:,:,1,L) = EMIS(:,:)
          end do
          FST         = 1.0
          TGT(:,:,1)  = TS
          TVT(:,:,1)  = TS
          EVT         = 0.0
          RVT         = 0.0

          if (NA == 0) then

             TAUA_IRT = 0.0
             SSAA_IRT = 0.0
             ASYA_IRT = 0.0
             TAUA_SOT = 0.0
             SSAA_SOT = 0.0
             ASYA_SOT = 0.0

          else

             !Place aerosol fields into single rank 4 array that can be looped over
             AEROT(:,:,:,1) = DU001T
             AEROT(:,:,:,2) = DU002T
             AEROT(:,:,:,3) = DU003T
             AEROT(:,:,:,4) = DU004T
             AEROT(:,:,:,5) = DU005T


             !Optical property calculation needs pressure as (1/g)dp/dz*q
             DO L = 1, LM
                DO J = 1, JM
                   DO I = 1, IM
                      X = ((PLE(I,J,L) - PLE(I,J,L-1))*0.01)*(100./MAPL_GRAV)
                      DO AI = 1, NA
                         SAEROT(I,J,L,AI) = X*AEROT(I,J,L,AI)
                      END DO
                   END DO
                END DO
             END DO

             TAUA_IRT = 0.0
             SSAA_IRT = 0.0
             ASYA_IRT = 0.0
             do J = 1,NBCHOU_IR
                do L = 1,NA
                   TAUA_IRT(:,:,:,J) = TAUA_IRT(:,:,:,J) + &
                                       SAEROT(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)
                   SSAA_IRT(:,:,:,J) = SSAA_IRT(:,:,:,J) + &
                                       SAEROT(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)*SSAA_IR_C(:,:,:,J,L)
                   ASYA_IRT(:,:,:,J) = ASYA_IRT(:,:,:,J) + &
                                       SAEROT(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)*SSAA_IR_C(:,:,:,J,L)*ASYA_IR_C(:,:,:,J,L)
                enddo
             enddo

             TAUA_SOT = 0.0
             SSAA_SOT = 0.0
             ASYA_SOT = 0.0
             do J = 1,NBCHOU_SO
                do L = 1,NA
                   TAUA_SOT(:,:,:,J) = TAUA_SOT(:,:,:,J) + &
                                       SAEROT(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)
                   SSAA_SOT(:,:,:,J) = SSAA_SOT(:,:,:,J) + &
                                       SAEROT(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)*SSAA_SO_C(:,:,:,J,L)
                   ASYA_SOT(:,:,:,J) = ASYA_SOT(:,:,:,J) + &
                                       SAEROT(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)*SSAA_SO_C(:,:,:,J,L)*ASYA_SO_C(:,:,:,J,L)
                enddo
             enddo

          endif

          do i = 1,IM
             do j = 1,JM
                COSZTMP(i,j) = max(.0001,COSZ(i,j))
             enddo
          enddo

          !Initialize the fluxes
          FLWUT = 0.0
          FLWDT = 0.0
          FLXNT = 0.0


          DO i = 1,IM
             DO j = 1,JM

            !SOLAR RADIATION (SHORT WAVE)
             call SORAD( 1,                          &
                           LM,                             &
                           NBCHOU_SO,                      &
                           COSZTMP,                        &
                           PLESO(I,J,:),                          &
                           TT(I,J,:),                          &
                           QVT(I,J,:),                        &
                           O3mmT(I,J,:),                    &
                           CO2,                            &
                           CWCT(I,J,:,:),                     &
                           RAD_CFT(I,J,:),                &
                           LCLDMH,                         &
                           LCLDLM,                         &
                           REFFT(I,J,:,:),                    &
                           HK_UV_TEMP,                     &
                           HK_IR_TEMP,                     &
                           TAUA_SOT(I,J,:,:),                       &
                           SSAA_SOT(I,J,:,:),                       &
                           ASYA_SOT(I,J,:,:),                       &
                           RGBUV,                          &
                           RGFUV,                          &
                           RGBIR,                          &
                           RGFIR,                          &
                           FLXNT(I,J,:),                    &
                           MAPL_GRAV,                      &
                           wk_uv, zk_uv, ry_uv,            &
                           xk_ir, ry_ir,                   &
                           cah, coa,                       &
                           aig_uv, awg_uv, arg_uv,         &
                           aib_uv, awb_uv, arb_uv,         &
                           aib_nir, awb_nir, arb_nir,      &
                           aia_nir, awa_nir, ara_nir,      &
                           aig_nir, awg_nir, arg_nir,      &
                           caib, caif                      )

   enddo
enddo


             FSWT_INT = FLXNT

DO i = 1,IM
   DO j = 1,JM
             !INFRARED RADIATION (LONG WAVE)
             call IRRAD(  LM,                             &
                           PLE(i,j,:),                            &
                           TT(i,j,:),                          &
                           QVT(i,j,:),                        &
                           O3mmT(i,j,:),                    &
                           T2MT(i,j),                      &
                           CO2,                            &
                           TRACE,                          &
                           N2OT(i,j,:),                           &
                           CH4T(i,j,:),                           &
                           CFC11T(i,j,:),                         &
                           CFC12T(i,j,:),                         &
                           CFC22T(i,j,:),                         &
                           CWCT(i,j,:,:),                      &
                           RAD_CFT(i,j,:),                &
                           LCLDMH,                         &
                           LCLDLM,                         &
                           REFFT(i,j,:,:),                    &
                           ns,                             &
                           FST(i,j,:),                            &
                           TGT(i,j,:),                            &
                           EGT(i,j,:,:),                            &
                           TVT(i,j,:),                            &
                           EVT(i,j,:,:),                            &
                           RVT(i,j,:,:),                            &
                           NA,                             &
                           NBCHOU_IR,                      &
                           TAUA_IRT(i,j,:,:),              &
                           SSAA_IRT(i,j,:,:),              &
                           ASYA_IRT(i,j,:,:),              &
                           FLWUT(i,j,:),                    &
                           FLWDT(i,j,:),                    &
                           DFDTST_INT(i,j,:),          &
                           aib_ir, awb_ir, aiw_ir,         &
                           aww_ir, aig_ir, awg_ir,         &
                           xkw, xke, mw, aw, bw, pm,       &
                           fkw, gkw, cb, dcb,              &
                           w11, w12, w13, p11, p12, p13,   &
                           dwe, dpe, c1, c2, c3,           &
                           oo1, oo2, oo3,                  &
                           h11, h12, h13, h21,             &
                           h22, h23, h81, h82, h83         )

   enddo
enddo


             FLWT_INT = FLWUT + FLWDT

endif

          do L = 0, LM
             FLWAVET(:,:,L) = FLWT_INT(:,:,L) + DFDTST_INT(:,:,L)*DELT
             FSWAVET(:,:,L) = FSWT_INT(:,:,L) * SLR
          end do
 
          !Tangent of fluxes to DTDt (mass-weighted)
          DTDt = ( (FLWAVET(:,:,0:LM-1) - FLWAVET(:,:,1:LM)) &
                  + (FSWAVET(:,:,0:LM-1) - FSWAVET(:,:,1:LM)) ) * (MAPL_GRAV/MAPL_CP)

          !Tangent of mass-weighted DTDt to updated temperature
          TT_out = TT + DT * DTDt / DELTAP

          !Tangent of temperature to potential temperature
          PTT1 = TT_out/PIf

          !Tangent of potential temperature to P0 = 1
          PTT = PTT1/(MAPL_P00**MAPL_KAPPA)



endsubroutine radiation_driver









subroutine RADCOUPLE(  TE,              & 
                       PL,              & 
                       CF,              & 
                       AF,              & 
                       QClLS,           & 
                       QCiLS,           & 
                       QClAN,           & 
                       QCiAN,           & 
                       RAD_QL,          &  
                       RAD_QI,          & 
                       RAD_CF,          & 
                       RAD_RL,          & 
                       RAD_RI,          & 
                       TEMPOR           )

 IMPLICIT NONE

 !Inputs
 real(8), intent(in ) :: TE, PL, TEMPOR
 real(8), intent(in ) :: AF, CF, QClAN, QCiAN, QClLS, QCiLS
! real(8), intent(in ) :: QRN_ALL, QSN_ALL

 !Outputs
 real(8), intent(out) :: RAD_QL,RAD_QI,RAD_CF,RAD_RL,RAD_RI
! real(8), intent(out) :: RAD_QR,RAD_QS

 !Locals
 real(8) :: ss, RAD_RI_AN, AFx, ALPH

 real(8), parameter :: MIN_RI = 20.e-6, MAX_RI = 40.e-6, RI_ANV = 30.e-6

 !Initialize outputs
 RAD_QL = 0.0
 RAD_QI = 0.0
 RAD_CF = 0.0
 RAD_RL = 0.0
 RAD_RI = 0.0
 !RAD_QR = 0.0
 !RAD_QS = 0.0

 ! Adjust Anvil fractions for warm clouds
 ALPH =  0.1
 SS   =  (280.-TE)/20.
 SS   =  MIN( 1.0 , SS )
 SS   =  MAX( 0.0 , SS )

 SS   =  ALPH + (SS**3) * ( 1.0 - ALPH )

 AFx  =  AF * SS * 0.5

 !Total cloud fraction
 RAD_CF = MIN( CF + AFx, 1.00 )

 !Total In-cloud liquid
 if ( RAD_CF > 10.0e-8 ) then  !0 -> 10e-8 FOR LINEARIZATION PROTECTION
    RAD_QL = ( QClLS + QClAN ) / RAD_CF
 else
    RAD_QL = 0.0
 end if
 RAD_QL = MIN( RAD_QL, 0.01 )

 ! Total In-cloud ice
 if (  RAD_CF > 10.0e-8 ) then !0 -> 10e-8 FOR LINEARIZATION PROTECTION
    RAD_QI = ( QCiLS + QCiAN ) / RAD_CF
 else
    RAD_QI = 0.0
 end if
 RAD_QI = MIN( RAD_QI, 0.01 )

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

 if (PL < 150. ) then
    RAD_RI = MAX_RI
 end if
 if (PL >= 150. ) then
    RAD_RI = MAX_RI*150./PL
 end if

 ! Weigh in a separate R_ice for Anvil Ice according to
 RAD_RI_AN  =  RAD_RI  

 if ( ( QCiLS + QCiAN ) > 0.0 ) then
    if (qcils/rad_ri+qcian/ri_anv .gt. 10e-8) then !LINEARIZATION PROTECTION
       RAD_RI_AN  = ( QCiLS + QCiAN ) / ( (QCiLS/RAD_RI) + (QCiAN/RI_ANV) )
    endif
 end if

 RAD_RI = MIN( RAD_RI, RAD_RI_AN )
 RAD_RI = MAX( RAD_RI, MIN_RI )

 ! Implement ramps for gradual change in effective radius
 if (PL < 300. ) then
    RAD_RL = 21.e-6
 end if
 if (PL >= 300. ) then
    RAD_RL = 21.e-6*300./PL
 end if
 RAD_RL = MAX( RAD_RL, 10.e-6 )

 ! Thicken low high lat clouds
 if ( PL .GE. 775.  .AND. TE .LE.  275. .AND. (tempor.eq.1.) ) then
    RAD_RL = max(min(-0.1 * PL + 87.5, 10.),5.)*1.e-6
 end if
 if ( PL .GE. 825.  .AND. TE .LE.  282. .AND. (tempor.eq.1.) ) then
    RAD_RL = max(0.71 * TE - 190.25, 5.)*1.e-6
 end if
 if ( PL .GE. 775.  .AND. PL .LT. 825. .AND. TE .LE.  282. .AND. TE .GT. 275. .AND. (tempor.eq.1.) ) then
    RAD_RL = min(-0.1*PL + 0.71 * TE - 107.75, 10.)*1.e-6
 end if
 if ( PL .GE. 825.  .AND. TE .LE.  275. .AND. (tempor.eq.1.) ) then
    RAD_RL = 5.*1.e-6
 end if

 ! Thin low tropical clouds
 if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
    RAD_RL = min(2.2 * TE - 617., 21.)*1.e-6
 end if
 if ( PL .GE. 925.  .AND. TE .GE.  290. ) then
    RAD_RL = min(0.44 * PL - 397., 21.)*1.e-6
 end if
 if ( PL .GE. 925.  .AND. PL .LT. 950. .AND. TE .GT.  285. .AND. TE .LT. 290.) then
    RAD_RL = max(min(0.44*PL + 2.2 * TE - 1035., 21.),10.)*1.e-6
 end if
 if ( PL .GE. 950.  .AND. TE .GE.  290. ) then
    RAD_RL = 21.*1.e-6
 end if

 if ( RAD_CF < 1.e-5 ) then
    RAD_QL = 0.
    RAD_QI = 0.
    RAD_CF = 0.
    !RAD_QR = 0.
    !RAD_QS = 0.
 end if

end subroutine RADCOUPLE


 subroutine irrad ( np          , & !Number of layers (LM)
                    ple_dev     , & !Pressure at level edges (Pa)
                    ta_dev      , & !d - Temperature (K)
                    wa_dev      , & !d - Specific humidity (g/g)
                    oa_dev      , & !d - Ozone (g/g)
                    tb_dev      , & !Surface air temperature (K)
                    co2         , & !* Carbon dioxide (pppv)
                    trace       , & !Option
                    n2o_dev     , & !* Nitrous oxide (pppv)
                    ch4_dev     , & !* Methane (pppv)
                    cfc11_dev   , & !* Trichlorofluoromethane (pppv)
                    cfc12_dev   , & !* Dichlorodifluoromethane (pppv)
                    cfc22_dev   , & !* Chlorodifluoromethane (pppv)
                    cwc_dev     , & !Cloud water mixing ratio (kg/kg) 
                    fcld_dev    , & !Cloud amount (fraction)
                    ict         , & !Level index separating high and middle clouds
                    icb         , & !Level index separating middle and low clouds
                    reff_dev    , & !Effective size of cloud particles (micron)
                    ns          , & !Number of sub-grid surface types
                    fs_dev      , & !Fractional cover of sub-grid regions
                    tg_dev      , & !Land or ocean surface temperature
                    eg_dev      , & !Land or ocean surface emissivity
                    tv_dev      , & !Vegetation temperature
                    ev_dev      , & !Vegetation emissivity
                    rv_dev      , & !Vegetation reflectivity 
                    na          , & !Number of bands
                    nb          , & !Number of bands in IRRAD calcs for Chou
                    taua_dev    , & !Aerosol optical thickness
                    ssaa_dev    , & !Aerosol single scattering albedo
                    asya_dev    , & !Aerosol asymmetry factor
                    flxu_dev    , & !Upwelling flux, all-sky
                    flxd_dev    , & !Downwelling flux, all-sky
                    dfdts_dev   , & !Sensitivity of net downward flux to surface temperature
                    aib_ir, awb_ir, aiw_ir,   &
                    aww_ir, aig_ir, awg_ir,   &
                    xkw, xke, mw, aw, bw, pm, &
                    fkw, gkw, cb, dcb,        &
                    w11, w12, w13, p11, p12,  &
                    p13, dwe, dpe,            &
                    c1,  c2,  c3,             &
                    oo1, oo2, oo3,            &
                    h11, h12, h13,            &
                    h21, h22, h23,            &
                    h81, h82, h83             )

 IMPLICIT NONE

 !Radiation constants, these need to be inputs for the autodiff tool
 real(8), intent(IN)    :: aib_ir(3,10), awb_ir(4,10), aiw_ir(4,10)
 real(8), intent(IN)    :: aww_ir(4,10), aig_ir(4,10), awg_ir(4,10)

 integer, intent(IN)    :: mw(9)
 real(8), intent(IN)    :: xkw(9), xke(9), aw(9), bw(9), pm(9)
 real(8), intent(IN)    :: fkw(6,9), gkw(6,3), cb(6,10), dcb(5,10)
 real(8), intent(IN)    :: w11, w12, w13, p11, p12
 real(8), intent(IN)    :: p13, dwe, dpe
 real(8), intent(IN)    :: c1(26,30),  c2(26,30),  c3(26,30)
 real(8), intent(IN)    :: oo1(26,21), oo2(26,21), oo3(26,21)
 real(8), intent(IN)    :: h11(26,31), h12(26,31), h13(26,31)
 real(8), intent(IN)    :: h21(26,31), h22(26,31), h23(26,31)
 real(8), intent(IN)    :: h81(26,31), h82(26,31), h83(26,31)

 !----- INPUTS -----
 INTEGER, INTENT(IN)                     :: np, ict, icb, ns, na, nb
 LOGICAL, INTENT(IN)                     :: trace 

 real(8), INTENT(IN)                        :: CO2

 !Rank 2 inputs
 real(8), INTENT(IN)          :: tb_dev

 !Rank 3 (Prognostic variables and tracers)
 real(8), DIMENSION(np),    INTENT(IN)    :: ta_dev, wa_dev, oa_dev, fcld_dev
 real(8), DIMENSION(np),    INTENT(IN)    :: n2o_dev, ch4_dev, cfc11_dev, cfc12_dev, cfc22_dev
 real(8), DIMENSION(np+1),  INTENT(IN)    :: ple_dev

 !Rank 3 (surface types)
 real(8), DIMENSION(ns),    INTENT(IN)    :: fs_dev, tg_dev, tv_dev
 real(8), DIMENSION(ns,10), INTENT(IN)    :: eg_dev, ev_dev, rv_dev

 !Rank 3 (diagnostic cloud parts)
 real(8), DIMENSION(np,4),  INTENT(IN)    :: cwc_dev, reff_dev

 !Rank 3 (aerosols)
 real(8), DIMENSION(np,nb), INTENT(INOUT) :: taua_dev, ssaa_dev, asya_dev


 !----- OUPUTS -----
 real(8), DIMENSION(np+1), INTENT(OUT)    :: flxu_dev
 real(8), DIMENSION(np+1), INTENT(OUT)    :: flxd_dev
 real(8), DIMENSION(np+1), INTENT(OUT)    :: dfdts_dev

 !----- LOCALS -----
 real(8), parameter :: CONS_GRAV = 9.80665

 integer, parameter :: nx1 = 26
 integer, parameter :: no1 = 21
 integer, parameter :: nc1 = 30
 integer, parameter :: nh1 = 31


 !Temporary arrays
 real(8) :: pa(0:np),dt(0:np)
 real(8) :: x1,x2,x3
 real(8) :: dh2o(0:np),dcont(0:np),dco2(0:np),do3(0:np)
 real(8) :: dn2o(0:np),dch4(0:np)
 real(8) :: df11(0:np),df12(0:np),df22(0:np)
 real(8) :: th2o(6),tcon(3),tco2(6)
 real(8) :: tn2o(4),tch4(4),tcom(6)
 real(8) :: tf11,tf12,tf22
 real(8) :: blayer(0:np+1),blevel(0:np+1)
 real(8) :: bd(0:np+1),bu(0:np+1)
 real(8) :: bs,dbs,rflxs
 real(8) :: dp(0:np)
 real(8) :: trant,tranal
 real(8) :: transfc(0:np+1)
 real(8) :: flxu(0:np+1),flxd(0:np+1)
 real(8) :: taerlyr(0:np)

 !OVERCAST
 integer :: ncld(3)
 integer :: icx(0:np)
 !OVERCAST

 integer :: idx, rc
 integer :: k,l,ip,iw,ibn,ik,iq,isb,k1,k2,ne

 real(8) :: enn(0:np)
 real(8) :: cldhi,cldmd,cldlw,tcldlyr(0:np),fclr
 real(8) :: x,xx,yy,p1,a1,b1,fk1,a2,b2,fk2
 real(8) :: w1,ff

 logical :: oznbnd,co2bnd,h2otable,conbnd,n2obnd
 logical :: ch4bnd,combnd,f11bnd,f12bnd,f22bnd,b10bnd
 logical :: do_aerosol 

 !Temp arrays and variables for consolidation of tables
 integer, parameter :: max_num_tables = 17
 real(8) :: exptbl(0:np,max_num_tables)
 type :: band_table
    integer :: start
    integer :: end
 end type band_table
 type(band_table) :: h2oexp
 type(band_table) :: conexp
 type(band_table) :: co2exp
 type(band_table) :: n2oexp
 type(band_table) :: ch4exp
 type(band_table) :: comexp
 type(band_table) :: f11exp
 type(band_table) :: f12exp
 type(band_table) :: f22exp

 !Variables for new getirtau routine
 real(8) :: dp_pa(np)
 real(8) :: fcld_col(np)
 real(8) :: reff_col(np,4)
 real(8) :: cwc_col(np,4)

 real(8) :: h2oexp_tmp(0:np,5),conexp_tmp(0:np),co2exp_tmp(0:np,6),n2oexp_tmp(0:np,2)

 !BEGIN CALCULATIONS ...

      !compute layer pressure (pa) and layer temperature minus 250K (dt) 
      do k=1,np
         pa(k) = 0.5*(ple_dev(k+1)+ple_dev(k))*0.01
         dp(k) =     (ple_dev(k+1)-ple_dev(k))*0.01
         dp_pa(k) =     (ple_dev(k+1)-ple_dev(k)) ! dp in Pascals for getirtau
         dt(k) = ta_dev(k)-250.0

         !compute layer absorber amount

         !dh2o : water vapor amount (g/cm**2)
         !dcont: scaled water vapor amount for continuum absorption
         !       (g/cm**2)
         !dco2 : co2 amount (cm-atm)stp
         !do3  : o3 amount (cm-atm)stp
         !dn2o : n2o amount (cm-atm)stp
         !dch4 : ch4 amount (cm-atm)stp
         !df11 : cfc11 amount (cm-atm)stp
         !df12 : cfc12 amount (cm-atm)stp
         !df22 : cfc22 amount (cm-atm)stp
         !the factor 1.02 is equal to 1000/980
         !factors 789 and 476 are for unit conversion
         !the factor 0.001618 is equal to 1.02/(.622*1013.25) 
         !the factor 6.081 is equal to 1800/296

         dh2o(k) = 1.02*wa_dev   (k)*dp(k)
         do3 (k) = 476.*oa_dev   (k)*dp(k)
         dco2(k) = 789.*co2           *dp(k)
         dch4(k) = 789.*ch4_dev  (k)*dp(k)
         dn2o(k) = 789.*n2o_dev  (k)*dp(k)
         df11(k) = 789.*cfc11_dev(k)*dp(k)
         df12(k) = 789.*cfc12_dev(k)*dp(k)
         df22(k) = 789.*cfc22_dev(k)*dp(k)

         dh2o(k) = max(dh2o(k),1.e-10)
         do3 (k) = max(do3 (k),1.e-6)
         dco2(k) = max(dco2(k),1.e-4)

         !Compute scaled water vapor amount for h2o continuum absorption
         !following eq. (4.21).

         xx=pa(k)*0.001618*wa_dev(k)*wa_dev(k)*dp(k)
         dcont(k) = xx*exp(1800./ta_dev(k)-6.081)

         !Fill the reff, cwc, and fcld for the column   
         fcld_col(k) = fcld_dev(k)
         do l = 1, 4
            reff_col(k,l) = reff_dev(k,l)
            cwc_col(k,l) = cwc_dev(k,l)
         end do

      end do

      !A layer is added above the top of the model atmosphere.
      !Index "0" is the layer above the top of the atmosphere.

      dp  (0) = max(ple_dev(1)*0.01,0.005)
      pa  (0) = 0.5*dp(0)
      dt  (0) = ta_dev(1)-250.0

      dh2o(0) = 1.02*wa_dev   (1)*dp(0)
      do3 (0) = 476.*oa_dev   (1)*dp(0)
      dco2(0) = 789.*co2           *dp(0)
      dch4(0) = 789.*ch4_dev  (1)*dp(0)
      dn2o(0) = 789.*n2o_dev  (1)*dp(0)
      df11(0) = 789.*cfc11_dev(1)*dp(0)
      df12(0) = 789.*cfc12_dev(1)*dp(0)
      df22(0) = 789.*cfc22_dev(1)*dp(0)

      dh2o(0) = max(dh2o(0),1.e-10)
      do3 (0) = max(do3(0),1.e-6)
      dco2(0) = max(dco2(0),1.e-4)

      xx=pa(0)*0.001618*wa_dev(1)*wa_dev(1)*dp(0)
      dcont(0) = xx*exp(1800./ta_dev(1)-6.081)

      !The surface (np+1) is treated as a layer filled with black clouds.
      !transfc is the transmittance between the surface and a pressure
      !level.
      transfc(np+1)=1.0

      !Initialize fluxes
      do k=1,np+1
         flxu_dev(k)  = 0.0
         flxd_dev(k)  = 0.0
         dfdts_dev(k)= 0.0
      end do

      !Integration over spectral bands
      do ibn=1,10

         if (ibn == 10 .and. .not. trace) return

         !if h2otable, compute h2o (line) transmittance using table look-up.
         !if conbnd,   compute h2o (continuum) transmittance in bands 2-7.
         !if co2bnd,   compute co2 transmittance in band 3.
         !if oznbnd,   compute  o3 transmittance in band 5.
         !if n2obnd,   compute n2o transmittance in bands 6 and 7.
         !if ch4bnd,   compute ch4 transmittance in bands 6 and 7.
         !if combnd,   compute co2-minor transmittance in bands 4 and 5.
         !if f11bnd,   compute cfc11 transmittance in bands 4 and 5.
         !if f12bnd,   compute cfc12 transmittance in bands 4 and 6.
         !if f22bnd,   compute cfc22 transmittance in bands 4 and 6.
         !if b10bnd,   compute flux reduction due to n2o in band 10.

         h2otable=ibn == 1 .or. ibn == 2 .or. ibn == 8
         conbnd  =ibn >= 2 .and. ibn <= 7
         co2bnd  =ibn == 3
         oznbnd  =ibn == 5
         n2obnd  =ibn == 6 .or. ibn == 7
         ch4bnd  =ibn == 6 .or. ibn == 7
         combnd  =ibn == 4 .or. ibn == 5
         f11bnd  =ibn == 4 .or. ibn == 5
         f12bnd  =ibn == 4 .or. ibn == 6
         f22bnd  =ibn == 4 .or. ibn == 6
         b10bnd  =ibn == 10

         do_aerosol = na > 0

         exptbl = 0.0

         !Control packing of the new exponential tables by band

         select case (ibn)
         case (2)
            conexp%start = 1
            conexp%end   = 1
         case (3)
            h2oexp%start = 1
            h2oexp%end   = 6
            conexp%start = 7
            conexp%end   = 9
         case (4)
            h2oexp%start = 1
            h2oexp%end   = 6
            conexp%start = 7
            conexp%end   = 7
            comexp%start = 8
            comexp%end   = 13
            f11exp%start = 14
            f11exp%end   = 14
            f12exp%start = 15
            f12exp%end   = 15
            f22exp%start = 16
            f22exp%end   = 16
         case (5)
            h2oexp%start = 1
            h2oexp%end   = 6
            conexp%start = 7
            conexp%end   = 7
            comexp%start = 8
            comexp%end   = 13
            f11exp%start = 14
            f11exp%end   = 14
         case (6)
            h2oexp%start = 1
            h2oexp%end   = 6
            conexp%start = 7
            conexp%end   = 7
            n2oexp%start = 8
            n2oexp%end   = 11
            ch4exp%start = 12
            ch4exp%end   = 15
            f12exp%start = 16
            f12exp%end   = 16
            f22exp%start = 17
            f22exp%end   = 17
         case (7)
            h2oexp%start = 1
            h2oexp%end   = 6
            conexp%start = 7
            conexp%end   = 7
            n2oexp%start = 8
            n2oexp%end   = 11
            ch4exp%start = 12
            ch4exp%end   = 15
         case (9)
            h2oexp%start = 1
            h2oexp%end   = 6
         case (10)
            h2oexp%start = 1
            h2oexp%end   = 5
            conexp%start = 6
            conexp%end   = 6
            co2exp%start = 7
            co2exp%end   = 12
            n2oexp%start = 13
            n2oexp%end   = 14
         end select

         !blayer is the spectrally integrated planck flux of the mean layer
         !temperature derived from eq. (3.11)
         !The fitting for the planck flux is valid for the range 160-345 K.

         do k=1,np
            call planck(ibn,cb,ta_dev(k),blayer(k))
         end do

         !Index "0" is the layer above the top of the atmosphere.
         blayer(0)=blayer(1)
         blevel(0)=blayer(1)

         !Surface emission and reflectivity. See Section 9.
         !bs and dbs include the effect of surface emissivity.
         call sfcflux (ibn,cb,dcb,ns,fs_dev,tg_dev,eg_dev,tv_dev,ev_dev,rv_dev,&
               bs,dbs,rflxs) 

         blayer(np+1)=bs

         !interpolate Planck function at model levels (linear in p)

         do k=2,np
            blevel(k)=(blayer(k-1)*dp(k)+blayer(k)*dp(k-1))/(dp(k-1)+dp(k))
         end do

         !Extrapolate blevel(1) from blayer(2) and blayer(1)

         blevel(1)=blayer(1)+(blayer(1)-blayer(2))*dp(1)/(dp(1)+dp(2))
         blevel(0)=blevel(1)

         !If the surface air temperature tb is known, compute blevel(np+1)
         call planck(ibn,cb,tb_dev,blevel(np+1))

         !Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
         !     NOTE: dp_pa is only dims(1:np) as the 0'th level isn't needed in getirtau.
         !           Plus, the pressures in getirtau *MUST* be in Pascals.
         !     Slots for reff, hydrometeors and tauall are as follows:
         !                 1         Cloud Ice
         !                 2         Cloud Liquid
         !                 3         Falling Liquid (Rain)
         !                 4         Falling Ice (Snow)

         call getirtau1(ibn,np,dp_pa,fcld_col,reff_col,cwc_col,&
                       tcldlyr,enn,aib_ir,awb_ir, &
                       aiw_ir, aww_ir, aig_ir, awg_ir, CONS_GRAV)

         do k=0,np
            icx(k) = k
         end do

         call mkicx(np,ict,icb,enn,icx,ncld)

         !Compute optical thickness, single-scattering albedo and asymmetry
         !factor for a mixture of "na" aerosol types. Eqs. (7.1)-(7.3)
         if (do_aerosol) then
            taerlyr(0)=1.0

            do k=1,np

              !-----taerlyr is the aerosol diffuse transmittance
               taerlyr(k)=1.0
               if (taua_dev(k,ibn) > 0.001) then 
                  if (ssaa_dev(k,ibn) > 0.001) then
                     asya_dev(k,ibn)=asya_dev(k,ibn)/ssaa_dev(k,ibn)
                     ssaa_dev(k,ibn)=ssaa_dev(k,ibn)/taua_dev(k,ibn)

                     !Parameterization of aerosol scattering following
                     ff=.5+(.3739+(0.0076+0.1185*asya_dev(k,ibn))*asya_dev(k,ibn))*asya_dev(k,ibn)
                     taua_dev(k,ibn)=taua_dev(k,ibn)*(1.-ssaa_dev(k,ibn)*ff)
                  end if
                  taerlyr(k)=exp(-1.66*taua_dev(k,ibn))
               end if
            end do
         end if

         !Compute the exponential terms (Eq. 8.21) at each layer due to
         !water vapor line absorption when k-distribution is used
         if (.not. h2otable .and. .not. b10bnd) then
            call h2oexps(ibn,np,dh2o,pa,dt,xkw,aw,bw,pm,mw,exptbl(:,h2oexp%start:h2oexp%end))
         end if

         !compute the exponential terms (Eq. 4.24) at each layer due to
         !water vapor continuum absorption.

         ne=0
         if (conbnd) then
            ne=1
            if (ibn == 3) ne=3
            call conexps(ibn,np,dcont,xke,exptbl(:,conexp%start:conexp%end))
         end if


         if (trace) then

            !compute the exponential terms at each layer due to n2o absorption
            if (n2obnd) then
               call n2oexps(ibn,np,dn2o,pa,dt,exptbl(:,n2oexp%start:n2oexp%end))
            end if

            !compute the exponential terms at each layer due to ch4 absorption
            if (ch4bnd) then
               call ch4exps(ibn,np,dch4,pa,dt,exptbl(:,ch4exp%start:ch4exp%end))
            end if

            !Compute the exponential terms due to co2 minor absorption
            if (combnd) then
               call comexps(ibn,np,dco2,dt,exptbl(:,comexp%start:comexp%end))
            end if

            !Compute the exponential terms due to cfc11 absorption.
            !The values of the parameters are given in Table 7.
            if (f11bnd) then
               a1  = 1.26610e-3
               b1  = 3.55940e-6
               fk1 = 1.89736e+1
               a2  = 8.19370e-4
               b2  = 4.67810e-6
               fk2 = 1.01487e+1
               call cfcexps(ibn,np,a1,b1,fk1,a2,b2,fk2,df11,dt,&
                     exptbl(:,f11exp%start:f11exp%end))
            end if

            !Compute the exponential terms due to cfc12 absorption.
            if (f12bnd) then
               a1  = 8.77370e-4
               b1  =-5.88440e-6
               fk1 = 1.58104e+1
               a2  = 8.62000e-4
               b2  =-4.22500e-6
               fk2 = 3.70107e+1
               call cfcexps(ibn,np,a1,b1,fk1,a2,b2,fk2,df12,dt,&
                     exptbl(:,f12exp%start:f12exp%end))
            end if

            !Compute the exponential terms due to cfc22 absorption.
            if (f22bnd) then
               a1  = 9.65130e-4
               b1  = 1.31280e-5
               fk1 = 6.18536e+0
               a2  =-3.00010e-5 
               b2  = 5.25010e-7
               fk2 = 3.27912e+1
               call cfcexps(ibn,np,a1,b1,fk1,a2,b2,fk2,df22,dt,&
                     exptbl(:,f22exp%start:f22exp%end))
            end if

            !Compute the exponential terms at each layer in band 10 due to
            !h2o line and continuum, co2, and n2o absorption

            if (b10bnd) then
               call b10exps(np,dh2o,dcont,dco2,dn2o,pa,dt, &
                     h2oexp_tmp,exptbl(:,conexp%start:conexp%end),co2exp_tmp,n2oexp_tmp)

               exptbl(:,h2oexp%start:h2oexp%end) = h2oexp_tmp
               exptbl(:,co2exp%start:co2exp%end) = co2exp_tmp
               exptbl(:,n2oexp%start:n2oexp%end) = n2oexp_tmp

            end if
         end if

         !blayer(np+1) includes the effect of surface emissivity.

         bu(0)=0.0 ! ALT: this was undefined, check with Max if 0.0 is good value
         bd(0)=blayer(1)
         bu(np+1)=blayer(np+1)

         !do-loop 1500 is for computing upward (bu) and downward (bd)
         !Here, trant is the transmittance of the layer k2-1.
         do k2=1,np+1

            !for h2o line transmission
            if (.not. h2otable) then
               th2o=1.0
            end if

            !for h2o continuum transmission
            tcon=1.0

            x1=0.0
            x2=0.0
            x3=0.0
            trant=1.0

            if (h2otable) then

               !Compute water vapor transmittance using table look-up.
               if (ibn == 1) then
                  call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                        w11,p11,dwe,dpe,h11,h12,h13,trant)
               end if
               if (ibn == 2) then
                  call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                        w11,p11,dwe,dpe,h21,h22,h23,trant)
               end if
               if (ibn == 8) then
                  call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                        w11,p11,dwe,dpe,h81,h82,h83,trant)
               end if

               !for water vapor continuum absorption
               if (conbnd) then
                  tcon(1)=tcon(1)*exptbl(k2-1,conexp%start) ! Only the first exp
                  trant=trant*tcon(1)
               end if

            else

               !compute water vapor transmittance using k-distribution
               if (.not. b10bnd) then
                  call h2okdis(ibn,np,k2-1,fkw,gkw,ne,&
                        exptbl(:,h2oexp%start:h2oexp%end), &
                        exptbl(:,conexp%start:conexp%end), &
                        th2o,tcon,trant)
               end if

            end if

            if (co2bnd) then
               !Compute co2 transmittance using table look-up method
               call tablup(nx1,nc1,dco2(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                     w12,p12,dwe,dpe,c1,c2,c3,trant)
            end if

            !Always use table look-up to compute o3 transmittance.
            if (oznbnd) then
               call tablup(nx1,no1,do3(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                     w13,p13,dwe,dpe,oo1,oo2,oo3,trant)
            end if

            !include aerosol effect
            if (do_aerosol) then
               trant=trant*taerlyr(k2-1)
            end if

            !Compute upward (bu) and downward (bd) emission of the layer k2-1
            xx=(1.-enn(k2-1))*trant
            yy=min(0.9999,xx)
            yy=max(0.00001,yy)
            xx=(blevel(k2-1)-blevel(k2))/log(yy)
            bd(k2-1)=(blevel(k2)-blevel(k2-1)*yy)/(1.0-yy)-xx
            bu(k2-1)=(blevel(k2-1)+blevel(k2))-bd(k2-1)

         end do

         !Initialize fluxes
         flxu = 0.0
         flxd = 0.0

         !Compute upward and downward fluxes for each spectral band, ibn.
         do k1=0,np

            !initialization
            cldlw = 0.0
            cldmd = 0.0
            cldhi = 0.0
            tranal= 1.0

            !for h2o line transmission
            if (.not. h2otable) then
               th2o=1.0
            end if

            !for h2o continuum transmission
            tcon=1.0

            if (trace) then !Add trace gases contribution

               if (n2obnd) then !n2o
                  tn2o=1.0
               end if

               if (ch4bnd) then !ch4
                  tch4=1.0
               end if

               if (combnd) then !co2-minor
                  tcom=1.0
               end if

               if (f11bnd) then !cfc-11
                  tf11=1.0
               end if

               if (f12bnd) then !cfc-12
                  tf12=1.0
               end if

               if (f22bnd) then !cfc-22
                  tf22=1.0
               end if

               if (b10bnd) then !
                  th2o=1.0

                  tco2=1.0

                  tcon(1)=1.0

                  tn2o=1.0
               end if

            end if

            x1=0.0
            x2=0.0
            x3=0.0

            do k2=k1+1,np+1

               trant=1.0
               fclr =1.0

               if (h2otable) then

                  !Compute water vapor transmittance using table look-up.
                  if (ibn == 1) then
                     call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                           w11,p11,dwe,dpe,h11,h12,h13,trant)
                  end if
                  if (ibn == 2) then
                     call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                           w11,p11,dwe,dpe,h21,h22,h23,trant)
                  end if
                  if (ibn == 8) then
                     call tablup(nx1,nh1,dh2o(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                           w11,p11,dwe,dpe,h81,h82,h83,trant)
                  end if

                  if (conbnd) then
                     tcon(1)=tcon(1)*exptbl(k2-1,conexp%start) ! Only the first exp
                     trant=trant*tcon(1)
                  end if

               else

                  !Compute water vapor transmittance using k-distribution
                  if (.not. b10bnd) then
                     call h2okdis(ibn,np,k2-1,fkw,gkw,ne,&
                           exptbl(:,h2oexp%start:h2oexp%end), &
                           exptbl(:,conexp%start:conexp%end), &
                           th2o,tcon,trant)
                  end if

               end if

               if (co2bnd) then

                  !Compute co2 transmittance using table look-up method.
                  call tablup(nx1,nc1,dco2(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                        w12,p12,dwe,dpe,c1,c2,c3,trant)
               end if

               if (oznbnd) then
                  !Always use table look-up to compute o3 transmittanc
                  call tablup(nx1,no1,do3(k2-1),pa(k2-1),dt(k2-1),x1,x2,x3, &
                        w13,p13,dwe,dpe,oo1,oo2,oo3,trant)
               end if

               if (trace) then !Add trace gas effects

                  if (n2obnd) then !n2o
                     call n2okdis(ibn,np,k2-1,exptbl(:,n2oexp%start:n2oexp%end),tn2o,trant)
                  end if

                  if (ch4bnd) then !ch4
                     call ch4kdis(ibn,np,k2-1,exptbl(:,ch4exp%start:ch4exp%end),tch4,trant)
                  end if

                  if (combnd) then !co2-minor
                     call comkdis(ibn,np,k2-1,exptbl(:,comexp%start:comexp%end),tcom,trant)
                  end if

                  if (f11bnd) then !cfc11
                     call cfckdis(np,k2-1,exptbl(:,f11exp%start:f11exp%end),tf11,trant)
                  end if

                  if (f12bnd) then !cfc12
                     call cfckdis(np,k2-1,exptbl(:,f12exp%start:f12exp%end),tf12,trant)
                  end if

                  if (f22bnd) then !CFC22
                     call cfckdis(np,k2-1,exptbl(:,f22exp%start:f22exp%end),tf22,trant)
                  end if

                  if (b10bnd) then
                     call b10kdis(np,k2-1,&
                           exptbl(:,h2oexp%start:h2oexp%end),&
                           exptbl(:,conexp%start:conexp%end),&
                           exptbl(:,co2exp%start:co2exp%end),&
                           exptbl(:,n2oexp%start:n2oexp%end),&
                           th2o,tcon,tco2,tn2o,trant)
                  end if

               end if


               if (do_aerosol) then
                  tranal=tranal*taerlyr(k2-1)
                  trant=trant *tranal
               end if

               if (enn(k2-1) >= 0.001) then
                  call cldovlp (np,k1,k2,ict,icb,icx,ncld,enn,tcldlyr,cldhi,cldmd,cldlw)
               end if

               fclr=(1.0-cldhi)*(1.0-cldmd)*(1.0-cldlw)

               if (k2 == k1+1 .and. ibn /= 10) then
                  flxu(k1)=flxu(k1)-bu(k1)
                  flxd(k2)=flxd(k2)+bd(k1)
               end if

               xx=trant*(bu(k2-1)-bu(k2))
               flxu(k1)=flxu(k1)+xx*fclr

               if (k1 == 0) then !mjs  bd(-1) is not defined
                  xx=-trant*bd(k1)
               else
                  xx= trant*(bd(k1-1)-bd(k1))
               end if
               flxd(k2)=flxd(k2)+xx*fclr

            end do

            transfc(k1) =trant*fclr

            if (k1 > 0)then
               dfdts_dev(k1) =dfdts_dev(k1)-dbs*transfc(k1)
            end if

         end do

         if (.not. b10bnd) then !surface emission

            flxu(np+1)       =             -blayer(np+1)
            dfdts_dev(np+1)= dfdts_dev(np+1)-dbs

            do k=1,np+1
               flxu(k) = flxu(k)-flxd(np+1)*transfc(k)*rflxs
            end do

         end if

         do k=1,np+1 !Summation of fluxes over spectral bands

            flxu_dev(k) = flxu_dev(k) + flxu(k)
            flxd_dev(k) = flxd_dev(k) + flxd(k)

         end do

      end do

   end subroutine irrad

!***********************************************************************
   subroutine planck(ibn,cb,t,xlayer)
!***********************************************************************
!
!-----Compute spectrally integrated Planck flux
!
   implicit none

   integer, intent(in) :: ibn                   ! spectral band index
   real(8), intent(in) :: cb(6,10)                 ! Planck table coefficients
   real(8), intent(in) :: t                        ! temperature (K)
   real(8), intent(out) :: xlayer                   ! planck flux (w/m2)

   xlayer=t*(t*(t*(t*(t*cb(6,ibn)+cb(5,ibn))+cb(4,ibn))+cb(3,ibn))+cb(2,ibn))+cb(1,ibn)

   end subroutine planck


!***********************************************************************
   subroutine plancd(ibn,dcb,t,dbdt) 
!***********************************************************************
!
!-----Compute the derivative of Planck flux wrt temperature
!
   implicit none

   integer, intent(in) :: ibn               ! spectral band index
   real(8), intent(in) :: dcb(5,10)            ! Planck function derivative coefficients
   real(8), intent(in) :: t                    ! temperature (K)
   real(8), intent(out) :: dbdt                 ! derivative of Planck flux wrt temperature

   dbdt=t*(t*(t*(t*dcb(5,ibn)+dcb(4,ibn))+dcb(3,ibn))+dcb(2,ibn))+dcb(1,ibn)

   end subroutine plancd


!**********************************************************************
   subroutine h2oexps(ib,np,dh2o,pa,dt,xkw,aw,bw,pm,mw,h2oexp)
!**********************************************************************
!   Compute exponentials for water vapor line absorption
!   in individual layers using Eqs. (8.21) and (8.22).
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  layer water vapor amount for line absorption (dh2o) 
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!  absorption coefficients for the first k-distribution
!     function due to h2o line absorption (xkw)
!  coefficients for the temperature and pressure scaling (aw,bw,pm)
!  ratios between neighboring absorption coefficients for
!     h2o line absorption (mw)
!
!---- output parameters
!  6 exponentials for each layer  (h2oexp)
!**********************************************************************
   implicit none

   integer :: ib,np,ik,k

!---- input parameters ------

   real(8) :: dh2o(0:np),pa(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: h2oexp(0:np,6)

!---- static data -----

   integer :: mw(9)
   real(8) :: xkw(9),aw(9),bw(9),pm(9)

!---- temporary arrays -----

   real(8) :: xh

!**********************************************************************
!    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
!    and bw,  therefore, h2oexp for these sub-bands are identical.
!**********************************************************************

   do k=0,np

!-----xh is the scaled water vapor amount for line absorption
!     computed from Eq. (4.4).

      xh = dh2o(k)*(pa(k)/500.)**pm(ib) &
            * ( 1.+(aw(ib)+bw(ib)* dt(k))*dt(k) )

!-----h2oexp is the water vapor transmittance of the layer k
!     due to line absorption

      h2oexp(k,1) = exp(-xh*xkw(ib))

!-----compute transmittances from Eq. (8.22)

      do ik=2,6
         if (mw(ib) == 6) then
            xh = h2oexp(k,ik-1)*h2oexp(k,ik-1)
            h2oexp(k,ik) = xh*xh*xh
         elseif (mw(ib) == 8) then
            xh = h2oexp(k,ik-1)*h2oexp(k,ik-1)
            xh = xh*xh
            h2oexp(k,ik) = xh*xh
         elseif (mw(ib) == 9) then
            xh=h2oexp(k,ik-1)*h2oexp(k,ik-1)*h2oexp(k,ik-1)
            h2oexp(k,ik) = xh*xh*xh
         else
            xh = h2oexp(k,ik-1)*h2oexp(k,ik-1)
            xh = xh*xh
            xh = xh*xh
            h2oexp(k,ik) = xh*xh
         end if
      end do
   end do

   end subroutine h2oexps


!**********************************************************************
   subroutine conexps(ib,np,dcont,xke,conexp)
!**********************************************************************
!   compute exponentials for continuum absorption in individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  layer scaled water vapor amount for continuum absorption (dcont) 
!  absorption coefficients for the first k-distribution function
!     due to water vapor continuum absorption (xke)
!
!---- output parameters
!  1 or 3 exponentials for each layer (conexp)
!**********************************************************************
   implicit none

   integer :: ib,np,k

!---- input parameters ------

   real(8) :: dcont(0:np)

!---- updated parameters -----

   real(8) :: conexp(0:np,3)

!---- static data -----

   real(8) :: xke(9)

!****************************************************************

   do k=0,np

      conexp(k,1) = exp(-dcont(k)*xke(ib))

!-----The absorption coefficients for sub-bands 3b and 3a are, respectively,
!     two and four times the absorption coefficient for sub-band 3c (Table 9).
!     Note that conexp(3) is for sub-band 3a. 

      if (ib  ==  3) then
         conexp(k,2) = conexp(k,1) *conexp(k,1)
         conexp(k,3) = conexp(k,2) *conexp(k,2)
      end if

   end do

   end subroutine conexps


!**********************************************************************
   subroutine n2oexps(ib,np,dn2o,pa,dt,n2oexp)
!**********************************************************************
!   Compute n2o exponentials for individual layers 
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  layer n2o amount (dn2o)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!  2 or 4 exponentials for each layer (n2oexp)
!**********************************************************************
   implicit none

   integer :: ib,k,np

!---- input parameters -----

   real(8) :: dn2o(0:np),pa(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: n2oexp(0:np,4)

!---- temporary arrays -----

   real(8) :: xc,xc1,xc2

!-----Scaling and absorption data are given in Table 5.
!     Transmittances are computed using Eqs. (8.21) and (8.22).

   do k=0,np

!-----four exponential by powers of 21 for band 6.

      if (ib == 6) then

         xc=dn2o(k)*(1.+(1.9297e-3+4.3750e-6*dt(k))*dt(k))
         n2oexp(k,1)=exp(-xc*6.31582e-2)

         xc=n2oexp(k,1)*n2oexp(k,1)*n2oexp(k,1)
         xc1=xc*xc
         xc2=xc1*xc1
         n2oexp(k,2)=xc*xc1*xc2

!-----four exponential by powers of 8 for band 7

      else

         xc=dn2o(k)*(pa(k)/500.0)**0.48 &
               *(1.+(1.3804e-3+7.4838e-6*dt(k))*dt(k))
         n2oexp(k,1)=exp(-xc*5.35779e-2)

         xc=n2oexp(k,1)*n2oexp(k,1)
         xc=xc*xc
         n2oexp(k,2)=xc*xc
         xc=n2oexp(k,2)*n2oexp(k,2)
         xc=xc*xc
         n2oexp(k,3)=xc*xc
         xc=n2oexp(k,3)*n2oexp(k,3)
         xc=xc*xc
         n2oexp(k,4)=xc*xc
      end if

   end do

   end subroutine n2oexps


!**********************************************************************
   subroutine ch4exps(ib,np,dch4,pa,dt,ch4exp)
!**********************************************************************
!   Compute ch4 exponentials for individual layers
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  layer ch4 amount (dch4)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!  1 or 4 exponentials for each layer (ch4exp)
!**********************************************************************
   implicit none

   integer :: ib,np,k

!---- input parameters -----

   real(8) :: dch4(0:np),pa(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: ch4exp(0:np,4)

!---- temporary arrays -----

   real(8) :: xc

!*****  Scaling and absorption data are given in Table 5  *****

   do k=0,np

!-----four exponentials for band 6

      if (ib == 6) then

         xc=dch4(k)*(1.+(1.7007e-2+1.5826e-4*dt(k))*dt(k))
         ch4exp(k,1)=exp(-xc*5.80708e-3)

!-----four exponentials by powers of 12 for band 7

      else

         xc=dch4(k)*(pa(k)/500.0)**0.65 &
               *(1.+(5.9590e-4-2.2931e-6*dt(k))*dt(k))
         ch4exp(k,1)=exp(-xc*6.29247e-2)

         xc=ch4exp(k,1)*ch4exp(k,1)*ch4exp(k,1)
         xc=xc*xc
         ch4exp(k,2)=xc*xc

         xc=ch4exp(k,2)*ch4exp(k,2)*ch4exp(k,2)
         xc=xc*xc
         ch4exp(k,3)=xc*xc

         xc=ch4exp(k,3)*ch4exp(k,3)*ch4exp(k,3)
         xc=xc*xc
         ch4exp(k,4)=xc*xc

      end if

   end do

   end subroutine ch4exps


!**********************************************************************
   subroutine comexps(ib,np,dcom,dt,comexp)
!**********************************************************************
!   Compute co2-minor exponentials for individual layers using 
!   Eqs. (8.21) and (8.22).
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  layer co2 amount (dcom)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!  6 exponentials for each layer (comexp)
!**********************************************************************
   implicit none

   integer :: ib,ik,np,k

!---- input parameters -----

   real(8) :: dcom(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: comexp(0:np,6)

!---- temporary arrays -----

   real(8) :: xc

!*****  Scaling and absorpton data are given in Table 6  *****

   do k=0,np

      if (ib == 4) then
         xc=dcom(k)*(1.+(3.5775e-2+4.0447e-4*dt(k))*dt(k))
      end if

      if (ib == 5) then
         xc=dcom(k)*(1.+(3.4268e-2+3.7401e-4*dt(k))*dt(k))
      end if

      comexp(k,1)=exp(-xc*1.922e-7)

      do ik=2,6
         xc=comexp(k,ik-1)*comexp(k,ik-1)
         xc=xc*xc
         comexp(k,ik)=xc*comexp(k,ik-1)
      end do

   end do

   end subroutine comexps


!**********************************************************************
   subroutine cfcexps(ib,np,a1,b1,fk1,a2,b2,fk2,dcfc,dt,cfcexp)
!**********************************************************************
!   compute cfc(-11, -12, -22) exponentials for individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of layers (np)
!  parameters for computing the scaled cfc amounts
!             for temperature scaling (a1,b1,a2,b2)
!  the absorption coefficients for the
!     first k-distribution function due to cfcs (fk1,fk2)
!  layer cfc amounts (dcfc)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!  1 exponential for each layer (cfcexp)
!**********************************************************************
   implicit none

   integer :: ib,np,k

!---- input parameters -----

   real(8) :: dcfc(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: cfcexp(0:np)

!---- static data -----

   real(8) :: a1,b1,fk1,a2,b2,fk2

!---- temporary arrays -----

   real(8) :: xf

!**********************************************************************

   do k=0,np

!-----compute the scaled cfc amount (xf) and exponential (cfcexp)

      if (ib == 4) then
         xf=dcfc(k)*(1.+(a1+b1*dt(k))*dt(k))
         cfcexp(k)=exp(-xf*fk1)
      else
         xf=dcfc(k)*(1.+(a2+b2*dt(k))*dt(k))
         cfcexp(k)=exp(-xf*fk2)
      end if

   end do

   end subroutine cfcexps


!**********************************************************************
   subroutine b10exps(np,dh2o,dcont,dco2,dn2o,pa,dt, &
      h2oexp,conexp,co2exp,n2oexp)
!**********************************************************************
!   Compute band3a exponentials for individual layers
!
!---- input parameters
!  number of layers (np)
!  layer h2o amount for line absorption (dh2o)
!  layer h2o amount for continuum absorption (dcont)
!  layer co2 amount (dco2)
!  layer n2o amount (dn2o)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!
!  exponentials for each layer (h2oexp,conexp,co2exp,n2oexp)
!**********************************************************************
   implicit none

   integer :: np,k

!---- input parameters -----

   real(8) :: dh2o(0:np),dcont(0:np),dn2o(0:np)
   real(8) :: dco2(0:np),pa(0:np),dt(0:np)

!---- output parameters -----

   real(8) :: h2oexp(0:np,5),conexp(0:np),co2exp(0:np,6),n2oexp(0:np,2)

!---- temporary arrays -----

   real(8) :: xx,xx1,xx2,xx3

!**********************************************************************

   do k=0,np

!-----Compute scaled h2o-line amount for Band 10 (Eq. 4.4 and Table 3).

      xx=dh2o(k)*(pa(k)/500.0) &
            *(1.+(0.0149+6.20e-5*dt(k))*dt(k))

!-----six exponentials by powers of 8

      h2oexp(k,1)=exp(-xx*0.10624)

      xx=h2oexp(k,1)*h2oexp(k,1)
      xx=xx*xx
      h2oexp(k,2)=xx*xx

      xx=h2oexp(k,2)*h2oexp(k,2)
      xx=xx*xx
      h2oexp(k,3)=xx*xx

      xx=h2oexp(k,3)*h2oexp(k,3)
      xx=xx*xx
      h2oexp(k,4)=xx*xx

      xx=h2oexp(k,4)*h2oexp(k,4)
      xx=xx*xx
      h2oexp(k,5)=xx*xx

!-----one exponential of h2o continuum for sub-band 3a (Table 9).

      conexp(k)=exp(-dcont(k)*109.0)

!-----Scaled co2 amount for the Band 10 (Eq. 4.4, Tables 3 and 6).

      xx=dco2(k)*(pa(k)/300.0)**0.5 &
            *(1.+(0.0179+1.02e-4*dt(k))*dt(k))

!-----six exponentials by powers of 8

      co2exp(k,1)=exp(-xx*2.656e-5)

      xx=co2exp(k,1)*co2exp(k,1)
      xx=xx*xx
      co2exp(k,2)=xx*xx

      xx=co2exp(k,2)*co2exp(k,2)
      xx=xx*xx
      co2exp(k,3)=xx*xx

      xx=co2exp(k,3)*co2exp(k,3)
      xx=xx*xx
      co2exp(k,4)=xx*xx

      xx=co2exp(k,4)*co2exp(k,4)
      xx=xx*xx
      co2exp(k,5)=xx*xx

      xx=co2exp(k,5)*co2exp(k,5)
      xx=xx*xx
      co2exp(k,6)=xx*xx

!-----Compute the scaled n2o amount for Band 10 (Table 5).

      xx=dn2o(k)*(1.+(1.4476e-3+3.6656e-6*dt(k))*dt(k))

!-----Two exponentials by powers of 58

      n2oexp(k,1)=exp(-xx*0.25238)

      xx=n2oexp(k,1)*n2oexp(k,1)
      xx1=xx*xx
      xx1=xx1*xx1
      xx2=xx1*xx1
      xx3=xx2*xx2
      n2oexp(k,2)=xx*xx1*xx2*xx3

   end do

   end subroutine b10exps


!**********************************************************************
   subroutine tablup(nx1,nh1,dw,p,dt,s1,s2,s3,w1,p1, &
      dwe,dpe,coef1,coef2,coef3,tran)
!**********************************************************************
!   Compute water vapor, co2 and o3 transmittances between level
!   k1 and and level k2 , using table look-up.
!
!   Calculations follow Eq. (4.16).
!
!---- input ---------------------
!
!  number of pressure intervals in the table (nx1)
!  number of absorber amount intervals in the table (nh1)
!  layer absorber amount (dw)
!  layer pressure in mb (p)
!  deviation of layer temperature from 250K (dt)
!  first value of absorber amount (log10) in the table (w1) 
!  first value of pressure (log10) in the table (p1) 
!  size of the interval of absorber amount (log10) in the table (dwe)
!  size of the interval of pressure (log10) in the table (dpe)
!  pre-computed coefficients (coef1, coef2, and coef3)
!
!---- updated ---------------------
!
!  column integrated absorber amount (s1)
!  absorber-weighted column pressure (s2)
!  absorber-weighted column temperature (s3)
!  transmittance (tran)
!
!  Note: Units of s1 are g/cm**2 for water vapor and
!       (cm-atm)stp for co2 and o3.
!   
!**********************************************************************
   implicit none

   integer :: nx1,nh1

!---- input parameters -----

   real(8) :: w1,p1,dwe,dpe
   real(8) :: dw,p,dt
   real(8) :: coef1(nx1,nh1),coef2(nx1,nh1),coef3(nx1,nh1)

!---- update parameter -----

   real(8) :: s1,s2,s3,tran

!---- temporary variables -----

   real(8) :: we,pe,fw,fp,pa,pb,pc,ax,ba,bb,t1,ca,cb,t2
   real(8) :: x1,x2,x3,xx, x1c
   integer :: iw,ip

!-----Compute effective pressure (x2) and temperature (x3) following 
!     Eqs. (8.28) and (8.29)

   s1=s1+dw
   s2=s2+p*dw
   s3=s3+dt*dw

   x1=s1
   x1c=1.0/s1
   x2=s2*x1c
   x3=s3*x1c

!-----normalize we and pe

!       we=(log10(x1)-w1)/dwe
!       pe=(log10(x2)-p1)/dpe
   we=(log10(x1)-w1)*dwe
   pe=(log10(x2)-p1)*dpe

!-----restrict the magnitudes of the normalized we and pe.

   we=min(we,real(nh1-1))
   pe=min(pe,real(nx1-1))

!-----assign iw and ip and compute the distance of we and pe 
!     from iw and ip.

   iw=int(we+1.0)
   iw=min(iw,nh1-1)
   iw=max(iw, 2)
   fw=we-real(iw-1)

   ip=int(pe+1.0)
   ip=min(ip,nx1-1)
   ip=max(ip, 1)
   fp=pe-real(ip-1)

!-----linear interpolation in pressure

   pa = coef1(ip,iw-1)+(coef1(ip+1,iw-1)-coef1(ip,iw-1))*fp
   pb = coef1(ip,  iw)+(coef1(ip+1,  iw)-coef1(ip,  iw))*fp
   pc = coef1(ip,iw+1)+(coef1(ip+1,iw+1)-coef1(ip,iw+1))*fp

!-----quadratic interpolation in absorber amount for coef1

   ax = ( (pc+pa)*fw + (pc-pa) ) *fw*0.5 + pb*(1.-fw*fw)

!-----linear interpolation in absorber amount for coef2 and coef3

   ba = coef2(ip,  iw)+(coef2(ip+1,  iw)-coef2(ip,  iw))*fp
   bb = coef2(ip,iw+1)+(coef2(ip+1,iw+1)-coef2(ip,iw+1))*fp
   t1 = ba+(bb-ba)*fw

   ca = coef3(ip,  iw)+(coef3(ip+1,  iw)-coef3(ip,  iw))*fp
   cb = coef3(ip,iw+1)+(coef3(ip+1,iw+1)-coef3(ip,iw+1))*fp
   t2 = ca + (cb-ca)*fw

!-----update the total transmittance between levels k1 and k2

   xx=(ax + (t1+t2*x3) * x3)
   xx= min(xx,0.9999999)
   xx= max(xx,0.0000001)
   tran= tran*xx

   end subroutine tablup


!**********************************************************************
   subroutine h2okdis(ib,np,k,fkw,gkw,ne,h2oexp,conexp, &
      th2o,tcon,tran)
!**********************************************************************
!   compute water vapor transmittance between levels k1 and k2, using the k-distribution method.
!
!---- input parameters
!  spectral band (ib)
!  number of levels (np)
!  current level (k)
!  planck-weighted k-distribution function due to
!    h2o line absorption (fkw)
!  planck-weighted k-distribution function due to
!    h2o continuum absorption (gkw)
!  number of terms used in each band to compute water vapor
!     continuum transmittance (ne)
!  exponentials for line absorption (h2oexp) 
!  exponentials for continuum absorption (conexp) 
!
!---- updated parameters
!  transmittance between levels k1 and k2 due to
!    water vapor line absorption (th2o)
!  transmittance between levels k1 and k2 due to
!    water vapor continuum absorption (tcon)
!  total transmittance (tran)
!
!**********************************************************************
   implicit none

!---- input parameters ------

   integer :: ib,ne,np,k
   real(8) :: h2oexp(0:np,6),conexp(0:np,3)
   real(8) ::  fkw(6,9),gkw(6,3)

!---- updated parameters -----

   real(8) :: th2o(6),tcon(3),tran

!---- temporary arrays -----

   real(8) :: trnth2o

!-----tco2 are the six exp factors between levels k1 and k2 
!     tran is the updated total transmittance between levels k1 and k2

!-----th2o is the 6 exp factors between levels k1 and k2 due to
!     h2o line absorption. 

!-----tcon is the 3 exp factors between levels k1 and k2 due to
!     h2o continuum absorption.

!-----trnth2o is the total transmittance between levels k1 and k2 due
!     to both line and continuum absorption.

!-----Compute th2o following Eq. (8.23).

   th2o(1) = th2o(1)*h2oexp(k,1)
   th2o(2) = th2o(2)*h2oexp(k,2)
   th2o(3) = th2o(3)*h2oexp(k,3)
   th2o(4) = th2o(4)*h2oexp(k,4)
   th2o(5) = th2o(5)*h2oexp(k,5)
   th2o(6) = th2o(6)*h2oexp(k,6)

   if (ne == 0) then

!-----Compute trnh2o following Eq. (8.25). fkw is given in Table 4.

      trnth2o =(fkw(1,ib)*th2o(1) &
              + fkw(2,ib)*th2o(2) &
              + fkw(3,ib)*th2o(3) &
              + fkw(4,ib)*th2o(4) &
              + fkw(5,ib)*th2o(5) &
              + fkw(6,ib)*th2o(6))


      tran    = tran*trnth2o

   elseif (ne == 1) then

!-----Compute trnh2o following Eqs. (8.25) and (4.27).

      tcon(1) = tcon(1)*conexp(k,1)

      trnth2o =(fkw(1,ib)*th2o(1) &
              + fkw(2,ib)*th2o(2) &
              + fkw(3,ib)*th2o(3) &
              + fkw(4,ib)*th2o(4) &
              + fkw(5,ib)*th2o(5) &
              + fkw(6,ib)*th2o(6))*tcon(1)

      tran    = tran*trnth2o

   else

!-----For band 3. This band is divided into 3 subbands.

      tcon(1)= tcon(1)*conexp(k,1)
      tcon(2)= tcon(2)*conexp(k,2)
      tcon(3)= tcon(3)*conexp(k,3)

!-----Compute trnh2o following Eqs. (4.29) and (8.25).

      trnth2o = (  gkw(1,1)*th2o(1) &
              + gkw(2,1)*th2o(2) &
              + gkw(3,1)*th2o(3) &
              + gkw(4,1)*th2o(4) &
              + gkw(5,1)*th2o(5) &
              + gkw(6,1)*th2o(6) ) * tcon(1) &
              + (  gkw(1,2)*th2o(1) &
              + gkw(2,2)*th2o(2) &
              + gkw(3,2)*th2o(3) &
              + gkw(4,2)*th2o(4) &
              + gkw(5,2)*th2o(5) &
              + gkw(6,2)*th2o(6) ) * tcon(2) &
              + (  gkw(1,3)*th2o(1) &
              + gkw(2,3)*th2o(2) &
              + gkw(3,3)*th2o(3) &
              + gkw(4,3)*th2o(4) &
              + gkw(5,3)*th2o(5) &
              + gkw(6,3)*th2o(6) ) * tcon(3)

      tran    = tran*trnth2o

   end if

   end subroutine h2okdis


!**********************************************************************
   subroutine n2okdis(ib,np,k,n2oexp,tn2o,tran)
!**********************************************************************
!   compute n2o transmittances between levels k1 and k2, using the k-distribution method with linear
!    pressure scaling.
!
!---- input parameters
!   spectral band (ib)
!   number of levels (np)
!   current level (k)
!   exponentials for n2o absorption (n2oexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to n2o absorption
!     for the various values of the absorption coefficient (tn2o)
!   total transmittance (tran)
!
!**********************************************************************
   implicit none
   integer :: ib,np,k

!---- input parameters -----

   real(8) :: n2oexp(0:np,4)

!---- updated parameters -----

   real(8) :: tn2o(4),tran

!---- temporary arrays -----

   real(8) :: xc

!-----tn2o is computed from Eq. (8.23). 
!     xc is the total n2o transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.

!-----band 6

   if (ib == 6) then

      tn2o(1)=tn2o(1)*n2oexp(k,1)
      xc=   0.940414*tn2o(1)

      tn2o(2)=tn2o(2)*n2oexp(k,2)
      xc=xc+0.059586*tn2o(2)

!-----band 7

   else

      tn2o(1)=tn2o(1)*n2oexp(k,1)
      xc=   0.561961*tn2o(1)

      tn2o(2)=tn2o(2)*n2oexp(k,2)
      xc=xc+0.138707*tn2o(2)

      tn2o(3)=tn2o(3)*n2oexp(k,3)
      xc=xc+0.240670*tn2o(3)

      tn2o(4)=tn2o(4)*n2oexp(k,4)
      xc=xc+0.058662*tn2o(4)

   end if

   tran=tran*xc

   end subroutine n2okdis


!**********************************************************************
   subroutine ch4kdis(ib,np,k,ch4exp,tch4,tran)
!**********************************************************************
!   compute ch4 transmittances between levels k1 and k2, using the k-distribution method with
!    linear pressure scaling.
!
!---- input parameters
!   spectral band (ib)
!   number of levels (np)
!   current level (k)
!   exponentials for ch4 absorption (ch4exp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to ch4 absorption
!     for the various values of the absorption coefficient (tch4)
!   total transmittance (tran)
!
!**********************************************************************
   implicit none
   integer :: ib,np,k

!---- input parameters -----

   real(8) :: ch4exp(0:np,4)

!---- updated parameters -----

   real(8) :: tch4(4),tran

!---- temporary arrays -----

   real(8) :: xc

!-----tch4 is computed from Eq. (8.23). 
!     xc is the total ch4 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.

!-----band 6

   if (ib == 6) then

      tch4(1)=tch4(1)*ch4exp(k,1)
      xc= tch4(1)

!-----band 7

   else

      tch4(1)=tch4(1)*ch4exp(k,1)
      xc=   0.610650*tch4(1)

      tch4(2)=tch4(2)*ch4exp(k,2)
      xc=xc+0.280212*tch4(2)

      tch4(3)=tch4(3)*ch4exp(k,3)
      xc=xc+0.107349*tch4(3)

      tch4(4)=tch4(4)*ch4exp(k,4)
      xc=xc+0.001789*tch4(4)

   end if

   tran=tran*xc

   end subroutine ch4kdis


!**********************************************************************
   subroutine comkdis(ib,np,k,comexp,tcom,tran)
!**********************************************************************
!  compute co2-minor transmittances between levels k1 and k2, using the k-distribution method
!   with linear pressure scaling.
!
!---- input parameters
!   spectral band (ib)
!   number of levels (np)
!   current level (k)
!   exponentials for co2-minor absorption (comexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to co2-minor absorption
!     for the various values of the absorption coefficient (tcom)
!   total transmittance (tran)
!
!**********************************************************************
   implicit none
   integer :: ib,np,k

!---- input parameters -----

   real(8) :: comexp(0:np,6)

!---- updated parameters -----

   real(8) :: tcom(6),tran

!---- temporary arrays -----

   real(8) :: xc

!-----tcom is computed from Eq. (8.23). 
!     xc is the total co2 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 6.

!-----band 4

   if (ib == 4) then

      tcom(1)=tcom(1)*comexp(k,1)
      xc=   0.12159*tcom(1)
      tcom(2)=tcom(2)*comexp(k,2)
      xc=xc+0.24359*tcom(2)
      tcom(3)=tcom(3)*comexp(k,3)
      xc=xc+0.24981*tcom(3)
      tcom(4)=tcom(4)*comexp(k,4)
      xc=xc+0.26427*tcom(4)
      tcom(5)=tcom(5)*comexp(k,5)
      xc=xc+0.07807*tcom(5)
      tcom(6)=tcom(6)*comexp(k,6)
      xc=xc+0.04267*tcom(6)

!-----band 5

   else

      tcom(1)=tcom(1)*comexp(k,1)
      xc=   0.06869*tcom(1)
      tcom(2)=tcom(2)*comexp(k,2)
      xc=xc+0.14795*tcom(2)
      tcom(3)=tcom(3)*comexp(k,3)
      xc=xc+0.19512*tcom(3)
      tcom(4)=tcom(4)*comexp(k,4)
      xc=xc+0.33446*tcom(4)
      tcom(5)=tcom(5)*comexp(k,5)
      xc=xc+0.17199*tcom(5)
      tcom(6)=tcom(6)*comexp(k,6)
      xc=xc+0.08179*tcom(6)

   end if

   tran=tran*xc

   end subroutine comkdis


!**********************************************************************
   subroutine cfckdis(np,k,cfcexp,tcfc,tran)
!**********************************************************************
!  compute cfc-(11,12,22) transmittances between levels k1 and k2, using the k-distribution method with
!   linear pressure scaling.
!
!---- input parameters
!   number of levels (np)
!   current level (k)
!   exponentials for cfc absorption (cfcexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to cfc absorption
!     for the various values of the absorption coefficient (tcfc)
!   total transmittance (tran)
!
!**********************************************************************
   implicit none

!---- input parameters -----

   integer :: k, np
   real(8) :: cfcexp(0:np)

!---- updated parameters -----

   real(8) :: tcfc,tran

!-----tcfc is the exp factors between levels k1 and k2. 

   tcfc=tcfc*cfcexp(k)
   tran=tran*tcfc

   end subroutine cfckdis


!**********************************************************************
   subroutine b10kdis(np,k,h2oexp,conexp,co2exp,n2oexp, &
      th2o,tcon,tco2,tn2o,tran)
!**********************************************************************
!
!   compute h2o (line and continuum),co2,n2o transmittances between
!   levels k1 and k2, using the k-distribution
!   method with linear pressure scaling.
!
!---- input parameters
!   number of levels (np)
!   current level (k)
!   exponentials for h2o line absorption (h2oexp)
!   exponentials for h2o continuum absorption (conexp)
!   exponentials for co2 absorption (co2exp)
!   exponentials for n2o absorption (n2oexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to h2o line absorption
!     for the various values of the absorption coefficient (th2o)
!   transmittance between levels k1 and k2 due to h2o continuum
!     absorption for the various values of the absorption
!     coefficient (tcon)
!   transmittance between levels k1 and k2 due to co2 absorption
!     for the various values of the absorption coefficient (tco2)
!   transmittance between levels k1 and k2 due to n2o absorption
!     for the various values of the absorption coefficient (tn2o)
!   total transmittance (tran)
!
!**********************************************************************
   implicit none

   integer :: np,k

!---- input parameters -----

   real(8) :: h2oexp(0:np,5),conexp(0:np),co2exp(0:np,6), &
         n2oexp(0:np,2)

!---- updated parameters -----

   real(8) :: th2o(6),tcon(3),tco2(6),tn2o(4), &
         tran

!---- temporary arrays -----

   real(8) :: xx

!-----For h2o line. The k-distribution functions are given in Table 4.

   th2o(1)=th2o(1)*h2oexp(k,1)
   xx=   0.3153*th2o(1)

   th2o(2)=th2o(2)*h2oexp(k,2)
   xx=xx+0.4604*th2o(2)

   th2o(3)=th2o(3)*h2oexp(k,3)
   xx=xx+0.1326*th2o(3)

   th2o(4)=th2o(4)*h2oexp(k,4)
   xx=xx+0.0798*th2o(4)

   th2o(5)=th2o(5)*h2oexp(k,5)
   xx=xx+0.0119*th2o(5)

   tran=xx

!-----For h2o continuum. Note that conexp(k,3) is for subband 3a.

   tcon(1)=tcon(1)*conexp(k)
   tran=tran*tcon(1)

!-----For co2 (Table 6)

   tco2(1)=tco2(1)*co2exp(k,1)
   xx=    0.2673*tco2(1)

   tco2(2)=tco2(2)*co2exp(k,2)
   xx=xx+ 0.2201*tco2(2)

   tco2(3)=tco2(3)*co2exp(k,3)
   xx=xx+ 0.2106*tco2(3)

   tco2(4)=tco2(4)*co2exp(k,4)
   xx=xx+ 0.2409*tco2(4)

   tco2(5)=tco2(5)*co2exp(k,5)
   xx=xx+ 0.0196*tco2(5)

   tco2(6)=tco2(6)*co2exp(k,6)
   xx=xx+ 0.0415*tco2(6)

   tran=tran*xx

!-----For n2o (Table 5)

   tn2o(1)=tn2o(1)*n2oexp(k,1)
   xx=   0.970831*tn2o(1)

   tn2o(2)=tn2o(2)*n2oexp(k,2)
   xx=xx+0.029169*tn2o(2)

   tran=tran*(xx-1.0)

   end subroutine b10kdis


!mjs

!***********************************************************************
   subroutine cldovlp (np,k1,k2,ict,icb,icx,ncld,enn,ett, &
      cldhi,cldmd,cldlw)
!***********************************************************************
!
!     update the effective superlayer cloud fractions between  levels k1
!     and k2 following Eqs.(6.18)-(6.21). This assumes that the input
!     are the fractions between k1 and k2-1.
!
! input parameters
!
!  np:      number of layers
!  k1:      index for top level
!  k2:      index for bottom level
!  ict:     the level separating high and middle clouds
!  icb:     the level separating middle and low clouds
!  icx:     indeces sorted by enn in the three superlayers
!  enn:     effective fractional cloud cover of a layer
!  ett:     transmittance of a cloud layer
!
! update paremeters
!
!  cldhi:   Effective high-level cloudiness
!  cldmd:   Effective middle-level cloudiness
!  cldlw:   Effective low-level cloudiness
!
!
!***********************************************************************

   implicit none

   integer, intent(IN   ) :: np,k1,k2,ict,icb,icx(0:np),ncld(3)
   real(8),    intent(IN   ) :: enn(0:np), ett(0:np)

   real(8),    intent(INOUT) :: cldhi,cldmd,cldlw

   integer :: j,k,km,kx

   km=k2-1

   if    (km <  ict               ) then ! do high clouds

      kx=ncld(1)

      if(kx==1 .or. cldhi==0.) then
         cldhi = enn(km)
      else
      !if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
         cldhi = 0.0
         if(kx/=0) then
            do k=ict-kx,ict-1
               j=icx(k)
               if(j>=k1 .and. j<=km) cldhi=enn(j)+ett(j)*cldhi
            end do
         end if
      end if

   elseif(km >= ict .and. km < icb) then ! do middle clouds

      kx=ncld(2)

      if(kx==1 .or. cldmd==0.) then
         cldmd = enn(km)
      else
      !if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
         cldmd = 0.0
         if(kx/=0) then
            do k=icb-kx,icb-1
               j=icx(k)
               if(j>=k1 .and. j<=km) cldmd=enn(j)+ett(j)*cldmd
            end do
         end if
      end if

   else                                  ! do low clouds

      kx=ncld(3)

      if(kx==1 .or. cldlw==0.) then
         cldlw = enn(km)
      else
      !if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
         cldlw = 0.0
         if(kx/=0) then
            do k=np+1-kx,np
               j=icx(k)
               if(j>=k1 .and. j<=km) cldlw=enn(j)+ett(j)*cldlw
            end do
         end if
      end if

   end if

   end subroutine cldovlp

!***********************************************************************
   subroutine sfcflux (ibn,cb,dcb,ns,fs,tg,eg,tv,ev,rv,bs,dbs,rflxs)
!***********************************************************************
! Compute emission and reflection by an homogeneous/inhomogeneous 
!  surface with vegetation cover.
!
!-----Input parameters
!  index for the spectral band (ibn)
!  number of sub-grid box (ns)
!  fractional cover of sub-grid box (fs)
!  sub-grid ground temperature (tg)
!  sub-grid ground emissivity (eg)
!  sub-grid vegetation temperature (tv)
!  sub-grid vegetation emissivity (ev)
!  sub-grid vegetation reflectivity (rv)
!                      
!-----Output parameters                                Unit
!  Emission by the surface (ground+vegetation) (bs)    W/m^2
!  Derivative of bs rwt to temperature (dbs)           W/m^2
!  Reflection by the surface (rflxs)                   W/m^2
!
!**********************************************************************
   implicit none

!---- input parameters -----
   integer :: ibn,ns
   real(8) :: cb(6,10)
   real(8) :: dcb(5,10)
   real(8) :: fs(ns),tg(ns),eg(ns,10)
   real(8) :: tv(ns),ev(ns,10),rv(ns,10)

!---- output parameters -----
   real(8) :: bs,dbs,rflxs

!---- temporary arrays -----
   integer :: j
   real(8) :: bg(ns),bv(ns),dbg(ns),dbv(ns)
   real(8) :: xx,yy,zz

!***************************************************************

!-----compute planck flux and the derivative wrt to temperature
!     c.f. eqs. (3.11) and (3.12).

   do j=1,ns
      call planck(ibn,cb,tg(j),bg(j))
      call planck(ibn,cb,tv(j),bv(j))
      call plancd(ibn,dcb,tg(j),dbg(j))
      call plancd(ibn,dcb,tv(j),dbv(j))
   end do

!-----

   if (fs(1) > 0.9999) then
      if (ev(1,ibn) < 0.0001 .and. rv(1,ibn) < 0.0001) then

!-----for homogeneous surface without vegetation
!     following Eqs. (9.4), (9.5), and (3.13)

         bs =eg(1,ibn)*bg(1)
         dbs=eg(1,ibn)*dbg(1)
         rflxs=1.0-eg(1,ibn)
      else

!-----With vegetation
!     following Eqs. (9.1), (9.3), and (9.13)

         xx=ev(1,ibn)*bv(1)
         yy=1.0-ev(1,ibn)-rv(1,ibn)
         zz=1.0-eg(1,ibn)
         bs=yy*(eg(1,ibn)*bg(1)+zz*xx)+xx

         xx=ev(1,ibn)*dbv(1)
         dbs=yy*(eg(1,ibn)*dbg(1)+zz*xx)+xx

         rflxs=rv(1,ibn)+zz*yy*yy/(1.0-rv(1,ibn)*zz)
      end if
   else

!-----for nonhomogeneous surface

      bs=0.0
      dbs=0.0
      rflxs=0.0

      if (ev(1,ibn) < 0.0001 .and. rv(1,ibn) < 0.0001) then

!-----No vegetation, following Eqs. (9.9), (9.10), and (9.13)

         do j=1,ns
            bs=bs+fs(j)*eg(j,ibn)*bg(j)
            dbs=dbs+fs(j)*eg(j,ibn)*dbg(j)
            rflxs=rflxs+fs(j)*(1.0-eg(j,ibn))
         end do
      else

!-----With vegetation, following Eqs. (9.6), (9.7), and (9.13)

         do j=1,ns
            xx=ev(j,ibn)*bv(j)
            yy=1.0-ev(j,ibn)-rv(j,ibn)
            zz=1.0-eg(j,ibn)
            bs=bs+fs(j)*(yy*(eg(j,ibn)*bg(j)+zz*xx)+xx)

            xx=ev(j,ibn)*dbv(j)
            dbs=dbs+fs(j)*(yy*(eg(j,ibn)*dbg(j)+zz*xx)+xx)

            rflxs=rflxs+fs(j)*(rv(j,ibn)+zz*yy*yy &
                  /(1.0-rv(j,ibn)*zz))
         end do
      end if
   end if

   end subroutine sfcflux

!mjs


!***********************************************************************
   subroutine mkicx(np,ict,icb,enn,icx,ncld)
!***********************************************************************

   implicit none

   integer, intent(IN   ) :: np, ict, icb
   real(8),    intent(IN   ) :: enn(0:np)
   integer, intent(INOUT) :: icx(0:np)
   integer, intent(  OUT) :: ncld(3)

   call sortit(enn(  0:ict-1),icx(  0:ict-1),ncld(1),  0, ict-1)
   call sortit(enn(ict:icb-1),icx(ict:icb-1),ncld(2),ict, icb-1)
   call sortit(enn(icb:np   ),icx(icb:   np),ncld(3),icb,    np)

   end subroutine mkicx

!***********************************************************************
   subroutine SORTIT(ENN,LST,NCL,ibg,iend)
!***********************************************************************

   implicit none
   
   integer, intent(  OUT) :: NCL
   integer, intent(IN   ) :: IBG, IEND
   real(8),    intent(IN   ) :: ENN(IBG:IEND)
   integer, intent(INOUT) :: LST(IBG:IEND)

   integer :: L, LL, I
   real(8) :: ENO

   if(ENN(IBG)>0) then
      NCL=1
   else
      NCL=0
   end if

   do L=IBG+1,IEND
      ENO = ENN(LST(L))

      if(ENO>0) NCL=NCL+1

      LL=LST(L)
      I=L-1
      do while(I>IBG-1)
         if (ENN(LST(I))<=ENO) exit
         LST(I+1)=LST(I)
         I=I-1
      end do
      LST(I+1)=LL
   end do

   end subroutine SORTIT


   subroutine getirtau1(ib,nlevs,dp,fcld,reff,hydromets,&
                       tcldlyr,enn, &
                       aib_ir1, awb_ir1, &
                       aiw_ir1, aww_ir1, &
                       aig_ir1, awg_ir1, CONS_GRAV)

! !USES:

      implicit none

! !INPUT PARAMETERS:
      integer, intent(IN ) :: ib                 !  Band number
      integer, intent(IN ) :: nlevs              !  Number of levels
      real(8),    intent(IN ) :: dp(nlevs)          !  Delta pressure in Pa
      real(8),    intent(IN ) :: fcld(nlevs)        !  Cloud fraction (used sometimes)
      real(8),    intent(IN ) :: reff(nlevs,4)      !  Effective radius (microns)
      real(8),    intent(IN ) :: hydromets(nlevs,4) !  Hydrometeors (kg/kg)
      real(8),    intent(IN ) :: aib_ir1(3,10), awb_ir1(4,10), aiw_ir1(4,10)
      real(8),    intent(IN ) :: aww_ir1(4,10), aig_ir1(4,10), awg_ir1(4,10)
      real(8),    intent(IN ) :: CONS_GRAV

! !OUTPUT PARAMETERS:
      real(8),    intent(OUT) :: tcldlyr(0:nlevs  ) !  Flux transmissivity?
      real(8),    intent(OUT) ::     enn(0:nlevs  ) !  Flux transmissivity of a cloud layer?

! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for infrared wavelengths
!  Slots for reff, hydrometeors and tauall are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
!
!  In the below calculations, the constants used in the tau calculation are in 
!  m$^2$ g$^{-1}$ and m$^2$ g$^{-1}$ $\mu$m. Thus, we must convert the kg contained in the 
!  pressure (Pa = kg m$^{-1}$ s$^{-2}$) to grams.
!
! !REVISION HISTORY: 
!    2011.11.18   MAT moved to Radiation_Shared and revised arg list, units
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer            :: k
      real(8)               :: taucld1,taucld2,taucld3,taucld4
      real(8)               :: g1,g2,g3,g4,gg
      real(8)               :: w1,w2,w3,w4,ww
      real(8)               :: ff,tauc

      real(8)               :: reff_snow


!-----Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
!     Rain optical thickness is set to 0.00307 /(gm/m**2).
!     It is for a specific drop size distribution provided by Q. Fu.

      tcldlyr(0) = 1.0
      enn    (0) = 0.0

      do k = 1, nlevs
         if(reff(k,1)<=0.0) then
            taucld1=0.0
         else
            taucld1=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,1))*(aib_ir1(1,ib)+aib_ir1(2,ib)/&
                  reff(k,1)**aib_ir1(3,ib))
         end if

            taucld2=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,2))*(awb_ir1(1,ib)+(awb_ir1(2,ib)+&
                  (awb_ir1(3,ib)+awb_ir1(4,ib)*reff(k,2))*reff(k,2))*reff(k,2))

            taucld3=0.00307*(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,3))

!-----Below, we use the table of coefficients tabulated for suspended
!     cloud ice particles (aib_ir1) for falling snow. These coefficients
!     lead to unphysical (negative) values of cloud optical thickness 
!     for effective radii greater than 113 microns. By restricting the 
!     effective radius of snow to 112 microns, we prevent unphysical 
!     optical thicknesses.

         reff_snow = min(reff(k,4),112.0)

         if(reff_snow<=0.0) then
            taucld4=0.0
         else
            taucld4=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,4))*(aib_ir1(1,ib)+aib_ir1(2,ib)/&
                  reff_snow**aib_ir1(3,ib))
         end if

!-----Compute cloud single-scattering albedo and asymmetry factor for
!     a mixture of ice particles and liquid drops following 
!     Eqs. (6.5), (6.6), (6.15) and (6.16).
!     Single-scattering albedo and asymmetry factor of rain are set
!     to 0.54 and 0.95, respectively, based on the information provided
!     by Prof. Qiang Fu.

         tauc=taucld1+taucld2+taucld3+taucld4

         if (tauc > 0.02 .and. fcld(k) > 0.01) then

            w1=taucld1*(aiw_ir1(1,ib)+(aiw_ir1(2,ib)+(aiw_ir1(3,ib) &
                  +aiw_ir1(4,ib)*reff(k,1))*reff(k,1))*reff(k,1))
            w2=taucld2*(aww_ir1(1,ib)+(aww_ir1(2,ib)+(aww_ir1(3,ib)&
                  +aww_ir1(4,ib)*reff(k,2))*reff(k,2))*reff(k,2))
            w3=taucld3*0.54
            w4=taucld4*(aiw_ir1(1,ib)+(aiw_ir1(2,ib)+(aiw_ir1(3,ib) &
                  +aiw_ir1(4,ib)*reff_snow)*reff_snow)*reff_snow)
            ww=(w1+w2+w3+w4)/tauc


            g1=w1*(aig_ir1(1,ib)+(aig_ir1(2,ib)+(aig_ir1(3,ib)&
                  +aig_ir1(4,ib)*reff(k,1))*reff(k,1))*reff(k,1))
            g2=w2*(awg_ir1(1,ib)+(awg_ir1(2,ib)+(awg_ir1(3,ib) &
                  +awg_ir1(4,ib)*reff(k,2))*reff(k,2))*reff(k,2))
            g3=w3*0.95
            g4=w4*(aig_ir1(1,ib)+(aig_ir1(2,ib)+(aig_ir1(3,ib)&
                  +aig_ir1(4,ib)*reff_snow)*reff_snow)*reff_snow)
       
           if (abs(w1+w2+w3+w4).gt.0.0) then 
            gg=(g1+g2+g3+g4)/(w1+w2+w3+w4)
	    else 
	    gg=0.5
	   end if 

!-----Parameterization of LW scattering following Eqs. (6.11)
!     and (6.12). 

            ff=0.5+(0.3739+(0.0076+0.1185*gg)*gg)*gg

!ALT: temporary protection against negative cloud optical thickness

            tauc=max((1.-ww*ff),0.0)*tauc

!-----compute cloud diffuse transmittance. It is approximated by using 
!     a diffusivity factor of 1.66.

            tcldlyr(k) = exp(-1.66*tauc)
            enn    (k) = fcld(k)*(1.0-tcldlyr(k)) ! N in the documentation (6.13)

         else

            tcldlyr(k) = 1.0
            enn    (k) = 0.0

         end if
      end do


   end subroutine getirtau1


 subroutine sorad ( m               , & !Number of soundings
                    np              , & !Number of model levels
                    nb              , & !Number of bands
                    cosz_dev        , & !Cosine of solar zenith angle
                    pl_dev          , & !Pressure (Pa)
                    ta_dev          , & !Temperature (K)
                    wa_dev          , & !Specific humidity (kgkg^-1)
                    oa_dev          , & !Ozone (kgkg^-1)
                    co2             , & !CO2 (pppv)
                    cwc_dev         , & !Cloud properties (kgkg^-1)
                    fcld_dev        , & !Cloud fractions (1)
                    ict             , & !Level index separating high and middle clouds
                    icb             , & !Level index separating middle and low clouds
                    reff_dev        , & !Moisture refflectivity properties
                    hk_uv           , & !Solar UV constant
                    hk_ir           , & !Solar IR constant
                    taua_dev        , & !Aerosol optical thickness 
                    ssaa_dev        , & !Aerosol single scattering albedo
                    asya_dev        , & !Aerosol asymmetry factor
                    rsuvbm_dev      , & !Surface reflectivity in the UV+par for beam insolation
                    rsuvdf_dev      , & !Surface reflectivity in the UV+par for diffuse insolation
                    rsirbm_dev      , & !Surface reflectivity in the near-ir region for beam insolation
                    rsirdf_dev      , & !Surface reflectivity in the near-ir region for diffuse insolation
!Outputs
                    flx_dev         , & !
!Constants            
                    CONS_GRAV,                   &
                    wk_uv, zk_uv, ry_uv,         &
                    xk_ir, ry_ir,                &
                    cah, coa,                    &
                    aig_uv, awg_uv, arg_uv,      &
                    aib_uv, awb_uv, arb_uv,      &
                    aib_nir, awb_nir, arb_nir,   &
                    aia_nir, awa_nir, ara_nir,   &
                    aig_nir, awg_nir, arg_nir,   &
                    caib, caif                   ) 

   IMPLICIT NONE

   ! Parameters
   ! ----------

   integer, parameter :: nu = 43
   integer, parameter :: nw = 37
   integer, parameter :: nx = 62
   integer, parameter :: ny = 101
   integer, parameter :: nband_uv = 5
   integer, parameter :: nk_ir = 10
   integer, parameter :: nband_ir = 3

   integer, parameter :: nband = nband_uv + nband_ir

   real(8),    parameter :: dsm = 0.602


!-----input values

!-----input parameters

      integer m,np,ict,icb,nb
      real(8) cosz_dev(m),pl_dev(m,np+1),ta_dev(m,np),wa_dev(m,np),oa_dev(m,np),co2
      real(8) cwc_dev(m,np,4),fcld_dev(m,np),reff_dev(m,np,4), hk_uv(5),hk_ir(3,10)
      real(8) rsuvbm_dev(m),rsuvdf_dev(m),rsirbm_dev(m),rsirdf_dev(m)

      real(8) taua_dev(m,np,nb)
      real(8) ssaa_dev(m,np,nb)
      real(8) asya_dev(m,np,nb)

      logical :: overcast

! Constants
real(8), intent(in) :: wk_uv(5), zk_uv(5), ry_uv(5)
real(8), intent(in) :: xk_ir(10), ry_ir(3)
real(8), intent(in) :: cah(43,37), coa(62,101)

real(8), intent(in) :: aig_uv(3), awg_uv(3), arg_uv(3)
real(8), intent(in) :: aib_uv, awb_uv(2), arb_uv(2)
real(8), intent(in) :: aib_nir, awb_nir(3,2), arb_nir(3,2)
real(8), intent(in) :: aia_nir(3,3), awa_nir(3,3), ara_nir(3,3)
real(8), intent(in) :: aig_nir(3,3), awg_nir(3,3), arg_nir(3,3)
real(8), intent(in) :: caib(11,9,11), caif(9,11)

real(8), intent(in) :: CONS_GRAV

!-----output parameters

      real(8) flx_dev(m,np+1),flc_dev(m,np+1)
      real(8) flxu_dev(m,np+1),flcu_dev(m,np+1)
      real(8) fdiruv_dev (m),fdifuv_dev (m)
      real(8) fdirpar_dev(m),fdifpar_dev(m)
      real(8) fdirir_dev (m),fdifir_dev (m)
      real(8) flx_sfc_band_dev(m,nband)

!-----temporary arrays

      integer :: i,j,k,l,in,ntop

      real(8) :: dp(np),wh(np),oh(np)
      real(8) :: scal(np)
      real(8) :: swh(np+1),so2(np+1),df(0:np+1)
      real(8) :: scal0, wvtoa, o3toa, pa
      real(8) :: snt,cnt,x,xx4,xtoa
      real(8) :: dp_pa(np)

!-----parameters for co2 transmission tables

      real(8) :: w1,dw,u1,du

      integer :: ib,rc
      real(8) :: tauclb(np),tauclf(np),asycl(np)
      real(8) :: taubeam(np,4),taudiff(np,4)
      real(8) :: fcld_col(np)
      real(8) :: cwc_col(np,4)
      real(8) :: reff_col(np,4)
      real(8) :: taurs,tauoz,tauwv
      real(8) :: tausto,ssatau,asysto
      real(8) :: tautob,ssatob,asytob
      real(8) :: tautof,ssatof,asytof
      real(8) :: rr(0:np+1,2),tt(0:np+1,2),td(0:np+1,2)
      real(8) :: rs(0:np+1,2),ts(0:np+1,2)
      real(8) :: fall(np+1),fclr(np+1),fsdir,fsdif
      real(8) :: fupa(np+1),fupc(np+1)
      real(8) :: cc1,cc2,cc3
      real(8) :: rrt,ttt,tdt,rst,tst

      integer :: iv,ik
      real(8) :: ssacl(np)

      integer :: im

      integer :: ic,iw 
      real(8) :: ulog,wlog,dc,dd,x0,x1,x2,y0,y1,y2,du2,dw2

      integer :: ih

      !if (overcast == true) then
      !real(8) :: rra(0:np+1),rxa(0:np+1)
      !real(8) :: ttaold,tdaold,rsaold
      !real(8) :: ttanew,tdanew,rsanew 
      !else
      real(8) :: rra(0:np+1,2,2),tta(0:np,2,2)
      real(8) :: tda(0:np,2,2)
      real(8) :: rsa(0:np,2,2),rxa(0:np+1,2,2)
      !endif

      real(8) :: flxdn
      real(8) :: fdndir,fdndif,fupdif
      real(8) :: denm,yy

      integer :: is
      real(8) :: ch,cm,ct

      integer :: foundtop
      real(8) :: dftop

!-----Variables for aerosols

      integer :: II, JJ, irhp1, an

      real(8) :: dum

      i = 1
      !RUN_LOOP: do i=1,m

         ntop = 0
         fdndir = 0.0
         fdndif = 0.0

!-----Beginning of sorad code

!-----wvtoa and o3toa are the water vapor and o3 amounts of the region 
!     above the pl(1) level.
!     snt is the secant of the solar zenith angle

         snt    = 1.0/cosz_dev(i)
         xtoa   = max(pl_dev(i,1),1.e-3)
         scal0  = xtoa*(0.5*xtoa/300.)**.8
         o3toa  = 1.02*oa_dev(i,1)*xtoa*466.7 + 1.0e-8
         wvtoa  = 1.02*wa_dev(i,1)*scal0 * (1.0+0.00135*(ta_dev(i,1)-240.)) + 1.0e-9
         swh(1)  = wvtoa

         do k=1,np

!-----compute layer thickness. indices for the surface level and
!     surface layer are np+1 and np, respectively.

            dp(k) = pl_dev(i,k+1)-pl_dev(i,k)
            dp_pa(k) = dp(k) * 100. ! dp in pascals
 
!-----compute scaled water vapor amount following Eqs. (3.3) and (3.5) 
!     unit is g/cm**2
!
            pa   = 0.5*(pl_dev(i,k)+pl_dev(i,k+1))
            scal(k) = dp(k)*(pa/300.)**.8
            wh(k)   = 1.02*wa_dev(i,k)*scal(k) * (1.+0.00135*(ta_dev(i,k)-240.)) + 1.e-9
            swh(k+1)= swh(k)+wh(k)

!-----compute ozone amount, unit is (cm-atm)stp
!     the number 466.7 is the unit conversion factor
!     from g/cm**2 to (cm-atm)stp

            oh(k)   = 1.02*oa_dev(i,k)*dp(k)*466.7 + 1.e-8

!-----Fill the reff, cwc, and fcld for the column

            fcld_col(k) = fcld_dev(i,k)
            do l = 1, 4
               reff_col(k,l) = reff_dev(i,k,l)
               cwc_col(k,l) = cwc_dev(i,k,l)
            end do

         end do

!-----Initialize temporary arrays to zero to avoid UMR

         rr = 0.0
         tt = 0.0
         td = 0.0
         rs = 0.0
         ts = 0.0

         rra = 0.0
         rxa = 0.0

!if( OVERCAST == .false. ) then
         tta = 0.0
         tda = 0.0
         rsa = 0.0
!endif

!-----initialize fluxes for all-sky (flx), clear-sky (flc), and
!     flux reduction (df)
!
         do k=1,np+1
            flx_dev(i,k)=0.
            flc_dev(i,k)=0.
            flxu_dev(i,k)=0.
            flcu_dev(i,k)=0.
         end do

!-----Initialize new per-band surface fluxes

         do ib = 1, nband
            flx_sfc_band_dev(i,ib) = 0.
         end do

!-----Begin inline of SOLUV

!-----compute solar uv and par fluxes

!-----initialize fdiruv, fdifuv, surface reflectances and transmittances.
!     the reflectance and transmittance of the clear and cloudy portions
!     of a layer are denoted by 1 and 2, respectively.
!     cc is the maximum cloud cover in each of the high, middle, and low
!     cloud groups.
!     1/dsm=1/cos(53) = 1.66

         fdiruv_dev(i)=0.0
         fdifuv_dev(i)=0.0
         rr(np+1,1)=rsuvbm_dev(i)
         rr(np+1,2)=rsuvbm_dev(i)
         rs(np+1,1)=rsuvdf_dev(i)
         rs(np+1,2)=rsuvdf_dev(i)
         td(np+1,1)=0.0
         td(np+1,2)=0.0
         tt(np+1,1)=0.0
         tt(np+1,2)=0.0
         ts(np+1,1)=0.0
         ts(np+1,2)=0.0
         rr(0,1)=0.0
         rr(0,2)=0.0
         rs(0,1)=0.0
         rs(0,2)=0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
         tt(0,1)=1.0
         tt(0,2)=1.0
         ts(0,1)=1.0
         ts(0,2)=1.0
         cc1=0.0
         cc2=0.0
         cc3=0.0

!-----options for scaling cloud optical thickness

!if ( OVERCAST == .true. ) then

!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud asymmetry factor
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

!         call getvistau1(np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,0,0,&
!                        taubeam,taudiff,asycl,                                  &
!                         aig_uv, awg_uv, arg_uv,                                 &
!                         aib_uv, awb_uv, arb_uv,                                 &
!                         aib_nir, awb_nir, arb_nir,                              &
!                         aia_nir, awa_nir, ara_nir,                              &
!                         aig_nir, awg_nir, arg_nir,                              &
!                         caib, caif,                                             &
!                         CONS_GRAV                                               )

!else

!-----Compute scaled cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud asymmetry factor
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

         call getvistau1(np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,ict,icb, &
                         taubeam,taudiff,asycl,                                  &
                         aig_uv, awg_uv, arg_uv,                                 &
                         aib_uv, awb_uv, arb_uv,                                 &
                         aib_nir, awb_nir, arb_nir,                              &
                         aia_nir, awa_nir, ara_nir,                              &
                         aig_nir, awg_nir, arg_nir,                              &
                         caib, caif,                                             &
                         CONS_GRAV                                               )

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
!     The cc1,2,3 are still needed in the flux calculations below

!MAT---DO NOT FUSE THIS LOOP
!MAT---Loop must run to completion so that cc[1,2,3] are correct.
         do k=1,np
            if (k.lt.ict) then
               cc1=max(cc1,fcld_dev(i,k))
            elseif (k.lt.icb) then
               cc2=max(cc2,fcld_dev(i,k))
            else
               cc3=max(cc3,fcld_dev(i,k))
            end if
         end do
!MAT---DO NOT FUSE THIS LOOP

!endif !overcast

         do k=1,np
            tauclb(k)=taubeam(k,1)+taubeam(k,2)+taubeam(k,3)+taubeam(k,4)
            tauclf(k)=taudiff(k,1)+taudiff(k,2)+taudiff(k,3)+taudiff(k,4)
         end do

!-----integration over spectral bands

!-----Compute optical thickness, single-scattering albedo and asymmetry
!     factor for a mixture of "na" aerosol types. [Eqs. (4.16)-(4.18)]

         do ib=1,nband_uv

!-----compute direct beam transmittances of the layer above pl(1)

            td(0,1)=exp(-(wvtoa*wk_uv(ib)+o3toa*zk_uv(ib))/cosz_dev(i))
            td(0,2)=td(0,1)

            do k=1,np

!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor (Eqs. 6.2-6.4)

               taurs=ry_uv(ib)*dp(k)
               tauoz=zk_uv(ib)*oh(k)
               tauwv=wk_uv(ib)*wh(k)

               tausto=taurs+tauoz+tauwv+taua_dev(i,k,ib)+1.0e-7
               ssatau=ssaa_dev(i,k,ib)+taurs
               asysto=asya_dev(i,k,ib)

               tautob=tausto
               asytob=asysto/ssatau
               ssatob=ssatau/tautob+1.0e-8
               ssatob=min(ssatob,0.999999)

!-----for direct incident radiation

               call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

               call deledd(tautob,ssatob,asytob,dsm,rst,tst,dum)

               rr(k,1)=rrt
               tt(k,1)=ttt
               td(k,1)=tdt
               rs(k,1)=rst
               ts(k,1)=tst

!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer

!-----for direct incident radiation
!     The effective layer optical properties. Eqs. (6.2)-(6.4)

               tautob=tausto+tauclb(k)
               ssatob=(ssatau+tauclb(k))/tautob+1.0e-8
               ssatob=min(ssatob,0.999999)
               asytob=(asysto+asycl(k)*tauclb(k))/(ssatob*tautob)

!-----for diffuse incident radiation

               tautof=tausto+tauclf(k)
               ssatof=(ssatau+tauclf(k))/tautof+1.0e-8
               ssatof=min(ssatof,0.999999)
               asytof=(asysto+asycl(k)*tauclf(k))/(ssatof*tautof)

!-----for direct incident radiation
!     note that the cloud optical thickness is scaled differently 
!     for direct and diffuse insolation, Eqs. (7.3) and (7.4).

               call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation 
!     with an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

               call deledd(tautof,ssatof,asytof,dsm,rst,tst,dum)

               rr(k,2)=rrt
               tt(k,2)=ttt
               td(k,2)=tdt
               rs(k,2)=rst
               ts(k,2)=tst
            end do

!-----flux calculations
!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)

            do k=1,np+1
               fclr(k)=0.0
               fall(k)=0.0
               fupa(k)=0.0
               fupc(k)=0.0
            end do

            fsdir=0.0
            fsdif=0.0

!if ( OVERCAST == .true. ) then

!-----Inline CLDFLXY

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is either 0 or 1.

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated by 
!         beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)

!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

!-----ih=1 for clear sky; ih=2 for cloudy sky.

!-----First set is ih = 1
!            rra(np+1)=rr(np+1,1)
!            rxa(np+1)=rs(np+1,1)
!
!            do k=np,0,-1
!               denm=ts(k,1)/(1.-rs(k,1)*rxa(k+1))
!               rra(k)=rr(k,1)+(td(k,1)*rra(k+1)+(tt(k,1)-td(k,1))*rxa(k+1))*denm
!               rxa(k)=rs(k,1)+ts(k,1)*rxa(k+1)*denm
!            end do
!
!            do k=1,np+1
!               if (k <= np) then
!                  if (k == 1) then
!                     tdaold = td(0,1)
!                     ttaold = tt(0,1)
!                     rsaold = rs(0,1)
!
!                     tdanew = 0.0
!                     ttanew = 0.0
!                     rsanew = 0.0
!                  end if
!                  denm=ts(k,1)/(1.-rsaold*rs(k,1))
!                  tdanew=tdaold*td(k,1)
!                  ttanew=tdaold*tt(k,1)+(tdaold*rsaold*rr(k,1)+ttaold-tdaold)*denm
!                  rsanew=rs(k,1)+ts(k,1)*rsaold*denm
!               end if
!
!               denm=1./(1.-rsaold*rxa(k))
!               fdndir=tdaold
!               xx4=tdaold*rra(k)
!               yy=ttaold-tdaold
!               fdndif=(xx4*rsaold+yy)*denm
!               fupdif=(xx4+yy*rxa(k))*denm
!               flxdn=fdndir+fdndif-fupdif
!               fupc(k)=fupdif 
!               fclr(k)=flxdn
!
!               tdaold = tdanew
!               ttaold = ttanew
!               rsaold = rsanew
!
!               tdanew = 0.0
!               ttanew = 0.0
!               rsanew = 0.0
!            end do
!
!!-----Second set is ih = 2
!
!            rra(np+1)=rr(np+1,2)
!            rxa(np+1)=rs(np+1,2)
!
!            do k=np,0,-1
!               denm=ts(k,2)/(1.-rs(k,2)*rxa(k+1))
!               rra(k)=rr(k,2)+(td(k,2)*rra(k+1)+(tt(k,2)-td(k,2))*rxa(k+1))*denm
!               rxa(k)=rs(k,2)+ts(k,2)*rxa(k+1)*denm
!            end do
!
!            do k=1,np+1
!               if (k <= np) then
!                  if (k == 1) then
!                     tdaold = td(0,2)
!                     ttaold = tt(0,2)
!                     rsaold = rs(0,2)

!                     tdanew = 0.0
!                     ttanew = 0.0
!                     rsanew = 0.0
!                  end if
!                  denm=ts(k,2)/(1.-rsaold*rs(k,2))
!                  tdanew=tdaold*td(k,2)
!                  ttanew=tdaold*tt(k,2)+(tdaold*rsaold*rr(k,2)+ttaold-tdaold)*denm
!                  rsanew=rs(k,2)+ts(k,2)*rsaold*denm
!               end if
!
!               denm=1./(1.-rsaold*rxa(k))
!               fdndir=tdaold
!               xx4=tdaold*rra(k)
!               yy=ttaold-tdaold
!               fdndif=(xx4*rsaold+yy)*denm
!               fupdif=(xx4+yy*rxa(k))*denm
!               flxdn=fdndir+fdndif-fupdif
!
!               fupa(k)=fupdif
!               fall(k)=flxdn
!
!               tdaold = tdanew
!               ttaold = ttanew
!               rsaold = rsanew
!
!               tdanew = 0.0
!               ttanew = 0.0
!               rsanew = 0.0
!            end do
!
!            fsdir=fdndir
!            fsdif=fdndif
!
!!-----End CLDFLXY inline
!
!else

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is allowed to be between 0 and 1.
!     the all-sky flux, fall is the summation inside the brackets
!     of Eq. (7.11)

!-----Inline CLDFLX

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated 
!         by beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----To save memory space, tda, tta, and rsa are pre-computed 
!     for k<icb. The dimension of these parameters is (m,np,2,2). 
!     It would have been (m,np,2,2,2) if these parameters were 
!     computed for all k's.

!-----for high clouds
!     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition

            do ih=1,2
               tda(0,ih,1)=td(0,ih)
               tta(0,ih,1)=tt(0,ih)
               rsa(0,ih,1)=rs(0,ih)
               tda(0,ih,2)=td(0,ih)
               tta(0,ih,2)=tt(0,ih)
               rsa(0,ih,2)=rs(0,ih)

               do k=1,ict-1
                  denm=ts(k,ih)/(1.-rsa(k-1,ih,1)*rs(k,ih))
                  tda(k,ih,1)=tda(k-1,ih,1)*td(k,ih)
                  tta(k,ih,1)=tda(k-1,ih,1)*tt(k,ih)+(tda(k-1,ih,1)*rsa(k-1,ih,1)&
                        *rr(k,ih)+tta(k-1,ih,1)-tda(k-1,ih,1))*denm
                  rsa(k,ih,1)=rs(k,ih)+ts(k,ih)*rsa(k-1,ih,1)*denm
                  tda(k,ih,2)=tda(k,ih,1)
                  tta(k,ih,2)=tta(k,ih,1)
                  rsa(k,ih,2)=rsa(k,ih,1)
               end do ! k loop

!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition

               do k=ict,icb-1
                  do im=1,2
                     denm=ts(k,im)/(1.-rsa(k-1,ih,im)*rs(k,im))
                     tda(k,ih,im)=tda(k-1,ih,im)*td(k,im)
                     tta(k,ih,im)=tda(k-1,ih,im)*tt(k,im)+(tda(k-1,ih,im)&
                           *rsa(k-1,ih,im)*rr(k,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                     rsa(k,ih,im)=rs(k,im)+ts(k,im)*rsa(k-1,ih,im)*denm
                  end do ! im loop
               end do ! k loop
            end do ! ih loop

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----To save memory space, rra and rxa are pre-computed for k>=icb.
!     the dimension of these parameters is (m,np,2,2). It would have
!     been (m,np,2,2,2) if these parameters were computed for all k's.

!-----for the low clouds
!     is=1 for clear-sky condition, is=2 for cloudy-sky condition

            do is=1,2
               rra(np+1,1,is)=rr(np+1,is)
               rxa(np+1,1,is)=rs(np+1,is)
               rra(np+1,2,is)=rr(np+1,is)
               rxa(np+1,2,is)=rs(np+1,is)

               do k=np,icb,-1
                  denm=ts(k,is)/(1.-rs(k,is)*rxa(k+1,1,is))
                  rra(k,1,is)=rr(k,is)+(td(k,is)*rra(k+1,1,is)+(tt(k,is)-td(k,is))&
                        *rxa(k+1,1,is))*denm
                  rxa(k,1,is)=rs(k,is)+ts(k,is)*rxa(k+1,1,is)*denm
                  rra(k,2,is)=rra(k,1,is)
                  rxa(k,2,is)=rxa(k,1,is)
               end do ! k loop

!-----for middle clouds

               do k=icb-1,ict,-1
                  do im=1,2
                     denm=ts(k,im)/(1.-rs(k,im)*rxa(k+1,im,is))
                     rra(k,im,is)=rr(k,im)+(td(k,im)*rra(k+1,im,is)+(tt(k,im)-td(k,im))&
                           *rxa(k+1,im,is))*denm
                     rxa(k,im,is)=rs(k,im)+ts(k,im)*rxa(k+1,im,is)*denm
                  end do ! im loop
               end do ! k loop
            end do ! is loop

!-----integration over eight sky situations.
!     ih, im, is denote high, middle and low cloud groups.

            do ih=1,2

!-----clear portion 
               if(ih.eq.1) then
                  ch=1.0-cc1
!-----cloudy portion
               else
                  ch=cc1
               end if

               do im=1,2
!-----clear portion 
                  if(im.eq.1) then
                     cm=ch*(1.0-cc2)
!-----cloudy portion
                  else
                     cm=ch*cc2 
                  end if

                  do is=1,2
!-----clear portion 
                     if(is.eq.1) then
                        ct=cm*(1.0-cc3) 
!-----cloudy portion
                     else
                        ct=cm*cc3
                     end if

!-----add one layer at a time, going down.

                     do k=icb,np
                        denm=ts(k,is)/(1.-rsa(k-1,ih,im)*rs(k,is))
                        tda(k,ih,im)=tda(k-1,ih,im)*td(k,is)
                        tta(k,ih,im)=tda(k-1,ih,im)*tt(k,is)+(tda(k-1,ih,im)*rr(k,is)&
                              *rsa(k-1,ih,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                        rsa(k,ih,im)=rs(k,is)+ts(k,is)*rsa(k-1,ih,im)*denm
                     end do ! k loop

!-----add one layer at a time, going up.

                     do k=ict-1,0,-1
                        denm=ts(k,ih)/(1.-rs(k,ih)*rxa(k+1,im,is))
                        rra(k,im,is)=rr(k,ih)+(td(k,ih)*rra(k+1,im,is)+(tt(k,ih)-td(k,ih))&
                              *rxa(k+1,im,is))*denm
                        rxa(k,im,is)=rs(k,ih)+ts(k,ih)*rxa(k+1,im,is)*denm
                     end do ! k loop

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)

!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

                     do k=1,np+1
                        denm=1./(1.-rsa(k-1,ih,im)*rxa(k,im,is))
                        fdndir=tda(k-1,ih,im)
                        xx4=tda(k-1,ih,im)*rra(k,im,is)
                        yy=tta(k-1,ih,im)-tda(k-1,ih,im)
                        fdndif=(xx4*rsa(k-1,ih,im)+yy)*denm
                        fupdif=(xx4+yy*rxa(k,im,is))*denm
                        flxdn=fdndir+fdndif-fupdif

!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)

                        if(ih.eq.1 .and. im.eq.1 .and. is.eq.1) then
                           fupc(k)=fupdif
                           fclr(k)=flxdn
                        end if
                        fupa(k)=fupa(k)+fupdif*ct
                        fall(k)=fall(k)+flxdn*ct
                     end do ! k loop
                     fsdir=fsdir+fdndir*ct
                     fsdif=fsdif+fdndif*ct
                  end do ! is loop
               end do ! im loop
            end do ! ih loop

!-----End CLDFLX inline

!endif !overcast

!-----flux integration, Eq. (6.1)

            do k=1,np+1
               flx_dev(i,k)=flx_dev(i,k)+fall(k)*hk_uv(ib)
               flc_dev(i,k)=flc_dev(i,k)+fclr(k)*hk_uv(ib)
               flxu_dev(i,k)=flxu_dev(i,k)+fupa(k)*hk_uv(ib)
               flcu_dev(i,k)=flcu_dev(i,k)+fupc(k)*hk_uv(ib)
            end do

!-----get surface flux for each band
            flx_sfc_band_dev(i,ib)=flx_sfc_band_dev(i,ib)+fall(np+1)*hk_uv(ib)

!-----compute direct and diffuse downward surface fluxes in the UV
!     and par regions

            if(ib.lt.5) then
               fdiruv_dev(i)=fdiruv_dev(i)+fsdir*hk_uv(ib)
               fdifuv_dev(i)=fdifuv_dev(i)+fsdif*hk_uv(ib)
            else
               fdirpar_dev(i)=fsdir*hk_uv(ib)
               fdifpar_dev(i)=fsdif*hk_uv(ib)
            end if
         end do

!-----Inline SOLIR

!-----compute and update solar ir fluxes

         fdirir_dev(i)=0.0
         fdifir_dev(i)=0.0
         rr(np+1,1)=rsirbm_dev(i)
         rr(np+1,2)=rsirbm_dev(i)
         rs(np+1,1)=rsirdf_dev(i)
         rs(np+1,2)=rsirdf_dev(i)
         td(np+1,1)=0.0
         td(np+1,2)=0.0
         tt(np+1,1)=0.0
         tt(np+1,2)=0.0
         ts(np+1,1)=0.0
         ts(np+1,2)=0.0
         rr(0,1)=0.0
         rr(0,2)=0.0
         rs(0,1)=0.0
         rs(0,2)=0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
         tt(0,1)=1.0
         tt(0,2)=1.0
         ts(0,1)=1.0
         ts(0,2)=1.0
         cc1=0.0
         cc2=0.0
         cc3=0.0

!-----integration over spectral bands

!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10)
!     The indices 1, 2, 3 are for ice, water, rain particles,
!     respectively.

         do ib=1,nband_ir
            iv=ib+5

!-----options for scaling cloud optical thickness

!if ( OVERCAST == .true. ) then

!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

!            call getnirtau1(ib,np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,0,0,&
!                           taubeam,taudiff,asycl,ssacl,                               &
!                            aig_uv, awg_uv, arg_uv,                                    &
!                            aib_uv, awb_uv, arb_uv,                                    &
!                            aib_nir, awb_nir, arb_nir,                                 &
!                            aia_nir, awa_nir, ara_nir,                                 &
!                            aig_nir, awg_nir, arg_nir,                                 &
!                            caib, caif,                                                &
!                            CONS_GRAV                                                  )

!else

!-----Compute scaled cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

            call getnirtau1(ib,np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,ict,icb, &
                            taubeam,taudiff,asycl,ssacl,                               &
                            aig_uv, awg_uv, arg_uv,                                    &
                            aib_uv, awb_uv, arb_uv,                                    &
                            aib_nir, awb_nir, arb_nir,                                 &
                            aia_nir, awa_nir, ara_nir,                                 &
                            aig_nir, awg_nir, arg_nir,                                 &
                            caib, caif,                                                &
                            CONS_GRAV                                                  )

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group

!MAT--DO NOT FUSE THIS LOOP
!MAT  Loop must run to completion so that cc[1,2,3] are correct.
            do k=1,np
               if (k.lt.ict) then
                  cc1=max(cc1,fcld_dev(i,k))
               elseif (k.lt.icb) then
                  cc2=max(cc2,fcld_dev(i,k))
               else
                  cc3=max(cc3,fcld_dev(i,k))
               end if
            end do
!MAT--DO NOT FUSE THIS LOOP

!endif !overcast

            do k=1,np
               tauclb(k)=taubeam(k,1)+taubeam(k,2)+taubeam(k,3)+taubeam(k,4)
               tauclf(k)=taudiff(k,1)+taudiff(k,2)+taudiff(k,3)+taudiff(k,4)
            end do


!-----integration over the k-distribution function

            do ik=1,nk_ir

!-----compute direct beam transmittances of the layer above pl(1)

               td(0,1)=exp(-wvtoa*xk_ir(ik)/cosz_dev(i))
               td(0,2)=td(0,1)

               do k=1,np
                  taurs=ry_ir(ib)*dp(k)
                  tauwv=xk_ir(ik)*wh(k)

!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor. Eqs.(6.2)-(6.4)

                  tausto=taurs+tauwv+taua_dev(i,k,iv)+1.0e-7
                  ssatau=ssaa_dev(i,k,iv)+taurs+1.0e-8
                  asysto=asya_dev(i,k,iv)
                  tautob=tausto
                  asytob=asysto/ssatau
                  ssatob=ssatau/tautob+1.0e-8
                  ssatob=min(ssatob,0.999999)

!-----Compute reflectance and transmittance of the clear portion 
!     of a layer

!-----for direct incident radiation

                  call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

                  call deledd(tautob,ssatob,asytob,dsm,rst,tst,dum)

                  rr(k,1)=rrt
                  tt(k,1)=ttt
                  td(k,1)=tdt
                  rs(k,1)=rst
                  ts(k,1)=tst

!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer

!-----for direct incident radiation. Eqs.(6.2)-(6.4)

                  tautob=tausto+tauclb(k)
                  ssatob=(ssatau+ssacl(k)*tauclb(k))/tautob+1.0e-8
                  ssatob=min(ssatob,0.999999)
                  asytob=(asysto+asycl(k)*ssacl(k)*tauclb(k))/(ssatob*tautob)

!-----for diffuse incident radiation

                  tautof=tausto+tauclf(k)
                  ssatof=(ssatau+ssacl(k)*tauclf(k))/tautof+1.0e-8
                  ssatof=min(ssatof,0.999999)
                  asytof=(asysto+asycl(k)*ssacl(k)*tauclf(k))/(ssatof*tautof)

!-----for direct incident radiation

                  call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs.(6.5) and (6.6)

                  call deledd(tautof,ssatof,asytof,dsm,rst,tst,dum)

                  rr(k,2)=rrt
                  tt(k,2)=ttt
                  td(k,2)=tdt
                  rs(k,2)=rst
                  ts(k,2)=tst
               end do

!-----FLUX CALCULATIONS

!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)

               do k=1,np+1
                  fclr(k)=0.0
                  fall(k)=0.0
                  fupc(k)=0.0
                  fupa(k)=0.0
               end do

               fsdir=0.0
               fsdif=0.0

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is either 0 or 1.

!if ( OVERCAST == .true. ) then

!-----Inline CLDFLXY

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated by 
!         beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)

!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

!-----ih=1 for clear sky; ih=2 for cloudy sky.

!-----First set is ih = 1
!               rra(np+1)=rr(np+1,1)
!               rxa(np+1)=rs(np+1,1)
!
!               do k=np,0,-1
!                  denm=ts(k,1)/(1.-rs(k,1)*rxa(k+1))
!                  rra(k)=rr(k,1)+(td(k,1)*rra(k+1)+(tt(k,1)-td(k,1))*rxa(k+1))*denm
!                  rxa(k)=rs(k,1)+ts(k,1)*rxa(k+1)*denm
!               end do
!
!               do k=1,np+1
!                  if (k <= np) then
!                     if (k == 1) then
!                        tdaold = td(0,1)
!                        ttaold = tt(0,1)
!                        rsaold = rs(0,1)
!
!                        tdanew = 0.0
!                        ttanew = 0.0
!                        rsanew = 0.0
!                     end if
!                     denm=ts(k,1)/(1.-rsaold*rs(k,1))
!                     tdanew=tdaold*td(k,1)
!                     ttanew=tdaold*tt(k,1)+(tdaold*rsaold*rr(k,1)+ttaold-tdaold)*denm
!                     rsanew=rs(k,1)+ts(k,1)*rsaold*denm
!                  end if
!
!                  denm=1./(1.-rsaold*rxa(k))
!                  fdndir=tdaold
!                  xx4=tdaold*rra(k)
!                  yy=ttaold-tdaold
!                  fdndif=(xx4*rsaold+yy)*denm
!                  fupdif=(xx4+yy*rxa(k))*denm
!                  flxdn=fdndir+fdndif-fupdif
!
!                  fupc(k)=fupdif
!                  fclr(k)=flxdn
!
!                  tdaold = tdanew
!                  ttaold = ttanew
!                  rsaold = rsanew
!
!                  tdanew = 0.0
!                  ttanew = 0.0
!                  rsanew = 0.0
!               end do
!
!!-----Second set is ih = 2
!
!               rra(np+1)=rr(np+1,2)
!               rxa(np+1)=rs(np+1,2)
!
!               do k=np,0,-1
!                  denm=ts(k,2)/(1.-rs(k,2)*rxa(k+1))
!                  rra(k)=rr(k,2)+(td(k,2)*rra(k+1)+(tt(k,2)-td(k,2))*rxa(k+1))*denm
!                  rxa(k)=rs(k,2)+ts(k,2)*rxa(k+1)*denm
!               end do
!
!               do k=1,np+1
!                  if (k <= np) then
!                     if (k == 1) then
!                        tdaold = td(0,2)
!                        ttaold = tt(0,2)
!                        rsaold = rs(0,2)
!
!                        tdanew = 0.0
!                        ttanew = 0.0
!                        rsanew = 0.0
!                     end if
!                     denm=ts(k,2)/(1.-rsaold*rs(k,2))
!                     tdanew=tdaold*td(k,2)
!                     ttanew=tdaold*tt(k,2)+(tdaold*rsaold*rr(k,2)+ttaold-tdaold)*denm
!                     rsanew=rs(k,2)+ts(k,2)*rsaold*denm
!                  end if
!
!                  denm=1./(1.-rsaold*rxa(k))
!                  fdndir=tdaold
!                  xx4=tdaold*rra(k)
!                  yy=ttaold-tdaold
!                  fdndif=(xx4*rsaold+yy)*denm
!                  fupdif=(xx4+yy*rxa(k))*denm
!                  flxdn=fdndir+fdndif-fupdif
!
!                  fupa(k)=fupdif
!                  fall(k)=flxdn
!
!                  tdaold = tdanew
!                  ttaold = ttanew
!                  rsaold = rsanew
!
!                  tdanew = 0.0
!                  ttanew = 0.0
!                  rsanew = 0.0
!               end do
!
!               fsdir=fdndir
!               fsdif=fdndif
!
!!-----End CLDFLXY inline
!
!else

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is allowed to be between 0 and 1.
!     the all-sky flux, fall is the summation inside the brackets
!     of Eq. (7.11)

!-----Inline CLDFLX

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated 
!         by beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----To save memory space, tda, tta, and rsa are pre-computed 
!     for k<icb. The dimension of these parameters is (m,np,2,2). 
!     It would have been (m,np,2,2,2) if these parameters were 
!     computed for all k's.

!-----for high clouds
!     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition

               do ih=1,2
                  tda(0,ih,1)=td(0,ih)
                  tta(0,ih,1)=tt(0,ih)
                  rsa(0,ih,1)=rs(0,ih)
                  tda(0,ih,2)=td(0,ih)
                  tta(0,ih,2)=tt(0,ih)
                  rsa(0,ih,2)=rs(0,ih)

                  do k=1,ict-1
                     denm=ts(k,ih)/(1.-rsa(k-1,ih,1)*rs(k,ih))
                     tda(k,ih,1)=tda(k-1,ih,1)*td(k,ih)
                     tta(k,ih,1)=tda(k-1,ih,1)*tt(k,ih)+(tda(k-1,ih,1)*rsa(k-1,ih,1)&
                           *rr(k,ih)+tta(k-1,ih,1)-tda(k-1,ih,1))*denm
                     rsa(k,ih,1)=rs(k,ih)+ts(k,ih)*rsa(k-1,ih,1)*denm
                     tda(k,ih,2)=tda(k,ih,1)
                     tta(k,ih,2)=tta(k,ih,1)
                     rsa(k,ih,2)=rsa(k,ih,1)
                  end do ! k loop

!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition

                  do k=ict,icb-1
                     do im=1,2
                        denm=ts(k,im)/(1.-rsa(k-1,ih,im)*rs(k,im))
                        tda(k,ih,im)=tda(k-1,ih,im)*td(k,im)
                        tta(k,ih,im)=tda(k-1,ih,im)*tt(k,im)+(tda(k-1,ih,im)&
                              *rsa(k-1,ih,im)*rr(k,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                        rsa(k,ih,im)=rs(k,im)+ts(k,im)*rsa(k-1,ih,im)*denm
                     end do ! im loop
                  end do ! k loop
               end do ! ih loop

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----To save memory space, rra and rxa are pre-computed for k>=icb.
!     the dimension of these parameters is (m,np,2,2). It would have
!     been (m,np,2,2,2) if these parameters were computed for all k's.

!-----for the low clouds
!     is=1 for clear-sky condition, is=2 for cloudy-sky condition

               do is=1,2
                  rra(np+1,1,is)=rr(np+1,is)
                  rxa(np+1,1,is)=rs(np+1,is)
                  rra(np+1,2,is)=rr(np+1,is)
                  rxa(np+1,2,is)=rs(np+1,is)

                  do k=np,icb,-1
                     denm=ts(k,is)/(1.-rs(k,is)*rxa(k+1,1,is))
                     rra(k,1,is)=rr(k,is)+(td(k,is)*rra(k+1,1,is)+(tt(k,is)-td(k,is))&
                           *rxa(k+1,1,is))*denm
                     rxa(k,1,is)=rs(k,is)+ts(k,is)*rxa(k+1,1,is)*denm
                     rra(k,2,is)=rra(k,1,is)
                     rxa(k,2,is)=rxa(k,1,is)
                  end do ! k loop

!-----for middle clouds

                  do k=icb-1,ict,-1
                     do im=1,2
                        denm=ts(k,im)/(1.-rs(k,im)*rxa(k+1,im,is))
                        rra(k,im,is)=rr(k,im)+(td(k,im)*rra(k+1,im,is)+(tt(k,im)-td(k,im))&
                              *rxa(k+1,im,is))*denm
                        rxa(k,im,is)=rs(k,im)+ts(k,im)*rxa(k+1,im,is)*denm
                     end do ! im loop
                  end do ! k loop
               end do ! is loop

!-----integration over eight sky situations.
!     ih, im, is denote high, middle and low cloud groups.

               do ih=1,2
!-----clear portion 
                  if(ih.eq.1) then
                     ch=1.0-cc1
!-----cloudy portion
                  else
                     ch=cc1
                  end if

                  do im=1,2
!-----clear portion 
                     if(im.eq.1) then
                        cm=ch*(1.0-cc2)
!-----cloudy portion
                     else
                        cm=ch*cc2 
                     end if

                     do is=1,2
!-----clear portion 
                        if(is.eq.1) then
                           ct=cm*(1.0-cc3) 
!-----cloudy portion
                        else
                           ct=cm*cc3
                        end if

!-----add one layer at a time, going down.

                        do k=icb,np
                           denm=ts(k,is)/(1.-rsa(k-1,ih,im)*rs(k,is))
                           tda(k,ih,im)=tda(k-1,ih,im)*td(k,is)
                           tta(k,ih,im)=tda(k-1,ih,im)*tt(k,is)+(tda(k-1,ih,im)*rr(k,is)&
                                 *rsa(k-1,ih,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                           rsa(k,ih,im)=rs(k,is)+ts(k,is)*rsa(k-1,ih,im)*denm
                        end do ! k loop

!-----add one layer at a time, going up.

                        do k=ict-1,0,-1
                           denm=ts(k,ih)/(1.-rs(k,ih)*rxa(k+1,im,is))
                           rra(k,im,is)=rr(k,ih)+(td(k,ih)*rra(k+1,im,is)+(tt(k,ih)-td(k,ih))&
                                 *rxa(k+1,im,is))*denm
                           rxa(k,im,is)=rs(k,ih)+ts(k,ih)*rxa(k+1,im,is)*denm
                        end do ! k loop

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

                        do k=1,np+1
                           denm=1./(1.-rsa(k-1,ih,im)*rxa(k,im,is))
                           fdndir=tda(k-1,ih,im)
                           xx4=tda(k-1,ih,im)*rra(k,im,is)
                           yy=tta(k-1,ih,im)-tda(k-1,ih,im)
                           fdndif=(xx4*rsa(k-1,ih,im)+yy)*denm
                           fupdif=(xx4+yy*rxa(k,im,is))*denm
                           flxdn=fdndir+fdndif-fupdif

!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)

                           if(ih.eq.1 .and. im.eq.1 .and. is.eq.1) then
                              fupc(k)=fupdif
                              fclr(k)=flxdn
                           end if
                           fupa(k)=fupa(k)+fupdif*ct
                           fall(k)=fall(k)+flxdn*ct
                        end do ! k loop
                        fsdir=fsdir+fdndir*ct
                        fsdif=fsdif+fdndif*ct
                     end do ! is loop
                  end do ! im loop
               end do ! ih loop

!-----End CLDFLX inline
!endif !overcast

!-----flux integration following Eq. (6.1)

               do k=1,np+1
                  flx_dev(i,k)=flx_dev(i,k)+fall(k)*hk_ir(ib,ik)
                  flc_dev(i,k)=flc_dev(i,k)+fclr(k)*hk_ir(ib,ik)
                  flxu_dev(i,k)=flxu_dev(i,k)+fupa(k)*hk_ir(ib,ik)
                  flcu_dev(i,k)=flcu_dev(i,k)+fupc(k)*hk_ir(ib,ik)
               end do

!-----compute downward surface fluxes in the ir region

               fdirir_dev(i)=fdirir_dev(i)+fsdir*hk_ir(ib,ik)
               fdifir_dev(i)=fdifir_dev(i)+fsdif*hk_ir(ib,ik)

!-----tabulate surface flux at ir bands
               flx_sfc_band_dev(i,iv)=flx_sfc_band_dev(i,iv)+fall(np+1)*hk_ir(ib,ik)

            end do ! ik loop
         end do

!-----compute pressure-scaled o2 amount following Eq. (3.5) with f=1.
!     unit is (cm-atm)stp. 165.22 = (1000/980)*23.14%*(22400/32)
!     compute flux reduction due to oxygen following Eq. (3.18). 0.0633 is the
!     fraction of insolation contained in the oxygen bands

         df(0) = 0.0
         cnt = 165.22*snt
         so2(1) = scal0*cnt
! LLT increased parameter 145 to 155 to enhance effect
         df(1) = 0.0633*(1.-exp(-0.000155*sqrt(so2(1))))

         do k=1,np
            so2(k+1) = so2(k) + scal(k)*cnt
! LLT increased parameter 145 to 155 to enhance effect
            df(k+1) = 0.0633*(1.0 - exp(-0.000155*sqrt(so2(k+1)))) 
         end do

!-----for solar heating due to co2 scaling follows Eq(3.5) with f=1.
!     unit is (cm-atm)stp. 789 = (1000/980)*(44/28.97)*(22400/44)

         so2(1) = (789.*co2)*scal0

         do k=1,np
            so2(k+1) = so2(k) + (789.*co2)*scal(k)
         end do

!-----The updated flux reduction for co2 absorption in Band 7 where absorption due to
!     water vapor and co2 are both moderate. df is given by the second term on the
!     right-hand-side of Eq. (3.24) divided by So. so2 and swh are the co2 and
!     water vapor amounts integrated from the top of the atmosphere

         u1 = -3.0
         du = 0.15
         w1 = -4.0
         dw = 0.15

!-----Inline RFLX
         du2=du*du
         dw2=dw*dw
         x0=u1+real(nu)*du
         y0=w1+real(nw)*dw

         x1=u1-0.5*du
         y1=w1-0.5*dw

         do k= 1, np+1
            ulog=min(log10(so2(k)*snt),x0)
            wlog=min(log10(swh(k)*snt),y0)
            ic=int((ulog-x1)/du+1.)
            iw=int((wlog-y1)/dw+1.)
            if(ic.lt.2)  ic=2
            if(iw.lt.2)  iw=2
            if(ic.gt.nu) ic=nu
            if(iw.gt.nw) iw=nw
            dc=ulog-real(ic-2)*du-u1
            dd=wlog-real(iw-2)*dw-w1   
            x2=cah(ic-1,iw-1)+(cah(ic-1,iw)-cah(ic-1,iw-1))/dw*dd
            y2=x2+(cah(ic,iw-1)-cah(ic-1,iw-1))/du*dc
            y2=max(y2,0.0)
            df(k)=df(k)+1.5*y2 ! LLT increase CO2 effect to help reduce cold tropopause bias
         end do

!-----df is the updated flux reduction for co2 absorption
!     in Band 8 where the co2 absorption has a large impact
!     on the heating of middle atmosphere. From the table
!     given by Eq. (3.19)

         u1 = 0.000250
         du = 0.000050
         w1 = -2.0
         dw = 0.05

!-----Inline RFLX
         du2=du*du
         dw2=dw*dw
         x0=u1+real(nx)*du
         y0=w1+real(ny)*dw

         x1=u1-0.5*du
         y1=w1-0.5*dw

         do k= 1, np+1
            ulog=min(co2*snt,x0)
            wlog=min(log10(pl_dev(i,k)),y0)
            ic=int((ulog-x1)/du+1.)
            iw=int((wlog-y1)/dw+1.)
            if(ic.lt.2)  ic=2
            if(iw.lt.2)  iw=2
            if(ic.gt.nx) ic=nx
            if(iw.gt.ny) iw=ny
            dc=ulog-real(ic-2)*du-u1
            dd=wlog-real(iw-2)*dw-w1   
            x2=coa(ic-1,iw-1)+(coa(ic-1,iw)-coa(ic-1,iw-1))/dw*dd
            y2=x2+(coa(ic,iw-1)-coa(ic-1,iw-1))/du*dc
            y2=max(y2,0.0)
            df(k)=df(k)+1.5*y2 ! LLT increase CO2 effect to help reduce cold tropopause bias
         end do

!-----adjust the o2-co2 reduction below cloud top following Eq. (6.18)

         foundtop = 0

         do k=1,np
            if (fcld_dev(i,k) > 0.02.and.foundtop.eq.0) then
               foundtop = 1
               ntop = k
            end if
         end do

         if (foundtop.eq.0) ntop=np+1

         dftop = df(ntop)

         do k=1,np+1
            if (k .gt. ntop) then
               xx4   = (flx_dev(i,k)/flx_dev(i,ntop))
               df(k) = dftop + xx4 * (df(k)-dftop)
            end if
         end do

!-----update the net fluxes

         do k=1,np+1
            df(k) = min(df(k),flx_dev(i,k)-1.0e-8)
!           df(k) = 0.0
            flx_dev(i,k) = flx_dev(i,k) - df(k)
            flc_dev(i,k) = flc_dev(i,k) - df(k)
         end do

!-----update the downward surface fluxes 

!        xx4 = fdirir (i) + fdifir (i) +&
!              fdiruv (i) + fdifuv (i) +&
!              fdirpar(i) + fdifpar(i)

         xx4 = flx_dev(i,np+1) + df(np+1)

         if ( abs(xx4) > epsilon(1.0) ) then
            xx4 = max(min(1.0 - df(np+1)/xx4,1.),0.)
         else
            xx4 = 0.0
         end if

         fdirir_dev(i)  = xx4*fdirir_dev(i) 
         fdifir_dev(i)  = xx4*fdifir_dev(i) 
         fdiruv_dev(i)  = xx4*fdiruv_dev(i) 
         fdifuv_dev(i)  = xx4*fdifuv_dev(i) 
         fdirpar_dev(i) = xx4*fdirpar_dev(i)
         fdifpar_dev(i) = xx4*fdifpar_dev(i)

         do ib = 1, nband
            flx_sfc_band_dev(i,ib) = xx4*flx_sfc_band_dev(i,ib)
         end do

      !end do RUN_LOOP

   end subroutine sorad


!*********************************************************************

   subroutine deledd(tau1,ssc1,g01,cza1,rr1,tt1,td1)

!*********************************************************************
!
!-----uses the delta-eddington approximation to compute the
!     bulk scattering properties of a single layer
!     coded following King and Harshvardhan (JAS, 1986)
!
!  inputs:
!     tau1:  optical thickness
!     ssc1:  single scattering albedo
!     g01:   asymmetry factor
!     cza1:  cosine o the zenith angle
!
!  outputs:
!
!     rr1:  reflection of the direct beam
!     tt1:  total (direct+diffuse) transmission of the direct beam
!     td1:  direct transmission of the direct beam
!
!*********************************************************************

      implicit none

      integer,parameter :: REAL_DE = 8 ! 8 byte real
      !integer,parameter :: REAL_SP = 4 ! 4 byte real

!-----input parameters

      real(8), intent(in) :: tau1,ssc1,g01,cza1

!-----output parameters

      real(8), intent(out) :: rr1, tt1, td1

!-----temporary parameters

      real*8, parameter :: zero = 0.0_REAL_DE
      real*8, parameter :: one = 1.0_REAL_DE
      real*8, parameter :: two = 2.0_REAL_DE
      real*8, parameter :: three = 3.0_REAL_DE
      real*8, parameter :: four = 4.0_REAL_DE
      real*8, parameter :: fourth = 0.25_REAL_DE
      real*8, parameter :: seven = 7.0_REAL_DE
      real*8, parameter :: thresh = 1.e-8_REAL_DE

      real*8 ::  tau,ssc,g0,rr,tt,td 
      real*8 ::  zth,ff,xx,taup,sscp,gp,gm1,gm2,gm3,akk,alf1,alf2
      real*8 ::  all,bll,st7,st8,cll,dll,fll,ell,st1,st2,st3,st4

      !zth = real(cza1,kind=REAL_DE)
      !g0  = real(g01 ,kind=REAL_DE)
      !tau = real(tau1,kind=REAL_DE)
      !ssc = real(ssc1,kind=REAL_DE)

      !zth = dble(cza1)
      !g0  = dble(g01)
      !tau = dble(tau1)
      !ssc = dble(ssc1)

      zth = cza1
      g0  = g01
      tau = tau1
      ssc = ssc1

      ff  = g0*g0
      xx  = one-ff*ssc
      taup= tau*xx
      sscp= ssc*(one-ff)/xx
      gp  = g0/(one+g0)

      xx  = three*gp 
      gm1 = (seven-sscp*(four+xx))*fourth
      gm2 =-(one  -sscp*(four-xx))*fourth

      akk = sqrt((gm1+gm2)*(gm1-gm2))

      xx  = akk*zth
      st7 = one-xx
      st8 = one+xx
      st3 = st7*st8

      if (abs(st3) .lt. thresh) then
         zth = zth+0.0010
         if(zth > 1.0) zth = zth-0.0020 
         xx  = akk*zth
         st7 = one-xx
         st8 = one+xx
         st3 = st7*st8
      end if

      td=exp(-taup/zth)

      gm3 = (two-zth*three*gp)*fourth
      xx  = gm1-gm2
      alf1= gm1-gm3*xx
      alf2= gm2+gm3*xx

      xx  = akk*two
      all = (gm3-alf2*zth    )*xx*td
      bll = (one-gm3+alf1*zth)*xx

      xx  = akk*gm3
      cll = (alf2+xx)*st7
      dll = (alf2-xx)*st8

      xx  = akk*(one-gm3)
      fll = (alf1+xx)*st8
      ell = (alf1-xx)*st7

      st2 = exp(-akk*taup)
      st4 = st2*st2

      st1 = sscp/((akk+gm1+(akk-gm1)*st4)*st3)

      rr = ( cll-dll*st4     - all*st2)*st1
      tt =-((fll-ell*st4)*td - bll*st2)*st1

      rr = max(rr,zero)
      tt = max(tt,zero)

      tt = tt+td

      !td1 = real(td,kind=REAL_SP)
      !rr1 = real(rr,kind=REAL_SP)
      !tt1 = real(tt,kind=REAL_SP)
      td1 = real(td)
      rr1 = real(rr)
      tt1 = real(tt)

   end subroutine deledd


   subroutine getvistau1(nlevs,cosz,dp,fcld,reff,hydromets,ict,icb,  &
                         taubeam,taudiff,asycl,                      &
                         aig_uv, awg_uv, arg_uv,                     &
                         aib_uv, awb_uv, arb_uv,                     &
                         aib_nir, awb_nir, arb_nir,                  &
                         aia_nir, awa_nir, ara_nir,                  &
                         aig_nir, awg_nir, arg_nir,                  &
                         caib, caif,                                 &
                         CONS_GRAV                                   )

! !USES:

      implicit none

! !INPUT PARAMETERS:
      integer, intent(IN ) :: nlevs              !  Number of levels
      real(8),    intent(IN ) :: cosz               !  Cosine of solar zenith angle
      real(8),    intent(IN ) :: dp(nlevs)          !  Delta pressure (Pa)
      real(8),    intent(IN ) :: fcld(nlevs)        !  Cloud fraction (used sometimes)
      real(8),    intent(IN ) :: reff(nlevs,4)      !  Effective radius (microns)
      real(8),    intent(IN ) :: hydromets(nlevs,4) !  Hydrometeors (kg/kg)
      integer, intent(IN ) :: ict, icb           !  Flags for various uses 
!                 ict  = 0   Indicates that in-cloud values have been given
!                            and are expected
!                 ict != 0   Indicates that overlap computation is needed, and:
!                               ict is the level of the mid-high boundary
!                               icb is the level of the low-mid  boundary
!                
! !OUTPUT PARAMETERS:
      real(8),    intent(OUT) :: taubeam(nlevs,4)    !  Optical Depth for Beam Radiation
      real(8),    intent(OUT) :: taudiff(nlevs,4)    !  Optical Depth for Diffuse Radiation
      real(8),    intent(OUT) ::   asycl(nlevs  )    !  Cloud Asymmetry Factor
! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for visible wavelengths
!  In general will compute in-cloud - will do grid mean when called
!  for diagnostic use only. ict flag will indicate which to do.
!  Slots for reff, hydrometeors, taubeam, taudiff, and asycl are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
!
!  In the below calculations, the constants used in the tau calculation are in 
!  m$^2$ g$^{-1}$ and m$^2$ g$^{-1}$ $\mu$m. Thus, we must convert the kg contained in the 
!  pressure (Pa = kg m$^{-1}$ s$^{-2}$) to grams.
!
! !REVISION HISTORY: 
!    2011.10.27   Molod moved to Radiation_Shared and revised arg list, units
!    2011.11.16   MAT: Generalized to a call that is per-column
!
!EOP
!------------------------------------------------------------------------------
!BOC

      integer            :: k,in,im,it,ia,kk
      real(8)               :: fm,ft,fa,xai,tauc,asyclt
      real(8)               :: cc(3)
      real(8)               :: taucld1,taucld2,taucld3,taucld4
      real(8)               :: g1,g2,g3,g4

      real(8)               :: reff_snow

      integer, parameter :: nm=11,nt=9,na=11
      real(8),    parameter :: dm=0.1,dt=0.30103,da=0.1,t1=-0.9031

      real(8), intent(in) :: aig_uv(3), awg_uv(3), arg_uv(3)
      real(8), intent(in) :: aib_uv, awb_uv(2), arb_uv(2)
      real(8), intent(in) :: aib_nir, awb_nir(3,2), arb_nir(3,2)
      real(8), intent(in) :: aia_nir(3,3), awa_nir(3,3), ara_nir(3,3)
      real(8), intent(in) :: aig_nir(3,3), awg_nir(3,3), arg_nir(3,3)
      real(8), intent(in) :: caib(11,9,11), caif(9,11)
      real(8), intent(in) :: CONS_GRAV

      taubeam = 0.0
      taudiff = 0.0

      if (ict.ne.0) then

!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group

         cc = 0.0

         do k = 1, ict-1
             cc(1)=max(cc(1),fcld(k))
         end do
         do k = ict, icb-1
             cc(2)=max(cc(2),fcld(k))
         end do
         do k = icb, nlevs
             cc(3)=max(cc(3),fcld(k))
         end do

      end if

!-----Compute cloud optical thickness.  Eqs. (4.6) and (4.10)
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     taucld1 is the optical thickness for ice particles
!     taucld2 is the optical thickness for liquid particles
!     taucld3 is the optical thickness for rain drops
!     taucld4 is the optical thickness for snow

      do k = 1, nlevs

         if (reff(k,1) <= 0.) then
            taucld1=0.
         else
            taucld1=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,1))*aib_uv/reff(k,1)
         end if

         if (reff(k,2) <= 0.) then
            taucld2=0.
         else
            taucld2=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,2))*(awb_uv(1)+awb_uv(2)/reff(k,2))
         end if

            taucld3=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,3))*arb_uv(1)

!-----In the IR optical thickness calculation (getirtau.code), it was
!     found that using the table of coefficients tabulated for suspended
!     cloud ice particles (aib_ir) for falling snow lead to unphysical 
!     (negative) values of cloud optical thickness for effective radii 
!     greater than 113 microns. By restricting the effective radius of 
!     snow to 112 microns, we prevent unphysical optical thicknesses. 
!     For consistency's sake, we limit snow effective radius similarly here.

         reff_snow = min(reff(k,4),112.0)

         if (reff_snow <= 0.) then
            taucld4=0.
         else
            taucld4=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,4))*aib_uv/reff_snow
         endif

         if ( ict .ne. 0 ) then

!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group

            if (k.lt.ict) then
               kk=1
            else if (k.ge.ict .and. k.lt.icb) then
               kk=2
            else
               kk=3
            end if
 
            tauc=taucld1+taucld2+taucld3+taucld4

            if (tauc.gt.0.02 .and. fcld(k).gt.0.01) then

!-----normalize cloud cover following Eq. (7.8)

               fa=fcld(k)/cc(kk)

!-----table look-up
 
               tauc=min(tauc,32.)

               fm=cosz/dm 
               ft=(log10(tauc)-t1)/dt
               fa=fa/da

               im=int(fm+1.5) 
               it=int(ft+1.5)
               ia=int(fa+1.5)
  
               im=max(im,2)
               it=max(it,2)
               ia=max(ia,2)
     
               im=min(im,nm-1)
               it=min(it,nt-1)
               ia=min(ia,na-1)
 
               fm=fm-real(im-1)
               ft=ft-real(it-1)
               fa=fa-real(ia-1)

!-----scale cloud optical thickness for beam radiation following 
!     Eq. (7.3).
!     the scaling factor, xai, is a function of the solar zenith
!     angle, optical thickness, and cloud cover.

               xai=    (-caib(im-1,it,ia)*(1.-fm)+&
                         caib(im+1,it,ia)*(1.+fm))*fm*.5+caib(im,it,ia)*(1.-fm*fm)

               xai=xai+(-caib(im,it-1,ia)*(1.-ft)+&
                         caib(im,it+1,ia)*(1.+ft))*ft*.5+caib(im,it,ia)*(1.-ft*ft)

               xai=xai+(-caib(im,it,ia-1)*(1.-fa)+&
                         caib(im,it,ia+1)*(1.+fa))*fa*.5+caib(im,it,ia)*(1.-fa*fa)

               xai=xai-2.*caib(im,it,ia)

               xai=max(xai,0.0)
               xai=min(xai,1.0)

               taubeam(k,1)=taucld1*xai
               taubeam(k,2)=taucld2*xai
               taubeam(k,3)=taucld3*xai
               taubeam(k,4)=taucld4*xai

!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
 
               xai=    (-caif(it-1,ia)*(1.-ft)+&
                         caif(it+1,ia)*(1.+ft))*ft*.5+caif(it,ia)*(1.-ft*ft)
 
               xai=xai+(-caif(it,ia-1)*(1.-fa)+&
                         caif(it,ia+1)*(1.+fa))*fa*.5+caif(it,ia)*(1.-fa*fa)
 
               xai=xai-caif(it,ia)
 
               xai=max(xai,0.0)
               xai=min(xai,1.0)
 
               taudiff(k,1)=taucld1*xai
               taudiff(k,2)=taucld2*xai
               taudiff(k,3)=taucld3*xai
               taudiff(k,4)=taucld4*xai
            end if
         else
         ! Overlap calculation scaling not needed
            taubeam(k,1)=taucld1
            taubeam(k,2)=taucld2
            taubeam(k,3)=taucld3
            taubeam(k,4)=taucld4

            taudiff(k,1)=taucld1
            taudiff(k,2)=taucld2
            taudiff(k,3)=taucld3
            taudiff(k,4)=taucld4
         end if

!-----cloud asymmetry factor for a mixture of liquid and ice particles.
!     unit of reff is micrometers. Eqs. (4.8) and (6.4)

         asyclt=1.0
         tauc=taucld1+taucld2+taucld3+taucld4

         if (tauc.gt.0.02 .and. fcld(k).gt.0.01) then
            g1=(aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff(k,1))*reff(k,1))*taucld1
            g2=(awg_uv(1)+(awg_uv(2)+awg_uv(3)*reff(k,2))*reff(k,2))*taucld2
            g3= arg_uv(1)                                           *taucld3
            g4=(aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff_snow)*reff_snow)*taucld4
            asyclt=(g1+g2+g3+g4)/tauc
         end if

         asycl(k)=asyclt

      end do

      return

!EOC
   end subroutine getvistau1



   subroutine getnirtau1(ib,nlevs,cosz,dp,fcld,reff,hydromets,ict,icb, &
                         taubeam,taudiff,asycl,ssacl,                  &
                         aig_uv, awg_uv, arg_uv,                       &
                         aib_uv, awb_uv, arb_uv,                       &
                         aib_nir, awb_nir, arb_nir,                    &
                         aia_nir, awa_nir, ara_nir,                    &
                         aig_nir, awg_nir, arg_nir,                    &
                         caib, caif,                                   &
                         CONS_GRAV                                     )


      implicit none

! !INPUT PARAMETERS:
      integer, intent(IN ) :: ib                 !  Band number
      integer, intent(IN ) :: nlevs              !  Number of levels
      real(8),    intent(IN ) :: cosz               !  Cosine of solar zenith angle
      real(8),    intent(IN ) :: dp(nlevs)          !  Delta pressure in Pa
      real(8),    intent(IN ) :: fcld(nlevs)        !  Cloud fraction (used sometimes)
      real(8),    intent(IN ) :: reff(nlevs,4)      !  Effective radius (microns)
      real(8),    intent(IN ) :: hydromets(nlevs,4) !  Hydrometeors (kg/kg)
      integer, intent(IN ) :: ict, icb           !  Flags for various uses 

      real(8), intent(in) :: aig_uv(3), awg_uv(3), arg_uv(3)
      real(8), intent(in) :: aib_uv, awb_uv(2), arb_uv(2)
      real(8), intent(in) :: aib_nir, awb_nir(3,2), arb_nir(3,2)
      real(8), intent(in) :: aia_nir(3,3), awa_nir(3,3), ara_nir(3,3)
      real(8), intent(in) :: aig_nir(3,3), awg_nir(3,3), arg_nir(3,3)
      real(8), intent(in) :: caib(11,9,11), caif(9,11)
      real(8), intent(in) :: CONS_GRAV

! !OUTPUT PARAMETERS:
      real(8),    intent(OUT) :: taubeam(nlevs,4)   !  Optical depth for beam radiation
      real(8),    intent(OUT) :: taudiff(nlevs,4)   !  Optical depth for diffuse radiation
      real(8),    intent(OUT) ::   ssacl(nlevs  )   !  Cloud single scattering albedo
      real(8),    intent(OUT) ::   asycl(nlevs  )   !  Cloud asymmetry factor

      integer            :: k,in,im,it,ia,kk
      real(8)               :: fm,ft,fa,xai,tauc,asyclt,ssaclt
      real(8)               :: cc(3)
      real(8)               :: taucld1,taucld2,taucld3,taucld4
      real(8)               :: g1,g2,g3,g4
      real(8)               :: w1,w2,w3,w4

      real(8)               :: reff_snow

      integer, parameter :: nm=11,nt=9,na=11
      real(8),    parameter :: dm=0.1,dt=0.30103,da=0.1,t1=-0.9031

      taubeam = 0.0
      taudiff = 0.0

      if (ict.ne.0) then

!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group

         cc = 0.0

         do k = 1, ict-1
             cc(1)=max(cc(1),fcld(k))
         end do
         do k = ict, icb-1
             cc(2)=max(cc(2),fcld(k))
         end do
         do k = icb, nlevs
             cc(3)=max(cc(3),fcld(k))
         end do

      end if

!-----Compute cloud optical thickness.  Eqs. (4.6) and (4.10)
!     taucld1 is the optical thickness for ice particles
!     taucld2 is the optical thickness for liquid particles
!     taucld3 is the optical thickness for rain drops
!     taucld4 is the optical thickness for snow

      do k = 1, nlevs

         if (reff(k,1) <= 0.) then
            taucld1=0.
         else
            taucld1=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,1))*aib_nir/reff(k,1)
         end if

         if (reff(k,2) <= 0.) then
            taucld2=0.
         else
            taucld2=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,2))*(awb_nir(ib,1)+awb_nir(ib,2)/reff(k,2))
         end if

            taucld3=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,3))*arb_nir(ib,1)

!-----In the IR optical thickness calculation (getirtau.code), it was 
!     found that using the table of coefficients tabulated for suspended
!     cloud ice particles (aib_ir) for falling snow lead to unphysical 
!     (negative) values of cloud optical thickness for effective radii 
!     greater than 113 microns. By restricting the effective radius of  
!     snow to 112 microns, we prevent unphysical optical thicknesses. 
!     For consistency's sake, we limit snow effective radius similarly here.

         reff_snow = min(reff(k,4),112.0)

         if (reff_snow <= 0.) then
            taucld4=0.
         else
            taucld4=(((dp(k)*1.0e3)/CONS_GRAV)*hydromets(k,4))*aib_nir/reff_snow
         endif

         if ( ict .ne. 0 ) then

!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group

            if (k.lt.ict) then
               kk=1
            else if (k.ge.ict .and. k.lt.icb) then
               kk=2
            else
               kk=3
            end if
 
            tauc=taucld1+taucld2+taucld3+taucld4

            if (tauc.gt.0.02 .and. fcld(k).gt.0.01) then

!-----normalize cloud cover following Eq. (7.8)
             if (cc(kk).ne.0.0) then
                 fa=fcld(k)/cc(kk)
	     else
	         fa=0.0
	     end if

!-----table look-up
 
               tauc=min(tauc,32.)

               fm=cosz/dm 
               ft=(log10(tauc)-t1)/dt
               fa=fa/da

               im=int(fm+1.5) 
               it=int(ft+1.5)
               ia=int(fa+1.5)
  
               im=max(im,2)
               it=max(it,2)
               ia=max(ia,2)
     
               im=min(im,nm-1)
               it=min(it,nt-1)
               ia=min(ia,na-1)
 
               fm=fm-real(im-1)
               ft=ft-real(it-1)
               fa=fa-real(ia-1)

!-----scale cloud optical thickness for beam radiation following 
!     Eq. (7.3).
!     the scaling factor, xai, is a function of the solar zenith
!     angle, optical thickness, and cloud cover.

               xai=    (-caib(im-1,it,ia)*(1.-fm)+&
                         caib(im+1,it,ia)*(1.+fm))*fm*.5+caib(im,it,ia)*(1.-fm*fm)

               xai=xai+(-caib(im,it-1,ia)*(1.-ft)+&
                         caib(im,it+1,ia)*(1.+ft))*ft*.5+caib(im,it,ia)*(1.-ft*ft)

               xai=xai+(-caib(im,it,ia-1)*(1.-fa)+&
                         caib(im,it,ia+1)*(1.+fa))*fa*.5+caib(im,it,ia)*(1.-fa*fa)

               xai=xai-2.*caib(im,it,ia)

               xai=max(xai,0.0)
               xai=min(xai,1.0)

               taubeam(k,1)=taucld1*xai
               taubeam(k,2)=taucld2*xai
               taubeam(k,3)=taucld3*xai
               taubeam(k,4)=taucld4*xai

!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
 
               xai=    (-caif(it-1,ia)*(1.-ft)+&
                         caif(it+1,ia)*(1.+ft))*ft*.5+caif(it,ia)*(1.-ft*ft)
 
               xai=xai+(-caif(it,ia-1)*(1.-fa)+&
                         caif(it,ia+1)*(1.+fa))*fa*.5+caif(it,ia)*(1.-fa*fa)
 
               xai=xai-caif(it,ia)
 
               xai=max(xai,0.0)
               xai=min(xai,1.0)
 
               taudiff(k,1)=taucld1*xai
               taudiff(k,2)=taucld2*xai
               taudiff(k,3)=taucld3*xai
               taudiff(k,4)=taucld4*xai
            end if
         else
         ! Overlap calculation scaling not needed
            taubeam(k,1)=taucld1
            taubeam(k,2)=taucld2
            taubeam(k,3)=taucld3
            taubeam(k,4)=taucld4

            taudiff(k,1)=taucld1
            taudiff(k,2)=taucld2
            taudiff(k,3)=taucld3
            taudiff(k,4)=taucld4
         end if

!-----compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)

         ssaclt=0.99999
         asyclt=1.0
         tauc=taucld1+taucld2+taucld3+taucld4

         if (tauc.gt.0.02 .and. fcld(k).gt.0.01) then

            w1=(1.-(aia_nir(ib,1)+(aia_nir(ib,2)+aia_nir(ib,3)*reff(k,1))*reff(k,1)))*taucld1
            w2=(1.-(awa_nir(ib,1)+(awa_nir(ib,2)+awa_nir(ib,3)*reff(k,2))*reff(k,2)))*taucld2
            w3=(1.- ara_nir(ib,1))                                                   *taucld3
            w4=(1.-(aia_nir(ib,1)+(aia_nir(ib,2)+aia_nir(ib,3)*reff_snow)*reff_snow))*taucld4
            ssaclt=(w1+w2+w3+w4)/tauc

            g1=(aig_nir(ib,1)+(aig_nir(ib,2)+aig_nir(ib,3)*reff(k,1))*reff(k,1))*w1
            g2=(awg_nir(ib,1)+(awg_nir(ib,2)+awg_nir(ib,3)*reff(k,2))*reff(k,2))*w2
            g3= arg_nir(ib,1)                                                   *w3

            g4=(aig_nir(ib,1)+(aig_nir(ib,2)+aig_nir(ib,3)*reff(k,4))*reff(k,4))*w4
	    
	    if ((w1+w2+w3+w4).ne.0.0) then
             asyclt=(g1+g2+g3+g4)/(w1+w2+w3+w4)
	    end if
	     

         end if

         ssacl(k)=ssaclt
         asycl(k)=asyclt

      end do

      return

   end subroutine getnirtau1

