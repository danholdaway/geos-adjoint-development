subroutine radiation_driver ( IM,JM,LM,RUNalarm,DT,UT,PTT,QVT,O3T,CFLST,CFCNT,QILST,QLLST,QICNT,QLCNT,PS,TS, &
                              EMIS, DELT, COSZ, SLR, RGBUV, RGFUV, RGBIR, RGFIR , SC, &
                              TAUA_IRT, SSAA_IRT, ASYA_IRT, TAUA_SOT, SSAA_SOT, ASYA_SOT, &
                              mapl_p00, mapl_rgas, mapl_cp, mapl_p00, mapl_kappa, mapl_grav, &
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
                              HK, HK_IR_TEMP, HK_UV_TEMP, &
                              HK_UV_OLD, HK_IR_OLD )


IMPLICIT NONE

!In
integer, intent(in) :: IM,JM,LM,RUNalarm
real, intent(in) :: SC, DT
real, intent(in), dimension(IM,JM,LM) :: UT, QVT, O3T, CFLST, CFCNT, QILST,QLLST,QICNT,QLCNT
real, intent(in), dimension(IM,JM) :: PS,TS
real, intent(in), dimension(IM,JM) :: EMIS, DELT, COSZ, SLR, RGBUV, RGFUV, RGBIR, RGFIR

real, intent(in), dimension(IM,JM,LM,10) :: TAUA_IRT, SSAA_IRT, ASYA_IRT
real, intent(in), dimension(IM,JM,LM,8)  :: TAUA_SOT, SSAA_SOT, ASYA_SOT

real, intent(in) :: mapl_p00, mapl_rgas, mapl_cp, mapl_p00, mapl_kappa, mapl_grav

 real, intent(IN)    :: aib_ir(3,10), awb_ir(4,10), aiw_ir(4,10)
 real, intent(IN)    :: aww_ir(4,10), aig_ir(4,10), awg_ir(4,10)

 integer, intent(IN) :: mw(9)
 real, intent(IN)    :: xkw(9), xke(9), aw(9), bw(9), pm(9)
 real, intent(IN)    :: fkw(6,9), gkw(6,3), cb(6,10), dcb(5,10)
 real, intent(IN)    :: w11, w12, w13, p11, p12
 real, intent(IN)    :: p13, dwe, dpe
 real, intent(IN)    :: c1(26,30),  c2(26,30),  c3(26,30)
 real, intent(IN)    :: oo1(26,21), oo2(26,21), oo3(26,21)
 real, intent(IN)    :: h11(26,31), h12(26,31), h13(26,31)
 real, intent(IN)    :: h21(26,31), h22(26,31), h23(26,31)
 real, intent(IN)    :: h81(26,31), h82(26,31), h83(26,31)

real, intent(in) :: wk_uv(5), zk_uv(5), ry_uv(5)
real, intent(in) :: xk_ir(10), ry_ir(3)
real, intent(in) :: cah(43,37), coa(62,101)

real, intent(in) :: aig_uv(3), awg_uv(3), arg_uv(3)
real, intent(in) :: aib_uv, awb_uv(2), arb_uv(2)
real, intent(in) :: aib_nir, awb_nir(3,2), arb_nir(3,2)
real, intent(in) :: aia_nir(3,3), awa_nir(3,3), ara_nir(3,3)
real, intent(in) :: aig_nir(3,3), awg_nir(3,3), arg_nir(3,3)
real, intent(in) :: caib(11,9,11), caif(9,11)

real, intent(in) ::  HK(8), HK_IR_TEMP(3,10), HK_UV_TEMP(5)

real, intent(in) ::  HK_UV_OLD(5), HK_IR_OLD(3,10)

!Inouts
real, intent(inout), dimension(IM,JM,LM) :: PTT

!Locals
integer :: i,j,k,l
real :: CO2
integer :: levs925, LCLDMH, LCLDLM
logical, parameter :: TRACE = .true., OVERCAST = .false. 
  integer, parameter :: NA = 0, NS = 1
  integer, parameter :: NBCHOU_IR = 10
  integer, parameter :: NBCHOU_SO = 8



real, dimension(IM,JM,LM) :: PTT1, TT, PLE, PLESO, PLOf, PIf, DELTAP, TT_OUT, DTDT
real, dimension(IM,JM,LM) :: O3mmT
real, dimension(IM,JM,LM) :: N2OT,CH4T,CFC11T,CFC12T,CFC22T
real, dimension(IM,JM,LM) :: RAD_QLT, RAD_QIT, RAD_CFT, RAD_RLT, RAD_RIT
real, dimension(IM,JM,LM,4) :: CWCT, REFFT
real, dimension(IM,JM,1)    :: FST, TGT, TVT
real, dimension(IM,JM,1,10) :: EGT, EVT, RVT

real, dimension(IM,JM,LM) :: FLWUT, FLWDT, FLXNT, FLWT_INT, FSWT_INT
real, dimension(IM,JM,LM) :: flwavet, fswavet, dfdtst_int


real, dimension(IM,JM) :: T2MT, TEMPOR, COSZTMP

real, dimension(0:LM) :: ak, bk, pref



       !Compute pref from ak and bk
       ak = 1.0
       bk = 1.0
       DO L = 0,LM
          pref(L) = AK(L+1) + BK(L+1)*MAPL_P00
       enddo

       !Pressure at the half levels from Ps
       DO L = 0,LM
          PLE(:,:,L) = AK(L+1) + BK(L+1)*PS(:,:)
       enddo
       DELTAP = ( PLE(:,:,1:LM) - PLE(:,:,0:LM-1) )

       !Pressure (hPa) and Exner pressure at the full levels
       PLOf(:,:,1:LM) = 0.01 * 0.5 * (PLE(:,:,0:LM-1) +  PLE(:,:,1:LM  ) )
       PIf(:,:,1:LM) = (PLOf(:,:,1:LM)/1000.)**(MAPL_RGAS/MAPL_CP)

       !Potential temperature with p0=10000
       PTT1 = PTT*(MAPL_P00**MAPL_KAPPA) !Place in holder so as not to overwrite

       !Temperature
       TT = PTT1*PIf !Save input for later.


if (RUNalarm == 1) then


          !Initialize the fluxes
          FLWUT = 0.0
          FLWDT = 0.0
          DFDTST_INT = 0.0

DO i = 1,IM
   DO j = 1,JM

 call irrad ( 1, LM, TT, FLWUT, FLWDT, DFDTST_INT)

   enddo
enddo


             FLWT_INT = FLWUT + FLWDT


endif


          do L = 0, LM
             FLWAVET(:,:,L) = FLWT_INT(:,:,L) + DFDTST_INT(:,:,L)*DELT
          end do
 
          !Tangent of fluxes to DTDt (mass-weighted)
          DTDt = ( (FLWAVET(:,:,0:LM-1) - FLWAVET(:,:,1:LM)) ) * (MAPL_GRAV/MAPL_CP)

          !Tangent of mass-weighted DTDt to updated temperature
          TT_out = TT + DT * DTDt / DELTAP

          !Tangent of temperature to potential temperature
          PTT1 = TT_out/PIf

          !Tangent of potential temperature to P0 = 1
          PTT = PTT1/(MAPL_P00**MAPL_KAPPA)

endsubroutine radiation_driver



 subroutine irrad ( m, np, ta_dev, flxu_dev, flxd_dev, dfdts_dev)

 IMPLICIT NONE

 !----- INPUTS -----
 INTEGER, INTENT(IN)                     :: m, np

 REAL, DIMENSION(m,np),    INTENT(IN)    :: ta_dev

 REAL, DIMENSION(m,np+1), INTENT(OUT)    :: flxu_dev
 REAL, DIMENSION(m,np+1), INTENT(OUT)    :: flxd_dev
 REAL, DIMENSION(m,np+1), INTENT(OUT)    :: dfdts_dev

 INTEGER :: i,j

 flxu_dev = 0.0
 flxd_dev = 0.0
 dfdts_dev = 0.0

 DO i = 1,m
    DO j = 1,np


       flxu_dev(i,j) = ta_dev(i,j) * 1e-6 + 1e-8*ta_dev(i,j)**2
       flxd_dev(i,j) = ta_dev(i,j) * 1e-4 + 1e-10*ta_dev(i,j)**3
       dfdts_dev(i,j) = 0.5*(flxu_dev(i,j) + flxd_dev(i,j))


    endDO
 endDO

   end subroutine irrad

