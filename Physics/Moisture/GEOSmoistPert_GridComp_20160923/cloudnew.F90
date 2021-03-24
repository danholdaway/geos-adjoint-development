! $Id: cloudnew.F90,v 1.40.2.7.2.2.10.1.4.1.4.8.2.1.2.3.6.10.24.1.4.1.2.1.2.1.12.3.4.5.12.1.2.42.2.18.2.12.2.16.2.4.8.1.42.3.10.1.46.2.2.2.16.1.6.4.12.1.2.1 2016-09-08 17:43:05 ltakacs Exp $
! $Name: drh-GEOSadas-5_16_0_H54 $

module cloudnew

#ifndef _CUDA
   use GEOS_UtilsMod,     only: QSAT=>GEOS_Qsat, DQSAT=>GEOS_DQsat
   use CLDPARAMS
#else
   use cudafor
   ! NOTE: GPUs use the QSAT and DQSAT at the end of this module
#endif

   use MAPL_ConstantsMod, only: MAPL_TICE , MAPL_CP   , &
                                MAPL_GRAV , MAPL_ALHS , &
                                MAPL_ALHL , MAPL_ALHF , &
                                MAPL_RGAS , MAPL_H2OMW, &
                                MAPL_AIRMW, MAPL_RVAP , &
                                MAPL_PI   , MAPL_R8   , &
                                MAPL_R4

   use MAPL_BaseMod,      only: MAPL_UNDEF

   implicit none

#ifndef _CUDA
   private

   !PUBLIC PROGNO_CLOUD
   !PUBLIC ICE_FRACTION
   PUBLIC T_CLOUD_CTL
#endif

   type T_CLOUD_CTL
      real  :: SCLMFDFR
      real  :: RSUB_RADIUS
   end type T_CLOUD_CTL

#ifdef _CUDA

   ! Inputs
   ! ------

   real, allocatable, dimension(:,:), device :: PP_dev
   real, allocatable, dimension(:,:), device :: EXNP_dev
   real, allocatable, dimension(:,:), device :: PPE_dev
   real, allocatable, dimension(:,:), device :: KH_dev
   real, allocatable, dimension(:  ), device :: FRLAND_dev
   real, allocatable, dimension(:,:), device :: RMFDTR_dev
   real, allocatable, dimension(:,:), device :: QLWDTR_dev
   real, allocatable, dimension(:,:), device :: U_dev
   real, allocatable, dimension(:,:), device :: V_dev
   real, allocatable, dimension(:,:), device :: QST3_dev
   real, allocatable, dimension(:,:), device :: DZET_dev
   real, allocatable, dimension(:,:), device :: QDDF3_dev
   real, allocatable, dimension(:  ), device :: TEMPOR_dev
   real, allocatable, dimension(:  ), device :: CNV_FRACTION_dev

   ! Inoutputs
   ! ---------

   real, allocatable, dimension(:,:), device :: TH_dev
   real, allocatable, dimension(:,:), device :: Q_dev
   real, allocatable, dimension(:,:), device :: QRN_CU_dev
   real, allocatable, dimension(:,:), device :: CNV_UPDFRC_dev ! on edges, but dims=1:LM
   real, allocatable, dimension(:,:), device :: QLW_LS_dev  
   real, allocatable, dimension(:,:), device :: QLW_AN_dev
   real, allocatable, dimension(:,:), device :: QIW_LS_dev  
   real, allocatable, dimension(:,:), device :: QIW_AN_dev
   real, allocatable, dimension(:,:), device :: ANVFRC_dev
   real, allocatable, dimension(:,:), device :: CLDFRC_dev

   ! Outputs
   ! -------

   real, allocatable, dimension(:,:), device :: RAD_CLDFRC_dev
   real, allocatable, dimension(:,:), device :: RAD_QL_dev
   real, allocatable, dimension(:,:), device :: RAD_QI_dev
   real, allocatable, dimension(:,:), device :: RAD_QR_dev
   real, allocatable, dimension(:,:), device :: RAD_QS_dev
   real, allocatable, dimension(:,:), device :: QPLS_dev
   real, allocatable, dimension(:,:), device :: CLDREFFL_dev
   real, allocatable, dimension(:,:), device :: CLDREFFI_dev
   real, allocatable, dimension(:  ), device :: PRELS_dev
   real, allocatable, dimension(:  ), device :: PRECU_dev
   real, allocatable, dimension(:  ), device :: PREAN_dev
   real, allocatable, dimension(:  ), device :: LSARF_dev
   real, allocatable, dimension(:  ), device :: CUARF_dev
   real, allocatable, dimension(:  ), device :: ANARF_dev
   real, allocatable, dimension(:  ), device :: SNRLS_dev
   real, allocatable, dimension(:  ), device :: SNRCU_dev
   real, allocatable, dimension(:  ), device :: SNRAN_dev

   real, allocatable, dimension(:,:), device :: PFL_CN_dev
   real, allocatable, dimension(:,:), device :: PFI_CN_dev
   real, allocatable, dimension(:,:), device :: PFL_AN_dev
   real, allocatable, dimension(:,:), device :: PFI_AN_dev
   real, allocatable, dimension(:,:), device :: PFL_LS_dev
   real, allocatable, dimension(:,:), device :: PFI_LS_dev

   real, allocatable, dimension(:,:), device :: RHX_dev
   real, allocatable, dimension(:,:), device :: REV_LS_dev
   real, allocatable, dimension(:,:), device :: REV_AN_dev
   real, allocatable, dimension(:,:), device :: REV_CN_dev
   real, allocatable, dimension(:,:), device :: RSU_LS_dev
   real, allocatable, dimension(:,:), device :: RSU_AN_dev
   real, allocatable, dimension(:,:), device :: RSU_CN_dev
   real, allocatable, dimension(:,:), device :: ACLL_CN_dev
   real, allocatable, dimension(:,:), device :: ACIL_CN_dev
   real, allocatable, dimension(:,:), device :: ACLL_AN_dev
   real, allocatable, dimension(:,:), device :: ACIL_AN_dev
   real, allocatable, dimension(:,:), device :: ACLL_LS_dev
   real, allocatable, dimension(:,:), device :: ACIL_LS_dev
   real, allocatable, dimension(:,:), device :: PDFL_dev
   real, allocatable, dimension(:,:), device :: PDFI_dev
   real, allocatable, dimension(:,:), device :: FIXL_dev
   real, allocatable, dimension(:,:), device :: FIXI_dev                          
   real, allocatable, dimension(:,:), device :: AUT_dev
   real, allocatable, dimension(:,:), device :: EVAPC_dev
   real, allocatable, dimension(:,:), device :: SDM_dev
   real, allocatable, dimension(:,:), device :: SUBLC_dev 
   real, allocatable, dimension(:,:), device :: FRZ_TT_dev
   real, allocatable, dimension(:,:), device :: FRZ_PP_dev
   real, allocatable, dimension(:,:), device :: DCNVL_dev
   real, allocatable, dimension(:,:), device :: DCNVi_dev
   real, allocatable, dimension(:,:), device :: ALPHT_dev
   real, allocatable, dimension(:,:), device :: ALPH1_dev
   real, allocatable, dimension(:,:), device :: ALPH2_dev
   real, allocatable, dimension(:,:), device :: CFPDF_dev
   real, allocatable, dimension(:,:), device :: RHCLR_dev
   real, allocatable, dimension(:,:), device :: DQRL_dev
   real, allocatable, dimension(:,:), device :: VFALLICE_AN_dev
   real, allocatable, dimension(:,:), device :: VFALLICE_LS_dev
   real, allocatable, dimension(:,:), device :: VFALLWAT_AN_dev
   real, allocatable, dimension(:,:), device :: VFALLWAT_LS_dev
   real, allocatable, dimension(:,:), device :: VFALLSN_AN_dev
   real, allocatable, dimension(:,:), device :: VFALLSN_LS_dev
   real, allocatable, dimension(:,:), device :: VFALLSN_CN_dev
   real, allocatable, dimension(:,:), device :: VFALLRN_AN_dev
   real, allocatable, dimension(:,:), device :: VFALLRN_LS_dev
   real, allocatable, dimension(:,:), device :: VFALLRN_CN_dev

   real, allocatable, dimension(:,:), device :: LIQANMOVE_dev
   real, allocatable, dimension(:,:), device :: ICEANMOVE_dev
   real, allocatable, dimension(:,:), device :: DANCLD_dev
   real, allocatable, dimension(:,:), device :: DLSCLD_dev
   real, allocatable, dimension(:,:), device :: CURAINMOVE_dev
   real, allocatable, dimension(:,:), device :: CUSNOWMOVE_dev

   ! Constants passed in from CLDPARAMS (from MAPL_GetResource)
   real,    constant :: CNV_BETA
   real,    constant :: ANV_BETA
   real,    constant :: LS_BETA
   real,    constant :: RH00
   real,    constant :: C_00
   real,    constant :: LWCRIT
   real,    constant :: C_ACC
   real,    constant :: C_EV_R
   real,    constant :: C_EV_S
   real,    constant :: CLDVOL2FRC
   real,    constant :: RHSUP_ICE
   real,    constant :: SHR_EVAP_FAC
   real,    constant :: MIN_CLD_WATER
   real,    constant :: CLD_EVP_EFF
   integer, constant :: NSMAX
   real,    constant :: LS_SDQV2
   real,    constant :: LS_SDQV3
   real,    constant :: LS_SDQVT1
   real,    constant :: ANV_SDQV2
   real,    constant :: ANV_SDQV3
   real,    constant :: ANV_SDQVT1
   real,    constant :: ANV_TO_LS
   real,    constant :: N_WARM
   real,    constant :: N_ICE
   real,    constant :: N_ANVIL
   real,    constant :: N_PBL
   integer, constant :: DISABLE_RAD
   integer, constant :: ICE_SETTLE
   real,    constant :: ANV_ICEFALL_C
   real,    constant :: LS_ICEFALL_C
   real,    constant :: REVAP_OFF_P
   real,    constant :: CNVENVFC
   real,    constant :: WRHODEP
   real,    constant :: T_ICE_ALL
   real,    constant :: CNVICEPARAM
   integer, constant :: ICEFRPWR
   real,    constant :: CNVDDRFC
   real,    constant :: ANVDDRFC
   real,    constant :: LSDDRFC
   integer, constant :: TANHRHCRIT
   real,    constant :: MINRHCRIT
   real,    constant :: MAXRHCRIT
   real,    constant :: TURNRHCRIT
   real,    constant :: MAXRHCRITLAND
   integer, constant :: FR_LS_WAT
   integer, constant :: FR_LS_ICE
   integer, constant :: FR_AN_WAT
   integer, constant :: FR_AN_ICE
   real,    constant :: MIN_RL
   real,    constant :: MIN_RI
   real,    constant :: MAX_RL
   real,    constant :: MAX_RI
   real,    constant :: RI_ANV
   integer, constant :: PDFFLAG

   ! Parameters for Internal DQSAT
   ! -----------------------------

   real, parameter :: ESFAC            = MAPL_H2OMW/MAPL_AIRMW
   real, parameter :: MAX_MIXING_RATIO = 1.
   real, parameter :: ZEROC            = MAPL_TICE

   real, parameter :: TMINTBL   =  150.0
   real, parameter :: TMAXTBL   =  333.0
   real, parameter :: DEGSUBS   =  100
   real, parameter :: ERFAC     = (DEGSUBS/ESFAC)
   real, parameter :: DELTA_T   =  1.0 / DEGSUBS
   real, parameter :: TABLESIZE =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1
   real, parameter :: TMIX      = -20.

   real, parameter :: TMINSTR = -95.
   real, parameter :: TSTARR1 = -75.
   real, parameter :: TSTARR2 = -65.
   real, parameter :: TSTARR3 = -50.
   real, parameter :: TSTARR4 = -40.
   real, parameter :: TMAXSTR = +60.

   real(kind=MAPL_R8), parameter :: B6 = 6.136820929E-11*100.0
   real(kind=MAPL_R8), parameter :: B5 = 2.034080948E-8 *100.0
   real(kind=MAPL_R8), parameter :: B4 = 3.031240396E-6 *100.0
   real(kind=MAPL_R8), parameter :: B3 = 2.650648471E-4 *100.0
   real(kind=MAPL_R8), parameter :: B2 = 1.428945805E-2 *100.0
   real(kind=MAPL_R8), parameter :: B1 = 4.436518521E-1 *100.0
   real(kind=MAPL_R8), parameter :: B0 = 6.107799961E+0 *100.0
   real(kind=MAPL_R8), parameter :: BI6= 1.838826904E-10*100.0
   real(kind=MAPL_R8), parameter :: BI5= 4.838803174E-8 *100.0
   real(kind=MAPL_R8), parameter :: BI4= 5.824720280E-6 *100.0
   real(kind=MAPL_R8), parameter :: BI3= 4.176223716E-4 *100.0
   real(kind=MAPL_R8), parameter :: BI2= 1.886013408E-2 *100.0
   real(kind=MAPL_R8), parameter :: BI1= 5.034698970E-1 *100.0
   real(kind=MAPL_R8), parameter :: BI0= 6.109177956E+0 *100.0
   real(kind=MAPL_R8), parameter :: S16= 0.516000335E-11*100.0
   real(kind=MAPL_R8), parameter :: S15= 0.276961083E-8 *100.0
   real(kind=MAPL_R8), parameter :: S14= 0.623439266E-6 *100.0
   real(kind=MAPL_R8), parameter :: S13= 0.754129933E-4 *100.0
   real(kind=MAPL_R8), parameter :: S12= 0.517609116E-2 *100.0
   real(kind=MAPL_R8), parameter :: S11= 0.191372282E+0 *100.0
   real(kind=MAPL_R8), parameter :: S10= 0.298152339E+1 *100.0
   real(kind=MAPL_R8), parameter :: S26= 0.314296723E-10*100.0
   real(kind=MAPL_R8), parameter :: S25= 0.132243858E-7 *100.0
   real(kind=MAPL_R8), parameter :: S24= 0.236279781E-5 *100.0
   real(kind=MAPL_R8), parameter :: S23= 0.230325039E-3 *100.0
   real(kind=MAPL_R8), parameter :: S22= 0.129690326E-1 *100.0
   real(kind=MAPL_R8), parameter :: S21= 0.401390832E+0 *100.0
   real(kind=MAPL_R8), parameter :: S20= 0.535098336E+1 *100.0

   real(kind=MAPL_R8), parameter :: DI(0:3) = (/ 57518.5606E08, 2.01889049, 3.56654, 20.947031 /)
   real(kind=MAPL_R8), parameter :: CI(0:3) = (/ 9.550426, -5723.265, 3.53068, -.00728332 /)
   real(kind=MAPL_R8), parameter :: DL(1:6) = (/ -7.902980, 5.02808, -1.3816, 11.344, 8.1328, -3.49149 /)
   real(kind=MAPL_R8), parameter :: LOGPS   = 3.005714898  ! log10(1013.246)
   real(kind=MAPL_R8), parameter :: TS      = 373.16
   real(kind=MAPL_R8), parameter :: CL(0:9) = (/54.842763, -6763.22, -4.21000, .000367, &
                                       .0415, 218.8,  53.878000, -1331.22, -9.44523, .014025  /)

   real, parameter :: TMINLQU = MAPL_TICE - 40.0
   real, parameter :: TMINICE = MAPL_TICE + -95.

#else

!! Some parameters set by CLDPARAMS 

   real    :: CNV_BETA
   real    :: ANV_BETA
   real    :: LS_BETA
   real    :: RH00
   real    :: C_00
   real    :: LWCRIT
   real    :: C_ACC
   real    :: C_EV_R
   real    :: C_EV_S
   real    :: CLDVOL2FRC
   real    :: RHSUP_ICE
   real    :: SHR_EVAP_FAC
   real    :: MIN_CLD_WATER
   real    :: CLD_EVP_EFF
   integer :: NSMAX
   real    :: LS_SDQV2
   real    :: LS_SDQV3
   real    :: LS_SDQVT1
   real    :: ANV_SDQV2
   real    :: ANV_SDQV3
   real    :: ANV_SDQVT1
   real    :: ANV_TO_LS
   real    :: N_WARM
   real    :: N_ICE
   real    :: N_ANVIL
   real    :: N_PBL
   integer :: DISABLE_RAD
   integer :: ICE_SETTLE
   real    :: ANV_ICEFALL_C
   real    :: LS_ICEFALL_C
   real    :: REVAP_OFF_P
   real    :: CNVENVFC
   real    :: WRHODEP
   real    :: T_ICE_ALL
   real    :: CNVICEPARAM
   integer :: ICEFRPWR
   real    :: CNVDDRFC
   real    :: ANVDDRFC
   real    :: LSDDRFC
   integer :: tanhrhcrit
   real    :: minrhcrit
   real    :: maxrhcrit
   real    :: turnrhcrit
   real    :: MIN_RI, MAX_RI, MIN_RL, MAX_RL, RI_ANV
   integer :: FR_LS_WAT, FR_LS_ICE, FR_AN_WAT, FR_AN_ICE
   real    :: maxrhcritland
   integer :: pdfflag

#endif

   real, parameter :: T_ICE_MAX    = MAPL_TICE  ! -7.0+MAPL_TICE
   real, parameter :: RHO_W        = 1.0e3      ! Density of liquid water in kg/m^3
   real, parameter :: MIN_CLD_FRAC = 1.0e-8

   ! There are two PIs in this routine: PI_0 and MAPL_PI
   real, parameter :: PI_0 = 4.*atan(1.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains


end module cloudnew
 
