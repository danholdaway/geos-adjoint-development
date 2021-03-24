! $Id: cloudnew.F90,v 1.40.2.7.2.2.10.1.4.1.4.8.2.1.2.3.6.10.24.1.4.1.2.1.2.1.12.3.4.5.12.1.2.42.2.18.2.12.2.16.2.4.8.1.42.3.10.1.46.2.2.2.16.1.6.4.12.1.2.1 2016-09-08 17:43:05 ltakacs Exp $
! $Name: drh-GEOSadas-5_16_0_H54 $

module cloudnew

!#ifndef _CUDA
!   use GEOS_UtilsMod,     only: QSAT=>GEOS_Qsat, DQSAT=>GEOS_DQsat
   use CLDPARAMS
!#else
!   use cudafor
!   ! NOTE: GPUs use the QSAT and DQSAT at the end of this module
!#endif

   use MAPL_ConstantsMod, only: MAPL8_TICE , MAPL8_CP   , &
                                MAPL8_GRAV , MAPL8_ALHS , &
                                MAPL8_ALHL , MAPL8_ALHF , &
                                MAPL8_RGAS , MAPL8_H2OMW, &
                                MAPL8_AIRMW, MAPL8_RVAP , &
                                MAPL8_PI   , MAPL8_R8   , &
                                MAPL8_R4

!   use MAPL_BaseMod,      only: MAPL_UNDEF

   use qsat_util

   implicit none

!#ifndef _CUDA
   private

   PUBLIC PROGNO_CLOUD, pre_progno_cloud
   !PUBLIC ICE_FRACTION
   PUBLIC T_CLOUD_CTL
!#endif

   type T_CLOUD_CTL
      real(8)  :: SCLMFDFR
      real(8)  :: RSUB_RADIUS
   end type T_CLOUD_CTL

!#ifdef _CUDA
!
!   ! Inputs
!   ! ------
!
!   real(8), allocatable, dimension(:,:), device :: PP_dev
!   real(8), allocatable, dimension(:,:), device :: EXNP_dev
!   real(8), allocatable, dimension(:,:), device :: PPE_dev
!   real(8), allocatable, dimension(:,:), device :: KH_dev
!   real(8), allocatable, dimension(:  ), device :: FRLAND_dev
!   real(8), allocatable, dimension(:,:), device :: RMFDTR_dev
!   real(8), allocatable, dimension(:,:), device :: QLWDTR_dev
!   real(8), allocatable, dimension(:,:), device :: U_dev
!   real(8), allocatable, dimension(:,:), device :: V_dev
!   real(8), allocatable, dimension(:,:), device :: QST3_dev
!   real(8), allocatable, dimension(:,:), device :: DZET_dev
!   real(8), allocatable, dimension(:,:), device :: QDDF3_dev
!   real(8), allocatable, dimension(:  ), device :: TEMPOR_dev
!   real(8), allocatable, dimension(:  ), device :: CNV_FRACTION_dev
!
!   ! Inoutputs
!   ! ---------
!
!   real(8), allocatable, dimension(:,:), device :: TH_dev
!   real(8), allocatable, dimension(:,:), device :: Q_dev
!   real(8), allocatable, dimension(:,:), device :: QRN_CU_dev
!   real(8), allocatable, dimension(:,:), device :: CNV_UPDFRC_dev ! on edges, but dims=1:LM
!   real(8), allocatable, dimension(:,:), device :: QLW_LS_dev  
!   real(8), allocatable, dimension(:,:), device :: QLW_AN_dev
!   real(8), allocatable, dimension(:,:), device :: QIW_LS_dev  
!   real(8), allocatable, dimension(:,:), device :: QIW_AN_dev
!   real(8), allocatable, dimension(:,:), device :: ANVFRC_dev
!   real(8), allocatable, dimension(:,:), device :: CLDFRC_dev
!
!   ! Outputs
!   ! -------
!
!   real(8), allocatable, dimension(:,:), device :: RAD_CLDFRC_dev
!   real(8), allocatable, dimension(:,:), device :: RAD_QL_dev
!   real(8), allocatable, dimension(:,:), device :: RAD_QI_dev
!   real(8), allocatable, dimension(:,:), device :: RAD_QR_dev
!   real(8), allocatable, dimension(:,:), device :: RAD_QS_dev
!   real(8), allocatable, dimension(:,:), device :: QPLS_dev
!   real(8), allocatable, dimension(:,:), device :: CLDREFFL_dev
!   real(8), allocatable, dimension(:,:), device :: CLDREFFI_dev
!   real(8), allocatable, dimension(:  ), device :: PRELS_dev
!   real(8), allocatable, dimension(:  ), device :: PRECU_dev
!   real(8), allocatable, dimension(:  ), device :: PREAN_dev
!   real(8), allocatable, dimension(:  ), device :: LSARF_dev
!   real(8), allocatable, dimension(:  ), device :: CUARF_dev
!   real(8), allocatable, dimension(:  ), device :: ANARF_dev
!   real(8), allocatable, dimension(:  ), device :: SNRLS_dev
!   real(8), allocatable, dimension(:  ), device :: SNRCU_dev
!   real(8), allocatable, dimension(:  ), device :: SNRAN_dev
!
!   real(8), allocatable, dimension(:,:), device :: PFL_CN_dev
!   real(8), allocatable, dimension(:,:), device :: PFI_CN_dev
!   real(8), allocatable, dimension(:,:), device :: PFL_AN_dev
!   real(8), allocatable, dimension(:,:), device :: PFI_AN_dev
!   real(8), allocatable, dimension(:,:), device :: PFL_LS_dev
!   real(8), allocatable, dimension(:,:), device :: PFI_LS_dev
!
!   real(8), allocatable, dimension(:,:), device :: RHX_dev
!   real(8), allocatable, dimension(:,:), device :: REV_LS_dev
!   real(8), allocatable, dimension(:,:), device :: REV_AN_dev
!   real(8), allocatable, dimension(:,:), device :: REV_CN_dev
!   real(8), allocatable, dimension(:,:), device :: RSU_LS_dev
!   real(8), allocatable, dimension(:,:), device :: RSU_AN_dev
!   real(8), allocatable, dimension(:,:), device :: RSU_CN_dev
!   real(8), allocatable, dimension(:,:), device :: ACLL_CN_dev
!   real(8), allocatable, dimension(:,:), device :: ACIL_CN_dev
!   real(8), allocatable, dimension(:,:), device :: ACLL_AN_dev
!   real(8), allocatable, dimension(:,:), device :: ACIL_AN_dev
!   real(8), allocatable, dimension(:,:), device :: ACLL_LS_dev
!   real(8), allocatable, dimension(:,:), device :: ACIL_LS_dev
!   real(8), allocatable, dimension(:,:), device :: PDFL_dev
!   real(8), allocatable, dimension(:,:), device :: PDFI_dev
!   real(8), allocatable, dimension(:,:), device :: FIXL_dev
!   real(8), allocatable, dimension(:,:), device :: FIXI_dev                          
!   real(8), allocatable, dimension(:,:), device :: AUT_dev
!   real(8), allocatable, dimension(:,:), device :: EVAPC_dev
!   real(8), allocatable, dimension(:,:), device :: SDM_dev
!   real(8), allocatable, dimension(:,:), device :: SUBLC_dev 
!   real(8), allocatable, dimension(:,:), device :: FRZ_TT_dev
!   real(8), allocatable, dimension(:,:), device :: FRZ_PP_dev
!   real(8), allocatable, dimension(:,:), device :: DCNVL_dev
!   real(8), allocatable, dimension(:,:), device :: DCNVi_dev
!   real(8), allocatable, dimension(:,:), device :: ALPHT_dev
!   real(8), allocatable, dimension(:,:), device :: ALPH1_dev
!   real(8), allocatable, dimension(:,:), device :: ALPH2_dev
!   real(8), allocatable, dimension(:,:), device :: CFPDF_dev
!   real(8), allocatable, dimension(:,:), device :: RHCLR_dev
!   real(8), allocatable, dimension(:,:), device :: DQRL_dev
!   real(8), allocatable, dimension(:,:), device :: VFALLICE_AN_dev
!   real(8), allocatable, dimension(:,:), device :: VFALLICE_LS_dev
!   real(8), allocatable, dimension(:,:), device :: VFALLWAT_AN_dev
!   real(8), allocatable, dimension(:,:), device :: VFALLWAT_LS_dev
!   real(8), allocatable, dimension(:,:), device :: VFALLSN_AN_dev
!   real(8), allocatable, dimension(:,:), device :: VFALLSN_LS_dev
!   real(8), allocatable, dimension(:,:), device :: VFALLSN_CN_dev
!   real(8), allocatable, dimension(:,:), device :: VFALLRN_AN_dev
!   real(8), allocatable, dimension(:,:), device :: VFALLRN_LS_dev
!   real(8), allocatable, dimension(:,:), device :: VFALLRN_CN_dev
!
!   real(8), allocatable, dimension(:,:), device :: LIQANMOVE_dev
!   real(8), allocatable, dimension(:,:), device :: ICEANMOVE_dev
!   real(8), allocatable, dimension(:,:), device :: DANCLD_dev
!   real(8), allocatable, dimension(:,:), device :: DLSCLD_dev
!   real(8), allocatable, dimension(:,:), device :: CURAINMOVE_dev
!   real(8), allocatable, dimension(:,:), device :: CUSNOWMOVE_dev
!
!   ! Constants passed in from CLDPARAMS (from MAPL_GetResource)
!   real(8),    constant :: CNV_BETA
!   real(8),    constant :: ANV_BETA
!   real(8),    constant :: LS_BETA
!   real(8),    constant :: RH00
!   real(8),    constant :: C_00
!   real(8),    constant :: LWCRIT
!   real(8),    constant :: C_ACC
!   real(8),    constant :: C_EV_R
!   real(8),    constant :: C_EV_S
!   real(8),    constant :: CLDVOL2FRC
!   real(8),    constant :: RHSUP_ICE
!   real(8),    constant :: SHR_EVAP_FAC
!   real(8),    constant :: MIN_CLD_WATER
!   real(8),    constant :: CLD_EVP_EFF
!   integer, constant :: NSMAX
!   real(8),    constant :: LS_SDQV2
!   real(8),    constant :: LS_SDQV3
!   real(8),    constant :: LS_SDQVT1
!   real(8),    constant :: ANV_SDQV2
!   real(8),    constant :: ANV_SDQV3
!   real(8),    constant :: ANV_SDQVT1
!   real(8),    constant :: ANV_TO_LS
!   real(8),    constant :: N_WARM
!   real(8),    constant :: N_ICE
!   real(8),    constant :: N_ANVIL
!   real(8),    constant :: N_PBL
!   integer, constant :: DISABLE_RAD
!   integer, constant :: ICE_SETTLE
!   real(8),    constant :: ANV_ICEFALL_C
!   real(8),    constant :: LS_ICEFALL_C
!   real(8),    constant :: REVAP_OFF_P
!   real(8),    constant :: CNVENVFC
!   real(8),    constant :: WRHODEP
!   real(8),    constant :: T_ICE_ALL
!   real(8),    constant :: CNVICEPARAM
!   integer, constant :: ICEFRPWR
!   real(8),    constant :: CNVDDRFC
!   real(8),    constant :: ANVDDRFC
!   real(8),    constant :: LSDDRFC
!   integer, constant :: TANHRHCRIT
!   real(8),    constant :: MINRHCRIT
!   real(8),    constant :: MAXRHCRIT
!   real(8),    constant :: TURNRHCRIT
!   real(8),    constant :: MAXRHCRITLAND
!   integer, constant :: FR_LS_WAT
!   integer, constant :: FR_LS_ICE
!   integer, constant :: FR_AN_WAT
!   integer, constant :: FR_AN_ICE
!   real(8),    constant :: MIN_RL
!   real(8),    constant :: MIN_RI
!   real(8),    constant :: MAX_RL
!   real(8),    constant :: MAX_RI
!   real(8),    constant :: RI_ANV
!   integer, constant :: PDFFLAG
!
!   ! Parameters for Internal DQSAT
!   ! -----------------------------
!
!   real(8), parameter :: ESFAC            = MAPL8_H2OMW/MAPL8_AIRMW
!   real(8), parameter :: MAX_MIXING_RATIO = 1.
!   real(8), parameter :: ZEROC            = MAPL8_TICE
!
!   real(8), parameter :: TMINTBL   =  150.0
!   real(8), parameter :: TMAXTBL   =  333.0
!   real(8), parameter :: DEGSUBS   =  100
!   real(8), parameter :: ERFAC     = (DEGSUBS/ESFAC)
!   real(8), parameter :: DELTA_T   =  1.0 / DEGSUBS
!   real(8), parameter :: TABLESIZE =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1
!   real(8), parameter :: TMIX      = -20.
!
!   real(8), parameter :: TMINSTR = -95.
!   real(8), parameter :: TSTARR1 = -75.
!   real(8), parameter :: TSTARR2 = -65.
!   real(8), parameter :: TSTARR3 = -50.
!   real(8), parameter :: TSTARR4 = -40.
!   real(8), parameter :: TMAXSTR = +60.
!
!   real(kind=MAPL_R8), parameter :: B6 = 6.136820929E-11*100.0
!   real(kind=MAPL_R8), parameter :: B5 = 2.034080948E-8 *100.0
!   real(kind=MAPL_R8), parameter :: B4 = 3.031240396E-6 *100.0
!   real(kind=MAPL_R8), parameter :: B3 = 2.650648471E-4 *100.0
!   real(kind=MAPL_R8), parameter :: B2 = 1.428945805E-2 *100.0
!   real(kind=MAPL_R8), parameter :: B1 = 4.436518521E-1 *100.0
!   real(kind=MAPL_R8), parameter :: B0 = 6.107799961E+0 *100.0
!   real(kind=MAPL_R8), parameter :: BI6= 1.838826904E-10*100.0
!   real(kind=MAPL_R8), parameter :: BI5= 4.838803174E-8 *100.0
!   real(kind=MAPL_R8), parameter :: BI4= 5.824720280E-6 *100.0
!   real(kind=MAPL_R8), parameter :: BI3= 4.176223716E-4 *100.0
!   real(kind=MAPL_R8), parameter :: BI2= 1.886013408E-2 *100.0
!   real(kind=MAPL_R8), parameter :: BI1= 5.034698970E-1 *100.0
!   real(kind=MAPL_R8), parameter :: BI0= 6.109177956E+0 *100.0
!   real(kind=MAPL_R8), parameter :: S16= 0.516000335E-11*100.0
!   real(kind=MAPL_R8), parameter :: S15= 0.276961083E-8 *100.0
!   real(kind=MAPL_R8), parameter :: S14= 0.623439266E-6 *100.0
!   real(kind=MAPL_R8), parameter :: S13= 0.754129933E-4 *100.0
!   real(kind=MAPL_R8), parameter :: S12= 0.517609116E-2 *100.0
!   real(kind=MAPL_R8), parameter :: S11= 0.191372282E+0 *100.0
!   real(kind=MAPL_R8), parameter :: S10= 0.298152339E+1 *100.0
!   real(kind=MAPL_R8), parameter :: S26= 0.314296723E-10*100.0
!   real(kind=MAPL_R8), parameter :: S25= 0.132243858E-7 *100.0
!   real(kind=MAPL_R8), parameter :: S24= 0.236279781E-5 *100.0
!   real(kind=MAPL_R8), parameter :: S23= 0.230325039E-3 *100.0
!   real(kind=MAPL_R8), parameter :: S22= 0.129690326E-1 *100.0
!   real(kind=MAPL_R8), parameter :: S21= 0.401390832E+0 *100.0
!   real(kind=MAPL_R8), parameter :: S20= 0.535098336E+1 *100.0
!
!   real(kind=MAPL_R8), parameter :: DI(0:3) = (/ 57518.5606E08, 2.01889049, 3.56654, 20.947031 /)
!   real(kind=MAPL_R8), parameter :: CI(0:3) = (/ 9.550426, -5723.265, 3.53068, -.00728332 /)
!   real(kind=MAPL_R8), parameter :: DL(1:6) = (/ -7.902980, 5.02808, -1.3816, 11.344, 8.1328, -3.49149 /)
!   real(kind=MAPL_R8), parameter :: LOGPS   = 3.005714898  ! log10(1013.246)
!   real(kind=MAPL_R8), parameter :: TS      = 373.16
!   real(kind=MAPL_R8), parameter :: CL(0:9) = (/54.842763, -6763.22, -4.21000, .000367, &
!                                       .0415, 218.8,  53.878000, -1331.22, -9.44523, .014025  /)
!
!   real(8), parameter :: TMINLQU = MAPL8_TICE - 40.0
!   real(8), parameter :: TMINICE = MAPL8_TICE + -95.
!
!#else

!! Some parameters set by CLDPARAMS 

   real(8) :: CNV_BETA
   real(8) :: ANV_BETA
   real(8) :: LS_BETA
   real(8) :: RH00
   real(8) :: C_00
   real(8) :: LWCRIT
   real(8) :: C_ACC
   real(8) :: C_EV_R
   real(8) :: C_EV_S
   real(8) :: CLDVOL2FRC
   real(8) :: RHSUP_ICE
   real(8) :: SHR_EVAP_FAC
   real(8) :: MIN_CLD_WATER
   real(8) :: CLD_EVP_EFF
   integer :: NSMAX
   real(8) :: LS_SDQV2
   real(8) :: LS_SDQV3
   real(8) :: LS_SDQVT1
   real(8) :: ANV_SDQV2
   real(8) :: ANV_SDQV3
   real(8) :: ANV_SDQVT1
   real(8) :: ANV_TO_LS
   real(8) :: N_WARM
   real(8) :: N_ICE
   real(8) :: N_ANVIL
   real(8) :: N_PBL
   integer :: DISABLE_RAD
   integer :: ICE_SETTLE
   real(8) :: ANV_ICEFALL_C
   real(8) :: LS_ICEFALL_C
   real(8) :: REVAP_OFF_P
   real(8) :: CNVENVFC
   real(8) :: WRHODEP
   real(8) :: T_ICE_ALL
   real(8) :: CNVICEPARAM
   integer :: ICEFRPWR
   real(8) :: CNVDDRFC
   real(8) :: ANVDDRFC
   real(8) :: LSDDRFC
   integer :: tanhrhcrit
   real(8) :: minrhcrit
   real(8) :: maxrhcrit
   real(8) :: turnrhcrit
   real(8) :: MIN_RI, MAX_RI, MIN_RL, MAX_RL, RI_ANV
   integer :: FR_LS_WAT, FR_LS_ICE, FR_AN_WAT, FR_AN_ICE
   real(8) :: maxrhcritland
   integer :: pdfflag
   integer :: KTOP

!#endif

   real(8), parameter :: T_ICE_MAX    = MAPL8_TICE  ! -7.0+MAPL8_TICE
   real(8), parameter :: RHO_W        = 1.0e3      ! Density of liquid water in kg/m^3
   real(8), parameter :: MIN_CLD_FRAC = 1.0e-8

   ! There are two PIs in this routine: PI_0 and MAPL8_PI
   real(8), parameter :: PI_0 = 4.*atan(1.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains 

! GPU The GPU main routine call is smaller due to CUDA limit on
!     number of arguments permitted in call. Most inputs and outputs
!     are USE-associated in the GridComp

!#ifdef _CUDA
!   attributes(global) subroutine progno_cloud(IRUN,LM,DT,SCLMFDFR)
!#else
   subroutine progno_cloud( &
!!! first vars are (in) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IRUN, LM         , &
         DT               , &
!         LATS_dev         , &
         PP_dev           , &
         PPE_dev          , &
         EXNP_dev         , &
         FRLAND_dev       , &
         KH_dev           , &
         RMFDTR_dev       , &
         QLWDTR_dev       , &              
         QRN_CU_dev       , &
         CNV_UPDFRC_dev   , &
         U_dev            , &
         V_dev            , & 
         TH_dev           , &
         Q_dev            , &
         QLW_LS_dev       , &
         QLW_AN_dev       , &
         QIW_LS_dev       , &
         QIW_AN_dev       , &
         ANVFRC_dev       , &
         CLDFRC_dev       , &
!         RAD_CLDFRC_dev   , &
!         RAD_QL_dev       , &
!         RAD_QI_dev       , &
!         RAD_QR_dev       , &
!         RAD_QS_dev       , &
!         QPLS_dev         , &
!         CLDREFFL_dev     , &
!         CLDREFFI_dev     , &
!         PRELS_dev        , &
!         PRECU_dev        , &
!         PREAN_dev        , &
!         LSARF_dev        , &
!         CUARF_dev        , &
!         ANARF_dev        , &
!         SNRLS_dev        , &
!         SNRCU_dev        , &
!         SNRAN_dev        , &
         CLDPARAMS        , &
         SCLMFDFR         , &
         QST3_dev         , &
         DZET_dev         , &
         QDDF3_dev        , &
         CNV_FRACTION_dev )!, &
!         RHX_dev          , &
!         REV_LS_dev       , &
!         REV_AN_dev       , &
!         REV_CN_dev       , &
!         RSU_LS_dev       , &
!         RSU_AN_dev       , &
!         RSU_CN_dev       , &
!         ACLL_CN_dev,ACIL_CN_dev, &
!         ACLL_AN_dev,ACIL_AN_dev, &
!         ACLL_LS_dev,ACIL_LS_dev, &
!         PFL_CN_dev,PFI_CN_dev, &
!         PFL_AN_dev,PFI_AN_dev, &
!         PFL_LS_dev,PFI_LS_dev, &
!         PDFL_dev,PDFI_dev,FIXL_dev,FIXI_dev, &
!         AUT_dev, EVAPC_dev, SDM_dev, SUBLC_dev,    &
!         FRZ_TT_dev, DCNVL_dev, DCNVi_dev,      &
!         ALPHT_dev, ALPH1_dev, ALPH2_dev,       &
!         CFPDF_dev, RHCLR_dev,       &
!         DQRL_dev,FRZ_PP_dev,               &
!         VFALLICE_AN_dev,VFALLICE_LS_dev,   &
!         VFALLWAT_AN_dev,VFALLWAT_LS_dev,   &
!         VFALLSN_AN_dev,VFALLSN_LS_dev,VFALLSN_CN_dev,  &
!         VFALLRN_AN_dev,VFALLRN_LS_dev,VFALLRN_CN_dev,  &
!         LIQANMOVE_dev, ICEANMOVE_dev, &
!         DANCLD_dev, DLSCLD_dev, &
!         CURAINMOVE_dev, CUSNOWMOVE_dev,   &
!         TEMPOR_dev)
!#endif

      implicit none

!#ifdef _CUDA
!      integer, intent(in   ), value :: IRUN
!      integer, intent(in   ), value :: LM
!      real(8), intent(in   ), value :: DT
!      real(8), intent(in   ), value :: SCLMFDFR   ! CLOUD_CTL%SCLMFDFR
!#else
      type (CLDPARAM_TYPE), intent(in)          :: CLDPARAMS

      integer, intent(in   )                    :: IRUN ! IM*JM
      integer, intent(in   )                    :: LM   ! LM
      real(8), intent(in   )                       :: DT   ! DT_MOIST
!      real(8), intent(in   ), dimension(IRUN)      :: LATS_dev    ! LATS
      real(8), intent(in   ), dimension(IRUN,  LM) :: PP_dev      ! PLO
      real(8), intent(in   ), dimension(IRUN,0:LM) :: PPE_dev     ! CNV_PLE
      real(8), intent(in   ), dimension(IRUN,  LM) :: EXNP_dev    ! PK
      real(8), intent(in   ), dimension(IRUN     ) :: FRLAND_dev  ! FRLAND
      real(8), intent(in   ), dimension(IRUN,0:LM) :: KH_dev      ! KH
      real(8), intent(in   ), dimension(IRUN,  LM) :: RMFDTR_dev  ! CNV_MFD
      real(8), intent(in   ), dimension(IRUN,  LM) :: QLWDTR_dev  ! CNV_DQLDT
      real(8), intent(inout), dimension(IRUN,  LM) :: QRN_CU_dev  ! CNV_PRC3 IS THIS INTENT IN?
      real(8), intent(inout), dimension(IRUN,  LM) :: CNV_UPDFRC_dev ! CNV_UPDF
      real(8), intent(in   ), dimension(IRUN,  LM) :: U_dev  ! U1
      real(8), intent(in   ), dimension(IRUN,  LM) :: V_dev  ! V1
      real(8), intent(inout), dimension(IRUN,  LM) :: TH_dev ! TH1
      real(8), intent(inout), dimension(IRUN,  LM) :: Q_dev  ! Q1
      real(8), intent(inout), dimension(IRUN,  LM) :: QLW_LS_dev ! QLLS
      real(8), intent(inout), dimension(IRUN,  LM) :: QLW_AN_dev ! QLCN
      real(8), intent(inout), dimension(IRUN,  LM) :: QIW_LS_dev ! QILS
      real(8), intent(inout), dimension(IRUN,  LM) :: QIW_AN_dev ! QICN
      real(8), intent(inout), dimension(IRUN,  LM) :: ANVFRC_dev ! CLCN
      real(8), intent(inout), dimension(IRUN,  LM) :: CLDFRC_dev ! CLLS
!      real(8), intent(  out), dimension(IRUN,  LM) :: RAD_CLDFRC_dev ! RAD_CF
!      real(8), intent(  out), dimension(IRUN,  LM) :: RAD_QL_dev ! RAD_QL
!      real(8), intent(  out), dimension(IRUN,  LM) :: RAD_QI_dev ! RAD_QI
!      real(8), intent(  out), dimension(IRUN,  LM) :: RAD_QR_dev ! QRAIN
!      real(8), intent(  out), dimension(IRUN,  LM) :: RAD_QS_dev ! QSNOW
!      real(8), intent(  out), dimension(IRUN,  LM) :: QPLS_dev ! QPLS
!      real(8), intent(  out), dimension(IRUN,  LM) :: CLDREFFL_dev ! CLDREFFL
!      real(8), intent(  out), dimension(IRUN,  LM) :: CLDREFFI_dev ! CLDREFFI
!      real(8), intent(  out), dimension(IRUN     ) :: PRELS_dev ! LS_PRC2
!      real(8), intent(  out), dimension(IRUN     ) :: PRECU_dev ! CN_PRC2
!      real(8), intent(  out), dimension(IRUN     ) :: PREAN_dev ! AN_PRC2
!      real(8), intent(  out), dimension(IRUN     ) :: LSARF_dev ! LS_ARFX
!      real(8), intent(  out), dimension(IRUN     ) :: CUARF_dev ! CN_ARFX
!      real(8), intent(  out), dimension(IRUN     ) :: ANARF_dev ! AN_ARFX
!      real(8), intent(  out), dimension(IRUN     ) :: SNRLS_dev ! LS_SNR
!      real(8), intent(  out), dimension(IRUN     ) :: SNRCU_dev ! CN_SNR
!      real(8), intent(  out), dimension(IRUN     ) :: SNRAN_dev ! AN_SNR
      real(8), intent(in   )                       :: SCLMFDFR   ! CLOUD_CTL%SCLMFDFR
      real(8), intent(in   ), dimension(IRUN,  LM) :: QST3_dev   ! QST3
      real(8), intent(in   ), dimension(IRUN,  LM) :: DZET_dev   ! DZET
      real(8), intent(in   ), dimension(IRUN,  LM) :: QDDF3_dev  ! QDDF3
      real(8), intent(in   ), dimension(IRUN)      :: CNV_FRACTION_dev   ! CNV_FRACTION

!      real(8), intent(  out), dimension(IRUN,  LM) :: RHX_dev    ! RHX
!      real(8), intent(  out), dimension(IRUN,  LM) :: REV_LS_dev ! REV_LS
!      real(8), intent(  out), dimension(IRUN,  LM) :: REV_AN_dev ! REV_AN
!      real(8), intent(  out), dimension(IRUN,  LM) :: REV_CN_dev ! REV_CN
!      real(8), intent(  out), dimension(IRUN,  LM) :: RSU_LS_dev ! RSU_LS
!      real(8), intent(  out), dimension(IRUN,  LM) :: RSU_AN_dev ! RSU_AN
!      real(8), intent(  out), dimension(IRUN,  LM) :: RSU_CN_dev ! RSU_CN
!      real(8), intent(  out), dimension(IRUN,  LM) :: ACLL_CN_dev ! ACLL_CN
!      real(8), intent(  out), dimension(IRUN,  LM) :: ACIL_CN_dev ! ACIL_CN
!      real(8), intent(  out), dimension(IRUN,  LM) :: ACLL_AN_dev ! ACLL_AN
!      real(8), intent(  out), dimension(IRUN,  LM) :: ACIL_AN_dev ! ACIL_AN
!      real(8), intent(  out), dimension(IRUN,  LM) :: ACLL_LS_dev ! ACLL_LS
!      real(8), intent(  out), dimension(IRUN,  LM) :: ACIL_LS_dev ! ACIL_LS
!      real(8), intent(  out), dimension(IRUN,0:LM) :: PFL_CN_dev ! PFL_CN
!      real(8), intent(  out), dimension(IRUN,0:LM) :: PFI_CN_dev ! PFI_CN
!      real(8), intent(  out), dimension(IRUN,0:LM) :: PFL_AN_dev ! PFL_AN
!      real(8), intent(  out), dimension(IRUN,0:LM) :: PFI_AN_dev ! PFI_AN
!      real(8), intent(  out), dimension(IRUN,0:LM) :: PFL_LS_dev ! PFL_LS
!      real(8), intent(  out), dimension(IRUN,0:LM) :: PFI_LS_dev ! PFI_LS
!      real(8), intent(  out), dimension(IRUN,  LM) :: PDFL_dev ! DlPDF
!      real(8), intent(  out), dimension(IRUN,  LM) :: PDFI_dev ! DiPDF
!      real(8), intent(  out), dimension(IRUN,  LM) :: FIXL_dev ! DlFIX
!      real(8), intent(  out), dimension(IRUN,  LM) :: FIXI_dev ! DiFIX                  
!      real(8), intent(  out), dimension(IRUN,  LM) :: AUT_dev   ! AUT
!      real(8), intent(  out), dimension(IRUN,  LM) :: EVAPC_dev ! EVAPC
!      real(8), intent(  out), dimension(IRUN,  LM) :: SDM_dev   ! SDM
!      real(8), intent(  out), dimension(IRUN,  LM) :: SUBLC_dev ! SUBLC
!      real(8), intent(  out), dimension(IRUN,  LM) :: FRZ_TT_dev ! FRZ_TT
!      real(8), intent(  out), dimension(IRUN,  LM) :: FRZ_PP_dev ! FRZ_PP
!      real(8), intent(  out), dimension(IRUN,  LM) :: DCNVL_dev ! DCNVL
!      real(8), intent(  out), dimension(IRUN,  LM) :: DCNVi_dev ! DCNVi
!      real(8), intent(  out), dimension(IRUN,  LM) :: ALPHT_dev ! ALPHT
!      real(8), intent(  out), dimension(IRUN,  LM) :: ALPH1_dev ! ALPH1
!      real(8), intent(  out), dimension(IRUN,  LM) :: ALPH2_dev ! ALPH2
!      real(8), intent(  out), dimension(IRUN,  LM) :: CFPDF_dev ! CFPDF
!      real(8), intent(  out), dimension(IRUN,  LM) :: RHCLR_dev ! RHCLR
!      real(8), intent(  out), dimension(IRUN,  LM) :: DQRL_dev ! DQRL
!      real(8), intent(  out), dimension(IRUN,  LM) :: VFALLICE_AN_dev ! VFALLICE_AN
!      real(8), intent(  out), dimension(IRUN,  LM) :: VFALLICE_LS_dev ! VFALLICE_LS
!      real(8), intent(  out), dimension(IRUN,  LM) :: VFALLWAT_AN_dev ! VFALLWAT_AN
!      real(8), intent(  out), dimension(IRUN,  LM) :: VFALLWAT_LS_dev ! VFALLWAT_LS
!      real(8), intent(  out), dimension(IRUN,  LM) :: VFALLSN_AN_dev ! VFALLSN_AN
!      real(8), intent(  out), dimension(IRUN,  LM) :: VFALLSN_LS_dev ! VFALLSN_LS
!      real(8), intent(  out), dimension(IRUN,  LM) :: VFALLSN_CN_dev ! VFALLSN_CN
!      real(8), intent(  out), dimension(IRUN,  LM) :: VFALLRN_AN_dev ! VFALLRN_AN
!      real(8), intent(  out), dimension(IRUN,  LM) :: VFALLRN_LS_dev ! VFALLRN_LS
!      real(8), intent(  out), dimension(IRUN,  LM) :: VFALLRN_CN_dev ! VFALLRN_CN

!      real(8), intent(  out), dimension(IRUN,  LM) :: LIQANMOVE_dev  ! LIQANMOVE
!      real(8), intent(  out), dimension(IRUN,  LM) :: ICEANMOVE_dev  ! ICEANMOVE
!      real(8), intent(  out), dimension(IRUN,  LM) :: DANCLD_dev     ! DANCLD
!      real(8), intent(  out), dimension(IRUN,  LM) :: DLSCLD_dev     ! DLSCLD
!      real(8), intent(  out), dimension(IRUN,  LM) :: CURAINMOVE_dev ! CURAINMOVE
!      real(8), intent(  out), dimension(IRUN,  LM) :: CUSNOWMOVE_dev ! CUSNOWMOVE

!      real(8), intent(in   ), dimension(IRUN     ) :: TEMPOR_dev  ! TEMPOR

!#endif

! GPU The GPUs need to know how big local arrays are during compile-time
!     as the GPUs cannot allocate memory themselves. This command resets
!     this a priori size to LM for the CPU.
!#ifndef GPU_MAXLEVS
!#define GPU_MAXLEVS LM
!#endif

      integer :: I , J , K , L

      integer :: FRACTION_REMOVAL1, FRACTION_REMOVAL2, FRACTION_REMOVAL3

      real(8) :: MASS, iMASS
      real(8) :: TOTFRC
      real(8) :: QRN_LS, QRN_AN, QRN_CU_1D
      real(8) :: QSN_LS, QSN_AN, QSN_CU
!      real(8) :: QRN_ALL, QSN_ALL
      real(8) :: QTMP1, QTMP2, QTMP3
      real(8) :: TEMP
      real(8) :: RHCRIT
      real(8) :: AA3, BB3, ALPHA
      real(8) :: VFALL, VFALLRN, VFALLSN

      real(8) :: TOT_PREC_UPD
      real(8) :: TOT_PREC_ANV
      real(8) :: TOT_PREC_LS
      real(8) :: AREA_UPD_PRC
      real(8) :: AREA_ANV_PRC
      real(8) :: AREA_LS_PRC

      real(8) :: AREA_UPD_PRC_tolayer
      real(8) :: AREA_ANV_PRC_tolayer
      real(8) :: AREA_LS_PRC_tolayer

      real(8) :: U_above,U_below
      real(8) :: V_above,V_below
      real(8) :: DZET_above,DZET_below

      real(8) :: PRN_CU_above, PSN_CU_above
      real(8) :: PRN_LS_above, PSN_LS_above
      real(8) :: PRN_AN_above, PSN_AN_above

      real(8) :: EVAP_DD_CU_above, SUBL_DD_CU_above
      real(8) :: EVAP_DD_LS_above, SUBL_DD_LS_above
      real(8) :: EVAP_DD_AN_above, SUBL_DD_AN_above

      logical :: use_autoconv_timescale

      real(8) :: TROPICAL, EXTRATROPICAL

      real(8) :: LSENVFC, ANVENVFC, ENVFC, DDRFC, BETA

      real(8) :: SDQV2, SDQV3, SDQVT1

      real(8) :: LSPDFLIQNEW, LSPDFICENEW, LSPDFFRACNEW

! These are in constant memory in CUDA and are set in the GridComp
!#ifndef _CUDA
         CNV_BETA      = CLDPARAMS%CNV_BETA  ! Area factor for convective rain showers (non-dim)
         ANV_BETA      = CLDPARAMS%ANV_BETA  ! Area factor for anvil rain showers (non-dim)
         LS_BETA       = CLDPARAMS%LS_BETA   ! Area factor for Large Scale rain showers (non-dim)
         RH00          = CLDPARAMS%RH_CRIT   ! Critical relative humidity
         C_00          = CLDPARAMS%AUTOC_LS
         LWCRIT        = CLDPARAMS%QC_CRIT_LS
         C_ACC         = CLDPARAMS%ACCRETION
         C_EV_R        = CLDPARAMS%RAIN_REVAP_FAC
         C_EV_S        = CLDPARAMS%SNOW_REVAP_FAC
         CLDVOL2FRC    = CLDPARAMS%VOL_TO_FRAC
         RHSUP_ICE     = CLDPARAMS%SUPERSAT
         SHR_EVAP_FAC  = CLDPARAMS%SHEAR_EVAP_FAC
         MIN_CLD_WATER = CLDPARAMS%MIN_ALLOW_CCW
         CLD_EVP_EFF   = CLDPARAMS%CCW_EVAP_EFF
         NSMAX         = INT( CLDPARAMS%NSUB_AUTOCONV  )
         LS_SDQV2      = CLDPARAMS%LS_SUND_INTER
         LS_SDQV3      = CLDPARAMS%LS_SUND_COLD
         LS_SDQVT1     = CLDPARAMS%LS_SUND_TEMP1
         ANV_SDQV2     = CLDPARAMS%ANV_SUND_INTER
         ANV_SDQV3     = CLDPARAMS%ANV_SUND_COLD
         ANV_SDQVT1    = CLDPARAMS%ANV_SUND_TEMP1
         ANV_TO_LS     = CLDPARAMS%ANV_TO_LS_TIME
         N_WARM        = CLDPARAMS%NCCN_WARM
         N_ICE         = CLDPARAMS%NCCN_ICE
         N_ANVIL       = CLDPARAMS%NCCN_ANVIL
         N_PBL         = CLDPARAMS%NCCN_PBL
         DISABLE_RAD   = INT( CLDPARAMS%DISABLE_RAD )
         ICE_SETTLE    = NINT( CLDPARAMS%ICE_SETTLE )
         ANV_ICEFALL_C = CLDPARAMS%ANV_ICEFALL
         LS_ICEFALL_C  = CLDPARAMS%LS_ICEFALL
         REVAP_OFF_P   = CLDPARAMS%REVAP_OFF_P
         CNVENVFC      = CLDPARAMS%CNV_ENVF
         WRHODEP       = CLDPARAMS%WRHODEP
         T_ICE_ALL     = CLDPARAMS%ICE_RAMP + MAPL8_TICE
         CNVICEPARAM   = CLDPARAMS%CNV_ICEPARAM
         ICEFRPWR      = INT( CLDPARAMS%CNV_ICEFRPWR + .001 )
         CNVDDRFC      = CLDPARAMS%CNV_DDRF
         ANVDDRFC      = CLDPARAMS%ANV_DDRF
         LSDDRFC       = CLDPARAMS%LS_DDRF
         TANHRHCRIT    = INT( CLDPARAMS%TANHRHCRIT )
         MINRHCRIT     = CLDPARAMS%MINRHCRIT
         MAXRHCRIT     = CLDPARAMS%MAXRHCRIT
         TURNRHCRIT    = CLDPARAMS%TURNRHCRIT
         MAXRHCRITLAND = CLDPARAMS%MAXRHCRITLAND
         FR_LS_WAT     = INT( CLDPARAMS%FR_LS_WAT )
         FR_LS_ICE     = INT( CLDPARAMS%FR_LS_ICE )
         FR_AN_WAT     = INT( CLDPARAMS%FR_AN_WAT )
         FR_AN_ICE     = INT( CLDPARAMS%FR_AN_ICE )
         MIN_RL        = CLDPARAMS%MIN_RL
         MIN_RI        = CLDPARAMS%MIN_RI
         MAX_RL        = CLDPARAMS%MAX_RL
         MAX_RI        = CLDPARAMS%MAX_RI
         RI_ANV        = CLDPARAMS%RI_ANV
         PDFFLAG       = INT(CLDPARAMS%PDFSHAPE)
         KTOP          = INT(CLDPARAMS%KTOP)
!#endif

      use_autoconv_timescale = .false.

!#ifdef _CUDA
!      i = (blockidx%x - 1) * blockdim%x + threadidx%x
!
!      RUN_LOOP: IF ( I <= IRUN ) THEN
!#else
      RUN_LOOP: DO I = 1, IRUN
!#endif
         K_LOOP: DO K = KTOP, LM         
            if (K == KTOP) then
               TOT_PREC_UPD = 0.
               TOT_PREC_ANV = 0.
               TOT_PREC_LS  = 0.

               AREA_UPD_PRC = 0.
               AREA_ANV_PRC = 0.
               AREA_LS_PRC  = 0.
            end if

            if (K == LM ) then
               !! ZERO DIAGNOSTIC OUTPUTS BEFORE SHOWERS !!
               
!               PRELS_dev(I) = 0.
!               PRECU_dev(I) = 0.
!               PREAN_dev(I) = 0.

!               SNRCU_dev(I) = 0. 
!               SNRLS_dev(I) = 0. 
!               SNRAN_dev(I) = 0. 

!               LSARF_dev(I) = 0.
!               CUARF_dev(I) = 0.
!               ANARF_dev(I) = 0.
            end if

            !Zero out/initialize precips, except QRN_CU which comes from RAS 
            QRN_LS    = 0.
            QRN_AN    = 0.
            QRN_CU_1D = 0.
            QSN_LS    = 0.
            QSN_AN    = 0.
            QSN_CU    = 0.
            VFALL     = 0.

!            RAD_QL_dev(I,K)     = 0.
!            RAD_QI_dev(I,K)     = 0.
!            RAD_QR_dev(I,K)     = 0.
!            RAD_QS_dev(I,K)     = 0.
!            QPLS_dev(I,K)       = 0.
!            RAD_CLDFRC_dev(I,K) = 0.
!            CLDREFFL_dev(I,K)   = 0.
!            CLDREFFI_dev(I,K)   = 0.

!            PFL_CN_dev(I,K) = 0.
!            PFI_CN_dev(I,K) = 0.
!            PFL_AN_dev(I,K) = 0.
!            PFI_AN_dev(I,K) = 0.
!            PFL_LS_dev(I,K) = 0.
!            PFI_LS_dev(I,K) = 0.

!            IF (K == KTOP) THEN
!               PFL_CN_dev(I,0) = 0.
!               PFI_CN_dev(I,0) = 0.
!               PFL_AN_dev(I,0) = 0.
!               PFI_AN_dev(I,0) = 0.
!               PFL_LS_dev(I,0) = 0.
!               PFI_LS_dev(I,0) = 0.
!            END IF

            ! Initialize other diagnostics 

!            RHX_dev(I,K) = MAPL8_UNDEF
!            REV_LS_dev(I,K) = MAPL8_UNDEF
!            REV_AN_dev(I,K) = MAPL8_UNDEF
!            REV_CN_dev(I,K) = MAPL8_UNDEF
!            RSU_LS_dev(I,K) = MAPL8_UNDEF
!            RSU_AN_dev(I,K) = MAPL8_UNDEF
!            RSU_CN_dev(I,K) = MAPL8_UNDEF
!            ACLL_CN_dev(I,K) = MAPL8_UNDEF
!            ACIL_CN_dev(I,K) = MAPL8_UNDEF
!            ACLL_AN_dev(I,K) = MAPL8_UNDEF
!            ACIL_AN_dev(I,K) = MAPL8_UNDEF
!            ACLL_LS_dev(I,K) = MAPL8_UNDEF
!            ACIL_LS_dev(I,K) = MAPL8_UNDEF
!            PDFL_dev(I,K) = MAPL8_UNDEF
!            PDFI_dev(I,K) = MAPL8_UNDEF
!            FIXL_dev(I,K) = MAPL8_UNDEF
!            FIXI_dev(I,K) = MAPL8_UNDEF
!            AUT_dev(I,K) = MAPL8_UNDEF
!            EVAPC_dev(I,K) = MAPL8_UNDEF
!            SDM_dev(I,K) = MAPL8_UNDEF
!            SUBLC_dev(I,K) = MAPL8_UNDEF
!            FRZ_TT_dev(I,K) = MAPL8_UNDEF
!            FRZ_PP_dev(I,K) = MAPL8_UNDEF
!            DCNVL_dev(I,K) = MAPL8_UNDEF
!            DCNVi_dev(I,K) = MAPL8_UNDEF
!            ALPHT_dev(I,K) = MAPL8_UNDEF
!            ALPH1_dev(I,K) = MAPL8_UNDEF
!            ALPH2_dev(I,K) = MAPL8_UNDEF
!            CFPDF_dev(I,K) = MAPL8_UNDEF
!            RHCLR_dev(I,K) = MAPL8_UNDEF
!            DQRL_dev(I,K) = MAPL8_UNDEF
!            VFALLICE_AN_dev(I,K) = MAPL8_UNDEF
!            VFALLICE_LS_dev(I,K) = MAPL8_UNDEF
!            VFALLWAT_AN_dev(I,K) = MAPL8_UNDEF
!            VFALLWAT_LS_dev(I,K) = MAPL8_UNDEF
!            VFALLSN_AN_dev(I,K) = MAPL8_UNDEF
!            VFALLSN_LS_dev(I,K) = MAPL8_UNDEF
!            VFALLSN_CN_dev(I,K) = MAPL8_UNDEF
!            VFALLRN_AN_dev(I,K) = MAPL8_UNDEF
!            VFALLRN_LS_dev(I,K) = MAPL8_UNDEF
!            VFALLRN_CN_dev(I,K) = MAPL8_UNDEF

            ! Copy QRN_CU into a temp scalar
            QRN_CU_1D = QRN_CU_dev(I,K)

            MASS =  ( PPE_dev(I,K) - PPE_dev(I,K-1) )*100./MAPL8_GRAV  ! layer-mass (kg/m**2)

            iMASS = 1.0 / MASS

            TEMP =  EXNP_dev(I,K) * TH_dev(I,K) 

!            FRZ_PP_dev(I,K) = 0.00

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Total Condensate Source
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!            FIXL_dev(I,K) = QLW_AN_dev(I,K) + QLW_LS_dev(I,K)
!            FIXI_dev(I,K) = QIW_AN_dev(I,K) + QIW_LS_dev(I,K)

            CALL fix_up_clouds(    &
                  Q_dev(I,K)     , &
                  TEMP           , &
                  QLW_LS_dev(I,K), & 
                  QIW_LS_dev(I,K), & 
                  CLDFRC_dev(I,K), &   
                  QLW_AN_dev(I,K), & 
                  QIW_AN_dev(I,K), & 
                  ANVFRC_dev(I,K))

!            FIXL_dev(I,K) = -( QLW_AN_dev(I,K) + QLW_LS_dev(I,K) - FIXL_dev(I,K) ) / DT 
!            FIXI_dev(I,K) = -( QIW_AN_dev(I,K) + QIW_LS_dev(I,K) - FIXI_dev(I,K) ) / DT 
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!            FRZ_TT_dev(I,K) = QIW_AN_dev(I,K) + QIW_LS_dev(I,K)

            CALL meltfrz (         &
                  DT             , &
                  TEMP           , &
                  QLW_LS_dev(I,K), & 
                  QIW_LS_dev(I,K))

            CALL meltfrz (         &
                  DT             , &
                  TEMP           , &
                  QLW_AN_dev(I,K), & 
                  QIW_AN_dev(I,K))

!            FRZ_TT_dev(I,K) = ( QIW_AN_dev(I,K) + QIW_LS_dev(I,K) - FRZ_TT_dev(I,K) ) / DT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!            DCNVi_dev(I,K) = QIW_AN_dev(I,K)
!            DCNVL_dev(I,K) = QLW_AN_dev(I,K)

            CALL cnvsrc (          &  
                  DT             , &
                  CNVICEPARAM    , &
                  SCLMFDFR       , &
                  MASS           , & 
                  iMASS          , &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  Q_dev(I,K)     , &
                  QLWDTR_dev(I,K), &
                  RMFDTR_dev(I,K), &
                  QLW_AN_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  CLDFRC_dev(I,K), & 
                  ANVFRC_dev(I,K), &
                  QST3_dev(I,K)  )

!            DCNVi_dev(I,K) = ( QIW_AN_dev(I,K) - DCNVi_dev(I,K) ) / DT
!            DCNVL_dev(I,K) = ( QLW_AN_dev(I,K) - DCNVL_dev(I,K) ) / DT
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!            PDFL_dev(I,K) = QLW_LS_dev(I,K)+QLW_AN_dev(I,K)
!            PDFI_dev(I,K) = QIW_LS_dev(I,K)+QIW_AN_dev(I,K)

            if (k == KTOP .or. k == lm) then
               U_above = 0.0
               U_below = 0.0
               V_above = 0.0
               V_below = 0.0
               DZET_above = 0.0
               DZET_below = 0.0
            else
               U_above = U_dev(i,k-1)
               U_below = U_dev(i,k+1)
               V_above = V_dev(i,k-1)
               V_below = V_dev(i,k+1)
               DZET_above = DZET_dev(i,k-1)
               DZET_below = DZET_dev(i,k+1)
            end if
  
            call pdf_spread (&
                  K,LM,&
                  U_dev(I,K),U_above,U_below,&
                  V_dev(I,K),V_above,V_below,&
                  KH_dev(I,K-1),DZET_above,DZET_below,&
                  CNV_UPDFRC_dev(I,K),PP_dev(I,K),ALPHA,&
!                  ALPHT_dev(I,K), ALPH1_dev(I,K),!ALPH2_dev(I,K), & 
                  FRLAND_dev(I) )  

            ! impose a minimum amount of variability
            ALPHA    = MAX(  ALPHA , 1.0 - RH00 )
   
            RHCRIT = 1.0 - ALPHA

            LSPDFLIQNEW = QLW_LS_dev(I,K)
            LSPDFICENEW = QIW_LS_dev(I,K)
            LSPDFFRACNEW= CLDFRC_dev(I,K)

            call hystpdf(          &
                  DT             , &
                  ALPHA          , &
                  PDFFLAG        , &
                  PP_dev(I,K)    , &
                  Q_dev(I,K)     , &
                  QLW_LS_dev(I,K), &
                  QLW_AN_dev(I,K), &
                  QIW_LS_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  TEMP           , &
                  CLDFRC_dev(I,K), & 
                  ANVFRC_dev(I,K)  ) 

            LSPDFLIQNEW = QLW_LS_dev(I,K) - LSPDFLIQNEW
            LSPDFICENEW = QIW_LS_dev(I,K) - LSPDFICENEW
            LSPDFFRACNEW= CLDFRC_dev(I,K) - LSPDFFRACNEW

!            RHX_dev(I,K)   = Q_dev(I,K)/QSAT( TEMP, PP_dev(I,K) )
!            CFPDF_dev(I,K) = CLDFRC_dev(I,K)
!            PDFL_dev(I,K)  = ( QLW_LS_dev(I,K) + QLW_AN_dev(I,K) - PDFL_dev(I,K) ) / DT 
!            PDFI_dev(I,K)  = ( QIW_LS_dev(I,K) + QIW_AN_dev(I,K) - PDFI_dev(I,K) ) / DT 

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            TOTFRC = CLDFRC_dev(I,K) + ANVFRC_dev(I,K)

            IF ( TOTFRC > 1.00 ) THEN
               CLDFRC_dev(I,k) = CLDFRC_dev(I,k)*(1.00 / TOTFRC )
               ANVFRC_dev(I,k) = ANVFRC_dev(I,k)*(1.00 / TOTFRC )
            END IF

            TOTFRC = CLDFRC_dev(I,K) + ANVFRC_dev(I,K)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! CONDENSATE/FRACTION SOURCES FINISHED. NOW LOSE CLOUD CONDENSATE !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Total Condensate Sink
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!       E  V  A  P  O  R  A  T  I  O  N
            !!                A  N  D 
            !!       S  U  B  L  I  M  A  T  I  O  N
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!            EVAPC_dev(I,K) = QLW_LS_dev(I,K)+QLW_AN_dev(I,K)
!            SUBLC_dev(I,K) = QIW_LS_dev(I,K)+QIW_AN_dev(I,K)

             ! 'Anvil' partition from RAS/Parameterized not done in hystpdf

            call evap3(            &
                  DT             , &
                  RHCRIT         , &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  Q_dev(I,K)     , &
                  QLW_AN_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  ANVFRC_dev(I,K), &
                  CLDFRC_dev(I,K), &
                  QST3_dev(I,K)  )  

            call subl3(            &
                  DT             , & 
                  RHCRIT         , &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  Q_dev(I,K)     , &
                  QLW_AN_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  ANVFRC_dev(I,K), &
                  CLDFRC_dev(I,K), &
                  QST3_dev(I,K)  ) 

!            EVAPC_dev(I,K) = ( EVAPC_dev(I,K) - (QLW_LS_dev(I,K)+QLW_AN_dev(I,K)) ) / DT
!            SUBLC_dev(I,K) = ( SUBLC_dev(I,K) - (QIW_LS_dev(I,K)+QIW_AN_dev(I,K)) ) / DT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!       A U T O C O N V E R S I  O  N
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !           FRACTION_REMOVAL = 0 -> none
            !                              1 -> constant in-cloud QC
            !                              2 -> trim high edge of PDF
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            AUT_dev(I,K) = QLW_LS_dev(I,K)+QLW_AN_dev(I,K)

            FRACTION_REMOVAL1 = fr_ls_wat
            call autocon3(         &
                  DT             , &
                  QLW_LS_dev(I,K), &
                  QRN_LS         , &
                  TEMP           , &
                  PP_dev(I,K)    , &
                  KH_dev(I,K-1)  , &
                  CLDFRC_dev(I,K), &
                  LS_SDQV2       , &
                  LS_SDQV3       , &   
                  LS_SDQVT1      , &
                  DZET_dev(I,K)  , &
                  VFALL          , &
                  FRACTION_REMOVAL1 )
         
!            VFALLWAT_LS_dev(I,K) = VFALL

            FRACTION_REMOVAL2 = fr_an_wat 

            call autocon3(         &
                  DT             , &
                  QLW_AN_dev(I,K), &
                  QRN_AN         , &
                  TEMP           , &
                  PP_dev(I,K)    , &
                  KH_dev(I,K)    , &
                  ANVFRC_dev(I,K), &
                  ANV_SDQV2      , &
                  ANV_SDQV3      , &
                  ANV_SDQVT1     , &
                  DZET_dev(I,K)  , &
                  VFALL          , &
                  FRACTION_REMOVAL2 )

!            VFALLWAT_AN_dev(I,K) = VFALL
!            AUT_dev(I,K) = ( AUT_dev(I,K) - ( QLW_AN_dev(I,K) + QLW_LS_dev(I,K) ) )/DT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !  Ice cloud settling
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!            SDM_dev(I,K) = QIW_AN_dev(I,K)+QIW_LS_dev(I,K)

            ! Parameterized (RAS) Ice Fall
            ! ----------------------------

            ! WMP: Adjustments to resolved scale ice fall speed options                            
            if (CNV_FRACTION_dev(I) >= 0.5) then
               FRACTION_REMOVAL3 = fr_an_ice
            else
               FRACTION_REMOVAL3 = fr_ls_ice
            endif

            SELECT CASE( ICE_SETTLE )
            CASE( 0 )
             ! MERRA-2 Formulation
              TROPICAL      = ANV_ICEFALL_C*1.0
              EXTRATROPICAL = ANV_ICEFALL_C*0.0
            CASE( 1 )
              TROPICAL      = CNV_FRACTION_dev(I) 
              EXTRATROPICAL = 1.0-TROPICAL
              TROPICAL      = ANV_ICEFALL_C*TROPICAL
              EXTRATROPICAL =  LS_ICEFALL_C*EXTRATROPICAL
            END SELECT

            CALL SETTLE_VEL(       &
                  WRHODEP        , &
                  QIW_AN_dev(I,K), &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  ANVFRC_dev(I,K), &
                  KH_dev(I,K-1)  , &
                  VFALL          , &
                  EXTRATROPICAL, TROPICAL )
!            VFALLICE_AN_dev(I,K) = VFALL

            CALL ICEFALL(          &
                  QIW_AN_dev(I,K), &
                  DZET_dev(I,K)  , &
                  QSN_AN         , &
                  VFALL          , &
                  ANVFRC_dev(I,K), &
                  DT             , &
                  FRACTION_REMOVAL3 )

            ! Resolved Scale Ice Fall
            ! -----------------------

            ! WMP: Adjustments to resolved scale ice fall speed options 
            SELECT CASE( ICE_SETTLE )
            CASE( 0 )
             ! MERRA-2 Formulation
              TROPICAL      =  LS_ICEFALL_C*0.0
              EXTRATROPICAL =  LS_ICEFALL_C*1.0
            CASE( 1 )
              TROPICAL      = CNV_FRACTION_dev(I)
              EXTRATROPICAL = 1.0-TROPICAL
              TROPICAL      = ANV_ICEFALL_C*TROPICAL
              EXTRATROPICAL =  LS_ICEFALL_C*EXTRATROPICAL
            END SELECT

            CALL SETTLE_VEL(       &
                  WRHODEP        , &
                  QIW_LS_dev(I,K), &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  CLDFRC_dev(I,K), &
                  KH_dev(I,K-1)  , &
                  VFALL          , &
                  EXTRATROPICAL, TROPICAL )
!            VFALLICE_LS_dev(I,K) = VFALL

            CALL ICEFALL(          &
                  QIW_LS_dev(I,K), &
                  DZET_dev(I,K)  , &
                  QSN_LS         , &
                  VFALL          , &
                  CLDFRC_dev(I,K), &
                  DT             , &
                  FRACTION_REMOVAL3 )

!            SDM_dev(I,K) = ( SDM_dev(I,K) - (QIW_LS_dev(I,K) + QIW_AN_dev(I,K)) )/DT
      
!            DQRL_dev(I,K) = ( QRN_LS + QRN_AN + QSN_LS + QSN_AN ) / DT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !  Add in convective rain 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! CU-FREEZE 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Also "freeze" out any conv. precip that needs
            ! to be since this isnt done in RAS. This is
            ! precip w/ large particles, so freezing is 
            ! strict. Check up on this!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            QTMP2 = 0.

            if ( TEMP < MAPL8_TICE ) then
               QTMP2     = QRN_CU_1D
               QSN_CU    = QRN_CU_1D
               QRN_CU_1D = 0.
               TEMP      = TEMP + QSN_CU*(MAPL8_ALHS-MAPL8_ALHL) / MAPL8_CP
            end if
      
!            FRZ_PP_dev(I,K) = FRZ_PP_dev(I,K) +  QTMP2/DT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !----------------------------------------------------------------------------------------------
            ! Column will now be swept from top-down for precip accumulation/accretion/re-evaporation

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            AREA_LS_PRC_tolayer  = 0.0
            AREA_UPD_PRC_tolayer = 0.0
            AREA_ANV_PRC_tolayer = 0.0

            TOT_PREC_UPD  = TOT_PREC_UPD + ( ( QRN_CU_1D + QSN_CU ) * MASS )
            AREA_UPD_PRC  = AREA_UPD_PRC + ( CNV_UPDFRC_dev(I,K)* ( QRN_CU_1D + QSN_CU )* MASS )

            TOT_PREC_ANV  = TOT_PREC_ANV + ( ( QRN_AN + QSN_AN ) * MASS )
            AREA_ANV_PRC  = AREA_ANV_PRC + ( ANVFRC_dev(I,K)* ( QRN_AN + QSN_AN )* MASS )

            TOT_PREC_LS   = TOT_PREC_LS  + ( ( QRN_LS + QSN_LS ) * MASS )
            AREA_LS_PRC   = AREA_LS_PRC  + ( CLDFRC_dev(I,K)* ( QRN_LS + QSN_LS )* MASS )

            if ( TOT_PREC_ANV > 0.0 ) AREA_ANV_PRC_tolayer = MAX( AREA_ANV_PRC/TOT_PREC_ANV, 1.E-6 )
            if ( TOT_PREC_UPD > 0.0 ) AREA_UPD_PRC_tolayer = MAX( AREA_UPD_PRC/TOT_PREC_UPD, 1.E-6 )
            if ( TOT_PREC_LS  > 0.0 ) AREA_LS_PRC_tolayer  = MAX( AREA_LS_PRC/TOT_PREC_LS,   1.E-6 )

            AREA_LS_PRC_tolayer  = LS_BETA  * AREA_LS_PRC_tolayer
            AREA_UPD_PRC_tolayer = CNV_BETA * AREA_UPD_PRC_tolayer
            AREA_ANV_PRC_tolayer = ANV_BETA * AREA_ANV_PRC_tolayer

            IF (K == LM) THEN ! Weve accumulated over the whole column
               if ( TOT_PREC_ANV > 0.0 ) AREA_ANV_PRC = MAX( AREA_ANV_PRC/TOT_PREC_ANV, 1.E-6 )
               if ( TOT_PREC_UPD > 0.0 ) AREA_UPD_PRC = MAX( AREA_UPD_PRC/TOT_PREC_UPD, 1.E-6 )
               if ( TOT_PREC_LS  > 0.0 ) AREA_LS_PRC  = MAX( AREA_LS_PRC/TOT_PREC_LS,   1.E-6 )

               AREA_LS_PRC  = LS_BETA  * AREA_LS_PRC
               AREA_UPD_PRC = CNV_BETA * AREA_UPD_PRC
               AREA_ANV_PRC = ANV_BETA * AREA_ANV_PRC

               !! "couple" to diagnostic areal(8) fraction output 
               !! Intensity factor in PRECIP3 is floored at
               !! 1.0. So this is fair.

!               LSARF_dev(I) = MIN( AREA_LS_PRC,  1.0 )
!               CUARF_dev(I) = MIN( AREA_UPD_PRC, 1.0 )
!               ANARF_dev(I) = MIN( AREA_ANV_PRC, 1.0 )
            END IF

!            QRN_ALL = 0.
!            QSN_ALL = 0.

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! GET SOME MICROPHYSICAL QUANTITIES 

            CALL MICRO_AA_BB_3( TEMP,PP_dev(I,K),QST3_dev(I,K),AA3,BB3 )

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            QTMP1 = QLW_LS_dev(I,K) + QLW_AN_dev(I,K)
            QTMP2 = QIW_LS_dev(I,K) + QIW_AN_dev(I,K)

          ! CNVENVFC passed in as argument (default=0.8)
            ANVENVFC = 1.00 ! Unlike CNVENVFC which is from PHYSPARAMS, this is set to 1.
            LSENVFC  = 1.00 ! Unlike CNVENVFC which is from PHYSPARAMS, this is set to 1.
          ! Apply convective fractions
            ENVFC  = CNVENVFC*CNV_FRACTION_dev(I) + LSENVFC*(1.0-CNV_FRACTION_dev(I))

            ! Convective
            ! ----------

            call  PRECIP3(          &
                  K,LM            , &
                  DT              , &
                  FRLAND_dev(I)   , &
                  RHCRIT          , &
                  QRN_CU_1D       , &
                  QSN_CU          , &
                  QTMP1           , &
                  QTMP2           , &
                  TEMP            , &
                  Q_dev(I,K)      , &
                  mass            , &
                  imass           , &
                  PP_dev(I,K)     , &
                  DZET_dev(I,K)   , &
                  QDDF3_dev(I,K)  , &
                  AA3             , &
                  BB3             , &
                  AREA_UPD_PRC_tolayer    , &
!                  PRECU_dev(I)    , & 
!                  SNRCU_dev(I)    , & 
                  PRN_CU_above    , &
                  PSN_CU_above    , &
                  EVAP_DD_CU_above, &
                  SUBL_DD_CU_above, &
!                  REV_CN_dev(I,K) , &
!                  RSU_CN_dev(I,K) , &
!                  ACLL_CN_dev(I,K), &
!                  ACIL_CN_dev(I,K), &
!                  PFL_CN_dev(I,K) , &
!                  PFI_CN_dev(I,K) , &
!                  VFALLRN         , &
!                  VFALLSN         , &
!                  FRZ_PP_dev(I,K) , &
                  ENVFC, CNVDDRFC )

!            VFALLSN_CN_dev(I,K) = VFALLSN
!            VFALLRN_CN_dev(I,K) = VFALLRN

!            if (.not.use_autoconv_timescale) then
!               if (VFALLSN.NE.0.) then
!                  QSN_ALL = QSN_ALL + PFI_CN_dev(I,K)/VFALLSN
!               end if
!               if (VFALLRN.NE.0.) then
!                  QRN_ALL = QRN_ALL + PFL_CN_dev(I,K)/VFALLRN
!               end if
!            end if

            ! Anvil
            ! -----

            call  PRECIP3(          &
                  K,LM            , &
                  DT              , & 
                  FRLAND_dev(I)   , &
                  RHCRIT          , &
                  QRN_AN          , &
                  QSN_AN          , &
                  QTMP1           , &
                  QTMP2           , &
                  TEMP            , &
                  Q_dev(I,K)      , &
                  mass            , & 
                  imass           , &
                  PP_dev(I,K)     , &
                  DZET_dev(I,K)   , &
                  QDDF3_dev(I,K)  , &
                  AA3             , &
                  BB3             , &
                  AREA_ANV_PRC_tolayer    , &
!                  PREAN_dev(I)    , & 
!                  SNRAN_dev(I)    , &
                  PRN_AN_above    , &
                  PSN_AN_above    , &
                  EVAP_DD_AN_above, &
                  SUBL_DD_AN_above, &
!                  REV_AN_dev(I,K) , &
!                  RSU_AN_dev(I,K) , &
!                  ACLL_AN_dev(I,K), &
!                  ACIL_AN_dev(I,K), &
!                  PFL_AN_dev(I,K) , &
!                  PFI_AN_dev(I,K) , &
!                  VFALLRN         , &
!                  VFALLSN         , &
!                  FRZ_PP_dev(I,K) , &
                  ENVFC, ANVDDRFC )

!            VFALLSN_AN_dev(I,K) = VFALLSN
!            VFALLRN_AN_dev(I,K) = VFALLRN

!            if (.not.use_autoconv_timescale) then
!               if (VFALLSN.NE.0.) then
!                  QSN_ALL = QSN_ALL + PFI_AN_dev(I,K)/VFALLSN
!               end if
!               if (VFALLRN.NE.0.) then
!                  QRN_ALL = QRN_ALL + PFL_AN_dev(I,K)/VFALLRN
!               end if
!            end if

            ! Largescale
            ! ----------

            call  PRECIP3(          &
                  K,LM            , &
                  DT              , &  
                  FRLAND_dev(I)   , &
                  RHCRIT          , &
                  QRN_LS          , &
                  QSN_LS          , &
                  QTMP1           , &
                  QTMP2           , &
                  TEMP            , &
                  Q_dev(I,K)      , &
                  mass            , &
                  imass           , &
                  PP_dev(I,K)     , &
                  DZET_dev(I,K)   , &
                  QDDF3_dev(I,K)  , &
                  AA3             , &
                  BB3             , &
                  AREA_LS_PRC_tolayer     , &
!                  PRELS_dev(I)    , & 
!                  SNRLS_dev(I)    , &
                  PRN_LS_above    , &
                  PSN_LS_above    , &
                  EVAP_DD_LS_above, &
                  SUBL_DD_LS_above, &
!                  REV_LS_dev(I,K) , &
!                  RSU_LS_dev(I,K) , &    
!                  ACLL_LS_dev(I,K), &
!                  ACIL_LS_dev(I,K), &
!                  PFL_LS_dev(I,K) , &
!                  PFI_LS_dev(I,K) , &
!                  VFALLRN         , &
!                  VFALLSN         , &
!                  FRZ_PP_dev(I,K) , &
                  ENVFC, LSDDRFC    )

!            VFALLSN_LS_dev(I,K) = VFALLSN
!            VFALLRN_LS_dev(I,K) = VFALLRN
      
!            if (.not.use_autoconv_timescale) then
!               if (VFALLSN.NE.0.) then
!                  QSN_ALL = QSN_ALL + PFI_LS_dev(I,K)/VFALLSN
!               end if
!               if (VFALLRN.NE.0.) then
!                  QRN_ALL = QRN_ALL + PFL_LS_dev(I,K)/VFALLRN
!               end if
!               if (VFALLRN.NE.0. .AND. VFALLSN.NE.0.) then
!                  QPLS_dev(I,K) = QPLS_dev(I,K) + PFL_LS_dev(I,K)/VFALLRN + PFI_LS_dev(I,K)/VFALLSN
!               end if 
!            end if

            IF ( (QLW_LS_dev(I,K)+QLW_AN_dev(I,K)) > 0.00 ) THEN
               QTMP3 = 1./(QLW_LS_dev(I,K)+QLW_AN_dev(I,K))
            ELSE
               QTMP3 = 0.0
            END IF
            QLW_LS_dev(I,K) = QLW_LS_dev(I,K) * QTMP1 * QTMP3
            QLW_AN_dev(I,K) = QLW_AN_dev(I,K) * QTMP1 * QTMP3
     
            IF ( (QIW_LS_dev(I,K)+QIW_AN_dev(I,K)) > 0.00 ) THEN
               QTMP3 = 1./(QIW_LS_dev(I,K)+QIW_AN_dev(I,K))
            ELSE
               QTMP3 = 0.0
            END IF
            QIW_LS_dev(I,K) = QIW_LS_dev(I,K) * QTMP2 * QTMP3
            QIW_AN_dev(I,K) = QIW_AN_dev(I,K) * QTMP2 * QTMP3

         
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            TH_dev(I,K)  =  TEMP / EXNP_dev(I,K) 

!            QRN_ALL = QRN_ALL / (100.*PP_dev(I,K) / (MAPL8_RGAS*TEMP ))
!            QSN_ALL = QSN_ALL / (100.*PP_dev(I,K) / (MAPL8_RGAS*TEMP ))

!            RHCLR_dev(I,K) = MIN( CLDFRC_dev(I,K) + ANVFRC_dev(I,K), 1.00 )
!            IF ( RHCLR_dev(I,K) < 1.00 ) THEN
!               RHCLR_dev(I,K) = ( MIN( Q_dev(I,K)/QSAT(TEMP,PP_dev(I,K)),1.0 ) - RHCLR_dev(I,K) ) / &
!                                     (1. - RHCLR_dev(I,K))
!               IF ( RHCLR_dev(I,K) < 0.00 ) THEN
!                  RHCLR_dev(I,K) = MAPL8_UNDEF
!               END IF
!            ELSE
!               RHCLR_dev(I,K) = MAPL8_UNDEF
!            END IF

!            IF (DISABLE_RAD==1) THEN
!               RAD_QL_dev(I,K)     = 0.
!               RAD_QI_dev(I,K)     = 0.
!               RAD_QR_dev(I,K)     = 0.
!               RAD_QS_dev(I,K)     = 0.
!               RAD_CLDFRC_dev(I,K) = 0.
!               CLDREFFL_dev(I,K)   = 0.
!               CLDREFFI_dev(I,K)   = 0.
!            ELSE
!               call RADCOUPLE ( TEMP, PP_dev(I,K), CLDFRC_dev(I,K), ANVFRC_dev(I,K), &
!                     QLW_LS_dev(I,K), QIW_LS_dev(I,K), QLW_AN_dev(I,K), QIW_AN_dev(I,K), QRN_ALL, QSN_ALL, & 
!                     RAD_QL_dev(I,K), RAD_QI_dev(I,K), RAD_QR_dev(I,K), RAD_QS_dev(I,K), RAD_CLDFRC_dev(I,K), & 
!                     CLDREFFL_dev(I,K), CLDREFFI_dev(I,K), CLDVOL2FRC,N_ANVIL*1.e6,N_ICE*1.e6,N_WARM*1.e6, &
!                     TEMPOR_dev(I))
!            END IF

            QRN_CU_dev(I,K) = QRN_CU_1D

         end do K_LOOP

!#ifndef _CUDA
      end do RUN_LOOP
!#else
!      end if RUN_LOOP
!#endif

   END SUBROUTINE PROGNO_CLOUD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!              P R O C E S S   S U B R O U T I N E S                 !!
!!                         * * * * *                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!              P R O C E S S   S U B R O U T I N E S                 !!
!!                         * * * * *                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!              P R O C E S S   S U B R O U T I N E S                 !!
!!                         * * * * *                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device) &
!#endif
   subroutine pdf_spread (K,LM,&
         U,U_above,U_below,&
         V,V_above,V_below,&
         KH,&
         DZ_above,DZ_below,&
         UPDF,PP,ALPHA,&
!         ALPHT_DIAG, ALPH1_DIAG, ALPH2_DIAG,&
         FRLAND )  

      integer, intent(in)  :: k,lm
      real(8),    intent(in)  :: U,U_above,U_below
      real(8),    intent(in)  :: V,V_above,V_below
      real(8),    intent(in)  :: DZ_above,DZ_below
      real(8),    intent(in)  :: UPDF,PP
      real(8),    intent(in)  :: KH
      real(8),    intent(out) :: ALPHA
!      real(8),    intent(out) :: ALPH1_DIAG, ALPH2_DIAG, ALPHT_DIAG
      real(8),    intent(in)  :: FRLAND

      real(8) :: A1,A2,A3
      real(8) :: tempmaxrh

      ! alpha is the 1/2*width so RH_crit=1.0-alpha

      if (tanhrhcrit.eq.1) then

         !  Use Slingo-Ritter (1985) formulation for critical relative humidity
         !  array a1 holds the critical rh, ranges from 0.8 to 1

         tempmaxrh = maxrhcrit
         if (frland > 0.05) tempmaxrh = maxrhcritland
         a1 = 1.0
         if (pp .le. turnrhcrit) then
            a1 = minrhcrit
         else
            a1 = minrhcrit + (tempmaxrh-minrhcrit)/(19.) * &
                  ((atan( (2.*(pp- turnrhcrit)/(1020.-turnrhcrit)-1.) * &
                  tan(20.*pi_0/21.-0.5*pi_0) ) + 0.5*pi_0) * 21./pi_0 - 1.)
         end if

         a1 = min(a1,1.)
         alpha = 1. - a1

      else

         alpha = 0.001 ! 0.1% RH SLOP

         !! DIRECTIONAL SHEAR == ABS( e_normal dot [U_z,V_z] ) 
   
         A1 = 0.

         A3 = 1./SQRT( U**2 + V**2 + 0.01 )  ! inverse of wind mag 

         A2 = V * A3 ! x-component of unit normal to (U,V) 

         if (k > 1 .and. k < lm) then
            A1 = ( A2 * ( U_above - U_below ) &
                  / ( DZ_above+DZ_below ) )
         end if

         A2 = -U * A3 ! y-component of unit normal to (U,V) 

         if (k > 1 .and. k < lm) then
            A1 = ( A2 * ( V_above - V_below )  & 
                  / ( DZ_above+DZ_below ) )  + A1
         end if

         A1 = ABS( A1 )  ! A1 is now magnitude of veering shear at layers in (m/s) /m.  Thus, A1=.001  ==> 1 m/s/km

         ALPHA = ALPHA  +  10.*A1  

!         ALPH1_DIAG = 10.*A1

         !! Total shear = SQRT( [U_z,V_z] dot [U_z,V_z] )

         A1  = 0.
         if (k > 1 .and. k < lm) then
            A1 = ( ( U_above - U_below )/ ( DZ_above+DZ_below ) )**2 & 
                  + ( ( V_above - V_below )/ ( DZ_above+DZ_below ) )**2  
         end if

         A1  = SQRT ( A1 )  ! A1 is now magnitude of TOTAL shear at layers in (m/s) /m.  Thus, A1=.001  ==> 1 m/s/km

         ALPHA = ALPHA  +  3.33*A1

         !! KH values ~100 m+2 s-1 typical of strong PBLs

         ALPHA = ALPHA  +  0.002*KH

!         ALPH2_DIAG = 0.002*KH

      end if               ! end of slingo ritter if-sequence

      ALPHA = MIN( ALPHA , 0.25 )  ! restrict RHcrit to > 75% 
!      ALPHT_DIAG = ALPHA

   end subroutine pdf_spread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device) &
!#endif
   subroutine fix_up_clouds( &
         QV, &
         TE, &
         QLC,&
         QIC,&
         CF, &
         QLA,&
         QIA,&
         AF  )

      real(8), intent(inout) :: TE,QV,QLC,CF,QLA,AF,QIC,QIA

      ! Fix if Anvil cloud fraction too small
      if (AF < 1.E-5) then
         QV  = QV + QLA + QIA
         TE  = TE - (MAPL8_ALHL/MAPL8_CP)*QLA - (MAPL8_ALHS/MAPL8_CP)*QIA
         AF  = 0.
         QLA = 0.
         QIA = 0.
      end if

      ! Fix if LS cloud fraction too small
      if ( CF < 1.E-5 ) then
         QV = QV + QLC + QIC
         TE = TE - (MAPL8_ALHL/MAPL8_CP)*QLC - (MAPL8_ALHS/MAPL8_CP)*QIC
         CF  = 0.
         QLC = 0.
         QIC = 0.
      end if
      
      ! LS LIQUID too small
      if ( QLC  < 1.E-8 ) then
         QV = QV + QLC 
         TE = TE - (MAPL8_ALHL/MAPL8_CP)*QLC
         QLC = 0.
      end if
      ! LS ICE too small
      if ( QIC  < 1.E-8 ) then
         QV = QV + QIC 
         TE = TE - (MAPL8_ALHS/MAPL8_CP)*QIC
         QIC = 0.
      end if

      ! Anvil LIQUID too small
      if ( QLA  < 1.E-8 ) then
         QV = QV + QLA 
         TE = TE - (MAPL8_ALHL/MAPL8_CP)*QLA
         QLA = 0.
      end if
      ! Anvil ICE too small
      if ( QIA  < 1.E-8 ) then
         QV = QV + QIA 
         TE = TE - (MAPL8_ALHS/MAPL8_CP)*QIA
         QIA = 0.
      end if

      ! Fix ALL cloud quants if Anvil cloud LIQUID+ICE too small
      if ( ( QLA + QIA ) < 1.E-8 ) then
         QV = QV + QLA + QIA
         TE = TE - (MAPL8_ALHL/MAPL8_CP)*QLA - (MAPL8_ALHS/MAPL8_CP)*QIA
         AF  = 0.
         QLA = 0.
         QIA = 0.
      end if
      ! Ditto if LS cloud LIQUID+ICE too small
      if ( ( QLC + QIC ) < 1.E-8 ) then
         QV = QV + QLC + QIC
         TE = TE - (MAPL8_ALHL/MAPL8_CP)*QLC - (MAPL8_ALHS/MAPL8_CP)*QIC
         CF  = 0.
         QLC = 0.
         QIC = 0.
      end if

   end subroutine fix_up_clouds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device) &
!#endif
   subroutine meltfrz( DT, TE, QL, QI )

      real(8), intent(in)    :: DT
      real(8), intent(inout) :: TE,QL,QI

      real(8)  :: fQi,dQil

      real(8)  ::  taufrz

      integer :: K

      taufrz = 1000.

      fQi  = ice_fraction( TE )
      dQil = 0.
      ! freeze liquid
      if ( TE <= T_ICE_MAX ) then
         dQil = Ql *(1.0 - EXP( -Dt * fQi / taufrz ) )
      end if
      dQil = max(  0., dQil )
      Qi   = Qi + dQil
      Ql   = Ql - dQil
      TE   = TE + (MAPL8_ALHS-MAPL8_ALHL)*dQil/MAPL8_CP

      dQil = 0.
      ! melt ice instantly above 0^C
      if ( TE > T_ICE_MAX ) then
         dQil = -Qi 
      end if
      !! Is this a bug?   dQil = max(  0., dQil )
      dQil = min(  0., dQil )
      Qi   = Qi + dQil
      Ql   = Ql - dQil
      TE   = TE + (MAPL8_ALHS-MAPL8_ALHL)*dQil/MAPL8_CP

   end subroutine meltfrz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device) &
!#endif
   subroutine hystpdf( &
         DT          , &
         ALPHA       , &
         PDFSHAPE    , &
         PL          , &
         QV          , &
         QCl         , &
         QAl         , &
         QCi         , &
         QAi         , &
         TE          , &
         CF          , &
         AF          )

      real(8), intent(in)    :: DT,ALPHA,PL
      integer, intent(in) :: pdfshape
      real(8), intent(inout) :: TE,QV,QCl,QCi,CF,QAl,QAi,AF

      ! internal arrays
      real(8) :: QCO, QVO, CFO, QAO, TAU
      real(8) :: QT, QMX, QMN, DQ, QVtop, sigmaqt1, sigmaqt2

      real(8) :: TEO,QSx,DQsx,QS,DQs

      real(8) :: TEp, QSp, CFp, QVp, QCp
      real(8) :: TEn, QSn, CFn, QVn, QCn

      real(8) :: QCx, QVx, CFx, QAx, QC, QA, fQi, fQi_A
      real(8) :: dQAi, dQAl, dQCi, dQCl 

      real(8) :: tmpARR
      real(8) :: ALHX
      ! internal scalars
      integer :: N

      QC = QCl + QCi
      QA = QAl + QAi

      if ( QA > 0.0 ) then
         fQi_A = QAi / QA 
      else
         fQi_A = 0.0
      end if

      TEo = TE

      call DQSATSCApert(DQSx,QSx,TEo,PL)

      if ( AF < 1.0 ) then
         tmpARR = 1./(1.-AF)
      else
         tmpARR = 0.0
      end if

      CFx = CF*tmpARR
      QCx = QC*tmpARR
      QVx = ( QV - QSx*AF )*tmpARR

      if ( AF >= 1.0 ) then
         QVx = QSx*1.e-4 
      end if


      if ( AF > 0. ) then
         QAx = QA/AF
      else
         QAx = 0.
      end if

      QT  = QCx + QVx

      TEp = TEo
      QSn = QSx
      TEn = TEo
      CFn = CFx
      QVn = QVx
      QCn = QCx

      DQS = DQSx

      do n=1,4

         QSp = QSn
         QVp = QVn
         QCp = QCn
         CFp = CFn

         call DQSATSCApert(DQS,QSn,TEn,PL)

         TEp = TEn
         fQi = ice_fraction( TEp )

         sigmaqt1  = ALPHA*QSn
         sigmaqt2  = ALPHA*QSn

         if(pdfflag.eq.2) then
! for triangular, symmetric: sigmaqt1 = sigmaqt2 = alpha*qsn (alpha is half width)
! for triangular, skewed r : sigmaqt1 < sigmaqt2
! try: skewed right below 500 mb
!!!       if(pl.lt.500.) then
          sigmaqt1  = ALPHA*QSn
          sigmaqt2  = ALPHA*QSn
!!!       else
!!!       sigmaqt1  = 2*ALPHA*QSn*0.4
!!!       sigmaqt2  = 2*ALPHA*QSn*0.6
!!!       endif
         endif

         call pdffrac(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,CFn)
         call pdfcondensate(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,QCn)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! These lines represent adjustments
         ! to anvil condensate due to the 
         ! assumption of a stationary TOTAL 
         ! water PDF subject to a varying 
         ! QSAT value during the iteration
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
         if ( AF > 0. ) then
            QAo = QAx  ! + QSx - QS 
         else
            QAo = 0.
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ALHX = (1.0-fQi)*MAPL8_ALHL + fQi*MAPL8_ALHS

         if(pdfflag.eq.1) then
          QCn = QCp + ( QCn - QCp ) / ( 1. - (CFn * (ALPHA-1.) - (QCn/QSn))*DQS*ALHX/MAPL8_CP)
         elseif(pdfflag.eq.2) then
! This next line needs correcting - need proper d(del qc)/dT derivative for triangular
! for now, just use relaxation of 1/2.
          if (n.ne.4) QCn = QCp + ( QCn - QCp ) *0.5
         endif
         QVn = QVp - (QCn - QCp)

         TEn = TEp + (1.0-fQi)*(MAPL8_ALHL/MAPL8_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF ) &
               +      fQi* (MAPL8_ALHS/MAPL8_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF )

      enddo ! qsat iteration

         CFo = CFn
         CF = CFn
         QCo = QCn
         QVo = QVn
         TEo = TEn

      ! Update prognostic variables.  Deal with special case of AF=1
      ! Temporary variables QCo, QAo become updated grid means.
      if ( AF < 1.0 ) then
         CF  = CFo * ( 1.-AF)
         QCo = QCo * ( 1.-AF)
         QAo = QAo *   AF  
      else

         ! Special case AF=1, i.e., box filled with anvil. 
         !   - Note: no guarantee QV_box > QS_box

         CF  = 0.          ! Remove any other cloud
         QAo = QA  + QC    ! Add any LS condensate to anvil type
         QCo = 0.          ! Remove same from LS   

         QT  = QAo + QV    ! Total water

         ! Now set anvil condensate to any excess of total water 
         ! over QSx (saturation value at top)
         QAo = MAX( QT - QSx, 0. )
      end if

      ! Now take {\em New} condensate and partition into ice and liquid
      ! taking care to keep both >=0 separately. New condensate can be
      ! less than old, so $\Delta$ can be < 0.

      QCx   = QCo - QC
      dQCl  = (1.0-fQi)*QCx
      dQCi  =    fQi  * QCx

      if ((QCl+dQCl)<0.) then
         dQCi = dQCi + (QCl+dQCl)
         dQCl = -QCl !== dQCl - (QCl+dQCl)
      end if
      if ((QCi+dQCi)<0.) then
         dQCl = dQCl + (QCi+dQCi)
         dQCi = -QCi !== dQCi - (QCi+dQCi)
      end if

      QAx   = QAo - QA
      dQAl  = QAx ! (1.0-fQi)*QAx
      dQAi  = 0.  !  fQi  * QAx

      if ((QAl+dQAl)<0.) then
         dQAi = dQAi + (QAl+dQAl)
         dQAl = -QAl
      end if
      if ((QAi+dQAi)<0.) then
         dQAl = dQAl + (QAi+dQAi)
         dQAi = -QAi 
      end if

      ! Clean-up cloud if fractions are too small
      if ( AF < 1.e-5 ) then
         dQAi = -QAi
         dQAl = -QAl
      end if
      if ( CF < 1.e-5 ) then
         dQCi = -QCi
         dQCl = -QCl
      end if

      QAi    = QAi + dQAi
      QAl    = QAl + dQAl
      QCi    = QCi + dQCi
      QCl    = QCl + dQCl
      QV     = QV  - ( dQAi+dQCi+dQAl+dQCl) 


      !!TE  = TE + (MAPL8_ALHS/MAPL8_CP)*(dQAi+dQCi) + (MAPL8_ALHL/MAPL8_CP)*(dQAl+dQCl)
      TE  = TE + (MAPL8_ALHL*( dQAi+dQCi+dQAl+dQCl)+MAPL8_ALHF*(dQAi+dQCi))/ MAPL8_CP

      ! We need to take care of situations where QS moves past QA
      ! during QSAT iteration. This should be only when QA/AF is small
      ! to begin with. Effect is to make QAo negative. So, we 
      ! "evaporate" offending QAs
      !
      ! We get rid of anvil fraction also, although strictly
      ! speaking, PDF-wise, we should not do this.
      if ( QAo <= 0. ) then
         QV  = QV + QAi + QAl
         TE  = TE - (MAPL8_ALHS/MAPL8_CP)*QAi - (MAPL8_ALHL/MAPL8_CP)*QAl
         QAi = 0.
         QAl = 0.
         AF  = 0.  
      end if

   end subroutine hystpdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device) &
!#endif
      subroutine pdffrac (flag,qtmean,sigmaqt1,sigmaqt2,qstar,clfrac)
      implicit none

      integer flag            ! flag to indicate shape of pdf
                              ! 1 for tophat, 2 for triangular, 3 for Gaussian
      real(8) qtmean             ! Grid box value of q total
      real(8) sigmaqt1           ! width of distribution (sigma)
      real(8) sigmaqt2           ! width of distribution (sigma)
      real(8) qstar              ! saturation q at grid box avg T
      real(8) clfrac             ! cloud fraction (area under pdf from qs)

      real(8) :: qtmode, qtmin, qtmax

      if(flag.eq.1) then
       if((qtmean+sigmaqt1).lt.qstar) then
        clfrac = 0.
       else
        if(sigmaqt1.gt.0.) then
        clfrac = min((qtmean + sigmaqt1 - qstar),2.*sigmaqt1)/(2.*sigmaqt1)
        else
        clfrac = 1.
        endif
       endif
      elseif(flag.eq.2) then
       qtmode =  qtmean + (sigmaqt1-sigmaqt2)/3.
       qtmin = min(qtmode-sigmaqt1,0.)
       qtmax = qtmode + sigmaqt2
       if(qtmax.lt.qstar) then
        clfrac = 0.
       elseif ( (qtmode.le.qstar).and.(qstar.lt.qtmax) ) then
        clfrac = (qtmax-qstar)*(qtmax-qstar) / ( (qtmax-qtmin)*(qtmax-qtmode) )
       elseif ( (qtmin.le.qstar).and.(qstar.lt.qtmode) ) then
        clfrac = 1. - ( (qstar-qtmin)*(qstar-qtmin) / ( (qtmax-qtmin)*(qtmode-qtmin) ) )
       elseif ( qstar.le.qtmin ) then
        clfrac = 1.
       endif
      endif

      return
      end subroutine pdffrac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef _CUDA
!   attributes(device) &
!#endif

      subroutine pdfcondensate (flag,qtmean4,sigmaqt14,sigmaqt24,qstar4,condensate4)
      implicit none

      integer flag            ! flag to indicate shape of pdf
                              ! 1 for tophat, 2 for triangular
      real(8) qtmean4            ! Grid box value of q total
      real(8) sigmaqt14          ! width of distribution (to left)
      real(8) sigmaqt24          ! width of distribution (to right)
      real(8) qstar4             ! saturation q at grid box avg T
      real(8) condensate4        ! condensate (area under (q*-qt)*pdf from qs)

      real(8) :: qtmode, qtmin, qtmax, constA, constB, cloudf
      real(8) :: term1, term2, term3
      real(8) :: qtmean, sigmaqt1, sigmaqt2, qstar, condensate

      qtmean = (qtmean4)
      sigmaqt1 = (sigmaqt14)
      sigmaqt2 = (sigmaqt24)
      qstar = (qstar4)

      if(flag.eq.1) then
       if(qtmean+sigmaqt1.lt.qstar) then
        condensate = 0.0
       elseif(qstar.gt.qtmean-sigmaqt1)then
        if(sigmaqt1.gt.0.0) then
         condensate = (min(qtmean + sigmaqt1 - qstar,2.0*sigmaqt1)**2)/ (4.0*sigmaqt1)
        else
         condensate = qtmean-qstar
        endif
       else
        condensate = qtmean-qstar
       endif
      elseif(flag.eq.2) then
       qtmode =  qtmean + (sigmaqt1-sigmaqt2)/3.0
       qtmin = min(qtmode-sigmaqt1,0.0)
       qtmax = qtmode + sigmaqt2
       if ( qtmax.lt.qstar ) then
        condensate = 0.0
       elseif ( (qtmode.le.qstar).and.(qstar.lt.qtmax) ) then
        constB = 2.0 / ( (qtmax - qtmin)*(qtmax-qtmode) )
        cloudf = (qtmax-qstar)*(qtmax-qstar) / ( (qtmax-qtmin)*(qtmax-qtmode) )
        term1 = (qstar*qstar*qstar)/3.0
        term2 = (qtmax*qstar*qstar)/2.0
        term3 = (qtmax*qtmax*qtmax)/6.0
        condensate = constB * (term1-term2+term3) - qstar*cloudf
       elseif ( (qtmin.le.qstar).and.(qstar.lt.qtmode) ) then
        constA = 2.0 / ( (qtmax - qtmin)*(qtmode-qtmin) )
        cloudf = 1.0 - ( (qstar-qtmin)*(qstar-qtmin) / ( (qtmax-qtmin)*(qtmode-qtmin) ) )
        term1 = (qstar*qstar*qstar)/3.0
        term2 = (qtmin*qstar*qstar)/2.0
        term3 = (qtmin*qtmin*qtmin)/6.0
        condensate = qtmean - ( constA * (term1-term2+term3) ) - qstar*cloudf
       elseif ( qstar.le.qtmin ) then
        condensate = qtmean-qstar
       endif
      endif
      condensate4 = (condensate)

      return
      end subroutine pdfcondensate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device) & 
!#endif
   subroutine cnvsrc( & 
         DT      , &
         ICEPARAM, &
         SCLMFDFR, &
         MASS    , &
         iMASS   , &
         PL      , &
         TE      , &
         QV      , &
         DCF     , &
         DMF     , &
         QLA     , & 
         QIA     , & 
         CF      , &
         AF      , &
         QS        )

      !INPUTS:
      !
      !       ICEPARAM: 0-1  controls how strongly new conv condensate is partitioned in ice-liquid
      !                 1 means partitioning follows ice_fraction(TE). 0 means all new condensate is
      !                 liquid 
      !
      !       SCLMFDFR: Scales detraining mass flux to a cloud fraction source - kludge. Thinly justified
      !                 by fuzziness of cloud boundaries and existence of PDF of condensates (for choices
      !                 0.-1.0) or by subgrid layering (for choices >1.0) 

      real(8), intent(in)    :: DT,ICEPARAM,SCLMFDFR
      real(8), intent(in)    :: MASS,iMASS,QS
      real(8), intent(in)    :: DMF,PL
      real(8), intent(in)    :: DCF,CF
      real(8), intent(inout) :: TE
      real(8), intent(inout) :: AF,QV
      real(8), intent(inout) :: QLA, QIA

      real(8) :: TEND,QVx,QCA,fQi

      integer  :: STRATEGY
      real(8)  :: minrhx 

      STRATEGY = 1

      !Minimum allowed env RH
      minrhx    = 0.001  

      !Addition of condensate from RAS
      TEND = DCF*iMASS
      fQi  = 0.0 + ICEPARAM*ice_fraction( TE )
      QLA  = QLA + (1.0-fQi)* TEND*DT
      QIA  = QIA +    fQi   * TEND*DT

      ! dont forget that conv cond has never frozen !!!!
      TE   = TE +   (MAPL8_ALHS-MAPL8_ALHL) * fQi * TEND*DT/MAPL8_CP

      QCA  = QLA + QIA

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Tiedtke-style anvil fraction !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TEND=DMF*iMASS * SCLMFDFR
      AF = AF + TEND*DT
      AF = MIN( AF , 0.99 )

      ! where ( (AF+CF) > 1.00 )
      !    AF=1.00-CF
      ! endwhere
      !!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Check for funny (tiny, negative)
      ! external QV s, resulting from assumed
      ! QV=QSAT within anvil.
      !
      ! Two strategies to fix 
      !   1) Simply constrain AF assume condensate
      !      just gets "packed" in
      !   2) Evaporate QCA to bring up QVx leave AF alone
      !      Should include QSAT iteration, but ...        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !   QS  = QSAT(         &
      !               TE    , &
      !               PL      )

      if ( AF < 1.0 ) then
         QVx  = ( QV - QS * AF )/(1.-AF)
      else
         QVx  = QS
      end if

      if (STRATEGY==1) then 
         if ( (( QVx - minrhx*QS ) < 0.0 ) .and. (AF > 0.) ) then
            AF  = (QV  - minrhx*QS )/( QS*(1.0-minrhx) )
         end if
         if ( AF < 0. ) then  ! If still cant make suitable env RH then destroy anvil
            AF  = 0.
            QV  = QV + QLA + QIA
            TE  = TE - (MAPL8_ALHL*QLA + MAPL8_ALHS*QIA)/MAPL8_CP        
            QLA = 0.
            QIA = 0.
         end if
      else if (STRATEGY==2) then
         if ( (( QVx - minrhx*QS ) < 0.0 ) .and. (AF > 0.) ) then
            QV  = QV  + (1.-AF)*( minrhx*QS - QVx )
            QCA = QCA - (1.-AF)*( minrhx*QS - QVx )
            TE  = TE  - (1.-AF)*( minrhx*QS - QVx )*MAPL8_ALHL/MAPL8_CP
         end if
      end if

   end subroutine cnvsrc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device) &
!#endif
   subroutine evap3(&
         DT      , &
         RHCR    , &
         PL      , &
         TE      , &
         QV      , &
         QL      , &
         QI      , &
         F       , &
         XF      , & 
         QS        )

      real(8), intent(in   ) :: DT 
      real(8), intent(in   ) :: RHCR
      real(8), intent(in   ) :: PL
      real(8), intent(inout) :: TE
      real(8), intent(inout) :: QV
      real(8), intent(inout) :: QL,QI
      real(8), intent(inout) :: F
      real(8), intent(in   ) :: XF
      real(8), intent(in   ) :: QS

      real(8) :: ES,NN,RADIUS,K1,K2,TEFF,QCm,EVAP,RHx,QC  !,QS

      real(8), parameter :: EPSILON =  MAPL8_H2OMW/MAPL8_AIRMW
      real(8), parameter :: K_COND  =  2.4e-2        ! J m**-1 s**-1 K**-1
      real(8), parameter :: DIFFU   =  2.2e-5        ! m**2 s**-1

      real(8) :: A_eff

      A_EFF = CLD_EVP_EFF

      NN = 50.*1.0e6

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!         EVAPORATION OF CLOUD WATER.             !!
      !!                                                 !!
      !!  DelGenio et al (1996, J. Clim., 9, 270-303)    !!
      !!  formulation  (Eq.s 15-17)                      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !   QS  = QSAT(         &
      !               TE    , &
      !               PL      )
      
      ES = 100.* PL * QS  / ( (EPSILON) + (1.0-(EPSILON))*QS )  ! (100's <-^ convert from mbar to Pa)

      RHx = MIN( QV/QS , 1.00 )


      K1 = (MAPL8_ALHL**2) * RHO_W / ( K_COND*MAPL8_RVAP*(TE**2))

      K2 = MAPL8_RVAP * TE * RHO_W / ( DIFFU * (1000./PL) * ES )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Here DIFFU is given for 1000 mb  !!
      !! so 1000./PR accounts for inc-    !!
      !! reased diffusivity at lower      !!
      !! pressure.                        !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      if ( ( F > 0.) .and. ( QL > 0. ) ) then
         QCm=QL/F
      else
         QCm=0.
      end if

      RADIUS = LDRADIUS3(PL,TE,QCm,NN)

      if ( (RHx < RHCR ) .and.(RADIUS > 0.0) ) then
         TEFF   =   (RHCR - RHx) / ((K1+K2)*RADIUS**2)  ! / (1.00 - RHx)
      else
         TEFF   = 0.0 ! -999.
      end if

      EVAP = a_eff*QL*DT*TEFF
      EVAP = MIN( EVAP , QL  )

      QC=QL+QI
      if (QC > 0.) then
         F    = F * ( QC - EVAP ) / QC
      end if

      QV   = QV   + EVAP
      QL   = QL   - EVAP
      TE   = TE   - (MAPL8_ALHL/MAPL8_CP)*EVAP

   end subroutine evap3



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device) &
!#endif
   subroutine subl3( &
         DT        , &
         RHCR      , &
         PL        , &
         TE        , &
         QV        , &
         QL        , &
         QI        , &
         F         , &
         XF        , & 
         QS        )

      real(8), intent(in   ) :: DT
      real(8), intent(in   ) :: RHCR
      real(8), intent(in   ) :: PL
      real(8), intent(inout) :: TE
      real(8), intent(inout) :: QV
      real(8), intent(inout) :: QL,QI
      real(8), intent(inout) :: F
      real(8), intent(in   ) :: XF
      real(8), intent(in   ) :: QS

      real(8) :: ES,NN,RADIUS,K1,K2,TEFF,QCm,SUBL,RHx,QC !, QS

      real(8), parameter :: EPSILON =  MAPL8_H2OMW/MAPL8_AIRMW
      real(8), parameter :: K_COND  =  2.4e-2        ! J m**-1 s**-1 K**-1
      real(8), parameter :: DIFFU   =  2.2e-5        ! m**2 s**-1

      real(8) :: A_eff

      A_EFF = CLD_EVP_EFF

      NN = 5.*1.0e6

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!         SUBLORATION OF CLOUD WATER.             !!
      !!                                                 !!
      !!  DelGenio et al (1996, J. Clim., 9, 270-303)    !!
      !!  formulation  (Eq.s 15-17)                      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !   QS  = QSAT(         &
      !               TE    , &
      !               PL      )

      ES = 100.* PL * QS  / ( (EPSILON) + (1.0-(EPSILON))*QS )  ! (100s <-^ convert from mbar to Pa)

      RHx = MIN( QV/QS , 1.00 )


      K1 = (MAPL8_ALHL**2) * RHO_W / ( K_COND*MAPL8_RVAP*(TE**2))

      K2 = MAPL8_RVAP * TE * RHO_W / ( DIFFU * (1000./PL) * ES )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Here DIFFU is given for 1000 mb  !!
      !! so 1000./PR accounts for inc-    !!
      !! reased diffusivity at lower      !!
      !! pressure.                        !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      if ( ( F > 0.) .and. ( QI > 0. ) ) then
         QCm=QI/F
      else
         QCm=0.
      end if

      RADIUS = LDRADIUS3(PL,TE,QCm,NN)

      if ( (RHx < RHCR ) .and.(RADIUS > 0.0) ) then
         TEFF   =   ( RHCR - RHx) / ((K1+K2)*RADIUS**2)  ! / (1.00 - RHx)
      else
         TEFF   = 0.0 ! -999.
      end if

      SUBL = a_eff*QI*DT*TEFF
      SUBL = MIN( SUBL , QI  )

      QC=QL+QI
      if (QC > 0.) then
         F    = F * ( QC - SUBL ) / QC
      end if

      QV   = QV   + SUBL
      QI   = QI   - SUBL
      TE   = TE   - (MAPL8_ALHS/MAPL8_CP)*SUBL

   end subroutine subl3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef _CUDA
!   attributes(device) &
!#endif
   subroutine autocon3( &
         DT       , &
         QC       , &
         QP       , &
         TE       , &
         PL       , &
         KH       , &
         F        , &
         SUNDQV2  , &
         SUNDQV3  , &
         SUNDQT1  , &
         DZET     , &
         VF       , &
         FRACTION_REMOVAL )

      integer, intent(in) :: FRACTION_REMOVAL

      real(8), intent(in   ) :: DT
      real(8), intent(in   ) :: TE
      real(8), intent(in   ) :: PL
      real(8), intent(in   ) :: KH
      real(8), intent(inout) :: QC
      real(8), intent(inout) :: QP
      real(8), intent(inout) :: F
      real(8), intent(in   ) :: DZET
      real(8), intent(inout) :: VF

      real(8), intent(in   ) :: SUNDQV2, SUNDQV3, SUNDQT1

      integer :: NSMX

      real(8) :: ACF0, ACF, DTX
      real(8) :: C00x, iQCcrx, F2, F3,RATE,dQP,QCm
      real(8) :: dqfac
      integer :: NS, K

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Precip. conversion from Smith (1990,    !
      !   QJRMS, 116, 435, Eq. 2.29). Similar    ! 
      !   to Del Genios Eq.(10).                 !
      !                                          !
      !   Coalesence term needs to be determined !
      !   through entire column and is done in   !
      !   subroutine "ACCRETE_EVAP_PRECIP"       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      NSMX = NSMAX
      DTX  = DT/NSMX

      CALL SUNDQ3_ICE3(TE, SUNDQV2, SUNDQV3, SUNDQT1, F2, F3 )

      C00x  = C_00 * F2 * F3
      !QCcrx = LWCRIT / ( F2 * F3 )
      iQCcrx = F2 * F3 / LWCRIT

      if ( ( F > 0.) .and. ( QC > 0. ) )then
         QCm = QC/F
      else
         QCm = 0.
      end if


      RATE = C00x * ( 1.0 - EXP( - ( QCm * iQCcrx )**2 ) )

      !!! Make up a fall velocity for liquid precipitation - in analogy to falling ice,
      !!!    think of the fall velocity as the ratio of autoconverted to existing condensate 
      !!!    multiplied by delta z / delta t  
      !!!   (ie, autoconversion/sec is related to residence time in layer)

      VF  =  (DZET/DT) * ( 1.0 - EXP( -RATE * DT ) )

      !! temporary kluge until we can figure a better to make
      !! thicker low clouds ( reuse arrays F2 and F3 )
      F2 = 1.0
      F3 = 1.0 

      ! Implement ramps for gradual change in autoconv
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Thicken low high lat clouds
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( PL .GE. 775.  .AND. TE .LE.  275. ) then
!!!      F3 = max(-0.016 * PL + 13.4, 0.2)
         F3 = 0.2
      end if
      if ( PL .GE. 825.  .AND. TE .LE.  282. ) then
!!!      F3 = max(0.11 * TE - 30.02, 0.2)
         F3 = 0.2
      end if
      if ( PL .GE. 775.  .AND. PL .LT. 825. .AND. TE .LE.  282. .AND. TE .GT. 275.) then
!!!      F3 = min(max(-0.016*PL + 0.11 * TE - 16.85, 0.2),1.)
         F3 = 0.2
      end if
      if ( PL .GE. 825.  .AND. TE .LE.  275. ) then
         F3 = 0.2
      end if
      if ( PL .LE. 775.  .OR. TE .GT.  282. ) then
         F3 = 1.
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Thin-out low tropical clouds
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
         F3 = min(0.2 * TE - 56, 2.)
      end if
      if ( PL .GE. 925.  .AND. TE .GE.  290. ) then
         F3 = min(0.04 * PL - 36., 2.)
      end if
      if ( PL .GE. 925.  .AND. PL .LT. 950. .AND. TE .GT.  285. .AND. TE .LT. 290.) then
         F3 = max(min(0.04*PL + 0.2 * TE - 94., 2.),1.)
      end if
      if ( PL .GE. 950.  .AND. TE .GE.  290. ) then
         F3 = 2.
      end if

      F3   = MAX( F3, 0.1 )

      RATE = F3 * RATE

      dQP  =  QC*( 1.0 - EXP( -RATE * DT ) )

      dQP  =  MAX( dQP , 0.0 )  ! Protects against floating point problems for tiny RATE

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !! Go ahead and totally wipe-out warm fogs
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !dqfac = 0.
     !if ( PL .GE. 975.  .AND. TE .GE.  280. ) then
     !   dqfac = max(min(0.2 * TE - 56., 1.),0.)
     !end if
     !if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
     !   dqfac = max(min(0.04 * PL - 38., 1.),0.)
     !end if
     !if ( PL .GE. 950.  .AND. PL .LT. 975. .AND. TE .GT.  280. .AND. TE .LT. 285.) then
     !   dqfac = max(min(0.04*PL + 0.2 * TE - 95., 1.),0.)
     !end if
     !if ( ( PL >= 975. ) .AND. (TE >= 285. ) ) then
     !   dqfac = 1.
     !end if
     !dQP = max(dQP, dqfac*QC)

      QC   = QC - dQP  
      QP   = QP + dQP


      !SELECT CASE( FRACTION_REMOVAL )

      !CASE( 0 )
      ! do NOTHING

      IF( FRACTION_REMOVAL == 1 ) THEN
         if ( (QC + dQP) > 0. ) then
            F = QC * F / (QC + dQP )
         end if

      ELSEIF( FRACTION_REMOVAL == 2 ) THEN
         if ( (QC + dQP) > 0. .and. QC > 0. ) then
            F = F * SQRT( QC / (QC + dQP ) )
         end if

      ELSEIF( FRACTION_REMOVAL == 3) THEN
         if ( (QC + dQP) > 0. .and. QC > 0. ) then
            F = F * ( QC / (QC + dQP ) )**0.333
         end if


      END IF

   END SUBROUTINE AUTOCON3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device)   &
!#endif
   subroutine PRECIP3(  &
         K,LM         , &
         DT           , & 
         FRLAND       , & 
         RHCR3        , &
         QPl          , &
         QPi          , &
         QCl          , &
         QCi          , &
         TE           , &
         QV           , &
         mass         , &
         imass        , &
         PL           , &
         dZE          , &
         QDDF3        , &
         AA           , &
         BB           , &
         AREA         , &
!         RAIN         , & 
!         SNOW         , &
         PFl_above    , &
         PFi_above    , &
         EVAP_DD_above, &
         SUBL_DD_above, &
!         REVAP_DIAG   , &
!         RSUBL_DIAG   , &
!         ACRLL_DIAG        , &
!         ACRIL_DIAG        , &
!         PFL_DIAG     , &
!         PFI_DIAG     , &
!         VFALLRN      , &
!         VFALLSN      , &
!         FRZ_DIAG     , &
         ENVFC,DDRFC  )


      integer, intent(in) :: K,LM

      real(8), intent(in   ) :: DT

      real(8), intent(inout) :: QV,QPl,QPi,QCl,QCi,TE

      real(8), intent(in   ) :: mass,imass
      real(8), intent(in   ) :: PL
      real(8), intent(in   ) :: AA,BB
      real(8), intent(in   ) :: RHCR3
      real(8), intent(in   ) :: dZE
      real(8), intent(in   ) :: QDDF3
!      real(8), intent(  out) :: RAIN,SNOW
      real(8), intent(in   ) :: AREA
      real(8), intent(in   ) :: FRLAND

      real(8), intent(inout) :: PFl_above, PFi_above
      real(8), intent(inout) :: EVAP_DD_above, SUBL_DD_above

!      real(8), intent(  out) :: REVAP_DIAG
!      real(8), intent(  out) :: RSUBL_DIAG
!      real(8), intent(  out) :: ACRLL_DIAG,ACRIL_DIAG
!      real(8), intent(  out) :: PFL_DIAG, PFI_DIAG
!      real(8), intent(inout) :: FRZ_DIAG
!      real(8), intent(  out) :: VFALLSN, VFALLRN

      real(8), intent(in   ) :: ENVFC,DDRFC


      real(8) :: PFi,PFl,QS,dQS,ENVFRAC
      real(8) :: TKo,QKo,QSTKo,DQSTKo,RH_BOX,T_ED,QPlKo,QPiKo
      real(8) :: Ifactor,RAINRAT0,SNOWRAT0
      real(8) :: FALLRN,FALLSN,VEsn,VErn,NRAIN,NSNOW,Efactor

      real(8) :: TinLAYERrn,DIAMrn,DROPRAD
      real(8) :: TinLAYERsn,DIAMsn,FLAKRAD

      real(8) :: EVAP,SUBL,ACCR,MLTFRZ,EVAPx,SUBLx
      real(8) :: EVAP_DD,SUBL_DD,DDFRACT
      real(8) :: LANDSEAF

      real(8), parameter :: TRMV_L = 1.0     ! m/s

      real(8) :: TAU_FRZ, TAU_MLT

      integer :: NS, NSMX, itr,L

      logical, parameter :: taneff = .false.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! fraction of precip falling through "environment" vs
      ! through cloud
      real(8), parameter :: B_SUB = 1.00

      if(taneff) then
         envfrac = 1.00

         if (pl .le. 600.) then
            envfrac = 0.25
         else
            envfrac = 0.25 + (1.-0.25)/(19.) *                    &
                  ((atan( (2.*(pl-600.)/(900.-600.)-1.) *       &
                  tan(20.*MAPL8_PI/21.-0.5*MAPL8_PI) ) + 0.5*MAPL8_PI) * 21./MAPL8_PI - 1.)
         end if

         envfrac = min(envfrac,1.)
      else
         ENVFRAC = ENVFC
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( AREA > 0. ) THEN
         Ifactor = 1./ ( AREA )
      ELSE
         Ifactor = 1.00
      END if

      Ifactor = MAX( Ifactor, 1.) !

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !   Start at top of precip column:
      !  
      !               a) Accrete                   
      !               b) Evaporate/Sublimate  
      !               c) Rain/Snow-out to next level down 
      !               d) return to (a)
      !
      !   ....................................................................
      !           
      !  Accretion formulated according to Smith (1990, Q.J.R.M.S., 116, 435
      !  Eq. 2.29)
      !  
      !  Evaporation (ibid. Eq. 2.32)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!! INITIALIZE DIAGNOSTIC ARRAYS !!!!!!!!!!!!!!!!!!!!!
!      PFL_DIAG =  0.
!      PFI_DIAG =  0.
!      ACRIL_DIAG    =  0.
!      ACRLL_DIAG    =  0.
!      REVAP_DIAG    =  0.
!      RSUBL_DIAG    =  0.

      !!!!!!!!!!!!!! UPDATE SATURATED HUMIDITY  !!!!!!!!!!!!!
      call DQSATSCApert(dQS,QS,TE,PL)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DDFRACT = DDRFC 

      IF (K == KTOP) THEN
         PFl=QPl*MASS
         PFi=QPi*MASS

         EVAP_DD = 0.
         SUBL_DD = 0.

!         VFALLRN = 0.0
!         VFALLSN = 0.0
      ELSE 
         QPl   = QPl + PFl_above * iMASS
         PFl = 0.00

         QPi   = QPi + PFi_above * iMASS
         PFi = 0.00


         ACCR = B_SUB * C_ACC * ( QPl*MASS ) *QCl   

         ACCR = MIN(  ACCR , QCl  )

         QPl     = QPl + ACCR
         QCl     = QCl - ACCR

!         ACRLL_DIAG = ACCR / DT

         !! Accretion of liquid condensate by falling ice/snow
         ACCR = B_SUB * C_ACC * ( QPi*MASS ) *QCl   

         ACCR = MIN(  ACCR , QCl  )

         QPi     = QPi + ACCR
         QCl     = QCl - ACCR
         !! Liquid freezes when accreted by snow
         TE      = TE + MAPL8_ALHF*ACCR/MAPL8_CP

 !        ACRIL_DIAG = ACCR / DT

         RAINRAT0 = Ifactor*QPl*MASS/DT
         SNOWRAT0 = Ifactor*QPi*MASS/DT

         call MARSHPALMQ2(RAINRAT0,PL,DIAMrn,NRAIN,FALLrn,VErn)
         call MARSHPALMQ2(SNOWRAT0,PL,DIAMsn,NSNOW,FALLsn,VEsn)

         IF ( FRLAND < 0.1 ) THEN
         !!      DIAMsn = MAX(  DIAMsn, 1.0e-3 )   ! Over Ocean
         END IF

!         VFALLRN = FALLrn
!         VFALLSN = FALLsn

         TinLAYERrn = dZE / ( FALLrn+0.01 )
         TinLAYERsn = dZE / ( FALLsn+0.01 )

         !*****************************************
         !  Melting of Frozen precipitation      
         !*****************************************
         TAU_FRZ = 5000.  ! time scale for freezing (s). 

         MLTFRZ = 0.0
         IF ( (TE > MAPL8_TICE ) .and.(TE <= MAPL8_TICE+5. ) ) THEN
            MLTFRZ= TinLAYERsn * QPi *( TE - MAPL8_TICE ) / TAU_FRZ 
            MLTFRZ= MIN( QPi , MLTFRZ )
            TE  = TE  - MAPL8_ALHF*MLTFRZ/MAPL8_CP
            QPl = QPl + MLTFRZ
            QPi = QPi - MLTFRZ
         END IF
!         FRZ_DIAG = FRZ_DIAG - MLTFRZ / DT

         MLTFRZ = 0.0
         IF ( TE > MAPL8_TICE+5.  ) THEN  ! Go Ahead and melt any snow/hail left above 5 C 
            MLTFRZ= QPi 
            TE  = TE  - MAPL8_ALHF*MLTFRZ/MAPL8_CP
            QPl = QPl + MLTFRZ
            QPi = QPi - MLTFRZ
         END IF
!         FRZ_DIAG = FRZ_DIAG - MLTFRZ / DT

         MLTFRZ = 0.0
         if ( K >= LM-1 ) THEN
            IF ( TE > MAPL8_TICE+0.  ) THEN ! Go Ahead and melt any snow/hail left above 0 C in lowest layers 
               MLTFRZ= QPi 
               TE  = TE  - MAPL8_ALHF*MLTFRZ/MAPL8_CP
               QPl = QPl + MLTFRZ
               QPi = QPi - MLTFRZ
            END IF
         endif
!         FRZ_DIAG = FRZ_DIAG - MLTFRZ / DT

         !*****************************************
         !  Freezing of liquid precipitation      
         !*****************************************
         MLTFRZ = 0.0
         IF ( TE <= MAPL8_TICE ) THEN
            TE  = TE + MAPL8_ALHF*QPl/MAPL8_CP
            QPi = QPl + QPi
            MLTFRZ= QPl
            QPl = 0.
         END IF
!         FRZ_DIAG = FRZ_DIAG + MLTFRZ / DT

         ! ******************************************
         !   In the exp below, evaporation time 
         !   scale is determined "microphysically"
         !   from temp, press, and drop size. In this
         !   context C_EV becomes a dimensionless 
         !   fudge-fraction.
         !   Also remember that these microphysics 
         !   are still only for liquid.
         ! ******************************************

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

            !!!!! RAin falling !!!!!!!!!!!!!!!!!!!!!!!
            if ( ( RH_BOX < RHCR3 ) .AND. ( DIAMrn > 0.00 ) .AND. &
                  ( PL > 100. ) .AND. ( PL < REVAP_OFF_P ) ) then
               DROPRAD=0.5*DIAMrn
               T_ED =  Efactor * DROPRAD**2 
               T_ED =  T_ED * ( 1.0 + DQSTKo*MAPL8_ALHL/MAPL8_CP )
               EVAP =  QPl*(1.0 - EXP( -C_EV_R * VErn * LANDSEAF * ENVFRAC * TinLAYERrn / T_ED ) )
            ELSE
               EVAP = 0.0
            END if

            !!!!! Snow falling !!!!!!!!!!!!!!!!!!!!!!!
            if ( ( RH_BOX < RHCR3 ) .AND. ( DIAMsn > 0.00 ) .AND. &
                  ( PL > 100. ) .AND. ( PL < REVAP_OFF_P ) ) then
               FLAKRAD=0.5*DIAMsn
               T_ED =  Efactor * FLAKRAD**2   
               T_ED =  T_ED * ( 1.0 + DQSTKo*MAPL8_ALHS/MAPL8_CP )
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

            QKo=QV + EVAP + SUBL
            TKo=TE - EVAP * MAPL8_ALHL / MAPL8_CP - SUBL * MAPL8_ALHS / MAPL8_CP

         enddo
      
         QPi  = QPi - SUBL
         QPl  = QPl - EVAP

         !! Put some re-evap/re-subl precip in to a \quote{downdraft} to be applied later
         EVAP_DD = EVAP_DD_above + DDFRACT*EVAP*MASS 
         EVAP    = EVAP          - DDFRACT*EVAP
         SUBL_DD = SUBL_DD_above + DDFRACT*SUBL*MASS 
         SUBL    = SUBL          - DDFRACT*SUBL
         ! -----

         QV   = QV  + EVAP + SUBL
         TE   = TE  - EVAP * MAPL8_ALHL / MAPL8_CP - SUBL * MAPL8_ALHS / MAPL8_CP

!         REVAP_DIAG = EVAP / DT
!         RSUBL_DIAG = SUBL / DT

         PFl  = QPl*MASS
         PFi  = QPi*MASS

!         PFL_DIAG =  PFl/DT
!         PFI_DIAG =  PFi/DT
      end if

      ! QDDF3 (<= QDDF3_dev) is calculated on the CPU in order to avoid
      ! the reverse loop on GPUs and thus save local memory use.
      EVAP = QDDF3*EVAP_DD/MASS
      SUBL = QDDF3*SUBL_DD/MASS
      QV   = QV  + EVAP + SUBL
      TE   = TE  - EVAP * MAPL8_ALHL / MAPL8_CP - SUBL * MAPL8_ALHS / MAPL8_CP
!      REVAP_DIAG = REVAP_DIAG + EVAP / DT
!      RSUBL_DIAG = RSUBL_DIAG + SUBL / DT

!      IF (K == LM) THEN
!         RAIN  = PFl/DT
!         SNOW  = PFi/DT
!      END IF

      QPi = 0.
      QPl = 0.

      PFl_above = PFl
      PFi_above = Pfi

      EVAP_DD_above = EVAP_DD
      SUBL_DD_above = SUBL_DD

   end subroutine precip3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef _CUDA
!   attributes(device) &
!#endif
   subroutine ICEFALL( QI, DZ, QP, VF, F, DT, FRACTION_REMOVAL )

      integer, intent(in) :: FRACTION_REMOVAL

      real(8), intent(inout) :: QI
      real(8), intent(in   ) :: DZ
      real(8), intent(inout) :: QP
      real(8), intent(in   ) :: VF
      real(8), intent(inout) :: F
      real(8), intent(in   ) :: DT

      real(8) :: QIxP

      QIxP = QI * ( VF * DT / DZ ) 
      QIxP = MIN( QIxP , QI )

      QIxP = MAX( QIxP, 0.0 ) ! protects against precision problem

      QP = QP + QIxP
      QI = QI - QIxP

      !SELECT CASE( FRACTION_REMOVAL )

      !CASE( 0 )
         ! do NOTHING

      IF( FRACTION_REMOVAL == 1 ) THEN
         if ( (QI + QIxP) > 0. ) then
            F = QI * F / (QI + QIxP )
         end if

      ELSEIF( FRACTION_REMOVAL == 2 ) THEN
         if ( (QI + QIxP) > 0. .and. QI > 0. ) then
            F = F * SQRT( QI / (QI + QIxP ) )
         end if

      ELSEIF( FRACTION_REMOVAL == 3 ) THEN
         if ( (QI + QIxP) > 0. .and. QI > 0. ) then
            F = F * ( QI / (QI + QIxP ) )**0.333
         end if

      END IF

   end subroutine ICEFALL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef _CUDA
!   attributes(device) &
!#endif
   subroutine SETTLE_VEL( WXR, QI, PL, TE, F, KH, VF, LARGESCALE, ANVIL )

      real(8), intent(in   ) :: WXR 
      real(8), intent(in   ) :: TE
      real(8), intent(in   ) :: QI, F, PL
      real(8), intent(in   ) :: KH
      real(8), intent(out  ) :: VF

      real(8), intent(in) :: ANVIL, LARGESCALE

      real(8) :: RHO, XIm,LXIm, VF_A, VF_L

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Uses Eq. 1 Lawrence and Crutzen (1998, Tellus 50B, 263-289) 
      ! Except midlat form is taken to be for LS cloud, and tropical
      ! form is taken to be for ANVIL cloud
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      RHO = 1000.*100.*PL/(MAPL8_RGAS*TE)  ! 1000 TAKES TO g m^-3 ; 100 takes mb TO Pa

      if ( ( F > 0.) .and. ( QI > 0. ) ) then
         XIm = (QI/F)*RHO
      else
         XIm = 0.
      end if

      if ( XIm > 0.) then
         LXIm = LOG10(XIm)
      else
         LXIm = 0.0
      end if

    ! Tropical ANVIL produced by resolved scale motions
       VF_A = 128.6 + 53.2*LXIm + 5.5*LXIm**2

    ! Mid-latitude cirrus
       VF_L = 109.0*(XIm**0.16)
       if (abs(XIm) .gt. 0.0) then !Linearisation security
          VF_L = 109.0*(XIm**0.16)
       else
          VF_L = 0.0
       endif
 
    ! Combine the two
       VF = ANVIL*VF_A + LARGESCALE*VF_L

      ! Reduce/increase fall speeds for high/low pressure (NOT in LC98!!! ) 
      ! Assume unmodified they represent situation at 100 mb
      if (WXR > 0.) then
         VF = VF * ( 100./MAX(PL,10.) )**WXR
      endif

      VF = VF/100.

      if ( KH > 2.0 ) then
         VF = 0.01 * VF
      end if

    ! Arctic stratus options
      !if ( KH > 2.0 ) then
      !   VF = 0.01 * VF
      !end if
      !where(PL > 700.)
      !  VF = 0.1*VF
      !endwhere
      !where( (TE >= 250.) .and. (TE<260.))
      !  VF =  ( ( -0.75/10.)*(TE-250.) + 1.0 ) * VF
      !endwhere
      !where(TE >= 260.)
      !  VF = 0.25*VF
      !endwhere

   end SUBROUTINE SETTLE_VEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef _CUDA
!   attributes(device) &
!#endif
   subroutine MARSHPALMQ2(RAIN,PR,DIAM3,NTOTAL,W,VE)

      real(8), intent(in ) :: RAIN,PR     ! in kg m**-2 s**-1, mbar
      real(8), intent(out) :: DIAM3,NTOTAL,W,VE

      real(8) :: RAIN_DAY,LAMBDA,A,B,SLOPR,DIAM1

      real(8), parameter  :: N0   = 0.08  ! # cm**-3

      INTEGER :: IQD

      real(8) :: RX(8) , D3X(8)

      ! Marshall-Palmer sizes at different rain-rates: avg(D^3)

      !RX = (/ 0.   , 5.   , 20.  , 80.  , 320. , 1280., 4*1280., 16*1280. /)  ! rain per in mm/day
      RX(1) = 0.
      RX(2) = 5.
      RX(3) = 20.
      RX(4) = 80.
      RX(5) = 320.
      RX(6) = 1280.
      RX(7) = 4*1280.
      RX(8) = 16*1280.

      !D3X= (/ 0.019, 0.032, 0.043, 0.057, 0.076, 0.102, 0.137  ,  0.183   /)
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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIAM3 = 0.664 * DIAM3  !!  DRYING/EVAP SHOULD PROBABLY GO AS          !!
      !!  D_1.5 == <<D^(3/2)>>^(2/3) NOT AS          !!
      !!  D_3   == <<D^3>>^(1/3)                     !!
      !!  RATIO D_1.5/D_3 =~ 0.66  (JTB 10/17/2002)  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      W      = (2483.8 * DIAM3 + 80.)*SQRT(1000./PR)
      !VE     = 1.0  + 28.0*DIAM3 
      VE     = MAX( 0.99*W/100. , 1.000 )

      DIAM1  = 3.0*DIAM3
      !  Change back to MKS units

      DIAM1  = DIAM1/100.
      DIAM3  = DIAM3/100.
      W      = W/100.
      NTOTAL = NTOTAL*1.0e6

   end subroutine MARSHPALMQ2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef _CUDA
!   attributes(device) &
!#endif
   subroutine MICRO_AA_BB_3(TEMP,PR,Q_SAT,AA,BB)

      real(8), intent(in ) :: TEMP,Q_SAT
      real(8), intent(in ) :: PR
      real(8), intent(out) :: AA,BB

      real(8) :: E_SAT

      real(8), parameter  :: EPSILON =  MAPL8_H2OMW/MAPL8_AIRMW
      real(8), parameter  :: K_COND  =  2.4e-2        ! J m**-1 s**-1 K**-1
      real(8), parameter  :: DIFFU   =  2.2e-5        ! m**2 s**-1

      E_SAT = 100.* PR * Q_SAT /( (EPSILON) + (1.0-(EPSILON))*Q_SAT )  ! (100 converts from mbar to Pa)
   
      AA  = ( GET_ALHX3(TEMP)**2 ) / ( K_COND*MAPL8_RVAP*(TEMP**2) )
      ! AA  = ( MAPL8_ALHL**2 ) / ( K_COND*MAPL8_RVAP*(TEMP**2) )

      BB  = MAPL8_RVAP*TEMP / ( DIFFU*(1000./PR)*E_SAT )

   end subroutine MICRO_AA_BB_3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef _CUDA
!   attributes(device) &
!#endif
   function LDRADIUS3(PL,TE,QCL,NN) RESULT(RADIUS)

      real(8), intent(in) :: TE,PL,NN,QCL
      real(8) :: RADIUS

      real(8) :: MUU,RHO


      RHO = 100.*PL / (MAPL8_RGAS*TE )
      MUU = QCL * RHO                       
      RADIUS = MUU/(NN*RHO_W*(4./3.)*MAPL8_PI)
      RADIUS = RADIUS**(1./3.)    ! Equiv. Spherical Cloud Particle Radius in m


   end function LDRADIUS3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef _CUDA
!   attributes(device) &
!#endif
   function ICE_FRACTION (TEMP) RESULT(ICEFRCT)
      real(8), intent(in) :: TEMP
      real(8)          :: ICEFRCT

      ICEFRCT  = 0.00
      if ( TEMP <= T_ICE_ALL ) then
         ICEFRCT = 1.000
      else if ( (TEMP > T_ICE_ALL) .AND. (TEMP <= T_ICE_MAX) ) then
         ICEFRCT = 1.00 -  ( TEMP - T_ICE_ALL ) / ( T_ICE_MAX - T_ICE_ALL ) 
      end if
      ICEFRCT = MIN(ICEFRCT,1.00)
      ICEFRCT = MAX(ICEFRCT,0.00)

      !!ICEFRCT = ICEFRCT**4
      ICEFRCT = ICEFRCT**ICEFRPWR

   end function ICE_FRACTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device) &
!#endif
   function GET_ALHX3(T) RESULT(ALHX3)

      real(8), intent(in) :: T
      real(8) :: ALHX3

      real(8) :: T_X

      T_X = T_ICE_MAX

      if ( T < T_ICE_ALL ) then
         ALHX3=MAPL8_ALHS
      end if

      if ( T > T_X ) then
         ALHX3=MAPL8_ALHL
      end if

      if ( (T <= T_X) .and. (T >= T_ICE_ALL) ) then
         ALHX3 = MAPL8_ALHS + (MAPL8_ALHL-MAPL8_ALHS)*( T - T_ICE_ALL ) /( T_X - T_ICE_ALL )
      end if

   end function GET_ALHX3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef _CUDA
!   attributes(device) &
!#endif
   real(8) function ICEFRAC(T,T_TRANS,T_FREEZ)

      real(8), intent(in) :: T
      real(8), intent(in),optional :: T_TRANS
      real(8), intent(in),optional :: T_FREEZ

      real(8) :: T_X,T_F

      if (present( T_TRANS )) then 
         T_X = T_TRANS
      else
         T_X = T_ICE_MAX
      endif
      if (present( T_FREEZ )) then 
         T_F = T_FREEZ
      else
         T_F = T_ICE_ALL
      endif


      if ( T < T_F ) ICEFRAC=1.000

      if ( T > T_X ) ICEFRAC=0.000

      if ( T <= T_X .and. T >= T_F ) then 
         ICEFRAC = 1.00 - ( T - T_F ) /( T_X - T_F )
      endif

   end function ICEFRAC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef _CUDA
!   attributes(device) &
!#endif
   SUBROUTINE SUNDQ3_ICE3(TEMP,RATE2,RATE3,TE1, F2, F3)

      real(8), INTENT( IN) :: RATE2,RATE3,TE1

      real(8), INTENT( IN) :: TEMP
      real(8), INTENT(OUT) :: F2, F3

      real(8) :: XX, YY,TE0,TE2,JUMP1  !,RATE2,RATE3,TE1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  Ice - phase treatment totally invented
      !!  Sharp increase in autoconversion in range
      !!  ~~TE1 K ~< T < TE0 K .
      !!  (JTB, 3/25/2003)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      TE0=273.
      !TE1=263.
      TE2=200.
      !RATE2=  10.
      JUMP1=  (RATE2-1.0) / ( ( TE0-TE1 )**0.333 ) 
      !RATE3=  25.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Ice - phase treatment  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      IF ( TEMP .GE. TE0 ) THEN
         F2   = 1.0
         F3   = 1.0
      END IF
      IF ( ( TEMP .GE. TE1 ) .AND. ( TEMP .LT. TE0 ) ) THEN
         F2   = 1.0 + JUMP1 * (( TE0 - TEMP )**0.3333)
         F3   = 1.0
      END IF
      IF ( TEMP .LT. TE1 ) THEN
         F2   = RATE2 + (RATE3-RATE2)*(TE1-TEMP)/(TE1-TE2)
         F3   = 1.0
      END IF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      F2 = MIN(F2,27.0)


   end  subroutine sundq3_ice3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device) &
!#endif
!   subroutine RADCOUPLE(  &
!         TE,              & 
!         PL,              & 
!         CF,              & 
!         AF,              & 
!         QClLS,           & 
!         QCiLS,           & 
!         QClAN,           & 
!         QCiAN,           & 
!         QRN_ALL,         & 
!         QSN_ALL,         & 
!         RAD_QL,          &  
!         RAD_QI,          & 
!         RAD_QR,          & 
!         RAD_QS,          & 
!         RAD_CF,          & 
!         RAD_RL,          & 
!         RAD_RI,          & 
!         CLDVOL2FRC,      &          
!         NN_ANVIL,NN_ICE,NN_WARM,&
!         TEMPOR)
!
!      real(8), intent(in ) :: NN_ANVIL,NN_ICE,NN_WARM 
!      real(8), intent(in ) :: CLDVOL2FRC 
!      real(8), intent(in ) :: TE
!      real(8), intent(in ) :: PL
!      real(8), intent(in ) :: AF,CF, QClAN, QCiAN, QClLS, QCiLS
!      real(8), intent(in ) :: QRN_ALL, QSN_ALL
!      real(8), intent(out) :: RAD_QL,RAD_QI,RAD_QR,RAD_QS,RAD_CF,RAD_RL,RAD_RI
!
!      real(8), intent(in )  :: tempor
!
!      real(8) :: RElAN, REiAN, RElLS, REiLS, QCm, NN, ss, RAD_RI_AN
!      real(8) :: QClANm, QCiANm, QClLSm, QCiLSm, QCtot, AFx
!      real(8) :: rampt, rampu, rampp
!
!      real(8) :: ALPH, POLAR_RL
!
!
!      ! Limits on Radii needed to ensure
!      ! correct behavior of cloud optical
!      ! properties currently calculated in 
!      ! sorad and irrad (1e-6 m = micron)
!
!      POLAR_RL=   5.0e-6  ! 11/09/2007 JTB - COLD low level clouds
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! Adjust Anvil fractions for
!      ! warm clouds
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   !  Needed for PRECIP Repartition
!   !  -----------------------------
!   !  ALPH =  0.05
!   !  SS   =  (260.-TE)/30.
!   !  SS   =  MIN( 1.0 , SS )
!   !  SS   =  MAX( 0.0 , SS )
!   !  SS   =  ALPH + (SS**3) * ( 1.0 - ALPH )
!   !  AFx  =  AF * SS 
!
!      ALPH =  0.1
!      SS   =  (280.-TE)/20.
!      SS   =  MIN( 1.0 , SS )
!      SS   =  MAX( 0.0 , SS )
!      SS   =  ALPH + (SS**3) * ( 1.0 - ALPH )
!      AFx  =  AF * SS * 0.5
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! Total cloud fraction
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      RAD_CF = MIN( CF + AFx, 1.00 )
!
!      ! Total In-cloud liquid
!      if ( RAD_CF > 0. ) then
!         RAD_QL = ( QClLS + QClAN ) / RAD_CF
!      else
!         RAD_QL = 0.0
!      end if
!      RAD_QL = MIN( RAD_QL, 0.01 )
!
!      ! Total In-cloud ice
!      if (  RAD_CF >0. ) then
!         RAD_QI = ( QCiLS + QCiAN ) / RAD_CF
!      else
!         RAD_QI = 0.0
!      end if
!      RAD_QI = MIN( RAD_QI, 0.01 )
!
!
!      ! Total In-cloud precipitation
!      if (  RAD_CF >0. ) then
!         RAD_QR = ( QRN_ALL ) / RAD_CF
!         RAD_QS = ( QSN_ALL ) / RAD_CF
!      else
!         RAD_QR = 0.0
!         RAD_QS = 0.0
!      end if
!      RAD_QR = MIN( RAD_QR, 0.01 )
!      RAD_QS = MIN( RAD_QS, 0.01 )
!
!      if (PL < 150. ) then
!         RAD_RI = MAX_RI
!      end if
!      if (PL >= 150. ) then
!         RAD_RI = MAX_RI*150./PL
!      end if
!
!      !! weigh in a separate R_ice for Anvil Ice according to
!      !
!      !       R_net_eff = (q_anv + q_ls) / ( q_anv/R_ice_anv + q_ls/R_ice_ls )
!      !-------------------------------------------------------------------------
!      RAD_RI_AN  =  RAD_RI ! 40.0e-6   ! MIN_RI 
!
!      if ( ( QCiLS + QCiAN ) > 0.0 ) then
!         RAD_RI_AN  = ( QCiLS + QCiAN ) / ( (QCiLS/RAD_RI) + (QCiAN/RI_ANV) )
!      end if
!
!      RAD_RI = MIN( RAD_RI, RAD_RI_AN )
!
!      RAD_RI = MAX( RAD_RI, MIN_RI )
!
!      ! Implement ramps for gradual change in effective radius
!      if (PL < 300. ) then
!         RAD_RL = 21.e-6
!      end if
!      if (PL >= 300. ) then
!         RAD_RL = 21.e-6*300./PL
!      end if
!      RAD_RL = MAX( RAD_RL, 10.e-6 )
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! Thicken low high lat clouds
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      ! NOTE: Due to how tempor is calculated, it is now calculated in the
!      ! GridComp and passed into progno_cloud
!
!      if ( PL .GE. 775.  .AND. TE .LE.  275. .AND. (tempor.eq.1.) ) then
!         RAD_RL = max(min(-0.1 * PL + 87.5, 10.),5.)*1.e-6
!      end if
!      if ( PL .GE. 825.  .AND. TE .LE.  282. .AND. (tempor.eq.1.) ) then
!         RAD_RL = max(0.71 * TE - 190.25, 5.)*1.e-6
!      end if
!      if ( PL .GE. 775.  .AND. PL .LT. 825. .AND. TE .LE.  282. .AND. TE .GT. 275. .AND. (tempor.eq.1.) ) then
!         RAD_RL = min(-0.1*PL + 0.71 * TE - 107.75, 10.)*1.e-6
!      end if
!      if ( PL .GE. 825.  .AND. TE .LE.  275. .AND. (tempor.eq.1.) ) then
!         RAD_RL = 5.*1.e-6
!      end if
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! Thin low tropical clouds
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
!         RAD_RL = min(2.2 * TE - 617., 21.)*1.e-6
!      end if
!      if ( PL .GE. 925.  .AND. TE .GE.  290. ) then
!         RAD_RL = min(0.44 * PL - 397., 21.)*1.e-6
!      end if
!      if ( PL .GE. 925.  .AND. PL .LT. 950. .AND. TE .GT.  285. .AND. TE .LT. 290.) then
!         RAD_RL = max(min(0.44*PL + 2.2 * TE - 1035., 21.),10.)*1.e-6
!      end if
!      if ( PL .GE. 950.  .AND. TE .GE.  290. ) then
!         RAD_RL = 21.*1.e-6
!      end if
!
!      if ( RAD_CF < 1.e-5 ) then
!         RAD_QL = 0.
!         RAD_QI = 0.
!         RAD_CF = 0.
!         RAD_QR = 0.
!         RAD_QS = 0.
!      end if
!
!   end subroutine RADCOUPLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef _CUDA
!   attributes(device) function QSAT(TL,PL,PASCALS)
!
!     real(8),              intent(IN) :: TL, PL
!     logical, optional, intent(IN) :: PASCALS
!     real(8) :: QSAT
!
!     real(8) :: URAMP, DD, QQ, TI, DQ, PP
!     integer :: IT
! 
!     URAMP = TMIX
!
!     if (present(PASCALS)) then 
!        if (PASCALS) then
!           PP = PL
!        else
!           PP = PL*100.
!        end if
!     else
!        PP = PL*100.
!     end if
! 
!     TI = TL - ZEROC
!
!     if    (TI <= URAMP) then
!        QSAT  =  QSATICE0(TL,PP,DQ)
!     elseif(TI >= 0.0  ) then
!        QSAT  =  QSATLQU0(TL,PP,DQ)
!     else
!        QSAT  =  QSATICE0(TL,PP,DQ)
!        QQ    =  QSATLQU0(TL,PP,DQ)
!        TI    =  TI/URAMP
!        QSAT  =  TI*(QSAT - QQ) +  QQ
!     end if
!
!   end function QSAT
!
!   attributes(device) function DQSAT(TL,PL,QSAT,PASCALS)
!
!      real(8),              intent(IN) :: TL, PL
!      real(8),              intent(OUT):: QSAT
!      logical, optional, intent(IN ):: PASCALS
!      real(8) :: DQSAT
!
!      real(8) :: URAMP, TT, WW, DD, DQQ, QQ, TI, DQI, QI, PP, DQ
!      integer :: IT
!
!      URAMP = TMIX
!
!      if (present(PASCALS)) then 
!         if (PASCALS) then
!            PP = PL
!         else
!            PP = PL*100.
!         end if
!      else
!         PP = PL*100.
!      end if
!
!      TI = TL - ZEROC
!
!      if    (TI <= URAMP) then
!         QQ  = QSATICE0(TL,PP,DQ)
!         QSAT  = QQ
!         DQSAT = DQ
!      elseif(TI >= 0.0  ) then
!         QQ  = QSATLQU0(TL,PP,DQ)
!         QSAT  = QQ
!         DQSAT = DQ
!      else
!         QQ  = QSATLQU0(TL,PP,DQQ)
!         QI  = QSATICE0(TL,PP,DQI)
!         TI  = TI/URAMP
!         DQSAT = TI*(DQI - DQQ) + DQQ
!         QSAT  = TI*(QI - QQ) +  QQ
!      end if
!
!   end function DQSAT
!
!   attributes(device) function QSATLQU0(TL,PL,DQ) result(QS)
!
!      real(8), intent(IN) :: TL
!      real(8), intent(IN) :: PL
!      real(8), intent(OUT):: DQ
!      real(8) :: QS
!
!      real(8) :: TI,W
!      real(8) :: DD
!      real(8) :: TT
!      real(8) :: DDQ
!      integer :: IT
!
!      integer, parameter :: TYPE = 1
!
!#define TX TL
!#define PX PL
!#define EX QS
!#define DX DQ
!
!
!   if    (TX<TMINLQU) then
!      TI = TMINLQU
!   elseif(TX>TMAXTBL) then
!      TI = TMAXTBL
!   else
!      TI = TX
!   end if
!
!#include "esatlqu.code"
!
!   if    (TX<TMINLQU) then
!      DDQ = 0.0
!   elseif(TX>TMAXTBL) then
!      DDQ = 0.0
!   else
!      if(PX>EX) then
!         DD = EX
!         TI = TX + DELTA_T
!#include "esatlqu.code"
!         DDQ = EX-DD
!         EX  = DD
!      endif
!   end if
!
!   if(PX > EX) then
!      DD = ESFAC/(PX - (1.0-ESFAC)*EX)
!      EX = EX*DD
!      DX = DDQ*ERFAC*PX*DD*DD
!   else
!      EX = MAX_MIXING_RATIO
!      DX = 0.0
!   end if
!
!#undef  DX
!#undef  TX
!#undef  EX
!#undef  PX
!
!      return
!   end function QSATLQU0
!
!   attributes(device) function QSATICE0(TL,PL,DQ) result(QS)
!
!      real(8), intent(IN) :: TL
!      real(8), intent(IN) :: PL
!      real(8), intent(OUT):: DQ
!      real(8) :: QS
!
!      real(8) :: TI,W
!      real(8) :: DD
!      real(8) :: TT
!      real(8) :: DDQ
!      integer :: IT
!
!      integer, parameter :: TYPE = 1
!
!#define TX TL
!#define PX PL
!#define EX QS
!#define DX DQ
!
!
!   if    (TX<TMINICE) then
!      TI = TMINICE
!   elseif(TX>ZEROC  ) then
!      TI = ZEROC
!   else
!      TI = TX
!   end if
!
!#include "esatice.code"
!
!   if    (TX<TMINICE) then
!      DDQ = 0.0
!   elseif(TX>ZEROC  ) then
!      DDQ = 0.0
!   else
!      if(PX>EX) then
!         DD = EX
!         TI = TX + DELTA_T
!#include "esatice.code"
!         DDQ = EX-DD
!         EX  = DD
!      endif
!   end if
!
!   if(PX > EX) then
!      DD = ESFAC/(PX - (1.0-ESFAC)*EX)
!      EX = EX*DD
!      DX = DDQ*ERFAC*PX*DD*DD
!   else
!      EX = MAX_MIXING_RATIO
!      DX = 0.0
!   end if
!
!#undef  DX
!#undef  TX
!#undef  EX
!#undef  PX
!
!         return
!   end function QSATICE0
!
!#endif

subroutine PRE_PROGNO_CLOUD (IM,JM,LM,TH1,PK,PLO,PKE,CNV_PLE,&
                             QST3,DZET,QDDF3,&
                             CNV_FRACTION,CLDPARAMS)

 implicit none

 !Inputs
 integer, intent(in) :: IM,JM,LM
 real(8), intent(in )  , dimension(IM,JM,0:LM) :: PKE,CNV_PLE
 real(8), intent(in )  , dimension(IM,JM,LM)   :: TH1, PK, PLO
 type(CLDPARAM_TYPE), intent(in)            :: CLDPARAMS

 !Outputs
 real(8), intent(out)  , dimension(IM,JM,LM)   :: QST3,DZET,QDDF3

 !Inouts
 real(8), intent(inout), dimension(IM,JM)      :: CNV_FRACTION
 
 !Locals
 integer :: i,j,k
 real(8), dimension(IM,JM,LM)   :: TEMP,MASS,DQST3
 real(8), dimension(IM,JM)      :: VMIP
 real(8), dimension(IM,JM,LM+1) :: ZET
 
 TEMP = TH1*PK

 DQST3 = 0.0
 QST3 = 0.0
 do I=1,IM
    do J=1,JM
       do k = 1,lm
          call DQSATSCApert(DQST3(i,j,k),QST3(i,j,k),TEMP(i,j,k),PLO(i,j,k))
       end do
    end do
 end do

 DZET = 0.0
 DZET(:,:,1:LM) = TH1(:,:,1:LM) * (PKE(:,:,1:LM) - PKE(:,:,0:LM-1)) * MAPL8_CP/MAPL8_GRAV
 MASS(:,:,1:LM) = ( CNV_PLE(:,:,1:LM) - CNV_PLE(:,:,0:LM-1) )*100./MAPL8_GRAV
 ZET(:,:,LM+1) = 0.0
 DO K = LM, 1, -1
    ZET(:,:,K) = ZET(:,:,K+1)+DZET(:,:,K)
 END DO

 QDDF3 = 0.0
 WHERE ( ZET(:,:,1:LM) < 3000. )
    QDDF3 = -( ZET(:,:,1:LM)-3000. ) * ZET(:,:,1:LM) * MASS
 ELSEWHERE
    QDDF3 = 0.
 END WHERE

 DO J=1,JM
    DO I=1,IM
       VMIP(I,J) = SUM(QDDF3(I,J,:))
    END DO
 END DO

 DO K = 1,LM
    QDDF3(:,:,K) = QDDF3(:,:,K) / VMIP
 END DO

 if(CLDPARAMS%MOVE2RAS.eq.0.) then
    CNV_FRACTION = 0.0 
 endif

end subroutine PRE_PROGNO_CLOUD

end module cloudnew
