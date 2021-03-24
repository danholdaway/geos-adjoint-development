! $Id: GEOS_MoistGridComp.F90,v 1.109.6.1 2014-05-15 16:31:13 drholdaw Exp $

! VERIFY_ and RETURN_ macros for error handling.

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_Moist -- A Module to compute moist processes, including convection,
!   large-scale condensation and precipitation and cloud parameters.

! !INTERFACE:

module GEOS_MoistGridCompMod

! !USES:

  use RAS       ! using module that contains ras code
 
#ifndef _CUDA
  use CLOUDNEW, only: PROGNO_CLOUD, ICE_FRACTION, T_CLOUD_CTL
#else
  use CLOUDNEW, only: &
        ! Subroutines
        PROGNO_CLOUD, ICE_FRACTION, &
        ! Derived Data Types
        T_CLOUD_CTL, &
        ! Inputs
        PP_DEV, EXNP_DEV, PPE_DEV, KH_DEV, FRLAND_DEV, &
        RMFDTR_DEV, QLWDTR_DEV, U_DEV, V_DEV, QST3_DEV, &
        DZET_DEV, QDDF3_DEV, TEMPOR_DEV, &
        ! Inoutputs
        TH_DEV, Q_DEV, QRN_CU_DEV, CNV_UPDFRC_DEV, QLW_LS_DEV, &
        QLW_AN_DEV, QIW_LS_DEV, QIW_AN_DEV, ANVFRC_DEV, CLDFRC_DEV, &
        ! Outputs
        RAD_CLDFRC_DEV, RAD_QL_DEV, RAD_QI_DEV, RAD_QR_DEV, RAD_QS_DEV, &
        CLDREFFL_DEV, CLDREFFI_DEV, PRELS_DEV, PRECU_DEV, PREAN_DEV, &
        LSARF_DEV, CUARF_DEV, ANARF_DEV, SNRLS_DEV, SNRCU_DEV, SNRAN_DEV, &
        ! Working arrays 
        PFL_CN_DEV, PFI_CN_DEV, PFL_AN_DEV, &
        PFI_AN_DEV, PFL_LS_DEV, PFI_LS_DEV, &
        ! Diagnostics
        RHX_DEV, &
        REV_LS_DEV, REV_AN_DEV, REV_CN_DEV, &
        RSU_LS_DEV, RSU_AN_DEV, RSU_CN_DEV, &
        ACLL_CN_DEV, ACIL_CN_DEV, ACLL_AN_DEV, &
        ACIL_AN_DEV, ACLL_LS_DEV, ACIL_LS_DEV, &
        PDFL_DEV, PDFI_DEV, FIXL_DEV, FIXI_DEV, &
        AUT_DEV, EVAPC_DEV, SDM_DEV, &
        SUBLC_DEV, FRZ_TT_DEV, DCNVL_DEV, DCNVI_DEV, &
        ALPHT_DEV, ALPH1_DEV, ALPH2_DEV, &
        CFPDF_DEV, RHCLR_DEV, DQRL_DEV, FRZ_PP_DEV, &
        VFALLICE_AN_DEV, VFALLICE_LS_DEV, &
        VFALLWAT_AN_DEV, VFALLWAT_LS_DEV, &
        VFALLSN_AN_DEV, VFALLSN_LS_DEV, &
        VFALLSN_CN_DEV, VFALLRN_AN_DEV, &
        VFALLRN_LS_DEV, VFALLRN_CN_DEV, &
        !CFPDFX is no longer calculated in PROGNO_CLOUD, but still an export
        !CFPDFX_DEV, &

        ! Constants
        ! PHYSPARAMS Constants are loaded into constant memory
        CNV_BETA, ANV_BETA, LS_BETA, RH00, C_00, LWCRIT, C_ACC, &
        C_EV_R, C_EV_S, CLDVOL2FRC, RHSUP_ICE, SHR_EVAP_FAC, MIN_CLD_WATER, &
        CLD_EVP_EFF, NSMAX, LS_SDQV2, LS_SDQV3, LS_SDQVT1, ANV_SDQV2, &
        ANV_SDQV3, ANV_SDQVT1, ANV_TO_LS, N_WARM, N_ICE, N_ANVIL, &
        N_PBL, DISABLE_RAD, ANV_ICEFALL_C, LS_ICEFALL_C, REVAP_OFF_P, CNVENVFC, &
        WRHODEP, T_ICE_ALL, CNVICEPARAM, ICEFRPWR, CNVDDRFC, ANVDDRFC, &
        LSDDRFC, TANHRHCRIT, MINRHCRIT, MAXRHCRIT, TURNRHCRIT, MAXRHCRITLAND, &
        FR_LS_WAT, FR_LS_ICE, FR_AN_WAT, FR_AN_ICE, MIN_RL, MIN_RI, MAX_RL, &
        MAX_RI, RI_ANV, PDFFLAG
  use cudafor
#endif

  use DDF
 
  use ESMF
  use MAPL_Mod
  use GEOS_UtilsMod
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !DESCRIPTION:
! 
!   {\tt GEOS\_MoistGridCompMod} implements moist processes in GEOS-5. These
!   include all processes that involve phase changes in the atmosphere, such
!   as large-scale condensation, convective clouds, and all rain and cloud
!   formation. It's state consists of water vapor, various types of condensate,
!   and fractions of various cloud types.
!

!EOP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code
    
! !DESCRIPTION:  {\tt GEOS\_MoistGridCompMod} uses the default Initialize and Finalize 
!                services, but registers its own Run method.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp    ), pointer   :: STATE 
    type (ESMF_Config          )            :: CF

    integer      :: RFRSHINT
    integer      :: AVRGNINT
    real         :: DT
    
!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,  &
                                      RC=STATUS)
    VERIFY_(STATUS)
    

! Get the configuration from the component
!-----------------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

    call ESMF_ConfigGetAttribute ( CF, DT, Label="RUN_DT:",                                   RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( CF, RFRSHINT, Label="REFRESH_INTERVAL:",  default=nint(DT), RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( CF, AVRGNINT, Label='AVERAGING_INTERVAL:',default=RFRSHINT, RC=STATUS)
    VERIFY_(STATUS)

!    call MAPL_GetResource ( STATE, RFRSHINT, Label="REFRESH_INTERVAL:",  default=nint(DT), RC=STATUS)
!    VERIFY_(STATUS)

!    call MAPL_GetResource ( STATE, AVRGNINT, Label='AVERAGING_INTERVAL:',default=RFRSHINT, RC=STATUS)
!    VERIFY_(STATUS)

! !INTERNAL STATE:

!BOS

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME = 'Q',                                          &
        LONG_NAME  = 'specific_humidity',                          &
        UNITS      = 'kg kg-1',                                    &
        FRIENDLYTO = 'DYNAMICS:TURBULENCE:CHEMISTRY:ANALYSIS',    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME = 'QLLS',                                       &
        LONG_NAME  = 'mass_fraction_of_large_scale_cloud_liquid_water', &
        UNITS      = 'kg kg-1',                                    &
        FRIENDLYTO = 'DYNAMICS:TURBULENCE',                       &
        !!FRIENDLYTO = 'TURBULENCE',                       &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 
     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME = 'QLCN',                                       &
        LONG_NAME  = 'mass_fraction_of_convective_cloud_liquid_water', &
        UNITS      = 'kg kg-1',                                    &
        FRIENDLYTO = 'DYNAMICS',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 
     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME = 'CLLS',                                       &
        LONG_NAME  = 'large_scale_cloud_area_fraction',            &
        UNITS      = '1',                                          &
        FRIENDLYTO = 'DYNAMICS',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME = 'CLCN',                                       &
        LONG_NAME  = 'convective_cloud_area_fraction',             &
        UNITS      = '1',                                          &
        FRIENDLYTO = 'DYNAMICS',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME = 'QILS',                                       &
        LONG_NAME  = 'mass_fraction_of_large_scale_cloud_ice_water', &
        UNITS      = 'kg kg-1',                                    &
        !!FRIENDLYTO = 'TURBULENCE',                       &
        FRIENDLYTO = 'DYNAMICS:TURBULENCE',                       &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME = 'QICN',                                       &
        LONG_NAME  = 'mass_fraction_of_convective_cloud_ice_water', &
        UNITS      = 'kg kg-1',                                    &
        FRIENDLYTO = 'DYNAMICS',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 

! !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME = 'PLE',                                         &
        LONG_NAME  = 'air_pressure',                                &
        UNITS      = 'Pa',                                          &
        DIMS       = MAPL_DimsHorzVert,                            &
        VLOCATION  = MAPL_VLocationEdge,                           &
        AVERAGING_INTERVAL = AVRGNINT,                             &
        REFRESH_INTERVAL   = RFRSHINT,                             &
                                                        RC=STATUS  )
     VERIFY_(STATUS)                                                                          

     call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME = 'PREF',                                       &
        LONG_NAME  = 'reference_air_pressure',                     &
        UNITS      = 'Pa',                                         &
        DIMS       = MAPL_DimsVertOnly,                            &
        VLOCATION  = MAPL_VLocationEdge,                           &
        AVERAGING_INTERVAL = AVRGNINT,                             &
        REFRESH_INTERVAL   = RFRSHINT,                             &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME = 'KH',                                         &
        LONG_NAME  = 'scalar_diffusivity',                         &
        UNITS      = 'm+2 s-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                            &
        VLOCATION  = MAPL_VLocationEdge,                           &
        AVERAGING_INTERVAL = AVRGNINT,                             &
        REFRESH_INTERVAL   = RFRSHINT,                             &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'TH',                                        &
        LONG_NAME  = 'potential_temperature',                     &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        AVERAGING_INTERVAL = AVRGNINT,                            &
        REFRESH_INTERVAL   = RFRSHINT,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'U',                                         &
        LONG_NAME  = 'eastward_wind',                             &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        AVERAGING_INTERVAL = AVRGNINT,                            &
        REFRESH_INTERVAL   = RFRSHINT,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'V',                                         &
        LONG_NAME  = 'northward_wind',                            &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        AVERAGING_INTERVAL = AVRGNINT,                            &
        REFRESH_INTERVAL   = RFRSHINT,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'TS',                                        &
        LONG_NAME  = 'surface temperature',                       &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
        AVERAGING_INTERVAL = AVRGNINT,                            &
        REFRESH_INTERVAL   = RFRSHINT,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'FRLAND',                                    &
        LONG_NAME  = 'areal_land_fraction',                       &
        UNITS      = '1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                        &
        AVERAGING_INTERVAL = AVRGNINT,                            &
        REFRESH_INTERVAL   = RFRSHINT,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'FROCEAN',                                   &
        LONG_NAME  = 'areal_ocean_fraction',                      &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
        AVERAGING_INTERVAL = AVRGNINT,                            &
        REFRESH_INTERVAL   = RFRSHINT,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

! These bundles should be changed when we merge w/ the head.

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'MTR',                                        &
        LONG_NAME  = 'tracers_for_moist',                          &
        UNITS      = 'X',                                          &
        DATATYPE   = MAPL_BundleItem,                             &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = .false.,                                     &
                                                        RC=STATUS )
     VERIFY_(STATUS)                                                                           

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'KPBL',                                       &
        LONG_NAME  = 'planetary_boundary_layer_level',             &
        UNITS      = '1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                            &
        VLOCATION  = MAPL_VLocationNone,                           &
        AVERAGING_INTERVAL = AVRGNINT,                             &
        REFRESH_INTERVAL   = RFRSHINT,                             &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

! !EXPORT STATE:
 
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'QCTOT',                                      &
        LONG_NAME  = 'mass_fraction_of_total_cloud_water',         &
        UNITS      = 'kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'QLTOT',                                      &
        LONG_NAME  = 'grid_box_mass_fraction_of_cloud_liquid_water',        &
        UNITS      = 'kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'QITOT',                                      &
        LONG_NAME  = 'grid_box_mass_fraction_of_cloud_ice_water',           &
        UNITS      = 'kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
                                                                               
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'QRTOT',                                      &
        LONG_NAME  = 'mass_fraction_of_falling_rain',              &
        UNITS      = 'kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
                                                                               
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'QSTOT',                                      &
        LONG_NAME  = 'mass_fraction_of_falling_snow',              &
        UNITS      = 'kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
                                                                               
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'MTRI',                                       &
        LONG_NAME  = 'tracer_tendencies_due_to_moist',             &
        UNITS      = 'X s-1',                                      &
!        TYPE       = MAPL_BUNDLE,                                &
        DATATYPE   = MAPL_BundleItem,                             &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                        RC=STATUS )
    VERIFY_(STATUS)                                                                           

    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'DTHDT ',                                     &
         LONG_NAME = 'pressure_weighted_potential_temperature_tendency_due_to_moist',&
         UNITS     = 'Pa K s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                           &
         VLOCATION = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'DTDTFRIC',                                   &
         LONG_NAME = 'pressure_weighted_temperature_tendency_due_to_moist_friction',&
         UNITS     = 'Pa K s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                           &
         VLOCATION = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'DQDT  ',                                     &
         LONG_NAME = 'specific_humidity_tendency_due_to_moist',    &
         UNITS     = 'kg kg-1 s-1',                                &
         DIMS      = MAPL_DimsHorzVert,                           &
         VLOCATION = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DUDT  ',                                      &
         LONG_NAME = 'zonal_wind_tendency_due_to_moist',            &
         UNITS     = 'm s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &                  
         SHORT_NAME = 'DVDT  ',                                      &
         LONG_NAME = 'meridional_wind_tendency_due_to_moist',       &
         UNITS     = 'm s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DTHDTCN',                                     &
         LONG_NAME = 'potential_temperature_tendency_due_to_convection',&
         UNITS     = 'K s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQDTCN',                                      &
         LONG_NAME = 'specific_humidity_tendency_due_to_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQLDT ',                                      &
         LONG_NAME = 'total_liq_water_tendency_due_to_moist',       &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'DQIDT ',                                      &
         LONG_NAME = 'total_ice_water_tendency_due_to_moist',       &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQCDTCN ',                                    &
         LONG_NAME = 'condensate_tendency_due_to_convection',       &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'CNV_DQLDT ',                                  &
         LONG_NAME = 'convective_condensate_source',                &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME= 'CNV_PRC3 ',                                   &
         LONG_NAME = 'convective_precipitation_from_RAS',           &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQRL   ',                                     &
         LONG_NAME = 'large_scale_rainwater_source',                &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'DQRC   ',                                     &
         LONG_NAME = 'convective_rainwater_source',                 &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
    

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'CNV_MF0',                                     &
         LONG_NAME = 'cloud_base_mass_flux',                        &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddExportSpec(GC,                               & 
         SHORT_NAME = 'CNV_MFD',                                     & 
         LONG_NAME = 'detraining_mass_flux',                        &
         UNITS     = 'kg m-2 s-1',                                  &    
         DIMS      = MAPL_DimsHorzVert,                            &  
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &                  
         SHORT_NAME = 'CNV_MFC',                                     & 
         LONG_NAME = 'cumulative_mass_flux',                        &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                           &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'CNV_FREQ',                                    & 
         LONG_NAME = 'convective_frequency',                        &
         UNITS     = 'fraction',                                    &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'CNV_BASEP',                                   & 
         LONG_NAME = 'pressure_at_convective_cloud_base',           &
         UNITS     = 'Pa',                                          &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'CNV_TOPP',                                    & 
         LONG_NAME = 'pressure_at_convective_cloud_top',            &
         UNITS     = 'Pa',                                          &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
                                                        RC=STATUS  )
    VERIFY_(STATUS)



    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CNV_UPDF',                                    &
         LONG_NAME = 'updraft_areal_fraction',                      &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )

    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CNV_CVW',                                     &
         LONG_NAME = 'updraft_vertical_velocity',                   &
         UNITS     = 'hPa s-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )

    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='CNV_QC',                                      &
         LONG_NAME ='grid_mean_convective_condensate',             &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='RL',                                          & 
         LONG_NAME ='liquid_cloud_particle_effective_radius',      &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RI',                                          & 
         LONG_NAME = 'ice_phase_cloud_particle_effective_radius',   &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RR',                                          & 
         LONG_NAME = 'falling_rain_particle_effective_radius',      &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RS',                                          & 
         LONG_NAME  = 'falling_ice_particle_effective_radius',       &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='CLDNCCN',                                     & 
         LONG_NAME ='number_concentration_of_cloud_particles',     &
         UNITS     ='m-3',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QSATI'  ,                                     & 
         LONG_NAME = 'saturation_spec_hum_over_ice',                &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QSATL'  ,                                     & 
         LONG_NAME = 'saturation_spec_hum_over_liquid',             &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ALPHT'  ,                                     & 
         LONG_NAME = 'pdf_spread_for_condensation_over_qsat_total', &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ALPH1'  ,                                     & 
         LONG_NAME ='pdf_spread_for_condensation_over_qsat_term1', &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ALPH2'  ,                                     & 
         LONG_NAME = 'pdf_spread_for_condensation_over_qsat_term2', &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CFPDFX'  ,                                    & 
         LONG_NAME = 'cloud_fraction_internal_in_PDF_scheme',       &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RHCLR'  ,                                     & 
         LONG_NAME = 'RH_clear_sky',                                &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CFPDF'  ,                                     & 
         LONG_NAME = 'cloud_fraction_after_PDF',                    &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'FCLD'  ,                                      & 
         LONG_NAME = 'cloud_fraction_for_radiation',                &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='QV',                                          & 
         LONG_NAME ='water_vapor_for_radiation',                   &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QL',                                          & 
         LONG_NAME = 'in_cloud_cloud_liquid_for_radiation',                  &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QI',                                          & 
         LONG_NAME = 'in_cloud_cloud_ice_for_radiation',                     &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QR',                                          & 
         LONG_NAME = 'Falling_rain_for_radiation',                  &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QS',                                          & 
         LONG_NAME = 'Falling_snow_for_radiation',                  &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='LS_PRCP',                                     & 
         LONG_NAME ='nonanvil_large_scale_precipitation',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'AN_PRCP',                                     & 
         LONG_NAME = 'anvil_precipitation',                         &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='CN_PRCP',                                     & 
         LONG_NAME ='convective_precipitation',                    &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ER_PRCP',                                     & 
         LONG_NAME = 'spurious_rain_from_RH_cleanup',          &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'FILLNQV',                                     & 
         LONG_NAME = 'filling_of_negative_Q',          &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PGENTOT',                                     & 
         LONG_NAME = 'Total_column_production_of_precipitation',    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PREVTOT',                                     & 
         LONG_NAME = 'Total_column_re-evap/subl_of_precipitation',    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'LS_ARF',                                      & 
         LONG_NAME = 'areal_fraction_of_nonanvil_large_scale_showers',&
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'AN_ARF',                                      & 
         LONG_NAME = 'areal_fraction_of_anvil_showers',             &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CN_ARF',                                      & 
         LONG_NAME = 'areal_fraction_of_convective_showers',        &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'SNO',                                         & 
         LONG_NAME = 'snowfall',                                    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PCU',                                         & 
         LONG_NAME = 'convective_rainfall',                         &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PLS',                                         & 
         LONG_NAME = 'large_scale_rainfall',                        &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TPREC',                                       & 
         LONG_NAME ='total_precipitation',                         &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='HOURNORAIN',                                  & 
         LONG_NAME ='time-during_an_hour_with_no_precipitation',   &
         UNITS     ='s',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TPW',                                         & 
         LONG_NAME ='total_precipitable_water',                    &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CCWP',                                        & 
         LONG_NAME ='grid_mean_conv_cond_water_path_diagnostic',   &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CWP',                                         & 
         LONG_NAME ='condensed_water_path',                        &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='LWP',                                         & 
         LONG_NAME ='liquid_water_path',                           &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='IWP',                                         & 
         LONG_NAME ='ice_water_path',                              &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='BYNCY',                                       & 
         LONG_NAME ='buoyancy_of surface_parcel',                  &
         UNITS     ='m s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            & 
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CAPE',                                        & 
         LONG_NAME ='cape_for_surface_parcel',                     &
         UNITS     ='J m-2',                                       &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='INHB',                                        & 
         LONG_NAME ='inhibition_for_surface_parcel',               &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)
     
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVQ0',                                        & 
         LONG_NAME ='Total_Water_Substance_Before',                &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVQ1',                                        & 
         LONG_NAME ='Total_Water_Substance_After',                 &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DCPTE',                                        & 
         LONG_NAME ='Total_VI_DcpT',                         &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVE0',                                        & 
         LONG_NAME ='Total_VI_MSE_Before',                         &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVE1',                                        & 
         LONG_NAME ='Total_VI_MSE_After',                          &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVEX',                                        & 
         LONG_NAME ='Total_VI_MSE_Somewhere',                      &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZPBLCN',                                      & 
         LONG_NAME ='boundary_layer_depth',                        &
         UNITS     ='m'   ,                                        &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZLCL',                                        & 
         LONG_NAME ='lifting_condensation_level',                  &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZLFC',                                        & 
         LONG_NAME ='level_of_free_convection',                    &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZCBL',                                        & 
         LONG_NAME ='height_of_cloud_base_layer',                  &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='MXDIAM',                                      & 
         LONG_NAME ='diameter_of_largest_RAS_plume',               &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RASTIME',                                     & 
         LONG_NAME ='timescale_for_deep_RAS_plumes',               &
         UNITS     ='s'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RASPBLQ',                                     & 
         LONG_NAME ='sqrt_of_integral_KH_dz',                      &
         UNITS     ='(m+3 s-1)+1/2',                               &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ENTLAM',                                      &
         LONG_NAME ='entrainment parameter',                       &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

#if 0 
! taken out since they are now friendly to dynamics
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME ='QLCN',                                       &
        LONG_NAME  ='mass_fraction_of_convective_cloud_liquid_water', &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME ='QICN',                                       &
        LONG_NAME  ='mass_fraction_of_convective_cloud_ice_water', &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME ='CLLS',                                       &
        LONG_NAME  ='large_scale_cloud_area_fraction',            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME ='CLCN',                                       &
        LONG_NAME  ='convective_cloud_area_fraction',             &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
#endif


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RH1',                                         & 
         LONG_NAME ='relative_humidity_before_moist',              &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RH2',                                         & 
         LONG_NAME ='relative_humidity_after_moist',               &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

!Outputs to give model trajectory in the moist TLM/ADJ
    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QILS_tlad',                                   & 
         LONG_NAME ='total_cloud_liquid_ice_before_moist',         &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QLLS_tlad',                                   & 
         LONG_NAME ='total_anvil_liquid_ice_before_moist',         &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QICN_tlad',                                   & 
         LONG_NAME ='total_cloud_liquid_water_before_moist',       &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QLCN_tlad',                                   & 
         LONG_NAME ='total_anvil_liquid_water_before_moist',       &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CFLS_tlad',                                   & 
         LONG_NAME ='large_scale_cloud_fraction_before_moist',     &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CFCN_tlad',                                   & 
         LONG_NAME ='convective_cloud_fraction_before_moist',      &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KCBL_tlad',                                   & 
         LONG_NAME ='KCBL_before_moist',                           &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='TSUR_tlad',                                   & 
         LONG_NAME ='surface_temp_before_moist',                   &
         UNITS     ='K',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QLS_tlad',                                   & 
         LONG_NAME ='total_cloud_liquid_ice_before_moist',         &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QCN_tlad',                                   & 
         LONG_NAME ='total_cloud_liquid_water_before_moist',       &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KHu_tlad',                                    & 
         LONG_NAME ='upper_index_where_Kh_greater_than_2',         &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KHl_tlad',                                    & 
         LONG_NAME ='lower_index_where_Kh_greater_than_2',         &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    !End outputs for trajectory of moist TLM/ADJ

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RHX',                                         & 
         LONG_NAME ='relative_humidity_after_PDF',               &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REVSU_CN',                                    & 
         LONG_NAME ='evap_subl_of_convective_precipitation',       &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REVSU_LSAN',                                    & 
         LONG_NAME ='evap_subl_of_non_convective_precipitation',       &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_CN',                                      & 
         LONG_NAME ='evaporation_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_AN',                                      & 
         LONG_NAME ='evaporation_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_LS',                                      & 
         LONG_NAME ='evaporation_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_CN',                                      & 
         LONG_NAME ='sublimation_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_AN',                                      & 
         LONG_NAME ='sublimation_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_LS',                                      & 
         LONG_NAME ='sublimation_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACR_TOT',                                      & 
         LONG_NAME ='total_accretion_of__precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRLL_CN',                                      & 
         LONG_NAME ='liq_liq_accretion_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRLL_AN',                                      & 
         LONG_NAME ='liq_liq_accretion_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRLL_LS',                                      & 
         LONG_NAME ='liq_liq_accretion_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRIL_CN',                                      & 
         LONG_NAME ='ice_liq_accretion_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRIL_AN',                                      & 
         LONG_NAME ='ice_liq_accretion_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRIL_LS',                                      & 
         LONG_NAME ='ice_liq_accretion_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_CN',                                      & 
         LONG_NAME ='3D_flux_of_ice_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_AN',                                      & 
         LONG_NAME ='3D_flux_of_ice_anvil_precipitation',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_LS',                                      & 
         LONG_NAME ='3D_flux_of_ice_nonanvil_large_scale_precipitation',&
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_LSAN',                                      & 
         LONG_NAME ='3D_flux_of_ice_nonconvective_precipitation'  ,&
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_CN',                                      & 
         LONG_NAME ='3D_flux_of_liquid_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_AN',                                      & 
         LONG_NAME ='3D_flux_of_liquid_anvil_precipitation',          &
         UNITS     ='kg m-2 s-1',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_LS',                                      & 
         LONG_NAME ='3D_flux_of_liquid_nonanvil_large_scale_precipitation',&
         UNITS     ='kg m-2 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_LSAN',                                      & 
         LONG_NAME ='3D_flux_of_liquid_nonconvective_precipitation'  ,&
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DCNVL',                                       & 
         LONG_NAME ='convective_source_of_cloud_liq',   &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DCNVI',                                       & 
         LONG_NAME ='convective_source_of_cloud_ice',        &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DLPDF',                                    & 
         LONG_NAME ='pdf_source_sink_of_cloud_liq',    &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DIPDF',                                    & 
         LONG_NAME ='pdf_source_sink_of_cloud_ice',    &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DLFIX',                                    & 
         LONG_NAME ='fix_source_sink_of_cloud_liq',          &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DIFIX',                                    & 
         LONG_NAME ='fix_source_sink_of_cloud_ice',          &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='AUT',                                      & 
         LONG_NAME ='autoconv_sink_of_cloud_liq',               &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='EVAPC',                                  & 
         LONG_NAME ='evaporation_of_cloud_liq',               &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SDM',                                    & 
         LONG_NAME ='sedimentation_sink_of_cloud_ice',        &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLICE_AN',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_anvil_snow',       &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLICE_LS',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_largescale_snow',  &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLWAT_AN',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_anvil_rain',     &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLWAT_LS',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_largescale_rain',&
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLRN_AN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_anvil_rain',              &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLRN_LS',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_largescale_rain',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLRN_CN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_convective_rain',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLSN_AN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_anvil_snow',              &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLSN_LS',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_largescale_snow',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLSN_CN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_convective_snow',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SUBLC',                                  & 
         LONG_NAME ='sublimation_of_cloud_ice',               &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='FRZ_TT',                                 & 
         LONG_NAME ='freezing_of_cloud_condensate',           &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='FRZ_PP',                                 & 
         LONG_NAME ='freezing_of_precip_condensate',          &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)




! Vertically integrated water substance conversions


    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='PDFLZ',                                        & 
         LONG_NAME ='statistical_source_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='PDFIZ',                                        & 
         LONG_NAME ='statistical_source_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CNVRNZ',                                        & 
         LONG_NAME ='convective_production_of_rain_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CNVLZ',                                        & 
         LONG_NAME ='convective_source_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CNVIZ',                                        & 
         LONG_NAME ='convective_source_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='EVPCZ',                                        & 
         LONG_NAME ='evaporation_loss_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SUBCZ',                                        & 
         LONG_NAME ='sumblimation_loss_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='EVPPZ',                                        & 
         LONG_NAME ='evaporation_loss_of_precip_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SUBPZ',                                        & 
         LONG_NAME ='sumblimation_loss_of_precip_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='AUTZ',                                        & 
         LONG_NAME ='autoconversion_loss_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SDMZ',                                        & 
         LONG_NAME ='sedimentation_loss_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='COLLLZ',                                        & 
         LONG_NAME ='accretion_loss_of_cloud_water_to_rain',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='COLLIZ',                                        & 
         LONG_NAME ='accretion_loss_of_cloud_water_to_snow',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='FRZCZ',                                        & 
         LONG_NAME ='net_freezing_of_cloud_condensate',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='FRZPZ',                                        & 
         LONG_NAME ='net_freezing_of_precip_condensate',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)








    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RCCODE',                                      & 
         LONG_NAME ='Convection_return_codes',                     &
         UNITS     ='codes',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TRIEDLV',                                     & 
         LONG_NAME ='Tested_for_convection_at_this_level',         &
         UNITS     ='0 or 1',                                      &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    ! MATMAT Exports for after-RAS inoutputs

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QVRAS',                                       & 
         LONG_NAME ='water_vapor_after_ras',                       &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='THRAS',                                       & 
         LONG_NAME ='potential_temperature_after_ras',             &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='URAS',                                        & 
         LONG_NAME ='eastward_wind_after_ras',                     &
         UNITS     ='m s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='VRAS',                                        & 
         LONG_NAME ='northward_wind_after_ras',                    &
         UNITS     ='m s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    ! MATMAT Exports for before-RAS inputs for RAStest

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='THOI',                                        & 
         LONG_NAME ='potential_temperature_before_ras',            &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QHOI',                                        & 
         LONG_NAME ='specific_humidity_before_ras',                &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QSSI',                                        & 
         LONG_NAME ='saturation_specific_humidity_before_ras',     &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='DQSI',                                        & 
         LONG_NAME ='deriv_sat_specific_humidity_wrt_t_before_ras',&
         UNITS     ='kg kg-1 K-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='PLEI',                                        & 
         LONG_NAME ='air_pressure_before_ras',                     &
         UNITS     ='Pa',                                          &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='TPERTI',                                      & 
         LONG_NAME ='temperature_perturbation_before_ras',         &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KCBLI',                                       & 
         LONG_NAME ='cloud_base_layer_before_ras',                 &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)


! !RECORD IMPORTS AND INTERNALS AT TOP OF MOIST:

     call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME ='QX0',                                          &
        LONG_NAME  ='specific_humidity',                          &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          

     call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME ='QLLSX0',                                       &
        LONG_NAME  ='initial_mass_fraction_of_large_scale_cloud_liquid_water', &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 
     call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME ='QLLSX1',                                       &
        LONG_NAME  ='final_mass_fraction_of_large_scale_cloud_liquid_water', &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 
     call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME ='QLCNX0',                                       &
        LONG_NAME  ='initial_mass_fraction_of_convective_cloud_liquid_water', &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 
     call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME ='QLCNX1',                                       &
        LONG_NAME  ='final_mass_fraction_of_convective_cloud_liquid_water', &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 
     call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME ='CLLSX0',                                       &
        LONG_NAME  ='large_scale_cloud_area_fraction',            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          

     call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME ='CLCNX0',                                       &
        LONG_NAME  ='convective_cloud_area_fraction',             &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          

     call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME ='QILSX0',                                       &
        LONG_NAME  ='initial_mass_fraction_of_large_scale_cloud_ice_water', &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 
     call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME ='QILSX1',                                       &
        LONG_NAME  ='final_mass_fraction_of_large_scale_cloud_ice_water', &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 

     call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME ='QICNX0',                                       &
        LONG_NAME  ='initial_mass_fraction_of_convective_cloud_ice_water', &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 
     call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME ='QICNX1',                                       &
        LONG_NAME  ='final_mass_fraction_of_convective_cloud_ice_water', &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          
 


     call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'KHX0',                                         &
        LONG_NAME  = 'scalar_diffusivity',                         &
        UNITS      = 'm+2 s-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                            &
        VLOCATION  = MAPL_VLocationEdge,                           &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'THX0',                                        &
        LONG_NAME  = 'potential_temperature',                     &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'UX0',                                         &
        LONG_NAME  = 'eastward_wind',                             &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'VX0',                                         &
        LONG_NAME  = 'northward_wind',                            &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TSX0',                                        &
        LONG_NAME  = 'surface temperature',                       &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'FRLANDX0',                                    &
        LONG_NAME  = 'areal_land_fraction',                       &
        UNITS      = '1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


!!! downdraft diags
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DDF_MFC',                                   &
        LONG_NAME  = 'Downdraft_mass_flux',                       &
        UNITS      = 'kg m-2 s-1',                                &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                          &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DDF_RH1',                                  &
        LONG_NAME  = 'Downdraft_in_cloud_RH_before',       &
        UNITS      = '1',                                       &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DDF_RH2',                                  &
        LONG_NAME  = 'Downdraft_in_cloud_RH_after',       &
        UNITS      = '1',                                       &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DDF_TC',                                    &
        LONG_NAME  = 'Temperature_excess_in_DDF',                 &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DDF_QVC',                                   &
        LONG_NAME  = 'Spec_hum_excess_in_DDF',                    &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DDF_BYNC',                                  &
        LONG_NAME  = 'Buoyancy_of_DDF',                           &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DDF_MUPH',                                  &
        LONG_NAME  = 'Downdraft_moistening_from_evap_subl',       &
        UNITS      = 'kg kg-1 s-1',                               &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DDF_DQDT',                                  &
        LONG_NAME  = 'Total_Downdraft_moistening',                      &
        UNITS      = 'kg kg-1 s-1',                               &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DDF_DTDT',                                  &
        LONG_NAME  = 'Total_Downdraft_heating',                         &
        UNITS      = 'K s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DDF_ZSCALE',                                &
        LONG_NAME  = 'vertical_scale_for_downdraft',              &
        UNITS      = 'm',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                            &
         SHORT_NAME = 'KEDISS',                                    &
         LONG_NAME  = 'kinetic_energy_diss_in_RAS',       &
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzVert,                          &
         VLOCATION  = MAPL_VLocationCenter,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                            &
         SHORT_NAME = 'KEMST',                                    &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_across_moist',       &
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                            &
         SHORT_NAME = 'KEMST2',                                    &
         LONG_NAME  = 'vertically_integrated_KE_dissipation_in_RAS',       &
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

! CAR 12/5/08
! Aerosol Scavenging diagnostics/export states
! ------------------------------
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DDUDT ',                                         & 
         LONG_NAME ='dust_tendency_due_to_conv_scav',                 &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DSSDT ',                                         &
         LONG_NAME ='sea_salt_tendency_due_to_conv_scav',             &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DOCDT ',                                         &
         LONG_NAME ='organic_carbon_tendency_due_to_conv_scav',       &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DBCDT ',                                         &
         LONG_NAME ='black_carbon_tendency_due_to_conv_scav',         &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DSUDT ',                                         &
         LONG_NAME ='sulfate_tendency_due_to_conv_scav',              &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DDUDTcarma ',                                    &
         LONG_NAME ='carma_dust_tendency_due_to_conv_scav',           &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DSSDTcarma ',                                    &
         LONG_NAME ='carma_seasalt_tendency_due_to_conv_scav',        &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
                                                       RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='LFR',                                            &
         LONG_NAME ='lightning_flash_rate',                           &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='A1X1',                                           &
         LONG_NAME ='LFR_Term_number_1',                              &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='A2X2',                                           &
         LONG_NAME ='LFR_Term_number_2',                              &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='A3X3',                                           &
         LONG_NAME ='LFR_Term_number_3',                              &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='A4X4',                                           &
         LONG_NAME ='LFR_Term_number_4',                              &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='A5X5',                                           &
         LONG_NAME ='LFR_Term_number_5',                              &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
                                                       RC=STATUS  )
    VERIFY_(STATUS)
!EOS


! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,name="DRIVER" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-PRE_RAS"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-RAS"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--RAS_RUN"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-POST_RAS"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-CLOUD"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--CLOUD_RUN"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--CLOUD_DATA"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="---CLOUD_DATA_DEVICE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="---CLOUD_DATA_CONST"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--CLOUD_ALLOC"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--CLOUD_DEALLOC"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-MISC"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--FLASH"    ,RC=STATUS)
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)
     
    RETURN_(ESMF_SUCCESS)
     
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!-!-!-!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: RUN -- Run method for the CONVECT component

! !INTERFACE:

  subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:
    
! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP


! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp), pointer   :: STATE
    type (ESMF_Config      )            :: CF
    type (ESMF_State       )            :: INTERNAL
    type (ESMF_Alarm)                   :: ALARM

! Local variables

    integer                             :: IM,JM,LM

    real, pointer, dimension(:,:)       :: LONS
    real, pointer, dimension(:,:)       :: LATS

!=============================================================================

! Begin... 

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'Run'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (STATE,"TOTAL")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get( STATE, IM=IM, JM=JM, LM=LM,   &
                               RUNALARM = ALARM,             &
                               CF       = CF,                &
                               LONS     = LONS,              &
                               LATS     = LATS,              &
                               INTERNAL_ESMF_STATE=INTERNAL, &
                                                   RC=STATUS )
    VERIFY_(STATUS)

! If its time, calculate convective tendencies
! --------------------------------------------

    if ( ESMF_AlarmIsRinging( ALARM, RC=status) ) then
       VERIFY_(STATUS)
       call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_TimerOn(STATE,"DRIVER")
         call MOIST_DRIVER(IM,JM,LM, RC=STATUS)
         VERIFY_(STATUS)
       call MAPL_TimerOff(STATE,"DRIVER")
    endif
    VERIFY_(STATUS)

    call MAPL_TimerOff(STATE,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MOIST_DRIVER(IM,JM,LM, RC)
      integer,           intent(IN ) :: IM, JM, LM
      integer, optional, intent(OUT) :: RC

      type RAS_Tracer_T
         real, pointer :: Q(:,:,:) => null()
      end type RAS_Tracer_T

!  Locals

      character(len=ESMF_MAXSTR)      :: IAm
      character(len=ESMF_MAXSTR)      :: NAME
      integer                         :: STATUS

      type (ESMF_TimeInterval)        :: TINT
      type (RAS_Tracer_T ), pointer   :: TRPtrs (:)
      type (RAS_Tracer_T ), pointer   :: TRIPtrs(:)

!  LOCAL COPY OF VARIABLES

      real, pointer, dimension(:,:,:) :: DQDTCN, DTHDTCN,DQCDTCN,DTDTFRIC
      real, pointer, dimension(:,:,:) :: CNV_DQLDT            , &
                                         CNV_MF0              , &
                                         CNV_MFD              , &
                                         CNV_MFC              , &
                                         CNV_UPDF             , &
                                         CNV_CVW              , &
                                         CNV_QC               , &
                                         RAD_CF               , &
                                         RAD_QL               , &
                                         RAD_QI               , &
                                         RAD_QR               , &
                                         RAD_QS               , &
                                         RAD_QV               , &
                                         CLDNCCN              , &                
                                         CLDREFFR             , &                
                                         CLDREFFS             , &                
                                         CLDREFFL             , &                
                                         CLDREFFI                             


      real, pointer, dimension(:,:,:) :: T, PLE, U, V, TH
      real, pointer, dimension(:,:,:) :: DQDT, UI, VI, TI, KH
      real, pointer, dimension(    :) :: PREF
      real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, BYNCY, QILS, QICN, QCTOT,QITOT,QLTOT
      real, pointer, dimension(:,:,:) :: QRTOT, QSTOT
      real, pointer, dimension(:,:  ) :: LS_PRCP,CN_PRCP,AN_PRCP,TT_PRCP,ER_PRCP,FILLNQV
      real, pointer, dimension(:,:  ) :: HOURNORAIN
      integer                         :: YEAR, MONTH, DAY, HR, SE, MN
      type (ESMF_Time)                :: CurrentTime
      real, pointer, dimension(:,:  ) :: LS_ARF, CN_ARF, AN_ARF
      real, pointer, dimension(:,:  ) :: SNR,PRECU,PRELS,TS,FRLAND,FROCEAN
      real, pointer, dimension(:,:  ) :: IWP,LWP,CWP,TPW,CAPE,ZPBLCN,INHB,ZLCL,ZLFC,ZCBL,CCWP,KPBLIN
      real, pointer, dimension(:,:  ) :: TVQ0,TVQ1,TVE0,TVE1,TVEX,DCPTE
      real, pointer, dimension(:,:,:,:) :: XHO, XHOI
      real, pointer, dimension(:,:  ) :: RASPBLQ, RASTIME, MXDIAM
 
      real, pointer, dimension(:,:  ) :: KEMST,KEMST2
      real, pointer, dimension(:,:,:) :: QSNOW, QRAIN

      logical, pointer                :: IS_FRIENDLY(:)
      logical, pointer                :: IS_WEIGHTED(:)

!     Tracer scavenging coefficients
!     FSCAV is the Fraction of tracer scavanged per km (=0 no scavenging, =1 full scavenging)
      real, pointer, dimension(:    ) :: FSCAV_, &  ! holding array for all tracers
                                         FSCAV      ! container for friendly to moist tracers

! Aerosol convective scavenging internal pointers (2D column-averages);  must deallocate!!!
! CAR 
      real, pointer, dimension(:,: )                           :: DDUDT, &
                           DSSDT, DOCDT, DBCDT, DSUDT, DDUDTcarma, DSSDTcarma
      character(len=ESMF_MAXSTR)                               :: QNAME,  CNAME, ENAME
      character(len=ESMF_MAXSTR), pointer, dimension(:)        :: QNAMES, CNAMES
      integer                                                  :: ind


      real, pointer, dimension(:,:)   :: PGENTOT , PREVTOT

 !Record vars at top pf moist
      real, pointer, dimension(:,:,:) :: Ux0, Vx0, THx0, KHx0
      real, pointer, dimension(:,:)   :: TSx0, FRLANDx0
      real, pointer, dimension(:,:,:) :: Qx0, QLLSx0, QLCNx0, CLLSx0, CLCNx0, QILSx0, QICNx0

      ! MATMAT Exports for pre-ras inputs for RAStest
      real, pointer, dimension(:,:,:) :: THOI,QHOI,QSSI,DQSI
      real, pointer, dimension(:,:,:) :: PLEI
      real, pointer, dimension(:,:  ) :: TPERTI,KCBLI

  !Extra outputs

      real, pointer, dimension(:,:,:) :: DQRL, DQRC, DQLDT, DQIDT, KEDISS
      real, pointer, dimension(:,:,:) :: RH1, RHX, RH2, XQLLS, XQLCN, XCLLS, XCLCN, XQILS, XQICN
      real, pointer, dimension(:,:,:) :: REV_CN, REV_AN, REV_LS
      real, pointer, dimension(:,:,:) :: RSU_CN, RSU_AN, RSU_LS,    ALPHT, ALPH1, ALPH2
      real, pointer, dimension(:,:,:) :: ENTLAM

      real, pointer, dimension(:,:,:) :: REVSU_CN, REVSU_LSAN
      real, pointer, dimension(:,:  ) :: CNV_FREQ, CNV_BASEP, CNV_TOPP


      real, pointer, dimension(:,:,:) :: ACLL_CN, ACLL_AN, ACLL_LS
      real, pointer, dimension(:,:,:) :: ACIL_CN, ACIL_AN, ACIL_LS
      real, pointer, dimension(:,:,:) :: ACR_TOT

      real, pointer, dimension(:,:,:) :: PFL_CN, PFL_AN, PFL_LS
      real, pointer, dimension(:,:,:) :: PFI_CN, PFI_AN, PFI_LS
      real, pointer, dimension(:,:,:) :: PFI_LSAN, PFL_LSAN

      real, pointer, dimension(:,:,:) :: DlPDF, DiPDF, DlFIX, DiFIX, FRZ_PP
      real, pointer, dimension(:,:,:) :: AUT, EVAPC, SDM, SUBLC, CFPDF, CFPDFX, RHCLR
      real, pointer, dimension(:,:,:) :: VFALLICE_AN, VFALLICE_LS
      real, pointer, dimension(:,:,:) :: VFALLWAT_AN,VFALLWAT_LS
      real, pointer, dimension(:,:,:) :: VFALLSN_AN,VFALLSN_LS,VFALLSN_CN
      real, pointer, dimension(:,:,:) :: VFALLRN_AN,VFALLRN_LS,VFALLRN_CN
      real, pointer, dimension(:,:,:) :: FRZ_TT, DCNVL, DCNVi,QSATi,QSATl,RCCODE,TRIEDLV,QVRAS

      real, pointer, dimension(:,:  ) :: LFR,A1X1,A2X2,A3X3,A4X4,A5X5

      !Trajectory for Moist TLM/ADJ
      real, pointer, dimension(:,:,:)  :: QILS_tlad, QLLS_tlad, QICN_tlad, QLCN_tlad, CFLS_tlad, CFCN_tlad
      real, pointer, dimension(:,:,:)  :: QLS_tlad, QCN_tlad
      real, pointer, dimension(:,:  )  :: KCBL_tlad, TSUR_tlad, KHu_tlad, KHl_tlad

      ! MATMAT Additional after-RAS exports
      real, pointer, dimension(:,:,:) :: THRAS,URAS,VRAS

    ! vapor-to-liquid  [ *ALHL ]
         real, pointer, dimension(:,:) :: PDFLZ,CNVLZ,CNVRNZ

    ! liquid-to-vapor  [ *ALHL ]
         real, pointer, dimension(:,:) :: EVPCZ, EVPPZ

    ! vapor-to-ice  [ *ALHS ]
         real, pointer, dimension(:,:) :: PDFIZ,CNVIZ 

    ! ice-to-vapor  [ *ALHS ]
         real, pointer, dimension(:,:) :: SUBCZ, SUBPZ

    ! liquid-to-ice  [ *ALHF ]
         real, pointer, dimension(:,:) :: FRZCZ, COLLIZ, FRZPZ

    ! liquid-to-liquid  [ * 0 ]
         real, pointer, dimension(:,:) :: AUTZ, COLLLZ

    ! ice-to-ice  [ * 0 ]
         real, pointer, dimension(:,:) :: SDMZ








  !Outputs for dpwndraft
      real, pointer, dimension(:,:,:) :: DDF_DQDT , DDF_DTDT , DDF_MFC , DDF_MUPH , DDF_RH1 , DDF_RH2 , DDF_BYNC
      real, pointer, dimension(:,:,:) :: DDF_QVC ,  DDF_TC
      real, pointer, dimension(:,:)   :: DDF_ZSCALE
      real, pointer, dimension(:)     :: DDF_DQDTz, DDF_DTDTz, DDF_MFCz, DDF_MUPHz, DDF_RH1z, DDF_RH2z, DDF_BYNCz
      real, pointer, dimension(:)     :: DDF_QVCz,  DDF_TCz
      real                            :: DDF_ZSCALEz

      type (ESMF_Field)               :: FIELD
      type (ESMF_FieldBundle)         :: TR
      type (ESMF_FieldBundle)         :: TRI

      real                            :: DT, DT_MOIST
      real(ESMF_KIND_R8)              :: DT_R8

      real                            :: PMIN_DET,AUTOC_CN_LAND,AUTOC_CN_OCN
      real                            :: PMIN_CBL
       
      integer, parameter              :: N_RAS_PARAMS     = 25
      integer, parameter              :: N_CLOUD_PARAMS   = 57

      real                            :: RASPARAMS   (N_RAS_PARAMS    )
      real                            :: CLOUDPARAMS (N_CLOUD_PARAMS  )
      real                            :: CBL_TPERT, CBL_QPERT, RASAL1, RASAL2, QUANX_0, QUANX_1
      real                            :: RHEXCESS,CBL_TPERT_MXOCN,CBL_TPERT_MXLND
      real                            :: HEARTBEAT

      integer                         :: L, K, I, J, KM, KK , Kii

      integer                         :: IDIM, IRUN, ICMIN, K0
      integer                         :: ITRCR,KSTRAP,CBL_METHOD,KCBLMIN,CLEANUP_RH

      real,    dimension(IM,JM,0:LM)  :: PKE
      real,    dimension(IM,JM,  LM)  :: DQS, QSS, PLO, ZLO, TEMP, PK, DM, DP, DQSDT
      real,    dimension(IM,JM,  LM)  :: KE0 , KEX, KE1, DKEX , DKE0, DKE1
      real,    dimension(IM,JM,  LM)  :: Q1, U1, V1, TH1, CNV_PRC3,fQi,CFPBL,CNV_HAIL

      real,    dimension(IM,JM,  LM)  :: Q1_dh, U1_dh, V1_dh, TH1_dh
      real,    dimension(IM,JM,  LM)  :: QILS_dh, QLLS_dh, QICN_dh, QLCN_dh, CFLS_dh, CFCN_dh

      real,    dimension(IM,JM,0:LM)  :: CNV_PLE,ZLE
      real,    dimension(      0:LM)  :: SIGE
      real,    dimension(IM,JM,  LM)  :: GZLO, HHO,HSO
      real,    dimension(IM,JM,0:LM)  :: GZLE
      real,    dimension(IM,JM)       :: RASPRCP
      real,    dimension(IM,JM)       :: CO_AUTO
      integer, dimension(IM,JM,2)     :: SEEDRAS
      integer, dimension(IM,JM)       :: KLCL, KLFC, KPBL

      real,    dimension(IM,JM)       :: LS_ARFX,CN_ARFX,AN_ARFX,QSSFC,IKEX,IKEX2
      real,    dimension(IM,JM)       :: LS_SNR, CN_SNR, AN_SNR,  ZCBLx, FILLQ, MXDIAMx
      real,    dimension(IM,JM)       :: LS_PRC2,CN_PRC2,AN_PRC2,ER_PRC2, TPREC, TVQX
      real,    dimension(IM,JM)       :: TPERT, QPERT, TSFCAIR, DTS, QUANX, QUANXTP
      integer, dimension(IM,JM)       :: IRAS, JRAS, KCBL
      real,    dimension(IM,JM,LM)    :: WGT0, WGT1
      real,    dimension(IM,JM,LM)    :: TRDLX
      integer, dimension(IM,JM,LM)    :: irccode

      type (T_CLOUD_CTL)              :: CLOUD_CTL
      type (T_DDF_CTL)                :: DDF_CTL

      !!real,    dimension(IM,JM,  LM)  :: QILS, QICN ! Soon to be moved into internal state

      logical ALLOC_CNV_DQLDT 
      logical ALLOC_CNV_MF0   
      logical ALLOC_CNV_MFD   
      logical ALLOC_CNV_MFC
      logical ALLOC_CNV_TOPP   
      logical ALLOC_CNV_UPDF  
      logical ALLOC_CNV_CVW  
      logical ALLOC_CNV_QC 
      logical ALLOC_RAD_CF    
      logical ALLOC_RAD_QV    
      logical ALLOC_RAD_QL    
      logical ALLOC_RAD_QI    
      logical ALLOC_RAD_QR    
      logical ALLOC_RAD_QS    
      logical ALLOC_CLDREFFL  
      logical ALLOC_CLDREFFI  
      logical ALLOC_CLDREFFR  
      logical ALLOC_CLDREFFS  
      logical ALLOC_CLDNCCN   
      logical ALLOC_BYNCY
      logical ALLOC_CAPE 
      logical ALLOC_INHB 

      logical ALLOC_DQRC

      logical ALLOC_ENTLAM

      character(len=ESMF_MAXSTR) :: GRIDNAME
      character(len=4)           :: imchar
      character(len=2)           :: dateline
      integer                    :: imsize,nn

      integer                    :: levs925
      real, dimension(IM,JM   )  :: tempor2d

! Manage diagnostic outputs for re-evaporation
!---------------------------------------------------
      real, dimension(IM,JM,LM) :: RSU_CN_X
      real, dimension(IM,JM,LM) :: RSU_AN_X
      real, dimension(IM,JM,LM) :: RSU_LS_X

      real, dimension(IM,JM,LM) :: REV_CN_X
      real, dimension(IM,JM,LM) :: REV_AN_X
      real, dimension(IM,JM,LM) :: REV_LS_X

! Manage diagnostic outputs for accretion 
!---------------------------------------------------
      real, dimension(IM,JM,LM) :: ACIL_CN_X
      real, dimension(IM,JM,LM) :: ACIL_AN_X
      real, dimension(IM,JM,LM) :: ACIL_LS_X

      real, dimension(IM,JM,LM) :: ACLL_CN_X
      real, dimension(IM,JM,LM) :: ACLL_AN_X
      real, dimension(IM,JM,LM) :: ACLL_LS_X

! Manage diagnostic outputs for 3D precip fluxes
!---------------------------------------------------
      real, dimension(IM,JM,0:LM) :: PFI_CN_X
      real, dimension(IM,JM,0:LM) :: PFI_AN_X
      real, dimension(IM,JM,0:LM) :: PFI_LS_X

      real, dimension(IM,JM,0:LM) :: PFL_CN_X
      real, dimension(IM,JM,0:LM) :: PFL_AN_X
      real, dimension(IM,JM,0:LM) :: PFL_LS_X

      real, dimension(IM,JM,LM) :: DLPDF_X
      real, dimension(IM,JM,LM) :: DIPDF_X
      real, dimension(IM,JM,LM) :: DLFIX_X
      real, dimension(IM,JM,LM) :: DIFIX_X

      real, dimension(IM,JM,LM) :: RHX_X
      real, dimension(IM,JM,LM) :: AUT_X
      real, dimension(IM,JM,LM) :: EVAPC_X
      real, dimension(IM,JM,LM) :: SDM_X
      real, dimension(IM,JM,LM) :: SUBLC_X
      real, dimension(IM,JM,LM) :: FRZ_TT_X
      real, dimension(IM,JM,LM) :: FRZ_PP_X
      real, dimension(IM,JM,LM) :: DCNVL_X
      real, dimension(IM,JM,LM) :: DCNVI_X

      real, dimension(IM,JM,LM) :: ALPHT_X
      real, dimension(IM,JM,LM) :: ALPH1_X
      real, dimension(IM,JM,LM) :: ALPH2_X

      real, dimension(IM,JM,LM) :: CFPDF_X
      real, dimension(IM,JM,LM) :: RHCLR_X
      real, dimension(IM,JM,LM) :: DQRL_X

      real, dimension(IM,JM,LM) :: VFALLICE_AN_X
      real, dimension(IM,JM,LM) :: VFALLICE_LS_X
      real, dimension(IM,JM,LM) :: VFALLWAT_AN_X
      real, dimension(IM,JM,LM) :: VFALLWAT_LS_X

      real, dimension(IM,JM,LM) :: VFALLRN_AN_X
      real, dimension(IM,JM,LM) :: VFALLRN_LS_X
      real, dimension(IM,JM,LM) :: VFALLRN_CN_X
      real, dimension(IM,JM,LM) :: VFALLSN_AN_X
      real, dimension(IM,JM,LM) :: VFALLSN_LS_X
      real, dimension(IM,JM,LM) :: VFALLSN_CN_X

! MATMAT Additional variables for GPUization
      real, dimension(IM,JM,LM)   :: DQST3,QST3,TEMP_0
      real, dimension(IM,JM,LM)   :: DZET, MASS, QDDF3
      real, dimension(IM,JM)      :: VMIP
      real, dimension(IM,JM,LM+1) :: ZET

! MATMAT CUDA Variables
#ifdef _CUDA
      type(dim3) :: Grid, Block
      integer :: blocksize

      logical :: COPY_RSU_CN
      logical :: COPY_RSU_AN
      logical :: COPY_RSU_LS

      logical :: COPY_REV_CN
      logical :: COPY_REV_AN
      logical :: COPY_REV_LS

      logical :: COPY_ACIL_CN
      logical :: COPY_ACIL_AN
      logical :: COPY_ACIL_LS

      logical :: COPY_ACLL_CN
      logical :: COPY_ACLL_AN
      logical :: COPY_ACLL_LS

      logical :: COPY_DLPDF
      logical :: COPY_DIPDF
      logical :: COPY_DLFIX
      logical :: COPY_DIFIX

      logical :: COPY_PFI_CN
      logical :: COPY_PFI_AN
      logical :: COPY_PFI_LS

      logical :: COPY_PFL_CN
      logical :: COPY_PFL_AN
      logical :: COPY_PFL_LS

      logical :: COPY_RHX
      logical :: COPY_AUT     
      logical :: COPY_EVAPC   
      logical :: COPY_SDM    
      logical :: COPY_SUBLC    
      logical :: COPY_FRZ_TT 
      logical :: COPY_FRZ_PP   
      logical :: COPY_DCNVL    
      logical :: COPY_DCNVI   

      logical :: COPY_ALPHT
      logical :: COPY_ALPH1
      logical :: COPY_ALPH2

      logical :: COPY_CFPDF
      logical :: COPY_RHCLR
      logical :: COPY_DQRL    

      logical :: COPY_VFALLICE_AN
      logical :: COPY_VFALLICE_LS
      logical :: COPY_VFALLWAT_AN
      logical :: COPY_VFALLWAT_LS

      logical :: COPY_VFALLRN_AN
      logical :: COPY_VFALLRN_LS
      logical :: COPY_VFALLRN_CN
      logical :: COPY_VFALLSN_AN
      logical :: COPY_VFALLSN_LS
      logical :: COPY_VFALLSN_CN

#endif

!  Begin...
!----------
      Iam = trim(COMP_NAME) // 'Convect_Driver'

      call MAPL_TimerOn(STATE,"-MISC")

! Get parameters from configuration
!----------------------------------
      call ESMF_ConfigGetAttribute (CF, HEARTBEAT, Label="RUN_DT:", RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_GetResource(STATE,CLEANUP_RH,   'CLEANUP_RH:',    DEFAULT= 1     ,RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS( 1),'CUFRICFAC:',     DEFAULT= 1.000, RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS( 2),'SHR_LAMBDA_FAC:',DEFAULT= 0.05,  RC=STATUS)
      call MAPL_GetResource(STATE,AUTOC_CN_OCN, 'AUTOC_CN:',      DEFAULT= 2.5e-3,RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS( 4),'QC_CRIT_CN:',    DEFAULT= 8.0e-4,RC=STATUS)
      call MAPL_GetResource(STATE,RASAL1       ,'RASAL1:',        DEFAULT=  1800.0,   RC=STATUS)
      call MAPL_GetResource(STATE,RASAL2       ,'RASAL2:',        DEFAULT= 43200.0,   RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS( 7),'RASNCL:',        DEFAULT= -300., RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS( 8),'LAMBDA_FAC:',    DEFAULT= 4.0 ,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS( 9),'LAMBMX_FAC:',    DEFAULT= 0.0 ,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(10),'MIN_DIAMETER:',  DEFAULT= 200.,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(11),'CUFRICLAMBDA:',  DEFAULT= 7.5e-4,RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(12),'RDTLEXPON:',     DEFAULT= 1.0 ,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(13),'STRAPPING:',     DEFAULT=-1.0 ,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(14),'SDQV2:',         DEFAULT= 1.3 ,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(15),'SDQV3:',         DEFAULT= 1.3 ,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(16),'SDQVT1:',        DEFAULT= 263.,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(17),'ACRITFAC:',      DEFAULT= 0.5 ,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(18),'HMINTRIGGER:',   DEFAULT= 1.0 ,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(19),'LLDISAGGXP:',    DEFAULT= 0.0 ,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(20),'PBLFRAC:',       DEFAULT= 0.1 ,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(21),'RASAUTORAMPB:',  DEFAULT= 0.8 ,  RC=STATUS)
      call MAPL_GetResource(STATE,PMIN_DET     ,'PMIN_DET',       DEFAULT= 3000.0,  RC=STATUS)
      call MAPL_GetResource(STATE,PMIN_CBL     ,'PMIN_CBL',       DEFAULT= 50000.0, RC=STATUS)
     
      call MAPL_GetResource(STATE,AUTOC_CN_LAND,'AUTOC_CN_LAND:', DEFAULT= AUTOC_CN_OCN, RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(22),'AUTOC_CN_ZDEP:', DEFAULT= 1.0   ,RC=STATUS)


      call MAPL_GetResource(STATE,GRIDNAME,'AGCM_GRIDNAME:', RC=STATUS)
      VERIFY_(STATUS)
      GRIDNAME =  AdjustL(GRIDNAME)
            nn = len_trim(GRIDNAME)
      dateline = GRIDNAME(nn-1:nn)
        imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
        read(imchar,*) imsize
      if(dateline.eq.'CF') imsize = imsize*4
      if( imsize.le.200       ) call MAPL_GetResource(STATE,RASPARAMS(23),'MAXDALLOWED:', DEFAULT= 4000.0 ,RC=STATUS)
      if( imsize.gt.200 .and. &
          imsize.le.400       ) call MAPL_GetResource(STATE,RASPARAMS(23),'MAXDALLOWED:', DEFAULT= 2000.0 ,RC=STATUS)
      if( imsize.gt.400 .and. &
          imsize.le.800       ) call MAPL_GetResource(STATE,RASPARAMS(23),'MAXDALLOWED:', DEFAULT=  700.0 ,RC=STATUS)
      if( imsize.gt.800 .and. &
          imsize.le.1600      ) call MAPL_GetResource(STATE,RASPARAMS(23),'MAXDALLOWED:', DEFAULT=  450.0 ,RC=STATUS)
      if( imsize.gt.1600      ) call MAPL_GetResource(STATE,RASPARAMS(23),'MAXDALLOWED:', DEFAULT=  450.0 ,RC=STATUS)

      call MAPL_GetResource(STATE,RASPARAMS(24),'RAS_RHMIN:',  DEFAULT= 0.5  ,  RC=STATUS)
      call MAPL_GetResource(STATE,RASPARAMS(25),'RAS_RHFULL:', DEFAULT= 0.65 ,  RC=STATUS)

      call MAPL_GetResource(STATE,CBL_METHOD,   'CBL_METHOD:',    DEFAULT= 6     , RC=STATUS)
      call MAPL_GetResource(STATE,CBL_TPERT,    'CBL_TPERT:',     DEFAULT= 1.0   , RC=STATUS)
      call MAPL_GetResource(STATE,CBL_QPERT,    'CBL_QPERT:',     DEFAULT= 0.0   , RC=STATUS)

      call MAPL_GetResource(STATE,CBL_TPERT_MXOCN, 'CBL_TPERT_MXOCN:',     DEFAULT= 2.0   , RC=STATUS)
      call MAPL_GetResource(STATE,CBL_TPERT_MXLND, 'CBL_TPERT_MXLND:',     DEFAULT= 4.0   , RC=STATUS)
      
      KSTRAP = INT( RASPARAMS(13) )

      RASPARAMS( 5) = RASAL1
      RASPARAMS( 6) = RASAL2

! Get the time step from the alarm
!---------------------------------

      call ESMF_AlarmGet( ALARM, RingInterval=TINT, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_TimeIntervalGet(TINT, S_R8=DT_R8, RC=STATUS)
      VERIFY_(STATUS)

      DT_MOIST = DT_R8

! Pointers to internals
!----------------------

      call MAPL_GetPointer(INTERNAL, Q,     'Q'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QLLS,  'QLLS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QLCN,  'QLCN'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, CLCN,  'CLCN'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, CLLS,  'CLLS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QILS,  'QILS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QICN,  'QICN'    , RC=STATUS); VERIFY_(STATUS)

! Pointers to imports
!--------------------
   
      call MAPL_GetPointer(IMPORT, PLE,     'PLE'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, PREF,    'PREF'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, KH,      'KH'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, TH,      'TH'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, U,       'U'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, V,       'V'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, FRLAND,  'FRLAND'  , RC=STATUS); VERIFY_(STATUS)
                                    ! frland =0.0
      call MAPL_GetPointer(IMPORT, FROCEAN, 'FROCEAN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, TS,      'TS'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, KPBLIN,  'KPBL'    , RC=STATUS); VERIFY_(STATUS)

      call ESMF_StateGet (IMPORT, 'MTR' ,   TR,        RC=STATUS); VERIFY_(STATUS)

! Pointers to exports
!--------------------

      call MAPL_GetPointer(EXPORT, TI,       'DTHDT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, UI,       'DUDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VI,       'DVDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQDT,     'DQDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SNR,      'SNO'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PRELS,    'PLS'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PRECU,    'PCU'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RH1   ,   'RH1'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RH2   ,   'RH2'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RHX   ,   'RHX'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, REV_CN,   'REV_CN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, REV_AN,   'REV_AN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, REV_LS,   'REV_LS'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RSU_CN,   'RSU_CN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RSU_AN,   'RSU_AN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RSU_LS,   'RSU_LS'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQRL,     'DQRL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQRC,     'DQRC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQDTCN,   'DQDTCN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQLDT,    'DQLDT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQIDT,    'DQIDT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQCDTCN,  'DQCDTCN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DTHDTCN,  'DTHDTCN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CCWP,     'CCWP'    , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, CWP,      'CWP'     , RC=STATUS); VERIFY_(STATUS)  
      call MAPL_GetPointer(EXPORT, IWP,      'IWP'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, LWP,      'LWP'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TPW,      'TPW'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, BYNCY,    'BYNCY'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CAPE,     'CAPE'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, INHB,     'INHB'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TVQ0,     'TVQ0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TVQ1,     'TVQ1'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DCPTE,    'DCPTE'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TVE0,     'TVE0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TVE1,     'TVE1'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TVEX,     'TVEX'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ZPBLCN,   'ZPBLCN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ZLCL,     'ZLCL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ZLFC,     'ZLFC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ZCBL,     'ZCBL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TT_PRCP,  'TPREC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, HOURNORAIN,  'HOURNORAIN'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, LS_PRCP,  'LS_PRCP' , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, CN_PRCP,  'CN_PRCP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, AN_PRCP,  'AN_PRCP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ER_PRCP,  'ER_PRCP' , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, FILLNQV,  'FILLNQV' , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, XQLLS,    'QLLSX1'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XQLCN,    'QLCNX1'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XQILS,    'QILSX1'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XQICN,    'QICNX1'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XCLCN,    'CLCN'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XCLLS,    'CLLS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QCTOT,    'QCTOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QITOT,    'QITOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QRTOT,    'QRTOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QSTOT,    'QSTOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLTOT,    'QLTOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, LS_ARF,   'LS_ARF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CN_ARF,   'CN_ARF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, AN_ARF,   'AN_ARF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_DQLDT,'CNV_DQLDT',RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_MF0,  'CNV_MF0' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_MFD,  'CNV_MFD' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_MFC,  'CNV_MFC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_FREQ, 'CNV_FREQ', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_BASEP,'CNV_BASEP', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_TOPP, 'CNV_TOPP', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_UPDF, 'CNV_UPDF', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_CVW,  'CNV_CVW' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_QC,   'CNV_QC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_CF,   'FCLD'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_QV,   'QV'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_QL,   'QL'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_QI,   'QI'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_QR,   'QR'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_QS,   'QS'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDREFFL, 'RL'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDREFFI, 'RI'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDREFFR, 'RR'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDREFFS, 'RS'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDNCCN,  'CLDNCCN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,DTDTFRIC, 'DTDTFRIC' , RC=STATUS); VERIFY_(STATUS)

      !Trajectory for the moist TLM/ADJ
      call MAPL_GetPointer(EXPORT, QILS_tlad, 'QILS_tlad', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLLS_tlad, 'QLLS_tlad', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QICN_tlad, 'QICN_tlad', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLCN_tlad, 'QLCN_tlad', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFLS_tlad, 'CFLS_tlad', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFCN_tlad, 'CFCN_tlad', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KCBL_tlad, 'KCBL_tlad', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TSUR_tlad, 'TSUR_tlad', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLS_tlad , 'QLS_tlad' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QCN_tlad , 'QCN_tlad' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KHu_tlad , 'KHu_tlad' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KHl_tlad , 'KHl_tlad' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, REVSU_LSAN, 'REVSU_LSAN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, REVSU_CN,   'REVSU_CN'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ACLL_CN,   'ACRLL_CN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ACLL_AN,   'ACRLL_AN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ACLL_LS,   'ACRLL_LS'  , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ACIL_CN,   'ACRIL_CN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ACIL_AN,   'ACRIL_AN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ACIL_LS,   'ACRIL_LS'  , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ACR_TOT,   'ACR_TOT'  , RC=STATUS); VERIFY_(STATUS)
      

      call MAPL_GetPointer(EXPORT, PFL_CN,   'PFL_CN'  , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFL_AN,   'PFL_AN'  , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFL_LS,   'PFL_LS'  , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFL_LSAN, 'PFL_LSAN' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, PFI_CN,   'PFI_CN'  , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFI_AN,   'PFI_AN'  , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFI_LS,   'PFI_LS'  , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFI_LSAN, 'PFI_LSAN'  , RC=STATUS); VERIFY_(STATUS)


      call MAPL_GetPointer(EXPORT, DlPDF,  'DLPDF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DiPDF,  'DIPDF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DlFIX,  'DLFIX'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DiFIX,  'DIFIX'  , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, AUT,    'AUT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, EVAPC,  'EVAPC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SDM,    'SDM'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLICE_AN, 'VFALLICE_AN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLICE_LS, 'VFALLICE_LS' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLWAT_AN, 'VFALLWAT_AN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLWAT_LS, 'VFALLWAT_LS' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLRN_AN, 'VFALLRN_AN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLRN_LS, 'VFALLRN_LS' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLRN_CN, 'VFALLRN_CN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLSN_AN, 'VFALLSN_AN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLSN_LS, 'VFALLSN_LS' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLSN_CN, 'VFALLSN_CN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SUBLC,  'SUBLC'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, FRZ_TT,    'FRZ_TT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, FRZ_PP,    'FRZ_PP'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DCNVL,     'DCNVL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DCNVI,     'DCNVI'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, QSATl,     'QSATL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QSATi,     'QSATI'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ALPHT,     'ALPHT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ALPH1,     'ALPH1'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ALPH2,     'ALPH2'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFPDF,     'CFPDF'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFPDFX,    'CFPDFX'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RHCLR,     'RHCLR'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, RCCODE,    'RCCODE'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TRIEDLV,   'TRIEDLV'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, QVRAS,     'QVRAS'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THRAS,     'THRAS'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, URAS,      'URAS '   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VRAS,      'VRAS '   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, RASPBLQ,   'RASPBLQ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RASTIME,   'RASTIME'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MXDIAM,    'MXDIAM'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, PREVTOT,   'PREVTOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PGENTOT,   'PGENTOT'   , RC=STATUS); VERIFY_(STATUS)


      call MAPL_GetPointer(EXPORT, FRZCZ,     'FRZCZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, FRZPZ,     'FRZPZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, COLLIZ,    'COLLIZ'  , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, CNVLZ,     'CNVLZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PDFLZ,     'PDFLZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNVRNZ,    'CNVRNZ'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, EVPCZ,     'EVPCZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, EVPPZ,     'EVPPZ'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, CNVIZ,     'CNVIZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PDFIZ,     'PDFIZ'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, SUBCZ,     'SUBCZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SUBPZ,     'SUBPZ'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, AUTZ,      'AUTZ'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, COLLLZ,    'COLLLZ'  , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, SDMZ,      'SDMZ'    , RC=STATUS); VERIFY_(STATUS)

      call ESMF_StateGet(EXPORT, 'MTRI',    TRI,        RC=STATUS); VERIFY_(STATUS)


      call MAPL_GetPointer(EXPORT, DDF_RH1,    'DDF_RH1'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_RH2,    'DDF_RH2'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_MUPH,   'DDF_MUPH'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_BYNC,   'DDF_BYNC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_TC,     'DDF_TC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_QVC,    'DDF_QVC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_DQDT,   'DDF_DQDT'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_DTDT,   'DDF_DTDT'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_MFC,    'DDF_MFC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_ZSCALE, 'DDF_ZSCALE', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KEMST,      'KEMST'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KEMST2,     'KEMST2'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KEDISS,     'KEDISS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ENTLAM,     'ENTLAM'    , RC=STATUS); VERIFY_(STATUS)
! Aerosol Scavenging
! CAR 12/5/08
      call MAPL_GetPointer(EXPORT, DDUDT,     'DDUDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DSSDT,     'DSSDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DBCDT,     'DBCDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DOCDT,     'DOCDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DSUDT,     'DSUDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDUDTcarma,  'DDUDTcarma'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DSSDTcarma,  'DSSDTcarma'    , RC=STATUS); VERIFY_(STATUS)




! Count the fields in TR...
!--------------------------

      call ESMF_FieldBundleGet(TR,FieldCount=KM , RC=STATUS)
      VERIFY_(STATUS)

! ...and make sure the other bundles are the same.
!-------------------------------------------------

      call ESMF_FieldBundleGet(TRI,FieldCount=K , RC=STATUS)
      VERIFY_(STATUS)
      ASSERT_(KM==K)

! Allocate tracer stuff
!----------------------

      allocate(IS_FRIENDLY(KM),stat=STATUS)
      VERIFY_(STATUS)
      allocate(IS_WEIGHTED(KM),stat=STATUS)
      VERIFY_(STATUS)
      allocate(TRPtrs     (KM),stat=STATUS)
      VERIFY_(STATUS)
      allocate(TRIPtrs    (KM),stat=STATUS)
      VERIFY_(STATUS)
      allocate(FSCAV_     (KM),stat=STATUS)
      VERIFY_(STATUS)
      allocate(QNAMES(KM), CNAMES(KM), stat=STATUS)
      VERIFY_(STATUS)
 
      QNAMES = " "
      CNAMES = " "
      QNAME = " "
      CNAME = " "

!CAR get the names of things
     call ESMF_FieldBundleGet(TR, fieldNameList=QNAMES, rc=STATUS )
     VERIFY_(STATUS) 

! Loop over all quantities to be diffused.
!----------------------------------------

      ITRCR = 0

      do K=1,KM

! Get the Kth Field from tracer bundle
!-------------------------------------

         call ESMF_FieldBundleGet(TR, K, FIELD, RC=STATUS)
         VERIFY_(STATUS)

         call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
         VERIFY_(STATUS)

! Get item's friendly status (default is not friendly)
!-----------------------------------------------------

         call ESMF_AttributeGet  (FIELD, NAME="FriendlyToMOIST",VALUE=IS_FRIENDLY(K), RC=STATUS)
         if(STATUS /= ESMF_SUCCESS) IS_FRIENDLY(K) = .false.

! Get item's weighting (default is unweighted tendencies)
!--------------------------------------------------------

         call ESMF_AttributeGet  (FIELD, NAME="WeightedTendency",VALUE=IS_WEIGHTED(K), RC=STATUS)
         if(STATUS /= ESMF_SUCCESS) IS_WEIGHTED(K) = .false.

! Get item's scavenging fraction
!-------------------------------

         call ESMF_AttributeGet  (FIELD, NAME="ScavengingFractionPerKm",VALUE=FSCAV_(K), RC=STATUS)
         if(STATUS /= ESMF_SUCCESS) FSCAV_(K) = 0.0 ! no scavenging

! Check aerosol names
! CAR 12/5/08

       !PRINT *, "*********NAME CHECKING:*******"
       !PRINT *, trim(QNAMES(K))
       
       ! Remove qualifier from variable name:  GOCART::du001->du001
       ! CAR 12/5/08
       !-------------------------------
       
       QNAME = QNAMES(K)
       ind= index(QNAME, '::')
       if (ind> 0) then
          CNAMES(K) = trim(QNAME(1:ind-1))  ! Component name (e.g., GOCART, CARMA)
          QNAMES(K) = trim(QNAME(ind+2:))
       end if
       !PRINT *, "******CROPPED NAME CHECKING*******"
       !PRINT *, trim(QNAMES(K)), FSCAV_(K)


! Get pointer to the quantity, its tendency, its surface value,
!   the surface flux, and the sensitivity of the surface flux.
! -------------------------------------------------------------

         if (IS_FRIENDLY(K)) then
            ITRCR = ITRCR + 1
            call ESMFL_BundleGetPointerToData(TR    ,      NAME,         TRPtrs (K)%Q, RC=STATUS)
            VERIFY_(STATUS)
            call ESMFL_BundleGetPointerToData(TRI   , trim(NAME)//'IM' , TRIPtrs(K)%Q, RC=STATUS)
            VERIFY_(STATUS)
         end if

      end do
! CAR END LOOP OVER K

! Allocate space for tracers
!---------------------------

      allocate(XHO (IM,JM,LM,ITRCR),stat=STATUS)
      VERIFY_(STATUS)
      allocate(XHOI(IM,JM,LM,ITRCR),stat=STATUS)
      VERIFY_(STATUS)
!     FSCAV changes dimensions of FSCAV_
      allocate(FSCAV(ITRCR),stat=STATUS)
      VERIFY_(STATUS)


      ALLOC_CNV_DQLDT = .not.associated(CNV_DQLDT )
      ALLOC_CNV_MF0   = .not.associated(CNV_MF0   )
      ALLOC_CNV_MFD   = .not.associated(CNV_MFD   )
      ALLOC_CNV_MFC   = .not.associated(CNV_MFC   )
      ALLOC_CNV_TOPP  = .not.associated(CNV_TOPP  )
      ALLOC_CNV_UPDF  = .not.associated(CNV_UPDF  )
      ALLOC_CNV_CVW   = .not.associated(CNV_CVW   )
      ALLOC_CNV_QC    = .not.associated(CNV_QC    )
      ALLOC_RAD_CF    = .not.associated(RAD_CF    )
      ALLOC_RAD_QV    = .not.associated(RAD_QV    )
      ALLOC_RAD_QL    = .not.associated(RAD_QL    )
      ALLOC_RAD_QI    = .not.associated(RAD_QI    )
      ALLOC_RAD_QR    = .not.associated(RAD_QR    )
      ALLOC_RAD_QS    = .not.associated(RAD_QS    )
      ALLOC_CLDREFFL  = .not.associated(CLDREFFL  )
      ALLOC_CLDREFFI  = .not.associated(CLDREFFI  )
      ALLOC_CLDREFFR  = .not.associated(CLDREFFR  )
      ALLOC_CLDREFFS  = .not.associated(CLDREFFS  )
      ALLOC_CLDNCCN   = .not.associated(CLDNCCN   )
      ALLOC_ENTLAM    = .not.associated(ENTLAM    )

      if(ALLOC_CNV_DQLDT) allocate(CNV_DQLDT(IM,JM,LM))
      if(ALLOC_CNV_MF0  ) allocate(CNV_MF0  (IM,JM,LM))
      if(ALLOC_CNV_MFD  ) allocate(CNV_MFD  (IM,JM,LM))
      if(ALLOC_CNV_MFC  ) allocate(CNV_MFC  (IM,JM,0:LM))
      if(ALLOC_CNV_TOPP ) allocate(CNV_TOPP (IM,JM))
      if(ALLOC_CNV_UPDF ) allocate(CNV_UPDF (IM,JM,LM))
      if(ALLOC_CNV_CVW  ) allocate(CNV_CVW  (IM,JM,LM))
      if(ALLOC_CNV_QC   ) allocate(CNV_QC   (IM,JM,LM))
      if(ALLOC_RAD_CF   ) allocate(RAD_CF   (IM,JM,LM))
      if(ALLOC_RAD_QV   ) allocate(RAD_QV   (IM,JM,LM))
      if(ALLOC_RAD_QL   ) allocate(RAD_QL   (IM,JM,LM))
      if(ALLOC_RAD_QI   ) allocate(RAD_QI   (IM,JM,LM))
      if(ALLOC_RAD_QR   ) allocate(RAD_QR   (IM,JM,LM))
      if(ALLOC_RAD_QS   ) allocate(RAD_QS   (IM,JM,LM))
      if(ALLOC_CLDREFFL ) allocate(CLDREFFL (IM,JM,LM))
      if(ALLOC_CLDREFFI ) allocate(CLDREFFI (IM,JM,LM))
      if(ALLOC_CLDREFFR ) allocate(CLDREFFR (IM,JM,LM))
      if(ALLOC_CLDREFFS ) allocate(CLDREFFS (IM,JM,LM))
      if(ALLOC_CLDNCCN  ) allocate(CLDNCCN  (IM,JM,LM))
      if(ALLOC_ENTLAM   ) allocate(ENTLAM   (IM,JM,LM))

      ALLOC_DQRC     = .not. associated( DQRC  )

      if(ALLOC_DQRC  ) allocate ( DQRC     (IM,JM,LM))

!  Allocate local space for some diagnostics that have to be in arg list to progno_cloud
      allocate(QSNOW(IM,JM,LM))
      allocate(QRAIN(IM,JM,LM))
      QSNOW = 0.
      QRAIN = 0.

! Recording of import/internal vars into export if desired
!---------------------------------------------------------

      call MAPL_GetPointer(EXPORT, Qx0,     'QX0'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLLSx0,  'QLLSX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLCNx0,  'QLCNX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLCNx0,  'CLCNX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLLSx0,  'CLLSX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QILSx0,  'QILSX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QICNx0,  'QICNX0'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, KHx0,    'KHX0'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THx0,    'THX0'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, Ux0,     'UX0'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, Vx0,     'VX0'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TSx0,    'TSX0'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, FRLANDx0,'FRLANDX0' , RC=STATUS); VERIFY_(STATUS)

      if(associated(Qx0       )) Qx0        = Q
      if(associated(QLLSx0    )) QLLSx0     = QLLS
      if(associated(QILSx0    )) QILSx0     = QILS
      if(associated(QLCNx0    )) QLCNx0     = QLCN
      if(associated(QICNx0    )) QICNx0     = QICN
      if(associated(CLLSx0    )) CLLSx0     = CLLS
      if(associated(CLCNx0    )) CLCNx0     = CLCN

      if(associated(KHx0      )) KHx0       = KH
      if(associated(THx0      )) THx0       = TH
      if(associated(Ux0       )) Ux0        = U
      if(associated(Vx0       )) Vx0        = V
      if(associated(TSx0      )) TSx0       = TS
      if(associated(FRLANDx0  )) FRLANDx0   = FRLAND

!  Copy incoming state vars to local arrays that will be adjusted
!  by physics.  Untouched state vars will later be used for 
!  post facto tendency calculations.
!----------------------------------------------------------------

      TH1      = TH
      Q1       = Q
      U1       = U
      V1       = V

      CNV_HAIL = 0.0
      CNV_PLE  = PLE*.01
      PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )

!      PKE      = (PLE/MAPL_P00)**(MAPL_KAPPA)
      PKE      = (CNV_PLE/1000.)**(MAPL_RGAS/MAPL_CP)
      DP       = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
!      PK       = (PLE(:,:,1:LM)*PKE(:,:,1:LM)-PLE(:,:,0:LM-1)*PKE(:,:,0:LM-1)) &
!               / ((1.+MAPL_KAPPA)*DP)
      PK       = (PLO/1000.)**(MAPL_RGAS/MAPL_CP)
      DM       = DP*(1./MAPL_GRAV)
      TEMP     = TH1*PK
      DQS      = GEOS_DQSAT(TEMP, PLO, qsat=QSS)

      QSSFC    = GEOS_QSAT( TS , CNV_PLE(:,:,LM) )


      ZLE(:,:,LM) = 0.
      do L=LM,1,-1
         ZLE(:,:,L-1) = TH (:,:,L) * (1.+MAPL_VIREPS*Q(:,:,L))
         ZLO(:,:,L  ) = ZLE(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
         ZLE(:,:,L-1) = ZLO(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
      end do

      GZLE  = MAPL_GRAV * ZLE
      GZLO  = MAPL_GRAV * ZLO
      RASPRCP  = 0.

      TSFCAIR  = TEMP(:,:,LM) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,LM)
      DTS      = TS - TSFCAIR

      IRAS       = nint(LONS*100)
      JRAS       = nint(LATS*100)


! << Quantity X >> to drive RAS time scales
! 4/25/07 - JTB
!--------------------------------------------      
      QUANX = SUM(  KH(:,:,0:lm-1)*( ZLE(:,:,0:lm-1) - ZLE(:,:,1:LM) ),3 )
      QUANX = SQRT( MAX( QUANX , 0.0 ) )
                                        ! Units = ( m+3 s-1 )+1/2 !!
                                        ! Reasonable dynamic range at 0.5x0.666
                                        ! < 50    (weak PBL), 
                                        ! > 150 > (moderate),
                                        ! > 250   (intense)

 
      QUANX_1 = 250.0
      QUANX_0 =  50.0

      QUANXTP = QUANX_0


      if(associated(RASPBLQ )) RASPBLQ    =  QUANX

      QUANX   = MIN( QUANX , QUANX_1 )
      QUANX   = MAX( QUANX , QUANX_0 )

                  ! Klugey enhancement of RASTIME over high elevation
                  ! jtb - 06/14/07
                  !--------------------------------------------------
      where( (CNV_PLE(:,:,LM) <= 900. ) .AND. (CNV_PLE(:,:,LM) > 800. ) )
        QUANXTP = QUANX_0 + 2.* ( 900. - CNV_PLE(:,:,LM) )
        QUANX   = MAX( QUANX , QUANXTP )
      endwhere

      where( (CNV_PLE(:,:,LM) <=800. ) )
        QUANX = QUANX_1
      endwhere



! Some export work prior to doing calculations
!---------------------------------------------
      if(associated(CNV_MFC)) CNV_MFC(:,:,LM) = 0.
      if(associated(RH1    )) RH1     = Q1/QSS
      if(associated(TVQ0   )) TVQ0    = SUM( (  Q +  QLLS + QLCN + QILS + QICN )*DM , 3 )
      if(associated(TVE0   )) TVE0    = SUM( (  MAPL_CP*TEMP + MAPL_ALHL*Q           & 
                                             -  MAPL_ALHF*(QILS+QICN) )*DM , 3 )
      if(associated(DCPTE  )) DCPTE   = SUM( MAPL_CP*TEMP*DM , 3 )
      if(associated(DQDTCN )) DQDTCN  = Q1
      if(associated(DTHDTCN)) DTHDTCN = TH1

      if(associated(DQLDT  )) DQLDT   = QLLS + QLCN
      if(associated(DQIDT  )) DQIDT   = QILS + QICN

      if(associated(DDF_RH1   )) allocate( DDF_RH1z(LM) )
      if(associated(DDF_RH2   )) allocate( DDF_RH2z(LM) )
      if(associated(DDF_MUPH  )) allocate( DDF_MUPHz(LM) )
      if(associated(DDF_BYNC  )) allocate( DDF_BYNCz(LM) )
      if(associated(DDF_QVC   )) allocate( DDF_QVCz(LM) )
      if(associated(DDF_TC    )) allocate( DDF_TCz(LM) )
      if(associated(DDF_DQDT  )) allocate( DDF_DQDTz(LM) )
      if(associated(DDF_DTDT  )) allocate( DDF_DTDTz(LM) )
      if(associated(DDF_MFC   )) allocate( DDF_MFCz(0:LM) )

      !Trajectory for Moist TLM/ADJ
      if(associated(QILS_tlad)) QILS_tlad = QILS
      if(associated(QLLS_tlad)) QLLS_tlad = QLLS
      if(associated(QICN_tlad)) QICN_tlad = QICN
      if(associated(QLCN_tlad)) QLCN_tlad = QLCN
      if(associated(CFLS_tlad)) CFLS_tlad = CLLS
      if(associated(CFCN_tlad)) CFCN_tlad = CLCN
      if(associated(QLS_tlad )) QLS_tlad  = QILS + QLLS
      if(associated(QCN_tlad )) QCN_tlad  = QICN + QLCN

      ALLOC_BYNCY     = .not.associated(BYNCY     )
      ALLOC_CAPE      = .not.associated(CAPE      )
      ALLOC_INHB      = .not.associated(INHB      )

      if(ALLOC_BYNCY    ) allocate(BYNCY    (IM,JM,LM))
      if(ALLOC_CAPE     ) allocate(CAPE     (IM,JM   ))
      if(ALLOC_INHB     ) allocate(INHB     (IM,JM   ))

      call BUOYANCY( TEMP, TH1, Q1, QSS, DQS, DM, ZLO, BYNCY, CAPE, INHB, IM, JM, LM )     

      K0 = LM
      ICMIN    = max(1,count(PREF < PMIN_DET))
      KCBLMIN  =       count(PREF < PMIN_CBL)
      levs925  = max(1,count(PREF < 92500.))


      KE0 = 0.5*( U1**2 + V1**2 )

      call MAPL_TimerOff(STATE,"-MISC")

! Do convection on IDIM soundings
!-------------------------------------------------

      call MAPL_TimerOn (STATE,"-PRE_RAS")

! Find Cloud base level (LCB), strap if desired
! and reset K0 accordingly
!-----------------------------------------------



! Copy tracers to local array
!----------------------------

      KK=0
      do K=1,KM
         if(IS_FRIENDLY(K)) then
            KK = KK+1
            !PRINT *, "*******TESTING: QNAME, FSCAV_, FSCAV********"
            FSCAV(KK) = FSCAV_(K)
            !PRINT *, QNAMES(K), FSCAV_(K), FSCAV(KK)
            XHO(:,:,:,KK) = TRPtrs(K)%Q(:,:,:)
            if(associated(TRIPtrs(K)%Q)) then
               XHOI(:,:,:,KK) = TRPtrs(K)%Q(:,:,:)
            end if
         end if
      end do

! Determine how to do cloud base 
!-------------------------------

      KLCL = FINDLCL( TH1, Q1, PLO, PK, IM, JM, LM )
      KLFC = FINDLFC( BYNCY, IM, JM, LM )
!!    KPBL = FINDPBL( KH, IM, JM, LM )
!! Set subcloud layer height to one level below PBL height level
!!   make sure subcloud layer is at least 2 levels thick
      do j = 1,jm
       do i = 1,im
        if(nint(kpblin(i,j)).ne.0) then
         kpbl(i,j) = min(nint(kpblin(i,j))+1,LM-1)
        else
         kpbl(i,j) = LM-1
        endif
       enddo
      enddo

      SIGE = PREF/PREF(LM) ! this should eventually change

      do I=1,IM
         do J=1,JM

            SELECT CASE( CBL_METHOD )

            CASE( 1 )
               KCBL(I,J)   =  K0 - KSTRAP
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J):K0) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J):K0) = 1.0 

            CASE( 2 )
               KCBL(I,J)   = KLCL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J):K0) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J):K0) = 1.0 

            CASE ( 3 )
               KCBL(I,J)   = KPBL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J)   ) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J)   ) = 1.0 

            CASE ( 4 )
               KCBL(I,J)   = KLCL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J)   ) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J)   ) = 1.0 
               do kii = kcbl(i,j)+1,k0
                  WGT0(I,J,kii  )      = wgt1(I,J,kii-1)*0.6
                  WGT1(I,J,kii  )      = wgt1(I,J,kii-1)*0.6  
               end do
 
            CASE ( 5 )
               KCBL(I,J)   = KPBL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J)   ) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J)   ) = 1.0 
               do kii = kcbl(i,j)+1,k0
                  WGT0(i,j,kii  )      = wgt1(i,j,kii-1)*0.6
                  WGT1(i,j,kii  )      = wgt1(i,j,kii-1)*0.6  
               end do 

            CASE( 6 )
               KCBL(I,J)   = KPBL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J):K0) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J):K0) = 1.0 

            END SELECT

            ZCBLx(I,J) = ZLE( I, J, KCBL(I,J)-1 )

            if(associated(ZLCL)) then
               if(KLCL(I,J)>0) then
                  ZLCL(I,J) = ZLE(I,J,KLCL(I,J)-1)
               else
                  ZLCL(I,J) = MAPL_UNDEF
               end if
            end if
 
            if(associated(ZLFC)) then
               if(KLFC(I,J)>0) then
                  ZLFC(I,J) = ZLE(I,J,KLFC(I,J)-1)
               else
                  ZLFC(I,J) = MAPL_UNDEF
               end if
            end if

         end do
      end do

! temp kluge to differentiate convective intensities according to PBL intensity 
! via Deep, slow  - RASPARAMS(6), and Shallow, fast - RASPARAMS(5) timescales (jtb 4/25/07).
! QUANX calculated above. MIN'ed MAX'ed at QUANX_{0,1}
! ---------------------------------------------------------------------

      RASPARAMS( 5)    = RASAL1 
      RASPARAMS( 6)    = RASAL2
! Remove dependance on kh
!!!         RASPARAMS( 6)    = ( RASAL1 * ( QUANX(i,j) - QUANX_0 ) + RASAL2 * ( QUANX_1 - QUANX(i,j) ) ) & 
!!!                                  /       ( QUANX_1 - QUANX_0 )

      if(associated(RASTIME )) RASTIME = RASPARAMS(6)


! Cheat by adding a kick to CB temp and q
! ---------------------------------------
      TPERT  = CBL_TPERT * ( TS - ( TEMP(:,:,LM)+ MAPL_GRAV*ZLO(:,:,LM)/MAPL_CP )  ) 
      QPERT  = CBL_QPERT * ( QSSFC - Q(:,:,LM) )          
      TPERT  = MAX( TPERT , 0.0 )
      QPERT  = MAX( QPERT , 0.0 )

      where (FRLAND<0.1) 
         TPERT = MIN( TPERT , CBL_TPERT_MXOCN ) ! ocean
      elsewhere
         TPERT = MIN( TPERT , CBL_TPERT_MXLND ) ! land
      end where


! Compute initial mass loading for aerosols; CAR 12/19/08
! -------------------------------------------------------
      !! First initialize everything to zero, just in case
      if(associated(DDUDT)) DDUDT =  0.0
      if(associated(DSSDT)) DSSDT =  0.0
      if(associated(DBCDT)) DBCDT =  0.0
      if(associated(DOCDT)) DOCDT =  0.0
      if(associated(DSUDT)) DSUDT =  0.0
      if(associated(DDUDTcarma)) DDUDTcarma =  0.0
      if(associated(DSSDTcarma)) DSSDTcarma =  0.0
      
      !! Now loop over tracers and accumulate initial column loading
      !! tendency  kg/m2/s CAR
      
      KK=0
      do K=1,KM
         if(IS_FRIENDLY(K)) then
            KK = KK + 1
            QNAME = trim(QNAMES(K))
            CNAME = trim(CNAMES(K))
            if(CNAME == 'GOCART') then   ! Diagnostics for GOCART tracers
               SELECT CASE (QNAME(1:3))
               CASE ('du0')
                  if(associated(DDUDT)) then
                     DDUDT = DDUDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST)
                  end if
               CASE ('ss0')
                  if(associated(DSSDT)) then
                     DSSDT = DSSDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST)
                  end if
               CASE ('BCp')
                  if(associated(DBCDT)) then
                     DBCDT = DBCDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST)
                  end if
               CASE ('OCp')
                  if(associated(DOCDT)) then
                     DOCDT = DOCDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST)
                  end if
               CASE ('SO4')
                  if(associated(DSUDT)) then
                     DSUDT = DSUDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST)
                  end if
               END SELECT
            endif
            if(CNAME == 'CARMA') then   ! Diagnostics for CARMA tracers
               ! Check name to see if it is a "pc" element
               ENAME = ''
               ind= index(QNAME, '::')
               if (ind> 0) then
                  ENAME = trim(QNAME(ind+2:ind+3))  ! Component name (e.g., GOCART, CARMA)
                  if(ENAME == 'pc') then
                     SELECT CASE (QNAME(1:4))
                     CASE ('dust') ! CARMA DUST
                        if(associated(DDUDTcarma)) then
                           DDUDTcarma = DDUDTcarma + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST)
                        end if
                     CASE ('seas') ! CARMA SEASALT
                        if(associated(DSSDTcarma)) then
                           DSSDTcarma = DSSDTcarma + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST)
                        end if
                     END SELECT
                  endif
               endif
            endif
         end if
      end do

      ! Myong-In I just   
      ! put these 100s    
      ! in. Think about  
      ! whether you want           |                         |
      ! to keep them               V                         V
      SEEDRAS(:,:,1) = 1000000 * ( 100*TEMP(:,:,LM)   - INT( 100*TEMP(:,:,LM) ) )
      SEEDRAS(:,:,2) = 1000000 * ( 100*TEMP(:,:,LM-1) - INT( 100*TEMP(:,:,LM-1) ) )

      call MAPL_TimerOff(STATE,"-PRE_RAS")

      call MAPL_TimerOn (STATE,"-RAS")

      ! MATMAT Export out the inputs before RAS needed for RAStest
      call MAPL_GetPointer(EXPORT, THOI,   'THOI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QHOI,   'QHOI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QSSI,   'QSSI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQSI,   'DQSI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PLEI,   'PLEI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TPERTI, 'TPERTI', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KCBLI,  'KCBLI' , RC=STATUS); VERIFY_(STATUS)

      if(associated(THOI  )) THOI   = TH1
      if(associated(QHOI  )) QHOI   = Q1
      if(associated(QSSI  )) QSSI   = QSS
      if(associated(DQSI  )) DQSI   = DQS
      if(associated(PLEI  )) PLEI   = CNV_PLE
      if(associated(TPERTI)) TPERTI = TPERT
      if(associated(KCBLI )) KCBLI  = KCBL

      ! temp kluge to differentiate ocean,land convective autoc (jtb 6/29/05)
      ! ---------------------------------------------------------------------
      where (FRLAND<0.1) 
         CO_AUTO = AUTOC_CN_OCN   ! ocean value
      elsewhere
         CO_AUTO = AUTOC_CN_LAND  ! land value
      end where

      IDIM       = IM*JM
      IRUN       = IM*JM

      call MAPL_TimerOn (STATE,"--RAS_RUN",RC=STATUS)
      VERIFY_(STATUS)

      Q1_dh = Q1
      U1_dh = U1
      V1_dh = V1
      TH1_dh = TH1
      QILS_dh = QILS
      QLLS_dh = QLLS
      QICN_dh = QICN
      QLCN_dh = QLCN
      CFLS_dh = CLLS
      CFCN_dh = CLCN

      call RASE(                        &                    
                 IDIM                 , &
                 IRUN                 , &                 
                 K0                   , &
                 ICMIN                , &
                 DT_MOIST             , &  !!where? see below.
                 MAPL_CP              , &
                 MAPL_ALHL            , &
                 MAPL_ALHS            , &
                 MAPL_TICE            , &
                 MAPL_GRAV            , &
                 SEEDRAS              , &
                 IRAS                 , &
                 JRAS                 , &
                 SIGE                 , &
! inputs for CBL
                 KCBL                 , &
                 WGT0                 , &
                 WGT1                 , &
                 ZCBLx                , &
                 MXDIAMx              , &
                 TPERT                , &
                 QPERT                , &
! inputs
                 TH1                  , &
                 Q1                   , & 
                 U1                   , &
                 V1                   , &
                 QSS                  , & 
                 DQS                  , &
! Pass in CO_AUTO
                 CO_AUTO              , &
! - new for ras 2
                 PK                   , &
                 PLO                  , &
                 GZLO                 , &   
                 GZLE                 , &   
                 QLCN                 , &   
                 QICN                 , &   
!
                 CNV_PLE              , &
                 PKE                  , &
! outputs
                 CNV_DQLDT            , &   ! -> progno_cloud
                 CNV_MF0              , &   ! -> diag
                 CNV_MFD              , &   ! -> progno_cloud
                 CNV_MFC              , &   ! -> diag
                 CNV_PRC3             , &   ! -> progno_cloud 
                 CNV_UPDF             , &   ! -> progno_cloud
                 CNV_CVW              , &   ! 
                 CNV_QC               , &   ! -> progno_cloud ???
                 ENTLAM               , &   ! -> diag
                 CLCN                 , &   ! -> upd if RAS-2 
                 HHO                  , &
                 HSO                  , &   ! -> diag
                 RASPRCP              , &
! params
                 RASPARAMS            , &
                 ITRCR                , &
                 irccode              , &
                 XHO = XHO            , &
                 TRIEDLEV_DIAG = trdlx, &
                 FSCAV  = FSCAV       , &
                 DISSKE = KEX           )

      Q1 = Q1_dh
      U1 = U1_dh
      V1 = V1_dh
      TH1 = TH1_dh
      QILS = QILS_dh
      QLLS = QLLS_dh
      QICN = QICN_dh
      QLCN = QLCN_dh
      CLLS = CFLS_dh
      CLCN = CFCN_dh

      call MAPL_TimerOff(STATE,"--RAS_RUN",RC=STATUS)
      VERIFY_(STATUS)

      if(associated(MXDIAM )) MXDIAM  = MXDIAMx
      if(associated(RCCODE )) RCCODE  = 1.0*IRCCODE
      if(associated(TRIEDLV)) TRIEDLV = TRDLX
      if(associated(ZCBL   )) ZCBL    = ZCBLx

      call MAPL_TimerOff(STATE,"-RAS")

      call MAPL_TimerOn (STATE,"-POST_RAS")

!     Compute new mass loading for aerosols; CAR 12/19/08
!     -----------------------------------------------------
!     Loop over tracers, accumulate, subtract new from initial, divide
!     time to get kg/m2/s
      
      KK=0
      do K=1,KM
         if(IS_FRIENDLY(K)) then
            KK = KK + 1
            QNAME = trim(QNAMES(K))
            CNAME = trim(CNAMES(K))
            if(CNAME == 'GOCART') then   ! Diagnostics for GOCART tracers
               SELECT CASE (QNAME(1:3))
               CASE ('du0')
                  if(associated(DDUDT)) then
                     DDUDT = DDUDT - (sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST))
                     where (DDUDT <= 0) DDUDT = 0.0
                  end if
               CASE ('ss0')
                  if(associated(DSSDT)) then
                     DSSDT = DSSDT - (sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST))
                     where (DSSDT <= 0) DSSDT = 0.0
                  end if
               CASE ('BCp')
                  if(associated(DBCDT)) then
                     DBCDT = DBCDT - (sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST))
                     where (DBCDT <= 0) DBCDT = 0.0
                  end if
               CASE ('OCp')
                  if(associated(DOCDT)) then
                     DOCDT = DOCDT - (sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST))
                     where (DOCDT <= 0) DOCDT = 0.0
                  end if
               CASE ('SO4')
                  if(associated(DSUDT)) then
                     DSUDT = DSUDT - (sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST))
                     where (DSUDT <= 0) DSUDT = 0.0
                  end if
               END SELECT
            endif
         end if
         if(CNAME == 'CARMA') then   ! Diagnostics for CARMA tracers
            ! Check name to see if it is a "pc" element
            ENAME = ''
            ind= index(QNAME, '::')
            if (ind> 0) then
               ENAME = trim(QNAME(ind+2:ind+3))  ! Component name (e.g., GOCART, CARMA)
               if(ENAME == 'pc') then
                  SELECT CASE (QNAME(1:4))
                  CASE ('dust') ! CARMA DUST
                     if(associated(DDUDTcarma)) then
                        DDUDTcarma = DDUDTcarma - (sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST))
                        where (DDUDTcarma <= 0) DDUDTcarma = 0.0
                     end if
                  CASE ('seas') ! CARMA SEASALT
                     if(associated(DSSDTcarma)) then
                        DSSDTcarma = DSSDTcarma - (sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)/(MAPL_GRAV*DT_MOIST))
                        where (DSSDTcarma <= 0) DSSDTcarma = 0.0
                     end if
                  END SELECT
               endif
            endif
         endif
      end do

#ifdef USEDDF

      DDF_CTL%ALPHA_DDF_UDF  = 0.5
      DDF_CTL%AREAL_FRACTION = 0.1
       
      DDF_CTL%PARTITION%TYPE  = "HEIGHT_DEP"
      DDF_CTL%PARTITION%MAX_RATIO = 0.9
      DDF_CTL%PARTITION%MIN_RATIO = 0.1
      DDF_CTL%THETA_IS_THVAR      = .TRUE.

      do I = 1, IM
         do J = 1, JM

            ! Call to ddf v 1.10.2.4 (tag=merra_dndrft_03)
            ! needs to change if more recent ddf used
            !----------------------------------------------
            call DDF1( DDF_CTL , DT_MOIST , TH1(I,J,:) , & 
                       Q1(I,J,:), QLCN(I,J,:), QICN(I,J,:), & 
                       CNV_PRC3(I,j,:), CNV_HAIL(I,j,:), &
                       ZLO(I,J,:), ZLE(I,J,:), PLO(I,J,:),  & 
                       CNV_PLE(I,J,:),  PK(I,J,:), CNV_MFC(I,J,:) ,  & 
                       DDF_ZSCALEz, DDF_DQDTz , DDF_DTDTz, DDF_MFCz )

            if(associated(DDF_DQDT   ))  DDF_DQDT   (I,J,:)  = DDF_DQDTz
            if(associated(DDF_DTDT   ))  DDF_DTDT   (I,J,:)  = DDF_DTDTz
            if(associated(DDF_MFC    ))  DDF_MFC    (I,J,:)  = DDF_MFCz
            if(associated(DDF_ZSCALE ))  DDF_ZSCALE (I,J)    = DDF_ZSCALEz

         end do
      end do

#endif

! Fill in tracer tendencies
!--------------------------

      KK=0
      do K=1,KM
         if(IS_FRIENDLY(K)) then
            KK = KK+1
            TRPtrs(K)%Q(:,:,:) =  XHO(:,:,:,KK)
            if(associated(TRIPtrs(K)%Q)) then
               !PRINT *, "TRIPtrs is associated"
               TRIPtrs(K)%Q(:,:,:) = (XHO(:,:,:,KK) - XHOI(:,:,:,KK)) / DT_MOIST
               if(IS_WEIGHTED(K)) then
                  TRIPtrs(K)%Q(:,:,:) = TRIPtrs(K)%Q(:,:,:)*DP(:,:,:)
               end if
            end if
         end if
      end do

      call MAPL_TimerOff(STATE,"-POST_RAS")

      call MAPL_TimerOn (STATE,"-MISC")

      !Trajectory for Moist TLM/ADJ
      if(associated(TSUR_tlad)) TSUR_tlad = TS
      if(associated(KCBL_tlad)) KCBL_tlad = KCBL
      if(associated(KHu_tlad)) then
         KHu_tlad = -1.0
         DO i = 1,IM
            DO j = 1,JM
               DO l = 0,LM
                  if (KH(i,j,l) .gt. 2.0) then 
                     KHu_tlad(i,j) = l * 1.0
                     exit
                  endif
               endDO
            endDO
         endDO
      endif
      if(associated(KHl_tlad)) then
         KHl_tlad = -1.0
         DO i = 1,IM
            DO j = 1,JM
               DO l = LM,0,-1
                  if (KH(i,j,l) .gt. 2.0) then 
                     KHl_tlad(i,j) = l * 1.0
                     exit
                  endif
               endDO
            endDO
         endDO
      endif

      if(associated(QVRAS  )) QVRAS    = Q1
      if(associated(THRAS  )) THRAS    = TH1
      if(associated(URAS   )) URAS     = U1
      if(associated(VRAS   )) VRAS     = V1

      if(associated(DQRC   )) DQRC     = CNV_PRC3           / DT_MOIST
      if(associated(DQDTCN )) DQDTCN   = ( Q1  -  DQDTCN  ) / DT_MOIST
      if(associated(DTHDTCN)) DTHDTCN  = ( TH1 -  DTHDTCN ) / DT_MOIST
      if(associated(TVEX   )) TVEX     = SUM( (MAPL_CP*TEMP + MAPL_ALHL*Q)*DM, 3 )
      if(associated(DQCDTCN)) DQCDTCN  = CNV_DQLDT * MAPL_GRAV / ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )

      !-- intermediate water column (all phases and reservoirs)------
      !! TVQX    = SUM( (  Q1 +  QLLS + QLCN + QILS + QICN + CNV_PRC3)*DM + CNV_DQLDT*DT_MOIST , 3 )

      deallocate(IS_FRIENDLY,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(IS_WEIGHTED,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(TRPtrs     ,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(TRIPtrs    ,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(XHO        ,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(XHOI       ,stat=STATUS)
      VERIFY_(STATUS)
! CAR 12/5/08
      deallocate(FSCAV      ,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(FSCAV_     ,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(QNAMES, CNAMES ,stat=STATUS)
      VERIFY_(STATUS)

      call MAPL_GetResource(STATE, CLOUDPARAMS( 1), 'CNV_BETA:',       DEFAULT= 10.0    )
      call MAPL_GetResource(STATE, CLOUDPARAMS( 2), 'ANV_BETA:',       DEFAULT= 4.0     )
      call MAPL_GetResource(STATE, CLOUDPARAMS( 3), 'LS_BETA:',        DEFAULT= 4.0     )
      call MAPL_GetResource(STATE, CLOUDPARAMS( 4), 'RH_CRIT:',        DEFAULT= 1.0     )
      call MAPL_GetResource(STATE, CLOUDPARAMS( 5), 'AUTOC_LS:',       DEFAULT= 2.0e-3  )
      call MAPL_GetResource(STATE, CLOUDPARAMS( 6), 'QC_CRIT_LS:',     DEFAULT= 8.0e-4  )
      call MAPL_GetResource(STATE, CLOUDPARAMS( 7), 'ACCRETION:',      DEFAULT= 2.0     )
      call MAPL_GetResource(STATE, CLOUDPARAMS( 8), 'RAIN_REVAP_FAC:', DEFAULT= 1.0     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(56), 'SNOW_REVAP_FAC:', DEFAULT= 1.0     )
      call MAPL_GetResource(STATE, CLOUDPARAMS( 9), 'VOL_TO_FRAC:',    DEFAULT= -1.0    )
      call MAPL_GetResource(STATE, CLOUDPARAMS(10), 'SUPERSAT:',       DEFAULT= 0.0     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(11), 'SHEAR_EVAP_FAC:', DEFAULT= 1.3     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(12), 'MIN_ALLOW_CCW:',  DEFAULT= 1.0e-9  )
      call MAPL_GetResource(STATE, CLOUDPARAMS(13), 'CCW_EVAP_EFF:',   DEFAULT= 3.3e-4  )
      call MAPL_GetResource(STATE, CLOUDPARAMS(14), 'NSUB_AUTOCONV:',  DEFAULT= 20.     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(15), 'LS_SUND_INTER:',  DEFAULT= 4.8     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(16), 'LS_SUND_COLD:',   DEFAULT= 4.8     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(17), 'LS_SUND_TEMP1:',  DEFAULT= 230.    )
      call MAPL_GetResource(STATE, CLOUDPARAMS(18), 'ANV_SUND_INTER:', DEFAULT= 1.0     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(19), 'ANV_SUND_COLD:',  DEFAULT= 1.0     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(20), 'ANV_SUND_TEMP1:', DEFAULT= 230.    )
      call MAPL_GetResource(STATE, CLOUDPARAMS(21), 'ANV_TO_LS_TIME:', DEFAULT= 14400.  )
      call MAPL_GetResource(STATE, CLOUDPARAMS(22), 'NCCN_WARM:',      DEFAULT= 50.     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(23), 'NCCN_ICE:',       DEFAULT= 0.01    )
      call MAPL_GetResource(STATE, CLOUDPARAMS(24), 'NCCN_ANVIL:',     DEFAULT= 0.1     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(25), 'NCCN_PBL:',       DEFAULT= 200.    )
      call MAPL_GetResource(STATE, CLOUDPARAMS(26), 'DISABLE_RAD:',    DEFAULT= 0.      )
      call MAPL_GetResource(STATE, CLOUDPARAMS(27), 'ICE_SETTLE:',     DEFAULT= 0.      )
      call MAPL_GetResource(STATE, CLOUDPARAMS(28), 'ANV_ICEFALL:',    DEFAULT= 0.5     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(29), 'LS_ICEFALL:',     DEFAULT= 0.5     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(30), 'REVAP_OFF_P:',    DEFAULT= 2000.   )
      call MAPL_GetResource(STATE, CLOUDPARAMS(31), 'CNV_ENVF:',       DEFAULT= 0.8     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(32), 'WRHODEP:',        DEFAULT= 0.5     )
      call MAPL_GetResource(STATE, CLOUDPARAMS(33), 'ICE_RAMP:',       DEFAULT= -40.0   )
      call MAPL_GetResource(STATE, CLOUDPARAMS(34), 'CNV_ICEPARAM:',   DEFAULT= 1.0   )
      call MAPL_GetResource(STATE, CLOUDPARAMS(35), 'CNV_ICEFRPWR:',   DEFAULT= 4.0   )
      call MAPL_GetResource(STATE, CLOUDPARAMS(36), 'CNV_DDRF:',       DEFAULT= 0.0   )
      call MAPL_GetResource(STATE, CLOUDPARAMS(37), 'ANV_DDRF:',       DEFAULT= 0.0   )
      call MAPL_GetResource(STATE, CLOUDPARAMS(38), 'LS_DDRF:',        DEFAULT= 0.0   )
      call MAPL_GetResource(STATE, CLOUDPARAMS(39), 'AUTOC_ANV:',      DEFAULT= 1.0e-3  )
      call MAPL_GetResource(STATE, CLOUDPARAMS(40), 'QC_CRIT_ANV:',    DEFAULT= 8.0e-4  )
      call MAPL_GetResource(STATE, CLOUDPARAMS(41), 'TANHRHCRIT:',     DEFAULT= 1.0     )

! Horizontal resolution dependant defaults for minimum RH crit
      if(imsize.le.200)call MAPL_GetResource(STATE,CLOUDPARAMS(42),'MINRHCRIT:',DEFAULT=0.80,RC=STATUS)
      if(imsize.gt.200 .and. imsize.le.400)  &
                       call MAPL_GetResource(STATE,CLOUDPARAMS(42),'MINRHCRIT:',DEFAULT=0.90,RC=STATUS)
      if(imsize.gt.400.and. imsize.le.800) &
                        call MAPL_GetResource(STATE,CLOUDPARAMS(42),'MINRHCRIT:',DEFAULT=0.93,RC=STATUS)
      if(imsize.gt.800 .and. imsize.le.1600) &
                       call MAPL_GetResource(STATE,CLOUDPARAMS(42),'MINRHCRIT:',DEFAULT=0.95,RC=STATUS)
      if(imsize.gt.1600)call MAPL_GetResource(STATE,CLOUDPARAMS(42),'MINRHCRIT:',DEFAULT=0.97,RC=STATUS)

      call MAPL_GetResource(STATE, CLOUDPARAMS(43), 'MAXRHCRIT:',      DEFAULT= 1.0 )
      call MAPL_GetResource(STATE, CLOUDPARAMS(44), 'PRECIPRAD:',      DEFAULT= 0.0 )
      call MAPL_GetResource(STATE, CLOUDPARAMS(45), 'TURNRHCRIT:',     DEFAULT= 750.0  )
      call MAPL_GetResource(STATE, CLOUDPARAMS(46), 'MAXRHCRITLAND:',DEFAULT= CLOUDPARAMS(42)+0.01  )

      call MAPL_GetResource(STATE, CLOUDPARAMS(47), 'FR_LS_WAT:',  DEFAULT= 1.0 )
      call MAPL_GetResource(STATE, CLOUDPARAMS(48), 'FR_LS_ICE:',  DEFAULT= 1.0 )
      call MAPL_GetResource(STATE, CLOUDPARAMS(49), 'FR_AN_WAT:',  DEFAULT= 0.0 )
      call MAPL_GetResource(STATE, CLOUDPARAMS(50), 'FR_AN_ICE:',  DEFAULT= 0.0 )

      call MAPL_GetResource(STATE, CLOUDPARAMS(51), 'MIN_RL:',  DEFAULT= 10.e-6 )
      call MAPL_GetResource(STATE, CLOUDPARAMS(52), 'MIN_RI:',  DEFAULT= 20.e-6 )
      call MAPL_GetResource(STATE, CLOUDPARAMS(53), 'MAX_RL:',  DEFAULT= 21.e-6 )
      call MAPL_GetResource(STATE, CLOUDPARAMS(54), 'MAX_RI:',  DEFAULT= 40.e-6 )
      call MAPL_GetResource(STATE, CLOUDPARAMS(55), 'RI_ANV:',  DEFAULT= 30.e-6 )
      call MAPL_GetResource(STATE, CLOUDPARAMS(57), 'PDFSHAPE:',  DEFAULT= 1.0 )

      call MAPL_GetResource(STATE, CLOUD_CTL%SCLMFDFR,    'SCLMFDFR:',      DEFAULT= 1.00  )
      call MAPL_GetResource(STATE, CLOUD_CTL%RSUB_RADIUS, 'RSUB_RADIUS:',   DEFAULT= 1.00e-3  )

      call MAPL_TimerOff(STATE,"-MISC")

      call MAPL_TimerOn (STATE,"-CLOUD")

      ! -----------------------------------------
      ! Preliminary calculations for progno_cloud
      ! -----------------------------------------

      ! GPU These are done here due to the backwards loop for ZET
      !     Done here, the progno_cloud code is completely local.

      ! Calculate QST3 and pass in to progno_cloud
      ! ------------------------------------------

      TEMP_0 = TH1 * PK
      DQST3 = GEOS_DQSAT(TEMP_0, PLO, QSAT=QST3)

      ! Calculate QDDF3 and pass in to progno_cloud=>precip3
      ! ----------------------------------------------------

      DZET(:,:,1:LM) = TH1(:,:,1:LM) * (PKE(:,:,1:LM) - PKE(:,:,0:LM-1)) * MAPL_CP/MAPL_GRAV

      MASS(:,:,1:LM) = ( CNV_PLE(:,:,1:LM) - CNV_PLE(:,:,0:LM-1) )*100./MAPL_GRAV

      ZET(:,:,LM+1) = 0.0
      DO K = LM, 1, -1
         ZET(:,:,K) = ZET(:,:,K+1)+DZET(:,:,K)
      END DO

      WHERE ( ZET(:,:,1:LM) < 3000. )
         QDDF3 = -( ZET(:,:,1:LM)-3000. ) * ZET(:,:,1:LM) * MASS
      ELSEWHERE
         QDDF3 = 0.
      END WHERE   

      VMIP = SUM(QDDF3, 3)
      DO K = 1,LM
         QDDF3(:,:,K) = QDDF3(:,:,K) / VMIP
      END DO

      ! Calculate tempor2d based off of levs925 => tempor in RADCOUPLE
      ! --------------------------------------------------------------

      tempor2d = 0.
      do l = levs925, lm
         where (u1(:,:,l).gt.4.) tempor2d(:,:) = 1.
      end do

#ifdef _CUDA

      call MAPL_GetResource(STATE,BLOCKSIZE,'BLOCKSIZE:',DEFAULT=128,RC=STATUS)
      VERIFY_(STATUS)

      Block = dim3(blocksize,1,1)
      Grid = dim3(ceiling(real(IM*JM)/real(blocksize)),1,1)

      ! GPU In order to save time on the GPU, we need to know the logic of all the 
      ! diagnostics that rely on these. When we integrate the code after PROGNO_CLOUD this will be redundant.
      COPY_RSU_CN = ASSOCIATED(RSU_CN) .OR. ASSOCIATED(REVSU_CN)   .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(SUBPZ)
      COPY_RSU_AN = ASSOCIATED(RSU_AN) .OR. ASSOCIATED(REVSU_LSAN) .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(SUBPZ)
      COPY_RSU_LS = ASSOCIATED(RSU_LS) .OR. ASSOCIATED(REVSU_LSAN) .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(SUBPZ)

      COPY_REV_CN = ASSOCIATED(REV_CN) .OR. ASSOCIATED(REVSU_CN)   .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(EVPPZ)
      COPY_REV_AN = ASSOCIATED(REV_AN) .OR. ASSOCIATED(REVSU_LSAN) .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(EVPPZ)
      COPY_REV_LS = ASSOCIATED(REV_LS) .OR. ASSOCIATED(REVSU_LSAN) .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(EVPPZ)

      COPY_ACIL_CN = ASSOCIATED(ACIL_CN) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLIZ)
      COPY_ACIL_AN = ASSOCIATED(ACIL_AN) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLIZ)
      COPY_ACIL_LS = ASSOCIATED(ACIL_LS) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLIZ)

      COPY_ACLL_CN = ASSOCIATED(ACLL_CN) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLLZ)
      COPY_ACLL_AN = ASSOCIATED(ACLL_AN) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLLZ)
      COPY_ACLL_LS = ASSOCIATED(ACLL_LS) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLLZ)

      COPY_PFI_CN = ASSOCIATED(PFI_CN)
      COPY_PFI_AN = ASSOCIATED(PFI_AN) .OR. ASSOCIATED(PFI_LSAN)
      COPY_PFI_LS = ASSOCIATED(PFI_LS) .OR. ASSOCIATED(PFI_LSAN)

      COPY_PFL_CN = ASSOCIATED(PFL_CN)
      COPY_PFL_AN = ASSOCIATED(PFL_AN) .OR. ASSOCIATED(PFL_LSAN)
      COPY_PFL_LS = ASSOCIATED(PFL_LS) .OR. ASSOCIATED(PFL_LSAN)

      COPY_DLPDF = ASSOCIATED(DLPDF) .OR. ASSOCIATED(PDFLZ)
      COPY_DIPDF = ASSOCIATED(DIPDF) .OR. ASSOCIATED(PDFIZ)
      COPY_DLFIX = ASSOCIATED(DLFIX) .OR. ASSOCIATED(EVPCZ)
      COPY_DIFIX = ASSOCIATED(DIFIX) .OR. ASSOCIATED(SUBCZ)

      COPY_RHX    = ASSOCIATED(RHX)
      COPY_AUT    = ASSOCIATED(AUT)    .OR. ASSOCIATED(AUTZ)
      COPY_EVAPC  = ASSOCIATED(EVAPC)  .OR. ASSOCIATED(EVPCZ)
      COPY_SDM    = ASSOCIATED(SDM)    .OR. ASSOCIATED(SDMZ)
      COPY_SUBLC  = ASSOCIATED(SUBLC)  .OR. ASSOCIATED(SUBCZ)
      COPY_FRZ_TT = ASSOCIATED(FRZ_TT) .OR. ASSOCIATED(FRZCZ)
      COPY_FRZ_PP = ASSOCIATED(FRZ_PP) .OR. ASSOCIATED(FRZPZ)
      COPY_DCNVL  = ASSOCIATED(DCNVL)  .OR. ASSOCIATED(CNVLZ)
      COPY_DCNVI  = ASSOCIATED(DCNVI)  .OR. ASSOCIATED(CNVIZ)

      COPY_ALPHT = ASSOCIATED(ALPHT)
      COPY_ALPH1 = ASSOCIATED(ALPH1)
      COPY_ALPH2 = ASSOCIATED(ALPH2)

      COPY_CFPDF = ASSOCIATED(CFPDF)
      COPY_RHCLR = ASSOCIATED(RHCLR)
      COPY_DQRL  = ASSOCIATED(DQRL) .OR. ASSOCIATED(PGENTOT)

      COPY_VFALLICE_AN = ASSOCIATED(VFALLICE_AN)
      COPY_VFALLICE_LS = ASSOCIATED(VFALLICE_LS)
      COPY_VFALLWAT_AN = ASSOCIATED(VFALLWAT_AN)
      COPY_VFALLWAT_LS = ASSOCIATED(VFALLWAT_LS)

      COPY_VFALLRN_AN = ASSOCIATED(VFALLRN_AN)
      COPY_VFALLRN_LS = ASSOCIATED(VFALLRN_LS)
      COPY_VFALLRN_CN = ASSOCIATED(VFALLRN_CN)
      COPY_VFALLSN_AN = ASSOCIATED(VFALLSN_AN)
      COPY_VFALLSN_LS = ASSOCIATED(VFALLSN_LS)
      COPY_VFALLSN_CN = ASSOCIATED(VFALLSN_CN)

      call MAPL_TimerOn (STATE,"--CLOUD_ALLOC",RC=STATUS)
      VERIFY_(STATUS)

      ! ----------------------
      ! Allocate device arrays
      ! ----------------------
  
      ! Inputs
      ! ------
  
      ALLOCATE(PP_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(EXNP_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PPE_DEV(IM*JM,0:LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(KH_DEV(IM*JM,0:LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(FRLAND_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(RMFDTR_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QLWDTR_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(U_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(V_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QST3_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(DZET_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QDDF3_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(TEMPOR_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
  
      ! Inoutputs
      ! ---------
  
      ALLOCATE(TH_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(Q_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QRN_CU_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CNV_UPDFRC_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QLW_LS_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QLW_AN_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QIW_LS_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QIW_AN_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(ANVFRC_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CLDFRC_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
  
      ! Outputs
      ! -------
  
      ALLOCATE(RAD_CLDFRC_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(RAD_QL_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(RAD_QI_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(RAD_QR_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(RAD_QS_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CLDREFFL_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CLDREFFI_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PRELS_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PRECU_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PREAN_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(LSARF_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CUARF_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(ANARF_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(SNRLS_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(SNRCU_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(SNRAN_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
  
      ! Because of use_autoconv_timescale in PROGNO_CLOUD, the PFL and PFI
      ! arrays can "silently" be used as working arrays
      ! They must always be allocated
      ALLOCATE(PFL_CN_DEV(IM*JM,0:LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PFI_CN_DEV(IM*JM,0:LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PFL_AN_DEV(IM*JM,0:LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PFI_AN_DEV(IM*JM,0:LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PFL_LS_DEV(IM*JM,0:LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PFI_LS_DEV(IM*JM,0:LM), STAT=STATUS)
      VERIFY_(STATUS)
  
      ! Diagnostics
      ! -----------

      ALLOCATE(RHX_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(REV_LS_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(REV_AN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(REV_CN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(RSU_LS_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(RSU_AN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(RSU_CN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(ACLL_CN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(ACIL_CN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(ACLL_AN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(ACIL_AN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(ACLL_LS_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(ACIL_LS_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(PDFL_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(PDFI_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(FIXL_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(FIXI_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(AUT_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(EVAPC_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(SDM_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(SUBLC_DEV (IM*JM,LM), STAT=STATUS)
      ALLOCATE(FRZ_TT_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(FRZ_PP_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(DCNVL_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(DCNVI_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(ALPHT_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(ALPH1_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(ALPH2_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(CFPDF_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(RHCLR_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(DQRL_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(VFALLICE_AN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(VFALLICE_LS_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(VFALLWAT_AN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(VFALLWAT_LS_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(VFALLSN_AN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(VFALLSN_LS_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(VFALLSN_CN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(VFALLRN_AN_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(VFALLRN_LS_DEV(IM*JM,LM), STAT=STATUS)
      ALLOCATE(VFALLRN_CN_DEV(IM*JM,LM), STAT=STATUS)

      call MAPL_TimerOff(STATE,"--CLOUD_ALLOC",RC=STATUS)
      VERIFY_(STATUS)

      ! Copy host inputs to device
      ! --------------------------

      call MAPL_TimerOn (STATE,"--CLOUD_DATA",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (STATE,"---CLOUD_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      STATUS = cudaMemcpy(PP_DEV,PLO,IM*JM*LM) 
      STATUS = cudaMemcpy(PPE_DEV,CNV_PLE,IM*JM*(LM+1))
      STATUS = cudaMemcpy(EXNP_DEV,PK,IM*JM*LM) 
      STATUS = cudaMemcpy(FRLAND_DEV,FRLAND,IM*JM) 
      STATUS = cudaMemcpy(KH_DEV,KH,IM*JM*(LM+1))
      STATUS = cudaMemcpy(RMFDTR_DEV,CNV_MFD,IM*JM*LM) 
      STATUS = cudaMemcpy(QLWDTR_DEV,CNV_DQLDT,IM*JM*LM) 
      STATUS = cudaMemcpy(U_DEV,U1,IM*JM*LM) 
      STATUS = cudaMemcpy(V_DEV,V1,IM*JM*LM) 
      STATUS = cudaMemcpy(QST3_DEV,QST3,IM*JM*LM) 
      STATUS = cudaMemcpy(DZET_DEV,DZET,IM*JM*LM) 
      STATUS = cudaMemcpy(QDDF3_DEV,QDDF3,IM*JM*LM) 
      STATUS = cudaMemcpy(TEMPOR_DEV,TEMPOR2D,IM*JM) 

      ! Inoutputs
      ! ---------

      STATUS = cudaMemcpy(QRN_CU_DEV,CNV_PRC3,IM*JM*LM) 
      STATUS = cudaMemcpy(CNV_UPDFRC_DEV,CNV_UPDF,IM*JM*LM) 
      STATUS = cudaMemcpy(TH_DEV,TH1,IM*JM*LM) 
      STATUS = cudaMemcpy(Q_DEV,Q1,IM*JM*LM) 
      STATUS = cudaMemcpy(QLW_LS_DEV,QLLS,IM*JM*LM) 
      STATUS = cudaMemcpy(QLW_AN_DEV,QLCN,IM*JM*LM) 
      STATUS = cudaMemcpy(QIW_LS_DEV,QILS,IM*JM*LM) 
      STATUS = cudaMemcpy(QIW_AN_DEV,QICN,IM*JM*LM) 
      STATUS = cudaMemcpy(ANVFRC_DEV,CLCN,IM*JM*LM) 
      STATUS = cudaMemcpy(CLDFRC_DEV,CLLS,IM*JM*LM) 

      call MAPL_TimerOff(STATE,"---CLOUD_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (STATE,"---CLOUD_DATA_CONST",RC=STATUS)
      VERIFY_(STATUS)
      
      ! ------------------------
      ! Load constants in device
      ! ------------------------

      ! PHYSPARAMS Constants are loaded into constant memory
      ! ----------------------------------------------------

      CNV_BETA      = CLOUDPARAMS(1)  ! Area factor for convective rain showers (non-dim)
      ANV_BETA      = CLOUDPARAMS(2)  ! Area factor for anvil rain showers (non-dim)
      LS_BETA       = CLOUDPARAMS(3)  ! Area factor for Large Scale rain showers (non-dim)
      RH00          = CLOUDPARAMS(4)  ! Critical relative humidity
      C_00          = CLOUDPARAMS(5)
      LWCRIT        = CLOUDPARAMS(6)
      C_ACC         = CLOUDPARAMS(7)
      C_EV_R        = CLOUDPARAMS(8)
      C_EV_S        = CLOUDPARAMS(56)
      CLDVOL2FRC    = CLOUDPARAMS(9)
      RHSUP_ICE     = CLOUDPARAMS(10)
      SHR_EVAP_FAC  = CLOUDPARAMS(11)
      MIN_CLD_WATER = CLOUDPARAMS(12)
      CLD_EVP_EFF   = CLOUDPARAMS(13)
      NSMAX         = INT( CLOUDPARAMS(14)  )
      LS_SDQV2      = CLOUDPARAMS(15)
      LS_SDQV3      = CLOUDPARAMS(16)
      LS_SDQVT1     = CLOUDPARAMS(17)
      ANV_SDQV2     = CLOUDPARAMS(18)
      ANV_SDQV3     = CLOUDPARAMS(19)
      ANV_SDQVT1    = CLOUDPARAMS(20)
      ANV_TO_LS     = CLOUDPARAMS(21)
      N_WARM        = CLOUDPARAMS(22)
      N_ICE         = CLOUDPARAMS(23)
      N_ANVIL       = CLOUDPARAMS(24)
      N_PBL         = CLOUDPARAMS(25)
      DISABLE_RAD   = INT( CLOUDPARAMS(26) )
      ANV_ICEFALL_C = CLOUDPARAMS(28)
      LS_ICEFALL_C  = CLOUDPARAMS(29)
      REVAP_OFF_P   = CLOUDPARAMS(30)
      CNVENVFC      = CLOUDPARAMS(31)
      WRHODEP       = CLOUDPARAMS(32) 
      T_ICE_ALL     = CLOUDPARAMS(33) + MAPL_TICE
      CNVICEPARAM   = CLOUDPARAMS(34)
      ICEFRPWR      = INT( CLOUDPARAMS(35) + .001 )
      CNVDDRFC      = CLOUDPARAMS(36)
      ANVDDRFC      = CLOUDPARAMS(37)
      LSDDRFC       = CLOUDPARAMS(38)
      TANHRHCRIT    = INT( CLOUDPARAMS(41) )
      MINRHCRIT     = CLOUDPARAMS(42)
      MAXRHCRIT     = CLOUDPARAMS(43)
      TURNRHCRIT    = CLOUDPARAMS(45)
      MAXRHCRITLAND = CLOUDPARAMS(46)
      FR_LS_WAT     = INT( CLOUDPARAMS(47) )
      FR_LS_ICE     = INT( CLOUDPARAMS(48) )
      FR_AN_WAT     = INT( CLOUDPARAMS(49) )
      FR_AN_ICE     = INT( CLOUDPARAMS(50) )
      MIN_RL        = CLOUDPARAMS(51)
      MIN_RI        = CLOUDPARAMS(52)
      MAX_RL        = CLOUDPARAMS(53)
      MAX_RI        = CLOUDPARAMS(54)
      RI_ANV        = CLOUDPARAMS(55)
      PDFFLAG       = INT(CLOUDPARAMS(57))

      STATUS = cudaDeviceSynchronize()

      call MAPL_TimerOff(STATE,"---CLOUD_DATA_CONST",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOff(STATE,"--CLOUD_DATA",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (STATE,"--CLOUD_RUN",RC=STATUS)
      VERIFY_(STATUS)

      STATUS = cudaFuncSetCacheConfig('cloudnew_progno_cloud',cudaFuncCachePreferL1)

      call PROGNO_CLOUD<<<Grid, Block>>>(IM*JM,LM,DT_MOIST,CLOUD_CTL%SCLMFDFR)

      STATUS = cudaGetLastError()
      if (STATUS /= 0) then 
         write (*,*) "Error code from PROGNO_CLOUD kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         ASSERT_(.FALSE.)
      end if

      STATUS = cudaDeviceSynchronize()

      call MAPL_TimerOff(STATE,"--CLOUD_RUN",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (STATE,"--CLOUD_DATA",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (STATE,"---CLOUD_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      ! Move device outputs to host
      ! ---------------------------

      ! Outputs
      ! -------

      STATUS = cudaMemcpy(RAD_CF,RAD_CLDFRC_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(RAD_QL,RAD_QL_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(RAD_QI,RAD_QI_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(QRAIN,RAD_QR_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(QSNOW,RAD_QS_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(CLDREFFL,CLDREFFL_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(CLDREFFI,CLDREFFI_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(LS_PRC2,PRELS_DEV,IM*JM)
      STATUS = cudaMemcpy(CN_PRC2,PRECU_DEV,IM*JM)
      STATUS = cudaMemcpy(AN_PRC2,PREAN_DEV,IM*JM)
      STATUS = cudaMemcpy(LS_ARFX,LSARF_DEV,IM*JM)
      STATUS = cudaMemcpy(CN_ARFX,CUARF_DEV,IM*JM)
      STATUS = cudaMemcpy(AN_ARFX,ANARF_DEV,IM*JM)
      STATUS = cudaMemcpy(LS_SNR,SNRLS_DEV,IM*JM)
      STATUS = cudaMemcpy(CN_SNR,SNRCU_DEV,IM*JM)
      STATUS = cudaMemcpy(AN_SNR,SNRAN_DEV,IM*JM)

      ! Inoutputs
      ! ---------

      STATUS = cudaMemcpy(CNV_PRC3,QRN_CU_DEV,IM*JM*LM) 
      STATUS = cudaMemcpy(CNV_UPDF,CNV_UPDFRC_DEV,IM*JM*LM) 
      STATUS = cudaMemcpy(TH1,TH_DEV,IM*JM*LM) 
      STATUS = cudaMemcpy(Q1,Q_DEV,IM*JM*LM) 
      STATUS = cudaMemcpy(QLLS,QLW_LS_DEV,IM*JM*LM) 
      STATUS = cudaMemcpy(QLCN,QLW_AN_DEV,IM*JM*LM) 
      STATUS = cudaMemcpy(QILS,QIW_LS_DEV,IM*JM*LM) 
      STATUS = cudaMemcpy(QICN,QIW_AN_DEV,IM*JM*LM) 
      STATUS = cudaMemcpy(CLCN,ANVFRC_DEV,IM*JM*LM) 
      STATUS = cudaMemcpy(CLLS,CLDFRC_DEV,IM*JM*LM) 

      ! Diagnostics
      ! -----------

      IF (COPY_RHX)     STATUS = cudaMemcpy(RHX_X,RHX_DEV,SIZE(RHX_X))
      IF (COPY_REV_LS)  STATUS = cudaMemcpy(REV_LS_X,REV_LS_DEV,SIZE(REV_LS_X))
      IF (COPY_REV_AN)  STATUS = cudaMemcpy(REV_AN_X,REV_AN_DEV,SIZE(REV_AN_X))
      IF (COPY_REV_CN)  STATUS = cudaMemcpy(REV_CN_X,REV_CN_DEV,SIZE(REV_CN_X))
      IF (COPY_RSU_LS)  STATUS = cudaMemcpy(RSU_LS_X,RSU_LS_DEV,SIZE(RSU_LS_X))
      IF (COPY_RSU_AN)  STATUS = cudaMemcpy(RSU_AN_X,RSU_AN_DEV,SIZE(RSU_AN_X))
      IF (COPY_RSU_CN)  STATUS = cudaMemcpy(RSU_CN_X,RSU_CN_DEV,SIZE(RSU_CN_X))
      IF (COPY_ACLL_CN) STATUS = cudaMemcpy(ACLL_CN_X,ACLL_CN_DEV,SIZE(ACLL_CN_X))
      IF (COPY_ACIL_CN) STATUS = cudaMemcpy(ACIL_CN_X,ACIL_CN_DEV,SIZE(ACIL_CN_X))
      IF (COPY_ACLL_AN) STATUS = cudaMemcpy(ACLL_AN_X,ACLL_AN_DEV,SIZE(ACLL_AN_X))
      IF (COPY_ACIL_AN) STATUS = cudaMemcpy(ACIL_AN_X,ACIL_AN_DEV,SIZE(ACIL_AN_X))
      IF (COPY_ACLL_LS) STATUS = cudaMemcpy(ACLL_LS_X,ACLL_LS_DEV,SIZE(ACLL_LS_X))
      IF (COPY_ACIL_LS) STATUS = cudaMemcpy(ACIL_LS_X,ACIL_LS_DEV,SIZE(ACIL_LS_X))
      IF (COPY_PFL_CN)  STATUS = cudaMemcpy(PFL_CN_X,PFL_CN_DEV,SIZE(PFL_CN_X))
      IF (COPY_PFI_CN)  STATUS = cudaMemcpy(PFI_CN_X,PFI_CN_DEV,SIZE(PFI_CN_X))
      IF (COPY_PFL_AN)  STATUS = cudaMemcpy(PFL_AN_X,PFL_AN_DEV,SIZE(PFL_AN_X))
      IF (COPY_PFI_AN)  STATUS = cudaMemcpy(PFI_AN_X,PFI_AN_DEV,SIZE(PFI_AN_X))
      IF (COPY_PFL_LS)  STATUS = cudaMemcpy(PFL_LS_X,PFL_LS_DEV,SIZE(PFL_LS_X))
      IF (COPY_PFI_LS)  STATUS = cudaMemcpy(PFI_LS_X,PFI_LS_DEV,SIZE(PFI_LS_X))
      IF (COPY_DLPDF)   STATUS = cudaMemcpy(DLPDF_X,PDFL_DEV,SIZE(DLPDF_X))
      IF (COPY_DIPDF)   STATUS = cudaMemcpy(DIPDF_X,PDFI_DEV,SIZE(DIPDF_X))
      IF (COPY_DLFIX)   STATUS = cudaMemcpy(DLFIX_X,FIXL_DEV,SIZE(DLFIX_X))
      IF (COPY_DIFIX)   STATUS = cudaMemcpy(DIFIX_X,FIXI_DEV,SIZE(DIFIX_X))
      IF (COPY_AUT)     STATUS = cudaMemcpy(AUT_X,AUT_DEV,SIZE(AUT_X))
      IF (COPY_EVAPC)   STATUS = cudaMemcpy(EVAPC_X,EVAPC_DEV,SIZE(EVAPC_X))
      IF (COPY_SDM)     STATUS = cudaMemcpy(SDM_X,SDM_DEV,SIZE(SDM_X))
      IF (COPY_SUBLC)   STATUS = cudaMemcpy(SUBLC_X,SUBLC_DEV,SIZE(SUBLC_X))
      IF (COPY_FRZ_TT)  STATUS = cudaMemcpy(FRZ_TT_X,FRZ_TT_DEV,SIZE(FRZ_TT_X))
      IF (COPY_FRZ_PP)  STATUS = cudaMemcpy(FRZ_PP_X,FRZ_PP_DEV,SIZE(FRZ_PP_X))
      IF (COPY_DCNVL)   STATUS = cudaMemcpy(DCNVL_X,DCNVL_DEV,SIZE(DCNVL_X))
      IF (COPY_DCNVi)   STATUS = cudaMemcpy(DCNVi_X,DCNVI_DEV,SIZE(DCNVi_X))
      IF (COPY_ALPHT)   STATUS = cudaMemcpy(ALPHT_X,ALPHT_DEV,SIZE(ALPHT_X))
      IF (COPY_ALPH1)   STATUS = cudaMemcpy(ALPH1_X,ALPH1_DEV,SIZE(ALPH1_X))
      IF (COPY_ALPH2)   STATUS = cudaMemcpy(ALPH2_X,ALPH2_DEV,SIZE(ALPH2_X))
      IF (COPY_CFPDF)   STATUS = cudaMemcpy(CFPDF_X,CFPDF_DEV,SIZE(CFPDF_X))
      IF (COPY_RHCLR)   STATUS = cudaMemcpy(RHCLR_X,RHCLR_DEV,SIZE(RHCLR_X))
      IF (COPY_DQRL)    STATUS = cudaMemcpy(DQRL_X,DQRL_DEV,SIZE(DQRL_X))
      IF (COPY_VFALLICE_AN) STATUS = cudaMemcpy(VFALLICE_AN_X,VFALLICE_AN_DEV,SIZE(VFALLICE_AN_X))
      IF (COPY_VFALLICE_LS) STATUS = cudaMemcpy(VFALLICE_LS_X,VFALLICE_LS_DEV,SIZE(VFALLICE_LS_X))
      IF (COPY_VFALLWAT_AN) STATUS = cudaMemcpy(VFALLWAT_AN_X,VFALLWAT_AN_DEV,SIZE(VFALLWAT_AN_X))
      IF (COPY_VFALLWAT_LS) STATUS = cudaMemcpy(VFALLWAT_LS_X,VFALLWAT_LS_DEV,SIZE(VFALLWAT_LS_X))
      IF (COPY_VFALLSN_AN)  STATUS = cudaMemcpy(VFALLSN_AN_X,VFALLSN_AN_DEV,SIZE(VFALLSN_AN_X))
      IF (COPY_VFALLSN_LS)  STATUS = cudaMemcpy(VFALLSN_LS_X,VFALLSN_LS_DEV,SIZE(VFALLSN_LS_X))
      IF (COPY_VFALLSN_CN)  STATUS = cudaMemcpy(VFALLSN_CN_X,VFALLSN_CN_DEV,SIZE(VFALLSN_CN_X))
      IF (COPY_VFALLRN_AN)  STATUS = cudaMemcpy(VFALLRN_AN_X,VFALLRN_AN_DEV,SIZE(VFALLRN_AN_X))
      IF (COPY_VFALLRN_LS)  STATUS = cudaMemcpy(VFALLRN_LS_X,VFALLRN_LS_DEV,SIZE(VFALLRN_LS_X))
      IF (COPY_VFALLRN_CN)  STATUS = cudaMemcpy(VFALLRN_CN_X,VFALLRN_CN_DEV,SIZE(VFALLRN_CN_X))

      STATUS = cudaDeviceSynchronize()

      call MAPL_TimerOff(STATE,"---CLOUD_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOff(STATE,"--CLOUD_DATA",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (STATE,"--CLOUD_DEALLOC",RC=STATUS)
      VERIFY_(STATUS)

      ! ------------------------
      ! Deallocate device arrays
      ! ------------------------
  
      ! Inputs
      ! ------
  
      DEALLOCATE(PP_DEV)
      DEALLOCATE(EXNP_DEV)
      DEALLOCATE(PPE_DEV)
      DEALLOCATE(KH_DEV)
      DEALLOCATE(FRLAND_DEV)
      DEALLOCATE(RMFDTR_DEV)
      DEALLOCATE(QLWDTR_DEV)
      DEALLOCATE(U_DEV)
      DEALLOCATE(V_DEV)
      DEALLOCATE(QST3_DEV)
      DEALLOCATE(DZET_DEV)
      DEALLOCATE(QDDF3_DEV)
      DEALLOCATE(TEMPOR_DEV)
  
      ! Inoutputs
      ! ---------
  
      DEALLOCATE(TH_DEV)
      DEALLOCATE(Q_DEV)
      DEALLOCATE(QRN_CU_DEV)
      DEALLOCATE(CNV_UPDFRC_DEV)
      DEALLOCATE(QLW_LS_DEV)
      DEALLOCATE(QLW_AN_DEV)
      DEALLOCATE(QIW_LS_DEV)
      DEALLOCATE(QIW_AN_DEV)
      DEALLOCATE(ANVFRC_DEV)
      DEALLOCATE(CLDFRC_DEV)
  
      ! Outputs
      ! -------
  
      DEALLOCATE(RAD_CLDFRC_DEV)
      DEALLOCATE(RAD_QL_DEV)
      DEALLOCATE(RAD_QI_DEV)
      DEALLOCATE(RAD_QR_DEV)
      DEALLOCATE(RAD_QS_DEV)
      DEALLOCATE(CLDREFFL_DEV)
      DEALLOCATE(CLDREFFI_DEV)
      DEALLOCATE(PRELS_DEV)
      DEALLOCATE(PRECU_DEV)
      DEALLOCATE(PREAN_DEV)
      DEALLOCATE(LSARF_DEV)
      DEALLOCATE(CUARF_DEV)
      DEALLOCATE(ANARF_DEV)
      DEALLOCATE(SNRLS_DEV)
      DEALLOCATE(SNRCU_DEV)
      DEALLOCATE(SNRAN_DEV)
  
      ! Because of use_autoconv_timescale in PROGNO_CLOUD, the PFL and PFI
      ! arrays can "silently" be used as working arrays
      ! They must always be allocated
      DEALLOCATE(PFL_CN_DEV)
      DEALLOCATE(PFI_CN_DEV)
      DEALLOCATE(PFL_AN_DEV)
      DEALLOCATE(PFI_AN_DEV)
      DEALLOCATE(PFL_LS_DEV)
      DEALLOCATE(PFI_LS_DEV)

      ! Diagnostics
      ! -----------

      DEALLOCATE(RHX_DEV)
      DEALLOCATE(REV_LS_DEV)
      DEALLOCATE(REV_AN_DEV)
      DEALLOCATE(REV_CN_DEV)
      DEALLOCATE(RSU_LS_DEV)
      DEALLOCATE(RSU_AN_DEV)
      DEALLOCATE(RSU_CN_DEV)
      DEALLOCATE(ACLL_CN_DEV)
      DEALLOCATE(ACIL_CN_DEV)
      DEALLOCATE(ACLL_AN_DEV)
      DEALLOCATE(ACIL_AN_DEV)
      DEALLOCATE(ACLL_LS_DEV)
      DEALLOCATE(ACIL_LS_DEV)
      DEALLOCATE(PDFL_DEV)
      DEALLOCATE(PDFI_DEV)
      DEALLOCATE(FIXL_DEV)
      DEALLOCATE(FIXI_DEV)
      DEALLOCATE(AUT_DEV)
      DEALLOCATE(EVAPC_DEV)
      DEALLOCATE(SDM_DEV)
      DEALLOCATE(SUBLC_DEV)
      DEALLOCATE(FRZ_TT_DEV)
      DEALLOCATE(FRZ_PP_DEV)
      DEALLOCATE(DCNVL_DEV)
      DEALLOCATE(DCNVI_DEV)
      DEALLOCATE(ALPHT_DEV)
      DEALLOCATE(ALPH1_DEV)
      DEALLOCATE(ALPH2_DEV)
      DEALLOCATE(CFPDF_DEV)
      DEALLOCATE(RHCLR_DEV)
      DEALLOCATE(DQRL_DEV)
      DEALLOCATE(VFALLICE_AN_DEV)
      DEALLOCATE(VFALLICE_LS_DEV)
      DEALLOCATE(VFALLWAT_AN_DEV)
      DEALLOCATE(VFALLWAT_LS_DEV)
      DEALLOCATE(VFALLSN_AN_DEV)
      DEALLOCATE(VFALLSN_LS_DEV)
      DEALLOCATE(VFALLSN_CN_DEV)
      DEALLOCATE(VFALLRN_AN_DEV)
      DEALLOCATE(VFALLRN_LS_DEV)
      DEALLOCATE(VFALLRN_CN_DEV)

      STATUS = cudaDeviceSynchronize()

      call MAPL_TimerOff(STATE,"--CLOUD_DEALLOC",RC=STATUS)
      VERIFY_(STATUS)

#else

      call MAPL_TimerOn (STATE,"--CLOUD_RUN",RC=STATUS)
      VERIFY_(STATUS)

      Q1_dh = Q1
      U1_dh = U1
      V1_dh = V1
      TH1_dh = TH1
      QILS_dh = QILS
      QLLS_dh = QLLS
      QICN_dh = QICN
      QLCN_dh = QLCN
      CFLS_dh = CLLS
      CFCN_dh = CLCN

      call  PROGNO_CLOUD (                    &
                          IM*JM, LM         , &
                          DT_MOIST          , &
                          PLO               , &
                          CNV_PLE           , &
                          PK                , &
                          FRLAND            , &   ! <- surf
                          KH                , &   ! <- turb
                          CNV_MFD           , &   ! <- ras
                          CNV_DQLDT         , &   ! <- ras              
                          CNV_PRC3          , &   ! <- ras   
                          CNV_UPDF          , &   ! <- ras
                          U1                , &
                          V1                , & 
                          TH1               , &              
                          Q1                , &
                          QLLS              , &
                          QLCN              , &
                          QILS              , &
                          QICN              , &
                          CLCN              , &
                          CLLS              , &
                          RAD_CF            , &
                          RAD_QL            , &
                          RAD_QI            , &
                          QRAIN             , &
                          QSNOW             , &
                          CLDREFFL          , &
                          CLDREFFI          , &
                          LS_PRC2           , &
                          CN_PRC2           , &
                          AN_PRC2           , &
                          LS_ARFX           , &
                          CN_ARFX           , &
                          AN_ARFX           , &
                          LS_SNR            , &
                          CN_SNR            , &
                          AN_SNR            , &
                          CLOUDPARAMS       , &
                          CLOUD_CTL%SCLMFDFR, &
                          QST3              , &
                          DZET              , &
                          QDDF3             , &
                          ! Diagnostics
                          RHX_X             , &
                          REV_LS_X          , &
                          REV_AN_X          , &
                          REV_CN_X          , &
                          RSU_LS_X          , &
                          RSU_AN_X          , &
                          RSU_CN_X          , &
                          ACLL_CN_X,ACIL_CN_X   , &
                          ACLL_AN_X,ACIL_AN_X   , &
                          ACLL_LS_X,ACIL_LS_X   , &
                          PFL_CN_X,PFI_CN_X     , &
                          PFL_AN_X,PFI_AN_X     , &
                          PFL_LS_X,PFI_LS_X     , &
                          DLPDF_X,DIPDF_X,DLFIX_X,DIFIX_X,    &
                          AUT_X, EVAPC_X , SDM_X , SUBLC_X ,  &
                          FRZ_TT_X, DCNVL_X, DCNVi_X,       &
                          ALPHT_X, ALPH1_X, ALPH2_X, CFPDF_X, &
                          RHCLR_X,                      &
                          DQRL_X, FRZ_PP_X,               &
                          VFALLICE_AN_X,VFALLICE_LS_X,    &
                          VFALLWAT_AN_X,VFALLWAT_LS_X,    &
                          VFALLSN_AN_X,VFALLSN_LS_X,VFALLSN_CN_X,  &
                          VFALLRN_AN_X,VFALLRN_LS_X,VFALLRN_CN_X,  &
                          ! End diagnostics
                          TEMPOR2D)

      Q1 = Q1_dh
      U1 = U1_dh
      V1 = V1_dh
      TH1 = TH1_dh
      QILS = QILS_dh
      QLLS = QLLS_dh
      QICN = QICN_dh
      QLCN = QLCN_dh
      CLLS = CFLS_dh
      CLCN = CFCN_dh


      VERIFY_(STATUS)

      call MAPL_TimerOff(STATE,"--CLOUD_RUN",RC=STATUS)
      VERIFY_(STATUS)

#endif

      IF (ASSOCIATED(RSU_CN)) RSU_CN = RSU_CN_X
      IF (ASSOCIATED(RSU_AN)) RSU_AN = RSU_AN_X
      IF (ASSOCIATED(RSU_LS)) RSU_LS = RSU_LS_X

      IF (ASSOCIATED(REV_CN)) REV_CN = REV_CN_X
      IF (ASSOCIATED(REV_AN)) REV_AN = REV_AN_X
      IF (ASSOCIATED(REV_LS)) REV_LS = REV_LS_X

      IF (ASSOCIATED(ACIL_CN)) ACIL_CN = ACIL_CN_X
      IF (ASSOCIATED(ACIL_AN)) ACIL_AN = ACIL_AN_X
      IF (ASSOCIATED(ACIL_LS)) ACIL_LS = ACIL_LS_X

      IF (ASSOCIATED(ACLL_CN)) ACLL_CN = ACLL_CN_X
      IF (ASSOCIATED(ACLL_AN)) ACLL_AN = ACLL_AN_X
      IF (ASSOCIATED(ACLL_LS)) ACLL_LS = ACLL_LS_X

      IF (ASSOCIATED(PFI_CN)) PFI_CN = PFI_CN_X
      IF (ASSOCIATED(PFI_AN)) PFI_AN = PFI_AN_X
      IF (ASSOCIATED(PFI_LS)) PFI_LS = PFI_LS_X

      IF (ASSOCIATED(PFL_CN)) PFL_CN = PFL_CN_X
      IF (ASSOCIATED(PFL_AN)) PFL_AN = PFL_AN_X
      IF (ASSOCIATED(PFL_LS)) PFL_LS = PFL_LS_X

      IF (ASSOCIATED(DLPDF)) DLPDF = DLPDF_X
      IF (ASSOCIATED(DIPDF)) DIPDF = DIPDF_X
      IF (ASSOCIATED(DLFIX)) DLFIX = DLFIX_X
      IF (ASSOCIATED(DIFIX)) DIFIX = DIFIX_X

      IF (ASSOCIATED(RHX))       RHX = RHX_X
      IF (ASSOCIATED(AUT))       AUT = AUT_X
      IF (ASSOCIATED(EVAPC))   EVAPC = EVAPC_X
      IF (ASSOCIATED(SDM))       SDM = SDM_X
      IF (ASSOCIATED(SUBLC))   SUBLC = SUBLC_X
      IF (ASSOCIATED(FRZ_TT)) FRZ_TT = FRZ_TT_X
      IF (ASSOCIATED(FRZ_PP)) FRZ_PP = FRZ_PP_X
      IF (ASSOCIATED(DCNVL))   DCNVL = DCNVL_X
      IF (ASSOCIATED(DCNVI))   DCNVI = DCNVI_X

      IF (ASSOCIATED(ALPHT)) ALPHT = ALPHT_X
      IF (ASSOCIATED(ALPH1)) ALPH1 = ALPH1_X
      IF (ASSOCIATED(ALPH2)) ALPH2 = ALPH2_X

      IF (ASSOCIATED(CFPDF)) CFPDF = CFPDF_X
      IF (ASSOCIATED(RHCLR)) RHCLR = RHCLR_X
      IF (ASSOCIATED(DQRL))   DQRL = DQRL_X

      IF (ASSOCIATED(VFALLICE_AN)) VFALLICE_AN = VFALLICE_AN_X
      IF (ASSOCIATED(VFALLICE_LS)) VFALLICE_LS = VFALLICE_LS_X
      IF (ASSOCIATED(VFALLWAT_AN)) VFALLWAT_AN = VFALLWAT_AN_X
      IF (ASSOCIATED(VFALLWAT_LS)) VFALLWAT_LS = VFALLWAT_LS_X

      IF (ASSOCIATED(VFALLRN_AN)) VFALLRN_AN = VFALLRN_AN_X
      IF (ASSOCIATED(VFALLRN_LS)) VFALLRN_LS = VFALLRN_LS_X
      IF (ASSOCIATED(VFALLRN_CN)) VFALLRN_CN = VFALLRN_CN_X
      IF (ASSOCIATED(VFALLSN_AN)) VFALLSN_AN = VFALLSN_AN_X
      IF (ASSOCIATED(VFALLSN_LS)) VFALLSN_LS = VFALLSN_LS_X
      IF (ASSOCIATED(VFALLSN_CN)) VFALLSN_CN = VFALLSN_CN_X

      call MAPL_TimerOff(STATE,"-CLOUD")

      call MAPL_TimerOn (STATE,"-MISC")

      RAD_QV   = max( Q1 , 0. )

      !-----------------------------------------
      !! kluge in some skewness for PBL clouds 
      !! where DTS:=TS-TSFCAIR > 0. i.e., unstable

      CFPBL = RAD_CF
      do L = 1,LM
         where( ( KH(:,:,L-1) > 5. ) .AND. (DTS > 0.) ) 
             CFPBL(:,:,L) = RAD_CF(:,:,L)**3
         endwhere
      enddo

      where( CFPBL > 0.0 ) 
       RAD_QL = RAD_QL*(RAD_CF/CFPBL)
       RAD_QI = RAD_QI*(RAD_CF/CFPBL)
       QRAIN  = QRAIN *(RAD_CF/CFPBL)
       QSNOW  = QSNOW *(RAD_CF/CFPBL)
      endwhere

      RAD_CF = CFPBL
      RAD_QL = MIN( RAD_QL , 0.001 )  ! Still a ridiculously large
      RAD_QI = MIN( RAD_QI , 0.001 )  ! value.
      QRAIN  = MIN( QRAIN , 0.01 )  ! value.
      QSNOW  = MIN( QSNOW , 0.01 )  ! value.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set rain water for radiation to 0 if preciprad flag is off (set to 0)
      if(cloudparams(44).eq.0.) then
       RAD_QR = 0.
       RAD_QS = 0.
      else
       RAD_QR = QRAIN
       RAD_QS = QSNOW
      endif

      if (associated(QRTOT)) QRTOT = QRAIN*RAD_CF
      if (associated(QSTOT)) QSTOT = QSNOW*RAD_CF



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CLDREFFR = 100.e-6
      CLDREFFS = 140.e-6
      CN_PRC2 = CN_PRC2 + RASPRCP

      CN_PRC2 = max(CN_PRC2 , 0.)
      LS_PRC2 = max(LS_PRC2 , 0.)
      AN_PRC2 = max(AN_PRC2 , 0.)
      CN_SNR  = max(CN_SNR  , 0.)
      LS_SNR  = max(LS_SNR  , 0.)
      AN_SNR  = max(AN_SNR  , 0.)

      TEMP    = TH1*PK

! Clean up Relative Humidity where RH > 110%
!---------------------------------------------
      if (CLEANUP_RH == 1)  then
          
           RHEXCESS = 1.1
           DQSDT = GEOS_DQSAT (TEMP , PLO, QSAT=QSS ) 
           where ( Q1 > RHEXCESS*QSS )
                   DQS = (Q1 - RHEXCESS*QSS)/( 1.0 + RHEXCESS*DQSDT*MAPL_ALHL/MAPL_CP )
           elsewhere
                   DQS = 0.0
           endwhere
           
           Q1      =  Q1    - DQS
           TEMP    =  TEMP  + (MAPL_ALHL/MAPL_CP)*DQS
           TH1     =  TEMP  / PK

           ER_PRC2 =  SUM( DQS * DM, 3)/DT_MOIST
           LS_PRC2 =  LS_PRC2 + ER_PRC2
      ELSE
           ER_PRC2 =  0.0
      ENDIF

      TPREC = CN_PRC2 + LS_PRC2 + AN_PRC2 + &
              CN_SNR  + LS_SNR  + AN_SNR 


! Clean up any negative specific humidity
!-----------------------------------------

     call FILLQ2ZERO( Q1, DM, FILLQ )

! Outputs for Spec
!----------------------------------------------

      if (associated(ACR_TOT   ))   ACR_TOT    = ACLL_CN_X + ACIL_CN_X + ACLL_AN_X + ACIL_AN_X + ACLL_LS_X + ACIL_LS_X 
      if (associated(REVSU_CN  ))   REVSU_CN   = REV_CN_X  + RSU_CN_X 
      if (associated(REVSU_LSAN))   REVSU_LSAN = REV_AN_X  + RSU_AN_X  + REV_LS_X  + RSU_LS_X 
      if (associated(PFL_LSAN  ))   then
                                    PFL_LSAN         = PFL_AN_X  + PFL_LS_X  
                                    PFL_LSAN(:,:,LM) = PFL_LSAN(:,:,LM) + ER_PRC2
      endif
      if (associated(PFI_LSAN  ))   PFI_LSAN   = PFI_AN_X  + PFI_LS_X  

      if (associated(PGENTOT   ))   PGENTOT    =  SUM( ( DQRL_X    + DQRC          &  
                                                       + ACLL_CN_X + ACIL_CN_X       & 
                                                       + ACLL_AN_X + ACIL_AN_X       & 
                                                       + ACLL_LS_X + ACIL_LS_X ) * DM , 3 ) + ER_PRC2
 
      if (associated(PREVTOT   ))   PREVTOT    =  SUM( ( REV_CN_X  + RSU_CN_X        & 
                                                       + REV_AN_X  + RSU_AN_X        & 
                                                       + REV_LS_X  + RSU_LS_X ) * DM , 3 )

      if(associated(cnv_mfc)) then
         if(associated(cnv_freq)) then
            cnv_freq(:,:) = 0.0
            do j=1,jm
               do i=1,im
                  do l=lm,1,-1
                     if(cnv_mfc(i,j,l)/=0.0) then
                        cnv_freq(i,j) = 1.0
                        exit
                     endif
                  enddo
               enddo
            enddo
         endif

         if(associated(cnv_basep)) then
            cnv_basep(:,:) = MAPL_UNDEF
            do j=1,jm
               do i=1,im
                  do l=lm,1,-1
                     if(cnv_mfc(i,j,l)/=0.0) then
                        cnv_basep(i,j) = ple(i,j,l-1)
                        exit
                     endif
                  enddo
               enddo
            enddo
         endif

         if(associated(cnv_topp)) then
             cnv_topp(:,:) = MAPL_UNDEF
             do j=1,jm
               do i=1,im
                  do l=1,lm
                     if(cnv_mfc(i,j,l)/=0.0) then
                        cnv_topp(i,j) = ple(i,j,l)
                        exit
                     endif
                  enddo
               enddo
            enddo
         endif
      endif




      if (associated(LS_ARF ))   LS_ARF  = LS_ARFX
      if (associated(AN_ARF ))   AN_ARF  = AN_ARFX
      if (associated(CN_ARF ))   CN_ARF  = CN_ARFX
      if (associated(XQLLS  ))   XQLLS   = QLLS
      if (associated(XQILS  ))   XQILS   = QILS
      if (associated(XCLLS  ))   XCLLS   = CLLS
      if (associated(XQLCN  ))   XQLCN   = QLCN
      if (associated(XQICN  ))   XQICN   = QICN
      if (associated(XCLCN  ))   XCLCN   = CLCN
      if (associated(QITOT  ))   QITOT   = QICN + QILS
      if (associated(QLTOT  ))   QLTOT   = QLCN + QLLS
      if (associated(QCTOT  ))   QCTOT   = QLCN + QLLS + QICN + QILS
      if (associated(TVQ1   ))   TVQ1    = SUM( ( Q1 +  QLLS + QLCN + QILS + QICN )*DM , 3 ) & 
                                                +  TPREC*DT_MOIST
      if (associated(TVE1   ))   TVE1    = SUM( (  MAPL_CP*TEMP + MAPL_ALHL*Q1             & 
                                                -  MAPL_ALHF*(QILS+QICN) )*DM , 3 )        &
                                                -  MAPL_ALHF*( CN_SNR  + LS_SNR  + AN_SNR )*DT_MOIST

      if (associated(DCPTE  ))   DCPTE   = (  SUM(  MAPL_CP*TEMP *DM , 3) - DCPTE )/DT_MOIST 
      if (associated(CWP    ))   CWP     = SUM( ( QLCN+QLLS+QICN+QILS )*DM , 3 )
      if (associated(LWP    ))   LWP     = SUM( ( QLCN+QLLS ) *DM , 3 )
      if (associated(IWP    ))   IWP     = SUM( ( QICN+QILS ) *DM , 3 )
      if (associated(CCWP   ))   CCWP    = SUM(   CNV_QC *DM , 3 )
      if (associated(TPW    ))   TPW     = SUM(   Q1         *DM , 3 )
      if (associated(RH2    ))   RH2     = max(MIN( Q1/GEOS_QSAT (TH1*PK, PLO) , 1.02 ),0.0)
      if (associated(PRECU  ))   PRECU   = CN_PRC2
      if (associated(PRELS  ))   PRELS   = LS_PRC2 + AN_PRC2
      if (associated(SNR    ))   SNR     = LS_SNR  + AN_SNR + CN_SNR
      if (associated(TT_PRCP))   TT_PRCP = TPREC

      if (associated(HOURNORAIN))   then
        call ESMF_ClockGet(CLOCK, currTime=CurrentTime, rc=STATUS)
        call ESMF_TimeGet (CurrentTime, YY=YEAR, MM=MONTH, DD=DAY, H=HR, M=MN, S=SE, RC=STATUS )
        if(SE==0 .and. MN==0) then
          HOURNORAIN  = 0.
        else
          where(TPREC.EQ.0.)
           HOURNORAIN  = HOURNORAIN + DT_MOIST
          endwhere
        endif
      endif

      if (associated(CN_PRCP))   CN_PRCP = CN_PRC2 + CN_SNR
      if (associated(AN_PRCP))   AN_PRCP = AN_PRC2 + AN_SNR
      if (associated(LS_PRCP))   LS_PRCP = LS_PRC2 + LS_SNR
      if (associated(ER_PRCP))   ER_PRCP = ER_PRC2 
      if (associated(FILLNQV))   FILLNQV = FILLQ / DT_MOIST 
      if (associated(DQDT   ))   DQDT    = (Q1  - Q )/DT_MOIST
      !!if(associated(QSATI   ))   DQS      = GEOS_DQSAT(TEMP, PLO, qsat=QSATi, OVER_ICE=.TRUE. )
      !!if(associated(QSATl   ))   DQS      = GEOS_DQSAT(TEMP, PLO, qsat=QSATl, OVER_LIQUID=.TRUE. )
      if (associated(TI     ))   TI      = (TH1 - TH)*(PLE(:,:,1:LM)-PLE(:,:,0:LM-1))/DT_MOIST
      if (associated(UI     ))   UI      = (U1  - U )/DT_MOIST
      if (associated(VI     ))   VI      = (V1  - V )/DT_MOIST
      if(associated(DQLDT   ))   then 
                    DQLDT   = (QLLS + QLCN - DQLDT)
                    where( ABS(DQLDT) < 1.e-20 )
                         DQLDT=0.0
                    endwhere
                    DQLDT   = DQLDT / DT_MOIST
      endif
      if(associated(DQIDT   ))   then 
                    DQIDT   = (QILS + QICN - DQIDT)
                    where( ABS(DQIDT) < 1.e-20 )
                         DQIDT=0.0
                    endwhere
                    DQIDT   = DQIDT / DT_MOIST
      endif

      temp = (0.5/DT_MOIST)*((V1**2+U1**2)  - (V**2+U**2))*DM
      !temp = (1.0/DT_MOIST)*( U *( U1 - U )  + V * ( V1 - V ) ) *DM
     
      IKEX  = SUM( TEMP , 3 )     
      IKEX2 = MAX(  SUM( KEX * DM , 3 ) ,  1.0e-6 ) ! floor at 1e-6 W m-2   

      if (associated(KEMST   ))  KEMST    = IKEX
      if (associated(KEMST2  ))  KEMST2   = IKEX2

                                   !scaled 3D kinetic energy dissipation
      if (associated(KEDISS  ))  then 
             do L=1,LM
               KEDISS(:,:,L)   = (IKEX/IKEX2) * KEX(:,:,L)
             enddo
      end if


      if(associated(DTDTFRIC)) then 
         do L=1,LM
            DTDTFRIC(:,:,l) = -(1./MAPL_CP)*(IKEX/IKEX2) * KEX(:,:,L) * (PLE(:,:,L)-PLE(:,:,L-1))
         end do
      end if

! Compute aerosol tendncies due to moist
! --------------------------------------
 ! if ( associated(ddudt) ) then
 !    do k = 1, km
 !       if it is dust: (based on qnames)
 !          ddudt = sum_k q(k) * delp(k) - ddudt
 !    else if is sea salt ...
!
 ! end do


      if (associated(CNVRNZ     ))   CNVRNZ    =  SUM( DQRC  * DM , 3 )
      if (associated(PDFLZ      ))   PDFLZ     =  SUM( DLPDF_X * DM , 3 )
      if (associated(CNVLZ      ))   CNVLZ     =  SUM( DCNVL_X * DM , 3 )

      if (associated(PDFIZ      ))   PDFIZ     =  SUM( DIPDF_X * DM , 3 )
      if (associated(CNVIZ      ))   CNVIZ     =  SUM( DCNVI_X * DM , 3 )

      if (associated(EVPCZ      ))   EVPCZ     =  SUM( ( EVAPC_X  + DLFIX_X ) * DM , 3 )
      if (associated(EVPPZ      ))   EVPPZ     =  SUM( ( REV_CN_X + REV_AN_X + REV_LS_X ) * DM , 3 )

      if (associated(SUBCZ      ))   SUBCZ     =  SUM( ( SUBLC_X  + DIFIX_X ) * DM , 3 )
      if (associated(SUBPZ      ))   SUBPZ     =  SUM( ( RSU_CN_X + RSU_AN_X + RSU_LS_X ) * DM , 3 )

      if (associated(FRZCZ      ))   FRZCZ     =  SUM( FRZ_TT_X * DM , 3 )
      if (associated(FRZPZ      ))   FRZPZ     =  SUM( FRZ_PP_X * DM , 3 )
      if (associated(COLLIZ     ))   COLLIZ    =  SUM( ( ACIL_CN_X + ACIL_AN_X + ACIL_LS_X ) * DM , 3 )

      if (associated(COLLLZ     ))   COLLLZ    =  SUM( ( ACLL_CN_X + ACLL_AN_X + ACLL_LS_X ) * DM , 3 )
      if (associated(AUTZ       ))   AUTZ      =  SUM( AUT_X * DM , 3 )

      if (associated(SDMZ       ))   SDMZ      =  SUM( SDM_X * DM , 3 )


! Replace the modified humidity
!------------------------------

      Q = Q1

      call MAPL_TimerOn (STATE,"--FLASH",RC=STATUS)
      VERIFY_(STATUS)

! Parameterized lightning flash rates [km^{-2} s^{-1}]
!-----------------------------------------------------

      CALL MAPL_GetPointer(EXPORT,  LFR,  'LFR', ALLOC=.TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(EXPORT, A1X1, 'A1X1', ALLOC=.TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(EXPORT, A2X2, 'A2X2', ALLOC=.TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(EXPORT, A3X3, 'A3X3', ALLOC=.TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(EXPORT, A4X4, 'A4X4', ALLOC=.TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(EXPORT, A5X5, 'A5X5', ALLOC=.TRUE., RC=STATUS)
      VERIFY_(STATUS)

      CALL flash_rate(STATE,    &
                      IM*JM,    &
                      LM,       &
                      TS,       &
                      CNV_TOPP, &
                      FROCEAN,  &
                      CN_PRCP,  &
                      CAPE,     &
                      CNV_MFC,  &
                      TH,       &
                      PLE,      &
                      ZLE,      &
                      LFR,      &
                      A1X1,     &
                      A2X2,     &
                      A3X3,     &
                      A4X4,     &
                      A5X5,     &
                      RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_TimerOff(STATE,"--FLASH",RC=STATUS)
      VERIFY_(STATUS)

! Deallocate temp space if necessary
!-----------------------------------

      if(associated(DDF_BYNC  )) deallocate( DDF_BYNCz )
      if(associated(DDF_MUPH  )) deallocate( DDF_MUPHz )
      if(associated(DDF_RH1   )) deallocate( DDF_RH1z )
      if(associated(DDF_RH2   )) deallocate( DDF_RH2z )
      if(associated(DDF_DQDT  )) deallocate( DDF_DQDTz )
      if(associated(DDF_DTDT  )) deallocate( DDF_DTDTz )
      if(associated(DDF_TC    )) deallocate( DDF_TCz  )
      if(associated(DDF_QVC   )) deallocate( DDF_QVCz )
      if(associated(DDF_MFC   )) deallocate( DDF_MFCz )

      if(ALLOC_BYNCY    ) deallocate(BYNCY    )
      if(ALLOC_CAPE     ) deallocate(CAPE     )
      if(ALLOC_INHB     ) deallocate(INHB     )
      if(ALLOC_CNV_DQLDT) deallocate(CNV_DQLDT)
      if(ALLOC_CNV_MF0  ) deallocate(CNV_MF0  )
      if(ALLOC_CNV_MFD  ) deallocate(CNV_MFD  )
      if(ALLOC_CNV_MFC  ) deallocate(CNV_MFC  )
      if(ALLOC_CNV_TOPP ) deallocate(CNV_TOPP )
      if(ALLOC_CNV_UPDF ) deallocate(CNV_UPDF )
      if(ALLOC_CNV_CVW  ) deallocate(CNV_CVW  )
      if(ALLOC_CNV_QC )   deallocate(CNV_QC   )
      if(ALLOC_RAD_CF   ) deallocate(RAD_CF   )
      if(ALLOC_RAD_QV   ) deallocate(RAD_QV   )
      if(ALLOC_RAD_QL   ) deallocate(RAD_QL   )
      if(ALLOC_RAD_QI   ) deallocate(RAD_QI   )
      if(ALLOC_RAD_QR   ) deallocate(RAD_QR   )
      if(ALLOC_RAD_QS   ) deallocate(RAD_QS   )

      if(ALLOC_CLDREFFL ) deallocate(CLDREFFL )
      if(ALLOC_CLDREFFI ) deallocate(CLDREFFI )
      if(ALLOC_CLDREFFR ) deallocate(CLDREFFR )
      if(ALLOC_CLDREFFS ) deallocate(CLDREFFS )
      if(ALLOC_CLDNCCN  ) deallocate(CLDNCCN  )

      if(ALLOC_DQRC )      deallocate ( DQRC  )

      if(ALLOC_ENTLAM )    deallocate ( ENTLAM )

      deallocate ( QRAIN )
      deallocate ( QSNOW )

      call MAPL_TimerOff(STATE,"-MISC")

!  All done
!-----------

      RETURN_(ESMF_SUCCESS)

    end subroutine MOIST_DRIVER

!!!!!!!!-!-!-!!!!!!!!
  end subroutine RUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FILLQ2ZERO( Q, DM, FILLQ  )

  real, dimension(:,:,:),   intent(inout)  :: Q
  real, dimension(:,:,:),   intent(in)     :: DM
  real, dimension(:,:),     intent(  out)  :: FILLQ
  integer                                  :: IM,JM,LM
  integer                                  :: I,J,K,L

  real                                     :: TPW, NEGTPW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fills in negative q values in a mass conserving way.
! Conservation of TPW was checked.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  IM = SIZE( Q, 1 )
  JM = SIZE( Q, 2 )
  LM = SIZE( Q, 3 )


  do j=1,JM
  do i=1,IM

     TPW = SUM( Q(i,j,:)*DM(i,j,:) )

     NEGTPW = 0.
     do l=1,LM
        if ( Q(i,j,l) < 0.0 ) then 
             NEGTPW   = NEGTPW + ( Q(i,j,l)*DM( i,j,l ) )
             Q(i,j,l) = 0.0
        endif
     enddo

     do l=1,LM
        if ( Q(i,j,l) >= 0.0 ) then 
             Q(i,j,l) = Q(i,j,l)*( 1.0+NEGTPW/(TPW-NEGTPW) )
        endif
     enddo

     FILLQ(i,j) = -NEGTPW

   end do
   end do


  

  end subroutine FILLQ2ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine BUOYANCY( T, TH, Q, QS, DQS, DM, ZLO, BUOY, CAPE, INHB, IM, JM, LM )


! !DESCRIPTION: Computes the buoyancy $ g \frac{T_c-T_e}{T_e} $ at each level
!  for a parcel raised from the surface. $T_c$ is the virtual temperature of
!  the parcel and $T_e$ is the virtual temperature of the environment.

  integer,                     intent(in)  :: IM,JM,LM
  !real, dimension(IM,JM,LM),   intent(in)  :: T, Q, QS, DQS, TH, DM, ZLO
  !real, dimension(IM,JM,LM),   intent(out) :: BUOY
  !real, dimension(IM,JM),      intent(out) :: CAPE, INHB

  real, dimension(:,:,:),   intent(in)  :: T, Q, QS, DQS, TH, DM, ZLO
  real, dimension(:,:,:),   intent(out) :: BUOY
  real, dimension(:,:),      intent(out) :: CAPE, INHB

  integer                         :: L
  real,    dimension(IM,JM     )  :: HC
  logical, dimension(IM,JM     )  :: UNSTABLE

       HC  =  T(:,:,LM) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,LM) + (MAPL_ALHL/MAPL_CP)*Q(:,:,LM)

       do L=LM-1,1,-1
          BUOY(:,:,L) = HC - (T(:,:,L) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,L) + (MAPL_ALHL/MAPL_CP)*QS(:,:,L))
          BUOY(:,:,L) = BUOY(:,:,L) / ( (1.+ (MAPL_ALHL/MAPL_CP)*DQS(:,:,L))*T(:,:,L) )
       enddo

       BUOY(:,:,LM) = 0.0

       UNSTABLE = .false.

       CAPE = 0.
       INHB = 0.

       do L=1,LM-1
          where(BUOY(:,:,L)>0.) UNSTABLE=.true.
          where(UNSTABLE) 
             CAPE = CAPE + BUOY(:,:,L)*DM(:,:,L)
          end where
          where(UNSTABLE.AND.BUOY(:,:,L) < 0.) 
             INHB = INHB - BUOY(:,:,L)*DM(:,:,L)
          end where
       end do

       UNSTABLE = CAPE > 0.0

       where(.not.UNSTABLE) 
          CAPE=MAPL_UNDEF
          INHB=MAPL_UNDEF
       elsewhere
          !!INHB=INHB/(CAPE-INHB)
       end where

     end subroutine BUOYANCY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function FINDLFC( BUOY, IM, JM, LM ) result( KLFC )
! !DESCRIPTION: 
     integer,                   intent(in)  :: IM, JM, LM
     real, dimension(IM,JM,LM), intent(in)  :: BUOY

     integer, dimension(IM,JM)              :: KLFC

     integer                                :: I, J, L

     do I = 1, IM
        do J = 1, JM

           KLFC(I,J) = 0    
           do L = LM,1,-1
              IF( BUOY(I,J,L) > 0. ) THEN 
                  KLFC(I,J) = L
                  EXIT
              ENDIF
           enddo

        end do
     end do
  end function FINDLFC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function FINDPBL( KH, IM, JM, LM ) result( KPBL )
! !DESCRIPTION: 
     integer                    , intent(in) :: IM,JM,LM
     real, dimension(IM,JM,0:LM), intent(in) :: KH

     integer, dimension(IM,JM)               :: KPBL

     integer                                 :: I, J, L
     real                                    :: KHCRIT

     KHCRIT = 2.0  ! m+2 s-1

     do I = 1, IM
        do J = 1, JM

           KPBL(I,J) = LM    
           do L = LM-1,1,-1
              IF( ( KH(I,J,L) >= KHCRIT ).AND.( KH(I,J,L-1) < KHCRIT ) ) THEN ! "top" is between L and L-1
                  KPBL(I,J) = L+1   ! returned index for CBL q,t etc. is just below PBL top
                  EXIT
              ENDIF
           enddo

           KPBL(I,J)=MIN( LM-1, KPBL(I,J) )      

        end do
     end do


  end function FINDPBL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function FINDLCL( THM, QM, PL, PK, IM, JM, LM ) result( KLCL )
! !DESCRIPTION: 
     integer,                      intent(in) :: IM, JM, LM
     real,    dimension(IM,JM,LM), intent(in) :: THM, QM
     real,    dimension(IM,JM,LM), intent(in) :: PL, PK

     integer, dimension(IM,JM)             :: KLCL

     real, dimension(LM) :: TPCL, QSPCL
     integer             :: I, J, L, KOFFSET

     do I = 1, IM
        do J = 1, JM

           TPCL  = THM(I,J,LM) * PK(I,J,:)
           QSPCL = GEOS_QSAT(TPCL, PL(I,J,:) )

           KLCL(I,J) = 0

           do L = LM,1,-1
              if( QM(I,J,LM) >= QSPCL(L) ) then
                  KLCL(I,J) = L
                  exit
              endif
           enddo

          
           !! ------------------------------------
           !!   Disabled for Daedalus (e0203) merge
           !! ------------------------------------
!!AMM      KOFFSET   = INT ( (LM - KLCL(I,J))/2 )   !! disable for Gan4_0
           KOFFSET   = 0
           KLCL(I,J) = MIN ( LM-1,  KLCL(I,J)+KOFFSET )

        end do
     end do

  end function FINDLCL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 SUBROUTINE flash_rate( STATE, nc, lm, TS, CCTP, FROCEAN, CN_PRCP, &
                        CAPE, CNV_MFC, TH, PLE, ZLE, strokeRate, &
                        A1X1, A2X2, A3X3, A4X4, A5X5, RC)

!=====================================================================================
!BOP
! !DESCRIPTION:
!  Generate lightning flash rates [km$^{-2}$ s$^{-1}$] using a six-variable polynomial fit.\\
!
!
!  ORIGIN AND CONTACT\\
!  Dr. Dale Allen, Associate Research Scientist\\
!  Dept. of Atmospheric and Oceanic Science\\
!  University of Maryland\\
!  College Park, MD 20742\\
!  301-405-7629 (ph); 301-314-9482 (fax)\\
!  http://www.meto.umd.edu/~allen\\
!  
!
!  FORMULATION NOTES\\
!  Predictor variables are set to zero where CN\_PRCP is zero or where the 
!   optical depth cloud top height is less than 5.5 km.
!  The fit returns flash rates in units km$^{-2}$ day$^{-1}$.  Convert to 
!   km$^{-2}$ s$^{-1}$ for the export state.\\
!
!
!  OTHER NOTES OF INTEREST\\
!  MOIST sets CNV\_TOPP to zero if there is an absence of convection.
!EOP
! !REVISION HISTORY
! 30 Nov 2011 Nielsen     First crack
! 29 Feb 2012 Nielsen     Accomodate CNV\_TOPP MAPL\_UNDEF for and after Fortuna-2\_5\_p4
!=====================================================================================

  TYPE(MAPL_MetaComp), POINTER :: STATE ! Internal MAPL_Generic state

  INTEGER, INTENT(IN) :: nc     ! Number of cells
  INTEGER, INTENT(IN) :: lm     ! Number of layers

  REAL, INTENT(IN), DIMENSION(nc) :: TS       ! Surface temperature [K]
  REAL, INTENT(IN), DIMENSION(nc) :: CCTP     ! Convective cloud top pressure [Pa] with MAPL_UNDEFs
  REAL, INTENT(IN), DIMENSION(nc) :: FROCEAN  ! Areal ocean fraction
  REAL, INTENT(IN), DIMENSION(nc) :: CN_PRCP  ! Convective precipitation [kg m^{-2} s^{-1}]
  REAL, INTENT(IN), DIMENSION(nc) :: CAPE     ! Convective available potential energy [J m^{-2}]

  REAL, INTENT(IN), DIMENSION(nc,lm) :: TH        ! Potential temperature [K]
  REAL, INTENT(IN), DIMENSION(nc,0:lm) :: CNV_MFC ! Convective mass flux [kg m^{-2} s^{-1}]
  REAL, INTENT(IN), DIMENSION(nc,0:lm) :: PLE     ! Layer interface pressures  [Pa]
  REAL, INTENT(IN), DIMENSION(nc,0:lm) :: ZLE     ! Layer depths [m]

  REAL, INTENT(OUT), DIMENSION(nc) :: strokeRate  ! Flashes per second
  REAL, INTENT(OUT), DIMENSION(nc) :: A1X1
  REAL, INTENT(OUT), DIMENSION(nc) :: A2X2
  REAL, INTENT(OUT), DIMENSION(nc) :: A3X3
  REAL, INTENT(OUT), DIMENSION(nc) :: A4X4
  REAL, INTENT(OUT), DIMENSION(nc) :: A5X5

! Error log variables
! -------------------
  INTEGER :: STATUS
  INTEGER, OPTIONAL, INTENT(OUT) :: RC
  CHARACTER(LEN=ESMF_MAXSTR) :: IAm

! Local variables
! ---------------
  INTEGER :: i            ! General-purpose integers
  INTEGER :: k
  INTEGER :: n

  REAL :: a0c,a0m         ! Coefficients at continental and marine locations
  REAL :: a1c,a1m
  REAL :: a2c,a2m
  REAL :: a3c,a3m
  REAL :: a4c,a4m
  REAL :: a5c,a5m

  REAL :: x1Divisor       ! Divisors for x1-x5.
  REAL :: x2Divisor
  REAL :: x3Divisor
  REAL :: x4Divisor
  REAL :: x5Divisor

  REAL :: x5Power         ! Exponent for the surface temperature deviation predictor

  REAL :: sfcTLimit       ! Temperature thresholds
  REAL :: airTLimit

  REAL :: hPaCldTop       ! Cloud top limiter for weak/no convection

  REAL, ALLOCATABLE, DIMENSION(:) :: x1         ! Five independent variables
  REAL, ALLOCATABLE, DIMENSION(:) :: x2
  REAL, ALLOCATABLE, DIMENSION(:) :: x3
  REAL, ALLOCATABLE, DIMENSION(:) :: x4
  REAL, ALLOCATABLE, DIMENSION(:) :: x5

  REAL, ALLOCATABLE, DIMENSION(:) :: cloudTopAG ! Cloud top height above ground
  REAL, ALLOCATABLE, DIMENSION(:) :: cnv_topp   ! Convective cloud top pressure with MAPL_UNDEFs
                                                ! changed to zero

  REAL, ALLOCATABLE, DIMENSION(:,:) :: dZ       ! Layer depths [m]
  REAL, ALLOCATABLE, DIMENSION(:,:) :: p        ! Pressure at middle of layer [Pa]
  REAL, ALLOCATABLE, DIMENSION(:,:) :: T        ! Air temperature at middle of layer [K]

  INTEGER, ALLOCATABLE, DIMENSION(:)   :: weakCnvMask   ! Weak or no convection mask
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask          ! Working mask
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: cloudTopMask  ! Mask is 1 below cloud top

! Preliminaries
! -------------
  RC = 0
  strokeRate = 0
  Iam = "flash_rate"

! Coefficients of the predictors, marine locations
! ------------------------------------------------
  CALL MAPL_GetResource(STATE,a0m,'MARINE_A0:',DEFAULT= 0.0139868,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,a1m,'MARINE_A1:',DEFAULT= 0.0358764,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,a2m,'MARINE_A2:',DEFAULT=-0.0610214,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,a3m,'MARINE_A3:',DEFAULT=-0.0102320,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,a4m,'MARINE_A4:',DEFAULT= 0.0031352,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,a5m,'MARINE_A5:',DEFAULT= 0.0346241,RC=STATUS)
  VERIFY_(STATUS)

! Coefficients of the predictors, continental locations
! -----------------------------------------------------
  CALL MAPL_GetResource(STATE,a0c,'CONTINENT_A0:',DEFAULT=-0.0183172,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,a1c,'CONTINENT_A1:',DEFAULT=-0.0562338,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,a2c,'CONTINENT_A2:',DEFAULT= 0.1862740,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,a3c,'CONTINENT_A3:',DEFAULT=-0.0023363,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,a4c,'CONTINENT_A4:',DEFAULT=-0.0013838,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,a5c,'CONTINENT_A5:',DEFAULT= 0.0114759,RC=STATUS)
  VERIFY_(STATUS)

! Divisors for nondimensionalization of the predictors
! ----------------------------------------------------
  CALL MAPL_GetResource(STATE,x1Divisor,'X1_DIVISOR:',DEFAULT=4.36,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,x2Divisor,'X2_DIVISOR:',DEFAULT=9.27,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,x3Divisor,'X3_DIVISOR:',DEFAULT=34.4,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,x4Divisor,'X4_DIVISOR:',DEFAULT=21.4,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,x5Divisor,'X5_DIVISOR:',DEFAULT=14600.,RC=STATUS)
  VERIFY_(STATUS)

! Exponent for the surface temperature deviation predictor
! --------------------------------------------------------
  CALL MAPL_GetResource(STATE,x5Power,'X5_EXPONENT:',DEFAULT=3.00,RC=STATUS)
  VERIFY_(STATUS)

! Threshold temperatures
! ----------------------
  CALL MAPL_GetResource(STATE,sfcTLimit,'SFC_T_LIMIT:',DEFAULT=273.0,RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_GetResource(STATE,airTLimit,'AIR_T_LIMIT:',DEFAULT=263.0,RC=STATUS)
  VERIFY_(STATUS)

! Cloud-top pressure limiter
! --------------------------
  CALL MAPL_GetResource(STATE,hPaCldTop,'CLOUD_TOP_LIMIT:',DEFAULT=500.,RC=STATUS)
  VERIFY_(STATUS)

! Layer depths [m]
! ----------------
  ALLOCATE(dZ(nc,lm),STAT=STATUS)
  VERIFY_(STATUS)
  dZ = zle(:,0:lm-1)-zle(:,1:lm)

! Pressure at mid-layer [Pa]
! --------------------------
  ALLOCATE(p(nc,lm),STAT=STATUS)
  VERIFY_(STATUS)
  p = (ple(:,1:lm)+ple(:,0:lm-1))*0.50

! Temperature at mid-layer [K]
! ----------------------------
  ALLOCATE(T(nc,lm),STAT=STATUS)
  VERIFY_(STATUS)
  T = TH*((p*1.00E-05)**(MAPL_RGAS/MAPL_CP))

! Reset CNV_TOPP's MAPL_UNDEFs to zeroes
! --------------------------------------
  ALLOCATE(cnv_topp(nc),STAT=STATUS)
  WHERE(CCTP == MAPL_UNDEF)
   cnv_topp = 0.00
  ELSEWHERE
   cnv_topp = CCTP
  END WHERE

! Set weak/no convection mask
! ---------------------------
  ALLOCATE(weakCnvMask(nc),STAT=STATUS)
  VERIFY_(STATUS)
  weakCnvMask = 0
  WHERE(cn_prcp == 0.00 .OR. cnv_topp >= hPaCldTop*100.00 .OR. cape >= MAPL_UNDEF) weakCnvMask = 1

! Convective cloud top mask
! -------------------------
  ALLOCATE(cloudTopMask(nc,lm),STAT=STATUS)
  VERIFY_(STATUS)
  cloudTopMask = 0
  DO k = 1,lm
   WHERE(ple(1:nc,k) > cnv_topp(1:nc) .AND. cnv_topp(1:nc) > 0.00) cloudTopMask(1:nc,k) = 1
  END DO
  
! Cloud top distance above ground [m]
! -----------------------------------
  ALLOCATE(cloudTopAG(nc),STAT=STATUS)
  VERIFY_(STATUS)
  cloudTopAG = 0.00
  DO i = 1,nc
   n = SUM(cloudTopMask(i,1:lm))
   IF(n > 0) cloudTopAG(i) = SUM(dZ(i,lm-n+1:lm))
  END DO

! X1: Cold cloud depth: Vertical extent [km] where T < airTLimit and p > cnv_topp
! -------------------------------------------------------------------------------
  ALLOCATE(x1(nc),STAT=STATUS)
  VERIFY_(STATUS)
  ALLOCATE(mask(nc,lm),STAT=STATUS)
  VERIFY_(STATUS)

  mask = 0
  WHERE(T < airTLimit .AND. cloudTopMask == 1) mask = 1

  x1 = 0.00
  DO i = 1,nc
   DO k = 1,lm
    IF(mask(i,k) == 1) x1(i) = x1(i)+dZ(i,k)*0.001
   END DO
  END DO
  WHERE(weakCnvMask == 1) x1 = 0.00
  x1 = x1/x1Divisor

! X4: Integrated convective mass flux
! -----------------------------------
  ALLOCATE(x4(nc),STAT=STATUS)
  VERIFY_(STATUS)
  x4 = 0.00
  DO i = 1,nc
   DO k = 1,lm
    IF(mask(i,k) == 1) x4(i) = x4(i)+cnv_mfc(i,k)*dZ(i,k)
   END DO
  END DO
  WHERE(weakCnvMask == 1) x4 = 0.00
  x4 = x4/x4Divisor

! X5: Surface temperature deviation from sfcTLimit, positive only.
! Note: UNDEF TS test retains the ability to boot-strap moist_import_rst.
! -----------------------------------------------------------------------
  ALLOCATE(x5(nc),STAT=STATUS)
  VERIFY_(STATUS)
  WHERE(TS == MAPL_UNDEF)
   x5 = 0.00
  ELSEWHERE
   x5 = TS-sfcTLimit
  END WHERE
  WHERE(weakCnvMask == 1) x5 = 0.00
  WHERE(x5 < 0.00) x5 = 0.00
  x5 = x5**x5Power/x5Divisor

! X2: Total cloud depth [km]
! --------------------------
  ALLOCATE(x2(nc),STAT=STATUS)
  VERIFY_(STATUS)
  x2 = cloudTopAG*0.001
  WHERE(weakCnvMask == 1) x2 = 0.00
  x2 = x2/x2Divisor

! X3: CAPE
! --------
  ALLOCATE(x3(nc),STAT=STATUS)
  VERIFY_(STATUS)
  x3 = cape
  WHERE(weakCnvMask == 1) x3 = 0.00
  x3 = x3/x3Divisor

! Polynomial fit [units: km^{-2} s^{-1}] and individual
! terms including marine and continental discrimination
! -----------------------------------------------------
  WHERE(frOcean >= 0.01)
   strokeRate = (a0m + a1m*x1 + a2m*x2 + a3m*x3 + a4m*x4 + a5m*x5)/86400.00
   A1X1 = a1m*x1/86400.00
   A2X2 = a2m*x2/86400.00
   A3X3 = a3m*x3/86400.00
   A4X4 = a4m*x4/86400.00
   A5X5 = a5m*x5/86400.00
  ELSEWHERE
   strokeRate = (a0c + a1c*x1 + a2c*x2 + a3c*x3 + a4c*x4 + a5c*x5)/86400.00
   A1X1 = a1c*x1/86400.00
   A2X2 = a2c*x2/86400.00
   A3X3 = a3c*x3/86400.00
   A4X4 = a4c*x4/86400.00
   A5X5 = a5c*x5/86400.00
  END WHERE

! Eliminate negatives
! -------------------
  WHERE(strokeRate < 0.00) strokeRate = 0.00

! Set rate to zero where any of x1 through x5 are zero
! ----------------------------------------------------
  WHERE(x1 == 0.00) strokeRate = 0.00
  WHERE(x2 == 0.00) strokeRate = 0.00
  WHERE(x3 == 0.00) strokeRate = 0.00
  WHERE(x4 == 0.00) strokeRate = 0.00
  WHERE(x5 == 0.00) strokeRate = 0.00

! Clean up
! --------
  DEALLOCATE(x1,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(x2,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(x3,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(x4,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(x5,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(cnv_topp,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(dZ,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(p,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(T,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(cloudTopAG,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(mask,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(cloudTopMask,STAT=STATUS)
  VERIFY_(STATUS)
  DEALLOCATE(weakCnvMask,STAT=STATUS)
  VERIFY_(STATUS)

 END SUBROUTINE flash_rate

end module GEOS_MoistGridCompMod

