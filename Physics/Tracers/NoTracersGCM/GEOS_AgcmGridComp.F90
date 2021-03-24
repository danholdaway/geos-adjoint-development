! $Id: GEOS_AgcmGridComp.F90,v 1.77 2014-01-24 21:06:29 ltakacs Exp $

#include "MAPL_Generic.h"

!#define PRINT_STATES
#define FULLPHYSICS
#define GCM
#define debug 0

! Held-Suarez is not in module GEOSGCM?????

#if defined HS
#define  GEOS_physicsGridCompMod GEOS_hsGridCompMod
#undef   FULLPHYSICS
#endif

#if defined SCM
#define GEOS_superdynGridCompMod GEOS_singcolGridCompMod
#undef  GCM
#endif

!=============================================================================
!BOP

! !MODULE: GEOS_AgcmGridCompMod -- A Module to combine Supedynamics and Physics Gridded Components

! !INTERFACE:

module GEOS_AgcmGridCompMod

! !USES:

  use ESMF
  use MAPL_Mod
  use GEOS_TopoGetMod

  use GEOS_superdynGridCompMod,  only:  SDYN_SetServices => SetServices
  use GEOS_physicsGridCompMod,   only:  PHYS_SetServices => SetServices
  use MAPL_OrbGridCompMod,       only:  ORB_SetServices => SetServices

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION: This gridded component (GC) combines the Superdynamics GC, 
!   and Physics GC into a new composite Agcm GC.

!\begin{verbatim}
!       DUDT .... Mass-Weighted U-Wind      Tendency (Pa m /s)
!       DVDT .... Mass-Weighted V-Wind      Tendency (Pa m /s)
!       DPEDT ... Edge-Pressure             Tendency (Pa   /s)
!       DTDT .... Mass-Weighted Temperature Tendency (Pa K /s)
!       TRACER .. Friendly Tracers                   (unknown)
!\end{verbatim}
 
!EOP
  
  integer :: SDYN
  integer :: PHYS
  integer :: ORB

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  The SetServices for the Physics GC needs to register its
!   Initialize and Run.  It uses the MAPL\_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs (SURF, CHEM, RADIATION, MOIST, TURBULENCE) and runs their
!   respective SetServices.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals

    integer                       :: I
    logical                       :: ANA_TS
    type (ESMF_Config)            :: CF
    character(len=ESMF_MAXSTR)    :: ReplayMode

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run       , RC=STATUS )
    VERIFY_(STATUS)

! Get the configuration from the component
!-----------------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

    call ESMF_ConfigGetAttribute ( CF, I, Label="ANALYZE_TS:", default=0, RC=STATUS)
    VERIFY_(STATUS)
    ANA_TS = I /= 0

    call ESMF_ConfigGetAttribute(CF, ReplayMode, Label='REPLAY_MODE:', default="NoReplay", RC=STATUS )
    VERIFY_(STATUS)
    if(ANA_TS) then
       ASSERT_(ReplayMode=="NoReplay")
    endif
 
!BOS

! !IMPORT STATE:

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_analysis_increment',          &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_analysis_increment',         &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'temperature_analysis_increment',            &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_analysis_increment',          &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DQVDT',                                     &
         LONG_NAME  = 'specific_humidity_analysis_increment',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DO3DT',                                     &
         LONG_NAME  = 'ozone_analysis_increment',                  &
         UNITS      = 'mol mol-1',                                 &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temperature_increment',                &
         UNITS      = 'K',                                         &
         RESTART    = ANA_TS,                                      &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

! !INTERNAL STATE:

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_bias_tendency',               &
         UNITS      = 'm s-2',                                     &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_bias_tendency',              &
         UNITS      = 'm s-2',                                     &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'temperature_bias_tendency',                 &
         UNITS      = 'K s-1',                                     &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_bias_tendency',               &
         UNITS      = 'Pa s-1',                                    &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DQVDT',                                     &
         LONG_NAME  = 'specific_humidity_bias_tendency',           &
         UNITS      = 'kg kg-1 s-1',                               &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DO3DT',                                     &
         LONG_NAME  = 'ozone_bias_tendency',                       &
         UNITS      = 'mol mol-1 s-1',                             &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temperature_tendency',                 &
         UNITS      = 'K s-1',                                     &
         FRIENDLYTO = trim(COMP_NAME),                             &
         RESTART    = ANA_TS,                                      &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

! !EXPORT STATE:


    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DUDT_ANA',                                  &
         LONG_NAME  = 'total_eastward_wind_analysis_tendency',     &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DVDT_ANA',                                  &
         LONG_NAME  = 'total_northward_wind_analysis_tendency',    &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DTDT_ANA',                                  &
         LONG_NAME  = 'total_temperature_analysis_tendency',       &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DPEDT_ANA',                                 &
         LONG_NAME  = 'total_edge_pressure_analysis_tendency',    &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DQVDT_ANA',                                 &
         LONG_NAME  = 'total_specific_humidity_analysis_tendency', &
         UNITS      = 'kg kg-1 s-1',                               &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DO3DT_ANA',                                 &
         LONG_NAME  = 'total_ozone_analysis_tendency',             &
         UNITS      = 'mol mol-1 s-1',                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DTSDT_ANA',                                 &
         LONG_NAME  = 'total_skin_temperature_tendency',           &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'DTHVDTFILINT',                                     &
         LONG_NAME  = 'vertically_integrated_thv_adjustment_from_filling',&
         UNITS      = 'K kg m-2 s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'PERES',                                            &
         LONG_NAME  = 'vertically_integrated_cpt_tendency_residual',      &
         UNITS      = 'W m-2',                                            &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'PEFILL',                                           &
         LONG_NAME  = 'vertically_integrated_cpt_adjustment_from_filling',&
         UNITS      = 'W m-2',                                            &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'QTFILL',                                           &
         LONG_NAME  = 'vertically_integrated_total_water_adjustment_from_filling', &
         UNITS      = 'kg m-2 s-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'QVFILL',                                           &
         LONG_NAME  = 'vertically_integrated_qv_adjustment_from_filling', &
         UNITS      = 'kg m-2 s-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'QLFILL',                                           &
         LONG_NAME  = 'vertically_integrated_ql_adjustment_from_filling', &
         UNITS      = 'kg m-2 s-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'QIFILL',                                           &
         LONG_NAME  = 'vertically_integrated_qi_adjustment_from_filling', &
         UNITS      = 'kg m-2 s-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'OXFILL',                                           &
         LONG_NAME  = 'vertically_integrated_ox_adjustment_from_filling', &
         UNITS      = 'kg m-2 s-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPP_EPV',                                                 &
       LONG_NAME          = 'tropopause_pressure_based_on_EPV_estimate',                 &
       UNITS              = 'Pa',                                                        &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPP_THERMAL',                                             &
       LONG_NAME          = 'tropopause_pressure_based_on_thermal_estimate',             &
       UNITS              = 'Pa',                                                        &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPP_BLENDED',                                             &
       LONG_NAME          = 'tropopause_pressure_based_on_blended_estimate',             &
       UNITS              = 'Pa',                                                        &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPT',                                                     &
       LONG_NAME          = 'tropopause_temperature_using_blended_TROPP_estimate',       &
       UNITS              = 'K',                                                         &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPQ',                                                     &
       LONG_NAME          = 'tropopause_specific_humidity_using_blended_TROPP_estimate', &
       UNITS              = 'kg kg-1',                                                   &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'TQV',                                        &
         LONG_NAME        = 'total_precipitable_water_vapor',             &
         UNITS            = 'kg m-2'  ,                                   &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'TQI',                                        &
         LONG_NAME        = 'total_precipitable_ice_water',               &
         UNITS            = 'kg m-2'  ,                                   &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'TQL',                                        &
         LONG_NAME        = 'total_precipitable_liquid_water',            &
         UNITS            = 'kg m-2'  ,                                   &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'TOX',                                        &
         LONG_NAME        = 'total_column_odd_oxygen',                    &
         UNITS            = 'kg m-2'  ,                                   &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'MASS',                                       &
         LONG_NAME        = 'atmospheric_mass',                           &
         UNITS            = 'kg m-2'  ,                                   &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'KE',                                         &
         LONG_NAME        = 'vertically_integrated_kinetic_energy',       &
         UNITS            = 'J m-2'  ,                                    &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'CPT',                                        &
         LONG_NAME        = 'vertically_integrated_enthalpy',             &
         UNITS            = 'J m-2'  ,                                    &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'THV',                                        &
         LONG_NAME        = 'vertically_integrated_virtual_potential_temperature',             &
         UNITS            = 'K'  ,                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'QLTOT',                                      &
         LONG_NAME        = 'mass_fraction_of_cloud_liquid_water',        &
         UNITS            = 'kg kg-1',                                    &
         DIMS             = MAPL_DimsHorzVert,                            &
         VLOCATION        = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'QITOT',                                      &
         LONG_NAME        = 'mass_fraction_of_cloud_ice_water',           &
         UNITS            = 'kg kg-1',                                    &
         DIMS             = MAPL_DimsHorzVert,                            &
         VLOCATION        = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'PHIS',                                       &
         LONG_NAME        = 'surface geopotential height',                &
         UNITS            = 'm+2 s-2',                                    &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'SGH',                                        &
         LONG_NAME        = 'isotropic stdv of GWD topography',           &
         UNITS            = 'm',                                          &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'GWDVARX',                                    &
         LONG_NAME        = 'east-west variance of GWD topography',       &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'GWDVARY',                                    &
         LONG_NAME        = 'north-south variance of GWD topography',     &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'GWDVARXY',                                   &
         LONG_NAME        = 'SW-NE variance of GWD topography',           &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'GWDVARYX',                                   &
         LONG_NAME        = 'NW-SE variance of GWD topography',           &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'TRBVAR',                                     &
         LONG_NAME        = 'isotropic variance of TRB topography',       &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'VARFLT',                                     &
         LONG_NAME        = 'isotropic variance of filtered topography',  &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)


! Create childrens gridded components and invoke their SetServices
! ----------------------------------------------------------------
#ifdef SCM
    SDYN = MAPL_AddChild(GC, NAME='SCMDYNAMICS', SS=SDYN_SetServices, RC=STATUS)
    VERIFY_(STATUS)
#else
    SDYN = MAPL_AddChild(GC, NAME='SUPERDYNAMICS', SS=SDYN_SetServices, RC=STATUS)
    VERIFY_(STATUS)
#endif
    PHYS = MAPL_AddChild(GC, NAME='PHYSICS', SS=PHYS_SetServices, RC=STATUS)
    VERIFY_(STATUS)

    ORB  = MAPL_AddChild(GC, NAME='ORBIT', SS=ORB_SetServices, RC=STATUS)
    VERIFY_(STATUS)

! Export for IAU or Analysis purposes
! -----------------------------------
!   call MAPL_AddExportSpec ( GC, &
!        SHORT_NAME = 'PHIS', &
!        CHILD_ID = SDYN, &
!        RC=STATUS )
!   VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'AREA', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'AK', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'BK', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'PLE', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec( GC, &
         SHORT_NAME = 'PS', &
         CHILD_ID = SDYN, &
         RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'PE', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'PT', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'TV', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'T', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'U', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'V', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'U_DGRID', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'V_DGRID', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'O3PPMV', &
         CHILD_ID   = PHYS,  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'OX',  &
         CHILD_ID   = PHYS,  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'Q', &
         CHILD_ID = PHYS, &
         RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'QCTOT', &
         CHILD_ID = PHYS, &
         RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'U10N', &
         CHILD_ID = PHYS, &
         RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'V10N', &
         CHILD_ID = PHYS, &
         RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'SNOMAS', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'WET1',  &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'TSOIL1', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'LWI', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'Z0', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'TS', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'TRANA', &
         CHILD_ID = PHYS, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'FRLAND', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'FRLANDICE', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'FRLAKE', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'FROCEAN', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'FRACI', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)
!EOS

! Set internal connections between the childrens IMPORTS and EXPORTS
! ------------------------------------------------------------------

    call MAPL_AddConnectivity ( GC,                                                        &
         SRC_NAME  = (/'U            ','V            ','TH           ','T            ',    &
                       'ZLE          ','PS           ','TA           ','QA           ',    &
                       'SPEED        ','DZ           ','PLE          ',                    &
                       'PREF         ','TROPP_BLENDED','S            ',                    &
                       'PV           '/),                                                  &
         DST_NAME  = (/'U     ','V     ','TH    ','T     ',                                &
                       'ZLE   ','PS    ','TA    ','QA    ',                                &
                       'SPEED ','DZ    ','PLE   ',                                         &
                       'PREF  ','TROPP ','S     ',                                         &
                       'PV    '/),                                                         &
         DST_ID = PHYS,                                                                    &
         SRC_ID = SDYN,                                                                    &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC,              &
         SRC_NAME    = 'PLE',                    &
         DST_NAME    = 'PLEINST',                &
         SRC_ID      = SDYN,                     &
         DST_ID      = PHYS,                     &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC, SRC_NAME = 'AREA', DST_NAME = 'AREA', &
         SRC_ID = SDYN, DST_ID = PHYS, RC=STATUS  )
     VERIFY_(STATUS)

! Bundle of quantities to be advected
!------------------------------------

     call MAPL_AddConnectivity ( GC,                               &
         SRC_NAME  = 'TRADV',                                      &
         DST_NAME  = 'TRADV',                                      &
         SRC_ID      = PHYS,                                       &
         DST_ID      = SDYN,                                       &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

! Orbital Component Bundle
     call MAPL_AddConnectivity( GC,                                &
         SRC_NAME='SATORB',                                        &
         DST_NAME='SATORB',                                        &
         SRC_ID = ORB,                                             &
         DST_ID = PHYS,                                            &
                                                        RC=STATUS  )
     VERIFY_(STATUS) 

#ifdef SCMSURF
    call MAPL_AddConnectivity ( GC,    &
         SHORT_NAME  = (/'TSKINOBS','QSKINOBS','LHOBS   ','SHOBS   '/), &
         DST_ID = PHYS,         &
         SRC_ID = SDYN,         &
         RC=STATUS  )
     VERIFY_(STATUS)
#endif

! We Terminate these IMPORTS which are manually filled
!-----------------------------------------------------

     call MAPL_TerminateImport    ( GC,                                                     &
          SHORT_NAME = (/'DUDT  ','DVDT  ','DTDT  ','DPEDT ','DQVANA','DOXANA','PHIS  '/),  &
          CHILD      = SDYN,                                                                &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_TerminateImport    ( GC,                          &
          SHORT_NAME = (/'VARFLT','PHIS  ','SGH   ', 'DTSDT '/), &
          CHILD      = PHYS,                                     &
          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_GenericSetServices    ( GC, RC=STATUS )
    VERIFY_(STATUS)

! Clocks
!-------

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
    VERIFY_(STATUS)

! All done
!---------

    RETURN_(ESMF_SUCCESS)  
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Initialize -- Initialize method for the composite Agcm Gridded Component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)           :: IAm 
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),  pointer  :: STATE
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
   type (ESMF_Field)                   :: FIELD
   type (ESMF_Time)                    :: CurrTime, RingTime
   type (ESMF_TimeInterval)            :: TIMEINT
   type (ESMF_Alarm)                   :: ALARM
   type (ESMF_Config)                  :: cf
   integer                             :: I, NQ
   real                                :: POFFSET, DT             
   real, pointer, dimension(:,:)       :: PHIS,SGH,VARFLT,PTR
   real, pointer, dimension(:,:,:)     :: TEND!
   character(len=ESMF_MAXSTR)          :: replayMode
   real                                :: RPL_INTERVAL
   real                                :: RPL_SHUTOFF
   character(len=ESMF_MAXSTR), parameter :: INITIALIZED_EXPORTS(3) = &
        (/'PHIS  ', 'SGH   ', 'VARFLT' /)

! =============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, config=cf, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Initialize"

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)


    call MAPL_TimerOn(STATE,"INITIALIZE")

! Call Initialize for every Child

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(STATE,"TOTAL")

! Get children and their im/ex states from my generic state.
!----------------------------------------------------------

    call MAPL_Get ( STATE, GIM=GIM, GEX=GEX, RC=STATUS )
    VERIFY_(STATUS)

! Make sure that the physics tendencies are allocated
!----------------------------------------------------

    call MAPL_GetPointer(GEX(PHYS), TEND, 'DUDT' , ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GEX(PHYS), TEND, 'DVDT' , ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GEX(PHYS), TEND, 'DTDT' , ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GEX(PHYS), TEND, 'DPEDT', ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)

! Fill Childrens TOPO variables and Diagnostics
!----------------------------------------------
    call MAPL_GetPointer(EXPORT, PHIS,   'PHIS',   ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SGH,    'SGH',    ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VARFLT, 'VARFLT', ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)

! PHIS ...
!---------
    call ESMF_StateGet( GIM(SDYN), 'PHIS', FIELD, rc=STATUS )
    VERIFY_(STATUS)
    Call GEOS_TopoGet ( cf, MEAN=FIELD, rc=STATUS )
    VERIFY_(STATUS)

    call ESMF_StateGet( GIM(PHYS), 'PHIS', FIELD, rc=STATUS )
    VERIFY_(STATUS)
    Call GEOS_TopoGet ( cf, MEAN=FIELD, rc=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldGet (FIELD, localDE=0, farrayPtr=PTR, rc = status)
    PHIS = PTR

! GWDVAR ...
!-----------
    call ESMF_StateGet( GIM(PHYS), 'SGH', FIELD, rc=STATUS )
    VERIFY_(STATUS)
    Call GEOS_TopoGet ( cf, GWDVAR=FIELD, rc=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldGet (FIELD, localDE=0, farrayPtr=PTR, rc = status)
    SGH = PTR

! TRBVAR ...
!-----------
    call ESMF_StateGet( GIM(PHYS), 'VARFLT', FIELD, rc=STATUS )
    VERIFY_(STATUS)
    Call GEOS_TopoGet ( cf, TRBVAR=FIELD, rc=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldGet (FIELD, localDE=0, farrayPtr=PTR, rc = status)
    VARFLT = PTR

! ======================================================================
!ALT: the next section addresses the problem when export variables have been
!     assigned values during Initialize. To prevent "connected" exports
!     being overwritten by DEFAULT in the Import spec in the other component
!     we label them as being "initailized by restart". A better solution
!     would be to move the computation to phase 2 of Initialize and
!     eliminate this section alltogether
! ======================================================================
    DO I = 1, size(INITIALIZED_EXPORTS)
       call ESMF_StateGet(EXPORT,INITIALIZED_EXPORTS(I), FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", &
                              VALUE=MAPL_InitialRestart, RC=STATUS)
       VERIFY_(STATUS)      
    END DO

! Initialize Predictor Alarm
!---------------------------

   call ESMF_ClockGet(clock, currTime=currTime, rc=status)
   VERIFY_(STATUS)

   call MAPL_GetResource( STATE, DT, Label="RUN_DT:", RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetResource( STATE, POFFSET, Label="PREDICTOR_OFFSET:", default=21600. , RC=STATUS)
   VERIFY_(STATUS)

   call ESMF_TimeIntervalSet(TIMEINT,  S=nint (POFFSET), RC=STATUS)
   VERIFY_(STATUS)

   ringTime = currTime+TIMEINT

   call ESMF_TimeIntervalSet(TIMEINT,  S=nint(DT) , RC=STATUS)
   VERIFY_(STATUS)

   ALARM = ESMF_AlarmCreate( name='PredictorAlarm', &
                             CLOCK = CLOCK, &
                             RingInterval = TIMEINT  ,  &
                             RingTime     = ringTime,  & 
                             RC           = STATUS      )
   VERIFY_(STATUS)
   if(ringTime == currTime) then
      call ESMF_AlarmRingerOn(Alarm, rc=status)
      VERIFY_(STATUS)
   end if

   call MAPL_StateAlarmAdd(STATE,ALARM,RC=status)
   VERIFY_(STATUS)

   call MAPL_GetResource(STATE, ReplayMode, 'REPLAY_MODE:', default="NoReplay", RC=STATUS )
   VERIFY_(STATUS)

   if(adjustl(ReplayMode)=="Exact") then
      call MAPL_GetResource(STATE, RPL_SHUTOFF, 'REPLAY_SHUTOFF:', default=4000*21600., RC=STATUS )
      VERIFY_(STATUS)
      call ESMF_TimeIntervalSet(TIMEINT, S=nint(RPL_SHUTOFF), RC=STATUS)
      VERIFY_(STATUS)

      ALARM = ESMF_AlarmCreate(name='ReplayShutOff', clock=CLOCK,      &
           ringInterval=TIMEINT, sticky=.false.,    &
           RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_GetResource(STATE, RPL_INTERVAL, 'REPLAY_INTERVAL:', default=21600., RC=STATUS )
      VERIFY_(STATUS)
      call ESMF_TimeIntervalSet(TIMEINT, S=nint(RPL_INTERVAL), RC=STATUS)
      VERIFY_(STATUS)

      ALARM = ESMF_AlarmCreate(name='ExactReplay', clock=CLOCK,      &
           ringInterval=TIMEINT, sticky=.false.,    &
           RC=STATUS )
      VERIFY_(STATUS)
      call ESMF_AlarmRingerOn(ALARM, rc=status)
      VERIFY_(STATUS)
      ASSERT_(POFFSET == RPL_INTERVAL)
   end if

    call MAPL_TimerOff(STATE,"TOTAL")
    call MAPL_TimerOff(STATE,"INITIALIZE")

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": IMPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( IMPORT, rc=STATUS )
    call WRITE_PARALLEL ( trim(Iam)//": EXPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( EXPORT, rc=STATUS )
#endif


    RETURN_(ESMF_SUCCESS)
 end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: Run -- Run method for the composite Agcm Gridded Component

! !INTERFACE:

  subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)           :: IAm 
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),      pointer  :: STATE
   type (ESMF_GridComp),      pointer  :: GCS(:)
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
   type (ESMF_State)                   :: INTERNAL
   type (ESMF_Alarm)                   :: ALARM
   type (ESMF_TimeInterval)            :: TINT
   type (ESMF_FieldBundle)             :: Bundle
   type (ESMF_FieldBundle)             :: Advect_Bundle
   type (ESMF_Grid)                    :: grid

   real, pointer, dimension(:,:,:)     :: U      => null()
   real, pointer, dimension(:,:,:)     :: V      => null()
   real, pointer, dimension(:,:,:)     :: T      => null()
   real, pointer, dimension(:,:,:)     :: Q      => null()
   real, pointer, dimension(:,:,:)     :: QLLS   => null()
   real, pointer, dimension(:,:,:)     :: QLCN   => null()
   real, pointer, dimension(:,:,:)     :: QILS   => null()
   real, pointer, dimension(:,:,:)     :: QICN   => null()
   real, pointer, dimension(:,:,:)     :: PLE    => null()
   real, pointer, dimension(:,:,:)     :: EPV    => null()
   real, pointer, dimension(:,:,:)     :: QLTOT  => null()
   real, pointer, dimension(:,:,:)     :: QITOT  => null()
   real, pointer, dimension(:,:,:)     :: DPEDT  => null()
   real, pointer, dimension(:,:,:)     :: TENDAN => null()

   real,   allocatable, dimension(:,:)   :: QFILL
   real,   allocatable, dimension(:,:)   :: QINT 
   real*8, allocatable, dimension(:,:)   :: SUMKE
   real*8, allocatable, dimension(:,:)   :: SUMCPT1, SUMCPT2
   real*8, allocatable, dimension(:,:)   :: SUMTHV1, SUMTHV2
   real*8, allocatable, dimension(:,:,:) :: PKE
   real*8, allocatable, dimension(:,:,:) :: PKZ

   real, pointer, dimension(:,:)       :: AREA   => null()
   real, pointer, dimension(:,:)       :: OXFILL => null()
   real, pointer, dimension(:,:)       :: QTFILL => null()
   real, pointer, dimension(:,:)       :: QVFILL => null()
   real, pointer, dimension(:,:)       :: QIFILL => null()
   real, pointer, dimension(:,:)       :: QLFILL => null()
   real, pointer, dimension(:,:)       :: TQV    => null()
   real, pointer, dimension(:,:)       :: TQI    => null()
   real, pointer, dimension(:,:)       :: TQL    => null()
   real, pointer, dimension(:,:)       :: TOX    => null()
   real, pointer, dimension(:,:)       :: TROPP1 => null()
   real, pointer, dimension(:,:)       :: TROPP2 => null()
   real, pointer, dimension(:,:)       :: TROPP3 => null()
   real, pointer, dimension(:,:)       :: TROPT  => null()
   real, pointer, dimension(:,:)       :: TROPQ  => null()
   real, pointer, dimension(:,:)       :: MASS   => null()
   real, pointer, dimension(:,:)       :: KE     => null()
   real, pointer, dimension(:,:)       :: CPT    => null()
   real, pointer, dimension(:,:)       :: THV    => null()

   real, pointer, dimension(:,:)       :: PERES        => null()
   real, pointer, dimension(:,:)       :: PEFILL       => null()
   real, pointer, dimension(:,:)       :: PEPHY_SDYN   => null()  ! D(CpT)DT from ADD_INCS (SuperDYNamics)
   real, pointer, dimension(:,:)       :: PEPHY_PHYS   => null()  ! D(CpT)DT from PHYSics 

   real, pointer, dimension(:,:)       :: DTHVDTFILINT => null()
   real, pointer, dimension(:,:)       :: DTHVDTPHYINT => null()
   real, pointer, dimension(:,:)       :: DQVDTPHYINT  => null()
   real, pointer, dimension(:,:)       :: DQLDTPHYINT  => null()
   real, pointer, dimension(:,:)       :: DQIDTPHYINT  => null()
   real, pointer, dimension(:,:)       :: DOXDTPHYINT  => null()

   real, pointer, dimension(:,:,:)     :: DP
   real, pointer, dimension(:,:,:)     :: PL
   real, pointer, dimension(:,:,:)     :: FC
   real, pointer, dimension(:,:,:)     :: TROP
   real, pointer, dimension(:,:,:)     :: XXINC
   real, pointer, dimension(:,:,:)     :: DQVANA
   real, pointer, dimension(:,:,:)     :: DOXANA
   real, pointer, dimension(:,:,:)     :: DQIMPORT => null()

   real, pointer, dimension(:,:)       :: ptr2d
   real, pointer, dimension(:,:,:)     :: ptr3d

   real, pointer, dimension(:)         :: AK
   real, pointer, dimension(:)         :: BK

   real*8, allocatable, dimension(:,:)   :: sumq
   real,   allocatable, dimension(:,:,:) ::  ple_ana
   real,   allocatable, dimension(:,:,:) ::   dp_ana

   real*8                                ::   gamma
   real*8                                ::   ps_ana_ave
   real*8                                :: pdry_bkg_ave
   real*8                                :: qint_bkg_ave
   real*8                                :: qint_ana_ave

   real*8                                ::   qvint_ana_ave
   real*8                                :: qllsint_ana_ave, qlcnint_ana_ave
   real*8                                :: qilsint_ana_ave, qicnint_ana_ave

   real*8                                ::   qvint_bkg_ave
   real*8                                :: qllsint_bkg_ave, qlcnint_bkg_ave
   real*8                                :: qilsint_bkg_ave, qicnint_bkg_ave

   real*8                                :: dqana2_ave

   real                                :: DT
   integer                             :: IM, JM, LM, L, TYPE, ISFCST
   integer                             :: NumFriendly
   integer                             :: K
   integer                             :: I
   logical                             :: DasMode
   logical                             :: DO_PREDICTOR
   logical                             :: LAST_CORRECTOR
   integer                             :: CONSTRAIN_DAS
   real                                :: ALPHA, BETA, TAUANL, DTX
   real                                :: ALPHAQ, BETAQ
   real                                :: ALPHAO, BETAO
   real                                :: ALF, BET
   real                                :: EPS

   character(len=ESMF_MAXSTR), pointer :: Names(:)
   character(len=ESMF_MAXSTR)          :: STRING

   integer, parameter                  :: FREERUN   = 0
   integer, parameter                  :: PREDICTOR = 1
   integer, parameter                  :: CORRECTOR = 2
   integer, parameter                  :: FORECAST  = 3

   integer                             :: unit
   logical                             :: is_ringing
   logical,save                        :: is_shutoff=.false.
   character(len=ESMF_MAXSTR)          :: FILENAME
   character(len=ESMF_MAXSTR)          :: FileTmpl
   character(len=ESMF_MAXSTR)          :: replayFile
   character(len=ESMF_MAXSTR)          :: replayMode
   character(len=ESMF_MAXSTR)          :: rplMode
   type(ESMF_Time)                     :: currTime

   CHARACTER(LEN=ESMF_MAXSTR)          :: fieldName

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (STATE,"TOTAL")
    call MAPL_TimerOn (STATE,"RUN"  )

! Get children and their im/ex states from my generic state.
!----------------------------------------------------------

    call MAPL_Get ( STATE, GCS=GCS, GIM=GIM, GEX=GEX,  &
                                INTERNAL_ESMF_STATE=INTERNAL,      &
                                IM=IM, JM=JM, LM=LM,               & 
                                RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, grid=grid, rc=status)
    VERIFY_(STATUS)

! Get the specific IAU alarm
!---------------------------

    call MAPL_StateAlarmGet(STATE, ALARM, NAME='PredictorAlarm', RC=STATUS)
    VERIFY_(STATUS)

! Set the various time scales
!----------------------------

    call MAPL_GetResource( STATE,     DT,          Label="RUN_DT:",                        RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE,  ALPHA,          Label="ALPHA:",         default=0.0,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE,   BETA,          Label="BETA:",          default=1.0,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE, ALPHAQ,          Label="ALPHAQ:",        default=ALPHA,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE,  BETAQ,          Label="BETAQ:",         default=BETA,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE, ALPHAO,          Label="ALPHAO:",        default=ALPHA,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE,  BETAO,          Label="BETAO:",         default=BETA,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE, TAUANL,          Label="TAUANL:",        default=21600., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE, ISFCST,          Label="IS_FCST:",       default=0,      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE, CONSTRAIN_DAS,   Label="CONSTRAIN_DAS:", default=1,      RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( STATE, STRING,                &
                           LABEL="IMPORT_RESTART_FILE:", &
                           RC=STATUS)
    IF (STATUS == ESMF_SUCCESS) THEN
       DasMode = .true.
    ELSE
       DasMode = .false.
    END IF

    
    call MAPL_GetResource(STATE, ReplayMode, 'REPLAY_MODE:', &
         default="NoReplay", RC=STATUS )
    VERIFY_(STATUS)

! NoReplay, Exact, Intermittent, Regular
    rplMode = adjustl(ReplayMode)
    if(rplMode=="Regular" .or. rplMode == "Exact") then
       DasMode = .true.
    end if

! Set type of update
!-------------------

    if     (ISFCST/=0  ) then
       TYPE = FORECAST
    else if(.not. DasMode) then
       TYPE = FREERUN
    else
       DO_PREDICTOR = ESMF_AlarmIsRinging( ALARM, rc=status)
       VERIFY_(STATUS)
       LAST_CORRECTOR = ESMF_AlarmWillRingNext( ALARM, rc=status)
       VERIFY_(STATUS)
       REPLAYING: if (rplMode == "Regular") then
          call ESMF_ClockGetAlarm(clock, 'startReplay', alarm, rc=status)
          VERIFY_(STATUS)
          LAST_CORRECTOR = ESMF_AlarmWillRingNext( ALARM, rc=status)
          VERIFY_(STATUS)
       else if (rplMode == "Exact") then
          call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
          VERIFY_(STATUS)
          DO_PREDICTOR = .FALSE.
          call MAPL_GetResource ( STATE, FileTmpl,'REPLAY_FILE:', RC=STATUS )
          VERIFY_(STATUS)

! If replay alarm is ringing, we need to reset state
!---------------------------------------------------

          if (is_shutoff) then ! once this alarm rings, is_shutoff will remain true for the rest of the run 
             call MAPL_GetPointer(IMPORT,ptr3d,'DUDT',RC=STATUS)
             ptr3d=0.0
             call MAPL_GetPointer(IMPORT,ptr3d,'DVDT',RC=STATUS)
             ptr3d=0.0
             call MAPL_GetPointer(IMPORT,ptr3d,'DTDT',RC=STATUS)
             ptr3d=0.0
             call MAPL_GetPointer(IMPORT,ptr3d,'DPEDT',RC=STATUS)
             ptr3d=0.0
             call MAPL_GetPointer(IMPORT,ptr3d,'DQVDT',RC=STATUS)
             ptr3d=0.0
             call MAPL_GetPointer(IMPORT,ptr3d,'DO3DT',RC=STATUS)
             ptr3d=0.0
             call MAPL_GetPointer(IMPORT,ptr2d,'DTSDT',RC=STATUS)
             if(associated(ptr2d))then
                ptr2d=0.0
             endif
          else
             call ESMF_ClockGetAlarm(Clock,'ReplayShutOff',Alarm,rc=Status)
             VERIFY_(status) 
             is_shutoff = ESMF_AlarmWillRingNext( Alarm,rc=status )
             VERIFY_(status)
          endif

          call ESMF_ClockGetAlarm(Clock,'ExactReplay',Alarm,rc=Status)
          VERIFY_(status) 

	  LAST_CORRECTOR = ESMF_AlarmWillRingNext( ALARM, rc=status)
          VERIFY_(STATUS)

          is_ringing = ESMF_AlarmIsRinging( Alarm,rc=status )
          VERIFY_(status) 

          is_ringing = is_ringing .and. (.not. is_shutoff)
          TIME_TO_REPLAY: if(is_ringing) then
             call ESMF_ClockGet(Clock, CurrTime=currTime, rc=Status)
             VERIFY_(status) 
             ! when testing check  alarms, predictor/corrector/last corrector
             ! template to file
             call MAPL_GetCurrentFile(FILETMPL=filetmpl, TIME=currTime, FILENAME=ReplayFile, &
                  RC=STATUS)
             VERIFY_(status) 
             unit = getfile(ReplayFile, FORM="unformatted", rc=status)
             VERIFY_(STATUS) 
             ! read import from file
             call MAPL_VarRead(UNIT=UNIT, STATE=IMPORT, RC=STATUS)
             VERIFY_(STATUS) 
             call FREE_FILE(unit, rc=status)
             VERIFY_(STATUS) 
          end if TIME_TO_REPLAY
       end if REPLAYING

       if(DO_PREDICTOR) then
          TYPE = PREDICTOR
       else
          TYPE = CORRECTOR
       end if
    end if

! Get Names Associated with Friendly Analysis Bundle
!---------------------------------------------------

    call ESMF_StateGet(GEX(PHYS), 'TRANA', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet(BUNDLE,FieldCount=NumFriendly,   RC=STATUS)
    VERIFY_(STATUS)

    ASSERT_(NumFriendly==2)

    allocate(Names(NumFriendly), stat=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet(BUNDLE, fieldNameList=Names, RC=STATUS)
    VERIFY_(STATUS)

! Prepare for update
!-------------------

    call MAPL_GetPointer( GEX(SDYN), PLE, 'PLE', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), U,   'U',   rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), V,   'V',   rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), T,   'T',   rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), AK,  'AK',  rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), BK,  'BK',  rc=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetPointer( GEX(SDYN), PEPHY_SDYN, 'PEPHY', alloc=.true., rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(PHYS), PEPHY_PHYS, 'PEPHY', alloc=.true., rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), DTHVDTPHYINT, 'DTHVDTPHYINT', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(PHYS), DQVDTPHYINT,  'DQVDTPHYINT',  rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(PHYS), DQLDTPHYINT,  'DQLDTPHYINT',  rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(PHYS), DQIDTPHYINT,  'DQIDTPHYINT',  rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(PHYS), DOXDTPHYINT,  'DOXDTPHYINT',  rc=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetPointer( GEX(SDYN), AREA, 'AREA', rc=STATUS )
    VERIFY_(STATUS)

    allocate( PL(IM,JM,LM),STAT=STATUS )
    VERIFY_(STATUS)
    allocate( DP(IM,JM,LM),STAT=STATUS )
    VERIFY_(STATUS)

    PL  = 0.5*(PLE(:,:,1:LM)+PLE(:,:,0:LM-1))
    DP  = PLE(:,:,1:LM)-PLE(:,:,0:LM-1)

    ALF = ALPHA
    BET = BETA

! Load Analysis Increments into Imports for Physics and Dynamics State Variables
!-------------------------------------------------------------------------------

    call DO_UPDATE_ANA2D ('DTSDT', PHYS)

! Note: Mass-Weighting is done for DTDT, No Mass-Weighting for DUDT,DVDT,DPEDT
! ----------------------------------------------------------------------------

    nullify(FC)

    call DO_UPDATE_ANA3D ('DUDT' , SDYN)
    call DO_UPDATE_ANA3D ('DVDT' , SDYN)
    call DO_UPDATE_ANA3D ('DPEDT', SDYN, CONSTRAIN_DAS = CONSTRAIN_DAS)

    call MAPL_GetPointer(GIM(SDYN), DPEDT, 'DPEDT', rc=STATUS)
    VERIFY_(STATUS)

    FC => DP

    call DO_UPDATE_ANA3D ('DTDT' , SDYN )

! Add Analysis Increment Directly to Friendlies
!----------------------------------------------

    if(TYPE /= FREERUN) then

       allocate(FC (IM,JM,LM),STAT=STATUS)
       VERIFY_(STATUS)

       do K=1,NumFriendly

          FC = 1.0

          NULLIFY(Q)  !ALT: ESMF requires that the data pointer is not associated
          call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
          VERIFY_(STATUS)

          STRING = TRIM(Names(K))
	  fieldName = MAPL_RmQualifier(STRING)

          if(TRIM(fieldName) == 'OX') then ! PCHEM OX or GOCART::OX, for example.

             ALF = ALPHAO
             BET = BETAO
             DTX = DT*1.0E-6

! Uncomment damping if problem with ozone
! ---------------------------------------
!            do L=1,LM
!               where(PL(:,:,L) < 100.0 .and. PL(:,:,L) > 0.0 )
!                     FC(:,:,L) = exp(-1.5*(log10(PL(:,:,L))-2.0)**2)
!               end where
!            end do

             call MAPL_GetPointer(GIM(SDYN), DOXANA, 'DOXANA', rc=STATUS)
             VERIFY_(STATUS)
             DOXANA = Q
             call DO_Friendly (Q,'DO3DT')
             DOXANA = Q - DOXANA

          else

             ALF = ALPHAQ
             BET = BETAQ
             DTX = DT

             if(NAMES(K)=='Q') then ! Q

                ! Initialize DQVANA Diagnostics with Background QV
                !-------------------------------------------------
                call MAPL_GetPointer(GIM(SDYN), DQVANA, 'DQVANA', rc=STATUS)
                VERIFY_(STATUS)

                DQVANA = Q

                if(TYPE  == CORRECTOR) then

                   allocate(     qint( IM,JM) )
                   allocate(     sumq( IM,JM) )
                   allocate(   dp_ana( IM,JM,  LM) )
                   allocate(  ple_ana( IM,JM,0:LM) )

                   ! Get Pointers to QL & QI from Friendly Advection Bundle
                   !-------------------------------------------------------
                   call ESMF_StateGet(GEX(PHYS), 'TRADV', Advect_BUNDLE, RC=STATUS )
                   VERIFY_(STATUS)

                   call ESMFL_BundleGetPointerToData( Advect_BUNDLE, 'QLLS', QLLS, RC=STATUS )
                   VERIFY_(STATUS)
                   call ESMFL_BundleGetPointerToData( Advect_BUNDLE, 'QLCN', QLCN, RC=STATUS )
                   VERIFY_(STATUS)
                   call ESMFL_BundleGetPointerToData( Advect_BUNDLE, 'QILS', QILS, RC=STATUS )
                   VERIFY_(STATUS)
                   call ESMFL_BundleGetPointerToData( Advect_BUNDLE, 'QICN', QICN, RC=STATUS )
                   VERIFY_(STATUS)

                   ! Ensure Global Dry Mass Conservation
                   !------------------------------------
                   IF( CONSTRAIN_DAS == 0 ) then
                       sumq = 0.0_8
                       do L=1,lm
                       sumq = sumq + ( q(:,:,L)+qlls(:,:,L)+qlcn(:,:,L)+qils(:,:,L)+qicn(:,:,L) )*dp(:,:,L)
                       enddo
                       qint = ple(:,:,LM)-sumq
                       call MAPL_AreaMean( pdry_bkg_ave, qint, area, grid, rc=STATUS ); VERIFY_(STATUS)
                   ENDIF
   
                   IF( CONSTRAIN_DAS == 1 ) then  ! BETA = SQRT( QT_ANA )
                   ! ----------------------------------------------------
                       sumq = 0.0_8
                       do L=1,lm
                       sumq = sumq + ( q(:,:,L)+qlls(:,:,L)+qlcn(:,:,L)+qils(:,:,L)+qicn(:,:,L) )*dp(:,:,L)
                       enddo
                       qint = sumq
                       call MAPL_AreaMean( qint_bkg_ave, qint, area, grid, rc=STATUS ); VERIFY_(STATUS)
                   ENDIF
   
                   IF( CONSTRAIN_DAS == 2 ) then  ! BETA = QT_ANA - QT_BKG
                   ! -----------------------------------------------------
#if debug
                       sumq = 0.0_8
                       do L=1,lm
                       sumq = sumq + ( q(:,:,L)+qlls(:,:,L)+qlcn(:,:,L)+qils(:,:,L)+qicn(:,:,L) )*dp(:,:,L)
                       enddo
                       qint = sumq
                       call MAPL_AreaMean( qint_bkg_ave, qint, area, grid, rc=STATUS ); VERIFY_(STATUS)
#endif
                       call DO_GLOBAL_MEAN (q   ,dp,area,  qvint_bkg_ave)
                       call DO_GLOBAL_MEAN (qlls,dp,area,qllsint_bkg_ave)
                       call DO_GLOBAL_MEAN (qlcn,dp,area,qlcnint_bkg_ave)
                       call DO_GLOBAL_MEAN (qils,dp,area,qilsint_bkg_ave)
                       call DO_GLOBAL_MEAN (qicn,dp,area,qicnint_bkg_ave)
                   ENDIF

                   IF( CONSTRAIN_DAS == 3 ) then  ! BETA = QT_ANA
                   ! --------------------------------------------
#if debug
                       sumq = 0.0_8
                       do L=1,lm
                       sumq = sumq + ( q(:,:,L)+qlls(:,:,L)+qlcn(:,:,L)+qils(:,:,L)+qicn(:,:,L) )*dp(:,:,L)
                       enddo
                       qint = sumq
                       call MAPL_AreaMean( qint_bkg_ave, qint, area, grid, rc=STATUS ); VERIFY_(STATUS)
#endif
                       call DO_GLOBAL_MEAN (q   ,dp,area,  qvint_bkg_ave)
                       call DO_GLOBAL_MEAN (qlls,dp,area,qllsint_bkg_ave)
                       call DO_GLOBAL_MEAN (qlcn,dp,area,qlcnint_bkg_ave)
                       call DO_GLOBAL_MEAN (qils,dp,area,qilsint_bkg_ave)
                       call DO_GLOBAL_MEAN (qicn,dp,area,qicnint_bkg_ave)
                   ENDIF

                ENDIF  ! End CORRECTOR Test
   
                ! Update Water Vapor from Analysis
                !---------------------------------

                 call DO_Friendly (Q,'DQVDT')

                if(TYPE  == CORRECTOR) then

                   ! Create Proxy for Updated Pressure due to Analysis Increment
                   !------------------------------------------------------------
                   ple_ana = ple + dt*dpedt
                    dp_ana = ple_ana(:,:,1:LM)-ple_ana(:,:,0:LM-1)

                   ! Constrain Global Analysis Increments of Water Variables to Vanish
                   !------------------------------------------------------------------

                   ! Ensure Global Dry Mass Conservation
                   !------------------------------------
                   IF( CONSTRAIN_DAS == 0 ) then
                       sumq = 0.0_8
                       do L=1,lm
                       sumq = sumq + ( q(:,:,L)+qlls(:,:,L)+qlcn(:,:,L)+qils(:,:,L)+qicn(:,:,L) )*dp_ana(:,:,L)
                       enddo
                       qint = sumq
                       call MAPL_AreaMean( qint_ana_ave, qint,            area, grid, rc=STATUS ); VERIFY_(STATUS)
                       call MAPL_AreaMean(   ps_ana_ave, ple_ana(:,:,LM), area, grid, rc=STATUS ); VERIFY_(STATUS)

                       ! Modify Water Variables to Conserve Dry Mass
                       !--------------------------------------------
                       gamma = ( ps_ana_ave - pdry_bkg_ave ) / qint_ana_ave
                          q  =    q*gamma
                       qlls  = qlls*gamma
                       qlcn  = qlcn*gamma
                       qils  = qils*gamma
                       qicn  = qicn*gamma
                   ENDIF
   
                   IF( CONSTRAIN_DAS == 1 ) then  ! BETA = SQRT( QT_ANA )
                   ! ----------------------------------------------------
                       sumq = 0.0_8
                       do L=1,lm
                       sumq = sumq + ( q(:,:,L)+qlls(:,:,L)+qlcn(:,:,L)+qils(:,:,L)+qicn(:,:,L) )*dp_ana(:,:,L)
                       enddo
                       qint = sumq
                       call MAPL_AreaMean( qint_ana_ave, qint, area, grid, rc=STATUS ); VERIFY_(STATUS)

                       ! Modify Water Variables to Conserve Dry Mass
                       !--------------------------------------------
                       gamma = qint_bkg_ave / qint_ana_ave
                          q  =    q*gamma
                       qlls  = qlls*gamma
                       qlcn  = qlcn*gamma
                       qils  = qils*gamma
                       qicn  = qicn*gamma
                   ENDIF
   
                   IF( CONSTRAIN_DAS == 2 ) then  ! BETA = QT_ANA - QT_BKG
                   ! -----------------------------------------------------
                       call DO_GLOBAL_MEAN (q   ,dp_ana,area,  qvint_ana_ave)
                       call DO_GLOBAL_MEAN (qlls,dp_ana,area,qllsint_ana_ave)
                       call DO_GLOBAL_MEAN (qlcn,dp_ana,area,qlcnint_ana_ave)
                       call DO_GLOBAL_MEAN (qils,dp_ana,area,qilsint_ana_ave)
                       call DO_GLOBAL_MEAN (qicn,dp_ana,area,qicnint_ana_ave)

                       sumq = 0.0_8
                       do L=1,lm
                       sumq = sumq + ( q(:,:,L)-dqvana(:,:,L) )**2 *dp_ana(:,:,L)
                       enddo
                       qint = sumq
                       call MAPL_AreaMean( dqana2_ave, qint, area, grid, rc=STATUS )
                       VERIFY_(STATUS)

                       ! Modify Water Variables to Conserve Dry Mass
                       !--------------------------------------------
                       qlls = qlls - ( q-dqvana )**2 *( qllsint_ana_ave - qllsint_bkg_ave )/dqana2_ave
                       qlcn = qlcn - ( q-dqvana )**2 *( qlcnint_ana_ave - qlcnint_bkg_ave )/dqana2_ave
                       qils = qils - ( q-dqvana )**2 *( qilsint_ana_ave - qilsint_bkg_ave )/dqana2_ave
                       qicn = qicn - ( q-dqvana )**2 *( qicnint_ana_ave - qicnint_bkg_ave )/dqana2_ave
                       q    = q    - ( q-dqvana )**2 *(   qvint_ana_ave -   qvint_bkg_ave )/dqana2_ave
                   ENDIF

                   IF( CONSTRAIN_DAS == 3 ) then  ! BETA = QT_ANA
                   ! --------------------------------------------
                       call DO_GLOBAL_MEAN (q   ,dp_ana,area,  qvint_ana_ave)
                       call DO_GLOBAL_MEAN (qlls,dp_ana,area,qllsint_ana_ave)
                       call DO_GLOBAL_MEAN (qlcn,dp_ana,area,qlcnint_ana_ave)
                       call DO_GLOBAL_MEAN (qils,dp_ana,area,qilsint_ana_ave)
                       call DO_GLOBAL_MEAN (qicn,dp_ana,area,qicnint_ana_ave)

                       sumq = 0.0_8
                       do L=1,lm
                       sumq = sumq + ( q(:,:,L)+qlls(:,:,L)+qlcn(:,:,L)+qils(:,:,L)+qicn(:,:,L) )**2 *dp_ana(:,:,L)
                       enddo
                       qint = sumq
                       call MAPL_AreaMean( dqana2_ave, qint, area, grid, rc=STATUS )
                       VERIFY_(STATUS)

                       ! Modify Water Variables to Conserve Dry Mass
                       !--------------------------------------------
                       do L=1,lm
                             sumq  =  ( q(:,:,L)+qlls(:,:,L)+qlcn(:,:,L)+qils(:,:,L)+qicn(:,:,L) )**2
                       qlls(:,:,L) = qlls(:,:,L) - sumq*( qllsint_ana_ave - qllsint_bkg_ave )/dqana2_ave
                       qlcn(:,:,L) = qlcn(:,:,L) - sumq*( qlcnint_ana_ave - qlcnint_bkg_ave )/dqana2_ave
                       qils(:,:,L) = qils(:,:,L) - sumq*( qilsint_ana_ave - qilsint_bkg_ave )/dqana2_ave
                       qicn(:,:,L) = qicn(:,:,L) - sumq*( qicnint_ana_ave - qicnint_bkg_ave )/dqana2_ave
                       q   (:,:,L) = q   (:,:,L) - sumq*(   qvint_ana_ave -   qvint_bkg_ave )/dqana2_ave
                       enddo
                   ENDIF

#if debug
                   ! Diagnostic Print
                   !-----------------
                   IF( CONSTRAIN_DAS == 0 ) then
                       sumq = 0.0_8
                       do L=1,lm
                       sumq = sumq + ( q(:,:,L)+qlls(:,:,L)+qlcn(:,:,L)+qils(:,:,L)+qicn(:,:,L) )*dp_ana(:,:,L)
                       enddo
                       qint = ple_ana(:,:,LM)-sumq
                       call MAPL_AreaMean(   ps_ana_ave, qint, area, grid, rc=STATUS ); VERIFY_(STATUS)

                       if(MAPL_AM_I_ROOT() ) then
                          write(6,1000) ps_ana_ave*0.01, pdry_bkg_ave*0.01, (ps_ana_ave-pdry_bkg_ave)*0.01
  1000                    format(5x,'PDRY_ANA: ',g,'  PDRY_BKG: ',g,'  DIFF: ',g)
                       endif
                   else
                       sumq = 0.0_8
                       do L=1,lm
                       sumq = sumq + ( q(:,:,L)+qlls(:,:,L)+qlcn(:,:,L)+qils(:,:,L)+qicn(:,:,L) )*dp_ana(:,:,L)
                       enddo
                       qint = sumq
                       call MAPL_AreaMean( qint_ana_ave, qint, area, grid, rc=STATUS )
                       VERIFY_(STATUS)

                       if(MAPL_AM_I_ROOT() ) then
                          write(6,1001) qint_ana_ave,qint_bkg_ave,qint_ana_ave-qint_bkg_ave
  1001                    format(5x,'GLO_QANA_DP: ',g,'  GLO_QBKG_DP: ',g,'  DIFF: ',g)
                       endif
                   endif
#endif
                   deallocate(    qint )
                   deallocate(    sumq )
                   deallocate(  dp_ana )
                   deallocate( ple_ana )

                ENDIF  ! End CORRECTOR Test

                DQVANA = Q - DQVANA

                ! Update Tendency Diagnostic due to CONSTRAINTS
                ! ---------------------------------------------
                call MAPL_GetPointer ( EXPORT, TENDAN, 'DQVDT_ANA', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(TENDAN)) TENDAN = DQVANA/DT

             else

                call DO_Friendly (Q,'D'//trim(Names(K))//'DT')

             end if

          end if

       end do
       deallocate(FC)

    end if ! not free-running

! Make Sure EPV is Allocated for TROPOPAUSE Diagnostics
!------------------------------------------------------
    call MAPL_GetPointer ( EXPORT, TROPP1, 'TROPP_THERMAL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TROPP2, 'TROPP_EPV'    , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TROPP3, 'TROPP_BLENDED', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TROPT, 'TROPT', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TROPQ, 'TROPQ', rc=STATUS )
    VERIFY_(STATUS)

    if( associated(TROPP1) .or. &
        associated(TROPP2) .or. &
        associated(TROPP3) .or. &
        associated(TROPT)  .or. &
        associated(TROPQ)       ) then
        call MAPL_GetPointer( GEX(SDYN),EPV,'EPV',ALLOC=.true.,rc=STATUS )
        VERIFY_(STATUS)
    endif

! Call basic run phase for both Child
!-------------------------------------

! Call run for satellite orbits
!------------------------------
    call MAPL_TimerOn (STATE,"ORBIT"  )
    call ESMF_GridCompRun(GCS(ORB), importState=GIM(ORB), exportState=GEX(ORB), clock=CLOCK, userRC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,"ORBIT"  )


    call MAPL_TimerOn (STATE,"SUPERDYNAMICS"  )
    call ESMF_GridCompRun(GCS(SDYN), importState=GIM(SDYN), exportState=GEX(SDYN), clock=CLOCK, PHASE=1, userRC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerOFF (STATE,"SUPERDYNAMICS"  )

    call MAPL_TimerOn (STATE,"PHYSICS"  )
    call ESMF_GridCompRun(GCS(PHYS), importState=GIM(PHYS), exportState=GEX(PHYS), clock=CLOCK, PHASE=1, userRC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,"PHYSICS"  )



! Load Physics Tendencies into Imports for RUN2 of Dynamics (ADD_INCS)
!---------------------------------------------------------------------

    call DO_UPDATE_PHY ('DUDT' )
    call DO_UPDATE_PHY ('DVDT' )
    call DO_UPDATE_PHY ('DPEDT')
    call DO_UPDATE_PHY ('DTDT' )


! Run RUN2 of SuperDynamics (ADD_INCS) to add Physics Diabatic Tendencies
!------------------------------------------------------------------------

    call MAPL_TimerOn (STATE,"SUPERDYNAMICS"  )

    call ESMF_GridCompRun(GCS(SDYN), importState=GIM(SDYN), exportState=GEX(SDYN), clock=CLOCK, PHASE=2, userRC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOff(STATE,"SUPERDYNAMICS"  )

! Get Names Associated with Friendly Advection Bundle for Final Check for Negative Tracers
!-----------------------------------------------------------------------------------------

    call ESMF_StateGet(GEX(PHYS), 'TRADV', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet(BUNDLE,FieldCount=NumFriendly,   RC=STATUS)
    VERIFY_(STATUS)

    deallocate(Names)
      allocate(Names(NumFriendly), stat=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet(BUNDLE, fieldNameList=Names, RC=STATUS)
    VERIFY_(STATUS)

! Get Pointers to Exports
!------------------------
    call MAPL_GetPointer ( EXPORT, QTFILL, 'QTFILL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, QVFILL, 'QVFILL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, QIFILL, 'QIFILL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, QLFILL, 'QLFILL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, OXFILL, 'OXFILL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TOX   , 'TOX'   , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TQV   , 'TQV'   , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TQI   , 'TQI'   , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TQL   , 'TQL'   , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, QLTOT , 'QLTOT' , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, QITOT , 'QITOT' , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, PERES       , 'PERES'        , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, PEFILL      , 'PEFILL'       , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, DTHVDTFILINT, 'DTHVDTFILINT' , rc=STATUS )
    VERIFY_(STATUS)

    if(associated(QTFILL)) QTFILL = 0.0
    if(associated(QIFILL)) QIFILL = 0.0
    if(associated(QLFILL)) QLFILL = 0.0
    if(associated(TQI)   ) TQI    = 0.0
    if(associated(TQL)   ) TQL    = 0.0
    if(associated(QLTOT) ) QLTOT  = 0.0
    if(associated(QITOT) ) QITOT  = 0.0

    allocate(QFILL(IM,JM)    ,STAT=STATUS )
    VERIFY_(STATUS)
    allocate(QINT (IM,JM)    ,STAT=STATUS )
    VERIFY_(STATUS)
    allocate( PKE(IM,JM,0:LM),STAT=STATUS )
    VERIFY_(STATUS)
    allocate( PKZ(IM,JM,1:LM),STAT=STATUS )
    VERIFY_(STATUS)

    PL  = 0.5*(PLE(:,:,1:LM)+PLE(:,:,0:LM-1))  ! Recompute Updated Pressure
    DP  =      PLE(:,:,1:LM)-PLE(:,:,0:LM-1)   ! Recompute Updated Pressure Thickness

    PKE = PLE**MAPL_KAPPA
    do L=1,LM
    PKZ(:,:,L) = ( PKE(:,:,L)-PKE(:,:,L-1) ) / ( MAPL_KAPPA*( log(PLE(:,:,L))-log(PLE(:,:,L-1)) ) )
    enddo

! Initialize Vertically Integrated Values of CPT and THV (before QFILL updates)
! -----------------------------------------------------------------------------
    if( associated(PEPHY_SDYN) .or. associated(PEFILL) ) then
        EPS = MAPL_RVAP/MAPL_RGAS-1.0
        do K=1,NumFriendly
           NULLIFY(Q)
           call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
           VERIFY_(STATUS)
           if(NAMES(K)=='Q') then
              allocate( SUMCPT1(IM,JM),STAT=STATUS )
              VERIFY_(STATUS)
              SUMCPT1 = 0.0
              do L=1,LM
              SUMCPT1 = SUMCPT1 + MAPL_CP*T(:,:,L)*(1.0+EPS*Q(:,:,L))*DP(:,:,L)
              enddo
              exit
           endif
        enddo
    endif

    if( associated(DTHVDTPHYINT) .or. associated(DTHVDTFILINT) ) then
        EPS = MAPL_RVAP/MAPL_RGAS-1.0
        do K=1,NumFriendly
           NULLIFY(Q)
           call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
           VERIFY_(STATUS)
           if(NAMES(K)=='Q') then
              allocate( SUMTHV1(IM,JM),STAT=STATUS )
              VERIFY_(STATUS)
              SUMTHV1 = 0.0
              do L=1,LM
              SUMTHV1 = SUMTHV1 + T(:,:,L)/PKZ(:,:,L)*(1.0+EPS*Q(:,:,L))*DP(:,:,L)
              enddo
              exit
           endif
        enddo
    endif

! Perform Final Check for Negative Friendlies
! -------------------------------------------

    do K=1,NumFriendly
       NULLIFY(Q)
       call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
       VERIFY_(STATUS)

       STRING = TRIM(Names(K))
       fieldName = MAPL_RmQualifier(STRING)

! Water Vapor
! -----------
       if(NAMES(K)=='Q') then
          call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(QVFILL))           QVFILL =               QFILL
          if(associated(QTFILL))           QTFILL = QTFILL      + QFILL
          if(associated(DQVDTPHYINT)) DQVDTPHYINT = DQVDTPHYINT + QFILL
          if(associated(TQV))                 TQV = QINT
       endif

! Ice Water
! ---------
       if(NAMES(K)=='QICN') then
          !call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(QIFILL))           QIFILL = QIFILL      + QFILL
          if(associated(QTFILL))           QTFILL = QTFILL      + QFILL
          if(associated(DQIDTPHYINT)) DQIDTPHYINT = DQIDTPHYINT + QFILL
          if(associated(TQI))                 TQI = TQI         + QINT
          if(associated(QITOT))             QITOT = QITOT       + Q
       endif

       if(NAMES(K)=='QILS') then
          !call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(QIFILL))           QIFILL = QIFILL      + QFILL
          if(associated(QTFILL))           QTFILL = QTFILL      + QFILL
          if(associated(DQIDTPHYINT)) DQIDTPHYINT = DQIDTPHYINT + QFILL
          if(associated(TQI))                 TQI = TQI         + QINT
          if(associated(QITOT))             QITOT = QITOT       + Q
       endif

! Liquid Water
! ------------
       if(NAMES(K)=='QLCN') then
          !call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(QLFILL))           QLFILL = QLFILL      + QFILL
          if(associated(QTFILL))           QTFILL = QTFILL      + QFILL
          if(associated(DQLDTPHYINT)) DQLDTPHYINT = DQLDTPHYINT + QFILL
          if(associated(TQL))                 TQL = TQL         + QINT
          if(associated(QLTOT))             QLTOT = QLTOT       + Q
       endif

       if(NAMES(K)=='QLLS') then
          !call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(QLFILL))           QLFILL = QLFILL      + QFILL
          if(associated(QTFILL))           QTFILL = QTFILL      + QFILL
          if(associated(DQLDTPHYINT)) DQLDTPHYINT = DQLDTPHYINT + QFILL
          if(associated(TQL))                 TQL = TQL         + QINT
          if(associated(QLTOT))             QLTOT = QLTOT       + Q
       endif

! Total Odd-Oxygen
! ----------------
       if(TRIM(fieldName) == 'OX') then
          call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(OXFILL))           OXFILL = QFILL              *(MAPL_O3MW/MAPL_AIRMW)
          if(associated(DOXDTPHYINT)) DOXDTPHYINT = DOXDTPHYINT + QFILL*(MAPL_O3MW/MAPL_AIRMW)
          if(associated(TOX))                 TOX = QINT               *(MAPL_O3MW/MAPL_AIRMW)
       endif

    enddo
    deallocate(QFILL)
    deallocate(QINT )


! Compute Additional Diagnostics
!-------------------------------

    call MAPL_GetPointer ( EXPORT, MASS , 'MASS' , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, KE   , 'KE'   , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, CPT  , 'CPT'  , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, THV  , 'THV'  , rc=STATUS )
    VERIFY_(STATUS)

    if( associated(MASS) ) MASS = (PLE(:,:,LM)-PLE(:,:,0)) * (1.0/MAPL_GRAV)

    if( associated(KE) ) then
        allocate( SUMKE(IM,JM),STAT=STATUS )
        VERIFY_(STATUS)
        SUMKE = 0.0
        do L=1,LM
        SUMKE = SUMKE + 0.5*( U(:,:,L)**2 + V(:,:,L)**2 )*DP(:,:,L)
        enddo
           KE = SUMKE * (1.0/MAPL_GRAV)
        deallocate(SUMKE)
    endif

! Instantaneous Values of CPT and THV are done here to include possible QFILL updates
! -----------------------------------------------------------------------------------
    if( associated(CPT) .or. associated(PEPHY_SDYN) .or. associated(PEFILL) ) then
        EPS = MAPL_RVAP/MAPL_RGAS-1.0
        do K=1,NumFriendly
           NULLIFY(Q)
           call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
           VERIFY_(STATUS)
           if(NAMES(K)=='Q') then
              allocate( SUMCPT2(IM,JM),STAT=STATUS )
              VERIFY_(STATUS)
              SUMCPT2 = 0.0
              do L=1,LM
              SUMCPT2 = SUMCPT2 + MAPL_CP*T(:,:,L)*(1.0+EPS*Q(:,:,L))*DP(:,:,L)
              enddo
              if( associated(CPT)         )       CPT = SUMCPT2 * (1.0/MAPL_GRAV)
              if( associated(PEFILL) .or. &
                  associated(PEPHY_SDYN)  )   then
                                              SUMCPT1 = ( SUMCPT2-SUMCPT1 ) / (DT*MAPL_GRAV)
              if( associated(PEFILL) )         PEFILL = SUMCPT1
              if( associated(PEPHY_SDYN) ) PEPHY_SDYN = PEPHY_SDYN + SUMCPT1
              deallocate(SUMCPT1)
              endif
              deallocate(SUMCPT2)
              exit
           endif
        enddo
    endif

    if( associated(PERES)       .and.  &
        associated(PEPHY_SDYN)  .and.  &
        associated(PEPHY_PHYS) ) PERES = PEPHY_SDYN - PEPHY_PHYS

    if( associated(THV) .or. associated(DTHVDTPHYINT) .or. associated(DTHVDTFILINT) ) then
        EPS = MAPL_RVAP/MAPL_RGAS-1.0
        do K=1,NumFriendly
           NULLIFY(Q)
           call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
           VERIFY_(STATUS)
           if(NAMES(K)=='Q') then
              allocate( SUMTHV2(IM,JM),STAT=STATUS )
              VERIFY_(STATUS)
              SUMTHV2 = 0.0
              do L=1,LM
              SUMTHV2 = SUMTHV2 + T(:,:,L)/PKZ(:,:,L)*(1.0+EPS*Q(:,:,L))*DP(:,:,L)
              enddo
              if( associated(THV)          )         THV  = SUMTHV2 * (MAPL_P00**MAPL_KAPPA) * (1.0/MAPL_GRAV)
              if( associated(DTHVDTFILINT)  .or. &
                  associated(DTHVDTPHYINT) ) then
                                                  SUMTHV1 = ( SUMTHV2-SUMTHV1 ) * (MAPL_P00**MAPL_KAPPA) / (DT*MAPL_GRAV)
              if( associated(DTHVDTFILINT) ) DTHVDTFILINT = SUMTHV1
              if( associated(DTHVDTPHYINT) ) DTHVDTPHYINT = DTHVDTPHYINT + SUMTHV1
              deallocate(SUMTHV1)
              endif
              deallocate(SUMTHV2)
              exit
           endif
        enddo
    endif


! Compute Tropopause Diagnostics
! ------------------------------
    if( associated(TROPP1) .or. &
        associated(TROPP2) .or. &
        associated(TROPP3) .or. &
        associated(TROPT)  .or. &
        associated(TROPQ)       ) then
        do K=1,NumFriendly
           NULLIFY(Q)
           call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
           VERIFY_(STATUS)
           if(NAMES(K)=='Q') then
              allocate( TROP(IM,JM,5),STAT=STATUS )
              VERIFY_(STATUS)
              call tropovars ( IM,JM,LM,PLE,PL,T,Q,EPV,TROP(:,:,1),TROP(:,:,2),TROP(:,:,3),TROP(:,:,4),TROP(:,:,5) )
               if( associated(TROPP1) )  TROPP1(:,:) = TROP(:,:,1)
               if( associated(TROPP2) )  TROPP2(:,:) = TROP(:,:,2)
               if( associated(TROPP3) )  TROPP3(:,:) = TROP(:,:,3)
               if( associated(TROPT ) )  TROPT (:,:) = TROP(:,:,4)
               if( associated(TROPQ ) )  TROPQ (:,:) = TROP(:,:,5)
                   deallocate(TROP )
               exit
           endif
        enddo
    endif

! Done
!-----

    deallocate(Names)
    deallocate(PL)
    deallocate(DP)
    deallocate(PKE)
    deallocate(PKZ)

    call MAPL_TimerOff(STATE,"RUN"  )
    call MAPL_TimerOff(STATE,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  contains

    subroutine DO_GLOBAL_MEAN (Q,DP,AREA,QAVE)
      real,   intent(IN)    ::    Q(:,:,:)
      real,   intent(IN)    ::   DP(:,:,:)
      real,   intent(IN)    :: AREA(:,:)
      real*8, intent(OUT)   :: QAVE

      real,   allocatable   :: qint(:,:)
      real*8, allocatable   :: sumq(:,:)

      integer               :: L,LM

         LM = size(Q,3)

         allocate( qint(size(q,1),size(q,2)), stat=STATUS)
         VERIFY_(STATUS)
         allocate( sumq(size(q,1),size(q,2)), stat=STATUS)
         VERIFY_(STATUS)

         sumq = 0.0_8
         do L=1,LM
         sumq = sumq + q(:,:,L)*dp(:,:,L)
         enddo
         qint = sumq
         call MAPL_AreaMean( qave, qint, area, grid, rc=STATUS )
         VERIFY_(STATUS)

         deallocate( qint )
         deallocate( sumq )

    end subroutine DO_GLOBAL_MEAN

    subroutine DO_UPDATE_ANA3D(NAME, COMP, CONSTRAIN_DAS)
      character*(*),     intent(IN) :: NAME
      integer,           intent(IN) :: COMP
      integer, optional, intent(IN) :: CONSTRAIN_DAS


      real,   pointer,     dimension(:,:,:) :: TENDSD   => null()
      real,   pointer,     dimension(:,:,:) :: TENDPH   => null()
      real,   pointer,     dimension(:,:,:) :: TENDBS   => null()
      real,   pointer,     dimension(:,:,:) :: TENDAN   => null()
      real,   pointer,     dimension(:,:,:) :: ANAINC   => null()
      real,   allocatable, dimension(:,:,:) :: TENDANAL
      real,   allocatable, dimension(:,:)   :: dummy
      real*8                                :: qave1
      real*8                                :: qave2
      real*8                                :: qave3
      integer                               :: L,LL,LU

      call MAPL_GetPointer(GIM(COMP), TENDSD, trim(NAME)        , rc=STATUS)
      VERIFY_(STATUS)
      
      select case (TYPE)
         case (FREERUN)

            TENDSD = 0.0

      case (PREDICTOR)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)
         if(associated(FC)) then
            TENDSD = TENDBS*FC
         else
            TENDSD = TENDBS
         end if

      case (FORECAST)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)

         TENDBS = BET * TENDBS
         if(associated(FC)) then
            TENDSD = TENDBS*FC
         else
            TENDSD = TENDBS
         end if

      case (CORRECTOR)

         LL = lbound(TENDSD,3)
         LU = ubound(TENDSD,3)

         allocate(TENDANAL(size(TENDSD,1),size(TENDSD,2),LL:LU), stat=STATUS)
         VERIFY_(STATUS)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT   , ANAINC, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)

         TENDANAL = ANAINC*(1.0/TAUANL)  ! No Constraints

         if( present(CONSTRAIN_DAS) ) then
                  if(CONSTRAIN_DAS == 1 ) then  ! ALPHA = SQRT( PS_ANA )
                  ! ----------------------------------------------------
                     allocate(dummy(size(TENDSD,1),size(TENDSD,2)), stat=STATUS)
                     VERIFY_(STATUS)
                     dummy = ANAINC(:,:,LU)                                     !                   ANAINC = PS_ANA - PS_BKG
                     call MAPL_AreaMean( qave1, dummy, area, grid, rc=STATUS )  ! qave1 = AreaMean( ANAINC )
                     VERIFY_(STATUS)
                     dummy = PLE(:,:,LU)                                        !                   P_n
                     call MAPL_AreaMean( qave2, dummy, area, grid, rc=STATUS )  ! qave2 = AreaMean( P_n    )
                     qave3 = qave2 + qave1*dt/TAUANL                            ! qave3 = AreaMean( P_n+1 = P_n + ANAINC*dt/tau )
                     VERIFY_(STATUS)
                     dummy = ANAINC(:,:,LU)*(qave2/qave3) - dummy * (qave1/qave3)
                     DO L=LL,LU
                     TENDANAL(:,:,L) = dummy(:,:)*BK(L)
                     ENDDO
                     TENDANAL = TENDANAL*(1.0/TAUANL)
                     deallocate(dummy)
                  endif

                  if(CONSTRAIN_DAS == 2 ) then  ! ALPHA = PS_ANA - PS_BKG
                  ! -----------------------------------------------------
                     allocate(dummy(size(TENDSD,1),size(TENDSD,2)), stat=STATUS)
                     VERIFY_(STATUS)
                     dummy = ANAINC(:,:,LU)
                     call MAPL_AreaMean( qave1, dummy, area, grid, rc=STATUS )
                     VERIFY_(STATUS)
                     dummy = ANAINC(:,:,LU)**2
                     call MAPL_AreaMean( qave2, dummy, area, grid, rc=STATUS )
                     VERIFY_(STATUS)
                     dummy = ANAINC(:,:,LU) * ( 1.0-ANAINC(:,:,LU)*(qave1/qave2) )
                     DO L=LL,LU
                     TENDANAL(:,:,L) = dummy(:,:)*BK(L)
                     ENDDO
                     TENDANAL = TENDANAL*(1.0/TAUANL)
                     deallocate(dummy)
                  endif

                  if(CONSTRAIN_DAS == 3 ) then  ! ALPHA = PS_ANA
                  ! --------------------------------------------
                     allocate(dummy(size(TENDSD,1),size(TENDSD,2)), stat=STATUS)
                     VERIFY_(STATUS)
                     dummy = ANAINC(:,:,LU)  ! ANAINC = PS_ANA - PS_BKG
                     call MAPL_AreaMean( qave1, dummy, area, grid, rc=STATUS )
                     VERIFY_(STATUS)
                     dummy = ( ANAINC(:,:,LU) + PLE(:,:,LU) )**2  ! PS_ANA**2
                     call MAPL_AreaMean( qave2, dummy, area, grid, rc=STATUS )
                     VERIFY_(STATUS)
                     dummy = ANAINC(:,:,LU) - dummy * (qave1/qave2)
                     DO L=LL,LU
                     TENDANAL(:,:,L) = dummy(:,:)*BK(L)
                     ENDDO
                     TENDANAL = TENDANAL*(1.0/TAUANL)
                     deallocate(dummy)
                  endif
         endif

         if(associated(FC)) then
            TENDSD = (TENDBS + TENDANAL)*FC
         else
            TENDSD = (TENDBS + TENDANAL)
         end if

         if (LAST_CORRECTOR) then
            TENDBS = BET*TENDBS + ALF*TENDANAL
         end if

         deallocate(TENDANAL)

      case default

         ASSERT_(.false.)

      end select

! Fill Total Increment Tendency (Current + Bias) Diagnostic
! ---------------------------------------------------------

      call MAPL_GetPointer ( EXPORT, TENDAN, trim(NAME)//'_ANA', rc=STATUS )
      VERIFY_(STATUS)

      if(associated(TENDAN)) then
          if(associated(FC)) then
                    TENDAN = TENDSD/FC
          else
                    TENDAN = TENDSD
          endif
      endif

    end subroutine DO_UPDATE_ANA3D

    subroutine DO_UPDATE_ANA2D(NAME, COMP)
      character*(*), intent(IN) :: NAME
      integer,       intent(IN) :: COMP

      real, pointer,     dimension(:,:)     :: TENDSD   => null()
      real, pointer,     dimension(:,:)     :: TENDPH   => null()
      real, pointer,     dimension(:,:)     :: TENDBS   => null()
      real, pointer,     dimension(:,:)     :: TENDAN   => null()
      real, pointer,     dimension(:,:)     :: ANAINC   => null()
      real, allocatable, dimension(:,:)     :: TENDANAL

      call MAPL_GetPointer(GIM(COMP), TENDSD, trim(NAME)        , rc=STATUS)
      VERIFY_(STATUS)
      
      select case (TYPE)
         case (FREERUN)

            TENDSD = 0.0

      case (PREDICTOR)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)
         TENDSD = TENDBS

      case (FORECAST)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)

         TENDBS = BET * TENDBS
         TENDSD = TENDBS

      case (CORRECTOR)

         allocate( TENDANAL(size(TENDSD,1),size(TENDSD,2)), stat=STATUS )
         VERIFY_(STATUS)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT   , ANAINC, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)

         TENDANAL = ANAINC*(1.0/TAUANL)
         TENDSD   = (TENDBS + TENDANAL)

         if (LAST_CORRECTOR) then
            TENDBS = BET*TENDBS + ALF*TENDANAL
         end if

         deallocate(TENDANAL)

      case default

         ASSERT_(.false.)

      end select

! Fill Total Increment Tendency (Current + Bias) Diagnostic
! ---------------------------------------------------------

      call MAPL_GetPointer ( EXPORT, TENDAN, trim(NAME)//'_ANA', rc=STATUS )
      VERIFY_(STATUS)

      if(associated(TENDAN)) then
                    TENDAN = TENDSD
      endif

    end subroutine DO_UPDATE_ANA2D


    subroutine DO_UPDATE_PHY (NAME)
      character*(*), intent(IN) :: NAME

      real, pointer,     dimension(:,:,:)     :: TENDSD   => null()
      real, pointer,     dimension(:,:,:)     :: TENDPH   => null()

      call MAPL_GetPointer(GIM(SDYN), TENDSD, trim(NAME)        , rc=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(PHYS), TENDPH, trim(NAME)        , rc=STATUS)
      VERIFY_(STATUS)
      
      TENDSD = TENDPH

    end subroutine DO_UPDATE_PHY


    subroutine DO_Friendly(Q, NAME)
      character*(*), intent(IN   ) :: NAME
      real,          intent(INOUT) :: Q(:,:,:)

      real, pointer,     dimension(:,:,:)     :: TENDBS   => null()
      real, pointer,     dimension(:,:,:)     :: TENDAN   => null()
      real, pointer,     dimension(:,:,:)     :: ANAINC   => null()
      real, allocatable, dimension(:,:,:)     :: TENDANAL
      real, allocatable, dimension(:,:,:)     :: QOLD

      allocate( QOLD(IM,JM,LM), stat=STATUS)
      VERIFY_(STATUS)

      QOLD = Q   ! Initialize Old Value for Total Tendency Diagnostic
      
      select case (TYPE)
      case (PREDICTOR)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)
         
         Q = Q + max( DTX*TENDBS*FC, -Q )  ! Prevent Negative Q

      case (FORECAST)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)
         
         TENDBS = BET * TENDBS

         Q = Q + max( DTX*TENDBS*FC, -Q )  ! Prevent Negative Q

      case (CORRECTOR)

         allocate(TENDANAL(IM,JM,LM), stat=STATUS)
         VERIFY_(STATUS)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT   , ANAINC, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)

         TENDANAL = ANAINC*(1.0/TAUANL)

         Q = Q + max( DTX*(TENDBS + TENDANAL)*FC, -Q )  ! Prevent Negative Q

         if (LAST_CORRECTOR) then
            TENDBS = BET*TENDBS + ALF*TENDANAL
         end if

         deallocate(TENDANAL)

      case default

         ASSERT_(.false.)

      end select

! Fill Total Increment Tendency (Current + Bias) Diagnostic
! ---------------------------------------------------------
      call MAPL_GetPointer ( EXPORT, TENDAN, trim(NAME)//'_ANA', rc=STATUS )
      VERIFY_(STATUS)
      if(associated(TENDAN)) TENDAN = (Q-QOLD)/DT

      deallocate(QOLD)

    end subroutine DO_FRIENDLY

    subroutine FILL_Friendly ( Q,DP,QFILL,QINT )
      real, intent(INOUT) ::   Q(:,:,:)
      real, intent(IN   ) ::  DP(:,:,:)
      real, intent(OUT  ) :: QFILL(:,:)
      real, intent(OUT  ) :: QINT (:,:)

      real*8, allocatable, dimension(:,:) :: QTEMP1
      real*8, allocatable, dimension(:,:) :: QTEMP2

      allocate(QTEMP1(IM,JM), stat=STATUS)
      VERIFY_(STATUS)
      allocate(QTEMP2(IM,JM), stat=STATUS)
      VERIFY_(STATUS)

      QTEMP1 = 0.0
      do L=1,LM
      QTEMP1(:,:) = QTEMP1(:,:) + Q(:,:,L)*DP(:,:,L)
      enddo 

      where( Q < 0.0 ) Q = 0.0

      QTEMP2 = 0.0
      do L=1,LM
      QTEMP2(:,:) = QTEMP2(:,:) + Q(:,:,L)*DP(:,:,L)
      enddo 

      where( qtemp2.ne.0.0_8 )
             qtemp2 = max( qtemp1/qtemp2, 0.0_8 )
      end where

      do L=1,LM
      Q(:,:,L) = Q(:,:,L)*qtemp2(:,:)
      enddo 

      QTEMP2 = 0.0
      do L=1,LM
      QTEMP2(:,:) = QTEMP2(:,:) + Q(:,:,L)*DP(:,:,L)
      enddo 

      WHERE( QTEMP1 >= 0.0 ) 
              QFILL  = 0.0
      ELSEWHERE
              QFILL = -QTEMP1 / (DT*MAPL_GRAV)
      END WHERE
      QINT  =  QTEMP2         /     MAPL_GRAV

      deallocate(QTEMP1)
      deallocate(QTEMP2)

    end subroutine FILL_FRIENDLY

  end subroutine Run

end module GEOS_AgcmGridCompMod
