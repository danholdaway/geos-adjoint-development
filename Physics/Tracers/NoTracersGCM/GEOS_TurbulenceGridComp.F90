!   $Id: GEOS_TurbulenceGridComp.F90,v 1.94.2.9.2.7.2.4.2.1.2.3.16.5 2013-06-07 17:30:09 ltakacs Exp $

#include "MAPL_Generic.h"

!=============================================================================

module GEOS_TurbulenceGridCompMod

!BOP

!  !MODULE: GEOS_Turbulence --- An GEOS generic atmospheric turbulence component

! !USES:

  use ESMF
  use GEOS_Mod
  use MAPL_Mod

#ifdef _CUDA
  use cudafor
  use LockEntrain, only: &
        ! Subroutines
        ENTRAIN, &
        ! Working Arrays
        ZFULL_DEV, TV_DEV, PV_DEV, RDZ_DEV, &
        DMI_DEV, PFULL_DEV,        &
        ! Inputs - Prelims
        T_DEV, QV_DEV, PHALF_DEV, TH_DEV,       &
        QLCN_DEV, QLLS_DEV, QICN_DEV, QILS_DEV, &
        ! Inputs - Louis
        U_DEV, V_DEV, ZPBL_DEV, &
        ! Inputs - Lock
        TDTLW_IN_DEV, U_STAR_DEV, B_STAR_DEV, FRLAND_DEV, &
        ! Inputs - Postlock
        CT_DEV, CQ_DEV, CU_DEV, &
        ! Inputs - Beljaars
        VARFLT_DEV, &
        ! Outputs - Prelims
        ZHALF_DEV, &
        ! Outputs - Louis
        DIFF_M_DEV, DIFF_T_DEV, &
        RI_DEV, DU_DEV,         &
        ! Outputs - Lock
        K_M_ENTR_DEV, K_T_ENTR_DEV,     &
        K_SFC_DIAG_DEV, K_RAD_DIAG_DEV, &
        ZCLOUD_DEV, ZRADML_DEV, &
        ZRADBASE_DEV, ZSML_DEV, &
        ! Outputs - Postlock
        AKQ_DEV, AKS_DEV, AKV_DEV, &
        BKQ_DEV, BKS_DEV, BKV_DEV, &
        CKQ_DEV, CKS_DEV, CKV_DEV, &
                          EKV_DEV, &
        ! Outputs - Beljaars
        FKV_DEV, &
        ! Outputs - Decomp
        DKQ_DEV, DKS_DEV, DKV_DEV, &
        ! Diagnostics - Louis
        ALH_DIAG_DEV, KMLS_DIAG_DEV, KHLS_DIAG_DEV, &
        ! Diagnostics - Lock
        ZCLDTOP_DIAG_DEV, WENTR_SFC_DIAG_DEV, WENTR_RAD_DIAG_DEV, &
        DEL_BUOY_DIAG_DEV, VSFC_DIAG_DEV, VRAD_DIAG_DEV,          &
        KENTRAD_DIAG_DEV, VBRV_DIAG_DEV, WENTR_BRV_DIAG_DEV,      &
        DSIEMS_DIAG_DEV, CHIS_DIAG_DEV, DELSINV_DIAG_DEV,         &
        SLMIXTURE_DIAG_DEV, CLDRADF_DIAG_DEV, RADRCODE_DIAG_DEV,  &
        ! Diagnostics - Postlock
        AKQODT_DIAG_DEV, AKSODT_DIAG_DEV, AKVODT_DIAG_DEV, &
        CKQODT_DIAG_DEV, CKSODT_DIAG_DEV, CKVODT_DIAG_DEV, &
        PPBL_DIAG_DEV,   TCZPBL_DIAG_DEV, &
        KPBL_DIAG_DEV,   &
        ZPBL2_DIAG_DEV, ZPBL10p_DIAG_DEV, ZPBLHTKE_DIAG_DEV, &
        ZPBLRI_DIAG_DEV, ZPBLTHV_DIAG_DEV, ZPBLRI2_DIAG_DEV, &
        ! Constants from MAPL_GetResource
        LOUIS_CONST, MINSHEAR_CONST, MINTHICK_CONST, AKHMMAX_CONST,     &
        LAMBDAM_CONST, LAMBDAM2_CONST, LAMBDAH_CONST, LAMBDAH2_CONST,   &
        ZKMENV_CONST, ZKHENV_CONST, PRANDTLSFC_CONST, PRANDTLRAD_CONST, &
        BETA_SURF_CONST, BETA_RAD_CONST, TPFAC_SFC_CONST,   &
        ENTRATE_SFC_CONST, PCEFF_SFC_CONST, KHRADFAC_CONST, &
        KHSFCFAC_CONST, LAMBDA_B_CONST, C_B_CONST, KPBLMIN_CONST
#else
  use LockEntrain, only: ENTRAIN
#endif
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !DESCRIPTION:
! 
!   {\tt GEOS\_TurbulenceGridComp} computes atmospheric tendencies due to turbulence.
!   Its physics is a combination of the first-order scheme of Louis---for stable PBLs
!   and free atmospheric turbulence---with a modified version of the non-local-K
!   scheme proposed by Lock for unstable and cloud-topped boundary layers.
!   In addition to diffusive tendencies, it adds the effects orographic form drag
!   for features with horizontal scales of 2 to 20 km following Beljaars et al. (2003,
!   ECMWF Tech. Memo. 427).  
!
!\vspace{12 pt}
!\noindent
!{\bf Grid Considerations}
!
!   Like all GEOS\_Generic-based components, it works on an inherited 
!   3-dimensional ESMF grid. It assumes that the first two (inner) dimensions span the
!   horizontal and the third (outer) dimension is the vertical. In the horizontal,
!   one or both dimensions can be degenerate, effectively supporting 
!   single-columns (1-D), and slices (2-D). No horizontal dimension needs to be
!   aligned with a particular coordinate. In the vertical, the only assumption
!   is that columns are indexed from top to bottom.
!
!\vspace{12 pt}
!\noindent
!{\bf Methods}
!
!   {\tt GEOS\_TurbulenceGridComp} uses the default Initialize and Finalize methods
!   of GEOS\_Generic. It has a 2-stage Run method that can be used in conjunction with
!   two-stage surface calculations to implement semi-implicit time differencing.
!
!\vspace{12 pt}
!\noindent
!{\bf Time Behavior}
!
!   {\tt GEOS\_TurbulenceGridComp} assumes both run stages will be invoked every 
!   RUN\_DT seconds, where RUN\_DT is required in the configuration. On this interval
!   both run stages will perform diffusion updates using diffusivities found in the
!   internal state.  The diffusivities in the internal state may be refreshed intermitently
!   by specifying MY\_STEP and ACCUMINT in the configuration. Accumulated imports used
!   in the intermittent refreshing are valid only on MY\_STEP intervals. Currently the
!   origin of these intervals is the beginning of the run. Accumulation of these imports
!   is done for a period ACCUMINT prior to the valid time. Both ACCUMINT and MY\_STEP are
!   in seconds.
!
!\vspace{12 pt}
!\noindent
!{\bf Working with Bundles and Friendlies}
!
!   {\tt GEOS\_TurbulenceGridComp} works on bundles of quantities to be diffused
!   and with corresponding bundles of their tendencies, surface values, etc.
!   These bundles may contain an arbitrary number of conservative quantities and
!   no requirements or restrictions are placed on what quantities they contain.
!   Quantities required for the calculation, such as pressures, stability, etc
!   are passed separately from the diffused quantities. Little distinction is made
!   of what is in the bundle, except that needed to decide what diffusivity applies
!   to the quantity and in what form its effects are implemented.
!
!   Quantities to be diffused can be marked as "Friendly-for-diffusion". In that case,
!   {\tt GEOS\_TurbulenceGridComp} directly updates the quantity; otherwise it 
!   merely computes its tendency, placing it in the appropriate bundle and treating
!   the quantity itself as read-only.
!
!   In working with bundled quantities, corresponding fields must appear in the 
!   same order in all bundles. Some of these fields, however, 
!   may be ``empty'' in the sense that the data pointer has not been allocated.
!   
!   {\tt GEOS\_TurbulenceGridComp} works with six bundles; three in the import
!   state and three in the export state. The import bundles are:
! \begin{itemize}
!   \item[]
!   \makebox[1in][l]{\bf TR} 
!   \parbox[t]{4in}{The quantity being diffused.}
!   \item[]
!   \makebox[1in][l]{\bf TRG} 
!   \parbox[t]{4in}{The surface (ground) value of the quantity being diffused.
!                   (Used only by Run2)}
!   \item[]
!   \makebox[1in][l]{\bf DTG} 
!   \parbox[t]{4in}{The change of TRG during the time step. (Used only by Run2)}
! \end{itemize}
!
!   The export bundles are:
! \begin{itemize}
!   \item[]
!   \makebox[1in][l]{\bf TRI} 
!   \parbox[t]{4in}{The tendency of the quantity being diffused.
!                   (Produced by Run1, updated by Run2.)  }
!   \item[]
!   \makebox[1in][l]{\bf FSTAR} 
!   \parbox[t]{4in}{After Run1, the ``preliminary'' (i.e., at the original surface
!    value) surface flux of the diffused quantity; after Run2, its final value.
!    (Produced by Run1, updated by Run2)}
!   \item[]
!   \makebox[1in][l]{\bf DFSTAR} 
!   \parbox[t]{4in}{The change of preliminary FSTAR per unit change in the 
!                   surface value. (Produced by Run1)}
! \end{itemize}
!
!   All fields in the export bundles are checked for associated pointers before being
!   updated.
!
!   Fields in the TR bundle can have four attributes:
! \begin{itemize}
! \item FriendlyTo[{\it Component Name}]: default=false --- If true, TR field is updated.
! \item WeightedTendency: default=true --- If true, tendencies (TRI) are pressure-weighted
! \item DiffuseLike: ('S','Q','M') default='S' --- Use mixing coefficients for either
!          heat, moisture or momentum.
! \end{itemize}
!  
!   Only fields in the TR bundle are checked for friendly status. Non-friendly
!   fields in TR and all other bundles are treated with the usual Import/Export
!   rules.
!
!\vspace{12 pt}
!\noindent
!{\bf Other imports and exports}
!
!   In addition to the updates of these bundles, {\tt GEOS\_TurbulenceGridComp} produces
!   a number of diagnostic exports, as well as frictional heating contributions. The latter 
!   are NOT added by {\tt GEOS\_TurbulenceGridComp}, but merely exported to be added
!   elsewhere in the GCM.
!
!\vspace{12 pt}
!\noindent
!{\bf Two-Stage Interactions with the Surface}
!
!   The two-stage scheme for interacting with the surface module is as follows:
! \begin{itemize}
!   \item  The first run stage takes the surface values of the diffused quantities
!      and the surface exchange coefficients as input. These are, of course, on the 
!      grid turbulence is working on.
!   \item  It then does the full diffusion calculation assuming the surface values are
!      fixed, i.e., the explicit surface case. In addition, it also computes derivatives of the
!      tendencies wrt surface values. These are to be used in the second stage.
!   \item The second run stage takes the increments of the surface values as inputs
!      and produces the final results, adding the implicit surface contributions. 
!   \item It also computes the frictional heating due to both implicit and explicit
!       surface contributions.
! \end{itemize}
!
!\vspace{12 pt}
!\noindent
!{\bf GEOS-5 Specific Aspects}
!
!   In GEOS-5, {\tt GEOS\_TurbulenceGridComp} works on the atmosphere's lat-lon grid,
!   while surface quantities are computed during the first run stage of the each of
!   the tiled surface components.  The tiled quantities are properly aggregated to
!   the GEOS-5 lat-lon grid by the first stage of {\tt GEOS\_SurfaceGridComp}, which
!   is called immediately before the first run stage of {\tt GEOS\_TurbulenceGridComp}.
!
!EOP

    logical                             :: dflt_false = .false.
    character(len=ESMF_MAXSTR)          :: dflt_q     = 'Q'
contains

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !DESCRIPTION: This version uses the {\tt GEOS\_GenericSetServices}, which sets
!               the Initialize and Finalize services to generic versions. It also
!   allocates our instance of a generic state and puts it in the 
!   gridded component (GC). Here we only set the two-stage run method and
!   declare the data services.
! \newline
! !REVISION HISTORY: 
!   ??Jul2006 E.Novak./Todling - Added output defining TLM/ADM trajectory

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code
!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Set the Run entry points
! ------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

!BOS

! !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'PLE',                                       &
        LONG_NAME  = 'air_pressure',                              &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                          &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

      call MAPL_AddImportSpec(GC,                                 &
         SHORT_NAME = 'T',                                        &
         LONG_NAME  = 'air_temperature',                          &
         UNITS      = 'K',                                        &
         DIMS       = MAPL_DimsHorzVert,                          &
         VLOCATION  = MAPL_VLocationCenter,                       &
         RESTART = .false.,                                       &
                                                        RC=STATUS  )
      VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'TH',                                        &
        LONG_NAME  = 'potential_temperature',                     &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'QV',                                        &
        LONG_NAME  = 'specific_humidity',                         &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'QLLS',                                      &
        LONG_NAME  = 'liquid_condensate_mixing_ratio',            &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'QILS',                                      &
        LONG_NAME  = 'frozen_condensate_mixing_ratio',            &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'CLLS',                                      &
        LONG_NAME  = 'cloud_fraction',                            &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'QLCN',                                      &
        LONG_NAME  = 'liquid_condensate_mixing_ratio',            &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'QICN',                                      &
        LONG_NAME  = 'frozen_condensate_mixing_ratio',            &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'CLCN',                                      &
        LONG_NAME  = 'cloud_fraction',                            &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'U',                                         &
        LONG_NAME  = 'eastward_wind',                             &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'V',                                         &
        LONG_NAME  = 'northward_wind',                            &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'CT',                                &
        LONG_NAME          = 'surface_heat_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'CQ',                                &
        LONG_NAME          = 'surface_moisture_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'CM',                                &
        LONG_NAME          = 'surface_momentum_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'BSTAR',                             &
        LONG_NAME          = 'surface_bouyancy_scale',            &
        UNITS              = 'm s-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'USTAR',                             &
        LONG_NAME          = 'surface_velocity_scale',            &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'FRLAND',                            &
        LONG_NAME          = 'land_fraction',                     &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'RADLW',                             &
        LONG_NAME          = 'air_temperature_tendency_due_to_longwave',&
        UNITS              = 'K s-1',                             &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'RADLWC',                            &
        LONG_NAME          = 'clearsky_air_temperature_tendency_lw',&
        UNITS              = 'K s-1',                             &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'PREF',                              &
        LONG_NAME          = 'reference_air_pressure',            &
        UNITS              = 'Pa',                                &
        DIMS               = MAPL_DimsVertOnly,                   &
        VLOCATION          = MAPL_VLocationEdge,                  &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VARFLT',                            &
        LONG_NAME          = 'variance_of_filtered_topography',   &
        UNITS              = 'm+2',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TR',                                &
        LONG_NAME          = 'diffused_quantities',               &
        UNITS              = 'X',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        DATATYPE           = MAPL_BundleItem,                     &
        RESTART = .false.,                                        &

                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TRG',                               &
        LONG_NAME          = 'surface_values_of_diffused_quantity',&
        UNITS              = 'X',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DATATYPE           = MAPL_BundleItem,                     &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'DTG',                               &
        LONG_NAME          = 'change_of_surface_values_of_diffused_quantity',&
        UNITS              = 'X',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DATATYPE           = MAPL_BundleItem,                     &
        RESTART = .false.,                                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

! !EXPORT STATE:

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME = 'TRI',                                       &
        LONG_NAME  = 'diffusion_tendencies',                      &
        UNITS      = 'X kg m-2 s-1',                              &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        DATATYPE   = MAPL_BundleItem,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME = 'FSTAR',                                     &
        LONG_NAME  = 'surface_fluxes',                            &
        UNITS      = 'X kg m-2 s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
        DATATYPE   = MAPL_BundleItem,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME = 'DFSTAR',                                    &
        LONG_NAME  = 'change_of_surface_fluxes_for_unit_change_of_surface_value',&
        UNITS      = 'kg m-2 s-1',                                &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
        DATATYPE   = MAPL_BundleItem,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'air_temperature',                                       &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'T',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'eastward_wind',                                         &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'U',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'northward_wind',                                        &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'V',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'specific_humidity',                                     &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'QV',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'total_momentum_diffusivity',                            &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KM',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'total_scalar_diffusivity',                              &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KH',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Richardson_number_from_Louis',                          &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'RI',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'bulk_shear_from_Louis',                                 &
       UNITS      = 's-1',                                                   &
       SHORT_NAME = 'DU',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'scalar_diffusivity_from_Louis',                         &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KHLS',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'momentum_diffusivity_from_Louis',                       &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KMLS',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_driven_scalar_diffusivity_from_Lock_scheme',    &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KHSFC',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'radiation_driven_scalar_diffusivity_from_Lock_scheme',  &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KHRAD',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'cloudy_LW_radiation_tendency_used_by_Lock_scheme',      &
       UNITS      = 'K s-1',                                                 &
       SHORT_NAME = 'LWCRT',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'entrainment_heat_diffusivity_from_Lock',                &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'EKH',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'entrainment_momentum_diffusivity_from_Lock',            &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'EKM',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Blackadar_length_scale_for_scalars',                    &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ALH',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'p-weighted_frictional_heating_rate_from_diffusion',     &
       UNITS      = 'K s-1 Pa',                                              &
       SHORT_NAME = 'INTDIS',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'p-weighted_frictional_heating_rate_from_orographic_drag',&
       UNITS      = 'K s-1 Pa',                                              &
       SHORT_NAME = 'TOPDIS',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'p-weighted_frictional_heating_rate_from_surface_drag',  &
       UNITS      = 'K s-1 Pa',                                              &
       SHORT_NAME = 'SRFDIS',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'KETRB',                                    &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_across_turbulence',&
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'KESRF',                                    &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_surface_friction',&
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'KEINT',                                    &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_diffusion',&
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'KETOP',                                    &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_topographic_friction',&
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'entrainment_velocity_from_surface_plume',               &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'WESFC',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'entrainment_velocity_from_radiation',                   &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'WERAD',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'entrainment_velocity_from_buoy_rev',                    &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'WEBRV',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Buoyancy_jump_across_inversion',                        &
       UNITS      = 'm s-2',                                                 &
       SHORT_NAME = 'DBUOY',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'turbulent_velocity_scale_for_sfc',                      &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'VSCSFC',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'turbulent_velocity_scale_for_cooling',                  &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'VSCRAD',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'turbulent_velocity_scale_for_buoy_rev',                 &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'VSCBRV',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'turbulent_entrainment_diff_from_cooling',               &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KERAD',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'cloud_top_radiative_forcing',                           &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'CLDRF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'pbltop_pressure',                                       &
       UNITS      = 'Pa',                                                    &
       SHORT_NAME = 'PPBL',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'pbltop_height_for_sfc_plume_LOCK',                      &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ZSML',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'depth_for_rad/brv_plume_LOCK',                          &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ZRADML',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'hght_of_base_for_rad/brv_plume_LOCK',                   &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ZRADBS',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'pbltop_cloud_depth_LOCK',                               &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ZCLD',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'pbltop_cloud_top_height_LOCK',                          &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ZCLDTOP',                                               &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'optimal_mixture_fraction_for_BRV',                      &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'CHIS',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 's_of_optimal_mixture_for_BRV',                          &
       UNITS      = 'J kg-1',                                                &
       SHORT_NAME = 'SMIXT',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Scaled_Del_s_at_Cloud_top',                             &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'DELSINV',                                               &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Siems_buoy_rev_parameter',                              &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DSIEMS',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Return_codes_for_Lock_top_driven_plume',                &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'RADRCODE',                                              &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ak_for_scalars_over_dt',                &
       SHORT_NAME = 'AKSODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ck_for_scalars_over_dt',                &
       SHORT_NAME = 'CKSODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ak_for_moisture_over_dt',               &
       SHORT_NAME = 'AKQODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ck_for_moisture_over_dt',               &
       SHORT_NAME = 'CKQODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ak_for_winds_over_dt',                  &
       SHORT_NAME = 'AKVODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ck_for_winds_over_dt',                  &
       SHORT_NAME = 'CKVODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'transcom_planetary_boundary_layer_height',              &
       SHORT_NAME = 'TCZPBL',                                                &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_threshold_2',           &
       SHORT_NAME = 'ZPBL2',                                                 &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_threshold_10p',         &
       SHORT_NAME = 'ZPBL10p',                                               &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_horiz_tke',             &
       SHORT_NAME = 'ZPBLHTKE',                                              &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_rich_0',                &
       SHORT_NAME = 'ZPBLRI',                                                &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_rich_02',               &
       SHORT_NAME = 'ZPBLRI2',                                               &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_thetav',                &
       SHORT_NAME = 'ZPBLTHV',                                               &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'pbltop_level',                                          &
       SHORT_NAME = 'KPBL',                                                  &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)



! !INTERNAL STATE:

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_ahat_for_scalars',                      &
       SHORT_NAME = 'AKS',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_bhat_for_scalars',                      &
       SHORT_NAME = 'BKS',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_c_for_scalars',                         &
       SHORT_NAME = 'CKS',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'sensitivity_of_tendency_to_surface_value_for_scalars',  &
       SHORT_NAME = 'DKS',                                                   &
       UNITS      = 's-1',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_ahat_for_moisture',                     &
       SHORT_NAME = 'AKQ',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_bhat_for_moisture',                     &
       SHORT_NAME = 'BKQ',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_c_for_moisture',                        &
       SHORT_NAME = 'CKQ',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'sensitivity_of_tendency_to_surface_value_for_moisture', &
       SHORT_NAME = 'DKQ',                                                   &
       UNITS      = 's-1',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_ahat_for_winds',                        &
       SHORT_NAME = 'AKV',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_bhat_for_winds',                        &
       SHORT_NAME = 'BKV',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_c_for_winds',                           &
       SHORT_NAME = 'CKV',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'sensitivity_of_tendency_to_surface_value_for_winds',    &
       SHORT_NAME = 'DKV',                                                   &
       UNITS      = 's-1',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'momentum_mixing_factor',                                &
       SHORT_NAME = 'EKV',                                                   &
       UNITS      = 'Pa s-1',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'topographic_roughness_factor',                          &
       SHORT_NAME = 'FKV',                                                   &
       UNITS      = 'Pa s-1',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'geopotential_height_above_surface',                     &
       SHORT_NAME = 'ZLE',                                                   &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'turbulence_tendency_for_dry_static_energy',             &
       SHORT_NAME = 'SINC',                                                  &
       UNITS      = 'm+2 s-3',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = .false.,                                                 &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'ZPBL',                                       &
       LONG_NAME  = 'planetary_boundary_layer_height',            &
       UNITS      = 'm',                                          &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzOnly,                           &
       VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)
!EOS

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,   name="-RUN1"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="--DIFFUSE"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="--REFRESHKS" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---REFRESHKS_RUN",RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---REFRESHKS_DATA",RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="----REFRESHKS_DATA_DEVICE",RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="----REFRESHKS_DATA_CONST",RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---REFRESHKS_ALLOC",RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---REFRESHKS_DEALLOC",RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---PRELIMS"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---LOUIS"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---LOCK"     ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---POSTLOCK" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---BELJAARS" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---DECOMP"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="-RUN2"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="--UPDATE"    ,RC=STATUS)
    VERIFY_(STATUS)
    
! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================


!BOP

! !IROUTINE: RUN1 -- First run stage for the {\tt MAPL_TurbulenceGridComp} component

! !INTERFACE:

  subroutine RUN1 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

! !DESCRIPTION: The first run stage of {\tt GEOS\_TurbulenceGridComp} computes the diffusivities,
!   sets-up the matrix for a backward-implicit computation of the surface fluxes,
!   and solves this system for a fixed surface value of the diffused quantity. Run1
!   takes as inputs the surface exchange coefficients (i.e., $\rho |U| C_{m,h,q}$) for
!   momentun, heat, and moisture, as well as the pressure, temperature, moisture, 
!   and winds for the sounding. These are used only for computing the diffusivities
!   and, as explained above,  are not the temperatures, moistures, etc. being diffused.
!
!   The computation of turbulence fluxes for fixed surface values is done at every
!   time step in the contained subroutine {\tt DIFFUSE}; but the computation of 
!   diffusivities and orographic drag coefficients, as well as the set-up of the
!   vertical difference matrix and its LU decomposition
!   can be done intermittently for economy in the contained subroutine  {\tt REFRESH}.
!   The results of this calculation are stored in an internal state. 
!   Run1 also computes the sensitivity of the 
!   atmospheric tendencies and the surface flux to changes in the surface value.
!
!   The diffusivities are computed by calls to {\tt LOUIS\_KS} and {\tt ENTRAIN}, which
!   compute the Louis et al. (1983) and Lock (2000) diffusivities. The Louis 
!   diffusivities are computed for all conditions, and {\tt ENTRAIN} overrides them 
!   where appropriate. Lock can be turned off from the resource file.


!

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp), pointer   :: MAPL
    type (ESMF_Config      )            :: CF
    type (ESMF_State       )            :: INTERNAL 
    type (ESMF_Alarm       )            :: ALARM   

! Local variables

    real, dimension(:,:,:), pointer     :: AKS, BKS, CKS, DKS
    real, dimension(:,:,:), pointer     :: AKQ, BKQ, CKQ, DKQ
    real, dimension(:,:,:), pointer     :: AKV, BKV, CKV, DKV, EKV, FKV
    real, dimension(:,:,:), pointer     :: PLE, ZLE, SINC
    real, dimension(:,:  ), pointer     :: CU, CT, CQ, ZPBL
    integer                             :: IM, JM, LM
    real                                :: DT

! Begin... 
!---------

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'Run1'

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"-RUN1")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,        &
         IM=IM, JM=JM, LM=LM,               &
         RUNALARM=ALARM,                    &
         INTERNAL_ESMF_STATE=INTERNAL,      &
                                  RC=STATUS )
    VERIFY_(STATUS)

! Get configuration from component
!---------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Get all pointers that are needed by both REFRESH and DIFFUSE
!-------------------------------------------------------------

! Get pressure structure; this is instantaneous.
!-----------------------------------------------

     call MAPL_GetPointer(IMPORT,  PLE,   'PLE',     RC=STATUS)
     VERIFY_(STATUS)

! Get surface exchange coefficients
!----------------------------------

     call MAPL_GetPointer(IMPORT,  CU,     'CM',     RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  CT,     'CT',     RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  CQ,     'CQ',     RC=STATUS)
     VERIFY_(STATUS)

! Get pointers from internal state
!---------------------------------

    call MAPL_GetPointer(INTERNAL, AKS,   'AKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, BKS,   'BKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CKS,   'CKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DKS,   'DKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, AKQ,   'AKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, BKQ,   'BKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CKQ,   'CKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DKQ,   'DKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, AKV,   'AKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, BKV,   'BKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CKV,   'CKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DKV,   'DKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, EKV,   'EKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, FKV,   'FKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, ZLE,   'ZLE',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, SINC,  'SINC',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, ZPBL,  'ZPBL',    RC=STATUS)
    VERIFY_(STATUS)

! Get application's timestep from configuration
!----------------------------------------------

    call ESMF_ConfigGetAttribute(CF, DT, Label="RUN_DT:" , RC=STATUS)
    VERIFY_(STATUS)

! If its time, do the refresh
! ---------------------------

    if ( ESMF_AlarmIsRinging(ALARM, rc=status) ) then
       VERIFY_(STATUS)
       call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
       VERIFY_(STATUS)

       call MAPL_TimerOn (MAPL,"--REFRESHKS")
        call REFRESH(IM,JM,LM,RC=STATUS)
        VERIFY_(STATUS)
       call MAPL_TimerOff(MAPL,"--REFRESHKS")
    endif

! Solve the free atmosphere problem
! ---------------------------------

    call MAPL_TimerOn (MAPL,"--DIFFUSE")
     call DIFFUSE(IM,JM,LM,RC=STATUS)
     VERIFY_(STATUS)
    call MAPL_TimerOff(MAPL,"--DIFFUSE")

!  All done with RUN1
!--------------------

    call MAPL_TimerOff(MAPL,"-RUN1")
    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  contains

!=============================================================================
!=============================================================================

!BOP

! !CROUTINE: REFRESH -- Refreshes diffusivities.

! !INTERFACE:

   subroutine REFRESH(IM,JM,LM,RC)

! !ARGUMENTS:

     integer,           intent(IN)       :: IM,JM,LM
     integer, optional, intent(OUT)      :: RC

! !DESCRIPTION: 
!   {\tt REFRESH} can be called intermittently to compute new values of the 
!   diffusivities. In addition it does all possible calculations that depend
!   only on these. In particular, it sets up the semi-implicit tridiagonal
!   solver in the vertical and does the LU decomposition. It also includes the
!   local effects of orographic drag, so that it to is done implicitly.
!
!   Diffusivities are first computed with the Louis scheme ({\tt LOUIS\_KS}),
!   and then, where appropriate,
!   they are overridden by the Lock values ({\tt ENTRAIN}).
!   Once diffusivities are computed, {\tt REFRESH} sets-up the tridiagonal
!   matrices for the semi-implicit vertical diffusion calculation and performs
!   their $LU$ decomposition. 
!
!   {\tt REFRESH} requires surface exchange coefficients for heat, moisture, and
!   momentum,  The calculations in the interior are also
!   done for momentum, heat, and water diffusion. Heat and water mixing
!   coefficients differ only at the surface, but these affect the entire $LU$
!   decomposition, and so all three decompositions are saved in the internal state. 
!
!   For a conservatively diffused quantity $q$, we have
!   $$
!   \frac{\partial q}{\partial t} = -g \frac{\partial }{\partial p} 
!       \left(\rho K_q \frac{\partial q}{\partial z} \right)
!   $$
!   In finite difference form, using backward time differencing, this becomes
!   $$
!   \begin{array}{rcl}
!   {q^{n+1}_l - q^{n}_l} & = & - \frac{g}{\delta_l p}^*
!     \delta_l \left[
!      \left( \frac{\Delta t \rho K_q}{\delta_l z} \right)^* (\delta_l q)^{n+1} \right]   \\
!   &&\\
!                         & = & - \alpha_l ( \beta_{l+\frac{1}{2}}(q_{l+1}-q_l)^{n+1} - 
!                                            \beta_{l-\frac{1}{2}}(q_l-q_{l-1})^{n+1} ) \\
!   &&\\
!   \alpha_l & = & \frac{g \Delta t}{(p_{l+\frac{1}{2}}-p_{l-\frac{1}{2}})^*} \\
!   &&\\
!   \beta_{l+\frac{1}{2}} & = & \left( \frac{ (\rho K_q)^*_{l+\frac{1}{2}}}{(z_{l+1}-z_{l})^*} \right) \\
!   \end{array}
!   $$
!   where the subscripts denote levels, superscripts denote times, and the $*$ superscript
!   denotes evaluation at the refresh time.
!   The following tridiagonal set is then solved for $q^{n+1}_l$:
!   $$
!   a_l q_{l-1} + b_l q_l + c_l q_{l+1} = q_l
!   $$
!   where
!   $$
!   \begin{array}{rcl}
!   a_l & = & \alpha_l \beta_{l-\frac{1}{2}} \\
!   c_l & = & \alpha_l \beta_{l+\frac{1}{2}} \\
!   b_l & = & 1 - a_l - c_l.
!   \end{array}
!   $$
!   At the top boundary, we assume $K_q=0$, so  $ \beta_{\frac{1}{2}}=0$ and $a_1=0$.
!   At the surface, $ \beta_{L+\frac{1}{2}}= \rho_s |U|_s C_{m,h,q}$, the surface exchange coefficient.
!   

!EOP
 
     character(len=ESMF_MAXSTR)          :: IAm='Refresh'
     integer                             :: STATUS

     real, dimension(:,:,:), pointer     :: TH, U, V, Q, T, RI, DU, RADLW, RADLWC, LWCRT
     real, dimension(:,:  ), pointer     :: VARFLT
     real, dimension(:,:,:), pointer     :: KH, KM, QLLS, QILS, CLLS, QLCN, QICN, CLCN
     real, dimension(:,:,:), pointer     :: ALH
     real, dimension(:    ), pointer     :: PREF

     real, dimension(IM,JM,1:LM-1)       :: TVE, RDZ
     real, dimension(IM,JM,LM)           :: THV, TV, Z, DMI, PLO, QL, QI, QA
     real, dimension(IM,JM,0:LM)         :: PKE

     real, dimension(:,:,:), pointer     :: EKH, EKM, KHLS, KMLS, KHRAD, KHSFC
     real, dimension(:,:  ), pointer     :: BSTAR, USTAR, PPBL, WERAD, WESFC,VSCRAD,KERAD,DBUOY,ZSML,ZCLD,ZRADML,FRLAND
     real, dimension(:,:  ), pointer     :: TCZPBL => null()
     real, dimension(:,:  ), pointer     :: ZPBL2 => null()
     real, dimension(:,:  ), pointer     :: ZPBL10P => null()
     real, dimension(:,:  ), pointer     :: ZPBLHTKE => null()
     real, dimension(:,:  ), pointer     :: ZPBLRI => null()
     real, dimension(:,:  ), pointer     :: ZPBLRI2 => null()
     real, dimension(:,:  ), pointer     :: ZPBLTHV => null()
     real, dimension(:,:  ), pointer     :: KPBL
     real, dimension(:,:  ), pointer     :: WEBRV,VSCBRV,DSIEMS,CHIS,ZCLDTOP,DELSINV,SMIXT,ZRADBS,CLDRF,VSCSFC,RADRCODE

     real, dimension(:,:,:), pointer     :: AKSODT, CKSODT
     real, dimension(:,:,:), pointer     :: AKQODT, CKQODT
     real, dimension(:,:,:), pointer     :: AKVODT, CKVODT


     logical, dimension(IM,JM     )      :: CONVECT
     logical                             :: ALLOC_TCZPBL, CALC_TCZPBL
     logical                             :: ALLOC_ZPBL2, CALC_ZPBL2
     logical                             :: ALLOC_ZPBL10p, CALC_ZPBL10p

     real                                :: LOUIS
     real                                :: LAMBDAM, LAMBDAM2
     real                                :: LAMBDAH, LAMBDAH2
     real                                :: ZKMENV, ZKHENV 
     real                                :: MINTHICK
     real                                :: MINSHEAR
     real                                :: AKHMMAX
     real                                :: C_B, LAMBDA_B, ZMAX_B,LOUIS_MEMORY
     real                                :: PRANDTLSFC,PRANDTLRAD,BETA_RAD,BETA_SURF,KHRADFAC,TPFAC_SURF,ENTRATE_SURF
     real                                :: PCEFF_SURF, KHSFCFAC,ZCHOKE

     integer                             :: I,J,L,LOCK_ON
     integer                             :: KPBLMIN,PBLHT_OPTION

     ! Below are automatic arrays corresponding to diagnostic 
     ! export pointers, e.g.:
     !    WESFC_X is the automatic for WESFC
     !
     ! To add a new diagnostic, add here and then use
     ! IF(ASSOCIATED(GGG)) GGG = GGG_X

     REAL, DIMENSION(IM,JM) :: ZCLDTOP_X
     REAL, DIMENSION(IM,JM) :: WESFC_X
     REAL, DIMENSION(IM,JM) :: WERAD_X
     REAL, DIMENSION(IM,JM) :: DBUOY_X
     REAL, DIMENSION(IM,JM) :: VSCRAD_X
     REAL, DIMENSION(IM,JM) :: VSCSFC_X
     REAL, DIMENSION(IM,JM) :: KERAD_X
     REAL, DIMENSION(IM,JM) :: VSCBRV_X
     REAL, DIMENSION(IM,JM) :: WEBRV_X
     REAL, DIMENSION(IM,JM) :: DSIEMS_X
     REAL, DIMENSION(IM,JM) :: CHIS_X
     REAL, DIMENSION(IM,JM) :: DELSINV_X
     REAL, DIMENSION(IM,JM) :: SMIXT_X
     REAL, DIMENSION(IM,JM) :: CLDRF_X
     REAL, DIMENSION(IM,JM) :: RADRCODE_X

     REAL, DIMENSION(IM,JM,0:LM) :: ALH_X
     REAL, DIMENSION(IM,JM,0:LM) :: KMLS_X
     REAL, DIMENSION(IM,JM,0:LM) :: KHLS_X

     REAL, DIMENSION(IM,JM,LM) :: AKQODT_X
     REAL, DIMENSION(IM,JM,LM) :: AKSODT_X
     REAL, DIMENSION(IM,JM,LM) :: AKVODT_X
     REAL, DIMENSION(IM,JM,LM) :: CKQODT_X
     REAL, DIMENSION(IM,JM,LM) :: CKSODT_X
     REAL, DIMENSION(IM,JM,LM) :: CKVODT_X
     REAL, DIMENSION(IM,JM   ) :: PPBL_X
     REAL, DIMENSION(IM,JM   ) :: KPBL_X
     REAL, DIMENSION(IM,JM   ) :: ZPBLHTKE_X
     REAL, DIMENSION(IM,JM   ) :: ZPBLRI_X
     REAL, DIMENSION(IM,JM   ) :: ZPBLRI2_X
     REAL, DIMENSION(IM,JM   ) :: ZPBLTHV_X
  
#ifdef _CUDA
     type(dim3) :: Grid, Block
     integer :: blocksize
#endif

     call MAPL_TimerOn(MAPL,"---PRELIMS")
     
! Get Sounding from the import state
!-----------------------------------

     call MAPL_GetPointer(IMPORT,     T,       'T', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,     Q,      'QV', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,    TH,      'TH', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,     U,       'U', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,     V,       'V', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,VARFLT,  'VARFLT', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  PREF,    'PREF', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT, RADLW,   'RADLW', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,RADLWC,  'RADLWC', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  QLLS,    'QLLS', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  QILS,    'QILS', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  QLCN,    'QLCN', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  QICN,    'QICN', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  CLLS,    'CLLS', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  CLCN,    'CLCN', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT, BSTAR,   'BSTAR', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT, USTAR,   'USTAR', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,FRLAND,  'FRLAND', RC=STATUS); VERIFY_(STATUS)

! Get turbulence parameters from configuration
!---------------------------------------------

     call MAPL_GetResource (MAPL, LOUIS,        trim(COMP_NAME)//"_LOUIS:",        default=5.0,          RC=STATUS)
     call MAPL_GetResource (MAPL, LAMBDAM,      trim(COMP_NAME)//"_LAMBDAM:",      default=160.0,        RC=STATUS)
     call MAPL_GetResource (MAPL, LAMBDAM2,     trim(COMP_NAME)//"_LAMBDAM2:",     default=1.0,          RC=STATUS)
     call MAPL_GetResource (MAPL, LAMBDAH,      trim(COMP_NAME)//"_LAMBDAH:",      default=160.0,        RC=STATUS)
     call MAPL_GetResource (MAPL, LAMBDAH2,     trim(COMP_NAME)//"_LAMBDAH2:",     default=1.0,          RC=STATUS)
     call MAPL_GetResource (MAPL, ZKMENV,       trim(COMP_NAME)//"_ZKMENV:",       default=3000.,        RC=STATUS)
     call MAPL_GetResource (MAPL, ZKHENV,       trim(COMP_NAME)//"_ZKHENV:",       default=3000.,        RC=STATUS)
     call MAPL_GetResource (MAPL, MINTHICK,     trim(COMP_NAME)//"_MINTHICK:",     default=0.1,          RC=STATUS)
     call MAPL_GetResource (MAPL, MINSHEAR,     trim(COMP_NAME)//"_MINSHEAR:",     default=0.0030,       RC=STATUS)
     call MAPL_GetResource (MAPL, C_B,          trim(COMP_NAME)//"_C_B:",          default=2.5101471e-8, RC=STATUS)
     call MAPL_GetResource (MAPL, LAMBDA_B,     trim(COMP_NAME)//"_LAMBDA_B:",     default=1500.,        RC=STATUS)
     call MAPL_GetResource (MAPL, AKHMMAX,      trim(COMP_NAME)//"_AKHMMAX:",      default=500.,         RC=STATUS)
     call MAPL_GetResource (MAPL, LOCK_ON,      trim(COMP_NAME)//"_LOCK_ON:",      default=1,            RC=STATUS)
     call MAPL_GetResource (MAPL, PRANDTLSFC,   trim(COMP_NAME)//"_PRANDTLSFC:",   default=1.0,          RC=STATUS)
     call MAPL_GetResource (MAPL, PRANDTLRAD,   trim(COMP_NAME)//"_PRANDTLRAD:",   default=0.75,         RC=STATUS)
     call MAPL_GetResource (MAPL, BETA_RAD,     trim(COMP_NAME)//"_BETA_RAD:",     default=0.50,         RC=STATUS)
     call MAPL_GetResource (MAPL, BETA_SURF,    trim(COMP_NAME)//"_BETA_SURF:",    default=0.25,         RC=STATUS)
     call MAPL_GetResource (MAPL, KHRADFAC,     trim(COMP_NAME)//"_KHRADFAC:",     default=0.85,         RC=STATUS)
     call MAPL_GetResource (MAPL, KHSFCFAC,     trim(COMP_NAME)//"_KHSFCFAC:",     default=0.45,         RC=STATUS)
     call MAPL_GetResource (MAPL, TPFAC_SURF,   trim(COMP_NAME)//"_TPFAC_SURF:",   default=20.0,         RC=STATUS)
     call MAPL_GetResource (MAPL, ENTRATE_SURF, trim(COMP_NAME)//"_ENTRATE_SURF:", default=1.5e-3,       RC=STATUS)
     call MAPL_GetResource (MAPL, PCEFF_SURF,   trim(COMP_NAME)//"_PCEFF_SURF:",   default=0.5,          RC=STATUS)
     call MAPL_GetResource (MAPL, LOUIS_MEMORY, trim(COMP_NAME)//"_LOUIS_MEMORY:", default=-999.,        RC=STATUS)
     call MAPL_GetResource (MAPL, PBLHT_OPTION, trim(COMP_NAME)//"_PBLHT_OPTION:", default=1,            RC=STATUS)


! Get pointers from export state...
!-----------------------------------

     call MAPL_GetPointer(EXPORT,      KH,      'KH', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,      KM,      'KM', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,      RI,      'RI', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,      DU,      'DU', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,     EKH,     'EKH', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,     EKM,     'EKM', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    KHLS,    'KHLS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    KMLS,    'KMLS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   KHSFC,   'KHSFC', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   KHRAD,   'KHRAD', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    PPBL,    'PPBL',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    KPBL,    'KPBL',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    TCZPBL,  'TCZPBL',             RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBL2,  'ZPBL2',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBL10p,  'ZPBL10p',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBLHTKE,  'ZPBLHTKE',         RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBLRI,  'ZPBLRI',             RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBLRI2,  'ZPBLRI2',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBLTHV,  'ZPBLTHV',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   LWCRT,   'LWCRT', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   WERAD,   'WERAD',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   WESFC,   'WESFC',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   DBUOY,   'DBUOY',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  VSCRAD,  'VSCRAD',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  VSCsfc,  'VSCSFC',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   KERAD,   'KERAD',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  VSCBRV,  'VSCBRV',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   WEBRV,   'WEBRV',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    CHIS,    'CHIS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  DSIEMS,  'DSIEMS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZCLD,    'ZCLD', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZSML,    'ZSML', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  ZRADML,  'ZRADML', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  ZRADBS,  'ZRADBS', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT, ZCLDTOP, 'ZCLDTOP',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT, DELSINV, 'DELSINV',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,RADRCODE,'RADRCODE',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   SMIXT,   'SMIXT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   CLDRF,   'CLDRF',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,     ALH,     'ALH',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  AKSODT,  'AKSODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  CKSODT,  'CKSODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  AKQODT,  'AKQODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  CKQODT,  'CKQODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  AKVODT,  'AKVODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  CKVODT,  'CKVODT',               RC=STATUS)
     VERIFY_(STATUS)

! Initialize some arrays

      LWCRT = RADLW - RADLWC

      KH    = 0.0
      KM    = 0.0
      RI    = 0.0
      DU    = 0.0
      EKH   = 0.0
      EKM   = 0.0
      KHSFC = 0.0
      KHRAD = 0.0
      if(associated(KHLS)) KHLS = 0.0
      if(associated(KMLS)) KMLS = 0.0

      call MAPL_TimerOff(MAPL,"---PRELIMS")

      ! MATMAT
      ! Note: There appears to be a bug in Intel 13 wherein we
      !       cannot use something like:
      !          if(ALLOC_TCZPBL.AND.PBLHT_OPTION==3) then
      !       in the tests below for allocating. This seems
      !       to be due to sending a pointer that hasn't
      !       been allocated to POSTLOCK which has an 
      !       explicit-shape interface to, say, TCZPBL_DIAG.
      !       My guess is this is a bug with the optimizer 
      !       in Intel 13 as both Intel 11 and PGI 12 allow
      !       this code and Intel 13 does not flag it as
      !       against the Fortran 2003 standard (-std). 
      !
      !       If fixed in the future, the more complex test
      !       should be restored.

      ALLOC_ZPBL2 = .FALSE.
      CALC_ZPBL2 = .FALSE.
      if(associated(ZPBL2).OR.PBLHT_OPTION==1) CALC_ZPBL2 = .TRUE.
      if(.not.associated(ZPBL2 )) then
         allocate(ZPBL2(IM,JM))
         ALLOC_ZPBL2 = .TRUE.
      endif

      ALLOC_ZPBL10p = .FALSE.
      CALC_ZPBL10p = .FALSE.
      if(associated(ZPBL10p).OR.PBLHT_OPTION==2) CALC_ZPBL10p = .TRUE.
      if(.not.associated(ZPBL10p )) then
         allocate(ZPBL10p(IM,JM))
         ALLOC_ZPBL10p = .TRUE.
      endif

      ALLOC_TCZPBL = .FALSE.
      CALC_TCZPBL = .FALSE.
      if(associated(TCZPBL).OR.PBLHT_OPTION==3) CALC_TCZPBL = .TRUE.
      if(.not.associated(TCZPBL)) then
                allocate(TCZPBL(IM,JM))
                   ALLOC_TCZPBL = .TRUE.
      endif

#ifdef _CUDA

      ASSERT_(LM <= GPU_MAXLEVS) !If this is tripped, GNUmakefile
                                 !must be changed

      call MAPL_GetResource(MAPL,BLOCKSIZE,'BLOCKSIZE:',DEFAULT=128,RC=STATUS)
      VERIFY_(STATUS)

      Block = dim3(blocksize,1,1)
      Grid = dim3(ceiling(real(IM*JM)/real(blocksize)),1,1)

      call MAPL_TimerOn (MAPL,name="---REFRESHKS_ALLOC" ,RC=STATUS)
      VERIFY_(STATUS)

      ! ----------------------
      ! Allocate device arrays
      ! ----------------------

      ! Working Arrays
      ! --------------

      ALLOCATE(ZFULL_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(TV_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PV_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(RDZ_DEV(IM*JM,(LM-1)), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(DMI_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PFULL_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      !ALLOCATE(QA_DEV(IM*JM,LM), STAT=STATUS) ! NOT CURRENTLY USED
      !VERIFY_(STATUS)

      ! Inputs - Prelims
      ! ----------------

      ALLOCATE(T_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QV_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PHALF_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(TH_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QLCN_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QLLS_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QICN_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(QILS_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      !ALLOCATE(CLCN_DEV(IM*JM,LM), STAT=STATUS) ! NOT CURRENTLY USED
      !VERIFY_(STATUS)
      !ALLOCATE(CLLS_DEV(IM*JM,LM), STAT=STATUS) ! NOT CURRENTLY USED
      !VERIFY_(STATUS)

      ! Inputs - Louis
      ! --------------

      ALLOCATE(U_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(V_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(ZPBL_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)

      ! Inputs - Lock
      ! -------------
   
      ALLOCATE(TDTLW_IN_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(U_STAR_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(B_STAR_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(FRLAND_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)

      ! Inputs - Postlock
      ! -----------------
   
      ALLOCATE(CT_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CQ_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CU_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
   
      ! Inputs - Beljaars
      ! -----------------

      ALLOCATE(VARFLT_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)

      ! Outputs - Prelims
      ! -----------------

      ALLOCATE(ZHALF_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)

      ! Outputs - Louis
      ! ---------------

      ALLOCATE(DIFF_M_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(DIFF_T_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(RI_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(DU_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)

      ! Outputs - Lock
      ! --------------

      ALLOCATE(K_M_ENTR_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(K_T_ENTR_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(K_SFC_DIAG_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(K_RAD_DIAG_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(ZCLOUD_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(ZRADML_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(ZRADBASE_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(ZSML_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)

      ! Outputs - Postlock
      ! ------------------

      ALLOCATE(AKQ_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(AKS_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(AKV_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(BKQ_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(BKS_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(BKV_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CKQ_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CKS_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CKV_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(EKV_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)

      ! Outputs - Beljaars
      ! ------------------

      ALLOCATE(FKV_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)

      ! Outputs - Decomp
      ! ----------------

      ALLOCATE(DKQ_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(DKS_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(DKV_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)

      ! Diagnostics - Louis
      ! -------------------

      ALLOCATE(ALH_DIAG_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(KMLS_DIAG_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(KHLS_DIAG_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
  
      ! Diagnostics - Lock
      ! ------------------

      ALLOCATE(ZCLDTOP_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(WENTR_SFC_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(WENTR_RAD_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(DEL_BUOY_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(VSFC_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(VRAD_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(KENTRAD_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(VBRV_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(WENTR_BRV_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(DSIEMS_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CHIS_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(DELSINV_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(SLMIXTURE_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CLDRADF_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(RADRCODE_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)

      ! Diagnostics - Postlock
      ! ----------------------
   
      ALLOCATE(AKQODT_DIAG_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(AKSODT_DIAG_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(AKVODT_DIAG_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CKQODT_DIAG_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CKSODT_DIAG_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CKVODT_DIAG_DEV(IM*JM,LM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(PPBL_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(KPBL_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      if(CALC_TCZPBL) then
         ALLOCATE(TCZPBL_DIAG_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
      endif
      if(CALC_ZPBL2) then
         ALLOCATE(ZPBL2_DIAG_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
      endif
      if(CALC_ZPBL10p) then
         ALLOCATE(ZPBL10p_DIAG_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
      endif
      ALLOCATE(ZPBLHTKE_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(ZPBLRI_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(ZPBLRI2_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(ZPBLTHV_DIAG_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOff(MAPL,name="---REFRESHKS_ALLOC" ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (MAPL,name="---REFRESHKS_DATA" ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (MAPL,name="----REFRESHKS_DATA_DEVICE" ,RC=STATUS)
      VERIFY_(STATUS)

      ! ---------------------
      ! Copy inputs to device
      ! ---------------------

      ! Inputs
      ! ------

      STATUS = cudaMemcpy(T_DEV,T,IM*JM*LM)
      STATUS = cudaMemcpy(QV_DEV,Q,IM*JM*LM)
      STATUS = cudaMemcpy(PHALF_DEV,PLE,IM*JM*(LM+1))
      STATUS = cudaMemcpy(TH_DEV,TH,IM*JM*LM)
      STATUS = cudaMemcpy(QLCN_DEV,QLCN,IM*JM*LM)
      STATUS = cudaMemcpy(QLLS_DEV,QLLS,IM*JM*LM)
      STATUS = cudaMemcpy(QICN_DEV,QICN,IM*JM*LM)
      STATUS = cudaMemcpy(QILS_DEV,QILS,IM*JM*LM)
      !STATUS = cudaMemcpy(CLCN_DEV,CLCN,IM*JM*LM) ! NOT CURRENTLY USED
      !STATUS = cudaMemcpy(CLLS_DEV,CLLS,IM*JM*LM) ! NOT CURRENTLY USED

      STATUS = cudaMemcpy(U_DEV,U,IM*JM*LM)
      STATUS = cudaMemcpy(V_DEV,V,IM*JM*LM)
      STATUS = cudaMemcpy(ZPBL_DEV,ZPBL,IM*JM)

      STATUS = cudaMemcpy(TDTLW_IN_DEV,RADLW,IM*JM*LM)
      STATUS = cudaMemcpy(U_STAR_DEV,USTAR,IM*JM)
      STATUS = cudaMemcpy(B_STAR_DEV,BSTAR,IM*JM)
      STATUS = cudaMemcpy(FRLAND_DEV,FRLAND,IM*JM)

      STATUS = cudaMemcpy(CT_DEV,CT,IM*JM)
      STATUS = cudaMemcpy(CQ_DEV,CQ,IM*JM)
      STATUS = cudaMemcpy(CU_DEV,CU,IM*JM)

      STATUS = cudaMemcpy(VARFLT_DEV,VARFLT,IM*JM)

      call MAPL_TimerOff(MAPL,name="----REFRESHKS_DATA_DEVICE" ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (MAPL,name="----REFRESHKS_DATA_CONST" ,RC=STATUS)
      VERIFY_(STATUS)

      ! ------------------------
      ! Load constants on device
      ! ------------------------

      ! Constants from MAPL_GetResource
      ! -------------------------------

      LOUIS_CONST       = LOUIS
      MINSHEAR_CONST    = MINSHEAR
      MINTHICK_CONST    = MINTHICK
      AKHMMAX_CONST     = AKHMMAX
      LAMBDAM_CONST     = LAMBDAM
      LAMBDAM2_CONST    = LAMBDAM2
      LAMBDAH_CONST     = LAMBDAH
      LAMBDAH2_CONST    = LAMBDAH2
      ZKMENV_CONST      = ZKMENV
      ZKHENV_CONST      = ZKHENV
      PRANDTLSFC_CONST  = PRANDTLSFC
      PRANDTLRAD_CONST  = PRANDTLRAD
      BETA_SURF_CONST   = BETA_SURF
      BETA_RAD_CONST    = BETA_RAD
      TPFAC_SFC_CONST   = TPFAC_SURF
      ENTRATE_SFC_CONST = ENTRATE_SURF
      PCEFF_SFC_CONST   = PCEFF_SURF
      KHRADFAC_CONST    = KHRADFAC
      KHSFCFAC_CONST    = KHSFCFAC
      LAMBDA_B_CONST    = LAMBDA_B
      C_B_CONST         = C_B

      KPBLMIN_CONST     = count(PREF < 50000.)

      STATUS = cudaDeviceSynchronize()

      call MAPL_TimerOff(MAPL,name="----REFRESHKS_DATA_CONST" ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOff(MAPL,name="---REFRESHKS_DATA" ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (MAPL,name="---REFRESHKS_RUN" ,RC=STATUS)
      VERIFY_(STATUS)

      STATUS = cudaFuncSetCacheConfig('geos_turbulencegridcompmod_prelims',cudaFuncCachePreferL1)

      call PRELIMS<<<Grid,Block>>>(IM*JM, LM, DT, &
                                   ! Inputs
                                   T_DEV,     &
                                   QV_DEV,    &
                                   PHALF_DEV, &
                                   TH_DEV,    &
                                   QLCN_DEV,  &
                                   QLLS_DEV,  &
                                   QICN_DEV,  &
                                   QILS_DEV,  &
                                   !CLCN_DEV, & ! CURRENTLY NOT USED
                                   !CLLS_DEV, & ! CURRENTLY NOT USED
                                   ! Outputs
                                   ZFULL_DEV, &
                                   ZHALF_DEV, &
                                   TV_DEV,    &
                                   PV_DEV,    &
                                   RDZ_DEV,   &
                                   DMI_DEV,   &
                                   !QA_DEV,   & ! CURRENTLY NOT USED
                                   PFULL_DEV  )


      STATUS = cudaGetLastError()
      if (STATUS /= 0) then 
         write (*,*) "Error code from PRELIMS kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         ASSERT_(.FALSE.)
      end if

      STATUS = cudaDeviceSynchronize()

      !   Refresh diffusivities: First compute Louis...
      !   ---------------------------------------------

      STATUS = cudaFuncSetCacheConfig('geos_turbulencegridcompmod_louis_ks',cudaFuncCachePreferL1)

      call LOUIS_KS<<<Grid,Block>>>(IM*JM, LM,            &
            DIFF_T_DEV,DIFF_M_DEV,RI_DEV,DU_DEV,ZPBL_DEV, &    
            ZFULL_DEV,ZHALF_DEV,PV_DEV,U_DEV,V_DEV,       &    
            LOUIS_CONST, MINSHEAR_CONST, MINTHICK_CONST,  &
            LAMBDAM_CONST, LAMBDAM2_CONST, &
            LAMBDAH_CONST, LAMBDAH2_CONST, & 
            ZKMENV_CONST, ZKHENV_CONST,    &
            AKHMMAX_CONST,                 &
            ALH_DIAG_DEV,KMLS_DIAG_DEV,KHLS_DIAG_DEV)

      STATUS = cudaGetLastError()
      if (STATUS /= 0) then 
         write (*,*) "Error code from LOUIS kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         ASSERT_(.FALSE.)
      end if

      STATUS = cudaDeviceSynchronize()

      !   ...then add Lock.
      !--------------------

      if (LOCK_ON==1) then

         STATUS = cudaFuncSetCacheConfig('lockentrain_entrain',cudaFuncCachePreferL1)

         call ENTRAIN<<<Grid,Block>>>(IM*JM,LM)

         STATUS = cudaGetLastError()
         if (STATUS /= 0) then 
            write (*,*) "Error code from ENTRAIN kernel call: ", STATUS
            write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
            ASSERT_(.FALSE.)
         end if

         STATUS = cudaDeviceSynchronize()

      endif

      STATUS = cudaFuncSetCacheConfig('geos_turbulencegridcompmod_postlock',cudaFuncCachePreferL1)

      call POSTLOCK<<<Grid,Block>>>(IM*JM, LM, DT, &
                                    ! Constants 
                                    KPBLMIN_CONST, &
                                    CALC_TCZPBL  , &
                                    CALC_ZPBL2  , &
                                    CALC_ZPBL10p  , &
                                    ! Inputs
                                    ZFULL_DEV, ZHALF_DEV, &
                                    PFULL_DEV, &
                                    RDZ_DEV,   &
                                    DMI_DEV,   &
                                    PHALF_DEV, &
                                    TV_DEV, &
                                    CT_DEV, &
                                    CQ_DEV, &
                                    CU_DEV, &
                                    T_DEV, &
                                    QV_DEV, &
                                    TH_DEV, &
                                    U_DEV, &
                                    V_DEV, &
                                    DU_DEV, &
                                    RI_DEV, &
                                    PBLHT_OPTION, &
                                    ! Inoutputs
                                    DIFF_T_DEV, &
                                    DIFF_M_DEV, &
                                    ! Outputs
                                    AKQ_DEV, AKS_DEV, AKV_DEV, &
                                    BKQ_DEV, BKS_DEV, BKV_DEV, &
                                    CKQ_DEV, CKS_DEV, CKV_DEV, &
                                                      EKV_DEV, &
                                    ZPBL_DEV,                  &
                                    ! Diagnostics
                                    AKQODT_DIAG_DEV, AKSODT_DIAG_DEV, AKVODT_DIAG_DEV, &
                                    CKQODT_DIAG_DEV, CKSODT_DIAG_DEV, CKVODT_DIAG_DEV, &
                                    PPBL_DIAG_DEV,   KPBL_DIAG_DEV, TCZPBL_DIAG_DEV,ZPBL2_DIAG_DEV,  &
                                    ZPBL10p_DIAG_DEV, ZPBLHTKE_DIAG_DEV,               &
                                    ZPBLRI_DIAG_DEV, ZPBLRI2_DIAG_DEV, ZPBLTHV_DIAG_DEV)

      STATUS = cudaGetLastError()
      if (STATUS /= 0) then 
         write (*,*) "Error code from POSTLOCK kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         ASSERT_(.FALSE.)
      end if

      STATUS = cudaDeviceSynchronize()

      STATUS = cudaFuncSetCacheConfig('geos_turbulencegridcompmod_beljaars',cudaFuncCachePreferL1)

      call BELJAARS<<<Grid,Block>>>(IM*JM, LM, DT,            &
                                    LAMBDA_B_CONST, C_B_CONST,&
                                    FKV_DEV, BKV_DEV,         &
                                    U_DEV, V_DEV, ZFULL_DEV,  &
                                    VARFLT_DEV, PHALF_DEV     )

      STATUS = cudaGetLastError()
      if (STATUS /= 0) then 
         write (*,*) "Error code from BELJAARS kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         ASSERT_(.FALSE.)
      end if

      STATUS = cudaDeviceSynchronize()

      ! Do LU decomposition; C is not modified.
      ! On exit, B is the main diagonals of the LU
      ! decomposition, and A is the r.h.s multiplier.
      !----------------------------------------------

      STATUS = cudaFuncSetCacheConfig('geos_turbulencegridcompmod_vtrilu',cudaFuncCachePreferL1)

      call VTRILU<<<Grid,Block>>>(IM*JM,LM,AKS_DEV,BKS_DEV,CKS_DEV)
      call VTRILU<<<Grid,Block>>>(IM*JM,LM,AKQ_DEV,BKQ_DEV,CKQ_DEV)
      call VTRILU<<<Grid,Block>>>(IM*JM,LM,AKV_DEV,BKV_DEV,CKV_DEV)

      STATUS = cudaGetLastError()
      if (STATUS /= 0) then 
         write (*,*) "Error code from VTRILU kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         ASSERT_(.FALSE.)
      end if

      STATUS = cudaDeviceSynchronize()

      ! Get the sensitivity of solution to a unit
      ! change in the surface value. B and C are
      ! not modified.
      !------------------------------------------

      STATUS = cudaFuncSetCacheConfig('geos_turbulencegridcompmod_vtrisolvesurf',cudaFuncCachePreferL1)

      call VTRISOLVESURF<<<Grid,Block>>>(IM*JM,LM,BKS_DEV,CKS_DEV,DKS_DEV)
      call VTRISOLVESURF<<<Grid,Block>>>(IM*JM,LM,BKQ_DEV,CKQ_DEV,DKQ_DEV)
      call VTRISOLVESURF<<<Grid,Block>>>(IM*JM,LM,BKV_DEV,CKV_DEV,DKV_DEV)

      STATUS = cudaGetLastError()
      if (STATUS /= 0) then 
         write (*,*) "Error code from VTRISOLVESURF kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         ASSERT_(.FALSE.)
      end if

      STATUS = cudaDeviceSynchronize()

      ! --------------
      ! Kernel is done
      ! --------------

      call MAPL_TimerOff(MAPL,name="---REFRESHKS_RUN" ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (MAPL,name="---REFRESHKS_DATA" ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (MAPL,name="----REFRESHKS_DATA_DEVICE" ,RC=STATUS)
      VERIFY_(STATUS)

      ! ------------------------
      ! Copy outputs to the host
      ! ------------------------

      ! Outputs - Prelims
      ! -----------------
   
      STATUS = cudaMemcpy(ZLE,ZHALF_DEV,IM*JM*(LM+1))

      ! Outputs - Louis
      ! ---------------
   
      STATUS = cudaMemcpy(RI,RI_DEV,IM*JM*(LM+1))
      STATUS = cudaMemcpy(DU,DU_DEV,IM*JM*(LM+1))

      ! Diagnostics - Louis
      ! -------------------

      IF (ASSOCIATED(ALH))  STATUS = cudaMemcpy(ALH, ALH_DIAG_DEV,IM*JM*(LM+1))
      IF (ASSOCIATED(KMLS)) STATUS = cudaMemcpy(KMLS,KMLS_DIAG_DEV,IM*JM*(LM+1))
      IF (ASSOCIATED(KHLS)) STATUS = cudaMemcpy(KHLS,KHLS_DIAG_DEV,IM*JM*(LM+1))

      if (LOCK_ON==1) then
         ! Outputs - Lock
         ! --------------
      
         STATUS = cudaMemcpy(EKM,K_M_ENTR_DEV,IM*JM*(LM+1))
         STATUS = cudaMemcpy(EKH,K_T_ENTR_DEV,IM*JM*(LM+1))
         STATUS = cudaMemcpy(KHSFC,K_SFC_DIAG_DEV,IM*JM*(LM+1))
         STATUS = cudaMemcpy(KHRAD,K_RAD_DIAG_DEV,IM*JM*(LM+1))
         STATUS = cudaMemcpy(ZCLD,ZCLOUD_DEV,IM*JM)
         STATUS = cudaMemcpy(ZRADML,ZRADML_DEV,IM*JM)
         STATUS = cudaMemcpy(ZRADBS,ZRADBASE_DEV,IM*JM)
         STATUS = cudaMemcpy(ZSML,ZSML_DEV,IM*JM)
      
         ! Diagnostics - Lock
         ! ------------------
      
         IF (ASSOCIATED(ZCLDTOP))  STATUS = cudaMemcpy(ZCLDTOP,ZCLDTOP_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(WESFC))    STATUS = cudaMemcpy(WESFC,WENTR_SFC_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(WERAD))    STATUS = cudaMemcpy(WERAD,WENTR_RAD_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(DBUOY))    STATUS = cudaMemcpy(DBUOY,DEL_BUOY_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(VSCSFC))   STATUS = cudaMemcpy(VSCSFC,VSFC_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(VSCRAD))   STATUS = cudaMemcpy(VSCRAD,VRAD_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(KERAD))    STATUS = cudaMemcpy(KERAD,KENTRAD_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(VSCBRV))   STATUS = cudaMemcpy(VSCBRV,VBRV_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(WEBRV))    STATUS = cudaMemcpy(WEBRV,WENTR_BRV_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(DSIEMS))   STATUS = cudaMemcpy(DSIEMS,DSIEMS_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(CHIS))     STATUS = cudaMemcpy(CHIS,CHIS_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(DELSINV))  STATUS = cudaMemcpy(DELSINV,DELSINV_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(SMIXT))    STATUS = cudaMemcpy(SMIXT,SLMIXTURE_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(CLDRF))    STATUS = cudaMemcpy(CLDRF,CLDRADF_DIAG_DEV,IM*JM)
         IF (ASSOCIATED(RADRCODE)) STATUS = cudaMemcpy(RADRCODE,RADRCODE_DIAG_DEV,IM*JM)
      else
         EKM  (:,:,0:LM-1) = 0.0
         EKH  (:,:,0:LM-1) = 0.0
         KHSFC(:,:,0:LM-1) = 0.0
         KHRAD(:,:,0:LM-1) = 0.0
      endif

      ! Outputs - Postlock
      ! ------------------

      STATUS = cudaMemcpy(KH,DIFF_T_DEV,IM*JM*(LM+1))
      STATUS = cudaMemcpy(KM,DIFF_M_DEV,IM*JM*(LM+1))

      STATUS = cudaMemcpy(CKQ,CKQ_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(CKS,CKS_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(CKV,CKV_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(EKV,EKV_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(ZPBL,ZPBL_DEV,IM*JM)

      ! Diagnostics - Postlock
      ! ----------------------
   
      IF (ASSOCIATED(AKSODT)) STATUS = cudaMemcpy(AKSODT,AKSODT_DIAG_DEV,IM*JM*LM)
      IF (ASSOCIATED(CKSODT)) STATUS = cudaMemcpy(CKSODT,CKSODT_DIAG_DEV,IM*JM*LM)
      IF (ASSOCIATED(AKQODT)) STATUS = cudaMemcpy(AKQODT,AKQODT_DIAG_DEV,IM*JM*LM)
      IF (ASSOCIATED(CKQODT)) STATUS = cudaMemcpy(CKQODT,CKQODT_DIAG_DEV,IM*JM*LM)
      IF (ASSOCIATED(AKVODT)) STATUS = cudaMemcpy(AKVODT,AKVODT_DIAG_DEV,IM*JM*LM)
      IF (ASSOCIATED(CKVODT)) STATUS = cudaMemcpy(CKVODT,CKVODT_DIAG_DEV,IM*JM*LM)
      IF (ASSOCIATED(PPBL  )) STATUS = cudaMemcpy(PPBL,PPBL_DIAG_DEV,IM*JM)
      IF (ASSOCIATED(KPBL  )) STATUS = cudaMemcpy(KPBL,KPBL_DIAG_DEV,IM*JM)
      IF (ASSOCIATED(TCZPBL)) STATUS = cudaMemcpy(TCZPBL,TCZPBL_DIAG_DEV,IM*JM)
      IF (ASSOCIATED(ZPBL2 )) STATUS = cudaMemcpy(ZPBL2,ZPBL2_DIAG_DEV,IM*JM)
      IF (ASSOCIATED(ZPBL10p )) STATUS = cudaMemcpy(ZPBL10p,ZPBL10p_DIAG_DEV,IM*JM)
      IF (ASSOCIATED(ZPBLHTKE )) STATUS = cudaMemcpy(ZPBLHTKE,ZPBLHTKE_DIAG_DEV,IM*JM)
      IF (ASSOCIATED(ZPBLRI )) STATUS = cudaMemcpy(ZPBLRI,ZPBLRI_DIAG_DEV,IM*JM)
      IF (ASSOCIATED(ZPBLRI2 )) STATUS = cudaMemcpy(ZPBLRI2,ZPBLRI2_DIAG_DEV,IM*JM)
      IF (ASSOCIATED(ZPBLTHV )) STATUS = cudaMemcpy(ZPBLTHV,ZPBLTHV_DIAG_DEV,IM*JM)

      ! Outputs - Beljaars
      ! ------------------

      STATUS = cudaMemcpy(FKV,FKV_DEV,IM*JM*LM)

      ! Outputs - Decomp
      ! ----------------

      STATUS = cudaMemcpy(AKQ,AKQ_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(AKS,AKS_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(AKV,AKV_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(BKQ,BKQ_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(BKS,BKS_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(BKV,BKV_DEV,IM*JM*LM)

      STATUS = cudaMemcpy(DKQ,DKQ_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(DKS,DKS_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(DKV,DKV_DEV,IM*JM*LM)

      STATUS = cudaDeviceSynchronize()

      call MAPL_TimerOff(MAPL,name="----REFRESHKS_DATA_DEVICE" ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOff(MAPL,name="---REFRESHKS_DATA" ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (MAPL,name="---REFRESHKS_DEALLOC" ,RC=STATUS)
      VERIFY_(STATUS)

      ! ------------------------
      ! Deallocate device arrays
      ! ------------------------
  
      ! Working Arrays
      ! --------------
  
      DEALLOCATE(ZFULL_DEV)
      DEALLOCATE(TV_DEV)
      DEALLOCATE(PV_DEV)
      DEALLOCATE(RDZ_DEV)
      DEALLOCATE(DMI_DEV)
      DEALLOCATE(PFULL_DEV)
      !DEALLOCATE(QA_DEV) ! NOT CURRENTLY USED
  
      ! Inputs - Prelims
      ! ----------------
  
      DEALLOCATE(T_DEV)
      DEALLOCATE(QV_DEV)
      DEALLOCATE(PHALF_DEV)
      DEALLOCATE(TH_DEV)
      DEALLOCATE(QLCN_DEV)
      DEALLOCATE(QLLS_DEV)
      DEALLOCATE(QICN_DEV)
      DEALLOCATE(QILS_DEV)
      !DEALLOCATE(CLCN_DEV) ! NOT CURRENTLY USED
      !DEALLOCATE(CLLS_DEV) ! NOT CURRENTLY USED
  
      ! Inputs - Louis
      ! --------------
  
      DEALLOCATE(U_DEV)
      DEALLOCATE(V_DEV)
      DEALLOCATE(ZPBL_DEV)
  
      ! Inputs - Lock
      ! -------------
  
      DEALLOCATE(TDTLW_IN_DEV)
      DEALLOCATE(U_STAR_DEV)
      DEALLOCATE(B_STAR_DEV)
      DEALLOCATE(FRLAND_DEV)
   
      ! Inputs - Postlock
      ! -----------------
  
      DEALLOCATE(CT_DEV)
      DEALLOCATE(CQ_DEV)
      DEALLOCATE(CU_DEV)
  
      ! Inputs - Beljaars
      ! -----------------
  
      DEALLOCATE(VARFLT_DEV)
  
      ! Outputs - Prelims
      ! -----------------

      DEALLOCATE(ZHALF_DEV)

      ! Outputs - Louis
      ! ---------------
  
      DEALLOCATE(DIFF_M_DEV)
      DEALLOCATE(DIFF_T_DEV)
      DEALLOCATE(RI_DEV)
      DEALLOCATE(DU_DEV)
  
      ! Outputs - Lock
      ! --------------
   
      DEALLOCATE(K_M_ENTR_DEV)
      DEALLOCATE(K_T_ENTR_DEV)
      DEALLOCATE(K_SFC_DIAG_DEV)
      DEALLOCATE(K_RAD_DIAG_DEV)
      DEALLOCATE(ZCLOUD_DEV)
      DEALLOCATE(ZRADML_DEV)
      DEALLOCATE(ZRADBASE_DEV)
      DEALLOCATE(ZSML_DEV)
     
      ! Outputs - Postlock
      ! ------------------
  
      DEALLOCATE(AKQ_DEV)
      DEALLOCATE(AKS_DEV)
      DEALLOCATE(AKV_DEV)
      DEALLOCATE(BKQ_DEV)
      DEALLOCATE(BKS_DEV)
      DEALLOCATE(BKV_DEV)
      DEALLOCATE(CKQ_DEV)
      DEALLOCATE(CKS_DEV)
      DEALLOCATE(CKV_DEV)
      DEALLOCATE(EKV_DEV)
  
      ! Outputs - Beljaars
      ! ------------------
  
      DEALLOCATE(FKV_DEV)
  
      ! Outputs - Decomp
      ! ----------------
  
      DEALLOCATE(DKQ_DEV)
      DEALLOCATE(DKS_DEV)
      DEALLOCATE(DKV_DEV)

      ! Diagnostics - Louis
      ! -------------------
  
      DEALLOCATE(ALH_DIAG_DEV)
      DEALLOCATE(KMLS_DIAG_DEV)
      DEALLOCATE(KHLS_DIAG_DEV)
  
      ! Diagnostics - Lock
      ! ------------------
   
      DEALLOCATE(ZCLDTOP_DIAG_DEV)
      DEALLOCATE(WENTR_SFC_DIAG_DEV)
      DEALLOCATE(WENTR_RAD_DIAG_DEV)
      DEALLOCATE(DEL_BUOY_DIAG_DEV)
      DEALLOCATE(VSFC_DIAG_DEV)
      DEALLOCATE(VRAD_DIAG_DEV)
      DEALLOCATE(KENTRAD_DIAG_DEV)
      DEALLOCATE(VBRV_DIAG_DEV)
      DEALLOCATE(WENTR_BRV_DIAG_DEV)
      DEALLOCATE(DSIEMS_DIAG_DEV)
      DEALLOCATE(CHIS_DIAG_DEV)
      DEALLOCATE(DELSINV_DIAG_DEV)
      DEALLOCATE(SLMIXTURE_DIAG_DEV)
      DEALLOCATE(CLDRADF_DIAG_DEV)
      DEALLOCATE(RADRCODE_DIAG_DEV)
  
      ! Diagnostics - Postlock
      ! ----------------------
   
      DEALLOCATE(AKQODT_DIAG_DEV)
      DEALLOCATE(AKSODT_DIAG_DEV)
      DEALLOCATE(AKVODT_DIAG_DEV)
      DEALLOCATE(CKQODT_DIAG_DEV)
      DEALLOCATE(CKSODT_DIAG_DEV)
      DEALLOCATE(CKVODT_DIAG_DEV)
      DEALLOCATE(PPBL_DIAG_DEV)
      DEALLOCATE(KPBL_DIAG_DEV)
      if(CALC_TCZPBL) then
         DEALLOCATE(TCZPBL_DIAG_DEV)
      end if
      if(CALC_ZPBL2) then
         DEALLOCATE(ZPBL2_DIAG_DEV)
      end if
      if(CALC_ZPBL10p) then
         DEALLOCATE(ZPBL10p_DIAG_DEV)
      end if
      DEALLOCATE(ZPBLHTKE_DIAG_DEV)
      DEALLOCATE(ZPBLRI_DIAG_DEV)
      DEALLOCATE(ZPBLRI2_DIAG_DEV)
      DEALLOCATE(ZPBLTHV_DIAG_DEV)

      STATUS = cudaDeviceSynchronize()

      call MAPL_TimerOff(MAPL,name="---REFRESHKS_DEALLOC" ,RC=STATUS)
      VERIFY_(STATUS)

#else

      call MAPL_TimerOn(MAPL,"---PRELIMS")

      call PRELIMS(IM*JM, LM, DT, &
                   ! Inputs
                   T,    &
                   Q,    &
                   PLE,  &
                   TH,   &
                   QLCN, &
                   QLLS, &
                   QICN, &
                   QILS, &
                   !CLCN, & ! CURRENTLY NOT USED
                   !CLLS, & ! CURRENTLY NOT USED
                   ! Outputs
                   Z,   &
                   ZLE, &
                   TV,  &
                   THV, &
                   RDZ, &
                   DMI, &
                   !QA, & ! CURRENTLY NOT USED
                   PLO  )

      call MAPL_TimerOff(MAPL,"---PRELIMS")

!   Refresh diffusivities: First compute Louis...
!   ---------------------------------------------

      call MAPL_TimerOn (MAPL,name="---LOUIS" ,RC=STATUS)
      VERIFY_(STATUS)

      call LOUIS_KS(IM*JM, LM,                    &
            KH,KM,RI,DU,ZPBL,                     &    
            Z,ZLE,THV,U,V,                        &    
            LOUIS, MINSHEAR, MINTHICK,            &
            LAMBDAM, LAMBDAM2, LAMBDAH, LAMBDAH2, & 
            ZKMENV, ZKHENV,                       &
            AKHMMAX,                              &
            ALH_X, KMLS_X, KHLS_X                 )

      if(associated(ALH )) ALH  = ALH_X
      if(associated(KHLS)) KHLS = KHLS_X
      if(associated(KMLS)) KMLS = KMLS_X

      call MAPL_TimerOff(MAPL,name="---LOUIS" ,RC=STATUS)
      VERIFY_(STATUS)

!   ...then add Lock.
!--------------------

      if (LOCK_ON==1) then
         call MAPL_TimerOn(MAPL,"---LOCK")

         CALL ENTRAIN(IM*JM,LM,                 &
                      RADLW,                    &
                      USTAR,                    &
                      BSTAR,                    &
                      FRLAND,                   &
                      T,                        &
                      Q,                        &
                      QLLS,                     &
                      QILS,                     &
                      !QA,                      & ! Not currently used in entrain
                      U,                        &
                      V,                        &
                      Z,                        &
                      PLO,                      &
                      ZLE,                      &
                      PLE,                      &
                      KM,                       &
                      KH,                       &
                      EKM,                      &
                      EKH,                      &
                      KHSFC,                    &
                      KHRAD,                    &
                      ZCLD,                     &
                      ZRADML,                   &
                      ZRADBS,                   &
                      ZSML,                     &
                      ! these are diagnostics
                      ZCLDTOP_X,                &
                      WESFC_X,                  &
                      WERAD_X,                  &
                      DBUOY_X,                  &
                      VSCSFC_X,                 &
                      VSCRAD_X,                 &
                      KERAD_X,                  &
                      VSCBRV_X,                 &
                      WEBRV_X,                  &
                      DSIEMS_X,                 &
                      CHIS_X,                   &
                      DELSINV_X,                &
                      SMIXT_X,                  &
                      CLDRF_X,                  &
                      RADRCODE_X,               &
                      ! these are input params
                      PRANDTLSFC, PRANDTLRAD,   &
                      BETA_SURF, BETA_RAD,      &
                      TPFAC_SURF, ENTRATE_SURF, &
                      PCEFF_SURF, KHRADFAC, KHSFCFAC )

         call MAPL_TimerOff(MAPL,"---LOCK")

         ! If we want an export diagnostic, fill it
         if (ASSOCIATED(ZCLDTOP))  ZCLDTOP  = ZCLDTOP_X
         if (ASSOCIATED(WESFC))    WESFC    = WESFC_X
         if (ASSOCIATED(WERAD))    WERAD    = WERAD_X
         if (ASSOCIATED(DBUOY))    DBUOY    = DBUOY_X
         if (ASSOCIATED(VSCRAD))   VSCRAD   = VSCRAD_X
         if (ASSOCIATED(VSCSFC))   VSCSFC   = VSCSFC_X
         if (ASSOCIATED(KERAD))    KERAD    = KERAD_X
         if (ASSOCIATED(VSCBRV))   VSCBRV   = VSCBRV_X
         if (ASSOCIATED(WEBRV))    WEBRV    = WEBRV_X
         if (ASSOCIATED(DSIEMS))   DSIEMS   = DSIEMS_X
         if (ASSOCIATED(CHIS))     CHIS     = CHIS_X
         if (ASSOCIATED(DELSINV))  DELSINV  = DELSINV_X
         if (ASSOCIATED(SMIXT))    SMIXT    = SMIXT_X
         if (ASSOCIATED(CLDRF))    CLDRF    = CLDRF_X
         if (ASSOCIATED(RADRCODE)) RADRCODE = RADRCODE_X

      else
         EKM  (:,:,0:LM-1) = 0.0
         EKH  (:,:,0:LM-1) = 0.0
         KHSFC(:,:,0:LM-1) = 0.0
         KHRAD(:,:,0:LM-1) = 0.0
      endif

      call MAPL_TimerOn (MAPL,"---POSTLOCK")

      KPBLMIN  = count(PREF < 50000.)

      call POSTLOCK(IM*JM,LM,DT,&
                    ! Constants
                    KPBLMIN,    &
                    CALC_TCZPBL,&
                    CALC_ZPBL2,&
                    CALC_ZPBL10p,&
                    ! Inputs
                    Z, ZLE, &
                    PLO,&
                    RDZ,&
                    DMI,&
                    PLE,&
                    TV, &
                    CT, &
                    CQ, &
                    CU, &
                    T,  &
                    Q,  &
                    TH, &
                    U,  &
                    V,  & 
                    DU,  &
                    RI,  &
                    PBLHT_OPTION, &
                    ! Inoutputs
                    KH, &
                    KM, &
                    ! Outputs
                    AKQ, AKS, AKV, &
                    BKQ, BKS, BKV, &
                    CKQ, CKS, CKV, &
                              EKV, &
                    ZPBL,          &
                    ! Diagnostics
                    AKQODT_X, AKSODT_X, AKVODT_X, &
                    CKQODT_X, CKSODT_X, CKVODT_X, &
                    PPBL_X,   KPBL_X, TCZPBL, ZPBL2, ZPBL10p, &
                    ZPBLHTKE_X, ZPBLRI_X, ZPBLRI2_X, ZPBLTHV_X)

      ! If we want an export diagnostic, fill it
      IF (ASSOCIATED(AKQODT)) AKQODT = AKQODT_X
      IF (ASSOCIATED(AKSODT)) AKSODT = AKSODT_X
      IF (ASSOCIATED(AKVODT)) AKVODT = AKVODT_X
      IF (ASSOCIATED(CKQODT)) CKQODT = CKQODT_X
      IF (ASSOCIATED(CKSODT)) CKSODT = CKSODT_X
      IF (ASSOCIATED(CKVODT)) CKVODT = CKVODT_X
      IF (ASSOCIATED(PPBL  )) PPBL   = PPBL_X
      IF (ASSOCIATED(KPBL  )) KPBL   = KPBL_X
      IF (ASSOCIATED(ZPBLHTKE )) ZPBLHTKE   = ZPBLHTKE_X
      IF (ASSOCIATED(ZPBLRI )) ZPBLRI   = ZPBLRI_X
      IF (ASSOCIATED(ZPBLRI2 )) ZPBLRI2   = ZPBLRI2_X
      IF (ASSOCIATED(ZPBLTHV )) ZPBLTHV   = ZPBLTHV_X

      call MAPL_TimerOff(MAPL,"---POSTLOCK")

      call MAPL_TimerOn(MAPL,"---BELJAARS")

      call BELJAARS(IM*JM, LM, DT, &
                    LAMBDA_B, C_B, &
                    FKV, BKV,      &
                    U, V, Z,       &
                    VARFLT, PLE    )

      call MAPL_TimerOff(MAPL,"---BELJAARS")

      call MAPL_TimerOn(MAPL,"---DECOMP")

! Do LU decomposition; C is not modified.
! On exit, B is the main diagonals of the LU
! decomposition, and A is the r.h.s multiplier.
!----------------------------------------------

      call VTRILU(IM*JM,LM,AKS,BKS,CKS)
      call VTRILU(IM*JM,LM,AKQ,BKQ,CKQ)
      call VTRILU(IM*JM,LM,AKV,BKV,CKV)

! Get the sensitivity of solution to a unit
! change in the surface value. B and C are
! not modified.
!------------------------------------------

      call VTRISOLVESURF(IM*JM,LM,BKS,CKS,DKS)
      call VTRISOLVESURF(IM*JM,LM,BKQ,CKQ,DKQ)
      call VTRISOLVESURF(IM*JM,LM,BKV,CKV,DKV)

      call MAPL_TimerOff(MAPL,"---DECOMP")


#endif
      if(ALLOC_TCZPBL) deallocate(TCZPBL)
      if(ALLOC_ZPBL2) deallocate(ZPBL2)
      if(ALLOC_ZPBL10p) deallocate(ZPBL10p)
      RETURN_(ESMF_SUCCESS)
     end subroutine REFRESH

!=============================================================================
!=============================================================================

!BOP

! !CROUTINE: DIFFUSE -- Solves for semi-implicit diffusive tendencies assuming fixed surface conditions.  

! !INTERFACE:

  subroutine DIFFUSE(IM,JM,LM,RC)

! !ARGUMENTS:

    integer,           intent(IN)       :: IM,JM,LM
    integer, optional, intent(OUT)      :: RC

! !DESCRIPTION: {\tt DIFFUSE} computes semi-implicit tendencies of all fields in
!  the TR bundle. Each field is examined for three attributes: {\tt DiffuseLike},
!  {\tt FriendlyToTURBULENCE}, and {\tt WeightedTendency}. These determine the behavior of
!  {\tt DIFFUSE} for that field. {\tt DiffuseLike} can be either 'U', 'Q', or 'S'; the default is 'Q'.
!  {\tt FriendlyToTURBULENCE}, and {\tt WeightedTendency} are ESMF logicals.
!  If {\tt FriendlyToTURBULENCE} is true, the field in TR is updated directly; otherwise
!  it is left untouched. In either case, If the corresponding pointer TRI bundle is associated, the
!  tendencies are returned there. If {\tt WeightedTendency} is true, the tendency in TRI, if any,
!  is pressure weighted.

!EOP

    character(len=ESMF_MAXSTR)          :: IAm='Diffuse'
    integer                             :: STATUS

    character(len=ESMF_MAXSTR)          :: TYPE
    character(len=ESMF_MAXSTR)          :: NAME
    type (ESMF_Field)                   :: FIELD
    type (ESMF_Array)                   :: ARRAY
    type (ESMF_FieldBundle)             :: TR
    type (ESMF_FieldBundle)             :: TRI
    type (ESMF_FieldBundle)             :: TRG
    type (ESMF_FieldBundle)             :: FSTAR
    type (ESMF_FieldBundle)             :: DFSTAR
    real, dimension(:,:,:), pointer     :: S, SOI, SOD
    real, dimension(:,:),   pointer     :: SG, SF, SDF, CX, SRG
    real, dimension(:,:,:), pointer     :: DX
    real, dimension(:,:,:), pointer     :: AK, BK, CK

    integer                             :: KM, K,L
    logical                             :: FRIENDLY
    logical                             :: WEIGHTED
    real, dimension(IM,JM,LM)           :: DP, SX

! Get the bundles containing the quantities to be diffused, 
!     their tendencies, their surface values, their surface
!     fluxes, and the derivatives of their surface fluxes
!     wrt the surface values. 
!----------------------------------------------------------

    call ESMF_StateGet(IMPORT, 'TR' ,    TR,     RC=STATUS); VERIFY_(STATUS)
    call ESMF_StateGet(IMPORT, 'TRG',    TRG,    RC=STATUS); VERIFY_(STATUS)

    call ESMF_StateGet(EXPORT, 'TRI',    TRI,    RC=STATUS); VERIFY_(STATUS)
    call ESMF_StateGet(EXPORT, 'FSTAR',  FSTAR,  RC=STATUS); VERIFY_(STATUS)
    call ESMF_StateGet(EXPORT, 'DFSTAR', DFSTAR, RC=STATUS); VERIFY_(STATUS)

! Count the firlds in TR...
!--------------------------

    call ESMF_FieldBundleGet(TR, fieldCOUNT=KM, RC=STATUS)
    VERIFY_(STATUS)

! ...and make sure the other bundles are the same.
!-------------------------------------------------

    call ESMF_FieldBundleGet(TRI,    FieldCount=K , RC=STATUS)
    VERIFY_(STATUS)
    ASSERT_(KM==K)
    call ESMF_FieldBundleGet(TRG,    FieldCount=K , RC=STATUS)
    VERIFY_(STATUS)
    ASSERT_(KM==K)
    call ESMF_FieldBundleGet(FSTAR,  FieldCount=K , RC=STATUS)
    VERIFY_(STATUS)
    ASSERT_(KM==K)
    call ESMF_FieldBundleGet(DFSTAR, FieldCount=K , RC=STATUS)
    VERIFY_(STATUS)
    ASSERT_(KM==K)

! Pressure thickness of layers
!-----------------------------

    DP = PLE(:,:,1:LM)-PLE(:,:,0:LM-1)

! Loop over all quantities to be diffused.
!----------------------------------------

    do K=1,KM

! Get the Kth Field and its name from tracer bundle
!--------------------------------------------------

       call ESMF_FieldBundleGet(TR, K, FIELD, RC=STATUS)
       VERIFY_(STATUS)

       call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
       VERIFY_(STATUS)

! Get item's diffusion type (U, S or Q; default is Q)
!----------------------------------------------------

       call ESMF_AttributeGet(FIELD, NAME="DiffuseLike",         &
            VALUE=TYPE, DEFAULTVALUE=dflt_q,    RC=STATUS)
       VERIFY_(STATUS)

! Get item's friendly status (default is not friendly)
!-----------------------------------------------------

       call ESMF_AttributeGet(FIELD, NAME="FriendlyToTURBULENCE", &
            VALUE=FRIENDLY, DEFAULTVALUE=dflt_false,    RC=STATUS)
       VERIFY_(STATUS)

! Get item's weighting (default is unweighted tendencies)
!--------------------------------------------------------

       call ESMF_AttributeGet(FIELD, NAME="WeightedTendency",   &
            VALUE=WEIGHTED, DEFAULTVALUE=dflt_false,    RC=STATUS)
       VERIFY_(STATUS)

! Get pointer to the quantity, its tendency, its surface value,
!   the surface flux, and the sensitivity of the surface flux.
! -------------------------------------------------------------

       call ESMFL_BundleGetPointerToData(TR    ,      NAME,         S  , RC=STATUS)
       VERIFY_(STATUS)
       call ESMFL_BundleGetPointerToData(TRI   , trim(NAME)//'IT' , SOI, RC=STATUS)
       VERIFY_(STATUS)
       call ESMFL_BundleGetPointerToData(TRG   , trim(NAME)//'HAT', SRG, RC=STATUS)
       VERIFY_(STATUS)
       call ESMFL_BundleGetPointerToData(FSTAR , trim(NAME)//'FLX', SF , RC=STATUS)
       VERIFY_(STATUS)
       call ESMFL_BundleGetPointerToData(DFSTAR, trim(NAME)//'DFL', SDF, RC=STATUS)
       VERIFY_(STATUS)

! The quantity must exist; others are optional.
!----------------------------------------------

       ASSERT_(associated(S ))

! If the surface values does not exists, we assume zero flux.
!------------------------------------------------------------
       
       if(associated(SRG)) then
          SG => SRG
       else
          allocate (SG(0,0), stat=STATUS)
          VERIFY_(STATUS)
       end if

! Pick the right exchange coefficients
!-------------------------------------

       if     ( TYPE=='U' ) then ! Momentum
          CX => CU
          DX => DKV
          AK => AKV; BK => BKV; CK => CKV
       else if( TYPE=='Q' ) then ! Water Vapor or other tracers
          CX => CQ
          DX => DKQ
          AK => AKQ; BK => BKQ; CK => CKQ
       else if( TYPE=='S' ) then ! Heat
          CX => CT
          DX => DKS
          AK => AKS; BK => BKS; CK => CKS
       else
          RETURN_(ESMF_FAILURE)
       endif

! Copy diffused quantity to temp buffer
! ------------------------------------------
       
       SX = S

! Solve for semi-implicit changes. This modifies SX
! -------------------------------------------------

if( TYPE .ne. 'Q' ) then
       call VTRISOLVE(AK,BK,CK,SX,SG)
else
   SG = 0.0
endif

! Compute the surface fluxes
!---------------------------

       if(associated(SF)) then
          if(size(SG)>0) then
             SF = CX*(SG - SX(:,:,LM))
          else
             SF = 0.0
          end if
       end if

! Create tendencies
!------------------

       if(associated(SOI)) then
          if( WEIGHTED ) then
             SOI = ( (SX - S)/DT )*DP
          else
             SOI = ( (SX - S)/DT )
          endif
       end if

       if( NAME=='S' ) then
          SINC = ( (SX - S)/DT )
       end if

! Update friendlies
!------------------

       if(FRIENDLY) then
          S = SX
       end if

! Compute the derivative of the surface flux wrt the surface value
!-----------------------------------------------------------------

       if(associated(SDF)) then
          SDF = CX * (1.0-DX(:,:,LM))
       endif

       if(.not.associated(SRG)) then
          deallocate (SG)
       end if

    enddo ! End loop over all quantities to be diffused
! -----------------------------------------------------

    RETURN_(ESMF_SUCCESS)
  end subroutine DIFFUSE

end subroutine RUN1


!*********************************************************************
!*********************************************************************
!*********************************************************************


!BOP

! !IROUTINE: RUN2 -- The second run stage for the TURBULENCE component

! !INTERFACE:

  subroutine RUN2 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Second run stage of {\tt GEOS\_TurbulenceGridComp} performs
!    the updates due to changes in surface quantities. Its input are the changes in 
!    surface quantities during the time step. It can also compute the frictional
!    dissipation terms as exports, but these are not added to the temperatures.


!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp), pointer       :: MAPL
    type (ESMF_Config      )            :: CF
    type (ESMF_State       )            :: INTERNAL 

! Local variables

    integer                             :: IM, JM, LM
    real                                :: DT

    real, pointer, dimension(:,:)       :: VARFLT
    real, pointer, dimension(:,:)       :: LATS

! Begin... 
!---------

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'Run2'

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"-RUN2")

! Get parameters from generic state.
!-----------------------------------

          call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,          &
                               LATS = LATS,                  &
                               INTERNAL_ESMF_STATE=INTERNAL, &
                                                   RC=STATUS )
      VERIFY_(STATUS)

! Get configuration from component
!---------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Get application's timestep from configuration
!----------------------------------------------

    call ESMF_ConfigGetAttribute( CF, DT, Label="RUN_DT:" , RC=STATUS)
    VERIFY_(STATUS)


    call MAPL_GetPointer(IMPORT,VARFLT,  'VARFLT', RC=STATUS)
    VERIFY_(STATUS)

! Solve the free atmosphere problem
! ---------------------------------

    call MAPL_TimerOn (MAPL,"--UPDATE")
      call UPDATE(IM,JM,LM,LATS,RC=STATUS)
      VERIFY_(STATUS)
    call MAPL_TimerOff(MAPL,"--UPDATE")

!  All done with RUN
!-------------------

    call MAPL_TimerOff(MAPL,"-RUN2")
    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  contains

!BOP

! !CROUTINE: UPDATE -- Updates diffusive effects for changes at surface.

! !INTERFACE:

    subroutine UPDATE(IM,JM,LM,LATS,RC)

! !ARGUMENTS:

      integer,           intent(IN)       :: IM,JM,LM
      integer, optional, intent(OUT)      :: RC

! !DESCRIPTION: 
!    Some description

!EOP
  
 
      character(len=ESMF_MAXSTR)          :: IAm='Update'
      integer                             :: STATUS

      character(len=ESMF_MAXSTR)          :: TYPE
      character(len=ESMF_MAXSTR)          :: NAME
      type (ESMF_Field)                   :: FIELD
      type (ESMF_FieldBundle)             :: TR
      type (ESMF_FieldBundle)             :: TRI
      type (ESMF_FieldBundle)             :: DTG
      type (ESMF_FieldBundle)             :: FSTAR
      type (ESMF_FieldBundle)             :: DFSTAR
      real, dimension(:,:,:), pointer     :: PLE
      real, dimension(:,:,:), pointer     :: ZLE
      real, dimension(:,:,:), pointer     :: S, SOI, SINC, INTDIS, TOPDIS
      real, dimension(:,:  ), pointer     :: DSG, SF, SDF, SRFDIS
      real, dimension(:,:  ), pointer     :: KETRB, KESRF, KETOP, KEINT
      real, dimension(:,:,:), pointer     :: DKS, DKV, DKQ, DKX, EKV, FKV
      integer                             :: KM, K, L, I, J
      logical                             :: FRIENDLY
      logical                             :: WEIGHTED
      real, dimension(IM,JM,LM)           :: DP, SX
      real, dimension(IM,JM,LM-1)         :: DF
      integer, allocatable                :: KK(:)

! The following variables are for SHVC parameterization

      real, dimension(IM,JM,LM)           :: SOIOFS, XINC
      real,    dimension(IM,JM)           :: z500, z1500, z7000, STDV
      integer, dimension(IM,JM)           :: L500, L1500, L7000
      integer, dimension(IM,JM)           :: LTOPS,LBOT,LTOPQ
      logical, dimension(IM,JM)           :: DidSHVC
      real                                :: REDUFAC, SUMSOI
      real                                :: SHVC_CRIT
      real                                :: SHVC_1500, SHVC_ZDEPTH
      real                                :: lat_in_degrees, lat_effect
      real,  dimension(IM,JM)             :: LATS
      real                                :: SHVC_ALPHA, SHVC_EFFECT, scaling
      logical                             :: DO_SHVC
      integer                             :: KS

      character(len=ESMF_MAXSTR) :: GRIDNAME
      character(len=4)           :: imchar
      character(len=2)           :: dateline
      integer                    :: imsize,nn

! Pressure-weighted dissipation heating rates
!--------------------------------------------

      call MAPL_GetPointer(EXPORT, KETRB ,  'KETRB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KESRF ,  'KESRF' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KETOP ,  'KETOP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KEINT ,  'KEINT' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, SRFDIS,  'SRFDIS',                  &
                       alloc=associated(KETRB) .or. associated(KESRF), &
                                                              RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, INTDIS,  'INTDIS',                  &
                       alloc=associated(KETRB) .or. associated(KEINT), &
                                                              RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TOPDIS,  'TOPDIS',                  &
                       alloc=associated(KETRB) .or. associated(KETOP), &
                                                              RC=STATUS)
      VERIFY_(STATUS)

! SHVC Resource parameters. SHVC_EFFECT can be set to zero to turn-off SHVC.
! SHVC_EFFECT = 1. is the tuned value for  2 degree horizontal resolution.
! It should be set to a lower number at higher resolution.

      call MAPL_GetResource( MAPL, SHVC_EFFECT, 'SHVC_EFFECT:', default=0.   , RC=STATUS )
      VERIFY_(STATUS)

      DO_SHVC = SHVC_EFFECT > 0.0

      if(DO_SHVC) then
         call MAPL_GetResource( MAPL, SHVC_CRIT,   'SHVC_CRIT:'  , default=300. , RC=STATUS )
         VERIFY_(STATUS)
         call MAPL_GetResource( MAPL, SHVC_ALPHA,  'SHVC_ALPHA:' , default=1.   , RC=STATUS )
         VERIFY_(STATUS)
         call MAPL_GetResource( MAPL, SHVC_1500,   'SHVC_1500:'  , default=2100., RC=STATUS )
         VERIFY_(STATUS)
         call MAPL_GetResource( MAPL, SHVC_ZDEPTH, 'SHVC_ZDEPTH:', default=3500., RC=STATUS )
         VERIFY_(STATUS)

         call MAPL_GetResource( MAPL, GRIDNAME, 'AGCM_GRIDNAME:', RC=STATUS )
         VERIFY_(STATUS)
         GRIDNAME =  AdjustL(GRIDNAME)
               nn = len_trim(GRIDNAME)
         dateline = GRIDNAME(nn-1:nn)
           imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
           read(imchar,*) imsize
         if(dateline.eq.'CF') imsize  = imsize*4

         if( imsize.le.200       ) scaling = 1.0  !              Resolution >= 2.000-deg
         if( imsize.gt.200 .and. &
             imsize.le.400       ) scaling = 1.0  !  2.000-deg > Resolution >= 1.000-deg
         if( imsize.gt.400 .and. &
             imsize.le.800       ) scaling = 7.0  !  1.000-deg > Resolution >= 0.500-deg
         if( imsize.gt.800 .and. &
             imsize.le.1600      ) scaling = 7.0  !  0.500-deg > Resolution >= 0.250-deg
         if( imsize.gt.1600      ) scaling = 7.0  !  0.250-deg > Resolution 
      end if

! Get imports
!------------

      call MAPL_GetPointer(IMPORT,    PLE,     'PLE', RC=STATUS); VERIFY_(STATUS)

! Get the tendecy sensitivities computed in RUN1
!-----------------------------------------------

      call MAPL_GetPointer(INTERNAL, DKS,   'DKS',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, DKV,   'DKV',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, DKQ,   'DKQ',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, EKV,   'EKV',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, FKV,   'FKV',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, ZLE,   'ZLE',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, SINC,  'SINC',    RC=STATUS)
      VERIFY_(STATUS)

! Get the bundles containing the quantities to be diffused, 
!     their tendencies, their surface values, their surface
!     fluxes, and the derivatives of their surface fluxes
!     wrt the surface values. 
!----------------------------------------------------------

      call ESMF_StateGet(IMPORT, 'TR' ,    TR,     RC=STATUS); VERIFY_(STATUS)
      call ESMF_StateGet(IMPORT, 'DTG',    DTG,    RC=STATUS); VERIFY_(STATUS)

      call ESMF_StateGet(EXPORT, 'TRI',    TRI,    RC=STATUS); VERIFY_(STATUS)
      call ESMF_StateGet(EXPORT, 'FSTAR' , FSTAR,  RC=STATUS); VERIFY_(STATUS)
      call ESMF_StateGet(EXPORT, 'DFSTAR', DFSTAR, RC=STATUS); VERIFY_(STATUS)

! Count them...
!--------------

      call ESMF_FieldBundleGet(TR , FieldCount=KM, RC=STATUS)
      VERIFY_(STATUS)

! and make sure the other bundles are the same.
!----------------------------------------------

      call ESMF_FieldBundleGet(DTG, FieldCount=K , RC=STATUS)
      VERIFY_(STATUS)

      ASSERT_(KM==K)

! KK gives the order in which quantities will be process.
!--------------------------------------------------------

      allocate(KK(KM), stat=STATUS)
      VERIFY_(STATUS)

      do K = 1,KM
         KK(K) = K
      end do

! Clear the accumulators for the dissipation.
!--------------------------------------------

      if(associated(SRFDIS)) SRFDIS = 0.0
      if(associated(INTDIS)) INTDIS = 0.0
      if(associated(TOPDIS)) TOPDIS = 0.0
      if(associated(KETRB )) KETRB  = 0.0
      if(associated(KESRF )) KESRF  = 0.0
      if(associated(KETOP )) KETOP  = 0.0
      if(associated(KEINT )) KEINT  = 0.0

! Pressure thickness of layers
!-----------------------------

      DP = PLE(:,:,1:LM)-PLE(:,:,0:LM-1)

! Section 1 of 2. SHVC parameterization (W. Chao, J. Atmos. Sci., May 2012, P.1547) 
!  Defining the top and bottom levels of the heat and moisture redistribution layer
!----------------------------------------------------------------------------------

      SHVC_INIT: if(DO_SHVC) then

! Ensure that S is processed first. This only matters for SHVC
!-------------------------------------------------------------

         KS = 0

         do K = 1,KM
            call ESMF_FieldBundleGet(TR, K, FIELD, RC=STATUS)
            VERIFY_(STATUS)

            call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
            VERIFY_(STATUS)

            if    (NAME == 'S') then
               KS=KK(1); KK(1)=K; KK(K)=KS
            end if
         end do

         ASSERT_(KS /= 0 )

! SHVC super-layers
!------------------

         z500  =  500.
         z1500 = 1500.
         z7000 = 7000.

         STDV = sqrt(varflt*scaling)   ! Scaling VARFLT based on resolution

         where (STDV >=700.)
            z1500 = SHVC_1500                   
         endwhere

         where ( (STDV >300.) .and. (STDV <700.) )
            z1500 = 1500.+ (SHVC_1500-1500.)* (STDV - 300.)/400.
         endwhere  

         z7000 = z1500 + SHVC_ZDEPTH



         L500=1.
         do L=LM,2,-1
            where (ZLE(:,:,L) <= z500 .and. ZLE(:,:,L-1) > z500)     
               L500=L-1    
            endwhere
         enddo

         L1500=1.
         do L=LM,2,-1
            where (ZLE(:,:,L) <= z1500 .and. ZLE(:,:,L-1) > z1500)    
               L1500=L-1
            endwhere
         enddo

         L7000=1.
         do L=LM,2,-1
            where (ZLE(:,:,L) <= z7000 .and. ZLE(:,:,L-1) > z7000)    
               L7000=L-1
            endwhere
         enddo

         LBOT  = L1500-1         
         LTOPS = L7000
         LTOPQ = L1500-(LM-L500)*2

         SOIOFS = 0.0

      end if SHVC_INIT

! Loop over all quantities to be diffused.
!-----------------------------------------

       TRACERS: do KS=1,KM

         K = KK(KS)

! Get Kth field from bundle
!--------------------------
          
         call ESMF_FieldBundleGet(TR, K, FIELD, RC=STATUS)
         VERIFY_(STATUS)

         call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
         VERIFY_(STATUS)

! Get item's diffusion type (U, S or Q; default is Q)
!----------------------------------------------------

         call ESMF_AttributeGet(FIELD, NAME="DiffuseLike",         &
              VALUE=TYPE, DEFAULTVALUE=dflt_Q,    RC=STATUS)
         VERIFY_(STATUS)

! Get item's friendly status (default is not friendly)
!-----------------------------------------------------

         call ESMF_AttributeGet(FIELD, NAME="FriendlyToTURBULENCE", &
              VALUE=Friendly, DEFAULTVALUE=dflt_false,    RC=STATUS)
         VERIFY_(STATUS)

! Get item's weighting (default is unweighted tendencies)
!--------------------------------------------------------

         call ESMF_AttributeGet(FIELD, NAME="WeightedTendency",    &
              VALUE=WEIGHTED, DEFAULTVALUE=dflt_false,    RC=STATUS)
         VERIFY_(STATUS)

! Get pointers to the quantity, its tendency, its surface increment,
!   the preliminary surface flux, and the sensitivity of the surface
!   flux to the surface value.
! ------------------------------------------------------------------

         call ESMFL_BundleGetPointerToData(TR    ,      NAME,         S  , RC=STATUS)
         VERIFY_(STATUS)
         call ESMFL_BundleGetPointerToData(TRI   , trim(NAME)//'IT' , SOI, RC=STATUS)
         VERIFY_(STATUS)
         call ESMFL_BundleGetPointerToData(DTG   , trim(NAME)//'DEL', DSG, RC=STATUS)
         VERIFY_(STATUS)
         call ESMFL_BundleGetPointerToData(FSTAR , trim(NAME)//'FLX', SF , RC=STATUS)
         VERIFY_(STATUS)
         call ESMFL_BundleGetPointerToData(DFSTAR, trim(NAME)//'DFL', SDF, RC=STATUS)
         VERIFY_(STATUS)

!  Point to the appropriate sensitivity
!--------------------------------------

         if      ( TYPE=='U' ) then
            DKX => DKV
         else if ( TYPE=='Q' ) then
            DKX => DKQ
            DKX = 0.0
         else if ( TYPE=='S' ) then
            DKX => DKS
         else
            RETURN_(ESMF_FAILURE)
         end if

! Update diffused quantity
!-------------------------

         SX = S

         if(associated(DSG)) then
            do L=1,LM 
               SX(:,:,L) = SX(:,:,L) + DKX(:,:,L)*DSG 
            end do
         end if

! Increment the dissipation
!-------------------------- 

         if( TYPE=='U' ) then
            if(associated(KETRB )) KETRB  = 0.0
            if(associated(KESRF )) KESRF  = 0.0
            if(associated(KETOP )) KETOP  = 0.0
            if(associated(KEINT )) KEINT  = 0.0
            if(associated(INTDIS)) then
               DF             = (0.5/(MAPL_CP))*EKV(:,:,1:LM-1)*(SX(:,:,1:LM-1)-SX(:,:,2:LM))**2
               INTDIS(:,:,1:LM-1) = INTDIS(:,:,1:LM-1) + DF
               INTDIS(:,:,2:LM  ) = INTDIS(:,:,2:LM  ) + DF
               if(associated(KETRB)) then
                  do L=1,LM
                     KETRB = KETRB - INTDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
                  end do
               end if
               if(associated(KEINT)) then
                  do L=1,LM
                     KEINT = KEINT - INTDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
                  end do
               end if
            endif
            if(associated(TOPDIS)) then
               TOPDIS = TOPDIS + (1.0/(MAPL_CP))*FKV*SX**2
               if(associated(KETRB)) then
                  do L=1,LM
                     KETRB = KETRB - TOPDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
                  end do
               end if
               if(associated(KETOP)) then
                  do L=1,LM
                     KETOP = KETOP - TOPDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
                  end do
               end if
             endif
            if(associated(SRFDIS)) then
               SRFDIS = SRFDIS + (1.0/(MAPL_CP))*EKV(:,:,LM)*SX(:,:,LM)**2
               if(associated(KETRB)) KETRB = KETRB - SRFDIS* (MAPL_CP/MAPL_GRAV)
               if(associated(KESRF)) KESRF = KESRF - SRFDIS* (MAPL_CP/MAPL_GRAV)
            endif
         end if

! Update tendencies
! -----------------

         if(associated(SOI).and.associated(DSG)) then
            if( WEIGHTED ) then
               do L=1,LM
                  SOI(:,:,L) = SOI(:,:,L) +  (DKX(:,:,L)*DSG/DT)*DP(:,:,L)
               end do
            else
               do L=1,LM
                  SOI(:,:,L) = SOI(:,:,L) +  (DKX(:,:,L)*DSG/DT)
               end do
            endif
         end if

! Section 2 of 2. SHVC parameterization   (W. Chao,  J. Atmos. Sci., 2012, p1547)   
!  To use SHVC set SHVC_EFFECT in AGCM.rc to > 0.0.
!--------------------------------------------------------------------------------

         RUN_SHVC: if (DO_SHVC) then

            XINC = 0.0

            S_or_Q: if (NAME=='S') then

               if(associated(DSG)) then
                  do L=1,LM
                     SINC(:,:,L) = SINC(:,:,L) +  (DKX(:,:,L)*DSG/DT)
                  end do
               end if

               do I=1,IM
                  do J=1,JM
                lat_effect = 1.
                lat_in_degrees= ABS(LATS(I,J)/(3.14159/2.)*90.)
                if (lat_in_degrees >=42.) lat_effect=0.
                if (lat_in_degrees >37. .and. lat_in_degrees < 42.)   &
                lat_effect = 1.0 - (lat_in_degrees-37.)/(42.-37.)
                     if (STDV(I,J) > SHVC_CRIT) then

                        SUMSOI = sum(SINC(I,J,L500(I,J):LM)*DP(I,J,L500(I,J):LM))
                        DidSHVC(I,J) = SUMSOI >= 0.0

                        if (DidSHVC(I,J)) then
                           if (STDV(I,J) >= 800.) then
                              REDUFAC = 1.0
                           elseif (STDV(i,j) >700. .and. STDV(I,J) <800.) then
                              REDUFAC = 0.95 + 0.05*(STDV(I,J)-700.)/100.
                           else
                              REDUFAC = max(min((STDV(I,J)-SHVC_CRIT)/100.,0.95),0.0)
                           end if

                           REDUFAC = REDUFAC * SHVC_EFFECT  *lat_effect       

                           SUMSOI = 0.
                           do L=L500(i,j),LM
                              SUMSOI        = SUMSOI + SINC(I,J,L)*REDUFAC*DP(I,J,L)
                              XINC  (I,J,L) = -SINC(I,J,L) * REDUFAC
                              SOIOFS(I,J,L) =  XINC(I,J,L) / SX(I,J,L)
                           enddo   !do L

                           XINC(I,J,LTOPS(I,J):LBOT(I,J)) = SUMSOI/SUM(DP(I,J,LTOPS(I,J):LBOT(I,J)))
                        endif
                     else
                        DidSHVC(I,J) = .false.
                     endif   ! end of if (STDV>SHVC_CRIT)
                  enddo   !do J
               enddo   !do I

            elseif (NAME == 'Q')  then

! SHVC_ALPHA below is the alpha factor mentioned on page 1552 of Chao (2012, cited above)
!----------------------------------------------------------------------------------------

               do J=1,JM
                  do I=1,IM
                     if (DidSHVC(I,J)) then
                        SUMSOI = 0.
                        do L=L500(I,J),LM
                           XINC(I,J,L) = SHVC_ALPHA*(SOIOFS(I,J,L)*SX(I,J,L))
                           SUMSOI      = SUMSOI +  XINC(I,J,L)*DP(I,J,L)
                        enddo

                        XINC(I,J,LTOPQ(I,J):LBOT(I,J)) = - SUMSOI/SUM(DP(I,J,LTOPQ(I,J):LBOT(I,J)))
                     endif
                  enddo
               enddo

            endif S_or_Q

            if (name == 'S' .or. name == 'Q') then
               SX  = SX + XINC * DT

               if(associated(SOI)) then
                  if(WEIGHTED) then
                     SOI = SOI + XINC*DP
                  else
                     SOI = SOI + XINC
                  end if
               end if
            end if


         end if RUN_SHVC

! Replace friendly
!-----------------

         if(FRIENDLY) then
            S = SX
         end if

! Update surface fluxes
! ---------------------

         if(associated(SF).and.associated(DSG)) then
            SF = SF + DSG*SDF
         end if

      enddo TRACERS

! End loop over all quantities to be diffused
!--------------------------------------------

      deallocate(KK)

      RETURN_(ESMF_SUCCESS)
    end subroutine UPDATE

  end subroutine RUN2


!*********************************************************************
!*********************************************************************
!*********************************************************************

!*********************************************************************

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine PRELIMS(IRUN,LM,DT,  &
                      ! Inputs
                      T,     &
                      QV,    &
                      PHALF, &
                      TH,    &
                      QLCN,  &
                      QLLS,  &
                      QICN,  &
                      QILS,  &
                      !CLCN,  & ! Currently not used in REFRESHKS
                      !CLLS,  & ! Currently not used in REFRESHKS
                      ! Outputs
                      ZFULL, &
                      ZHALF, &
                      TV,    &
                      PV,    &
                      RDZ,   &
                      DMI,   &
                      !QA,    & ! Currently not used in REFRESHKS
                      PFULL  )

#ifdef _CUDA
      integer, value :: IRUN,LM
      real,    value :: DT
#else
      integer, intent(IN   )                       :: IRUN,LM
      real,    intent(IN   )                       :: DT
#endif

      real,    intent(IN   ), dimension(IRUN,LM)   :: T
      real,    intent(IN   ), dimension(IRUN,LM)   :: QV
      real,    intent(IN   ), dimension(IRUN,LM+1) :: PHALF
      real,    intent(IN   ), dimension(IRUN,LM)   :: TH
      real,    intent(IN   ), dimension(IRUN,LM)   :: QLCN
      real,    intent(IN   ), dimension(IRUN,LM)   :: QLLS
      real,    intent(IN   ), dimension(IRUN,LM)   :: QICN
      real,    intent(IN   ), dimension(IRUN,LM)   :: QILS
      !real,    intent(IN   ), dimension(IRUN,LM)   :: CLCN ! Currently not used in REFRESHKS
      !real,    intent(IN   ), dimension(IRUN,LM)   :: CLLS ! Currently not used in REFRESHKS

      real,    intent(  OUT), dimension(IRUN,  LM  ) :: ZFULL
      real,    intent(  OUT), dimension(IRUN,  LM+1) :: ZHALF
      real,    intent(  OUT), dimension(IRUN,  LM  ) :: TV
      real,    intent(  OUT), dimension(IRUN,  LM  ) :: PV
      real,    intent(  OUT), dimension(IRUN,1:LM-1) :: RDZ
      real,    intent(  OUT), dimension(IRUN,  LM  ) :: DMI
      real,    intent(  OUT), dimension(IRUN,  LM  ) :: PFULL
      !real,    intent(  OUT), dimension(IRUN,  LM  ) :: QA ! Currently not used in REFRESHKS

      integer :: I,L

      real :: PKE_atL, PKE_atLp1
      real :: TVE
      real :: QS_dummy
      real :: QL_tot, QI_tot

#ifdef _CUDA
      I = (blockidx%x - 1) * blockdim%x + threadidx%x

      if ( I <= irun ) then
#else
      do I = 1, IRUN
#endif
         ! Compute the edge heights using Arakawa-Suarez hydrostatic equation
         !---------------------------------------------------------------------------

         ZHALF(I,LM+1) = 0.0
         do L = LM, 1, -1
            PKE_atLp1  = (PHALF(I,L+1)/MAPL_P00)**MAPL_KAPPA
            PKE_atL    = (PHALF(I,L  )/MAPL_P00)**MAPL_KAPPA

            ZHALF(I,L) = ZHALF(I,L+1) + (MAPL_CP/MAPL_GRAV)*TH(I,L)*(PKE_atLp1-PKE_atL)
         end do

         ! Layer height, pressure, and virtual temperatures
         !-------------------------------------------------

         do L = 1, LM
            QL_tot = QLCN(I,L) + QLLS(I,L)
            QI_tot = QICN(I,L) + QILS(I,L)
            !QA(I,J) = CLCN(I,L) + CLLS(I,L) ! Currently not used in REFRESHKS

            ZFULL(I,L) = 0.5*(ZHALF(I,L)+ZHALF(I,L+1))
            PFULL(I,L) = 0.5*(PHALF(I,L)+PHALF(I,L+1))

            TV(I,L)  = T(I,L) *( 1.0 + MAPL_VIREPS *QV(I,L) - QL_tot - QI_tot ) 
            PV(I,L) = TV(I,L)*(TH(I,L)/T(I,L))
         end do

         do L = 1, LM-1
            TVE = (TV(I,L) + TV(I,L+1))*0.5

            ! Miscellaneous factors
            !----------------------

            RDZ(I,L) = PHALF(I,L+1) / ( MAPL_RGAS * TVE )
            RDZ(I,L) = RDZ(I,L) / (ZFULL(I,L)-ZFULL(I,L+1))
         end do

         do L = 1, LM
            DMI(I,L) = (MAPL_GRAV*DT)/(PHALF(I,L+1)-PHALF(I,L))
         end do

         !===> Running 1-2-1 smooth of bottom 5 levels of Virtual Pot. Temp.
         PV(I,LM  ) = PV(I,LM-1)*0.25 + PV(I,LM  )*0.75
         PV(I,LM-1) = PV(I,LM-2)*0.25 + PV(I,LM-1)*0.50 + PV(I,LM  )*0.25 
         PV(I,LM-2) = PV(I,LM-3)*0.25 + PV(I,LM-2)*0.50 + PV(I,LM-1)*0.25 
         PV(I,LM-3) = PV(I,LM-4)*0.25 + PV(I,LM-3)*0.50 + PV(I,LM-2)*0.25 
         PV(I,LM-4) = PV(I,LM-5)*0.25 + PV(I,LM-4)*0.50 + PV(I,LM-3)*0.25 
         PV(I,LM-5) = PV(I,LM-6)*0.25 + PV(I,LM-5)*0.50 + PV(I,LM-4)*0.25 
#ifndef _CUDA
      end do 
#else
      end if
#endif

   end subroutine PRELIMS

!*********************************************************************

!BOP

! !IROUTINE:  GETKS -- Computes atmospheric diffusivities at interior levels

! !INTERFACE:

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine LOUIS_KS(IRUN,LM,     &
         KH,KM,RI,DU,ZPBL,          &
         ZZ,ZE,PV,UU,VV,            &
         LOUIS, MINSHEAR, MINTHICK, &
         LAMBDAM, LAMBDAM2,         &
         LAMBDAH, LAMBDAH2,         &
         ZKMENV, ZKHENV,            &
         AKHMMAX,                   &
         ALH_DIAG,KMLS_DIAG,KHLS_DIAG)

! !ARGUMENTS:

#ifdef _CUDA
      integer, value :: IRUN,LM
#else
      integer, intent(IN   ) :: IRUN,LM
#endif
      real,    intent(IN   ) :: LOUIS         ! Louis scheme parameters (usually 5).
      real,    intent(IN   ) :: MINSHEAR      ! Min shear allowed in Ri calculation (s-1).
      real,    intent(IN   ) :: MINTHICK      ! Min layer thickness (m).
      real,    intent(IN   ) :: AKHMMAX       ! Maximum allowe diffusivity (m+2 s-1).
      real,    intent(IN   ) :: LAMBDAM       ! Blackadar(1962) length scale parameter for momentum (m).
      real,    intent(IN   ) :: LAMBDAM2      ! Second Blackadar parameter for momentum (m).
      real,    intent(IN   ) :: LAMBDAH       ! Blackadar(1962) length scale parameter for heat (m).
      real,    intent(IN   ) :: LAMBDAH2      ! Second Blackadar parameter for heat (m).
      real,    intent(IN   ) :: ZKMENV        ! Transition height for Blackadar param for momentum (m)
      real,    intent(IN   ) :: ZKHENV        ! Transition height for Blackadar param for heat     (m)
      real,    intent(IN   ) :: ZZ(IRUN,LM)   ! Height of layer center above the surface (m).
      real,    intent(IN   ) :: PV(IRUN,LM)   ! Virtual potential temperature at layer center (K).
      real,    intent(IN   ) :: UU(IRUN,LM)   ! Eastward velocity at layer center (m s-1).
      real,    intent(IN   ) :: VV(IRUN,LM)   ! Northward velocity at layer center (m s-1).
      real,    intent(IN   ) :: ZE(IRUN,LM+1) ! Height of layer base above the surface (m).

      ! These are 1:LM+1 here but 0:LM in the GC
      ! Old code only passed in 1:LM-1 from GC which is 2:LM here.
      real,    intent(  OUT) :: KM(IRUN,LM+1) ! Momentum diffusivity at base of each layer (m+2 s-1).
      real,    intent(  OUT) :: KH(IRUN,LM+1) ! Heat diffusivity at base of each layer  (m+2 s-1).
      real,    intent(  OUT) :: RI(IRUN,LM+1) ! Richardson number
      real,    intent(  OUT) :: DU(IRUN,LM+1) ! Magnitude of wind shear (s-1).
      real,    intent(IN   ) :: ZPBL(IRUN)    ! PBL Depth (m)
   
      real,    intent(  OUT) :: ALH_DIAG(IRUN,LM+1)  ! Blackadar Length Scale diagnostic (m) [Optional] 
      real,    intent(  OUT) :: KMLS_DIAG(IRUN,LM+1) ! Momentum diffusivity at base of each layer (m+2 s-1).
      real,    intent(  OUT) :: KHLS_DIAG(IRUN,LM+1) ! Heat diffusivity at base of each layer  (m+2 s-1).

! !DESCRIPTION: Computes Louis et al.(1979) Richardson-number-based diffusivites,
!                as well as an additional ``entrainment'' diffusivity.
!                The Louis diffusivities for momentum, $K_m$, and for heat
!   and moisture, $K_h$, are defined at the interior layer edges. For LM layers,
!   we define diffusivities at the base of the top LM-1 layers. All indexing
!   is from top to bottom of the atmosphere. 
!
!
!  The Richardson number, Ri, is defined at the same edges as the diffusivities. 
!  $$
!  {\rm Ri}_l = \frac{ \frac{g}{\left(\overline{\theta_v}\right)_l}\left(\frac{\delta \theta_v}{\delta z}\right)_l }
!                    { \left(\frac{\delta {\bf |V|}}{\delta z}\right)^2_l             }, \, \,  l=1,LM-1
!  $$
!  where $\theta_v=\theta(1+\epsilon q)$ is the virtual potential temperature,
!  $\epsilon=\frac{M_a}{M_w}-1$, $M_a$ and $M_w$ are the molecular weights of
!  dry air and water, and $q$ is the specific humidity.
!  $\delta \theta_v$ is the difference of $\theta_v$ in the layers above and below the edge 
!  at which Ri$_l$ is defined; $\overline{\theta_v}$ is their average.
!
!  The diffusivities at the layer edges have the form:
!  $$
!  K^m_l = (\ell^2_m)_l \left(\frac{\delta {\bf |V|}}{\delta z}\right)_l f_m({\rm Ri}_l)
!  $$
!  and
!  $$
!  K^h_l = (\ell^2_h)_l \left(\frac{\delta {\bf |V|}}{\delta z}\right)_l f_h({\rm Ri}_l),
!  $$
!  where $k$ is the Von Karman constant, and $\ell$ is the 
!  Blackdar(1962) length scale, also defined at the layer edges.
!
!  Different turbulent length scales can be used for heat and momentum. 
!  in both cases, we use the  traditional formulation:
!  $$
!  (\ell_{(m,h)})_l  = \frac{kz_l}{1 + \frac{kz_l}{\lambda_{(m,h)}}},
!  $$
!  where, near the surface, the scale is proportional to $z_l$, the height above 
!  the surface of edge level $l$, and far from the surface it approaches $\lambda$.
!  The length scale $\lambda$ is usually taken to be a constant (order 150 m), assuming
!  the same scale for the outre boundary layer and the free atmosphere. We make it
!  a function of height, reducing its value in the free atmosphere. The momentum
!  length scale written as:
!  $$
!  \lambda_m = \max(\lambda_1 e^{\left(\frac{z_l}{z_T}\right)^2}, \lambda_2)
!  $$
!  where $\lambda_2 \le \lambda_1$ and $z_T$ is the top of the boundary layer.
!  The length scale for heat and other scalers is taken as: $\lambda_h =  \sqrt\frac{3d}{2} \lambda_m$,
!  following the scheme used at ECMWF.
!
!  The two universal functions of the Richardson number,  $f_m$ and $f_h$,
!  are taken from Louis et al (1982). For unstable conditions (Ri$\le 0$),
!  they are:
!  $$
!  f_m = (1 - 2b \psi)
!  $$
!  and
!  $$
!  f_h = (1 - 3b \psi),
!  $$
!  where
!  $$
!  \psi = \frac{ {\rm Ri} }{ 1 + 3bC(z)\sqrt{-{\rm Ri}} },
!  $$
!  and
!  $$
!  C(z)=
!  $$

!  For stable condition (Ri$\ge 0$), they are
!  $$
!  f_m = \frac{1}{1.0 + \frac{2b{\rm Ri}}{\psi}}
!  $$
!  and
!  $$
!  f_h = \frac{1}{1.0 + 3b{\rm Ri}\psi},
!  $$
!  where
!  $$
!  \psi =  \sqrt{1+d{\rm Ri}}.
!  $$
!  As in Louis et al (1982), the parameters appearing in these are taken  
!  as $b = c = d = 5$. 


!EOP

! Locals

      real    :: ALH, ALM, DZ, DT, TM, PS, LAMBDAM_X, LAMBDAH_X
      integer :: L,Levs,i,j,k
      !real    :: Zchoke 
      real    :: pbllocal
      real    :: alhfac,almfac

! Begin...
!===>   Number of layers; edge levels will be one less (LM-1).

      !LM = size(ZZ,3)

      almfac = 1.2
      alhfac = 1.2

#ifdef _CUDA
      i = (blockidx%x - 1) * blockdim%x + threadidx%x

      if ( i <= irun ) then
#else
      do i = 1, irun
#endif
         ALH_DIAG(i,   1) = 0.0
         ALH_DIAG(i,LM+1) = 0.0

         KH(i,   1) = 0.0
         KH(i,LM+1) = 0.0
         KM(i,   1) = 0.0
         KM(i,LM+1) = 0.0
         DU(i,   1) = 0.0
         DU(i,LM+1) = 0.0
         RI(i,   1) = 0.0
         RI(i,LM+1) = 0.0
         KMLS_DIAG(i,   1) = 0.0
         KMLS_DIAG(i,LM+1) = 0.0
         KHLS_DIAG(i,   1) = 0.0
         KHLS_DIAG(i,LM+1) = 0.0

!MAT>   Initialize pbllocal

         pbllocal = ZPBL(i)

         if ( pbllocal .LE. ZZ(i,LM) ) pbllocal = ZZ(i,LM)

!===>   Quantities needed for Richardson number

         do k = 2, lm ! This is equivalent to 1->LM-1 but we move all 0:72 elements
            KH(i,k) = 0.0
            KM(i,k) = 0.0
            DU(i,k) = 0.0
            RI(i,k) = 0.0

            DZ      = (ZZ(i,k-1) - ZZ(i,k))
            TM      = (PV(i,k-1) + PV(i,k))*0.5
            DT      = (PV(i,k-1) - PV(i,k))
            DU(i,k) = (UU(i,k-1) - UU(i,k))**2 + &
                      (VV(i,k-1) - VV(i,k))**2

!===>   Limits on distance between layer centers and vertical shear at edges.

            DZ      = max(DZ, MINTHICK)
            DU(i,k) = sqrt(DU(i,k))/DZ

!===>   Richardson number  ( RI = G*(DTheta_v/DZ) / (Theta_v*|DV/DZ|^2) )

            RI(i,k) = MAPL_GRAV*(DT/DZ)/(TM*( max(DU(i,k), MINSHEAR)**2))

!===>   Blackadar(1962) length scale: $1/l = 1/(kz) + 1/\lambda$

!!! LAMBDAM_X = MAX( LAMBDAM * EXP( -(ZE / ZKMENV )**2 ) , LAMBDAM2 )
!!! LAMBDAH_X = MAX( LAMBDAH * EXP( -(ZE / ZKHENV )**2 ) , LAMBDAH2 )

            LAMBDAM_X = MAX( 0.1 * PBLLOCAL * EXP( -(ZE(i,k) / ZKMENV )**2 ) , LAMBDAM2 )
            LAMBDAH_X = MAX( 0.1 * PBLLOCAL * EXP( -(ZE(i,k) / ZKHENV )**2 ) , LAMBDAH2 )

            ALM = almfac * ( MAPL_KARMAN*ZE(i,k)/( 1.0 + MAPL_KARMAN*(ZE(i,k)/LAMBDAM_X) ) )**2
            ALH = alhfac * ( MAPL_KARMAN*ZE(i,k)/( 1.0 + MAPL_KARMAN*(ZE(i,k)/LAMBDAH_X) ) )**2

            ALH_DIAG(i,k) = SQRT( ALH )

!===>   Unstable case: Uses (3.14, 3.18, 3.27) in Louis-scheme
!                      should approach (3.13) for small -RI.

            if ( RI(i,k)  < 0.0 ) then
               PS = ( (ZZ(i,k-1)/ZZ(i,k))**(1./3.) - 1.0 ) ** 3
               PS = ALH*sqrt( PS/(ZE(i,k)*(DZ**3)) )
               PS = RI(i,k)/(1.0 + (3.0*LOUIS*LOUIS)*PS*sqrt(-RI(i,k)))
      
               KH(i,k) = 1.0 - (LOUIS*3.0)*PS
               KM(i,k) = 1.0 - (LOUIS*2.0)*PS
            end if

!===>   Choke off unstable KH below Zchoke (m). JTB 2/2/06    
!!! Zchoke = 500.
!!! where( (RI < 0.0) .and. (ZE < Zchoke ) )  
!!!    KH = KH * (( ZE / Zchoke )**3)            
!!! endwhere

!===>   Stable case

            if ( RI(i,k) >= 0.0 ) then
               PS      = sqrt  (1.0 +  LOUIS     *RI(i,k)   )
         
               KH(i,k) = 1.0 / (1.0 + (LOUIS*3.0)*RI(i,k)*PS)
               KM(i,k) = PS  / (PS  + (LOUIS*2.0)*RI(i,k)   )
            end if

!===>   DIMENSIONALIZE Kz and  LIMIT DIFFUSIVITY

            ALM = DU(i,k)*ALM
            ALH = DU(i,k)*ALH

            KM(i,k) = min(KM(i,k)*ALM, AKHMMAX)
            KH(i,k) = min(KH(i,k)*ALH, AKHMMAX)
            KMLS_DIAG(i,k) = KM(i,k)
            KHLS_DIAG(i,k) = KH(i,k)
         end do
#ifndef _CUDA
      end do 
#else
      end if
#endif

  end subroutine LOUIS_KS

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine BELJAARS(IRUN,LM,DT,     &
                       LAMBDA_B, &
                       C_B,      &
                       FKV,            &
                       BKV,            &
                       U,              &
                       V,              &
                       ZFULL,          &
                       VARFLT,         &
                       PHALF           )

!BOP
!
!   Orograpghic drag follows  Beljaars (2003):
!   $$
!   \frac{\partial}{\partial z}\frac{\tau}{\rho} = \frac{C_B}{\lambda_B} |U(z)| U(z) 
!          e^{-\tilde{z}^\frac{3}{2}}\tilde{z}^{-1.2},
!   $$
!   where $z$ is the height above the surface in meters, 
!   $\tilde{z}=\frac{z}{\lambda_B}$, $\tau$ is the orographic stress at $z$,
!   $\rho$ is the air density, $U(z)$ is the wind velocity, and $\lambda_B$ is a vertical length scale.
!   Beljaars uses $\lambda_B = 1500$m, for which the non-dimensional parameter $C_B = 2.5101471 \times 10^{-8}$.
!   These are the default values, but both can be modified from the configuration. To avoid underflow.
!   the tendency is set to zero once $\tilde{z}$ exceeds 4 (i.e., 6 km from the surface for default values). 
!
!EOP

#ifdef _CUDA
      integer, value :: IRUN,LM
      real,    value :: DT
#else
      integer, intent(IN   )                       :: IRUN,LM
      real,    intent(IN   )                       :: DT
#endif
      real,    intent(IN   )                       :: LAMBDA_B
      real,    intent(IN   )                       :: C_B

      real,    intent(  OUT), dimension(IRUN,  LM) :: FKV
      real,    intent(INOUT), dimension(IRUN,  LM) :: BKV
      real,    intent(IN   ), dimension(IRUN,  LM) :: U
      real,    intent(IN   ), dimension(IRUN,  LM) :: V
      real,    intent(IN   ), dimension(IRUN,  LM) :: ZFULL
      real,    intent(IN   ), dimension(IRUN     ) :: VARFLT
      real,    intent(IN   ), dimension(IRUN,LM+1) :: PHALF

      integer :: I,L
      real    :: FKV_temp

#ifdef _CUDA
      I = (blockidx%x - 1) * blockdim%x + threadidx%x

      if ( I <= irun ) then
#else
      do I = 1, IRUN
#endif
         do L=LM,1,-1
            FKV(I,L) = 0.0

            if (ZFULL(I,L) < 4.0*LAMBDA_B) then
               FKV_temp = ZFULL(I,L)*(1.0/LAMBDA_B)
               FKV_temp = VARFLT(I) * exp(-FKV_temp*sqrt(FKV_temp))*(FKV_temp**(-1.2))
               FKV_temp = (C_B/LAMBDA_B)*min( sqrt(U(I,L)**2+V(I,L)**2),5.0 )*FKV_temp

               BKV(I,L) = BKV(I,L) + DT*FKV_temp
               FKV(I,L) = FKV_temp * (PHALF(I,L+1)-PHALF(I,L))
            end if
         end do
#ifndef _CUDA
      end do 
#else
      end if
#endif

   end subroutine BELJAARS

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine POSTLOCK(IRUN,LM,DT,    &
                       ! Constants
                       KPBLMIN,       &
                       CALC_TCZPBL,   &
                       CALC_ZPBL2,   &
                       CALC_ZPBL10p,   &
                       ! Inputs
                       ZFULL, ZEDGE,  &
                       PFULL,         &
                       RDZ,           &
                       DMI,           &
                       PHALF,         &
                       TV,            &
                       CT,            &
                       CQ,            &
                       CU,            &
                       T,             &
                       Q,             &
                       TH,            &
                       U,             &
                       V,             &
                       DU,            &
                       RI,            &
                       PBLHT_OPTION,  &
                       ! Inoutputs
                       DIFF_T,        &
                       DIFF_M,        &
                       ! Outputs
                       AKQ, AKS, AKV, &
                       BKQ, BKS, BKV, &
                       CKQ, CKS, CKV, &
                                 EKV, &
                       ZPBL,          &
                       ! Diagnostics
                       AKQODT_DIAG, AKSODT_DIAG, AKVODT_DIAG, &
                       CKQODT_DIAG, CKSODT_DIAG, CKVODT_DIAG, &
                       PPBL_DIAG,   KPBL, TCZPBL_DIAG, ZPBL2_DIAG,  &
                       ZPBL10p_DIAG, ZPBLHTKE_DIAG, ZPBLRI_DIAG, &
                       ZPBLRI2_DIAG, ZPBLTHV_DIAG)

#ifdef _CUDA
      integer, value :: IRUN,LM
      real,    value :: DT
      integer, value :: PBLHT_OPTION
      logical, value :: CALC_TCZPBL
      logical, value :: CALC_ZPBL2
      logical, value :: CALC_ZPBL10p
#else
      integer, intent(IN   )                   :: IRUN,LM
      real,    intent(IN   )                   :: DT
      integer, intent(IN   )                   :: PBLHT_OPTION
      logical, intent(IN   )                   :: CALC_TCZPBL
      logical, intent(IN   )                   :: CALC_ZPBL2
      logical, intent(IN   )                   :: CALC_ZPBL10p
#endif

      integer, intent(IN   )                   :: KPBLMIN

      real,    intent(IN   ), dimension(IRUN,LM) :: ZFULL
      real,    intent(IN   ), dimension(IRUN,LM) :: ZEDGE
      real,    intent(IN   ), dimension(IRUN,LM) :: PFULL
      real,    intent(IN   ), dimension(IRUN,LM-1) :: RDZ
      real,    intent(IN   ), dimension(IRUN,LM) :: DMI
      real,    intent(IN   ), dimension(IRUN,LM+1) :: PHALF
      real,    intent(IN   ), dimension(IRUN,LM) :: TV
      real,    intent(IN   ), dimension(IRUN  ) :: CT
      real,    intent(IN   ), dimension(IRUN  ) :: CQ
      real,    intent(IN   ), dimension(IRUN  ) :: CU
      real,    intent(IN   ), dimension(IRUN,LM) :: T,Q,TH,U,V
      real,    intent(IN   ), dimension(IRUN,LM+1) :: DU, RI

      real,    intent(INOUT), dimension(IRUN,LM+1) :: DIFF_T
      real,    intent(INOUT), dimension(IRUN,LM+1) :: DIFF_M

      real,    intent(  OUT), dimension(IRUN,LM) :: AKQ, AKS, AKV
      real,    intent(  OUT), dimension(IRUN,LM) :: BKQ, BKS, BKV
      real,    intent(  OUT), dimension(IRUN,LM) :: CKQ, CKS, CKV
      real,    intent(  OUT), dimension(IRUN,LM) ::           EKV

      real,    intent(  OUT), dimension(IRUN,LM) :: AKQODT_DIAG, AKSODT_DIAG, AKVODT_DIAG
      real,    intent(  OUT), dimension(IRUN,LM) :: CKQODT_DIAG, CKSODT_DIAG, CKVODT_DIAG

      real,    intent(  OUT), dimension(IRUN   ) :: TCZPBL_DIAG
      real,    intent(  OUT), dimension(IRUN   ) :: ZPBL2_DIAG, ZPBL10p_DIAG, ZPBLRI2_DIAG
      real,    intent(  OUT), dimension(IRUN   ) :: ZPBLHTKE_DIAG, ZPBLRI_DIAG, ZPBLTHV_DIAG
      real,    intent(  OUT), dimension(IRUN   ) :: ZPBL, PPBL_DIAG, KPBL

! CUDA must know the size of all GPU arrays a priori
#ifndef _CUDA
#define GPU_MAXLEVS LM
#endif

      real,             dimension(GPU_MAXLEVS+1) :: temparray, htke
      real,             dimension(GPU_MAXLEVS  ) :: tcrib !TransCom bulk Ri
      real,             dimension(GPU_MAXLEVS+1) :: thetav

      integer :: I,L,locmax
      real    :: maxkh,maxhtke,minlval,thetavs,thetavh,uv2h,kpbltc,kpbl2,kpbl10p
      real    :: maxdthvdz,dthvdz

! PBL-top diagnostic
! -----------------------------------------

      real, parameter :: tcri_crit = 0.25
      real, parameter :: ri_crit = 0.00
      real, parameter :: ri_crit2 = 0.20

#ifdef _CUDA
      I = (blockidx%x - 1) * blockdim%x + threadidx%x

      if ( I <= irun ) then
#else
      do I = 1, IRUN
#endif
         ZPBL(I) = MAPL_UNDEF
         PPBL_DIAG(I) = MAPL_UNDEF

         if(CALC_TCZPBL) then
            TCZPBL_DIAG(I) = MAPL_UNDEF
            thetavs = T(I,LM)*(1.0+MAPL_VIREPS*Q(I,LM)/(1.0-Q(I,LM)))*(TH(I,LM)/T(I,LM))
            do L=LM-1,1,-1
               thetavh = T(I,L)*(1.0+MAPL_VIREPS*Q(I,L)/(1.0-Q(I,L)))*(TH(I,L)/T(I,L))
               uv2h = max(U(I,L)**2+V(I,L)**2,1.0E-8)
               tcrib(L) = MAPL_GRAV*(thetavh-thetavs)*ZFULL(I,L)/(thetavs*uv2h)
               if (tcrib(L) >= tcri_crit) then
                 TCZPBL_DIAG(I) = ZFULL(I,L+1)+(tcri_crit-tcrib(L+1))/(tcrib(L)-tcrib(L+1))*(ZFULL(I,L)-ZFULL(I,L+1))
                 KPBLTC = float(L)
                 exit
               endif
            enddo
            if(TCZPBL_DIAG(I)<0.) then
             TCZPBL_DIAG(I) = ZFULL(I,LM)
             KPBLTC = float(LM)
            endif
         endif

         if(CALC_ZPBL2) then
            ZPBL2_DIAG(I) = MAPL_UNDEF
            do L=LM,2,-1
               if ((DIFF_T(I,L) < 2.).and.(DIFF_T(I,L+1) >= 2.).and.(ZPBL2_DIAG(I)==MAPL_UNDEF)) then
                  ZPBL2_DIAG(I) = ZFULL(I,L)
                  KPBL2 = float(L)
               end if
            end do
            if (  ZPBL2_DIAG(I) .eq. MAPL_UNDEF ) then
             ZPBL2_DIAG(I) = ZFULL(I,LM)
             KPBL2 = float(LM)
            endif
            ZPBL2_DIAG(I) = MIN(ZPBL2_DIAG(I),ZFULL(I,KPBLMIN))
         endif

         if(CALC_ZPBL10p) then
            ZPBL10p_DIAG(I) = MAPL_UNDEF
            temparray = DIFF_T(I,:)
            do L = LM,2,-1
               ! PGI has an ICE with CUDA and maxloc, so we call our own 
               ! function (which is probably inefficient, but is 
               ! equivalent to the intrinsic, math-wise)
#ifdef _CUDA
               locmax = gpu_maxloc(temparray,LM+1)
#else
               locmax = maxloc(temparray,1)
#endif
               minlval = max(0.001,0.0001*maxval(temparray))
               if(temparray(locmax-1)<minlval.and.temparray(locmax+1)<minlval) temparray(locmax) = minlval
            enddo
            maxkh = temparray(LM)
            do L = LM-1,2,-1
               if(temparray(L)>maxkh) maxkh = temparray(L)
               if(temparray(L-1)<minlval) exit
            end do
            do L=LM-1,2,-1
               if ( (temparray(L) < 0.1*maxkh) .and. (temparray(L+1) >= 0.1*maxkh)  &
                 .and. (ZPBL10p_DIAG(I) == MAPL_UNDEF ) ) then
                  ZPBL10p_DIAG(I) = ZEDGE(I,L+1)+ &
                ((ZEDGE(I,L)-ZEDGE(I,L+1))/(temparray(L)-temparray(L+1))) * (0.1*maxkh-temparray(L+1))
                  KPBL10p = float(L)
               end if
            end do
            if (  ZPBL10p_DIAG(I) .eq. MAPL_UNDEF .or. (maxkh.lt.1.)) then
             ZPBL10p_DIAG(I) = ZFULL(I,LM)
             KPBL10p = float(LM)
            endif
            ZPBL10p_DIAG(I) = MIN(ZPBL10p_DIAG(I),ZFULL(I,KPBLMIN))
         end if

! HTKE pbl height
            ZPBLHTKE_DIAG(I) = MAPL_UNDEF
            do L=LM+1,1,-1
               HTKE(L) = DU(I,L)*DIFF_M(I,L)
            end do
            temparray = HTKE(:)
            maxhtke = temparray(LM)
            do L = LM,1,-1
#ifdef _CUDA
               locmax = gpu_maxloc(temparray,LM+1)
#else
               locmax = maxloc(temparray,1)
#endif
              minlval = max(1E-6,0.001*maxval(temparray))
              if(temparray(locmax-1)<minlval.and.temparray(locmax+1)<minlval) temparray(locmax) = minlval
              if(temparray(L)>maxhtke) maxhtke = HTKE(L)
              if(temparray(L-1)<minlval) exit
            enddo
            do L=LM-1,1,-1
              if((HTKE(L)<0.1*maxhtke).and.(HTKE(L+1)>=0.1*maxhtke) &
                .and.(ZPBLHTKE_DIAG(I)==MAPL_UNDEF)) then
                ZPBLHTKE_DIAG(I) = ZFULL(I,L+1)+(0.1*maxhtke-temparray(L+1))/ &
                    (temparray(L)-temparray(L+1))*(ZFULL(I,L)-ZFULL(I,L+1))
               end if
            end do
            if ( ZPBLHTKE_DIAG(I) .eq. MAPL_UNDEF ) ZPBLHTKE_DIAG(I) = ZFULL(I,LM)
            if ( ZPBLHTKE_DIAG(I) < 0 ) ZPBLHTKE_DIAG = ZFULL(I,LM)
            ZPBLHTKE_DIAG(I) = MIN(ZPBLHTKE_DIAG(I),ZFULL(I,KPBLMIN))

! RI local diagnostic for pbl height thresh 0.
            ZPBLRI_DIAG(I) = MAPL_UNDEF
            if(RI(I,LM)>ri_crit) ZPBLRI_DIAG(I) = ZFULL(I,LM)
            do L=LM-1,1,-1
               if( (RI(I,L)>ri_crit) .and. (ZPBLRI_DIAG(I) == MAPL_UNDEF) ) then
                ZPBLRI_DIAG(I) = ZFULL(I,L+1)+(ri_crit-RI(I,L+1))/(RI(I,L)-RI(I,L+1))*(ZFULL(I,L)-ZFULL(I,L+1))
              endif
            enddo
            if (  ZPBLRI_DIAG(I) .eq. MAPL_UNDEF ) ZPBLRI_DIAG(I) = ZFULL(I,LM)
            ZPBLRI_DIAG(I) = MIN(ZPBLRI_DIAG(I),ZFULL(I,KPBLMIN))
            if(ZPBLRI_DIAG(I)<0) ZPBLRI_DIAG(I) = ZFULL(I,LM)

! RI local diagnostic for pbl height thresh 0.2
            ZPBLRI2_DIAG(I) = MAPL_UNDEF
            if(RI(I,LM)>ri_crit2) ZPBLRI2_DIAG(I) = ZFULL(I,LM)
            do L=LM-1,1,-1
               if( (RI(I,L)>ri_crit2) .and. (ZPBLRI2_DIAG(I) == MAPL_UNDEF) ) then
                ZPBLRI2_DIAG(I) = ZFULL(I,L+1)+(ri_crit2-RI(I,L+1))/(RI(I,L)-RI(I,L+1))*(ZFULL(I,L)-ZFULL(I,L+1))
              endif
            enddo
            if (  ZPBLRI2_DIAG(I) .eq. MAPL_UNDEF ) ZPBLRI2_DIAG(I) = ZFULL(I,LM)
            ZPBLRI2_DIAG(I) = MIN(ZPBLRI2_DIAG(I),ZFULL(I,KPBLMIN))
            if(ZPBLRI2_DIAG(I)<0) ZPBLRI2_DIAG(I) = ZFULL(I,LM)

! thetav gradient based pbl height diagnostic
            ZPBLTHV_DIAG(I) = MAPL_UNDEF
            do L=LM,1,-1
               thetav(L) = T(I,L)*(1.0*MAPL_VIREPS*Q(I,L)/(1.0-Q(I,L)))*(TH(I,L)/T(I,L))
            end do
            maxdthvdz = 0
            do L=LM-1,1,-1
              if(ZFULL(I,L)<=ZFULL(I,KPBLMIN)) then
                 dthvdz = (thetav(L+1)-thetav(L))/(ZFULL(I,L+1)-ZFULL(I,L))
                 if(dthvdz>maxdthvdz) then
                    maxdthvdz = dthvdz
                    ZPBLTHV_DIAG(I) = 0.5*(ZFULL(I,L+1)+ZFULL(I,L))
                 endif
              endif
            end do

         SELECT CASE( PBLHT_OPTION)

         CASE( 1 )
            ZPBL(I) = ZPBL2_DIAG(I)
            KPBL(I) = kpbl2
            ZPBL(I) = MIN(ZPBL(I),ZFULL(I,KPBLMIN))
            KPBL(I) = MAX(KPBL(I),float(KPBLMIN))
            PPBL_DIAG(I) = PFULL(I,KPBL(I))
            PPBL_DIAG(I) = MAX(PPBL_DIAG(I),PFULL(I,KPBLMIN))

         CASE( 2 )
            ZPBL(I) = ZPBL10p_DIAG(I)
            KPBL(I) = kpbl10p
            ZPBL(I) = MIN(ZPBL(I),ZFULL(I,KPBLMIN))
            KPBL(I) = MAX(KPBL(I),float(KPBLMIN))
            PPBL_DIAG(I) = PFULL(I,KPBL(I))
            PPBL_DIAG(I) = MAX(PPBL_DIAG(I),PFULL(I,KPBLMIN))

         CASE( 3 )
            ZPBL(I) = TCZPBL_DIAG(I)
            KPBL(I) = kpbltc
            ZPBL(I) = MIN(ZPBL(I),ZFULL(I,KPBLMIN))
            KPBL(I) = MAX(KPBL(I),float(KPBLMIN))
            PPBL_DIAG(I) = PFULL(I,KPBL(I))
            PPBL_DIAG(I) = MAX(PPBL_DIAG(I),PFULL(I,KPBLMIN))

         END SELECT


! Second difference coefficients for scalars; RDZ is RHO/DZ, DMI is (G DT)/DP
! ---------------------------------------------------------------------------

         do L = 1, LM
            if (L <= LM-1) then
               CKS(I,L) = -DIFF_T(I,L+1) * RDZ(I,L)
               if (L == 1) AKS(I,L) = 0.0
               AKS(I,L+1) = CKS(I,L) * DMI(I,L+1)
               CKS(I,L) = CKS(I,L) * DMI(I,L)
            end if
            if (L == LM) CKS(I,L) = -CT(I) * DMI(I,L)

! Fill DIFF_T at level LM+1 with CT * RDZ for diagnostic output
            if (L == LM) DIFF_T(I,L+1) = CT(I) * (PHALF(I,L+1)/(MAPL_RGAS * TV(I,L))) / ZFULL(I,L)

! Water vapor can differ at the surface
!--------------------------------------

            CKQ(I,L) = CKS(I,L)
            AKQ(I,L) = AKS(I,L)
            if (L == LM) CKQ(I,L) = -CQ(I) * DMI(I,L)

! Second difference coefficients for winds
! EKV is saved to use in the frictional heating calc.
! ---------------------------------------------------
   
            if (L <= LM-1) then
               EKV(I,L) = -DIFF_M(I,L+1) * RDZ(I,L)
               if (L == 1) AKV(I,L) = 0.0
               AKV(I,L+1) = EKV(I,L) * DMI(I,L+1)
               CKV(I,L) = EKV(I,L) * DMI(I,L)
               EKV(I,L) = -MAPL_GRAV * EKV(I,L)
            end if

            if (L == LM) then
               CKV(I,L) = -  CU(I) * DMI(I,L)
               EKV(I,L) =  MAPL_GRAV *  CU(I)
            end if

! Fill DIFF_M at level LM+1 with CU * RDZ for diagnostic output
            if (L == LM) DIFF_M(I,L+1) = CU(I) * (PHALF(I,L+1)/(MAPL_RGAS * TV(I,L))) / ZFULL(I,L)

! Setup the tridiagonal matrix
! ----------------------------

            BKS(I,L) = 1.00 - (AKS(I,L)+CKS(I,L))
            BKQ(I,L) = 1.00 - (AKQ(I,L)+CKQ(I,L))
            BKV(I,L) = 1.00 - (AKV(I,L)+CKV(I,L))

! Add the topographic roughness term
! ----------------------------------

            if (L == 1) then
               AKSODT_DIAG(I,L) = 0.0
            else
               AKSODT_DIAG(I,L) = -AKS(I,L)/DT
            end if

            if (L == LM) then
               CKSODT_DIAG(I,L) = 0.0
            else
               CKSODT_DIAG(I,L) = -CKS(I,L)/DT
            end if

            if (L == 1) then
               AKQODT_DIAG(I,L) = 0.0
            else
               AKQODT_DIAG(I,L) = -AKQ(I,L)/DT
            end if

            if (L == LM) then
               CKQODT_DIAG(I,L) = 0.0
            else
               CKQODT_DIAG(I,L) = -CKQ(I,L)/DT
            end if

            AKVODT_DIAG(I,L) = -AKV(I,L)/DT

            CKVODT_DIAG(I,L) = -CKV(I,L)/DT
         end do
#ifndef _CUDA
      end do 
#else
      end if
#endif

   end subroutine POSTLOCK


!*********************************************************************

!BOP

! !IROUTINE:  VTRILU --  Does LU decomposition of tridiagonal matrix.

! !INTERFACE:

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine VTRILU(IRUN,LM,A,B,C)

! !ARGUMENTS:

#ifdef _CUDA
      integer, value :: IRUN, LM
#else
      integer,                     intent(IN   ) :: IRUN, LM
#endif
      real,    dimension(IRUN,LM), intent(IN   ) :: C
      real,    dimension(IRUN,LM), intent(INOUT) :: A, B

! !DESCRIPTION: {\tt VTRILU} performs an $LU$ decomposition on
! a tridiagonal matrix $M=LU$.
!
! $$
! M = \left( \begin{array}{ccccccc}
!      b_1 & c_1 & & & & & \\
!      a_2 & b_2 & c_2 & & & &  \\
!      &  \cdot& \cdot & \cdot & & &  \\
!      & & \cdot& \cdot & \cdot & &  \\
!      &&  & \cdot& \cdot & \cdot &  \\
!      &&&& a_{K-1} & b_{K-1} & c_{K-1}   \\
!      &&&&& a_{K} & b_{K}
!    \end{array} \right)
! $$
!
!
! $$
! \begin{array}{lr}
! L = \left( \begin{array}{ccccccc}
!      1 &&&&&& \\
!      \hat{a}_2 & 1 & &&&&  \\
!      &  \cdot& \cdot &  & & &  \\
!      & & \cdot& \cdot &  &&  \\
!      &&  & \cdot& \cdot &  &  \\
!      &&&& \hat{a}_{K-1} & 1 &   \\
!      &&&&& \hat{a}_{K} & 1
!    \end{array} \right)
! &
! U = \left( \begin{array}{ccccccc}
!      \hat{b}_1 & c_1 &&&&& \\
!       & \hat{b}_2 & c_2 &&&&  \\
!      &  & \cdot & \cdot & & &  \\
!      & & & \cdot & \cdot &&  \\
!      &&  & & \cdot & \cdot &  \\
!      &&&&  & \hat{b}_{K-1} & c_{K-1}   \\
!      &&&&&  & \hat{b}_{K}
!    \end{array} \right)
! \end{array}
! $$
!
!
! On input, A, B, and C contain, $a_k$, $b_k$, and $c_k$
! the lower, main, and upper diagonals of the matrix, respectively.
! On output, B contains $1/\hat{b}_k$, the inverse of the main diagonal of $U$,
! and A contains $\hat{a}_k$,
! the lower diagonal of $L$. C contains the upper diagonal of the original matrix and of $U$.
!
! The new diagonals $\hat{a}_k$ and $\hat{b}_k$ are:
! $$
! \begin{array}{rcl}
! \hat{b}_1 & = & b_1, \\
! \hat{a}_k & = & \makebox[2 in][l]{$a_k / \hat{b}_{k-1}$,}  k=2, K, \\
! \hat{b}_k & = & \makebox[2 in][l]{$b_k - c_{k-1} \hat{a}_k$,} k=2, K. 
! \end{array}
! $$
!EOP

      integer :: I, L

#ifdef _CUDA
      I = (blockidx%x - 1) * blockdim%x + threadidx%x

      if ( I <= irun ) then
#else
      do I = 1, IRUN
#endif
         B(I,1) = 1. / B(I,1)

         do L = 2,LM
            A(I,L) = A(I,L) * B(I,L-1)
            B(I,L) = 1. / ( B(I,L) - C(I,L-1) * A(I,L) )
         end do
#ifndef _CUDA
      end do 
#else
      end if
#endif

   end subroutine VTRILU

!*********************************************************************

!BOP

! !IROUTINE:  VTRISOLVESURF -- Solves for sensitivity to surface value


! !INTERFACE:

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine VTRISOLVESURF(IRUN,LM,B,C,Y)

! !ARGUMENTS:

#ifdef _CUDA
      integer, value :: IRUN, LM
#else
      integer,                     intent(IN   ) :: IRUN, LM
#endif
      real,    dimension(IRUN,LM), intent(IN   ) :: B, C
      real,    dimension(IRUN,LM), intent(  OUT) :: Y

! !DESCRIPTION: Solves tridiagonal system that has been LU decomposed
!   for the special case
!   where the surface Y (YG) is 1 and the rest of the input Ys are 0.
!   Everything else is as in {\tt VTRISOLVE}. This gives the sensitivity of the
!   solution to a unit change in the surface values.

!EOP

      integer :: I, L

#ifdef _CUDA
      I = (blockidx%x - 1) * blockdim%x + threadidx%x

      if ( I <= irun ) then
#else
      do I = 1, IRUN
#endif
         Y(I,LM) = -C(I,LM) * B(I,LM)

         do L = LM-1,1,-1
            Y(I,L) = -C(I,L) * Y(I,L+1) * B(I,L)
         end do
#ifndef _CUDA
      end do 
#else
      end if
#endif

   end subroutine VTRISOLVESURF

!BOP

! !IROUTINE:  VTRISOLVE -- Solves for tridiagonal system that has been decomposed by VTRILU


! !INTERFACE:

  subroutine VTRISOLVE ( A,B,C,Y,YG )

! !ARGUMENTS:

    real, dimension(:,:,:),  intent(IN   ) ::  A, B, C
    real, dimension(:,:,:),  intent(INOUT) ::  Y
    real, dimension(:,:),    intent(IN)    ::  YG

! !DESCRIPTION: Solves tridiagonal system that has been LU decomposed
!   $LU x = f$. This is done by first solving $L g = f$ for $g$, and 
!   then solving $U x = g$ for $x$. The solutions are:
! $$
! \begin{array}{rcl}
! g_1 & = & f_1, \\
! g_k & = & \makebox[2 in][l]{$f_k - g_{k-1} \hat{a}_{k}$,}  k=2, K, \\
! \end{array}
! $$
! and  
! $$
! \begin{array}{rcl}
! x_K & = & g_K /\hat{b}_K, \\
! x_k & = & \makebox[2 in][l]{($g_k - c_k g_{k+1}) / \hat{b}_{k}$,}  k=K-1, 1 \\
! \end{array}
! $$
!  
!  On input A contains the $\hat{a}_k$, the lower diagonal of $L$,
!   B contains the $1/\hat{b}_k$, inverse of the  main diagonal of $U$,
!   C contains the $c_k$, the upper diagonal of $U$. The forcing, $f_k$ is
!   
!   It returns the
!   solution in the r.h.s input vector, Y. A has the multiplier from the
!   decomposition, B the 
!   matrix (U), and C the upper diagonal of the original matrix and of U.
!   YG is the LM+1 (Ground) value of Y.

!EOP

    integer :: LM, L

    LM = size(Y,3)

! Sweep down, modifying rhs with multiplier A

    do L = 2,LM
       Y(:,:,L) = Y(:,:,L) - Y(:,:,L-1) * A(:,:,L)
    enddo

! Sweep up, solving for updated value. Note B has the inverse of the main diagonal

    if(size(YG)>0) then
       Y(:,:,LM)   = (Y(:,:,LM) - C(:,:,LM) * YG        )*B(:,:,LM)
    else
       Y(:,:,LM)   =  Y(:,:,LM)*B(:,:,LM-1)/(B(:,:,LM-1) - A(:,:,LM)*(1.0+C(:,:,LM-1)*B(:,:,LM-1) ))
    endif

    do L = LM-1,1,-1
       Y(:,:,L) = (Y(:,:,L ) - C(:,:,L ) * Y(:,:,L+1))*B(:,:,L )
    enddo

    return
  end subroutine VTRISOLVE

!*********************************************************************

#ifdef _CUDA
   attributes(device) &
#endif
   function gpu_maxloc(array,len)

#ifdef _CUDA
      integer, value      :: len
#else
      integer, intent(in) :: len
#endif
      real,    intent(in) :: array(len)

      real    :: gpu_maxval
      integer :: gpu_maxloc

      integer :: i

      gpu_maxloc = 1
      gpu_maxval = array(1)

      do i = 1, len
         if (array(i) > gpu_maxval) then
            gpu_maxloc = i
            gpu_maxval = array(i)
         end if
      end do

   end function gpu_maxloc

end module GEOS_TurbulenceGridCompMod

