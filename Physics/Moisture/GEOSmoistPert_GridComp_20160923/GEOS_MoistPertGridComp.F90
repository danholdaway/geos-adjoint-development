! $Id: GEOS_MoistPertGridComp.F90,v 1.3.2.1.2.3.6.1 2014-12-18 16:06:51 drholdaw Exp $

! VERIFY_ and RETURN_ macros for error handling.

#include "MAPL_Generic.h"

#define MAPL_FieldBundleGetPointer ESMFL_BundleGetPointerToData

!=============================================================================
!BOP

! !MODULE: GEOS_MoistPert -- A Module to compute linearised moist processes, 
! including convection, large-scale condensation and precipitation and cloud.
! Although the model works with real*4 variables it is found that the moist 
! physics schemes do not satisfy the dot product test of equivalence between 
! the tangent and adjoint in real*4. Precision issues lead to differences 
! between the tlm and adjoint and the tlm at real*4 and the tlm at real*8. 
! By working in real*8 within this module avoids these issues from occuring
! and is found to have little impact on efficiency.

! !INTERFACE:

module GEOS_MoistPertGridCompMod

! GEOS Modules
  use ESMF
  use MAPL_MOD
  use GEOS_Mod
  use GEOS_PertSharedMod

! Saturation table
  use qsat_util

! Scheme constants
  use RASPARAMS
  use CLDPARAMS

! Moist Schemes
  use ras    , only: rase, rase_nowinds
  use ras_tlm, only: rase_tlm, rase_fast_tlm
  use ras_adm, only: rase_fwd, rase_bwd

  use cloudnew

!  use CONVECTION
!  use CONVECTION_AD
!  use CONVECTION_TL
!  use CLOUD
!  use CLOUD_AD
!  use CLOUD_TL

! Pert timers
  use mpp_mod, only: mpp_pe
  use fv_timing_mod,  only: timing_on, timing_off, timing_prt

  implicit none
  private

  logical, save :: comp_qsat_table = .true.

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION:
! 
!   {\tt GEOS_MoistPertGridCompMod} implements linearised moist processes in GEOS-5. 
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
!EOP

!=============================================================================
!
! ErrLog Variables

    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: IAm
    character(len=ESMF_MAXSTR)              :: COMP_NAME
    
!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,   name="RUNTLM"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="RUNADJ"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="-MOIST"       ,RC=STATUS)
    VERIFY_(STATUS)

! The same run method is registered twice
!  It queries the phase to do TLM or ADJ
! ---------------------------------------

    call MAPL_GridCompSetEntryPoint (gc, ESMF_METHOD_RUN,  Run,        rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint (gc, ESMF_METHOD_RUN,  Run,        rc=status)
    VERIFY_(STATUS)


! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)
     
    RETURN_(ESMF_SUCCESS)
     
  end subroutine SetServices


! Main run component
! ------------------

subroutine Run(gc, import, export, clock, rc)

  !Run arguments
  type(ESMF_GridComp), intent(inout) :: gc
  type (ESMF_State),   intent(inout) :: import
  type (ESMF_State),   intent(inout) :: export
  type (ESMF_Clock),   intent(inout) :: clock
  integer,  optional,  intent(  out) :: rc 

  !Mapl variables
  type (MAPL_MetaComp), pointer      :: MAPL
  integer                            :: STATUS
  character(len=ESMF_MAXSTR)         :: IAm
  character(len=ESMF_MAXSTR)         :: COMP_NAME
  character(len=ESMF_MAXSTR)         :: PHASE_NAME

  !Grid lookup for params
  character(len=ESMF_MAXSTR)         :: GRIDNAME
  character(len=4)                   :: imchar
  character(len=2)                   :: dateline
  integer                            :: imsize, nn
  
  !Control variables
  real(8)                            :: DT
  integer                            :: IM, JM, LM, IDIM, IRUN
  integer                            :: DO_MOIST_PHYS
  integer                            :: PHASE
  integer                            :: I, J, L
  type (ESMF_FieldBundle)            :: Traj3, Traj2, Pert

  !Reference and perturbation trajectories (io single precision)
  real(4), pointer, dimension(:,:,:) :: UT4, VT4, PTT4, QVT4            , QLST4, QCNT4, CFCNT4
  real(4), pointer, dimension(:,:,:) :: UP4, VP4, PTP4, QVP4, QIP4, QLP4              , CFCNP4
  real(4), pointer, dimension(:,:)   :: PS4, KCBL4, TS4, FRLAND4, KHu4, KHl4

  !Reference and perturbation trajectories (working double precision)
  real(8), allocatable, dimension(:,:,:) :: UT, VT, PTT, QVT
  real(8), allocatable, dimension(:,:,:) :: UP, VP, PTP, QVP            
  real(8), allocatable, dimension(:,:)   :: PS, TS, FRLAND
  integer, allocatable, dimension(:,:)   :: KCBL, KHu, KHl

  !Pressure and temperature variables
  real(8), pointer                       :: ak(:), bk(:)
  real(8), allocatable, dimension(:)     :: pref
  real(8), allocatable, dimension(:,:,:) :: PLE, CNV_PLE, PLO, PK, PKE
  real(8), allocatable, dimension(:,:,:) :: TEMP

  !Convection scheme variables
  integer :: ICMIN, MAXCONDEP, ITRCR
  real(8), parameter :: PMIN_DET = 3000.0, AUTOC_CN_OCN  = 2.5e-3, AUTOC_CN_LAND = AUTOC_CN_OCN
  real(8), allocatable, dimension(:,:,:) :: DQST, QSST, DQSP, QSSP
  real(8), allocatable, dimension(:,:,:) :: CNV_DQLDTT, CNV_MFDT, CNV_PRC3T, CNV_UPDFT
  real(8), allocatable, dimension(:,:,:) :: CNV_DQLDTP, CNV_MFDP, CNV_PRC3P, CNV_UPDFP
  real(8), allocatable, dimension(:,:,:) :: PTT_C, QVT_C
  real(8), allocatable, dimension(:,:,:) :: CNV_DQLDTT_C, CNV_MFDT_C, CNV_PRC3T_C, CNV_UPDFT_C
  real(8), allocatable, dimension(:,:,:) :: PTT_F, QVT_F
  real(8), allocatable, dimension(:,:,:) :: PTT_L, QVT_L
  integer, allocatable, dimension(:,:,:) :: SEEDRAS
  real(8), allocatable, dimension(:,:)   :: CO_AUTO, RASAL2_2d
  real(8), allocatable, dimension(:,:)   :: TPERTT, QPERTT, TPERTP, QPERTP
  real(8), allocatable, dimension(:)     :: SIGE
  real(8), allocatable, dimension(:,:,:) :: WGT0, WGT1
  integer, parameter :: n_rasparams = 25
  type(RASPARAM_TYPE)                    :: RASPARAMS
  real(8)                                :: RASAL1, RASAL2, CNV_Q600_MIN, CNV_Q600_MAX, tmp
  real(8), allocatable, dimension(:,:)   :: CNV_FRACTION, QV600
  integer                                :: levs600
  real(8)                                :: CBL_TPERT,CBL_QPERT,CBL_TPERT_MXOCN,CBL_TPERT_MXLND

  !Convection filtering
  real(8), allocatable, dimension(:,:,:) :: HEAT
  integer, allocatable, dimension(:,:)   :: DOCONVEC, CTOP
  real(8), allocatable, dimension(:,:)   :: sumHEAT
  real(8), allocatable, dimension(:,:) :: JACOBIAN 
  real(8), allocatable, dimension(:)   :: PT_pert, QV_pert
  real(8), allocatable, dimension(:)   :: PT_pert_in, QV_pert_in
  real(8), allocatable, dimension(:)   :: H_pert, M_pert
  real(8), allocatable, dimension(:)   :: TPERTTj, TPERTPj, QPERTTj, DQSTj, DQSPj, QSSTj, QSSPj

  !Cloud scheme variables
  type(CLDPARAM_TYPE)                    :: CLDPARAMS
  type (T_CLOUD_CTL)              :: CLOUD_CTL
  real(8), allocatable, dimension(:,:,:) :: QILST, QLLST, QICNT, QLCNT
  real(8), allocatable, dimension(:,:,:) :: QILSP, QLLSP, QICNP, QLCNP
  real(8), allocatable, dimension(:,:,:) :: CFLST, CFCNT
  real(8), allocatable, dimension(:,:,:) :: CFLSP, CFCNP
  real(8), allocatable, dimension(:,:,:) :: fQi
  real(8), allocatable, dimension(:,:,:) :: ILSF, ICNF, LLSF, LCNF

  !R8 versions of the required MAPL_Constants
  real(8) :: MAPL8_CP, MAPL8_ALHL, MAPL8_GRAV, MAPL8_P00, MAPL8_KAPPA
  real(8) :: MAPL8_RGAS, MAPL8_H2OMW, MAPL8_AIRMW, MAPL8_VIREPS
  real(8) :: MAPL8_RUNIV, MAPL8_ALHF, MAPL8_PI, MAPL8_ALHS
  real(8) :: MAPL8_TICE, MAPL8_RVAP
  real(8) :: MAPL8_UNDEF

  !One/two moment microphysics switch
  character(len=ESMF_MAXSTR) :: CLDMICRO

! Begin
!------
    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, currentPhase=Phase, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the generic state
! -----------------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

! Turn on MAPL timers
! -------------------
    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"-MOIST")
    call MAPL_TimerOn(MAPL,PHASE_NAME)

! Get time step of the linear model
! ---------------------------------
    call MAPL_GetResource(MAPL, DT, Label="RUN_DT:", RC=STATUS)
    VERIFY_(STATUS)

! Get moist physics switch
! ------------------------
    call MAPL_GetResource(MAPL, DO_MOIST_PHYS, Label="DO_MOIST_PHYSICS:", default=0,  RC=STATUS)
    VERIFY_(STATUS)

! Set phase
! ---------
    select case(phase)
    case(TLMphase)
       PHASE_NAME = "RUNTLM"
    case(ADJphase)
       PHASE_NAME = "RUNADJ"
    case default
       ASSERT_(.false.)
    end select

! Begin linearised moist physics calculations if asked for in AGCM_apert
! ----------------------------------------------------------------------
    IF (DO_MOIST_PHYS == 1 .or. DO_MOIST_PHYS == 2) THEN

! One or two-moment microphysics
! ------------------------------
       CLDMICRO = "1MOMENT"

! Get the phase
! -------------
       if (PHASE==ADJphase) then         
          call timing_on('MOIST_ADM')             
       elseif (PHASE==TLMphase) then         
          call timing_on('MOIST_TLM')            
       endif
 
! Read in reference and perturbation trajectories
! -----------------------------------------------
       call ESMF_StateGet(Import, "TRAJWRK3", Traj3, rc=status)
       VERIFY_(STATUS)
       call ESMF_StateGet(Import, 'TRAJWRK2', Traj2, rc=STATUS)
       VERIFY_(STATUS)
       call ESMF_StateGet(Import, "PERTWRK" , Pert , rc=STATUS)
       VERIFY_(STATUS)

       !U wind speed
       call MAPL_FieldBundleGetPointer(Traj3, "U"     , UT4    , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert , "U"     , UP4    , rc=STATUS)
       VERIFY_(STATUS)

       !V wind speed
       call MAPL_FieldBundleGetPointer(Traj3, "V"     , VT4    , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert , "V"     , VP4    , rc=STATUS)
       VERIFY_(STATUS)
       
       !Potential temperature
       call MAPL_FieldBundleGetPointer(Traj3, "PT"    , PTT4   , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert , "PT"    , PTP4   , rc=STATUS)
       VERIFY_(STATUS)
       
       !Specific humidity
       call MAPL_FieldBundleGetPointer(Traj3, "QV"    , QVT4   , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert , "QV"    , QVP4   , rc=STATUS)
       VERIFY_(STATUS)
  
       !Total cloud liquid ice and water (pert only)
       call MAPL_FieldBundleGetPointer(Pert , "QI"    , QIP4   , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert , "QL"    , QLP4   , rc=STATUS)
       VERIFY_(STATUS)

       !Cloud totals for large scale and convective (traj only)
       call MAPL_FieldBundleGetPointer(Traj3, "QLS"   , QLST4  , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Traj3, "QCN"   , QCNT4  , rc=STATUS)
       VERIFY_(STATUS)

       !Anvil fraction
       call MAPL_FieldBundleGetPointer(Traj3, "CFCN"  , CFCNT4 , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert , "CFCN"  , CFCNP4 , rc=STATUS)
       VERIFY_(STATUS)

       !Rank 2 variables (traj only)
       call MAPL_FieldBundleGetPointer(Traj2, "PS"     , PS4    , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Traj2, "KCBL"   , KCBL4  , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Traj2, "TS"     , TS4    , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Traj2, "FRLAND" , FRLAND4, rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Traj2, "KHL"    , KHL4   , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Traj2, "KHU"    , KHU4   , rc=STATUS)
       VERIFY_(STATUS)

! Find grid dimensions for processor
! ----------------------------------
       IM = size(UP4,1)
       JM = size(UP4,2)
       LM = size(UP4,3)

! Allocate memory for all the local variables
! -------------------------------------------
       !Double precision trajectories
       allocate(UT(im,jm,lm),           stat=status); VERIFY_(STATUS)
       allocate(VT(im,jm,lm),           stat=status); VERIFY_(STATUS)
       allocate(PTT(im,jm,lm),          stat=status); VERIFY_(STATUS)
       allocate(QVT(im,jm,lm),          stat=status); VERIFY_(STATUS)
       allocate(UP(im,jm,lm),           stat=status); VERIFY_(STATUS)
       allocate(VP(im,jm,lm),           stat=status); VERIFY_(STATUS)
       allocate(PTP(im,jm,lm),          stat=status); VERIFY_(STATUS)
       allocate(QVP(im,jm,lm),          stat=status); VERIFY_(STATUS)
       allocate(PS(im,jm),              stat=status); VERIFY_(STATUS)
       allocate(TS(im,jm),              stat=status); VERIFY_(STATUS)
       allocate(FRLAND(im,jm),          stat=status); VERIFY_(STATUS)
       allocate(KCBL(im,jm),            stat=status); VERIFY_(STATUS)
       allocate(KHl(IM,JM),             stat=status); VERIFY_(STATUS)
       allocate(KHu(IM,JM),             stat=status); VERIFY_(STATUS)

       !Pressure and temperature variables
       allocate(ak(size(pert_ak)),      stat=status); VERIFY_(STATUS)
       allocate(bk(size(pert_ak)),      stat=status); VERIFY_(STATUS)
       allocate(pref(0:size(AK)-1),     stat=status); VERIFY_(STATUS)
       allocate(PLE(im,jm,0:lm),        stat=status); VERIFY_(STATUS)
       allocate(CNV_PLE(im,jm,0:lm),    stat=status); VERIFY_(STATUS)
       allocate(PKE(im,jm,0:lm),        stat=status); VERIFY_(STATUS)
       allocate(PLO(im,jm,lm),          stat=status); VERIFY_(STATUS)
       allocate(PK(im,jm,lm),           stat=status); VERIFY_(STATUS)
       allocate(TEMP(IM,JM,LM),         stat=status); VERIFY_(STATUS)

       !Convection varaibles
       allocate(PTT_F(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QVT_F(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(PTT_L(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QVT_L(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(PTT_C(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QVT_C(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(DQST(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(QSST(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(DQSP(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(QSSP(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(CNV_DQLDTT_C(im,jm,lm), stat=status); VERIFY_(STATUS)
       allocate(CNV_MFDT_C(im,jm,lm),   stat=status); VERIFY_(STATUS)
       allocate(CNV_PRC3T_C(im,jm,lm),  stat=status); VERIFY_(STATUS)
       allocate(CNV_UPDFT_C(im,jm,lm),  stat=status); VERIFY_(STATUS)
       allocate(CNV_DQLDTT(im,jm,lm),   stat=status); VERIFY_(STATUS)
       allocate(CNV_DQLDTP(im,jm,lm),   stat=status); VERIFY_(STATUS)
       allocate(CNV_MFDT(im,jm,lm),     stat=status); VERIFY_(STATUS)
       allocate(CNV_MFDP(im,jm,lm),     stat=status); VERIFY_(STATUS)
       allocate(CNV_PRC3T(im,jm,lm),    stat=status); VERIFY_(STATUS)
       allocate(CNV_PRC3P(im,jm,lm),    stat=status); VERIFY_(STATUS)
       allocate(CNV_UPDFT(im,jm,lm),    stat=status); VERIFY_(STATUS)
       allocate(CNV_UPDFP(im,jm,lm),    stat=status); VERIFY_(STATUS)
       allocate(SIGE(0:lm),             stat=status); VERIFY_(STATUS)
       allocate(WGT0(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(WGT1(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(CO_AUTO(im,jm),         stat=status); VERIFY_(STATUS)
       allocate(RASAL2_2d(im,jm),       stat=status); VERIFY_(STATUS)
       allocate(TPERTT(im,jm),          stat=status); VERIFY_(STATUS)
       allocate(QPERTT(im,jm),          stat=status); VERIFY_(STATUS)
       allocate(TPERTP(im,jm),          stat=status); VERIFY_(STATUS)
       allocate(QPERTP(im,jm),          stat=status); VERIFY_(STATUS)
       allocate(SEEDRAS(im,jm,2),       stat=status); VERIFY_(STATUS)
       allocate(CNV_FRACTION(im,jm),    stat=status); VERIFY_(STATUS)
       allocate(QV600(im,jm),           stat=status); VERIFY_(STATUS)
       allocate(DOCONVEC(im,jm),        stat=status); VERIFY_(STATUS)
       allocate(HEAT(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(CTOP(im,jm),            stat=status); VERIFY_(STATUS)
       allocate(sumHEAT(im,jm),         stat=status); VERIFY_(STATUS)
       allocate(JACOBIAN(2*lm,2),       stat=status); VERIFY_(STATUS)
       allocate(H_pert(lm),             stat=status); VERIFY_(STATUS)
       allocate(M_pert(lm),             stat=status); VERIFY_(STATUS) 
       allocate(PT_pert(lm),            stat=status); VERIFY_(STATUS)
       allocate(QV_pert(lm),            stat=status); VERIFY_(STATUS)
       allocate(PT_pert_in(lm),         stat=status); VERIFY_(STATUS)
       allocate(QV_pert_in(lm),         stat=status); VERIFY_(STATUS)
       allocate(TPERTTj(lm),         stat=status); VERIFY_(STATUS)
       allocate(TPERTPj(lm),         stat=status); VERIFY_(STATUS)
       allocate(QPERTTj(lm),         stat=status); VERIFY_(STATUS)
       allocate(DQSTj(lm),         stat=status); VERIFY_(STATUS)
       allocate(DQSPj(lm),         stat=status); VERIFY_(STATUS)
       allocate(QSSTj(lm),         stat=status); VERIFY_(STATUS)
       allocate(QSSPj(lm),         stat=status); VERIFY_(STATUS)

       !Cloud variables
       allocate(QILST(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QLLST(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QICNT(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QLCNT(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QILSP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QLLSP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QICNP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QLCNP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(CFLST(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(CFCNT(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(CFLSP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(CFCNP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(fQi(IM,JM,LM),          stat=status); VERIFY_(STATUS)
       allocate(ILSF(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(ICNF(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(LLSF(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(LCNF(im,jm,lm),         stat=status); VERIFY_(STATUS)

! Create double precision versions of required constants and varibles
! -------------------------------------------------------------------
       !MAPL_Constants
       MAPL8_CP     = dble(MAPL_CP)
       MAPL8_ALHL   = dble(MAPL_ALHL)
       MAPL8_GRAV   = dble(MAPL_GRAV)
       MAPL8_P00    = dble(MAPL_P00)
       MAPL8_KAPPA  = dble(MAPL_KAPPA)
       MAPL8_RGAS   = dble(MAPL_RGAS)
       MAPL8_H2OMW  = dble(MAPL_H2OMW)
       MAPL8_AIRMW  = dble(MAPL_AIRMW)
       MAPL8_VIREPS = dble(MAPL_VIREPS)
       MAPL8_RUNIV  = dble(MAPL_RUNIV)
       MAPL8_ALHF   = dble(MAPL_ALHF)
       MAPL8_PI     = dble(MAPL_PI)
       MAPL8_ALHS   = dble(MAPL_ALHS)
       MAPL8_TICE   = dble(MAPL_TICE)
       MAPL8_RVAP   = dble(MAPL_RVAP)
       MAPL8_UNDEF  = dble(MAPL_UNDEF)

       !Copy reference trajectory to double precision
       UT     = dble(UT4)
       VT     = dble(VT4)
       PTT    = dble(PTT4)*(MAPL8_P00**MAPL8_KAPPA)
       QVT    = dble(QVT4)
       CFLST  = 0.0
       CFCNT  = dble(CFCNT4)
       PS     = dble(PS4)
       TS     = dble(TS4)
       FRLAND = dble(FRLAND4)
       KCBL   = nint(KCBL4)
       KHL    = nint(KHL4)
       KHU    = nint(KHU4)

       !Copy perturbation trajectory to double precision
       UP = dble(UP4)
       VP = dble(VP4)
       if (PHASE==ADJphase) then
          PTP = dble(PTP4)/(MAPL8_P00**MAPL8_KAPPA)
       elseif (PHASE==TLMphase) then 
          PTP = dble(PTP4)*(MAPL8_P00**MAPL8_KAPPA)
       endif
       QVP = dble(QVP4)
       CFLSP = 0.0
       CFCNP = dble(CFCNP4)

! Create additional versions of the trajectory to be overwritten
! --------------------------------------------------------------
       PTT_C = PTT
       QVT_C = QVT
       CNV_DQLDTT_C  = 0.0
       CNV_MFDT_C    = 0.0
       CNV_PRC3T_C   = 0.0
       CNV_UPDFT_C   = 0.0

       PTT_F = PTT
       QVT_F = QVT

       PTT_L = PTT
       QVT_L = QVT

! Define up the constants used by the convection scheme
! -----------------------------------------------------
       !RASPARAMS%MAXDALLOWED_S and CLDPARAMS%MINRHCRIT Set based on resolution.
       !The following reads PHYSICS_PERT_GRIDNAME: from AGCM_APERT. E.G. "PHYSICS_PERT_GRIDNAME: PE90x540-CF"
       !90X540 means face of cube is 90 BY 90 (1 degree), 540/90 = 6 sides of the cube)
       !Equivalnet format for lat-long is "GRIDNAME: PC288x181-DC"
       call MAPL_GetResource(MAPL,GRIDNAME,'PHYSICS_PERT_GRIDNAME:', RC=STATUS)
       VERIFY_(STATUS)
       GRIDNAME =  AdjustL(GRIDNAME)
             nn = len_trim(GRIDNAME)
       dateline = GRIDNAME(nn-1:nn)
       imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
       read(imchar,*) imsize
       if(dateline.eq.'CF') imsize = imsize*4

       RASPARAMS%CUFRICFAC      = 1.0
       RASPARAMS%SHR_LAMBDA_FAC = 0.05
       RASPARAMS%QC_CRIT_CN     = 8.0e-4
       RASPARAMS%RASAL1         = 1800.0
      if( imsize.le.200                      ) RASPARAMS%RASAL2 = -43200.0
      if( imsize.gt.200 .and. imsize.le.400  ) RASPARAMS%RASAL2 = -43200.0
      if( imsize.gt.400 .and. imsize.le.800  ) RASPARAMS%RASAL2 = -43200.0
      if( imsize.gt.800 .and. imsize.le.1600 ) RASPARAMS%RASAL2 = -64800.0
      if( imsize.gt.1600                     ) RASPARAMS%RASAL2 = -86400.0
       RASPARAMS%RASNCL         = -300.
       RASPARAMS%LAMBDA_FAC     = 4.0 
       RASPARAMS%LAMBMX_FAC     = 0.0 
       RASPARAMS%MIN_DIAMETER   = 200.
       RASPARAMS%CUFRICLAMBDA   = 7.5e-4
       RASPARAMS%RDTLEXPON      = 1.0 
       RASPARAMS%STRAPPING      = -1.0 
       RASPARAMS%SDQV2          = 1.3 
       RASPARAMS%SDQV3          = 1.3 
       RASPARAMS%SDQVT1         = 263.
       RASPARAMS%ACRITFAC       = 0.5 
       RASPARAMS%HMINTRIGGER    = 1.0 
       RASPARAMS%LLDISAGGXP     = 0.0 
       RASPARAMS%PBLFRAC        = 0.1 
       RASPARAMS%RASAUTORAMPB   = 0.8 
       RASPARAMS%AUTOC_CN_ZDEP  = 1.0
       if( imsize.le.200                      ) RASPARAMS%MAXDALLOWED_S = 4000.0
       if( imsize.gt.200 .and. imsize.le.400  ) RASPARAMS%MAXDALLOWED_S = 2000.0
       if( imsize.gt.400 .and. imsize.le.800  ) RASPARAMS%MAXDALLOWED_S =  700.0
       if( imsize.gt.800 .and. imsize.le.1600 ) RASPARAMS%MAXDALLOWED_S =  450.0
       if( imsize.gt.1600                     ) RASPARAMS%MAXDALLOWED_S =  300.0
       RASPARAMS%MAXDALLOWED_D  = RASPARAMS%MAXDALLOWED_S
       RASPARAMS%RASAL_EXP      = 1
       RASPARAMS%RAS_RHMIN      = 0.5
       RASPARAMS%RAS_RHFULL     = 0.65
       RASPARAMS%CLDMICRO       = 0.0

       CBL_QPERT  = 0.0
       if( imsize.le.800                      ) CBL_TPERT = -1.0
       if( imsize.gt.800 .and. imsize.le.1600 ) CBL_TPERT = -0.5
       if( imsize.gt.1600                     ) CBL_TPERT =  0.0
       CBL_TPERT_MXOCN = 2.0
       CBL_TPERT_MXLND = 0.0

! Define constants that are used by the cloud scheme
! --------------------------------------------------
       !SET SBAC PARAMETERS
       CLDPARAMS%CNV_BETA        = 10.0
       CLDPARAMS%ANV_BETA        = 4.0
       CLDPARAMS%LS_BETA         = 4.0
       CLDPARAMS%RH_CRIT         = 1.0
       CLDPARAMS%AUTOC_LS        = 2.0e-3
       CLDPARAMS%QC_CRIT_LS      = 8.0e-4
       CLDPARAMS%ACCRETION       = 2.0
       if( imsize.le.200                      ) CLDPARAMS%RAIN_REVAP_FAC = 1.0
       if( imsize.gt.200 .and. imsize.le.400  ) CLDPARAMS%RAIN_REVAP_FAC = 1.0
       if( imsize.gt.400 .and. imsize.le.800  ) CLDPARAMS%RAIN_REVAP_FAC = 1.0
       if( imsize.gt.800 .and. imsize.le.1600 ) CLDPARAMS%RAIN_REVAP_FAC = 2.0
       if( imsize.gt.1600                     ) CLDPARAMS%RAIN_REVAP_FAC = 3.0
       CLDPARAMS%VOL_TO_FRAC     = -1.0
       CLDPARAMS%SUPERSAT        = 0.0
       CLDPARAMS%SHEAR_EVAP_FAC  = 1.3
       CLDPARAMS%MIN_ALLOW_CCW   = 1.0e-9
       CLDPARAMS%CCW_EVAP_EFF    = 3.3e-4
       CLDPARAMS%NSUB_AUTOCONV   = 20.
       CLDPARAMS%LS_SUND_INTER   = 4.8
       CLDPARAMS%LS_SUND_COLD    = 4.8
       CLDPARAMS%LS_SUND_TEMP1   = 230.
       CLDPARAMS%ANV_SUND_INTER  = 1.0
       CLDPARAMS%ANV_SUND_COLD   = 1.0
       CLDPARAMS%ANV_SUND_TEMP1  = 230.
       CLDPARAMS%ANV_TO_LS_TIME  = 14400.
       CLDPARAMS%NCCN_WARM       = 50.
       CLDPARAMS%NCCN_ICE        = 0.01
       CLDPARAMS%NCCN_ANVIL      = 0.1
       CLDPARAMS%NCCN_PBL        = 200.
       CLDPARAMS%DISABLE_RAD     = 0.
       CLDPARAMS%ICE_SETTLE      = 1.0
       CLDPARAMS%ANV_ICEFALL     = 0.5
       CLDPARAMS%LS_ICEFALL      = 1.0
       CLDPARAMS%REVAP_OFF_P     = 2000.
!      CLDPARAMS%CNV_ENVF        = 0.8
       CLDPARAMS%WRHODEP         = 0.5
       CLDPARAMS%ICE_RAMP        = -40.0
       CLDPARAMS%CNV_ICEPARAM    = 1.0
       CLDPARAMS%CNV_ICEFRPWR    = 4.0
       CLDPARAMS%CNV_DDRF        = 0.0
       CLDPARAMS%ANV_DDRF        = 0.0
       CLDPARAMS%LS_DDRF         = 0.0
       CLDPARAMS%AUTOC_ANV       = 1.0e-3
       CLDPARAMS%QC_CRIT_ANV     = 8.0e-4
       CLDPARAMS%TANHRHCRIT      = 1.0
       if( imsize.le.200                      ) CLDPARAMS%MINRHCRIT = 0.80
       if( imsize.gt.200 .and. imsize.le.400  ) CLDPARAMS%MINRHCRIT = 0.90
       if( imsize.gt.400 .and. imsize.le.800  ) CLDPARAMS%MINRHCRIT = 0.93
       if( imsize.gt.800 .and. imsize.le.1600 ) CLDPARAMS%MINRHCRIT = 0.95
       if( imsize.gt.1600                     ) CLDPARAMS%MINRHCRIT = 0.97
       CLDPARAMS%MAXRHCRIT = 1.0
       if(adjustl(CLDMICRO) =="2MOMENT") then
          CLDPARAMS%PRECIPRAD       = 1.0
          CLDPARAMS%MAXRHCRITLAND   = 1.0
          CLDPARAMS%SNOW_REVAP_FAC  = 0.5
	  CLDPARAMS%CNV_ENVF        = 1.0
   	  CLOUD_CTL%SCLMFDFR        = 1.0
  	  CLDPARAMS%TURNRHCRIT      = 884.0
          CLDPARAMS%MOVE2RAS        = 0.0
       else   
          CLDPARAMS%PRECIPRAD       = 0.0
          CLDPARAMS%MAXRHCRITLAND   = CLDPARAMS%MINRHCRIT+0.01
          CLDPARAMS%SNOW_REVAP_FAC  = 1.0
	  CLDPARAMS%CNV_ENVF        = 0.8
	  CLOUD_CTL%SCLMFDFR        = 1.00
	  CLDPARAMS%TURNRHCRIT      = 750.0
          CLDPARAMS%MOVE2RAS        = 1.0
       end if
       CLDPARAMS%FR_LS_WAT       = 1.0
       CLDPARAMS%FR_AN_WAT       = 0.0
       CLDPARAMS%FR_LS_ICE       = 0.0
       CLDPARAMS%FR_AN_ICE       = 0.0
       CLDPARAMS%MIN_RL          = 10.e-6
       CLDPARAMS%MIN_RI          = 20.e-6
       CLDPARAMS%MAX_RL          = 21.e-6
       CLDPARAMS%MAX_RI          = 40.e-6
       CLDPARAMS%RI_ANV          = 30.e-6
       CLDPARAMS%PDFSHAPE        = 1.0
       CLDPARAMS%TURNRHCRIT_UP   = 300.0

! Generate saturation vapour pressure loopup table
! ------------------------------------------------
       if (comp_qsat_table == .true.) then
          call ESINIT
          comp_qsat_table = .false.
       endif

! Compute reference pressure, pressure, Exner pressure and temperature
! -------------------------------------------------------------------
       !Compute pref from ak and bk
       ak = pert_ak
       bk = pert_bk
       DO L = 0,lm
          pref(L) = AK(L+1) + BK(L+1)*MAPL8_P00
       enddo

       !Pressure at the half levels from Ps
       DO L = 0,lm
          PLE(:,:,L) = AK(L+1) + BK(L+1)*PS(:,:)
       enddo

! Prepare convection inputs and call nonlinear convection scheme 
! --------------------------------------------------------------

       CNV_PLE  = PLE*.01
       PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
       PKE      = (CNV_PLE/1000.)**(MAPL8_RGAS/MAPL8_CP)
       PK       = (PLO/1000.)**(MAPL8_RGAS/MAPL8_CP)
       TEMP     = PTT*PK !Keep a higher level version for non-linearized quantities

       ICMIN = max(1,count(PREF < PMIN_DET))
       SIGE = PREF/PREF(LM)

       ITRCR = 0

       !Strapping levels
       DO I = 1,IM
          DO J = 1,JM
             WGT0(I,J,:)            = 0.0
             WGT0(I,J,KCBL(I,J):LM) = 1.0 
             WGT1(I,J,:)            = 0.0
             WGT1(I,J,KCBL(I,J):LM) = 1.0 
          ENDDO
       ENDDO

       where (FRLAND<0.1) 
          CO_AUTO = AUTOC_CN_OCN   ! ocean value
       elsewhere
          CO_AUTO = AUTOC_CN_LAND  ! land value
       end where

       IDIM = IM*JM
       IRUN = IM*JM

       tmp = log(60000.)
       call VertInterp(QV600,QVT,log(PLE),tmp,MAPL8_UNDEF,STATUS)
       VERIFY_(STATUS)
       !Fill undefs (600mb below the surface) with surface QV values L=LM
       levs600  = max(1,count(PREF < 60000.))
       WHERE (QV600 == MAPL8_UNDEF)
          QV600 = QVT(:,:,levs600)
       END WHERE

       if ( imsize.le.800  ) then
         CNV_Q600_MIN = 0.00250
         CNV_Q600_MAX = 0.00600
       endif
       if ( imsize.gt.800 .and. imsize.le.1600 ) then
         CNV_Q600_MIN = 0.00375
         CNV_Q600_MAX = 0.00675
       endif
       if ( imsize.gt.1600 ) then
          CNV_Q600_MIN = 0.00500
          CNV_Q600_MAX = 0.00750
       endif

       CNV_FRACTION = 0.0 !Dont linearize for now
       if ( CNV_Q600_MAX > CNV_Q600_MIN ) then
          DO J = 1,JM
             DO I = 1,IM
                CNV_FRACTION(I,J) = MAX(0.0,MIN(1.0,(QV600(I,J)-CNV_Q600_MIN)/(CNV_Q600_MAX-CNV_Q600_MIN)))
             END DO
          END DO
       endif

       RASAL1 = RASPARAMS%RASAL1
       RASAL2 = RASPARAMS%RASAL2

       if (RASAL2 > 0.0) then
          RASAL2_2d(:,:) = RASAL2
       else
          !if (0) then 
          !   ! include KH dependence
          !   DO J=1, JM
          !      DO I=1, IM
          !         RASAL2_2d(I,J) = ( RASAL1 * ( QUANX(I,J) - QUANX_0 ) + RASAL2 * ( QUANX_1 - QUANX(I,J) ) ) &
          !                          /       ( QUANX_1    - QUANX_0 )
          !      END DO
          !   END DO
          !else
             ! include CNV dependence
             DO J=1, JM
                DO I=1, IM
                   RASAL2_2d(I,J) = CNV_FRACTION(I,J)*ABS(RASAL2) + (1-CNV_FRACTION(I,J))*RASAL1
                END DO
             END DO
          !endif
       endif

       SEEDRAS(:,:,1) = 1000000 * ( 100*TEMP(:,:,LM)   - INT( 100*TEMP(:,:,LM) ) )
       SEEDRAS(:,:,2) = 1000000 * ( 100*TEMP(:,:,LM-1) - INT( 100*TEMP(:,:,LM-1) ) )

       !Call nonlinear convection scheme. We need to do this because the
       !cloud adjoint scheme needs the outputs from the convection.
       !We will also only call the convective linearizations for profiles
       !where convection is occuring, which is a big time saver.
       CALL PRE_RASE( IM,JM,LM,PLE,PTT_C,QVT_C,TS,CNV_FRACTION,FRLAND, &
                      PLO,PKE,PK, &
                      CBL_TPERT,CBL_QPERT,CBL_TPERT_MXOCN,CBL_TPERT_MXLND, &
                      TPERTT,QPERTT,DQST,QSST )
 
       call RASE_NOWINDS( IDIM                 , &
                          IRUN                 , &                 
                          LM                   , &
                          ICMIN                , &
                          DT                   , &
                          MAPL8_CP             , &
                          MAPL8_ALHL           , &
                          MAPL8_ALHS           , &
                          MAPL8_TICE           , &
                          MAPL8_GRAV           , &
                          SEEDRAS              , &
                          SIGE                 , &
                          KCBL                 , &
                          WGT0                 , &
                          WGT1                 , &
                          TPERTT               , &
                          QPERTT               , &
                          PTT_C                , &
                          QVT_C                , & 
                          QSST                 , & 
                          DQST                 , &
                          CNV_FRACTION         , &
                          RASAL2_2d            , &
                          CO_AUTO              , &
                          CNV_PLE              , &
                          PKE                  , &
                          CNV_DQLDTT_C         , &
                          CNV_MFDT_C           , &
                          CNV_PRC3T_C          , &
                          CNV_UPDFT_C          , &
                          RASPARAMS            , &
                          ITRCR                  )

! Do the filtering to determine whether linear convection should be called
! ------------------------------------------------------------------------
       !Figure out whether or not each convective profile should be linearized.
       !Conditions  - Convection is happening (heating rate is nonzero)
       !            - Convection is deep enough (>= 10 levels)
       !            - The heating rate profile is not just a spike at one level
       !            - Gradients are not too steep (Jacobian filtering).
       DOCONVEC = 0
       HEAT = 0.0
       CTOP = LM
       sumHEAT = 0.0

       !Be more lenient on profiles let in for 4DVAR, less likely to encounter problems
       !at shorter lead times and gives more realistic low level cloud perturbations.
       if (DO_MOIST_PHYS == 1) then
          MAXCONDEP = 1
       elseif (DO_MOIST_PHYS == 2) then
          MAXCONDEP = 10
       endif

       DO I = 1,IM
          DO J = 1,JM

             !Compute the heating rate.
             HEAT(I,J,:) = (PTT_C(I,J,:) - PTT(I,J,:))/DT

             !Starting at the top scan downwards to look for nonzero heating rate,
             !record index of highest level convection reaches (ignoring v.small heating rate)
             DO L = 1,LM
                IF ( abs(HEAT(I,J,L)) .gt. 0.01*maxval(abs(HEAT(I,J,:)),1) ) THEN
                   CTOP(I,J) = L
                   exit
                ENDIF
             ENDDO
             
             !Compute sort of integral of the heating rate.
             if ( (CTOP(I,J) .ne. LM) .and. (KCBL(I,J) - CTOP(I,J) > 0) ) then
                sumHEAT(I,J) = (   sum(abs(HEAT(I,J,CTOP(I,J):KCBL(I,J)-1))) &
                                       - maxval(abs(HEAT(I,J,CTOP(I,J):KCBL(I,J)-1)),1) ) &
                                       / ( KCBL(I,J) - CTOP(I,J)  )
             endif

             !Compare `integral` to maximum absolute heating rate
             IF ( KCBL(I,J) - CTOP(I,J) >= MAXCONDEP ) THEN
                IF (sumHEAT(I,J) / maxval(abs(HEAT(I,J,1:KCBL(I,J)-1)),1) > 0.125 ) then
                   DOCONVEC(I,J) = 1
                endif
             endif

             !Compute two columns of the Jacobian to check for steep gradients
             !This prevents instability and floating point issues that can cause failure of the TLM/ADJ dot prod test
             IF ( DOCONVEC(I,J) == 1 ) THEN
                call jacobian_filter_tlm
             ENDIF

          enddo
       enddo

! Prepare the inputs for the cloud scheme 
! ---------------------------------------
       !Compute the ice fraction for each grid box
       DO i = 1,IM
          DO j = 1,JM
             DO l = 1,LM
                call IceFraction8( TEMP(i,j,l), fQi(i,j,l) )
             enddo
          enddo
       enddo

       !Split the input large scale and convective cloud into ice and liquid parts.
       QILST = QLST4 * fQi
       QLLST = QLST4 * (1-fQi)
       QICNT = QCNT4 * fQi
       QLCNT = QCNT4 * (1-fQi)

       !Split the perturbations for total cloud water and ice into the convective and large scale parts.
       !Spilitting is based on fraction of total cloud ice/water that is attributed to large scale and anvil
       !in the trajectory at this time step. Note that we don't really need to protect for small values of
       !total cloud sum, if the denominator is small then so is the numerator, though compiler may complain.
       ILSF = 0.0
       ICNF = 0.0
       LLSF = 0.0
       LCNF = 0.0
       DO I = 1,IM
          DO J = 1,JM
             DO L = 1,LM
 
                if ( QILST(i,j,l) + QICNT(i,j,l) .gt. 0.0 ) then
                   ILSF(i,j,l) = QILST(i,j,l) / ( QILST(i,j,l) + QICNT(i,j,l) )
                   ICNF(i,j,l) = QICNT(i,j,l) / ( QILST(i,j,l) + QICNT(i,j,l) )
                endif
                if ( QLLST(i,j,l) + QLCNT(i,j,l) .gt. 0.0 ) then
                   LLSF(i,j,l) = QLLST(i,j,l) / ( QLLST(i,j,l) + QLCNT(i,j,l) )
                   LCNF(i,j,l) = QLCNT(i,j,l) / ( QLLST(i,j,l) + QLCNT(i,j,l) )
                endif

             enddo
          enddo
       enddo

       !Use fraction to split the perturbation total liquid ice and water.
       if (PHASE==ADJphase) then
          !Adjoint of summing to QIT/QLT
          QILSP = QIP4
          QICNP = QIP4
          QLLSP = QLP4
          QLCNP = QLP4
       elseif (PHASE==TLMphase) then
          !Splitting to large scale and convective parts
          QILSP = QIP4 * ILSF
          QICNP = QIP4 * ICNF
          QLLSP = QLP4 * LLSF
          QLCNP = QLP4 * LCNF
       endif

! Reset/initialise the RAS outputs to zero
! -------------------------------------------
       CNV_DQLDTT  = 0.0
       CNV_DQLDTP  = 0.0
       CNV_MFDT    = 0.0
       CNV_MFDP    = 0.0
       CNV_PRC3T   = 0.0
       CNV_PRC3P   = 0.0
       CNV_UPDFT   = 0.0
       CNV_UPDFP   = 0.0

! Perform tangent linear and adjoint computations of the main convection and cloud schemes
! ----------------------------------------------------------------------------------------
       if (PHASE==ADJphase) then ! Do adjoint


          !Call the adjoint of the cloud scheme
          !CALL CLOUD_DRIVER_B ( DT, IM, JM, LM, PTT_C, PTP, QVT_C, QVP, PLE,        &
          !                      CNV_DQLDTT_C, CNV_DQLDTP, CNV_MFDT_C, CNV_MFDP,     &
          !                      CNV_PRC3T_C, CNV_PRC3P, CNV_UPDFT_C, CNV_UPDFP,     &
          !                      QILST, QILSP, QLLST, QLLSP,                         &
          !                      QICNT, QICNP, QLCNT, QLCNP,                         &
          !                      CFLST, CFLSP, CFCNT, CFCNP,                         &
          !                      FRLAND, CLDPARAMS, ESTBLX, KHu, KHl,              &
          !                      MAPL8_RUNIV, MAPL8_KAPPA, MAPL8_airmw, MAPL8_h2omw, &
          !                      MAPL8_GRAV, MAPL8_ALHL, MAPL8_ALHF, MAPL8_PI,       &
          !                      MAPL8_RGAS, MAPL8_CP, MAPL8_VIREPS, MAPL8_ALHS,     &
          !                      MAPL8_TICE, MAPL8_RVAP, MAPL8_P00, DO_MOIST_PHYS    )

          !Call the adjoint of the convection scheme
          DO I = 1,IM
             DO J = 1,JM
                if ( DOCONVEC(I,J) == 1) then

                   !call RASE_B (1, 1, LM, ICMIN, DT,                              &
                   !             MAPL8_CP, MAPL8_ALHL, MAPL8_GRAV, MAPL8_RGAS,     &
                   !             MAPL8_H2OMW, MAPL8_AIRMW, MAPL8_VIREPS,           &
                   !             SEEDRAS(I,J), SIGE,                               &
                   !             KCBL(I,J),                                        &
                   !             WGT0(I,J,:), WGT1(I,J,:),                         &
                   !             FRLAND(I,J), Ts(I,J),                             &
                   !             PTT_L(I,J,:), PTP(I,J,:),                         &
                   !             QVT_L(I,J,:), QVP(I,J,:),                         &
                   !             UT(I,J,:), UP(I,J,:),                             &
                   !             VT(I,J,:), VP(I,J,:),                             &
                   !             CO_AUTO(I,J), CNV_PLE(I,J,:),                     &
                   !             CNV_DQLDTT(I,J,:), CNV_DQLDTP(I,J,:),             &
                   !             CNV_MFDT(I,J,:),   CNV_MFDP(I,J,:),               &
                   !             CNV_PRC3T(I,J,:),  CNV_PRC3P(I,J,:),              &
                   !             CNV_UPDFT(I,J,:),  CNV_UPDFP(I,J,:),              &
                   !             RASPARAMS, ESTBLX                                 )

                endif
             enddo
          enddo
       

          !CALL PRE_RASE_FWD( IM, JM, LM, PLE, PTT, QVT, TS,                &
          !                   CNV_FRACTION, FRLAND, PLO, PKE, PK,                     &
          !                   CBL_TPERT, CBL_QPERT, CBL_TPERT_MXOCN, CBL_TPERT_MXLND, &
          !                   TPERTT, QPERTT, DQST, QSST         )
          !CALL RASE_FWD( IDIM, IRUN, LM, ICMIN, DT, DOCONVEC, MAPL8_CP,                 &
          !               MAPL8_ALHL, MAPL8_ALHS, MAPL8_TICE, MAPL8_GRAV, SEEDRAS, SIGE, &
          !               KCBL, WGT0, WGT1, TPERTT, QPERTT, PTT, QVT,       &
          !               UT, VT, QSST, DQST,                   &
          !               CNV_FRACTION, RASAL2_2d, CO_AUTO, CNV_PLE, PKE,                &
          !               CNV_DQLDTT,                                         &
          !               CNV_MFDT,                                             &
          !               CNV_PRC3T,                                           &
          !               CNV_UPDFT,                                           &
          !               RASPARAMS, ITRCR                                               )
          !CALL RASE_BWD( IDIM, IRUN, LM, ICMIN, DT, DOCONVEC, MAPL8_CP,                 &
          !               MAPL8_ALHL, MAPL8_ALHS, MAPL8_TICE, MAPL8_GRAV, SEEDRAS, SIGE, &
          !               KCBL, WGT0, WGT1, TPERTT, TPERTP, QPERTT, PTT, PTP, QVT,       &
          !               QVP, UT, UP, VT, VP, QSST, QSSP, DQST, DQSP,                   &
          !               CNV_FRACTION, RASAL2_2d, CO_AUTO, CNV_PLE, PKE,                &
          !               CNV_DQLDTT, CNV_DQLDTP,                                        &
          !               CNV_MFDT,   CNV_MFDP,                                          &
          !               CNV_PRC3T,  CNV_PRC3P,                                         &
          !               CNV_UPDFT,  CNV_UPDFP,                                         &
          !               RASPARAMS, ITRCR                                               )
          !CALL PRE_RASE_BWD( IM, JM, LM, PLE, PTT, PTP, QVT, QVP, TS,                &
          !                   CNV_FRACTION, FRLAND, PLO, PKE, PK,                     &
          !                   CBL_TPERT, CBL_QPERT, CBL_TPERT_MXOCN, CBL_TPERT_MXLND, &
          !                   TPERTT, TPERTP, QPERTT, DQST, DQSP, QSST, QSSP          )

       else

!          CALL PRE_RASE_TLM( IM, JM, LM, PLE, PTT, PTP, QVT, QVP, TS,                &
!                             CNV_FRACTION, FRLAND, PLO, PKE, PK,                     &
!                             CBL_TPERT, CBL_QPERT, CBL_TPERT_MXOCN, CBL_TPERT_MXLND, &
!                             TPERTT, TPERTP, QPERTT, DQST, DQSP, QSST, QSSP          )
!
!          CALL RASE_TLM( IDIM, IRUN, LM, ICMIN, DT, DOCONVEC, MAPL8_CP,                 &
!                         MAPL8_ALHL, MAPL8_ALHS, MAPL8_TICE, MAPL8_GRAV, SEEDRAS, SIGE, &
!                         KCBL, WGT0, WGT1, TPERTT, TPERTP, QPERTT, PTT, PTP, QVT,       &
!                         QVP, UT, UP, VT, VP, QSST, QSSP, DQST, DQSP,                   &
!                         CNV_FRACTION, RASAL2_2d, CO_AUTO, CNV_PLE, PKE,                &
!                         CNV_DQLDTT, CNV_DQLDTP,                                        &
!                         CNV_MFDT,   CNV_MFDP,                                          &
!                         CNV_PRC3T,  CNV_PRC3P,                                         &
!                         CNV_UPDFT,  CNV_UPDFP,                                         &
!                         RASPARAMS, ITRCR                                               )

          DO I = 1,IM
             DO J = 1,JM
                if ( DOCONVEC(I,J) == 1) then

          CALL PRE_RASE_TLM( 1, 1, LM, PLE(I,J,:), PTT(I,J,:), PTP(I,J,:), QVT(I,J,:), QVP(I,J,:), TS,                &
                             CNV_FRACTION(I,J), FRLAND(I,J), PLO(I,J,:), PKE(I,J,:), PK(I,J,:),                     &
                             CBL_TPERT, CBL_QPERT, CBL_TPERT_MXOCN, CBL_TPERT_MXLND, &
                             TPERTT(I,J), TPERTP(I,J), QPERTT(I,J), DQST(I,J,:), DQSP(I,J,:), QSST(I,J,:), QSSP(I,J,:)          )

          CALL RASE_TLM( 1, 1, LM, ICMIN, DT, MAPL8_CP,                 &
                         MAPL8_ALHL, MAPL8_ALHS, MAPL8_TICE, MAPL8_GRAV, SEEDRAS(I,J,:), SIGE, &
                         KCBL(I,J), WGT0(I,J,:), WGT1(I,J,:), TPERTT(I,J), TPERTP(I,J), QPERTT(I,J), PTT(I,J,:), PTP(I,J,:), QVT(I,J,:),       &
                         QVP(I,J,:), UT(I,J,:), UP(I,J,:), VT(I,J,:), VP(I,J,:), QSST(I,J,:), QSSP(I,J,:), DQST(I,J,:), DQSP(I,J,:),                   &
                         CNV_FRACTION(I,J), RASAL2_2d(I,J), CO_AUTO(I,J), CNV_PLE(I,J,:), PKE(I,J,:),                &
                         CNV_DQLDTT(I,J,:), CNV_DQLDTP(I,J,:),                                        &
                         CNV_MFDT(I,J,:),   CNV_MFDP(I,J,:),                                          &
                         CNV_PRC3T(I,J,:),  CNV_PRC3P(I,J,:),                                         &
                         CNV_UPDFT(I,J,:),  CNV_UPDFP(I,J,:),                                         &
                         RASPARAMS, ITRCR                                               )

                endif
             enddo
          enddo


          !Call the tangent linear convection scheme
          !DO I = 1,IM
           !  DO J = 1,JM
            !    if ( DOCONVEC(I,J) == 1) then

                   !call RASE_D (1, 1, LM, ICMIN, DT,                              &
                   !             MAPL8_CP, MAPL8_ALHL, MAPL8_GRAV, MAPL8_RGAS,     &
                   !             MAPL8_H2OMW, MAPL8_AIRMW, MAPL8_VIREPS,           &
                   !             SEEDRAS(I,J), SIGE,                               &
                   !             KCBL(I,J),                                        &
                   !             WGT0(I,J,:), WGT1(I,J,:),                         &
                   !             FRLAND(I,J), Ts(I,J),                             &
                   !             PTT(I,J,:), PTP(I,J,:),                           &
                   !             QVT(I,J,:), QVP(I,J,:),                           &
                   !             UT(I,J,:), UP(I,J,:),                             &
                   !             VT(I,J,:), VP(I,J,:),                             &
                   !             CO_AUTO(I,J), CNV_PLE(I,J,:),                     &
                   !             CNV_DQLDTT(I,J,:), CNV_DQLDTP(I,J,:),             &
                   !             CNV_MFDT(I,J,:),   CNV_MFDP(I,J,:),               &
                   !             CNV_PRC3T(I,J,:),  CNV_PRC3P(I,J,:),              &
                   !             CNV_UPDFT(I,J,:),  CNV_UPDFP(I,J,:),              &
                   !             RASPARAMS, ESTBLX                                 )

            !    endif
           !  enddo
          !enddo

          !Call the tangent linear cloud scheme.
          !CALL CLOUD_DRIVER_D ( DT, IM, JM, LM, PTT_C, PTP, QVT_C, QVP, PLE,        &
          !                      CNV_DQLDTT_C, CNV_DQLDTP, CNV_MFDT_C, CNV_MFDP,     &
          !                      CNV_PRC3T_C, CNV_PRC3P, CNV_UPDFT_C, CNV_UPDFP,     &
          !                      QILST, QILSP, QLLST, QLLSP,                         &
          !                      QICNT, QICNP, QLCNT, QLCNP,                         &
          !                      CFLST, CFLSP, CFCNT, CFCNP,                         &
          !                      FRLAND, CLDPARAMS, ESTBLX, KHu, KHl,              &
          !                      MAPL8_RUNIV, MAPL8_KAPPA, MAPL8_airmw, MAPL8_h2omw, &
          !                      MAPL8_GRAV, MAPL8_ALHL, MAPL8_ALHF, MAPL8_PI,       &
          !                      MAPL8_RGAS, MAPL8_CP, MAPL8_VIREPS, MAPL8_ALHS,     &
          !                      MAPL8_TICE, MAPL8_RVAP, MAPL8_P00, DO_MOIST_PHYS    )

       endif 

! Move perturbation tractory back into the original single precision containers
! -----------------------------------------------------------------------------
       UP4 = real(UP,4)
       VP4 = real(VP,4)
       if (PHASE==ADJphase) then
          PTP4 = real(PTP*(MAPL8_P00**MAPL8_KAPPA),4)
       elseif (PHASE==TLMphase) then 
          PTP4 = real(PTP/(MAPL8_P00**MAPL8_KAPPA),4)
       endif
       QVP4 = real(QVP,4)
       if (PHASE==ADJphase) then
          QIP4 = real(QILSP*ILSF + QICNP*ICNF,4)
          QLP4 = real(QLLSP*LLSF + QLCNP*LCNF,4)
       elseif (PHASE==TLMphase) then
          QIP4 = real(QILSP + QICNP,4)
          QLP4 = real(QLLSP + QLCNP,4)
       endif
       CFCNP4 = real(CFCNP,4)

! Deallocate memory used to store all the local arrays
! ----------------------------------------------------
       deallocate(UT,VT,PTT,QVT,UP,VP,PTP,QVP,PS,TS,FRLAND,KCBL,KHl,KHu)
       deallocate(ak,bk,pref,PLE,CNV_PLE,PLO,PK,PKE,TEMP)
       deallocate(PTT_F,QVT_F,PTT_L,QVT_L,PTT_C,QVT_C,DQST,QSST,DQSP,QSSP)
       deallocate(CNV_DQLDTT_C,CNV_MFDT_C,CNV_PRC3T_C,CNV_UPDFT_C)
       deallocate(CNV_DQLDTT,CNV_DQLDTP,CNV_MFDT,CNV_MFDP,CNV_PRC3T,CNV_PRC3P,CNV_UPDFT,CNV_UPDFP)
       deallocate(SIGE,WGT0,WGT1,CO_AUTO,RASAL2_2d,TPERTT,QPERTT,TPERTP,QPERTP)
       deallocate(SEEDRAS,DOCONVEC,HEAT,CTOP,sumHEAT,CNV_FRACTION,QV600)
       deallocate(JACOBIAN,H_pert,M_pert,PT_pert,QV_pert,PT_pert_in,QV_pert_in)
       deallocate(TPERTTj, TPERTPj, QPERTTj, DQSTj, DQSPj, QSSTj, QSSPj)
       deallocate(QILST,QLLST,QICNT,QLCNT,QILSP,QLLSP,QICNP,QLCNP)
       deallocate(CFLST,CFCNT,CFLSP,CFCNP)
       deallocate(fQi,ILSF,ICNF,LLSF,LCNF)

! Turn of the timers
! ------------------
       if (PHASE==ADJphase) then         
          call timing_off('MOIST_ADM')             
       elseif (PHASE==TLMphase) then         
          call timing_off('MOIST_TLM')            
       endif

    ENDIF !DO_MOIST_PHYS

! All done
! --------
    call MAPL_TimerOff(MAPL,PHASE_NAME) 
    call MAPL_TimerOff(MAPL,"TOTAL")
    call MAPL_TimerOff(MAPL,"-MOIST")

    RETURN_(ESMF_SUCCESS)

 contains


! Subroutines
! -----------

 subroutine jacobian_filter_tlm

 !Compute two specific columns of the Jacobian for a single profile. Then filter if values in those
 !columns are over a certain value, implying too steep gradients.

  JACOBIAN = 0.0
  
  DO L = 1,1 !Perturb level by level
 
     PT_pert = 0.0
     QV_pert = 0.0
 
     if (L == 1) then
 
        PT_pert(KCBL(I,J)) = 1.0 
 
     elseif (L == 2) then
 
        if (KCBL(I,J) == lm) then
          QV_pert(KCBL(I,J)) = 1.0 
        else
           QV_pert(KCBL(I,J) + 1) = 1.0 
        endif
 
     endif
 
     !SAVE PRECALL PROGNOSTIC VARIABLES
     PT_pert_in = PT_pert
     QV_pert_in = QV_pert
         
     !DO NONLINEAR CONVECTION TO CREATE INPUTS TRAJECTORY FOR LARGE SCALE ADJOINT
     CALL PRE_RASE_TLM( 1, 1, LM, PLE(I,J,:), PTT_F(I,J,:), PT_PERT, QVT_F(I,J,:), QV_PERT, TS(I,J), &
                        CNV_FRACTION(I,J), FRLAND(I,J), PLO(I,J,:), PKE(I,J,:), PK(I,J,:),           &
                        CBL_TPERT, CBL_QPERT, CBL_TPERT_MXOCN, CBL_TPERT_MXLND, &
                        TPERTTj, TPERTPj, QPERTTj, DQSTj, DQSPj, QSSTj, QSSPj )

      call RASE_FAST_TLM( 1                    , &
                          1                    , &                 
                          LM                   , &
                          ICMIN                , &
                          DT                   , &
                          MAPL8_CP             , &
                          MAPL8_ALHL           , &
                          MAPL8_ALHS           , &
                          MAPL8_TICE           , &
                          MAPL8_GRAV           , &
                          SEEDRAS(I,J,:)       , &
                          SIGE                 , &
                          KCBL(I,J)            , &
                          WGT0(I,J,:)          , &
                          WGT1(I,J,:)          , &
                          TPERTTj              , &
                          TPERTPj              , &
                          QPERTTj              , &
                          PTT_F(I,J,:)         , &
                          PT_PERT              , &
                          QVT_F(I,J,:)         , & 
                          QV_PERT              , &
                          QSSTj                , & 
                          QSSPj                , & 
                          DQSTj                , &
                          DQSPj                , &
                          CNV_FRACTION(I,J)    , &
                          RASAL2_2d(I,J)       , &
                          CO_AUTO(I,J)         , &
                          CNV_PLE(I,J,:)       , &
                          PKE(I,J,:)           , &
                          RASPARAMS              )
 
     !COMPUTE PERTURBATION HEATING AND MOISTENING RATES
     H_pert = (PT_pert - PT_pert_in)/DT
     M_pert = (QV_pert - QV_pert_in)/DT
       
     !Uncomment here if just doing two columns of the Jacobian
     if (L == 1) then
        JACOBIAN(0*LM+1:1*LM,1) = H_pert
        JACOBIAN(1*LM+1:2*LM,1) = M_pert
     elseif (L == 2) then
        JACOBIAN(0*LM+1:1*LM,2) = H_pert
        JACOBIAN(1*LM+1:2*LM,2) = M_pert
     endif
 
  endDO

  !Constants here determined so as to remove as many of the problematic points as possible.
  !The constants used in this if loop are tuned from looking at many Jacobians for many time steps. Values are choosen
  !so as to balance between keeping the natural behahiour for as many points as possible without 
  IF ( (maxval(abs(Jacobian(1:lm     ,1))) .gt. 0.00010  ) .or. &
       (maxval(abs(Jacobian(1:lm     ,2))) .gt. 0.25     ) .or. & 
       (maxval(abs(Jacobian(lm+1:2*lm,1))) .gt. 1.0e-07  ) .or. &
       (maxval(abs(Jacobian(lm+1:2*lm,2))) .gt. 0.000250 ) ) then

     DOCONVEC(I,J) = 0

  else

     DOCONVEC(I,J) = 1

  endIF 

 endsubroutine jacobian_filter_tlm

 end subroutine Run

 subroutine PRE_RASE( IM,JM,LM,PLE,TH1,Q1,TS,CNV_FRACTION,FRLAND, &
                      PLO,PKE,PK, &
                      CBL_TPERT,CBL_QPERT,CBL_TPERT_MXOCN,CBL_TPERT_MXLND, &
                      TPERT,QPERT,DQS,QSS )

  implicit none

  integer, intent(in) :: IM,JM,LM
  real(8), intent(in ), dimension(IM,JM,0:LM) :: PLE
  real(8), intent(in ), dimension(IM,JM,LM)   :: TH1, Q1
  real(8), intent(in ), dimension(IM,JM)      :: TS, CNV_FRACTION, FRLAND
  real(8), intent(in )                        :: CBL_TPERT,CBL_QPERT, CBL_TPERT_MXOCN, CBL_TPERT_MXLND
  real(8), intent(in ), dimension(IM,JM,0:LM) :: PKE
  real(8), intent(in ), dimension(IM,JM,LM)   :: PLO, PK

  real(8), intent(out), dimension(IM,JM,LM)   :: DQS, QSS
  real(8), intent(out), dimension(IM,JM)      :: TPERT, QPERT

  integer                                     :: i,j,l
  real(8),              dimension(IM,JM,0:LM) :: ZLE
  real(8),              dimension(IM,JM,LM)   :: TEMP1, ZLO

   TEMP1     = TH1*PK

   ZLE(:,:,LM) = 0.
   do L=LM,1,-1
      ZLE(:,:,L-1) = TH1(:,:,L) * (1.+MAPL_VIREPS*Q1(:,:,L))
      ZLO(:,:,L  ) = ZLE(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
      ZLE(:,:,L-1) = ZLO(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
   end do

   do J=1,JM
      do I=1,IM
         call DQSATpert(DQS(i,j,:),QSS(i,j,:),TEMP1(i,j,:),PLO(i,j,:),LM)
      end do
   end do

   TPERT  = ABS(CBL_TPERT) * ( TS - ( TEMP1(:,:,LM)+ MAPL_GRAV*ZLO(:,:,LM)/MAPL_CP )  ) 
   if (CBL_TPERT < 0) then
      ! Make TPERT 0 in areas of deep convection
      TPERT = TPERT*(1.0-CNV_FRACTION)
   endif
   QPERT  = 0.0 !CBL_QPERT * ( QSSFC - Q(:,:,LM) )      !dh: CBL_QPERT = 0.0
   TPERT  = MAX( TPERT , 0.0 )
   QPERT  = MAX( QPERT , 0.0 )

   where (FRLAND<0.1) 
      TPERT = MIN( TPERT , CBL_TPERT_MXOCN ) ! ocean
   elsewhere
      TPERT = MIN( TPERT , CBL_TPERT_MXLND ) ! land
   end where

 end subroutine PRE_RASE

  SUBROUTINE PRE_RASE_TLM(im, jm, lm, ple, th1, th1_tl, q1, q1_tl, ts, &
&   cnv_fraction, frland, plo, pke, pk, cbl_tpert, cbl_qpert, &
&   cbl_tpert_mxocn, cbl_tpert_mxlnd, tpert, tpert_tl, qpert, dqs, &
&   dqs_tl, qss, qss_tl)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im, jm, lm
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: ple
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: th1, q1
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: th1_tl, q1_tl
    REAL*8, DIMENSION(im, jm), INTENT(IN) :: ts, cnv_fraction, frland
    REAL*8, INTENT(IN) :: cbl_tpert, cbl_qpert, cbl_tpert_mxocn, &
&   cbl_tpert_mxlnd
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: pke
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: plo, pk
    REAL*8, DIMENSION(im, jm, lm), INTENT(OUT) :: dqs, qss
    REAL*8, DIMENSION(im, jm, lm), INTENT(OUT) :: dqs_tl, qss_tl
    REAL*8, DIMENSION(im, jm), INTENT(OUT) :: tpert, qpert
    REAL*8, DIMENSION(im, jm), INTENT(OUT) :: tpert_tl
    INTEGER :: i, j, l
    REAL*8, DIMENSION(im, jm, 0:lm) :: zle
    REAL*8, DIMENSION(im, jm, 0:lm) :: zle_tl
    REAL*8, DIMENSION(im, jm, lm) :: temp1, zlo
    REAL*8, DIMENSION(im, jm, lm) :: temp1_tl, zlo_tl
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    REAL*8 :: abs0
    temp1_tl = pk*th1_tl
    temp1 = th1*pk
    zle(:, :, lm) = 0.
    zlo_tl = 0.0_8
    zle_tl = 0.0_8
    DO l=lm,1,-1
      zle_tl(:, :, l-1) = th1_tl(:, :, l)*(1.+mapl_vireps*q1(:, :, l)) +&
&       th1(:, :, l)*mapl_vireps*q1_tl(:, :, l)
      zle(:, :, l-1) = th1(:, :, l)*(1.+mapl_vireps*q1(:, :, l))
      zlo_tl(:, :, l) = zle_tl(:, :, l) + mapl_cp*(pke(:, :, l)-pk(:, :&
&       , l))*zle_tl(:, :, l-1)/mapl_grav
      zlo(:, :, l) = zle(:, :, l) + mapl_cp/mapl_grav*(pke(:, :, l)-pk(:&
&       , :, l))*zle(:, :, l-1)
      zle_tl(:, :, l-1) = zlo_tl(:, :, l) + mapl_cp*(pk(:, :, l)-pke(:, &
&       :, l-1))*zle_tl(:, :, l-1)/mapl_grav
      zle(:, :, l-1) = zlo(:, :, l) + mapl_cp/mapl_grav*(pk(:, :, l)-pke&
&       (:, :, l-1))*zle(:, :, l-1)
    END DO
    dqs_tl = 0.0_8
    qss_tl = 0.0_8
    DO j=1,jm
      DO i=1,im
        CALL DQSATPERT_TLM(dqs(i, j, :), dqs_tl(i, j, :), qss(i, j, :), &
&                    qss_tl(i, j, :), temp1(i, j, :), temp1_tl(i, j, :)&
&                    , plo(i, j, :), lm)
      END DO
    END DO
    IF (cbl_tpert .GE. 0.) THEN
      abs0 = cbl_tpert
    ELSE
      abs0 = -cbl_tpert
    END IF
    tpert_tl = abs0*(-temp1_tl(:, :, lm)-mapl_grav*zlo_tl(:, :, lm)/&
&     mapl_cp)
    tpert = abs0*(ts-(temp1(:, :, lm)+mapl_grav*zlo(:, :, lm)/mapl_cp))
    IF (cbl_tpert .LT. 0) THEN
! Make TPERT 0 in areas of deep convection
      tpert_tl = (1.0-cnv_fraction)*tpert_tl
      tpert = tpert*(1.0-cnv_fraction)
    END IF
!CBL_QPERT * ( QSSFC - Q(:,:,LM) )      !dh: CBL_QPERT = 0.0
    qpert = 0.0
    WHERE (tpert .LT. 0.0) 
      tpert_tl = 0.0_8
      tpert = 0.0
    ELSEWHERE
      tpert = tpert
    END WHERE
    WHERE (qpert .LT. 0.0) 
      qpert = 0.0
    ELSEWHERE
      qpert = qpert
    END WHERE
    WHERE (tpert .GT. cbl_tpert_mxocn) 
      tpert_tl = 0.0_8
      tpert = cbl_tpert_mxocn
    ELSEWHERE
      tpert = tpert
    END WHERE
    WHERE (tpert .GT. cbl_tpert_mxlnd) 
      tpert_tl = 0.0_8
      tpert = cbl_tpert_mxlnd
    ELSEWHERE
      tpert = tpert
    END WHERE
  END SUBROUTINE PRE_RASE_TLM
  SUBROUTINE PRE_RASE_FWD(im, jm, lm, ple, th1, q1, ts, cnv_fraction, &
&   frland, plo, pke, pk, cbl_tpert, cbl_qpert, cbl_tpert_mxocn, &
&   cbl_tpert_mxlnd, tpert, qpert, dqs, qss)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im, jm, lm
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: ple
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: th1, q1
    REAL*8, DIMENSION(im, jm), INTENT(IN) :: ts, cnv_fraction, frland
    REAL*8, INTENT(IN) :: cbl_tpert, cbl_qpert, cbl_tpert_mxocn, &
&   cbl_tpert_mxlnd
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: pke
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: plo, pk
    REAL*8, DIMENSION(im, jm, lm) :: dqs, qss
    REAL*8, DIMENSION(im, jm) :: tpert, qpert
    INTEGER :: i, j, l
    REAL*8, DIMENSION(im, jm, 0:lm) :: zle
    REAL*8, DIMENSION(im, jm, lm) :: temp1, zlo
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    LOGICAL :: mask3(im, jm)
    LOGICAL :: mask2(im, jm)
    LOGICAL :: mask1(im, jm)
    LOGICAL :: mask0(im, jm)
    REAL*8 :: abs0
    LOGICAL :: mask(im, jm)
    temp1 = th1*pk
    zle(:, :, lm) = 0.
    DO l=lm,1,-1
      zle(:, :, l-1) = th1(:, :, l)*(1.+mapl_vireps*q1(:, :, l))
      zlo(:, :, l) = zle(:, :, l) + mapl_cp/mapl_grav*(pke(:, :, l)-pk(:&
&       , :, l))*zle(:, :, l-1)
      zle(:, :, l-1) = zlo(:, :, l) + mapl_cp/mapl_grav*(pk(:, :, l)-pke&
&       (:, :, l-1))*zle(:, :, l-1)
    END DO
    DO j=1,jm
      DO i=1,im
        CALL DQSATPERT(dqs(i, j, :), qss(i, j, :), temp1(i, j, :), plo(i&
&                , j, :), lm)
      END DO
    END DO
    IF (cbl_tpert .GE. 0.) THEN
      abs0 = cbl_tpert
    ELSE
      abs0 = -cbl_tpert
    END IF
    tpert = abs0*(ts-(temp1(:, :, lm)+mapl_grav*zlo(:, :, lm)/mapl_cp))
    IF (cbl_tpert .LT. 0) THEN
! Make TPERT 0 in areas of deep convection
      tpert = tpert*(1.0-cnv_fraction)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
!CBL_QPERT * ( QSSFC - Q(:,:,LM) )      !dh: CBL_QPERT = 0.0
    qpert = 0.0
    mask = tpert .LT. 0.0
    WHERE (mask) 
      tpert = 0.0
    ELSEWHERE
      tpert = tpert
    END WHERE
    mask0 = qpert .LT. 0.0
    WHERE (mask0) 
      qpert = 0.0
    ELSEWHERE
      qpert = qpert
    END WHERE
    mask1 = frland .LT. 0.1
    mask2 = tpert .GT. cbl_tpert_mxocn
    WHERE (mask2) 
      tpert = cbl_tpert_mxocn
    ELSEWHERE
      tpert = tpert
    END WHERE
    mask3 = tpert .GT. cbl_tpert_mxlnd
    WHERE (mask3) 
      tpert = cbl_tpert_mxlnd
    ELSEWHERE
      tpert = tpert
    END WHERE
    CALL PUSHBOOLEANARRAY(mask, im*jm)
    CALL PUSHREAL8(abs0)
    CALL PUSHBOOLEANARRAY(mask0, im*jm)
    CALL PUSHBOOLEANARRAY(mask1, im*jm)
    CALL PUSHBOOLEANARRAY(mask2, im*jm)
    CALL PUSHBOOLEANARRAY(mask3, im*jm)
  END SUBROUTINE PRE_RASE_FWD
!  Differentiation of pre_rase in reverse (adjoint) mode, backward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS_
!MoistGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE qs
!at_util.DQSAT cloudnew.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensa
!te cloudnew.cnvsrc cloudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudne
!w.MARSHPALMQ2 cloudnew.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ
!3_ICE3)):
!   gradient     of useful results: dqs tpert th1 q1 qss
!   with respect to varying inputs: th1 q1
  SUBROUTINE PRE_RASE_BWD(im, jm, lm, ple, th1, th1_ad, q1, q1_ad, ts, &
&   cnv_fraction, frland, plo, pke, pk, cbl_tpert, cbl_qpert, &
&   cbl_tpert_mxocn, cbl_tpert_mxlnd, tpert, tpert_ad, qpert, dqs, &
&   dqs_ad, qss, qss_ad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im, jm, lm
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: ple
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: th1, q1
    REAL*8, DIMENSION(im, jm, lm) :: th1_ad, q1_ad
    REAL*8, DIMENSION(im, jm), INTENT(IN) :: ts, cnv_fraction, frland
    REAL*8, INTENT(IN) :: cbl_tpert, cbl_qpert, cbl_tpert_mxocn, &
&   cbl_tpert_mxlnd
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: pke
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: plo, pk
    REAL*8, DIMENSION(im, jm, lm) :: dqs, qss
    REAL*8, DIMENSION(im, jm, lm) :: dqs_ad, qss_ad
    REAL*8, DIMENSION(im, jm) :: tpert, qpert
    REAL*8, DIMENSION(im, jm) :: tpert_ad
    INTEGER :: i, j, l
    REAL*8, DIMENSION(im, jm, 0:lm) :: zle
    REAL*8, DIMENSION(im, jm, 0:lm) :: zle_ad
    REAL*8, DIMENSION(im, jm, lm) :: temp1, zlo
    REAL*8, DIMENSION(im, jm, lm) :: temp1_ad, zlo_ad
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: branch
    LOGICAL :: mask3(im, jm)
    LOGICAL :: mask2(im, jm)
    LOGICAL :: mask1(im, jm)
    LOGICAL :: mask0(im, jm)
    REAL*8 :: abs0
    LOGICAL :: mask(im, jm)
    CALL POPBOOLEANARRAY(mask3, im*jm)
    CALL POPBOOLEANARRAY(mask2, im*jm)
    CALL POPBOOLEANARRAY(mask1, im*jm)
    CALL POPBOOLEANARRAY(mask0, im*jm)
    CALL POPREAL8(abs0)
    CALL POPBOOLEANARRAY(mask, im*jm)
    WHERE (mask3) tpert_ad = 0.0_8
    WHERE (mask2) tpert_ad = 0.0_8
    WHERE (mask) tpert_ad = 0.0_8
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) tpert_ad = (1.0-cnv_fraction)*tpert_ad
    temp1 = th1*pk
    zlo_ad = 0.0_8
    temp1_ad = 0.0_8
    temp1_ad(:, :, lm) = temp1_ad(:, :, lm) - abs0*tpert_ad
    zlo_ad(:, :, lm) = zlo_ad(:, :, lm) - mapl_grav*abs0*tpert_ad/&
&     mapl_cp
    DO j=jm,1,-1
      DO i=im,1,-1
        CALL DQSATPERT_ADM(dqs(i, j, :), dqs_ad(i, j, :), qss(i, j, :), &
&                    qss_ad(i, j, :), temp1(i, j, :), temp1_ad(i, j, :)&
&                    , plo(i, j, :), lm)
      END DO
    END DO
    zle_ad = 0.0_8
    DO l=1,lm,1
      zlo_ad(:, :, l) = zlo_ad(:, :, l) + zle_ad(:, :, l-1)
      zle_ad(:, :, l-1) = mapl_cp*(pk(:, :, l)-pke(:, :, l-1))*zle_ad(:&
&       , :, l-1)/mapl_grav
      zle_ad(:, :, l) = zle_ad(:, :, l) + zlo_ad(:, :, l)
      zle_ad(:, :, l-1) = zle_ad(:, :, l-1) + mapl_cp*(pke(:, :, l)-pk(:&
&       , :, l))*zlo_ad(:, :, l)/mapl_grav
      zlo_ad(:, :, l) = 0.0_8
      th1_ad(:, :, l) = th1_ad(:, :, l) + (mapl_vireps*q1(:, :, l)+1.)*&
&       zle_ad(:, :, l-1)
      q1_ad(:, :, l) = q1_ad(:, :, l) + mapl_vireps*th1(:, :, l)*zle_ad(&
&       :, :, l-1)
      zle_ad(:, :, l-1) = 0.0_8
    END DO
    th1_ad = th1_ad + pk*temp1_ad
  END SUBROUTINE PRE_RASE_BWD

  subroutine VertInterp(v2,v3,ple,pp,MAPL8_UNDEF,rc)

    real(8)    , intent(OUT) :: v2(:,:)
    real(8)    , intent(IN ) :: v3(:,:,:)
    real(8)    , intent(IN ) :: ple(:,:,:)
    real(8)    , intent(IN ) :: pp, MAPL8_UNDEF
    integer, optional, intent(OUT) :: rc

    real(8), dimension(size(v2,1),size(v2,2)) :: al,PT,PB
    integer k,km
    logical edge

    character*(10) :: Iam='VertInterp'

    km   = size(ple,3)-1
    edge = size(v3,3)==km+1

    ASSERT_(edge .or. size(v3,3)==km)

    v2   = MAPL8_UNDEF

    if(EDGE) then
       pb   = ple(:,:,km+1)
       do k=km,1,-1
          pt = ple(:,:,k)
          if(all(pb<pp)) exit
          where(pp>pt .and. pp<=pb)
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k)*al + v3(:,:,k+1)*(1.0-al)
          end where
          pb = pt
       end do
    else
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
       do k=km,2,-1
          pt = 0.5*(ple(:,:,k-1)+ple(:,:,k))
          if(all(pb<pp)) exit
          where( (pp>pt.and.pp<=pb) )
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k-1)*al + v3(:,:,k)*(1.0-al)
          end where
          pb = pt
       end do
       pt = 0.5*(ple(:,:,km)+ple(:,:,km-1))
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
          where( (pp>pb.and.pp<=ple(:,:,km+1)) )
             v2 = v3(:,:,km)
          end where
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine VertInterp

end module GEOS_MoistPertGridCompMod


