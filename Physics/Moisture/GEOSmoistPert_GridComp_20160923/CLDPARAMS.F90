      module CLDPARAMS

      implicit none

! Define CLDPARAM
! ---------------
      type CLDPARAM_TYPE
           real(8)            :: CNV_BETA              ! 1
           real(8)            :: ANV_BETA              ! 2
           real(8)            :: LS_BETA               ! 3
           real(8)            :: RH_CRIT               ! 4
           real(8)            :: AUTOC_LS              ! 5
           real(8)            :: QC_CRIT_LS            ! 6
           real(8)            :: ACCRETION             ! 7
           real(8)            :: RAIN_REVAP_FAC        ! 8
           real(8)            :: VOL_TO_FRAC           ! 9
           real(8)            :: SUPERSAT              ! 10
           real(8)            :: SHEAR_EVAP_FAC        ! 11
           real(8)            :: MIN_ALLOW_CCW         ! 12
           real(8)            :: CCW_EVAP_EFF          ! 13
           real(8)            :: NSUB_AUTOCONV         ! 14
           real(8)            :: LS_SUND_INTER         ! 15
           real(8)            :: LS_SUND_COLD          ! 16
           real(8)            :: LS_SUND_TEMP1         ! 17
           real(8)            :: ANV_SUND_INTER        ! 18
           real(8)            :: ANV_SUND_COLD         ! 19
           real(8)            :: ANV_SUND_TEMP1        ! 20
           real(8)            :: ANV_TO_LS_TIME        ! 21
           real(8)            :: NCCN_WARM             ! 22
           real(8)            :: NCCN_ICE              ! 23
           real(8)            :: NCCN_ANVIL            ! 24
           real(8)            :: NCCN_PBL              ! 25
           real(8)            :: DISABLE_RAD           ! 26
           real(8)            :: ICE_SETTLE            ! 27
           real(8)            :: ANV_ICEFALL           ! 28
           real(8)            :: LS_ICEFALL            ! 29
           real(8)            :: REVAP_OFF_P           ! 30
           real(8)            :: CNV_ENVF              ! 31
           real(8)            :: WRHODEP               ! 32
           real(8)            :: ICE_RAMP              ! 33
           real(8)            :: CNV_ICEPARAM          ! 34
           real(8)            :: CNV_ICEFRPWR          ! 35
           real(8)            :: CNV_DDRF              ! 36
           real(8)            :: ANV_DDRF              ! 37
           real(8)            :: LS_DDRF               ! 38
           real(8)            :: AUTOC_ANV             ! 39
           real(8)            :: QC_CRIT_ANV           ! 40
           real(8)            :: TANHRHCRIT            ! 41
           real(8)            :: MINRHCRIT             ! 42
           real(8)            :: MAXRHCRIT             ! 43
           real(8)            :: PRECIPRAD             ! 44
           real(8)            :: TURNRHCRIT            ! 45
           real(8)            :: MAXRHCRITLAND         ! 46
           real(8)            :: FR_LS_WAT             ! 47
           real(8)            :: FR_LS_ICE             ! 48
           real(8)            :: FR_AN_WAT             ! 49
           real(8)            :: FR_AN_ICE             ! 50
           real(8)            :: MIN_RL                ! 51
           real(8)            :: MIN_RI                ! 52
           real(8)            :: MAX_RL                ! 53
           real(8)            :: MAX_RI                ! 54
           real(8)            :: RI_ANV                ! 55
           real(8)            :: SNOW_REVAP_FAC        ! 56
           real(8)            :: PDFSHAPE              ! 57
           real(8)            :: TURNRHCRIT_UP         ! 58
           real(8)            :: MOVE2RAS              ! 59
      endtype CLDPARAM_TYPE

  end module CLDPARAMS
