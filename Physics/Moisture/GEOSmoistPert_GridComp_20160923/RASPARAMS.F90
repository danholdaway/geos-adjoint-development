      module RASPARAMS

      implicit none

! Define RASPARAM
! ---------------
      type RASPARAM_TYPE
           real(8)            :: CUFRICFAC             ! 1
           real(8)            :: SHR_LAMBDA_FAC        ! 2

           real(8)            :: QC_CRIT_CN            ! 4
           real(8)            :: RASAL1                ! 5
           real(8)            :: RASAL2                ! 6
           real(8)            :: RASNCL                ! 7
           real(8)            :: LAMBDA_FAC            ! 8
           real(8)            :: LAMBMX_FAC            ! 9
           real(8)            :: MIN_DIAMETER          ! 10
           real(8)            :: CUFRICLAMBDA          ! 11
           real(8)            :: RDTLEXPON             ! 12
           real(8)            :: STRAPPING             ! 13
           real(8)            :: SDQV2                 ! 14
           real(8)            :: SDQV3                 ! 15
           real(8)            :: SDQVT1                ! 16
           real(8)            :: ACRITFAC              ! 17
           real(8)            :: HMINTRIGGER           ! 18
           real(8)            :: LLDISAGGXP            ! 19
           real(8)            :: PBLFRAC               ! 20
           real(8)            :: RASAUTORAMPB          ! 21
           real(8)            :: AUTOC_CN_ZDEP         ! 22
           real(8)            :: MAXDALLOWED_S         ! 23
           real(8)            :: MAXDALLOWED_D         ! 24
           integer            :: RASAL_EXP
           real(8)            :: RAS_RHMIN             ! 25
           real(8)            :: RAS_RHFULL            ! 26
           real(8)            :: CLDMICRO              ! 27
           real(8)            :: FDROP_DUST            ! 28
           real(8)            :: FDROP_SOOT            ! 29
      endtype RASPARAM_TYPE

  end module RASPARAMS
