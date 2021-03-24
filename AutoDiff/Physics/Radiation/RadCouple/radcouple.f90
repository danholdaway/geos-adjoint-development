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
 real, intent(in ) :: TE, PL, TEMPOR
 real, intent(in ) :: AF, CF, QClAN, QCiAN, QClLS, QCiLS
! real, intent(in ) :: QRN_ALL, QSN_ALL

 !Outputs
 real, intent(out) :: RAD_QL,RAD_QI,RAD_CF,RAD_RL,RAD_RI
! real, intent(out) :: RAD_QR,RAD_QS

 !Locals
 real :: ss, RAD_RI_AN, AFx, ALPH

 real, parameter :: MIN_RI = 20.e-6, MAX_RI = 40.e-6, RI_ANV = 30.e-6

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

