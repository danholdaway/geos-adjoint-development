subroutine aeroopt ( IM,JM,LM,NA,NBCHOU_IR,NBCHOU_SO,AEROT,PLE,                   &
                     TAUA_IR_C,SSAA_IR_C,ASYA_IR_C,TAUA_IRT,SSAA_IRT,ASYA_IRT, &
                     TAUA_SO_C,SSAA_SO_C,ASYA_SO_C,TAUA_SOT,SSAA_SOT,ASYA_SOT  )  

implicit none

integer, intent(in) :: IM, JM, LM
integer, intent(in) :: NA, NBCHOU_IR, NBCHOU_SO
real, parameter :: MAPL_GRAV = 9.81

!Inputs
real, intent(in), dimension(IM,JM,LM,NA) :: AEROT
real, intent(in), dimension(IM,JM,0:LM)  :: PLE

real, intent(in), dimension(IM,JM,LM,NBCHOU_IR,NA) :: TAUA_IR_C,SSAA_IR_C,ASYA_IR_C
real, intent(in), dimension(IM,JM,LM,NBCHOU_SO,NA) :: TAUA_SO_C,SSAA_SO_C,ASYA_SO_C

!Outputs
real, intent(out), dimension(IM,JM,LM,NBCHOU_IR) :: TAUA_IRT,SSAA_IRT,ASYA_IRT
real, intent(out), dimension(IM,JM,LM,NBCHOU_SO) :: TAUA_SOT,SSAA_SOT,ASYA_SOT

!Locals
integer :: I, J, L, AI
real :: x
real, dimension(IM,JM,LM,NA) :: SAEROT

 DO L = 1, LM
    DO J = 1, JM
       DO I = 1, IM
          X = ((PLE(I,J,L) - PLE(I,J,L-1))*0.01)*(100./MAPL_GRAV)
          DO AI = 1, NA
             SAEROT(I,J,L,AI) = X*AEROT(I,J,L,AI)
          END DO
       END DO
    END DO
 END DO

 !Compute aerosol optical properties for IR
 TAUA_IRT = 0.0
 SSAA_IRT = 0.0
 ASYA_IRT = 0.0
 do J = 1,NBCHOU_IR
    do L = 1,NA
       TAUA_IRT(:,:,:,J) = TAUA_IRT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)
       SSAA_IRT(:,:,:,J) = SSAA_IRT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)*SSAA_IR_C(:,:,:,J,L)
       ASYA_IRT(:,:,:,J) = ASYA_IRT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)*SSAA_IR_C(:,:,:,J,L)*ASYA_IR_C(:,:,:,J,L)
    enddo
 enddo

 !Compute aerosol optical properties for SO
 TAUA_SOT = 0.0
 SSAA_SOT = 0.0
 ASYA_SOT = 0.0
 do J = 1,NBCHOU_SO
    do L = 1,NA
       TAUA_SOT(:,:,:,J) = TAUA_SOT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)
       SSAA_SOT(:,:,:,J) = SSAA_SOT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)*SSAA_SO_C(:,:,:,J,L)
       ASYA_SOT(:,:,:,J) = ASYA_SOT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)*SSAA_SO_C(:,:,:,J,L)*ASYA_IR_C(:,:,:,J,L)
    enddo
 enddo


end subroutine aeroopt
