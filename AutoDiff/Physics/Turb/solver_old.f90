subroutine G5TDSOLVER(IM,JM,LM,A,B,C,Y,PHASE)

   implicit none

   !Arguments
   integer,                      intent(IN   ) :: IM, JM, LM, PHASE
   real,    dimension(IM,JM,LM), intent(IN   ) :: A, B, C
   real,    dimension(IM,JM,LM), intent(INOUT) :: Y

   !Locals
   integer :: L
   real, dimension(IM,JM,LM) :: A1, B1, C1


    !Make temporary copies and dont overwrite
    A1 = A
    B1 = B
    C1 = C

    !VTRILU
    !------
    B1(:,:,1) = 1. / B1(:,:,1)
    do L = 2,LM
       A1(:,:,L) = A1(:,:,L) * B1(:,:,L-1)
       B1(:,:,L) = 1. / ( B1(:,:,L) - C1(:,:,L-1) * A1(:,:,L) )
    end do


    !VTRISOLVE
    !---------
    if (PHASE == TLMPhase) then

       !Sweep 1
       do L = 2,LM
          Y(:,:,L) = Y(:,:,L) - Y(:,:,L-1) * A1(:,:,L)
       enddo

       !Sweep 2
       Y(:,:,LM) = Y(:,:,LM)*B1(:,:,LM-1)/(B1(:,:,LM-1) - A1(:,:,LM)*(1.0+C1(:,:,LM-1)*B1(:,:,LM-1) ))
       do L = LM-1,1,-1
          Y(:,:,L) = (Y(:,:,L ) - C1(:,:,L ) * Y(:,:,L+1))*B1(:,:,L )
       enddo

    elseif (PHASE == ADJPhase) then)

       !Adjoint of sweep 2
       do l=1,lm-1,1
          TMP = B1(:,:,L)*Y(:,:,L)
          Y(:,:,L+1) = Y(:,:,L+1) - C1(:,:,L)*TMP
          Y(:,:,L)   = TMP
       enddo
       Y(:,:,lm) = B1(:,:,LM-1)*Y(:,:,LM)/(B1(:,:,LM-1)-A1(:,:,LM)*(C1(:,:,LM-1)*B1(:,:,LM-1)+1.0))

       !Adjoint of sweep 1
       do L=lm,2,-1
          Y(:,:,L-1) = Y(:,:,L-1) - A1(:,:,L)*Y(:,:,L)
       enddo

    endif

end subroutine G5TDSOLVER
