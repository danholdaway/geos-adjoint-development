subroutine G5TDSOLVER(IM,JM,LM,A,B,C,Y)

   implicit none

   !Arguments
   integer,                      intent(IN   ) :: IM, JM, LM
   real(8),    dimension(IM,JM,LM), intent(IN   ) :: A, B, C
   real(8),    dimension(IM,JM,LM), intent(INOUT) :: Y

   !Locals
   integer :: L


! Sweep down, modifying rhs with multiplier A

    do L = 2,LM
       Y(:,:,L) = Y(:,:,L) - Y(:,:,L-1) * A(:,:,L)
    enddo

! Sweep up, solving for updated value. Note B has the inverse of the main diagonal

    Y(:,:,LM)   =  Y(:,:,LM)*B(:,:,LM-1)/(B(:,:,LM-1) - A(:,:,LM)*(1.0+C(:,:,LM-1)*B(:,:,LM-1) ))

    do L = LM-1,1,-1
       Y(:,:,L) = (Y(:,:,L ) - C(:,:,L ) * Y(:,:,L+1))*B(:,:,L )
    enddo
    
    

end subroutine G5TDSOLVER
