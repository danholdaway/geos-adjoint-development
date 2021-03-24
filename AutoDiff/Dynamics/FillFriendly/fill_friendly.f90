subroutine FILL_Friendly ( IM,JM,LM,Q,DP )

 IMPLICIT NONE

 !Inputs
 INTEGER, INTENT(IN) :: IM, JM, LM
 REAL(8), INTENT(IN) :: DP(IM,JM,LM)

 !Prognostic
 REAL(8), INTENT(INOUT) :: Q(IM,JM,LM)

 !Locals
 INTEGER :: I, J, L
 REAL(8) :: QTEMP1(IM,JM), QTEMP2(IM,JM)


 QTEMP1 = 0.0
 do L=1,LM
    QTEMP1(:,:) = QTEMP1(:,:) + Q(:,:,L)*DP(:,:,L)
 enddo 

 do i = 1,IM
    do j = 1,JM
       do l = 1,LM
          if (Q(i,j,l) < 0.0) then

             Q(i,j,l) = 0.0

          endif
       enddo
    enddo
 enddo

 QTEMP2 = 0.0
 do L=1,LM
    QTEMP2(:,:) = QTEMP2(:,:) + Q(:,:,L)*DP(:,:,L)
 enddo 

 do i = 1,IM
    do j = 1,JM
          if (abs(qtemp2(i,j)) > 0.0) then

             qtemp2(i,j) = max( qtemp1(i,j)/qtemp2(i,j), 0.0 )

          endif
    enddo
 enddo

 do L=1,LM
    Q(:,:,L) = Q(:,:,L)*qtemp2(:,:)
 enddo 

endsubroutine FILL_FRIENDLY
