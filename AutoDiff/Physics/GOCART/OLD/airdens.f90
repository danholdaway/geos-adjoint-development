subroutine dummy ( im, jm, lm, rho, P, PT, QV, &
                     CONS_RVAP, CONS_RGAS, CONS_P00, CONS_KAPPA, CONS_CP )

 implicit none

 !Inputs
 integer, intent(in)  :: im, jm, lm        ! grid
 real(8), intent(in)  :: P(im,jm,0:lm)    ! pressure edges
 real(8), intent(in)  :: PT(im,jm,lm)      ! (dry) potential temperature
 real(8), intent(in)  :: QV(im,jm,lm)       ! specific humidity
 real(8), intent(in)  :: CONS_RVAP, CONS_RGAS, CONS_P00, CONS_KAPPA, CONS_CP

 !Outputs
 real(8), intent(out) :: rho(im,jm,lm)     ! air density [kg/m3]

PT = QV*P*PT
QV = QV**2*P

call airdens ( im, jm, lm, rho, P, PT, QV, &
                     CONS_RVAP, CONS_RGAS, CONS_P00, CONS_KAPPA, CONS_CP )


rho = rho*QV*PT*P

endsubroutine dummy

subroutine airdens ( im, jm, lm, rho, P, PT, QV, &
                     CONS_RVAP, CONS_RGAS, CONS_P00, CONS_KAPPA, CONS_CP )

 implicit none

 !Inputs
 integer, intent(in)  :: im, jm, lm        ! grid
 real(8), intent(in)  :: P(im,jm,0:lm)    ! pressure edges
 real(8), intent(in)  :: PT(im,jm,lm)      ! (dry) potential temperature
 real(8), intent(in)  :: QV(im,jm,lm)       ! specific humidity
 real(8), intent(in)  :: CONS_RVAP, CONS_RGAS, CONS_P00, CONS_KAPPA, CONS_CP

 !Outputs
 real(8), intent(out) :: rho(im,jm,lm)     ! air density [kg/m3]

 !Locals
 integer :: l
 real(8) :: eps, tmp(im,jm)
 real(8) :: npk(im,jm,0:lm)                 ! normalized pk = (P/p0)^kappa

 eps = CONS_RVAP / CONS_RGAS - 1.0

 !Compute normalized P**Kappa
 npk = (P/CONS_P00)**CONS_KAPPA

 !Compute rho from hydrostatic equation
 do l = 1, lm
    tmp = PT(:,:,l)*(1.+eps*QV(:,:,l))
    rho(:,:,l) = (P(:,:,l)-P(:,:,l-1))/( CONS_CP*tmp*(npk(:,:,l)-npk(:,:,l-1)))
 end do

end subroutine airdens
