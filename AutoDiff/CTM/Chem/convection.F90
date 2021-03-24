  SUBROUTINE convection(i1, i2, j1, j2, km, n1, n2, dt_30m, aero_type, kin, &
                        tc, cldmas, dtrain, area, delz, delp, vud, &
                        airmass, airmol, tmpu, ple, bcnv, h2o2 )

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: i1, i2, j1, j2, km, n1, n2, dt_30m
  character(len=*)     :: aero_type
  REAL*8,    INTENT(INOUT) :: tc(i1:i2,j1:j2,km,n1:n2)
!  REAL*8,    INTENT(INOUT) :: cldscv(i1:i2,j1:j2,km,n1:n2), cldso2(i1:i2,j1:j2,km)
!  REAL*8,    INTENT(INOUT) :: cldso4(i1:i2,j1:j2,km), cldmsa(i1:i2,j1:j2,km)
  REAL*8                    :: cldscv(i1:i2,j1:j2,km,n1:n2), cldso2(i1:i2,j1:j2,km)
  REAL*8                    :: cldso4(i1:i2,j1:j2,km), cldmsa(i1:i2,j1:j2,km)
!  REAL*8,    INTENT(INOUT) :: tcnv(i1:i2,j1:j2,n1:n2)
!  REAL*8,    INTENT(INOUT) :: wet_conv_in(i1:i2,j1:j2,km,n1:n2)
  REAL*8,    INTENT(INOUT) :: airmass(i1:i2,j1:j2,km)
  REAL*8,    DIMENSION(i1:i2,j1:j2,km),   INTENT(IN) :: vud
  REAL*8,    DIMENSION(i1:i2,j1:j2,km+1), INTENT(IN) :: cldmas
  REAL*8,    DIMENSION(i1:i2,j1:j2,km),   INTENT(IN) :: dtrain
  REAL*8,    DIMENSION(i1:i2,j1:j2),      INTENT(IN) :: area
  REAL*8,    DIMENSION(i1:i2,j1:j2,km),   INTENT(IN) :: delz, delp
  REAL*8,    DIMENSION(i1:i2,j1:j2,km),   INTENT(IN) :: airmol
  REAL*8,    DIMENSION(i1:i2,j1:j2,km),   INTENT(IN) :: tmpu
  REAL*8,    DIMENSION(i1:i2,j1:j2,km+1), INTENT(IN) :: ple
  REAL*8,    INTENT(OUT)   :: bcnv(i1:i2,j1:j2,n1:n2)
  REAL*8,    INTENT(INOUT), optional :: h2o2(i1:i2,j1:j2,km)
  LOGICAL, INTENT(INOUT)   :: KIN  ! true for aerosol

  REAL*8    :: tc1(i1:i2,j1:j2,km,n1:n2), f(i1:i2,j1:j2,km,n1:n2)
  REAL*8    :: cldmas_tmp(i1:i2,j1:j2,km), so2loss
! epsilon: A very small positive number   [unitless]
  REAL*8,  PARAMETER   :: EPSILON = 1.0E-32
  REAL*8,  PARAMETER   :: R = 8.2057d-2  ! universal gas constant [L*atm/moles/K]
  REAL*8,  PARAMETER   :: INV_T0 = 1d0 / 298d0
  REAL*8,  PARAMETER   :: conv_NH3 = 5.69209978831d-1 ! 0.6*SQRT(0.9) for ice to gas ratio
  REAL*8  :: kg, Kstar298, H298_R, I2G, L2G, C_TOT, F_L, F_I
  INTEGER :: n, i, j, l, NSO2, NSO4, NMSA
  REAL*8,    DIMENSION(i1:i2,j1:j2,km) :: c_h2o
  REAL*8,    DIMENSION(i1:i2,j1:j2,km) :: cldliq
  REAL*8,    DIMENSION(i1:i2,j1:j2,km) :: cldice

!  Initialize local variables
!  --------------------------
!  c_h2o, cldliq, and cldice are respectively intended to be the 
!  water mixing ratio (liquid or vapor?, in or out of cloud?)
!  cloud liquid water mixing ratio
!  cloud ice water mixing ratio
   c_h2o  = (10d0**(-2663.5d0/tmpu(:,:,:) + 12.537d0 ) ) /  &
                   (ple(:,:,1:km)+ple(:,:,2:km+1)) /2d0   
   cldliq = 0.d0
   where(tmpu >  248.) cldliq = 1.d-6 * ( ( tmpu - 248.d0) / 20.d0 )
   where(tmpu >= 268.) cldliq = 1.d-6
   cldice = 1.d-6 - cldliq

  ! executable statements

  tc1(:,:,:,:) = tc(:,:,:,:)
  !if (MAPL_AM_I_ROOT()) print *, 'hbian convection tmpu =', tmpu(i1,j1,1), tmpu(i1,j1,km)
  
! compute the fraction of tracer scavenged in convective cloud updrafts
  f = 0.0
  kg = 0d0

  DO n = n1, n2
     if (TRIM(aero_type) .eq. 'nitrate' .and. n .eq. n1 )  kin = .false.  ! treat NH3 as a gas tracer
     if (TRIM(aero_type) .eq. 'nitrate' .and. n .gt. n1 )  kin = .true.   ! treat others as aerosol

     if (kin) then
        CALL f_aerosol(i1, i2, j1, j2, km, kc, f(:,:,:,n), delz, vud )
     else
        ! gas tracer NH3
        if (TRIM(aero_type) .eq. 'nitrate' .and. n .eq. n1 )  then
        ! values adopted in Umich/IMPACT and GMI, effective Herry's law coefficient at pH = 5
           Kstar298 = 1.05d6
           H298_R = -4.2d3
        endif
           DO L = 2, KM
           DO J = j1, j2
           DO I = i1, i2
              ! ice to gas ratio
              if ( c_h2o(i,j,l) > 0.d0) then
                 I2G = (cldice(i,j,l) / c_h2o(i,j,l)) * conv_NH3
              else
                 I2G = 0.d0
              endif
              L2G = cldliq(i,j,l) * R * tmpu(i,j,l) * &
                      Kstar298 * EXP( -H298_R * ( ( 1d0 / tmpu(i,j,l) ) - INV_T0 ) )
              ! fraction of NH3 in liquid & ice phases
              C_TOT = 1d0 + L2G + I2G
              F_L = L2G / C_TOT
              F_I = I2G / C_TOT
              ! compute kg, the retention factor for liquid NH3 is 0 at T < 248K and
              ! 0.05 at 248K < T < 268K
              if (tmpu(i,j,l) >=268d0) then
                 kg = kc * ( F_L+F_I )
              elseif ( (248d0 < tmpu(i,j,l)) .and. (tmpu(i,j,l) < 268d0) ) then
                 kg = kc * ( (0.05*F_L)+F_I )
              else
                 kg = kc * F_I
              endif
              if(kg > 0.d0 .and. vud(i,j,l) > 1.e-14) &
               f(i,j,l,n) = 1.0 - EXP( -kg * delz(i,j,l) / vud(i,j,l) )
           ENDDO
           ENDDO
           ENDDO
     endif

!    Special treatment for DMS and SO2 if aero_type is "sulfur"
     if(trim(aero_type) .eq. 'sulfur') then

        if(.not.present(h2o2)) call die ('GOCARTConvectionMod.F90', &
                                         'missing required H2O2 for sulfur') 

        if(n .eq. n1)    f(:,:,:,n1) = 0.0    ! DMS
        !if(n .eq. n1+1)  f(:,:,:,n1+1) = 0.0  ! SO2 for now is not scavenged   
#undef PRC
!#ifdef PRC
        if(n .eq. n1+1) then                  ! SO2 requires special handling

        !==============================================================
        ! Coupled full chemistry/aerosol simulation:
        ! Use the wet scavenging formula of Chin et al [1996], 
        ! such that a soluble fraction of SO2 is limited by the
        ! availability of H2O2 in the precipitating grid box. 
        ! Scavenge the soluble SO2 at the same rate as the sulfate.
        ! Update H2O2_sav and SO2_sav for use in RAINOUT, WASHOUT
        !==============================================================
        DO L = 2, KM
           DO J = j1, j2
              DO I = i1, i2

                 ! Make sure to deplete H2O2s the same as SO2s.
                 ! (dkh, rjp, bmy, 11/17/05)
                 ! based on GEOS-Chem. tq, 01/09
                 
                 IF ( tc1(i,j,l,n) > epsilon ) THEN
                    
                    ! limit f
                    so2loss  = MIN( h2o2(i,j,l), tc1(i,j,l,n) )
                    f(i,j,l,n) = f(i,j,l,n) * so2loss / tc1(i,j,l,n)
                    f(i,j,l,n) = MAX(f(i,j,l,n), 0.0)
                    
                    ! update saved h2o2 concentration
                    h2o2(i,j,l) = h2o2(i,j,l) - ( tc1(i,j,l,n) * f(i,j,l,n) )
                    h2o2(i,j,l) = MAX( h2o2(i,j,l), epsilon )
                    
                 ELSE
                    
                    ! set f = 0 if so2 < epsilon (dkh, rjp, bmy, 11/17/05)
                    f(i,j,l,n) = 0.d0
                    
                 END IF
                 
              ENDDO
           ENDDO
        ENDDO
        endif                            ! SO2
!#endif
     endif                               ! sulfur

  ENDDO  ! n

! if tracer is type "carbon" then set coefficient to 0 for hydrophobic
! implementing QQ Wang's change by Huisheng Bian (4/24/2015)
! not scavenging  BCn1 (hydrophobic) when T > 258 K
  if(trim(aero_type) .eq. 'OC') f(:,:,:,n1) = 0d0

! suppress scavenging most aerosols at cold T except BCn1 (hydrophobic), dust, and HNO3
  if (trim(aero_type) .eq. 'BC') then
     where (tmpu >= 258.d0)
        f(:,:,:,n1) = 0.d0
     end where
  end if

  if (trim(aero_type) .eq. 'BC') then
     where (tmpu < 258.d0)
        f(:,:,:,n2) = 0.d0
     end where
  end if

  if (trim(aero_type) .eq. 'OC'       .or. &
      trim(aero_type) .eq. 'sea_salt' .or. &
      trim(aero_type) .eq. 'sulfur'   .or. &
      trim(aero_type) .eq. 'seasalt'  .or. &
      trim(aero_type) .eq. 'sulfate'  .or. &
      trim(aero_type) .eq. 'nitrate'  .or. &
      trim(aero_type) .eq. 'NH3'      .or. &
      trim(aero_type) .eq. 'NH4a') then

      do n = n1, n2
         where (tmpu < 258.d0 )
            f(:,:,:,n) = 0.d0
         endwhere
      end do

  end if


  ! re-index for routine cldcnv
  cldmas_tmp(:,:,1:km) = cldmas(:,:,2:km+1)

  ! internal time step for the convection routine is 300s
  CALL cldcnv(i1, i2, j1, j2, km, n1, n2, dt_30m, aero_type, &
              tc, f, airmass, area, cldmas_tmp, dtrain, delz, delp)

  ! -- Mass balance
  SELECT CASE(TRIM(aero_type))

  CASE('sulfur')
     NSO2 = n1+1
     NSO4 = n1+2
     NMSA = n1+3

     cldso2 = 0.0d0
     cldso4 = 0.0d0
     cldmsa = 0.0d0

     DO l = 1,km
        DO j = j1, j2
           DO i = i1, i2
              cldso2(i,j,l) = cldso2(i,j,l) + &
                  (tc(i,j,l,NSO2) - tc1(i,j,l,NSO2)) * airmass(i,j,l)
              cldso4(i,j,l) = cldso4(i,j,l) + &
                  (tc(i,j,l,NSO4) - tc1(i,j,l,NSO4)) * airmass(i,j,l)
              cldmsa(i,j,l) = cldmsa(i,j,l) + &
                  (tc(i,j,l,NMSA) - tc1(i,j,l,NMSA)) * airmass(i,j,l)
           END DO
        END DO
     END DO
     
     DO n = n1, n2
        DO i = i1, i2
           DO j = j1, j2
              bcnv(i,j,n) = 0.0
              DO l = 1,km
                 IF (tc(i,j,l,n) < 0.0) tc(i,j,l,n) = 1.0E-32
                 ! kg tracer
                 bcnv(i,j,n) = bcnv(i,j,n) + &
                     (tc(i,j,l,n) - tc1(i,j,l,n)) *airmass(i,j,l)
              END DO
           END DO
        END DO
     END DO

  CASE('co')
     cldscv = 0.0d0
     DO n = n1, n2
        DO l = 1,km
           DO j = j1, j2
              DO i = i1, i2
                 cldscv(i,j,l,n) = cldscv(i,j,l,n) + &
 !                wet_conv_in(i,j,l,n) = wet_conv_in(i,j,l,n) + &
                      (tc(i,j,l,n) - tc1(i,j,l,n)) * airmol(i,j,l)
              END DO
           END DO
        END DO
     END DO

     bcnv(:,:,:) = 0.0
     DO n = n1, n2
        DO i = i1, i2
           DO j = j1, j2
              DO l = 1,km
                 IF (tc(i,j,l,n) < 0.0) tc(i,j,l,n) = 1.0E-32
                 bcnv(i,j,n) = bcnv(i,j,n) + &
                      (tc(i,j,l,n) - tc1(i,j,l,n)) * airmol(i,j,l)
              END DO
           END DO
        END DO
     END DO

  CASE DEFAULT

!!$     DO n = n1, n2
!!$        DO l = 1,km
!!$           DO j = j1, j2
!!$              DO i = i1, i2
!!$                 wet_conv_in(i,j,l,n) = wet_conv_in(i,j,l,n) + &
!!$                      (tc(i,j,l,n) - tc1(i,j,l,n)) * airmass(i,j,l)
!!$              END DO
!!$           END DO
!!$        END DO
!!$     END DO

     DO n = n1, n2
        DO i = i1, i2
           DO j = j1, j2
              bcnv(i,j,n) = 0.0
              DO l = 1,km
                 IF (tc(i,j,l,n) < 0.0) tc(i,j,l,n) = 1.0E-32
                 bcnv(i,j,n) = bcnv(i,j,n) + &
                      (tc(i,j,l,n) - tc1(i,j,l,n)) * airmass(i,j,l)
              END DO
           END DO
        END DO
     END DO

  END SELECT

!  tcnv(:,:,:) = tcnv(:,:,:) + bcnv(:,:,:)

END SUBROUTINE convection
