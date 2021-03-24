subroutine largescalemoist_driver(kmax,DT,TH,Q,P,PRC3_LSC)

IMPLICIT NONE

!MAPL CONSTANTS
real, parameter :: GRAV   = 9.80                   ! m^2/s
real, parameter :: ALHL   = 2.4665E6               ! J/kg @15C
real, parameter :: AIRMW  = 28.97                  ! kg/Kmole
real, parameter :: RUNIV  = 8314.3                 ! J/(Kmole K)
real, parameter :: KAPPA  = 2.0/7.0                ! --
real, parameter :: RGAS   = RUNIV/AIRMW  ! J/(kg K)
real, parameter :: CP     = RGAS/KAPPA   ! J/(kg K)

!!INPUTS!!
INTEGER :: kmax !Vertical levels
REAL :: DT !Moist time step.

REAL, DIMENSION(kmax+1) :: PREF !Reference vertical pressure.

!!INPUT/OUTPUT!!
REAL, DIMENSION(kmax) :: TH, Q !Input vertical Profiles of Potential tempearture and Humidity
REAL, DIMENSION(kmax+1) :: P   !Input of edge pressure in hPa

!!OUTPUTS!!
REAL, DIMENSION(kmax) :: PRC3_LSC !Diagnosed column precipitation.

!!LOCALS!!
INTEGER :: k !Loop index.

REAL, DIMENSION(kmax) :: TH_lsc, Q_lsc !Temporary storing of variables
REAL, DIMENSION(kmax+1) :: P_lsc

REAL, DIMENSION(kmax) :: T_lsc, Ph_lsc, Pih_lsc, cpm, Qs_lsc, DqsDT_lsc, deltaq, deltaTH, col_rain

   TH_lsc = TH
   q_lsc = Q
   P_lsc = P

   col_rain = 0.0

   !Compute pressure and Exner pressure pi at mid point
   Ph_lsc(1:kmax) = 0.5*( P_lsc(1:kmax) +  P_lsc(2:kmax+1) ) 
   Pih_lsc = (Ph_lsc/1000.0)**KAPPA

   !Compute moist cp
   cpm = CP*(1+0.887*q_lsc)

   !Compute temperature
   T_lsc = TH_lsc*Pih_lsc

   !Compute sat vap pressure (qs) and its gradient (dqsdT) from T and Ph.
   call MRQSAT1(T_lsc,Ph_lsc,qs_lsc,DqsdT_lsc,kmax)

   deltaq = 0.0
   DO k = 40,kmax
      if (q_lsc(k) - qs_lsc(k) .gt. 0) THEN !If supersaturated

         !Compute excess moisture per grid cell due to supersaturation.
         deltaq(k) = (q_lsc(k)-qs_lsc(k))/(1+dqsdT_lsc(k)*ALHL/cpm(k))

         !Remove excess moisture, i.e. rain out.
         q_lsc(k) = q_lsc(k) - deltaq(k)

         !Compute temperture change caused by latent heating.
         deltaTH(k) = (ALHL*deltaq(k))/(Pih_lsc(k)*cpm(k))
   
         !Update temeprature.
         TH_lsc(k) = TH_lsc(k) + deltaTH(k)

         !Compute column precipitation rate in mm/s.
         col_rain(k) = (1/DT) * 100 * deltaq(k) * (P_lsc(k+1)-P_lsc(k))/GRAV
      
      endif
   endDo

   !Return updated TH and q as well as the column precip.
   TH = TH_lsc
   Q = Q_lsc
   PRC3_lsc = col_rain

end subroutine largescalemoist_driver


SUBROUTINE MRQSAT1(TT,P,Qs,DQsDT,kmax)

   implicit none

   !Inputs
   integer :: kmax
   real, dimension(kmax) :: TT, P, Qs, DQsDT
   logical :: LDQDT

   !Constants
   real, parameter :: svpt0=273.15, svp1=.6112, svp2=17.67, svp3=29.65, airmw = 28.97, h2omw = 18.01, one = 1.0, ep2 = 0.622
   real, parameter :: esfac = h2omw/airmw, erfac = (one-esfac)/esfac

   !Internals
   real :: dqdta, qsatvp, qsatvp_exc, d
   integer :: kk 

   do kk=40,kmax

      qsatvp_exc = svp2*(TT(kk) - svpt0)/(TT(kk) - svp3)
      qsatvp_exc = min(qsatvp_exc, 10.0)

      qsatvp = svp1*exp(qsatvp_exc)

      d = 1./(.1*p(kk)-ep2*erfac*qsatvp)
      Qs(kk)=ep2*qsatvp*d

      dqdta =  (svp2*(svpt0-svp3)/(svp3-TT(kk))**2 )*qsatvp
      dQsdt(kk) = ep2*d*dqdta*(1. + erfac*Qs(kk))

   enddo

end subroutine MRQSAT1


