subroutine dummy ( im, jm, km, nbins, flag, &
                                radiusInp, rhopInp, cdt, &
                                delp, RH, AERO, &
                                tmpu, rhoa, &
                                dz, &
                                correctionMaring, grav )

  !INPUTS PARAMETERS:
  integer, intent(in) :: im,jm,km, nbins
  integer, intent(in) :: flag 
  real(8), intent(in) :: cdt 
  real(8), intent(in), dimension(nbins)    :: radiusInp, rhopInp
  real(8), intent(in), dimension(im,jm,km) :: RH, delp
  real(8), intent(in), dimension(im,jm,km) :: tmpu, rhoa, dz
  real(8) :: grav
  logical, intent(in)    :: correctionMaring
  

! !OUTPUT PARAMETERS:
  real(8), intent(inout), dimension(im,jm,km,nbins) :: AERO


RH = RH**2*RHOA
RHOA = RHOA**2*RH
tmpu = tmpu**2*RH*RHOA

call Chem_Settling_Pert ( im, jm, km, nbins, flag, &
                          radiusInp, rhopInp, cdt, &
                               delp, RH, AERO, &
                                tmpu, rhoa, &
                                dz, &
                                correctionMaring, grav )

AERO = AERO**2
RH = RH**2*RHOA
RHOA = RHOA**2*RH
tmpu = tmpu**2*RH*RHOA


endsubroutine dummy




subroutine Chem_Settling_Pert ( im, jm, km, nbins, flag, &
                                radiusInp, rhopInp, cdt, &
                                delp, RH, AERO, &
                                tmpu, rhoa, &
                                dz, &
                                correctionMaring, grav )

  implicit NONE

  !radiusInp - Particle radius
  !rhopInp   - Soil density
  !cdt       - Time step
  !delp      - Level thickness (pressure)
  !RH        - Relative humidity     (potential pert)
  !AERO      - Aerosol                         (pert)
  !tmpu      - Temperature                     (pert)
  !rhoa      - Airdens
  !dz        - Level thickness (m)


  !INPUTS PARAMETERS:
  integer, intent(in) :: im,jm,km, nbins
  integer, intent(in) :: flag 
  real(8), intent(in) :: cdt 
  real(8), intent(in), dimension(nbins)    :: radiusInp, rhopInp
  real(8), intent(in), dimension(im,jm,km) :: RH, delp
  real(8), intent(in), dimension(im,jm,km) :: tmpu, rhoa, dz
  real(8) :: grav
  logical, intent(in)    :: correctionMaring
  

! !OUTPUT PARAMETERS:
  real(8), intent(inout), dimension(im,jm,km,nbins) :: AERO


! !Local Variables
   integer  ::  i, j, k, iit, n
   real(8), parameter ::  rhow = 1000.  ! Density of water [kg m-3]
   real(8) :: pdog(im,jm)               ! air mass factor dp/g [kg m-2]
   real(8) :: pdog_m1(im,jm)            ! air mass factor dp/g [kg m-2]
   real(8) :: vsettle(im,jm,km)         ! fall speed [m s-1]
   real(8) :: q_save(im,jm)             ! save the dust mmr [kg kg-1]
   real(8) :: q_before(im,jm)           ! save the dust mmr [kg kg-1]
   real(8) :: diff_coef                 ! Brownian diffusion coefficient [m2 s-1]

   !The following parameters relate to the swelling of seasalt like particles
   !following Fitzgerald, Journal of Applied Meteorology, 1975.
   real(8), parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
   real(8), parameter :: alphaNaCl = 1.35
   real(8) :: alpha, alpha1, alpharat, beta, theta, f1, f2

   !Parameter from Gerber 1985 (units require radius in cm, see rcm)
   real(8) :: rcm
   real(8), parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424
   !Parameters for ammonium sulfate
   real(8), parameter :: SU_c1=0.4809, SU_c2=3.082, SU_c3=3.110e-11, SU_c4=-1.428

   !Parameters from Maring et al, 2003
   real(8), parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real(8), parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]

   real(8) :: sat, rrat
   real(8) :: radius, rhop   ! particle radius and density passed to
   real(8) :: dt_settle, minTime, qmin, qmax
   integer :: nSubSteps, ijl


   ijl = im*jm

!  Loop over the number of dust bins
   do n = 1, nbins

    radius = radiusInp(n)
    rhop = rhopInp(n)

!   Reset a (large) minimum time to cross a grid cell in settling
    minTime = cdt


!   If radius le 0 then get out of loop
    if(radius .le. 0.) cycle

    do k = 1, km
     do j = 1, jm
      do i = 1, im

!      Adjust the particle size for relative humidity effects
       sat = max(rh(i,j,k),1.0) ! to avoid zero FPE

!      Fitzgerald
       if(flag .eq. 1 .and. sat .ge. 0.80) then
!       parameterization blows up for RH > 0.995, so set that as max
!       rh needs to be scaled 0 - 1
        sat = min(0.995,sat)
!       Calculate the alpha and beta parameters for the wet particle
!       relative to amonium sulfate
        beta = exp( (0.00077*sat) / (1.009-sat) )
        if(sat .le. 0.97) then
         theta = 1.058
        else
         theta = 1.058 - (0.0155*(sat-0.97)) /(1.02-sat**1.4)
        endif
        alpha1 = 1.2*exp( (0.066*sat) / (theta-sat) )
        f1 = 10.2 - 23.7*sat + 14.5*sat**2.
        f2 = -6.7 + 15.5*sat - 9.2*sat**2.
        alpharat = 1. - f1*(1.-epsilon) - f2*(1.-epsilon**2.)
        alpha = alphaNaCl * (alpha1*alpharat)
!       radius is the radius of the wet particle
        radius = alpha * radiusInp(n)**beta
        rrat = (radiusInp(n)/radius)**3.
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       elseif(flag .eq. 2) then   ! Gerber
        sat = min(0.995,sat)
        rcm = radiusInp(n)*100.
        radius = 0.01 * (   c1*rcm**c2 / (c3*rcm**c4-log(sat)) &
                          + rcm**3.)**(1./3.)
        rrat = (radiusInp(n)/radius)**3.
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       elseif(flag .eq. 3) then   
!       Gerber parameterization for Ammonium Sulfate
        sat = min(0.995,sat)
        rcm = radiusInp(n)*100.
        radius = 0.01 * (   SU_c1*rcm**SU_c2 / (SU_c3*rcm**SU_c4-log(sat)) &
                      + rcm**3.)**(1./3.)
        rrat = (radiusInp(n)/radius)**3.
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       elseif(flag .eq. 4) then
!       Petters and Kreidenweis (ACP2007) parameterization
        sat = min(0.99,sat)
        radius = (radiusInp(n)**3 * (1+1.19*sat/(1-sat)))**(1./3.)
        rrat = (radiusInp(n)/radius)**3
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       endif

!      Calculate the settling velocity
       call Chem_CalcVsettle_Pert( radius, rhop, rhoa(i,j,k), &
                                   tmpu(i,j,k), diff_coef, vsettle(i,j,k), grav)
      end do
     end do
    end do

    if ((correctionMaring) .and. (radiusInp(n) .le. (0.5*diameterMaring))) then
        vsettle = max(1.0e-9, vsettle - v_upwardMaring)
    endif

!   Determine global max/min time to cross grid cell
    call pmaxmin ( dz/vsettle, qmin, qmax, im, jm, km)
    minTime = min(minTime,qmin)


!   Now, how many iterations do we need to do?
    if ( minTime < 0 ) then
         nSubSteps = 0
    else if(minTime .ge. cdt) then
     nSubSteps = 1
     dt_settle = cdt
    else
     nSubSteps = int(cdt/minTime+1)
     dt_settle = cdt/nSubSteps
    endif

!   Loop over sub-timestep
    do iit = 1, nSubSteps

     q_save = AERO(:,:,1,n)
     AERO(:,:,1,n) = AERO(:,:,1,n) / (1.+dt_settle*vsettle(:,:,1)/dz(:,:,1))

     do k = 2, km

      !Air mass factors of layers k, k-1
      pdog = delp(:,:,k)/grav
      pdog_m1 = delp(:,:,k-1)/grav
      q_before = AERO(:,:,k,n)
      AERO(:,:,k,n) = 1./(1.+dt_settle*vsettle(:,:,k)/dz(:,:,k)) &
                    * ( AERO(:,:,k,n) &
                    + ( dt_settle*vsettle(:,:,k-1)/dz(:,:,k-1) &
                    * q_save*pdog_m1(:,:)/pdog(:,:) ) )
      q_save = q_before
     end do ! k

    end do  ! iit


   end do   ! n

endsubroutine Chem_Settling_Pert 


subroutine Chem_CalcVsettle_Pert ( radius, rhop, rhoa, tmpu, diff_coef, vsettle, grav )

  implicit NONE

  !INPUTS
  real(8), intent(in)    :: radius              ! Particle radius [m]
  real(8), intent(in)    :: rhop                ! Particle density [kg m-3]
  real(8), intent(in)    :: rhoa                ! Layer air density [kg m-3]
  real(8), intent(in)    :: tmpu                ! Layer temperature [K]
  real(8), intent(in)    :: grav                ! Layer temperature [K]

  !OUTPUTS
  real(8), intent(out)   :: diff_coef               ! Brownian diffusion coefficient [m2 s-1]
  real(8), intent(out)   :: vsettle                 ! Layer fall speed [m s-1]

  !Locals
  real(8) :: rmu                         ! Dynamic viscosity [kg m-1 s-1]
  real(8) :: vt                          ! Thermal velocity of air molecule [m s-1]
  real(8) :: rmfp                        ! Air molecule mean free path [m]
  real(8) :: bpm                         ! Cunningham slip correction factor
  real(8) :: rkn                         ! Knudsen number
  real(8) :: re, x, y                    ! reynolds number and parameters

  real(8), parameter :: kb = 1.3807e-23    ! Boltzmann constant [kg m2 s-1 K-1 mol-1]
  real(8), parameter :: m_air = 4.8096e-26 ! Mass of <avg> air molecule [kg]
  real(8), parameter :: pi = 3.141529265


!  Dynamic viscosity from corrected Sutherlands Equation
   rmu = 1.8325e-5*(416.16/(tmpu+120.))*(tmpu/296.16)**1.5

!  Thermal velocity of air molecule
   vt = sqrt(8.*kb*tmpu/pi/m_air)

!  Air molecule mean free path
   rmfp = 2.*rmu/rhoa/vt

!  Knudsen number
   rkn = rmfp/radius

!  Cunningham slip correction factor
   bpm = 1. + 1.246*rkn + 0.42*rkn*exp(-0.87/rkn)

!  Brownian diffusion coefficient
   diff_coef = kb*tmpu*bpm/3./pi/rmu/(2.*radius)

!  Fall speed (assumes Reynolds # < 0.01)
   vsettle = 2./9.*rhop*radius**2.*grav*bpm/rmu

!  Check the Reynolds number to see if we need a drag correction
!  First guess at Reynolds number using Stokes calculation
   re = 2.*rhoa*radius*vsettle/rmu

!  If Re > 0.01 then apply drag correction following Pruppacher and
!  Klett regime 2 (eq. 10-142).  Assuming reasonable aerosols we
!  do not consider that particle Re may exceed 300.
   if(re .gt. 0.01) then
    x = log(24.*re/bpm)
    y = -3.18657 + 0.992696   *x     - .00153193   *x**2. &
                 - 0.000987059*x**3. - .000578878  *x**4. &
                 + 8.55176E-05*x**5. -  3.27815E-06*x**6.
    re = exp(y)*bpm
    vsettle = rmu*re/2./rhoa/radius
   endif

   end subroutine Chem_CalcVsettle_Pert


subroutine pmaxmin(a, qmin, qmax, im, jm, km)

implicit none

integer, intent(in) :: im, jm, km
real(8), intent(in), dimension(im,jm,km) :: a
real(8), intent(out) :: qmin, qmax


qmin = min(abs(a))
qmax = max(abs(a))



endsubroutine pmaxmin
