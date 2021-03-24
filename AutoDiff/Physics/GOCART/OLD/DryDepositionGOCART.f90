subroutine dummy ( im, jm, km, &
                                      tmpu, rhoa, dz1, oro, ustar, &
                                      pblh, shflux, z0h, drydepf, &
                                      radius, rhop, u10m, v10m, fraclake, gwettop, &
                                      dustflag, von_karman, cpd, grav )

  implicit none

  !INPUTS
  integer, intent(in)                      :: im, jm, km            ! grid information
  integer, intent(in)                      :: dustflag              ! flag for if doing dust
  real(8), intent(in)                      :: von_karman, cpd, grav !constants
  real(8), intent(in), dimension(im,jm,km) :: tmpu                  ! temperature [K]                 (pert)
  real(8), intent(in), dimension(im,jm,km) :: rhoa                  ! air density [kg m-3]
  real(8), intent(in), dimension(im,jm)    :: dz1                   ! lowest layer thickness [m]
  real(8), intent(in), dimension(im,jm)    :: oro                   ! orography flag
  real(8), intent(in), dimension(im,jm)    :: ustar                 ! friction speed                  (should be)
  real(8), intent(in), dimension(im,jm)    :: pblh                  ! PBL height [m]
  real(8), intent(in), dimension(im,jm)    :: shflux                ! sfc. sens. heat flux [W m-2]    (should be)
  real(8), intent(in), dimension(im,jm)    :: z0h                   ! rough height, sens. heat [m]
  real(8), intent(in)                      :: radius                ! particle radius [m]
  real(8), intent(in)                      :: rhop                  ! particle density [kg m-3]
  real(8), intent(in), dimension(im,jm)    :: u10m                  ! 10-m u-wind component [m s-1]   (pert)
  real(8), intent(in), dimension(im,jm)    :: v10m                  ! 10-m v-wind component [m s-1]   (pert)
  real(8), intent(in), dimension(im,jm)    :: fraclake              ! fraction covered by water
  real(8), intent(in), dimension(im,jm)    :: gwettop               ! fraction soil moisture

  !OUTPUTS
  real(8), intent(out), dimension(im,jm)   :: drydepf      ! Deposition frequency [s-1]


tmpu = tmpu**2*rhoa
rhoa = tmpu*rhoa**2
u10m = u10m**2*v10m
v10m = u10m*v10m**2

call DryDepositionGOCART_Pert ( im, jm, km, &
                                      tmpu, rhoa, dz1, oro, ustar, &
                                      pblh, shflux, z0h, drydepf, &
                                      radius, rhop, u10m, v10m, fraclake, gwettop, &
                                      dustflag, von_karman, cpd, grav )

drydepf = drydepf**2
tmpu = tmpu**2*rhoa
rhoa = tmpu*rhoa**2
u10m = u10m**2*v10m
v10m = u10m*v10m**2


endsubroutine dummy



subroutine DryDepositionGOCART_Pert ( im, jm, km, &
                                      tmpu, rhoa, dz1, oro, ustar, &
                                      pblh, shflux, z0h, drydepf, &
                                      radius, rhop, u10m, v10m, fraclake, gwettop, &
                                      dustflag, von_karman, cpd, grav )

  implicit none

  !INPUTS
  integer, intent(in)                      :: im, jm, km            ! grid information
  integer, intent(in)                      :: dustflag              ! flag for if doing dust
  real(8), intent(in)                      :: von_karman, cpd, grav !constants
  real(8), intent(in), dimension(im,jm,km) :: tmpu                  ! temperature [K]                 (pert)
  real(8), intent(in), dimension(im,jm,km) :: rhoa                  ! air density [kg m-3]
  real(8), intent(in), dimension(im,jm)    :: dz1                   ! lowest layer thickness [m]
  real(8), intent(in), dimension(im,jm)    :: oro                   ! orography flag
  real(8), intent(in), dimension(im,jm)    :: ustar                 ! friction speed                  (should be)
  real(8), intent(in), dimension(im,jm)    :: pblh                  ! PBL height [m]
  real(8), intent(in), dimension(im,jm)    :: shflux                ! sfc. sens. heat flux [W m-2]    (should be)
  real(8), intent(in), dimension(im,jm)    :: z0h                   ! rough height, sens. heat [m]
  real(8), intent(in)                      :: radius                ! particle radius [m]
  real(8), intent(in)                      :: rhop                  ! particle density [kg m-3]
  real(8), intent(in), dimension(im,jm)    :: u10m                  ! 10-m u-wind component [m s-1]   (pert)
  real(8), intent(in), dimension(im,jm)    :: v10m                  ! 10-m v-wind component [m s-1]   (pert)
  real(8), intent(in), dimension(im,jm)    :: fraclake              ! fraction covered by water
  real(8), intent(in), dimension(im,jm)    :: gwettop               ! fraction soil moisture

  !OUTPUTS
  real(8), intent(out), dimension(im,jm)   :: drydepf      ! Deposition frequency [s-1]

  !Locals
  integer :: i, j, k, n
  real(8), parameter :: rhow = 1000.      ! density of water [kg m-3]
  real(8), parameter :: coll_size = 0.002 ! collector size [m]
  real(8) :: rmu(im,jm)                   ! dynamic viscosity [kg m-1 s-1]
  real(8) :: Ra(im,jm)                    ! aerodynamic resistance
  real(8) :: Rs(im,jm)                    ! surface resistance
  real(8) :: vdep(im,jm)                  ! Deposition speed [m s-1]
  real(8) :: obk(im,jm)                   ! Obukhov Length [m]

  real(8) :: qmin, qmax

  real(8) :: Rttl                         ! total surface resistance
  real(8) :: dc
  real(8) :: diff_coef, vsettle
  real(8) :: Sc                           ! Schmidt number
  real(8) :: Eb                           ! Brownian diffusion collection efficiency
  real(8) :: St                           ! Stokes number
  real(8) :: Ein                          ! Interception collection efficiency
  real(8) :: Eim                          ! Impaction collection efficiency
  real(8) :: alpha, gamma

  real(8) :: R1, R2, w10m, u_thresh0
  real(8) :: vds, vdsmax, czh, factor
  real(8) :: frac, cz, psi_h, eps, logmfrac, z0h_min, z0h_, r8_cdt
  real(8) :: one = 1.0, zero = 0.0

  real(8), parameter :: OCEAN=0.0

  !Initialize output to zeros
  drydepf = 0.0

  !Calculate the viscosity and thickness of the surface level
  rmu = 1.8325e-5*(416.16/(tmpu(:,:,km)+120.)) *(tmpu(:,:,km)/296.16)**1.5

  z0h_min = 100.0   ! because sometimes we may get z0h=0.

  !Calculate the Obukhov length scale
  call ObukhovLength_Pert ( im, jm, tmpu(:,:,km), rhoa(:,:,km), shflux, ustar, obk, von_karman, cpd, grav )

  !Aerodynamic Resistance
  do j = 1, jm
     do i = 1, im

        cz = dz1(i,j) / 2.
        frac = cz / obk(i,j)

        if (frac .gt. 1.) then
           frac = 1.
        endif

        if (frac .gt. 0. .and. frac .le. 1.) then
           psi_h = -5.0*frac
        else if (frac .lt. 0.) then
           eps = min(one,-frac)
           logmfrac = log(eps)
           psi_h = exp(0.598 + 0.39*logmfrac - 0.09*(logmfrac)**2.)
        endif

        z0h_ = max ( z0h(i,j), z0h_min )

        Ra(i,j) = (log(cz/z0h_) - psi_h) / (von_karman*ustar(i,j))

     enddo
  enddo

  !Surface Resistance term for aerosols
  do j = 1, jm
     do i = 1, im

        !Calculate the surface resistance term
        vds = 0.002*ustar(i,j)

        !Set to small value of vds if ustar too small
        vds = max(vds, 0.002 * 0.00001)
 
        if (obk(i,j) .lt. 0.) then
           vds = vds*(1.+(-300./obk(i,j))**0.6667)
        endif

        czh = pblh(i,j)/obk(i,j)

        if (czh .lt. -30.) then
           vds = 0.0009*ustar(i,j)*(-czh)**0.6667
        endif

        vdsMax = 0.01

        Rs(i,j) = 1./min(vds,vdsmax)

        if (Rs(i,j) .gt. 9999.) then
           Rs(i,j) = 9999.
        endif

        if (Rs(i,j) .lt. 1.) then
           Rs(i,j) = 1.
        endif

        !If doing dust over land, possibly re-emit. Logic is to check on optional provided parameter and modify R2
        R2 = 1.
        if ( dustflag == 1 ) then

           !Calculate the threshold velocity for dust emissions
           u_thresh0 = 0.13 * sqrt(rhop*grav*2.*radius/rhoa(i,j,km)) &
                            * sqrt(1.+6.e-7/(rhop*grav*(2.*radius)**2.5)) &
                            / sqrt(1.928*(1331.*(100.*2.*radius)**1.56+0.38)**0.092 - 1.)
           w10m = sqrt(u10m(i,j)**2. + v10m(i,j)**2.)

           !Calculate the coefficient for resuspension
           if (abs(oro(i,j)) .le. 1e-4) then
              R2 = 1.
           else
              R2 = fraclake(i,j)+(1.-fraclake(i,j))*( gwettop(i,j)+(1.-gwettop(i,j))*exp(-max(zero,(w10m-u_thresh0))))
           endif

        endif

        !Now what is the deposition velocity
        Rttl = Ra(i,j) + Rs(i,j)

        vdep(i,j) = 1./Rttl*R2

        !Set a minimum value of deposition velocity
        vdep(i,j) = max(vdep(i,j),1.e-4)

        !Save the dry deposition frequency for the chemical removal terms in units of s-1
        drydepf(i,j) = max(0.,vdep(i,j) / dz1(i,j))

     enddo  ! i
  enddo   ! j

end subroutine DryDepositionGOCART_Pert


subroutine ObukhovLength_Pert ( im, jm, t, rhoa, shflux, ustar, obk, von_karman, cpd, grav )

 implicit none

 !INPUTS
 integer, intent(in)                    :: im, jm                 ! grid information
 real(8), intent(in)                    :: von_karman, cpd, grav  ! constants 
 real(8), intent(in), dimension(im,jm)  :: t                      ! temperature lowest level [K]
 real(8), intent(in), dimension(im,jm)  :: rhoa                   ! air density [kg m-3]
 real(8), intent(in), dimension(im,jm)  :: ustar                  ! friction speed [m s-1]
 real(8), intent(in), dimension(im,jm)  :: shflux                 ! sfc. sens. heat flux [W m-2]

 !OUTPUTS
 real(8), intent(out), dimension(im,jm) :: obk                    ! Obukhov length [m]

 !Locals
 integer :: i, j

 obk = 1.e5

 do i = 1,im
    do j = 1,jm

       if (abs(shflux(i,j)) > 1.e-32) then

          obk(i,j) =  - rhoa(i,j) * cpd * t(i,j) * ustar(i,j)**3. / (von_karman * grav * shflux(i,j))

       endif

    enddo
 enddo

end subroutine ObukhovLength_Pert
