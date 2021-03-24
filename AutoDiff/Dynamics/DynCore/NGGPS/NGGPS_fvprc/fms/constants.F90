module constants_mod

implicit none
private

integer, parameter :: FVPRC = 8

real(FVPRC) :: realnumber
real(FVPRC), public, parameter :: RADIUS = 6371.0e3   
real(FVPRC), public, parameter :: OMEGA  = 7.292e-5 
real(FVPRC), public, parameter :: GRAV   = 9.80    
real(FVPRC), public, parameter :: RDGAS  = 287.04 
real(FVPRC), public, parameter :: KAPPA  = 2./7.  
real(FVPRC), public, parameter :: CP_AIR = RDGAS/KAPPA 
real(FVPRC), public, parameter :: CP_OCEAN = 3989.24495292815
real(FVPRC), public, parameter :: RHO0    = 1.035e3
real(FVPRC), public, parameter :: RHO0R   = 1.0/RHO0
real(FVPRC), public, parameter :: RHO_CP  = RHO0*CP_OCEAN
real(FVPRC), public, parameter :: ES0 = 1.0 
real(FVPRC), public, parameter :: RVGAS = 461.50 
real(FVPRC), public, parameter :: CP_VAPOR = 4.0*RVGAS
real(FVPRC), public, parameter :: DENS_H2O = 1000. 
real(FVPRC), public, parameter :: HLV = 2.500e6   
real(FVPRC), public, parameter :: HLF = 3.34e5   
real(FVPRC), public, parameter :: HLS = HLV + HLF
real(FVPRC), public, parameter :: TFREEZE = 273.16    
real(FVPRC), public, parameter :: WTMAIR = 2.896440E+01
real(FVPRC), public, parameter :: WTMH2O = WTMAIR*(RDGAS/RVGAS)
real(FVPRC), public, parameter :: WTMOZONE =  47.99820
real(FVPRC), public, parameter :: WTMC     =  12.00000
real(FVPRC), public, parameter :: WTMCO2   =  44.00995
real(FVPRC), public, parameter :: WTMO2    =  31.9988
real(FVPRC), public, parameter :: WTMCFC11 = 137.3681
real(FVPRC), public, parameter :: WTMCFC12 = 120.9135
real(FVPRC), public, parameter :: DIFFAC = 1.660000E+00
real(FVPRC), public, parameter :: SECONDS_PER_DAY  = 8.640000E+04, SECONDS_PER_HOUR = 3600., SECONDS_PER_MINUTE=60.
real(FVPRC), public, parameter :: AVOGNO = 6.023000E+23
real(FVPRC), public, parameter :: PSTD   = 1.013250E+06
real(FVPRC), public, parameter :: PSTD_MKS    = 101325.0
real(FVPRC), public, parameter :: RADCON = ((1.0E+02*GRAV)/(1.0E+04*CP_AIR))*SECONDS_PER_DAY
real(FVPRC), public, parameter :: RADCON_MKS  = (GRAV/CP_AIR)*SECONDS_PER_DAY
real(FVPRC), public, parameter :: O2MIXRAT    = 2.0953E-01
real(FVPRC), public, parameter :: RHOAIR      = 1.292269
real(FVPRC), public, parameter :: ALOGMIN     = -50.0
real(FVPRC), public, parameter :: STEFAN  = 5.6734e-8 
real(FVPRC), public, parameter :: VONKARM = 0.40     
real(FVPRC), public, parameter :: PI      = 3.14159265358979323846
real(FVPRC), public, parameter :: RAD_TO_DEG=180./PI
real(FVPRC), public, parameter :: DEG_TO_RAD=PI/180.
real(FVPRC), public, parameter :: RADIAN  = RAD_TO_DEG
real(FVPRC), public, parameter :: C2DBARS = 1.e-4
real(FVPRC), public, parameter :: KELVIN  = 273.15
real(FVPRC), public, parameter :: EPSLN   = 1.0e-30

public :: constants_init

contains

subroutine constants_init

! dummy routine.

end subroutine constants_init

end module constants_mod


