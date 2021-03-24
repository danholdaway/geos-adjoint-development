subroutine dummy( im, jm, radius, &
                                    fraclake, gwettop, oro, &
                                    u10m, v10m, &
                                    emissions, grav )

   implicit none

   !Inputs
   integer, intent(in)                   :: im, jm
   real(8), intent(in)                   :: grav                   
   real(8), intent(in)                   :: radius                   ! particle radius [m]
   real(8), intent(in), dimension(im,jm) :: fraclake, gwettop, oro
   real(8), intent(in), dimension(im,jm) :: u10m
   real(8), intent(in), dimension(im,jm) :: v10m

   !Outputs
   real(8), intent(out)     ::  emissions(im,jm)             ! Local emission


u10m = v10m*u10m
v10m = v10m**2*u10m

call DustEmissionGOCART_Pert( im, jm, radius, &
                                    fraclake, gwettop, oro, &
                                    u10m, v10m, &
                                    emissions, grav )

u10m = v10m*u10m
v10m = v10m**2*u10m
emissions = emissions*u10m*v10m

  end subroutine dummy


subroutine DustEmissionGOCART_Pert( im, jm, radius, &
                                    fraclake, gwettop, oro, &
                                    u10m, v10m, &
                                    emissions, grav )

   implicit none

   !Inputs
   integer, intent(in)                   :: im, jm
   real(8), intent(in)                   :: grav                   
   real(8), intent(in)                   :: radius                   ! particle radius [m]
   real(8), intent(in), dimension(im,jm) :: fraclake, gwettop, oro
   real(8), intent(in), dimension(im,jm) :: u10m
   real(8), intent(in), dimension(im,jm) :: v10m

   !Outputs
   real(8), intent(out)     ::  emissions(im,jm)             ! Local emission

   !Local
   integer            ::  i, j
   integer            ::  n1
   real(8), parameter ::  air_dens = 1.25  ! Air density = 1.25 kg m-3
   real(8), parameter ::  soil_density  = 2650.  ! km m-3
   real(8)            ::  diameter         ! dust effective diameter [m]
   real(8)            ::  u_thresh0
   real(8)            ::  u_thresh
   real(8)            ::  w10m

!  Initialize local variables
!  --------------------------
   emissions(:,:) = 0.

   diameter = 2. * radius
   u_thresh0 = 0.13 * sqrt(soil_density*grav*diameter/air_dens) &
                    * sqrt(1.+6.e-7/(soil_density*grav*diameter**2.5)) &
                    / sqrt(1.928*(1331.*(100.*diameter)**1.56+0.38)**0.092 - 1.)


!  Spatially dependent part of calculation
!  ---------------------------------------
   do j = 1, jm
      do i = 1, im

         if ( nint(oro(i,j)) == 1 ) then ! only over LAND gridpoints

            w10m = sqrt(u10m(i,j)**2.+v10m(i,j)**2.)

            !Modify the threshold depending on soil moisture as in Ginoux et al. [2001]
            if (gwettop(i,j) .lt. 0.5) then

               u_thresh = max(0.,u_thresh0* (1.2+0.2*log(max(1.e-3,gwettop(i,j)))))

               if (w10m .gt. u_thresh) then     

                  !Emission of dust [kg m-2 s-1]
                  emissions(i,j) = (1.-fraclake(i,j)) * w10m**2. * (w10m-u_thresh)

               endif
            endif

         endif

      enddo   ! i
   enddo    ! j

  end subroutine DustEmissionGOCART_Pert
