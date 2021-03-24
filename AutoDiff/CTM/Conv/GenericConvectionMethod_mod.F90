!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GenericConvectionMethod_mod 

module GenericConvectionMethod_mod

use ESMF
USE convectiveTransport_mod
use GmiArrayBundlePointer_mod

implicit none
private

public doGenericConvectiveTransport

  CONTAINS

      subroutine doGenericConvectiveTransport (det_ent, do_downdraft, pbl, cldmas, &
                      dtrn, eu, ed, md, grid_height, mass, kel, press3e,       &
                      concentration, isFixedConcentration, mcor, tdt,          &
                      i1, i2, j1, j2, k1, k2, ilong, ivert, num_species)
!
! !INPUT PARAMETER:
      integer, intent(in) :: i1, i2, j1, j2, k1, k2
      integer, intent(in) :: ilong, ivert, num_species
      logical, intent(in) :: isFixedConcentration(:)
      logical, intent(in) :: det_ent ! flag for doing detrainment then entrainment
      logical, intent(in) :: do_downdraft ! flag for doing downdrafts
      real*8 , intent(in) :: tdt          ! model time step  (s)
      real*8 , intent(in) :: mcor       (i1:i2, j1:j2) ! area of each grid box (m^2)
      real*8 , intent(in) :: pbl        (i1:i2,j1:j2) ! planetary boundary layer thickness (m)
      real*8 , intent(in) :: cldmas     (i1:i2,j1:j2,k1:k2) ! convective mass flux in     updraft (kg/m^2/s)
      real*8 , intent(in) :: dtrn       (i1:i2,j1:j2,k1:k2) ! detrainment rate (DAO:kg/m^2*s, NCAR:s^-1)
      real*8 , intent(in) :: eu         (i1:i2,j1:j2,k1:k2) ! ntrainment into convective updraft (s^-1)
      real*8 , intent(in) :: ed         (i1:i2,j1:j2,k1:k2) ! entrainment into convective downdraft (s^-1)
      real*8 , intent(in) :: md         (i1:i2,j1:j2,k1:k2) ! convective mass flux in downdraft (kg/m^2/s)
      real*8 , intent(in) :: grid_height(i1:i2,j1:j2,k1:k2) ! grid box height  (m)
      real*8 , intent(in) :: mass       (i1:i2,j1:j2,k1:k2) ! mass of air in each grid box (kg)
      real*8 , intent(in) :: kel        (i1:i2,j1:j2,k1:k2) ! temperature      (degK)
      real*8 , intent(in) :: press3e    (i1:i2,j1:j2,k1-1:k2) ! atmospheric pressure at the edge of each grid box (mb)
!
! !INPUT/OUTPUT PARAMETERS:
                             ! species concentration, known at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!
! !DEFINED PARAMETERS:
      real*8, parameter :: MBSTH = 1.0d-15 ! threshold below which we treat 
                                           ! mass fluxes as zero (mb/s)
      real*8, parameter :: GMI_G  =   9.81d0 ! mean surface gravity accel. (m/s^2)
      real*8, parameter :: PASPMB = 100.00d0 ! pascals  per millibar
!
! !DESCRIPTION:
! This is the interface routine to Convective Transport.  
! It formats the gem variables to satisfy the Convective Transport routine.
!
!EOP
!------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      integer :: iku
      integer :: il, ij, ik, ic
      integer :: il2g                      ! gathered index to operate over
      integer :: itdt_conv
      integer :: num_conv_steps
      real*8  :: rnum_conv_steps
      real*8  :: tdt_conv
      real*8  :: xmbsth
      integer :: ideep(ilong)              ! gathering array
      integer :: pbli (ilong)              ! index of pbl height
      real*8  :: updraft_velocity(i1:i2)   ! velocity in convective updraft
                                           ! (m/s)
      real*8  :: dpi(ilong,  k1:k2)        ! delta pressure between interfaces
      real*8  :: dui(ilong,  k1:k2)        ! mass detraining from updraft
      real*8  :: eui(ilong,  k1:k2)        ! mass entraining into updraft
      real*8  :: mui(ilong,  k1:k2)        ! mass flux up
      real*8  :: mdi(ilong,  k1:k2)        ! mass flux down
      real*8  :: fracis(i1:i2, k1:k2, num_species)  ! insoluble fraction of tracer
      real*8  :: qq    (i1:i2, k1:k2, num_species)  ! tracer array including moisture

      ideep(:) = 0
      pbli (:) = 0

      dpi(:,:) = 0.0d0
      dui(:,:) = 0.0d0
      eui(:,:) = 0.0d0
      mui(:,:) = 0.0d0
      mdi(:,:) = 0.0d0

      fracis(:,:,:) = 0.0d0

      xmbsth = MBSTH

      updraft_velocity(:) = 0.0d0

      ! -----------------------------------------------------------
      ! Calculate any needed sub-cycling of the convection operator
      ! by comparing the mass flux over the full time step to the
      ! mass within the grid box.
      ! -----------------------------------------------------------

      num_conv_steps = int(Maxval (tdt * cldmas(:,:,:) * Spread (mcor(:,:), 3, k2-k1+1) /  &
                               mass(:,:,:)) + 1.0d0)

      rnum_conv_steps = real(num_conv_steps,8)
      tdt_conv        = tdt / rnum_conv_steps

      IJLOOP: do ij = j1, j2
         il2g = 0

         ILLOOP: do il = i1, i2
            ! ----------------------------------------------------
            ! Verify that there is convection in the current cell.
            ! ----------------------------------------------------
            if (Maxval (cldmas(il,ij,:)) >= 0.0d0) then
               il2g = il2g + 1
               ideep(il2g) = il

               dpi(il2g,:) = (press3e(il,ij,k2-1:k1-1:-1) -  &
                              press3e(il,ij,k2:k1:-1)) * PASPMB

               mui(il2g,:) = cldmas(il,ij,k2:k1:-1) * GMI_G

               eui(il2g,k1+1:k2-1) = Max (0.0d0, cldmas(il,ij,k2-1:k1+1:-1) -  &
                                                 cldmas(il,ij,k2-2:k1  :-1) +  &
                                                 dtrn  (il,ij,k2-1:k1+1:-1))

               eui(il2g,:) = eui(il2g,:) * GMI_G

               if (det_ent) dui(il2g,:) = dtrn(il,ij,k2:k1:-1) * GMI_G
            end if
         end do ILLOOP

         IL2GIF: if (il2g /= 0) then
            fracis(:,:,:) = 1.0d0

            ! ----------------------------------------------------------------
            ! Find the index of the top of the planetary boundary layer (pbl).
            ! Convection will assume well mixed tracers below that level.
            ! ----------------------------------------------------------------
            do il = 1, il2g
               pbli(il) = 0

               IKLOOP: do ik = k1, k2
                  if (pbl(ideep(il),ij) < Sum (grid_height(ideep(il),ij,k1:ik))) then
                     pbli(il) = ik
                     exit IKLOOP
                  end if
               end do IKLOOP

               if (pbli(il) == 0) then
                  err_msg = 'Could not find pbl in doGenericConvectiveTransport.'
                  PRINT*, err_msg
                  stop
               end if

               pbli(il) = k2 - pbli(il)
            end do

            ITDTCLOOP: do itdt_conv = 1, num_conv_steps
               do ic = 1, num_species
                  qq(:,k2:k1:-1,ic) = concentration(ic)%pArray3D(:,ij,k1:k2)
               end do

               call convectiveTransport (il2g, tdt_conv, xmbsth, ideep,  &
                           pbli, dui, eui, mui, mdi, dpi, fracis, qq, &
                           isFixedConcentration, i1, i2, k1, k2, ilong, num_species)

               do ic = 1, num_species
                  concentration(ic)%pArray3D(:,ij,k1:k2) = qq(:,k2:k1:-1,ic)
               end do
            end do ITDTCLOOP
         end if IL2GIF
      end do IJLOOP

      return

      end subroutine doGenericConvectiveTransport

  end module GenericConvectionMethod_mod
