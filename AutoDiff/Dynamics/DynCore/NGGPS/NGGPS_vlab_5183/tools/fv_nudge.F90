module fv_nwp_nudge_mod

 use fv_arrays_mod,     only: fv_grid_bounds_type, fv_grid_type
 use mpp_domains_mod,   only: domain2d
 use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed

implicit none
private

public do_adiabatic_init, breed_slp_inline

logical, parameter :: do_adiabatic_init = .false.

contains

 subroutine breed_slp_inline(nstep, dt, npz, ak, bk, phis, pe, pk, peln, pkz, delp, u, v, pt, q, nwat,   &
                             zvir, gridstruct, ks, domain_local, bd, hydrostatic)
!------------------------------------------------------------------------------------------
! Purpose:  Vortex-breeding by nudging sea-level-pressure towards single point observations
! Note: conserve water mass, geopotential, and momentum at the expense of dry air mass
!------------------------------------------------------------------------------------------
! Input
      integer, intent(in):: nstep, npz, nwat, ks
      real, intent(in):: dt       ! (small) time step in seconds
      real, intent(in):: zvir
      real, intent(in), dimension(npz+1):: ak, bk
      logical, intent(in):: hydrostatic
      type(fv_grid_bounds_type), intent(IN) :: bd
      real, intent(in):: phis(isd:ied,jsd:jed)
      type(domain2d), intent(INOUT) :: domain_local
! Input/Output
      real, intent(inout):: u(isd:ied,jsd:jed+1,npz)
      real, intent(inout):: v(isd:ied+1,jsd:jed,npz)
      real, intent(inout), dimension(isd:ied,jsd:jed,npz):: delp, pt
      real, intent(inout)::q(isd:ied,jsd:jed,npz,*)

      real, intent(inout):: pk(is:ie,js:je, npz+1)          ! pe**kappa
      real, intent(inout):: pe(is-1:ie+1, npz+1,js-1:je+1)  ! edge pressure (pascal)
      real, intent(inout):: pkz(is:ie,js:je,npz) 
      real, intent(out):: peln(is:ie,npz+1,js:je)           ! ln(pe)

      type(fv_grid_type), target :: gridstruct

      !Stubbed, nonlinear
    
  end subroutine breed_slp_inline

end module fv_nwp_nudge_mod

