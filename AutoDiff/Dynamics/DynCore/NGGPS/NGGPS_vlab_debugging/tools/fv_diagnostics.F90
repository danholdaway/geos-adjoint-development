module fv_diagnostics_mod

use fv_arrays_mod, only: r_grid
use mpp_domains_mod, only: domain2d

implicit none
private

public fv_time, prt_mxm, range_check, prt_minmax, prt_maxmin, sphum_ll_fix

logical, parameter :: prt_minmax = .false.
integer :: fv_time !dummy
real :: sphum_ll_fix

contains


 subroutine range_check(qname, q, is, ie, js, je, n_g, km, pos, q_low, q_hi, bad_range)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je
      integer, intent(in):: n_g, km
      real, intent(in)::    q(is-n_g:ie+n_g, js-n_g:je+n_g, km)
      real, intent(in):: pos(is-n_g:ie+n_g, js-n_g:je+n_g,2)
      real, intent(in):: q_low, q_hi
      logical, optional, intent(out):: bad_range

      !Stubbed

 end subroutine range_check

 subroutine prt_mxm(qname, q, is, ie, js, je, n_g, km, fac, area, domain)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je
      integer, intent(in):: n_g, km
      real, intent(in)::    q(is-n_g:ie+n_g, js-n_g:je+n_g, km)
      real, intent(in)::    fac
      real(kind=R_GRID), intent(IN)::    area(is-3:ie+3, js-3:je+3)
      type(domain2d), intent(INOUT) :: domain

      !Stubbed

 end subroutine prt_mxm

 subroutine prt_maxmin(qname, q, is, ie, js, je, n_g, km, fac)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je
      integer, intent(in):: n_g, km
      real, intent(in)::    q(is-n_g:ie+n_g, js-n_g:je+n_g, km)
      real, intent(in)::    fac

      !Stubbed

 end subroutine prt_maxmin


end module fv_diagnostics_mod

