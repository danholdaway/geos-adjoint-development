module nh_core_mod

   use fv_arrays_mod,  only: g_precision, p_precision, pg_precision
   use mpp_mod,        only: FATAL, mpp_error

   implicit none
   private

   public Riem_Solver, Riem_Solver_C, update_dz_c, update_dz_d

CONTAINS 

  subroutine update_dz_c(is, ie, js, je, km, ng, area,   &
                         zh, ut, vt, dz_in, dz_out, wk)
! !INPUT PARAMETERS:
  integer, intent(in):: is, ie, js, je, ng, km
  real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: ut, vt, zh
  real(g_precision), intent(in ):: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(in ):: dz_in (is:ie,js:je,km) 
  real, intent(out):: dz_out(is:ie,js:je,km) 
  real(pg_precision), intent(out):: wk(is-ng:ie+ng,js-ng:je+ng,km+1)  ! work array

  call mpp_error(FATAL,'nh core stubbed')

  end subroutine update_dz_c



  subroutine update_dz_d(hord, is, ie, js, je, km, ng, npx, npy, area,    &
                         zh, crx, cry, xfx, yfx, delz, wk, delp, n_sponge)

  integer, intent(in):: is, ie, js, je, ng, km, npx, npy
  integer, intent(in):: hord, n_sponge
  real(g_precision), intent(in)   :: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(inout) ::  zh(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout) ::delz(is:ie,js:je,km)
  real, intent(inout) ::delp(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout), dimension(is:ie+1,js-ng:je+ng,km):: crx, xfx
  real, intent(inout), dimension(is-ng:ie+ng,js:je+1,km):: cry, yfx
  real, intent(  out) ::   wk(is:ie,js:je,km)

  call mpp_error(FATAL,'nh core stubbed')

  end subroutine update_dz_d


  subroutine Riem_Solver_C(dt,   is,  ie,   js, je, km,   ng,  &
                           akap, rdgas, grav, cp,  ptop, hs, w,  delz, pt,  &
                           delp, gz,  pk,   ip)

   integer, intent(in):: is, ie, js, je, ng, km
   integer, intent(in):: ip
   real, intent(in):: dt
   real(p_precision) :: ptop
   real, intent(in)::  akap, rdgas, grav, cp
   real, intent(in):: delz(is:ie,js:je,km)
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: pt, delp
   real, intent(in)::       hs(is-ng:ie+ng,js-ng:je+ng)
   real, intent(inout):: w(is-ng:ie+ng,js-ng:je+ng,km)
! OUTPUT PARAMETERS 
   real(pg_precision), intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz
   real(pg_precision), intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: pk

  call mpp_error(FATAL,'nh core stubbed')

  end subroutine Riem_Solver_C


  subroutine Riem_Solver(dt,   is,   ie,   js, je, km, ng,    &
                         akap, rdgas, grav, cp,   ptop, hs, peln, w,  delz, pt,  &
                         delp, gz,   pkc, pk, pe, last_call, ip)

   integer, intent(in):: is, ie, js, je, km, ng
   integer, intent(in):: ip      
   real, intent(in):: dt         
   real(p_precision) :: ptop
   real, intent(in):: akap, rdgas, grav, cp
   real, intent(in):: hs(is-ng:ie+ng,js-ng:je+ng)
   logical, intent(in):: last_call
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km):: w, delp, pt
   real, intent(inout):: delz(is:ie,js:je,km)
   real(pg_precision), intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz
   real(pg_precision), intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: pkc
   real(p_precision), intent(out):: pk(is:ie,js:je,km+1)
   real(p_precision), intent(out):: pe(is-1:ie+1,km+1,js-1:je+1)
   real(p_precision), intent(out):: peln(is:ie,km+1,js:je)

   call mpp_error(FATAL,'nh core stubbed')

  end subroutine Riem_Solver

end module nh_core_mod
