module fv_cmp_mod

use fv_arrays_mod, only : r_grid

  implicit none
 private

 public fv_sat_adj, qs_init

contains

 subroutine fv_sat_adj(mdt, zvir, is, ie, js, je, ng, hydrostatic, consv_te, &
                       te0, qv, ql, qi, qr, qs, qg, dpln, delz, pt, dp,  &
                       q_con, cappa, area, dtdt, out_dt, last_step, do_qa, qa)
! This is designed for 6-class micro-physics schemes; handles the heat release
! due to in situ phase changes
! input pt is T_vir
 integer, intent(in):: is, ie, js, je, ng
 real, intent(in):: mdt ! remapping time step
 real, intent(in):: zvir
 logical, intent(in):: hydrostatic, consv_te, out_dt
 logical, intent(in):: last_step
 logical, intent(in):: do_qa
 real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng):: dp, delz
 real, intent(in):: dpln(is:ie,js:je)
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng):: pt, qv, ql, qi, qr, qs, qg
 real, intent(out):: qa(is-ng:ie+ng,js-ng:je+ng)
 real(kind=R_GRID), intent(in), dimension(is-ng:ie+ng,js-ng:je+ng):: area
 real, intent(inout), dimension(is-ng:,js-ng:):: q_con
 real, intent(inout), dimension(is-ng:,js-ng:):: cappa
 real, intent(inout)::dtdt(is:ie,js:je)
 real, intent(out):: te0(is-ng:ie+ng,js-ng:je+ng)
!---
 real, dimension(is:ie):: wqsat, dq2dt, qpz, cvm, t0, pt1, icp2, lcp2, tcp2, tcp3,    &
                          den, q_liq, q_sol, src, hvar
 real, dimension(is:ie):: mc_air, lhl, lhi  ! latent heat
 real:: sink, qsw, rh, fac_v2l, fac_l2v
 real:: tc, qsi, dqsdt, dq, dq0, pidep, qi_crt, tmp, dtmp
 real:: condensates, tin, qstar, rqi, q_plus, q_minus
 real:: sdt, dt_Bigg, adj_fac, fac_s, fac_r, fac_i2s, fac_mlt, fac_l2r
 real:: factor, qim, tice0, c_air, c_vap
 integer i,j


end subroutine fv_sat_adj


 subroutine qs_init(kmp)
 integer, intent(in):: kmp
 integer, parameter:: length=2621
 real, parameter:: rhor  = 1.0e3  ! LFO83
 real, parameter:: vdifu = 2.11e-5
 real, parameter:: tcond = 2.36e-2
 real, parameter:: visk = 1.259e-5
 real, parameter:: hltc = 2.5e6
 real, parameter:: gam290 = 1.827363
 real, parameter:: gam380 = 4.694155
 real, parameter:: alin = 842.0
!Intercept parameters
 real, parameter:: rnzr = 8.0e6
 real, parameter:: c_cracw = 0.9      ! rain accretion efficiency
 real:: scm3, act2
 integer i

!Stubbed, nonlinear?

 end subroutine qs_init

end module fv_cmp_mod
