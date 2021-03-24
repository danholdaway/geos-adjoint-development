module fv_mapz_tlm_mod

  use fv_arrays_mod,      only: p_precision
  use fv_grid_tools_mod,  only: area, dx, dy, rdxa, rdya
  use fv_grid_utils_mod,  only: g_sum, ptop, ptop_min
  use fv_grid_utils_mod,  only: cosa_s, rsin2
  use fv_fill_tlm_mod,    only: fillz_tlm
  use fv_mp_mod,          only: domain
  use mpp_domains_mod,    only: mpp_update_domains

  implicit none

  real(p_precision), parameter::  D0_5 = 0.5, D0_0 = 0.0
  real(p_precision), parameter::  r3 = 1./3., r23 = 2./3., r12 = 1./12.
  real(kind=4) :: E_Flux
  private

  public compute_total_energy_tlm, Lagrangian_to_Eulerian_tlm

CONTAINS

  SUBROUTINE LAGRANGIAN_TO_EULERIAN_TLM(do_consv, consv, ps, pe, pe_tl, &
&   delp, delp_tl, pkz, pkz_tl, pk, pk_tl, pdt, km, is, ie, js, je, isd&
&   , ied, jsd, jed, nq, sphum, u, u_tl, v, v_tl, w, w_tl, delz, delz_tl&
&   , pt, pt_tl, q, q_tl, hs, r_vir, cp, akap, pi, radius, grav, kord_mt&
&   , kord_tr, kord_tm, kord_mt_pert, kord_tr_pert, kord_tm_pert, peln&
&   , peln_tl, te0_2d, te0_2d_tl, ng, ua, ua_tl&
&   , va, omga, omga_tl, te, te_tl, pem, fill, reproduce_sum, ak, bk, ks&
&   , ze0, ze0_tl, te_method, remap_t, hydrostatic, hybrid_z, do_omega, &
&   ktop, ncnst, mfx, mfy)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: do_consv
! phys time step
    REAL, INTENT(IN) :: pdt
    INTEGER, INTENT(IN) :: km
! number of tracers (including h2o)
    INTEGER, INTENT(IN) :: nq
! index for water vapor (specific humidity)
    INTEGER, INTENT(IN) :: sphum
    INTEGER, INTENT(IN) :: ng
! starting & ending X-Dir index
    INTEGER, INTENT(IN) :: is, ie, isd, ied
! starting & ending Y-Dir index
    INTEGER, INTENT(IN) :: js, je, jsd, jed
    INTEGER, INTENT(IN) :: ks, ktop
! Mapping oder for the vector winds
    INTEGER, INTENT(IN) :: kord_mt
! Mapping oder for tracers
    INTEGER, INTENT(IN) :: kord_tr
! Mapping oder for thermodynamics
    INTEGER, INTENT(IN) :: kord_tm
    INTEGER, INTENT(IN) :: kord_mt_pert, kord_tr_pert, kord_tm_pert
! Added 4TAF
    INTEGER, INTENT(IN) :: ncnst
! factor for TE conservation
    REAL, INTENT(IN) :: consv
    REAL, INTENT(IN) :: r_vir
    REAL, INTENT(IN) :: cp
    REAL, INTENT(IN) :: akap
    REAL, INTENT(IN) :: pi, radius, grav
! surface geopotential
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: te0_2d(is:ie, js:je)
    REAL, INTENT(IN) :: te0_2d_tl(is:ie, js:je)
! fill negative tracers
    LOGICAL, INTENT(IN) :: fill
    LOGICAL, INTENT(IN) :: reproduce_sum
    LOGICAL, INTENT(IN) :: do_omega
    REAL, INTENT(IN) :: ak(km+1)
    REAL, INTENT(IN) :: bk(km+1)
! !INPUT/OUTPUT
! pe to the kappa
    REAL(p_precision), INTENT(INOUT) :: pk(is:ie, js:je, km+1)
    REAL(p_precision), INTENT(INOUT) :: pk_tl(is:ie, js:je, km+1)
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, km, *)
    REAL, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed, km, *)
! pressure thickness
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: delp_tl(isd:ied, jsd:jed, km)
! pressure at layer edges
    REAL(p_precision), INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL(p_precision), INTENT(INOUT) :: pe_tl(is-1:ie+1, km+1, js-1:je+1&
&   )
    REAL, INTENT(INOUT) :: pem(is-1:ie+1, km+1, js-1:je+1)
! surface pressure
    REAL(p_precision), INTENT(INOUT) :: ps(isd:ied, jsd:jed)
! Specified height at edges (m)
    REAL, INTENT(INOUT) :: ze0(is:ie, js:je, km+1)
    REAL, INTENT(INOUT) :: ze0_tl(is:ie, js:je, km+1)
! u-wind will be ghosted one latitude to the north upon exit
! u-wind (m/s)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: u_tl(isd:ied, jsd:jed+1, km)
! v-wind (m/s)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_tl(isd:ied+1, jsd:jed, km)
! vertical velocity (m/s)
    REAL, INTENT(INOUT) :: w(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: w_tl(isd:ied, jsd:jed, km)
! cp*virtual potential temperature 
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: pt_tl(isd:ied, jsd:jed, km)
! as input; output: temperature
! delta-height (m)
    REAL, INTENT(INOUT) :: delz(is:ie, js:je, km)
    REAL, INTENT(INOUT) :: delz_tl(is:ie, js:je, km)
    INTEGER, INTENT(IN) :: te_method
    LOGICAL, INTENT(IN) :: remap_t
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: hybrid_z
! u-wind (m/s) on physics grid
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: ua_tl(isd:ied, jsd:jed, km)
! v-wind (m/s) on physics grid
    REAL, INTENT(INOUT) :: va(isd:ied, jsd:jed, km)
! vertical press. velocity (pascal/sec)
    REAL, INTENT(INOUT) :: omga(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: omga_tl(isd:ied, jsd:jed, km)
! log(pe)
    REAL(p_precision), INTENT(INOUT) :: peln(is:ie, km+1, js:je)
    REAL(p_precision), INTENT(INOUT) :: peln_tl(is:ie, km+1, js:je)
! layer-mean pk for converting t to pt
    REAL(p_precision), INTENT(OUT) :: pkz(is:ie, js:je, km)
    REAL(p_precision), INTENT(OUT) :: pkz_tl(is:ie, js:je, km)
    REAL, INTENT(OUT) :: te(is:ie, js:je, km)
    REAL, INTENT(OUT) :: te_tl(is:ie, js:je, km)
! Mass fluxes
! X-dir Mass Flux
    REAL, OPTIONAL, INTENT(INOUT) :: mfx(is:ie+1, js:je, km)
! Y-dir Mass Flux
    REAL, OPTIONAL, INTENT(INOUT) :: mfy(is:ie, js:je+1, km)
! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!-----------------------------------------------------------------------
    INTEGER :: i, j, k
! numerical tracer source from surface
    REAL(p_precision) :: q_source(is:ie, js:je, nq)
! in case fillz is not sufficient
    REAL(p_precision) :: te_2d(is:ie, js:je)
    REAL(p_precision) :: te_2d_tl(is:ie, js:je)
    REAL(p_precision) :: zsum0(is:ie, js:je)
    REAL(p_precision) :: zsum0_tl(is:ie, js:je)
    REAL(p_precision) :: zsum1(is:ie, js:je)
    REAL(p_precision) :: zsum1_tl(is:ie, js:je)
    REAL(p_precision) :: q2(is:ie, km)
    REAL(p_precision) :: q2_tl(is:ie, km)
    REAL(p_precision) :: dp2(is:ie, km)
    REAL(p_precision) :: dp2_tl(is:ie, km)
    REAL(p_precision) :: pe1(is:ie, km+1)
    REAL(p_precision) :: pe1_tl(is:ie, km+1)
    REAL(p_precision) :: pe2(is:ie, km+1)
    REAL(p_precision) :: pe2_tl(is:ie, km+1)
    REAL(p_precision) :: pk1(is:ie, km+1)
    REAL(p_precision) :: pk1_tl(is:ie, km+1)
    REAL(p_precision) :: pk2(is:ie, km+1)
    REAL(p_precision) :: pk2_tl(is:ie, km+1)
    REAL(p_precision) :: pn1(is:ie, km+1)
    REAL(p_precision) :: pn2(is:ie, km+1)
    REAL(p_precision) :: pn2_tl(is:ie, km+1)
    REAL(p_precision) :: pe0(is:ie+1, km+1)
    REAL(p_precision) :: pe0_tl(is:ie+1, km+1)
    REAL(p_precision) :: pe3(is:ie+1, km+1)
    REAL(p_precision) :: pe3_tl(is:ie+1, km+1)
    REAL(p_precision) :: phis(is:ie, km+1)
    REAL(p_precision) :: phis_tl(is:ie, km+1)
    REAL(p_precision) :: gz(is:ie)
    REAL(p_precision) :: gz_tl(is:ie)
! for nonhydrostatic option with hybrid_z coordinate
    REAL(p_precision) :: ze1(is:ie, km+1), ze2(is:ie, km+1), deng(is:ie&
&   , km)
    REAL(p_precision) :: ze1_tl(is:ie, km+1), ze2_tl(is:ie, km+1), &
&   deng_tl(is:ie, km)
    REAL :: dz1(km), ztop, z_rat
    REAL :: z_rat_tl
    REAL(p_precision) :: rcp, rg, ak1, tmp, tpe, cv, rgama, rrg
    REAL(p_precision) :: tmp_tl, tpe_tl
    REAL(p_precision) :: bkh
    REAL(p_precision) :: dtmp
    REAL(p_precision) :: dtmp_tl
    REAL(p_precision) :: dlnp
    REAL(p_precision) :: dlnp_tl
    INTEGER :: iq, n, kp, k_next
    LOGICAL :: te_map
    REAL(p_precision) :: k1k, kapag
    REAL(p_precision) :: tmpsum
    REAL(p_precision) :: tmpsum_tl
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC ABS
    INTRINSIC PRESENT
    INTRINSIC REAL
    REAL*8 :: arg1
    REAL*8 :: arg1_tl
    REAL(p_precision) :: arg2
    REAL(p_precision) :: arg2_tl
    INTEGER :: abs4
    INTEGER :: abs3
    INTEGER :: abs2
    INTEGER :: abs1
    INTEGER :: abs0
! rg/Cv=0.4
    k1k = akap/(1.-akap)
    kapag = -(akap/grav)
    rg = akap*cp
    cv = cp - rg
! cv/cp
    rgama = 1. - akap
    rcp = 1./cp
    ak1 = (akap+1.)/akap
    rrg = -(rg/grav)
    IF (kord_tm .LT. 0) THEN
      te_map = .false.
! remap_t test
      IF (remap_t) THEN
! hydro test
! Note: pt at this stage is cp*Theta_v
! Transform virtual pt to virtual Temp
        IF (hydrostatic) THEN
          DO k=1,km
            DO j=js,je
              DO i=is,ie
                pt_tl(i, j, k) = ((pt_tl(i, j, k)*(pk(i, j, k+1)-pk(i, j&
&                 , k))+pt(i, j, k)*(pk_tl(i, j, k+1)-pk_tl(i, j, k)))*&
&                 rg*(peln(i, k+1, j)-peln(i, k, j))-pt(i, j, k)*(pk(i, &
&                 j, k+1)-pk(i, j, k))*rg*(peln_tl(i, k+1, j)-peln_tl(i&
&                 , k, j)))/(rg*(peln(i, k+1, j)-peln(i, k, j)))**2
                pt(i, j, k) = pt(i, j, k)*(pk(i, j, k+1)-pk(i, j, k))/(&
&                 rg*(peln(i, k+1, j)-peln(i, k, j)))
              END DO
            END DO
          END DO
        ELSE
          IF (ktop .GT. 1) THEN
            DO k=1,ktop-1
              DO j=js,je
                DO i=is,ie
                  pt_tl(i, j, k) = ((pt_tl(i, j, k)*(pk(i, j, k+1)-pk(i&
&                   , j, k))+pt(i, j, k)*(pk_tl(i, j, k+1)-pk_tl(i, j, k&
&                   )))*rg*(peln(i, k+1, j)-peln(i, k, j))-pt(i, j, k)*(&
&                   pk(i, j, k+1)-pk(i, j, k))*rg*(peln_tl(i, k+1, j)-&
&                   peln_tl(i, k, j)))/(rg*(peln(i, k+1, j)-peln(i, k, j&
&                   )))**2
                  pt(i, j, k) = pt(i, j, k)*(pk(i, j, k+1)-pk(i, j, k))/&
&                   (rg*(peln(i, k+1, j)-peln(i, k, j)))
                END DO
              END DO
            END DO
          END IF
          DO k=ktop,km
            DO j=js,je
              DO i=is,ie
                arg1_tl = (kapag*delp_tl(i, j, k)*delz(i, j, k)-kapag*&
&                 delp(i, j, k)*delz_tl(i, j, k))*pt(i, j, k)/delz(i, j&
&                 , k)**2 + kapag*delp(i, j, k)*pt_tl(i, j, k)/delz(i, j&
&                 , k)
                arg1 = kapag*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
                arg2_tl = k1k*arg1_tl/arg1
                arg2 = k1k*LOG(arg1)
                pt_tl(i, j, k) = rcp*(pt_tl(i, j, k)*EXP(arg2)+pt(i, j, &
&                 k)*arg2_tl*EXP(arg2))
                pt(i, j, k) = rcp*pt(i, j, k)*EXP(arg2)
              END DO
            END DO
          END DO
        END IF
      END IF
    ELSE
      te_map = .true.
      CALL PKEZ_TLM(km, is, ie, js, je, pe, pk, pk_tl, akap, peln, &
&             peln_tl, pkz, pkz_tl)
!           call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km, 1)
! Compute cp*T + KE
      DO k=1,km
        DO j=js,je
          DO i=is,ie
            te_tl(i, j, k) = 0.25*rsin2(i, j)*(2*u(i, j, k)*u_tl(i, j, k&
&             )+2*u(i, j+1, k)*u_tl(i, j+1, k)+2*v(i, j, k)*v_tl(i, j, k&
&             )+2*v(i+1, j, k)*v_tl(i+1, j, k)-cosa_s(i, j)*((u_tl(i, j&
&             , k)+u_tl(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k&
&             )+u(i, j+1, k))*(v_tl(i, j, k)+v_tl(i+1, j, k)))) + pt_tl(&
&             i, j, k)*pkz(i, j, k) + pt(i, j, k)*pkz_tl(i, j, k)
            te(i, j, k) = 0.25*rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**&
&             2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*&
&             (v(i, j, k)+v(i+1, j, k))*cosa_s(i, j)) + pt(i, j, k)*pkz(&
&             i, j, k)
          END DO
        END DO
      END DO
    END IF
    IF (.NOT.hydrostatic .AND. (.NOT.hybrid_z)) THEN
      DO k=1,km
        DO j=js,je
          DO i=is,ie
! ="specific volume"/grav
            delz_tl(i, j, k) = -((delz_tl(i, j, k)*delp(i, j, k)-delz(i&
&             , j, k)*delp_tl(i, j, k))/delp(i, j, k)**2)
            delz(i, j, k) = -(delz(i, j, k)/delp(i, j, k))
          END DO
        END DO
      END DO
      phis_tl = 0.0_8
      pe0_tl = 0.0_8
      pe1_tl = 0.0_8
      pe2_tl = 0.0_8
      pe3_tl = 0.0_8
      ze1_tl = 0.0_8
      dp2_tl = 0.0_8
      ze2_tl = 0.0_8
      q2_tl = 0.0_8
      deng_tl = 0.0_8
      pn2_tl = 0.0_8
      pk1_tl = 0.0_8
      pk2_tl = 0.0_8
    ELSE
      phis_tl = 0.0_8
      pe0_tl = 0.0_8
      pe1_tl = 0.0_8
      pe2_tl = 0.0_8
      pe3_tl = 0.0_8
      ze1_tl = 0.0_8
      dp2_tl = 0.0_8
      ze2_tl = 0.0_8
      q2_tl = 0.0_8
      deng_tl = 0.0_8
      pn2_tl = 0.0_8
      pk1_tl = 0.0_8
      pk2_tl = 0.0_8
    END IF
    DO j=js,je+1
      DO k=1,km+1
        DO i=is,ie
          pe1_tl(i, k) = pe_tl(i, k, j)
          pe1(i, k) = pe(i, k, j)
        END DO
      END DO
      DO i=is,ie
        pe2_tl(i, 1) = 0.0_8
        pe2(i, 1) = ptop
        pe2_tl(i, km+1) = pe_tl(i, km+1, j)
        pe2(i, km+1) = pe(i, km+1, j)
      END DO
!(j < je+1)
      IF (j .LT. je + 1) THEN
! update ps
        DO i=is,ie
          ps(i, j) = pe1(i, km+1)
        END DO
        IF (hybrid_z) THEN
!--------------------------
! hybrid z_p coordinate
!--------------------------
          DO i=is,ie
            ze1_tl(i, km+1) = ze0_tl(i, j, km+1)
            ze1(i, km+1) = ze0(i, j, km+1)
          END DO
          DO k=km,1,-1
            DO i=is,ie
! current height
              ze1_tl(i, k) = ze1_tl(i, k+1) - delz_tl(i, j, k)
              ze1(i, k) = ze1(i, k+1) - delz(i, j, k)
            END DO
          END DO
!
! Copy ztop; the top layer must be thick enough to prevent numerical problems.
!
          DO i=is,ie
            ze2_tl(i, 1) = ze1_tl(i, 1)
            ze2(i, 1) = ze1(i, 1)
! Note: ze0 (top) updated
            ze0_tl(i, j, 1) = ze1_tl(i, 1)
            ze0(i, j, 1) = ze1(i, 1)
          END DO
          DO k=2,km+1
            DO i=is,ie
! specified height
              ze2_tl(i, k) = ze0_tl(i, j, k)
              ze2(i, k) = ze0(i, j, k)
            END DO
          END DO
! Check/fix monotonicity of top 3 layers thickness
!-----------------------------------------------
          dz1 = 1.
          DO i=is,ie
            DO k=1,3
!          if ( (ze2(i,k)-ze2(i,k+1)) < (ze2(i,k+1)-ze2(i,k+2)) ) then
!                ze2(i,k+1) = 0.5*(ze2(i,k) + ze2(i,k+2))
!          endif
              z_rat_tl = (ze2_tl(i, 1)-ze2_tl(i, 4))/(dz1(1)+dz1(2)+dz1(&
&               3))
              z_rat = (ze2(i, 1)-ze2(i, 4))/(dz1(1)+dz1(2)+dz1(3))
              ze2_tl(i, 3) = ze2_tl(i, 4) + dz1(3)*z_rat_tl
              ze2(i, 3) = ze2(i, 4) + dz1(3)*z_rat
              ze2_tl(i, 2) = ze2_tl(i, 3) + dz1(2)*z_rat_tl
              ze2(i, 2) = ze2(i, 3) + dz1(2)*z_rat
            END DO
          END DO
!-----------------------------------------------
          DO k=1,km
            DO i=is,ie
! density * grav
              deng_tl(i, k) = -((delp_tl(i, j, k)*delz(i, j, k)-delp(i, &
&               j, k)*delz_tl(i, j, k))/delz(i, j, k)**2)
              deng(i, k) = -(delp(i, j, k)/delz(i, j, k))
            END DO
          END DO
          IF (kord_tm .GE. 0.) THEN
            abs0 = kord_tm
          ELSE
            abs0 = -kord_tm
          END IF
          CALL REMAP_Z_NOLIMITERS_TLM(km, ze1, ze1_tl, ze2, ze2_tl, deng, deng_tl, &
&                    is, ie, 1, kord_tm_pert)
          CALL REMAP_Z(km, ze1, ze2, deng, is, ie, 1, abs0)
!-------------
! Update delz
!-------------
          DO k=1,km
            DO i=is,ie
              delz_tl(i, j, k) = ze2_tl(i, k+1) - ze2_tl(i, k)
              delz(i, j, k) = ze2(i, k+1) - ze2(i, k)
            END DO
          END DO
!------------
! update delp
!------------
          DO k=1,km-1
            DO i=is,ie
              dp2_tl(i, k) = -(delz_tl(i, j, k)*deng(i, k)+delz(i, j, k)&
&               *deng_tl(i, k))
              dp2(i, k) = -(delz(i, j, k)*deng(i, k))
              pe2_tl(i, k+1) = pe2_tl(i, k) + dp2_tl(i, k)
              pe2(i, k+1) = pe2(i, k) + dp2(i, k)
            END DO
          END DO
          DO i=is,ie
! to reduce rounding error
            dp2_tl(i, km) = pe2_tl(i, km+1) - pe2_tl(i, km)
            dp2(i, km) = pe2(i, km+1) - pe2(i, km)
          END DO
        ELSE
!
! Hybrid sigma-P coordinate:
!
          DO k=2,ks+1
            DO i=is,ie
              pe2_tl(i, k) = 0.0_8
              pe2(i, k) = ak(k)
            END DO
          END DO
          DO k=ks+2,km
            DO i=is,ie
              pe2_tl(i, k) = bk(k)*pe_tl(i, km+1, j)
              pe2(i, k) = ak(k) + bk(k)*pe(i, km+1, j)
            END DO
          END DO
          DO k=1,km
            DO i=is,ie
              dp2_tl(i, k) = pe2_tl(i, k+1) - pe2_tl(i, k)
              dp2(i, k) = pe2(i, k+1) - pe2(i, k)
            END DO
          END DO
        END IF
!------------
! update delp
!------------
        DO k=1,km
          DO i=is,ie
            delp_tl(i, j, k) = dp2_tl(i, k)
            delp(i, j, k) = dp2(i, k)
          END DO
        END DO
!----------------
! Map constituents
!----------------
        IF (nq .NE. 0) THEN
!------------------------------------------------------------------
! Do remapping one tracer at a time; seems to be faster on the SGI
! It requires less memory than mapn_ppm
!------------------------------------------------------------------
          DO iq=1,nq
            CALL MAP1_Q2_NOLIMITERS_TLM(km, pe1, pe1_tl, q(isd, jsd, 1, iq), q_tl(&
&                      isd, jsd, 1, iq), km, pe2, pe2_tl, q2, q2_tl, dp2&
&                      , dp2_tl, is, ie, 0, kord_tr_pert, j, isd, ied, jsd, &
&                      jed)
            CALL MAP1_Q2(km, pe1, q(isd, jsd, 1, iq), km, pe2, q2, dp2&
&                      , is, ie, 0, kord_tr, j, isd, ied, jsd, jed)
!           if (fill) call fillz(ie-is+1, km, 1, q2, dp2, q_source(is,j,iq))
            IF (fill) CALL FILLZ_TLM(ie - is + 1, km, 1, q2, q2_tl, dp2&
&                              , dp2_tl)
            DO k=1,km
              DO i=is,ie
                q_tl(i, j, k, iq) = q2_tl(i, k)
                q(i, j, k, iq) = q2(i, k)
              END DO
            END DO
          END DO
        END IF
!------------------
! Compute p**cappa
!------------------
        DO k=1,km+1
          DO i=is,ie
            pk1_tl(i, k) = pk_tl(i, j, k)
            pk1(i, k) = pk(i, j, k)
          END DO
        END DO
        DO i=is,ie
          pn1(i, :) = peln(i, :, j)
          pn2_tl(i, 1) = peln_tl(i, 1, j)
          pn2(i, 1) = peln(i, 1, j)
          pn2_tl(i, km+1) = peln_tl(i, km+1, j)
          pn2(i, km+1) = peln(i, km+1, j)
          pk2_tl(i, 1) = pk1_tl(i, 1)
          pk2(i, 1) = pk1(i, 1)
          pk2_tl(i, km+1) = pk1_tl(i, km+1)
          pk2(i, km+1) = pk1(i, km+1)
        END DO
        DO k=2,km
          DO i=is,ie
!        pk2(i,k) = pe2(i,k) ** akap
            pn2_tl(i, k) = pe2_tl(i, k)/pe2(i, k)
            pn2(i, k) = LOG(pe2(i, k))
            pk2_tl(i, k) = akap*pn2_tl(i, k)*EXP(akap*pn2(i, k))
            pk2(i, k) = EXP(akap*pn2(i, k))
          END DO
        END DO
        IF (te_map) THEN
!---------------------
! Compute Total Energy
!---------------------
          DO i=is,ie
            phis_tl(i, km+1) = 0.0_8
            phis(i, km+1) = hs(i, j)
          END DO
          DO k=km,1,-1
            DO i=is,ie
              phis_tl(i, k) = phis_tl(i, k+1) + pt_tl(i, j, k)*(pk1(i, k&
&               +1)-pk1(i, k)) + pt(i, j, k)*(pk1_tl(i, k+1)-pk1_tl(i, k&
&               ))
              phis(i, k) = phis(i, k+1) + pt(i, j, k)*(pk1(i, k+1)-pk1(i&
&               , k))
            END DO
          END DO
          DO k=1,km+1
            DO i=is,ie
              phis_tl(i, k) = phis_tl(i, k)*pe1(i, k) + phis(i, k)*&
&               pe1_tl(i, k)
              phis(i, k) = phis(i, k)*pe1(i, k)
            END DO
          END DO
          DO k=1,km
            DO i=is,ie
              te_tl(i, j, k) = te_tl(i, j, k) + ((phis_tl(i, k+1)-&
&               phis_tl(i, k))*(pe1(i, k+1)-pe1(i, k))-(phis(i, k+1)-&
&               phis(i, k))*(pe1_tl(i, k+1)-pe1_tl(i, k)))/(pe1(i, k+1)-&
&               pe1(i, k))**2
              te(i, j, k) = te(i, j, k) + (phis(i, k+1)-phis(i, k))/(pe1&
&               (i, k+1)-pe1(i, k))
            END DO
          END DO
!----------------
! Map Total Energy
!----------------
          SELECT CASE  (te_method) 
          CASE (1) 
            CALL MAP1_CUBIC_TE_TLM(km, pe1, pe1_tl, pe2, pe2_tl, te, &
&                            te_tl, is, ie, j, is, ie, js, je, 1, &
&                            kord_tm)
          CASE DEFAULT
            CALL MAP1_PPM_NOLIMITERS_TLM(km, pe1, pe1_tl, pe2, pe2_tl, te, te_tl, &
&                       is, ie, j, is, ie, js, je, 1, kord_tm_pert)
            CALL MAP1_PPM(km, pe1, pe2, te, is, ie, j, is, ie, js, je, 1, kord_tm)
          END SELECT
        ELSE IF (remap_t) THEN
!----------------------------------
! Map t using ze1 (hybrid_z) or logp 
!----------------------------------
          IF (hybrid_z) THEN
            DO k=1,km
              DO i=is,ie
                deng_tl(i, k) = pt_tl(i, j, k)
                deng(i, k) = pt(i, j, k)
              END DO
            END DO
            IF (kord_tm .GE. 0.) THEN
              abs1 = kord_tm
            ELSE
              abs1 = -kord_tm
            END IF
            CALL REMAP_Z_NOLIMITERS_TLM(km, ze1, ze1_tl, ze2, ze2_tl, deng, deng_tl&
&                      , is, ie, 2, kord_tm_pert)
            CALL REMAP_Z(km, ze1, ze2, deng, is, ie, 2, abs1)
            DO k=1,km
              DO i=is,ie
                pt_tl(i, j, k) = deng_tl(i, k)
                pt(i, j, k) = deng(i, k)
              END DO
            END DO
          ELSE
            IF (kord_tm .GE. 0.) THEN
              abs2 = kord_tm
            ELSE
              abs2 = -kord_tm
            END IF
            CALL MAP1_PPM_NOLIMITERS_TLM(km, peln(is, 1, j), peln_tl(is, 1, j), pn2&
&                       , pn2_tl, pt, pt_tl, is, ie, j, isd, ied, jsd, &
&                       jed, 2, kord_tm_pert)
            CALL MAP1_PPM(km, peln(is, 1, j), pn2, pt, is, ie, j, isd, ied, jsd, jed, 2, abs2)
          END IF
        ELSE
          IF (kord_tm .GE. 0.) THEN
            abs3 = kord_tm
          ELSE
            abs3 = -kord_tm
          END IF
!----------------
! Map pt using pk
!----------------
          CALL MAP1_PPM_NOLIMITERS_TLM(km, pk1, pk1_tl, pk2, pk2_tl, pt, pt_tl, is&
&                     , ie, j, isd, ied, jsd, jed, 1, kord_tm_pert)
          CALL MAP1_PPM(km, pk1, pk2, pt, is, ie, j, isd, ied, jsd, jed, 1, abs3)
        END IF
        IF (.NOT.hydrostatic) THEN
! Remap vertical wind:
          CALL MAP1_PPM_NOLIMITERS_TLM(km, pe1, pe1_tl, pe2, pe2_tl, w, w_tl, is, &
&                     ie, j, isd, ied, jsd, jed, -2, kord_mt)
          CALL MAP1_PPM(km, pe1, pe2, w, is, ie, j, isd, ied, jsd, jed, -2, kord_mt)
          IF (.NOT.hybrid_z) THEN
            IF (kord_tm .GE. 0.) THEN
              abs4 = kord_tm
            ELSE
              abs4 = -kord_tm
            END IF
! Remap delz for hybrid sigma-p coordinate
            CALL MAP1_PPM_NOLIMITERS_TLM(km, pe1, pe1_tl, pe2, pe2_tl, delz, &
&                       delz_tl, is, ie, j, is, ie, js, je, 1, kord_tm_pert)
            CALL MAP1_PPM(km, pe1, pe2, delz, is, ie, j, is, ie, js, je, 1, abs4)
            DO k=1,km
              DO i=is,ie
                delz_tl(i, j, k) = -(delz_tl(i, j, k)*dp2(i, k)+delz(i, &
&                 j, k)*dp2_tl(i, k))
                delz(i, j, k) = -(delz(i, j, k)*dp2(i, k))
              END DO
            END DO
          END IF
        END IF
!----------
! Update pk
!----------
        DO k=2,km
          DO i=is,ie
            pk_tl(i, j, k) = pk2_tl(i, k)
            pk(i, j, k) = pk2(i, k)
          END DO
        END DO
!----------------
        IF (do_omega) THEN
! Start do_omega
! Copy omega field to pe3
          DO i=is,ie
            pe3_tl(i, 1) = 0.0_8
            pe3(i, 1) = 0.
          END DO
          DO k=2,km+1
            DO i=is,ie
              pe3_tl(i, k) = omga_tl(i, j, k-1)
              pe3(i, k) = omga(i, j, k-1)
            END DO
          END DO
        END IF
        DO k=1,km+1
          DO i=is,ie
            pe0_tl(i, k) = peln_tl(i, k, j)
            pe0(i, k) = peln(i, k, j)
            peln_tl(i, k, j) = pn2_tl(i, k)
            peln(i, k, j) = pn2(i, k)
          END DO
        END DO
!------------
! Compute pkz
!------------
        IF (hydrostatic) THEN
          DO k=1,km
            DO i=is,ie
              pkz_tl(i, j, k) = ((pk2_tl(i, k+1)-pk2_tl(i, k))*akap*(&
&               peln(i, k+1, j)-peln(i, k, j))-(pk2(i, k+1)-pk2(i, k))*&
&               akap*(peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(akap*(peln(&
&               i, k+1, j)-peln(i, k, j)))**2
              pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1&
&               , j)-peln(i, k, j)))
            END DO
          END DO
        ELSE
          IF (ktop .GT. 1) THEN
            DO k=1,ktop-1
              DO i=is,ie
                pkz_tl(i, j, k) = ((pk2_tl(i, k+1)-pk2_tl(i, k))*akap*(&
&                 peln(i, k+1, j)-peln(i, k, j))-(pk2(i, k+1)-pk2(i, k))&
&                 *akap*(peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(akap*(&
&                 peln(i, k+1, j)-peln(i, k, j)))**2
                pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+&
&                 1, j)-peln(i, k, j)))
              END DO
            END DO
          END IF
          IF (remap_t) THEN
! Note: pt at this stage is T_v
            DO k=ktop,km
              DO i=is,ie
                arg1_tl = (rrg*delp_tl(i, j, k)*delz(i, j, k)-rrg*delp(i&
&                 , j, k)*delz_tl(i, j, k))*pt(i, j, k)/delz(i, j, k)**2&
&                 + rrg*delp(i, j, k)*pt_tl(i, j, k)/delz(i, j, k)
                arg1 = rrg*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
                arg2_tl = akap*arg1_tl/arg1
                arg2 = akap*LOG(arg1)
                pkz_tl(i, j, k) = arg2_tl*EXP(arg2)
                pkz(i, j, k) = EXP(arg2)
              END DO
            END DO
          ELSE
! Note: pt at this stage is cp*Theta_v
            DO k=ktop,km
              DO i=is,ie
                arg1_tl = (kapag*delp_tl(i, j, k)*delz(i, j, k)-kapag*&
&                 delp(i, j, k)*delz_tl(i, j, k))*pt(i, j, k)/delz(i, j&
&                 , k)**2 + kapag*delp(i, j, k)*pt_tl(i, j, k)/delz(i, j&
&                 , k)
                arg1 = kapag*delp(i, j, k)/delz(i, j, k)*pt(i, j, k)
                arg2_tl = k1k*arg1_tl/arg1
                arg2 = k1k*LOG(arg1)
                pkz_tl(i, j, k) = arg2_tl*EXP(arg2)
                pkz(i, j, k) = EXP(arg2)
              END DO
            END DO
          END IF
        END IF
! end do_omega
! Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
        IF (do_omega) THEN
          DO k=1,km
            DO i=is,ie
              dp2_tl(i, k) = 0.5*(peln_tl(i, k, j)+peln_tl(i, k+1, j))
              dp2(i, k) = 0.5*(peln(i, k, j)+peln(i, k+1, j))
            END DO
          END DO
          DO i=is,ie
            k_next = 1
            DO 110 n=1,km
              kp = k_next
              DO k=kp,km
                IF (dp2(i, n) .LE. pe0(i, k+1) .AND. dp2(i, n) .GE. pe0(&
&                   i, k)) GOTO 100
              END DO
              GOTO 110
 100          omga_tl(i, j, n) = pe3_tl(i, k) + (((pe3_tl(i, k+1)-pe3_tl&
&               (i, k))*(dp2(i, n)-pe0(i, k))+(pe3(i, k+1)-pe3(i, k))*(&
&               dp2_tl(i, n)-pe0_tl(i, k)))*(pe0(i, k+1)-pe0(i, k))-(pe3&
&               (i, k+1)-pe3(i, k))*(dp2(i, n)-pe0(i, k))*(pe0_tl(i, k+1&
&               )-pe0_tl(i, k)))/(pe0(i, k+1)-pe0(i, k))**2
              omga(i, j, n) = pe3(i, k) + (pe3(i, k+1)-pe3(i, k))*(dp2(i&
&               , n)-pe0(i, k))/(pe0(i, k+1)-pe0(i, k))
              k_next = k
 110        CONTINUE
          END DO
        END IF
      END IF
! end hybrid_z check
      IF (.NOT.hybrid_z) THEN
        DO i=is,ie+1
          pe0_tl(i, 1) = pe_tl(i, 1, j)
          pe0(i, 1) = pe(i, 1, j)
        END DO
!------
! map u
!------
        DO k=2,km+1
          DO i=is,ie
            pe0_tl(i, k) = 0.5*(pe_tl(i, k, j-1)+pe1_tl(i, k))
            pe0(i, k) = 0.5*(pe(i, k, j-1)+pe1(i, k))
          END DO
        END DO
        DO k=1,ks+1
          DO i=is,ie+1
            pe3_tl(i, k) = 0.0_8
            pe3(i, k) = ak(k)
          END DO
        END DO
        DO k=ks+2,km+1
          bkh = 0.5*bk(k)
          DO i=is,ie
            pe3_tl(i, k) = bkh*(pe_tl(i, km+1, j-1)+pe1_tl(i, km+1))
            pe3(i, k) = ak(k) + bkh*(pe(i, km+1, j-1)+pe1(i, km+1))
          END DO
        END DO
        CALL MAP1_PPM_NOLIMITERS_TLM(km, pe0(is:ie, :), pe0_tl(is:ie, :), pe3(is:ie&
&                   , :), pe3_tl(is:ie, :), u, u_tl, is, ie, j, isd, ied&
&                   , jsd, jed + 1, -1, kord_mt_pert)
        CALL MAP1_PPM(km, pe0(is:ie, :), pe3(is:ie, :), u, is, ie, j, isd, ied&
&                   , jsd, jed + 1, -1, kord_mt)
! (j < je+1)
        IF (j .LT. je + 1) THEN
!------
! map v
!------
          DO k=2,km+1
            DO i=is,ie+1
              pe0_tl(i, k) = 0.5*(pe_tl(i-1, k, j)+pe_tl(i, k, j))
              pe0(i, k) = 0.5*(pe(i-1, k, j)+pe(i, k, j))
            END DO
          END DO
          DO k=ks+2,km+1
            bkh = 0.5*bk(k)
            DO i=is,ie+1
              pe3_tl(i, k) = bkh*(pe_tl(i-1, km+1, j)+pe_tl(i, km+1, j))
              pe3(i, k) = ak(k) + bkh*(pe(i-1, km+1, j)+pe(i, km+1, j))
            END DO
          END DO
          CALL MAP1_PPM_NOLIMITERS_TLM(km, pe0, pe0_tl, pe3, pe3_tl, v, v_tl, is, &
&                     ie + 1, j, isd, ied + 1, jsd, jed, -1, kord_mt_pert)
          CALL MAP1_PPM(km, pe0, pe3, v, is, ie + 1, j, isd, ied + 1, jsd, jed, -1, kord_mt)
        END IF
      END IF
      DO k=1,km
        DO i=is,ie
          ua_tl(i, j, k) = pe2_tl(i, k+1)
          ua(i, j, k) = pe2(i, k+1)
        END DO
      END DO
    END DO
    IF (hybrid_z) THEN
!------- Hybrid_z section ---------------
! if ( square_domain ) then
!   call mpp_update_domains(ua , domain,  whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
! else
!      CALL MPP_UPDATE_DOMAINS_DUM3_TLM(ua, ua_tl, isd, ied, jsd, jed, km&
!&                               )
      call mpp_update_domains(ua    , domain, complete=.true.)
      call mpp_update_domains(ua_tl , domain, complete=.true.)
! endif
! u-wind
      DO j=js,je+1
        DO i=is,ie
          pe1_tl(i, 1) = 0.0_8
          pe1(i, 1) = ptop
          pe2_tl(i, 1) = 0.0_8
          pe2(i, 1) = ptop
        END DO
        DO k=2,km+1
          DO i=is,ie
            pe1_tl(i, k) = 0.5*(pe_tl(i, k, j-1)+pe_tl(i, k, j))
            pe1(i, k) = 0.5*(pe(i, k, j-1)+pe(i, k, j))
            pe2_tl(i, k) = 0.5*(ua_tl(i, j-1, k-1)+ua_tl(i, j, k-1))
            pe2(i, k) = 0.5*(ua(i, j-1, k-1)+ua(i, j, k-1))
          END DO
        END DO
        CALL MAP1_PPM_NOLIMITERS_TLM(km, pe1, pe1_tl, pe2, pe2_tl, u, u_tl, is, ie&
&                   , j, isd, ied, jsd, jed + 1, -1, kord_mt_pert)
        CALL MAP1_PPM(km, pe1, pe2, u, is, ie, j, isd, ied, jsd, jed + 1, -1, kord_mt)
      END DO
! v-wind
      DO j=js,je
        DO i=is,ie+1
          pe0_tl(i, 1) = 0.0_8
          pe0(i, 1) = ptop
          pe3_tl(i, 1) = 0.0_8
          pe3(i, 1) = ptop
        END DO
        DO k=2,km+1
          DO i=is,ie+1
            pe0_tl(i, k) = 0.5*(pe_tl(i-1, k, j)+pe_tl(i, k, j))
            pe0(i, k) = 0.5*(pe(i-1, k, j)+pe(i, k, j))
            pe3_tl(i, k) = 0.5*(ua_tl(i-1, j, k-1)+ua_tl(i, j, k-1))
            pe3(i, k) = 0.5*(ua(i-1, j, k-1)+ua(i, j, k-1))
          END DO
        END DO
        CALL MAP1_PPM_NOLIMITERS_TLM(km, pe0, pe0_tl, pe3, pe3_tl, v, v_tl, is, ie &
&                   + 1, j, isd, ied + 1, jsd, jed, -1, kord_mt_pert)
        CALL MAP1_PPM(km, pe0, pe3, v, is, ie + 1, j, isd, ied + 1, jsd, jed, -1, kord_mt)
      END DO
    END IF
!------------- Hybrid_z section ----------------------
    DO k=2,km
      DO j=js,je
        DO i=is,ie
          pe_tl(i, k, j) = ua_tl(i, j, k-1)
          pe(i, k, j) = ua(i, j, k-1)
        END DO
      END DO
    END DO
! end consv check
!  call cubed_to_latlon(u,  v, ua, va, dx, dy, rdxa, rdya, km, 1)
    IF (do_consv .AND. consv .GT. 0.) THEN
      IF (te_map) THEN
        te_2d_tl = 0.0_8
        DO j=js,je
          DO i=is,ie
            te_2d_tl(i, j) = te_tl(i, j, 1)*delp(i, j, 1) + te(i, j, 1)*&
&             delp_tl(i, j, 1)
            te_2d(i, j) = te(i, j, 1)*delp(i, j, 1)
          END DO
          DO k=2,km
            DO i=is,ie
              te_2d_tl(i, j) = te_2d_tl(i, j) + te_tl(i, j, k)*delp(i, j&
&               , k) + te(i, j, k)*delp_tl(i, j, k)
              te_2d(i, j) = te_2d(i, j) + te(i, j, k)*delp(i, j, k)
            END DO
          END DO
        END DO
        gz_tl = 0.0_8
        zsum0_tl = 0.0_8
        zsum1_tl = 0.0_8
      ELSE
        gz_tl = 0.0_8
        te_2d_tl = 0.0_8
        DO j=js,je
          IF (remap_t) THEN
            IF (hydrostatic) THEN
              DO i=is,ie
                gz_tl(i) = 0.0_8
                gz(i) = hs(i, j)
                DO k=1,km
                  gz_tl(i) = gz_tl(i) + rg*(pt_tl(i, j, k)*(peln(i, k+1&
&                   , j)-peln(i, k, j))+pt(i, j, k)*(peln_tl(i, k+1, j)-&
&                   peln_tl(i, k, j)))
                  gz(i) = gz(i) + rg*pt(i, j, k)*(peln(i, k+1, j)-peln(i&
&                   , k, j))
                END DO
              END DO
              DO i=is,ie
                te_2d_tl(i, j) = hs(i, j)*pe_tl(i, km+1, j) - pe_tl(i, 1&
&                 , j)*gz(i) - pe(i, 1, j)*gz_tl(i)
                te_2d(i, j) = pe(i, km+1, j)*hs(i, j) - pe(i, 1, j)*gz(i&
&                 )
              END DO
              DO k=1,km
                DO i=is,ie
                  te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(cp&
&                   *pt(i, j, k)+0.25*rsin2(i, j)*(u(i, j, k)**2+u(i, j+&
&                   1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u&
&                   (i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*cosa_s(i, j))&
&                   ) + delp(i, j, k)*(cp*pt_tl(i, j, k)+0.25*rsin2(i, j&
&                   )*(2*u(i, j, k)*u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i&
&                   , j+1, k)+2*v(i, j, k)*v_tl(i, j, k)+2*v(i+1, j, k)*&
&                   v_tl(i+1, j, k)-cosa_s(i, j)*((u_tl(i, j, k)+u_tl(i&
&                   , j+1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k)+u(i&
&                   , j+1, k))*(v_tl(i, j, k)+v_tl(i+1, j, k)))))
                  te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cp*pt(i, j&
&                   , k)+0.25*rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2&
&                   +v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1&
&                   , k))*(v(i, j, k)+v(i+1, j, k))*cosa_s(i, j)))
                END DO
              END DO
            ELSE
! non-hydrostatic & remap_t
              DO i=is,ie
                phis_tl(i, km+1) = 0.0_8
                phis(i, km+1) = hs(i, j)
                DO k=km,1,-1
                  phis_tl(i, k) = phis_tl(i, k+1) - grav*delz_tl(i, j, k&
&                   )
                  phis(i, k) = phis(i, k+1) - grav*delz(i, j, k)
                END DO
              END DO
              DO i=is,ie
                te_2d_tl(i, j) = 0.0_8
                te_2d(i, j) = 0.
              END DO
              DO k=1,km
                DO i=is,ie
! KE using 3D winds:
                  te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(cv&
&                   *pt(i, j, k)+0.5*(phis(i, k)+phis(i, k+1)+w(i, j, k)&
&                   **2+0.5*rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v&
&                   (i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k&
&                   ))*(v(i, j, k)+v(i+1, j, k))*cosa_s(i, j)))) + delp(&
&                   i, j, k)*(cv*pt_tl(i, j, k)+0.5*(phis_tl(i, k)+&
&                   phis_tl(i, k+1)+2*w(i, j, k)*w_tl(i, j, k)+0.5*rsin2&
&                   (i, j)*(2*u(i, j, k)*u_tl(i, j, k)+2*u(i, j+1, k)*&
&                   u_tl(i, j+1, k)+2*v(i, j, k)*v_tl(i, j, k)+2*v(i+1, &
&                   j, k)*v_tl(i+1, j, k)-cosa_s(i, j)*((u_tl(i, j, k)+&
&                   u_tl(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, &
&                   k)+u(i, j+1, k))*(v_tl(i, j, k)+v_tl(i+1, j, k))))))
                  te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cv*pt(i, j&
&                   , k)+0.5*(phis(i, k)+phis(i, k+1)+w(i, j, k)**2+0.5*&
&                   rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k&
&                   )**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i&
&                   , j, k)+v(i+1, j, k))*cosa_s(i, j))))
                END DO
              END DO
            END IF
          ELSE IF (hydrostatic) THEN
            DO i=is,ie
              gz_tl(i) = 0.0_8
              gz(i) = hs(i, j)
              DO k=1,km
                gz_tl(i) = gz_tl(i) + pt_tl(i, j, k)*(pk(i, j, k+1)-pk(i&
&                 , j, k)) + pt(i, j, k)*(pk_tl(i, j, k+1)-pk_tl(i, j, k&
&                 ))
                gz(i) = gz(i) + pt(i, j, k)*(pk(i, j, k+1)-pk(i, j, k))
              END DO
            END DO
            DO i=is,ie
              te_2d_tl(i, j) = hs(i, j)*pe_tl(i, km+1, j) - pe_tl(i, 1, &
&               j)*gz(i) - pe(i, 1, j)*gz_tl(i)
              te_2d(i, j) = pe(i, km+1, j)*hs(i, j) - pe(i, 1, j)*gz(i)
            END DO
            DO k=1,km
              DO i=is,ie
                te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(pt(i&
&                 , j, k)*pkz(i, j, k)+0.25*rsin2(i, j)*(u(i, j, k)**2+u&
&                 (i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, &
&                 k)+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*cosa_s(i, j&
&                 ))) + delp(i, j, k)*(pt_tl(i, j, k)*pkz(i, j, k)+pt(i&
&                 , j, k)*pkz_tl(i, j, k)+0.25*rsin2(i, j)*(2*u(i, j, k)&
&                 *u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i, j+1, k)+2*v(i, j&
&                 , k)*v_tl(i, j, k)+2*v(i+1, j, k)*v_tl(i+1, j, k)-&
&                 cosa_s(i, j)*((u_tl(i, j, k)+u_tl(i, j+1, k))*(v(i, j&
&                 , k)+v(i+1, j, k))+(u(i, j, k)+u(i, j+1, k))*(v_tl(i, &
&                 j, k)+v_tl(i+1, j, k)))))
                te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(pt(i, j, k)*&
&                 pkz(i, j, k)+0.25*rsin2(i, j)*(u(i, j, k)**2+u(i, j+1&
&                 , k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i&
&                 , j+1, k))*(v(i, j, k)+v(i+1, j, k))*cosa_s(i, j)))
              END DO
            END DO
          ELSE
!-----------------
! Non-hydrostatic:
!-----------------
            DO i=is,ie
              phis_tl(i, km+1) = 0.0_8
              phis(i, km+1) = hs(i, j)
              DO k=km,1,-1
                phis_tl(i, k) = phis_tl(i, k+1) - grav*delz_tl(i, j, k)
                phis(i, k) = phis(i, k+1) - grav*delz(i, j, k)
              END DO
            END DO
            DO i=is,ie
              te_2d_tl(i, j) = 0.0_8
              te_2d(i, j) = 0.
            END DO
            DO k=1,km
              DO i=is,ie
! KE using 3D winds:
                te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(&
&                 rgama*pt(i, j, k)*pkz(i, j, k)+0.5*(phis(i, k)+phis(i&
&                 , k+1)+w(i, j, k)**2+0.5*rsin2(i, j)*(u(i, j, k)**2+u(&
&                 i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k&
&                 )+u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*cosa_s(i, j)&
&                 ))) + delp(i, j, k)*(rgama*(pt_tl(i, j, k)*pkz(i, j, k&
&                 )+pt(i, j, k)*pkz_tl(i, j, k))+0.5*(phis_tl(i, k)+&
&                 phis_tl(i, k+1)+2*w(i, j, k)*w_tl(i, j, k)+0.5*rsin2(i&
&                 , j)*(2*u(i, j, k)*u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i&
&                 , j+1, k)+2*v(i, j, k)*v_tl(i, j, k)+2*v(i+1, j, k)*&
&                 v_tl(i+1, j, k)-cosa_s(i, j)*((u_tl(i, j, k)+u_tl(i, j&
&                 +1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k)+u(i, j+1&
&                 , k))*(v_tl(i, j, k)+v_tl(i+1, j, k))))))
                te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(rgama*pt(i, j&
&                 , k)*pkz(i, j, k)+0.5*(phis(i, k)+phis(i, k+1)+w(i, j&
&                 , k)**2+0.5*rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2&
&                 +v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k&
&                 ))*(v(i, j, k)+v(i+1, j, k))*cosa_s(i, j))))
              END DO
            END DO
          END IF
        END DO
        zsum0_tl = 0.0_8
        zsum1_tl = 0.0_8
      END IF
      DO j=js,je
        DO i=is,ie
          zsum1_tl(i, j) = pkz_tl(i, j, 1)*delp(i, j, 1) + pkz(i, j, 1)*&
&           delp_tl(i, j, 1)
          zsum1(i, j) = pkz(i, j, 1)*delp(i, j, 1)
        END DO
        DO k=2,km
          DO i=is,ie
            zsum1_tl(i, j) = zsum1_tl(i, j) + pkz_tl(i, j, k)*delp(i, j&
&             , k) + pkz(i, j, k)*delp_tl(i, j, k)
            zsum1(i, j) = zsum1(i, j) + pkz(i, j, k)*delp(i, j, k)
          END DO
        END DO
        DO i=is,ie
          zsum0_tl(i, j) = ptop*(pk_tl(i, j, 1)-pk_tl(i, j, km+1)) + &
&           zsum1_tl(i, j)
          zsum0(i, j) = ptop*(pk(i, j, 1)-pk(i, j, km+1)) + zsum1(i, j)
          te_2d_tl(i, j) = te0_2d_tl(i, j) - te_2d_tl(i, j)
          te_2d(i, j) = te0_2d(i, j) - te_2d(i, j)
        END DO
      END DO
      tpe_tl = 0.0
      tpe    = consv*g_sum(te_2d, is, ie, js, je, ng, area, 0)
      if ( hydrostatic ) then
           dtmp_tl = 0.0
           dtmp    = tpe / (cp*g_sum(zsum0,  is, ie, js, je, ng, area, 0))
      else
           dtmp_tl = 0.0
           dtmp    = tpe / (cv*g_sum(zsum1,  is, ie, js, je, ng, area, 0))
      endif

!-------------------------------------------------------------------------------
! One may use this quick fix to ensure reproducibility at the expense of a lower
! floating precision; this is fine for the TE correction
!-------------------------------------------------------------------------------
! convert to 4-byte real
      IF (reproduce_sum) dtmp = REAL(dtmp, 4)
    ELSE
      dtmp = 0.
      gz_tl = 0.0_8
      dtmp_tl = 0.0_8
    END IF
    IF (te_map) THEN
      DO j=js,je
        DO i=is,ie
          gz_tl(i) = 0.0_8
          gz(i) = hs(i, j)
        END DO
        DO k=km,1,-1
          DO i=is,ie
            tpe_tl = te_tl(i, j, k) - gz_tl(i) - 0.25*rsin2(i, j)*(2*u(i&
&             , j, k)*u_tl(i, j, k)+2*u(i, j+1, k)*u_tl(i, j+1, k)+2*v(i&
&             , j, k)*v_tl(i, j, k)+2*v(i+1, j, k)*v_tl(i+1, j, k)-&
&             cosa_s(i, j)*((u_tl(i, j, k)+u_tl(i, j+1, k))*(v(i, j, k)+&
&             v(i+1, j, k))+(u(i, j, k)+u(i, j+1, k))*(v_tl(i, j, k)+&
&             v_tl(i+1, j, k))))
            tpe = te(i, j, k) - gz(i) - 0.25*rsin2(i, j)*(u(i, j, k)**2+&
&             u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+&
&             u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*cosa_s(i, j))
            dlnp_tl = rg*(peln_tl(i, k+1, j)-peln_tl(i, k, j))
            dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
            tmp_tl = (tpe_tl*(cp-pe(i, k, j)*dlnp/delp(i, j, k))+tpe*((&
&             pe_tl(i, k, j)*dlnp+pe(i, k, j)*dlnp_tl)*delp(i, j, k)-pe(&
&             i, k, j)*dlnp*delp_tl(i, j, k))/delp(i, j, k)**2)/(cp-pe(i&
&             , k, j)*dlnp/delp(i, j, k))**2
            tmp = tpe/(cp-pe(i, k, j)*dlnp/delp(i, j, k))
            pt_tl(i, j, k) = cp*((tmp_tl*pkz(i, j, k)-tmp*pkz_tl(i, j, k&
&             ))/pkz(i, j, k)**2+dtmp_tl)
            pt(i, j, k) = cp*(tmp/pkz(i, j, k)+dtmp)
            gz_tl(i) = gz_tl(i) + dlnp_tl*tmp + dlnp*tmp_tl
            gz(i) = gz(i) + dlnp*tmp
          END DO
        END DO
      END DO
    ELSE IF (remap_t) THEN
! end k-loop
      DO k=1,km
        DO j=js,je
          DO i=is,ie
            pt_tl(i, j, k) = cp*((pt_tl(i, j, k)*pkz(i, j, k)-pt(i, j, k&
&             )*pkz_tl(i, j, k))/pkz(i, j, k)**2+dtmp_tl)
            pt(i, j, k) = cp*(pt(i, j, k)/pkz(i, j, k)+dtmp)
          END DO
        END DO
      END DO
    ELSE
      DO k=1,km
        DO j=js,je
          DO i=is,ie
            pt_tl(i, j, k) = pt_tl(i, j, k) + cp*dtmp_tl
            pt(i, j, k) = pt(i, j, k) + cp*dtmp
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE LAGRANGIAN_TO_EULERIAN_TLM

  SUBROUTINE COMPUTE_TOTAL_ENERGY_TLM(is, ie, js, je, isd, ied, jsd, jed&
&   , km, u, u_tl, v, v_tl, w, delz, delz_tl, pt, pt_tl, delp, delp_tl, &
&   q, q_tl, pe, pe_tl, peln, peln_tl, hs, grav, r_vir, cp, rg, hlv, &
&   te_2d, te_2d_tl, ua, va, teq, moist_phys, sphum, hydrostatic, id_te)
    IMPLICIT NONE
!     do j=js,je
!        do i=is,ie
!           teq(i,j) = teq(i,j) / (pe(i,km,j) - pe(i,1,j))
!        enddo
!     enddo
!------------------------------------------------------
! Compute vertically integrated total energy per column
!------------------------------------------------------
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: km, is, ie, js, je, isd, ied, jsd, jed, id_te&
&   , sphum
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: pt, delp
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: pt_tl, delp_tl
    REAL, DIMENSION(isd:ied, jsd:jed, km, sphum), INTENT(IN) :: q
    REAL, DIMENSION(isd:ied, jsd:jed, km, sphum), INTENT(IN) :: q_tl
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: u_tl(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_tl(isd:ied+1, jsd:jed, km)
! vertical velocity (m/s)
    REAL, INTENT(IN) :: w(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: delz(is:ie, js:je, km)
    REAL, INTENT(IN) :: delz_tl(is:ie, js:je, km)
! surface geopotential
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
! pressure at layer edges
    REAL(p_precision), INTENT(IN) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL(p_precision), INTENT(IN) :: pe_tl(is-1:ie+1, km+1, js-1:je+1)
! log(pe)
    REAL(p_precision), INTENT(IN) :: peln(is:ie, km+1, js:je)
    REAL(p_precision), INTENT(IN) :: peln_tl(is:ie, km+1, js:je)
    REAL, INTENT(IN) :: cp, rg, r_vir, hlv
    REAL, INTENT(IN) :: grav
    LOGICAL, INTENT(IN) :: moist_phys, hydrostatic
! Output:
! vertically integrated TE
    REAL, INTENT(OUT) :: te_2d(is:ie, js:je)
    REAL, INTENT(OUT) :: te_2d_tl(is:ie, js:je)
! Moist TE
    REAL, INTENT(OUT) :: teq(is:ie, js:je)
! Local
    REAL(p_precision), DIMENSION(is:ie, km) :: tv
    REAL(p_precision), DIMENSION(is:ie, km) :: tv_tl
    REAL(p_precision) :: phiz(is:ie, km+1)
    REAL(p_precision) :: phiz_tl(is:ie, km+1)
    REAL(p_precision) :: cv
    INTEGER :: i, j, k
    cv = cp - rg
    te_2d_tl = 0.0
    phiz_tl = 0.0_8
    tv_tl = 0.0_8
!----------------------
! Output lat-lon winds:
!----------------------
!  call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km)
    DO j=js,je
      IF (hydrostatic) THEN
        DO i=is,ie
          phiz_tl(i, km+1) = 0.0_8
          phiz(i, km+1) = hs(i, j)
        END DO
        DO k=km,1,-1
          DO i=is,ie
            tv_tl(i, k) = pt_tl(i, j, k)*(1.+r_vir*q(i, j, k, sphum)) + &
&             pt(i, j, k)*r_vir*q_tl(i, j, k, sphum)
            tv(i, k) = pt(i, j, k)*(1.+r_vir*q(i, j, k, sphum))
            phiz_tl(i, k) = phiz_tl(i, k+1) + rg*(tv_tl(i, k)*(peln(i, k&
&             +1, j)-peln(i, k, j))+tv(i, k)*(peln_tl(i, k+1, j)-peln_tl&
&             (i, k, j)))
            phiz(i, k) = phiz(i, k+1) + rg*tv(i, k)*(peln(i, k+1, j)-&
&             peln(i, k, j))
          END DO
        END DO
        DO i=is,ie
          te_2d_tl(i, j) = pe_tl(i, km+1, j)*phiz(i, km+1) + pe(i, km+1&
&           , j)*phiz_tl(i, km+1) - pe_tl(i, 1, j)*phiz(i, 1) - pe(i, 1&
&           , j)*phiz_tl(i, 1)
          te_2d(i, j) = pe(i, km+1, j)*phiz(i, km+1) - pe(i, 1, j)*phiz(&
&           i, 1)
        END DO
        DO k=1,km
          DO i=is,ie
            te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(cp*tv(i&
&             , k)+0.25*rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, &
&             j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j&
&             , k)+v(i+1, j, k))*cosa_s(i, j))) + delp(i, j, k)*(cp*&
&             tv_tl(i, k)+0.25*rsin2(i, j)*(2*u(i, j, k)*u_tl(i, j, k)+2&
&             *u(i, j+1, k)*u_tl(i, j+1, k)+2*v(i, j, k)*v_tl(i, j, k)+2&
&             *v(i+1, j, k)*v_tl(i+1, j, k)-cosa_s(i, j)*((u_tl(i, j, k)&
&             +u_tl(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k)+u(&
&             i, j+1, k))*(v_tl(i, j, k)+v_tl(i+1, j, k)))))
            te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cp*tv(i, k)+0.25*&
&             rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v&
&             (i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+v(i+1&
&             , j, k))*cosa_s(i, j)))
          END DO
        END DO
      ELSE
!-----------------
! Non-hydrostatic:
!-----------------
        DO i=is,ie
          phiz_tl(i, km+1) = 0.0_8
          phiz(i, km+1) = hs(i, j)
          DO k=km,1,-1
            phiz_tl(i, k) = phiz_tl(i, k+1) - grav*delz_tl(i, j, k)
            phiz(i, k) = phiz(i, k+1) - grav*delz(i, j, k)
          END DO
        END DO
        DO i=is,ie
          te_2d_tl(i, j) = 0.0
          te_2d(i, j) = 0.
        END DO
        DO k=1,km
          DO i=is,ie
!          te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv*pt(i,j,k)*(1.+r_vir*q(i,j,k,sphum)) +  &
!                       0.5*(phiz(i,k)+phiz(i,k+1)+ua(i,j,k)**2+va(i,j,k)**2+w(i,j,k)**2) )
            te_2d_tl(i, j) = te_2d_tl(i, j) + delp_tl(i, j, k)*(cv*pt(i&
&             , j, k)*(1.+r_vir*q(i, j, k, sphum))+0.5*(phiz(i, k)+phiz(&
&             i, k+1)+0.5*rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i&
&             , j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i&
&             , j, k)+v(i+1, j, k))*cosa_s(i, j)))) + delp(i, j, k)*(cv*&
&             (pt_tl(i, j, k)*(1.+r_vir*q(i, j, k, sphum))+pt(i, j, k)*&
&             r_vir*q_tl(i, j, k, sphum))+0.5*(phiz_tl(i, k)+phiz_tl(i, &
&             k+1)+0.5*rsin2(i, j)*(2*u(i, j, k)*u_tl(i, j, k)+2*u(i, j+&
&             1, k)*u_tl(i, j+1, k)+2*v(i, j, k)*v_tl(i, j, k)+2*v(i+1, &
&             j, k)*v_tl(i+1, j, k)-cosa_s(i, j)*((u_tl(i, j, k)+u_tl(i&
&             , j+1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k)+u(i, j+1&
&             , k))*(v_tl(i, j, k)+v_tl(i+1, j, k))))))
            te_2d(i, j) = te_2d(i, j) + delp(i, j, k)*(cv*pt(i, j, k)*(&
&             1.+r_vir*q(i, j, k, sphum))+0.5*(phiz(i, k)+phiz(i, k+1)+&
&             0.5*rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+v(i, j, k)&
&             **2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i, j, k)+&
&             v(i+1, j, k))*cosa_s(i, j))))
          END DO
        END DO
      END IF
    END DO
!-------------------------------------
! Doganostics computation for moist TE
!-------------------------------------
    IF (id_te .GT. 0) THEN
      DO j=js,je
        DO i=is,ie
          teq(i, j) = te_2d(i, j)
        END DO
      END DO
      IF (moist_phys) THEN
        DO k=1,km
          DO j=js,je
            DO i=is,ie
              teq(i, j) = teq(i, j) + hlv*q(i, j, k, sphum)*delp(i, j, k&
&               )
            END DO
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE COMPUTE_TOTAL_ENERGY_TLM

  SUBROUTINE PKEZ_TLM(km, ifirst, ilast, jfirst, jlast, pe, pk, pk_tl, &
&   akap, peln, peln_tl, pkz, pkz_tl)
    IMPLICIT NONE
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: km
! Latitude strip
    INTEGER, INTENT(IN) :: ifirst, ilast
! Latitude strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    REAL, INTENT(IN) :: akap
    REAL(p_precision), INTENT(IN) :: pe(ifirst-1:ilast+1, km+1, jfirst-1&
&   :jlast+1)
    REAL(p_precision), INTENT(IN) :: pk(ifirst:ilast, jfirst:jlast, km+1&
&   )
    REAL(p_precision), INTENT(IN) :: pk_tl(ifirst:ilast, jfirst:jlast, &
&   km+1)
! !OUTPUT
    REAL(p_precision), INTENT(OUT) :: pkz(ifirst:ilast, jfirst:jlast, km&
&   )
    REAL(p_precision), INTENT(OUT) :: pkz_tl(ifirst:ilast, jfirst:jlast&
&   , km)
! log (pe)
    REAL(p_precision), INTENT(INOUT) :: peln(ifirst:ilast, km+1, jfirst:&
&   jlast)
    REAL(p_precision), INTENT(INOUT) :: peln_tl(ifirst:ilast, km+1, &
&   jfirst:jlast)
! Local
    REAL(p_precision) :: pk2(ifirst:ilast, km+1)
    REAL(p_precision) :: pk2_tl(ifirst:ilast, km+1)
    REAL(p_precision) :: pek
    REAL(p_precision) :: pek_tl
    REAL(p_precision) :: lnp
    REAL(p_precision) :: ak1
    INTEGER :: i, j, k
    INTRINSIC LOG
    ak1 = (akap+1.)/akap
    pk2_tl = 0.0_8
    DO j=jfirst,jlast
      pek_tl = pk_tl(ifirst, j, 1)
      pek = pk(ifirst, j, 1)
      DO i=ifirst,ilast
        pk2_tl(i, 1) = pek_tl
        pk2(i, 1) = pek
      END DO
      DO k=2,km+1
        DO i=ifirst,ilast
!             peln(i,k,j) =  log(pe(i,k,j))
          pk2_tl(i, k) = pk_tl(i, j, k)
          pk2(i, k) = pk(i, j, k)
        END DO
      END DO
!---- GFDL modification
      IF (ptop .LT. ptop_min) THEN
        DO i=ifirst,ilast
          peln_tl(i, 1, j) = peln_tl(i, 2, j)
          peln(i, 1, j) = peln(i, 2, j) - ak1
        END DO
      ELSE
        lnp = LOG(ptop)
        DO i=ifirst,ilast
          peln_tl(i, 1, j) = 0.0_8
          peln(i, 1, j) = lnp
        END DO
      END IF
!---- GFDL modification
      DO k=1,km
        DO i=ifirst,ilast
          pkz_tl(i, j, k) = ((pk2_tl(i, k+1)-pk2_tl(i, k))*akap*(peln(i&
&           , k+1, j)-peln(i, k, j))-(pk2(i, k+1)-pk2(i, k))*akap*(&
&           peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(akap*(peln(i, k+1, j)&
&           -peln(i, k, j)))**2
          pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1, j)-&
&           peln(i, k, j)))
        END DO
      END DO
    END DO
  END SUBROUTINE PKEZ_TLM

  SUBROUTINE REMAP_Z_TLM(km, pe1, pe1_tl, pe2, pe2_tl, q2, q2_tl, i1, i2&
&   , iv, kord)
    IMPLICIT NONE
! !INPUT PARAMETERS:
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! Method order
    INTEGER, INTENT(IN) :: kord
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: iv
! height at layer edges 
    REAL(p_precision), INTENT(IN) :: pe1(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! hieght at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe2_tl(i1:i2, km+1)
! (from model top to bottom surface)
!      real(p_precision), intent(in) ::  q1(i1:i2,km)        ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL(p_precision), INTENT(INOUT) :: q2(i1:i2, km)
    REAL(p_precision), INTENT(INOUT) :: q2_tl(i1:i2, km)
! !LOCAL VARIABLES:
    REAL(p_precision) :: dp1(i1:i2, km)
    REAL(p_precision) :: dp1_tl(i1:i2, km)
    REAL(p_precision) :: q4(4, i1:i2, km)
    REAL(p_precision) :: q4_tl(4, i1:i2, km)
    REAL(p_precision) :: pl, pr, qsum, delp, esl
    REAL(p_precision) :: pl_tl, pr_tl, qsum_tl, delp_tl, esl_tl
    INTEGER :: i, k, l, m, k0
    dp1_tl = 0.0_8
    q4_tl = 0.0_8
    DO k=1,km
      DO i=i1,i2
! negative
        dp1_tl(i, k) = pe1_tl(i, k+1) - pe1_tl(i, k)
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4_tl(1, i, k) = q2_tl(i, k)
        q4(1, i, k) = q2(i, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv, kord)
      qsum_tl = 0.0_8
    ELSE
      CALL PPM_PROFILE_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv, kord)
      qsum_tl = 0.0_8
    END IF
! Mapping
    DO i=i1,i2
      k0 = 1
      DO 555 k=1,km
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .LE. pe1(i, l) .AND. pe2(i, k) .GE. pe1(i, l+1)&
&         ) GOTO 110
        END DO
        GOTO 123
 110    pl_tl = ((pe2_tl(i, k)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k)-pe1(i&
&         , l))*dp1_tl(i, l))/dp1(i, l)**2
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .GE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr_tl = ((pe2_tl(i, k+1)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k+1)-&
&           pe1(i, l))*dp1_tl(i, l))/dp1(i, l)**2
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          q2_tl(i, k) = q4_tl(2, i, l) + 0.5*((q4_tl(4, i, l)+q4_tl(3, i&
&           , l)-q4_tl(2, i, l))*(pr+pl)+(q4(4, i, l)+q4(3, i, l)-q4(2, &
&           i, l))*(pr_tl+pl_tl)) - r3*(q4_tl(4, i, l)*(pr*(pr+pl)+pl**2&
&           )+q4(4, i, l)*(pr_tl*(pr+pl)+pr*(pr_tl+pl_tl)+2*pl*pl_tl))
          q2(i, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i&
&           , l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
          k0 = l
          GOTO 555
        ELSE
! Fractional area...
          qsum_tl = (pe1_tl(i, l+1)-pe2_tl(i, k))*(q4(2, i, l)+0.5*(q4(4&
&           , i, l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.&
&           +pl*(1.+pl)))) + (pe1(i, l+1)-pe2(i, k))*(q4_tl(2, i, l)+0.5&
&           *((q4_tl(4, i, l)+q4_tl(3, i, l)-q4_tl(2, i, l))*(1.+pl)+(q4&
&           (4, i, l)+q4(3, i, l)-q4(2, i, l))*pl_tl)-r3*(q4_tl(4, i, l)&
&           *(1.+pl*(1.+pl))+q4(4, i, l)*(pl_tl*(1.+pl)+pl*pl_tl)))
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .LT. pe1(i, m+1)) THEN
! Whole layer..
              qsum_tl = qsum_tl + dp1_tl(i, m)*q4(1, i, m) + dp1(i, m)*&
&               q4_tl(1, i, m)
              qsum = qsum + dp1(i, m)*q4(1, i, m)
            ELSE
              GOTO 120
            END IF
          END DO
          GOTO 123
 120      delp_tl = pe2_tl(i, k+1) - pe1_tl(i, m)
          delp = pe2(i, k+1) - pe1(i, m)
          esl_tl = (delp_tl*dp1(i, m)-delp*dp1_tl(i, m))/dp1(i, m)**2
          esl = delp/dp1(i, m)
          qsum_tl = qsum_tl + delp_tl*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-&
&           q4(2, i, m)+q4(4, i, m)*(1.-r23*esl))) + delp*(q4_tl(2, i, m&
&           )+0.5*(esl_tl*(q4(3, i, m)-q4(2, i, m)+q4(4, i, m)*(1.-r23*&
&           esl))+esl*(q4_tl(3, i, m)-q4_tl(2, i, m)+q4_tl(4, i, m)*(1.-&
&           r23*esl)-q4(4, i, m)*r23*esl_tl)))
          qsum = qsum + delp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, &
&           m)+q4(4, i, m)*(1.-r23*esl)))
          k0 = m
        END IF
 123    q2_tl(i, k) = (qsum_tl*(pe2(i, k+1)-pe2(i, k))-qsum*(pe2_tl(i, k&
&         +1)-pe2_tl(i, k)))/(pe2(i, k+1)-pe2(i, k))**2
        q2(i, k) = qsum/(pe2(i, k+1)-pe2(i, k))
 555  CONTINUE
    END DO
  END SUBROUTINE REMAP_Z_TLM

  SUBROUTINE MAP1_PPM_TLM(km, pe1, pe1_tl, pe2, pe2_tl, q2, q2_tl, i1, &
&   i2, j, ibeg, iend, jbeg, jend, iv, kord)
    IMPLICIT NONE
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! Mode: 0 == constituents  1 == ???
    INTEGER, INTENT(IN) :: iv
!       2 == remap temp with cs scheme
! Method order
    INTEGER, INTENT(IN) :: kord
! Current latitude
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe1(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe2_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(INOUT) :: q2_tl(ibeg:iend, jbeg:jend, km)
! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
    REAL(p_precision) :: dp1(i1:i2, km)
    REAL(p_precision) :: dp1_tl(i1:i2, km)
    REAL(p_precision) :: q4(4, i1:i2, km)
    REAL(p_precision) :: q4_tl(4, i1:i2, km)
    REAL(p_precision) :: pl, pr, qsum, dp, esl
    REAL(p_precision) :: pl_tl, pr_tl, qsum_tl, dp_tl, esl_tl
    INTEGER :: i, k, l, m, k0
    dp1_tl = 0.0_8
    q4_tl = 0.0_8
    DO k=1,km
      DO i=i1,i2
        dp1_tl(i, k) = pe1_tl(i, k+1) - pe1_tl(i, k)
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4_tl(1, i, k) = q2_tl(i, j, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv, kord)
      qsum_tl = 0.0_8
    ELSE
      CALL PPM_PROFILE_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv, kord)
      qsum_tl = 0.0_8
    END IF
    DO i=i1,i2
      k0 = 1
      DO 555 k=1,km
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) GOTO 100
        END DO
        GOTO 123
 100    pl_tl = ((pe2_tl(i, k)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k)-pe1(i&
&         , l))*dp1_tl(i, l))/dp1(i, l)**2
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr_tl = ((pe2_tl(i, k+1)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k+1)-&
&           pe1(i, l))*dp1_tl(i, l))/dp1(i, l)**2
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          q2_tl(i, j, k) = q4_tl(2, i, l) + 0.5*((q4_tl(4, i, l)+q4_tl(3&
&           , i, l)-q4_tl(2, i, l))*(pr+pl)+(q4(4, i, l)+q4(3, i, l)-q4(&
&           2, i, l))*(pr_tl+pl_tl)) - r3*(q4_tl(4, i, l)*(pr*(pr+pl)+pl&
&           **2)+q4(4, i, l)*(pr_tl*(pr+pl)+pr*(pr_tl+pl_tl)+2*pl*pl_tl)&
&           )
          q2(i, j, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2&
&           , i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
          k0 = l
          GOTO 555
        ELSE
! Fractional area...
          qsum_tl = (pe1_tl(i, l+1)-pe2_tl(i, k))*(q4(2, i, l)+0.5*(q4(4&
&           , i, l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.&
&           +pl*(1.+pl)))) + (pe1(i, l+1)-pe2(i, k))*(q4_tl(2, i, l)+0.5&
&           *((q4_tl(4, i, l)+q4_tl(3, i, l)-q4_tl(2, i, l))*(1.+pl)+(q4&
&           (4, i, l)+q4(3, i, l)-q4(2, i, l))*pl_tl)-r3*(q4_tl(4, i, l)&
&           *(1.+pl*(1.+pl))+q4(4, i, l)*(pl_tl*(1.+pl)+pl*pl_tl)))
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
              qsum_tl = qsum_tl + dp1_tl(i, m)*q4(1, i, m) + dp1(i, m)*&
&               q4_tl(1, i, m)
              qsum = qsum + dp1(i, m)*q4(1, i, m)
            ELSE
              GOTO 110
            END IF
          END DO
          GOTO 123
 110      dp_tl = pe2_tl(i, k+1) - pe1_tl(i, m)
          dp = pe2(i, k+1) - pe1(i, m)
          esl_tl = (dp_tl*dp1(i, m)-dp*dp1_tl(i, m))/dp1(i, m)**2
          esl = dp/dp1(i, m)
          qsum_tl = qsum_tl + dp_tl*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4&
&           (2, i, m)+q4(4, i, m)*(1.-r23*esl))) + dp*(q4_tl(2, i, m)+&
&           0.5*(esl_tl*(q4(3, i, m)-q4(2, i, m)+q4(4, i, m)*(1.-r23*esl&
&           ))+esl*(q4_tl(3, i, m)-q4_tl(2, i, m)+q4_tl(4, i, m)*(1.-r23&
&           *esl)-q4(4, i, m)*r23*esl_tl)))
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
        END IF
 123    q2_tl(i, j, k) = (qsum_tl*(pe2(i, k+1)-pe2(i, k))-qsum*(pe2_tl(i&
&         , k+1)-pe2_tl(i, k)))/(pe2(i, k+1)-pe2(i, k))**2
        q2(i, j, k) = qsum/(pe2(i, k+1)-pe2(i, k))
 555  CONTINUE
    END DO
  END SUBROUTINE MAP1_PPM_TLM

  SUBROUTINE MAP1_Q2_TLM(km, pe1, pe1_tl, q1, q1_tl, kn, pe2, pe2_tl, q2&
&   , q2_tl, dp2, dp2_tl, i1, i2, iv, kord, j, ibeg, iend, jbeg, jend)
    IMPLICIT NONE
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: i1, i2
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Mode: 0 ==  constituents 1 == ???
    INTEGER, INTENT(IN) :: iv
    INTEGER, INTENT(IN) :: kord
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe1(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL(p_precision), INTENT(IN) :: pe2_tl(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! Field input
    REAL, INTENT(IN) :: q1(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(IN) :: q1_tl(ibeg:iend, jbeg:jend, km)
    REAL(p_precision), INTENT(IN) :: dp2(i1:i2, kn)
    REAL(p_precision), INTENT(IN) :: dp2_tl(i1:i2, kn)
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL(p_precision), INTENT(INOUT) :: q2(i1:i2, kn)
    REAL(p_precision), INTENT(INOUT) :: q2_tl(i1:i2, kn)
! !LOCAL VARIABLES:
    REAL(p_precision) :: dp1(i1:i2, km)
    REAL(p_precision) :: dp1_tl(i1:i2, km)
    REAL(p_precision) :: q4(4, i1:i2, km)
    REAL(p_precision) :: q4_tl(4, i1:i2, km)
    REAL(p_precision) :: pl, pr, qsum, dp, esl
    REAL(p_precision) :: pl_tl, pr_tl, qsum_tl, dp_tl, esl_tl
    INTEGER :: i, k, l, m, k0
    dp1_tl = 0.0_8
    q4_tl = 0.0_8
    DO k=1,km
      DO i=i1,i2
        dp1_tl(i, k) = pe1_tl(i, k+1) - pe1_tl(i, k)
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4_tl(1, i, k) = q1_tl(i, j, k)
        q4(1, i, k) = q1(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv, kord)
      qsum_tl = 0.0_8
    ELSE
      CALL PPM_PROFILE_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv, kord)
      qsum_tl = 0.0_8
    END IF
! Mapping
    DO i=i1,i2
      k0 = 1
      DO 555 k=1,kn
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) GOTO 110
        END DO
        GOTO 123
 110    pl_tl = ((pe2_tl(i, k)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k)-pe1(i&
&         , l))*dp1_tl(i, l))/dp1(i, l)**2
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr_tl = ((pe2_tl(i, k+1)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k+1)-&
&           pe1(i, l))*dp1_tl(i, l))/dp1(i, l)**2
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          q2_tl(i, k) = q4_tl(2, i, l) + 0.5*((q4_tl(4, i, l)+q4_tl(3, i&
&           , l)-q4_tl(2, i, l))*(pr+pl)+(q4(4, i, l)+q4(3, i, l)-q4(2, &
&           i, l))*(pr_tl+pl_tl)) - r3*(q4_tl(4, i, l)*(pr*(pr+pl)+pl**2&
&           )+q4(4, i, l)*(pr_tl*(pr+pl)+pr*(pr_tl+pl_tl)+2*pl*pl_tl))
          q2(i, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i&
&           , l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
          k0 = l
          GOTO 555
        ELSE
! Fractional area...
          qsum_tl = (pe1_tl(i, l+1)-pe2_tl(i, k))*(q4(2, i, l)+0.5*(q4(4&
&           , i, l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.&
&           +pl*(1.+pl)))) + (pe1(i, l+1)-pe2(i, k))*(q4_tl(2, i, l)+0.5&
&           *((q4_tl(4, i, l)+q4_tl(3, i, l)-q4_tl(2, i, l))*(1.+pl)+(q4&
&           (4, i, l)+q4(3, i, l)-q4(2, i, l))*pl_tl)-r3*(q4_tl(4, i, l)&
&           *(1.+pl*(1.+pl))+q4(4, i, l)*(pl_tl*(1.+pl)+pl*pl_tl)))
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
              qsum_tl = qsum_tl + dp1_tl(i, m)*q4(1, i, m) + dp1(i, m)*&
&               q4_tl(1, i, m)
              qsum = qsum + dp1(i, m)*q4(1, i, m)
            ELSE
              GOTO 120
            END IF
          END DO
          GOTO 123
 120      dp_tl = pe2_tl(i, k+1) - pe1_tl(i, m)
          dp = pe2(i, k+1) - pe1(i, m)
          esl_tl = (dp_tl*dp1(i, m)-dp*dp1_tl(i, m))/dp1(i, m)**2
          esl = dp/dp1(i, m)
          qsum_tl = qsum_tl + dp_tl*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4&
&           (2, i, m)+q4(4, i, m)*(1.-r23*esl))) + dp*(q4_tl(2, i, m)+&
&           0.5*(esl_tl*(q4(3, i, m)-q4(2, i, m)+q4(4, i, m)*(1.-r23*esl&
&           ))+esl*(q4_tl(3, i, m)-q4_tl(2, i, m)+q4_tl(4, i, m)*(1.-r23&
&           *esl)-q4(4, i, m)*r23*esl_tl)))
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
        END IF
 123    q2_tl(i, k) = (qsum_tl*dp2(i, k)-qsum*dp2_tl(i, k))/dp2(i, k)**2
        q2(i, k) = qsum/dp2(i, k)
 555  CONTINUE
    END DO
  END SUBROUTINE MAP1_Q2_TLM

  SUBROUTINE MAP1_CUBIC_TE_TLM(km, pe1, pe1_tl, pe2, pe2_tl, q2, q2_tl, &
&   i1, i2, j, ibeg, iend, jbeg, jend, iv, kord)
    IMPLICIT NONE
!EOC
! !INPUT PARAMETERS:
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! Mode: 0 ==  constituents  1 == ???
    INTEGER, INTENT(IN) :: iv
! Method order
    INTEGER, INTENT(IN) :: kord
! Current latitude
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe1(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe2_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the new vertical coordinate
!      real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(INOUT) :: q2_tl(ibeg:iend, jbeg:jend, km)
! !DESCRIPTION:
!
!     Perform Cubic Interpolation a given latitude
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY:
!    05.11.14   Takacs    Initial Code
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    REAL(p_precision) :: qx(i1:i2, km)
    REAL(p_precision) :: qx_tl(i1:i2, km)
    REAL(p_precision) :: logpl1(i1:i2, km)
    REAL(p_precision) :: logpl1_tl(i1:i2, km)
    REAL(p_precision) :: logpl2(i1:i2, km)
    REAL(p_precision) :: logpl2_tl(i1:i2, km)
    REAL(p_precision) :: dlogp1(i1:i2, km)
    REAL(p_precision) :: dlogp1_tl(i1:i2, km)
    REAL(p_precision) :: vsum1(i1:i2)
    REAL(p_precision) :: vsum1_tl(i1:i2)
    REAL(p_precision) :: vsum2(i1:i2)
    REAL(p_precision) :: vsum2_tl(i1:i2)
    REAL(p_precision) :: am2, am1, ap0, ap1, p, plp1, plp0, plm1, plm2, &
&   dlp0, dlm1, dlm2
    REAL(p_precision) :: am2_tl, am1_tl, ap0_tl, ap1_tl, p_tl, plp1_tl, &
&   plp0_tl, plm1_tl, plm2_tl, dlp0_tl, dlm1_tl, dlm2_tl
    INTEGER :: i, k, lm2, lm1, lp0, lp1
    INTRINSIC LOG
    INTRINSIC MAX
    INTRINSIC MIN
    qx_tl = 0.0_8
    logpl1_tl = 0.0_8
! Initialization
! --------------
    DO k=1,km
      qx_tl(:, k) = q2_tl(:, j, k)
      qx(:, k) = q2(:, j, k)
      logpl1_tl(:, k) = (pe1_tl(:, k)+pe1_tl(:, k+1))/(pe1(:, k)+pe1(:, &
&       k+1))
      logpl1(:, k) = LOG(d0_5*(pe1(:, k)+pe1(:, k+1)))
    END DO
    logpl2_tl = 0.0_8
    DO k=1,km
      logpl2_tl(:, k) = (pe2_tl(:, k)+pe2_tl(:, k+1))/(pe2(:, k)+pe2(:, &
&       k+1))
      logpl2(:, k) = LOG(d0_5*(pe2(:, k)+pe2(:, k+1)))
    END DO
    dlogp1_tl = 0.0_8
    DO k=1,km-1
      dlogp1_tl(:, k) = logpl1_tl(:, k+1) - logpl1_tl(:, k)
      dlogp1(:, k) = logpl1(:, k+1) - logpl1(:, k)
    END DO
! Compute vertical integral of Input TE
! -------------------------------------
    vsum1(:) = d0_0
    vsum1_tl = 0.0_8
    DO i=i1,i2
      DO k=1,km
        vsum1_tl(i) = vsum1_tl(i) + qx_tl(i, k)*(pe1(i, k+1)-pe1(i, k)) &
&         + qx(i, k)*(pe1_tl(i, k+1)-pe1_tl(i, k))
        vsum1(i) = vsum1(i) + qx(i, k)*(pe1(i, k+1)-pe1(i, k))
      END DO
      vsum1_tl(i) = (vsum1_tl(i)*(pe1(i, km+1)-pe1(i, 1))-vsum1(i)*(&
&       pe1_tl(i, km+1)-pe1_tl(i, 1)))/(pe1(i, km+1)-pe1(i, 1))**2
      vsum1(i) = vsum1(i)/(pe1(i, km+1)-pe1(i, 1))
    END DO
! Interpolate TE onto target Pressures
! ------------------------------------
    DO i=i1,i2
      DO k=1,km
        lm1 = 1
        lp0 = 1
        DO WHILE (lp0 .LE. km)
          IF (logpl1(i, lp0) .LT. logpl2(i, k)) THEN
            lp0 = lp0 + 1
          ELSE
            GOTO 100
          END IF
        END DO
 100    IF (lp0 - 1 .LT. 1) THEN
          lm1 = 1
        ELSE
          lm1 = lp0 - 1
        END IF
        IF (lp0 .GT. km) THEN
          lp0 = km
        ELSE
          lp0 = lp0
        END IF
! Extrapolate Linearly in LogP above first model level
! ----------------------------------------------------
        IF (lm1 .EQ. 1 .AND. lp0 .EQ. 1) THEN
          q2_tl(i, j, k) = qx_tl(i, 1) + (((qx_tl(i, 2)-qx_tl(i, 1))*(&
&           logpl2(i, k)-logpl1(i, 1))+(qx(i, 2)-qx(i, 1))*(logpl2_tl(i&
&           , k)-logpl1_tl(i, 1)))*(logpl1(i, 2)-logpl1(i, 1))-(qx(i, 2)&
&           -qx(i, 1))*(logpl2(i, k)-logpl1(i, 1))*(logpl1_tl(i, 2)-&
&           logpl1_tl(i, 1)))/(logpl1(i, 2)-logpl1(i, 1))**2
          q2(i, j, k) = qx(i, 1) + (qx(i, 2)-qx(i, 1))*(logpl2(i, k)-&
&           logpl1(i, 1))/(logpl1(i, 2)-logpl1(i, 1))
! Extrapolate Linearly in LogP below last model level
! ---------------------------------------------------
        ELSE IF (lm1 .EQ. km .AND. lp0 .EQ. km) THEN
          q2_tl(i, j, k) = qx_tl(i, km) + (((qx_tl(i, km)-qx_tl(i, km-1)&
&           )*(logpl2(i, k)-logpl1(i, km))+(qx(i, km)-qx(i, km-1))*(&
&           logpl2_tl(i, k)-logpl1_tl(i, km)))*(logpl1(i, km)-logpl1(i, &
&           km-1))-(qx(i, km)-qx(i, km-1))*(logpl2(i, k)-logpl1(i, km))*&
&           (logpl1_tl(i, km)-logpl1_tl(i, km-1)))/(logpl1(i, km)-logpl1&
&           (i, km-1))**2
          q2(i, j, k) = qx(i, km) + (qx(i, km)-qx(i, km-1))*(logpl2(i, k&
&           )-logpl1(i, km))/(logpl1(i, km)-logpl1(i, km-1))
! Interpolate Linearly in LogP between levels 1 => 2 and km-1 => km
! -----------------------------------------------------------------
        ELSE IF (lm1 .EQ. 1 .OR. lp0 .EQ. km) THEN
          q2_tl(i, j, k) = qx_tl(i, lp0) + (((qx_tl(i, lm1)-qx_tl(i, lp0&
&           ))*(logpl2(i, k)-logpl1(i, lp0))+(qx(i, lm1)-qx(i, lp0))*(&
&           logpl2_tl(i, k)-logpl1_tl(i, lp0)))*(logpl1(i, lm1)-logpl1(i&
&           , lp0))-(qx(i, lm1)-qx(i, lp0))*(logpl2(i, k)-logpl1(i, lp0)&
&           )*(logpl1_tl(i, lm1)-logpl1_tl(i, lp0)))/(logpl1(i, lm1)-&
&           logpl1(i, lp0))**2
          q2(i, j, k) = qx(i, lp0) + (qx(i, lm1)-qx(i, lp0))*(logpl2(i, &
&           k)-logpl1(i, lp0))/(logpl1(i, lm1)-logpl1(i, lp0))
! Interpolate Cubicly in LogP between other model levels
! ------------------------------------------------------
        ELSE
          lp1 = lp0 + 1
          lm2 = lm1 - 1
          p_tl = logpl2_tl(i, k)
          p = logpl2(i, k)
          plp1_tl = logpl1_tl(i, lp1)
          plp1 = logpl1(i, lp1)
          plp0_tl = logpl1_tl(i, lp0)
          plp0 = logpl1(i, lp0)
          plm1_tl = logpl1_tl(i, lm1)
          plm1 = logpl1(i, lm1)
          plm2_tl = logpl1_tl(i, lm2)
          plm2 = logpl1(i, lm2)
          dlp0_tl = dlogp1_tl(i, lp0)
          dlp0 = dlogp1(i, lp0)
          dlm1_tl = dlogp1_tl(i, lm1)
          dlm1 = dlogp1(i, lm1)
          dlm2_tl = dlogp1_tl(i, lm2)
          dlm2 = dlogp1(i, lm2)
          ap1_tl = ((((p_tl-plp0_tl)*(p-plm1)+(p-plp0)*(p_tl-plm1_tl))*(&
&           p-plm2)+(p-plp0)*(p-plm1)*(p_tl-plm2_tl))*dlp0*(dlp0+dlm1)*(&
&           dlp0+dlm1+dlm2)-(p-plp0)*(p-plm1)*(p-plm2)*((dlp0_tl*(dlp0+&
&           dlm1)+dlp0*(dlp0_tl+dlm1_tl))*(dlp0+dlm1+dlm2)+dlp0*(dlp0+&
&           dlm1)*(dlp0_tl+dlm1_tl+dlm2_tl)))/(dlp0*(dlp0+dlm1)*(dlp0+&
&           dlm1+dlm2))**2
          ap1 = (p-plp0)*(p-plm1)*(p-plm2)/(dlp0*(dlp0+dlm1)*(dlp0+dlm1+&
&           dlm2))
          ap0_tl = ((((plp1_tl-p_tl)*(p-plm1)+(plp1-p)*(p_tl-plm1_tl))*(&
&           p-plm2)+(plp1-p)*(p-plm1)*(p_tl-plm2_tl))*dlp0*dlm1*(dlm1+&
&           dlm2)-(plp1-p)*(p-plm1)*(p-plm2)*((dlp0_tl*dlm1+dlp0*dlm1_tl&
&           )*(dlm1+dlm2)+dlp0*dlm1*(dlm1_tl+dlm2_tl)))/(dlp0*dlm1*(dlm1&
&           +dlm2))**2
          ap0 = (plp1-p)*(p-plm1)*(p-plm2)/(dlp0*dlm1*(dlm1+dlm2))
          am1_tl = ((((plp1_tl-p_tl)*(plp0-p)+(plp1-p)*(plp0_tl-p_tl))*(&
&           p-plm2)+(plp1-p)*(plp0-p)*(p_tl-plm2_tl))*dlm1*dlm2*(dlp0+&
&           dlm1)-(plp1-p)*(plp0-p)*(p-plm2)*((dlm1_tl*dlm2+dlm1*dlm2_tl&
&           )*(dlp0+dlm1)+dlm1*dlm2*(dlp0_tl+dlm1_tl)))/(dlm1*dlm2*(dlp0&
&           +dlm1))**2
          am1 = (plp1-p)*(plp0-p)*(p-plm2)/(dlm1*dlm2*(dlp0+dlm1))
          am2_tl = ((((plp1_tl-p_tl)*(plp0-p)+(plp1-p)*(plp0_tl-p_tl))*(&
&           plm1-p)+(plp1-p)*(plp0-p)*(plm1_tl-p_tl))*dlm2*(dlm1+dlm2)*(&
&           dlp0+dlm1+dlm2)-(plp1-p)*(plp0-p)*(plm1-p)*((dlm2_tl*(dlm1+&
&           dlm2)+dlm2*(dlm1_tl+dlm2_tl))*(dlp0+dlm1+dlm2)+dlm2*(dlm1+&
&           dlm2)*(dlp0_tl+dlm1_tl+dlm2_tl)))/(dlm2*(dlm1+dlm2)*(dlp0+&
&           dlm1+dlm2))**2
          am2 = (plp1-p)*(plp0-p)*(plm1-p)/(dlm2*(dlm1+dlm2)*(dlp0+dlm1+&
&           dlm2))
          q2_tl(i, j, k) = ap1_tl*qx(i, lp1) + ap1*qx_tl(i, lp1) + &
&           ap0_tl*qx(i, lp0) + ap0*qx_tl(i, lp0) + am1_tl*qx(i, lm1) + &
&           am1*qx_tl(i, lm1) + am2_tl*qx(i, lm2) + am2*qx_tl(i, lm2)
          q2(i, j, k) = ap1*qx(i, lp1) + ap0*qx(i, lp0) + am1*qx(i, lm1)&
&           + am2*qx(i, lm2)
        END IF
      END DO
    END DO
! Compute vertical integral of Output TE
! --------------------------------------
    vsum2(:) = d0_0
    vsum2_tl = 0.0_8
    DO i=i1,i2
      DO k=1,km
        vsum2_tl(i) = vsum2_tl(i) + q2_tl(i, j, k)*(pe2(i, k+1)-pe2(i, k&
&         )) + q2(i, j, k)*(pe2_tl(i, k+1)-pe2_tl(i, k))
        vsum2(i) = vsum2(i) + q2(i, j, k)*(pe2(i, k+1)-pe2(i, k))
      END DO
      vsum2_tl(i) = (vsum2_tl(i)*(pe2(i, km+1)-pe2(i, 1))-vsum2(i)*(&
&       pe2_tl(i, km+1)-pe2_tl(i, 1)))/(pe2(i, km+1)-pe2(i, 1))**2
      vsum2(i) = vsum2(i)/(pe2(i, km+1)-pe2(i, 1))
    END DO
! Adjust Final TE to conserve
! ---------------------------
    DO i=i1,i2
      DO k=1,km
        q2_tl(i, j, k) = q2_tl(i, j, k) + vsum1_tl(i) - vsum2_tl(i)
        q2(i, j, k) = q2(i, j, k) + vsum1(i) - vsum2(i)
      END DO
    END DO
!        q2(i,j,k) = q2(i,j,k) * vsum1(i)/vsum2(i)
    RETURN
  END SUBROUTINE MAP1_CUBIC_TE_TLM

  SUBROUTINE CS_PROFILE_TLM(a4, a4_tl, delp, delp_tl, km, i1, i2, iv, &
&   kord)
    IMPLICIT NONE
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
    INTEGER, INTENT(IN) :: i1, i2
! vertical dimension
    INTEGER, INTENT(IN) :: km
! iv =-1: winds
    INTEGER, INTENT(IN) :: iv
! iv = 0: positive definite scalars
! iv = 1: others
! layer pressure thickness
    REAL(p_precision), INTENT(IN) :: delp(i1:i2, km)
    REAL(p_precision), INTENT(IN) :: delp_tl(i1:i2, km)
! Interpolated values
    REAL(p_precision), INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL(p_precision), INTENT(INOUT) :: a4_tl(4, i1:i2, km)
    INTEGER, INTENT(IN) :: kord
!-----------------------------------------------------------------------
    LOGICAL :: extm(i1:i2, km)
    REAL(p_precision) :: gam(i1:i2, km)
    REAL(p_precision) :: gam_tl(i1:i2, km)
    REAL(p_precision) :: q(i1:i2, km+1)
    REAL(p_precision) :: q_tl(i1:i2, km+1)
    REAL(p_precision) :: d4(i1:i2)
    REAL(p_precision) :: d4_tl(i1:i2)
    REAL(p_precision) :: bet, a_bot, grat, pmp, lac
    REAL(p_precision) :: bet_tl, a_bot_tl, grat_tl
    REAL(p_precision) :: pmp_1, lac_1, pmp_2, lac_2
    REAL(p_precision) :: pmp_1_tl, lac_1_tl, pmp_2_tl, lac_2_tl
    INTEGER :: i, k, im
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    REAL*8 :: y1_tl
    REAL*8 :: y18_tl
    REAL*8 :: y20_tl
    REAL*8 :: y2_tl
    REAL*8 :: y19_tl
    REAL*8 :: y29
    REAL*8 :: y28
    REAL*8 :: y27
    REAL*8 :: y26
    REAL*8 :: y25
    REAL*8 :: y24
    REAL*8 :: y21_tl
    REAL*8 :: y23
    REAL*8 :: y22
    REAL*8 :: y3_tl
    REAL*8 :: y21
    REAL*8 :: y20
    REAL*8 :: x1_tl
    REAL*8 :: x6
    REAL*8 :: y22_tl
    REAL*8 :: x5
    REAL*8 :: x4
    REAL*8 :: y4_tl
    REAL*8 :: x3
    REAL*8 :: y10_tl
    REAL*8 :: x2
    REAL*8 :: x1
    REAL*8 :: x2_tl
    REAL*8 :: y23_tl
    REAL*8 :: y5_tl
    REAL*8 :: y11_tl
    REAL*8 :: x3_tl
    REAL*8 :: y19
    REAL*8 :: y18
    REAL*8 :: y24_tl
    REAL*8 :: y17
    REAL*8 :: y16
    REAL*8 :: y6_tl
    REAL*8 :: y15
    REAL*8 :: y12_tl
    REAL*8 :: y14
    REAL*8 :: y13
    REAL*8 :: x4_tl
    REAL*8 :: y12
    REAL*8 :: y11
    REAL*8 :: y10
    REAL*8 :: y25_tl
    REAL*8 :: y7_tl
    REAL*8 :: y13_tl
    REAL*8 :: x5_tl
    REAL*8 :: y26_tl
    REAL*8 :: y8_tl
    REAL*8 :: y14_tl
    REAL*8 :: x6_tl
    REAL*8 :: y27_tl
    REAL*8 :: abs9
    REAL*8 :: y9_tl
    REAL*8 :: abs8
    REAL*8 :: y15_tl
    REAL*8 :: abs7
    INTEGER :: abs6
    REAL*8 :: abs5
    INTEGER :: abs4
    INTEGER :: abs3
    INTEGER :: abs2
    REAL*8 :: y28_tl
    INTEGER :: abs1
    INTEGER :: abs0
    REAL*8 :: y16_tl
    REAL*8 :: y29_tl
    REAL*8 :: y9
    REAL*8 :: y17_tl
    REAL*8 :: y8
    REAL*8 :: y7
    REAL*8 :: y6
    REAL*8 :: y5
    REAL*8 :: y4
    REAL*8 :: y3
    REAL*8 :: y2
    REAL*8 :: y1
    q_tl = 0.0_8
    gam_tl = 0.0_8
    DO i=i1,i2
! grid ratio
      grat_tl = (delp_tl(i, 2)*delp(i, 1)-delp(i, 2)*delp_tl(i, 1))/delp&
&       (i, 1)**2
      grat = delp(i, 2)/delp(i, 1)
      bet_tl = grat_tl*(grat+0.5) + grat*grat_tl
      bet = grat*(grat+0.5)
      q_tl(i, 1) = (((2*grat_tl*(grat+1.)+(grat+grat)*grat_tl)*a4(1, i, &
&       1)+(grat+grat)*(grat+1.)*a4_tl(1, i, 1)+a4_tl(1, i, 2))*bet-((&
&       grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))*bet_tl)/bet**2
      q(i, 1) = ((grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))/bet
      gam_tl(i, 1) = ((grat_tl*(grat+1.5)+grat*grat_tl)*bet-(1.+grat*(&
&       grat+1.5))*bet_tl)/bet**2
      gam(i, 1) = (1.+grat*(grat+1.5))/bet
    END DO
    d4_tl = 0.0_8
    DO k=2,km
      DO i=i1,i2
        d4_tl(i) = (delp_tl(i, k-1)*delp(i, k)-delp(i, k-1)*delp_tl(i, k&
&         ))/delp(i, k)**2
        d4(i) = delp(i, k-1)/delp(i, k)
        bet_tl = 2*d4_tl(i) - gam_tl(i, k-1)
        bet = 2. + d4(i) + d4(i) - gam(i, k-1)
        q_tl(i, k) = ((3.*(a4_tl(1, i, k-1)+d4_tl(i)*a4(1, i, k)+d4(i)*&
&         a4_tl(1, i, k))-q_tl(i, k-1))*bet-(3.*(a4(1, i, k-1)+d4(i)*a4(&
&         1, i, k))-q(i, k-1))*bet_tl)/bet**2
        q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
        gam_tl(i, k) = (d4_tl(i)*bet-d4(i)*bet_tl)/bet**2
        gam(i, k) = d4(i)/bet
      END DO
    END DO
    DO i=i1,i2
      a_bot_tl = d4_tl(i)*(d4(i)+1.5) + d4(i)*d4_tl(i)
      a_bot = 1. + d4(i)*(d4(i)+1.5)
      q_tl(i, km+1) = ((2.*((d4_tl(i)*(d4(i)+1.)+d4(i)*d4_tl(i))*a4(1, i&
&       , km)+d4(i)*(d4(i)+1.)*a4_tl(1, i, km))+a4_tl(1, i, km-1)-&
&       a_bot_tl*q(i, km)-a_bot*q_tl(i, km))*(d4(i)*(d4(i)+0.5)-a_bot*&
&       gam(i, km))-(2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&       a_bot*q(i, km))*(d4_tl(i)*(d4(i)+0.5)+d4(i)*d4_tl(i)-a_bot_tl*&
&       gam(i, km)-a_bot*gam_tl(i, km)))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i&
&       , km))**2
      q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&       a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
    END DO
    DO k=km,1,-1
      DO i=i1,i2
        q_tl(i, k) = q_tl(i, k) - gam_tl(i, k)*q(i, k+1) - gam(i, k)*&
&         q_tl(i, k+1)
        q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
      END DO
    END DO
!------------------
! Apply constraints
!------------------
    im = i2 - i1 + 1
! Apply *large-scale* constraints 
    DO i=i1,i2
      IF (a4(1, i, 1) .LT. a4(1, i, 2)) THEN
        y1_tl = a4_tl(1, i, 2)
        y1 = a4(1, i, 2)
      ELSE
        y1_tl = a4_tl(1, i, 1)
        y1 = a4(1, i, 1)
      END IF
      IF (q(i, 2) .GT. y1) THEN
        q_tl(i, 2) = y1_tl
        q(i, 2) = y1
      ELSE
        q(i, 2) = q(i, 2)
      END IF
      IF (a4(1, i, 1) .GT. a4(1, i, 2)) THEN
        y2_tl = a4_tl(1, i, 2)
        y2 = a4(1, i, 2)
      ELSE
        y2_tl = a4_tl(1, i, 1)
        y2 = a4(1, i, 1)
      END IF
      IF (q(i, 2) .LT. y2) THEN
        q_tl(i, 2) = y2_tl
        q(i, 2) = y2
      ELSE
        q(i, 2) = q(i, 2)
      END IF
    END DO
    DO k=2,km
      DO i=i1,i2
        gam_tl(i, k) = a4_tl(1, i, k) - a4_tl(1, i, k-1)
        gam(i, k) = a4(1, i, k) - a4(1, i, k-1)
      END DO
    END DO
    IF (kord .GE. 0.) THEN
      abs0 = kord
    ELSE
      abs0 = -kord
    END IF
! Interior:
    IF (abs0 .LT. 11) THEN
      DO k=3,km-1
        DO i=i1,i2
          IF (gam(i, k-1)*gam(i, k+1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y3_tl = a4_tl(1, i, k)
              y3 = a4(1, i, k)
            ELSE
              y3_tl = a4_tl(1, i, k-1)
              y3 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .GT. y3) THEN
              q_tl(i, k) = y3_tl
              q(i, k) = y3
            ELSE
              q(i, k) = q(i, k)
            END IF
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y4_tl = a4_tl(1, i, k)
              y4 = a4(1, i, k)
            ELSE
              y4_tl = a4_tl(1, i, k-1)
              y4 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .LT. y4) THEN
              q_tl(i, k) = y4_tl
              q(i, k) = y4
            ELSE
              q(i, k) = q(i, k)
            END IF
          ELSE IF (gam(i, k-1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y5_tl = a4_tl(1, i, k)
              y5 = a4(1, i, k)
            ELSE
              y5_tl = a4_tl(1, i, k-1)
              y5 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .LT. y5) THEN
              q_tl(i, k) = y5_tl
              q(i, k) = y5
            ELSE
              q(i, k) = q(i, k)
            END IF
          ELSE
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y6_tl = a4_tl(1, i, k)
              y6 = a4(1, i, k)
            ELSE
              y6_tl = a4_tl(1, i, k-1)
              y6 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .GT. y6) THEN
              q_tl(i, k) = y6_tl
              q(i, k) = y6
            ELSE
              q(i, k) = q(i, k)
            END IF
            IF (iv .EQ. 0) THEN
              IF (0. .LT. q(i, k)) THEN
                q(i, k) = q(i, k)
              ELSE
                q_tl(i, k) = 0.0_8
                q(i, k) = 0.
              END IF
            END IF
          END IF
        END DO
      END DO
    ELSE
! abs(kord) >=11
      DO k=3,km-1
        DO i=i1,i2
          IF (gam(i, k-1)*gam(i, k+1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y7_tl = a4_tl(1, i, k)
              y7 = a4(1, i, k)
            ELSE
              y7_tl = a4_tl(1, i, k-1)
              y7 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .GT. y7) THEN
              q_tl(i, k) = y7_tl
              q(i, k) = y7
            ELSE
              q(i, k) = q(i, k)
            END IF
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y8_tl = a4_tl(1, i, k)
              y8 = a4(1, i, k)
            ELSE
              y8_tl = a4_tl(1, i, k-1)
              y8 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .LT. y8) THEN
              q_tl(i, k) = y8_tl
              q(i, k) = y8
            ELSE
              q(i, k) = q(i, k)
            END IF
          ELSE IF (gam(i, k-1) .GT. 0.) THEN
            IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
              y9_tl = a4_tl(1, i, k)
              y9 = a4(1, i, k)
            ELSE
              y9_tl = a4_tl(1, i, k-1)
              y9 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .LT. y9) THEN
              q_tl(i, k) = y9_tl
              q(i, k) = y9
            ELSE
              q(i, k) = q(i, k)
            END IF
            IF (a4(1, i, k-1) + 0.5*gam(i, k-1) .GT. a4(1, i, k) - 0.5*&
&               gam(i, k+1)) THEN
              y10_tl = a4_tl(1, i, k) - 0.5*gam_tl(i, k+1)
              y10 = a4(1, i, k) - 0.5*gam(i, k+1)
            ELSE
              y10_tl = a4_tl(1, i, k-1) + 0.5*gam_tl(i, k-1)
              y10 = a4(1, i, k-1) + 0.5*gam(i, k-1)
            END IF
            IF (q(i, k) .GT. y10) THEN
              q_tl(i, k) = y10_tl
              q(i, k) = y10
            ELSE
              q(i, k) = q(i, k)
            END IF
          ELSE IF (iv .EQ. 0) THEN
! There exists a local min
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y11_tl = a4_tl(1, i, k)
              y11 = a4(1, i, k)
            ELSE
              y11_tl = a4_tl(1, i, k-1)
              y11 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .GT. y11) THEN
              q_tl(i, k) = y11_tl
              q(i, k) = y11
            ELSE
              q(i, k) = q(i, k)
            END IF
            IF (0. .GT. a4(1, i, k) - 0.5*gam(i, k+1)) THEN
              y23_tl = a4_tl(1, i, k) - 0.5*gam_tl(i, k+1)
              y23 = a4(1, i, k) - 0.5*gam(i, k+1)
            ELSE
              y23 = 0.
              y23_tl = 0.0_8
            END IF
            IF (a4(1, i, k-1) + 0.5*gam(i, k-1) .GT. y23) THEN
              y12_tl = y23_tl
              y12 = y23
            ELSE
              y12_tl = a4_tl(1, i, k-1) + 0.5*gam_tl(i, k-1)
              y12 = a4(1, i, k-1) + 0.5*gam(i, k-1)
            END IF
            IF (q(i, k) .LT. y12) THEN
              q_tl(i, k) = y12_tl
              q(i, k) = y12
            ELSE
              q(i, k) = q(i, k)
            END IF
          ELSE
            IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
              y13_tl = a4_tl(1, i, k)
              y13 = a4(1, i, k)
            ELSE
              y13_tl = a4_tl(1, i, k-1)
              y13 = a4(1, i, k-1)
            END IF
            IF (q(i, k) .GT. y13) THEN
              q_tl(i, k) = y13_tl
              q(i, k) = y13
            ELSE
              q(i, k) = q(i, k)
            END IF
            IF (a4(1, i, k-1) + 0.5*gam(i, k-1) .LT. a4(1, i, k) - 0.5*&
&               gam(i, k+1)) THEN
              y14_tl = a4_tl(1, i, k) - 0.5*gam_tl(i, k+1)
              y14 = a4(1, i, k) - 0.5*gam(i, k+1)
            ELSE
              y14_tl = a4_tl(1, i, k-1) + 0.5*gam_tl(i, k-1)
              y14 = a4(1, i, k-1) + 0.5*gam(i, k-1)
            END IF
            IF (q(i, k) .LT. y14) THEN
              q_tl(i, k) = y14_tl
              q(i, k) = y14
            ELSE
              q(i, k) = q(i, k)
            END IF
          END IF
        END DO
      END DO
    END IF
! Bottom:
    DO i=i1,i2
      IF (a4(1, i, km-1) .LT. a4(1, i, km)) THEN
        y15_tl = a4_tl(1, i, km)
        y15 = a4(1, i, km)
      ELSE
        y15_tl = a4_tl(1, i, km-1)
        y15 = a4(1, i, km-1)
      END IF
      IF (q(i, km) .GT. y15) THEN
        q_tl(i, km) = y15_tl
        q(i, km) = y15
      ELSE
        q(i, km) = q(i, km)
      END IF
      IF (a4(1, i, km-1) .GT. a4(1, i, km)) THEN
        y16_tl = a4_tl(1, i, km)
        y16 = a4(1, i, km)
      ELSE
        y16_tl = a4_tl(1, i, km-1)
        y16 = a4(1, i, km-1)
      END IF
      IF (q(i, km) .LT. y16) THEN
        q_tl(i, km) = y16_tl
        q(i, km) = y16
      ELSE
        q(i, km) = q(i, km)
      END IF
    END DO
    DO k=1,km
      DO i=i1,i2
        a4_tl(2, i, k) = q_tl(i, k)
        a4(2, i, k) = q(i, k)
        a4_tl(3, i, k) = q_tl(i, k+1)
        a4(3, i, k) = q(i, k+1)
      END DO
    END DO
    DO k=2,km-1
      DO i=i1,i2
        IF (gam(i, k)*gam(i, k+1) .GT. 0.0) THEN
          extm(i, k) = .false.
        ELSE
          extm(i, k) = .true.
        END IF
      END DO
    END DO
!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping
    IF (iv .EQ. 0) THEN
      DO i=i1,i2
        IF (0. .LT. a4(2, i, 1)) THEN
          a4(2, i, 1) = a4(2, i, 1)
        ELSE
          a4_tl(2, i, 1) = 0.0_8
          a4(2, i, 1) = 0.
        END IF
      END DO
    ELSE IF (iv .EQ. -1) THEN
      DO i=i1,i2
        IF (a4(2, i, 1)*a4(1, i, 1) .LE. 0.) THEN
          a4_tl(2, i, 1) = 0.0_8
          a4(2, i, 1) = 0.
        END IF
      END DO
    ELSE
      IF (iv .GE. 0.) THEN
        abs1 = iv
      ELSE
        abs1 = -iv
      END IF
      IF (abs1 .EQ. 2) THEN
        DO i=i1,i2
          a4_tl(2, i, 1) = a4_tl(1, i, 1)
          a4(2, i, 1) = a4(1, i, 1)
          a4_tl(3, i, 1) = a4_tl(1, i, 1)
          a4(3, i, 1) = a4(1, i, 1)
          a4_tl(4, i, 1) = 0.0_8
          a4(4, i, 1) = 0.
        END DO
      END IF
    END IF
    IF (iv .GE. 0.) THEN
      abs2 = iv
    ELSE
      abs2 = -iv
    END IF
    IF (abs2 .NE. 2) THEN
      DO i=i1,i2
        a4_tl(4, i, 1) = 3.*(2.*a4_tl(1, i, 1)-a4_tl(2, i, 1)-a4_tl(3, i&
&         , 1))
        a4(4, i, 1) = 3.*(2.*a4(1, i, 1)-(a4(2, i, 1)+a4(3, i, 1)))
      END DO
      CALL CS_LIMITERS_TLM(im, extm(i1, 1), a4(1, i1, 1), a4_tl(1, i1, 1&
&                    ), 1)
    END IF
! k=2
    DO i=i1,i2
      a4_tl(4, i, 2) = 3.*(2.*a4_tl(1, i, 2)-a4_tl(2, i, 2)-a4_tl(3, i, &
&       2))
      a4(4, i, 2) = 3.*(2.*a4(1, i, 2)-(a4(2, i, 2)+a4(3, i, 2)))
    END DO
    CALL CS_LIMITERS_TLM(im, extm(i1, 2), a4(1, i1, 2), a4_tl(1, i1, 2)&
&                  , 2)
!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
    DO k=3,km-2
      IF (kord .GE. 0.) THEN
        abs3 = kord
      ELSE
        abs3 = -kord
      END IF
      IF (abs3 .LT. 9) THEN
        DO i=i1,i2
! Left  edges
          pmp_1_tl = a4_tl(1, i, k) - 2.*gam_tl(i, k+1)
          pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
          lac_1_tl = pmp_1_tl + 1.5*gam_tl(i, k+2)
          lac_1 = pmp_1 + 1.5*gam(i, k+2)
          IF (a4(1, i, k) .GT. pmp_1) THEN
            IF (pmp_1 .GT. lac_1) THEN
              y24_tl = lac_1_tl
              y24 = lac_1
            ELSE
              y24_tl = pmp_1_tl
              y24 = pmp_1
            END IF
          ELSE IF (a4(1, i, k) .GT. lac_1) THEN
            y24_tl = lac_1_tl
            y24 = lac_1
          ELSE
            y24_tl = a4_tl(1, i, k)
            y24 = a4(1, i, k)
          END IF
          IF (a4(2, i, k) .LT. y24) THEN
            x1_tl = y24_tl
            x1 = y24
          ELSE
            x1_tl = a4_tl(2, i, k)
            x1 = a4(2, i, k)
          END IF
          IF (a4(1, i, k) .LT. pmp_1) THEN
            IF (pmp_1 .LT. lac_1) THEN
              y17_tl = lac_1_tl
              y17 = lac_1
            ELSE
              y17_tl = pmp_1_tl
              y17 = pmp_1
            END IF
          ELSE IF (a4(1, i, k) .LT. lac_1) THEN
            y17_tl = lac_1_tl
            y17 = lac_1
          ELSE
            y17_tl = a4_tl(1, i, k)
            y17 = a4(1, i, k)
          END IF
          IF (x1 .GT. y17) THEN
            a4_tl(2, i, k) = y17_tl
            a4(2, i, k) = y17
          ELSE
            a4_tl(2, i, k) = x1_tl
            a4(2, i, k) = x1
          END IF
! Right edges
          pmp_2_tl = a4_tl(1, i, k) + 2.*gam_tl(i, k)
          pmp_2 = a4(1, i, k) + 2.*gam(i, k)
          lac_2_tl = pmp_2_tl - 1.5*gam_tl(i, k-1)
          lac_2 = pmp_2 - 1.5*gam(i, k-1)
          IF (a4(1, i, k) .GT. pmp_2) THEN
            IF (pmp_2 .GT. lac_2) THEN
              y25_tl = lac_2_tl
              y25 = lac_2
            ELSE
              y25_tl = pmp_2_tl
              y25 = pmp_2
            END IF
          ELSE IF (a4(1, i, k) .GT. lac_2) THEN
            y25_tl = lac_2_tl
            y25 = lac_2
          ELSE
            y25_tl = a4_tl(1, i, k)
            y25 = a4(1, i, k)
          END IF
          IF (a4(3, i, k) .LT. y25) THEN
            x2_tl = y25_tl
            x2 = y25
          ELSE
            x2_tl = a4_tl(3, i, k)
            x2 = a4(3, i, k)
          END IF
          IF (a4(1, i, k) .LT. pmp_2) THEN
            IF (pmp_2 .LT. lac_2) THEN
              y18_tl = lac_2_tl
              y18 = lac_2
            ELSE
              y18_tl = pmp_2_tl
              y18 = pmp_2
            END IF
          ELSE IF (a4(1, i, k) .LT. lac_2) THEN
            y18_tl = lac_2_tl
            y18 = lac_2
          ELSE
            y18_tl = a4_tl(1, i, k)
            y18 = a4(1, i, k)
          END IF
          IF (x2 .GT. y18) THEN
            a4_tl(3, i, k) = y18_tl
            a4(3, i, k) = y18
          ELSE
            a4_tl(3, i, k) = x2_tl
            a4(3, i, k) = x2
          END IF
          a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(3&
&           , i, k))
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
      ELSE
        IF (kord .GE. 0.) THEN
          abs4 = kord
        ELSE
          abs4 = -kord
        END IF
        IF (abs4 .EQ. 9) THEN
          DO i=i1,i2
            IF (extm(i, k) .AND. (extm(i, k-1) .OR. extm(i, k+1))) THEN
! c90_mp122
! grid-scale 2-delta-z wave detected
              a4_tl(2, i, k) = a4_tl(1, i, k)
              a4(2, i, k) = a4(1, i, k)
              a4_tl(3, i, k) = a4_tl(1, i, k)
              a4(3, i, k) = a4(1, i, k)
              a4_tl(4, i, k) = 0.0_8
              a4(4, i, k) = 0.
            ELSE
              a4_tl(4, i, k) = 6.*a4_tl(1, i, k) - 3.*(a4_tl(2, i, k)+&
&               a4_tl(3, i, k))
              a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3, i, k)&
&               )
              IF (a4(4, i, k) .GE. 0.) THEN
                abs5 = a4(4, i, k)
              ELSE
                abs5 = -a4(4, i, k)
              END IF
              IF (a4(2, i, k) - a4(3, i, k) .GE. 0.) THEN
                abs8 = a4(2, i, k) - a4(3, i, k)
              ELSE
                abs8 = -(a4(2, i, k)-a4(3, i, k))
              END IF
! Check within the smooth region if subgrid profile is non-monotonic
              IF (abs5 .GT. abs8) THEN
                pmp_1_tl = a4_tl(1, i, k) - 2.*gam_tl(i, k+1)
                pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                lac_1_tl = pmp_1_tl + 1.5*gam_tl(i, k+2)
                lac_1 = pmp_1 + 1.5*gam(i, k+2)
                IF (a4(1, i, k) .GT. pmp_1) THEN
                  IF (pmp_1 .GT. lac_1) THEN
                    y26_tl = lac_1_tl
                    y26 = lac_1
                  ELSE
                    y26_tl = pmp_1_tl
                    y26 = pmp_1
                  END IF
                ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                  y26_tl = lac_1_tl
                  y26 = lac_1
                ELSE
                  y26_tl = a4_tl(1, i, k)
                  y26 = a4(1, i, k)
                END IF
                IF (a4(2, i, k) .LT. y26) THEN
                  x3_tl = y26_tl
                  x3 = y26
                ELSE
                  x3_tl = a4_tl(2, i, k)
                  x3 = a4(2, i, k)
                END IF
                IF (a4(1, i, k) .LT. pmp_1) THEN
                  IF (pmp_1 .LT. lac_1) THEN
                    y19_tl = lac_1_tl
                    y19 = lac_1
                  ELSE
                    y19_tl = pmp_1_tl
                    y19 = pmp_1
                  END IF
                ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                  y19_tl = lac_1_tl
                  y19 = lac_1
                ELSE
                  y19_tl = a4_tl(1, i, k)
                  y19 = a4(1, i, k)
                END IF
                IF (x3 .GT. y19) THEN
                  a4_tl(2, i, k) = y19_tl
                  a4(2, i, k) = y19
                ELSE
                  a4_tl(2, i, k) = x3_tl
                  a4(2, i, k) = x3
                END IF
                pmp_2_tl = a4_tl(1, i, k) + 2.*gam_tl(i, k)
                pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                lac_2_tl = pmp_2_tl - 1.5*gam_tl(i, k-1)
                lac_2 = pmp_2 - 1.5*gam(i, k-1)
                IF (a4(1, i, k) .GT. pmp_2) THEN
                  IF (pmp_2 .GT. lac_2) THEN
                    y27_tl = lac_2_tl
                    y27 = lac_2
                  ELSE
                    y27_tl = pmp_2_tl
                    y27 = pmp_2
                  END IF
                ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                  y27_tl = lac_2_tl
                  y27 = lac_2
                ELSE
                  y27_tl = a4_tl(1, i, k)
                  y27 = a4(1, i, k)
                END IF
                IF (a4(3, i, k) .LT. y27) THEN
                  x4_tl = y27_tl
                  x4 = y27
                ELSE
                  x4_tl = a4_tl(3, i, k)
                  x4 = a4(3, i, k)
                END IF
                IF (a4(1, i, k) .LT. pmp_2) THEN
                  IF (pmp_2 .LT. lac_2) THEN
                    y20_tl = lac_2_tl
                    y20 = lac_2
                  ELSE
                    y20_tl = pmp_2_tl
                    y20 = pmp_2
                  END IF
                ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                  y20_tl = lac_2_tl
                  y20 = lac_2
                ELSE
                  y20_tl = a4_tl(1, i, k)
                  y20 = a4(1, i, k)
                END IF
                IF (x4 .GT. y20) THEN
                  a4_tl(3, i, k) = y20_tl
                  a4(3, i, k) = y20
                ELSE
                  a4_tl(3, i, k) = x4_tl
                  a4(3, i, k) = x4
                END IF
                a4_tl(4, i, k) = 6.*a4_tl(1, i, k) - 3.*(a4_tl(2, i, k)+&
&                 a4_tl(3, i, k))
                a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3, i, &
&                 k))
              END IF
            END IF
          END DO
        ELSE
          IF (kord .GE. 0.) THEN
            abs6 = kord
          ELSE
            abs6 = -kord
          END IF
          IF (abs6 .EQ. 10) THEN
            DO i=i1,i2
              IF (extm(i, k)) THEN
                IF (extm(i, k-1) .OR. extm(i, k+1)) THEN
! grid-scale 2-delta-z wave detected
                  a4_tl(2, i, k) = a4_tl(1, i, k)
                  a4(2, i, k) = a4(1, i, k)
                  a4_tl(3, i, k) = a4_tl(1, i, k)
                  a4(3, i, k) = a4(1, i, k)
                  a4_tl(4, i, k) = 0.0_8
                  a4(4, i, k) = 0.
                ELSE
! True local extremum
                  a4_tl(4, i, k) = 6.*a4_tl(1, i, k) - 3.*(a4_tl(2, i, k&
&                   )+a4_tl(3, i, k))
                  a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3, i&
&                   , k))
                END IF
              ELSE
! not a local extremum
                a4_tl(4, i, k) = 6.*a4_tl(1, i, k) - 3.*(a4_tl(2, i, k)+&
&                 a4_tl(3, i, k))
                a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3, i, &
&                 k))
                IF (a4(4, i, k) .GE. 0.) THEN
                  abs7 = a4(4, i, k)
                ELSE
                  abs7 = -a4(4, i, k)
                END IF
                IF (a4(2, i, k) - a4(3, i, k) .GE. 0.) THEN
                  abs9 = a4(2, i, k) - a4(3, i, k)
                ELSE
                  abs9 = -(a4(2, i, k)-a4(3, i, k))
                END IF
! Check within the smooth region if subgrid profile is non-monotonic
                IF (abs7 .GT. abs9) THEN
                  pmp_1_tl = a4_tl(1, i, k) - 2.*gam_tl(i, k+1)
                  pmp_1 = a4(1, i, k) - 2.*gam(i, k+1)
                  lac_1_tl = pmp_1_tl + 1.5*gam_tl(i, k+2)
                  lac_1 = pmp_1 + 1.5*gam(i, k+2)
                  IF (a4(1, i, k) .GT. pmp_1) THEN
                    IF (pmp_1 .GT. lac_1) THEN
                      y28_tl = lac_1_tl
                      y28 = lac_1
                    ELSE
                      y28_tl = pmp_1_tl
                      y28 = pmp_1
                    END IF
                  ELSE IF (a4(1, i, k) .GT. lac_1) THEN
                    y28_tl = lac_1_tl
                    y28 = lac_1
                  ELSE
                    y28_tl = a4_tl(1, i, k)
                    y28 = a4(1, i, k)
                  END IF
                  IF (a4(2, i, k) .LT. y28) THEN
                    x5_tl = y28_tl
                    x5 = y28
                  ELSE
                    x5_tl = a4_tl(2, i, k)
                    x5 = a4(2, i, k)
                  END IF
                  IF (a4(1, i, k) .LT. pmp_1) THEN
                    IF (pmp_1 .LT. lac_1) THEN
                      y21_tl = lac_1_tl
                      y21 = lac_1
                    ELSE
                      y21_tl = pmp_1_tl
                      y21 = pmp_1
                    END IF
                  ELSE IF (a4(1, i, k) .LT. lac_1) THEN
                    y21_tl = lac_1_tl
                    y21 = lac_1
                  ELSE
                    y21_tl = a4_tl(1, i, k)
                    y21 = a4(1, i, k)
                  END IF
                  IF (x5 .GT. y21) THEN
                    a4_tl(2, i, k) = y21_tl
                    a4(2, i, k) = y21
                  ELSE
                    a4_tl(2, i, k) = x5_tl
                    a4(2, i, k) = x5
                  END IF
                  pmp_2_tl = a4_tl(1, i, k) + 2.*gam_tl(i, k)
                  pmp_2 = a4(1, i, k) + 2.*gam(i, k)
                  lac_2_tl = pmp_2_tl - 1.5*gam_tl(i, k-1)
                  lac_2 = pmp_2 - 1.5*gam(i, k-1)
                  IF (a4(1, i, k) .GT. pmp_2) THEN
                    IF (pmp_2 .GT. lac_2) THEN
                      y29_tl = lac_2_tl
                      y29 = lac_2
                    ELSE
                      y29_tl = pmp_2_tl
                      y29 = pmp_2
                    END IF
                  ELSE IF (a4(1, i, k) .GT. lac_2) THEN
                    y29_tl = lac_2_tl
                    y29 = lac_2
                  ELSE
                    y29_tl = a4_tl(1, i, k)
                    y29 = a4(1, i, k)
                  END IF
                  IF (a4(3, i, k) .LT. y29) THEN
                    x6_tl = y29_tl
                    x6 = y29
                  ELSE
                    x6_tl = a4_tl(3, i, k)
                    x6 = a4(3, i, k)
                  END IF
                  IF (a4(1, i, k) .LT. pmp_2) THEN
                    IF (pmp_2 .LT. lac_2) THEN
                      y22_tl = lac_2_tl
                      y22 = lac_2
                    ELSE
                      y22_tl = pmp_2_tl
                      y22 = pmp_2
                    END IF
                  ELSE IF (a4(1, i, k) .LT. lac_2) THEN
                    y22_tl = lac_2_tl
                    y22 = lac_2
                  ELSE
                    y22_tl = a4_tl(1, i, k)
                    y22 = a4(1, i, k)
                  END IF
                  IF (x6 .GT. y22) THEN
                    a4_tl(3, i, k) = y22_tl
                    a4(3, i, k) = y22
                  ELSE
                    a4_tl(3, i, k) = x6_tl
                    a4(3, i, k) = x6
                  END IF
                  a4_tl(4, i, k) = 6.*a4_tl(1, i, k) - 3.*(a4_tl(2, i, k&
&                   )+a4_tl(3, i, k))
                  a4(4, i, k) = 6.*a4(1, i, k) - 3.*(a4(2, i, k)+a4(3, i&
&                   , k))
                END IF
              END IF
            END DO
          ELSE
! kord = 11, ...
            DO i=i1,i2
              IF (extm(i, k) .AND. (extm(i, k-1) .OR. extm(i, k+1))) &
&             THEN
! Noisy region:
                a4_tl(2, i, k) = a4_tl(1, i, k)
                a4(2, i, k) = a4(1, i, k)
                a4_tl(3, i, k) = a4_tl(1, i, k)
                a4(3, i, k) = a4(1, i, k)
                a4_tl(4, i, k) = 0.0_8
                a4(4, i, k) = 0.
              ELSE
                a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-&
&                 a4_tl(3, i, k))
                a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k&
&                 )))
              END IF
            END DO
          END IF
        END IF
      END IF
! Additional constraint to ensure positivity
      IF (iv .EQ. 0) CALL CS_LIMITERS_TLM(im, extm(i1, k), a4(1, i1, k)&
&                                   , a4_tl(1, i1, k), 0)
    END DO
! k-loop
!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
    IF (iv .EQ. 0) THEN
      DO i=i1,i2
        IF (0. .LT. a4(3, i, km)) THEN
          a4(3, i, km) = a4(3, i, km)
        ELSE
          a4_tl(3, i, km) = 0.0_8
          a4(3, i, km) = 0.
        END IF
      END DO
    ELSE IF (iv .LT. 0) THEN
      DO i=i1,i2
        IF (a4(3, i, km)*a4(1, i, km) .LE. 0.) THEN
          a4_tl(3, i, km) = 0.0_8
          a4(3, i, km) = 0.
        END IF
      END DO
    END IF
    DO k=km-1,km
      DO i=i1,i2
        a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(3, i&
&         , k))
        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
      END DO
      IF (k .EQ. km - 1) CALL CS_LIMITERS_TLM(im, extm(i1, k), a4(1, i1&
&                                       , k), a4_tl(1, i1, k), 2)
      IF (k .EQ. km) CALL CS_LIMITERS_TLM(im, extm(i1, k), a4(1, i1, k)&
&                                   , a4_tl(1, i1, k), 1)
    END DO
  END SUBROUTINE CS_PROFILE_TLM

  SUBROUTINE CS_LIMITERS_TLM(im, extm, a4, a4_tl, iv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    INTEGER, INTENT(IN) :: iv
    LOGICAL, INTENT(IN) :: extm(im)
! PPM array
    REAL(p_precision), INTENT(INOUT) :: a4(4, im)
    REAL(p_precision), INTENT(INOUT) :: a4_tl(4, im)
! !LOCAL VARIABLES:
    REAL(p_precision) :: da1, da2, a6da
    REAL(p_precision) :: fmin
    INTEGER :: i
    INTRINSIC ABS
    REAL*8 :: abs0
    IF (iv .EQ. 0) THEN
! Positive definite constraint
      DO i=1,im
        IF (a4(1, i) .LE. 0.) THEN
          a4_tl(2, i) = a4_tl(1, i)
          a4(2, i) = a4(1, i)
          a4_tl(3, i) = a4_tl(1, i)
          a4(3, i) = a4(1, i)
          a4_tl(4, i) = 0.0_8
          a4(4, i) = 0.
        ELSE
          IF (a4(3, i) - a4(2, i) .GE. 0.) THEN
            abs0 = a4(3, i) - a4(2, i)
          ELSE
            abs0 = -(a4(3, i)-a4(2, i))
          END IF
          IF (abs0 .LT. -a4(4, i)) THEN
            IF (a4(1, i) + 0.25*(a4(3, i)-a4(2, i))**2/a4(4, i) + a4(4, &
&               i)*r12 .LT. 0.) THEN
! local minimum is negative
              IF (a4(1, i) .LT. a4(3, i) .AND. a4(1, i) .LT. a4(2, i)) &
&             THEN
                a4_tl(3, i) = a4_tl(1, i)
                a4(3, i) = a4(1, i)
                a4_tl(2, i) = a4_tl(1, i)
                a4(2, i) = a4(1, i)
                a4_tl(4, i) = 0.0_8
                a4(4, i) = 0.
              ELSE IF (a4(3, i) .GT. a4(2, i)) THEN
                a4_tl(4, i) = 3.*(a4_tl(2, i)-a4_tl(1, i))
                a4(4, i) = 3.*(a4(2, i)-a4(1, i))
                a4_tl(3, i) = a4_tl(2, i) - a4_tl(4, i)
                a4(3, i) = a4(2, i) - a4(4, i)
              ELSE
                a4_tl(4, i) = 3.*(a4_tl(3, i)-a4_tl(1, i))
                a4(4, i) = 3.*(a4(3, i)-a4(1, i))
                a4_tl(2, i) = a4_tl(3, i) - a4_tl(4, i)
                a4(2, i) = a4(3, i) - a4(4, i)
              END IF
            END IF
          END IF
        END IF
      END DO
    ELSE IF (iv .EQ. 1) THEN
      DO i=1,im
        IF ((a4(1, i)-a4(2, i))*(a4(1, i)-a4(3, i)) .GE. 0.) THEN
          a4_tl(2, i) = a4_tl(1, i)
          a4(2, i) = a4(1, i)
          a4_tl(3, i) = a4_tl(1, i)
          a4(3, i) = a4(1, i)
          a4_tl(4, i) = 0.0_8
          a4(4, i) = 0.
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(2, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            a4_tl(3, i) = a4_tl(2, i) - a4_tl(4, i)
            a4(3, i) = a4(2, i) - a4(4, i)
          ELSE IF (a6da .GT. da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(3, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            a4_tl(2, i) = a4_tl(3, i) - a4_tl(4, i)
            a4(2, i) = a4(3, i) - a4(4, i)
          END IF
        END IF
      END DO
    ELSE
! Standard PPM constraint
      DO i=1,im
        IF (extm(i)) THEN
          a4_tl(2, i) = a4_tl(1, i)
          a4(2, i) = a4(1, i)
          a4_tl(3, i) = a4_tl(1, i)
          a4(3, i) = a4(1, i)
          a4_tl(4, i) = 0.0_8
          a4(4, i) = 0.
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(2, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            a4_tl(3, i) = a4_tl(2, i) - a4_tl(4, i)
            a4(3, i) = a4(2, i) - a4(4, i)
          ELSE IF (a6da .GT. da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(3, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            a4_tl(2, i) = a4_tl(3, i) - a4_tl(4, i)
            a4(2, i) = a4(3, i) - a4(4, i)
          END IF
        END IF
      END DO
    END IF
  END SUBROUTINE CS_LIMITERS_TLM

  SUBROUTINE PPM_PROFILE_TLM(a4, a4_tl, delp, delp_tl, km, i1, i2, iv, &
&   kord)
    IMPLICIT NONE
! !INPUT PARAMETERS:
! iv =-1: winds
    INTEGER, INTENT(IN) :: iv
! iv = 0: positive definite scalars
! iv = 1: others
! iv = 2: temp (if remap_t) and w (iv=-2)
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! vertical dimension
    INTEGER, INTENT(IN) :: km
! Order (or more accurately method no.):
    INTEGER, INTENT(IN) :: kord
! 
! layer pressure thickness
    REAL(p_precision), INTENT(IN) :: delp(i1:i2, km)
    REAL(p_precision), INTENT(IN) :: delp_tl(i1:i2, km)
! !INPUT/OUTPUT PARAMETERS:
! Interpolated values
    REAL(p_precision), INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL(p_precision), INTENT(INOUT) :: a4_tl(4, i1:i2, km)
! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
! 
! !REVISION HISTORY: 
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
    REAL(p_precision) :: dc(i1:i2, km)
    REAL(p_precision) :: dc_tl(i1:i2, km)
    REAL(p_precision) :: h2(i1:i2, km)
    REAL(p_precision) :: h2_tl(i1:i2, km)
    REAL(p_precision) :: delq(i1:i2, km)
    REAL(p_precision) :: delq_tl(i1:i2, km)
    REAL(p_precision) :: df2(i1:i2, km)
    REAL(p_precision) :: df2_tl(i1:i2, km)
    REAL(p_precision) :: d4(i1:i2, km)
    REAL(p_precision) :: d4_tl(i1:i2, km)
! local scalars:
    INTEGER :: i, k, km1, lmt, it
    REAL(p_precision) :: fac
    REAL(p_precision) :: a1, a2, c1, c2, c3, d1, d2
    REAL(p_precision) :: a1_tl, a2_tl, c1_tl, c2_tl, c3_tl, d1_tl, d2_tl
    REAL(p_precision) :: qm, dq, lac, qmp, pmp
    REAL(p_precision) :: qm_tl, dq_tl, lac_tl, qmp_tl, pmp_tl
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL*8 :: y1_tl
    REAL*8 :: y2_tl
    REAL*8 :: min2
    REAL*8 :: min1
    REAL*8 :: y3_tl
    REAL*8 :: x1_tl
    REAL*8 :: y4_tl
    REAL*8 :: x3
    REAL*8 :: x2
    REAL*8 :: x1
    REAL*8 :: x2_tl
    REAL*8 :: y5_tl
    REAL*8 :: x3_tl
    REAL*8 :: y6_tl
    REAL*8 :: max1_tl
    REAL*8 :: y7_tl
    REAL*8 :: min1_tl
    REAL*8 :: y8_tl
    REAL*8 :: min2_tl
    REAL*8 :: z1
    REAL*8 :: y9_tl
    REAL*8 :: z1_tl
    INTEGER :: abs0
    REAL*8 :: y9
    REAL*8 :: y8
    REAL*8 :: y7
    REAL*8 :: y6
    REAL*8 :: max1
    REAL*8 :: y5
    REAL*8 :: y4
    REAL*8 :: y3
    REAL*8 :: y2
    REAL*8 :: y1
    km1 = km - 1
    it = i2 - i1 + 1
    delq_tl = 0.0_8
    d4_tl = 0.0_8
    DO k=2,km
      DO i=i1,i2
        delq_tl(i, k-1) = a4_tl(1, i, k) - a4_tl(1, i, k-1)
        delq(i, k-1) = a4(1, i, k) - a4(1, i, k-1)
        d4_tl(i, k) = delp_tl(i, k-1) + delp_tl(i, k)
        d4(i, k) = delp(i, k-1) + delp(i, k)
      END DO
    END DO
    df2_tl = 0.0_8
    dc_tl = 0.0_8
    DO k=2,km1
      DO i=i1,i2
        c1_tl = ((delp_tl(i, k-1)+0.5*delp_tl(i, k))*d4(i, k+1)-(delp(i&
&         , k-1)+0.5*delp(i, k))*d4_tl(i, k+1))/d4(i, k+1)**2
        c1 = (delp(i, k-1)+0.5*delp(i, k))/d4(i, k+1)
        c2_tl = ((delp_tl(i, k+1)+0.5*delp_tl(i, k))*d4(i, k)-(delp(i, k&
&         +1)+0.5*delp(i, k))*d4_tl(i, k))/d4(i, k)**2
        c2 = (delp(i, k+1)+0.5*delp(i, k))/d4(i, k)
        df2_tl(i, k) = ((delp_tl(i, k)*(c1*delq(i, k)+c2*delq(i, k-1))+&
&         delp(i, k)*(c1_tl*delq(i, k)+c1*delq_tl(i, k)+c2_tl*delq(i, k-&
&         1)+c2*delq_tl(i, k-1)))*(d4(i, k)+delp(i, k+1))-delp(i, k)*(c1&
&         *delq(i, k)+c2*delq(i, k-1))*(d4_tl(i, k)+delp_tl(i, k+1)))/(&
&         d4(i, k)+delp(i, k+1))**2
        df2(i, k) = delp(i, k)*(c1*delq(i, k)+c2*delq(i, k-1))/(d4(i, k)&
&         +delp(i, k+1))
        IF (df2(i, k) .GE. 0.) THEN
          x1_tl = df2_tl(i, k)
          x1 = df2(i, k)
        ELSE
          x1_tl = -df2_tl(i, k)
          x1 = -df2(i, k)
        END IF
        IF (a4(1, i, k-1) .LT. a4(1, i, k)) THEN
          IF (a4(1, i, k) .LT. a4(1, i, k+1)) THEN
            max1_tl = a4_tl(1, i, k+1)
            max1 = a4(1, i, k+1)
          ELSE
            max1_tl = a4_tl(1, i, k)
            max1 = a4(1, i, k)
          END IF
        ELSE IF (a4(1, i, k-1) .LT. a4(1, i, k+1)) THEN
          max1_tl = a4_tl(1, i, k+1)
          max1 = a4(1, i, k+1)
        ELSE
          max1_tl = a4_tl(1, i, k-1)
          max1 = a4(1, i, k-1)
        END IF
        y1_tl = max1_tl - a4_tl(1, i, k)
        y1 = max1 - a4(1, i, k)
        IF (a4(1, i, k-1) .GT. a4(1, i, k)) THEN
          IF (a4(1, i, k) .GT. a4(1, i, k+1)) THEN
            min2_tl = a4_tl(1, i, k+1)
            min2 = a4(1, i, k+1)
          ELSE
            min2_tl = a4_tl(1, i, k)
            min2 = a4(1, i, k)
          END IF
        ELSE IF (a4(1, i, k-1) .GT. a4(1, i, k+1)) THEN
          min2_tl = a4_tl(1, i, k+1)
          min2 = a4(1, i, k+1)
        ELSE
          min2_tl = a4_tl(1, i, k-1)
          min2 = a4(1, i, k-1)
        END IF
        z1_tl = a4_tl(1, i, k) - min2_tl
        z1 = a4(1, i, k) - min2
        IF (x1 .GT. y1) THEN
          IF (y1 .GT. z1) THEN
            min1_tl = z1_tl
            min1 = z1
          ELSE
            min1_tl = y1_tl
            min1 = y1
          END IF
        ELSE IF (x1 .GT. z1) THEN
          min1_tl = z1_tl
          min1 = z1
        ELSE
          min1_tl = x1_tl
          min1 = x1
        END IF
        dc_tl(i, k) = min1_tl*SIGN(1.d0, min1*df2(i, k))
        dc(i, k) = SIGN(min1, df2(i, k))
      END DO
    END DO
!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------
    DO k=3,km1
      DO i=i1,i2
        c1_tl = ((delq_tl(i, k-1)*delp(i, k-1)+delq(i, k-1)*delp_tl(i, k&
&         -1))*d4(i, k)-delq(i, k-1)*delp(i, k-1)*d4_tl(i, k))/d4(i, k)&
&         **2
        c1 = delq(i, k-1)*delp(i, k-1)/d4(i, k)
        a1_tl = (d4_tl(i, k-1)*(d4(i, k)+delp(i, k-1))-d4(i, k-1)*(d4_tl&
&         (i, k)+delp_tl(i, k-1)))/(d4(i, k)+delp(i, k-1))**2
        a1 = d4(i, k-1)/(d4(i, k)+delp(i, k-1))
        a2_tl = (d4_tl(i, k+1)*(d4(i, k)+delp(i, k))-d4(i, k+1)*(d4_tl(i&
&         , k)+delp_tl(i, k)))/(d4(i, k)+delp(i, k))**2
        a2 = d4(i, k+1)/(d4(i, k)+delp(i, k))
        a4_tl(2, i, k) = a4_tl(1, i, k-1) + c1_tl + 2.*(delp_tl(i, k)*(&
&         c1*(a1-a2)+a2*dc(i, k-1))+delp(i, k)*(c1_tl*(a1-a2)+c1*(a1_tl-&
&         a2_tl)+a2_tl*dc(i, k-1)+a2*dc_tl(i, k-1))-delp_tl(i, k-1)*a1*&
&         dc(i, k)-delp(i, k-1)*(a1_tl*dc(i, k)+a1*dc_tl(i, k)))/(d4(i, &
&         k-1)+d4(i, k+1)) - 2.*(d4_tl(i, k-1)+d4_tl(i, k+1))*(delp(i, k&
&         )*(c1*(a1-a2)+a2*dc(i, k-1))-delp(i, k-1)*a1*dc(i, k))/(d4(i, &
&         k-1)+d4(i, k+1))**2
        a4(2, i, k) = a4(1, i, k-1) + c1 + 2./(d4(i, k-1)+d4(i, k+1))*(&
&         delp(i, k)*(c1*(a1-a2)+a2*dc(i, k-1))-delp(i, k-1)*a1*dc(i, k)&
&         )
      END DO
    END DO
    IF (km .GT. 8 .AND. kord .GT. 4) CALL STEEPZ_TLM(i1, i2, km, a4, &
&                                              a4_tl, df2, df2_tl, dc, &
&                                              dc_tl, delq, delq_tl, &
&                                              delp, delp_tl, d4, d4_tl)
! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
    DO i=i1,i2
      d1_tl = delp_tl(i, 1)
      d1 = delp(i, 1)
      d2_tl = delp_tl(i, 2)
      d2 = delp(i, 2)
      qm_tl = ((d2_tl*a4(1, i, 1)+d2*a4_tl(1, i, 1)+d1_tl*a4(1, i, 2)+d1&
&       *a4_tl(1, i, 2))*(d1+d2)-(d2*a4(1, i, 1)+d1*a4(1, i, 2))*(d1_tl+&
&       d2_tl))/(d1+d2)**2
      qm = (d2*a4(1, i, 1)+d1*a4(1, i, 2))/(d1+d2)
      dq_tl = (2.*(a4_tl(1, i, 2)-a4_tl(1, i, 1))*(d1+d2)-2.*(a4(1, i, 2&
&       )-a4(1, i, 1))*(d1_tl+d2_tl))/(d1+d2)**2
      dq = 2.*(a4(1, i, 2)-a4(1, i, 1))/(d1+d2)
      c1_tl = (4.*(a4_tl(2, i, 3)-qm_tl-d2_tl*dq-d2*dq_tl)*d2*(2.*d2*d2+&
&       d1*(d2+3.*d1))-4.*(a4(2, i, 3)-qm-d2*dq)*(d2_tl*(2.*d2*d2+d1*(d2&
&       +3.*d1))+d2*(2.*(d2_tl*d2+d2*d2_tl)+d1_tl*(d2+3.*d1)+d1*(d2_tl+&
&       3.*d1_tl))))/(d2*(2.*d2*d2+d1*(d2+3.*d1)))**2
      c1 = 4.*(a4(2, i, 3)-qm-d2*dq)/(d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3_tl = dq_tl - 0.5*(c1_tl*(d2*(5.*d1+d2)-3.*d1*d1)+c1*(d2_tl*(5.*&
&       d1+d2)+d2*(5.*d1_tl+d2_tl)-3.*(d1_tl*d1+d1*d1_tl)))
      c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      a4_tl(2, i, 2) = qm_tl - 0.25*(((c1_tl*d1+c1*d1_tl)*d2+c1*d1*d2_tl&
&       )*(d2+3.*d1)+c1*d1*d2*(d2_tl+3.*d1_tl))
      a4(2, i, 2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
      a4_tl(2, i, 1) = d1_tl*(2.*c1*d1**2-c3) + d1*(2.*(c1_tl*d1**2+c1*2&
&       *d1*d1_tl)-c3_tl) + a4_tl(2, i, 2)
      a4(2, i, 1) = d1*(2.*c1*d1**2-c3) + a4(2, i, 2)
      IF (a4(1, i, 1) .GT. a4(1, i, 2)) THEN
        y2_tl = a4_tl(1, i, 2)
        y2 = a4(1, i, 2)
      ELSE
        y2_tl = a4_tl(1, i, 1)
        y2 = a4(1, i, 1)
      END IF
      IF (a4(2, i, 2) .LT. y2) THEN
        a4_tl(2, i, 2) = y2_tl
        a4(2, i, 2) = y2
      ELSE
        a4(2, i, 2) = a4(2, i, 2)
      END IF
      IF (a4(1, i, 1) .LT. a4(1, i, 2)) THEN
        y3_tl = a4_tl(1, i, 2)
        y3 = a4(1, i, 2)
      ELSE
        y3_tl = a4_tl(1, i, 1)
        y3 = a4(1, i, 1)
      END IF
      IF (a4(2, i, 2) .GT. y3) THEN
        a4_tl(2, i, 2) = y3_tl
        a4(2, i, 2) = y3
      ELSE
        a4(2, i, 2) = a4(2, i, 2)
      END IF
      dc_tl(i, 1) = 0.5*(a4_tl(2, i, 2)-a4_tl(1, i, 1))
      dc(i, 1) = 0.5*(a4(2, i, 2)-a4(1, i, 1))
    END DO
! Enforce monotonicity  within the top layer
    IF (iv .EQ. 0) THEN
      DO i=i1,i2
        IF (0. .LT. a4(2, i, 1)) THEN
          a4(2, i, 1) = a4(2, i, 1)
        ELSE
          a4_tl(2, i, 1) = 0.0_8
          a4(2, i, 1) = 0.
        END IF
        IF (0. .LT. a4(2, i, 2)) THEN
          a4(2, i, 2) = a4(2, i, 2)
        ELSE
          a4_tl(2, i, 2) = 0.0_8
          a4(2, i, 2) = 0.
        END IF
      END DO
    ELSE IF (iv .EQ. -1) THEN
      DO i=i1,i2
        IF (a4(2, i, 1)*a4(1, i, 1) .LE. 0.) THEN
          a4_tl(2, i, 1) = 0.0_8
          a4(2, i, 1) = 0.
        END IF
      END DO
    ELSE
      IF (iv .GE. 0.) THEN
        abs0 = iv
      ELSE
        abs0 = -iv
      END IF
      IF (abs0 .EQ. 2) THEN
        DO i=i1,i2
          a4_tl(2, i, 1) = a4_tl(1, i, 1)
          a4(2, i, 1) = a4(1, i, 1)
          a4_tl(3, i, 1) = a4_tl(1, i, 1)
          a4(3, i, 1) = a4(1, i, 1)
        END DO
      END IF
    END IF
! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
    DO i=i1,i2
      d1_tl = delp_tl(i, km)
      d1 = delp(i, km)
      d2_tl = delp_tl(i, km1)
      d2 = delp(i, km1)
      qm_tl = ((d2_tl*a4(1, i, km)+d2*a4_tl(1, i, km)+d1_tl*a4(1, i, km1&
&       )+d1*a4_tl(1, i, km1))*(d1+d2)-(d2*a4(1, i, km)+d1*a4(1, i, km1)&
&       )*(d1_tl+d2_tl))/(d1+d2)**2
      qm = (d2*a4(1, i, km)+d1*a4(1, i, km1))/(d1+d2)
      dq_tl = (2.*(a4_tl(1, i, km1)-a4_tl(1, i, km))*(d1+d2)-2.*(a4(1, i&
&       , km1)-a4(1, i, km))*(d1_tl+d2_tl))/(d1+d2)**2
      dq = 2.*(a4(1, i, km1)-a4(1, i, km))/(d1+d2)
      c1_tl = ((a4_tl(2, i, km1)-qm_tl-d2_tl*dq-d2*dq_tl)*d2*(2.*d2*d2+&
&       d1*(d2+3.*d1))-(a4(2, i, km1)-qm-d2*dq)*(d2_tl*(2.*d2*d2+d1*(d2+&
&       3.*d1))+d2*(2.*(d2_tl*d2+d2*d2_tl)+d1_tl*(d2+3.*d1)+d1*(d2_tl+3.&
&       *d1_tl))))/(d2*(2.*d2*d2+d1*(d2+3.*d1)))**2
      c1 = (a4(2, i, km1)-qm-d2*dq)/(d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3_tl = dq_tl - 2.0*(c1_tl*(d2*(5.*d1+d2)-3.*d1*d1)+c1*(d2_tl*(5.*&
&       d1+d2)+d2*(5.*d1_tl+d2_tl)-3.*(d1_tl*d1+d1*d1_tl)))
      c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
      a4_tl(2, i, km) = qm_tl - ((c1_tl*d1+c1*d1_tl)*d2+c1*d1*d2_tl)*(d2&
&       +3.*d1) - c1*d1*d2*(d2_tl+3.*d1_tl)
      a4(2, i, km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
      a4_tl(3, i, km) = d1_tl*(8.*c1*d1**2-c3) + d1*(8.*(c1_tl*d1**2+c1*&
&       2*d1*d1_tl)-c3_tl) + a4_tl(2, i, km)
      a4(3, i, km) = d1*(8.*c1*d1**2-c3) + a4(2, i, km)
      IF (a4(1, i, km) .GT. a4(1, i, km1)) THEN
        y4_tl = a4_tl(1, i, km1)
        y4 = a4(1, i, km1)
      ELSE
        y4_tl = a4_tl(1, i, km)
        y4 = a4(1, i, km)
      END IF
      IF (a4(2, i, km) .LT. y4) THEN
        a4_tl(2, i, km) = y4_tl
        a4(2, i, km) = y4
      ELSE
        a4(2, i, km) = a4(2, i, km)
      END IF
      IF (a4(1, i, km) .LT. a4(1, i, km1)) THEN
        y5_tl = a4_tl(1, i, km1)
        y5 = a4(1, i, km1)
      ELSE
        y5_tl = a4_tl(1, i, km)
        y5 = a4(1, i, km)
      END IF
      IF (a4(2, i, km) .GT. y5) THEN
        a4_tl(2, i, km) = y5_tl
        a4(2, i, km) = y5
      ELSE
        a4(2, i, km) = a4(2, i, km)
      END IF
      dc_tl(i, km) = 0.5*(a4_tl(1, i, km)-a4_tl(2, i, km))
      dc(i, km) = 0.5*(a4(1, i, km)-a4(2, i, km))
    END DO
! Enforce constraint on the "slope" at the surface
    IF (iv .EQ. 0) THEN
      DO i=i1,i2
        IF (0. .LT. a4(2, i, km)) THEN
          a4(2, i, km) = a4(2, i, km)
        ELSE
          a4_tl(2, i, km) = 0.0_8
          a4(2, i, km) = 0.
        END IF
        IF (0. .LT. a4(3, i, km)) THEN
          a4(3, i, km) = a4(3, i, km)
        ELSE
          a4_tl(3, i, km) = 0.0_8
          a4(3, i, km) = 0.
        END IF
      END DO
    ELSE IF (iv .LT. 0) THEN
      DO i=i1,i2
        IF (a4(1, i, km)*a4(3, i, km) .LE. 0.) THEN
          a4_tl(3, i, km) = 0.0_8
          a4(3, i, km) = 0.
        END IF
      END DO
    END IF
    DO k=1,km1
      DO i=i1,i2
        a4_tl(3, i, k) = a4_tl(2, i, k+1)
        a4(3, i, k) = a4(2, i, k+1)
      END DO
    END DO
!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
    DO k=1,2
      DO i=i1,i2
        a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(3, i&
&         , k))
        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
      END DO
      CALL PPM_LIMITERS_TLM(dc(i1, k), dc_tl(i1, k), a4(1, i1, k), a4_tl&
&                     (1, i1, k), it, 0)
    END DO
    IF (kord .GE. 7) THEN
      h2_tl = 0.0_8
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      DO k=2,km1
        DO i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
          h2_tl(i, k) = (2.*((dc_tl(i, k+1)*delp(i, k+1)-dc(i, k+1)*&
&           delp_tl(i, k+1))/delp(i, k+1)**2-(dc_tl(i, k-1)*delp(i, k-1)&
&           -dc(i, k-1)*delp_tl(i, k-1))/delp(i, k-1)**2)*(delp(i, k)+&
&           0.5*(delp(i, k-1)+delp(i, k+1)))-2.*(dc(i, k+1)/delp(i, k+1)&
&           -dc(i, k-1)/delp(i, k-1))*(delp_tl(i, k)+0.5*(delp_tl(i, k-1&
&           )+delp_tl(i, k+1))))*delp(i, k)**2/(delp(i, k)+0.5*(delp(i, &
&           k-1)+delp(i, k+1)))**2 + 2.*(dc(i, k+1)/delp(i, k+1)-dc(i, k&
&           -1)/delp(i, k-1))*2*delp(i, k)*delp_tl(i, k)/(delp(i, k)+0.5&
&           *(delp(i, k-1)+delp(i, k+1)))
          h2(i, k) = 2.*(dc(i, k+1)/delp(i, k+1)-dc(i, k-1)/delp(i, k-1)&
&           )/(delp(i, k)+0.5*(delp(i, k-1)+delp(i, k+1)))*delp(i, k)**2
        END DO
      END DO
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
! original quasi-monotone
      fac = 1.5
      DO k=3,km-2
        DO i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
          pmp_tl = 2.*dc_tl(i, k)
          pmp = 2.*dc(i, k)
          qmp_tl = a4_tl(1, i, k) + pmp_tl
          qmp = a4(1, i, k) + pmp
          lac_tl = a4_tl(1, i, k) + fac*h2_tl(i, k-1) + dc_tl(i, k)
          lac = a4(1, i, k) + fac*h2(i, k-1) + dc(i, k)
          IF (a4(1, i, k) .GT. qmp) THEN
            IF (qmp .GT. lac) THEN
              y8_tl = lac_tl
              y8 = lac
            ELSE
              y8_tl = qmp_tl
              y8 = qmp
            END IF
          ELSE IF (a4(1, i, k) .GT. lac) THEN
            y8_tl = lac_tl
            y8 = lac
          ELSE
            y8_tl = a4_tl(1, i, k)
            y8 = a4(1, i, k)
          END IF
          IF (a4(3, i, k) .LT. y8) THEN
            x2_tl = y8_tl
            x2 = y8
          ELSE
            x2_tl = a4_tl(3, i, k)
            x2 = a4(3, i, k)
          END IF
          IF (a4(1, i, k) .LT. qmp) THEN
            IF (qmp .LT. lac) THEN
              y6_tl = lac_tl
              y6 = lac
            ELSE
              y6_tl = qmp_tl
              y6 = qmp
            END IF
          ELSE IF (a4(1, i, k) .LT. lac) THEN
            y6_tl = lac_tl
            y6 = lac
          ELSE
            y6_tl = a4_tl(1, i, k)
            y6 = a4(1, i, k)
          END IF
          IF (x2 .GT. y6) THEN
            a4_tl(3, i, k) = y6_tl
            a4(3, i, k) = y6
          ELSE
            a4_tl(3, i, k) = x2_tl
            a4(3, i, k) = x2
          END IF
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
          qmp_tl = a4_tl(1, i, k) - pmp_tl
          qmp = a4(1, i, k) - pmp
          lac_tl = a4_tl(1, i, k) + fac*h2_tl(i, k+1) - dc_tl(i, k)
          lac = a4(1, i, k) + fac*h2(i, k+1) - dc(i, k)
          IF (a4(1, i, k) .GT. qmp) THEN
            IF (qmp .GT. lac) THEN
              y9_tl = lac_tl
              y9 = lac
            ELSE
              y9_tl = qmp_tl
              y9 = qmp
            END IF
          ELSE IF (a4(1, i, k) .GT. lac) THEN
            y9_tl = lac_tl
            y9 = lac
          ELSE
            y9_tl = a4_tl(1, i, k)
            y9 = a4(1, i, k)
          END IF
          IF (a4(2, i, k) .LT. y9) THEN
            x3_tl = y9_tl
            x3 = y9
          ELSE
            x3_tl = a4_tl(2, i, k)
            x3 = a4(2, i, k)
          END IF
          IF (a4(1, i, k) .LT. qmp) THEN
            IF (qmp .LT. lac) THEN
              y7_tl = lac_tl
              y7 = lac
            ELSE
              y7_tl = qmp_tl
              y7 = qmp
            END IF
          ELSE IF (a4(1, i, k) .LT. lac) THEN
            y7_tl = lac_tl
            y7 = lac
          ELSE
            y7_tl = a4_tl(1, i, k)
            y7 = a4(1, i, k)
          END IF
          IF (x3 .GT. y7) THEN
            a4_tl(2, i, k) = y7_tl
            a4(2, i, k) = y7
          ELSE
            a4_tl(2, i, k) = x3_tl
            a4(2, i, k) = x3
          END IF
!-------------
! Recompute A6
!-------------
          a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(3&
&           , i, k))
          a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
        END DO
! Additional constraint to ensure positivity when kord=7
        IF (iv .EQ. 0 .AND. kord .GE. 6) CALL PPM_LIMITERS_TLM(dc(i1, k)&
&                                                        , dc_tl(i1, k)&
&                                                        , a4(1, i1, k)&
&                                                        , a4_tl(1, i1, &
&                                                        k), it, 2)
      END DO
    ELSE
      lmt = kord - 3
      IF (0 .LT. lmt) THEN
        lmt = lmt
      ELSE
        lmt = 0
      END IF
      IF (iv .EQ. 0) THEN
        IF (2 .GT. lmt) THEN
          lmt = lmt
        ELSE
          lmt = 2
        END IF
      END IF
      DO k=3,km-2
        IF (kord .NE. 4) THEN
          DO i=i1,i2
            a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(&
&             3, i, k))
            a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
          END DO
        END IF
        IF (kord .NE. 6) CALL PPM_LIMITERS_TLM(dc(i1, k), dc_tl(i1, k), &
&                                        a4(1, i1, k), a4_tl(1, i1, k), &
&                                        it, lmt)
      END DO
    END IF
    DO k=km1,km
      DO i=i1,i2
        a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(3, i&
&         , k))
        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
      END DO
      CALL PPM_LIMITERS_TLM(dc(i1, k), dc_tl(i1, k), a4(1, i1, k), a4_tl&
&                     (1, i1, k), it, 0)
    END DO
  END SUBROUTINE PPM_PROFILE_TLM

  SUBROUTINE PPM_LIMITERS_TLM(dm, dm_tl, a4, a4_tl, itot, lmt)
    IMPLICIT NONE
! !INPUT PARAMETERS:
! the linear slope
    REAL(p_precision), INTENT(IN) :: dm(*)
    REAL(p_precision), INTENT(IN) :: dm_tl(*)
! Total Longitudes
    INTEGER, INTENT(IN) :: itot
! 0: Standard PPM constraint
    INTEGER, INTENT(IN) :: lmt
! 1: Improved full monotonicity constraint (Lin)
! 2: Positive definite constraint
! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
! PPM array
    REAL(p_precision), INTENT(INOUT) :: a4(4, *)
    REAL(p_precision), INTENT(INOUT) :: a4_tl(4, *)
! AA <-- a4(1,i)
! AL <-- a4(2,i)
! AR <-- a4(3,i)
! A6 <-- a4(4,i)
! !LOCAL VARIABLES:
    REAL(p_precision) :: qmp
    REAL(p_precision) :: qmp_tl
    REAL(p_precision) :: da1, da2, a6da
    REAL(p_precision) :: fmin
    INTEGER :: i
    INTRINSIC ABS
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL*8 :: y1_tl
    REAL*8 :: y2_tl
    REAL*8 :: min2
    REAL*8 :: min1
    REAL*8 :: x1_tl
    REAL*8 :: x2
    REAL*8 :: x1
    REAL*8 :: x2_tl
    REAL*8 :: min1_tl
    REAL*8 :: min2_tl
    REAL*8 :: abs0
    REAL*8 :: y2
    REAL*8 :: y1
! Developer: S.-J. Lin
    IF (lmt .EQ. 3) THEN
      RETURN
    ELSE IF (lmt .EQ. 0) THEN
! Standard PPM constraint
      DO i=1,itot
        IF (dm(i) .EQ. 0.) THEN
          a4_tl(2, i) = a4_tl(1, i)
          a4(2, i) = a4(1, i)
          a4_tl(3, i) = a4_tl(1, i)
          a4(3, i) = a4(1, i)
          a4_tl(4, i) = 0.0_8
          a4(4, i) = 0.
        ELSE
          da1 = a4(3, i) - a4(2, i)
          da2 = da1**2
          a6da = a4(4, i)*da1
          IF (a6da .LT. -da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(2, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(2, i)-a4(1, i))
            a4_tl(3, i) = a4_tl(2, i) - a4_tl(4, i)
            a4(3, i) = a4(2, i) - a4(4, i)
          ELSE IF (a6da .GT. da2) THEN
            a4_tl(4, i) = 3.*(a4_tl(3, i)-a4_tl(1, i))
            a4(4, i) = 3.*(a4(3, i)-a4(1, i))
            a4_tl(2, i) = a4_tl(3, i) - a4_tl(4, i)
            a4(2, i) = a4(3, i) - a4(4, i)
          END IF
        END IF
      END DO
    ELSE IF (lmt .EQ. 1) THEN
! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      DO i=1,itot
        qmp_tl = 2.*dm_tl(i)
        qmp = 2.*dm(i)
        IF (qmp .GE. 0.) THEN
          x1_tl = qmp_tl
          x1 = qmp
        ELSE
          x1_tl = -qmp_tl
          x1 = -qmp
        END IF
        IF (a4(2, i) - a4(1, i) .GE. 0.) THEN
          y1_tl = a4_tl(2, i) - a4_tl(1, i)
          y1 = a4(2, i) - a4(1, i)
        ELSE
          y1_tl = -(a4_tl(2, i)-a4_tl(1, i))
          y1 = -(a4(2, i)-a4(1, i))
        END IF
        IF (x1 .GT. y1) THEN
          min1_tl = y1_tl
          min1 = y1
        ELSE
          min1_tl = x1_tl
          min1 = x1
        END IF
        a4_tl(2, i) = a4_tl(1, i) - min1_tl*SIGN(1.d0, min1*qmp)
        a4(2, i) = a4(1, i) - SIGN(min1, qmp)
        IF (qmp .GE. 0.) THEN
          x2_tl = qmp_tl
          x2 = qmp
        ELSE
          x2_tl = -qmp_tl
          x2 = -qmp
        END IF
        IF (a4(3, i) - a4(1, i) .GE. 0.) THEN
          y2_tl = a4_tl(3, i) - a4_tl(1, i)
          y2 = a4(3, i) - a4(1, i)
        ELSE
          y2_tl = -(a4_tl(3, i)-a4_tl(1, i))
          y2 = -(a4(3, i)-a4(1, i))
        END IF
        IF (x2 .GT. y2) THEN
          min2_tl = y2_tl
          min2 = y2
        ELSE
          min2_tl = x2_tl
          min2 = x2
        END IF
        a4_tl(3, i) = a4_tl(1, i) + min2_tl*SIGN(1.d0, min2*qmp)
        a4(3, i) = a4(1, i) + SIGN(min2, qmp)
        a4_tl(4, i) = 3.*(2.*a4_tl(1, i)-a4_tl(2, i)-a4_tl(3, i))
        a4(4, i) = 3.*(2.*a4(1, i)-(a4(2, i)+a4(3, i)))
      END DO
    ELSE IF (lmt .EQ. 2) THEN
! Positive definite constraint
      DO i=1,itot
        IF (a4(3, i) - a4(2, i) .GE. 0.) THEN
          abs0 = a4(3, i) - a4(2, i)
        ELSE
          abs0 = -(a4(3, i)-a4(2, i))
        END IF
        IF (abs0 .LT. -a4(4, i)) THEN
          fmin = a4(1, i) + 0.25*(a4(3, i)-a4(2, i))**2/a4(4, i) + a4(4&
&           , i)*r12
          IF (fmin .LT. 0.) THEN
            IF (a4(1, i) .LT. a4(3, i) .AND. a4(1, i) .LT. a4(2, i)) &
&           THEN
              a4_tl(3, i) = a4_tl(1, i)
              a4(3, i) = a4(1, i)
              a4_tl(2, i) = a4_tl(1, i)
              a4(2, i) = a4(1, i)
              a4_tl(4, i) = 0.0_8
              a4(4, i) = 0.
            ELSE IF (a4(3, i) .GT. a4(2, i)) THEN
              a4_tl(4, i) = 3.*(a4_tl(2, i)-a4_tl(1, i))
              a4(4, i) = 3.*(a4(2, i)-a4(1, i))
              a4_tl(3, i) = a4_tl(2, i) - a4_tl(4, i)
              a4(3, i) = a4(2, i) - a4(4, i)
            ELSE
              a4_tl(4, i) = 3.*(a4_tl(3, i)-a4_tl(1, i))
              a4(4, i) = 3.*(a4(3, i)-a4(1, i))
              a4_tl(2, i) = a4_tl(3, i) - a4_tl(4, i)
              a4(2, i) = a4(3, i) - a4(4, i)
            END IF
          END IF
        END IF
      END DO
    END IF
  END SUBROUTINE PPM_LIMITERS_TLM

  SUBROUTINE STEEPZ_TLM(i1, i2, km, a4, a4_tl, df2, df2_tl, dm, dm_tl, &
&   dq, dq_tl, dp, dp_tl, d4, d4_tl)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km, i1, i2
! grid size
    REAL(p_precision), INTENT(IN) :: dp(i1:i2, km)
    REAL(p_precision), INTENT(IN) :: dp_tl(i1:i2, km)
! backward diff of q
    REAL(p_precision), INTENT(IN) :: dq(i1:i2, km)
    REAL(p_precision), INTENT(IN) :: dq_tl(i1:i2, km)
! backward sum:  dp(k)+ dp(k-1) 
    REAL(p_precision), INTENT(IN) :: d4(i1:i2, km)
    REAL(p_precision), INTENT(IN) :: d4_tl(i1:i2, km)
! first guess mismatch
    REAL(p_precision), INTENT(IN) :: df2(i1:i2, km)
    REAL(p_precision), INTENT(IN) :: df2_tl(i1:i2, km)
! monotonic mismatch
    REAL(p_precision), INTENT(IN) :: dm(i1:i2, km)
    REAL(p_precision), INTENT(IN) :: dm_tl(i1:i2, km)
! !INPUT/OUTPUT PARAMETERS:
! first guess/steepened
    REAL(p_precision), INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL(p_precision), INTENT(INOUT) :: a4_tl(4, i1:i2, km)
! !LOCAL VARIABLES:
    INTEGER :: i, k
    REAL(p_precision) :: alfa(i1:i2, km)
    REAL(p_precision) :: alfa_tl(i1:i2, km)
    REAL(p_precision) :: f(i1:i2, km)
    REAL(p_precision) :: f_tl(i1:i2, km)
    REAL(p_precision) :: rat(i1:i2, km)
    REAL(p_precision) :: rat_tl(i1:i2, km)
    REAL(p_precision) :: dg2
    REAL(p_precision) :: dg2_tl
    INTRINSIC MIN
    INTRINSIC MAX
    REAL*8 :: y1_tl
    REAL*8 :: y1
    rat_tl = 0.0_8
! Compute ratio of dq/dp
    DO k=2,km
      DO i=i1,i2
        rat_tl(i, k) = (dq_tl(i, k-1)*d4(i, k)-dq(i, k-1)*d4_tl(i, k))/&
&         d4(i, k)**2
        rat(i, k) = dq(i, k-1)/d4(i, k)
      END DO
    END DO
    f_tl = 0.0_8
! Compute F
    DO k=2,km-1
      DO i=i1,i2
        f_tl(i, k) = ((rat_tl(i, k+1)-rat_tl(i, k))*(dp(i, k-1)+dp(i, k)&
&         +dp(i, k+1))-(rat(i, k+1)-rat(i, k))*(dp_tl(i, k-1)+dp_tl(i, k&
&         )+dp_tl(i, k+1)))/(dp(i, k-1)+dp(i, k)+dp(i, k+1))**2
        f(i, k) = (rat(i, k+1)-rat(i, k))/(dp(i, k-1)+dp(i, k)+dp(i, k+1&
&         ))
      END DO
    END DO
    alfa_tl = 0.0_8
    DO k=3,km-2
      DO i=i1,i2
        IF (f(i, k+1)*f(i, k-1) .LT. 0. .AND. df2(i, k) .NE. 0.) THEN
          dg2_tl = (f_tl(i, k+1)-f_tl(i, k-1))*((dp(i, k+1)-dp(i, k-1))&
&           **2+d4(i, k)*d4(i, k+1)) + (f(i, k+1)-f(i, k-1))*(2*(dp(i, k&
&           +1)-dp(i, k-1))*(dp_tl(i, k+1)-dp_tl(i, k-1))+d4_tl(i, k)*d4&
&           (i, k+1)+d4(i, k)*d4_tl(i, k+1))
          dg2 = (f(i, k+1)-f(i, k-1))*((dp(i, k+1)-dp(i, k-1))**2+d4(i, &
&           k)*d4(i, k+1))
          IF (0.5 .GT. -(0.1875*dg2/df2(i, k))) THEN
            y1_tl = -((0.1875*dg2_tl*df2(i, k)-0.1875*dg2*df2_tl(i, k))/&
&             df2(i, k)**2)
            y1 = -(0.1875*dg2/df2(i, k))
          ELSE
            y1 = 0.5
            y1_tl = 0.0_8
          END IF
          IF (0. .LT. y1) THEN
            alfa_tl(i, k) = y1_tl
            alfa(i, k) = y1
          ELSE
            alfa_tl(i, k) = 0.0_8
            alfa(i, k) = 0.
          END IF
        ELSE
          alfa_tl(i, k) = 0.0_8
          alfa(i, k) = 0.
        END IF
      END DO
    END DO
    DO k=4,km-2
      DO i=i1,i2
        a4_tl(2, i, k) = (-alfa_tl(i, k-1)-alfa_tl(i, k))*a4(2, i, k) + &
&         (1.-alfa(i, k-1)-alfa(i, k))*a4_tl(2, i, k) + alfa_tl(i, k-1)*&
&         (a4(1, i, k)-dm(i, k)) + alfa(i, k-1)*(a4_tl(1, i, k)-dm_tl(i&
&         , k)) + alfa_tl(i, k)*(a4(1, i, k-1)+dm(i, k-1)) + alfa(i, k)*&
&         (a4_tl(1, i, k-1)+dm_tl(i, k-1))
        a4(2, i, k) = (1.-alfa(i, k-1)-alfa(i, k))*a4(2, i, k) + alfa(i&
&         , k-1)*(a4(1, i, k)-dm(i, k)) + alfa(i, k)*(a4(1, i, k-1)+dm(i&
&         , k-1))
      END DO
    END DO
  END SUBROUTINE STEEPZ_TLM

  SUBROUTINE MAP1_Q2_NOLIMITERS_TLM(km, pe1, pe1_tl, q1, q1_tl, kn, pe2, &
&   pe2_tl, q2, q2_tl, dp2, dp2_tl, i1, i2, iv, kord, j, ibeg, iend, &
&   jbeg, jend)
    IMPLICIT NONE
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: i1, i2
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Mode: 0 ==  constituents 1 == ???
    INTEGER, INTENT(IN) :: iv
    INTEGER, INTENT(IN) :: kord
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! Target vertical dimension
    INTEGER, INTENT(IN) :: kn
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe1(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL(p_precision), INTENT(IN) :: pe2_tl(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! Field input
    REAL, INTENT(IN) :: q1(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(IN) :: q1_tl(ibeg:iend, jbeg:jend, km)
    REAL(p_precision), INTENT(IN) :: dp2(i1:i2, kn)
    REAL(p_precision), INTENT(IN) :: dp2_tl(i1:i2, kn)
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL(p_precision), INTENT(INOUT) :: q2(i1:i2, kn)
    REAL(p_precision), INTENT(INOUT) :: q2_tl(i1:i2, kn)
! !LOCAL VARIABLES:
    REAL(p_precision) :: dp1(i1:i2, km)
    REAL(p_precision) :: dp1_tl(i1:i2, km)
    REAL(p_precision) :: q4(4, i1:i2, km)
    REAL(p_precision) :: q4_tl(4, i1:i2, km)
    REAL(p_precision) :: pl, pr, qsum, dp, esl
    REAL(p_precision) :: pl_tl, pr_tl, qsum_tl, dp_tl, esl_tl
    INTEGER :: i, k, l, m, k0
    dp1_tl = 0.0_8
    q4_tl = 0.0_8
    DO k=1,km
      DO i=i1,i2
        dp1_tl(i, k) = pe1_tl(i, k+1) - pe1_tl(i, k)
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4_tl(1, i, k) = q1_tl(i, j, k)
        q4(1, i, k) = q1(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE_NOLIMITERS_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv&
&                          , kord)
      qsum_tl = 0.0_8
    ELSE
      CALL PPM_PROFILE_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv, kord)
      qsum_tl = 0.0_8
    END IF
! Mapping
    DO i=i1,i2
      k0 = 1
      DO 555 k=1,kn
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) GOTO 110
        END DO
        GOTO 123
 110    pl_tl = ((pe2_tl(i, k)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k)-pe1(i&
&         , l))*dp1_tl(i, l))/dp1(i, l)**2
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr_tl = ((pe2_tl(i, k+1)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k+1)-&
&           pe1(i, l))*dp1_tl(i, l))/dp1(i, l)**2
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          q2_tl(i, k) = q4_tl(2, i, l) + 0.5*((q4_tl(4, i, l)+q4_tl(3, i&
&           , l)-q4_tl(2, i, l))*(pr+pl)+(q4(4, i, l)+q4(3, i, l)-q4(2, &
&           i, l))*(pr_tl+pl_tl)) - r3*(q4_tl(4, i, l)*(pr*(pr+pl)+pl**2&
&           )+q4(4, i, l)*(pr_tl*(pr+pl)+pr*(pr_tl+pl_tl)+2*pl*pl_tl))
          k0 = l
          GOTO 555
        ELSE
! Fractional area...
          qsum_tl = (pe1_tl(i, l+1)-pe2_tl(i, k))*(q4(2, i, l)+0.5*(q4(4&
&           , i, l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.&
&           +pl*(1.+pl)))) + (pe1(i, l+1)-pe2(i, k))*(q4_tl(2, i, l)+0.5&
&           *((q4_tl(4, i, l)+q4_tl(3, i, l)-q4_tl(2, i, l))*(1.+pl)+(q4&
&           (4, i, l)+q4(3, i, l)-q4(2, i, l))*pl_tl)-r3*(q4_tl(4, i, l)&
&           *(1.+pl*(1.+pl))+q4(4, i, l)*(pl_tl*(1.+pl)+pl*pl_tl)))
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
              qsum_tl = qsum_tl + dp1_tl(i, m)*q4(1, i, m) + dp1(i, m)*&
&               q4_tl(1, i, m)
              qsum = qsum + dp1(i, m)*q4(1, i, m)
            ELSE
              GOTO 120
            END IF
          END DO
          GOTO 123
 120      dp_tl = pe2_tl(i, k+1) - pe1_tl(i, m)
          dp = pe2(i, k+1) - pe1(i, m)
          esl_tl = (dp_tl*dp1(i, m)-dp*dp1_tl(i, m))/dp1(i, m)**2
          esl = dp/dp1(i, m)
          qsum_tl = qsum_tl + dp_tl*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4&
&           (2, i, m)+q4(4, i, m)*(1.-r23*esl))) + dp*(q4_tl(2, i, m)+&
&           0.5*(esl_tl*(q4(3, i, m)-q4(2, i, m)+q4(4, i, m)*(1.-r23*esl&
&           ))+esl*(q4_tl(3, i, m)-q4_tl(2, i, m)+q4_tl(4, i, m)*(1.-r23&
&           *esl)-q4(4, i, m)*r23*esl_tl)))
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
        END IF
 123    q2_tl(i, k) = (qsum_tl*dp2(i, k)-qsum*dp2_tl(i, k))/dp2(i, k)**2
 555  CONTINUE
    END DO
  END SUBROUTINE MAP1_Q2_NOLIMITERS_TLM

  SUBROUTINE REMAP_Z_NOLIMITERS_TLM(km, pe1, pe1_tl, pe2, pe2_tl, q2, q2_tl&
&   , i1, i2, iv, kord)
    IMPLICIT NONE
! !INPUT PARAMETERS:
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! Method order
    INTEGER, INTENT(IN) :: kord
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: iv
! height at layer edges 
    REAL(p_precision), INTENT(IN) :: pe1(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! hieght at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe2_tl(i1:i2, km+1)
! (from model top to bottom surface)
!      real(p_precision), intent(in) ::  q1(i1:i2,km)        ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL(p_precision), INTENT(IN) :: q2(i1:i2, km)
    REAL(p_precision), INTENT(INOUT) :: q2_tl(i1:i2, km)
! !LOCAL VARIABLES:
    REAL(p_precision) :: dp1(i1:i2, km)
    REAL(p_precision) :: dp1_tl(i1:i2, km)
    REAL(p_precision) :: q4(4, i1:i2, km)
    REAL(p_precision) :: q4_tl(4, i1:i2, km)
    REAL(p_precision) :: pl, pr, qsum, delp, esl
    REAL(p_precision) :: pl_tl, pr_tl, qsum_tl, delp_tl, esl_tl
    INTEGER :: i, k, l, m, k0
    dp1_tl = 0.0_8
    q4_tl = 0.0_8
    DO k=1,km
      DO i=i1,i2
! negative
        dp1_tl(i, k) = pe1_tl(i, k+1) - pe1_tl(i, k)
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4_tl(1, i, k) = q2_tl(i, k)
        q4(1, i, k) = q2(i, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE_NOLIMITERS_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv&
&                          , kord)
      qsum_tl = 0.0_8
    ELSE
      CALL PPM_PROFILE_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv, kord)
      qsum_tl = 0.0_8
    END IF
! Mapping
    DO i=i1,i2
      k0 = 1
      DO 555 k=1,km
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .LE. pe1(i, l) .AND. pe2(i, k) .GE. pe1(i, l+1)&
&         ) GOTO 110
        END DO
        GOTO 123
 110    pl_tl = ((pe2_tl(i, k)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k)-pe1(i&
&         , l))*dp1_tl(i, l))/dp1(i, l)**2
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .GE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr_tl = ((pe2_tl(i, k+1)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k+1)-&
&           pe1(i, l))*dp1_tl(i, l))/dp1(i, l)**2
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          q2_tl(i, k) = q4_tl(2, i, l) + 0.5*((q4_tl(4, i, l)+q4_tl(3, i&
&           , l)-q4_tl(2, i, l))*(pr+pl)+(q4(4, i, l)+q4(3, i, l)-q4(2, &
&           i, l))*(pr_tl+pl_tl)) - r3*(q4_tl(4, i, l)*(pr*(pr+pl)+pl**2&
&           )+q4(4, i, l)*(pr_tl*(pr+pl)+pr*(pr_tl+pl_tl)+2*pl*pl_tl))
          k0 = l
          GOTO 555
        ELSE
! Fractional area...
          qsum_tl = (pe1_tl(i, l+1)-pe2_tl(i, k))*(q4(2, i, l)+0.5*(q4(4&
&           , i, l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.&
&           +pl*(1.+pl)))) + (pe1(i, l+1)-pe2(i, k))*(q4_tl(2, i, l)+0.5&
&           *((q4_tl(4, i, l)+q4_tl(3, i, l)-q4_tl(2, i, l))*(1.+pl)+(q4&
&           (4, i, l)+q4(3, i, l)-q4(2, i, l))*pl_tl)-r3*(q4_tl(4, i, l)&
&           *(1.+pl*(1.+pl))+q4(4, i, l)*(pl_tl*(1.+pl)+pl*pl_tl)))
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .LT. pe1(i, m+1)) THEN
! Whole layer..
              qsum_tl = qsum_tl + dp1_tl(i, m)*q4(1, i, m) + dp1(i, m)*&
&               q4_tl(1, i, m)
              qsum = qsum + dp1(i, m)*q4(1, i, m)
            ELSE
              GOTO 120
            END IF
          END DO
          GOTO 123
 120      delp_tl = pe2_tl(i, k+1) - pe1_tl(i, m)
          delp = pe2(i, k+1) - pe1(i, m)
          esl_tl = (delp_tl*dp1(i, m)-delp*dp1_tl(i, m))/dp1(i, m)**2
          esl = delp/dp1(i, m)
          qsum_tl = qsum_tl + delp_tl*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-&
&           q4(2, i, m)+q4(4, i, m)*(1.-r23*esl))) + delp*(q4_tl(2, i, m&
&           )+0.5*(esl_tl*(q4(3, i, m)-q4(2, i, m)+q4(4, i, m)*(1.-r23*&
&           esl))+esl*(q4_tl(3, i, m)-q4_tl(2, i, m)+q4_tl(4, i, m)*(1.-&
&           r23*esl)-q4(4, i, m)*r23*esl_tl)))
          qsum = qsum + delp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, &
&           m)+q4(4, i, m)*(1.-r23*esl)))
          k0 = m
        END IF
 123    q2_tl(i, k) = (qsum_tl*(pe2(i, k+1)-pe2(i, k))-qsum*(pe2_tl(i, k&
&         +1)-pe2_tl(i, k)))/(pe2(i, k+1)-pe2(i, k))**2
 555  CONTINUE
    END DO
  END SUBROUTINE REMAP_Z_NOLIMITERS_TLM

  SUBROUTINE MAP1_PPM_NOLIMITERS_TLM(km, pe1, pe1_tl, pe2, pe2_tl, q2, q2_tl&
&   , i1, i2, j, ibeg, iend, jbeg, jend, iv, kord)
    IMPLICIT NONE
! Starting longitude
    INTEGER, INTENT(IN) :: i1
! Finishing longitude
    INTEGER, INTENT(IN) :: i2
! Mode: 0 == constituents  1 == ???
    INTEGER, INTENT(IN) :: iv
!       2 == remap temp with cs scheme
! Method order
    INTEGER, INTENT(IN) :: kord
! Current latitude
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: ibeg, iend, jbeg, jend
! Original vertical dimension
    INTEGER, INTENT(IN) :: km
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe1(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe1_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, km+1)
    REAL(p_precision), INTENT(IN) :: pe2_tl(i1:i2, km+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(IN) :: q2(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(INOUT) :: q2_tl(ibeg:iend, jbeg:jend, km)
! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
    REAL(p_precision) :: dp1(i1:i2, km)
    REAL(p_precision) :: dp1_tl(i1:i2, km)
    REAL(p_precision) :: q4(4, i1:i2, km)
    REAL(p_precision) :: q4_tl(4, i1:i2, km)
    REAL(p_precision) :: pl, pr, qsum, dp, esl
    REAL(p_precision) :: pl_tl, pr_tl, qsum_tl, dp_tl, esl_tl
    INTEGER :: i, k, l, m, k0
    dp1_tl = 0.0_8
    q4_tl = 0.0_8
    DO k=1,km
      DO i=i1,i2
        dp1_tl(i, k) = pe1_tl(i, k+1) - pe1_tl(i, k)
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4_tl(1, i, k) = q2_tl(i, j, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL CS_PROFILE_NOLIMITERS_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv&
&                          , kord)
      qsum_tl = 0.0_8
    ELSE
      CALL PPM_PROFILE_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2, iv, kord)
      qsum_tl = 0.0_8
    END IF
    DO i=i1,i2
      k0 = 1
      DO 555 k=1,km
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) GOTO 100
        END DO
        GOTO 123
 100    pl_tl = ((pe2_tl(i, k)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k)-pe1(i&
&         , l))*dp1_tl(i, l))/dp1(i, l)**2
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          pr_tl = ((pe2_tl(i, k+1)-pe1_tl(i, l))*dp1(i, l)-(pe2(i, k+1)-&
&           pe1(i, l))*dp1_tl(i, l))/dp1(i, l)**2
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          q2_tl(i, j, k) = q4_tl(2, i, l) + 0.5*((q4_tl(4, i, l)+q4_tl(3&
&           , i, l)-q4_tl(2, i, l))*(pr+pl)+(q4(4, i, l)+q4(3, i, l)-q4(&
&           2, i, l))*(pr_tl+pl_tl)) - r3*(q4_tl(4, i, l)*(pr*(pr+pl)+pl&
&           **2)+q4(4, i, l)*(pr_tl*(pr+pl)+pr*(pr_tl+pl_tl)+2*pl*pl_tl)&
&           )
          k0 = l
          GOTO 555
        ELSE
! Fractional area...
          qsum_tl = (pe1_tl(i, l+1)-pe2_tl(i, k))*(q4(2, i, l)+0.5*(q4(4&
&           , i, l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.&
&           +pl*(1.+pl)))) + (pe1(i, l+1)-pe2(i, k))*(q4_tl(2, i, l)+0.5&
&           *((q4_tl(4, i, l)+q4_tl(3, i, l)-q4_tl(2, i, l))*(1.+pl)+(q4&
&           (4, i, l)+q4(3, i, l)-q4(2, i, l))*pl_tl)-r3*(q4_tl(4, i, l)&
&           *(1.+pl*(1.+pl))+q4(4, i, l)*(pl_tl*(1.+pl)+pl*pl_tl)))
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
              qsum_tl = qsum_tl + dp1_tl(i, m)*q4(1, i, m) + dp1(i, m)*&
&               q4_tl(1, i, m)
              qsum = qsum + dp1(i, m)*q4(1, i, m)
            ELSE
              GOTO 110
            END IF
          END DO
          GOTO 123
 110      dp_tl = pe2_tl(i, k+1) - pe1_tl(i, m)
          dp = pe2(i, k+1) - pe1(i, m)
          esl_tl = (dp_tl*dp1(i, m)-dp*dp1_tl(i, m))/dp1(i, m)**2
          esl = dp/dp1(i, m)
          qsum_tl = qsum_tl + dp_tl*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4&
&           (2, i, m)+q4(4, i, m)*(1.-r23*esl))) + dp*(q4_tl(2, i, m)+&
&           0.5*(esl_tl*(q4(3, i, m)-q4(2, i, m)+q4(4, i, m)*(1.-r23*esl&
&           ))+esl*(q4_tl(3, i, m)-q4_tl(2, i, m)+q4_tl(4, i, m)*(1.-r23&
&           *esl)-q4(4, i, m)*r23*esl_tl)))
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
        END IF
 123    q2_tl(i, j, k) = (qsum_tl*(pe2(i, k+1)-pe2(i, k))-qsum*(pe2_tl(i&
&         , k+1)-pe2_tl(i, k)))/(pe2(i, k+1)-pe2(i, k))**2
 555  CONTINUE
    END DO
  END SUBROUTINE MAP1_PPM_NOLIMITERS_TLM

  SUBROUTINE CS_PROFILE_NOLIMITERS_TLM(a4, a4_tl, delp, delp_tl, km, i1, i2&
&   , iv, kord)
    IMPLICIT NONE
! A perfectly linear version of the vertical remapping algorithm
! Warning: not for general purpose; non-monotonic
! Latest: 20141211 S.-J. Lin, NOAA/GFDL
    INTEGER, INTENT(IN) :: i1, i2
! dummy argument
    INTEGER, INTENT(IN) :: kord
! vertical dimension
    INTEGER, INTENT(IN) :: km
! iv =-1: winds
    INTEGER, INTENT(IN) :: iv
! iv = 0: positive definite scalars
! iv = 1: others
! layer pressure thickness
    REAL, INTENT(IN) :: delp(i1:i2, km)
    REAL, INTENT(IN) :: delp_tl(i1:i2, km)
! Interpolated values
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL, INTENT(INOUT) :: a4_tl(4, i1:i2, km)
!-----------------------------------------------------------------------
    REAL :: gam(i1:i2, km)
    REAL :: gam_tl(i1:i2, km)
    REAL :: q(i1:i2, km+1)
    REAL :: q_tl(i1:i2, km+1)
    REAL :: d4(i1:i2)
    REAL :: d4_tl(i1:i2)
    REAL :: bet, a_bot, grat
    REAL :: bet_tl, a_bot_tl, grat_tl
    INTEGER :: i, k
    q_tl = 0.0
    gam_tl = 0.0
    DO i=i1,i2
! grid ratio
      grat_tl = (delp_tl(i, 2)*delp(i, 1)-delp(i, 2)*delp_tl(i, 1))/delp&
&       (i, 1)**2
      grat = delp(i, 2)/delp(i, 1)
      bet_tl = grat_tl*(grat+0.5) + grat*grat_tl
      bet = grat*(grat+0.5)
      q_tl(i, 1) = (((2*grat_tl*(grat+1.)+(grat+grat)*grat_tl)*a4(1, i, &
&       1)+(grat+grat)*(grat+1.)*a4_tl(1, i, 1)+a4_tl(1, i, 2))*bet-((&
&       grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))*bet_tl)/bet**2
      q(i, 1) = ((grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))/bet
      gam_tl(i, 1) = ((grat_tl*(grat+1.5)+grat*grat_tl)*bet-(1.+grat*(&
&       grat+1.5))*bet_tl)/bet**2
      gam(i, 1) = (1.+grat*(grat+1.5))/bet
    END DO
    d4_tl = 0.0
    DO k=2,km
      DO i=i1,i2
        d4_tl(i) = (delp_tl(i, k-1)*delp(i, k)-delp(i, k-1)*delp_tl(i, k&
&         ))/delp(i, k)**2
        d4(i) = delp(i, k-1)/delp(i, k)
        bet_tl = 2*d4_tl(i) - gam_tl(i, k-1)
        bet = 2. + d4(i) + d4(i) - gam(i, k-1)
        q_tl(i, k) = ((3.*(a4_tl(1, i, k-1)+d4_tl(i)*a4(1, i, k)+d4(i)*&
&         a4_tl(1, i, k))-q_tl(i, k-1))*bet-(3.*(a4(1, i, k-1)+d4(i)*a4(&
&         1, i, k))-q(i, k-1))*bet_tl)/bet**2
        q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
        gam_tl(i, k) = (d4_tl(i)*bet-d4(i)*bet_tl)/bet**2
        gam(i, k) = d4(i)/bet
      END DO
    END DO
    DO i=i1,i2
      a_bot_tl = d4_tl(i)*(d4(i)+1.5) + d4(i)*d4_tl(i)
      a_bot = 1. + d4(i)*(d4(i)+1.5)
      q_tl(i, km+1) = ((2.*((d4_tl(i)*(d4(i)+1.)+d4(i)*d4_tl(i))*a4(1, i&
&       , km)+d4(i)*(d4(i)+1.)*a4_tl(1, i, km))+a4_tl(1, i, km-1)-&
&       a_bot_tl*q(i, km)-a_bot*q_tl(i, km))*(d4(i)*(d4(i)+0.5)-a_bot*&
&       gam(i, km))-(2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&       a_bot*q(i, km))*(d4_tl(i)*(d4(i)+0.5)+d4(i)*d4_tl(i)-a_bot_tl*&
&       gam(i, km)-a_bot*gam_tl(i, km)))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i&
&       , km))**2
      q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&       a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
    END DO
    DO k=km,1,-1
      DO i=i1,i2
        q_tl(i, k) = q_tl(i, k) - gam_tl(i, k)*q(i, k+1) - gam(i, k)*&
&         q_tl(i, k+1)
        q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
      END DO
    END DO
    DO k=1,km
      DO i=i1,i2
        a4_tl(2, i, k) = q_tl(i, k)
        a4(2, i, k) = q(i, k)
        a4_tl(3, i, k) = q_tl(i, k+1)
        a4(3, i, k) = q(i, k+1)
        a4_tl(4, i, k) = 3.*(2.*a4_tl(1, i, k)-a4_tl(2, i, k)-a4_tl(3, i&
&         , k))
        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
      END DO
    END DO
  END SUBROUTINE CS_PROFILE_NOLIMITERS_TLM

 subroutine remap_z(km, pe1, pe2, q2, i1, i2, iv, kord)

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: iv

      real(p_precision), intent(in) ::  pe1(i1:i2,km+1)     ! height at layer edges 
                                               ! (from model top to bottom surface)
      real(p_precision), intent(in) ::  pe2(i1:i2,km+1)     ! hieght at layer edges 
                                               ! (from model top to bottom surface)
!      real(p_precision), intent(in) ::  q1(i1:i2,km)        ! Field input

! !INPUT/OUTPUT PARAMETERS:
      real(p_precision), intent(inout)::  q2(i1:i2,km)      ! Field output

! !LOCAL VARIABLES:
      real(p_precision)  dp1(  i1:i2,km)
      real(p_precision)   q4(4,i1:i2,km)
      real(p_precision)   pl, pr, qsum, delp, esl
      integer i, k, l, m, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)      ! negative
            q4(1,i,k) = q2(i,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,km
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) <= pe1(i,l) .and. pe2(i,k) >= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) >= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) < pe1(i,m+1) ) then
! Whole layer..
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                    delp = pe2(i,k+1)-pe1(i,m)
                    esl = delp / dp1(i,m)
                    qsum = qsum + delp*(q4(2,i,m)+0.5*esl*               &
                         (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                    k0 = m
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue

 end subroutine remap_z

 subroutine map1_ppm( km,   pe1,                        &
                            pe2,    q2,   i1, i2,       &
                      j,    ibeg, iend, jbeg, jend, iv,  kord)
 integer, intent(in) :: i1                ! Starting longitude
 integer, intent(in) :: i2                ! Finishing longitude
 integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == ???
                                          !       2 == remap temp with cs scheme
 integer, intent(in) :: kord              ! Method order
 integer, intent(in) :: j                 ! Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                ! Original vertical dimension
 real(p_precision), intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the original vertical coordinate
 real(p_precision), intent(in) ::  pe2(i1:i2,km+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the new vertical coordinate
! real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout)::  q2(ibeg:iend,jbeg:jend,km) ! Field output

! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
   real(p_precision)    dp1(i1:i2,km)
   real(p_precision)   q4(4,i1:i2,km)
   real(p_precision)    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q2(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,km
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map1_ppm

 subroutine map1_q2(km,   pe1,   q1,            &
                    kn,   pe2,   q2,   dp2,     &
                    i1,   i2,    iv,   kord, j, &
                    ibeg, iend, jbeg, jend )


! !INPUT PARAMETERS:
      integer, intent(in) :: j
      integer, intent(in) :: i1, i2
      integer, intent(in) :: ibeg, iend, jbeg, jend
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents 1 == ???
      integer, intent(in) :: kord
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real(p_precision), intent(in) ::  pe1(i1:i2,km+1)     ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real(p_precision), intent(in) ::  pe2(i1:i2,kn+1)     ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real             , intent(in) ::  q1(ibeg:iend,jbeg:jend,km) ! Field input
      real(p_precision), intent(in) ::  dp2(i1:i2,kn)
! !INPUT/OUTPUT PARAMETERS:
      real(p_precision), intent(inout):: q2(i1:i2,kn) ! Field output
! !LOCAL VARIABLES:
      real(p_precision)   dp1(i1:i2,km)
      real(p_precision)   q4(4,i1:i2,km)
      real(p_precision)   pl, pr, qsum, dp, esl

      integer i, k, l, m, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                    esl = dp / dp1(i,m)
                   qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                       (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                   k0 = m
                   goto 123
                 endif
              enddo
              goto 123
          endif
      endif
100   continue
123   q2(i,k) = qsum / dp2(i,k)
555   continue
1000  continue

 end subroutine map1_q2

 subroutine cs_profile(a4, delp, km, i1, i2, iv, kord)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 real(p_precision), intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real(p_precision), intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
 integer, intent(in):: kord
!-----------------------------------------------------------------------
 logical:: extm(i1:i2,km)
 real(p_precision)  gam(i1:i2,km)
 real(p_precision)    q(i1:i2,km+1)
 real(p_precision)   d4(i1:i2)
 real(p_precision)   bet, a_bot, grat, pmp, lac
 real(p_precision)   pmp_1, lac_1, pmp_2, lac_2
 integer i, k, im

  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo
 
  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo

!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints 
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
     enddo
  enddo

! Interior:
  if ( abs(kord)<11 ) then
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo
  else
! abs(kord) >=11
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1) > 0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
               q(i,k) = min(q(i,k), min(a4(1,i,k-1)+0.5*gam(i,k-1),a4(1,i,k  )-0.5*gam(i,k+1) ) )
          else
! There exists a local min
            if ( iv==0 ) then
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
                 q(i,k) = max(q(i,k), min(a4(1,i,k-1)+0.5*gam(i,k-1),     &
                                  min(0., a4(1,i,k  )-0.5*gam(i,k+1)) ) )
            else
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
                 q(i,k) = max(q(i,k), max(a4(1,i,k-1)+0.5*gam(i,k-1), a4(1,i,k  )-0.5*gam(i,k+1) ))
            endif
          endif
        endif
     enddo
  enddo
  endif

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=2,km-1
     do i=i1,i2
        if ( gam(i,k)*gam(i,k+1) > 0.0 ) then
             extm(i,k) = .false. 
        else
             extm(i,k) = .true.
        endif
     enddo
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then 
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( abs(iv)==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( abs(iv)/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. (extm(i,k-1).or.extm(i,k+1)) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( extm(i,k) ) then
              if( extm(i,k-1) .or. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
                   a4(4,i,k) = 0.
              else
! True local extremum
                a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
              endif
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     else      ! kord = 11, ...
       do i=i1,i2
         if ( extm(i,k) .and. (extm(i,k-1) .or. extm(i,k+1)) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv<0 ) then 
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine cs_profile

 subroutine cs_limiters(im, extm, a4, iv)
 integer, intent(in) :: im
 integer, intent(in) :: iv
 logical, intent(in) :: extm(im)
 real(p_precision) , intent(inout) :: a4(4,im)   ! PPM array
! !LOCAL VARIABLES:
 real(p_precision)  da1, da2, a6da
 real(p_precision)  fmin
 integer i

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
    if( a4(1,i)<=0.) then
        a4(2,i) = a4(1,i)
        a4(3,i) = a4(1,i)
        a4(4,i) = 0.
    else
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
         if( (a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12) < 0. ) then
! local minimum is negative
             if( a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i) ) then
                 a4(3,i) = a4(1,i)
                 a4(2,i) = a4(1,i)
                 a4(4,i) = 0.
             elseif( a4(3,i) > a4(2,i) ) then
                 a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                 a4(3,i) = a4(2,i) - a4(4,i)
             else
                 a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                 a4(2,i) = a4(3,i) - a4(4,i)
             endif
         endif
      endif
    endif
    enddo
 elseif ( iv==1 ) then
    do i=1,im
      if( (a4(1,i)-a4(2,i))*(a4(1,i)-a4(3,i))>=0. ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
      if( extm(i) ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 endif
 end subroutine cs_limiters

 subroutine ppm_profile(a4, delp, km, i1, i2, iv, kord)

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
                               ! iv = 2: temp (if remap_t) and w (iv=-2)
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real(p_precision) , intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real(p_precision) , intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
! 
! !REVISION HISTORY: 
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
      real(p_precision)    dc(i1:i2,km)
      real(p_precision)    h2(i1:i2,km)
      real(p_precision)  delq(i1:i2,km)
      real(p_precision)   df2(i1:i2,km)
      real(p_precision)    d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt, it
      real(p_precision)  fac
      real(p_precision)  a1, a2, c1, c2, c3, d1, d2
      real(p_precision)  qm, dq, lac, qmp, pmp

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
                 c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
                 c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            df2(i,k) = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            dc(i,k) = sign( min(abs(df2(i,k)),              &
                            max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))-a4(1,i,k),  &
                  a4(1,i,k)-min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))), df2(i,k) )
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=3,km1
         do i=i1,i2
            c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
            a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
            a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
            a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                      ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                        delp(i,k-1)*a1*dc(i,k  ) )
         enddo
      enddo

      if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
!-------------------------------------------------------
!        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
!-------------------------------------------------------
! No over- and undershoot condition
         a4(2,i,2) = max( a4(2,i,2), min(a4(1,i,1), a4(1,i,2)) )
         a4(2,i,2) = min( a4(2,i,2), max(a4(1,i,1), a4(1,i,2)) )
         dc(i,1) =  0.5*(a4(2,i,2) - a4(1,i,1))
      enddo

! Enforce monotonicity  within the top layer

      if( iv==0 ) then
         do i=i1,i2
            a4(2,i,1) = max(0., a4(2,i,1))
            a4(2,i,2) = max(0., a4(2,i,2))
         enddo 
      elseif( iv==-1 ) then
         do i=i1,i2
            if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
         enddo
      elseif( abs(iv)==2 ) then
         do i=i1,i2
            a4(2,i,1) = a4(1,i,1)
            a4(3,i,1) = a4(1,i,1)
         enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
!        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
!-----------------------------------------------------
!        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
! No over- and under-shoot condition
         a4(2,i,km) = max( a4(2,i,km), min(a4(1,i,km), a4(1,i,km1)) )
         a4(2,i,km) = min( a4(2,i,km), max(a4(1,i,km), a4(1,i,km1)) )
         dc(i,km) = 0.5*(a4(1,i,km) - a4(2,i,km))
      enddo


! Enforce constraint on the "slope" at the surface

      if( iv==0 ) then
          do i=i1,i2
             a4(2,i,km) = max(0.,a4(2,i,km))
             a4(3,i,km) = max(0.,a4(3,i,km))
          enddo
      elseif( iv<0 ) then
          do i=i1,i2
             if( a4(1,i,km)*a4(3,i,km) <= 0. )  a4(3,i,km) = 0.
          enddo
      endif

   do k=1,km1
      do i=i1,i2
         a4(3,i,k) = a4(2,i,k+1)
      enddo
   enddo

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=2,km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
            h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))  &
                     / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )        &
                     * delp(i,k)**2 
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      fac = 1.5           ! original quasi-monotone

      do k=3,km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), qmp, lac)),    &
                                        max(a4(1,i,k), qmp, lac) )
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         a4(2,i,k) = min(max(a4(2,i,k),  min(a4(1,i,k), qmp, lac)),   &
                     max(a4(1,i,k), qmp, lac))
!-------------
! Recompute A6
!-------------
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord >= 6 )                      &
             call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 2)
      enddo

      else

         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

         do k=3,km-2
            if( kord /= 4) then
              do i=i1,i2
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              enddo
            endif
            if(kord/=6) call ppm_limiters(dc(i1,k), a4(1,i1,k), it, lmt)
         enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

 end subroutine ppm_profile


 subroutine ppm_limiters(dm, a4, itot, lmt)

! !INPUT PARAMETERS:
      real(p_precision) , intent(in):: dm(*)     ! the linear slope
      integer, intent(in) :: itot      ! Total Longitudes
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real(p_precision) , intent(inout) :: a4(4,*)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)
! !LOCAL VARIABLES:
      real(p_precision)  qmp
      real(p_precision)  da1, da2, a6da
      real(p_precision)  fmin
      integer i

! Developer: S.-J. Lin

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

 end subroutine ppm_limiters

 subroutine steepz(i1, i2, km, a4, df2, dm, dq, dp, d4)
 integer, intent(in) :: km, i1, i2
   real(p_precision) , intent(in) ::  dp(i1:i2,km)       ! grid size
   real(p_precision) , intent(in) ::  dq(i1:i2,km)       ! backward diff of q
   real(p_precision) , intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
   real(p_precision) , intent(in) :: df2(i1:i2,km)       ! first guess mismatch
   real(p_precision) , intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch
! !INPUT/OUTPUT PARAMETERS:
      real(p_precision) , intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened
! !LOCAL VARIABLES:
      integer i, k
      real(p_precision)  alfa(i1:i2,km)
      real(p_precision)     f(i1:i2,km)
      real(p_precision)   rat(i1:i2,km)
      real(p_precision)   dg2

! Compute ratio of dq/dp
      do k=2,km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=2,km-1
         do i=i1,i2
            f(i,k) =   (rat(i,k+1) - rat(i,k))                          &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=3,km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<0. .and. df2(i,k)/=0.) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(0., min(0.5, -0.1875*dg2/df2(i,k)))
         else
            alfa(i,k) = 0.
         endif
         enddo
      enddo

      do k=4,km-2
         do i=i1,i2
            a4(2,i,k) = (1.-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

 end subroutine steepz

end module fv_mapz_tlm_mod
