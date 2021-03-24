module fv_mapz_tlm_mod

  use fv_arrays_mod,      only: p_precision
!#ifndef MAPL_MODE
!  use constants_mod, only: radius, pi, rvgas, rdgas, grav
!#endif
  use fv_grid_tools_mod,  only: area, dx, dy, rdxa, rdya
  use fv_grid_utils_mod,  only: g_sum, ptop, ptop_min
  use fv_grid_utils_mod,  only: cosa_s, rsin2
  use fv_fill_mod,        only: fillz
  use fv_fill_tlm_mod,    only: fillz_tlm
  use fv_mp_mod,          only: gid, domain
  use mpp_domains_mod,    only: mpp_update_domains

  implicit none

  real(p_precision), parameter::  D0_5 = 0.5, D0_0 = 0.0
  real(p_precision), parameter::  r3 = 1./3., r23 = 2./3., r12 = 1./12.
  real(kind=4) :: E_Flux
  private

  public compute_total_energy_tlm, Lagrangian_to_Eulerian_tlm

!      INTERFACE mappm
!        MODULE PROCEDURE mappm_r4
!        MODULE PROCEDURE mappm_r8
!      END INTERFACE
!
!      INTERFACE ana_remap
!        MODULE PROCEDURE ana_remap_
!      END INTERFACE

CONTAINS
!  Differentiation of lagrangian_to_eulerian in forward (tangent) mode (with options r8):
!   variations   of useful results: peln q u v delp pkz pe pt
!   with respect to varying inputs: peln q u v delp ua te0_2d pkz
!                pe pk pt te
!#ifdef MAPL_MODE
!#endif
  SUBROUTINE LAGRANGIAN_TO_EULERIAN_TLM(do_consv, consv, pe, pe_tl, delp&
&   , delp_tl, pkz, pkz_tl, pk, pk_tl, pdt, km, is, ie, js, je, isd, ied&
&   , jsd, jed, nq, sphum, u, u_tl, v, v_tl, pt, pt_tl, q, q_tl, hs, &
&   r_vir, cp, akap, pi, radius, grav, kord_mt, kord_tr, kord_tm, peln, &
&   peln_tl, te0_2d, te0_2d_tl, ng, ua, ua_tl, te, te_tl, pem, fill, &
&   reproduce_sum, ak, bk, ks, te_method, remap_t, ncnst)
    IMPLICIT NONE
! end k-loop
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
    INTEGER, INTENT(IN) :: ks
! Mapping oder for the vector winds
    INTEGER, INTENT(IN) :: kord_mt
! Mapping oder for tracers
    INTEGER, INTENT(IN) :: kord_tr
! Mapping oder for thermodynamics
    INTEGER, INTENT(IN) :: kord_tm
! Added 4TAF
    INTEGER, INTENT(IN) :: ncnst
! factor for TE conservation
    REAL, INTENT(IN) :: consv
    REAL, INTENT(IN) :: r_vir
    REAL, INTENT(IN) :: cp
    REAL, INTENT(IN) :: akap
!#ifdef MAPL_MODE
    REAL, INTENT(IN) :: pi, radius, grav
!#endif
! surface geopotential
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: te0_2d(is:ie, js:je)
    REAL, INTENT(IN) :: te0_2d_tl(is:ie, js:je)
! fill negative tracers
    LOGICAL, INTENT(IN) :: fill
    LOGICAL, INTENT(IN) :: reproduce_sum
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
! u-wind will be ghosted one latitude to the north upon exit
! u-wind (m/s)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: u_tl(isd:ied, jsd:jed+1, km)
! v-wind (m/s)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_tl(isd:ied+1, jsd:jed, km)
! cp*virtual potential temperature 
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: pt_tl(isd:ied, jsd:jed, km)
! as input; output: temperature
    INTEGER, INTENT(IN) :: te_method
    LOGICAL, INTENT(IN) :: remap_t
! u-wind (m/s) on physics grid
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: ua_tl(isd:ied, jsd:jed, km)
! log(pe)
    REAL(p_precision), INTENT(INOUT) :: peln(is:ie, km+1, js:je)
    REAL(p_precision), INTENT(INOUT) :: peln_tl(is:ie, km+1, js:je)
! layer-mean pk for converting t to pt
    REAL(p_precision), INTENT(OUT) :: pkz(is:ie, js:je, km)
    REAL(p_precision), INTENT(OUT) :: pkz_tl(is:ie, js:je, km)
    REAL, INTENT(OUT) :: te(is:ie, js:je, km)
    REAL, INTENT(OUT) :: te_tl(is:ie, js:je, km)
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
    INTRINSIC REAL
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
    te_map = .true.
    CALL PKEZ_TLM(km, is, ie, js, je, pe, pk, pk_tl, akap, peln, peln_tl&
&           , pkz, pkz_tl)
    DO k=1,km
      DO j=js,je
        DO i=is,ie
          te_tl(i, j, k) = 0.25*rsin2(i, j)*(2*u(i, j, k)*u_tl(i, j, k)+&
&           2*u(i, j+1, k)*u_tl(i, j+1, k)+2*v(i, j, k)*v_tl(i, j, k)+2*&
&           v(i+1, j, k)*v_tl(i+1, j, k)-cosa_s(i, j)*((u_tl(i, j, k)+&
&           u_tl(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))+(u(i, j, k)+u(i, &
&           j+1, k))*(v_tl(i, j, k)+v_tl(i+1, j, k)))) + pt_tl(i, j, k)*&
&           pkz(i, j, k) + pt(i, j, k)*pkz_tl(i, j, k)
          te(i, j, k) = 0.25*rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+&
&           v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i&
&           , j, k)+v(i+1, j, k))*cosa_s(i, j)) + pt(i, j, k)*pkz(i, j, &
&           k)
        END DO
      END DO
    END DO
    phis_tl = 0.0_8
    pe0_tl = 0.0_8
    pe1_tl = 0.0_8
    pe2_tl = 0.0_8
    pe3_tl = 0.0_8
    dp2_tl = 0.0_8
    q2_tl = 0.0_8
    pn2_tl = 0.0_8
    pk1_tl = 0.0_8
    pk2_tl = 0.0_8
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
!------------
! update delp
!------------
        DO k=1,km
          DO i=is,ie
            delp_tl(i, j, k) = dp2_tl(i, k)
            delp(i, j, k) = dp2(i, k)
          END DO
        END DO
!-----------------
! Map constituents
!-----------------
        IF (nq .NE. 0) THEN
          DO iq=1,nq
            CALL MAP1_Q2_NOLIMITERS_TLM(km, pe1, pe1_tl, q(isd, jsd, 1, &
&                                 iq), q_tl(isd, jsd, 1, iq), km, pe2, &
&                                 pe2_tl, q2, q2_tl, dp2, dp2_tl, is, ie&
&                                 , 0, kord_tr, j, isd, ied, jsd, jed)
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
          CALL MAP1_CUBIC_TE_TLM(km, pe1, pe1_tl, pe2, pe2_tl, te, te_tl&
&                          , is, ie, j, is, ie, js, je, 1, kord_tm)
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
        DO k=1,km
          DO i=is,ie
            pkz_tl(i, j, k) = ((pk2_tl(i, k+1)-pk2_tl(i, k))*akap*(peln(&
&             i, k+1, j)-peln(i, k, j))-(pk2(i, k+1)-pk2(i, k))*akap*(&
&             peln_tl(i, k+1, j)-peln_tl(i, k, j)))/(akap*(peln(i, k+1, &
&             j)-peln(i, k, j)))**2
            pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1, j&
&             )-peln(i, k, j)))
          END DO
        END DO
      END IF
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
      CALL MAP1_PPM_NOLIMITERS_TLM(km, pe0(is:ie, :), pe0_tl(is:ie, :), &
&                            pe3(is:ie, :), pe3_tl(is:ie, :), u, u_tl, &
&                            is, ie, j, isd, ied, jsd, jed + 1, -1, &
&                            kord_mt)
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
        CALL MAP1_PPM_NOLIMITERS_TLM(km, pe0, pe0_tl, pe3, pe3_tl, v, &
&                              v_tl, is, ie + 1, j, isd, ied + 1, jsd, &
&                              jed, -1, kord_mt)
      END IF
      DO k=1,km
        DO i=is,ie
          ua_tl(i, j, k) = pe2_tl(i, k+1)
          ua(i, j, k) = pe2(i, k+1)
        END DO
      END DO
    END DO
    DO k=2,km
      DO j=js,je
        DO i=is,ie
          pe_tl(i, k, j) = ua_tl(i, j, k-1)
          pe(i, k, j) = ua(i, j, k-1)
        END DO
      END DO
    END DO
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
        zsum0_tl = 0.0_8
        zsum1_tl = 0.0_8
      ELSE
        te_2d_tl = 0.0_8
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
!on tpe = consv*g_sum(te_2d, is, ie, js, je, ng, area, 0)
!on tmpsum = g_sum(zsum0,  is, ie, js, je, ng, area, 0)
!on dtmp = tpe / (cp*tmpsum)

      tpe = consv*g_sum(te_2d, is, ie, js, je, ng, area, 0)
      tmpsum = g_sum(zsum0,  is, ie, js, je, ng, area, 0)
      tpe_tl = consv*g_sum(te_2d_tl, is, ie, js, je, ng, area, 0)
      tmpsum_tl = g_sum(zsum0_tl,  is, ie, js, je, ng, area, 0)

      dtmp_tl = (tpe_tl*cp*tmpsum-tpe*cp*tmpsum_tl)/(cp*tmpsum)**2
      dtmp = tpe/(cp*tmpsum)
! convert to 4-byte real
      IF (reproduce_sum) dtmp = REAL(dtmp, 4)
    ELSE
      dtmp = 0.
      dtmp_tl = 0.0_8
    END IF
    IF (te_map) THEN
      gz_tl = 0.0_8
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
    END IF
  END SUBROUTINE LAGRANGIAN_TO_EULERIAN_TLM
  SUBROUTINE COMPUTE_TOTAL_ENERGY_TLM(is, ie, js, je, isd, ied, jsd, jed&
&   , km, u, u_tl, v, v_tl, w, delz, delz_tl, pt, pt_tl, delp, delp_tl, &
&   q, q_tl, pe, pe_tl, peln, peln_tl, hs, grav, r_vir, cp, rg, hlv, &
&   te_2d, te_2d_tl, ua, teq, moist_phys, sphum, hydrostatic, id_te)
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
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(INOUT) :: ua
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
!#ifdef MAPL_MODE                 
    REAL, INTENT(IN) :: grav
!#endif
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
    te_2d_tl = 0.0_8
    phiz_tl = 0.0_8
    tv_tl = 0.0_8
!----------------------
! Output lat-lon winds:
!----------------------
!  call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km)
!$omp parallel do default(shared) private(phiz, tv)
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
          te_2d_tl(i, j) = 0.0_8
          te_2d(i, j) = 0.
        END DO
        DO k=1,km
          DO i=is,ie
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
!$omp parallel do default(shared) 
      DO j=js,je
        DO i=is,ie
          teq(i, j) = te_2d(i, j)
        END DO
      END DO
      IF (moist_phys) THEN
!$omp parallel do default(shared) 
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
!$omp parallel do default(shared) private(lnp, pek, pk2)
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
  SUBROUTINE MAP1_PPM_NOLIMITERS_TLM(km, pe1, pe1_tl, pe2, pe2_tl, q2, &
&   q2_tl, i1, i2, j, ibeg, iend, jbeg, jend, iv, kord)
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
      CALL CS_PROFILE_NOLIMITERS_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2&
&                              , iv, kord)
      qsum_tl = 0.0_8
    ELSE
      qsum_tl = 0.0_8
    END IF
!call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
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
  END SUBROUTINE MAP1_PPM_NOLIMITERS_TLM
  SUBROUTINE MAP1_Q2_NOLIMITERS_TLM(km, pe1, pe1_tl, q1, q1_tl, kn, pe2&
&   , pe2_tl, q2, q2_tl, dp2, dp2_tl, i1, i2, iv, kord, j, ibeg, iend, &
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
      CALL CS_PROFILE_NOLIMITERS_TLM(q4, q4_tl, dp1, dp1_tl, km, i1, i2&
&                              , iv, kord)
      qsum_tl = 0.0_8
    ELSE
      qsum_tl = 0.0_8
    END IF
!call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
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
  END SUBROUTINE MAP1_Q2_NOLIMITERS_TLM
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
  SUBROUTINE CS_PROFILE_NOLIMITERS_TLM(a4, a4_tl, delp, delp_tl, km, i1&
&   , i2, iv, kord)
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

end module fv_mapz_tlm_mod
