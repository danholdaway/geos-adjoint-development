module fv_mapz_adm_mod

  use fv_arrays_mod,      only: p_precision
!#ifndef MAPL_MODE
!  use constants_mod, only: radius, pi, rvgas, rdgas, grav
!#endif
  use fv_grid_tools_mod,  only: area, dx, dy, rdxa, rdya
  use fv_grid_utils_mod,  only: g_sum, ptop, ptop_min
  use fv_grid_utils_mod,  only: cosa_s, rsin2
  use fv_fill_adm_mod,    only: fillz_adm
  use fv_mp_mod,          only: gid, domain
  use mpp_domains_mod,    only: mpp_update_domains

  implicit none

  real(p_precision), parameter::  D0_5 = 0.5, D0_0 = 0.0
  real(p_precision), parameter::  r3 = 1./3., r23 = 2./3., r12 = 1./12.
  real(kind=4) :: E_Flux
  private

  public LAGRANGIAN_TO_EULERIAN, LAGRANGIAN_TO_EULERIAN_ADM, &
         compute_total_energy, compute_total_energy_adm

!      INTERFACE mappm
!        MODULE PROCEDURE mappm_r4
!        MODULE PROCEDURE mappm_r8
!      END INTERFACE
!
!      INTERFACE ana_remap
!        MODULE PROCEDURE ana_remap_
!      END INTERFACE

CONTAINS

  SUBROUTINE LAGRANGIAN_TO_EULERIAN_ADM(do_consv, consv, pe, pe_ad, delp&
&   , delp_ad, pkz, pkz_ad, pk, pk_ad, pdt, km, is, ie, js, je, isd, ied&
&   , jsd, jed, nq, sphum, u, u_ad, v, v_ad, pt, pt_ad, q, q_ad, hs, &
&   r_vir, cp, akap, pi, radius, grav, kord_mt, kord_tr, kord_tm, peln, &
&   peln_ad, te0_2d, te0_2d_ad, ng, ua, ua_ad, te, te_ad, pem, fill, &
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
    REAL :: te0_2d_ad(is:ie, js:je)
! fill negative tracers
    LOGICAL, INTENT(IN) :: fill
    LOGICAL, INTENT(IN) :: reproduce_sum
    REAL, INTENT(IN) :: ak(km+1)
    REAL, INTENT(IN) :: bk(km+1)
! !INPUT/OUTPUT
! pe to the kappa
    REAL(p_precision), INTENT(INOUT) :: pk(is:ie, js:je, km+1)
    REAL(p_precision) :: pk_ad(is:ie, js:je, km+1)
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, km, *)
    REAL, INTENT(INOUT) :: q_ad(isd:ied, jsd:jed, km, *)
! pressure thickness
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: delp_ad(isd:ied, jsd:jed, km)
! pressure at layer edges
    REAL(p_precision), INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL(p_precision) :: pe_ad(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(INOUT) :: pem(is-1:ie+1, km+1, js-1:je+1)
! u-wind will be ghosted one latitude to the north upon exit
! u-wind (m/s)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: u_ad(isd:ied, jsd:jed+1, km)
! v-wind (m/s)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_ad(isd:ied+1, jsd:jed, km)
! cp*virtual potential temperature 
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: pt_ad(isd:ied, jsd:jed, km)
! as input; output: temperature
    INTEGER, INTENT(IN) :: te_method
    LOGICAL, INTENT(IN) :: remap_t
! u-wind (m/s) on physics grid
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, km)
    REAL, INTENT(INOUT) :: ua_ad(isd:ied, jsd:jed, km)
! log(pe)
    REAL(p_precision), INTENT(INOUT) :: peln(is:ie, km+1, js:je)
    REAL(p_precision) :: peln_ad(is:ie, km+1, js:je)
! layer-mean pk for converting t to pt
    REAL(p_precision) :: pkz(is:ie, js:je, km)
    REAL(p_precision) :: pkz_ad(is:ie, js:je, km)
    REAL :: te(is:ie, js:je, km)
    REAL :: te_ad(is:ie, js:je, km)
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
    REAL(p_precision) :: te_2d_ad(is:ie, js:je)
    REAL(p_precision) :: zsum0(is:ie, js:je)
    REAL(p_precision) :: zsum0_ad(is:ie, js:je)
    REAL(p_precision) :: zsum1(is:ie, js:je)
    REAL(p_precision) :: zsum1_ad(is:ie, js:je)
    REAL(p_precision) :: q2(is:ie, km)
    REAL(p_precision) :: q2_ad(is:ie, km)
    REAL(p_precision) :: dp2(is:ie, km)
    REAL(p_precision) :: dp2_ad(is:ie, km)
    REAL(p_precision) :: pe1(is:ie, km+1)
    REAL(p_precision) :: pe1_ad(is:ie, km+1)
    REAL(p_precision) :: pe2(is:ie, km+1)
    REAL(p_precision) :: pe2_ad(is:ie, km+1)
    REAL(p_precision) :: pk1(is:ie, km+1)
    REAL(p_precision) :: pk1_ad(is:ie, km+1)
    REAL(p_precision) :: pk2(is:ie, km+1)
    REAL(p_precision) :: pk2_ad(is:ie, km+1)
    REAL(p_precision) :: pn1(is:ie, km+1)
    REAL(p_precision) :: pn2(is:ie, km+1)
    REAL(p_precision) :: pn2_ad(is:ie, km+1)
    REAL(p_precision) :: pe0(is:ie+1, km+1)
    REAL(p_precision) :: pe0_ad(is:ie+1, km+1)
    REAL(p_precision) :: pe3(is:ie+1, km+1)
    REAL(p_precision) :: pe3_ad(is:ie+1, km+1)
    REAL(p_precision) :: phis(is:ie, km+1)
    REAL(p_precision) :: phis_ad(is:ie, km+1)
    REAL(p_precision) :: gz(is:ie)
    REAL(p_precision) :: gz_ad(is:ie)
    REAL(p_precision) :: rcp, rg, ak1, tmp, tpe, cv, rgama, rrg
    REAL(p_precision) :: tmp_ad, tpe_ad
    REAL(p_precision) :: bkh
    REAL(p_precision) :: dtmp
    REAL(p_precision) :: dtmp_ad
    REAL(p_precision) :: dlnp
    REAL(p_precision) :: dlnp_ad
    INTEGER :: iq, n, kp, k_next
    LOGICAL :: te_map
    REAL(p_precision) :: k1k, kapag
    REAL(p_precision) :: tmpsum
    REAL(p_precision) :: tmpsum_ad
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC REAL
    INTEGER :: arg1
    INTEGER :: arg2
    INTEGER :: branch
    REAL :: temp3
    REAL(p_precision) :: temp2
    REAL*8 :: temp1
    REAL*8 :: temp0
    REAL*8 :: temp_ad9
    REAL(p_precision) :: temp_ad8
    REAL(p_precision) :: temp_ad7
    REAL(p_precision) :: temp_ad6
    REAL(p_precision) :: temp_ad5
    REAL(p_precision) :: temp_ad4
    REAL(p_precision) :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    REAL :: temp_ad15
    REAL :: temp_ad14
    REAL :: temp_ad13
    REAL :: temp_ad12
    REAL(p_precision) :: temp_ad11
    REAL(p_precision) :: temp_ad10
    REAL(p_precision) :: temp
! rg/Cv=0.4
    rg = akap*cp
! cv/cp
    te_map = .true.
    CALL PUSHREAL8ARRAY(peln, (ie-is+1)*(km+1)*(je-js+1))
    CALL PKEZ(km, is, ie, js, je, pe, pk, akap, peln, pkz)
    DO k=1,km
      DO j=js,je
        DO i=is,ie
          te(i, j, k) = 0.25*rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+&
&           v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i&
&           , j, k)+v(i+1, j, k))*cosa_s(i, j)) + pt(i, j, k)*pkz(i, j, &
&           k)
        END DO
      END DO
    END DO
    DO j=js,je+1
      DO k=1,km+1
        DO i=is,ie
          CALL PUSHREAL8(pe1(i, k))
          pe1(i, k) = pe(i, k, j)
        END DO
      END DO
      DO i=is,ie
        CALL PUSHREAL8(pe2(i, 1))
        pe2(i, 1) = ptop
        CALL PUSHREAL8(pe2(i, km+1))
        pe2(i, km+1) = pe(i, km+1, j)
      END DO
!(j < je+1)
      IF (j .LT. je + 1) THEN
!
! Hybrid sigma-P coordinate:
!
        DO k=2,ks+1
          DO i=is,ie
            CALL PUSHREAL8(pe2(i, k))
            pe2(i, k) = ak(k)
          END DO
        END DO
        DO k=ks+2,km
          DO i=is,ie
            CALL PUSHREAL8(pe2(i, k))
            pe2(i, k) = ak(k) + bk(k)*pe(i, km+1, j)
          END DO
        END DO
        DO k=1,km
          DO i=is,ie
            CALL PUSHREAL8(dp2(i, k))
            dp2(i, k) = pe2(i, k+1) - pe2(i, k)
          END DO
        END DO
!------------
! update delp
!------------
        DO k=1,km
          DO i=is,ie
            delp(i, j, k) = dp2(i, k)
          END DO
        END DO
!-----------------
! Map constituents
!-----------------
        IF (nq .NE. 0) THEN
          DO iq=1,nq
            CALL MAP1_Q2_NOLIMITERS(km, pe1, q(isd, jsd, 1, iq), km, pe2&
&                             , q2, dp2, is, ie, 0, kord_tr, j, isd, ied&
&                             , jsd, jed)
            DO k=1,km
              DO i=is,ie
                CALL PUSHREAL8(q(i, j, k, iq))
                q(i, j, k, iq) = q2(i, k)
              END DO
            END DO
          END DO
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
!------------------
! Compute p**cappa
!------------------
        DO k=1,km+1
          DO i=is,ie
            CALL PUSHREAL8(pk1(i, k))
            pk1(i, k) = pk(i, j, k)
          END DO
        END DO
        DO i=is,ie
          CALL PUSHREAL8(pn2(i, 1))
          pn2(i, 1) = peln(i, 1, j)
          CALL PUSHREAL8(pn2(i, km+1))
          pn2(i, km+1) = peln(i, km+1, j)
          CALL PUSHREAL8(pk2(i, 1))
          pk2(i, 1) = pk1(i, 1)
          CALL PUSHREAL8(pk2(i, km+1))
          pk2(i, km+1) = pk1(i, km+1)
        END DO
        DO k=2,km
          DO i=is,ie
!        pk2(i,k) = pe2(i,k) ** akap
            CALL PUSHREAL8(pn2(i, k))
            pn2(i, k) = LOG(pe2(i, k))
            CALL PUSHREAL8(pk2(i, k))
            pk2(i, k) = EXP(akap*pn2(i, k))
          END DO
        END DO
        IF (te_map) THEN
!---------------------
! Compute Total Energy
!---------------------
          DO i=is,ie
            CALL PUSHREAL8(phis(i, km+1))
            phis(i, km+1) = hs(i, j)
          END DO
          DO k=km,1,-1
            DO i=is,ie
              CALL PUSHREAL8(phis(i, k))
              phis(i, k) = phis(i, k+1) + pt(i, j, k)*(pk1(i, k+1)-pk1(i&
&               , k))
            END DO
          END DO
          DO k=1,km+1
            DO i=is,ie
              CALL PUSHREAL8(phis(i, k))
              phis(i, k) = phis(i, k)*pe1(i, k)
            END DO
          END DO
          DO k=1,km
            DO i=is,ie
              te(i, j, k) = te(i, j, k) + (phis(i, k+1)-phis(i, k))/(pe1&
&               (i, k+1)-pe1(i, k))
            END DO
          END DO
!----------------
! Map Total Energy
!----------------
          CALL PUSHREAL8ARRAY(te, (ie-is+1)*(je-js+1)*km)
          CALL MAP1_CUBIC_TE(km, pe1, pe2, te, is, ie, j, is, ie, js, je&
&                      , 1, kord_tm)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
!----------
! Update pk
!----------
        DO k=2,km
          DO i=is,ie
            CALL PUSHREAL8(pk(i, j, k))
            pk(i, j, k) = pk2(i, k)
          END DO
        END DO
        DO k=1,km+1
          DO i=is,ie
            CALL PUSHREAL8(pe0(i, k))
            pe0(i, k) = peln(i, k, j)
            CALL PUSHREAL8(peln(i, k, j))
            peln(i, k, j) = pn2(i, k)
          END DO
        END DO
!------------
! Compute pkz
!------------
        DO k=1,km
          DO i=is,ie
            CALL PUSHREAL8(pkz(i, j, k))
            pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1, j&
&             )-peln(i, k, j)))
          END DO
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      DO i=is,ie+1
        CALL PUSHREAL8(pe0(i, 1))
        pe0(i, 1) = pe(i, 1, j)
      END DO
!------
! map u
!------
      DO k=2,km+1
        DO i=is,ie
          CALL PUSHREAL8(pe0(i, k))
          pe0(i, k) = 0.5*(pe(i, k, j-1)+pe1(i, k))
        END DO
      END DO
      DO k=1,ks+1
        DO i=is,ie+1
          CALL PUSHREAL8(pe3(i, k))
          pe3(i, k) = ak(k)
        END DO
      END DO
      DO k=ks+2,km+1
        CALL PUSHREAL8(bkh)
        bkh = 0.5*bk(k)
        DO i=is,ie
          CALL PUSHREAL8(pe3(i, k))
          pe3(i, k) = ak(k) + bkh*(pe(i, km+1, j-1)+pe1(i, km+1))
        END DO
      END DO
      arg1 = jed + 1
      CALL PUSHREAL8ARRAY(u, (ied-isd+1)*(jed-jsd+2)*km)
      CALL MAP1_PPM_NOLIMITERS(km, pe0(is:ie, :), pe3(is:ie, :), u, is, &
&                        ie, j, isd, ied, jsd, arg1, -1, kord_mt)
! (j < je+1)
      IF (j .LT. je + 1) THEN
!------
! map v
!------
        DO k=2,km+1
          DO i=is,ie+1
            CALL PUSHREAL8(pe0(i, k))
            pe0(i, k) = 0.5*(pe(i-1, k, j)+pe(i, k, j))
          END DO
        END DO
        DO k=ks+2,km+1
          CALL PUSHREAL8(bkh)
          bkh = 0.5*bk(k)
          DO i=is,ie+1
            CALL PUSHREAL8(pe3(i, k))
            pe3(i, k) = ak(k) + bkh*(pe(i-1, km+1, j)+pe(i, km+1, j))
          END DO
        END DO
        arg1 = ie + 1
        arg2 = ied + 1
        CALL PUSHREAL8ARRAY(v, (ied-isd+2)*(jed-jsd+1)*km)
        CALL MAP1_PPM_NOLIMITERS(km, pe0, pe3, v, is, arg1, j, isd, arg2&
&                          , jsd, jed, -1, kord_mt)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      DO k=1,km
        DO i=is,ie
          ua(i, j, k) = pe2(i, k+1)
        END DO
      END DO
    END DO
    DO k=2,km
      DO j=js,je
        DO i=is,ie
          pe(i, k, j) = ua(i, j, k-1)
        END DO
      END DO
    END DO
    IF (do_consv .AND. consv .GT. 0.) THEN
      IF (te_map) THEN
        DO j=js,je
          DO i=is,ie
            te_2d(i, j) = te(i, j, 1)*delp(i, j, 1)
          END DO
          DO k=2,km
            DO i=is,ie
              te_2d(i, j) = te_2d(i, j) + te(i, j, k)*delp(i, j, k)
            END DO
          END DO
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      DO j=js,je
        DO i=is,ie
          zsum1(i, j) = pkz(i, j, 1)*delp(i, j, 1)
        END DO
        DO k=2,km
          DO i=is,ie
            zsum1(i, j) = zsum1(i, j) + pkz(i, j, k)*delp(i, j, k)
          END DO
        END DO
        DO i=is,ie
          zsum0(i, j) = ptop*(pk(i, j, 1)-pk(i, j, km+1)) + zsum1(i, j)
          te_2d(i, j) = te0_2d(i, j) - te_2d(i, j)
        END DO
      END DO
      tpe = consv*g_sum(te_2d, is, ie, js, je, ng, area, 0)
      tmpsum = g_sum(zsum0,  is, ie, js, je, ng, area, 0)
      dtmp = tpe/(cp*tmpsum)
! convert to 4-byte real
      IF (reproduce_sum) THEN
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
    IF (te_map) THEN
      DO j=js,je
        DO i=is,ie
          gz(i) = hs(i, j)
        END DO
        DO k=km,1,-1
          DO i=is,ie
            CALL PUSHREAL8(tpe)
            tpe = te(i, j, k) - gz(i) - 0.25*rsin2(i, j)*(u(i, j, k)**2+&
&             u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+&
&             u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*cosa_s(i, j))
            dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
            tmp = tpe/(cp-pe(i, k, j)*dlnp/delp(i, j, k))
            gz(i) = gz(i) + dlnp*tmp
          END DO
        END DO
      END DO
      te_ad = 0.0_8
      gz_ad = 0.0_8
      dtmp_ad = 0.0_8
      DO j=je,js,-1
        DO k=1,km,1
          DO i=ie,is,-1
            temp_ad10 = cp*pt_ad(i, j, k)/pkz(i, j, k)
            dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
            tmp = tpe/(cp-pe(i, k, j)*dlnp/delp(i, j, k))
            tmp_ad = temp_ad10 + dlnp*gz_ad(i)
            pkz_ad(i, j, k) = pkz_ad(i, j, k) - tmp*temp_ad10/pkz(i, j, &
&             k)
            dtmp_ad = dtmp_ad + cp*pt_ad(i, j, k)
            pt_ad(i, j, k) = 0.0_8
            temp3 = delp(i, j, k)
            temp2 = pe(i, k, j)
            temp1 = temp2*dlnp/temp3
            temp_ad11 = tmp_ad/(cp-temp1)
            temp_ad9 = tpe*temp_ad11/((cp-temp1)*temp3)
            dlnp_ad = temp2*temp_ad9 + tmp*gz_ad(i)
            tpe_ad = temp_ad11
            pe_ad(i, k, j) = pe_ad(i, k, j) + dlnp*temp_ad9
            delp_ad(i, j, k) = delp_ad(i, j, k) - temp1*temp_ad9
            peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + rg*dlnp_ad
            peln_ad(i, k, j) = peln_ad(i, k, j) - rg*dlnp_ad
            CALL POPREAL8(tpe)
            temp_ad12 = -(rsin2(i, j)*0.25*tpe_ad)
            temp_ad13 = -(cosa_s(i, j)*temp_ad12)
            temp_ad14 = (v(i, j, k)+v(i+1, j, k))*temp_ad13
            temp_ad15 = (u(i, j, k)+u(i, j+1, k))*temp_ad13
            te_ad(i, j, k) = te_ad(i, j, k) + tpe_ad
            gz_ad(i) = gz_ad(i) - tpe_ad
            u_ad(i, j, k) = u_ad(i, j, k) + temp_ad14 + 2*u(i, j, k)*&
&             temp_ad12
            u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp_ad14 + 2*u(i, j+1, &
&             k)*temp_ad12
            v_ad(i, j, k) = v_ad(i, j, k) + temp_ad15 + 2*v(i, j, k)*&
&             temp_ad12
            v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp_ad15 + 2*v(i+1, j, &
&             k)*temp_ad12
          END DO
        END DO
        DO i=ie,is,-1
          gz_ad(i) = 0.0_8
        END DO
      END DO
    ELSE
      te_ad = 0.0_8
      dtmp_ad = 0.0_8
    END IF
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      tpe = consv*g_sum(te_2d, is, ie, js, je, ng, area, 0)
    ELSE IF (branch .NE. 1) THEN
      te0_2d_ad = 0.0_8
      pk_ad = 0.0_8
      GOTO 100
    END IF

    !Forward TLM code
    !tpe = consv*g_sum(te_2d, is, ie, js, je, ng, area, 0)
    !tmpsum = g_sum(zsum0,  is, ie, js, je, ng, area, 0)
    !tpe_tl = consv*g_sum(te_2d_tl, is, ie, js, je, ng, area, 0)
    !tmpsum_tl = g_sum(zsum0_tl,  is, ie, js, je, ng, area, 0)
    !dtmp_tl = (tpe_tl*cp*tmpsum-tpe*cp*tmpsum_tl)/(cp*tmpsum)**2
    !dtmp = tpe/(cp*tmpsum)

    !Tapenade linearized from arbitrary self adjoint function
    !tpe = 2*te_2d(is,ie)
    !tmpsum = 2*zsum0(is,ie)

    temp_ad8 = dtmp_ad/(cp*tmpsum)
    tpe_ad = temp_ad8
    tmpsum_ad = -(tpe*temp_ad8/tmpsum)

    zsum0_ad = 0.0_8
    zsum0_ad = zsum0_ad + tmpsum_ad
    te_2d_ad = 0.0_8
    te_2d_ad = te_2d_ad + tpe_ad

    te0_2d_ad = 0.0_8
    pk_ad = 0.0_8
    zsum1_ad = 0.0_8
    DO j=je,js,-1
      DO i=ie,is,-1
        te0_2d_ad(i, j) = te0_2d_ad(i, j) + te_2d_ad(i, j)
        te_2d_ad(i, j) = -te_2d_ad(i, j)
        pk_ad(i, j, 1) = pk_ad(i, j, 1) + ptop*zsum0_ad(i, j)
        pk_ad(i, j, km+1) = pk_ad(i, j, km+1) - ptop*zsum0_ad(i, j)
        zsum1_ad(i, j) = zsum1_ad(i, j) + zsum0_ad(i, j)
        zsum0_ad(i, j) = 0.0_8
      END DO
      DO k=km,2,-1
        DO i=ie,is,-1
          pkz_ad(i, j, k) = pkz_ad(i, j, k) + delp(i, j, k)*zsum1_ad(i, &
&           j)
          delp_ad(i, j, k) = delp_ad(i, j, k) + pkz(i, j, k)*zsum1_ad(i&
&           , j)
        END DO
      END DO
      DO i=ie,is,-1
        pkz_ad(i, j, 1) = pkz_ad(i, j, 1) + delp(i, j, 1)*zsum1_ad(i, j)
        delp_ad(i, j, 1) = delp_ad(i, j, 1) + pkz(i, j, 1)*zsum1_ad(i, j&
&         )
        zsum1_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      DO j=je,js,-1
        DO k=km,2,-1
          DO i=ie,is,-1
            te_ad(i, j, k) = te_ad(i, j, k) + delp(i, j, k)*te_2d_ad(i, &
&             j)
            delp_ad(i, j, k) = delp_ad(i, j, k) + te(i, j, k)*te_2d_ad(i&
&             , j)
          END DO
        END DO
        DO i=ie,is,-1
          te_ad(i, j, 1) = te_ad(i, j, 1) + delp(i, j, 1)*te_2d_ad(i, j)
          delp_ad(i, j, 1) = delp_ad(i, j, 1) + te(i, j, 1)*te_2d_ad(i, &
&           j)
          te_2d_ad(i, j) = 0.0_8
        END DO
      END DO
    END IF
 100 ua_ad = 0.0_8
    DO k=km,2,-1
      DO j=je,js,-1
        DO i=ie,is,-1
          ua_ad(i, j, k-1) = ua_ad(i, j, k-1) + pe_ad(i, k, j)
          pe_ad(i, k, j) = 0.0_8
        END DO
      END DO
    END DO
    phis_ad = 0.0_8
    pe0_ad = 0.0_8
    pe1_ad = 0.0_8
    pe2_ad = 0.0_8
    pe3_ad = 0.0_8
    dp2_ad = 0.0_8
    q2_ad = 0.0_8
    pn2_ad = 0.0_8
    pk1_ad = 0.0_8
    pk2_ad = 0.0_8
    DO j=je+1,js,-1
      DO k=km,1,-1
        DO i=ie,is,-1
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + ua_ad(i, j, k)
          ua_ad(i, j, k) = 0.0_8
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        arg1 = ie + 1
        CALL POPREAL8ARRAY(v, (ied-isd+2)*(jed-jsd+1)*km)
        CALL MAP1_PPM_NOLIMITERS_ADM(km, pe0, pe0_ad, pe3, pe3_ad, v, &
&                              v_ad, is, arg1, j, isd, arg2, jsd, jed, -&
&                              1, kord_mt)
        DO k=km+1,ks+2,-1
          DO i=ie+1,is,-1
            CALL POPREAL8(pe3(i, k))
            pe_ad(i-1, km+1, j) = pe_ad(i-1, km+1, j) + bkh*pe3_ad(i, k)
            pe_ad(i, km+1, j) = pe_ad(i, km+1, j) + bkh*pe3_ad(i, k)
            pe3_ad(i, k) = 0.0_8
          END DO
          CALL POPREAL8(bkh)
        END DO
        DO k=km+1,2,-1
          DO i=ie+1,is,-1
            CALL POPREAL8(pe0(i, k))
            pe_ad(i-1, k, j) = pe_ad(i-1, k, j) + 0.5*pe0_ad(i, k)
            pe_ad(i, k, j) = pe_ad(i, k, j) + 0.5*pe0_ad(i, k)
            pe0_ad(i, k) = 0.0_8
          END DO
        END DO
      END IF
      CALL POPREAL8ARRAY(u, (ied-isd+1)*(jed-jsd+2)*km)
      CALL MAP1_PPM_NOLIMITERS_ADM(km, pe0(is:ie, :), pe0_ad(is:ie, :), &
&                            pe3(is:ie, :), pe3_ad(is:ie, :), u, u_ad, &
&                            is, ie, j, isd, ied, jsd, arg1, -1, kord_mt&
&                           )
      DO k=km+1,ks+2,-1
        DO i=ie,is,-1
          CALL POPREAL8(pe3(i, k))
          pe_ad(i, km+1, j-1) = pe_ad(i, km+1, j-1) + bkh*pe3_ad(i, k)
          pe1_ad(i, km+1) = pe1_ad(i, km+1) + bkh*pe3_ad(i, k)
          pe3_ad(i, k) = 0.0_8
        END DO
        CALL POPREAL8(bkh)
      END DO
      DO k=ks+1,1,-1
        DO i=ie+1,is,-1
          CALL POPREAL8(pe3(i, k))
          pe3_ad(i, k) = 0.0_8
        END DO
      END DO
      DO k=km+1,2,-1
        DO i=ie,is,-1
          CALL POPREAL8(pe0(i, k))
          pe_ad(i, k, j-1) = pe_ad(i, k, j-1) + 0.5*pe0_ad(i, k)
          pe1_ad(i, k) = pe1_ad(i, k) + 0.5*pe0_ad(i, k)
          pe0_ad(i, k) = 0.0_8
        END DO
      END DO
      DO i=ie+1,is,-1
        CALL POPREAL8(pe0(i, 1))
        pe_ad(i, 1, j) = pe_ad(i, 1, j) + pe0_ad(i, 1)
        pe0_ad(i, 1) = 0.0_8
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        DO k=km,1,-1
          DO i=ie,is,-1
            CALL POPREAL8(pkz(i, j, k))
            temp0 = akap*(peln(i, k+1, j)-peln(i, k, j))
            temp_ad6 = pkz_ad(i, j, k)/temp0
            temp_ad7 = -((pk2(i, k+1)-pk2(i, k))*akap*temp_ad6/temp0)
            pk2_ad(i, k+1) = pk2_ad(i, k+1) + temp_ad6
            pk2_ad(i, k) = pk2_ad(i, k) - temp_ad6
            peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad7
            peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad7
            pkz_ad(i, j, k) = 0.0_8
          END DO
        END DO
        DO k=km+1,1,-1
          DO i=ie,is,-1
            CALL POPREAL8(peln(i, k, j))
            pn2_ad(i, k) = pn2_ad(i, k) + peln_ad(i, k, j)
            peln_ad(i, k, j) = pe0_ad(i, k)
            CALL POPREAL8(pe0(i, k))
            pe0_ad(i, k) = 0.0_8
          END DO
        END DO
        DO k=km,2,-1
          DO i=ie,is,-1
            CALL POPREAL8(pk(i, j, k))
            pk2_ad(i, k) = pk2_ad(i, k) + pk_ad(i, j, k)
            pk_ad(i, j, k) = 0.0_8
          END DO
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          CALL POPREAL8ARRAY(te, (ie-is+1)*(je-js+1)*km)
          CALL MAP1_CUBIC_TE_ADM(km, pe1, pe1_ad, pe2, pe2_ad, te, te_ad&
&                          , is, ie, j, is, ie, js, je, 1, kord_tm)
          DO k=km,1,-1
            DO i=ie,is,-1
              temp = pe1(i, k+1) - pe1(i, k)
              temp_ad4 = te_ad(i, j, k)/temp
              temp_ad5 = -((phis(i, k+1)-phis(i, k))*temp_ad4/temp)
              phis_ad(i, k+1) = phis_ad(i, k+1) + temp_ad4
              phis_ad(i, k) = phis_ad(i, k) - temp_ad4
              pe1_ad(i, k+1) = pe1_ad(i, k+1) + temp_ad5
              pe1_ad(i, k) = pe1_ad(i, k) - temp_ad5
            END DO
          END DO
          DO k=km+1,1,-1
            DO i=ie,is,-1
              CALL POPREAL8(phis(i, k))
              pe1_ad(i, k) = pe1_ad(i, k) + phis(i, k)*phis_ad(i, k)
              phis_ad(i, k) = pe1(i, k)*phis_ad(i, k)
            END DO
          END DO
          DO k=1,km,1
            DO i=ie,is,-1
              CALL POPREAL8(phis(i, k))
              temp_ad3 = pt(i, j, k)*phis_ad(i, k)
              phis_ad(i, k+1) = phis_ad(i, k+1) + phis_ad(i, k)
              pt_ad(i, j, k) = pt_ad(i, j, k) + (pk1(i, k+1)-pk1(i, k))*&
&               phis_ad(i, k)
              pk1_ad(i, k+1) = pk1_ad(i, k+1) + temp_ad3
              pk1_ad(i, k) = pk1_ad(i, k) - temp_ad3
              phis_ad(i, k) = 0.0_8
            END DO
          END DO
          DO i=ie,is,-1
            CALL POPREAL8(phis(i, km+1))
            phis_ad(i, km+1) = 0.0_8
          END DO
        END IF
        DO k=km,2,-1
          DO i=ie,is,-1
            CALL POPREAL8(pk2(i, k))
            pn2_ad(i, k) = pn2_ad(i, k) + EXP(akap*pn2(i, k))*akap*&
&             pk2_ad(i, k)
            pk2_ad(i, k) = 0.0_8
            CALL POPREAL8(pn2(i, k))
            pe2_ad(i, k) = pe2_ad(i, k) + pn2_ad(i, k)/pe2(i, k)
            pn2_ad(i, k) = 0.0_8
          END DO
        END DO
        DO i=ie,is,-1
          CALL POPREAL8(pk2(i, km+1))
          pk1_ad(i, km+1) = pk1_ad(i, km+1) + pk2_ad(i, km+1)
          pk2_ad(i, km+1) = 0.0_8
          CALL POPREAL8(pk2(i, 1))
          pk1_ad(i, 1) = pk1_ad(i, 1) + pk2_ad(i, 1)
          pk2_ad(i, 1) = 0.0_8
          CALL POPREAL8(pn2(i, km+1))
          peln_ad(i, km+1, j) = peln_ad(i, km+1, j) + pn2_ad(i, km+1)
          pn2_ad(i, km+1) = 0.0_8
          CALL POPREAL8(pn2(i, 1))
          peln_ad(i, 1, j) = peln_ad(i, 1, j) + pn2_ad(i, 1)
          pn2_ad(i, 1) = 0.0_8
        END DO
        DO k=km+1,1,-1
          DO i=ie,is,-1
            CALL POPREAL8(pk1(i, k))
            pk_ad(i, j, k) = pk_ad(i, j, k) + pk1_ad(i, k)
            pk1_ad(i, k) = 0.0_8
          END DO
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          DO iq=nq,1,-1
            DO k=km,1,-1
              DO i=ie,is,-1
                CALL POPREAL8(q(i, j, k, iq))
                q2_ad(i, k) = q2_ad(i, k) + q_ad(i, j, k, iq)
                q_ad(i, j, k, iq) = 0.0_8
              END DO
            END DO
            CALL MAP1_Q2_NOLIMITERS_ADM(km, pe1, pe1_ad, q(isd, jsd, 1, &
&                                 iq), q_ad(isd, jsd, 1, iq), km, pe2, &
&                                 pe2_ad, q2, q2_ad, dp2, dp2_ad, is, ie&
&                                 , 0, kord_tr, j, isd, ied, jsd, jed)
          END DO
        END IF
        DO k=km,1,-1
          DO i=ie,is,-1
            dp2_ad(i, k) = dp2_ad(i, k) + delp_ad(i, j, k)
            delp_ad(i, j, k) = 0.0_8
          END DO
        END DO
        DO k=km,1,-1
          DO i=ie,is,-1
            CALL POPREAL8(dp2(i, k))
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp2_ad(i, k)
            pe2_ad(i, k) = pe2_ad(i, k) - dp2_ad(i, k)
            dp2_ad(i, k) = 0.0_8
          END DO
        END DO
        DO k=km,ks+2,-1
          DO i=ie,is,-1
            CALL POPREAL8(pe2(i, k))
            pe_ad(i, km+1, j) = pe_ad(i, km+1, j) + bk(k)*pe2_ad(i, k)
            pe2_ad(i, k) = 0.0_8
          END DO
        END DO
        DO k=ks+1,2,-1
          DO i=ie,is,-1
            CALL POPREAL8(pe2(i, k))
            pe2_ad(i, k) = 0.0_8
          END DO
        END DO
      END IF
      DO i=ie,is,-1
        CALL POPREAL8(pe2(i, km+1))
        pe_ad(i, km+1, j) = pe_ad(i, km+1, j) + pe2_ad(i, km+1)
        pe2_ad(i, km+1) = 0.0_8
        CALL POPREAL8(pe2(i, 1))
        pe2_ad(i, 1) = 0.0_8
      END DO
      DO k=km+1,1,-1
        DO i=ie,is,-1
          CALL POPREAL8(pe1(i, k))
          pe_ad(i, k, j) = pe_ad(i, k, j) + pe1_ad(i, k)
          pe1_ad(i, k) = 0.0_8
        END DO
      END DO
    END DO
    DO k=km,1,-1
      DO j=je,js,-1
        DO i=ie,is,-1
          temp_ad = rsin2(i, j)*0.25*te_ad(i, j, k)
          temp_ad0 = -(cosa_s(i, j)*temp_ad)
          temp_ad1 = (v(i, j, k)+v(i+1, j, k))*temp_ad0
          temp_ad2 = (u(i, j, k)+u(i, j+1, k))*temp_ad0
          u_ad(i, j, k) = u_ad(i, j, k) + temp_ad1 + 2*u(i, j, k)*&
&           temp_ad
          u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp_ad1 + 2*u(i, j+1, k)*&
&           temp_ad
          v_ad(i, j, k) = v_ad(i, j, k) + temp_ad2 + 2*v(i, j, k)*&
&           temp_ad
          v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp_ad2 + 2*v(i+1, j, k)*&
&           temp_ad
          pt_ad(i, j, k) = pt_ad(i, j, k) + pkz(i, j, k)*te_ad(i, j, k)
          pkz_ad(i, j, k) = pkz_ad(i, j, k) + pt(i, j, k)*te_ad(i, j, k)
          te_ad(i, j, k) = 0.0_8
        END DO
      END DO
    END DO
    CALL POPREAL8ARRAY(peln, (ie-is+1)*(km+1)*(je-js+1))
    CALL PKEZ_ADM(km, is, ie, js, je, pe, pk, pk_ad, akap, peln, peln_ad&
&           , pkz, pkz_ad)
  END SUBROUTINE LAGRANGIAN_TO_EULERIAN_ADM
!#ifdef MAPL_MODE
!#endif
  SUBROUTINE LAGRANGIAN_TO_EULERIAN(do_consv, consv, pe, delp, pkz, pk, &
&   pdt, km, is, ie, js, je, isd, ied, jsd, jed, nq, sphum, u, v, pt, q&
&   , hs, r_vir, cp, akap, pi, radius, grav, kord_mt, kord_tr, kord_tm, &
&   peln, te0_2d, ng, ua, te, pem, fill, reproduce_sum, ak, bk, ks, &
&   te_method, remap_t, ncnst)
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
! fill negative tracers
    LOGICAL, INTENT(IN) :: fill
    LOGICAL, INTENT(IN) :: reproduce_sum
    REAL, INTENT(IN) :: ak(km+1)
    REAL, INTENT(IN) :: bk(km+1)
! !INPUT/OUTPUT
! pe to the kappa
    REAL(p_precision), INTENT(INOUT) :: pk(is:ie, js:je, km+1)
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, km, *)
! pressure thickness
    REAL, INTENT(INOUT) :: delp(isd:ied, jsd:jed, km)
! pressure at layer edges
    REAL(p_precision), INTENT(INOUT) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL, INTENT(INOUT) :: pem(is-1:ie+1, km+1, js-1:je+1)
! u-wind will be ghosted one latitude to the north upon exit
! u-wind (m/s)
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
! v-wind (m/s)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
! cp*virtual potential temperature 
    REAL, INTENT(INOUT) :: pt(isd:ied, jsd:jed, km)
! as input; output: temperature
    INTEGER, INTENT(IN) :: te_method
    LOGICAL, INTENT(IN) :: remap_t
! u-wind (m/s) on physics grid
    REAL, INTENT(INOUT) :: ua(isd:ied, jsd:jed, km)
! log(pe)
    REAL(p_precision), INTENT(INOUT) :: peln(is:ie, km+1, js:je)
! layer-mean pk for converting t to pt
    REAL(p_precision), INTENT(OUT) :: pkz(is:ie, js:je, km)
    REAL, INTENT(OUT) :: te(is:ie, js:je, km)
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
    REAL(p_precision) :: zsum0(is:ie, js:je)
    REAL(p_precision) :: zsum1(is:ie, js:je)
    REAL(p_precision) :: q2(is:ie, km)
    REAL(p_precision) :: dp2(is:ie, km)
    REAL(p_precision) :: pe1(is:ie, km+1)
    REAL(p_precision) :: pe2(is:ie, km+1)
    REAL(p_precision) :: pk1(is:ie, km+1)
    REAL(p_precision) :: pk2(is:ie, km+1)
    REAL(p_precision) :: pn1(is:ie, km+1)
    REAL(p_precision) :: pn2(is:ie, km+1)
    REAL(p_precision) :: pe0(is:ie+1, km+1)
    REAL(p_precision) :: pe3(is:ie+1, km+1)
    REAL(p_precision) :: phis(is:ie, km+1)
    REAL(p_precision) :: gz(is:ie)
    REAL(p_precision) :: rcp, rg, ak1, tmp, tpe, cv, rgama, rrg
    REAL(p_precision) :: bkh
    REAL(p_precision) :: dtmp
    REAL(p_precision) :: dlnp
    INTEGER :: iq, n, kp, k_next
    LOGICAL :: te_map
    REAL(p_precision) :: k1k, kapag
    REAL(p_precision) :: tmpsum
    INTRINSIC LOG
    INTRINSIC EXP
    INTRINSIC REAL
    INTEGER :: arg1
    INTEGER :: arg2
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
    CALL PKEZ(km, is, ie, js, je, pe, pk, akap, peln, pkz)
    DO k=1,km
      DO j=js,je
        DO i=is,ie
          te(i, j, k) = 0.25*rsin2(i, j)*(u(i, j, k)**2+u(i, j+1, k)**2+&
&           v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+u(i, j+1, k))*(v(i&
&           , j, k)+v(i+1, j, k))*cosa_s(i, j)) + pt(i, j, k)*pkz(i, j, &
&           k)
        END DO
      END DO
    END DO
    DO j=js,je+1
      DO k=1,km+1
        DO i=is,ie
          pe1(i, k) = pe(i, k, j)
        END DO
      END DO
      DO i=is,ie
        pe2(i, 1) = ptop
        pe2(i, km+1) = pe(i, km+1, j)
      END DO
!(j < je+1)
      IF (j .LT. je + 1) THEN
!
! Hybrid sigma-P coordinate:
!
        DO k=2,ks+1
          DO i=is,ie
            pe2(i, k) = ak(k)
          END DO
        END DO
        DO k=ks+2,km
          DO i=is,ie
            pe2(i, k) = ak(k) + bk(k)*pe(i, km+1, j)
          END DO
        END DO
        DO k=1,km
          DO i=is,ie
            dp2(i, k) = pe2(i, k+1) - pe2(i, k)
          END DO
        END DO
!------------
! update delp
!------------
        DO k=1,km
          DO i=is,ie
            delp(i, j, k) = dp2(i, k)
          END DO
        END DO
!-----------------
! Map constituents
!-----------------
        IF (nq .NE. 0) THEN
          DO iq=1,nq
            CALL MAP1_Q2_NOLIMITERS(km, pe1, q(isd, jsd, 1, iq), km, pe2&
&                             , q2, dp2, is, ie, 0, kord_tr, j, isd, ied&
&                             , jsd, jed)
            DO k=1,km
              DO i=is,ie
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
            pk1(i, k) = pk(i, j, k)
          END DO
        END DO
        DO i=is,ie
          pn1(i, :) = peln(i, :, j)
          pn2(i, 1) = peln(i, 1, j)
          pn2(i, km+1) = peln(i, km+1, j)
          pk2(i, 1) = pk1(i, 1)
          pk2(i, km+1) = pk1(i, km+1)
        END DO
        DO k=2,km
          DO i=is,ie
!        pk2(i,k) = pe2(i,k) ** akap
            pn2(i, k) = LOG(pe2(i, k))
            pk2(i, k) = EXP(akap*pn2(i, k))
          END DO
        END DO
        IF (te_map) THEN
!---------------------
! Compute Total Energy
!---------------------
          DO i=is,ie
            phis(i, km+1) = hs(i, j)
          END DO
          DO k=km,1,-1
            DO i=is,ie
              phis(i, k) = phis(i, k+1) + pt(i, j, k)*(pk1(i, k+1)-pk1(i&
&               , k))
            END DO
          END DO
          DO k=1,km+1
            DO i=is,ie
              phis(i, k) = phis(i, k)*pe1(i, k)
            END DO
          END DO
          DO k=1,km
            DO i=is,ie
              te(i, j, k) = te(i, j, k) + (phis(i, k+1)-phis(i, k))/(pe1&
&               (i, k+1)-pe1(i, k))
            END DO
          END DO
!----------------
! Map Total Energy
!----------------
          CALL MAP1_CUBIC_TE(km, pe1, pe2, te, is, ie, j, is, ie, js, je&
&                      , 1, kord_tm)
        END IF
!----------
! Update pk
!----------
        DO k=2,km
          DO i=is,ie
            pk(i, j, k) = pk2(i, k)
          END DO
        END DO
        DO k=1,km+1
          DO i=is,ie
            pe0(i, k) = peln(i, k, j)
            peln(i, k, j) = pn2(i, k)
          END DO
        END DO
!------------
! Compute pkz
!------------
        DO k=1,km
          DO i=is,ie
            pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1, j&
&             )-peln(i, k, j)))
          END DO
        END DO
      END IF
      DO i=is,ie+1
        pe0(i, 1) = pe(i, 1, j)
      END DO
!------
! map u
!------
      DO k=2,km+1
        DO i=is,ie
          pe0(i, k) = 0.5*(pe(i, k, j-1)+pe1(i, k))
        END DO
      END DO
      DO k=1,ks+1
        DO i=is,ie+1
          pe3(i, k) = ak(k)
        END DO
      END DO
      DO k=ks+2,km+1
        bkh = 0.5*bk(k)
        DO i=is,ie
          pe3(i, k) = ak(k) + bkh*(pe(i, km+1, j-1)+pe1(i, km+1))
        END DO
      END DO
      arg1 = jed + 1
      CALL MAP1_PPM_NOLIMITERS(km, pe0(is:ie, :), pe3(is:ie, :), u, is, &
&                        ie, j, isd, ied, jsd, arg1, -1, kord_mt)
! (j < je+1)
      IF (j .LT. je + 1) THEN
!------
! map v
!------
        DO k=2,km+1
          DO i=is,ie+1
            pe0(i, k) = 0.5*(pe(i-1, k, j)+pe(i, k, j))
          END DO
        END DO
        DO k=ks+2,km+1
          bkh = 0.5*bk(k)
          DO i=is,ie+1
            pe3(i, k) = ak(k) + bkh*(pe(i-1, km+1, j)+pe(i, km+1, j))
          END DO
        END DO
        arg1 = ie + 1
        arg2 = ied + 1
        CALL MAP1_PPM_NOLIMITERS(km, pe0, pe3, v, is, arg1, j, isd, arg2&
&                          , jsd, jed, -1, kord_mt)
      END IF
      DO k=1,km
        DO i=is,ie
          ua(i, j, k) = pe2(i, k+1)
        END DO
      END DO
    END DO
    DO k=2,km
      DO j=js,je
        DO i=is,ie
          pe(i, k, j) = ua(i, j, k-1)
        END DO
      END DO
    END DO
    IF (do_consv .AND. consv .GT. 0.) THEN
      IF (te_map) THEN
        DO j=js,je
          DO i=is,ie
            te_2d(i, j) = te(i, j, 1)*delp(i, j, 1)
          END DO
          DO k=2,km
            DO i=is,ie
              te_2d(i, j) = te_2d(i, j) + te(i, j, k)*delp(i, j, k)
            END DO
          END DO
        END DO
      END IF
      DO j=js,je
        DO i=is,ie
          zsum1(i, j) = pkz(i, j, 1)*delp(i, j, 1)
        END DO
        DO k=2,km
          DO i=is,ie
            zsum1(i, j) = zsum1(i, j) + pkz(i, j, k)*delp(i, j, k)
          END DO
        END DO
        DO i=is,ie
          zsum0(i, j) = ptop*(pk(i, j, 1)-pk(i, j, km+1)) + zsum1(i, j)
          te_2d(i, j) = te0_2d(i, j) - te_2d(i, j)
        END DO
      END DO
      tpe = consv*g_sum(te_2d, is, ie, js, je, ng, area, 0)
      tmpsum = g_sum(zsum0,  is, ie, js, je, ng, area, 0)
      dtmp = tpe/(cp*tmpsum)
! convert to 4-byte real
      IF (reproduce_sum) dtmp = REAL(dtmp, 4)
    ELSE
      dtmp = 0.
    END IF
    IF (te_map) THEN
      DO j=js,je
        DO i=is,ie
          gz(i) = hs(i, j)
        END DO
        DO k=km,1,-1
          DO i=is,ie
            tpe = te(i, j, k) - gz(i) - 0.25*rsin2(i, j)*(u(i, j, k)**2+&
&             u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-(u(i, j, k)+&
&             u(i, j+1, k))*(v(i, j, k)+v(i+1, j, k))*cosa_s(i, j))
            dlnp = rg*(peln(i, k+1, j)-peln(i, k, j))
            tmp = tpe/(cp-pe(i, k, j)*dlnp/delp(i, j, k))
            pt(i, j, k) = cp*(tmp/pkz(i, j, k)+dtmp)
            gz(i) = gz(i) + dlnp*tmp
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE LAGRANGIAN_TO_EULERIAN
!  Differentiation of compute_total_energy in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: peln q u v delp delz te_2d
!                pe pt
!   with respect to varying inputs: peln q u v delp delz pe pt
!#ifdef MAPL_MODE
!#endif
  SUBROUTINE COMPUTE_TOTAL_ENERGY_ADM(is, ie, js, je, isd, ied, jsd, jed&
&   , km, u, u_ad, v, v_ad, w, delz, delz_ad, pt, pt_ad, delp, delp_ad, &
&   q, q_ad, pe, pe_ad, peln, peln_ad, hs, grav, r_vir, cp, rg, hlv, &
&   te_2d, te_2d_ad, ua, teq, moist_phys, sphum, hydrostatic, id_te)
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
    REAL, DIMENSION(isd:ied, jsd:jed, km) :: pt_ad, delp_ad
    REAL, DIMENSION(isd:ied, jsd:jed, km, sphum), INTENT(IN) :: q
    REAL, DIMENSION(isd:ied, jsd:jed, km, sphum) :: q_ad
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: u_ad(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
    REAL, INTENT(INOUT) :: v_ad(isd:ied+1, jsd:jed, km)
! vertical velocity (m/s)
    REAL, INTENT(IN) :: w(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: delz(is:ie, js:je, km)
    REAL :: delz_ad(is:ie, js:je, km)
! surface geopotential
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
! pressure at layer edges
    REAL(p_precision), INTENT(IN) :: pe(is-1:ie+1, km+1, js-1:je+1)
    REAL(p_precision) :: pe_ad(is-1:ie+1, km+1, js-1:je+1)
! log(pe)
    REAL(p_precision), INTENT(IN) :: peln(is:ie, km+1, js:je)
    REAL(p_precision) :: peln_ad(is:ie, km+1, js:je)
    REAL, INTENT(IN) :: cp, rg, r_vir, hlv
!#ifdef MAPL_MODE                 
    REAL, INTENT(IN) :: grav
!#endif
    LOGICAL, INTENT(IN) :: moist_phys, hydrostatic
! Output:
! vertically integrated TE
    REAL :: te_2d(is:ie, js:je)
    REAL :: te_2d_ad(is:ie, js:je)
! Moist TE
    REAL, INTENT(OUT) :: teq(is:ie, js:je)
! Local
    REAL(p_precision), DIMENSION(is:ie, km) :: tv
    REAL(p_precision), DIMENSION(is:ie, km) :: tv_ad
    REAL(p_precision) :: phiz(is:ie, km+1)
    REAL(p_precision) :: phiz_ad(is:ie, km+1)
    REAL(p_precision) :: cv
    INTEGER :: i, j, k
    INTEGER :: branch
    REAL(p_precision) :: temp3
    REAL(p_precision) :: temp2
    REAL :: temp1
    REAL :: temp0
    REAL :: temp_ad6
    REAL :: temp_ad5
    REAL(p_precision) :: temp_ad4
    REAL(p_precision) :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL*8 :: temp_ad0
    REAL(p_precision) :: temp_ad
    REAL*8 :: temp
    REAL :: temp5
    REAL :: temp4
    cv = cp - rg
!----------------------
! Output lat-lon winds:
!----------------------
!  call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km)
!$omp parallel do default(shared) private(phiz, tv)
    DO j=js,je
      IF (hydrostatic) THEN
        DO i=is,ie
          CALL PUSHREAL8(phiz(i, km+1))
          phiz(i, km+1) = hs(i, j)
        END DO
        DO k=km,1,-1
          DO i=is,ie
            CALL PUSHREAL8(tv(i, k))
            tv(i, k) = pt(i, j, k)*(1.+r_vir*q(i, j, k, sphum))
            CALL PUSHREAL8(phiz(i, k))
            phiz(i, k) = phiz(i, k+1) + rg*tv(i, k)*(peln(i, k+1, j)-&
&             peln(i, k, j))
          END DO
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
!-----------------
! Non-hydrostatic:
!-----------------
        DO i=is,ie
          CALL PUSHREAL8(phiz(i, km+1))
          phiz(i, km+1) = hs(i, j)
          DO k=km,1,-1
            CALL PUSHREAL8(phiz(i, k))
            phiz(i, k) = phiz(i, k+1) - grav*delz(i, j, k)
          END DO
        END DO
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    phiz_ad = 0.0_8
    tv_ad = 0.0_8
    DO j=je,js,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO k=km,1,-1
          DO i=ie,is,-1
            temp5 = v(i, j, k) + v(i+1, j, k)
            temp4 = u(i, j, k) + u(i, j+1, k)
            temp3 = 0.5*rsin2(i, j)
            temp2 = r_vir*q(i, j, k, sphum) + 1.
            temp_ad3 = delp(i, j, k)*te_2d_ad(i, j)
            temp_ad4 = 0.5*temp_ad3
            temp_ad5 = temp3*temp_ad4
            temp_ad6 = -(cosa_s(i, j)*temp_ad5)
            delp_ad(i, j, k) = delp_ad(i, j, k) + (cv*(pt(i, j, k)*temp2&
&             )+0.5*(phiz(i, k)+phiz(i, k+1)+temp3*(u(i, j, k)**2+u(i, j&
&             +1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-cosa_s(i, j)*(&
&             temp4*temp5))))*te_2d_ad(i, j)
            pt_ad(i, j, k) = pt_ad(i, j, k) + cv*temp2*temp_ad3
            q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) + pt(i, j, k)*cv&
&             *r_vir*temp_ad3
            phiz_ad(i, k) = phiz_ad(i, k) + temp_ad4
            phiz_ad(i, k+1) = phiz_ad(i, k+1) + temp_ad4
            u_ad(i, j, k) = u_ad(i, j, k) + temp5*temp_ad6 + 2*u(i, j, k&
&             )*temp_ad5
            u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp5*temp_ad6 + 2*u(i, &
&             j+1, k)*temp_ad5
            v_ad(i, j, k) = v_ad(i, j, k) + temp4*temp_ad6 + 2*v(i, j, k&
&             )*temp_ad5
            v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp4*temp_ad6 + 2*v(i+1&
&             , j, k)*temp_ad5
          END DO
        END DO
        DO i=ie,is,-1
          te_2d_ad(i, j) = 0.0_8
        END DO
        DO i=ie,is,-1
          DO k=1,km,1
            CALL POPREAL8(phiz(i, k))
            phiz_ad(i, k+1) = phiz_ad(i, k+1) + phiz_ad(i, k)
            delz_ad(i, j, k) = delz_ad(i, j, k) - grav*phiz_ad(i, k)
            phiz_ad(i, k) = 0.0_8
          END DO
          CALL POPREAL8(phiz(i, km+1))
          phiz_ad(i, km+1) = 0.0_8
        END DO
      ELSE
        DO k=km,1,-1
          DO i=ie,is,-1
            temp1 = v(i, j, k) + v(i+1, j, k)
            temp0 = u(i, j, k) + u(i, j+1, k)
            temp = 0.25*rsin2(i, j)
            temp_ad0 = delp(i, j, k)*te_2d_ad(i, j)
            temp_ad1 = temp*temp_ad0
            temp_ad2 = -(cosa_s(i, j)*temp_ad1)
            delp_ad(i, j, k) = delp_ad(i, j, k) + (cp*tv(i, k)+temp*(u(i&
&             , j, k)**2+u(i, j+1, k)**2+v(i, j, k)**2+v(i+1, j, k)**2-&
&             cosa_s(i, j)*(temp0*temp1)))*te_2d_ad(i, j)
            tv_ad(i, k) = tv_ad(i, k) + cp*temp_ad0
            u_ad(i, j, k) = u_ad(i, j, k) + temp1*temp_ad2 + 2*u(i, j, k&
&             )*temp_ad1
            u_ad(i, j+1, k) = u_ad(i, j+1, k) + temp1*temp_ad2 + 2*u(i, &
&             j+1, k)*temp_ad1
            v_ad(i, j, k) = v_ad(i, j, k) + temp0*temp_ad2 + 2*v(i, j, k&
&             )*temp_ad1
            v_ad(i+1, j, k) = v_ad(i+1, j, k) + temp0*temp_ad2 + 2*v(i+1&
&             , j, k)*temp_ad1
          END DO
        END DO
        DO i=ie,is,-1
          pe_ad(i, km+1, j) = pe_ad(i, km+1, j) + phiz(i, km+1)*te_2d_ad&
&           (i, j)
          phiz_ad(i, km+1) = phiz_ad(i, km+1) + pe(i, km+1, j)*te_2d_ad(&
&           i, j)
          pe_ad(i, 1, j) = pe_ad(i, 1, j) - phiz(i, 1)*te_2d_ad(i, j)
          phiz_ad(i, 1) = phiz_ad(i, 1) - pe(i, 1, j)*te_2d_ad(i, j)
          te_2d_ad(i, j) = 0.0_8
        END DO
        DO k=1,km,1
          DO i=ie,is,-1
            CALL POPREAL8(phiz(i, k))
            temp_ad = rg*tv(i, k)*phiz_ad(i, k)
            phiz_ad(i, k+1) = phiz_ad(i, k+1) + phiz_ad(i, k)
            tv_ad(i, k) = tv_ad(i, k) + rg*(peln(i, k+1, j)-peln(i, k, j&
&             ))*phiz_ad(i, k)
            peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad
            peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad
            phiz_ad(i, k) = 0.0_8
            CALL POPREAL8(tv(i, k))
            pt_ad(i, j, k) = pt_ad(i, j, k) + (r_vir*q(i, j, k, sphum)+&
&             1.)*tv_ad(i, k)
            q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) + pt(i, j, k)*&
&             r_vir*tv_ad(i, k)
            tv_ad(i, k) = 0.0_8
          END DO
        END DO
        DO i=ie,is,-1
          CALL POPREAL8(phiz(i, km+1))
          phiz_ad(i, km+1) = 0.0_8
        END DO
      END IF
    END DO
  END SUBROUTINE COMPUTE_TOTAL_ENERGY_ADM
!#ifdef MAPL_MODE
!#endif
  SUBROUTINE COMPUTE_TOTAL_ENERGY(is, ie, js, je, isd, ied, jsd, jed, km&
&   , u, v, w, delz, pt, delp, q, pe, peln, hs, grav, r_vir, cp, rg, hlv&
&   , te_2d, ua, teq, moist_phys, sphum, hydrostatic, id_te)
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
    REAL, DIMENSION(isd:ied, jsd:jed, km, sphum), INTENT(IN) :: q
    REAL, INTENT(INOUT) :: u(isd:ied, jsd:jed+1, km)
    REAL, INTENT(INOUT) :: v(isd:ied+1, jsd:jed, km)
! vertical velocity (m/s)
    REAL, INTENT(IN) :: w(isd:ied, jsd:jed, km)
    REAL, INTENT(IN) :: delz(is:ie, js:je, km)
! surface geopotential
    REAL, INTENT(IN) :: hs(isd:ied, jsd:jed)
! pressure at layer edges
    REAL(p_precision), INTENT(IN) :: pe(is-1:ie+1, km+1, js-1:je+1)
! log(pe)
    REAL(p_precision), INTENT(IN) :: peln(is:ie, km+1, js:je)
    REAL, INTENT(IN) :: cp, rg, r_vir, hlv
!#ifdef MAPL_MODE                 
    REAL, INTENT(IN) :: grav
!#endif
    LOGICAL, INTENT(IN) :: moist_phys, hydrostatic
! Output:
! vertically integrated TE
    REAL, INTENT(OUT) :: te_2d(is:ie, js:je)
! Moist TE
    REAL, INTENT(OUT) :: teq(is:ie, js:je)
! Local
    REAL(p_precision), DIMENSION(is:ie, km) :: tv
    REAL(p_precision) :: phiz(is:ie, km+1)
    REAL(p_precision) :: cv
    INTEGER :: i, j, k
    cv = cp - rg
!----------------------
! Output lat-lon winds:
!----------------------
!  call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km)
!$omp parallel do default(shared) private(phiz, tv)
    DO j=js,je
      IF (hydrostatic) THEN
        DO i=is,ie
          phiz(i, km+1) = hs(i, j)
        END DO
        DO k=km,1,-1
          DO i=is,ie
            tv(i, k) = pt(i, j, k)*(1.+r_vir*q(i, j, k, sphum))
            phiz(i, k) = phiz(i, k+1) + rg*tv(i, k)*(peln(i, k+1, j)-&
&             peln(i, k, j))
          END DO
        END DO
        DO i=is,ie
          te_2d(i, j) = pe(i, km+1, j)*phiz(i, km+1) - pe(i, 1, j)*phiz(&
&           i, 1)
        END DO
        DO k=1,km
          DO i=is,ie
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
          phiz(i, km+1) = hs(i, j)
          DO k=km,1,-1
            phiz(i, k) = phiz(i, k+1) - grav*delz(i, j, k)
          END DO
        END DO
        DO i=is,ie
          te_2d(i, j) = 0.
        END DO
        DO k=1,km
          DO i=is,ie
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
  END SUBROUTINE COMPUTE_TOTAL_ENERGY
!  Differentiation of pkez in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: peln pkz pk
!   with respect to varying inputs: peln pkz pk
  SUBROUTINE PKEZ_ADM(km, ifirst, ilast, jfirst, jlast, pe, pk, pk_ad, &
&   akap, peln, peln_ad, pkz, pkz_ad)
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
    REAL(p_precision) :: pk_ad(ifirst:ilast, jfirst:jlast, km+1)
! !OUTPUT
    REAL(p_precision) :: pkz(ifirst:ilast, jfirst:jlast, km)
    REAL(p_precision) :: pkz_ad(ifirst:ilast, jfirst:jlast, km)
! log (pe)
    REAL(p_precision), INTENT(INOUT) :: peln(ifirst:ilast, km+1, jfirst:&
&   jlast)
    REAL(p_precision) :: peln_ad(ifirst:ilast, km+1, jfirst:jlast)
! Local
    REAL(p_precision) :: pk2(ifirst:ilast, km+1)
    REAL(p_precision) :: pk2_ad(ifirst:ilast, km+1)
    REAL(p_precision) :: pek
    REAL(p_precision) :: pek_ad
    REAL(p_precision) :: lnp
    REAL(p_precision) :: ak1
    INTEGER :: i, j, k
    INTRINSIC LOG
    INTEGER :: branch
    REAL(p_precision) :: temp_ad0
    REAL(p_precision) :: temp_ad
    REAL*8 :: temp
    ak1 = (akap+1.)/akap
!$omp parallel do default(shared) private(lnp, pek, pk2)
    DO j=jfirst,jlast
      pek = pk(ifirst, j, 1)
      DO i=ifirst,ilast
        CALL PUSHREAL8(pk2(i, 1))
        pk2(i, 1) = pek
      END DO
      DO k=2,km+1
        DO i=ifirst,ilast
!             peln(i,k,j) =  log(pe(i,k,j))
          CALL PUSHREAL8(pk2(i, k))
          pk2(i, k) = pk(i, j, k)
        END DO
      END DO
!---- GFDL modification
      IF (ptop .LT. ptop_min) THEN
        DO i=ifirst,ilast
          peln(i, 1, j) = peln(i, 2, j) - ak1
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        lnp = LOG(ptop)
        DO i=ifirst,ilast
          peln(i, 1, j) = lnp
        END DO
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    pk2_ad = 0.0_8
    DO j=jlast,jfirst,-1
      DO k=km,1,-1
        DO i=ilast,ifirst,-1
          temp = akap*(peln(i, k+1, j)-peln(i, k, j))
          temp_ad = pkz_ad(i, j, k)/temp
          temp_ad0 = -((pk2(i, k+1)-pk2(i, k))*akap*temp_ad/temp)
          pk2_ad(i, k+1) = pk2_ad(i, k+1) + temp_ad
          pk2_ad(i, k) = pk2_ad(i, k) - temp_ad
          peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad0
          peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad0
          pkz_ad(i, j, k) = 0.0_8
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO i=ilast,ifirst,-1
          peln_ad(i, 1, j) = 0.0_8
        END DO
      ELSE
        DO i=ilast,ifirst,-1
          peln_ad(i, 2, j) = peln_ad(i, 2, j) + peln_ad(i, 1, j)
          peln_ad(i, 1, j) = 0.0_8
        END DO
      END IF
      DO k=km+1,2,-1
        DO i=ilast,ifirst,-1
          CALL POPREAL8(pk2(i, k))
          pk_ad(i, j, k) = pk_ad(i, j, k) + pk2_ad(i, k)
          pk2_ad(i, k) = 0.0_8
        END DO
      END DO
      pek_ad = 0.0_8
      DO i=ilast,ifirst,-1
        CALL POPREAL8(pk2(i, 1))
        pek_ad = pek_ad + pk2_ad(i, 1)
        pk2_ad(i, 1) = 0.0_8
      END DO
      pk_ad(ifirst, j, 1) = pk_ad(ifirst, j, 1) + pek_ad
    END DO
  END SUBROUTINE PKEZ_ADM
  SUBROUTINE PKEZ(km, ifirst, ilast, jfirst, jlast, pe, pk, akap, peln, &
&   pkz)
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
! !OUTPUT
    REAL(p_precision), INTENT(OUT) :: pkz(ifirst:ilast, jfirst:jlast, km&
&   )
! log (pe)
    REAL(p_precision), INTENT(INOUT) :: peln(ifirst:ilast, km+1, jfirst:&
&   jlast)
! Local
    REAL(p_precision) :: pk2(ifirst:ilast, km+1)
    REAL(p_precision) :: pek
    REAL(p_precision) :: lnp
    REAL(p_precision) :: ak1
    INTEGER :: i, j, k
    INTRINSIC LOG
    ak1 = (akap+1.)/akap
!$omp parallel do default(shared) private(lnp, pek, pk2)
    DO j=jfirst,jlast
      pek = pk(ifirst, j, 1)
      DO i=ifirst,ilast
        pk2(i, 1) = pek
      END DO
      DO k=2,km+1
        DO i=ifirst,ilast
!             peln(i,k,j) =  log(pe(i,k,j))
          pk2(i, k) = pk(i, j, k)
        END DO
      END DO
!---- GFDL modification
      IF (ptop .LT. ptop_min) THEN
        DO i=ifirst,ilast
          peln(i, 1, j) = peln(i, 2, j) - ak1
        END DO
      ELSE
        lnp = LOG(ptop)
        DO i=ifirst,ilast
          peln(i, 1, j) = lnp
        END DO
      END IF
!---- GFDL modification
      DO k=1,km
        DO i=ifirst,ilast
          pkz(i, j, k) = (pk2(i, k+1)-pk2(i, k))/(akap*(peln(i, k+1, j)-&
&           peln(i, k, j)))
        END DO
      END DO
    END DO
  END SUBROUTINE PKEZ
!  Differentiation of map1_ppm_nolimiters in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: pe1 pe2 q2
!   with respect to varying inputs: pe1 pe2 q2
  SUBROUTINE MAP1_PPM_NOLIMITERS_ADM(km, pe1, pe1_ad, pe2, pe2_ad, q2, &
&   q2_ad, i1, i2, j, ibeg, iend, jbeg, jend, iv, kord)
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
    REAL(p_precision) :: pe1_ad(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, km+1)
    REAL(p_precision) :: pe2_ad(i1:i2, km+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(INOUT) :: q2_ad(ibeg:iend, jbeg:jend, km)
! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
    REAL(p_precision) :: dp1(i1:i2, km)
    REAL(p_precision) :: dp1_ad(i1:i2, km)
    REAL(p_precision) :: q4(4, i1:i2, km)
    REAL(p_precision) :: q4_ad(4, i1:i2, km)
    REAL(p_precision) :: pl, pr, qsum, dp, esl
    REAL(p_precision) :: pl_ad, pr_ad, qsum_ad, dp_ad, esl_ad
    INTEGER :: i, k, l, m, k0
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_count0
    INTEGER :: i3
    REAL(p_precision) :: temp1
    REAL(p_precision) :: temp0
    REAL(p_precision) :: temp_ad9
    REAL(p_precision) :: temp_ad8
    REAL(p_precision) :: temp_ad7
    REAL(p_precision) :: temp_ad6
    REAL(p_precision) :: temp_ad5
    REAL(p_precision) :: temp_ad4
    REAL(p_precision) :: temp_ad3
    REAL(p_precision) :: temp_ad2
    REAL(p_precision) :: temp_ad1
    REAL(p_precision) :: temp_ad0
    REAL(p_precision) :: temp_ad
    REAL(p_precision) :: temp_ad11
    REAL(p_precision) :: temp_ad10
    REAL(p_precision) :: temp
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL PUSHREAL8ARRAY(q4, 4*(i2-i1+1)*km)
      CALL CS_PROFILE_NOLIMITERS(q4, dp1, km, i1, i2, iv, kord)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
!call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
    DO i=i1,i2
      k0 = 1
      DO 120 k=1,km
        CALL PUSHINTEGER4(l)
        ad_count = 1
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            GOTO 100
          ELSE
            CALL PUSHINTEGER4(l)
            ad_count = ad_count + 1
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
        CALL PUSHINTEGER4(ad_count)
        CALL PUSHCONTROL2B(2)
        GOTO 123
 100    CALL PUSHCONTROL1B(1)
        CALL PUSHINTEGER4(ad_count)
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          k0 = l
          CALL PUSHCONTROL1B(0)
          GOTO 120
        ELSE
! Fractional area...
          CALL PUSHREAL8(qsum)
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          CALL PUSHINTEGER4(m)
          ad_count0 = 1
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
              qsum = qsum + dp1(i, m)*q4(1, i, m)
              CALL PUSHINTEGER4(m)
              ad_count0 = ad_count0 + 1
            ELSE
              GOTO 110
            END IF
          END DO
          CALL PUSHCONTROL1B(0)
          CALL PUSHINTEGER4(ad_count0)
          CALL PUSHCONTROL2B(1)
          GOTO 123
 110      CALL PUSHCONTROL1B(1)
          CALL PUSHINTEGER4(ad_count0)
          dp = pe2(i, k+1) - pe1(i, m)
          esl = dp/dp1(i, m)
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
          CALL PUSHCONTROL2B(0)
        END IF
 123    CALL PUSHCONTROL1B(1)
 120  CONTINUE
    END DO
    dp1_ad = 0.0_8
    qsum_ad = 0.0_8
    q4_ad = 0.0_8
    DO i=i2,i1,-1
      DO k=km,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          temp_ad0 = 0.5*(pr+pl)*q2_ad(i, j, k)
          temp_ad1 = 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i, l))*q2_ad(i, &
&           j, k)
          temp_ad2 = -(r3*q4(4, i, l)*q2_ad(i, j, k))
          q4_ad(2, i, l) = q4_ad(2, i, l) + q2_ad(i, j, k) - temp_ad0
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad0 - r3*(pr*(pr+pl)+pl&
&           **2)*q2_ad(i, j, k)
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad0
          pr_ad = (2*pr+pl)*temp_ad2 + temp_ad1
          pl_ad = (2*pl+pr)*temp_ad2 + temp_ad1
          q2_ad(i, j, k) = 0.0_8
          temp_ad3 = pr_ad/dp1(i, l)
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad3
          pe1_ad(i, l) = pe1_ad(i, l) - temp_ad3
          dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k+1)-pe1(i, l))*temp_ad3&
&           /dp1(i, l)
        ELSE
          temp1 = pe2(i, k+1) - pe2(i, k)
          temp_ad11 = -(qsum*q2_ad(i, j, k)/temp1**2)
          qsum_ad = qsum_ad + q2_ad(i, j, k)/temp1
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad11
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad11
          q2_ad(i, j, k) = 0.0_8
          CALL POPCONTROL2B(branch)
          IF (branch .EQ. 0) THEN
            dp = pe2(i, k+1) - pe1(i, m)
            esl = dp/dp1(i, m)
            temp0 = q4(3, i, m) - q4(2, i, m) + q4(4, i, m)*(-(r23*esl)+&
&             1.)
            temp_ad8 = dp*qsum_ad
            temp_ad9 = 0.5*esl*temp_ad8
            q4_ad(2, i, m) = q4_ad(2, i, m) + temp_ad8 - temp_ad9
            esl_ad = 0.5*temp0*temp_ad8 - q4(4, i, m)*r23*temp_ad9
            q4_ad(3, i, m) = q4_ad(3, i, m) + temp_ad9
            q4_ad(4, i, m) = q4_ad(4, i, m) + (1.-r23*esl)*temp_ad9
            temp_ad10 = esl_ad/dp1(i, m)
            dp_ad = temp_ad10 + (q4(2, i, m)+0.5*(esl*temp0))*qsum_ad
            dp1_ad(i, m) = dp1_ad(i, m) - dp*temp_ad10/dp1(i, m)
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp_ad
            pe1_ad(i, m) = pe1_ad(i, m) - dp_ad
          ELSE IF (branch .NE. 1) THEN
            GOTO 130
          END IF
          CALL POPINTEGER4(ad_count0)
          DO i3=1,ad_count0
            IF (i3 .EQ. 1) THEN
              CALL POPCONTROL1B(branch)
            ELSE
              dp1_ad(i, m) = dp1_ad(i, m) + q4(1, i, m)*qsum_ad
              q4_ad(1, i, m) = q4_ad(1, i, m) + dp1(i, m)*qsum_ad
            END IF
            CALL POPINTEGER4(m)
          END DO
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          CALL POPREAL8(qsum)
          temp = q4(4, i, l) + q4(3, i, l) - q4(2, i, l)
          temp_ad4 = (q4(2, i, l)+0.5*(temp*(pl+1.))-r3*(q4(4, i, l)*(pl&
&           *(pl+1.)+1.)))*qsum_ad
          temp_ad5 = (pe1(i, l+1)-pe2(i, k))*qsum_ad
          temp_ad6 = 0.5*(pl+1.)*temp_ad5
          temp_ad7 = -(r3*q4(4, i, l)*temp_ad5)
          pe1_ad(i, l+1) = pe1_ad(i, l+1) + temp_ad4
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad4
          q4_ad(2, i, l) = q4_ad(2, i, l) + temp_ad5 - temp_ad6
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad6 - r3*(pl*(pl+1.)+1.&
&           )*temp_ad5
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad6
          pl_ad = (2*pl+1.)*temp_ad7 + 0.5*temp*temp_ad5
          qsum_ad = 0.0_8
        END IF
        temp_ad = pl_ad/dp1(i, l)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad
        pe1_ad(i, l) = pe1_ad(i, l) - temp_ad
        dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k)-pe1(i, l))*temp_ad/dp1(&
&         i, l)
 130    CALL POPINTEGER4(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL1B(branch)
          CALL POPINTEGER4(l)
        END DO
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      CALL POPREAL8ARRAY(q4, 4*(i2-i1+1)*km)
      CALL CS_PROFILE_NOLIMITERS_ADM(q4, q4_ad, dp1, dp1_ad, km, i1, i2&
&                              , iv, kord)
    END IF
    DO k=km,1,-1
      DO i=i2,i1,-1
        q2_ad(i, j, k) = q2_ad(i, j, k) + q4_ad(1, i, k)
        q4_ad(1, i, k) = 0.0_8
        pe1_ad(i, k+1) = pe1_ad(i, k+1) + dp1_ad(i, k)
        pe1_ad(i, k) = pe1_ad(i, k) - dp1_ad(i, k)
        dp1_ad(i, k) = 0.0_8
      END DO
    END DO
  END SUBROUTINE MAP1_PPM_NOLIMITERS_ADM
  SUBROUTINE MAP1_PPM_NOLIMITERS(km, pe1, pe2, q2, i1, i2, j, ibeg, iend&
&   , jbeg, jend, iv, kord)
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
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, km+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, km)
! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
    REAL(p_precision) :: dp1(i1:i2, km)
    REAL(p_precision) :: q4(4, i1:i2, km)
    REAL(p_precision) :: pl, pr, qsum, dp, esl
    INTEGER :: i, k, l, m, k0
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q2(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) CALL CS_PROFILE_NOLIMITERS(q4, dp1, km, i1, i2, iv&
&                                         , kord)
!call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
    DO i=i1,i2
      k0 = 1
      DO k=1,km
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
            IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
              pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
              q2(i, j, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-&
&               q4(2, i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
              k0 = l
              GOTO 555
            ELSE
! Fractional area...
              qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, &
&               l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+&
&               pl*(1.+pl))))
              DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer
                  qsum = qsum + dp1(i, m)*q4(1, i, m)
                ELSE
                  dp = pe2(i, k+1) - pe1(i, m)
                  esl = dp/dp1(i, m)
                  qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(&
&                   2, i, m)+q4(4, i, m)*(1.-r23*esl)))
                  k0 = m
                  GOTO 123
                END IF
              END DO
              GOTO 123
            END IF
          END IF
        END DO
 123    q2(i, j, k) = qsum/(pe2(i, k+1)-pe2(i, k))
 555    CONTINUE
      END DO
    END DO
  END SUBROUTINE MAP1_PPM_NOLIMITERS
!  Differentiation of map1_q2_nolimiters in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: pe1 pe2 dp2 q1 q2
!   with respect to varying inputs: pe1 pe2 dp2 q1 q2
  SUBROUTINE MAP1_Q2_NOLIMITERS_ADM(km, pe1, pe1_ad, q1, q1_ad, kn, pe2&
&   , pe2_ad, q2, q2_ad, dp2, dp2_ad, i1, i2, iv, kord, j, ibeg, iend, &
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
    REAL(p_precision) :: pe1_ad(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, kn+1)
    REAL(p_precision) :: pe2_ad(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! Field input
    REAL, INTENT(IN) :: q1(ibeg:iend, jbeg:jend, km)
    REAL :: q1_ad(ibeg:iend, jbeg:jend, km)
    REAL(p_precision), INTENT(IN) :: dp2(i1:i2, kn)
    REAL(p_precision) :: dp2_ad(i1:i2, kn)
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL(p_precision), INTENT(INOUT) :: q2(i1:i2, kn)
    REAL(p_precision) :: q2_ad(i1:i2, kn)
! !LOCAL VARIABLES:
    REAL(p_precision) :: dp1(i1:i2, km)
    REAL(p_precision) :: dp1_ad(i1:i2, km)
    REAL(p_precision) :: q4(4, i1:i2, km)
    REAL(p_precision) :: q4_ad(4, i1:i2, km)
    REAL(p_precision) :: pl, pr, qsum, dp, esl
    REAL(p_precision) :: pl_ad, pr_ad, qsum_ad, dp_ad, esl_ad
    INTEGER :: i, k, l, m, k0
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    INTEGER :: ad_count0
    INTEGER :: i3
    REAL(p_precision) :: temp0
    REAL(p_precision) :: temp_ad9
    REAL(p_precision) :: temp_ad8
    REAL(p_precision) :: temp_ad7
    REAL(p_precision) :: temp_ad6
    REAL(p_precision) :: temp_ad5
    REAL(p_precision) :: temp_ad4
    REAL(p_precision) :: temp_ad3
    REAL(p_precision) :: temp_ad2
    REAL(p_precision) :: temp_ad1
    REAL(p_precision) :: temp_ad0
    REAL(p_precision) :: temp_ad
    REAL(p_precision) :: temp_ad11
    REAL(p_precision) :: temp_ad10
    REAL(p_precision) :: temp
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q1(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) THEN
      CALL PUSHREAL8ARRAY(q4, 4*(i2-i1+1)*km)
      CALL CS_PROFILE_NOLIMITERS(q4, dp1, km, i1, i2, iv, kord)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
!call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
! Mapping
    DO i=i1,i2
      k0 = 1
      DO 120 k=1,kn
        CALL PUSHINTEGER4(l)
        ad_count = 1
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            GOTO 100
          ELSE
            CALL PUSHINTEGER4(l)
            ad_count = ad_count + 1
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
        CALL PUSHINTEGER4(ad_count)
        CALL PUSHCONTROL2B(2)
        GOTO 123
 100    CALL PUSHCONTROL1B(1)
        CALL PUSHINTEGER4(ad_count)
        pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
        IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
          k0 = l
          CALL PUSHCONTROL1B(0)
          GOTO 120
        ELSE
! Fractional area...
          CALL PUSHREAL8(qsum)
          qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+&
&           q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+pl*(1.+&
&           pl))))
          CALL PUSHINTEGER4(m)
          ad_count0 = 1
          DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
            IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
              qsum = qsum + dp1(i, m)*q4(1, i, m)
              CALL PUSHINTEGER4(m)
              ad_count0 = ad_count0 + 1
            ELSE
              GOTO 110
            END IF
          END DO
          CALL PUSHCONTROL1B(0)
          CALL PUSHINTEGER4(ad_count0)
          CALL PUSHCONTROL2B(1)
          GOTO 123
 110      CALL PUSHCONTROL1B(1)
          CALL PUSHINTEGER4(ad_count0)
          dp = pe2(i, k+1) - pe1(i, m)
          esl = dp/dp1(i, m)
          qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(2, i, m)&
&           +q4(4, i, m)*(1.-r23*esl)))
          k0 = m
          CALL PUSHCONTROL2B(0)
        END IF
 123    CALL PUSHCONTROL1B(1)
 120  CONTINUE
    END DO
    dp1_ad = 0.0_8
    qsum_ad = 0.0_8
    q4_ad = 0.0_8
    DO i=i2,i1,-1
      DO k=kn,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
          temp_ad0 = 0.5*(pr+pl)*q2_ad(i, k)
          temp_ad1 = 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i, l))*q2_ad(i, &
&           k)
          temp_ad2 = -(r3*q4(4, i, l)*q2_ad(i, k))
          q4_ad(2, i, l) = q4_ad(2, i, l) + q2_ad(i, k) - temp_ad0
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad0 - r3*(pr*(pr+pl)+pl&
&           **2)*q2_ad(i, k)
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad0
          pr_ad = (2*pr+pl)*temp_ad2 + temp_ad1
          pl_ad = (2*pl+pr)*temp_ad2 + temp_ad1
          q2_ad(i, k) = 0.0_8
          temp_ad3 = pr_ad/dp1(i, l)
          pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad3
          pe1_ad(i, l) = pe1_ad(i, l) - temp_ad3
          dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k+1)-pe1(i, l))*temp_ad3&
&           /dp1(i, l)
        ELSE
          temp_ad11 = q2_ad(i, k)/dp2(i, k)
          qsum_ad = qsum_ad + temp_ad11
          dp2_ad(i, k) = dp2_ad(i, k) - qsum*temp_ad11/dp2(i, k)
          q2_ad(i, k) = 0.0_8
          CALL POPCONTROL2B(branch)
          IF (branch .EQ. 0) THEN
            dp = pe2(i, k+1) - pe1(i, m)
            esl = dp/dp1(i, m)
            temp0 = q4(3, i, m) - q4(2, i, m) + q4(4, i, m)*(-(r23*esl)+&
&             1.)
            temp_ad8 = dp*qsum_ad
            temp_ad9 = 0.5*esl*temp_ad8
            q4_ad(2, i, m) = q4_ad(2, i, m) + temp_ad8 - temp_ad9
            esl_ad = 0.5*temp0*temp_ad8 - q4(4, i, m)*r23*temp_ad9
            q4_ad(3, i, m) = q4_ad(3, i, m) + temp_ad9
            q4_ad(4, i, m) = q4_ad(4, i, m) + (1.-r23*esl)*temp_ad9
            temp_ad10 = esl_ad/dp1(i, m)
            dp_ad = temp_ad10 + (q4(2, i, m)+0.5*(esl*temp0))*qsum_ad
            dp1_ad(i, m) = dp1_ad(i, m) - dp*temp_ad10/dp1(i, m)
            pe2_ad(i, k+1) = pe2_ad(i, k+1) + dp_ad
            pe1_ad(i, m) = pe1_ad(i, m) - dp_ad
          ELSE IF (branch .NE. 1) THEN
            GOTO 130
          END IF
          CALL POPINTEGER4(ad_count0)
          DO i3=1,ad_count0
            IF (i3 .EQ. 1) THEN
              CALL POPCONTROL1B(branch)
            ELSE
              dp1_ad(i, m) = dp1_ad(i, m) + q4(1, i, m)*qsum_ad
              q4_ad(1, i, m) = q4_ad(1, i, m) + dp1(i, m)*qsum_ad
            END IF
            CALL POPINTEGER4(m)
          END DO
          pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
          CALL POPREAL8(qsum)
          temp = q4(4, i, l) + q4(3, i, l) - q4(2, i, l)
          temp_ad4 = (q4(2, i, l)+0.5*(temp*(pl+1.))-r3*(q4(4, i, l)*(pl&
&           *(pl+1.)+1.)))*qsum_ad
          temp_ad5 = (pe1(i, l+1)-pe2(i, k))*qsum_ad
          temp_ad6 = 0.5*(pl+1.)*temp_ad5
          temp_ad7 = -(r3*q4(4, i, l)*temp_ad5)
          pe1_ad(i, l+1) = pe1_ad(i, l+1) + temp_ad4
          pe2_ad(i, k) = pe2_ad(i, k) - temp_ad4
          q4_ad(2, i, l) = q4_ad(2, i, l) + temp_ad5 - temp_ad6
          q4_ad(4, i, l) = q4_ad(4, i, l) + temp_ad6 - r3*(pl*(pl+1.)+1.&
&           )*temp_ad5
          q4_ad(3, i, l) = q4_ad(3, i, l) + temp_ad6
          pl_ad = (2*pl+1.)*temp_ad7 + 0.5*temp*temp_ad5
          qsum_ad = 0.0_8
        END IF
        temp_ad = pl_ad/dp1(i, l)
        pe2_ad(i, k) = pe2_ad(i, k) + temp_ad
        pe1_ad(i, l) = pe1_ad(i, l) - temp_ad
        dp1_ad(i, l) = dp1_ad(i, l) - (pe2(i, k)-pe1(i, l))*temp_ad/dp1(&
&         i, l)
 130    CALL POPINTEGER4(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL1B(branch)
          CALL POPINTEGER4(l)
        END DO
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      CALL POPREAL8ARRAY(q4, 4*(i2-i1+1)*km)
      CALL CS_PROFILE_NOLIMITERS_ADM(q4, q4_ad, dp1, dp1_ad, km, i1, i2&
&                              , iv, kord)
    END IF
    DO k=km,1,-1
      DO i=i2,i1,-1
        q1_ad(i, j, k) = q1_ad(i, j, k) + q4_ad(1, i, k)
        q4_ad(1, i, k) = 0.0_8
        pe1_ad(i, k+1) = pe1_ad(i, k+1) + dp1_ad(i, k)
        pe1_ad(i, k) = pe1_ad(i, k) - dp1_ad(i, k)
        dp1_ad(i, k) = 0.0_8
      END DO
    END DO
  END SUBROUTINE MAP1_Q2_NOLIMITERS_ADM
  SUBROUTINE MAP1_Q2_NOLIMITERS(km, pe1, q1, kn, pe2, q2, dp2, i1, i2, &
&   iv, kord, j, ibeg, iend, jbeg, jend)
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
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, kn+1)
! (from model top to bottom surface)
! in the new vertical coordinate
! Field input
    REAL, INTENT(IN) :: q1(ibeg:iend, jbeg:jend, km)
    REAL(p_precision), INTENT(IN) :: dp2(i1:i2, kn)
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL(p_precision), INTENT(INOUT) :: q2(i1:i2, kn)
! !LOCAL VARIABLES:
    REAL(p_precision) :: dp1(i1:i2, km)
    REAL(p_precision) :: q4(4, i1:i2, km)
    REAL(p_precision) :: pl, pr, qsum, dp, esl
    INTEGER :: i, k, l, m, k0
    DO k=1,km
      DO i=i1,i2
        dp1(i, k) = pe1(i, k+1) - pe1(i, k)
        q4(1, i, k) = q1(i, j, k)
      END DO
    END DO
! Compute vertical subgrid distribution
    IF (kord .GT. 7) CALL CS_PROFILE_NOLIMITERS(q4, dp1, km, i1, i2, iv&
&                                         , kord)
!call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
! Mapping
    DO i=i1,i2
      k0 = 1
      DO k=1,kn
        DO l=k0,km
! locate the top edge: pe2(i,k)
          IF (pe2(i, k) .GE. pe1(i, l) .AND. pe2(i, k) .LE. pe1(i, l+1)&
&         ) THEN
            pl = (pe2(i, k)-pe1(i, l))/dp1(i, l)
            IF (pe2(i, k+1) .LE. pe1(i, l+1)) THEN
! entire new grid is within the original grid
              pr = (pe2(i, k+1)-pe1(i, l))/dp1(i, l)
              q2(i, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2&
&               , i, l))*(pr+pl) - q4(4, i, l)*r3*(pr*(pr+pl)+pl**2)
              k0 = l
              GOTO 555
            ELSE
! Fractional area...
              qsum = (pe1(i, l+1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, &
&               l)+q4(3, i, l)-q4(2, i, l))*(1.+pl)-q4(4, i, l)*(r3*(1.+&
&               pl*(1.+pl))))
              DO m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                IF (pe2(i, k+1) .GT. pe1(i, m+1)) THEN
! Whole layer..
                  qsum = qsum + dp1(i, m)*q4(1, i, m)
                ELSE
                  dp = pe2(i, k+1) - pe1(i, m)
                  esl = dp/dp1(i, m)
                  qsum = qsum + dp*(q4(2, i, m)+0.5*esl*(q4(3, i, m)-q4(&
&                   2, i, m)+q4(4, i, m)*(1.-r23*esl)))
                  k0 = m
                  GOTO 123
                END IF
              END DO
              GOTO 123
            END IF
          END IF
        END DO
 123    q2(i, k) = qsum/dp2(i, k)
 555    CONTINUE
      END DO
    END DO
  END SUBROUTINE MAP1_Q2_NOLIMITERS
!  Differentiation of map1_cubic_te in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: pe1 pe2 q2
!   with respect to varying inputs: pe1 pe2 q2
!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  map1_cubic_te --- Cubic Interpolation for TE mapping
!
! !INTERFACE:
  SUBROUTINE MAP1_CUBIC_TE_ADM(km, pe1, pe1_ad, pe2, pe2_ad, q2, q2_ad, &
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
    REAL(p_precision) :: pe1_ad(i1:i2, km+1)
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, km+1)
    REAL(p_precision) :: pe2_ad(i1:i2, km+1)
! (from model top to bottom surface)
! in the new vertical coordinate
!      real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, km)
    REAL, INTENT(INOUT) :: q2_ad(ibeg:iend, jbeg:jend, km)
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
    REAL(p_precision) :: qx_ad(i1:i2, km)
    REAL(p_precision) :: logpl1(i1:i2, km)
    REAL(p_precision) :: logpl1_ad(i1:i2, km)
    REAL(p_precision) :: logpl2(i1:i2, km)
    REAL(p_precision) :: logpl2_ad(i1:i2, km)
    REAL(p_precision) :: dlogp1(i1:i2, km)
    REAL(p_precision) :: dlogp1_ad(i1:i2, km)
    REAL(p_precision) :: vsum1(i1:i2)
    REAL(p_precision) :: vsum1_ad(i1:i2)
    REAL(p_precision) :: vsum2(i1:i2)
    REAL(p_precision) :: vsum2_ad(i1:i2)
    REAL(p_precision) :: am2, am1, ap0, ap1, p, plp1, plp0, plm1, plm2, &
&   dlp0, dlm1, dlm2
    REAL(p_precision) :: am2_ad, am1_ad, ap0_ad, ap1_ad, p_ad, plp1_ad, &
&   plp0_ad, plm1_ad, plm2_ad, dlp0_ad, dlm1_ad, dlm2_ad
    INTEGER :: i, k, lm2, lm1, lp0, lp1
    INTRINSIC LOG
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    REAL(p_precision) :: temp3
    REAL(p_precision) :: temp2
    REAL(p_precision) :: temp1
    REAL(p_precision) :: temp0
    REAL(p_precision) :: temp19
    REAL(p_precision) :: temp18
    REAL(p_precision) :: temp17
    REAL(p_precision) :: temp16
    REAL(p_precision) :: temp15
    REAL(p_precision) :: temp14
    REAL(p_precision) :: temp13
    REAL(p_precision) :: temp12
    REAL(p_precision) :: temp11
    REAL(p_precision) :: temp10
    REAL(p_precision) :: temp_ad26
    REAL(p_precision) :: temp_ad25
    REAL(p_precision) :: temp_ad24
    REAL(p_precision) :: temp_ad23
    REAL(p_precision) :: temp_ad22
    REAL(p_precision) :: temp_ad21
    REAL(p_precision) :: temp_ad20
    REAL(p_precision) :: temp_ad9
    REAL(p_precision) :: temp_ad8
    REAL(p_precision) :: temp_ad7
    REAL(p_precision) :: temp_ad6
    REAL(p_precision) :: temp_ad5
    REAL(p_precision) :: temp_ad4
    REAL(p_precision) :: temp_ad3
    REAL(p_precision) :: temp_ad2
    REAL(p_precision) :: temp_ad1
    REAL(p_precision) :: temp_ad0(i2-i1+1)
    REAL(p_precision) :: temp_ad(i2-i1+1)
    REAL(p_precision) :: temp_ad19
    REAL(p_precision) :: temp_ad18
    REAL(p_precision) :: temp_ad17
    REAL(p_precision) :: temp_ad16
    REAL(p_precision) :: temp_ad15
    REAL(p_precision) :: temp_ad14
    REAL(p_precision) :: temp_ad13
    REAL(p_precision) :: temp_ad12
    REAL(p_precision) :: temp_ad11
    REAL(p_precision) :: temp_ad10
    REAL(p_precision) :: temp
    REAL(p_precision) :: temp9
    REAL(p_precision) :: temp8
    REAL(p_precision) :: temp7
    REAL(p_precision) :: temp6
    REAL(p_precision) :: temp5
    REAL(p_precision) :: temp4
! Initialization
! --------------
    DO k=1,km
      qx(:, k) = q2(:, j, k)
      logpl1(:, k) = LOG(d0_5*(pe1(:, k)+pe1(:, k+1)))
    END DO
    DO k=1,km
      logpl2(:, k) = LOG(d0_5*(pe2(:, k)+pe2(:, k+1)))
    END DO
    DO k=1,km-1
      dlogp1(:, k) = logpl1(:, k+1) - logpl1(:, k)
    END DO
! Compute vertical integral of Input TE
! -------------------------------------
    vsum1(:) = d0_0
    DO i=i1,i2
      DO k=1,km
        vsum1(i) = vsum1(i) + qx(i, k)*(pe1(i, k+1)-pe1(i, k))
      END DO
    END DO
! Interpolate TE onto target Pressures
! ------------------------------------
    DO i=i1,i2
      DO k=1,km
        CALL PUSHINTEGER4(lp0)
        lp0 = 1
        ad_count = 1
        DO WHILE (lp0 .LE. km)
          IF (logpl1(i, lp0) .LT. logpl2(i, k)) THEN
            lp0 = lp0 + 1
            ad_count = ad_count + 1
          ELSE
            GOTO 100
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
        CALL PUSHINTEGER4(ad_count)
        GOTO 110
 100    CALL PUSHCONTROL1B(1)
        CALL PUSHINTEGER4(ad_count)
 110    IF (lp0 - 1 .LT. 1) THEN
          CALL PUSHINTEGER4(lm1)
          lm1 = 1
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHINTEGER4(lm1)
          lm1 = lp0 - 1
          CALL PUSHCONTROL1B(1)
        END IF
        IF (lp0 .GT. km) THEN
          lp0 = km
        ELSE
          lp0 = lp0
        END IF
! Extrapolate Linearly in LogP above first model level
! ----------------------------------------------------
        IF (lm1 .EQ. 1 .AND. lp0 .EQ. 1) THEN
          q2(i, j, k) = qx(i, 1) + (qx(i, 2)-qx(i, 1))*(logpl2(i, k)-&
&           logpl1(i, 1))/(logpl1(i, 2)-logpl1(i, 1))
! Extrapolate Linearly in LogP below last model level
! ---------------------------------------------------
          CALL PUSHCONTROL2B(3)
        ELSE IF (lm1 .EQ. km .AND. lp0 .EQ. km) THEN
          q2(i, j, k) = qx(i, km) + (qx(i, km)-qx(i, km-1))*(logpl2(i, k&
&           )-logpl1(i, km))/(logpl1(i, km)-logpl1(i, km-1))
! Interpolate Linearly in LogP between levels 1 => 2 and km-1 => km
! -----------------------------------------------------------------
          CALL PUSHCONTROL2B(2)
        ELSE IF (lm1 .EQ. 1 .OR. lp0 .EQ. km) THEN
          q2(i, j, k) = qx(i, lp0) + (qx(i, lm1)-qx(i, lp0))*(logpl2(i, &
&           k)-logpl1(i, lp0))/(logpl1(i, lm1)-logpl1(i, lp0))
! Interpolate Cubicly in LogP between other model levels
! ------------------------------------------------------
          CALL PUSHCONTROL2B(1)
        ELSE
          CALL PUSHINTEGER4(lp1)
          lp1 = lp0 + 1
          CALL PUSHINTEGER4(lm2)
          lm2 = lm1 - 1
          p = logpl2(i, k)
          plp1 = logpl1(i, lp1)
          plp0 = logpl1(i, lp0)
          plm1 = logpl1(i, lm1)
          plm2 = logpl1(i, lm2)
          dlp0 = dlogp1(i, lp0)
          dlm1 = dlogp1(i, lm1)
          dlm2 = dlogp1(i, lm2)
          CALL PUSHREAL8(ap1)
          ap1 = (p-plp0)*(p-plm1)*(p-plm2)/(dlp0*(dlp0+dlm1)*(dlp0+dlm1+&
&           dlm2))
          CALL PUSHREAL8(ap0)
          ap0 = (plp1-p)*(p-plm1)*(p-plm2)/(dlp0*dlm1*(dlm1+dlm2))
          CALL PUSHREAL8(am1)
          am1 = (plp1-p)*(plp0-p)*(p-plm2)/(dlm1*dlm2*(dlp0+dlm1))
          CALL PUSHREAL8(am2)
          am2 = (plp1-p)*(plp0-p)*(plm1-p)/(dlm2*(dlm1+dlm2)*(dlp0+dlm1+&
&           dlm2))
          q2(i, j, k) = ap1*qx(i, lp1) + ap0*qx(i, lp0) + am1*qx(i, lm1)&
&           + am2*qx(i, lm2)
          CALL PUSHCONTROL2B(0)
        END IF
      END DO
    END DO
! Compute vertical integral of Output TE
! --------------------------------------
    vsum2(:) = d0_0
    DO i=i1,i2
      DO k=1,km
        vsum2(i) = vsum2(i) + q2(i, j, k)*(pe2(i, k+1)-pe2(i, k))
      END DO
    END DO
    vsum1_ad = 0.0_8
    vsum2_ad = 0.0_8
    DO i=i2,i1,-1
      DO k=km,1,-1
        vsum1_ad(i) = vsum1_ad(i) + q2_ad(i, j, k)
        vsum2_ad(i) = vsum2_ad(i) - q2_ad(i, j, k)
      END DO
    END DO
    DO i=i2,i1,-1
      temp19 = pe2(i, km+1) - pe2(i, 1)
      temp_ad26 = -(vsum2(i)*vsum2_ad(i)/temp19**2)
      pe2_ad(i, km+1) = pe2_ad(i, km+1) + temp_ad26
      pe2_ad(i, 1) = pe2_ad(i, 1) - temp_ad26
      vsum2_ad(i) = vsum2_ad(i)/temp19
      DO k=km,1,-1
        temp_ad25 = q2(i, j, k)*vsum2_ad(i)
        q2_ad(i, j, k) = q2_ad(i, j, k) + (pe2(i, k+1)-pe2(i, k))*&
&         vsum2_ad(i)
        pe2_ad(i, k+1) = pe2_ad(i, k+1) + temp_ad25
        pe2_ad(i, k) = pe2_ad(i, k) - temp_ad25
      END DO
    END DO
    qx_ad = 0.0_8
    logpl1_ad = 0.0_8
    logpl2_ad = 0.0_8
    dlogp1_ad = 0.0_8
    DO i=i2,i1,-1
      DO k=km,1,-1
        CALL POPCONTROL2B(branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            ap1_ad = qx(i, lp1)*q2_ad(i, j, k)
            qx_ad(i, lp1) = qx_ad(i, lp1) + ap1*q2_ad(i, j, k)
            ap0_ad = qx(i, lp0)*q2_ad(i, j, k)
            qx_ad(i, lp0) = qx_ad(i, lp0) + ap0*q2_ad(i, j, k)
            am1_ad = qx(i, lm1)*q2_ad(i, j, k)
            qx_ad(i, lm1) = qx_ad(i, lm1) + am1*q2_ad(i, j, k)
            am2_ad = qx(i, lm2)*q2_ad(i, j, k)
            qx_ad(i, lm2) = qx_ad(i, lm2) + am2*q2_ad(i, j, k)
            q2_ad(i, j, k) = 0.0_8
            plp0 = logpl1(i, lp0)
            plp1 = logpl1(i, lp1)
            dlp0 = dlogp1(i, lp0)
            plm1 = logpl1(i, lm1)
            p = logpl2(i, k)
            dlm1 = dlogp1(i, lm1)
            dlm2 = dlogp1(i, lm2)
            CALL POPREAL8(am2)
            temp18 = dlm2*(dlm1+dlm2)
            temp17 = temp18*(dlp0+dlm1+dlm2)
            temp_ad9 = am2_ad/temp17
            temp_ad10 = (plm1-p)*temp_ad9
            temp16 = (plp1-p)*(plp0-p)
            temp_ad11 = -(temp16*(plm1-p)*temp_ad9/temp17)
            temp_ad12 = (dlp0+dlm1+dlm2)*temp_ad11
            temp_ad13 = temp18*temp_ad11
            plm2 = logpl1(i, lm2)
            CALL POPREAL8(am1)
            temp15 = dlm1*dlm2
            temp_ad16 = am1_ad/(temp15*(dlp0+dlm1))
            temp_ad14 = (p-plm2)*temp_ad16
            temp14 = (plp1-p)*(plp0-p)
            temp_ad20 = -(temp14*(p-plm2)*temp_ad16/(temp15*(dlp0+dlm1))&
&             )
            CALL POPREAL8(ap0)
            temp13 = dlp0*dlm1
            temp_ad19 = ap0_ad/(temp13*(dlm1+dlm2))
            temp_ad15 = (p-plm2)*temp_ad19
            plp1_ad = (plp0-p)*temp_ad14 + (p-plm1)*temp_ad15 + (plp0-p)&
&             *temp_ad10
            temp12 = (plp1-p)*(p-plm1)
            temp_ad22 = -(temp12*(p-plm2)*temp_ad19/(temp13*(dlm1+dlm2))&
&             )
            CALL POPREAL8(ap1)
            temp11 = dlp0*(dlp0+dlm1)
            temp10 = temp11*(dlp0+dlm1+dlm2)
            temp_ad18 = ap1_ad/temp10
            temp_ad17 = (p-plm2)*temp_ad18
            plp0_ad = (plp1-p)*temp_ad14 - (p-plm1)*temp_ad17 + (plp1-p)&
&             *temp_ad10
            plm1_ad = temp16*temp_ad9 - (p-plp0)*temp_ad17 - (plp1-p)*&
&             temp_ad15
            temp9 = (p-plp0)*(p-plm1)
            p_ad = (2*p-plp1-plp0)*temp_ad14 + temp14*temp_ad16 + (2*p-&
&             plp0-plm1)*temp_ad17 + temp9*temp_ad18 + temp12*temp_ad19 &
&             + (plp1-2*p+plm1)*temp_ad15 - temp16*temp_ad9 + (2*p-plp1-&
&             plp0)*temp_ad10
            plm2_ad = -(temp12*temp_ad19) - temp9*temp_ad18 - temp14*&
&             temp_ad16
            temp_ad24 = -(temp9*(p-plm2)*temp_ad18/temp10)
            temp_ad23 = (dlp0+dlm1+dlm2)*temp_ad24
            temp_ad21 = temp11*temp_ad24
            dlm2_ad = (dlp0+dlm1)*dlm1*temp_ad20 + temp_ad21 + temp13*&
&             temp_ad22 + temp_ad13 + (2*dlm2+dlm1)*temp_ad12
            dlm1_ad = (temp15+(dlp0+dlm1)*dlm2)*temp_ad20 + dlp0*&
&             temp_ad23 + temp_ad21 + (temp13+(dlm1+dlm2)*dlp0)*&
&             temp_ad22 + temp_ad13 + dlm2*temp_ad12
            dlp0_ad = temp15*temp_ad20 + (2*dlp0+dlm1)*temp_ad23 + &
&             temp_ad21 + (dlm1+dlm2)*dlm1*temp_ad22 + temp_ad13
            dlogp1_ad(i, lm2) = dlogp1_ad(i, lm2) + dlm2_ad
            dlogp1_ad(i, lm1) = dlogp1_ad(i, lm1) + dlm1_ad
            dlogp1_ad(i, lp0) = dlogp1_ad(i, lp0) + dlp0_ad
            logpl1_ad(i, lm2) = logpl1_ad(i, lm2) + plm2_ad
            logpl1_ad(i, lm1) = logpl1_ad(i, lm1) + plm1_ad
            logpl1_ad(i, lp0) = logpl1_ad(i, lp0) + plp0_ad
            logpl1_ad(i, lp1) = logpl1_ad(i, lp1) + plp1_ad
            logpl2_ad(i, k) = logpl2_ad(i, k) + p_ad
            CALL POPINTEGER4(lm2)
            CALL POPINTEGER4(lp1)
          ELSE
            temp6 = logpl1(i, lm1) - logpl1(i, lp0)
            temp_ad7 = q2_ad(i, j, k)/temp6
            temp8 = logpl2(i, k) - logpl1(i, lp0)
            temp7 = qx(i, lm1) - qx(i, lp0)
            temp_ad8 = -(temp7*temp8*temp_ad7/temp6)
            qx_ad(i, lp0) = qx_ad(i, lp0) + q2_ad(i, j, k) - temp8*&
&             temp_ad7
            qx_ad(i, lm1) = qx_ad(i, lm1) + temp8*temp_ad7
            logpl2_ad(i, k) = logpl2_ad(i, k) + temp7*temp_ad7
            logpl1_ad(i, lp0) = logpl1_ad(i, lp0) - temp_ad8 - temp7*&
&             temp_ad7
            logpl1_ad(i, lm1) = logpl1_ad(i, lm1) + temp_ad8
            q2_ad(i, j, k) = 0.0_8
          END IF
        ELSE IF (branch .EQ. 2) THEN
          temp3 = logpl1(i, km) - logpl1(i, km-1)
          temp_ad5 = q2_ad(i, j, k)/temp3
          temp5 = logpl2(i, k) - logpl1(i, km)
          temp4 = qx(i, km) - qx(i, km-1)
          temp_ad6 = -(temp4*temp5*temp_ad5/temp3)
          qx_ad(i, km) = qx_ad(i, km) + temp5*temp_ad5 + q2_ad(i, j, k)
          qx_ad(i, km-1) = qx_ad(i, km-1) - temp5*temp_ad5
          logpl2_ad(i, k) = logpl2_ad(i, k) + temp4*temp_ad5
          logpl1_ad(i, km) = logpl1_ad(i, km) + temp_ad6 - temp4*&
&           temp_ad5
          logpl1_ad(i, km-1) = logpl1_ad(i, km-1) - temp_ad6
          q2_ad(i, j, k) = 0.0_8
        ELSE
          temp1 = logpl1(i, 2) - logpl1(i, 1)
          temp2 = qx(i, 2) - qx(i, 1)
          temp0 = temp2/temp1
          temp_ad3 = (logpl2(i, k)-logpl1(i, 1))*q2_ad(i, j, k)/temp1
          temp_ad4 = -(temp0*temp_ad3)
          qx_ad(i, 1) = qx_ad(i, 1) + q2_ad(i, j, k) - temp_ad3
          logpl2_ad(i, k) = logpl2_ad(i, k) + temp0*q2_ad(i, j, k)
          logpl1_ad(i, 1) = logpl1_ad(i, 1) - temp_ad4 - temp0*q2_ad(i, &
&           j, k)
          qx_ad(i, 2) = qx_ad(i, 2) + temp_ad3
          logpl1_ad(i, 2) = logpl1_ad(i, 2) + temp_ad4
          q2_ad(i, j, k) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPINTEGER4(lm1)
        ELSE
          CALL POPINTEGER4(lm1)
        END IF
        CALL POPINTEGER4(ad_count)
        DO i0=1,ad_count
          IF (i0 .EQ. 1) CALL POPCONTROL1B(branch)
        END DO
        CALL POPINTEGER4(lp0)
      END DO
    END DO
    DO i=i2,i1,-1
      temp = pe1(i, km+1) - pe1(i, 1)
      temp_ad2 = -(vsum1(i)*vsum1_ad(i)/temp**2)
      pe1_ad(i, km+1) = pe1_ad(i, km+1) + temp_ad2
      pe1_ad(i, 1) = pe1_ad(i, 1) - temp_ad2
      vsum1_ad(i) = vsum1_ad(i)/temp
      DO k=km,1,-1
        temp_ad1 = qx(i, k)*vsum1_ad(i)
        qx_ad(i, k) = qx_ad(i, k) + (pe1(i, k+1)-pe1(i, k))*vsum1_ad(i)
        pe1_ad(i, k+1) = pe1_ad(i, k+1) + temp_ad1
        pe1_ad(i, k) = pe1_ad(i, k) - temp_ad1
      END DO
    END DO
    DO k=km-1,1,-1
      logpl1_ad(:, k+1) = logpl1_ad(:, k+1) + dlogp1_ad(:, k)
      logpl1_ad(:, k) = logpl1_ad(:, k) - dlogp1_ad(:, k)
      dlogp1_ad(:, k) = 0.0_8
    END DO
    DO k=km,1,-1
      temp_ad0 = logpl2_ad(:, k)/(pe2(:, k)+pe2(:, k+1))
      pe2_ad(:, k) = pe2_ad(:, k) + temp_ad0
      pe2_ad(:, k+1) = pe2_ad(:, k+1) + temp_ad0
      logpl2_ad(:, k) = 0.0_8
    END DO
    DO k=km,1,-1
      temp_ad = logpl1_ad(:, k)/(pe1(:, k)+pe1(:, k+1))
      pe1_ad(:, k) = pe1_ad(:, k) + temp_ad
      pe1_ad(:, k+1) = pe1_ad(:, k+1) + temp_ad
      logpl1_ad(:, k) = 0.0_8
      q2_ad(:, j, k) = q2_ad(:, j, k) + qx_ad(:, k)
      qx_ad(:, k) = 0.0_8
    END DO
  END SUBROUTINE MAP1_CUBIC_TE_ADM
!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  map1_cubic_te --- Cubic Interpolation for TE mapping
!
! !INTERFACE:
  SUBROUTINE MAP1_CUBIC_TE(km, pe1, pe2, q2, i1, i2, j, ibeg, iend, jbeg&
&   , jend, iv, kord)
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
! (from model top to bottom surface)
! in the original vertical coordinate
! pressure at layer edges 
    REAL(p_precision), INTENT(IN) :: pe2(i1:i2, km+1)
! (from model top to bottom surface)
! in the new vertical coordinate
!      real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
! Field output
    REAL, INTENT(INOUT) :: q2(ibeg:iend, jbeg:jend, km)
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
    REAL(p_precision) :: logpl1(i1:i2, km)
    REAL(p_precision) :: logpl2(i1:i2, km)
    REAL(p_precision) :: dlogp1(i1:i2, km)
    REAL(p_precision) :: vsum1(i1:i2)
    REAL(p_precision) :: vsum2(i1:i2)
    REAL(p_precision) :: am2, am1, ap0, ap1, p, plp1, plp0, plm1, plm2, &
&   dlp0, dlm1, dlm2
    INTEGER :: i, k, lm2, lm1, lp0, lp1
    INTRINSIC LOG
    INTRINSIC MAX
    INTRINSIC MIN
! Initialization
! --------------
    DO k=1,km
      qx(:, k) = q2(:, j, k)
      logpl1(:, k) = LOG(d0_5*(pe1(:, k)+pe1(:, k+1)))
    END DO
    DO k=1,km
      logpl2(:, k) = LOG(d0_5*(pe2(:, k)+pe2(:, k+1)))
    END DO
    DO k=1,km-1
      dlogp1(:, k) = logpl1(:, k+1) - logpl1(:, k)
    END DO
! Compute vertical integral of Input TE
! -------------------------------------
    vsum1(:) = d0_0
    DO i=i1,i2
      DO k=1,km
        vsum1(i) = vsum1(i) + qx(i, k)*(pe1(i, k+1)-pe1(i, k))
      END DO
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
          q2(i, j, k) = qx(i, 1) + (qx(i, 2)-qx(i, 1))*(logpl2(i, k)-&
&           logpl1(i, 1))/(logpl1(i, 2)-logpl1(i, 1))
! Extrapolate Linearly in LogP below last model level
! ---------------------------------------------------
        ELSE IF (lm1 .EQ. km .AND. lp0 .EQ. km) THEN
          q2(i, j, k) = qx(i, km) + (qx(i, km)-qx(i, km-1))*(logpl2(i, k&
&           )-logpl1(i, km))/(logpl1(i, km)-logpl1(i, km-1))
! Interpolate Linearly in LogP between levels 1 => 2 and km-1 => km
! -----------------------------------------------------------------
        ELSE IF (lm1 .EQ. 1 .OR. lp0 .EQ. km) THEN
          q2(i, j, k) = qx(i, lp0) + (qx(i, lm1)-qx(i, lp0))*(logpl2(i, &
&           k)-logpl1(i, lp0))/(logpl1(i, lm1)-logpl1(i, lp0))
! Interpolate Cubicly in LogP between other model levels
! ------------------------------------------------------
        ELSE
          lp1 = lp0 + 1
          lm2 = lm1 - 1
          p = logpl2(i, k)
          plp1 = logpl1(i, lp1)
          plp0 = logpl1(i, lp0)
          plm1 = logpl1(i, lm1)
          plm2 = logpl1(i, lm2)
          dlp0 = dlogp1(i, lp0)
          dlm1 = dlogp1(i, lm1)
          dlm2 = dlogp1(i, lm2)
          ap1 = (p-plp0)*(p-plm1)*(p-plm2)/(dlp0*(dlp0+dlm1)*(dlp0+dlm1+&
&           dlm2))
          ap0 = (plp1-p)*(p-plm1)*(p-plm2)/(dlp0*dlm1*(dlm1+dlm2))
          am1 = (plp1-p)*(plp0-p)*(p-plm2)/(dlm1*dlm2*(dlp0+dlm1))
          am2 = (plp1-p)*(plp0-p)*(plm1-p)/(dlm2*(dlm1+dlm2)*(dlp0+dlm1+&
&           dlm2))
          q2(i, j, k) = ap1*qx(i, lp1) + ap0*qx(i, lp0) + am1*qx(i, lm1)&
&           + am2*qx(i, lm2)
        END IF
      END DO
    END DO
! Compute vertical integral of Output TE
! --------------------------------------
    vsum2(:) = d0_0
    DO i=i1,i2
      DO k=1,km
        vsum2(i) = vsum2(i) + q2(i, j, k)*(pe2(i, k+1)-pe2(i, k))
      END DO
      vsum2(i) = vsum2(i)/(pe2(i, km+1)-pe2(i, 1))
    END DO
! Adjust Final TE to conserve
! ---------------------------
    DO i=i1,i2
      DO k=1,km
        q2(i, j, k) = q2(i, j, k) + vsum1(i) - vsum2(i)
      END DO
    END DO
!        q2(i,j,k) = q2(i,j,k) * vsum1(i)/vsum2(i)
    RETURN
  END SUBROUTINE MAP1_CUBIC_TE
!  Differentiation of cs_profile_nolimiters in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: delp a4
!   with respect to varying inputs: delp a4
  SUBROUTINE CS_PROFILE_NOLIMITERS_ADM(a4, a4_ad, delp, delp_ad, km, i1&
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
    REAL :: delp_ad(i1:i2, km)
! Interpolated values
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
    REAL, INTENT(INOUT) :: a4_ad(4, i1:i2, km)
!-----------------------------------------------------------------------
    REAL :: gam(i1:i2, km)
    REAL :: gam_ad(i1:i2, km)
    REAL :: q(i1:i2, km+1)
    REAL :: q_ad(i1:i2, km+1)
    REAL :: d4(i1:i2)
    REAL :: d4_ad(i1:i2)
    REAL :: bet, a_bot, grat
    REAL :: bet_ad, a_bot_ad, grat_ad
    INTEGER :: i, k
    REAL :: temp1
    REAL :: temp0
    REAL :: temp_ad9
    REAL :: temp_ad8
    REAL :: temp_ad7
    REAL :: temp_ad6
    REAL :: temp_ad5
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    REAL :: temp
    DO i=i1,i2
! grid ratio
      grat = delp(i, 2)/delp(i, 1)
      bet = grat*(grat+0.5)
      q(i, 1) = ((grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))/bet
      gam(i, 1) = (1.+grat*(grat+1.5))/bet
    END DO
    DO k=2,km
      DO i=i1,i2
        CALL PUSHREAL8(d4(i))
        d4(i) = delp(i, k-1)/delp(i, k)
        CALL PUSHREAL8(bet)
        bet = 2. + d4(i) + d4(i) - gam(i, k-1)
        CALL PUSHREAL8(q(i, k))
        q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
        gam(i, k) = d4(i)/bet
      END DO
    END DO
    DO i=i1,i2
      a_bot = 1. + d4(i)*(d4(i)+1.5)
      CALL PUSHREAL8(q(i, km+1))
      q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&       a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
    END DO
    DO k=km,1,-1
      DO i=i1,i2
        CALL PUSHREAL8(q(i, k))
        q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
      END DO
    END DO
    q_ad = 0.0_8
    DO k=km,1,-1
      DO i=i2,i1,-1
        temp_ad9 = 3.*a4_ad(4, i, k)
        a4_ad(1, i, k) = a4_ad(1, i, k) + 2.*temp_ad9
        a4_ad(2, i, k) = a4_ad(2, i, k) - temp_ad9
        a4_ad(3, i, k) = a4_ad(3, i, k) - temp_ad9
        a4_ad(4, i, k) = 0.0_8
        q_ad(i, k+1) = q_ad(i, k+1) + a4_ad(3, i, k)
        a4_ad(3, i, k) = 0.0_8
        q_ad(i, k) = q_ad(i, k) + a4_ad(2, i, k)
        a4_ad(2, i, k) = 0.0_8
      END DO
    END DO
    gam_ad = 0.0_8
    DO k=1,km,1
      DO i=i2,i1,-1
        CALL POPREAL8(q(i, k))
        gam_ad(i, k) = gam_ad(i, k) - q(i, k+1)*q_ad(i, k)
        q_ad(i, k+1) = q_ad(i, k+1) - gam(i, k)*q_ad(i, k)
      END DO
    END DO
    d4_ad = 0.0_8
    DO i=i2,i1,-1
      a_bot = 1. + d4(i)*(d4(i)+1.5)
      CALL POPREAL8(q(i, km+1))
      temp1 = d4(i)*(d4(i)+0.5) - a_bot*gam(i, km)
      temp_ad6 = q_ad(i, km+1)/temp1
      temp0 = d4(i)*(d4(i)+1.)
      temp_ad7 = 2.*a4(1, i, km)*temp_ad6
      temp_ad8 = -((2.*(temp0*a4(1, i, km))+a4(1, i, km-1)-a_bot*q(i, km&
&       ))*temp_ad6/temp1)
      a4_ad(1, i, km) = a4_ad(1, i, km) + 2.*temp0*temp_ad6
      a4_ad(1, i, km-1) = a4_ad(1, i, km-1) + temp_ad6
      a_bot_ad = -(gam(i, km)*temp_ad8) - q(i, km)*temp_ad6
      d4_ad(i) = d4_ad(i) + (2*d4(i)+1.5)*a_bot_ad + (2*d4(i)+0.5)*&
&       temp_ad8 + (2*d4(i)+1.)*temp_ad7
      q_ad(i, km) = q_ad(i, km) - a_bot*temp_ad6
      gam_ad(i, km) = gam_ad(i, km) - a_bot*temp_ad8
      q_ad(i, km+1) = 0.0_8
    END DO
    DO k=km,2,-1
      DO i=i2,i1,-1
        temp_ad4 = q_ad(i, k)/bet
        temp_ad3 = 3.*temp_ad4
        CALL POPREAL8(q(i, k))
        bet_ad = -((3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))*&
&         temp_ad4/bet) - d4(i)*gam_ad(i, k)/bet**2
        d4_ad(i) = d4_ad(i) + a4(1, i, k)*temp_ad3 + 2*bet_ad + gam_ad(i&
&         , k)/bet
        gam_ad(i, k) = 0.0_8
        a4_ad(1, i, k-1) = a4_ad(1, i, k-1) + temp_ad3
        a4_ad(1, i, k) = a4_ad(1, i, k) + d4(i)*temp_ad3
        q_ad(i, k-1) = q_ad(i, k-1) - temp_ad4
        q_ad(i, k) = 0.0_8
        CALL POPREAL8(bet)
        gam_ad(i, k-1) = gam_ad(i, k-1) - bet_ad
        CALL POPREAL8(d4(i))
        temp_ad5 = d4_ad(i)/delp(i, k)
        delp_ad(i, k-1) = delp_ad(i, k-1) + temp_ad5
        delp_ad(i, k) = delp_ad(i, k) - delp(i, k-1)*temp_ad5/delp(i, k)
        d4_ad(i) = 0.0_8
      END DO
    END DO
    DO i=i2,i1,-1
      grat = delp(i, 2)/delp(i, 1)
      bet = grat*(grat+0.5)
      temp_ad = gam_ad(i, 1)/bet
      gam_ad(i, 1) = 0.0_8
      temp_ad1 = q_ad(i, 1)/bet
      temp_ad0 = a4(1, i, 1)*temp_ad1
      temp = 2*grat*(grat+1.)
      bet_ad = -((temp*a4(1, i, 1)+a4(1, i, 2))*temp_ad1/bet) - (grat*(&
&       grat+1.5)+1.)*temp_ad/bet
      grat_ad = (4*grat+2*1.)*temp_ad0 + (2*grat+0.5)*bet_ad + (2*grat+&
&       1.5)*temp_ad
      a4_ad(1, i, 1) = a4_ad(1, i, 1) + temp*temp_ad1
      a4_ad(1, i, 2) = a4_ad(1, i, 2) + temp_ad1
      q_ad(i, 1) = 0.0_8
      temp_ad2 = grat_ad/delp(i, 1)
      delp_ad(i, 2) = delp_ad(i, 2) + temp_ad2
      delp_ad(i, 1) = delp_ad(i, 1) - delp(i, 2)*temp_ad2/delp(i, 1)
    END DO
  END SUBROUTINE CS_PROFILE_NOLIMITERS_ADM
  SUBROUTINE CS_PROFILE_NOLIMITERS(a4, delp, km, i1, i2, iv, kord)
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
! Interpolated values
    REAL, INTENT(INOUT) :: a4(4, i1:i2, km)
!-----------------------------------------------------------------------
    REAL :: gam(i1:i2, km)
    REAL :: q(i1:i2, km+1)
    REAL :: d4(i1:i2)
    REAL :: bet, a_bot, grat
    INTEGER :: i, k
    DO i=i1,i2
! grid ratio
      grat = delp(i, 2)/delp(i, 1)
      bet = grat*(grat+0.5)
      q(i, 1) = ((grat+grat)*(grat+1.)*a4(1, i, 1)+a4(1, i, 2))/bet
      gam(i, 1) = (1.+grat*(grat+1.5))/bet
    END DO
    DO k=2,km
      DO i=i1,i2
        d4(i) = delp(i, k-1)/delp(i, k)
        bet = 2. + d4(i) + d4(i) - gam(i, k-1)
        q(i, k) = (3.*(a4(1, i, k-1)+d4(i)*a4(1, i, k))-q(i, k-1))/bet
        gam(i, k) = d4(i)/bet
      END DO
    END DO
    DO i=i1,i2
      a_bot = 1. + d4(i)*(d4(i)+1.5)
      q(i, km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1, i, km)+a4(1, i, km-1)-&
&       a_bot*q(i, km))/(d4(i)*(d4(i)+0.5)-a_bot*gam(i, km))
    END DO
    DO k=km,1,-1
      DO i=i1,i2
        q(i, k) = q(i, k) - gam(i, k)*q(i, k+1)
      END DO
    END DO
    DO k=1,km
      DO i=i1,i2
        a4(2, i, k) = q(i, k)
        a4(3, i, k) = q(i, k+1)
        a4(4, i, k) = 3.*(2.*a4(1, i, k)-(a4(2, i, k)+a4(3, i, k)))
      END DO
    END DO
  END SUBROUTINE CS_PROFILE_NOLIMITERS

end module fv_mapz_adm_mod
