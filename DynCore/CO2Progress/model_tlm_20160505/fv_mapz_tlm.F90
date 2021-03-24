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

!  public compute_total_energy_tlm, Lagrangian_to_Eulerian_tlm
  public Lagrangian_to_Eulerian_tlm

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
!   variations   of useful results: peln q delp pkz
!   with respect to varying inputs: peln q delp pkz pe pk
!#ifdef MAPL_MODE
!#endif
  SUBROUTINE LAGRANGIAN_TO_EULERIAN_TLM(do_consv, consv, pe, pe_tl, delp&
&   , delp_tl, pkz, pkz_tl, pk, pk_tl, pdt, km, is, ie, js, je, isd, ied&
&   , jsd, jed, nq, sphum, u, v, pt, q, q_tl, hs, r_vir, cp, akap, pi, &
&   radius, grav, kord_mt, kord_tr, kord_tm, peln, peln_tl, te0_2d, ng, &
&   ua, te, pem, fill, reproduce_sum, ak, bk, ks, te_method, remap_t, &
&   ncnst)
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
    REAL(p_precision), INTENT(INOUT) :: pk_tl(is:ie, js:je, km+1)
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, km, ncnst)
    REAL, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed, km, ncnst)
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
    REAL(p_precision), INTENT(INOUT) :: peln_tl(is:ie, km+1, js:je)
! layer-mean pk for converting t to pt
    REAL(p_precision), INTENT(OUT) :: pkz(is:ie, js:je, km)
    REAL(p_precision), INTENT(OUT) :: pkz_tl(is:ie, js:je, km)
    REAL, INTENT(OUT) :: te(is:ie, js:je, km)
! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!-----------------------------------------------------------------------
    INTEGER :: i, j, k
    REAL(p_precision) :: te_2d(is:ie, js:je)
    REAL(p_precision) :: zsum0(is:ie, js:je)
    REAL(p_precision) :: zsum1(is:ie, js:je)
    REAL(p_precision) :: q2(is:ie, km)
    REAL(p_precision) :: q2_tl(is:ie, km)
    REAL(p_precision) :: dp2(is:ie, km)
    REAL(p_precision) :: dp2_tl(is:ie, km)
    REAL(p_precision) :: pe1(is:ie, km+1)
    REAL(p_precision) :: pe1_tl(is:ie, km+1)
    REAL(p_precision) :: pe2(is:ie, km+1)
    REAL(p_precision) :: pe2_tl(is:ie, km+1)
    REAL(p_precision) :: pk1(is:ie, km+1)
    REAL(p_precision) :: pk2(is:ie, km+1)
    REAL(p_precision) :: pn1(is:ie, km+1)
    REAL(p_precision) :: pn2(is:ie, km+1)
    REAL(p_precision) :: pe0(is:ie+1, km+1)
    REAL(p_precision) :: pe3(is:ie+1, km+1)
    REAL(p_precision) :: phis(is:ie, km+1)
    REAL(p_precision) :: gz(is:ie)
    REAL(p_precision) :: rg, tmp, tpe
    REAL(p_precision) :: bkh
    REAL(p_precision) :: dtmp
    REAL(p_precision) :: dlnp
    INTEGER :: iq
    LOGICAL :: te_map
    REAL(p_precision) :: tmpsum
    rg = akap*cp
    te_map = .true.
    CALL PKEZ_TLM(km, is, ie, js, je, pe, pk, pk_tl, akap, peln, peln_tl&
&           , pkz, pkz_tl)
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
    pe1_tl = 0.0_8
    pe2_tl = 0.0_8
    dp2_tl = 0.0_8
    q2_tl = 0.0_8
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
        DO k=1,km
          DO i=is,ie
            delp_tl(i, j, k) = dp2_tl(i, k)
            delp(i, j, k) = dp2(i, k)
          END DO
        END DO
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
      END IF
    END DO
  END SUBROUTINE LAGRANGIAN_TO_EULERIAN_TLM
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
