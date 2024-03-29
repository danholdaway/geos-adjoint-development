!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5903) - 14 Dec 2015 10:32
!
MODULE GEOS_MOISTGRIDCOMP_B
  USE QSAT_UTIL_B
  USE RAS_B
  USE CLOUDNEW_B
  USE RASPARAMS
  USE CLDPARAMS
  USE MAPL_CONSTANTSMOD
  IMPLICIT NONE
!      call RASE_FAST( IDIM                 , &
!                 IRUN                 , &                 
!                 LM                   , &
!                 ICMIN                , &
!                 DT_MOIST             , &
!                 MAPL8_CP              , &
!                 MAPL8_ALHL            , &
!                 MAPL8_ALHS            , &
!                 MAPL8_TICE            , &
!                 MAPL8_GRAV            , &
!                 SEEDRAS              , &
!                 SIGE                 , &
!                 KCBL                 , &
!                 WGT0                 , &
!                 WGT1                 , &
!                 TPERT                , &
!                 QPERT                , &
!                 TH1                  , &
!                 Q1                   , & 
!                 QSS                  , & 
!                 DQS                  , &
!                 CNV_FRACTION         , &
!                 RASAL2_2d            , &
!                 CO_AUTO              , &
!                 CNV_PLE              , &
!                 PKE                  , &
!                 RASPARAMS              )
  INCLUDE 'preamble'
  PRIVATE 
  PUBLIC moist_run
  PUBLIC moist_run_adm

CONTAINS
!  Differentiation of moist_run in reverse (adjoint) mode (with options r8 split(GEOS_MoistGridComp.PRE_PROGNO_CLOUD)):
!   gradient     of useful results: v1 qlls th1 qlcn qils q1 clls
!                qicn clcn u1
!   with respect to varying inputs: v1 qlls th1 qlcn qils q1 clls
!                qicn clcn u1
!   RW status of diff variables: v1:in-out qlls:in-out th1:in-out
!                qlcn:in-out qils:in-out q1:in-out clls:in-out
!                qicn:in-out clcn:in-out u1:in-out
  SUBROUTINE MOIST_RUN_ADM(im, jm, lm, dt_moist, idim, irun, icmin, &
&   doconvec, cnv_ple, pke, plo, pk, sige, iras, jras, wgt0, wgt1, &
&   co_auto, rasparams, cldparams, sclmfdfr, frland, kcbl, ts, kh, &
&   seedras, rasal2_2d, cnv_fraction, itrcr, irccode, cbl_tpert, &
&   cbl_qpert, cbl_tpert_mxocn, cbl_tpert_mxlnd, th1, th1_ad, q1, q1_ad&
&   , u1, u1_ad, v1, v1_ad, ple, qlls, qlls_ad, qlcn, qlcn_ad, qils, &
&   qils_ad, qicn, qicn_ad, clcn, clcn_ad, clls, clls_ad)
    IMPLICIT NONE
!Inputs - not linearized
    INTEGER, INTENT(IN) :: im, jm, lm, idim, irun, icmin, itrcr
    REAL*8, INTENT(IN) :: dt_moist, sige(0:lm)
    INTEGER, DIMENSION(im, jm), INTENT(IN) :: iras, jras
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: wgt0, wgt1
    INTEGER, DIMENSION(im, jm), INTENT(IN) :: kcbl, seedras, doconvec
    REAL*8, DIMENSION(im, jm), INTENT(IN) :: ts, co_auto, frland
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: ple, kh
    TYPE(RASPARAM_TYPE), INTENT(IN) :: rasparams
    TYPE(CLDPARAM_TYPE), INTENT(IN) :: cldparams
    REAL*8, INTENT(IN) :: sclmfdfr, cbl_tpert, cbl_qpert, &
&   cbl_tpert_mxocn, cbl_tpert_mxlnd
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: cnv_ple, pke
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: plo, pk
!Inputs - could linearze but not be sensible, e.g. because stochasitc
    REAL*8, DIMENSION(im, jm), INTENT(IN) :: rasal2_2d, cnv_fraction
!Inouts
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: th1, q1, u1, v1
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: th1_ad
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: qlls, qlcn, qils, &
&   qicn, clcn, clls
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: qlls_ad
!Outs
    INTEGER, INTENT(OUT) :: irccode
    REAL*8, DIMENSION(im, jm, lm) :: dqs, qss
    REAL*8, DIMENSION(im, jm, lm) :: dqs_ad, qss_ad
    REAL*8, DIMENSION(im, jm, lm) :: qst3
    REAL*8, DIMENSION(im, jm, lm) :: qst3_ad
    REAL*8, DIMENSION(im, jm, lm) :: dzet, qddf3
    REAL*8, DIMENSION(im, jm, lm) :: dzet_ad, qddf3_ad
    REAL*8, DIMENSION(im, jm) :: tpert, qpert
    REAL*8, DIMENSION(im, jm) :: tpert_ad
    REAL*8, DIMENSION(im, jm) :: mxdiamx
    REAL*8, DIMENSION(im, jm, lm) :: cnv_dqldt, cnv_mfd, cnv_prc3, &
&   cnv_updf, cnv_cvw, entlam
    REAL*8, DIMENSION(im, jm, lm) :: cnv_dqldt_ad, cnv_mfd_ad, &
&   cnv_prc3_ad, cnv_updf_ad
    INTEGER :: i, j, l, k
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: clcn_ad
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: u1_ad
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: q1_ad
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: v1_ad
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: qils_ad
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: clls_ad
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: qlcn_ad
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: qicn_ad
    CALL PRE_RASE(im, jm, lm, th1, q1, ts, cnv_fraction, frland, plo, &
&           pke, pk, cbl_tpert, cbl_qpert, cbl_tpert_mxocn, &
&           cbl_tpert_mxlnd, tpert, qpert, dqs, qss)
    CALL PUSHREAL8ARRAY(v1, im*jm*lm)
    CALL PUSHREAL8ARRAY(u1, im*jm*lm)
    CALL PUSHREAL8ARRAY(q1, im*jm*lm)
    CALL PUSHREAL8ARRAY(th1, im*jm*lm)
    CALL RASE(idim, irun, lm, icmin, dt_moist, doconvec, mapl8_cp, &
&       mapl8_alhl, mapl8_alhs, mapl8_tice, mapl8_grav, seedras, sige, &
&       kcbl, wgt0, wgt1, tpert, qpert, th1, q1, u1, v1, qss, dqs, &
&       cnv_fraction, rasal2_2d, co_auto, cnv_ple, pke, cnv_dqldt, &
&       cnv_mfd, cnv_prc3, cnv_updf, rasparams, itrcr)
    CALL PRE_PROGNO_CLOUD_FWD(im, jm, lm, th1, pk, plo, pke, cnv_ple, &
&                       qst3, dzet, qddf3, cnv_fraction, cldparams)
    CALL PROGNO_CLOUD_ADM(idim, lm, dt_moist, plo, cnv_ple, pk, frland, &
&                   kh, cnv_mfd, cnv_mfd_ad, cnv_dqldt, cnv_dqldt_ad, &
&                   cnv_prc3, cnv_prc3_ad, cnv_updf, cnv_updf_ad, u1, &
&                   u1_ad, v1, v1_ad, th1, th1_ad, q1, q1_ad, qlls, &
&                   qlls_ad, qlcn, qlcn_ad, qils, qils_ad, qicn, qicn_ad&
&                   , clcn, clcn_ad, clls, clls_ad, cldparams, sclmfdfr&
&                   , qst3, qst3_ad, dzet, dzet_ad, qddf3, qddf3_ad, &
&                   cnv_fraction)
    CALL PRE_PROGNO_CLOUD_BWD(im, jm, lm, th1, th1_ad, pk, plo, pke, &
&                       cnv_ple, qst3, qst3_ad, dzet, dzet_ad, qddf3, &
&                       qddf3_ad, cnv_fraction, cldparams)
    CALL POPREAL8ARRAY(th1, im*jm*lm)
    CALL POPREAL8ARRAY(q1, im*jm*lm)
    CALL POPREAL8ARRAY(u1, im*jm*lm)
    CALL POPREAL8ARRAY(v1, im*jm*lm)
    CALL RASE_ADM(idim, irun, lm, icmin, dt_moist, doconvec, mapl8_cp, &
&           mapl8_alhl, mapl8_alhs, mapl8_tice, mapl8_grav, seedras, &
&           sige, kcbl, wgt0, wgt1, tpert, tpert_ad, qpert, th1, th1_ad&
&           , q1, q1_ad, u1, u1_ad, v1, v1_ad, qss, qss_ad, dqs, dqs_ad&
&           , cnv_fraction, rasal2_2d, co_auto, cnv_ple, pke, cnv_dqldt&
&           , cnv_dqldt_ad, cnv_mfd, cnv_mfd_ad, cnv_prc3, cnv_prc3_ad, &
&           cnv_updf, cnv_updf_ad, rasparams, itrcr)
    CALL PRE_RASE_ADM(im, jm, lm, th1, th1_ad, q1, q1_ad, ts, &
&               cnv_fraction, frland, plo, pke, pk, cbl_tpert, cbl_qpert&
&               , cbl_tpert_mxocn, cbl_tpert_mxlnd, tpert, tpert_ad, &
&               qpert, dqs, dqs_ad, qss, qss_ad)
  END SUBROUTINE MOIST_RUN_ADM
  SUBROUTINE MOIST_RUN(im, jm, lm, dt_moist, idim, irun, icmin, doconvec&
&   , cnv_ple, pke, plo, pk, sige, iras, jras, wgt0, wgt1, co_auto, &
&   rasparams, cldparams, sclmfdfr, frland, kcbl, ts, kh, seedras, &
&   rasal2_2d, cnv_fraction, itrcr, irccode, cbl_tpert, cbl_qpert, &
&   cbl_tpert_mxocn, cbl_tpert_mxlnd, th1, q1, u1, v1, ple, qlls, qlcn, &
&   qils, qicn, clcn, clls)
    IMPLICIT NONE
!Inputs - not linearized
    INTEGER, INTENT(IN) :: im, jm, lm, idim, irun, icmin, itrcr
    REAL*8, INTENT(IN) :: dt_moist, sige(0:lm)
    INTEGER, DIMENSION(im, jm), INTENT(IN) :: iras, jras
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: wgt0, wgt1
    INTEGER, DIMENSION(im, jm), INTENT(IN) :: kcbl, seedras, doconvec
    REAL*8, DIMENSION(im, jm), INTENT(IN) :: ts, co_auto, frland
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: ple, kh
    TYPE(RASPARAM_TYPE), INTENT(IN) :: rasparams
    TYPE(CLDPARAM_TYPE), INTENT(IN) :: cldparams
    REAL*8, INTENT(IN) :: sclmfdfr, cbl_tpert, cbl_qpert, &
&   cbl_tpert_mxocn, cbl_tpert_mxlnd
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: cnv_ple, pke
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: plo, pk
!Inputs - could linearze but not be sensible, e.g. because stochasitc
    REAL*8, DIMENSION(im, jm), INTENT(IN) :: rasal2_2d, cnv_fraction
!Inouts
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: th1, q1, u1, v1
    REAL*8, DIMENSION(im, jm, lm), INTENT(INOUT) :: qlls, qlcn, qils, &
&   qicn, clcn, clls
!Outs
    INTEGER, INTENT(OUT) :: irccode
    REAL*8, DIMENSION(im, jm, lm) :: dqs, qss
    REAL*8, DIMENSION(im, jm, lm) :: qst3
    REAL*8, DIMENSION(im, jm, lm) :: dzet, qddf3
    REAL*8, DIMENSION(im, jm) :: tpert, qpert
    REAL*8, DIMENSION(im, jm) :: mxdiamx
    REAL*8, DIMENSION(im, jm, lm) :: cnv_dqldt, cnv_mfd, cnv_prc3, &
&   cnv_updf, cnv_cvw, entlam
    INTEGER :: i, j, l, k
    CALL PRE_RASE(im, jm, lm, th1, q1, ts, cnv_fraction, frland, plo, &
&           pke, pk, cbl_tpert, cbl_qpert, cbl_tpert_mxocn, &
&           cbl_tpert_mxlnd, tpert, qpert, dqs, qss)
    CALL RASE(idim, irun, lm, icmin, dt_moist, doconvec, mapl8_cp, &
&       mapl8_alhl, mapl8_alhs, mapl8_tice, mapl8_grav, seedras, sige, &
&       kcbl, wgt0, wgt1, tpert, qpert, th1, q1, u1, v1, qss, dqs, &
&       cnv_fraction, rasal2_2d, co_auto, cnv_ple, pke, cnv_dqldt, &
&       cnv_mfd, cnv_prc3, cnv_updf, rasparams, itrcr)
    CALL PRE_PROGNO_CLOUD(im, jm, lm, th1, pk, plo, pke, cnv_ple, qst3, &
&                   dzet, qddf3, cnv_fraction, cldparams)
    CALL PROGNO_CLOUD(idim, lm, dt_moist, plo, cnv_ple, pk, frland, kh, &
&               cnv_mfd, cnv_dqldt, cnv_prc3, cnv_updf, u1, v1, th1, q1&
&               , qlls, qlcn, qils, qicn, clcn, clls, cldparams, &
&               sclmfdfr, qst3, dzet, qddf3, cnv_fraction)
  END SUBROUTINE MOIST_RUN
!  Differentiation of pre_rase in reverse (adjoint) mode (with options r8 split(GEOS_MoistGridComp.PRE_PROGNO_CLOUD)):
!   gradient     of useful results: dqs tpert th1 q1 qss
!   with respect to varying inputs: th1 q1
  SUBROUTINE PRE_RASE_ADM(im, jm, lm, th1, th1_ad, q1, q1_ad, ts, &
&   cnv_fraction, frland, plo, pke, pk, cbl_tpert, cbl_qpert, &
&   cbl_tpert_mxocn, cbl_tpert_mxlnd, tpert, tpert_ad, qpert, dqs, &
&   dqs_ad, qss, qss_ad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im, jm, lm
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: th1, q1
    REAL*8, DIMENSION(im, jm, lm) :: th1_ad, q1_ad
    REAL*8, DIMENSION(im, jm), INTENT(IN) :: ts, cnv_fraction, frland
    REAL*8, INTENT(IN) :: cbl_tpert, cbl_qpert, cbl_tpert_mxocn, &
&   cbl_tpert_mxlnd
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: pke
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: plo, pk
    REAL*8, DIMENSION(im, jm, lm) :: dqs, qss
    REAL*8, DIMENSION(im, jm, lm) :: dqs_ad, qss_ad
    REAL*8, DIMENSION(im, jm) :: tpert, qpert
    REAL*8, DIMENSION(im, jm) :: tpert_ad
    INTEGER :: i, j, l
    REAL*8, DIMENSION(im, jm, 0:lm) :: zle
    REAL*8, DIMENSION(im, jm, 0:lm) :: zle_ad
    REAL*8, DIMENSION(im, jm, lm) :: temp1, zlo
    REAL*8, DIMENSION(im, jm, lm) :: temp1_ad, zlo_ad
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: branch
    LOGICAL :: mask3(im, jm)
    LOGICAL :: mask2(im, jm)
    LOGICAL :: mask1(im, jm)
    LOGICAL :: mask0(im, jm)
    REAL*8 :: abs0
    LOGICAL :: mask(im, jm)
    temp1 = th1*pk
    zle(:, :, lm) = 0.
    DO l=lm,1,-1
      zle(:, :, l-1) = th1(:, :, l)*(1.+mapl8_vireps*q1(:, :, l))
      zlo(:, :, l) = zle(:, :, l) + mapl8_cp/mapl8_grav*(pke(:, :, l)-pk&
&       (:, :, l))*zle(:, :, l-1)
      zle(:, :, l-1) = zlo(:, :, l) + mapl8_cp/mapl8_grav*(pk(:, :, l)-&
&       pke(:, :, l-1))*zle(:, :, l-1)
    END DO
    IF (cbl_tpert .GE. 0.) THEN
      abs0 = cbl_tpert
    ELSE
      abs0 = -cbl_tpert
    END IF
    tpert = abs0*(ts-(temp1(:, :, lm)+mapl8_grav*zlo(:, :, lm)/mapl8_cp)&
&     )
    IF (cbl_tpert .LT. 0) THEN
! Make TPERT 0 in areas of deep convection
      tpert = tpert*(1.0-cnv_fraction)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
!CBL_QPERT * ( QSSFC - Q(:,:,LM) )      !dh: CBL_QPERT = 0.0
    mask = tpert .LT. 0.0
    WHERE (mask) 
      tpert = 0.0
    ELSEWHERE
      tpert = tpert
    END WHERE
    mask0 = qpert .LT. 0.0
    mask1 = frland .LT. 0.1
    mask2 = tpert .GT. cbl_tpert_mxocn
    WHERE (mask2) 
      tpert = cbl_tpert_mxocn
    ELSEWHERE
      tpert = tpert
    END WHERE
    mask3 = tpert .GT. cbl_tpert_mxlnd
    WHERE (mask3) tpert_ad = 0.0_8
    WHERE (mask2) tpert_ad = 0.0_8
    WHERE (mask) tpert_ad = 0.0_8
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) tpert_ad = (1.0-cnv_fraction)*tpert_ad
    zlo_ad = 0.0_8
    temp1_ad = 0.0_8
    temp1_ad(:, :, lm) = temp1_ad(:, :, lm) - abs0*tpert_ad
    zlo_ad(:, :, lm) = zlo_ad(:, :, lm) - mapl8_grav*abs0*tpert_ad/&
&     mapl8_cp
    DO j=jm,1,-1
      DO i=im,1,-1
        CALL DQSATPERT_ADM(dqs(i, j, :), dqs_ad(i, j, :), qss(i, j, :), &
&                    qss_ad(i, j, :), temp1(i, j, :), temp1_ad(i, j, :)&
&                    , plo(i, j, :), lm)
      END DO
    END DO
    zle_ad = 0.0_8
    DO l=1,lm,1
      zlo_ad(:, :, l) = zlo_ad(:, :, l) + zle_ad(:, :, l-1)
      zle_ad(:, :, l-1) = mapl8_cp*(pk(:, :, l)-pke(:, :, l-1))*zle_ad(:&
&       , :, l-1)/mapl8_grav
      zle_ad(:, :, l) = zle_ad(:, :, l) + zlo_ad(:, :, l)
      zle_ad(:, :, l-1) = zle_ad(:, :, l-1) + mapl8_cp*(pke(:, :, l)-pk(&
&       :, :, l))*zlo_ad(:, :, l)/mapl8_grav
      zlo_ad(:, :, l) = 0.0_8
      th1_ad(:, :, l) = th1_ad(:, :, l) + (mapl8_vireps*q1(:, :, l)+1.)*&
&       zle_ad(:, :, l-1)
      q1_ad(:, :, l) = q1_ad(:, :, l) + mapl8_vireps*th1(:, :, l)*zle_ad&
&       (:, :, l-1)
      zle_ad(:, :, l-1) = 0.0_8
    END DO
    th1_ad = th1_ad + pk*temp1_ad
  END SUBROUTINE PRE_RASE_ADM
  SUBROUTINE PRE_RASE(im, jm, lm, th1, q1, ts, cnv_fraction, frland, plo&
&   , pke, pk, cbl_tpert, cbl_qpert, cbl_tpert_mxocn, cbl_tpert_mxlnd, &
&   tpert, qpert, dqs, qss)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im, jm, lm
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: th1, q1
    REAL*8, DIMENSION(im, jm), INTENT(IN) :: ts, cnv_fraction, frland
    REAL*8, INTENT(IN) :: cbl_tpert, cbl_qpert, cbl_tpert_mxocn, &
&   cbl_tpert_mxlnd
    REAL*8, DIMENSION(im, jm, 0:lm), INTENT(IN) :: pke
    REAL*8, DIMENSION(im, jm, lm), INTENT(IN) :: plo, pk
    REAL*8, DIMENSION(im, jm, lm), INTENT(OUT) :: dqs, qss
    REAL*8, DIMENSION(im, jm), INTENT(OUT) :: tpert, qpert
    INTEGER :: i, j, l
    REAL*8, DIMENSION(im, jm, 0:lm) :: zle
    REAL*8, DIMENSION(im, jm, lm) :: temp1, zlo
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    LOGICAL :: mask3(im, jm)
    LOGICAL :: mask2(im, jm)
    LOGICAL :: mask1(im, jm)
    LOGICAL :: mask0(im, jm)
    REAL*8 :: abs0
    LOGICAL :: mask(im, jm)
    temp1 = th1*pk
    zle(:, :, lm) = 0.
    DO l=lm,1,-1
      zle(:, :, l-1) = th1(:, :, l)*(1.+mapl8_vireps*q1(:, :, l))
      zlo(:, :, l) = zle(:, :, l) + mapl8_cp/mapl8_grav*(pke(:, :, l)-pk&
&       (:, :, l))*zle(:, :, l-1)
      zle(:, :, l-1) = zlo(:, :, l) + mapl8_cp/mapl8_grav*(pk(:, :, l)-&
&       pke(:, :, l-1))*zle(:, :, l-1)
    END DO
    DO j=1,jm
      DO i=1,im
        CALL DQSATPERT(dqs(i, j, :), qss(i, j, :), temp1(i, j, :), plo(i&
&                , j, :), lm)
      END DO
    END DO
    IF (cbl_tpert .GE. 0.) THEN
      abs0 = cbl_tpert
    ELSE
      abs0 = -cbl_tpert
    END IF
    tpert = abs0*(ts-(temp1(:, :, lm)+mapl8_grav*zlo(:, :, lm)/mapl8_cp)&
&     )
    IF (cbl_tpert .LT. 0) tpert = tpert*(1.0-cnv_fraction)
! Make TPERT 0 in areas of deep convection
!CBL_QPERT * ( QSSFC - Q(:,:,LM) )      !dh: CBL_QPERT = 0.0
    qpert = 0.0
    mask = tpert .LT. 0.0
    WHERE (mask) 
      tpert = 0.0
    ELSEWHERE
      tpert = tpert
    END WHERE
    mask0 = qpert .LT. 0.0
    WHERE (mask0) 
      qpert = 0.0
    ELSEWHERE
      qpert = qpert
    END WHERE
    mask1 = frland .LT. 0.1
    mask2 = tpert .GT. cbl_tpert_mxocn
    WHERE (mask2) 
      tpert = cbl_tpert_mxocn
    ELSEWHERE
      tpert = tpert
    END WHERE
    mask3 = tpert .GT. cbl_tpert_mxlnd
    WHERE (mask3) 
      tpert = cbl_tpert_mxlnd
    ELSEWHERE
      tpert = tpert
    END WHERE
  END SUBROUTINE PRE_RASE
END MODULE GEOS_MOISTGRIDCOMP_B
