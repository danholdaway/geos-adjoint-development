!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.9 (r5096) - 24 Feb 2014 16:53
!
!  Differentiation of rase0 in forward (tangent) mode:
!   variations   of useful results: tho qho
!   with respect to varying inputs: tho qho
!   RW status of diff variables: tho:in-out qho:in-out
SUBROUTINE RASE0_D(idim, irun, k0, icmin, dt, cons_cp, cons_alhl, &
& cons_grav, cons_rgas, cons_h2omw, cons_airmw, cons_vireps, seedras, &
& sige, kcbl, wgt0, wgt1, frland, ts, tho, thod, qho, qhod, co_auto, ple&
& , rasparams, estblx)
  IMPLICIT NONE
!INPUTS
  INTEGER, INTENT(IN) :: idim, irun, k0, icmin
  REAL*8, DIMENSION(idim, k0 + 1), INTENT(IN) :: ple
  REAL*8, DIMENSION(k0 + 1), INTENT(IN) :: sige
  REAL*8, INTENT(IN) :: dt, cons_cp, cons_alhl, cons_grav, cons_rgas
  REAL*8, INTENT(IN) :: cons_h2omw, cons_airmw, cons_vireps
  INTEGER, DIMENSION(idim), INTENT(IN) :: seedras
  INTEGER, DIMENSION(idim), INTENT(IN) :: kcbl
  REAL*8, DIMENSION(idim), INTENT(IN) :: ts, frland
  REAL*8, DIMENSION(idim), INTENT(IN) :: co_auto
  REAL*8, DIMENSION(idim, k0), INTENT(IN) :: wgt0, wgt1
  REAL*8, DIMENSION(:), INTENT(IN) :: rasparams
  REAL*8, DIMENSION(:), INTENT(IN) :: estblx
!PROGNOSTIC
  REAL*8, DIMENSION(idim, k0), INTENT(INOUT) :: tho, qho
  REAL*8, DIMENSION(idim, k0), INTENT(INOUT) :: thod, qhod
!LOCALS
  INTEGER :: i, ic, l, kk, k
!Parameters
  REAL*8, PARAMETER :: onepkap=1.+2./7., daylen=86400.0
  REAL*8, PARAMETER :: rhmax=0.9999
  REAL*8, PARAMETER :: cbl_qpert=0.0, cbl_tpert=1.0
  REAL*8, PARAMETER :: cbl_tpert_mxocn=2.0, cbl_tpert_mxlnd=4.0
!Constants
  REAL*8 :: grav, cp, alhl, cpbg, alhi, cpi, gravi, ddt, lbcp
!Rasparams
  REAL*8 :: fricfac, cli_crit, rasal1, rasal2
  REAL*8 :: friclambda
  REAL*8 :: sdqv2, sdqv3, sdqvt1
  REAL*8 :: acritfac, pblfrac, autorampb
  REAL*8 :: maxdallowed, rhmn, rhmx
  REAL*8 :: mxdiam
  REAL*8 :: tx2, tx3, akm, acr, alm, tth, qqh, dqx
  REAL*8 :: tx2d, akmd, almd
  REAL*8 :: wfn, tem, trg, trgexp, evp, wlq, qcc
  REAL*8 :: wfnd, temd, trgd
  REAL*8 :: cli, te_a, c00_x, cli_crit_x, toki
  REAL*8 :: tokid
  REAL*8 :: dt_lyr, rate, cvw_x, closs, f2, f3, f4
  REAL*8 :: wght0, prcbl, rndu
  REAL*8 :: lambda_min, lambda_max
  REAL*8 :: tpert, qpert
  REAL*8 :: tpertd
  REAL*8, DIMENSION(k0) :: poi_sv, qoi_sv
  REAL*8, DIMENSION(k0) :: poi_svd, qoi_svd
  REAL*8, DIMENSION(k0) :: poi, qoi, dqq, bet, gam, cll
  REAL*8, DIMENSION(k0) :: poid, qoid, dqqd, betd, gamd
  REAL*8, DIMENSION(k0) :: poi_c, qoi_c
  REAL*8, DIMENSION(k0) :: poi_cd, qoi_cd
  REAL*8, DIMENSION(k0) :: prh, pri, ght, dpt, dpb, pki
  REAL*8, DIMENSION(k0) :: prhd, prid, ghtd, dptd, dpbd, pkid
  REAL*8, DIMENSION(k0) :: cln, rns, pol, dm
  REAL*8, DIMENSION(k0) :: pold
  REAL*8, DIMENSION(k0) :: qst, ssl, rmf, rnn, rn1, rmfc, rmfp
  REAL*8, DIMENSION(k0) :: qstd, ssld
  REAL*8, DIMENSION(k0) :: gms, eta, gmh, eht, gm1, hcc, rmfd
  REAL*8, DIMENSION(k0) :: gmsd, etad, gmhd, ehtd, hccd
  REAL*8, DIMENSION(k0) :: hol, hst, qol, zol, hcld, cll0, cllx, clli
  REAL*8, DIMENSION(k0) :: hold, hstd, qold, zold
  REAL*8, DIMENSION(k0) :: bke, cvw, updfrc
  REAL*8, DIMENSION(k0) :: rasal, updfrp, bk2, dll0, dllx
  REAL*8, DIMENSION(k0) :: rasald
  REAL*8, DIMENSION(k0) :: wght, massf
  REAL*8, DIMENSION(k0) :: wghtd
  REAL*8, DIMENSION(k0) :: qss, dqs, pf, pk, tempf, zlo
  REAL*8, DIMENSION(k0) :: qssd, dqsd, tempfd, zlod
  REAL*8, DIMENSION(k0 + 1) :: prj, prs, qht, sht, zet, zle, pke
  REAL*8, DIMENSION(k0+1) :: prjd, prsd, qhtd, shtd, zetd, zled
  INTRINSIC MAX
  INTRINSIC MIN
  INTRINSIC SQRT
  INTRINSIC EXP
  INTRINSIC SUM
  REAL*8, DIMENSION(k0+1) :: pwx1
  REAL*8 :: pwy1
  REAL*8, DIMENSION(k0) :: pwx10
  REAL*8 :: pwx11
  REAL*8 :: pwr1
  REAL*8 :: arg1
  REAL*8 :: x1
  REAL*8 :: x1d
  REAL*8 :: max2
  REAL*8 :: max1
  REAL*8 :: y1
!  ---  1
  fricfac = rasparams(1)
!  ---  4
  cli_crit = rasparams(4)
!  ---  5
  rasal1 = rasparams(5)
!  ---  6
  rasal2 = rasparams(6)
!  --- 11
  friclambda = rasparams(11)
!  --- 14
  sdqv2 = rasparams(14)
!  --- 15
  sdqv3 = rasparams(15)
!  --- 16
  sdqvt1 = rasparams(16)
!  --- 17
  acritfac = rasparams(17)
!  --- 20
  pblfrac = rasparams(20)
!  --- 21
  autorampb = rasparams(21)
!  --- 24
  rhmn = rasparams(24)
!  --- 24
  maxdallowed = rasparams(23)
!  --- 25
  rhmx = rasparams(25)
  grav = cons_grav
  alhl = cons_alhl
  cp = cons_cp
  cpi = 1.0/cp
  alhi = 1.0/alhl
  gravi = 1.0/grav
  cpbg = cp*gravi
  ddt = daylen/dt
  lbcp = alhl*cpi
  ehtd = 0.0_8
  hccd = 0.0_8
  dqqd = 0.0_8
  dqsd = 0.0_8
  ghtd = 0.0_8
  hstd = 0.0_8
  betd = 0.0_8
  qhtd = 0.0_8
  qold = 0.0_8
  rasald = 0.0_8
  etad = 0.0_8
  shtd = 0.0_8
  gmhd = 0.0_8
  qssd = 0.0_8
  qstd = 0.0_8
  gmsd = 0.0_8
  ssld = 0.0_8
  zetd = 0.0_8
  zold = 0.0_8
  gamd = 0.0_8
  DO i=1,irun
!CALL FINDBASE
    k = kcbl(i)
    IF (k .GT. 0) THEN
!Get saturation specific humidity and gradient wrt to T
      pwx1 = ple(i, :)/1000.
      pwy1 = cons_rgas/cons_cp
      pke = pwx1**pwy1
      pf = 0.5*(ple(i, 1:k0)+ple(i, 2:k0+1))
      pwx10 = pf/1000.
      pwy1 = cons_rgas/cons_cp
      pk = pwx10**pwy1
      tempfd = pk*thod(i, :)
      tempf = tho(i, :)*pk
      zle = 0.0
      zlo = 0.0
      zled(k0+1) = 0.0_8
      zle(k0+1) = 0.
      zlod = 0.0_8
      zled = 0.0_8
      DO l=k0,1,-1
        zled(l) = thod(i, l)*(1.+cons_vireps*qho(i, l)) + tho(i, l)*&
&         cons_vireps*qhod(i, l)
        zle(l) = tho(i, l)*(1.+cons_vireps*qho(i, l))
        zlod(l) = zled(l+1) + cons_cp*(pke(l+1)-pk(l))*zled(l)/cons_grav
        zlo(l) = zle(l+1) + cons_cp/cons_grav*(pke(l+1)-pk(l))*zle(l)
        zled(l) = zlod(l) + cons_cp*(pk(l)-pke(l))*zled(l)/cons_grav
        zle(l) = zlo(l) + cons_cp/cons_grav*(pk(l)-pke(l))*zle(l)
      END DO
      tpertd = cbl_tpert*(-tempfd(k0)-cons_grav*zlod(k0)/cons_cp)
      tpert = cbl_tpert*(ts(i)-(tempf(k0)+cons_grav*zlo(k0)/cons_cp))
!* ( QSSFC - Q(:,:,K0) ) [CBL_QPERT = 0.0]
      qpert = cbl_qpert
      IF (tpert .LT. 0.0) THEN
        tpert = 0.0
        tpertd = 0.0_8
      ELSE
        tpert = tpert
      END IF
      IF (qpert .LT. 0.0) THEN
        qpert = 0.0
      ELSE
        qpert = qpert
      END IF
      IF (frland(i) .LT. 0.1) THEN
        IF (tpert .GT. cbl_tpert_mxocn) THEN
          tpert = cbl_tpert_mxocn
          tpertd = 0.0_8
        ELSE
          tpert = tpert
        END IF
      ELSE IF (tpert .GT. cbl_tpert_mxlnd) THEN
        tpert = cbl_tpert_mxlnd
        tpertd = 0.0_8
      ELSE
        tpert = tpert
      END IF
      CALL DQSAT_RAS_D(dqs, dqsd, qss, qssd, tempf, tempfd, pf, k0, &
&                estblx, cons_h2omw, cons_airmw)
      DO kk=icmin,k+1
        prjd(kk) = 0.0_8
        prj(kk) = pke(kk)
      END DO
! These initialized here in order not to confuse Valgrind debugger
      poi = 0.
! Do not believe it actually makes any difference.
      qoi = 0.
      prsd(icmin:k0+1) = 0.0_8
      prs(icmin:k0+1) = ple(i, icmin:k0+1)
      poid = 0.0_8
      poid(icmin:k) = thod(i, icmin:k)
      poi(icmin:k) = tho(i, icmin:k)
      qoid = 0.0_8
      qoid(icmin:k) = qhod(i, icmin:k)
      qoi(icmin:k) = qho(i, icmin:k)
      qstd(icmin:k) = qssd(icmin:k)
      qst(icmin:k) = qss(icmin:k)
      dqqd(icmin:k) = dqsd(icmin:k)
      dqq(icmin:k) = dqs(icmin:k)
!Mass fraction of each layer below cloud base
      massf(:) = wgt0(i, :)
!RESET PRESSURE at bottom edge of CBL 
      prcbl = prs(k)
      DO l=k,k0
        prcbl = prcbl + massf(l)*(prs(l+1)-prs(l))
      END DO
      prsd(k+1) = 0.0_8
      prs(k+1) = prcbl
      pwx11 = prs(k+1)/1000.
      pwy1 = cons_rgas/cons_cp
      prjd(k+1) = 0.0_8
      prj(k+1) = pwx11**pwy1
      DO l=k,icmin,-1
        pold(l) = 0.0_8
        pol(l) = 0.5*(prs(l)+prs(l+1))
        prhd(l) = 0.0_8
        prh(l) = (prs(l+1)*prj(l+1)-prs(l)*prj(l))/(onepkap*(prs(l+1)-&
&         prs(l)))
        pkid(l) = 0.0_8
        pki(l) = 1.0/prh(l)
        dptd(l) = 0.0_8
        dpt(l) = prh(l) - prj(l)
        dpbd(l) = 0.0_8
        dpb(l) = prj(l+1) - prh(l)
        prid(l) = 0.0_8
        pri(l) = .01/(prs(l+1)-prs(l))
      END DO
!RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
      IF (k .LE. k0) THEN
        poid(k) = 0.0_8
        poi(k) = 0.
        qoid(k) = 0.0_8
        qoi(k) = 0.
!SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
        wght = 0.
        DO l=k,k0
          wghtd(l) = 0.0_8
          wght(l) = massf(l)*(ple(i, l+1)-ple(i, l))/(prs(k+1)-prs(k))
        END DO
        DO l=k,k0
          poid(k) = poid(k) + wght(l)*thod(i, l)
          poi(k) = poi(k) + wght(l)*tho(i, l)
          qoid(k) = qoid(k) + wght(l)*qhod(i, l)
          qoi(k) = qoi(k) + wght(l)*qho(i, l)
        END DO
        CALL DQSATS_RAS_D(dqq(k), dqqd(k), qst(k), qstd(k), poi(k)*prh(k&
&                   ), prh(k)*poid(k), pol(k), estblx, cons_h2omw, &
&                   cons_airmw)
      END IF
      IF (seedras(i)/1000000. .LT. 1e-6) THEN
        rndu = 1e-6
      ELSE
        rndu = seedras(i)/1000000.
      END IF
      pwr1 = rndu**(-(1./2.))
      mxdiam = maxdallowed*pwr1
      DO l=k,icmin,-1
!*
        betd(l) = pki(l)*dqqd(l)
        bet(l) = dqq(l)*pki(l)
!*
        gamd(l) = -(pki(l)*lbcp*dqqd(l)/(1.0+lbcp*dqq(l))**2)
        gam(l) = pki(l)/(1.0+lbcp*dqq(l))
        IF (l .LT. k) THEN
          ghtd(l+1) = dpb(l)*gamd(l) + dpt(l+1)*gamd(l+1)
          ght(l+1) = gam(l)*dpb(l) + gam(l+1)*dpt(l+1)
          gm1(l+1) = 0.5*lbcp*(dqq(l)/(alhl*(1.0+lbcp*dqq(l)))+dqq(l+1)/&
&           (alhl*(1.0+lbcp*dqq(l+1))))
        END IF
      END DO
      rns = 0.
      cll = 0.
      rmf = 0.
      rmfd = 0.
      rmfc = 0.
      rmfp = 0.
      cll0 = 0.
      dll0 = 0.
      cllx = 0.
      dllx = 0.
      clli = 0.
      poi_svd = poid
      poi_sv = poi
      qoi_svd = qoid
      qoi_sv = qoi
      cvw = 0.0
      updfrc = 0.0
      updfrp = 0.0
! HOL initialized here in order not to confuse Valgrind debugger
      hol = 0.
      zetd(k+1) = 0.0_8
      zet(k+1) = 0
      shtd(k+1) = cp*prj(k+1)*poid(k)
      sht(k+1) = cp*poi(k)*prj(k+1)
      hold = 0.0_8
      DO l=k,icmin,-1
        IF (qst(l)*rhmax .GT. qoi(l)) THEN
          qold(l) = qoid(l)
          qol(l) = qoi(l)
        ELSE
          qold(l) = rhmax*qstd(l)
          qol(l) = qst(l)*rhmax
        END IF
        IF (0.000 .LT. qol(l)) THEN
          qol(l) = qol(l)
        ELSE
          qold(l) = 0.0_8
          qol(l) = 0.000
        END IF
        ssld(l) = cp*prj(l+1)*poid(l) + grav*zetd(l+1)
        ssl(l) = cp*prj(l+1)*poi(l) + grav*zet(l+1)
        hold(l) = ssld(l) + alhl*qold(l)
        hol(l) = ssl(l) + qol(l)*alhl
        hstd(l) = ssld(l) + alhl*qstd(l)
        hst(l) = ssl(l) + qst(l)*alhl
        temd = (prj(l+1)-prj(l))*cpbg*poid(l)
        tem = poi(l)*(prj(l+1)-prj(l))*cpbg
        zetd(l) = zetd(l+1) + temd
        zet(l) = zet(l+1) + tem
        zold(l) = zetd(l+1) + (prj(l+1)-prh(l))*cpbg*poid(l)
        zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi(l)*cpbg
      END DO
      DO ic=k+1,icmin+1,-1
        alm = 0.
        IF (1. .GT. (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)) THEN
          trgd = (qoid(k)*qst(k)-qoi(k)*qstd(k))/qst(k)**2/(rhmx-rhmn)
          trg = (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)
        ELSE
          trg = 1.
          trgd = 0.0_8
        END IF
        IF (0.0 .LT. (autorampb-sige(ic))/0.2) THEN
          y1 = (autorampb-sige(ic))/0.2
        ELSE
          y1 = 0.0
        END IF
        IF (1.0 .GT. y1) THEN
          f4 = y1
        ELSE
          f4 = 1.0
        END IF
        IF (trg .GT. 1.0e-5) THEN
!================>>
!RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL
          poi_cd = poid
          poi_c = poi
          qoi_cd = qoid
          qoi_c = qoi
          poi_cd(k) = poi_cd(k) + tpertd
          poi_c(k) = poi_c(k) + tpert
          qoi_c(k) = qoi_c(k) + qpert
          zetd(k+1) = 0.0_8
          zet(k+1) = 0.
          shtd(k+1) = cp*prj(k+1)*poi_cd(k)
          sht(k+1) = cp*poi_c(k)*prj(k+1)
          DO l=k,ic,-1
            IF (qst(l)*rhmax .GT. qoi_c(l)) THEN
              qold(l) = qoi_cd(l)
              qol(l) = qoi_c(l)
            ELSE
              qold(l) = rhmax*qstd(l)
              qol(l) = qst(l)*rhmax
            END IF
            IF (0.000 .LT. qol(l)) THEN
              qol(l) = qol(l)
            ELSE
              qold(l) = 0.0_8
              qol(l) = 0.000
            END IF
            ssld(l) = cp*prj(l+1)*poi_cd(l) + grav*zetd(l+1)
            ssl(l) = cp*prj(l+1)*poi_c(l) + grav*zet(l+1)
            hold(l) = ssld(l) + alhl*qold(l)
            hol(l) = ssl(l) + qol(l)*alhl
            hstd(l) = ssld(l) + alhl*qstd(l)
            hst(l) = ssl(l) + qst(l)*alhl
            temd = (prj(l+1)-prj(l))*cpbg*poi_cd(l)
            tem = poi_c(l)*(prj(l+1)-prj(l))*cpbg
            zetd(l) = zetd(l+1) + temd
            zet(l) = zet(l+1) + tem
            zold(l) = zetd(l+1) + (prj(l+1)-prh(l))*cpbg*poi_cd(l)
            zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi_c(l)*cpbg
          END DO
          DO l=ic+1,k
            tem = (prj(l)-prh(l-1))/(prh(l)-prh(l-1))
            shtd(l) = ssld(l-1) + tem*(ssld(l)-ssld(l-1))
            sht(l) = ssl(l-1) + tem*(ssl(l)-ssl(l-1))
            qhtd(l) = .5*(qold(l)+qold(l-1))
            qht(l) = .5*(qol(l)+qol(l-1))
          END DO
!CALCULATE LAMBDA, ETA, AND WORKFUNCTION
          lambda_min = .2/mxdiam
          lambda_max = .2/200.
          IF (hol(k) .GT. hst(ic)) THEN
!================>>
!LAMBDA CALCULATION: MS-A18
            temd = (hstd(ic)-hold(ic))*(zol(ic)-zet(ic+1)) + (hst(ic)-&
&             hol(ic))*(zold(ic)-zetd(ic+1))
            tem = (hst(ic)-hol(ic))*(zol(ic)-zet(ic+1))
            DO l=ic+1,k-1
              temd = temd + (hstd(ic)-hold(l))*(zet(l)-zet(l+1)) + (hst(&
&               ic)-hol(l))*(zetd(l)-zetd(l+1))
              tem = tem + (hst(ic)-hol(l))*(zet(l)-zet(l+1))
            END DO
            IF (tem .GT. 0.0) THEN
!================>>
              almd = ((hold(k)-hstd(ic))*tem-(hol(k)-hst(ic))*temd)/tem&
&               **2
              alm = (hol(k)-hst(ic))/tem
              IF (alm .LE. lambda_max) THEN
!================>>
                toki = 1.0
                IF (alm .LT. lambda_min) THEN
                  tokid = 2*alm*almd/lambda_min**2
                  toki = (alm/lambda_min)**2
                ELSE
                  tokid = 0.0_8
                END IF
!ETA CALCULATION: MS-A2
                DO l=ic+1,k
                  etad(l) = almd*(zet(l)-zet(k)) + alm*(zetd(l)-zetd(k))
                  eta(l) = 1.0 + alm*(zet(l)-zet(k))
                END DO
                etad(ic) = almd*(zol(ic)-zet(k)) + alm*(zold(ic)-zetd(k)&
&                 )
                eta(ic) = 1.0 + alm*(zol(ic)-zet(k))
!WORKFUNCTION CALCULATION:  MS-A22
                wfn = 0.0
                hccd(k) = hold(k)
                hcc(k) = hol(k)
                wfnd = 0.0_8
                DO l=k-1,ic+1,-1
                  hccd(l) = hccd(l+1) + (etad(l)-etad(l+1))*hol(l) + (&
&                   eta(l)-eta(l+1))*hold(l)
                  hcc(l) = hcc(l+1) + (eta(l)-eta(l+1))*hol(l)
                  temd = dpb(l)*hccd(l+1) + dpt(l)*hccd(l)
                  tem = hcc(l+1)*dpb(l) + hcc(l)*dpt(l)
                  ehtd(l) = dpb(l)*etad(l+1) + dpt(l)*etad(l)
                  eht(l) = eta(l+1)*dpb(l) + eta(l)*dpt(l)
                  wfnd = wfnd + (temd-ehtd(l)*hst(l)-eht(l)*hstd(l))*gam&
&                   (l) + (tem-eht(l)*hst(l))*gamd(l)
                  wfn = wfn + (tem-eht(l)*hst(l))*gam(l)
                END DO
                hccd(ic) = hstd(ic)*eta(ic) + hst(ic)*etad(ic)
                hcc(ic) = hst(ic)*eta(ic)
                wfnd = wfnd + dpb(ic)*((hccd(ic+1)-hstd(ic)*eta(ic+1)-&
&                 hst(ic)*etad(ic+1))*gam(ic)+(hcc(ic+1)-hst(ic)*eta(ic+&
&                 1))*gamd(ic))
                wfn = wfn + (hcc(ic+1)-hst(ic)*eta(ic+1))*gam(ic)*dpb(ic&
&                 )
!VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
                bk2(k) = 0.0
                bke(k) = 0.0
                hcld(k) = hol(k)
                DO l=k-1,ic,-1
                  hcld(l) = (eta(l+1)*hcld(l+1)+(eta(l)-eta(l+1))*hol(l)&
&                   )/eta(l)
                  tem = (hcld(l)-hst(l))*(zet(l)-zet(l+1))/(1.0+lbcp*dqq&
&                   (l))
                  bke(l) = bke(l+1) + grav*tem/(cp*prj(l+1)*poi(l))
                  IF (tem .LT. 0.0) THEN
                    max1 = 0.0
                  ELSE
                    max1 = tem
                  END IF
                  bk2(l) = bk2(l+1) + grav*max1/(cp*prj(l+1)*poi(l))
                  IF (bk2(l) .LT. 0.0) THEN
                    max2 = 0.0
                  ELSE
                    max2 = bk2(l)
                  END IF
                  cvw(l) = SQRT(2.0*max2)
                END DO
!ALPHA CALCULATION 
                IF (zet(ic) .LT. 2000.) THEN
                  rasald(ic) = 0.0_8
                  rasal(ic) = rasal1
                END IF
                IF (zet(ic) .GE. 2000.) THEN
                  rasald(ic) = (rasal2-rasal1)*zetd(ic)/8000.
                  rasal(ic) = rasal1 + (rasal2-rasal1)*(zet(ic)-2000.)/&
&                   8000.
                END IF
                IF (rasal(ic) .GT. 1.0e5) THEN
                  rasald(ic) = 0.0_8
                  rasal(ic) = 1.0e5
                ELSE
                  rasal(ic) = rasal(ic)
                END IF
                rasald(ic) = -(dt*rasald(ic)/rasal(ic)**2)
                rasal(ic) = dt/rasal(ic)
                WHERE (cvw(ic:k) .LT. 1.00) 
                  cvw(ic:k) = 1.00
                ELSEWHERE
                  cvw(ic:k) = cvw(ic:k)
                END WHERE
                CALL ACRITN(pol(ic), prs(k), acr, acritfac)
                IF (wfn .GT. acr) THEN
!================>>
                  wlq = qol(k)
                  rnn(k) = 0.
                  cll0(k) = 0.
                  DO l=k-1,ic,-1
                    tem = eta(l) - eta(l+1)
                    wlq = wlq + tem*qol(l)
                    IF (l .GT. ic) THEN
                      tx2 = 0.5*(qst(l)+qst(l-1))*eta(l)
                      tx3 = 0.5*(hst(l)+hst(l-1))*eta(l)
                      qcc = tx2 + gm1(l)*(hcc(l)-tx3)
                      cll0(l) = wlq - qcc
                    ELSE
                      cll0(l) = wlq - qst(ic)*eta(ic)
                    END IF
                    IF (cll0(l) .LT. 0.00) THEN
                      cll0(l) = 0.00
                    ELSE
                      cll0(l) = cll0(l)
                    END IF
                    cli = cll0(l)/eta(l)
                    te_a = poi(l)*prh(l)
                    CALL SUNDQ3_ICE(te_a, sdqv2, sdqv3, sdqvt1, f2, f3)
                    c00_x = co_auto(i)*f2*f3*f4
                    cli_crit_x = cli_crit/(f2*f3)
                    arg1 = -(cli**2/cli_crit_x**2)
                    rate = c00_x*(1.0-EXP(arg1))
                    IF (cvw(l) .LT. 1.00) THEN
                      cvw_x = 1.00
                    ELSE
                      cvw_x = cvw(l)
                    END IF
                    dt_lyr = (zet(l)-zet(l+1))/cvw_x
                    closs = cll0(l)*rate*dt_lyr
                    IF (closs .GT. cll0(l)) THEN
                      closs = cll0(l)
                    ELSE
                      closs = closs
                    END IF
                    cll0(l) = cll0(l) - closs
                    dll0(l) = closs
                    IF (closs .GT. 0.) THEN
                      wlq = wlq - closs
                      rnn(l) = closs
                    ELSE
                      rnn(l) = 0.
                    END IF
                  END DO
                  wlq = wlq - qst(ic)*eta(ic)
!CALCULATE GAMMAS AND KERNEL
                  gmsd(k) = pri(k)*(shtd(k)-ssld(k))
                  gms(k) = (sht(k)-ssl(k))*pri(k)
                  gmhd(k) = gmsd(k) + pri(k)*alhl*(qhtd(k)-qold(k))
                  gmh(k) = gms(k) + (qht(k)-qol(k))*pri(k)*alhl
                  akmd = dpb(k-1)*(gmhd(k)*gam(k-1)+gmh(k)*gamd(k-1))
                  akm = gmh(k)*gam(k-1)*dpb(k-1)
                  tx2d = gmhd(k)
                  tx2 = gmh(k)
                  DO l=k-1,ic+1,-1
                    gmsd(l) = pri(l)*(etad(l)*(sht(l)-ssl(l))+eta(l)*(&
&                     shtd(l)-ssld(l))+etad(l+1)*(ssl(l)-sht(l+1))+eta(l&
&                     +1)*(ssld(l)-shtd(l+1)))
                    gms(l) = (eta(l)*(sht(l)-ssl(l))+eta(l+1)*(ssl(l)-&
&                     sht(l+1)))*pri(l)
                    gmhd(l) = gmsd(l) + alhl*pri(l)*(etad(l)*(qht(l)-qol&
&                     (l))+eta(l)*(qhtd(l)-qold(l))+etad(l+1)*(qol(l)-&
&                     qht(l+1))+eta(l+1)*(qold(l)-qhtd(l+1)))
                    gmh(l) = gms(l) + (eta(l)*(qht(l)-qol(l))+eta(l+1)*(&
&                     qol(l)-qht(l+1)))*alhl*pri(l)
                    tx2d = tx2d + (etad(l)-etad(l+1))*gmh(l) + (eta(l)-&
&                     eta(l+1))*gmhd(l)
                    tx2 = tx2 + (eta(l)-eta(l+1))*gmh(l)
                    akmd = akmd - pki(l)*(gmsd(l)*eht(l)+gms(l)*ehtd(l))&
&                     + tx2d*ght(l) + tx2*ghtd(l)
                    akm = akm - gms(l)*eht(l)*pki(l) + tx2*ght(l)
                  END DO
                  gmsd(ic) = pri(ic)*(etad(ic+1)*(ssl(ic)-sht(ic+1))+eta&
&                   (ic+1)*(ssld(ic)-shtd(ic+1)))
                  gms(ic) = eta(ic+1)*(ssl(ic)-sht(ic+1))*pri(ic)
                  akmd = akmd - dpb(ic)*pki(ic)*(gmsd(ic)*eta(ic+1)+gms(&
&                   ic)*etad(ic+1))
                  akm = akm - gms(ic)*eta(ic+1)*dpb(ic)*pki(ic)
                  gmhd(ic) = gmsd(ic) + pri(ic)*(alhl*(etad(ic+1)*(qol(&
&                   ic)-qht(ic+1))+eta(ic+1)*(qold(ic)-qhtd(ic+1)))+etad&
&                   (ic)*(hst(ic)-hol(ic))+eta(ic)*(hstd(ic)-hold(ic)))
                  gmh(ic) = gms(ic) + (eta(ic+1)*(qol(ic)-qht(ic+1))*&
&                   alhl+eta(ic)*(hst(ic)-hol(ic)))*pri(ic)
!CLOUD BASE MASS FLUX
                  IF (.NOT.(akm .GE. 0.0 .OR. wlq .LT. 0.0)) THEN
!================>>
                    wfnd = -((wfnd*akm-(wfn-acr)*akmd)/akm**2)
                    wfn = -((wfn-acr)/akm)
                    x1d = (rasald(ic)*wfn+rasal(ic)*wfnd)*trg*toki + &
&                     rasal(ic)*wfn*(trgd*toki+trg*tokid)
                    x1 = rasal(ic)*trg*toki*wfn
                    IF (x1 .GT. (prs(k+1)-prs(k))*(100.*pblfrac)) THEN
                      wfn = (prs(k+1)-prs(k))*(100.*pblfrac)
                      wfnd = 0.0_8
                    ELSE
                      wfnd = x1d
                      wfn = x1
                    END IF
!CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
                    tem = wfn*gravi
                    cll(ic) = cll(ic) + wlq*tem
                    rmf(ic) = rmf(ic) + tem
                    rmfd(ic) = rmfd(ic) + tem*eta(ic)
                    DO l=ic+1,k
                      rmfp(l) = tem*eta(l)
                      rmfc(l) = rmfc(l) + rmfp(l)
                      dllx(l) = dllx(l) + tem*dll0(l)
                      IF (cvw(l) .GT. 0.0) THEN
                        updfrp(l) = rmfp(l)*(ddt/daylen)*1000./(cvw(l)*&
&                         prs(l))
                      ELSE
                        updfrp(l) = 0.0
                      END IF
                      clli(l) = cll0(l)/eta(l)
                      updfrc(l) = updfrc(l) + updfrp(l)
                    END DO
!THETA AND Q CHANGE DUE TO CLOUD TYPE IC
                    DO l=ic,k
                      rns(l) = rns(l) + rnn(l)*tem
                      gmhd(l) = gmhd(l)*wfn + gmh(l)*wfnd
                      gmh(l) = gmh(l)*wfn
                      gmsd(l) = gmsd(l)*wfn + gms(l)*wfnd
                      gms(l) = gms(l)*wfn
                      qoid(l) = qoid(l) + alhi*(gmhd(l)-gmsd(l))
                      qoi(l) = qoi(l) + (gmh(l)-gms(l))*alhi
                      poid(l) = poid(l) + pki(l)*cpi*gmsd(l)
                      poi(l) = poi(l) + gms(l)*pki(l)*cpi
                      qstd(l) = qstd(l) + cpi*(gmsd(l)*bet(l)+gms(l)*&
&                       betd(l))
                      qst(l) = qst(l) + gms(l)*bet(l)*cpi
                    END DO
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
!CLOUD LOOP
      IF (SUM(rmf(icmin:k)) .GT. 0.0) THEN
        thod(i, icmin:k-1) = poid(icmin:k-1)
        tho(i, icmin:k-1) = poi(icmin:k-1)
        qhod(i, icmin:k-1) = qoid(icmin:k-1)
        qho(i, icmin:k-1) = qoi(icmin:k-1)
!De-strap tendencies from RAS
        wght = wgt1(i, :)
!Scale properly by layer masses
        wght0 = 0.
        DO l=k,k0
          wght0 = wght0 + wght(l)*(ple(i, l+1)-ple(i, l))
        END DO
        wght0 = (prs(k+1)-prs(k))/wght0
        wght = wght0*wght
        DO l=k,k0
          thod(i, l) = thod(i, l) + wght(l)*(poid(k)-poi_svd(k))
          tho(i, l) = tho(i, l) + wght(l)*(poi(k)-poi_sv(k))
          qhod(i, l) = qhod(i, l) + wght(l)*(qoid(k)-qoi_svd(k))
          qho(i, l) = qho(i, l) + wght(l)*(qoi(k)-qoi_sv(k))
        END DO
      END IF
    END IF
  END DO
END SUBROUTINE RASE0_D

!  Differentiation of dqsat_ras in forward (tangent) mode:
!   variations   of useful results: dqsi qssi
!   with respect to varying inputs: temp dqsi qssi
SUBROUTINE DQSAT_RAS_D(dqsi, dqsid, qssi, qssid, temp, tempd, plo, lm, &
& estblx, cons_h2omw, cons_airmw)
  IMPLICIT NONE
!Inputs
  INTEGER :: lm
  REAL*8, DIMENSION(lm) :: temp, plo
  REAL*8, DIMENSION(lm) :: tempd
  REAL*8 :: estblx(:)
  REAL*8 :: cons_h2omw, cons_airmw
!Outputs
  REAL*8, DIMENSION(lm) :: dqsi, qssi
  REAL*8, DIMENSION(lm) :: dqsid, qssid
!Locals
  REAL*8, PARAMETER :: max_mixing_ratio=1.0
  REAL*8 :: esfac
  INTEGER :: k
  REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
  REAL*8 :: tld, ttd, tid, dqsatd, qsatd, qqd, ddd
  INTEGER :: it
  INTEGER, PARAMETER :: degsubs=100
  REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
  INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
  INTRINSIC NINT
  INTRINSIC INT
  esfac = cons_h2omw/cons_airmw
  DO k=1,lm
    tld = tempd(k)
    tl = temp(k)
    pl = plo(k)
    pp = pl*100.0
    IF (tl .LE. tmintbl) THEN
      ti = tmintbl
      tid = 0.0_8
    ELSE IF (tl .GE. tmaxtbl - .001) THEN
      tid = 0.0_8
      ti = tmaxtbl - .001
    ELSE
      tid = tld
      ti = tl
    END IF
    ttd = degsubs*tid
    tt = (ti-tmintbl)*degsubs + 1
    it = INT(tt)
    dqq = estblx(it+1) - estblx(it)
    qqd = dqq*ttd
    qq = (tt-it)*dqq + estblx(it)
    IF (pp .LE. qq) THEN
      qsat = max_mixing_ratio
      dqsat = 0.0
      qsatd = 0.0_8
      dqsatd = 0.0_8
    ELSE
      ddd = -((-((1.0-esfac)*qqd))/(pp-(1.0-esfac)*qq)**2)
      dd = 1.0/(pp-(1.0-esfac)*qq)
      qsatd = esfac*(qqd*dd+qq*ddd)
      qsat = esfac*qq*dd
      dqsatd = esfac*degsubs*dqq*pp*(ddd*dd+dd*ddd)
      dqsat = esfac*degsubs*dqq*pp*(dd*dd)
    END IF
    dqsid(k) = dqsatd
    dqsi(k) = dqsat
    qssid(k) = qsatd
    qssi(k) = qsat
  END DO
END SUBROUTINE DQSAT_RAS_D

!  Differentiation of dqsats_ras in forward (tangent) mode:
!   variations   of useful results: dqsi qssi
!   with respect to varying inputs: temp
SUBROUTINE DQSATS_RAS_D(dqsi, dqsid, qssi, qssid, temp, tempd, plo, &
& estblx, cons_h2omw, cons_airmw)
  IMPLICIT NONE
!Inputs
  REAL*8 :: temp, plo
  REAL*8 :: tempd
  REAL*8 :: estblx(:)
  REAL*8 :: cons_h2omw, cons_airmw
!Outputs
  REAL*8 :: dqsi, qssi
  REAL*8 :: dqsid, qssid
!Locals
  REAL*8, PARAMETER :: max_mixing_ratio=1.0
  REAL*8 :: esfac
  REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
  REAL*8 :: tld, ttd, tid, dqsatd, qsatd, qqd, ddd
  INTEGER :: it
  INTEGER, PARAMETER :: degsubs=100
  REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
  INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
  INTRINSIC NINT
  INTRINSIC INT
  esfac = cons_h2omw/cons_airmw
  tld = tempd
  tl = temp
  pl = plo
  pp = pl*100.0
  IF (tl .LE. tmintbl) THEN
    ti = tmintbl
    tid = 0.0_8
  ELSE IF (tl .GE. tmaxtbl - .001) THEN
    ti = tmaxtbl - .001
    tid = 0.0_8
  ELSE
    tid = tld
    ti = tl
  END IF
  ttd = degsubs*tid
  tt = (ti-tmintbl)*degsubs + 1
  it = INT(tt)
  dqq = estblx(it+1) - estblx(it)
  qqd = dqq*ttd
  qq = (tt-it)*dqq + estblx(it)
  IF (pp .LE. qq) THEN
    qsat = max_mixing_ratio
    dqsat = 0.0
    qsatd = 0.0_8
    dqsatd = 0.0_8
  ELSE
    ddd = -((-((1.0-esfac)*qqd))/(pp-(1.0-esfac)*qq)**2)
    dd = 1.0/(pp-(1.0-esfac)*qq)
    qsatd = esfac*(qqd*dd+qq*ddd)
    qsat = esfac*qq*dd
    dqsatd = esfac*degsubs*dqq*pp*(ddd*dd+dd*ddd)
    dqsat = esfac*degsubs*dqq*pp*(dd*dd)
  END IF
  dqsid = dqsatd
  dqsi = dqsat
  qssid = qsatd
  qssi = qsat
END SUBROUTINE DQSATS_RAS_D

