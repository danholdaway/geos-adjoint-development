!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.7 (r4786) - 21 Feb 2013 15:53
!
!CNV_DQLDT
!CNV_MFD
!CNV_PRC3
!CNV_UPDFRC
SUBROUTINE RASE(k0, icmin, dt, seedras, sige, kcbl, wgt0, wgt1, frland, &
&  ts, tho, qho, uho, vho, co_auto, ple, clw, flxd, cnv_prc3, cnv_updfrc&
&  , rasparams, cons_cp, cons_alhl, cons_alhs, cons_grav, cons_undef, &
&  cons_rgas, cons_h2omw, cons_airmw, cons_vireps, estblx, tap_dummy, &
&  itrcr, xho, fscav, disske)
  IMPLICIT NONE
!INPUTS
  INTEGER, INTENT(IN) :: k0, icmin, kcbl, seedras(2)
  REAL*8, INTENT(IN) :: dt, frland, ts, co_auto, tap_dummy
  REAL*8, DIMENSION(k0), INTENT(IN) :: wgt0, wgt1
  REAL*8, DIMENSION(k0 + 1), INTENT(IN) :: ple, sige
  REAL*8, DIMENSION(:), INTENT(IN) :: rasparams, estblx
  REAL*4, INTENT(IN) :: cons_alhs, cons_cp, cons_alhl
  REAL*4, INTENT(IN) :: cons_grav, cons_undef, cons_rgas
  REAL*4, INTENT(IN) :: cons_h2omw, cons_airmw, cons_vireps
!OUTPUTS
! real(8), DIMENSION (K0+1), INTENT(INOUT) ::  FLXC
! real(8), DIMENSION (K0  ), INTENT(INOUT) ::  FLX, CNV_CVW, CNV_QC
! real(8), DIMENSION (K0+1)              ::  FLXC
! real(8), DIMENSION (K0  )              ::  FLX, CNV_CVW, CNV_QC
!PROGNOSTIC
  REAL*8, DIMENSION(k0), INTENT(INOUT) :: tho, qho, uho, vho
  REAL*8, DIMENSION(k0), INTENT(INOUT) :: clw, flxd, cnv_prc3, &
&  cnv_updfrc
!SCAVENGING
  INTEGER, INTENT(IN) :: itrcr
! = 0 (no scav), = 1 (full scav)
  REAL*8, OPTIONAL, INTENT(IN) :: fscav(itrcr)
  REAL*8, OPTIONAL, INTENT(INOUT) :: disske(k0)
  REAL*8, OPTIONAL, INTENT(INOUT) :: xho(k0, itrcr)
!DIAGNOSTICS
!LOCALS
  REAL*8, PARAMETER :: onepkap=1.+2./7., daylen=86400.0
  REAL*8, PARAMETER :: rhmax=0.9999
  REAL*8, PARAMETER :: cbl_qpert=0.0, cbl_tpert=1.0
  REAL*8, PARAMETER :: cbl_tpert_mxocn=2.0, cbl_tpert_mxlnd=4.0
!Constants
  REAL*8 :: cpbg, alhi, cpi, gravi, ddt, lbcp
!Rasparams
  REAL*8 :: fricfac, cli_crit, rasal1, rasal2, rasncl
  REAL*8 :: lambda_fac, lambmx_fac, diammn_min, friclambda
  REAL*8 :: rdtlexpon, strapping
  REAL*8 :: sdqv2, sdqv3, sdqvt1
  REAL*8 :: acritfac, hmintrigger, lldisaggxp, pblfrac, autorampb
  REAL*8 :: co_zdep, maxdallowed, rhmn, rhmx
  INTEGER :: k, kk, ic, l, itr, rc(k0-1)
  INTEGER :: icl_v(k0)
  REAL*8 :: lambda_min, lambda_max
  REAL*8 :: tpert, qpert, mxdiam
  REAL*8 :: tx2, tx3, uht, vht, akm, acr, alm, tth, qqh, dqx
  REAL*8 :: wfn, tem, trg, wlq, qcc
  REAL*8 :: wfnog
  REAL*8 :: prcbl, rndu
  REAL*8 :: cli, te_a, c00_x, cli_crit_x, toki, gmhx, hstx
  REAL*8 :: dt_lyr, rate, cvw_x, closs, f2, f3, f4, f5
  REAL*8, DIMENSION(k0 + 1) :: prj, prs, qht, sht, zet
  REAL*8, DIMENSION(k0 + 1) :: pke, zle
  REAL*8, DIMENSION(k0) :: qss, dqs, tempf, pf, pk, zlo
  REAL*8, DIMENSION(k0) :: poi_sv, qoi_sv, uoi_sv, voi_sv
  REAL*8, DIMENSION(k0) :: poi, qoi, uoi, voi, dqq, bet, gam, cll
  REAL*8, DIMENSION(k0) :: poi_c, qoi_c
  REAL*8, DIMENSION(k0) :: prh, pri, ght, dpt, dpb, pki, dissk0
  REAL*8, DIMENSION(k0) :: ucu, vcu, rns, pol
  REAL*8, DIMENSION(k0) :: qst, ssl, rmf, rnn, rn1, rmfc, rmfp
  REAL*8, DIMENSION(k0) :: gms, eta, gmh, eht, gm1, hcc, rmfd
  REAL*8, DIMENSION(k0) :: hol, hst, qol, zol, hcld, cll0, clli, cllb
  REAL*8, DIMENSION(k0) :: bke, cvw, updfrc
  REAL*8, DIMENSION(k0) :: rasal, updfrp, bk2
  REAL*8, DIMENSION(k0) :: wght, wght0, wght1, massf
! real(8), DIMENSION(K0+1) :: FLXCTMP
! real(8), DIMENSION(K0)   :: FLXTMP, CNV_CVWTMP, CNV_QCTMP
  REAL*8, DIMENSION(k0) :: clwtmp, flxdtmp, cnv_prc3tmp, cnv_updfrctmp, &
&  dissketmp
!SCAVANGING RELATED PARAMETERS
  LOGICAL :: do_tracers
! layer thickness in km
  REAL*8 :: delzkm
! fraction of tracer *not* scavenged
  REAL*8 :: fnoscav
! Fraction scavenged per km
  REAL*8 :: fscav_(itrcr)
  REAL*8, DIMENSION(itrcr) :: xht
  REAL*8, DIMENSION(k0, itrcr) :: xoi, xcu, xoi_sv
  INTRINSIC EXP
  INTRINSIC MAX
  REAL*8 :: x5
  REAL*8 :: x4
  REAL*8 :: x3
  REAL*8 :: x2
  REAL*8 :: x1
  INTRINSIC PRESENT
  INTRINSIC SUM
  INTRINSIC MIN
  INTRINSIC SQRT
  REAL*8 :: max2
  REAL*8 :: max1
! BEGIN CALCULATIONS
!-------------------
!Initialize local arrays
  prj = 0.0
  prs = 0.0
  qht = 0.0
  sht = 0.0
  zet = 0.0
  zle = 0.0
  qss = 0.0
  dqs = 0.0
  tempf = 0.0
  pf = 0.0
  pk = 0.0
  zlo = 0.0
  poi_sv = 0.0
  qoi_sv = 0.0
  uoi_sv = 0.0
  voi_sv = 0.0
  poi = 0.0
  qoi = 0.0
  uoi = 0.0
  voi = 0.0
  dqq = 0.0
  bet = 0.0
  gam = 0.0
  cll = 0.0
  poi_c = 0.0
  qoi_c = 0.0
  prh = 0.0
  pri = 0.0
  ght = 0.0
  dpt = 0.0
  dpb = 0.0
  pki = 0.0
  dissk0 = 0.0
  ucu = 0.0
  vcu = 0.0
  rns = 0.0
  pol = 0.0
  qst = 0.0
  ssl = 0.0
  rmf = 0.0
  rnn = 0.0
  rn1 = 0.0
  rmfc = 0.0
  rmfp = 0.0
  gms = 0.0
  eta = 0.0
  gmh = 0.0
  eht = 0.0
  gm1 = 0.0
  hcc = 0.0
  rmfd = 0.0
  hol = 0.0
  hst = 0.0
  qol = 0.0
  zol = 0.0
  hcld = 0.0
  cll0 = 0.0
  clli = 0.0
  cllb = 0.0
  bke = 0.0
  cvw = 0.0
  updfrc = 0.0
  rasal = 0.0
  updfrp = 0.0
  bk2 = 0.0
  wght = 0.0
  wght1 = 0.0
  massf = 0.0
  clwtmp = 0.0
  flxdtmp = 0.0
  cnv_prc3tmp = 0.0
  cnv_updfrctmp = 0.0
  dissketmp = 0.0
  xht = 0.0
  xoi = 0.0
  xcu = 0.0
  xoi_sv = 0.0
  IF (PRESENT(fscav)) THEN
    fscav_ = fscav
  ELSE
! NO SCAVENGING BY DEFAULT
    fscav_ = 0.0
  END IF
!Make outputs prognostic to trick the autodiff
  clwtmp = clw*tap_dummy
  flxdtmp = flxd*tap_dummy
  cnv_prc3tmp = cnv_prc3*tap_dummy
  cnv_updfrctmp = cnv_updfrc*tap_dummy
! FLXCTMP       = FLXC       * TAP_DUMMY
! FLXTMP        = FLX        * TAP_DUMMY
! CNV_CVWTMP    = CNV_CVW    * TAP_DUMMY
! CNV_QCTMP     = CNV_QC     * TAP_DUMMY
  IF (PRESENT(disske)) dissketmp = disske*tap_dummy
  IF (PRESENT(disske)) disske = 0.
!  ---  1
  fricfac = rasparams(1)
!  ---  4
  cli_crit = rasparams(4)
!  ---  5
  rasal1 = rasparams(5)
!  ---  6
  rasal2 = rasparams(6)
!  ---  7
  rasncl = rasparams(7)
!  ---  8
  lambda_fac = rasparams(8)
!  ---  9
  lambmx_fac = rasparams(9)
!  --- 10
  diammn_min = rasparams(10)
!  --- 11
  friclambda = rasparams(11)
!  --- 12
  rdtlexpon = rasparams(12)
!  --- 13
  strapping = rasparams(13)
!  --- 14
  sdqv2 = rasparams(14)
!  --- 15
  sdqv3 = rasparams(15)
!  --- 16
  sdqvt1 = rasparams(16)
!  --- 17
  acritfac = rasparams(17)
!  --- 18
  hmintrigger = rasparams(18)
!  --- 19
  lldisaggxp = rasparams(19)
!  --- 20
  pblfrac = rasparams(20)
!  --- 21
  autorampb = rasparams(21)
!  --- 22
  co_zdep = rasparams(22)
!  --- 23
  maxdallowed = rasparams(23)
!  --- 24
  rhmn = rasparams(24)
!  --- 25
  rhmx = rasparams(25)
  do_tracers = PRESENT(xho) .AND. itrcr .GT. 0
  cpi = 1.0/cons_cp
  alhi = 1.0/cons_alhl
  gravi = 1.0/cons_grav
  cpbg = cons_cp*gravi
  ddt = daylen/dt
  lbcp = cons_alhl*cpi
!Get saturation specific humidity and gradient wrt to T
  pke = (ple/1000.)**(cons_rgas/cons_cp)
  pf = 0.5*(ple(1:k0)+ple(2:k0+1))
  pk = (pf/1000.)**(cons_rgas/cons_cp)
  tempf = tho(:)*pk
  zle(k0) = 0.
  DO l=k0,1,-1
    zle(l-1) = tho(l)*(1.+cons_vireps*qho(l))
    zlo(l) = zle(l) + cons_cp/cons_grav*(pke(l+1)-pk(l))*zle(l-1)
    zle(l-1) = zlo(l) + cons_cp/cons_grav*(pk(l)-pke(l))*zle(l-1)
  END DO
  tpert = cbl_tpert*(ts-(tempf(k0)+cons_grav*zlo(k0)/cons_cp))
!* ( QSSFC - Q(:,:,K0) ) [CBL_QPERT = 0.0]
  qpert = cbl_qpert
  IF (tpert .LT. 0.0) THEN
    tpert = 0.0
  ELSE
    tpert = tpert
  END IF
  IF (qpert .LT. 0.0) THEN
    qpert = 0.0
  ELSE
    qpert = qpert
  END IF
  IF (frland .LT. 0.1) THEN
    IF (tpert .GT. cbl_tpert_mxocn) THEN
      tpert = cbl_tpert_mxocn
    ELSE
      tpert = tpert
    END IF
  ELSE IF (tpert .GT. cbl_tpert_mxlnd) THEN
    tpert = cbl_tpert_mxlnd
  ELSE
    tpert = tpert
  END IF
  CALL DQSAT_RAS(dqs, qss, tempf, pf, k0, estblx, cons_h2omw, cons_airmw&
&          )
!CALL FINDBASE
  k = kcbl
  rc(icmin) = 0
!---------
!STRAP FINAL=0 SUBROUTINE
!------------------
  DO kk=icmin,k+1
    prj(kk) = pke(kk)
  END DO
! These initialized here in order not to confuse Valgrind debugger
  poi = 0.
! Do not believe it actually makes any difference.
  qoi = 0.
  uoi = 0.
  voi = 0.
  prs(icmin:k0+1) = ple(icmin:k0+1)
  poi(icmin:k) = tho(icmin:k)
  qoi(icmin:k) = qho(icmin:k)
  uoi(icmin:k) = uho(icmin:k)
  voi(icmin:k) = vho(icmin:k)
  qst(icmin:k) = qss(icmin:k)
  dqq(icmin:k) = dqs(icmin:k)
  IF (do_tracers) THEN
    DO itr=1,itrcr
      xoi(icmin:k, itr) = xho(icmin:k, itr)
    END DO
  END IF
!Mass fraction of each layer below cloud base contributed to aggregate cloudbase layer 
  massf(:) = wgt0(:)
!RESET PRESSURE at bottom edge of CBL 
  prcbl = prs(k)
  DO l=k,k0
    prcbl = prcbl + massf(l)*(prs(l+1)-prs(l))
  END DO
  prs(k+1) = prcbl
  prj(k+1) = (prs(k+1)/1000.)**(cons_rgas/cons_cp)
  DO l=k,icmin,-1
    pol(l) = 0.5*(prs(l)+prs(l+1))
    prh(l) = (prs(l+1)*prj(l+1)-prs(l)*prj(l))/(onepkap*(prs(l+1)-prs(l)&
&      ))
    pki(l) = 1.0/prh(l)
    dpt(l) = prh(l) - prj(l)
    dpb(l) = prj(l+1) - prh(l)
    pri(l) = .01/(prs(l+1)-prs(l))
  END DO
!RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
  IF (k .LE. k0) THEN
    poi(k) = 0.
    qoi(k) = 0.
    uoi(k) = 0.
    voi(k) = 0.
!SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
    wght = 0.
    DO l=k,k0
      wght(l) = massf(l)*(ple(l+1)-ple(l))/(prs(k+1)-prs(k))
    END DO
    DO l=k,k0
      poi(k) = poi(k) + wght(l)*tho(l)
      qoi(k) = qoi(k) + wght(l)*qho(l)
      uoi(k) = uoi(k) + wght(l)*uho(l)
      voi(k) = voi(k) + wght(l)*vho(l)
    END DO
    IF (do_tracers) THEN
      xoi(k, :) = 0.
      DO itr=1,itrcr
        DO l=k,k0
          xoi(k, itr) = xoi(k, itr) + wght(l)*xho(l, itr)
        END DO
      END DO
    END IF
    CALL DQSATS_RAS(dqq(k), qst(k), poi(k)*prh(k), pol(k), estblx, &
&              cons_h2omw, cons_airmw)
  END IF
  IF (seedras(1)/1000000. .LT. 1e-6) THEN
    rndu = 1e-6
  ELSE
    rndu = seedras(1)/1000000.
  END IF
  mxdiam = maxdallowed*rndu**(-(1./2.))
  DO l=k,icmin,-1
!*
    bet(l) = dqq(l)*pki(l)
!*
    gam(l) = pki(l)/(1.0+lbcp*dqq(l))
    IF (l .LT. k) THEN
      ght(l+1) = gam(l)*dpb(l) + gam(l+1)*dpt(l+1)
      gm1(l+1) = 0.5*lbcp*(dqq(l)/(cons_alhl*(1.0+lbcp*dqq(l)))+dqq(l+1)&
&        /(cons_alhl*(1.0+lbcp*dqq(l+1))))
    END IF
  END DO
  poi_sv = poi
  qoi_sv = qoi
  uoi_sv = uoi
  voi_sv = voi
  IF (do_tracers) xoi_sv = xoi
!END STRAP (FINAL=0) SUBROUTINE
!------------------------------
!HTEST 
!-----
  hol = 0.0
  zet(k+1) = 0
  sht(k+1) = cons_cp*poi(k)*prj(k+1)
  DO l=k,icmin,-1
    IF (qst(l)*rhmax .GT. qoi(l)) THEN
      qol(l) = qoi(l)
    ELSE
      qol(l) = qst(l)*rhmax
    END IF
    IF (0.000 .LT. qol(l)) THEN
      qol(l) = qol(l)
    ELSE
      qol(l) = 0.000
    END IF
    ssl(l) = cons_cp*prj(l+1)*poi(l) + cons_grav*zet(l+1)
    hol(l) = ssl(l) + qol(l)*cons_alhl
    hst(l) = ssl(l) + qst(l)*cons_alhl
    tem = poi(l)*(prj(l+1)-prj(l))*cpbg
    zet(l) = zet(l+1) + tem
    zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi(l)*cpbg
  END DO
!HTEST
  DO ic=kcbl+1,icmin+1,-1
!MAIN CLOUD SUBROUTINE
!---------------------
!CALL CLOUDE
!INTEGER, INTENT(IN ) :: IC
    alm = 0.0
    lambda_min = 0.0
    f4 = 0.0
    trg = 0.0
    toki = 1.0
    IF (alm .LT. lambda_min) toki = (alm/lambda_min)**2
!ETA CALCULATION: MS-A2
    DO l=ic+1,k
      eta(l) = 1.0 + alm*(zet(l)-zet(k))
    END DO
    eta(ic) = 1.0 + alm*(zol(ic)-zet(k))
!WORKFUNCTION CALCULATION:  MS-A22
    wfn = 0.0
    hcc(k) = hol(k)
    DO l=k-1,ic+1,-1
      hcc(l) = hcc(l+1) + (eta(l)-eta(l+1))*hol(l)
      tem = hcc(l+1)*dpb(l) + hcc(l)*dpt(l)
      eht(l) = eta(l+1)*dpb(l) + eta(l)*dpt(l)
      wfn = wfn + (tem-eht(l)*hst(l))*gam(l)
    END DO
    hcc(ic) = hst(ic)*eta(ic)
    wfn = wfn + (hcc(ic+1)-hst(ic)*eta(ic+1))*gam(ic)*dpb(ic)
!VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
    bk2(k) = 0.0
    bke(k) = 0.0
    hcld(k) = hol(k)
    DO l=k-1,ic,-1
      hcld(l) = (eta(l+1)*hcld(l+1)+(eta(l)-eta(l+1))*hol(l))/eta(l)
      tem = (hcld(l)-hst(l))*(zet(l)-zet(l+1))/(1.0+lbcp*dqq(l))
      bke(l) = bke(l+1) + cons_grav*tem/(cons_cp*prj(l+1)*poi(l))
      IF (tem .LT. 0.0) THEN
        max1 = 0.0
      ELSE
        max1 = tem
      END IF
      bk2(l) = bk2(l+1) + cons_grav*max1/(cons_cp*prj(l+1)*poi(l))
      IF (bk2(l) .LT. 0.0) THEN
        max2 = 0.0
      ELSE
        max2 = bk2(l)
      END IF
      cvw(l) = SQRT(2.0*max2)
    END DO
!ALPHA CALCULATION 
    IF (zet(ic) .LT. 2000.) rasal(ic) = rasal1
    IF (zet(ic) .GE. 2000.) rasal(ic) = rasal1 + (rasal2-rasal1)*(zet(ic&
&        )-2000.)/8000.
    IF (rasal(ic) .GT. 1.0e5) THEN
      rasal(ic) = 1.0e5
    ELSE
      rasal(ic) = rasal(ic)
    END IF
    rasal(ic) = dt/rasal(ic)
    DO kk=ic,k
      IF (cvw(kk) .LT. 1.00) THEN
        cvw(kk) = 1.00
      ELSE
        cvw(kk) = cvw(kk)
      END IF
    END DO
!TEST FOR CRITICAL WORK FUNCTION
    CALL ACRITN(pol(ic), prs(k), acr, acritfac)
    IF (wfn .LE. acr) rc(ic) = 4
! SUB-CRITICAL WORK FUNCTION ======>>
    IF (rc(ic) .EQ. 0) THEN
      IF (do_tracers) THEN
        DO itr=1,itrcr
!Scavenging of the below cloud tracer
          delzkm = (zet(ic)-zet(k))/1000.
          x4 = EXP(-(fscav_(itr)*delzkm))
          IF (x4 .GT. 1.) THEN
            x1 = 1.
          ELSE
            x1 = x4
          END IF
          IF (x1 .LT. 0.) THEN
            fnoscav = 0.
          ELSE
            fnoscav = x1
          END IF
          xht(itr) = xoi(k, itr)*fnoscav
        END DO
      END IF
      wlq = qol(k)
      uht = uoi(k)
      vht = voi(k)
      rnn(k) = 0.
      cll0(k) = 0.
      DO l=k-1,ic,-1
        tem = eta(l) - eta(l+1)
        wlq = wlq + tem*qol(l)
        uht = uht + tem*uoi(l)
        vht = vht + tem*voi(l)
        IF (do_tracers) THEN
          DO itr=1,itrcr
!Scavenging of the entrained tracer.  Updates transported tracer mass.
            delzkm = (zet(ic)-zet(l+1))/1000.
            x5 = EXP(-(fscav_(itr)*delzkm))
            IF (x5 .GT. 1.) THEN
              x2 = 1.
            ELSE
              x2 = x5
            END IF
            IF (x2 .LT. 0.) THEN
              fnoscav = 0.
            ELSE
              fnoscav = x2
            END IF
            xht(itr) = xht(itr) + tem*xoi(l, itr)*fnoscav
          END DO
        END IF
!How much condensate (CLI) is present here? 
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
! condensate (kg/kg)
        cli = cll0(l)/eta(l)
! Temperature (K)
        te_a = poi(l)*prh(l)
        CALL SUNDQ3_ICE(te_a, sdqv2, sdqv3, sdqvt1, f2, f3)
!F4 reduces AUTO for shallow clouds, F5 modifies auto for deep clouds
!* F5
        c00_x = co_auto*f2*f3*f4
        cli_crit_x = cli_crit/(f2*f3)
        rate = c00_x*(1.0-EXP(-(cli**2/cli_crit_x**2)))
        IF (cvw(l) .LT. 1.00) THEN
          cvw_x = 1.00
        ELSE
          cvw_x = cvw(l)
        END IF
! l.h.s. DT_LYR => time in layer (L,L+1)
        dt_lyr = (zet(l)-zet(l+1))/cvw_x
        closs = cll0(l)*rate*dt_lyr
        IF (closs .GT. cll0(l)) THEN
          closs = cll0(l)
        ELSE
          closs = closs
        END IF
        cll0(l) = cll0(l) - closs
        IF (closs .GT. 0.) THEN
          wlq = wlq - closs
          rnn(l) = closs
        ELSE
          rnn(l) = 0.
        END IF
      END DO
      wlq = wlq - qst(ic)*eta(ic)
!CALCULATE GAMMAS AND KERNEL
! MS-A30 (W/O GRAV)
      gms(k) = (sht(k)-ssl(k))*pri(k)
! MS-A31 (W/O GRAV)
      gmh(k) = gms(k) + (qht(k)-qol(k))*pri(k)*cons_alhl
! MS-A37 (W/O GRAV)
      akm = gmh(k)*gam(k-1)*dpb(k-1)
      tx2 = gmh(k)
      DO l=k-1,ic+1,-1
        gms(l) = (eta(l)*(sht(l)-ssl(l))+eta(l+1)*(ssl(l)-sht(l+1)))*pri&
&          (l)
        gmh(l) = gms(l) + (eta(l)*(qht(l)-qol(l))+eta(l+1)*(qol(l)-qht(l&
&          +1)))*cons_alhl*pri(l)
        tx2 = tx2 + (eta(l)-eta(l+1))*gmh(l)
        akm = akm - gms(l)*eht(l)*pki(l) + tx2*ght(l)
      END DO
      gms(ic) = eta(ic+1)*(ssl(ic)-sht(ic+1))*pri(ic)
      akm = akm - gms(ic)*eta(ic+1)*dpb(ic)*pki(ic)
      gmh(ic) = gms(ic) + (eta(ic+1)*(qol(ic)-qht(ic+1))*cons_alhl+eta(&
&        ic)*(hst(ic)-hol(ic)))*pri(ic)
!CLOUD BASE MASS FLUX
      IF (akm .GE. 0.0 .OR. wlq .LT. 0.0) rc(ic) = 5
!  =========>
      IF (rc(ic) .EQ. 0) THEN
! MS-A39 MASS-FLUX IN Pa/step
        wfn = -((wfn-acr)/akm)
        x3 = rasal(ic)*trg*toki*wfn
        IF (x3 .GT. (prs(k+1)-prs(k))*(100.*pblfrac)) THEN
          wfn = (prs(k+1)-prs(k))*(100.*pblfrac)
        ELSE
          wfn = x3
        END IF
!CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
        IF (PRESENT(disske)) wfnog = wfn*gravi
        tem = wfn*gravi
! (kg/m^2/step)
        cll(ic) = cll(ic) + wlq*tem
! (kg/m^2/step)
        rmf(ic) = rmf(ic) + tem
! (kg/m^2/step)
        rmfd(ic) = rmfd(ic) + tem*eta(ic)
        DO l=ic+1,k
! (kg/m^2/step)
          rmfp(l) = tem*eta(l)
! (kg/m^2/step)
          rmfc(l) = rmfc(l) + rmfp(l)
          IF (cvw(l) .GT. 0.0) THEN
            updfrp(l) = rmfp(l)*(ddt/daylen)*1000./(cvw(l)*prs(l))
          ELSE
            updfrp(l) = 0.0
          END IF
! current cloud; incloud condensate        
          clli(l) = cll0(l)/eta(l)
!  cumulative grid mean convective condensate        
          cllb(l) = cllb(l) + updfrp(l)*clli(l)
          updfrc(l) = updfrc(l) + updfrp(l)
        END DO
!THETA AND Q CHANGE DUE TO CLOUD TYPE IC
        DO l=ic,k
! (kg/m^2/step)
          rns(l) = rns(l) + rnn(l)*tem
          gmh(l) = gmh(l)*wfn
          gms(l) = gms(l)*wfn
          qoi(l) = qoi(l) + (gmh(l)-gms(l))*alhi
          poi(l) = poi(l) + gms(l)*pki(l)*cpi
          qst(l) = qst(l) + gms(l)*bet(l)*cpi
        END DO
        IF (do_tracers) THEN
!*FRICFAC*0.5
          wfn = wfn*0.5*1.0
          tem = wfn*pri(k)
          xcu(icmin:, :) = 0.
          DO itr=1,itrcr
            xcu(k, itr) = xcu(k, itr) + tem*(xoi(k-1, itr)-xoi(k, itr))
          END DO
          DO itr=1,itrcr
            DO l=k-1,ic+1,-1
              tem = wfn*pri(l)
              xcu(l, itr) = xcu(l, itr) + tem*((xoi(l-1, itr)-xoi(l, itr&
&                ))*eta(l)+(xoi(l, itr)-xoi(l+1, itr))*eta(l+1))
            END DO
          END DO
          tem = wfn*pri(ic)
          DO itr=1,itrcr
            xcu(ic, itr) = xcu(ic, itr) + (2.*(xht(itr)-xoi(ic, itr)*(&
&              eta(ic)-eta(ic+1)))-(xoi(ic, itr)+xoi(ic+1, itr))*eta(ic+1&
&              ))*tem
          END DO
          DO itr=1,itrcr
            DO l=ic,k
              xoi(l, itr) = xoi(l, itr) + xcu(l, itr)
            END DO
          END DO
        ELSE
!*FRICFAC*0.5
          wfn = wfn*0.5*1.0
        END IF
!Cumulus friction
        wfn = wfn*fricfac*EXP(-(alm/friclambda))
        tem = wfn*pri(k)
! This change makes cumulus friction
        ucu(icmin:) = 0.
! correct.
        vcu(icmin:) = 0.
        ucu(k) = ucu(k) + tem*(uoi(k-1)-uoi(k))
        vcu(k) = vcu(k) + tem*(voi(k-1)-voi(k))
        DO l=k-1,ic+1,-1
          tem = wfn*pri(l)
          ucu(l) = ucu(l) + tem*((uoi(l-1)-uoi(l))*eta(l)+(uoi(l)-uoi(l+&
&            1))*eta(l+1))
          vcu(l) = vcu(l) + tem*((voi(l-1)-voi(l))*eta(l)+(voi(l)-voi(l+&
&            1))*eta(l+1))
        END DO
        tem = wfn*pri(ic)
        ucu(ic) = ucu(ic) + (2.*(uht-uoi(ic)*(eta(ic)-eta(ic+1)))-(uoi(&
&          ic)+uoi(ic+1))*eta(ic+1))*tem
        vcu(ic) = vcu(ic) + (2.*(vht-voi(ic)*(eta(ic)-eta(ic+1)))-(voi(&
&          ic)+voi(ic+1))*eta(ic+1))*tem
        IF (PRESENT(disske)) dissk0(ic) = eta(ic)*cons_grav*wfnog*pri(ic&
&            )*0.5*((uht/eta(ic)-uoi(ic))**2+(vht/eta(ic)-voi(ic))**2)
        DO l=ic,k
          uoi(l) = uoi(l) + ucu(l)
          voi(l) = voi(l) + vcu(l)
        END DO
        rc(ic) = 10
      END IF
    END IF
  END DO
!---------------------
!MAIN CLOUD SUBROUTINE
  IF (SUM(rmf(icmin:k)) .GT. 0.0) THEN
!STRAP (FINAL = 1) SUBROUTINE
!----------------------------
!Scale properly by layer masses
    wght0 = 0.
    DO l=k,k0
      wght0(l) = wght0(l-1) + wgt1(l)*(ple(l+1)-ple(l))
    END DO
    wght1 = wgt1*(prs(k+1)-prs(k))/wght0(k0)
    tho(icmin:k-1) = poi(icmin:k-1)
    qho(icmin:k-1) = qoi(icmin:k-1)
    uho(icmin:k-1) = uoi(icmin:k-1)
    vho(icmin:k-1) = voi(icmin:k-1)
    DO l=k,k0
      tho(l) = tho(l) + wght1(l)*(poi(k)-poi_sv(k))
      qho(l) = qho(l) + wght1(l)*(qoi(k)-qoi_sv(k))
      uho(l) = uho(l) + wght1(l)*(uoi(k)-uoi_sv(k))
      vho(l) = vho(l) + wght1(l)*(voi(k)-voi_sv(k))
    END DO
!Outputs -> Cloud
!  (KG/m^2/s )
    clw(icmin:k) = cll(icmin:k)*ddt/daylen + clwtmp(icmin:k)
!  (KG/m^2/s @ CLOUD TOP)
    flxd(icmin:k) = rmfd(icmin:k)*ddt/daylen + flxdtmp(icmin:k)
    DO l=icmin,k
      tem = pri(l)*cons_grav
      cnv_prc3(l) = rns(l)*tem + cnv_prc3tmp(l)
    END DO
    cnv_updfrc(icmin:k-1) = updfrc(icmin:k-1) + cnv_updfrctmp(icmin:k-1)
    IF (k .LT. k0) THEN
      clw(k:k0) = 0.
      flxd(k:k0) = 0.
    END IF
!Other outputs
!    FLX (ICMIN:K) = RMF (ICMIN:K) * DDT/DAYLEN + FLXTMP (ICMIN:K) !  (KG/m^2/s @ CLOUD BASE)
!    FLXC(ICMIN:K) = RMFC(ICMIN:K) * DDT/DAYLEN + FLXCTMP(ICMIN:K) !  (KG/m^2/s @ CLOUD TOP)
!    CNV_CVW   (ICMIN:K-1)   =  CVW(ICMIN:K-1)    + CNV_CVWTMP   (ICMIN:K-1)
!    CNV_QC    (ICMIN:K-1)   =  CLLB(ICMIN:K-1)   + CNV_QCTMP    (ICMIN:K-1)
!    IF ( K < K0 ) THEN 
!       FLX (K:K0) = 0.
!       FLXC(K:K0) = 0.
!    END IF
    IF (do_tracers) THEN
      xho(icmin:k-1, :) = xoi(icmin:k-1, :)
      DO itr=1,itrcr
        DO l=k,k0
          xho(l, itr) = xho(l, itr) + wght1(l)*(xoi(k, itr)-xoi_sv(k, &
&            itr))
        END DO
      END DO
    END IF
    IF (PRESENT(disske)) disske(icmin:k-1) = dissk0(icmin:k-1)*ddt/&
&        daylen + dissketmp(icmin:k-1)
  END IF
END SUBROUTINE RASE

! Subroutines
!------------
SUBROUTINE ACRITN(pl, plb, acr, acritfac)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: pl, plb, acritfac
  REAL*8, INTENT(OUT) :: acr
  INTEGER :: iwk
  REAL*8, PARAMETER :: ph(15)=(/150.0, 200.0, 250.0, 300.0, 350.0, 400.0&
&    , 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/)
  REAL*8, PARAMETER :: a(15)=(/1.6851, 1.1686, 0.7663, 0.5255, 0.4100, &
&    0.3677, 0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664, 0.0553, &
&    0.0445, 0.0633/)
  INTRINSIC INT
  iwk = INT(pl*0.02 - 0.999999999)
  IF (iwk .GT. 1 .AND. iwk .LE. 15) THEN
    acr = a(iwk-1) + (pl-ph(iwk-1))*.02*(a(iwk)-a(iwk-1))
  ELSE IF (iwk .GT. 15) THEN
    acr = a(15)
  ELSE
    acr = a(1)
  END IF
  acr = acritfac*acr*(plb-pl)
END SUBROUTINE ACRITN

SUBROUTINE SUNDQ3_ICE(temp, rate2, rate3, te1, f2, f3)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: temp, rate2, rate3, te1
  REAL*8, INTENT(OUT) :: f2, f3
!,RATE2,RATE3,TE1
  REAL*8 :: xx, yy, te0, te2, jump1
  te0 = 273.
  te2 = 200.
  jump1 = (rate2-1.0)/(te0-te1)**0.333
! Ice - phase treatment
  IF (temp .GE. te0) THEN
    f2 = 1.0
    f3 = 1.0
  END IF
  IF (temp .GE. te1 .AND. temp .LT. te0) THEN
    f2 = 1.0 + jump1*(te0-temp)**0.3333
    f3 = 1.0
  END IF
  IF (temp .LT. te1) THEN
    f2 = rate2 + (rate3-rate2)*(te1-temp)/(te1-te2)
    f3 = 1.0
  END IF
  IF (f2 .GT. 27.0) f2 = 27.0
END SUBROUTINE SUNDQ3_ICE

SUBROUTINE DQSAT_RAS(dqsi, qssi, temp, plo, lm, estblx, cons_h2omw, &
&  cons_airmw)
  IMPLICIT NONE
!Inputs
  INTEGER :: lm
  REAL*8, DIMENSION(lm) :: temp, plo
  REAL*8 :: estblx(:)
  REAL*4 :: cons_h2omw, cons_airmw
!Outputs
  REAL*8, DIMENSION(lm) :: dqsi, qssi
!Locals
  REAL*8, PARAMETER :: max_mixing_ratio=1.0
  REAL*8 :: esfac
  INTEGER :: k
  REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
  INTEGER :: it
  INTEGER, PARAMETER :: degsubs=100
  REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
  INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
  INTRINSIC NINT
  INTRINSIC INT
  esfac = cons_h2omw/cons_airmw
  DO k=1,lm
    tl = temp(k)
    pl = plo(k)
    pp = pl*100.0
    IF (tl .LE. tmintbl) THEN
      ti = tmintbl
    ELSE IF (tl .GE. tmaxtbl - .001) THEN
      ti = tmaxtbl - .001
    ELSE
      ti = tl
    END IF
    tt = (ti-tmintbl)*degsubs + 1
    it = INT(tt)
    dqq = estblx(it+1) - estblx(it)
    qq = (tt-it)*dqq + estblx(it)
    IF (pp .LE. qq) THEN
      qsat = max_mixing_ratio
      dqsat = 0.0
    ELSE
      dd = 1.0/(pp-(1.0-esfac)*qq)
      qsat = esfac*qq*dd
      dqsat = esfac*degsubs*dqq*pp*(dd*dd)
    END IF
    dqsi(k) = dqsat
    qssi(k) = qsat
  END DO
END SUBROUTINE DQSAT_RAS

SUBROUTINE DQSATS_RAS(dqsi, qssi, temp, plo, estblx, cons_h2omw, &
&  cons_airmw)
  IMPLICIT NONE
!Inputs
  REAL*8 :: temp, plo
  REAL*8 :: estblx(:)
  REAL*4 :: cons_h2omw, cons_airmw
!Outputs
  REAL*8 :: dqsi, qssi
!Locals
  REAL*8, PARAMETER :: max_mixing_ratio=1.0
  REAL*8 :: esfac
  REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
  INTEGER :: it
  INTEGER, PARAMETER :: degsubs=100
  REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
  INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
  INTRINSIC NINT
  INTRINSIC INT
  esfac = cons_h2omw/cons_airmw
  tl = temp
  pl = plo
  pp = pl*100.0
  IF (tl .LE. tmintbl) THEN
    ti = tmintbl
  ELSE IF (tl .GE. tmaxtbl - .001) THEN
    ti = tmaxtbl - .001
  ELSE
    ti = tl
  END IF
  tt = (ti-tmintbl)*degsubs + 1
  it = INT(tt)
  dqq = estblx(it+1) - estblx(it)
  qq = (tt-it)*dqq + estblx(it)
  IF (pp .LE. qq) THEN
    qsat = max_mixing_ratio
    dqsat = 0.0
  ELSE
    dd = 1.0/(pp-(1.0-esfac)*qq)
    qsat = esfac*qq*dd
    dqsat = esfac*degsubs*dqq*pp*(dd*dd)
  END IF
  dqsi = dqsat
  qssi = qsat
END SUBROUTINE DQSATS_RAS

