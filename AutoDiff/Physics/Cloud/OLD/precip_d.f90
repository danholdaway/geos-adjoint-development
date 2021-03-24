!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.7 (r4786) - 21 Feb 2013 15:53
!
!  Differentiation of precipandevap in forward (tangent) mode:
!   variations   of useful results: evap_dd_above_out qv subl_dd_above_out
!                qcl qpi pfi_above_out qpl pfl_above_out te
!   with respect to varying inputs: aa area pfl_above_in mass imass
!                qv pfi_above_in bb evap_dd_above_in qcl qpi qpl
!                subl_dd_above_in pl rhcr3 dze qddf3 te
!   RW status of diff variables: evap_dd_above_out:out aa:in area:in
!                pfl_above_in:in mass:in imass:in qv:in-out pfi_above_in:in
!                bb:in subl_dd_above_out:out evap_dd_above_in:in
!                qcl:in-out qpi:in-zero pfi_above_out:out qpl:in-zero
!                subl_dd_above_in:in pl:in rhcr3:in dze:in qddf3:in
!                pfl_above_out:out te:in-out
SUBROUTINE PRECIPANDEVAP_D(k, lm, dt, frland, rhcr3, rhcr3d, qpl, qpld, &
&  qpi, qpid, qcl, qcld, qci, te, ted, qv, qvd, mass, massd, imass, &
&  imassd, pl, pld, dze, dzed, qddf3, qddf3d, aa, aad, bb, bbd, area, &
&  aread, pfl_above_in, pfl_above_ind, pfl_above_out, pfl_above_outd, &
&  pfi_above_in, pfi_above_ind, pfi_above_out, pfi_above_outd, &
&  evap_dd_above_in, evap_dd_above_ind, evap_dd_above_out, &
&  evap_dd_above_outd, subl_dd_above_in, subl_dd_above_ind, &
&  subl_dd_above_out, subl_dd_above_outd, frz_diag, envfc, ddrfc, &
&  cons_alhf, cons_alhs, cons_alhl, cons_cp, cons_tice, revap_off_p, &
&  c_acc, c_ev_r, c_ev_s, rho_w, estblx)
  IMPLICIT NONE
!REVAP_DIAG = REVAP_DIAG + EVAP / DT
!RSUBL_DIAG = RSUBL_DIAG + SUBL / DT
!Inputs
  INTEGER, INTENT(IN) :: k, lm
  REAL*8, INTENT(IN) :: dt, mass, imass, pl, aa, bb, rhcr3, dze, qddf3, &
&  area, frland, envfc, ddrfc
  REAL*8, INTENT(IN) :: massd, imassd, pld, aad, bbd, rhcr3d, dzed, &
&  qddf3d, aread
  REAL*8, INTENT(IN) :: cons_alhf, cons_alhs, cons_alhl, cons_cp, &
&  cons_tice, revap_off_p
  REAL*8, INTENT(IN) :: c_acc, c_ev_r, c_ev_s, rho_w
  REAL*8, INTENT(IN) :: estblx(50)
!Prognostics
  REAL*8, INTENT(INOUT) :: qv, qpl, qpi, qcl, qci, te, frz_diag
  REAL*8, INTENT(INOUT) :: qvd, qpld, qpid, qcld, ted
  REAL*8, INTENT(INOUT) :: pfl_above_in, pfl_above_out, pfi_above_in, &
&  pfi_above_out
  REAL*8, INTENT(INOUT) :: pfl_above_ind, pfl_above_outd, pfi_above_ind&
&  , pfi_above_outd
  REAL*8, INTENT(INOUT) :: evap_dd_above_in, evap_dd_above_out, &
&  subl_dd_above_in, subl_dd_above_out
  REAL*8, INTENT(INOUT) :: evap_dd_above_ind, evap_dd_above_outd, &
&  subl_dd_above_ind, subl_dd_above_outd
!Outputs - not used for now
!real(8), intent(  out) :: RAIN, SNOW, REVAP_DIAG, RSUBL_DIAG, ACRLL_DIAG, ACRIL_DIAG, PFL_DIAG, PFI_DIAG, VFALLSN, VFALLRN
!Locals
  INTEGER :: ns, nsmx, itr, l
  REAL*8 :: pfi, pfl, qs, dqs, envfrac, tko, qko, qstko, dqstko, rh_box&
&  , t_ed, qplko, qpiko
  REAL*8 :: pfid, pfld, qsd, dqsd, tkod, qkod, qstkod, dqstkod, rh_boxd&
&  , t_edd
  REAL*8 :: ifactor, rainrat0, snowrat0, fallrn, fallsn, vesn, vern, &
&  nrain, nsnow, efactor
  REAL*8 :: ifactord, rainrat0d, snowrat0d, fallrnd, fallsnd, vesnd, &
&  vernd, efactord
  REAL*8 :: tinlayerrn, diamrn, droprad, tinlayersn, diamsn, flakrad
  REAL*8 :: tinlayerrnd, diamrnd, dropradd, tinlayersnd, diamsnd, &
&  flakradd
  REAL*8 :: evap, subl, accr, mltfrz, evapx, sublx, evap_dd, subl_dd, &
&  ddfract, landseaf
  REAL*8 :: evapd, subld, accrd, mltfrzd, evapxd, sublxd, evap_ddd, &
&  subl_ddd
  REAL*8 :: tau_frz, tau_mlt
!m/s
  REAL*8, PARAMETER :: trmv_l=1.0
  LOGICAL, PARAMETER :: taneff=.false.
!Fraction of precip falling through "environment" vs through cloud
  REAL*8, PARAMETER :: b_sub=1.00
  REAL*8 :: arg1
  REAL*8 :: arg1d
  INTRINSIC EXP
  INTRINSIC MAX
  INTRINSIC MIN
!Initialize diagnostics
!rain = 0.0
!snow = 0.0
!revap_diag = 0.0
!rsubl_diag = 0.0
!acrll_diag = 0.0
!acril_diag = 0.0
!pfl_diag = 0.0
!pfi_diag = 0.0
!vfallsn = 0.0
!vfallrn = 0.0
  IF (.NOT.taneff) envfrac = envfc
!envfrac = 1.00
!if (pl .le. 600.) then
!   envfrac = 0.25
!else
!   envfrac = 0.25 + (1.-0.25)/(19.) * ((atan( (2.*(pl-600.)/(900.-600.)-1.) *       &
!                                    tan(20.*CONS_PI/21.-0.5*CONS_PI) ) + 0.5*CONS_PI) * 21./CONS_PI - 1.)
!end if
!
!envfrac = min(envfrac,1.)
  IF (area .GT. 0.) THEN
    ifactord = -(aread/area**2)
    ifactor = 1./area
  ELSE
    ifactor = 1.00
    ifactord = 0.0_8
  END IF
  IF (ifactor .LT. 1.) THEN
    ifactor = 1.
    ifactord = 0.0_8
  ELSE
    ifactor = ifactor
  END IF
!Start at top of precip column:
!
!   a) Accrete                   
!   b) Evaporate/Sublimate  
!   c) Rain/Snow-out to next level down 
!   d) return to (a)
!Update saturated humidity
  CALL DQSAT_SUB_SCA_D(dqs, dqsd, qs, qsd, te, ted, pl, pld, estblx)
  ddfract = ddrfc
  IF (k .EQ. 1) THEN
    pfld = qpld*mass + qpl*massd
    pfl = qpl*mass
    pfid = qpid*mass + qpi*massd
    pfi = qpi*mass
    evap_dd = 0.
    subl_dd = 0.
!VFALLRN = 0.0
!VFALLSN = 0.0
    evap_ddd = 0.0_8
    subl_ddd = 0.0_8
  ELSE
    qpld = qpld + pfl_above_ind*imass + pfl_above_in*imassd
    qpl = qpl + pfl_above_in*imass
    pfl = 0.00
    qpid = qpid + pfi_above_ind*imass + pfi_above_in*imassd
    qpi = qpi + pfi_above_in*imass
    pfi = 0.00
    accrd = b_sub*c_acc*((qpld*mass+qpl*massd)*qcl+qpl*mass*qcld)
    accr = b_sub*c_acc*(qpl*mass)*qcl
    IF (accr .GT. qcl) THEN
      accrd = qcld
      accr = qcl
    ELSE
      accr = accr
    END IF
    qpld = qpld + accrd
    qpl = qpl + accr
    qcld = qcld - accrd
    qcl = qcl - accr
!ACRLL_DIAG = ACCR / DT
!Accretion of liquid condensate by falling ice/snow
    accrd = b_sub*c_acc*((qpid*mass+qpi*massd)*qcl+qpi*mass*qcld)
    accr = b_sub*c_acc*(qpi*mass)*qcl
    IF (accr .GT. qcl) THEN
      accrd = qcld
      accr = qcl
    ELSE
      accr = accr
    END IF
    qpid = qpid + accrd
    qpi = qpi + accr
    qcld = qcld - accrd
    qcl = qcl - accr
!! Liquid freezes when accreted by snow
    ted = ted + cons_alhf*accrd/cons_cp
    te = te + cons_alhf*accr/cons_cp
!ACRIL_DIAG = ACCR / DT
    rainrat0d = ((ifactord*qpl+ifactor*qpld)*mass+ifactor*qpl*massd)/dt
    rainrat0 = ifactor*qpl*mass/dt
    snowrat0d = ((ifactord*qpi+ifactor*qpid)*mass+ifactor*qpi*massd)/dt
    snowrat0 = ifactor*qpi*mass/dt
    CALL MARSHPALM_D(rainrat0, rainrat0d, pl, pld, diamrn, diamrnd, &
&               nrain, fallrn, fallrnd, vern, vernd)
    CALL MARSHPALM_D(snowrat0, snowrat0d, pl, pld, diamsn, diamsnd, &
&               nsnow, fallsn, fallsnd, vesn, vesnd)
!VFALLRN = FALLrn
!VFALLSN = FALLsn
    tinlayerrnd = (dzed*(fallrn+0.01)-dze*fallrnd)/(fallrn+0.01)**2
    tinlayerrn = dze/(fallrn+0.01)
    tinlayersnd = (dzed*(fallsn+0.01)-dze*fallsnd)/(fallsn+0.01)**2
    tinlayersn = dze/(fallsn+0.01)
!Melting of Frozen precipitation      
! time scale for freezing (s). 
    tau_frz = 5000.
    mltfrz = 0.0
    IF (te .GT. cons_tice .AND. te .LE. cons_tice + 5.) THEN
      mltfrzd = ((tinlayersnd*qpi+tinlayersn*qpid)*(te-cons_tice)+&
&        tinlayersn*qpi*ted)/tau_frz
      mltfrz = tinlayersn*qpi*(te-cons_tice)/tau_frz
      IF (qpi .GT. mltfrz) THEN
        mltfrz = mltfrz
      ELSE
        mltfrzd = qpid
        mltfrz = qpi
      END IF
      ted = ted - cons_alhf*mltfrzd/cons_cp
      te = te - cons_alhf*mltfrz/cons_cp
      qpld = qpld + mltfrzd
      qpl = qpl + mltfrz
      qpid = qpid - mltfrzd
      qpi = qpi - mltfrz
    END IF
    frz_diag = frz_diag - mltfrz/dt
    mltfrz = 0.0
    IF (te .GT. cons_tice + 5.) THEN
! Go Ahead and melt any snow/hail left above 5 C 
      mltfrzd = qpid
      mltfrz = qpi
      ted = ted - cons_alhf*mltfrzd/cons_cp
      te = te - cons_alhf*mltfrz/cons_cp
      qpld = qpld + mltfrzd
      qpl = qpl + mltfrz
      qpid = qpid - mltfrzd
      qpi = qpi - mltfrz
    END IF
    frz_diag = frz_diag - mltfrz/dt
    mltfrz = 0.0
    IF (k .GE. lm - 1) THEN
      IF (te .GT. cons_tice + 0.) THEN
! Go Ahead and melt any snow/hail left above 0 C in lowest layers 
        mltfrzd = qpid
        mltfrz = qpi
        ted = ted - cons_alhf*mltfrzd/cons_cp
        te = te - cons_alhf*mltfrz/cons_cp
        qpld = qpld + mltfrzd
        qpl = qpl + mltfrz
        qpid = qpid - mltfrzd
        qpi = qpi - mltfrz
      END IF
    END IF
    frz_diag = frz_diag - mltfrz/dt
!Freezing of liquid precipitation      
    mltfrz = 0.0
    IF (te .LE. cons_tice) THEN
      ted = ted + cons_alhf*qpld/cons_cp
      te = te + cons_alhf*qpl/cons_cp
      qpid = qpld + qpid
      qpi = qpl + qpi
      mltfrz = qpl
      qpl = 0.
      qpld = 0.0_8
    END IF
    frz_diag = frz_diag + mltfrz/dt
!In the exp below, evaporation time scale is determined "microphysically" from temp, 
!press, and drop size. In this context C_EV becomes a dimensionless fudge-fraction. 
!Also remember that these microphysics are still only for liquid.
    qkod = qvd
    qko = qv
    tkod = ted
    tko = te
    qplko = qpl
    qpiko = qpi
    sublxd = 0.0_8
    evapd = 0.0_8
    subld = 0.0_8
    evapxd = 0.0_8
    DO itr=1,3
      dqstkod = dqsd
      dqstko = dqs
      qstkod = qsd + dqstkod*(tko-te) + dqstko*(tkod-ted)
      qstko = qs + dqstko*(tko-te)
      IF (qstko .LT. 1.0e-7) THEN
        qstko = 1.0e-7
        qstkod = 0.0_8
      ELSE
        qstko = qstko
      END IF
      rh_boxd = (qkod*qstko-qko*qstkod)/qstko**2
      rh_box = qko/qstko
      qko = qv
      tko = te
      IF (rh_box .LT. rhcr3) THEN
        efactord = (rho_w*(aad+bbd)*(rhcr3-rh_box)-rho_w*(aa+bb)*(rhcr3d&
&          -rh_boxd))/(rhcr3-rh_box)**2
        efactor = rho_w*(aa+bb)/(rhcr3-rh_box)
      ELSE
        efactor = 9.99e9
        efactord = 0.0_8
      END IF
      IF (frland .LT. 0.1) THEN
! Over Ocean
        landseaf = 0.5
      ELSE
! Over Land
        landseaf = 0.5
      END IF
      landseaf = 1.00
!Rain falling
      IF (rh_box .LT. rhcr3 .AND. diamrn .GT. 0.00 .AND. pl .GT. 100. &
&          .AND. pl .LT. revap_off_p) THEN
        dropradd = 0.5*diamrnd
        droprad = 0.5*diamrn
        t_edd = efactord*droprad**2 + efactor*2*droprad*dropradd
        t_ed = efactor*droprad**2
        t_edd = t_edd*(1.0+dqstko*cons_alhl/cons_cp) + t_ed*cons_alhl*&
&          dqstkod/cons_cp
        t_ed = t_ed*(1.0+dqstko*cons_alhl/cons_cp)
        arg1d = -((c_ev_r*landseaf*envfrac*(vernd*tinlayerrn+vern*&
&          tinlayerrnd)*t_ed-c_ev_r*vern*landseaf*envfrac*tinlayerrn*&
&          t_edd)/t_ed**2)
        arg1 = -(c_ev_r*vern*landseaf*envfrac*tinlayerrn/t_ed)
        evapd = qpld*(1.0-EXP(arg1)) - qpl*arg1d*EXP(arg1)
        evap = qpl*(1.0-EXP(arg1))
      ELSE
        evap = 0.0
        evapd = 0.0_8
      END IF
!Snow falling
      IF (rh_box .LT. rhcr3 .AND. diamsn .GT. 0.00 .AND. pl .GT. 100. &
&          .AND. pl .LT. revap_off_p) THEN
        flakradd = 0.5*diamsnd
        flakrad = 0.5*diamsn
        t_edd = efactord*flakrad**2 + efactor*2*flakrad*flakradd
        t_ed = efactor*flakrad**2
        t_edd = t_edd*(1.0+dqstko*cons_alhs/cons_cp) + t_ed*cons_alhs*&
&          dqstkod/cons_cp
        t_ed = t_ed*(1.0+dqstko*cons_alhs/cons_cp)
        arg1d = -((c_ev_s*landseaf*envfrac*(vesnd*tinlayersn+vesn*&
&          tinlayersnd)*t_ed-c_ev_s*vesn*landseaf*envfrac*tinlayersn*&
&          t_edd)/t_ed**2)
        arg1 = -(c_ev_s*vesn*landseaf*envfrac*tinlayersn/t_ed)
        subld = qpid*(1.0-EXP(arg1)) - qpi*arg1d*EXP(arg1)
        subl = qpi*(1.0-EXP(arg1))
      ELSE
        subl = 0.0
        subld = 0.0_8
      END IF
      IF (itr .EQ. 1) THEN
        evapxd = evapd
        evapx = evap
        sublxd = subld
        sublx = subl
      ELSE
        evapd = (evapd+evapxd)/2.0
        evap = (evap+evapx)/2.0
        subld = (subld+sublxd)/2.0
        subl = (subl+sublx)/2.0
      END IF
      qkod = qvd + evapd + subld
      qko = qv + evap + subl
      tkod = ted - cons_alhl*evapd/cons_cp - cons_alhs*subld/cons_cp
      tko = te - evap*cons_alhl/cons_cp - subl*cons_alhs/cons_cp
    END DO
    qpid = qpid - subld
    qpi = qpi - subl
    qpld = qpld - evapd
    qpl = qpl - evap
!Put some re-evap/re-subl precip in to a \quote{downdraft} to be applied later
    evap_ddd = evap_dd_above_ind + ddfract*(evapd*mass+evap*massd)
    evap_dd = evap_dd_above_in + ddfract*evap*mass
    evapd = evapd - ddfract*evapd
    evap = evap - ddfract*evap
    subl_ddd = subl_dd_above_ind + ddfract*(subld*mass+subl*massd)
    subl_dd = subl_dd_above_in + ddfract*subl*mass
    subld = subld - ddfract*subld
    subl = subl - ddfract*subl
    qvd = qvd + evapd + subld
    qv = qv + evap + subl
    ted = ted - cons_alhl*evapd/cons_cp - cons_alhs*subld/cons_cp
    te = te - evap*cons_alhl/cons_cp - subl*cons_alhs/cons_cp
!REVAP_DIAG = EVAP / DT
!RSUBL_DIAG = SUBL / DT
    pfld = qpld*mass + qpl*massd
    pfl = qpl*mass
    pfid = qpid*mass + qpi*massd
    pfi = qpi*mass
!PFL_DIAG =  PFl/DT
!PFI_DIAG =  PFi/DT
  END IF
  evapd = ((qddf3d*evap_dd+qddf3*evap_ddd)*mass-qddf3*evap_dd*massd)/&
&    mass**2
  evap = qddf3*evap_dd/mass
  subld = ((qddf3d*subl_dd+qddf3*subl_ddd)*mass-qddf3*subl_dd*massd)/&
&    mass**2
  subl = qddf3*subl_dd/mass
  qvd = qvd + evapd + subld
  qv = qv + evap + subl
  ted = ted - cons_alhl*evapd/cons_cp - cons_alhs*subld/cons_cp
  te = te - evap*cons_alhl/cons_cp - subl*cons_alhs/cons_cp
!RAIN  = PFl/DT
!SNOW  = PFi/DT
  qpi = 0.
  qpl = 0.
  pfl_above_outd = pfld
  pfl_above_out = pfl
  pfi_above_outd = pfid
  pfi_above_out = pfi
  evap_dd_above_outd = evap_ddd
  evap_dd_above_out = evap_dd
  subl_dd_above_outd = subl_ddd
  subl_dd_above_out = subl_dd
  qpid = 0.0_8
  qpld = 0.0_8
END SUBROUTINE PRECIPANDEVAP_D

!  Differentiation of dqsat_sub_sca in forward (tangent) mode:
!   variations   of useful results: dqsi qssi
!   with respect to varying inputs: temp plo
SUBROUTINE DQSAT_SUB_SCA_D(dqsi, dqsid, qssi, qssid, temp, tempd, plo, &
&  plod, estblx)
  IMPLICIT NONE
!Inputs
  REAL*8 :: temp, plo
  REAL*8 :: tempd, plod
  REAL*8 :: estblx(:)
!Outputs
  REAL*8 :: dqsi, qssi
  REAL*8 :: dqsid, qssid
!Locals
  REAL*8, PARAMETER :: max_mixing_ratio=1.0
  REAL*8, PARAMETER :: cons_h2omw=18.01, cons_airmw=28.97
  REAL*8, PARAMETER :: esfac=cons_h2omw/cons_airmw
  REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
  REAL*8 :: tld, ttd, tid, dqsatd, qsatd, qqd, pld, ppd, ddd
  INTEGER :: it
  INTEGER, PARAMETER :: degsubs=100
  REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
  INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
  INTRINSIC NINT
  INTRINSIC INT
  tld = tempd
  tl = temp
  pld = plod
  pl = plo
  ppd = 100.0*pld
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
    ddd = -((ppd-(1.0-esfac)*qqd)/(pp-(1.0-esfac)*qq)**2)
    dd = 1.0/(pp-(1.0-esfac)*qq)
    qsatd = esfac*(qqd*dd+qq*ddd)
    qsat = esfac*qq*dd
    dqsatd = esfac*degsubs*dqq*((ppd*dd+pp*ddd)*dd+pp*dd*ddd)
    dqsat = esfac*degsubs*dqq*pp*(dd*dd)
  END IF
  dqsid = dqsatd
  dqsi = dqsat
  qssid = qsatd
  qssi = qsat
END SUBROUTINE DQSAT_SUB_SCA_D

!  Differentiation of marshpalm in forward (tangent) mode:
!   variations   of useful results: diam3 w ve
!   with respect to varying inputs: rain pr
SUBROUTINE MARSHPALM_D(rain, raind, pr, prd, diam3, diam3d, ntotal, w, &
&  wd, ve, ved)
  IMPLICIT NONE
!Inputs
! in kg m^-2 s^-1, mbar
  REAL*8, INTENT(IN) :: rain, pr
  REAL*8, INTENT(IN) :: raind, prd
!Outputs
  REAL*8, INTENT(OUT) :: diam3, ntotal, w, ve
  REAL*8, INTENT(OUT) :: diam3d, wd, ved
!Locals
  INTEGER :: iqd
!cm^-3
  REAL*8, PARAMETER :: n0=0.08
  REAL*8 :: rain_day, slopr, diam1
  REAL*8 :: rain_dayd
  REAL*8 :: rx(8), d3x(8)
  REAL*8 :: rxd(8), d3xd(8)
  REAL*8 :: result1
  REAL*8 :: result1d
  INTRINSIC MAX
  INTRINSIC SQRT
  rain_day = 0.0
  slopr = 0.0
  diam1 = 0.0
!Marshall-Palmer sizes at different rain-rates: avg(D^3)
!RX = (/ 0.   , 5.   , 20.  , 80.  , 320. , 1280., 5120., 20480. /)  ! rain per in mm/day
  rxd(1) = 0.0_8
  rx(1) = 0.
  rxd(2) = 0.0_8
  rx(2) = 5.
  rxd(3) = 0.0_8
  rx(3) = 20.
  rxd(4) = 0.0_8
  rx(4) = 80.
  rxd(5) = 0.0_8
  rx(5) = 320.
  rxd(6) = 0.0_8
  rx(6) = 1280.
  rxd(7) = 0.0_8
  rx(7) = 5120.
  rxd(8) = 0.0_8
  rx(8) = 20480.
!D3X= (/ 0.019, 0.032, 0.043, 0.057, 0.076, 0.102, 0.137, 0.183  /)
  d3xd(1) = 0.0_8
  d3x(1) = 0.019
  d3xd(2) = 0.0_8
  d3x(2) = 0.032
  d3xd(3) = 0.0_8
  d3x(3) = 0.043
  d3xd(4) = 0.0_8
  d3x(4) = 0.057
  d3xd(5) = 0.0_8
  d3x(5) = 0.076
  d3xd(6) = 0.0_8
  d3x(6) = 0.102
  d3xd(7) = 0.0_8
  d3x(7) = 0.137
  d3xd(8) = 0.0_8
  d3x(8) = 0.183
  rain_dayd = 3600.*24.*raind
  rain_day = rain*3600.*24.
  IF (rain_day .LE. 0.00) THEN
    diam1 = 0.00
    diam3 = 0.00
    ntotal = 0.00
    w = 0.00
    diam3d = 0.0_8
  ELSE
    diam3d = 0.0_8
  END IF
  DO iqd=1,7
    IF (rain_day .LE. rx(iqd+1) .AND. rain_day .GT. rx(iqd)) THEN
      slopr = (d3x(iqd+1)-d3x(iqd))/(rx(iqd+1)-rx(iqd))
      diam3d = slopr*rain_dayd
      diam3 = d3x(iqd) + (rain_day-rx(iqd))*slopr
    END IF
  END DO
  IF (rain_day .GE. rx(8)) THEN
    diam3 = d3x(8)
    diam3d = 0.0_8
  END IF
  ntotal = 0.019*diam3
  diam3d = 0.664*diam3d
  diam3 = 0.664*diam3
  result1d = -(1000.*prd/(pr**2*2.0*SQRT(1000./pr)))
  result1 = SQRT(1000./pr)
  wd = 2483.8*diam3d*result1 + (2483.8*diam3+80.)*result1d
  w = (2483.8*diam3+80.)*result1
  IF (0.99*w/100. .LT. 1.000) THEN
    ve = 1.000
    ved = 0.0_8
  ELSE
    ved = 0.99*wd/100.
    ve = 0.99*w/100.
  END IF
  diam1 = 3.0*diam3
  diam1 = diam1/100.
  diam3d = diam3d/100.
  diam3 = diam3/100.
  wd = wd/100.
  w = w/100.
  ntotal = ntotal*1.0e6
END SUBROUTINE MARSHPALM_D
