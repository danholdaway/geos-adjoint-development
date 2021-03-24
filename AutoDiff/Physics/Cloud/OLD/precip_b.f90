!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.7 (r4786) - 21 Feb 2013 15:53
!
!  Differentiation of precipandevap in reverse (adjoint) mode:
!   gradient     of useful results: evap_dd_above_out aa area pfl_above_in
!                mass imass qv pfi_above_in bb subl_dd_above_out
!                evap_dd_above_in qcl pfi_above_out subl_dd_above_in
!                pl rhcr3 dze qddf3 pfl_above_out te
!   with respect to varying inputs: evap_dd_above_out aa area pfl_above_in
!                mass imass qv pfi_above_in bb subl_dd_above_out
!                evap_dd_above_in qcl qpi pfi_above_out qpl subl_dd_above_in
!                pl rhcr3 dze qddf3 pfl_above_out te
!   RW status of diff variables: evap_dd_above_out:in-zero aa:incr
!                area:incr pfl_above_in:incr mass:incr imass:incr
!                qv:in-out pfi_above_in:incr bb:incr subl_dd_above_out:in-zero
!                evap_dd_above_in:incr qcl:in-out qpi:out pfi_above_out:in-zero
!                qpl:out subl_dd_above_in:incr pl:incr rhcr3:incr
!                dze:incr qddf3:incr pfl_above_out:in-zero te:in-out
SUBROUTINE PRECIPANDEVAP_B(k, lm, dt, frland, rhcr3, rhcr3b, qpl, qplb, &
&  qpi, qpib, qcl, qclb, qci, te, teb, qv, qvb, mass, massb, imass, &
&  imassb, pl, plb, dze, dzeb, qddf3, qddf3b, aa, aab, bb, bbb, area, &
&  areab, pfl_above_in, pfl_above_inb, pfl_above_out, pfl_above_outb, &
&  pfi_above_in, pfi_above_inb, pfi_above_out, pfi_above_outb, &
&  evap_dd_above_in, evap_dd_above_inb, evap_dd_above_out, &
&  evap_dd_above_outb, subl_dd_above_in, subl_dd_above_inb, &
&  subl_dd_above_out, subl_dd_above_outb, frz_diag, envfc, ddrfc, &
&  cons_alhf, cons_alhs, cons_alhl, cons_cp, cons_tice, revap_off_p, &
&  c_acc, c_ev_r, c_ev_s, rho_w, estblx)
  IMPLICIT NONE
!REVAP_DIAG = REVAP_DIAG + EVAP / DT
!RSUBL_DIAG = RSUBL_DIAG + SUBL / DT
!Inputs
  INTEGER, INTENT(IN) :: k, lm
  REAL*8, INTENT(IN) :: dt, mass, imass, pl, aa, bb, rhcr3, dze, qddf3, &
&  area, frland, envfc, ddrfc
  REAL*8 :: massb, imassb, plb, aab, bbb, rhcr3b, dzeb, qddf3b, areab
  REAL*8, INTENT(IN) :: cons_alhf, cons_alhs, cons_alhl, cons_cp, &
&  cons_tice, revap_off_p
  REAL*8, INTENT(IN) :: c_acc, c_ev_r, c_ev_s, rho_w
  REAL*8, INTENT(IN) :: estblx(50)
!Prognostics
  REAL*8, INTENT(INOUT) :: qv, qpl, qpi, qcl, qci, te, frz_diag
  REAL*8 :: qvb, qplb, qpib, qclb, teb
  REAL*8, INTENT(INOUT) :: pfl_above_in, pfl_above_out, pfi_above_in, &
&  pfi_above_out
  REAL*8 :: pfl_above_inb, pfl_above_outb, pfi_above_inb, pfi_above_outb
  REAL*8, INTENT(INOUT) :: evap_dd_above_in, evap_dd_above_out, &
&  subl_dd_above_in, subl_dd_above_out
  REAL*8 :: evap_dd_above_inb, evap_dd_above_outb, subl_dd_above_inb, &
&  subl_dd_above_outb
!Outputs - not used for now
!real(8), intent(  out) :: RAIN, SNOW, REVAP_DIAG, RSUBL_DIAG, ACRLL_DIAG, ACRIL_DIAG, PFL_DIAG, PFI_DIAG, VFALLSN, VFALLRN
!Locals
  INTEGER :: ns, nsmx, itr, l
  REAL*8 :: pfi, pfl, qs, dqs, envfrac, tko, qko, qstko, dqstko, rh_box&
&  , t_ed, qplko, qpiko
  REAL*8 :: pfib, pflb, qsb, dqsb, tkob, qkob, qstkob, dqstkob, rh_boxb&
&  , t_edb
  REAL*8 :: ifactor, rainrat0, snowrat0, fallrn, fallsn, vesn, vern, &
&  nrain, nsnow, efactor
  REAL*8 :: ifactorb, rainrat0b, snowrat0b, fallrnb, fallsnb, vesnb, &
&  vernb, efactorb
  REAL*8 :: tinlayerrn, diamrn, droprad, tinlayersn, diamsn, flakrad
  REAL*8 :: tinlayerrnb, diamrnb, dropradb, tinlayersnb, diamsnb, &
&  flakradb
  REAL*8 :: evap, subl, accr, mltfrz, evapx, sublx, evap_dd, subl_dd, &
&  ddfract, landseaf
  REAL*8 :: evapb, sublb, accrb, mltfrzb, evapxb, sublxb, evap_ddb, &
&  subl_ddb
  REAL*8 :: tau_frz, tau_mlt
!m/s
  REAL*8, PARAMETER :: trmv_l=1.0
  LOGICAL, PARAMETER :: taneff=.false.
!Fraction of precip falling through "environment" vs through cloud
  REAL*8, PARAMETER :: b_sub=1.00
  INTEGER :: branch
  REAL*8 :: temp3
  REAL*8 :: temp2
  REAL*8 :: temp1
  REAL*8 :: temp0
  INTRINSIC EXP
  REAL*8 :: tempb8
  REAL*8 :: tempb7
  REAL*8 :: tempb6
  REAL*8 :: tempb5
  REAL*8 :: tempb4
  REAL*8 :: tempb3
  REAL*8 :: tempb2
  REAL*8 :: tempb1
  REAL*8 :: tempb0
  INTRINSIC MAX
  REAL*8 :: temp5b0
  REAL*8 :: tempb
  REAL*8 :: temp2b
  REAL*8 :: temp5b
  INTRINSIC MIN
  REAL*8 :: temp
  REAL*8 :: temp4
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
    ifactor = 1./area
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
    ifactor = 1.00
  END IF
  IF (ifactor .LT. 1.) THEN
    ifactor = 1.
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    ifactor = ifactor
  END IF
!Start at top of precip column:
!
!   a) Accrete                   
!   b) Evaporate/Sublimate  
!   c) Rain/Snow-out to next level down 
!   d) return to (a)
!Update saturated humidity
  CALL DQSAT_SUB_SCA(dqs, qs, te, pl, estblx)
  ddfract = ddrfc
  IF (k .EQ. 1) THEN
    evap_dd = 0.
    subl_dd = 0.
!VFALLRN = 0.0
!VFALLSN = 0.0
    CALL PUSHCONTROL1B(0)
  ELSE
    qpl = qpl + pfl_above_in*imass
    qpi = qpi + pfi_above_in*imass
    accr = b_sub*c_acc*(qpl*mass)*qcl
    IF (accr .GT. qcl) THEN
      accr = qcl
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
      accr = accr
    END IF
    CALL PUSHREAL8(qpl)
    qpl = qpl + accr
    CALL PUSHREAL8(qcl)
    qcl = qcl - accr
!ACRLL_DIAG = ACCR / DT
!Accretion of liquid condensate by falling ice/snow
    accr = b_sub*c_acc*(qpi*mass)*qcl
    IF (accr .GT. qcl) THEN
      accr = qcl
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
      accr = accr
    END IF
    CALL PUSHREAL8(qpi)
    qpi = qpi + accr
    CALL PUSHREAL8(te)
!! Liquid freezes when accreted by snow
    te = te + cons_alhf*accr/cons_cp
!ACRIL_DIAG = ACCR / DT
    rainrat0 = ifactor*qpl*mass/dt
    snowrat0 = ifactor*qpi*mass/dt
    CALL PUSHREAL8(diamrn)
    CALL MARSHPALM(rainrat0, pl, diamrn, nrain, fallrn, vern)
    CALL PUSHREAL8(diamsn)
    CALL MARSHPALM(snowrat0, pl, diamsn, nsnow, fallsn, vesn)
!VFALLRN = FALLrn
!VFALLSN = FALLsn
    tinlayerrn = dze/(fallrn+0.01)
    tinlayersn = dze/(fallsn+0.01)
!Melting of Frozen precipitation      
! time scale for freezing (s). 
    tau_frz = 5000.
    IF (te .GT. cons_tice .AND. te .LE. cons_tice + 5.) THEN
      mltfrz = tinlayersn*qpi*(te-cons_tice)/tau_frz
      IF (qpi .GT. mltfrz) THEN
        CALL PUSHCONTROL1B(0)
        mltfrz = mltfrz
      ELSE
        mltfrz = qpi
        CALL PUSHCONTROL1B(1)
      END IF
      CALL PUSHREAL8(te)
      te = te - cons_alhf*mltfrz/cons_cp
      CALL PUSHREAL8(qpl)
      qpl = qpl + mltfrz
      CALL PUSHREAL8(qpi)
      qpi = qpi - mltfrz
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (te .GT. cons_tice + 5.) THEN
! Go Ahead and melt any snow/hail left above 5 C 
      mltfrz = qpi
      te = te - cons_alhf*mltfrz/cons_cp
      CALL PUSHREAL8(qpl)
      qpl = qpl + mltfrz
      CALL PUSHREAL8(qpi)
      qpi = qpi - mltfrz
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (k .GE. lm - 1) THEN
      IF (te .GT. cons_tice + 0.) THEN
! Go Ahead and melt any snow/hail left above 0 C in lowest layers 
        mltfrz = qpi
        te = te - cons_alhf*mltfrz/cons_cp
        CALL PUSHREAL8(qpl)
        qpl = qpl + mltfrz
        CALL PUSHREAL8(qpi)
        qpi = qpi - mltfrz
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
!Freezing of liquid precipitation      
    IF (te .LE. cons_tice) THEN
      te = te + cons_alhf*qpl/cons_cp
      CALL PUSHREAL8(qpi)
      qpi = qpl + qpi
      CALL PUSHREAL8(qpl)
      qpl = 0.
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
!In the exp below, evaporation time scale is determined "microphysically" from temp, 
!press, and drop size. In this context C_EV becomes a dimensionless fudge-fraction. 
!Also remember that these microphysics are still only for liquid.
    qko = qv
    tko = te
    DO itr=1,3
      dqstko = dqs
      CALL PUSHREAL8(qstko)
      qstko = qs + dqstko*(tko-te)
      IF (qstko .LT. 1.0e-7) THEN
        qstko = 1.0e-7
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        qstko = qstko
      END IF
      CALL PUSHREAL8(rh_box)
      rh_box = qko/qstko
      IF (rh_box .LT. rhcr3) THEN
        CALL PUSHREAL8(efactor)
        efactor = rho_w*(aa+bb)/(rhcr3-rh_box)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHREAL8(efactor)
        efactor = 9.99e9
        CALL PUSHCONTROL1B(1)
      END IF
      landseaf = 1.00
!Rain falling
      IF (rh_box .LT. rhcr3 .AND. diamrn .GT. 0.00 .AND. pl .GT. 100. &
&          .AND. pl .LT. revap_off_p) THEN
        droprad = 0.5*diamrn
        CALL PUSHREAL8(t_ed)
        t_ed = efactor*droprad**2
        CALL PUSHREAL8(t_ed)
        t_ed = t_ed*(1.0+dqstko*cons_alhl/cons_cp)
        evap = qpl*(1.0-EXP(-(c_ev_r*vern*landseaf*envfrac*tinlayerrn/&
&          t_ed)))
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        evap = 0.0
      END IF
!Snow falling
      IF (rh_box .LT. rhcr3 .AND. diamsn .GT. 0.00 .AND. pl .GT. 100. &
&          .AND. pl .LT. revap_off_p) THEN
        flakrad = 0.5*diamsn
        CALL PUSHREAL8(t_ed)
        t_ed = efactor*flakrad**2
        CALL PUSHREAL8(t_ed)
        t_ed = t_ed*(1.0+dqstko*cons_alhs/cons_cp)
        subl = qpi*(1.0-EXP(-(c_ev_s*vesn*landseaf*envfrac*tinlayersn/&
&          t_ed)))
        CALL PUSHCONTROL1B(0)
      ELSE
        subl = 0.0
        CALL PUSHCONTROL1B(1)
      END IF
      IF (itr .EQ. 1) THEN
        evapx = evap
        sublx = subl
        CALL PUSHCONTROL1B(0)
      ELSE
        evap = (evap+evapx)/2.0
        subl = (subl+sublx)/2.0
        CALL PUSHCONTROL1B(1)
      END IF
      CALL PUSHREAL8(qko)
      qko = qv + evap + subl
      CALL PUSHREAL8(tko)
      tko = te - evap*cons_alhl/cons_cp - subl*cons_alhs/cons_cp
    END DO
    CALL PUSHREAL8(qpi)
    qpi = qpi - subl
    CALL PUSHREAL8(qpl)
    qpl = qpl - evap
!Put some re-evap/re-subl precip in to a \quote{downdraft} to be applied later
    evap_dd = evap_dd_above_in + ddfract*evap*mass
    subl_dd = subl_dd_above_in + ddfract*subl*mass
!REVAP_DIAG = EVAP / DT
!RSUBL_DIAG = SUBL / DT
!PFL_DIAG =  PFl/DT
!PFI_DIAG =  PFi/DT
    CALL PUSHCONTROL1B(1)
  END IF
  evapb = qvb - cons_alhl*teb/cons_cp
  sublb = qvb - cons_alhs*teb/cons_cp
  temp5b = sublb/mass
  temp5b0 = evapb/mass
  subl_ddb = qddf3*temp5b + subl_dd_above_outb
  evap_ddb = qddf3*temp5b0 + evap_dd_above_outb
  pfib = pfi_above_outb
  pflb = pfl_above_outb
  qddf3b = qddf3b + evap_dd*temp5b0 + subl_dd*temp5b
  massb = massb - qddf3*evap_dd*temp5b0/mass - qddf3*subl_dd*temp5b/mass
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    qpib = mass*pfib
    massb = massb + qpl*pflb + qpi*pfib
    qplb = mass*pflb
    dqsb = 0.0_8
    qsb = 0.0_8
    ifactorb = 0.0_8
  ELSE
    qpib = mass*pfib
    massb = massb + qpl*pflb + ddfract*evap*evap_ddb + ddfract*subl*&
&      subl_ddb + qpi*pfib
    qplb = mass*pflb
    evapb = qvb - cons_alhl*teb/cons_cp
    sublb = qvb - cons_alhs*teb/cons_cp
    sublb = ddfract*mass*subl_ddb - qpib + (1.0_8-ddfract)*sublb
    subl_dd_above_inb = subl_dd_above_inb + subl_ddb
    evapb = ddfract*mass*evap_ddb - qplb + (1.0_8-ddfract)*evapb
    evap_dd_above_inb = evap_dd_above_inb + evap_ddb
    CALL POPREAL8(qpl)
    CALL POPREAL8(qpi)
    tkob = 0.0_8
    sublxb = 0.0_8
    dqsb = 0.0_8
    qsb = 0.0_8
    qkob = 0.0_8
    diamrnb = 0.0_8
    vernb = 0.0_8
    tinlayerrnb = 0.0_8
    diamsnb = 0.0_8
    evapxb = 0.0_8
    vesnb = 0.0_8
    tinlayersnb = 0.0_8
    DO itr=3,1,-1
      CALL POPREAL8(tko)
      teb = teb + tkob
      evapb = evapb + qkob - cons_alhl*tkob/cons_cp
      sublb = sublb + qkob - cons_alhs*tkob/cons_cp
      CALL POPREAL8(qko)
      qvb = qvb + qkob
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        sublb = sublb + sublxb
        evapb = evapb + evapxb
        sublxb = 0.0_8
        evapxb = 0.0_8
      ELSE
        sublxb = sublxb + sublb/2.0
        sublb = sublb/2.0
        evapxb = evapxb + evapb/2.0
        evapb = evapb/2.0
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        landseaf = 1.00
        temp2 = vesn*tinlayersn/t_ed
        temp4 = c_ev_s*landseaf*envfrac
        temp3 = -(temp4*temp2)
        temp2b = temp4*EXP(temp3)*qpi*sublb/t_ed
        qpib = qpib + (1.0-EXP(temp3))*sublb
        vesnb = vesnb + tinlayersn*temp2b
        tinlayersnb = tinlayersnb + vesn*temp2b
        t_edb = -(temp2*temp2b)
        dqstko = dqs
        flakrad = 0.5*diamsn
        CALL POPREAL8(t_ed)
        dqstkob = cons_alhs*t_ed*t_edb/cons_cp
        t_edb = (cons_alhs*(dqstko/cons_cp)+1.0)*t_edb
        CALL POPREAL8(t_ed)
        efactorb = flakrad**2*t_edb
        flakradb = efactor*2*flakrad*t_edb
        diamsnb = diamsnb + 0.5*flakradb
      ELSE
        landseaf = 1.00
        dqstko = dqs
        dqstkob = 0.0_8
        efactorb = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        temp = vern*tinlayerrn/t_ed
        temp1 = c_ev_r*landseaf*envfrac
        temp0 = -(temp1*temp)
        tempb8 = temp1*EXP(temp0)*qpl*evapb/t_ed
        qplb = qplb + (1.0-EXP(temp0))*evapb
        vernb = vernb + tinlayerrn*tempb8
        tinlayerrnb = tinlayerrnb + vern*tempb8
        t_edb = -(temp*tempb8)
        droprad = 0.5*diamrn
        CALL POPREAL8(t_ed)
        dqstkob = dqstkob + cons_alhl*t_ed*t_edb/cons_cp
        t_edb = (cons_alhl*(dqstko/cons_cp)+1.0)*t_edb
        CALL POPREAL8(t_ed)
        efactorb = efactorb + droprad**2*t_edb
        dropradb = efactor*2*droprad*t_edb
        diamrnb = diamrnb + 0.5*dropradb
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8(efactor)
        tempb6 = rho_w*efactorb/(rhcr3-rh_box)
        tempb7 = -((aa+bb)*tempb6/(rhcr3-rh_box))
        aab = aab + tempb6
        bbb = bbb + tempb6
        rhcr3b = rhcr3b + tempb7
        rh_boxb = -tempb7
      ELSE
        CALL POPREAL8(efactor)
        rh_boxb = 0.0_8
      END IF
      CALL POPREAL8(rh_box)
      qkob = rh_boxb/qstko
      qstkob = -(qko*rh_boxb/qstko**2)
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) qstkob = 0.0_8
      CALL POPREAL8(qstko)
      qsb = qsb + qstkob
      dqstkob = dqstkob + (tko-te)*qstkob
      tkob = dqstko*qstkob
      teb = teb - dqstko*qstkob
      dqsb = dqsb + dqstkob
      evapb = 0.0_8
      sublb = 0.0_8
    END DO
    teb = teb + tkob
    qvb = qvb + qkob
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(qpl)
      CALL POPREAL8(qpi)
      qplb = cons_alhf*teb/cons_cp + qpib
    END IF
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(qpi)
      mltfrzb = qplb - cons_alhf*teb/cons_cp - qpib
      CALL POPREAL8(qpl)
      qpib = qpib + mltfrzb
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(qpi)
      mltfrzb = qplb - cons_alhf*teb/cons_cp - qpib
      CALL POPREAL8(qpl)
      qpib = qpib + mltfrzb
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(qpi)
      mltfrzb = qplb - cons_alhf*teb/cons_cp - qpib
      CALL POPREAL8(qpl)
      CALL POPREAL8(te)
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        qpib = qpib + mltfrzb
        mltfrzb = 0.0_8
      END IF
      tempb5 = (te-cons_tice)*mltfrzb/tau_frz
      tinlayersnb = tinlayersnb + qpi*tempb5
      qpib = qpib + tinlayersn*tempb5
      teb = teb + tinlayersn*qpi*mltfrzb/tau_frz
    END IF
    tempb2 = tinlayerrnb/(fallrn+0.01)
    tempb1 = tinlayersnb/(fallsn+0.01)
    dzeb = dzeb + tempb2 + tempb1
    fallsnb = -(dze*tempb1/(fallsn+0.01))
    fallrnb = -(dze*tempb2/(fallrn+0.01))
    CALL POPREAL8(diamsn)
    CALL MARSHPALM_B(snowrat0, snowrat0b, pl, plb, diamsn, diamsnb, &
&               nsnow, fallsn, fallsnb, vesn, vesnb)
    CALL POPREAL8(diamrn)
    CALL MARSHPALM_B(rainrat0, rainrat0b, pl, plb, diamrn, diamrnb, &
&               nrain, fallrn, fallrnb, vern, vernb)
    tempb3 = mass*snowrat0b/dt
    qpib = qpib + ifactor*tempb3
    massb = massb + ifactor*qpl*rainrat0b/dt + ifactor*qpi*snowrat0b/dt
    tempb4 = mass*rainrat0b/dt
    ifactorb = qpl*tempb4 + qpi*tempb3
    qplb = qplb + ifactor*tempb4
    CALL POPREAL8(te)
    accrb = qpib - qclb + cons_alhf*teb/cons_cp
    CALL POPREAL8(qpi)
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      qclb = qclb + accrb
      accrb = 0.0_8
    END IF
    tempb0 = b_sub*c_acc*accrb
    qpib = qpib + qcl*mass*tempb0
    massb = massb + qcl*qpi*tempb0
    qclb = qclb + qpi*mass*tempb0
    CALL POPREAL8(qcl)
    accrb = qplb - qclb
    CALL POPREAL8(qpl)
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      qclb = qclb + accrb
      accrb = 0.0_8
    END IF
    tempb = b_sub*c_acc*accrb
    qplb = qplb + qcl*mass*tempb
    massb = massb + qcl*qpl*tempb
    qclb = qclb + qpl*mass*tempb
    pfi_above_inb = pfi_above_inb + imass*qpib
    imassb = imassb + pfl_above_in*qplb + pfi_above_in*qpib
    pfl_above_inb = pfl_above_inb + imass*qplb
  END IF
  CALL DQSAT_SUB_SCA_B(dqs, dqsb, qs, qsb, te, teb, pl, plb, estblx)
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) ifactorb = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) areab = areab - ifactorb/area**2
  evap_dd_above_outb = 0.0_8
  subl_dd_above_outb = 0.0_8
  pfi_above_outb = 0.0_8
  pfl_above_outb = 0.0_8
END SUBROUTINE PRECIPANDEVAP_B

!  Differentiation of dqsat_sub_sca in reverse (adjoint) mode:
!   gradient     of useful results: temp dqsi qssi plo
!   with respect to varying inputs: temp plo
SUBROUTINE DQSAT_SUB_SCA_B(dqsi, dqsib, qssi, qssib, temp, tempb, plo, &
&  plob, estblx)
  IMPLICIT NONE
!Inputs
  REAL*8 :: temp, plo
  REAL*8 :: tempb, plob
  REAL*8 :: estblx(:)
!Outputs
  REAL*8 :: dqsi, qssi
  REAL*8 :: dqsib, qssib
!Locals
  REAL*8, PARAMETER :: max_mixing_ratio=1.0
  REAL*8, PARAMETER :: cons_h2omw=18.01, cons_airmw=28.97
  REAL*8, PARAMETER :: esfac=cons_h2omw/cons_airmw
  REAL*8 :: tl, tt, ti, dqsat, qsat, dqq, qq, pl, pp, dd
  REAL*8 :: tlb, ttb, tib, dqsatb, qsatb, qqb, plb, ppb, ddb
  INTEGER :: it
  INTEGER, PARAMETER :: degsubs=100
  REAL*8, PARAMETER :: tmintbl=150.0, tmaxtbl=333.0
  INTEGER, PARAMETER :: tablesize=NINT(tmaxtbl-tmintbl)*degsubs+1
  INTEGER :: branch
  REAL*8 :: temp0
  REAL*8 :: temp0b
  INTRINSIC NINT
  INTRINSIC INT
  REAL*8 :: temp1b
  tl = temp
  pl = plo
  pp = pl*100.0
  IF (tl .LE. tmintbl) THEN
    ti = tmintbl
    CALL PUSHCONTROL2B(0)
  ELSE IF (tl .GE. tmaxtbl - .001) THEN
    ti = tmaxtbl - .001
    CALL PUSHCONTROL2B(1)
  ELSE
    ti = tl
    CALL PUSHCONTROL2B(2)
  END IF
  tt = (ti-tmintbl)*degsubs + 1
  it = INT(tt)
  dqq = estblx(it+1) - estblx(it)
  qq = (tt-it)*dqq + estblx(it)
  IF (pp .LE. qq) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    dd = 1.0/(pp-(1.0-esfac)*qq)
    CALL PUSHCONTROL1B(1)
  END IF
  qsatb = qssib
  dqsatb = dqsib
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    qqb = 0.0_8
    ppb = 0.0_8
  ELSE
    temp1b = esfac*degsubs*dqq*dqsatb
    ddb = esfac*qq*qsatb + pp*2*dd*temp1b
    temp0 = pp - (-esfac+1.0)*qq
    temp0b = -(ddb/temp0**2)
    ppb = temp0b + dd**2*temp1b
    qqb = esfac*dd*qsatb - (1.0-esfac)*temp0b
  END IF
  ttb = dqq*qqb
  tib = degsubs*ttb
  CALL POPCONTROL2B(branch)
  IF (branch .EQ. 0) THEN
    tlb = 0.0_8
  ELSE IF (branch .EQ. 1) THEN
    tlb = 0.0_8
  ELSE
    tlb = tib
  END IF
  plb = 100.0*ppb
  plob = plob + plb
  tempb = tempb + tlb
END SUBROUTINE DQSAT_SUB_SCA_B

!  Differentiation of marshpalm in reverse (adjoint) mode:
!   gradient     of useful results: diam3 w ve pr
!   with respect to varying inputs: rain pr
SUBROUTINE MARSHPALM_B(rain, rainb, pr, prb, diam3, diam3b, ntotal, w, &
&  wb, ve, veb)
  IMPLICIT NONE
!Inputs
! in kg m^-2 s^-1, mbar
  REAL*8, INTENT(IN) :: rain, pr
  REAL*8 :: rainb, prb
!Outputs
  REAL*8 :: diam3, ntotal, w, ve
  REAL*8 :: diam3b, wb, veb
!Locals
  INTEGER :: iqd
!cm^-3
  REAL*8, PARAMETER :: n0=0.08
  REAL*8 :: rain_day, slopr, diam1
  REAL*8 :: rain_dayb
  REAL*8 :: rx(8), d3x(8)
  INTEGER :: branch
  REAL*8 :: temp0
  INTRINSIC MAX
  INTRINSIC SQRT
  REAL*8 :: temp
!Marshall-Palmer sizes at different rain-rates: avg(D^3)
!RX = (/ 0.   , 5.   , 20.  , 80.  , 320. , 1280., 5120., 20480. /)  ! rain per in mm/day
  rx(1) = 0.
  rx(2) = 5.
  rx(3) = 20.
  rx(4) = 80.
  rx(5) = 320.
  rx(6) = 1280.
  rx(7) = 5120.
  rx(8) = 20480.
!D3X= (/ 0.019, 0.032, 0.043, 0.057, 0.076, 0.102, 0.137, 0.183  /)
  d3x(1) = 0.019
  d3x(2) = 0.032
  d3x(3) = 0.043
  d3x(4) = 0.057
  d3x(5) = 0.076
  d3x(6) = 0.102
  d3x(7) = 0.137
  d3x(8) = 0.183
  rain_day = rain*3600.*24.
  IF (rain_day .LE. 0.00) diam3 = 0.00
  DO iqd=1,7
    IF (rain_day .LE. rx(iqd+1) .AND. rain_day .GT. rx(iqd)) THEN
      CALL PUSHREAL8(slopr)
      slopr = (d3x(iqd+1)-d3x(iqd))/(rx(iqd+1)-rx(iqd))
      diam3 = d3x(iqd) + (rain_day-rx(iqd))*slopr
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  IF (rain_day .GE. rx(8)) THEN
    diam3 = d3x(8)
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
  END IF
  diam3 = 0.664*diam3
  w = (2483.8*diam3+80.)*SQRT(1000./pr)
  IF (0.99*w/100. .LT. 1.000) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  wb = wb/100.
  diam3b = diam3b/100.
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) wb = wb + 0.99*veb/100.
  temp0 = 1000./pr
  temp = SQRT(temp0)
  diam3b = diam3b + temp*2483.8*wb
  IF (.NOT.temp0 .EQ. 0.0) prb = prb - (2483.8*diam3+80.)*temp0*wb/(pr*&
&      2.0*temp)
  diam3b = 0.664*diam3b
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) diam3b = 0.0_8
  rain_dayb = 0.0_8
  DO iqd=7,1,-1
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      rain_dayb = rain_dayb + slopr*diam3b
      CALL POPREAL8(slopr)
      diam3b = 0.0_8
    END IF
  END DO
  rainb = 24.*3600.*rain_dayb
END SUBROUTINE MARSHPALM_B
