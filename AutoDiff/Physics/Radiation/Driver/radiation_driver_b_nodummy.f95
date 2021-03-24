SUBROUTINE RADIATION_DRIVER_B(im, jm, lm, runalarm, dt, ut, ptt, pttb, &
& qvt, qvtb, o3t, o3tb, cflst, cflstb, cfcnt, cfcntb, qit, qitb, qlt, &
& qltb, ps, ts, emis, delt, cosz, slr, rgbuv, rgfuv, rgbir, rgfir, sc, &
& du001t, du001tb, du002t, du002tb, du003t, du003tb, du004t, du004tb, &
& du005t, du005tb, taua_ir_c, ssaa_ir_c, asya_ir_c, taua_so_c, ssaa_so_c&
& , asya_so_c, ilsf, icnf, llsf, lcnf, mapl_p00, mapl_rgas, mapl_cp, &
& mapl_p00, mapl_kappa, mapl_grav, aib_ir, awb_ir, aiw_ir, aww_ir, &
& aig_ir, awg_ir, xkw, xke, mw, aw, bw, pm, fkw, gkw, cb, dcb, w11, w12&
& , w13, p11, p12, p13, dwe, dpe, c1, c2, c3, oo1, oo2, oo3, h11, h12, &
& h13, h21, h22, h23, h81, h82, h83, wk_uv, zk_uv, ry_uv, xk_ir, ry_ir, &
& cah, coa, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv, arb_uv, aib_nir, &
& awb_nir, arb_nir, aia_nir, awa_nir, ara_nir, aig_nir, awg_nir, arg_nir&
& , caib, caif, hk, hk_ir_temp, hk_uv_temp, hk_uv_old, hk_ir_old)
  IMPLICIT NONE
!In
  INTEGER, INTENT(IN) :: im, jm, lm, runalarm
  INTEGER, INTENT(IN) :: na, nbchou_ir, nbchou_so
  REAL, INTENT(IN) :: sc, dt
  REAL, DIMENSION(im, jm, lm), INTENT(IN) :: ut, qvt, o3t, cflst, cfcnt&
& , qit, qlt
  REAL, DIMENSION(im, jm, lm) :: qvtb, o3tb, cflstb, cfcntb, qitb, qltb
  REAL, DIMENSION(im, jm), INTENT(IN) :: ps, ts
  REAL, DIMENSION(im, jm), INTENT(IN) :: emis, delt, cosz, slr, rgbuv, &
& rgfuv, rgbir, rgfir
  REAL, DIMENSION(im, jm, lm), INTENT(IN) :: du001t, du002t, du003t, &
& du004t, du005t
  REAL, DIMENSION(im, jm, lm) :: du001tb, du002tb, du003tb, du004tb, &
& du005tb
  REAL, DIMENSION(im, jm, lm), INTENT(IN) :: ilsf, icnf, llsf, lcnf
  REAL, DIMENSION(im, jm, lm, nbchou_ir, na), INTENT(IN) :: taua_ir_c, &
& ssaa_ir_c, asya_ir_c
  REAL, DIMENSION(im, jm, lm, nbchou_so, na), INTENT(IN) :: taua_so_c, &
& ssaa_so_c, asya_so_c
  REAL, INTENT(IN) :: mapl_p00, mapl_rgas, mapl_cp, mapl_p00, mapl_kappa&
& , mapl_grav
  REAL, INTENT(IN) :: aib_ir(3, 10), awb_ir(4, 10), aiw_ir(4, 10)
  REAL, INTENT(IN) :: aww_ir(4, 10), aig_ir(4, 10), awg_ir(4, 10)
  INTEGER, INTENT(IN) :: mw(9)
  REAL, INTENT(IN) :: xkw(9), xke(9), aw(9), bw(9), pm(9)
  REAL, INTENT(IN) :: fkw(6, 9), gkw(6, 3), cb(6, 10), dcb(5, 10)
  REAL, INTENT(IN) :: w11, w12, w13, p11, p12
  REAL, INTENT(IN) :: p13, dwe, dpe
  REAL, INTENT(IN) :: c1(26, 30), c2(26, 30), c3(26, 30)
  REAL, INTENT(IN) :: oo1(26, 21), oo2(26, 21), oo3(26, 21)
  REAL, INTENT(IN) :: h11(26, 31), h12(26, 31), h13(26, 31)
  REAL, INTENT(IN) :: h21(26, 31), h22(26, 31), h23(26, 31)
  REAL, INTENT(IN) :: h81(26, 31), h82(26, 31), h83(26, 31)
  REAL, INTENT(IN) :: wk_uv(5), zk_uv(5), ry_uv(5)
  REAL, INTENT(IN) :: xk_ir(10), ry_ir(3)
  REAL, INTENT(IN) :: cah(43, 37), coa(62, 101)
  REAL, INTENT(IN) :: aig_uv(3), awg_uv(3), arg_uv(3)
  REAL, INTENT(IN) :: aib_uv, awb_uv(2), arb_uv(2)
  REAL, INTENT(IN) :: aib_nir, awb_nir(3, 2), arb_nir(3, 2)
  REAL, INTENT(IN) :: aia_nir(3, 3), awa_nir(3, 3), ara_nir(3, 3)
  REAL, INTENT(IN) :: aig_nir(3, 3), awg_nir(3, 3), arg_nir(3, 3)
  REAL, INTENT(IN) :: caib(11, 9, 11), caif(9, 11)
  REAL, INTENT(IN) :: hk(8), hk_ir_temp(3, 10), hk_uv_temp(5)
  REAL, INTENT(IN) :: hk_uv_old(5), hk_ir_old(3, 10)
!Inouts
  REAL, DIMENSION(im, jm, lm), INTENT(INOUT) :: ptt
  REAL, DIMENSION(im, jm, lm), INTENT(INOUT) :: pttb
!Locals
  INTEGER :: i, j, k, l
  REAL :: co2
  INTEGER :: levs925, lcldmh, lcldlm
  LOGICAL, PARAMETER :: trace=.true., overcast=.false.
  INTEGER, PARAMETER :: ns=1
  REAL, DIMENSION(im, jm, lm, 10) :: taua_irt, ssaa_irt, asya_irt
  REAL, DIMENSION(im, jm, lm, 10) :: taua_irtb, ssaa_irtb, asya_irtb
  REAL, DIMENSION(im, jm, lm, 8) :: taua_sot, ssaa_sot, asya_sot
  REAL, DIMENSION(im, jm, lm, 8) :: taua_sotb, ssaa_sotb, asya_sotb
  REAL, DIMENSION(im, jm, lm, na) :: aerot, saerot
  REAL, DIMENSION(im, jm, lm, na) :: aerotb, saerotb
  CHARACTER(len=50) :: aerosols(na)
  REAL :: x
  INTEGER :: ai
  REAL, DIMENSION(im, jm, lm) :: ptt1, tt, ple, pleso, plof, pif, deltap&
& , tt_out, dtdt
  REAL, DIMENSION(im, jm, lm) :: ptt1b, ttb, tt_outb, dtdtb
  REAL, DIMENSION(im, jm, lm) :: o3mmt
  REAL, DIMENSION(im, jm, lm) :: o3mmtb
  REAL, DIMENSION(im, jm, lm) :: n2ot, ch4t, cfc11t, cfc12t, cfc22t
  REAL, DIMENSION(im, jm, lm) :: rad_qlt, rad_qit, rad_cft, rad_rlt, &
& rad_rit
  REAL, DIMENSION(im, jm, lm) :: rad_qltb, rad_qitb, rad_cftb, rad_rltb&
& , rad_ritb
  REAL, DIMENSION(im, jm, lm, 4) :: cwct, refft
  REAL, DIMENSION(im, jm, lm, 4) :: cwctb, refftb
  REAL, DIMENSION(im, jm, 1) :: fst, tgt, tvt
  REAL, DIMENSION(im, jm, 1, 10) :: egt, evt, rvt
  REAL, DIMENSION(im, jm, lm) :: flwut, flwdt, flxnt, flwt_int, fswt_int
  REAL, DIMENSION(im, jm, lm) :: flwutb, flwdtb, flxntb, flwt_intb, &
& fswt_intb
  REAL, DIMENSION(im, jm, lm) :: flwavet, fswavet, dfdtst_int
  REAL, DIMENSION(im, jm, lm) :: flwavetb, fswavetb, dfdtst_intb
  REAL, DIMENSION(im, jm, lm) :: qilst, qllst, qicnt, qlcnt
  REAL, DIMENSION(im, jm, lm) :: qilstb, qllstb, qicntb, qlcntb
  REAL, DIMENSION(im, jm) :: t2mt, tempor, cosztmp
  REAL, DIMENSION(im, jm) :: t2mtb
  REAL, DIMENSION(0:lm) :: ak, bk, pref
  INTRINSIC MAX
  INTRINSIC SUM
  INTEGER :: arg1
  INTEGER :: branch
  REAL :: tempb(im, jm, lm)
  LOGICAL :: mask(im, jm)
!Compute pref from ak and bk
  ak = 1.0
  bk = 1.0
  DO l=0,lm
    pref(l) = ak(l+1) + bk(l+1)*mapl_p00
  END DO
!Pressure at the half levels from Ps
  DO l=0,lm
    ple(:, :, l) = ak(l+1) + bk(l+1)*ps(:, :)
  END DO
  deltap = ple(:, :, 1:lm) - ple(:, :, 0:lm-1)
!Pressure (hPa) and Exner pressure at the full levels
  plof(:, :, 1:lm) = 0.01*0.5*(ple(:, :, 0:lm-1)+ple(:, :, 1:lm))
  pif(:, :, 1:lm) = (plof(:, :, 1:lm)/1000.)**(mapl_rgas/mapl_cp)
!Potential temperature with p0=10000
!Place in holder so as not to overwrite
  ptt1 = ptt*mapl_p00**mapl_kappa
!Temperature
!Save input for later.
  tt = ptt1*pif
  IF (runalarm .EQ. 1) THEN
!2m temperature
    t2mt = tt(:, :, lm)*(0.5*(1.0+ple(:, :, lm-1)/ple(:, :, lm)))**(-&
&     mapl_kappa)
!Pressure in mb for SORAD
    pleso(:, :, 1:lm+1) = 0.01*ple(:, :, 0:lm)
!Caclulate TEMPOR for RAD-Cloud coupling
    levs925 = 0
    DO l=0,lm
      IF (pref(l) .LT. 92500.) levs925 = levs925 + 1
    END DO
    IF (1 .LT. levs925) THEN
      levs925 = levs925
    ELSE
      levs925 = 1
    END IF
    tempor = 0.
    DO l=levs925,lm
      mask(:, :) = ut(:, :, l) .GT. 4.
      WHERE (mask(:, :)) tempor(:, :) = 1.
    END DO
!Convert ozone from ppmv to kgkg 
    o3mmt = o3t/1.e6
    co2 = 1.0
!Initialize other gases, not available right now.
    n2ot = 0.0
    ch4t = 0.0
    cfc11t = 0.0
    cfc12t = 0.0
    cfc22t = 0.0
    qilst = qit*ilsf
    qicnt = qit*icnf
    qllst = qlt*llsf
    qlcnt = qlt*lcnf
!Initialize RAD/CLOUD variables
    rad_cft = 0.0
    rad_qlt = 0.0
    rad_qit = 0.0
    rad_rlt = 0.0
    rad_rit = 0.0
!Call RADcouple to produce cloud-rad variables
    DO i=1,im
      DO j=1,jm
        DO l=1,lm
          CALL RADCOUPLE(tt(i, j, l), plof(i, j, l), cflst(i, j, l), &
&                  cfcnt(i, j, l), qllst(i, j, l), qilst(i, j, l), qlcnt&
&                  (i, j, l), qicnt(i, j, l), rad_qlt(i, j, l), rad_qit(&
&                  i, j, l), rad_cft(i, j, l), rad_rlt(i, j, l), rad_rit&
&                  (i, j, l), tempor(i, j))
        END DO
      END DO
    END DO
!QI
    cwct(:, :, :, 1) = rad_qit
!QL
    cwct(:, :, :, 2) = rad_qlt
!QR - this is not available right now
    cwct(:, :, :, 3) = 0.0
!QI - this is not available right now
    cwct(:, :, :, 4) = 0.0
!RI
    refft(:, :, :, 1) = rad_rit*1.0e6
!RL
    refft(:, :, :, 2) = rad_rlt*1.0e6
!RR
    refft(:, :, :, 3) = 100.e-6*1.0e6
!RS
    refft(:, :, :, 4) = 140.e-6*1.0e6
!Determine the model level seperating high and middle clouds
    lcldmh = 1
    DO l=1,lm
      IF (pref(l) .GE. 40000.) GOTO 100
    END DO
    GOTO 110
 100 lcldmh = l
!Determine the model level seperating low and middle clouds
 110 lcldlm = lm
    DO l=1,lm
      IF (pref(l) .GE. 70000.) GOTO 120
    END DO
    GOTO 130
 120 lcldlm = l
!Set surface quantities
 130 egt = 0.0
    tgt = 0.0
    tvt = 0.0
    DO l=1,10
      egt(:, :, 1, l) = emis(:, :)
    END DO
    fst = 1.0
    tgt(:, :, 1) = ts
    tvt(:, :, 1) = ts
    evt = 0.0
    rvt = 0.0
!Get IRRAD aerosols
    IF (na .EQ. 0) THEN
      taua_irt = 0.0
      ssaa_irt = 0.0
      asya_irt = 0.0
      taua_sot = 0.0
      ssaa_sot = 0.0
      asya_sot = 0.0
      CALL PUSHCONTROL2B(0)
    ELSE IF (na .GT. 0) THEN
!Place aerosol fields into single rank 4 array that can be looped over
      aerot(:, :, :, 1) = du001t
      aerot(:, :, :, 2) = du002t
      aerot(:, :, :, 3) = du003t
      aerot(:, :, :, 4) = du004t
      aerot(:, :, :, 5) = du005t
!Names of aerosols as named in NL model
!Optical property calculation needs pressure as (1/g)dp/dz*q
      DO l=1,lm
        DO j=1,jm
          DO i=1,im
            CALL PUSHREAL4(x)
            x = (ple(i, j, l)-ple(i, j, l-1))*0.01*(100./mapl_grav)
            DO ai=1,na
              saerot(i, j, l, ai) = x*aerot(i, j, l, ai)
            END DO
          END DO
        END DO
      END DO
      taua_irt = 0.0
      ssaa_irt = 0.0
      asya_irt = 0.0
      DO j=1,nbchou_ir
        DO l=1,na
          taua_irt(:, :, :, j) = taua_irt(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_ir_c(:, :, :, j, l)
          ssaa_irt(:, :, :, j) = ssaa_irt(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_ir_c(:, :, :, j, l)*ssaa_ir_c(:, :, :, j, l)
          asya_irt(:, :, :, j) = asya_irt(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_ir_c(:, :, :, j, l)*ssaa_ir_c(:, :, :, j, l)*&
&           asya_ir_c(:, :, :, j, l)
        END DO
      END DO
      taua_sot = 0.0
      ssaa_sot = 0.0
      asya_sot = 0.0
      DO j=1,nbchou_so
        DO l=1,na
          taua_sot(:, :, :, j) = taua_sot(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_so_c(:, :, :, j, l)
          ssaa_sot(:, :, :, j) = ssaa_sot(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_so_c(:, :, :, j, l)*ssaa_so_c(:, :, :, j, l)
          asya_sot(:, :, :, j) = asya_sot(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_so_c(:, :, :, j, l)*ssaa_so_c(:, :, :, j, l)*&
&           asya_so_c(:, :, :, j, l)
        END DO
      END DO
      CALL PUSHCONTROL2B(1)
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
    IF (sc .LT. 0.0) THEN
      hk_uv_temp = hk(1:5)
      DO l=1,3
        hk_ir_temp(l, :) = hk_ir_old(l, :)*(hk(5+l)/SUM(hk_ir_old(l, :))&
&         )
      END DO
    ELSE
      hk_uv_temp = hk_uv_old
      hk_ir_temp = hk_ir_old
    END IF
    DO i=1,im
      DO j=1,jm
        IF (.0001 .LT. cosz(i, j)) THEN
          cosztmp(i, j) = cosz(i, j)
        ELSE
          cosztmp(i, j) = .0001
        END IF
      END DO
    END DO
!Initialize the fluxes
    flxnt = 0.0
!SOLAR RADIATION (SHORT WAVE)
    arg1 = im*jm
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
  END IF
  ptt1b = 0.0
  ptt1b = pttb/mapl_p00**mapl_kappa
  tt_outb = 0.0
  tt_outb = ptt1b/pif
  ttb = 0.0
  dtdtb = 0.0
  ttb = tt_outb
  dtdtb = dt*tt_outb/deltap
  flwavetb = 0.0
  fswavetb = 0.0
  tempb = mapl_grav*dtdtb/mapl_cp
  flwavetb(:, :, 0:lm-1) = flwavetb(:, :, 0:lm-1) + tempb
  flwavetb(:, :, 1:lm) = -tempb
  fswavetb(:, :, 0:lm-1) = fswavetb(:, :, 0:lm-1) + tempb
  fswavetb(:, :, 1:lm) = -tempb
  dfdtst_intb = 0.0
  fswt_intb = 0.0
  flwt_intb = 0.0
  DO l=lm,0,-1
    fswt_intb(:, :, l) = fswt_intb(:, :, l) + slr*fswavetb(:, :, l)
    fswavetb(:, :, l) = 0.0
    flwt_intb(:, :, l) = flwt_intb(:, :, l) + flwavetb(:, :, l)
    dfdtst_intb(:, :, l) = dfdtst_intb(:, :, l) + delt*flwavetb(:, :, l)
    flwavetb(:, :, l) = 0.0
  END DO
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    du005tb = 0.0
    du002tb = 0.0
    cfcntb = 0.0
    cflstb = 0.0
    du003tb = 0.0
    qltb = 0.0
    qitb = 0.0
    qvtb = 0.0
    o3tb = 0.0
    du004tb = 0.0
    du001tb = 0.0
  ELSE
    flwdtb = 0.0
    flwutb = 0.0
    flwutb = flwt_intb
    flwdtb = flwt_intb
    qvtb = 0.0
    t2mtb = 0.0
    ssaa_irtb = 0.0
    o3mmtb = 0.0
    rad_cftb = 0.0
    cwctb = 0.0
    asya_irtb = 0.0
    refftb = 0.0
    taua_irtb = 0.0
    DO i=im,1,-1
      DO j=jm,1,-1
        CALL IRRAD_B(1, lm, ple(i, j, :), tt(i, j, :), ttb(i, j, :), qvt&
&              (i, j, :), qvtb(i, j, :), o3mmt(i, j, :), o3mmtb(i, j, :)&
&              , t2mt(i, j), t2mtb(i, j), co2, trace, n2ot(i, j, :), &
&              ch4t(i, j, :), cfc11t(i, j, :), cfc12t(i, j, :), cfc22t(i&
&              , j, :), cwct(i, j, :, :), cwctb(i, j, :, :), rad_cft(i, &
&              j, :), rad_cftb(i, j, :), lcldmh, lcldlm, refft(i, j, :, &
&              :), refftb(i, j, :, :), ns, fst(i, j, :), tgt(i, j, :), &
&              egt(i, j, :, :), tvt(i, j, :), evt(i, j, :, :), rvt(i, j&
&              , :, :), na, nbchou_ir, taua_irt(i, j, :, :), taua_irtb(i&
&              , j, :, :), ssaa_irt(i, j, :, :), ssaa_irtb(i, j, :, :), &
&              asya_irt(i, j, :, :), asya_irtb(i, j, :, :), flwut(i, j, &
&              :), flwutb(i, j, :), flwdt(i, j, :), flwdtb(i, j, :), &
&              dfdtst_int(i, j, :), dfdtst_intb(i, j, :), overcast, &
&              aib_ir, awb_ir, aiw_ir, aww_ir, aig_ir, awg_ir, xkw, xke&
&              , mw, aw, bw, pm, fkw, gkw, cb, dcb, w11, w12, w13, p11, &
&              p12, p13, dwe, dpe, c1, c2, c3, oo1, oo2, oo3, h11, h12, &
&              h13, h21, h22, h23, h81, h82, h83)
      END DO
    END DO
    flxntb = 0.0
    flxntb = fswt_intb
    CALL SORAD_B(arg1, lm, nbchou_so, cosztmp, pleso, tt, ttb, qvt, qvtb&
&          , o3mmt, o3mmtb, co2, cwct, cwctb, rad_cft, rad_cftb, lcldmh&
&          , lcldlm, refft, refftb, hk_uv_temp, hk_ir_temp, taua_sot, &
&          taua_sotb, ssaa_sot, ssaa_sotb, asya_sot, asya_sotb, rgbuv, &
&          rgfuv, rgbir, rgfir, flxnt, flxntb, mapl_grav, wk_uv, zk_uv, &
&          ry_uv, xk_ir, ry_ir, cah, coa, aig_uv, awg_uv, arg_uv, aib_uv&
&          , awb_uv, arb_uv, aib_nir, awb_nir, arb_nir, aia_nir, awa_nir&
&          , ara_nir, aig_nir, awg_nir, arg_nir, caib, caif)
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      du005tb = 0.0
      du002tb = 0.0
      du003tb = 0.0
      du004tb = 0.0
      du001tb = 0.0
    ELSE IF (branch .EQ. 1) THEN
      saerotb = 0.0
      DO j=nbchou_so,1,-1
        DO l=na,1,-1
          saerotb(:, :, :, l) = saerotb(:, :, :, l) + taua_so_c(:, :, :&
&           , j, l)*ssaa_so_c(:, :, :, j, l)*ssaa_sotb(:, :, :, j) + &
&           taua_so_c(:, :, :, j, l)*taua_sotb(:, :, :, j) + taua_so_c(:&
&           , :, :, j, l)*ssaa_so_c(:, :, :, j, l)*asya_so_c(:, :, :, j&
&           , l)*asya_sotb(:, :, :, j)
        END DO
      END DO
      DO j=nbchou_ir,1,-1
        DO l=na,1,-1
          saerotb(:, :, :, l) = saerotb(:, :, :, l) + taua_ir_c(:, :, :&
&           , j, l)*ssaa_ir_c(:, :, :, j, l)*ssaa_irtb(:, :, :, j) + &
&           taua_ir_c(:, :, :, j, l)*taua_irtb(:, :, :, j) + taua_ir_c(:&
&           , :, :, j, l)*ssaa_ir_c(:, :, :, j, l)*asya_ir_c(:, :, :, j&
&           , l)*asya_irtb(:, :, :, j)
        END DO
      END DO
      aerotb = 0.0
      DO l=lm,1,-1
        DO j=jm,1,-1
          DO i=im,1,-1
            DO ai=na,1,-1
              aerotb(i, j, l, ai) = aerotb(i, j, l, ai) + x*saerotb(i, j&
&               , l, ai)
              saerotb(i, j, l, ai) = 0.0
            END DO
            CALL POPREAL4(x)
          END DO
        END DO
      END DO
      du005tb = 0.0
      du005tb = aerotb(:, :, :, 5)
      aerotb(:, :, :, 5) = 0.0
      du004tb = 0.0
      du004tb = aerotb(:, :, :, 4)
      aerotb(:, :, :, 4) = 0.0
      du003tb = 0.0
      du003tb = aerotb(:, :, :, 3)
      aerotb(:, :, :, 3) = 0.0
      du002tb = 0.0
      du002tb = aerotb(:, :, :, 2)
      aerotb(:, :, :, 2) = 0.0
      du001tb = 0.0
      du001tb = aerotb(:, :, :, 1)
    ELSE
      du005tb = 0.0
      du002tb = 0.0
      du003tb = 0.0
      du004tb = 0.0
      du001tb = 0.0
    END IF
    refftb(:, :, :, 4) = 0.0
    refftb(:, :, :, 3) = 0.0
    rad_rltb = 0.0
    rad_rltb = 1.0e6*refftb(:, :, :, 2)
    refftb(:, :, :, 2) = 0.0
    rad_ritb = 0.0
    rad_ritb = 1.0e6*refftb(:, :, :, 1)
    cwctb(:, :, :, 4) = 0.0
    cwctb(:, :, :, 3) = 0.0
    rad_qltb = 0.0
    rad_qltb = cwctb(:, :, :, 2)
    cwctb(:, :, :, 2) = 0.0
    rad_qitb = 0.0
    rad_qitb = cwctb(:, :, :, 1)
    cfcntb = 0.0
    cflstb = 0.0
    qlcntb = 0.0
    qllstb = 0.0
    qicntb = 0.0
    qilstb = 0.0
    DO i=im,1,-1
      DO j=jm,1,-1
        DO l=lm,1,-1
          CALL RADCOUPLE_B(tt(i, j, l), ttb(i, j, l), plof(i, j, l), &
&                    cflst(i, j, l), cflstb(i, j, l), cfcnt(i, j, l), &
&                    cfcntb(i, j, l), qllst(i, j, l), qllstb(i, j, l), &
&                    qilst(i, j, l), qilstb(i, j, l), qlcnt(i, j, l), &
&                    qlcntb(i, j, l), qicnt(i, j, l), qicntb(i, j, l), &
&                    rad_qlt(i, j, l), rad_qltb(i, j, l), rad_qit(i, j, &
&                    l), rad_qitb(i, j, l), rad_cft(i, j, l), rad_cftb(i&
&                    , j, l), rad_rlt(i, j, l), rad_rltb(i, j, l), &
&                    rad_rit(i, j, l), rad_ritb(i, j, l), tempor(i, j))
          rad_qltb(i, j, l) = 0.0
          rad_qitb(i, j, l) = 0.0
          rad_cftb(i, j, l) = 0.0
          rad_rltb(i, j, l) = 0.0
          rad_ritb(i, j, l) = 0.0
        END DO
      END DO
    END DO
    qltb = 0.0
    qltb = llsf*qllstb + lcnf*qlcntb
    qitb = 0.0
    qitb = ilsf*qilstb + icnf*qicntb
    o3tb = 0.0
    o3tb = o3mmtb/1.e6
    DO l=lm,levs925,-1
      mask(:, :) = ut(:, :, l) .GT. 4.
    END DO
    ttb(:, :, lm) = ttb(:, :, lm) + ((ple(:, :, lm-1)/ple(:, :, lm)+1.0)&
&     *0.5)**(-mapl_kappa)*t2mtb
  END IF
  ptt1b = 0.0
  ptt1b = pif*ttb
  pttb = 0.0
  pttb = mapl_p00**mapl_kappa*ptt1b
END SUBROUTINE RADIATION_DRIVER_B

!  Differentiation of radcouple in reverse (adjoint) mode:
!   gradient     of useful results: qcian af rad_cf qclls qcils
!                rad_qi rad_ql cf qclan rad_ri rad_rl te
!   with respect to varying inputs: qcian af qclls qcils cf qclan
!                te
SUBROUTINE RADCOUPLE_B(te, teb, pl, cf, cfb, af, afb, qclls, qcllsb, &
& qcils, qcilsb, qclan, qclanb, qcian, qcianb, rad_ql, rad_qlb, rad_qi, &
& rad_qib, rad_cf, rad_cfb, rad_rl, rad_rlb, rad_ri, rad_rib, tempor)
  IMPLICIT NONE
!Inputs
  REAL, INTENT(IN) :: te, pl, tempor
  REAL :: teb
  REAL, INTENT(IN) :: af, cf, qclan, qcian, qclls, qcils
  REAL :: afb, cfb, qclanb, qcianb, qcllsb, qcilsb
! real, intent(in ) :: QRN_ALL, QSN_ALL
!Outputs
  REAL :: rad_ql, rad_qi, rad_cf, rad_rl, rad_ri
  REAL :: rad_qlb, rad_qib, rad_cfb, rad_rlb, rad_rib
! real, intent(out) :: RAD_QR,RAD_QS
!Locals
  REAL :: ss, rad_ri_an, afx, alph
  REAL :: ssb, rad_ri_anb, afxb
  REAL, PARAMETER :: min_ri=20.e-6, max_ri=40.e-6, ri_anv=30.e-6
  INTRINSIC MIN
  INTRINSIC MAX
  INTEGER :: branch
  REAL :: min3
  REAL :: min2
  REAL :: max2b
  REAL :: min1
  REAL :: tempb2
  REAL :: tempb1
  REAL :: tempb0
  REAL :: min1b
  REAL :: x2
  REAL :: x1
  REAL :: x2b
  REAL :: tempb
  REAL :: max3b
  REAL :: temp
  REAL :: max3
  REAL :: max2
  REAL :: max1
  REAL :: min2b
!Initialize outputs
  rad_ri = 0.0
!RAD_QR = 0.0
!RAD_QS = 0.0
! Adjust Anvil fractions for warm clouds
  alph = 0.1
  ss = (280.-te)/20.
  IF (1.0 .GT. ss) THEN
    CALL PUSHCONTROL1B(0)
    ss = ss
  ELSE
    ss = 1.0
    CALL PUSHCONTROL1B(1)
  END IF
  IF (0.0 .LT. ss) THEN
    CALL PUSHCONTROL1B(0)
    ss = ss
  ELSE
    ss = 0.0
    CALL PUSHCONTROL1B(1)
  END IF
  CALL PUSHREAL4(ss)
  ss = alph + ss**3*(1.0-alph)
  afx = af*ss*0.5
  IF (cf + afx .GT. 1.00) THEN
    rad_cf = 1.00
    CALL PUSHCONTROL1B(0)
  ELSE
    rad_cf = cf + afx
    CALL PUSHCONTROL1B(1)
  END IF
!Total In-cloud liquid
  IF (rad_cf .GT. 10.0e-8) THEN
!0 -> 10e-8 FOR LINEARIZATION PROTECTION
    rad_ql = (qclls+qclan)/rad_cf
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
    rad_ql = 0.0
  END IF
  IF (rad_ql .GT. 0.01) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
! Total In-cloud ice
  IF (rad_cf .GT. 10.0e-8) THEN
!0 -> 10e-8 FOR LINEARIZATION PROTECTION
    rad_qi = (qcils+qcian)/rad_cf
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
    rad_qi = 0.0
  END IF
  IF (rad_qi .GT. 0.01) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
! Total In-cloud precipitation
! if (  RAD_CF >0. ) then
!    RAD_QR = ( QRN_ALL ) / RAD_CF
!    RAD_QS = ( QSN_ALL ) / RAD_CF
! else
!    RAD_QR = 0.0
!    RAD_QS = 0.0
! end if
! RAD_QR = MIN( RAD_QR, 0.01 )
! RAD_QS = MIN( RAD_QS, 0.01 )
  IF (pl .LT. 150.) rad_ri = max_ri
  IF (pl .GE. 150.) rad_ri = max_ri*150./pl
! Weigh in a separate R_ice for Anvil Ice according to
  rad_ri_an = rad_ri
  IF (qcils + qcian .GT. 0.0) THEN
    IF (qcils/rad_ri + qcian/ri_anv .GT. 10e-8) THEN
!LINEARIZATION PROTECTION
      rad_ri_an = (qcils+qcian)/(qcils/rad_ri+qcian/ri_anv)
      CALL PUSHCONTROL2B(2)
    ELSE
      CALL PUSHCONTROL2B(1)
    END IF
  ELSE
    CALL PUSHCONTROL2B(0)
  END IF
  IF (rad_ri .GT. rad_ri_an) THEN
    CALL PUSHREAL4(rad_ri)
    rad_ri = rad_ri_an
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHREAL4(rad_ri)
    rad_ri = rad_ri
    CALL PUSHCONTROL1B(1)
  END IF
  IF (rad_ri .LT. min_ri) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 825. .AND. te .LE. 282. .AND. tempor .EQ. 1.) THEN
    IF (0.71*te - 190.25 .LT. 5.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 775. .AND. pl .LT. 825. .AND. te .LE. 282. .AND. te .GT. &
&     275. .AND. tempor .EQ. 1.) THEN
    IF (-(0.1*pl) + 0.71*te - 107.75 .GT. 10.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 825. .AND. te .LE. 275. .AND. tempor .EQ. 1.) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
! Thin low tropical clouds
  IF (pl .GE. 950. .AND. te .GE. 285.) THEN
    IF (2.2*te - 617. .GT. 21.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 925. .AND. te .GE. 290.) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 925. .AND. pl .LT. 950. .AND. te .GT. 285. .AND. te .LT. &
&     290.) THEN
    IF (0.44*pl + 2.2*te - 1035. .GT. 21.) THEN
      CALL PUSHCONTROL1B(0)
      x2 = 21.
    ELSE
      x2 = 0.44*pl + 2.2*te - 1035.
      CALL PUSHCONTROL1B(1)
    END IF
    IF (x2 .LT. 10.) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (pl .GE. 950. .AND. te .GE. 290.) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  IF (rad_cf .LT. 1.e-5) THEN
    rad_cfb = 0.0
    rad_qib = 0.0
    rad_qlb = 0.0
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rlb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    max3b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      x2b = 0.0
    ELSE
      x2b = max3b
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) teb = teb + 2.2*x2b
    rad_rlb = 0.0
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rlb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    min2b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) teb = teb + 2.2*min2b
    rad_rlb = 0.0
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rlb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    min1b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) teb = teb + 0.71*min1b
    rad_rlb = 0.0
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    max2b = 1.e-6*rad_rlb
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) teb = teb + 0.71*max2b
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_rib = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    CALL POPREAL4(rad_ri)
    rad_ri_anb = rad_rib
  ELSE
    CALL POPREAL4(rad_ri)
    rad_ri_anb = 0.0
  END IF
  CALL POPCONTROL2B(branch)
  IF (branch .NE. 0) THEN
    IF (branch .NE. 1) THEN
      temp = qcils/rad_ri + qcian/ri_anv
      tempb1 = rad_ri_anb/temp
      tempb2 = -((qcils+qcian)*tempb1/temp)
      qcilsb = qcilsb + tempb2/rad_ri + tempb1
      qcianb = qcianb + tempb2/ri_anv + tempb1
    END IF
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_qib = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) THEN
    tempb0 = rad_qib/rad_cf
    qcilsb = qcilsb + tempb0
    qcianb = qcianb + tempb0
    rad_cfb = rad_cfb - (qcils+qcian)*tempb0/rad_cf
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rad_qlb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) THEN
    tempb = rad_qlb/rad_cf
    qcllsb = qcllsb + tempb
    qclanb = qclanb + tempb
    rad_cfb = rad_cfb - (qclls+qclan)*tempb/rad_cf
  END IF
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    afxb = 0.0
  ELSE
    cfb = cfb + rad_cfb
    afxb = rad_cfb
  END IF
  afb = afb + 0.5*ss*afxb
  ssb = 0.5*af*afxb
  CALL POPREAL4(ss)
  ssb = (1.0-alph)*3*ss**2*ssb
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) ssb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) ssb = 0.0
  teb = teb - ssb/20.
END SUBROUTINE RADCOUPLE_B

!  Differentiation of irrad in reverse (adjoint) mode:
!   gradient     of useful results: ssaa_dev wa_dev fcld_dev cwc_dev
!                flxu_dev ta_dev asya_dev tb_dev dfdts_dev reff_dev
!                taua_dev flxd_dev oa_dev
!   with respect to varying inputs: ssaa_dev wa_dev fcld_dev cwc_dev
!                flxu_dev ta_dev asya_dev tb_dev dfdts_dev reff_dev
!                taua_dev flxd_dev oa_dev
!Number of soundings (IM*JM)
!Number of layers (LM)
!Pressure at level edges (Pa)
!d - Temperature (K)
!d - Specific humidity (g/g)
!d - Ozone (g/g)
!Surface air temperature (K)
!* Carbon dioxide (pppv)
!Option
!* Nitrous oxide (pppv)
!* Methane (pppv)
!* Trichlorofluoromethane (pppv)
!* Dichlorodifluoromethane (pppv)
!* Chlorodifluoromethane (pppv)
!Cloud water mixing ratio (kg/kg) 
!Cloud amount (fraction)
!Level index separating high and middle clouds
!Level index separating middle and low clouds
!Effective size of cloud particles (micron)
!Number of sub-grid surface types
!Fractional cover of sub-grid regions
!Land or ocean surface temperature
!Land or ocean surface emissivity
!Vegetation temperature
!Vegetation emissivity
!Vegetation reflectivity 
!Number of bands
!Number of bands in IRRAD calcs for Chou
!Aerosol optical thickness
!Aerosol single scattering albedo
!Aerosol asymmetry factor
!Outputs
!Upwelling flux, all-sky
!                    flcu_dev    , & !Upwelling flux, clear-sky
!                    flau_dev    , & !Upwelling flux, clear-sky no aerosol
!Downwelling flux, all-sky
!                    flcd_dev    , & !Downwelling flux, clear-sky
!                    flad_dev    , & !Downwelling flux, clear-sky no aerosol
!Sensitivity of net downward flux to surface temperature
!                    sfcem_dev   , & !Emission by the surface
!Switch for overcast simplification
SUBROUTINE IRRAD_B(m, np, ple_dev, ta_dev, ta_devb, wa_dev, wa_devb, &
& oa_dev, oa_devb, tb_dev, tb_devb, co2, trace, n2o_dev, ch4_dev, &
& cfc11_dev, cfc12_dev, cfc22_dev, cwc_dev, cwc_devb, fcld_dev, &
& fcld_devb, ict, icb, reff_dev, reff_devb, ns, fs_dev, tg_dev, eg_dev, &
& tv_dev, ev_dev, rv_dev, na, nb, taua_dev, taua_devb, ssaa_dev, &
& ssaa_devb, asya_dev, asya_devb, flxu_dev, flxu_devb, flxd_dev, &
& flxd_devb, dfdts_dev, dfdts_devb, overcastl, aib_ir, awb_ir, aiw_ir, &
& aww_ir, aig_ir, awg_ir, xkw, xke, mw, aw, bw, pm, fkw, gkw, cb, dcb, &
& w11, w12, w13, p11, p12, p13, dwe, dpe, c1, c2, c3, oo1, oo2, oo3, h11&
& , h12, h13, h21, h22, h23, h81, h82, h83)
  IMPLICIT NONE
!end do
!Radiation constants, these need to be inputs for the autodiff tool
  REAL, INTENT(IN) :: aib_ir(3, 10), awb_ir(4, 10), aiw_ir(4, 10)
  REAL, INTENT(IN) :: aww_ir(4, 10), aig_ir(4, 10), awg_ir(4, 10)
  INTEGER, INTENT(IN) :: mw(9)
  REAL, INTENT(IN) :: xkw(9), xke(9), aw(9), bw(9), pm(9)
  REAL, INTENT(IN) :: fkw(6, 9), gkw(6, 3), cb(6, 10), dcb(5, 10)
  REAL, INTENT(IN) :: w11, w12, w13, p11, p12
  REAL, INTENT(IN) :: p13, dwe, dpe
  REAL, INTENT(IN) :: c1(26, 30), c2(26, 30), c3(26, 30)
  REAL, INTENT(IN) :: oo1(26, 21), oo2(26, 21), oo3(26, 21)
  REAL, INTENT(IN) :: h11(26, 31), h12(26, 31), h13(26, 31)
  REAL, INTENT(IN) :: h21(26, 31), h22(26, 31), h23(26, 31)
  REAL, INTENT(IN) :: h81(26, 31), h82(26, 31), h83(26, 31)
!----- INPUTS -----
  INTEGER, INTENT(IN) :: m, np, ict, icb, ns, na, nb
  LOGICAL, INTENT(IN) :: trace, overcastl
  REAL, INTENT(IN) :: co2
!Rank 2 inputs
  REAL, DIMENSION(m), INTENT(IN) :: tb_dev
  REAL, DIMENSION(m) :: tb_devb
!Rank 3 (Prognostic variables and tracers)
  REAL, DIMENSION(m, np), INTENT(IN) :: ta_dev, wa_dev, oa_dev, fcld_dev
  REAL, DIMENSION(m, np) :: ta_devb, wa_devb, oa_devb, fcld_devb
  REAL, DIMENSION(m, np), INTENT(IN) :: n2o_dev, ch4_dev, cfc11_dev, &
& cfc12_dev, cfc22_dev
  REAL, DIMENSION(m, np + 1), INTENT(IN) :: ple_dev
!Rank 3 (surface types)
  REAL, DIMENSION(m, ns), INTENT(IN) :: fs_dev, tg_dev, tv_dev
  REAL, DIMENSION(m, ns, 10), INTENT(IN) :: eg_dev, ev_dev, rv_dev
!Rank 3 (diagnostic cloud parts)
  REAL, DIMENSION(m, np, 4), INTENT(IN) :: cwc_dev, reff_dev
  REAL, DIMENSION(m, np, 4) :: cwc_devb, reff_devb
!Rank 3 (aerosols)
  REAL, DIMENSION(m, np, nb), INTENT(IN) :: taua_dev, ssaa_dev, asya_dev
  REAL, DIMENSION(m, np, nb) :: taua_devb, ssaa_devb, asya_devb
!----- OUPUTS -----
!REAL, DIMENSION(m), INTENT(OUT)         :: sfcem_dev
!, flcu_dev, flau_dev
  REAL, DIMENSION(m, np + 1) :: flxu_dev
  REAL, DIMENSION(m, np+1) :: flxu_devb
!, flcd_dev, flad_dev
  REAL, DIMENSION(m, np + 1) :: flxd_dev
  REAL, DIMENSION(m, np+1) :: flxd_devb
  REAL, DIMENSION(m, np + 1) :: dfdts_dev
  REAL, DIMENSION(m, np+1) :: dfdts_devb
!----- LOCALS -----
  REAL, PARAMETER :: cons_grav=9.80665
  INTEGER, PARAMETER :: nx1=26
  INTEGER, PARAMETER :: no1=21
  INTEGER, PARAMETER :: nc1=30
  INTEGER, PARAMETER :: nh1=31
!Temporary arrays
  REAL :: pa(0:np), dt(0:np)
  REAL :: dtb(0:np)
  REAL :: x1, x2, x3
  REAL :: x1b, x2b, x3b
  REAL :: dh2o(0:np), dcont(0:np), dco2(0:np), do3(0:np)
  REAL :: dh2ob(0:np), dcontb(0:np), dco2b(0:np), do3b(0:np)
  REAL :: dn2o(0:np), dch4(0:np)
  REAL :: df11(0:np), df12(0:np), df22(0:np)
  REAL :: th2o(6), tcon(3), tco2(6)
  REAL :: th2ob(6), tconb(3), tco2b(6)
  REAL :: tn2o(4), tch4(4), tcom(6)
  REAL :: tn2ob(4), tch4b(4), tcomb(6)
  REAL :: tf11, tf12, tf22
  REAL :: tf11b, tf12b, tf22b
  REAL :: blayer(0:np+1), blevel(0:np+1)
  REAL :: blayerb(0:np+1), blevelb(0:np+1)
  REAL :: cd(0:np+1), cu(0:np+1)
  REAL :: bd(0:np+1), bu(0:np+1)
  REAL :: bdb(0:np+1), bub(0:np+1)
  REAL :: ad(0:np+1), au(0:np+1)
  REAL :: bs, dbs, rflxs
  REAL :: dp(0:np)
  REAL :: taant
  REAL :: trant, tranal
  REAL :: trantb, tranalb
  REAL :: transfc(0:np+1), trantcr(0:np+1), trantca(0:np+1)
  REAL :: transfcb(0:np+1)
  REAL :: flau(0:np+1), flad(0:np+1)
  REAL :: flcu(0:np+1), flcd(0:np+1)
  REAL :: flxu(0:np+1), flxd(0:np+1)
  REAL :: flxub(0:np+1), flxdb(0:np+1)
  REAL :: taerlyr(0:np)
  REAL :: taerlyrb(0:np)
!OVERCAST
  INTEGER :: ncld(3)
  INTEGER :: icx(0:np)
!OVERCAST
  INTEGER :: idx, rc
  INTEGER :: i, j, k, l, ip, iw, ibn, ik, iq, isb, k1, k2, ne
  REAL :: enn(0:np)
  REAL :: ennb(0:np)
  REAL :: cldhi, cldmd, cldlw, tcldlyr(0:np), fclr, fclr_above
  REAL :: cldhib, cldmdb, cldlwb, tcldlyrb(0:np), fclrb, fclr_aboveb
  REAL :: x, xx, yy, p1, a1, b1, fk1, a2, b2, fk2
  REAL :: xxb, yyb
  REAL :: w1, ff
  REAL :: ffb
  LOGICAL :: oznbnd, co2bnd, h2otable, conbnd, n2obnd
  LOGICAL :: ch4bnd, combnd, f11bnd, f12bnd, f22bnd, b10bnd
  LOGICAL :: do_aerosol
!Temp arrays and variables for consolidation of tables
  INTEGER, PARAMETER :: max_num_tables=17
  REAL :: exptbl(0:np, max_num_tables)
  REAL :: exptblb(0:np, max_num_tables)
  TYPE BAND_TABLE
      INTEGER :: start
      INTEGER :: end
  END TYPE BAND_TABLE
  TYPE(BAND_TABLE) :: h2oexp
  TYPE(BAND_TABLE) :: conexp
  TYPE(BAND_TABLE) :: co2exp
  TYPE(BAND_TABLE) :: n2oexp
  TYPE(BAND_TABLE) :: ch4exp
  TYPE(BAND_TABLE) :: comexp
  TYPE(BAND_TABLE) :: f11exp
  TYPE(BAND_TABLE) :: f12exp
  TYPE(BAND_TABLE) :: f22exp
!Variables for new getirtau routine
  REAL :: dp_pa(np)
  REAL :: taudiaglyr(np, 4)
  REAL :: fcld_col(np)
  REAL :: fcld_colb(np)
  REAL :: reff_col(np, 4)
  REAL :: reff_colb(np, 4)
  REAL :: cwc_col(np, 4)
  REAL :: cwc_colb(np, 4)
  REAL :: h2oexp_tmp(0:np, 5), conexp_tmp(0:np), co2exp_tmp(0:np, 6), &
& n2oexp_tmp(0:np, 2)
  REAL :: h2oexp_tmpb(0:np, 5), co2exp_tmpb(0:np, 6), n2oexp_tmpb(0:np, &
& 2)
  REAL, DIMENSION(m, np, nb) :: taua_dev_tmp, ssaa_dev_tmp, asya_dev_tmp
  REAL, DIMENSION(m, np, nb) :: taua_dev_tmpb, ssaa_dev_tmpb, &
& asya_dev_tmpb
  INTRINSIC MAX
  INTRINSIC EXP
  INTRINSIC MIN
  INTRINSIC ALOG
  INTEGER :: arg1
  INTEGER :: branch
  INTEGER :: ad_from
  INTEGER :: ad_count
  INTEGER :: i0
  REAL :: temp3
  REAL :: temp2
  REAL :: temp1
  REAL :: temp0
  REAL :: tempb6
  REAL :: tempb5
  REAL :: tempb4
  REAL :: tempb3
  REAL :: tempb2
  REAL :: tempb1
  REAL :: tempb0
  REAL :: tempb
  REAL :: temp
!BEGIN CALCULATIONS ...
  taua_dev_tmp = taua_dev
  ssaa_dev_tmp = ssaa_dev
  asya_dev_tmp = asya_dev
  i = 1
!do i=1,m
!-----compute layer pressure (pa) and layer temperature minus 250K (dt)
  DO k=1,np
    pa(k) = 0.5*(ple_dev(i, k+1)+ple_dev(i, k))*0.01
    dp(k) = (ple_dev(i, k+1)-ple_dev(i, k))*0.01
! dp in Pascals for getirtau
    dp_pa(k) = ple_dev(i, k+1) - ple_dev(i, k)
    dt(k) = ta_dev(i, k) - 250.0
!-----compute layer absorber amount
!     dh2o : water vapor amount (g/cm^2)
!     dcont: scaled water vapor amount for continuum absorption
!            (g/cm^2)
!     dco2 : co2 amount (cm-atm)stp
!     do3  : o3 amount (cm-atm)stp
!     dn2o : n2o amount (cm-atm)stp
!     dch4 : ch4 amount (cm-atm)stp
!     df11 : cfc11 amount (cm-atm)stp
!     df12 : cfc12 amount (cm-atm)stp
!     df22 : cfc22 amount (cm-atm)stp
!     the factor 1.02 is equal to 1000/980
!     factors 789 and 476 are for unit conversion
!     the factor 0.001618 is equal to 1.02/(.622*1013.25) 
!     the factor 6.081 is equal to 1800/296
    dh2o(k) = 1.02*wa_dev(i, k)*dp(k)
    do3(k) = 476.*oa_dev(i, k)*dp(k)
    dco2(k) = 789.*co2*dp(k)
    dch4(k) = 789.*ch4_dev(i, k)*dp(k)
    dn2o(k) = 789.*n2o_dev(i, k)*dp(k)
    df11(k) = 789.*cfc11_dev(i, k)*dp(k)
    df12(k) = 789.*cfc12_dev(i, k)*dp(k)
    df22(k) = 789.*cfc22_dev(i, k)*dp(k)
    IF (dh2o(k) .LT. 1.e-10) THEN
      dh2o(k) = 1.e-10
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
      dh2o(k) = dh2o(k)
    END IF
    IF (do3(k) .LT. 1.e-6) THEN
      do3(k) = 1.e-6
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
      do3(k) = do3(k)
    END IF
    IF (dco2(k) .LT. 1.e-4) THEN
      dco2(k) = 1.e-4
    ELSE
      dco2(k) = dco2(k)
    END IF
!-----compute scaled water vapor amount for h2o continuum absorption
!     following eq. (4.21).
    xx = pa(k)*0.001618*wa_dev(i, k)*wa_dev(i, k)*dp(k)
    dcont(k) = xx*EXP(1800./ta_dev(i, k)-6.081)
!-----Fill the reff, cwc, and fcld for the column
    fcld_col(k) = fcld_dev(i, k)
    DO l=1,4
      reff_col(k, l) = reff_dev(i, k, l)
      cwc_col(k, l) = cwc_dev(i, k, l)
    END DO
  END DO
  IF (ple_dev(i, 1)*0.01 .LT. 0.005) THEN
    CALL PUSHREAL4(dp(0))
    dp(0) = 0.005
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHREAL4(dp(0))
    dp(0) = ple_dev(i, 1)*0.01
    CALL PUSHCONTROL1B(1)
  END IF
  CALL PUSHREAL4(pa(0))
  pa(0) = 0.5*dp(0)
  dt(0) = ta_dev(i, 1) - 250.0
  dh2o(0) = 1.02*wa_dev(i, 1)*dp(0)
  do3(0) = 476.*oa_dev(i, 1)*dp(0)
  dco2(0) = 789.*co2*dp(0)
  dch4(0) = 789.*ch4_dev(i, 1)*dp(0)
  dn2o(0) = 789.*n2o_dev(i, 1)*dp(0)
  df11(0) = 789.*cfc11_dev(i, 1)*dp(0)
  df12(0) = 789.*cfc12_dev(i, 1)*dp(0)
  df22(0) = 789.*cfc22_dev(i, 1)*dp(0)
  IF (dh2o(0) .LT. 1.e-10) THEN
    dh2o(0) = 1.e-10
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    dh2o(0) = dh2o(0)
  END IF
  IF (do3(0) .LT. 1.e-6) THEN
    do3(0) = 1.e-6
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    do3(0) = do3(0)
  END IF
  IF (dco2(0) .LT. 1.e-4) THEN
    dco2(0) = 1.e-4
  ELSE
    dco2(0) = dco2(0)
  END IF
  xx = pa(0)*0.001618*wa_dev(i, 1)*wa_dev(i, 1)*dp(0)
  dcont(0) = xx*EXP(1800./ta_dev(i, 1)-6.081)
!-----the surface (np+1) is treated as a layer filled with black clouds.
!     transfc is the transmittance between the surface and a pressure
!     level.
!     trantcr is the clear-sky transmittance between the surface and a
!     pressure level.
!sfcem_dev(i) =0.0
  transfc(np+1) = 1.0
  CALL PUSHINTEGER4(ibn)
  ad_count = 1
!-----integration over spectral bands
  DO ibn=1,10
    IF (ibn .EQ. 10 .AND. (.NOT.trace)) THEN
      GOTO 100
    ELSE
!-----if h2otable, compute h2o (line) transmittance using table look-up.
!     if conbnd,   compute h2o (continuum) transmittance in bands 2-7.
!     if co2bnd,   compute co2 transmittance in band 3.
!     if oznbnd,   compute  o3 transmittance in band 5.
!     if n2obnd,   compute n2o transmittance in bands 6 and 7.
!     if ch4bnd,   compute ch4 transmittance in bands 6 and 7.
!     if combnd,   compute co2-minor transmittance in bands 4 and 5.
!     if f11bnd,   compute cfc11 transmittance in bands 4 and 5.
!     if f12bnd,   compute cfc12 transmittance in bands 4 and 6.
!     if f22bnd,   compute cfc22 transmittance in bands 4 and 6.
!     if b10bnd,   compute flux reduction due to n2o in band 10.
      h2otable = (ibn .EQ. 1 .OR. ibn .EQ. 2) .OR. ibn .EQ. 8
      conbnd = ibn .GE. 2 .AND. ibn .LE. 7
      co2bnd = ibn .EQ. 3
      oznbnd = ibn .EQ. 5
      n2obnd = ibn .EQ. 6 .OR. ibn .EQ. 7
      ch4bnd = ibn .EQ. 6 .OR. ibn .EQ. 7
      combnd = ibn .EQ. 4 .OR. ibn .EQ. 5
      f11bnd = ibn .EQ. 4 .OR. ibn .EQ. 5
      f12bnd = ibn .EQ. 4 .OR. ibn .EQ. 6
      f22bnd = ibn .EQ. 4 .OR. ibn .EQ. 6
      b10bnd = ibn .EQ. 10
      do_aerosol = na .GT. 0
      CALL PUSHREAL4ARRAY(exptbl, (np+1)*17)
      exptbl = 0.0
!-----Control packing of the new exponential tables by band
      SELECT CASE  (ibn) 
      CASE (2) 
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 1
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 1
        CALL PUSHCONTROL4B(7)
      CASE (3) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 7
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 9
        CALL PUSHCONTROL4B(6)
      CASE (4) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 7
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 7
        CALL PUSHINTEGER4(comexp%start)
        comexp%start = 8
        CALL PUSHINTEGER4(comexp%end)
        comexp%end = 13
        CALL PUSHINTEGER4(f11exp%start)
        f11exp%start = 14
        CALL PUSHINTEGER4(f11exp%end)
        f11exp%end = 14
        CALL PUSHINTEGER4(f12exp%start)
        f12exp%start = 15
        CALL PUSHINTEGER4(f12exp%end)
        f12exp%end = 15
        CALL PUSHINTEGER4(f22exp%start)
        f22exp%start = 16
        CALL PUSHINTEGER4(f22exp%end)
        f22exp%end = 16
        CALL PUSHCONTROL4B(5)
      CASE (5) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 7
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 7
        CALL PUSHINTEGER4(comexp%start)
        comexp%start = 8
        CALL PUSHINTEGER4(comexp%end)
        comexp%end = 13
        CALL PUSHINTEGER4(f11exp%start)
        f11exp%start = 14
        CALL PUSHINTEGER4(f11exp%end)
        f11exp%end = 14
        CALL PUSHCONTROL4B(4)
      CASE (6) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 7
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 7
        CALL PUSHINTEGER4(n2oexp%start)
        n2oexp%start = 8
        CALL PUSHINTEGER4(n2oexp%end)
        n2oexp%end = 11
        CALL PUSHINTEGER4(ch4exp%start)
        ch4exp%start = 12
        CALL PUSHINTEGER4(ch4exp%end)
        ch4exp%end = 15
        CALL PUSHINTEGER4(f12exp%start)
        f12exp%start = 16
        CALL PUSHINTEGER4(f12exp%end)
        f12exp%end = 16
        CALL PUSHINTEGER4(f22exp%start)
        f22exp%start = 17
        CALL PUSHINTEGER4(f22exp%end)
        f22exp%end = 17
        CALL PUSHCONTROL4B(3)
      CASE (7) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 7
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 7
        CALL PUSHINTEGER4(n2oexp%start)
        n2oexp%start = 8
        CALL PUSHINTEGER4(n2oexp%end)
        n2oexp%end = 11
        CALL PUSHINTEGER4(ch4exp%start)
        ch4exp%start = 12
        CALL PUSHINTEGER4(ch4exp%end)
        ch4exp%end = 15
        CALL PUSHCONTROL4B(2)
      CASE (9) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 6
        CALL PUSHCONTROL4B(1)
      CASE (10) 
        CALL PUSHINTEGER4(h2oexp%start)
        h2oexp%start = 1
        CALL PUSHINTEGER4(h2oexp%end)
        h2oexp%end = 5
        CALL PUSHINTEGER4(conexp%start)
        conexp%start = 6
        CALL PUSHINTEGER4(conexp%end)
        conexp%end = 6
        CALL PUSHINTEGER4(co2exp%start)
        co2exp%start = 7
        CALL PUSHINTEGER4(co2exp%end)
        co2exp%end = 12
        CALL PUSHINTEGER4(n2oexp%start)
        n2oexp%start = 13
        CALL PUSHINTEGER4(n2oexp%end)
        n2oexp%end = 14
        CALL PUSHCONTROL4B(0)
      CASE DEFAULT
        CALL PUSHCONTROL4B(8)
      END SELECT
!-----blayer is the spectrally integrated planck flux of the mean layer
!     temperature derived from eq. (3.11)
!     The fitting for the planck flux is valid for the range 160-345 K.
      DO k=1,np
        CALL PLANCK(ibn, cb, ta_dev(i, k), blayer(k))
      END DO
!-----Index "0" is the layer above the top of the atmosphere.
      blayer(0) = blayer(1)
      CALL PUSHREAL4(blevel(0))
      blevel(0) = blayer(1)
!-----Surface emission and reflectivity. See Section 9.
!     bs and dbs include the effect of surface emissivity.
      CALL PUSHREAL4(rflxs)
      CALL PUSHREAL4(dbs)
      CALL SFCFLUX(ibn, m, i, cb, dcb, ns, fs_dev, tg_dev, eg_dev, &
&            tv_dev, ev_dev, rv_dev, bs, dbs, rflxs)
      blayer(np+1) = bs
!------interpolate Planck function at model levels (linear in p)
      DO k=2,np
        CALL PUSHREAL4(blevel(k))
        blevel(k) = (blayer(k-1)*dp(k)+blayer(k)*dp(k-1))/(dp(k-1)+dp(k)&
&         )
      END DO
!-----Extrapolate blevel(1) from blayer(2) and blayer(1)
      CALL PUSHREAL4(blevel(1))
      blevel(1) = blayer(1) + (blayer(1)-blayer(2))*dp(1)/(dp(1)+dp(2))
      CALL PUSHREAL4(blevel(0))
      blevel(0) = blevel(1)
!-----If the surface air temperature tb is known, compute blevel(np+1)
      CALL PUSHREAL4(blevel(np+1))
      CALL PLANCK(ibn, cb, tb_dev(i), blevel(np+1))
!-----if not, extrapolate blevel(np+1) from blayer(np-1) and blayer(np)
!        blevel(np+1)=blayer(np)+(blayer(np)-blayer(np-1))&
!                      *dp(np)/(dp(np)+dp(np-1))
!-----Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
!     NOTE: dp_pa is only dims(1:np) as the 0'th level isn't needed in getirtau.
!           Plus, the pressures in getirtau *MUST* be in Pascals.
!     Slots for reff, hydrometeors and tauall are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
      CALL PUSHREAL4ARRAY(enn, np + 1)
      CALL PUSHREAL4ARRAY(tcldlyr, np + 1)
      CALL GETIRTAU1(ibn, np, dp_pa, fcld_col, reff_col, cwc_col, &
&              tcldlyr, enn, aib_ir, awb_ir, aiw_ir, aww_ir, aig_ir, &
&              awg_ir, cons_grav)
!     The getirtau call returns results for a single band and column. Thus
!     we need to transfer the taudiag for that band/column to the overall array.
!MAT-- icx and ncld only used when overcast=.false.
!Overcast
      IF (overcastl .EQ. .false.) THEN
        DO k=0,np
          CALL PUSHINTEGER4(icx(k))
          icx(k) = k
        END DO
        CALL PUSHINTEGER4ARRAY(ncld, 3)
        CALL PUSHINTEGER4ARRAY(icx, np + 1)
        CALL MKICX(np, ict, icb, enn, icx, ncld)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!-----Compute optical thickness, single-scattering albedo and asymmetry
!     factor for a mixture of "na" aerosol types. Eqs. (7.1)-(7.3)
      IF (do_aerosol) THEN
        CALL PUSHREAL4(taerlyr(0))
        taerlyr(0) = 1.0
        DO k=1,np
!-----taerlyr is the aerosol diffuse transmittance
          CALL PUSHREAL4(taerlyr(k))
          taerlyr(k) = 1.0
          IF (taua_dev_tmp(i, k, ibn) .GT. 0.001) THEN
            IF (ssaa_dev_tmp(i, k, ibn) .GT. 0.001) THEN
              CALL PUSHREAL4(asya_dev_tmp(i, k, ibn))
              asya_dev_tmp(i, k, ibn) = asya_dev_tmp(i, k, ibn)/&
&               ssaa_dev_tmp(i, k, ibn)
              CALL PUSHREAL4(ssaa_dev_tmp(i, k, ibn))
              ssaa_dev_tmp(i, k, ibn) = ssaa_dev_tmp(i, k, ibn)/&
&               taua_dev_tmp(i, k, ibn)
!-----Parameterization of aerosol scattering following Eqs. (6.11)
!     and (6.12). 
              CALL PUSHREAL4(ff)
              ff = .5 + (.3739+(0.0076+0.1185*asya_dev_tmp(i, k, ibn))*&
&               asya_dev_tmp(i, k, ibn))*asya_dev_tmp(i, k, ibn)
              CALL PUSHREAL4(taua_dev_tmp(i, k, ibn))
              taua_dev_tmp(i, k, ibn) = taua_dev_tmp(i, k, ibn)*(1.-&
&               ssaa_dev_tmp(i, k, ibn)*ff)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            taerlyr(k) = EXP(-(1.66*taua_dev_tmp(i, k, ibn)))
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!-----Compute the exponential terms (Eq. 8.21) at each layer due to
!     water vapor line absorption when k-distribution is used
      IF (.NOT.h2otable .AND. (.NOT.b10bnd)) THEN
        CALL PUSHREAL4ARRAY(exptbl(:, h2oexp%start:h2oexp%end), (np+1)*(&
&                     h2oexp%end-h2oexp%start+1))
        CALL H2OEXPS(ibn, np, dh2o, pa, dt, xkw, aw, bw, pm, mw, exptbl(&
&              :, h2oexp%start:h2oexp%end))
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!-----compute the exponential terms (Eq. 4.24) at each layer due to
!     water vapor continuum absorption.
!     ne is the number of terms used in each band to compute water 
!     vapor continuum transmittance (Table 9).
      CALL PUSHINTEGER4(ne)
      ne = 0
      IF (conbnd) THEN
        ne = 1
        IF (ibn .EQ. 3) ne = 3
        CALL PUSHREAL4ARRAY(exptbl(:, conexp%start:conexp%end), (np+1)*(&
&                     conexp%end-conexp%start+1))
        CALL CONEXPS(ibn, np, dcont, xke, exptbl(:, conexp%start:conexp%&
&              end))
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!----- for trace gases 
      IF (trace) THEN
!-----compute the exponential terms at each layer due to n2o absorption
        IF (n2obnd) THEN
          CALL PUSHREAL4ARRAY(exptbl(:, n2oexp%start:n2oexp%end), (np+1)&
&                       *(n2oexp%end-n2oexp%start+1))
          CALL N2OEXPS(ibn, np, dn2o, pa, dt, exptbl(:, n2oexp%start:&
&                n2oexp%end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!-----compute the exponential terms at each layer due to ch4 absorption
        IF (ch4bnd) THEN
          CALL PUSHREAL4ARRAY(exptbl(:, ch4exp%start:ch4exp%end), (np+1)&
&                       *(ch4exp%end-ch4exp%start+1))
          CALL CH4EXPS(ibn, np, dch4, pa, dt, exptbl(:, ch4exp%start:&
&                ch4exp%end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!-----Compute the exponential terms due to co2 minor absorption
        IF (combnd) THEN
          CALL PUSHREAL4ARRAY(exptbl(:, comexp%start:comexp%end), (np+1)&
&                       *(comexp%end-comexp%start+1))
          CALL COMEXPS(ibn, np, dco2, dt, exptbl(:, comexp%start:comexp%&
&                end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!-----Compute the exponential terms due to cfc11 absorption.
!     The values of the parameters are given in Table 7.
        IF (f11bnd) THEN
          a1 = 1.26610e-3
          b1 = 3.55940e-6
          fk1 = 1.89736e+1
          a2 = 8.19370e-4
          b2 = 4.67810e-6
          fk2 = 1.01487e+1
          CALL CFCEXPS(ibn, np, a1, b1, fk1, a2, b2, fk2, df11, dt, &
&                exptbl(:, f11exp%start:f11exp%end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!-----Compute the exponential terms due to cfc12 absorption.
        IF (f12bnd) THEN
          a1 = 8.77370e-4
          b1 = -5.88440e-6
          fk1 = 1.58104e+1
          a2 = 8.62000e-4
          b2 = -4.22500e-6
          fk2 = 3.70107e+1
          CALL CFCEXPS(ibn, np, a1, b1, fk1, a2, b2, fk2, df12, dt, &
&                exptbl(:, f12exp%start:f12exp%end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!-----Compute the exponential terms due to cfc22 absorption.
        IF (f22bnd) THEN
          a1 = 9.65130e-4
          b1 = 1.31280e-5
          fk1 = 6.18536e+0
          a2 = -3.00010e-5
          b2 = 5.25010e-7
          fk2 = 3.27912e+1
          CALL CFCEXPS(ibn, np, a1, b1, fk1, a2, b2, fk2, df22, dt, &
&                exptbl(:, f22exp%start:f22exp%end))
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!-----Compute the exponential terms at each layer in band 10 due to
!     h2o line and continuum, co2, and n2o absorption
        IF (b10bnd) THEN
          CALL PUSHREAL4ARRAY(n2oexp_tmp, (np+1)*2)
          CALL PUSHREAL4ARRAY(co2exp_tmp, (np+1)*6)
          CALL PUSHREAL4ARRAY(h2oexp_tmp, (np+1)*5)
          CALL B10EXPS(np, dh2o, dcont, dco2, dn2o, pa, dt, h2oexp_tmp, &
&                exptbl(:, conexp%start:conexp%end), co2exp_tmp, &
&                n2oexp_tmp)
          exptbl(:, h2oexp%start:h2oexp%end) = h2oexp_tmp
!              exptbl(:,conexp%start:conexp%end) = conexp_tmp
          exptbl(:, co2exp%start:co2exp%end) = co2exp_tmp
          exptbl(:, n2oexp%start:n2oexp%end) = n2oexp_tmp
          CALL PUSHCONTROL2B(0)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(2)
      END IF
!-----blayer(np+1) includes the effect of surface emissivity.
! ALT: this was undefined, check with Max if 0.0 is good value
      CALL PUSHREAL4(bu(0))
      bu(0) = 0.0
      CALL PUSHREAL4(bd(0))
      bd(0) = blayer(1)
      CALL PUSHREAL4(bu(np+1))
      bu(np+1) = blayer(np+1)
! ALT: this was undefined, check with Max if 0.0 is good value
! ALT: this was undefined, check with Max if 0.0 is good value
!-----do-loop 1500 is for computing upward (bu) and downward (bd)
!     emission of a layer following Eqs. (8.17), (8.18), (8.19).
!     Here, trant is the transmittance of the layer k2-1.
      DO k2=1,np+1
!-----for h2o line transmission
        IF (.NOT.h2otable) THEN
          th2o = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!-----for h2o continuum transmission
        CALL PUSHREAL4ARRAY(tcon, 3)
        tcon = 1.0
        x1 = 0.0
        x2 = 0.0
        x3 = 0.0
        CALL PUSHREAL4(trant)
        trant = 1.0
        IF (h2otable) THEN
!-----Compute water vapor transmittance using table look-up.
!     The following values are taken from Table 8.
!bdc
!              w1=-8.0
!              p1=-2.0
!              dwe=0.3
!              dpe=0.2
          IF (ibn .EQ. 1) THEN
            CALL PUSHREAL4(trant)
            CALL PUSHREAL4(x3)
            CALL PUSHREAL4(x2)
            CALL PUSHREAL4(x1)
            CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, x2&
&                 , x3, w11, p11, dwe, dpe, h11, h12, h13, trant)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ibn .EQ. 2) THEN
            CALL PUSHREAL4(trant)
            CALL PUSHREAL4(x3)
            CALL PUSHREAL4(x2)
            CALL PUSHREAL4(x1)
            CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, x2&
&                 , x3, w11, p11, dwe, dpe, h21, h22, h23, trant)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ibn .EQ. 8) THEN
            CALL PUSHREAL4(trant)
            CALL PUSHREAL4(x3)
            CALL PUSHREAL4(x2)
            CALL PUSHREAL4(x1)
            CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, x2&
&                 , x3, w11, p11, dwe, dpe, h81, h82, h83, trant)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!bdc
!-----for water vapor continuum absorption
          IF (conbnd) THEN
! Only the first exp
            CALL PUSHREAL4(tcon(1))
            tcon(1) = tcon(1)*exptbl(k2-1, conexp%start)
            CALL PUSHREAL4(trant)
            trant = trant*tcon(1)
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE IF (.NOT.b10bnd) THEN
!-----compute water vapor transmittance using k-distribution
          arg1 = k2 - 1
          CALL PUSHREAL4(trant)
          CALL PUSHREAL4ARRAY(tcon, 3)
          CALL PUSHREAL4ARRAY(th2o, 6)
          CALL H2OKDIS(ibn, np, arg1, fkw, gkw, ne, exptbl(:, h2oexp%&
&                start:h2oexp%end), exptbl(:, conexp%start:conexp%end), &
&                th2o, tcon, trant)
          CALL PUSHCONTROL2B(2)
        ELSE
          CALL PUSHCONTROL2B(3)
        END IF
        IF (co2bnd) THEN
!-----Compute co2 transmittance using table look-up method.
!     The following values are taken from Table 8.
!              w1=-4.0
!              p1=-2.0
!              dwe=0.3
!              dpe=0.2
          CALL PUSHREAL4(trant)
          CALL PUSHREAL4(x3)
          CALL PUSHREAL4(x2)
          CALL PUSHREAL4(x1)
          CALL TABLUP(nx1, nc1, dco2(k2-1), pa(k2-1), dt(k2-1), x1, x2, &
&               x3, w12, p12, dwe, dpe, c1, c2, c3, trant)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!-----Always use table look-up to compute o3 transmittance.
!     The following values are taken from Table 8.
        IF (oznbnd) THEN
!              w1=-6.0
!              p1=-2.0
!              dwe=0.3
!              dpe=0.2
          CALL PUSHREAL4(trant)
          CALL PUSHREAL4(x3)
          CALL PUSHREAL4(x2)
          CALL PUSHREAL4(x1)
          CALL TABLUP(nx1, no1, do3(k2-1), pa(k2-1), dt(k2-1), x1, x2, &
&               x3, w13, p13, dwe, dpe, oo1, oo2, oo3, trant)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!-----include aerosol effect
        IF (do_aerosol) THEN
          CALL PUSHREAL4(trant)
          trant = trant*taerlyr(k2-1)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
!-----Compute upward (bu) and downward (bd) emission of the layer k2-1
!     following Eqs.(8.17) and (8.18).
!     The effect of clouds on the transmission of a layer is taken
!     into account, following Eq. (8.19).
!     trant is the total transmittance of the layer k2-1.
        CALL PUSHREAL4(xx)
        xx = (1.-enn(k2-1))*trant
        IF (0.9999 .GT. xx) THEN
          CALL PUSHREAL4(yy)
          yy = xx
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHREAL4(yy)
          yy = 0.9999
          CALL PUSHCONTROL1B(1)
        END IF
        IF (0.00001 .LT. yy) THEN
          CALL PUSHCONTROL1B(0)
          yy = yy
        ELSE
          yy = 0.00001
          CALL PUSHCONTROL1B(1)
        END IF
        xx = (blevel(k2-1)-blevel(k2))/ALOG(yy)
        CALL PUSHREAL4(bd(k2-1))
        bd(k2-1) = (blevel(k2)-blevel(k2-1)*yy)/(1.0-yy) - xx
        CALL PUSHREAL4(bu(k2-1))
        bu(k2-1) = blevel(k2-1) + blevel(k2) - bd(k2-1)
      END DO
!-----initialize fluxes
      CALL PUSHREAL4ARRAY(flxd, np + 2)
      flxd = 0.0
!-----Compute upward and downward fluxes for each spectral band, ibn.
      DO k1=0,np
!-----initialization
!
!     cldlw, cldmd, and cldhi are the equivalent black-cloud fractions
!     of low, middle, and high troposphere.
!     tranal is the aerosol transmission function
        CALL PUSHREAL4(cldlw)
        cldlw = 0.0
        CALL PUSHREAL4(cldmd)
        cldmd = 0.0
        CALL PUSHREAL4(cldhi)
        cldhi = 0.0
        CALL PUSHREAL4(tranal)
        tranal = 1.0
!-----for h2o line transmission
        IF (.NOT.h2otable) THEN
          th2o = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!-----for h2o continuum transmission
        CALL PUSHREAL4ARRAY(tcon, 3)
        tcon = 1.0
!----- for trace gases
        IF (trace) THEN
!-----for n2o transmission using k-distribution method.
          IF (n2obnd) THEN
            tn2o = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!-----for ch4 transmission using k-distribution method.
          IF (ch4bnd) THEN
            tch4 = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!-----for co2-minor transmission using k-distribution method.
          IF (combnd) THEN
            tcom = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!-----for cfc-11 transmission using k-distribution method.
          IF (f11bnd) THEN
            tf11 = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!-----for cfc-12 transmission using k-distribution method.
          IF (f12bnd) THEN
            tf12 = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!-----for cfc-22 transmission when using k-distribution method.
          IF (f22bnd) THEN
            tf22 = 1.0
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!-----for the transmission in band 10 using k-distribution method.
          IF (b10bnd) THEN
            th2o = 1.0
            tco2 = 1.0
            tcon(1) = 1.0
            tn2o = 1.0
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
!----- end trace gases
        x1 = 0.0
        x2 = 0.0
        x3 = 0.0
!-----do-loop 3000 are for computing (a) transmittance, trant,
!     and (b) clear line-of-sight, fclr(k2), between levels k1 and k2.
        fclr_above = 1.0
        ad_from = k1 + 1
!MAT--Beginning of original 3000 loop
        DO k2=ad_from,np+1
          CALL PUSHREAL4(trant)
          trant = 1.0
          IF (h2otable) THEN
!-----Compute water vapor transmittance using table look-up.
!     The following values are taken from Table 8.
!                 w1=-8.0
!                 p1=-2.0
!                 dwe=0.3
!                 dpe=0.2
            IF (ibn .EQ. 1) THEN
              CALL PUSHREAL4(trant)
              CALL PUSHREAL4(x3)
              CALL PUSHREAL4(x2)
              CALL PUSHREAL4(x1)
              CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, &
&                   x2, x3, w11, p11, dwe, dpe, h11, h12, h13, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (ibn .EQ. 2) THEN
              CALL PUSHREAL4(trant)
              CALL PUSHREAL4(x3)
              CALL PUSHREAL4(x2)
              CALL PUSHREAL4(x1)
              CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, &
&                   x2, x3, w11, p11, dwe, dpe, h21, h22, h23, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (ibn .EQ. 8) THEN
              CALL PUSHREAL4(trant)
              CALL PUSHREAL4(x3)
              CALL PUSHREAL4(x2)
              CALL PUSHREAL4(x1)
              CALL TABLUP(nx1, nh1, dh2o(k2-1), pa(k2-1), dt(k2-1), x1, &
&                   x2, x3, w11, p11, dwe, dpe, h81, h82, h83, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            IF (conbnd) THEN
! Only the first exp
              CALL PUSHREAL4(tcon(1))
              tcon(1) = tcon(1)*exptbl(k2-1, conexp%start)
              CALL PUSHREAL4(trant)
              trant = trant*tcon(1)
              CALL PUSHCONTROL2B(0)
            ELSE
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE IF (.NOT.b10bnd) THEN
!-----compute water vapor transmittance using k-distribution
            arg1 = k2 - 1
            CALL PUSHREAL4(trant)
            CALL PUSHREAL4ARRAY(tcon, 3)
            CALL PUSHREAL4ARRAY(th2o, 6)
            CALL H2OKDIS(ibn, np, arg1, fkw, gkw, ne, exptbl(:, h2oexp%&
&                  start:h2oexp%end), exptbl(:, conexp%start:conexp%end)&
&                  , th2o, tcon, trant)
            CALL PUSHCONTROL2B(2)
          ELSE
            CALL PUSHCONTROL2B(3)
          END IF
          IF (co2bnd) THEN
!-----Compute co2 transmittance using table look-up method.
!     The following values are taken from Table 8.
!                 w1=-4.0
!                 p1=-2.0
!                 dwe=0.3
!                 dpe=0.2
            CALL PUSHREAL4(trant)
            CALL PUSHREAL4(x3)
            CALL PUSHREAL4(x2)
            CALL PUSHREAL4(x1)
            CALL TABLUP(nx1, nc1, dco2(k2-1), pa(k2-1), dt(k2-1), x1, x2&
&                 , x3, w12, p12, dwe, dpe, c1, c2, c3, trant)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!-----Always use table look-up to compute o3 transmittance.
!     The following values are taken from Table 8.
          IF (oznbnd) THEN
!                 w1=-6.0
!                 p1=-2.0
!                 dwe=0.3
!                 dpe=0.2
            CALL PUSHREAL4(trant)
            CALL PUSHREAL4(x3)
            CALL PUSHREAL4(x2)
            CALL PUSHREAL4(x1)
            CALL TABLUP(nx1, no1, do3(k2-1), pa(k2-1), dt(k2-1), x1, x2&
&                 , x3, w13, p13, dwe, dpe, oo1, oo2, oo3, trant)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!------ for trace gases 
          IF (trace) THEN
!-----compute n2o transmittance using k-distribution method
            IF (n2obnd) THEN
              arg1 = k2 - 1
              CALL PUSHREAL4(trant)
              CALL PUSHREAL4ARRAY(tn2o, 4)
              CALL N2OKDIS(ibn, np, arg1, exptbl(:, n2oexp%start:n2oexp%&
&                    end), tn2o, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
!-----compute ch4 transmittance using k-distribution method
            IF (ch4bnd) THEN
              arg1 = k2 - 1
              CALL PUSHREAL4(trant)
              CALL PUSHREAL4ARRAY(tch4, 4)
              CALL CH4KDIS(ibn, np, arg1, exptbl(:, ch4exp%start:ch4exp%&
&                    end), tch4, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
!-----compute co2-minor transmittance using k-distribution method
            IF (combnd) THEN
              arg1 = k2 - 1
              CALL PUSHREAL4(trant)
              CALL PUSHREAL4ARRAY(tcom, 6)
              CALL COMKDIS(ibn, np, arg1, exptbl(:, comexp%start:comexp%&
&                    end), tcom, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
!-----compute cfc11 transmittance using k-distribution method
            IF (f11bnd) THEN
              arg1 = k2 - 1
              CALL PUSHREAL4(trant)
              CALL PUSHREAL4(tf11)
              CALL CFCKDIS(np, arg1, exptbl(:, f11exp%start:f11exp%end)&
&                    , tf11, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
!-----compute cfc12 transmittance using k-distribution method
            IF (f12bnd) THEN
              arg1 = k2 - 1
              CALL PUSHREAL4(trant)
              CALL PUSHREAL4(tf12)
              CALL CFCKDIS(np, arg1, exptbl(:, f12exp%start:f12exp%end)&
&                    , tf12, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
!-----compute cfc22 transmittance using k-distribution method
            IF (f22bnd) THEN
              arg1 = k2 - 1
              CALL PUSHREAL4(trant)
              CALL PUSHREAL4(tf22)
              CALL CFCKDIS(np, arg1, exptbl(:, f22exp%start:f22exp%end)&
&                    , tf22, trant)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
!-----Compute transmittance in band 10 using k-distribution method.
!     For band 10, trant is the change in transmittance due to n2o 
!     absorption.
            IF (b10bnd) THEN
              arg1 = k2 - 1
              CALL PUSHREAL4ARRAY(tn2o, 4)
              CALL PUSHREAL4ARRAY(tco2, 6)
              CALL PUSHREAL4ARRAY(tcon, 3)
              CALL PUSHREAL4ARRAY(th2o, 6)
              CALL B10KDIS(np, arg1, exptbl(:, h2oexp%start:h2oexp%end)&
&                    , exptbl(:, conexp%start:conexp%end), exptbl(:, &
&                    co2exp%start:co2exp%end), exptbl(:, n2oexp%start:&
&                    n2oexp%end), th2o, tcon, tco2, tn2o, trant)
              CALL PUSHCONTROL2B(0)
            ELSE
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE
            CALL PUSHCONTROL2B(2)
          END IF
!-----   end trace gases
!-----include aerosol effect
          IF (do_aerosol) THEN
            CALL PUSHREAL4(tranal)
            tranal = tranal*taerlyr(k2-1)
            CALL PUSHREAL4(trant)
            trant = trant*tranal
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!----- cloud overlapping 
!OVERCAST
          IF (overcastl .EQ. .false.) THEN
            IF (enn(k2-1) .GE. 0.001) THEN
              CALL PUSHREAL4(cldlw)
              CALL PUSHREAL4(cldmd)
              CALL PUSHREAL4(cldhi)
              CALL CLDOVLP(np, k1, k2, ict, icb, icx, ncld, enn, tcldlyr&
&                    , cldhi, cldmd, cldlw)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREAL4(fclr)
            fclr = (1.0-cldhi)*(1.0-cldmd)*(1.0-cldlw)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL4(fclr)
            fclr = fclr_above*tcldlyr(k2-1)
            CALL PUSHCONTROL1B(1)
          END IF
!Overcast
!MAT--End of original 3000 loop
!-----do-loop 4000 is for computing upward and downward fluxes
!     for each spectral band
!     flau, flad: clear-sky aerosol-free upward and downward fluxes
!     flcu, flcd: clear-sky upward and downward fluxes
!     flxu, flxd: all-sky   upward and downward fluxes
!MAT--Beginning of original 4000 loop
!-----The first terms on the rhs of Eqs. (8.15) and (8.16)
          IF (k2 .EQ. k1 + 1 .AND. ibn .NE. 10) THEN
            flxd(k2) = flxd(k2) + bd(k1)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
!-----The summation terms on the rhs of Eqs. (8.15) and (8.16).
!     Also see Eqs. (5.4) and (5.5) for Band 10.
          CALL PUSHREAL4(xx)
          IF (k1 .EQ. 0) THEN
!mjs  bd(-1) is not defined
            xx = -(trant*bd(k1))
            CALL PUSHCONTROL1B(0)
          ELSE
            xx = trant*(bd(k1-1)-bd(k1))
            CALL PUSHCONTROL1B(1)
          END IF
          flxd(k2) = flxd(k2) + xx*fclr
!MAT--End of original 4000 loop
          CALL PUSHREAL4(fclr_above)
          fclr_above = fclr
        END DO
        CALL PUSHINTEGER4(ad_from)
!-----Here, fclr and trant are, respectively, the clear line-of-sight 
!     and the transmittance between k1 and the surface.
        CALL PUSHREAL4(transfc(k1))
        transfc(k1) = trant*fclr
!-----compute the partial derivative of fluxes with respect to
!     surface temperature (Eq. 3.12). 
!     Note: upward flux is negative, and so is dfdts.
        IF (k1 .GT. 0) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
!-----For surface emission.
      IF (.NOT.b10bnd) THEN
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      CALL PUSHINTEGER4(ibn)
      ad_count = ad_count + 1
    END IF
  END DO
  CALL PUSHCONTROL1B(0)
  CALL PUSHINTEGER4(ad_count)
  GOTO 110
 100 CALL PUSHCONTROL1B(1)
  CALL PUSHINTEGER4(ad_count)
 110 CALL POPINTEGER4(ad_count)
  DO i0=1,ad_count
    IF (i0 .EQ. 1) THEN
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        dh2ob = 0.0
        dcontb = 0.0
        n2oexp_tmpb = 0.0
        dtb = 0.0
        reff_colb = 0.0
        blayerb = 0.0
        tn2ob = 0.0
        transfcb = 0.0
        bdb = 0.0
        h2oexp_tmpb = 0.0
        tcomb = 0.0
        tf11b = 0.0
        tf12b = 0.0
        ssaa_dev_tmpb = 0.0
        trantb = 0.0
        bub = 0.0
        taua_dev_tmpb = 0.0
        th2ob = 0.0
        tcldlyrb = 0.0
        tch4b = 0.0
        tf22b = 0.0
        asya_dev_tmpb = 0.0
        ennb = 0.0
        fcld_colb = 0.0
        fclrb = 0.0
        cwc_colb = 0.0
        tco2b = 0.0
        co2exp_tmpb = 0.0
        taerlyrb = 0.0
        do3b = 0.0
        blevelb = 0.0
      ELSE
        dh2ob = 0.0
        dcontb = 0.0
        n2oexp_tmpb = 0.0
        dtb = 0.0
        reff_colb = 0.0
        blayerb = 0.0
        tn2ob = 0.0
        transfcb = 0.0
        bdb = 0.0
        h2oexp_tmpb = 0.0
        tcomb = 0.0
        tf11b = 0.0
        tf12b = 0.0
        ssaa_dev_tmpb = 0.0
        trantb = 0.0
        bub = 0.0
        taua_dev_tmpb = 0.0
        th2ob = 0.0
        tcldlyrb = 0.0
        tch4b = 0.0
        tf22b = 0.0
        asya_dev_tmpb = 0.0
        ennb = 0.0
        fcld_colb = 0.0
        fclrb = 0.0
        cwc_colb = 0.0
        tco2b = 0.0
        co2exp_tmpb = 0.0
        taerlyrb = 0.0
        do3b = 0.0
        blevelb = 0.0
      END IF
    ELSE
      flxub = 0.0
      flxdb = 0.0
      DO k=np+1,1,-1
        flxdb(k) = flxdb(k) + flxd_devb(i, k)
        flxub(k) = flxub(k) + flxu_devb(i, k)
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        DO k=np+1,1,-1
          flxdb(np+1) = flxdb(np+1) - rflxs*transfc(k)*flxub(k)
          transfcb(k) = transfcb(k) - rflxs*flxd(np+1)*flxub(k)
        END DO
        blayerb(np+1) = blayerb(np+1) - flxub(np+1)
        flxub(np+1) = 0.0
      END IF
      exptblb = 0.0
      DO k1=np,0,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) transfcb(k1) = transfcb(k1) - dbs*dfdts_devb(&
&           i, k1)
        CALL POPREAL4(transfc(k1))
        trantb = trantb + fclr*transfcb(k1)
        fclrb = fclrb + trant*transfcb(k1)
        transfcb(k1) = 0.0
        cldhib = 0.0
        tconb = 0.0
        tranalb = 0.0
        x1b = 0.0
        x2b = 0.0
        x3b = 0.0
        cldlwb = 0.0
        fclr_aboveb = 0.0
        cldmdb = 0.0
        CALL POPINTEGER4(ad_from)
        DO k2=np+1,ad_from,-1
          CALL POPREAL4(fclr_above)
          fclrb = fclrb + fclr_aboveb
          xxb = fclr*flxdb(k2)
          fclrb = fclrb + xx*flxdb(k2)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            trantb = trantb - bd(k1)*xxb
            bdb(k1) = bdb(k1) - trant*xxb
          ELSE
            trantb = trantb + (bd(k1-1)-bd(k1))*xxb
            bdb(k1-1) = bdb(k1-1) + trant*xxb
            bdb(k1) = bdb(k1) - trant*xxb
          END IF
          xx = trant*(bu(k2-1)-bu(k2))
          xxb = fclr*flxub(k1)
          fclrb = fclrb + xx*flxub(k1)
          CALL POPREAL4(xx)
          trantb = trantb + (bu(k2-1)-bu(k2))*xxb
          bub(k2-1) = bub(k2-1) + trant*xxb
          bub(k2) = bub(k2) - trant*xxb
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            bdb(k1) = bdb(k1) + flxdb(k2)
            bub(k1) = bub(k1) - flxub(k1)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4(fclr)
            cldhib = cldhib - (1.0-cldmd)*(1.0-cldlw)*fclrb
            cldmdb = cldmdb - (1.0-cldhi)*(1.0-cldlw)*fclrb
            cldlwb = cldlwb - (1.0-cldhi)*(1.0-cldmd)*fclrb
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL4(cldhi)
              CALL POPREAL4(cldmd)
              CALL POPREAL4(cldlw)
              CALL CLDOVLP_B(np, k1, k2, ict, icb, icx, ncld, enn, ennb&
&                      , tcldlyr, tcldlyrb, cldhi, cldhib, cldmd, cldmdb&
&                      , cldlw, cldlwb)
            END IF
            fclr_aboveb = 0.0
          ELSE
            CALL POPREAL4(fclr)
            fclr_aboveb = tcldlyr(k2-1)*fclrb
            tcldlyrb(k2-1) = tcldlyrb(k2-1) + fclr_above*fclrb
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4(trant)
            tranalb = tranalb + trant*trantb
            trantb = tranal*trantb
            CALL POPREAL4(tranal)
            taerlyrb(k2-1) = taerlyrb(k2-1) + tranal*tranalb
            tranalb = taerlyr(k2-1)*tranalb
          END IF
          CALL POPCONTROL2B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL4ARRAY(th2o, 6)
            CALL POPREAL4ARRAY(tcon, 3)
            CALL POPREAL4ARRAY(tco2, 6)
            CALL POPREAL4ARRAY(tn2o, 4)
            CALL B10KDIS_B(np, arg1, exptbl(:, h2oexp%start:h2oexp%end)&
&                    , exptblb(:, h2oexp%start:h2oexp%end), exptbl(:, &
&                    conexp%start:conexp%end), exptblb(:, conexp%start:&
&                    conexp%end), exptbl(:, co2exp%start:co2exp%end), &
&                    exptblb(:, co2exp%start:co2exp%end), exptbl(:, &
&                    n2oexp%start:n2oexp%end), exptblb(:, n2oexp%start:&
&                    n2oexp%end), th2o, th2ob, tcon, tconb, tco2, tco2b&
&                    , tn2o, tn2ob, trant, trantb)
            trantb = 0.0
          ELSE IF (branch .NE. 1) THEN
            GOTO 120
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL4(tf22)
            CALL POPREAL4(trant)
            CALL CFCKDIS_B(np, arg1, exptbl(:, f22exp%start:f22exp%end)&
&                    , exptblb(:, f22exp%start:f22exp%end), tf22, tf22b&
&                    , trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL4(tf12)
            CALL POPREAL4(trant)
            CALL CFCKDIS_B(np, arg1, exptbl(:, f12exp%start:f12exp%end)&
&                    , exptblb(:, f12exp%start:f12exp%end), tf12, tf12b&
&                    , trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL4(tf11)
            CALL POPREAL4(trant)
            CALL CFCKDIS_B(np, arg1, exptbl(:, f11exp%start:f11exp%end)&
&                    , exptblb(:, f11exp%start:f11exp%end), tf11, tf11b&
&                    , trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL4ARRAY(tcom, 6)
            CALL POPREAL4(trant)
            CALL COMKDIS_B(ibn, np, arg1, exptbl(:, comexp%start:comexp%&
&                    end), exptblb(:, comexp%start:comexp%end), tcom, &
&                    tcomb, trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL4ARRAY(tch4, 4)
            CALL POPREAL4(trant)
            CALL CH4KDIS_B(ibn, np, arg1, exptbl(:, ch4exp%start:ch4exp%&
&                    end), exptblb(:, ch4exp%start:ch4exp%end), tch4, &
&                    tch4b, trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            arg1 = k2 - 1
            CALL POPREAL4ARRAY(tn2o, 4)
            CALL POPREAL4(trant)
            CALL N2OKDIS_B(ibn, np, arg1, exptbl(:, n2oexp%start:n2oexp%&
&                    end), exptblb(:, n2oexp%start:n2oexp%end), tn2o, &
&                    tn2ob, trant, trantb)
          END IF
 120      CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4(x1)
            CALL POPREAL4(x2)
            CALL POPREAL4(x3)
            CALL POPREAL4(trant)
            CALL TABLUP_B(nx1, no1, do3(k2-1), do3b(k2-1), pa(k2-1), dt(&
&                   k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w13, &
&                   p13, dwe, dpe, oo1, oo2, oo3, trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4(x1)
            CALL POPREAL4(x2)
            CALL POPREAL4(x3)
            CALL POPREAL4(trant)
            dco2b = 0.0
            CALL TABLUP_B(nx1, nc1, dco2(k2-1), dco2b(k2-1), pa(k2-1), &
&                   dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w12&
&                   , p12, dwe, dpe, c1, c2, c3, trant, trantb)
          END IF
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              CALL POPREAL4(trant)
              tconb(1) = tconb(1) + trant*trantb
              trantb = tcon(1)*trantb
              CALL POPREAL4(tcon(1))
              exptblb(k2-1, conexp%start) = exptblb(k2-1, conexp%start) &
&               + tcon(1)*tconb(1)
              tconb(1) = exptbl(k2-1, conexp%start)*tconb(1)
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL4(x1)
              CALL POPREAL4(x2)
              CALL POPREAL4(x3)
              CALL POPREAL4(trant)
              CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1)&
&                     , dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, &
&                     w11, p11, dwe, dpe, h81, h82, h83, trant, trantb)
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL4(x1)
              CALL POPREAL4(x2)
              CALL POPREAL4(x3)
              CALL POPREAL4(trant)
              CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1)&
&                     , dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, &
&                     w11, p11, dwe, dpe, h21, h22, h23, trant, trantb)
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL4(x1)
              CALL POPREAL4(x2)
              CALL POPREAL4(x3)
              CALL POPREAL4(trant)
              CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1)&
&                     , dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, &
&                     w11, p11, dwe, dpe, h11, h12, h13, trant, trantb)
            END IF
          ELSE IF (branch .EQ. 2) THEN
            arg1 = k2 - 1
            CALL POPREAL4ARRAY(th2o, 6)
            CALL POPREAL4ARRAY(tcon, 3)
            CALL POPREAL4(trant)
            CALL H2OKDIS_B(ibn, np, arg1, fkw, gkw, ne, exptbl(:, h2oexp&
&                    %start:h2oexp%end), exptblb(:, h2oexp%start:h2oexp%&
&                    end), exptbl(:, conexp%start:conexp%end), exptblb(:&
&                    , conexp%start:conexp%end), th2o, th2ob, tcon, &
&                    tconb, trant, trantb)
          END IF
          CALL POPREAL4(trant)
          trantb = 0.0
          fclrb = 0.0
        END DO
        CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          tn2ob = 0.0
          th2ob = 0.0
          tco2b = 0.0
        ELSE IF (branch .NE. 1) THEN
          GOTO 130
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tf22b = 0.0
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tf12b = 0.0
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tf11b = 0.0
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tcomb = 0.0
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tch4b = 0.0
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) tn2ob = 0.0
 130    CALL POPREAL4ARRAY(tcon, 3)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) th2ob = 0.0
        CALL POPREAL4(tranal)
        CALL POPREAL4(cldhi)
        CALL POPREAL4(cldmd)
        CALL POPREAL4(cldlw)
      END DO
      CALL POPREAL4ARRAY(flxd, np + 2)
      DO k2=np+1,1,-1
        bdb(k2-1) = bdb(k2-1) - bub(k2-1)
        tempb5 = bdb(k2-1)/(1.0-yy)
        xxb = -bdb(k2-1)
        temp3 = ALOG(yy)
        tempb6 = xxb/temp3
        CALL POPREAL4(bu(k2-1))
        blevelb(k2-1) = blevelb(k2-1) + bub(k2-1)
        blevelb(k2) = blevelb(k2) + tempb5 + bub(k2-1)
        bub(k2-1) = 0.0
        CALL POPREAL4(bd(k2-1))
        blevelb(k2-1) = blevelb(k2-1) + tempb6 - yy*tempb5
        yyb = ((blevel(k2)-blevel(k2-1)*yy)/(1.0-yy)-blevel(k2-1))*&
&         tempb5 - (blevel(k2-1)-blevel(k2))*tempb6/(temp3*yy)
        bdb(k2-1) = 0.0
        blevelb(k2) = blevelb(k2) - tempb6
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) yyb = 0.0
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL4(yy)
          xxb = yyb
        ELSE
          CALL POPREAL4(yy)
          xxb = 0.0
        END IF
        CALL POPREAL4(xx)
        ennb(k2-1) = ennb(k2-1) - trant*xxb
        trantb = trantb + (1.-enn(k2-1))*xxb
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          CALL POPREAL4(trant)
          taerlyrb(k2-1) = taerlyrb(k2-1) + trant*trantb
          trantb = taerlyr(k2-1)*trantb
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL4(x1)
          CALL POPREAL4(x2)
          CALL POPREAL4(x3)
          CALL POPREAL4(trant)
          x1b = 0.0
          x2b = 0.0
          x3b = 0.0
          CALL TABLUP_B(nx1, no1, do3(k2-1), do3b(k2-1), pa(k2-1), dt(k2&
&                 -1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w13, p13, &
&                 dwe, dpe, oo1, oo2, oo3, trant, trantb)
        ELSE
          x1b = 0.0
          x2b = 0.0
          x3b = 0.0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL4(x1)
          CALL POPREAL4(x2)
          CALL POPREAL4(x3)
          CALL POPREAL4(trant)
          dco2b = 0.0
          CALL TABLUP_B(nx1, nc1, dco2(k2-1), dco2b(k2-1), pa(k2-1), dt(&
&                 k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w12, p12&
&                 , dwe, dpe, c1, c2, c3, trant, trantb)
        END IF
        CALL POPCONTROL2B(branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            tconb = 0.0
            CALL POPREAL4(trant)
            tconb(1) = tconb(1) + trant*trantb
            trantb = tcon(1)*trantb
            CALL POPREAL4(tcon(1))
            exptblb(k2-1, conexp%start) = exptblb(k2-1, conexp%start) + &
&             tcon(1)*tconb(1)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4(x1)
            CALL POPREAL4(x2)
            CALL POPREAL4(x3)
            CALL POPREAL4(trant)
            CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1), &
&                   dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w11&
&                   , p11, dwe, dpe, h81, h82, h83, trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4(x1)
            CALL POPREAL4(x2)
            CALL POPREAL4(x3)
            CALL POPREAL4(trant)
            CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1), &
&                   dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w11&
&                   , p11, dwe, dpe, h21, h22, h23, trant, trantb)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4(x1)
            CALL POPREAL4(x2)
            CALL POPREAL4(x3)
            CALL POPREAL4(trant)
            CALL TABLUP_B(nx1, nh1, dh2o(k2-1), dh2ob(k2-1), pa(k2-1), &
&                   dt(k2-1), dtb(k2-1), x1, x1b, x2, x2b, x3, x3b, w11&
&                   , p11, dwe, dpe, h11, h12, h13, trant, trantb)
          END IF
        ELSE IF (branch .EQ. 2) THEN
          arg1 = k2 - 1
          CALL POPREAL4ARRAY(th2o, 6)
          CALL POPREAL4ARRAY(tcon, 3)
          CALL POPREAL4(trant)
          tconb = 0.0
          CALL H2OKDIS_B(ibn, np, arg1, fkw, gkw, ne, exptbl(:, h2oexp%&
&                  start:h2oexp%end), exptblb(:, h2oexp%start:h2oexp%end&
&                  ), exptbl(:, conexp%start:conexp%end), exptblb(:, &
&                  conexp%start:conexp%end), th2o, th2ob, tcon, tconb, &
&                  trant, trantb)
        END IF
        CALL POPREAL4(trant)
        CALL POPREAL4ARRAY(tcon, 3)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) th2ob = 0.0
        trantb = 0.0
      END DO
      CALL POPREAL4(bu(np+1))
      blayerb(np+1) = blayerb(np+1) + bub(np+1)
      bub(np+1) = 0.0
      CALL POPREAL4(bd(0))
      blayerb(1) = blayerb(1) + bdb(0)
      bdb(0) = 0.0
      CALL POPREAL4(bu(0))
      bub(0) = 0.0
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        n2oexp_tmpb = n2oexp_tmpb + exptblb(:, n2oexp%start:n2oexp%end)
        exptblb(:, n2oexp%start:n2oexp%end) = 0.0
        co2exp_tmpb = co2exp_tmpb + exptblb(:, co2exp%start:co2exp%end)
        exptblb(:, co2exp%start:co2exp%end) = 0.0
        h2oexp_tmpb = h2oexp_tmpb + exptblb(:, h2oexp%start:h2oexp%end)
        exptblb(:, h2oexp%start:h2oexp%end) = 0.0
        CALL POPREAL4ARRAY(h2oexp_tmp, (np+1)*5)
        CALL POPREAL4ARRAY(co2exp_tmp, (np+1)*6)
        CALL POPREAL4ARRAY(n2oexp_tmp, (np+1)*2)
        CALL B10EXPS_B(np, dh2o, dh2ob, dcont, dcontb, dco2, dn2o, pa, &
&                dt, dtb, h2oexp_tmp, h2oexp_tmpb, exptbl(:, conexp%&
&                start:conexp%end), exptblb(:, conexp%start:conexp%end)&
&                , co2exp_tmp, co2exp_tmpb, n2oexp_tmp, n2oexp_tmpb)
      ELSE IF (branch .NE. 1) THEN
        GOTO 140
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        a1 = 9.65130e-4
        a2 = -3.00010e-5
        b1 = 1.31280e-5
        b2 = 5.25010e-7
        fk1 = 6.18536e+0
        fk2 = 3.27912e+1
        CALL CFCEXPS_B(ibn, np, a1, b1, fk1, a2, b2, fk2, df22, dt, dtb&
&                , exptbl(:, f22exp%start:f22exp%end), exptblb(:, f22exp&
&                %start:f22exp%end))
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        a1 = 8.77370e-4
        a2 = 8.62000e-4
        b1 = -5.88440e-6
        b2 = -4.22500e-6
        fk1 = 1.58104e+1
        fk2 = 3.70107e+1
        CALL CFCEXPS_B(ibn, np, a1, b1, fk1, a2, b2, fk2, df12, dt, dtb&
&                , exptbl(:, f12exp%start:f12exp%end), exptblb(:, f12exp&
&                %start:f12exp%end))
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        a1 = 1.26610e-3
        a2 = 8.19370e-4
        b1 = 3.55940e-6
        b2 = 4.67810e-6
        fk1 = 1.89736e+1
        fk2 = 1.01487e+1
        CALL CFCEXPS_B(ibn, np, a1, b1, fk1, a2, b2, fk2, df11, dt, dtb&
&                , exptbl(:, f11exp%start:f11exp%end), exptblb(:, f11exp&
&                %start:f11exp%end))
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4ARRAY(exptbl(:, comexp%start:comexp%end), (np+1)*(&
&                    comexp%end-comexp%start+1))
        CALL COMEXPS_B(ibn, np, dco2, dt, dtb, exptbl(:, comexp%start:&
&                comexp%end), exptblb(:, comexp%start:comexp%end))
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4ARRAY(exptbl(:, ch4exp%start:ch4exp%end), (np+1)*(&
&                    ch4exp%end-ch4exp%start+1))
        CALL CH4EXPS_B(ibn, np, dch4, pa, dt, dtb, exptbl(:, ch4exp%&
&                start:ch4exp%end), exptblb(:, ch4exp%start:ch4exp%end))
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4ARRAY(exptbl(:, n2oexp%start:n2oexp%end), (np+1)*(&
&                    n2oexp%end-n2oexp%start+1))
        CALL N2OEXPS_B(ibn, np, dn2o, pa, dt, dtb, exptbl(:, n2oexp%&
&                start:n2oexp%end), exptblb(:, n2oexp%start:n2oexp%end))
      END IF
 140  CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4ARRAY(exptbl(:, conexp%start:conexp%end), (np+1)*(&
&                    conexp%end-conexp%start+1))
        CALL CONEXPS_B(ibn, np, dcont, dcontb, xke, exptbl(:, conexp%&
&                start:conexp%end), exptblb(:, conexp%start:conexp%end))
      END IF
      CALL POPINTEGER4(ne)
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4ARRAY(exptbl(:, h2oexp%start:h2oexp%end), (np+1)*(&
&                    h2oexp%end-h2oexp%start+1))
        CALL H2OEXPS_B(ibn, np, dh2o, dh2ob, pa, dt, dtb, xkw, aw, bw, &
&                pm, mw, exptbl(:, h2oexp%start:h2oexp%end), exptblb(:, &
&                h2oexp%start:h2oexp%end))
        exptblb(:, h2oexp%start:h2oexp%end) = 0.0
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO k=np,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            taua_dev_tmpb(i, k, ibn) = taua_dev_tmpb(i, k, ibn) - EXP(-(&
&             1.66*taua_dev_tmp(i, k, ibn)))*1.66*taerlyrb(k)
            taerlyrb(k) = 0.0
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL4(taua_dev_tmp(i, k, ibn))
              tempb1 = taua_dev_tmp(i, k, ibn)*taua_dev_tmpb(i, k, ibn)
              ssaa_dev_tmpb(i, k, ibn) = ssaa_dev_tmpb(i, k, ibn) - ff*&
&               tempb1
              ffb = -(ssaa_dev_tmp(i, k, ibn)*tempb1)
              taua_dev_tmpb(i, k, ibn) = (1.-ssaa_dev_tmp(i, k, ibn)*ff)&
&               *taua_dev_tmpb(i, k, ibn)
              CALL POPREAL4(ff)
              tempb2 = asya_dev_tmp(i, k, ibn)*ffb
              temp2 = 0.1185*asya_dev_tmp(i, k, ibn) + 0.0076
              asya_dev_tmpb(i, k, ibn) = asya_dev_tmpb(i, k, ibn) + (&
&               temp2*asya_dev_tmp(i, k, ibn)+.3739)*ffb + (temp2+&
&               asya_dev_tmp(i, k, ibn)*0.1185)*tempb2
              CALL POPREAL4(ssaa_dev_tmp(i, k, ibn))
              tempb3 = ssaa_dev_tmpb(i, k, ibn)/taua_dev_tmp(i, k, ibn)
              taua_dev_tmpb(i, k, ibn) = taua_dev_tmpb(i, k, ibn) - &
&               ssaa_dev_tmp(i, k, ibn)*tempb3/taua_dev_tmp(i, k, ibn)
              CALL POPREAL4(asya_dev_tmp(i, k, ibn))
              tempb4 = asya_dev_tmpb(i, k, ibn)/ssaa_dev_tmp(i, k, ibn)
              ssaa_dev_tmpb(i, k, ibn) = tempb3 - asya_dev_tmp(i, k, ibn&
&               )*tempb4/ssaa_dev_tmp(i, k, ibn)
              asya_dev_tmpb(i, k, ibn) = tempb4
            END IF
          END IF
          CALL POPREAL4(taerlyr(k))
          taerlyrb(k) = 0.0
        END DO
        CALL POPREAL4(taerlyr(0))
        taerlyrb(0) = 0.0
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPINTEGER4ARRAY(icx, np + 1)
        CALL POPINTEGER4ARRAY(ncld, 3)
        DO k=np,0,-1
          CALL POPINTEGER4(icx(k))
        END DO
      END IF
      CALL POPREAL4ARRAY(tcldlyr, np + 1)
      CALL POPREAL4ARRAY(enn, np + 1)
      CALL GETIRTAU1_B(ibn, np, dp_pa, fcld_col, fcld_colb, reff_col, &
&                reff_colb, cwc_col, cwc_colb, tcldlyr, tcldlyrb, enn, &
&                ennb, aib_ir, awb_ir, aiw_ir, aww_ir, aig_ir, awg_ir, &
&                cons_grav)
      CALL POPREAL4(blevel(np+1))
      CALL PLANCK_B(ibn, cb, tb_dev(i), tb_devb(i), blevel(np+1), &
&             blevelb(np+1))
      blevelb(np+1) = 0.0
      CALL POPREAL4(blevel(0))
      blevelb(1) = blevelb(1) + blevelb(0)
      blevelb(0) = 0.0
      CALL POPREAL4(blevel(1))
      tempb0 = dp(1)*blevelb(1)/(dp(1)+dp(2))
      blayerb(1) = blayerb(1) + tempb0 + blevelb(1)
      blayerb(2) = blayerb(2) - tempb0
      blevelb(1) = 0.0
      DO k=np,2,-1
        CALL POPREAL4(blevel(k))
        tempb = blevelb(k)/(dp(k-1)+dp(k))
        blayerb(k-1) = blayerb(k-1) + dp(k)*tempb
        blayerb(k) = blayerb(k) + dp(k-1)*tempb
        blevelb(k) = 0.0
      END DO
      blayerb(np+1) = 0.0
      CALL POPREAL4(dbs)
      CALL POPREAL4(rflxs)
      CALL POPREAL4(blevel(0))
      blayerb(1) = blayerb(1) + blayerb(0) + blevelb(0)
      blevelb(0) = 0.0
      blayerb(0) = 0.0
      DO k=np,1,-1
        CALL PLANCK_B(ibn, cb, ta_dev(i, k), ta_devb(i, k), blayer(k), &
&               blayerb(k))
        blayerb(k) = 0.0
      END DO
      CALL POPCONTROL4B(branch)
      IF (branch .LT. 4) THEN
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            CALL POPINTEGER4(n2oexp%end)
            CALL POPINTEGER4(n2oexp%start)
            CALL POPINTEGER4(co2exp%end)
            CALL POPINTEGER4(co2exp%start)
            CALL POPINTEGER4(conexp%end)
            CALL POPINTEGER4(conexp%start)
            CALL POPINTEGER4(h2oexp%end)
            CALL POPINTEGER4(h2oexp%start)
          ELSE
            CALL POPINTEGER4(h2oexp%end)
            CALL POPINTEGER4(h2oexp%start)
          END IF
        ELSE IF (branch .EQ. 2) THEN
          CALL POPINTEGER4(ch4exp%end)
          CALL POPINTEGER4(ch4exp%start)
          CALL POPINTEGER4(n2oexp%end)
          CALL POPINTEGER4(n2oexp%start)
          CALL POPINTEGER4(conexp%end)
          CALL POPINTEGER4(conexp%start)
          CALL POPINTEGER4(h2oexp%end)
          CALL POPINTEGER4(h2oexp%start)
        ELSE
          CALL POPINTEGER4(f22exp%end)
          CALL POPINTEGER4(f22exp%start)
          CALL POPINTEGER4(f12exp%end)
          CALL POPINTEGER4(f12exp%start)
          CALL POPINTEGER4(ch4exp%end)
          CALL POPINTEGER4(ch4exp%start)
          CALL POPINTEGER4(n2oexp%end)
          CALL POPINTEGER4(n2oexp%start)
          CALL POPINTEGER4(conexp%end)
          CALL POPINTEGER4(conexp%start)
          CALL POPINTEGER4(h2oexp%end)
          CALL POPINTEGER4(h2oexp%start)
        END IF
      ELSE IF (branch .LT. 6) THEN
        IF (branch .EQ. 4) THEN
          CALL POPINTEGER4(f11exp%end)
          CALL POPINTEGER4(f11exp%start)
          CALL POPINTEGER4(comexp%end)
          CALL POPINTEGER4(comexp%start)
          CALL POPINTEGER4(conexp%end)
          CALL POPINTEGER4(conexp%start)
          CALL POPINTEGER4(h2oexp%end)
          CALL POPINTEGER4(h2oexp%start)
        ELSE
          CALL POPINTEGER4(f22exp%end)
          CALL POPINTEGER4(f22exp%start)
          CALL POPINTEGER4(f12exp%end)
          CALL POPINTEGER4(f12exp%start)
          CALL POPINTEGER4(f11exp%end)
          CALL POPINTEGER4(f11exp%start)
          CALL POPINTEGER4(comexp%end)
          CALL POPINTEGER4(comexp%start)
          CALL POPINTEGER4(conexp%end)
          CALL POPINTEGER4(conexp%start)
          CALL POPINTEGER4(h2oexp%end)
          CALL POPINTEGER4(h2oexp%start)
        END IF
      ELSE IF (branch .EQ. 6) THEN
        CALL POPINTEGER4(conexp%end)
        CALL POPINTEGER4(conexp%start)
        CALL POPINTEGER4(h2oexp%end)
        CALL POPINTEGER4(h2oexp%start)
      ELSE IF (branch .EQ. 7) THEN
        CALL POPINTEGER4(conexp%end)
        CALL POPINTEGER4(conexp%start)
      END IF
      CALL POPREAL4ARRAY(exptbl, (np+1)*17)
    END IF
    CALL POPINTEGER4(ibn)
  END DO
  DO k=np+1,1,-1
    dfdts_devb(i, k) = 0.0
    flxd_devb(i, k) = 0.0
    flxu_devb(i, k) = 0.0
  END DO
  xx = pa(0)*0.001618*wa_dev(i, 1)*wa_dev(i, 1)*dp(0)
  temp1 = 1800./ta_dev(i, 1)
  xxb = EXP(temp1-6.081)*dcontb(0)
  ta_devb(i, 1) = ta_devb(i, 1) - EXP(temp1-6.081)*xx*temp1*dcontb(0)/&
&   ta_dev(i, 1)
  dcontb(0) = 0.0
  wa_devb(i, 1) = wa_devb(i, 1) + pa(0)*0.001618*dp(0)*2*wa_dev(i, 1)*&
&   xxb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) do3b(0) = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) dh2ob(0) = 0.0
  oa_devb(i, 1) = oa_devb(i, 1) + dp(0)*476.*do3b(0)
  do3b(0) = 0.0
  wa_devb(i, 1) = wa_devb(i, 1) + dp(0)*1.02*dh2ob(0)
  dh2ob(0) = 0.0
  ta_devb(i, 1) = ta_devb(i, 1) + dtb(0)
  dtb(0) = 0.0
  CALL POPREAL4(pa(0))
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    CALL POPREAL4(dp(0))
  ELSE
    CALL POPREAL4(dp(0))
  END IF
  DO k=np,1,-1
    DO l=4,1,-1
      cwc_devb(i, k, l) = cwc_devb(i, k, l) + cwc_colb(k, l)
      cwc_colb(k, l) = 0.0
      reff_devb(i, k, l) = reff_devb(i, k, l) + reff_colb(k, l)
      reff_colb(k, l) = 0.0
    END DO
    fcld_devb(i, k) = fcld_devb(i, k) + fcld_colb(k)
    fcld_colb(k) = 0.0
    xx = pa(k)*0.001618*wa_dev(i, k)*wa_dev(i, k)*dp(k)
    temp0 = ta_dev(i, k)
    temp = 1800./temp0
    xxb = EXP(temp-6.081)*dcontb(k)
    ta_devb(i, k) = ta_devb(i, k) - EXP(temp-6.081)*xx*temp*dcontb(k)/&
&     temp0
    dcontb(k) = 0.0
    wa_devb(i, k) = wa_devb(i, k) + pa(k)*0.001618*dp(k)*2*wa_dev(i, k)*&
&     xxb
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) do3b(k) = 0.0
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) dh2ob(k) = 0.0
    oa_devb(i, k) = oa_devb(i, k) + dp(k)*476.*do3b(k)
    do3b(k) = 0.0
    wa_devb(i, k) = wa_devb(i, k) + dp(k)*1.02*dh2ob(k)
    dh2ob(k) = 0.0
    ta_devb(i, k) = ta_devb(i, k) + dtb(k)
    dtb(k) = 0.0
  END DO
  asya_devb = asya_devb + asya_dev_tmpb
  ssaa_devb = ssaa_devb + ssaa_dev_tmpb
  taua_devb = taua_devb + taua_dev_tmpb
END SUBROUTINE IRRAD_B

!  Differentiation of planck in reverse (adjoint) mode:
!   gradient     of useful results: t xlayer
!   with respect to varying inputs: t
!----------------------------------!
!      SUBROUTINES GO HERE         !
!----------------------------------!
SUBROUTINE PLANCK_B(ibn, cb, t, tb, xlayer, xlayerb)
  IMPLICIT NONE
! spectral band index
  INTEGER :: ibn
! Planck table coefficients
  REAL :: cb(6, 10)
! temperature (K)
  REAL :: t
  REAL :: tb
! planck flux (w/m2)
  REAL :: xlayer
  REAL :: xlayerb
  REAL :: temp1
  REAL :: temp0
  REAL :: tempb
  REAL :: temp
  temp1 = cb(5, ibn) + cb(6, ibn)*t
  temp0 = cb(4, ibn) + t*temp1
  temp = cb(3, ibn) + t*temp0
  tempb = t**2*xlayerb
  tb = tb + (t**2*cb(6, ibn)+t*temp1+temp0)*tempb + (2*(t*temp)+cb(2, &
&   ibn))*xlayerb
END SUBROUTINE PLANCK_B

!  Differentiation of h2oexps in reverse (adjoint) mode:
!   gradient     of useful results: dh2o dt h2oexp
!   with respect to varying inputs: dh2o dt
SUBROUTINE H2OEXPS_B(ib, np, dh2o, dh2ob, pa, dt, dtb, xkw, aw, bw, pm, &
& mw, h2oexp, h2oexpb)
  IMPLICIT NONE
  INTEGER :: ib, np, ik, k
!---- input parameters ------
  REAL :: dh2o(0:np), pa(0:np), dt(0:np)
  REAL :: dh2ob(0:np), dtb(0:np)
!---- output parameters -----
  REAL :: h2oexp(0:np, 6)
  REAL :: h2oexpb(0:np, 6)
!---- static data -----
  INTEGER :: mw(9)
  REAL :: xkw(9), aw(9), bw(9), pm(9)
!---- temporary arrays -----
  REAL :: xh
  REAL :: xhb
  INTRINSIC EXP
  INTEGER :: branch
  REAL :: tempb
  REAL :: temp
!    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
!    and bw,  therefore, h2oexp for these sub-bands are identical.
  DO k=0,np
!-----xh is the scaled water vapor amount for line absorption
!     computed from Eq. (4.4).
    CALL PUSHREAL4(xh)
    xh = dh2o(k)*(pa(k)/500.)**pm(ib)*(1.+(aw(ib)+bw(ib)*dt(k))*dt(k))
!-----h2oexp is the water vapor transmittance of the layer k
!     due to line absorption
    h2oexp(k, 1) = EXP(-(xh*xkw(ib)))
!-----compute transmittances from Eq. (8.22)
    DO ik=2,6
      IF (mw(ib) .EQ. 6) THEN
        CALL PUSHREAL4(xh)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        CALL PUSHREAL4(h2oexp(k, ik))
        h2oexp(k, ik) = xh*xh*xh
        CALL PUSHCONTROL2B(3)
      ELSE IF (mw(ib) .EQ. 8) THEN
        CALL PUSHREAL4(xh)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        CALL PUSHREAL4(xh)
        xh = xh*xh
        CALL PUSHREAL4(h2oexp(k, ik))
        h2oexp(k, ik) = xh*xh
        CALL PUSHCONTROL2B(2)
      ELSE IF (mw(ib) .EQ. 9) THEN
        CALL PUSHREAL4(xh)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)*h2oexp(k, ik-1)
        CALL PUSHREAL4(h2oexp(k, ik))
        h2oexp(k, ik) = xh*xh*xh
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHREAL4(xh)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        CALL PUSHREAL4(xh)
        xh = xh*xh
        CALL PUSHREAL4(xh)
        xh = xh*xh
        CALL PUSHREAL4(h2oexp(k, ik))
        h2oexp(k, ik) = xh*xh
        CALL PUSHCONTROL2B(0)
      END IF
    END DO
  END DO
  DO k=np,0,-1
    DO ik=6,2,-1
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          CALL POPREAL4(h2oexp(k, ik))
          xhb = 2*xh*h2oexpb(k, ik)
          h2oexpb(k, ik) = 0.0
          CALL POPREAL4(xh)
          xhb = 2*xh*xhb
          CALL POPREAL4(xh)
          xhb = 2*xh*xhb
          CALL POPREAL4(xh)
          h2oexpb(k, ik-1) = h2oexpb(k, ik-1) + 2*h2oexp(k, ik-1)*xhb
        ELSE
          CALL POPREAL4(h2oexp(k, ik))
          xhb = 3*xh**2*h2oexpb(k, ik)
          h2oexpb(k, ik) = 0.0
          CALL POPREAL4(xh)
          h2oexpb(k, ik-1) = h2oexpb(k, ik-1) + 3*h2oexp(k, ik-1)**2*xhb
        END IF
      ELSE IF (branch .EQ. 2) THEN
        CALL POPREAL4(h2oexp(k, ik))
        xhb = 2*xh*h2oexpb(k, ik)
        h2oexpb(k, ik) = 0.0
        CALL POPREAL4(xh)
        xhb = 2*xh*xhb
        CALL POPREAL4(xh)
        h2oexpb(k, ik-1) = h2oexpb(k, ik-1) + 2*h2oexp(k, ik-1)*xhb
      ELSE
        CALL POPREAL4(h2oexp(k, ik))
        xhb = 3*xh**2*h2oexpb(k, ik)
        h2oexpb(k, ik) = 0.0
        CALL POPREAL4(xh)
        h2oexpb(k, ik-1) = h2oexpb(k, ik-1) + 2*h2oexp(k, ik-1)*xhb
      END IF
    END DO
    xhb = -(EXP(-(xkw(ib)*xh))*xkw(ib)*h2oexpb(k, 1))
    h2oexpb(k, 1) = 0.0
    CALL POPREAL4(xh)
    temp = aw(ib) + bw(ib)*dt(k)
    tempb = (pa(k)/500.)**pm(ib)*xhb
    dh2ob(k) = dh2ob(k) + (temp*dt(k)+1.)*tempb
    dtb(k) = dtb(k) + (dh2o(k)*temp+dt(k)*dh2o(k)*bw(ib))*tempb
  END DO
END SUBROUTINE H2OEXPS_B

!  Differentiation of conexps in reverse (adjoint) mode:
!   gradient     of useful results: dcont conexp
!   with respect to varying inputs: dcont conexp
SUBROUTINE CONEXPS_B(ib, np, dcont, dcontb, xke, conexp, conexpb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters ------
  REAL :: dcont(0:np)
  REAL :: dcontb(0:np)
!---- updated parameters -----
  REAL :: conexp(0:np, 3)
  REAL :: conexpb(0:np, 3)
!---- static data -----
  REAL :: xke(9)
  INTRINSIC EXP
  INTEGER :: branch
  DO k=0,np
    conexp(k, 1) = EXP(-(dcont(k)*xke(ib)))
!-----The absorption coefficients for sub-bands 3b and 3a are, respectively,
!     two and four times the absorption coefficient for sub-band 3c (Table 9).
!     Note that conexp(3) is for sub-band 3a. 
    IF (ib .EQ. 3) THEN
      CALL PUSHREAL4(conexp(k, 2))
      conexp(k, 2) = conexp(k, 1)*conexp(k, 1)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO k=np,0,-1
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      conexpb(k, 2) = conexpb(k, 2) + 2*conexp(k, 2)*conexpb(k, 3)
      conexpb(k, 3) = 0.0
      CALL POPREAL4(conexp(k, 2))
      conexpb(k, 1) = conexpb(k, 1) + 2*conexp(k, 1)*conexpb(k, 2)
      conexpb(k, 2) = 0.0
    END IF
    dcontb(k) = dcontb(k) - EXP(-(xke(ib)*dcont(k)))*xke(ib)*conexpb(k, &
&     1)
    conexpb(k, 1) = 0.0
  END DO
END SUBROUTINE CONEXPS_B

!  Differentiation of n2oexps in reverse (adjoint) mode:
!   gradient     of useful results: dt n2oexp
!   with respect to varying inputs: dt n2oexp
SUBROUTINE N2OEXPS_B(ib, np, dn2o, pa, dt, dtb, n2oexp, n2oexpb)
  IMPLICIT NONE
  INTEGER :: ib, k, np
!---- input parameters -----
  REAL :: dn2o(0:np), pa(0:np), dt(0:np)
  REAL :: dtb(0:np)
!---- output parameters -----
  REAL :: n2oexp(0:np, 4)
  REAL :: n2oexpb(0:np, 4)
!---- temporary arrays -----
  REAL :: xc, xc1, xc2
  REAL :: xcb, xc1b, xc2b
  INTRINSIC EXP
  INTEGER :: branch
  REAL :: tempb
!-----Scaling and absorption data are given in Table 5.
!     Transmittances are computed using Eqs. (8.21) and (8.22).
  DO k=0,np
!-----four exponential by powers of 21 for band 6.
    IF (ib .EQ. 6) THEN
      CALL PUSHREAL4(xc)
      xc = dn2o(k)*(1.+(1.9297e-3+4.3750e-6*dt(k))*dt(k))
      n2oexp(k, 1) = EXP(-(xc*6.31582e-2))
      xc = n2oexp(k, 1)*n2oexp(k, 1)*n2oexp(k, 1)
!-----four exponential by powers of 8 for band 7
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL4(xc)
      xc = dn2o(k)*(pa(k)/500.0)**0.48*(1.+(1.3804e-3+7.4838e-6*dt(k))*&
&       dt(k))
      n2oexp(k, 1) = EXP(-(xc*5.35779e-2))
      xc = n2oexp(k, 1)*n2oexp(k, 1)
      CALL PUSHREAL4(xc)
      xc = xc*xc
      CALL PUSHREAL4(n2oexp(k, 2))
      n2oexp(k, 2) = xc*xc
      CALL PUSHREAL4(xc)
      xc = n2oexp(k, 2)*n2oexp(k, 2)
      CALL PUSHREAL4(xc)
      xc = xc*xc
      CALL PUSHREAL4(n2oexp(k, 3))
      n2oexp(k, 3) = xc*xc
      CALL PUSHREAL4(xc)
      xc = n2oexp(k, 3)*n2oexp(k, 3)
      CALL PUSHREAL4(xc)
      xc = xc*xc
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO k=np,0,-1
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      xcb = 2*xc*n2oexpb(k, 4)
      n2oexpb(k, 4) = 0.0
      CALL POPREAL4(xc)
      xcb = 2*xc*xcb
      CALL POPREAL4(xc)
      n2oexpb(k, 3) = n2oexpb(k, 3) + 2*n2oexp(k, 3)*xcb
      CALL POPREAL4(n2oexp(k, 3))
      xcb = 2*xc*n2oexpb(k, 3)
      n2oexpb(k, 3) = 0.0
      CALL POPREAL4(xc)
      xcb = 2*xc*xcb
      CALL POPREAL4(xc)
      n2oexpb(k, 2) = n2oexpb(k, 2) + 2*n2oexp(k, 2)*xcb
      CALL POPREAL4(n2oexp(k, 2))
      xcb = 2*xc*n2oexpb(k, 2)
      n2oexpb(k, 2) = 0.0
      CALL POPREAL4(xc)
      xcb = 2*xc*xcb
      n2oexpb(k, 1) = n2oexpb(k, 1) + 2*n2oexp(k, 1)*xcb
      xc = dn2o(k)*(pa(k)/500.0)**0.48*(1.+(1.3804e-3+7.4838e-6*dt(k))*&
&       dt(k))
      xcb = -(EXP(-(5.35779e-2*xc))*5.35779e-2*n2oexpb(k, 1))
      n2oexpb(k, 1) = 0.0
      CALL POPREAL4(xc)
      tempb = (pa(k)/500.0)**0.48*dn2o(k)*xcb
      dtb(k) = dtb(k) + (7.4838e-6*dt(k)+dt(k)*7.4838e-6+1.3804e-3)*&
&       tempb
    ELSE
      xc1 = xc*xc
      xc2 = xc1*xc1
      xc2b = xc*xc1*n2oexpb(k, 2)
      xc1b = 2*xc1*xc2b + xc2*xc*n2oexpb(k, 2)
      xcb = 2*xc*xc1b + xc2*xc1*n2oexpb(k, 2)
      n2oexpb(k, 2) = 0.0
      n2oexpb(k, 1) = n2oexpb(k, 1) + 3*n2oexp(k, 1)**2*xcb
      xc = dn2o(k)*(1.+(1.9297e-3+4.3750e-6*dt(k))*dt(k))
      xcb = -(EXP(-(6.31582e-2*xc))*6.31582e-2*n2oexpb(k, 1))
      n2oexpb(k, 1) = 0.0
      CALL POPREAL4(xc)
      dtb(k) = dtb(k) + (dn2o(k)*(4.3750e-6*dt(k)+1.9297e-3)+dt(k)*dn2o(&
&       k)*4.3750e-6)*xcb
    END IF
  END DO
END SUBROUTINE N2OEXPS_B

!  Differentiation of ch4exps in reverse (adjoint) mode:
!   gradient     of useful results: dt ch4exp
!   with respect to varying inputs: dt ch4exp
SUBROUTINE CH4EXPS_B(ib, np, dch4, pa, dt, dtb, ch4exp, ch4expb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL :: dch4(0:np), pa(0:np), dt(0:np)
  REAL :: dtb(0:np)
!---- output parameters -----
  REAL :: ch4exp(0:np, 4)
  REAL :: ch4expb(0:np, 4)
!---- temporary arrays -----
  REAL :: xc
  REAL :: xcb
  INTRINSIC EXP
  INTEGER :: branch
  REAL :: tempb
!-----  Scaling and absorption data are given in Table 5 
  DO k=0,np
!-----four exponentials for band 6
    IF (ib .EQ. 6) THEN
      CALL PUSHREAL4(xc)
!-----four exponentials by powers of 12 for band 7
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL4(xc)
      xc = dch4(k)*(pa(k)/500.0)**0.65*(1.+(5.9590e-4-2.2931e-6*dt(k))*&
&       dt(k))
      ch4exp(k, 1) = EXP(-(xc*6.29247e-2))
      xc = ch4exp(k, 1)*ch4exp(k, 1)*ch4exp(k, 1)
      CALL PUSHREAL4(xc)
      xc = xc*xc
      CALL PUSHREAL4(ch4exp(k, 2))
      ch4exp(k, 2) = xc*xc
      CALL PUSHREAL4(xc)
      xc = ch4exp(k, 2)*ch4exp(k, 2)*ch4exp(k, 2)
      CALL PUSHREAL4(xc)
      xc = xc*xc
      CALL PUSHREAL4(ch4exp(k, 3))
      ch4exp(k, 3) = xc*xc
      CALL PUSHREAL4(xc)
      xc = ch4exp(k, 3)*ch4exp(k, 3)*ch4exp(k, 3)
      CALL PUSHREAL4(xc)
      xc = xc*xc
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO k=np,0,-1
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      xcb = 2*xc*ch4expb(k, 4)
      ch4expb(k, 4) = 0.0
      CALL POPREAL4(xc)
      xcb = 2*xc*xcb
      CALL POPREAL4(xc)
      ch4expb(k, 3) = ch4expb(k, 3) + 3*ch4exp(k, 3)**2*xcb
      CALL POPREAL4(ch4exp(k, 3))
      xcb = 2*xc*ch4expb(k, 3)
      ch4expb(k, 3) = 0.0
      CALL POPREAL4(xc)
      xcb = 2*xc*xcb
      CALL POPREAL4(xc)
      ch4expb(k, 2) = ch4expb(k, 2) + 3*ch4exp(k, 2)**2*xcb
      CALL POPREAL4(ch4exp(k, 2))
      xcb = 2*xc*ch4expb(k, 2)
      ch4expb(k, 2) = 0.0
      CALL POPREAL4(xc)
      xcb = 2*xc*xcb
      ch4expb(k, 1) = ch4expb(k, 1) + 3*ch4exp(k, 1)**2*xcb
      xc = dch4(k)*(pa(k)/500.0)**0.65*(1.+(5.9590e-4-2.2931e-6*dt(k))*&
&       dt(k))
      xcb = -(EXP(-(6.29247e-2*xc))*6.29247e-2*ch4expb(k, 1))
      ch4expb(k, 1) = 0.0
      CALL POPREAL4(xc)
      tempb = (pa(k)/500.0)**0.65*dch4(k)*xcb
      dtb(k) = dtb(k) + (5.9590e-4-dt(k)*2.2931e-6-2.2931e-6*dt(k))*&
&       tempb
    ELSE
      xc = dch4(k)*(1.+(1.7007e-2+1.5826e-4*dt(k))*dt(k))
      xcb = -(EXP(-(5.80708e-3*xc))*5.80708e-3*ch4expb(k, 1))
      ch4expb(k, 1) = 0.0
      CALL POPREAL4(xc)
      dtb(k) = dtb(k) + (dch4(k)*(1.5826e-4*dt(k)+1.7007e-2)+dt(k)*dch4(&
&       k)*1.5826e-4)*xcb
    END IF
  END DO
END SUBROUTINE CH4EXPS_B

!  Differentiation of comexps in reverse (adjoint) mode:
!   gradient     of useful results: dt comexp
!   with respect to varying inputs: dt comexp
SUBROUTINE COMEXPS_B(ib, np, dcom, dt, dtb, comexp, comexpb)
  IMPLICIT NONE
  INTEGER :: ib, ik, np, k
!---- input parameters -----
  REAL :: dcom(0:np), dt(0:np)
  REAL :: dtb(0:np)
!---- output parameters -----
  REAL :: comexp(0:np, 6)
  REAL :: comexpb(0:np, 6)
!---- temporary arrays -----
  REAL :: xc
  REAL :: xcb
  INTRINSIC EXP
  INTEGER :: branch
!-----  Scaling and absorpton data are given in Table 6
  DO k=0,np
    IF (ib .EQ. 4) THEN
      CALL PUSHREAL4(xc)
      xc = dcom(k)*(1.+(3.5775e-2+4.0447e-4*dt(k))*dt(k))
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (ib .EQ. 5) THEN
      CALL PUSHREAL4(xc)
      xc = dcom(k)*(1.+(3.4268e-2+3.7401e-4*dt(k))*dt(k))
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    comexp(k, 1) = EXP(-(xc*1.922e-7))
    DO ik=2,6
      CALL PUSHREAL4(xc)
      xc = comexp(k, ik-1)*comexp(k, ik-1)
      CALL PUSHREAL4(xc)
      xc = xc*xc
      CALL PUSHREAL4(comexp(k, ik))
      comexp(k, ik) = xc*comexp(k, ik-1)
    END DO
  END DO
  xcb = 0.0
  DO k=np,0,-1
    DO ik=6,2,-1
      CALL POPREAL4(comexp(k, ik))
      xcb = xcb + comexp(k, ik-1)*comexpb(k, ik)
      comexpb(k, ik-1) = comexpb(k, ik-1) + xc*comexpb(k, ik)
      comexpb(k, ik) = 0.0
      CALL POPREAL4(xc)
      xcb = 2*xc*xcb
      CALL POPREAL4(xc)
      comexpb(k, ik-1) = comexpb(k, ik-1) + 2*comexp(k, ik-1)*xcb
      xcb = 0.0
    END DO
    xcb = xcb - EXP(-(1.922e-7*xc))*1.922e-7*comexpb(k, 1)
    comexpb(k, 1) = 0.0
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(xc)
      dtb(k) = dtb(k) + (dcom(k)*(3.7401e-4*dt(k)+3.4268e-2)+dt(k)*dcom(&
&       k)*3.7401e-4)*xcb
      xcb = 0.0
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(xc)
      dtb(k) = dtb(k) + (dcom(k)*(4.0447e-4*dt(k)+3.5775e-2)+dt(k)*dcom(&
&       k)*4.0447e-4)*xcb
      xcb = 0.0
    END IF
  END DO
END SUBROUTINE COMEXPS_B

!  Differentiation of cfcexps in reverse (adjoint) mode:
!   gradient     of useful results: dt cfcexp
!   with respect to varying inputs: dt cfcexp
SUBROUTINE CFCEXPS_B(ib, np, a1, b1, fk1, a2, b2, fk2, dcfc, dt, dtb, &
& cfcexp, cfcexpb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL :: dcfc(0:np), dt(0:np)
  REAL :: dtb(0:np)
!---- output parameters -----
  REAL :: cfcexp(0:np)
  REAL :: cfcexpb(0:np)
!---- static data -----
  REAL :: a1, b1, fk1, a2, b2, fk2
!---- temporary arrays -----
  REAL :: xf
  REAL :: xfb
  INTRINSIC EXP
  INTEGER :: branch
  DO k=0,np
!-----compute the scaled cfc amount (xf) and exponential (cfcexp)
    IF (ib .EQ. 4) THEN
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO k=np,0,-1
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      xf = dcfc(k)*(1.+(a2+b2*dt(k))*dt(k))
      xfb = -(EXP(-(fk2*xf))*fk2*cfcexpb(k))
      cfcexpb(k) = 0.0
      dtb(k) = dtb(k) + (dcfc(k)*(a2+b2*dt(k))+dt(k)*dcfc(k)*b2)*xfb
    ELSE
      xf = dcfc(k)*(1.+(a1+b1*dt(k))*dt(k))
      xfb = -(EXP(-(fk1*xf))*fk1*cfcexpb(k))
      cfcexpb(k) = 0.0
      dtb(k) = dtb(k) + (dcfc(k)*(a1+b1*dt(k))+dt(k)*dcfc(k)*b1)*xfb
    END IF
  END DO
END SUBROUTINE CFCEXPS_B

!  Differentiation of b10exps in reverse (adjoint) mode:
!   gradient     of useful results: dh2o dcont dt co2exp h2oexp
!                n2oexp conexp
!   with respect to varying inputs: dh2o dcont dt co2exp h2oexp
!                n2oexp conexp
SUBROUTINE B10EXPS_B(np, dh2o, dh2ob, dcont, dcontb, dco2, dn2o, pa, dt&
& , dtb, h2oexp, h2oexpb, conexp, conexpb, co2exp, co2expb, n2oexp, &
& n2oexpb)
  IMPLICIT NONE
  INTEGER :: np, k
!---- input parameters -----
  REAL :: dh2o(0:np), dcont(0:np), dn2o(0:np)
  REAL :: dh2ob(0:np), dcontb(0:np)
  REAL :: dco2(0:np), pa(0:np), dt(0:np)
  REAL :: dtb(0:np)
!---- output parameters -----
  REAL :: h2oexp(0:np, 5), conexp(0:np), co2exp(0:np, 6), n2oexp(0:np, 2&
& )
  REAL :: h2oexpb(0:np, 5), conexpb(0:np), co2expb(0:np, 6), n2oexpb(0:&
& np, 2)
!---- temporary arrays -----
  REAL :: xx, xx1, xx2, xx3
  REAL :: xxb, xx1b, xx2b, xx3b
  INTRINSIC EXP
  REAL :: tempb2
  REAL :: tempb1
  REAL :: tempb0
  REAL :: tempb
  DO k=0,np
!-----Compute scaled h2o-line amount for Band 10 (Eq. 4.4 and Table 3).
    CALL PUSHREAL4(xx)
    xx = dh2o(k)*(pa(k)/500.0)*(1.+(0.0149+6.20e-5*dt(k))*dt(k))
!-----six exponentials by powers of 8
    h2oexp(k, 1) = EXP(-(xx*0.10624))
    xx = h2oexp(k, 1)*h2oexp(k, 1)
    CALL PUSHREAL4(xx)
    xx = xx*xx
    CALL PUSHREAL4(h2oexp(k, 2))
    h2oexp(k, 2) = xx*xx
    CALL PUSHREAL4(xx)
    xx = h2oexp(k, 2)*h2oexp(k, 2)
    CALL PUSHREAL4(xx)
    xx = xx*xx
    CALL PUSHREAL4(h2oexp(k, 3))
    h2oexp(k, 3) = xx*xx
    CALL PUSHREAL4(xx)
    xx = h2oexp(k, 3)*h2oexp(k, 3)
    CALL PUSHREAL4(xx)
    xx = xx*xx
    CALL PUSHREAL4(h2oexp(k, 4))
    h2oexp(k, 4) = xx*xx
    CALL PUSHREAL4(xx)
    xx = h2oexp(k, 4)*h2oexp(k, 4)
    CALL PUSHREAL4(xx)
    xx = xx*xx
!-----one exponential of h2o continuum for sub-band 3a (Table 9).
!-----Scaled co2 amount for the Band 10 (Eq. 4.4, Tables 3 and 6).
    CALL PUSHREAL4(xx)
    xx = dco2(k)*(pa(k)/300.0)**0.5*(1.+(0.0179+1.02e-4*dt(k))*dt(k))
!-----six exponentials by powers of 8
    co2exp(k, 1) = EXP(-(xx*2.656e-5))
    xx = co2exp(k, 1)*co2exp(k, 1)
    CALL PUSHREAL4(xx)
    xx = xx*xx
    CALL PUSHREAL4(co2exp(k, 2))
    co2exp(k, 2) = xx*xx
    CALL PUSHREAL4(xx)
    xx = co2exp(k, 2)*co2exp(k, 2)
    CALL PUSHREAL4(xx)
    xx = xx*xx
    CALL PUSHREAL4(co2exp(k, 3))
    co2exp(k, 3) = xx*xx
    CALL PUSHREAL4(xx)
    xx = co2exp(k, 3)*co2exp(k, 3)
    CALL PUSHREAL4(xx)
    xx = xx*xx
    CALL PUSHREAL4(co2exp(k, 4))
    co2exp(k, 4) = xx*xx
    CALL PUSHREAL4(xx)
    xx = co2exp(k, 4)*co2exp(k, 4)
    CALL PUSHREAL4(xx)
    xx = xx*xx
    CALL PUSHREAL4(co2exp(k, 5))
    co2exp(k, 5) = xx*xx
    CALL PUSHREAL4(xx)
    xx = co2exp(k, 5)*co2exp(k, 5)
    CALL PUSHREAL4(xx)
    xx = xx*xx
!-----Compute the scaled n2o amount for Band 10 (Table 5).
    CALL PUSHREAL4(xx)
    xx = dn2o(k)*(1.+(1.4476e-3+3.6656e-6*dt(k))*dt(k))
!-----Two exponentials by powers of 58
    n2oexp(k, 1) = EXP(-(xx*0.25238))
    xx = n2oexp(k, 1)*n2oexp(k, 1)
    CALL PUSHREAL4(xx1)
    xx1 = xx*xx
    CALL PUSHREAL4(xx1)
    xx1 = xx1*xx1
  END DO
  DO k=np,0,-1
    xx2 = xx1*xx1
    xx3 = xx2*xx2
    tempb = xx2*xx3*n2oexpb(k, 2)
    tempb0 = xx*xx1*n2oexpb(k, 2)
    xxb = xx1*tempb
    xx3b = xx2*tempb0
    xx2b = 2*xx2*xx3b + xx3*tempb0
    xx1b = 2*xx1*xx2b + xx*tempb
    n2oexpb(k, 2) = 0.0
    CALL POPREAL4(xx1)
    xx1b = 2*xx1*xx1b
    CALL POPREAL4(xx1)
    xxb = xxb + 2*xx*xx1b
    n2oexpb(k, 1) = n2oexpb(k, 1) + 2*n2oexp(k, 1)*xxb
    xx = dn2o(k)*(1.+(1.4476e-3+3.6656e-6*dt(k))*dt(k))
    xxb = -(EXP(-(0.25238*xx))*0.25238*n2oexpb(k, 1))
    n2oexpb(k, 1) = 0.0
    CALL POPREAL4(xx)
    dtb(k) = dtb(k) + (dn2o(k)*(3.6656e-6*dt(k)+1.4476e-3)+dt(k)*dn2o(k)&
&     *3.6656e-6)*xxb
    xxb = 2*xx*co2expb(k, 6)
    co2expb(k, 6) = 0.0
    CALL POPREAL4(xx)
    xxb = 2*xx*xxb
    CALL POPREAL4(xx)
    co2expb(k, 5) = co2expb(k, 5) + 2*co2exp(k, 5)*xxb
    CALL POPREAL4(co2exp(k, 5))
    xxb = 2*xx*co2expb(k, 5)
    co2expb(k, 5) = 0.0
    CALL POPREAL4(xx)
    xxb = 2*xx*xxb
    CALL POPREAL4(xx)
    co2expb(k, 4) = co2expb(k, 4) + 2*co2exp(k, 4)*xxb
    CALL POPREAL4(co2exp(k, 4))
    xxb = 2*xx*co2expb(k, 4)
    co2expb(k, 4) = 0.0
    CALL POPREAL4(xx)
    xxb = 2*xx*xxb
    CALL POPREAL4(xx)
    co2expb(k, 3) = co2expb(k, 3) + 2*co2exp(k, 3)*xxb
    CALL POPREAL4(co2exp(k, 3))
    xxb = 2*xx*co2expb(k, 3)
    co2expb(k, 3) = 0.0
    CALL POPREAL4(xx)
    xxb = 2*xx*xxb
    CALL POPREAL4(xx)
    co2expb(k, 2) = co2expb(k, 2) + 2*co2exp(k, 2)*xxb
    CALL POPREAL4(co2exp(k, 2))
    xxb = 2*xx*co2expb(k, 2)
    co2expb(k, 2) = 0.0
    CALL POPREAL4(xx)
    xxb = 2*xx*xxb
    co2expb(k, 1) = co2expb(k, 1) + 2*co2exp(k, 1)*xxb
    xx = dco2(k)*(pa(k)/300.0)**0.5*(1.+(0.0179+1.02e-4*dt(k))*dt(k))
    xxb = -(EXP(-(2.656e-5*xx))*2.656e-5*co2expb(k, 1))
    co2expb(k, 1) = 0.0
    CALL POPREAL4(xx)
    tempb1 = (pa(k)/300.0)**0.5*dco2(k)*xxb
    dcontb(k) = dcontb(k) - EXP(-(109.0*dcont(k)))*109.0*conexpb(k)
    conexpb(k) = 0.0
    xxb = 2*xx*h2oexpb(k, 5)
    h2oexpb(k, 5) = 0.0
    CALL POPREAL4(xx)
    xxb = 2*xx*xxb
    CALL POPREAL4(xx)
    h2oexpb(k, 4) = h2oexpb(k, 4) + 2*h2oexp(k, 4)*xxb
    CALL POPREAL4(h2oexp(k, 4))
    xxb = 2*xx*h2oexpb(k, 4)
    h2oexpb(k, 4) = 0.0
    CALL POPREAL4(xx)
    xxb = 2*xx*xxb
    CALL POPREAL4(xx)
    h2oexpb(k, 3) = h2oexpb(k, 3) + 2*h2oexp(k, 3)*xxb
    CALL POPREAL4(h2oexp(k, 3))
    xxb = 2*xx*h2oexpb(k, 3)
    h2oexpb(k, 3) = 0.0
    CALL POPREAL4(xx)
    xxb = 2*xx*xxb
    CALL POPREAL4(xx)
    h2oexpb(k, 2) = h2oexpb(k, 2) + 2*h2oexp(k, 2)*xxb
    CALL POPREAL4(h2oexp(k, 2))
    xxb = 2*xx*h2oexpb(k, 2)
    h2oexpb(k, 2) = 0.0
    CALL POPREAL4(xx)
    xxb = 2*xx*xxb
    h2oexpb(k, 1) = h2oexpb(k, 1) + 2*h2oexp(k, 1)*xxb
    xx = dh2o(k)*(pa(k)/500.0)*(1.+(0.0149+6.20e-5*dt(k))*dt(k))
    xxb = -(EXP(-(0.10624*xx))*0.10624*h2oexpb(k, 1))
    h2oexpb(k, 1) = 0.0
    CALL POPREAL4(xx)
    tempb2 = pa(k)*dh2o(k)*xxb/500.0
    dtb(k) = dtb(k) + (6.20e-5*dt(k)+dt(k)*6.20e-5+0.0149)*tempb2 + (&
&     1.02e-4*dt(k)+dt(k)*1.02e-4+0.0179)*tempb1
    dh2ob(k) = dh2ob(k) + ((6.20e-5*dt(k)+0.0149)*dt(k)+1.)*pa(k)*xxb/&
&     500.0
  END DO
END SUBROUTINE B10EXPS_B

!  Differentiation of tablup in reverse (adjoint) mode:
!   gradient     of useful results: s1 s2 s3 dt dw tran
!   with respect to varying inputs: s1 s2 s3 dt dw tran
SUBROUTINE TABLUP_B(nx1, nh1, dw, dwb, p, dt, dtb, s1, s1b, s2, s2b, s3&
& , s3b, w1, p1, dwe, dpe, coef1, coef2, coef3, tran, tranb)
  IMPLICIT NONE
  INTEGER :: nx1, nh1
!---- input parameters -----
  REAL :: w1, p1, dwe, dpe
  REAL :: dw, p, dt
  REAL :: dwb, dtb
  REAL :: coef1(nx1, nh1), coef2(nx1, nh1), coef3(nx1, nh1)
!---- update parameter -----
  REAL :: s1, s2, s3, tran
  REAL :: s1b, s2b, s3b, tranb
!---- temporary variables -----
  REAL :: we, pe, fw, fp, pa, pb, pc, ax, ba, bb, t1, ca, cb, t2
  REAL :: web, peb, fwb, fpb, pab, pbb, pcb, axb, bab, bbb, t1b, cab, &
& cbb, t2b
  REAL :: x1, x2, x3, xx, x1c
  REAL :: x1b, x2b, x3b, xxb, x1cb
  INTEGER :: iw, ip
  INTRINSIC LOG10
  INTRINSIC REAL
  INTRINSIC MIN
  INTRINSIC INT
  INTRINSIC MAX
  INTEGER :: branch
  REAL :: tempb0
  REAL :: tempb
  REAL :: y2
  REAL :: y1
!-----Compute effective pressure (x2) and temperature (x3) following 
!     Eqs. (8.28) and (8.29)
  s1 = s1 + dw
  s2 = s2 + p*dw
  s3 = s3 + dt*dw
  x1 = s1
  x1c = 1.0/s1
  x2 = s2*x1c
  x3 = s3*x1c
!-----normalize we and pe
!       we=(log10(x1)-w1)/dwe
!       pe=(log10(x2)-p1)/dpe
  we = (LOG10(x1)-w1)*dwe
  pe = (LOG10(x2)-p1)*dpe
  y1 = REAL(nh1 - 1)
  IF (we .GT. y1) THEN
    we = y1
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    we = we
  END IF
  y2 = REAL(nx1 - 1)
  IF (pe .GT. y2) THEN
    pe = y2
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    pe = pe
  END IF
!-----assign iw and ip and compute the distance of we and pe 
!     from iw and ip.
  iw = INT(we + 1.0)
  IF (iw .GT. nh1 - 1) THEN
    iw = nh1 - 1
  ELSE
    iw = iw
  END IF
  IF (iw .LT. 2) THEN
    iw = 2
  ELSE
    iw = iw
  END IF
  fw = we - REAL(iw - 1)
  ip = INT(pe + 1.0)
  IF (ip .GT. nx1 - 1) THEN
    ip = nx1 - 1
  ELSE
    ip = ip
  END IF
  IF (ip .LT. 1) THEN
    ip = 1
  ELSE
    ip = ip
  END IF
  fp = pe - REAL(ip - 1)
!-----linear interpolation in pressure
!       pa = coef1(ip,iw-1)*(1.-fp)+coef1(ip+1,iw-1)*fp
!       pb = coef1(ip,  iw)*(1.-fp)+coef1(ip+1,  iw)*fp
!       pc = coef1(ip,iw+1)*(1.-fp)+coef1(ip+1,iw+1)*fp
  pa = coef1(ip, iw-1) + (coef1(ip+1, iw-1)-coef1(ip, iw-1))*fp
  pb = coef1(ip, iw) + (coef1(ip+1, iw)-coef1(ip, iw))*fp
  pc = coef1(ip, iw+1) + (coef1(ip+1, iw+1)-coef1(ip, iw+1))*fp
!-----quadratic interpolation in absorber amount for coef1
!       ax = (-pa*(1.-fw)+pc*(1.+fw)) *fw*0.5 + pb*(1.-fw*fw)
  ax = ((pc+pa)*fw+(pc-pa))*fw*0.5 + pb*(1.-fw*fw)
!-----linear interpolation in absorber amount for coef2 and coef3
!       ba = coef2(ip,  iw)*(1.-fp)+coef2(ip+1,  iw)*fp
!       bb = coef2(ip,iw+1)*(1.-fp)+coef2(ip+1,iw+1)*fp
!       t1 = ba*(1.-fw) + bb*fw
  ba = coef2(ip, iw) + (coef2(ip+1, iw)-coef2(ip, iw))*fp
  bb = coef2(ip, iw+1) + (coef2(ip+1, iw+1)-coef2(ip, iw+1))*fp
  t1 = ba + (bb-ba)*fw
!       ca = coef3(ip,  iw)*(1.-fp)+coef3(ip+1,  iw)*fp
!       cb = coef3(ip,iw+1)*(1.-fp)+coef3(ip+1,iw+1)*fp
!       t2 = ca*(1.-fw) + cb*fw
  ca = coef3(ip, iw) + (coef3(ip+1, iw)-coef3(ip, iw))*fp
  cb = coef3(ip, iw+1) + (coef3(ip+1, iw+1)-coef3(ip, iw+1))*fp
  t2 = ca + (cb-ca)*fw
!-----update the total transmittance between levels k1 and k2
  xx = ax + (t1+t2*x3)*x3
  IF (xx .GT. 0.9999999) THEN
    xx = 0.9999999
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    xx = xx
  END IF
  IF (xx .LT. 0.0000001) THEN
    xx = 0.0000001
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    xx = xx
  END IF
  xxb = tran*tranb
  tranb = xx*tranb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) xxb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) xxb = 0.0
  tempb = x3*xxb
  axb = xxb
  t1b = tempb
  t2b = x3*tempb
  x3b = (t1+t2*x3)*xxb + t2*tempb
  cab = (1.0-fw)*t2b
  cbb = fw*t2b
  bab = (1.0-fw)*t1b
  bbb = fw*t1b
  tempb0 = 0.5*fw*axb
  fwb = (bb-ba)*t1b + (0.5*((pc+pa)*fw+pc-pa)-pb*2*fw)*axb + (pc+pa)*&
&   tempb0 + (cb-ca)*t2b
  pcb = (fw+1.0)*tempb0
  pab = (fw-1.0)*tempb0
  pbb = (1.-fw**2)*axb
  fpb = (coef3(ip+1, iw)-coef3(ip, iw))*cab + (coef2(ip+1, iw)-coef2(ip&
&   , iw))*bab + (coef1(ip+1, iw)-coef1(ip, iw))*pbb + (coef1(ip+1, iw-1&
&   )-coef1(ip, iw-1))*pab + (coef1(ip+1, iw+1)-coef1(ip, iw+1))*pcb + (&
&   coef2(ip+1, iw+1)-coef2(ip, iw+1))*bbb + (coef3(ip+1, iw+1)-coef3(ip&
&   , iw+1))*cbb
  peb = fpb
  web = fwb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) peb = 0.0
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) web = 0.0
  x2b = dpe*peb/(x2*LOG(10.0))
  x1b = dwe*web/(x1*LOG(10.0))
  s3b = s3b + x1c*x3b
  x1cb = s2*x2b + s3*x3b
  s2b = s2b + x1c*x2b
  s1b = s1b + x1b - x1cb/s1**2
  dtb = dtb + dw*s3b
  dwb = dwb + p*s2b + s1b + dt*s3b
END SUBROUTINE TABLUP_B

!  Differentiation of h2okdis in reverse (adjoint) mode:
!   gradient     of useful results: tran tcon h2oexp th2o conexp
!   with respect to varying inputs: tcon h2oexp th2o conexp
SUBROUTINE H2OKDIS_B(ib, np, k, fkw, gkw, ne, h2oexp, h2oexpb, conexp, &
& conexpb, th2o, th2ob, tcon, tconb, tran, tranb)
  IMPLICIT NONE
!---- input parameters ------
  INTEGER :: ib, ne, np, k
  REAL :: h2oexp(0:np, 6), conexp(0:np, 3)
  REAL :: h2oexpb(0:np, 6), conexpb(0:np, 3)
  REAL :: fkw(6, 9), gkw(6, 3)
!---- updated parameters -----
  REAL :: th2o(6), tcon(3), tran
  REAL :: th2ob(6), tconb(3), tranb
!---- temporary arrays -----
  REAL :: trnth2o
  REAL :: trnth2ob
  INTEGER :: i
  REAL :: tempb2
  REAL :: tempb1
  REAL :: tempb0
  REAL :: tempb
!-----tco2 are the six exp factors between levels k1 and k2 
!     tran is the updated total transmittance between levels k1 and k2
!-----th2o is the 6 exp factors between levels k1 and k2 due to
!     h2o line absorption. 
!-----tcon is the 3 exp factors between levels k1 and k2 due to
!     h2o continuum absorption.
!-----trnth2o is the total transmittance between levels k1 and k2 due
!     to both line and continuum absorption.
!-----Compute th2o following Eq. (8.23).
  CALL PUSHREAL4(th2o(1))
  th2o(1) = th2o(1)*h2oexp(k, 1)
  CALL PUSHREAL4(th2o(2))
  th2o(2) = th2o(2)*h2oexp(k, 2)
  CALL PUSHREAL4(th2o(3))
  th2o(3) = th2o(3)*h2oexp(k, 3)
  CALL PUSHREAL4(th2o(4))
  th2o(4) = th2o(4)*h2oexp(k, 4)
  CALL PUSHREAL4(th2o(5))
  th2o(5) = th2o(5)*h2oexp(k, 5)
  CALL PUSHREAL4(th2o(6))
  th2o(6) = th2o(6)*h2oexp(k, 6)
  IF (ne .EQ. 0) THEN
    trnth2ob = tran*tranb
    th2ob(1) = th2ob(1) + fkw(1, ib)*trnth2ob
    th2ob(2) = th2ob(2) + fkw(2, ib)*trnth2ob
    th2ob(3) = th2ob(3) + fkw(3, ib)*trnth2ob
    th2ob(4) = th2ob(4) + fkw(4, ib)*trnth2ob
    th2ob(5) = th2ob(5) + fkw(5, ib)*trnth2ob
    th2ob(6) = th2ob(6) + fkw(6, ib)*trnth2ob
  ELSE IF (ne .EQ. 1) THEN
!-----Compute trnh2o following Eqs. (8.25) and (4.27).
    CALL PUSHREAL4(tcon(1))
    tcon(1) = tcon(1)*conexp(k, 1)
    trnth2ob = tran*tranb
    tempb = tcon(1)*trnth2ob
    th2ob(1) = th2ob(1) + fkw(1, ib)*tempb
    th2ob(2) = th2ob(2) + fkw(2, ib)*tempb
    th2ob(3) = th2ob(3) + fkw(3, ib)*tempb
    th2ob(4) = th2ob(4) + fkw(4, ib)*tempb
    th2ob(5) = th2ob(5) + fkw(5, ib)*tempb
    th2ob(6) = th2ob(6) + fkw(6, ib)*tempb
    tconb(1) = tconb(1) + (fkw(1, ib)*th2o(1)+fkw(2, ib)*th2o(2)+fkw(3, &
&     ib)*th2o(3)+fkw(4, ib)*th2o(4)+fkw(5, ib)*th2o(5)+fkw(6, ib)*th2o(&
&     6))*trnth2ob
    CALL POPREAL4(tcon(1))
    conexpb(k, 1) = conexpb(k, 1) + tcon(1)*tconb(1)
    tconb(1) = conexp(k, 1)*tconb(1)
  ELSE
!-----For band 3. This band is divided into 3 subbands.
    CALL PUSHREAL4(tcon(1))
    tcon(1) = tcon(1)*conexp(k, 1)
    CALL PUSHREAL4(tcon(2))
    tcon(2) = tcon(2)*conexp(k, 2)
    CALL PUSHREAL4(tcon(3))
    tcon(3) = tcon(3)*conexp(k, 3)
!-----Compute trnh2o following Eqs. (4.29) and (8.25).
    trnth2ob = tran*tranb
    tempb0 = tcon(1)*trnth2ob
    tempb1 = tcon(2)*trnth2ob
    tempb2 = tcon(3)*trnth2ob
    th2ob(1) = th2ob(1) + gkw(1, 3)*tempb2 + gkw(1, 2)*tempb1 + gkw(1, 1&
&     )*tempb0
    th2ob(2) = th2ob(2) + gkw(2, 3)*tempb2 + gkw(2, 2)*tempb1 + gkw(2, 1&
&     )*tempb0
    th2ob(3) = th2ob(3) + gkw(3, 3)*tempb2 + gkw(3, 2)*tempb1 + gkw(3, 1&
&     )*tempb0
    th2ob(4) = th2ob(4) + gkw(4, 3)*tempb2 + gkw(4, 2)*tempb1 + gkw(4, 1&
&     )*tempb0
    th2ob(5) = th2ob(5) + gkw(5, 3)*tempb2 + gkw(5, 2)*tempb1 + gkw(5, 1&
&     )*tempb0
    th2ob(6) = th2ob(6) + gkw(6, 3)*tempb2 + gkw(6, 2)*tempb1 + gkw(6, 1&
&     )*tempb0
    tconb(1) = tconb(1) + (gkw(1, 1)*th2o(1)+gkw(2, 1)*th2o(2)+gkw(3, 1)&
&     *th2o(3)+gkw(4, 1)*th2o(4)+gkw(5, 1)*th2o(5)+gkw(6, 1)*th2o(6))*&
&     trnth2ob
    tconb(2) = tconb(2) + (gkw(1, 2)*th2o(1)+gkw(2, 2)*th2o(2)+gkw(3, 2)&
&     *th2o(3)+gkw(4, 2)*th2o(4)+gkw(5, 2)*th2o(5)+gkw(6, 2)*th2o(6))*&
&     trnth2ob
    tconb(3) = tconb(3) + (gkw(1, 3)*th2o(1)+gkw(2, 3)*th2o(2)+gkw(3, 3)&
&     *th2o(3)+gkw(4, 3)*th2o(4)+gkw(5, 3)*th2o(5)+gkw(6, 3)*th2o(6))*&
&     trnth2ob
    CALL POPREAL4(tcon(3))
    conexpb(k, 3) = conexpb(k, 3) + tcon(3)*tconb(3)
    tconb(3) = conexp(k, 3)*tconb(3)
    CALL POPREAL4(tcon(2))
    conexpb(k, 2) = conexpb(k, 2) + tcon(2)*tconb(2)
    tconb(2) = conexp(k, 2)*tconb(2)
    CALL POPREAL4(tcon(1))
    conexpb(k, 1) = conexpb(k, 1) + tcon(1)*tconb(1)
    tconb(1) = conexp(k, 1)*tconb(1)
  END IF
  CALL POPREAL4(th2o(6))
  h2oexpb(k, 6) = h2oexpb(k, 6) + th2o(6)*th2ob(6)
  th2ob(6) = h2oexp(k, 6)*th2ob(6)
  CALL POPREAL4(th2o(5))
  h2oexpb(k, 5) = h2oexpb(k, 5) + th2o(5)*th2ob(5)
  th2ob(5) = h2oexp(k, 5)*th2ob(5)
  CALL POPREAL4(th2o(4))
  h2oexpb(k, 4) = h2oexpb(k, 4) + th2o(4)*th2ob(4)
  th2ob(4) = h2oexp(k, 4)*th2ob(4)
  CALL POPREAL4(th2o(3))
  h2oexpb(k, 3) = h2oexpb(k, 3) + th2o(3)*th2ob(3)
  th2ob(3) = h2oexp(k, 3)*th2ob(3)
  CALL POPREAL4(th2o(2))
  h2oexpb(k, 2) = h2oexpb(k, 2) + th2o(2)*th2ob(2)
  th2ob(2) = h2oexp(k, 2)*th2ob(2)
  CALL POPREAL4(th2o(1))
  h2oexpb(k, 1) = h2oexpb(k, 1) + th2o(1)*th2ob(1)
  th2ob(1) = h2oexp(k, 1)*th2ob(1)
END SUBROUTINE H2OKDIS_B

!  Differentiation of n2okdis in reverse (adjoint) mode:
!   gradient     of useful results: tran tn2o n2oexp
!   with respect to varying inputs: tran tn2o n2oexp
SUBROUTINE N2OKDIS_B(ib, np, k, n2oexp, n2oexpb, tn2o, tn2ob, tran, &
& tranb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL :: n2oexp(0:np, 4)
  REAL :: n2oexpb(0:np, 4)
!---- updated parameters -----
  REAL :: tn2o(4), tran
  REAL :: tn2ob(4), tranb
!---- temporary arrays -----
  REAL :: xc
  REAL :: xcb
  INTEGER :: branch
!-----tn2o is computed from Eq. (8.23). 
!     xc is the total n2o transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.
!-----band 6
  IF (ib .EQ. 6) THEN
    CALL PUSHREAL4(tn2o(1))
    tn2o(1) = tn2o(1)*n2oexp(k, 1)
    xc = 0.940414*tn2o(1)
    CALL PUSHREAL4(tn2o(2))
    tn2o(2) = tn2o(2)*n2oexp(k, 2)
    xc = xc + 0.059586*tn2o(2)
!-----band 7
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHREAL4(tn2o(1))
    tn2o(1) = tn2o(1)*n2oexp(k, 1)
    xc = 0.561961*tn2o(1)
    CALL PUSHREAL4(tn2o(2))
    tn2o(2) = tn2o(2)*n2oexp(k, 2)
    xc = xc + 0.138707*tn2o(2)
    CALL PUSHREAL4(tn2o(3))
    tn2o(3) = tn2o(3)*n2oexp(k, 3)
    xc = xc + 0.240670*tn2o(3)
    CALL PUSHREAL4(tn2o(4))
    tn2o(4) = tn2o(4)*n2oexp(k, 4)
    xc = xc + 0.058662*tn2o(4)
    CALL PUSHCONTROL1B(1)
  END IF
  xcb = tran*tranb
  tranb = xc*tranb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    tn2ob(2) = tn2ob(2) + 0.059586*xcb
    CALL POPREAL4(tn2o(2))
    n2oexpb(k, 2) = n2oexpb(k, 2) + tn2o(2)*tn2ob(2)
    tn2ob(2) = n2oexp(k, 2)*tn2ob(2)
    tn2ob(1) = tn2ob(1) + 0.940414*xcb
    CALL POPREAL4(tn2o(1))
    n2oexpb(k, 1) = n2oexpb(k, 1) + tn2o(1)*tn2ob(1)
    tn2ob(1) = n2oexp(k, 1)*tn2ob(1)
  ELSE
    tn2ob(4) = tn2ob(4) + 0.058662*xcb
    CALL POPREAL4(tn2o(4))
    n2oexpb(k, 4) = n2oexpb(k, 4) + tn2o(4)*tn2ob(4)
    tn2ob(4) = n2oexp(k, 4)*tn2ob(4)
    tn2ob(3) = tn2ob(3) + 0.240670*xcb
    CALL POPREAL4(tn2o(3))
    n2oexpb(k, 3) = n2oexpb(k, 3) + tn2o(3)*tn2ob(3)
    tn2ob(3) = n2oexp(k, 3)*tn2ob(3)
    tn2ob(2) = tn2ob(2) + 0.138707*xcb
    CALL POPREAL4(tn2o(2))
    n2oexpb(k, 2) = n2oexpb(k, 2) + tn2o(2)*tn2ob(2)
    tn2ob(2) = n2oexp(k, 2)*tn2ob(2)
    tn2ob(1) = tn2ob(1) + 0.561961*xcb
    CALL POPREAL4(tn2o(1))
    n2oexpb(k, 1) = n2oexpb(k, 1) + tn2o(1)*tn2ob(1)
    tn2ob(1) = n2oexp(k, 1)*tn2ob(1)
  END IF
END SUBROUTINE N2OKDIS_B

!  Differentiation of ch4kdis in reverse (adjoint) mode:
!   gradient     of useful results: tran ch4exp tch4
!   with respect to varying inputs: tran ch4exp tch4
SUBROUTINE CH4KDIS_B(ib, np, k, ch4exp, ch4expb, tch4, tch4b, tran, &
& tranb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL :: ch4exp(0:np, 4)
  REAL :: ch4expb(0:np, 4)
!---- updated parameters -----
  REAL :: tch4(4), tran
  REAL :: tch4b(4), tranb
!---- temporary arrays -----
  REAL :: xc
  REAL :: xcb
  INTEGER :: branch
!-----tch4 is computed from Eq. (8.23). 
!     xc is the total ch4 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.
!-----band 6
  IF (ib .EQ. 6) THEN
    CALL PUSHREAL4(tch4(1))
    tch4(1) = tch4(1)*ch4exp(k, 1)
    xc = tch4(1)
!-----band 7
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHREAL4(tch4(1))
    tch4(1) = tch4(1)*ch4exp(k, 1)
    xc = 0.610650*tch4(1)
    CALL PUSHREAL4(tch4(2))
    tch4(2) = tch4(2)*ch4exp(k, 2)
    xc = xc + 0.280212*tch4(2)
    CALL PUSHREAL4(tch4(3))
    tch4(3) = tch4(3)*ch4exp(k, 3)
    xc = xc + 0.107349*tch4(3)
    CALL PUSHREAL4(tch4(4))
    tch4(4) = tch4(4)*ch4exp(k, 4)
    xc = xc + 0.001789*tch4(4)
    CALL PUSHCONTROL1B(1)
  END IF
  xcb = tran*tranb
  tranb = xc*tranb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    tch4b(1) = tch4b(1) + xcb
    CALL POPREAL4(tch4(1))
    ch4expb(k, 1) = ch4expb(k, 1) + tch4(1)*tch4b(1)
    tch4b(1) = ch4exp(k, 1)*tch4b(1)
  ELSE
    tch4b(4) = tch4b(4) + 0.001789*xcb
    CALL POPREAL4(tch4(4))
    ch4expb(k, 4) = ch4expb(k, 4) + tch4(4)*tch4b(4)
    tch4b(4) = ch4exp(k, 4)*tch4b(4)
    tch4b(3) = tch4b(3) + 0.107349*xcb
    CALL POPREAL4(tch4(3))
    ch4expb(k, 3) = ch4expb(k, 3) + tch4(3)*tch4b(3)
    tch4b(3) = ch4exp(k, 3)*tch4b(3)
    tch4b(2) = tch4b(2) + 0.280212*xcb
    CALL POPREAL4(tch4(2))
    ch4expb(k, 2) = ch4expb(k, 2) + tch4(2)*tch4b(2)
    tch4b(2) = ch4exp(k, 2)*tch4b(2)
    tch4b(1) = tch4b(1) + 0.610650*xcb
    CALL POPREAL4(tch4(1))
    ch4expb(k, 1) = ch4expb(k, 1) + tch4(1)*tch4b(1)
    tch4b(1) = ch4exp(k, 1)*tch4b(1)
  END IF
END SUBROUTINE CH4KDIS_B

!  Differentiation of comkdis in reverse (adjoint) mode:
!   gradient     of useful results: tran tcom comexp
!   with respect to varying inputs: tran tcom comexp
SUBROUTINE COMKDIS_B(ib, np, k, comexp, comexpb, tcom, tcomb, tran, &
& tranb)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL :: comexp(0:np, 6)
  REAL :: comexpb(0:np, 6)
!---- updated parameters -----
  REAL :: tcom(6), tran
  REAL :: tcomb(6), tranb
!---- temporary arrays -----
  REAL :: xc
  REAL :: xcb
  INTEGER :: branch
!-----tcom is computed from Eq. (8.23). 
!     xc is the total co2 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 6.
!-----band 4
  IF (ib .EQ. 4) THEN
    CALL PUSHREAL4(tcom(1))
    tcom(1) = tcom(1)*comexp(k, 1)
    xc = 0.12159*tcom(1)
    CALL PUSHREAL4(tcom(2))
    tcom(2) = tcom(2)*comexp(k, 2)
    xc = xc + 0.24359*tcom(2)
    CALL PUSHREAL4(tcom(3))
    tcom(3) = tcom(3)*comexp(k, 3)
    xc = xc + 0.24981*tcom(3)
    CALL PUSHREAL4(tcom(4))
    tcom(4) = tcom(4)*comexp(k, 4)
    xc = xc + 0.26427*tcom(4)
    CALL PUSHREAL4(tcom(5))
    tcom(5) = tcom(5)*comexp(k, 5)
    xc = xc + 0.07807*tcom(5)
    CALL PUSHREAL4(tcom(6))
    tcom(6) = tcom(6)*comexp(k, 6)
    xc = xc + 0.04267*tcom(6)
!-----band 5
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHREAL4(tcom(1))
    tcom(1) = tcom(1)*comexp(k, 1)
    xc = 0.06869*tcom(1)
    CALL PUSHREAL4(tcom(2))
    tcom(2) = tcom(2)*comexp(k, 2)
    xc = xc + 0.14795*tcom(2)
    CALL PUSHREAL4(tcom(3))
    tcom(3) = tcom(3)*comexp(k, 3)
    xc = xc + 0.19512*tcom(3)
    CALL PUSHREAL4(tcom(4))
    tcom(4) = tcom(4)*comexp(k, 4)
    xc = xc + 0.33446*tcom(4)
    CALL PUSHREAL4(tcom(5))
    tcom(5) = tcom(5)*comexp(k, 5)
    xc = xc + 0.17199*tcom(5)
    CALL PUSHREAL4(tcom(6))
    tcom(6) = tcom(6)*comexp(k, 6)
    xc = xc + 0.08179*tcom(6)
    CALL PUSHCONTROL1B(1)
  END IF
  xcb = tran*tranb
  tranb = xc*tranb
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    tcomb(6) = tcomb(6) + 0.04267*xcb
    CALL POPREAL4(tcom(6))
    comexpb(k, 6) = comexpb(k, 6) + tcom(6)*tcomb(6)
    tcomb(6) = comexp(k, 6)*tcomb(6)
    tcomb(5) = tcomb(5) + 0.07807*xcb
    CALL POPREAL4(tcom(5))
    comexpb(k, 5) = comexpb(k, 5) + tcom(5)*tcomb(5)
    tcomb(5) = comexp(k, 5)*tcomb(5)
    tcomb(4) = tcomb(4) + 0.26427*xcb
    CALL POPREAL4(tcom(4))
    comexpb(k, 4) = comexpb(k, 4) + tcom(4)*tcomb(4)
    tcomb(4) = comexp(k, 4)*tcomb(4)
    tcomb(3) = tcomb(3) + 0.24981*xcb
    CALL POPREAL4(tcom(3))
    comexpb(k, 3) = comexpb(k, 3) + tcom(3)*tcomb(3)
    tcomb(3) = comexp(k, 3)*tcomb(3)
    tcomb(2) = tcomb(2) + 0.24359*xcb
    CALL POPREAL4(tcom(2))
    comexpb(k, 2) = comexpb(k, 2) + tcom(2)*tcomb(2)
    tcomb(2) = comexp(k, 2)*tcomb(2)
    tcomb(1) = tcomb(1) + 0.12159*xcb
    CALL POPREAL4(tcom(1))
    comexpb(k, 1) = comexpb(k, 1) + tcom(1)*tcomb(1)
    tcomb(1) = comexp(k, 1)*tcomb(1)
  ELSE
    tcomb(6) = tcomb(6) + 0.08179*xcb
    CALL POPREAL4(tcom(6))
    comexpb(k, 6) = comexpb(k, 6) + tcom(6)*tcomb(6)
    tcomb(6) = comexp(k, 6)*tcomb(6)
    tcomb(5) = tcomb(5) + 0.17199*xcb
    CALL POPREAL4(tcom(5))
    comexpb(k, 5) = comexpb(k, 5) + tcom(5)*tcomb(5)
    tcomb(5) = comexp(k, 5)*tcomb(5)
    tcomb(4) = tcomb(4) + 0.33446*xcb
    CALL POPREAL4(tcom(4))
    comexpb(k, 4) = comexpb(k, 4) + tcom(4)*tcomb(4)
    tcomb(4) = comexp(k, 4)*tcomb(4)
    tcomb(3) = tcomb(3) + 0.19512*xcb
    CALL POPREAL4(tcom(3))
    comexpb(k, 3) = comexpb(k, 3) + tcom(3)*tcomb(3)
    tcomb(3) = comexp(k, 3)*tcomb(3)
    tcomb(2) = tcomb(2) + 0.14795*xcb
    CALL POPREAL4(tcom(2))
    comexpb(k, 2) = comexpb(k, 2) + tcom(2)*tcomb(2)
    tcomb(2) = comexp(k, 2)*tcomb(2)
    tcomb(1) = tcomb(1) + 0.06869*xcb
    CALL POPREAL4(tcom(1))
    comexpb(k, 1) = comexpb(k, 1) + tcom(1)*tcomb(1)
    tcomb(1) = comexp(k, 1)*tcomb(1)
  END IF
END SUBROUTINE COMKDIS_B

!  Differentiation of cfckdis in reverse (adjoint) mode:
!   gradient     of useful results: tran tcfc cfcexp
!   with respect to varying inputs: tran tcfc cfcexp
SUBROUTINE CFCKDIS_B(np, k, cfcexp, cfcexpb, tcfc, tcfcb, tran, tranb)
  IMPLICIT NONE
!---- input parameters -----
  INTEGER :: k, np
  REAL :: cfcexp(0:np)
  REAL :: cfcexpb(0:np)
!---- updated parameters -----
  REAL :: tcfc, tran
  REAL :: tcfcb, tranb
!-----tcfc is the exp factors between levels k1 and k2. 
  CALL PUSHREAL4(tcfc)
  tcfc = tcfc*cfcexp(k)
  tcfcb = tcfcb + tran*tranb
  tranb = tcfc*tranb
  CALL POPREAL4(tcfc)
  cfcexpb(k) = cfcexpb(k) + tcfc*tcfcb
  tcfcb = cfcexp(k)*tcfcb
END SUBROUTINE CFCKDIS_B

!  Differentiation of b10kdis in reverse (adjoint) mode:
!   gradient     of useful results: tran tn2o co2exp tcon h2oexp
!                n2oexp th2o tco2 conexp
!   with respect to varying inputs: tn2o co2exp tcon h2oexp n2oexp
!                th2o tco2 conexp
SUBROUTINE B10KDIS_B(np, k, h2oexp, h2oexpb, conexp, conexpb, co2exp, &
& co2expb, n2oexp, n2oexpb, th2o, th2ob, tcon, tconb, tco2, tco2b, tn2o&
& , tn2ob, tran, tranb)
  IMPLICIT NONE
  INTEGER :: np, k
!---- input parameters -----
  REAL :: h2oexp(0:np, 5), conexp(0:np), co2exp(0:np, 6), n2oexp(0:np, 2&
& )
  REAL :: h2oexpb(0:np, 5), conexpb(0:np), co2expb(0:np, 6), n2oexpb(0:&
& np, 2)
!---- updated parameters -----
  REAL :: th2o(6), tcon(3), tco2(6), tn2o(4), tran
  REAL :: th2ob(6), tconb(3), tco2b(6), tn2ob(4), tranb
!---- temporary arrays -----
  REAL :: xx
  REAL :: xxb
!-----For h2o line. The k-distribution functions are given in Table 4.
  CALL PUSHREAL4(th2o(1))
  th2o(1) = th2o(1)*h2oexp(k, 1)
  xx = 0.3153*th2o(1)
  CALL PUSHREAL4(th2o(2))
  th2o(2) = th2o(2)*h2oexp(k, 2)
  xx = xx + 0.4604*th2o(2)
  CALL PUSHREAL4(th2o(3))
  th2o(3) = th2o(3)*h2oexp(k, 3)
  xx = xx + 0.1326*th2o(3)
  CALL PUSHREAL4(th2o(4))
  th2o(4) = th2o(4)*h2oexp(k, 4)
  xx = xx + 0.0798*th2o(4)
  CALL PUSHREAL4(th2o(5))
  th2o(5) = th2o(5)*h2oexp(k, 5)
  xx = xx + 0.0119*th2o(5)
  tran = xx
!-----For h2o continuum. Note that conexp(k,3) is for subband 3a.
  CALL PUSHREAL4(tcon(1))
  tcon(1) = tcon(1)*conexp(k)
  CALL PUSHREAL4(tran)
  tran = tran*tcon(1)
!-----For co2 (Table 6)
  CALL PUSHREAL4(tco2(1))
  tco2(1) = tco2(1)*co2exp(k, 1)
  xx = 0.2673*tco2(1)
  CALL PUSHREAL4(tco2(2))
  tco2(2) = tco2(2)*co2exp(k, 2)
  xx = xx + 0.2201*tco2(2)
  CALL PUSHREAL4(tco2(3))
  tco2(3) = tco2(3)*co2exp(k, 3)
  xx = xx + 0.2106*tco2(3)
  CALL PUSHREAL4(tco2(4))
  tco2(4) = tco2(4)*co2exp(k, 4)
  xx = xx + 0.2409*tco2(4)
  CALL PUSHREAL4(tco2(5))
  tco2(5) = tco2(5)*co2exp(k, 5)
  xx = xx + 0.0196*tco2(5)
  CALL PUSHREAL4(tco2(6))
  tco2(6) = tco2(6)*co2exp(k, 6)
  xx = xx + 0.0415*tco2(6)
  CALL PUSHREAL4(tran)
  tran = tran*xx
!-----For n2o (Table 5)
  CALL PUSHREAL4(tn2o(1))
  tn2o(1) = tn2o(1)*n2oexp(k, 1)
  CALL PUSHREAL4(xx)
  xx = 0.970831*tn2o(1)
  CALL PUSHREAL4(tn2o(2))
  tn2o(2) = tn2o(2)*n2oexp(k, 2)
  xx = xx + 0.029169*tn2o(2)
  xxb = tran*tranb
  tranb = (xx-1.0)*tranb
  tn2ob(2) = tn2ob(2) + 0.029169*xxb
  CALL POPREAL4(tn2o(2))
  n2oexpb(k, 2) = n2oexpb(k, 2) + tn2o(2)*tn2ob(2)
  tn2ob(2) = n2oexp(k, 2)*tn2ob(2)
  CALL POPREAL4(xx)
  tn2ob(1) = tn2ob(1) + 0.970831*xxb
  CALL POPREAL4(tn2o(1))
  n2oexpb(k, 1) = n2oexpb(k, 1) + tn2o(1)*tn2ob(1)
  tn2ob(1) = n2oexp(k, 1)*tn2ob(1)
  CALL POPREAL4(tran)
  xxb = tran*tranb
  tranb = xx*tranb
  tco2b(6) = tco2b(6) + 0.0415*xxb
  CALL POPREAL4(tco2(6))
  co2expb(k, 6) = co2expb(k, 6) + tco2(6)*tco2b(6)
  tco2b(6) = co2exp(k, 6)*tco2b(6)
  tco2b(5) = tco2b(5) + 0.0196*xxb
  CALL POPREAL4(tco2(5))
  co2expb(k, 5) = co2expb(k, 5) + tco2(5)*tco2b(5)
  tco2b(5) = co2exp(k, 5)*tco2b(5)
  tco2b(4) = tco2b(4) + 0.2409*xxb
  CALL POPREAL4(tco2(4))
  co2expb(k, 4) = co2expb(k, 4) + tco2(4)*tco2b(4)
  tco2b(4) = co2exp(k, 4)*tco2b(4)
  tco2b(3) = tco2b(3) + 0.2106*xxb
  CALL POPREAL4(tco2(3))
  co2expb(k, 3) = co2expb(k, 3) + tco2(3)*tco2b(3)
  tco2b(3) = co2exp(k, 3)*tco2b(3)
  tco2b(2) = tco2b(2) + 0.2201*xxb
  CALL POPREAL4(tco2(2))
  co2expb(k, 2) = co2expb(k, 2) + tco2(2)*tco2b(2)
  tco2b(2) = co2exp(k, 2)*tco2b(2)
  tco2b(1) = tco2b(1) + 0.2673*xxb
  CALL POPREAL4(tco2(1))
  co2expb(k, 1) = co2expb(k, 1) + tco2(1)*tco2b(1)
  tco2b(1) = co2exp(k, 1)*tco2b(1)
  CALL POPREAL4(tran)
  tconb(1) = tconb(1) + tran*tranb
  tranb = tcon(1)*tranb
  CALL POPREAL4(tcon(1))
  conexpb(k) = conexpb(k) + tcon(1)*tconb(1)
  tconb(1) = conexp(k)*tconb(1)
  xxb = tranb
  th2ob(5) = th2ob(5) + 0.0119*xxb
  CALL POPREAL4(th2o(5))
  h2oexpb(k, 5) = h2oexpb(k, 5) + th2o(5)*th2ob(5)
  th2ob(5) = h2oexp(k, 5)*th2ob(5)
  th2ob(4) = th2ob(4) + 0.0798*xxb
  CALL POPREAL4(th2o(4))
  h2oexpb(k, 4) = h2oexpb(k, 4) + th2o(4)*th2ob(4)
  th2ob(4) = h2oexp(k, 4)*th2ob(4)
  th2ob(3) = th2ob(3) + 0.1326*xxb
  CALL POPREAL4(th2o(3))
  h2oexpb(k, 3) = h2oexpb(k, 3) + th2o(3)*th2ob(3)
  th2ob(3) = h2oexp(k, 3)*th2ob(3)
  th2ob(2) = th2ob(2) + 0.4604*xxb
  CALL POPREAL4(th2o(2))
  h2oexpb(k, 2) = h2oexpb(k, 2) + th2o(2)*th2ob(2)
  th2ob(2) = h2oexp(k, 2)*th2ob(2)
  th2ob(1) = th2ob(1) + 0.3153*xxb
  CALL POPREAL4(th2o(1))
  h2oexpb(k, 1) = h2oexpb(k, 1) + th2o(1)*th2ob(1)
  th2ob(1) = h2oexp(k, 1)*th2ob(1)
END SUBROUTINE B10KDIS_B

!  Differentiation of cldovlp in reverse (adjoint) mode:
!   gradient     of useful results: cldhi ett cldlw enn cldmd
!   with respect to varying inputs: cldhi ett cldlw enn cldmd
!mjs
SUBROUTINE CLDOVLP_B(np, k1, k2, ict, icb, icx, ncld, enn, ennb, ett, &
& ettb, cldhi, cldhib, cldmd, cldmdb, cldlw, cldlwb)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: np, k1, k2, ict, icb, icx(0:np), ncld(3)
  REAL, INTENT(IN) :: enn(0:np), ett(0:np)
  REAL :: ennb(0:np), ettb(0:np)
  REAL, INTENT(INOUT) :: cldhi, cldmd, cldlw
  REAL, INTENT(INOUT) :: cldhib
  INTEGER :: i, j, k, km, kx
  INTEGER :: branch
  REAL, INTENT(INOUT) :: cldmdb
  REAL, INTENT(INOUT) :: cldlwb
  km = k2 - 1
  IF (km .LT. ict) THEN
! do high clouds
    kx = ncld(1)
    IF (kx .EQ. 1 .OR. cldhi .EQ. 0.) THEN
      ennb(km) = ennb(km) + cldhib
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      CALL PUSHREAL4(cldhi)
      cldhi = 0.0
      IF (kx .NE. 0) THEN
        DO k=ict-kx,ict-1
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            CALL PUSHREAL4(cldhi)
            cldhi = enn(j) + ett(j)*cldhi
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        DO k=ict-1,ict-kx,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            j = icx(k)
            CALL POPREAL4(cldhi)
            ennb(j) = ennb(j) + cldhib
            ettb(j) = ettb(j) + cldhi*cldhib
            cldhib = ett(j)*cldhib
          END IF
        END DO
      END IF
      CALL POPREAL4(cldhi)
    END IF
    cldhib = 0.0
  ELSE IF (km .GE. ict .AND. km .LT. icb) THEN
! do middle clouds
    kx = ncld(2)
    IF (kx .EQ. 1 .OR. cldmd .EQ. 0.) THEN
      ennb(km) = ennb(km) + cldmdb
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      CALL PUSHREAL4(cldmd)
      cldmd = 0.0
      IF (kx .NE. 0) THEN
        DO k=icb-kx,icb-1
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            CALL PUSHREAL4(cldmd)
            cldmd = enn(j) + ett(j)*cldmd
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        DO k=icb-1,icb-kx,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            j = icx(k)
            CALL POPREAL4(cldmd)
            ennb(j) = ennb(j) + cldmdb
            ettb(j) = ettb(j) + cldmd*cldmdb
            cldmdb = ett(j)*cldmdb
          END IF
        END DO
      END IF
      CALL POPREAL4(cldmd)
    END IF
    cldmdb = 0.0
  ELSE
! do low clouds
    kx = ncld(3)
    IF (kx .EQ. 1 .OR. cldlw .EQ. 0.) THEN
      ennb(km) = ennb(km) + cldlwb
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      CALL PUSHREAL4(cldlw)
      cldlw = 0.0
      IF (kx .NE. 0) THEN
        DO k=np+1-kx,np
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            CALL PUSHREAL4(cldlw)
            cldlw = enn(j) + ett(j)*cldlw
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        DO k=np,np+1-kx,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            j = icx(k)
            CALL POPREAL4(cldlw)
            ennb(j) = ennb(j) + cldlwb
            ettb(j) = ettb(j) + cldlw*cldlwb
            cldlwb = ett(j)*cldlwb
          END IF
        END DO
      END IF
      CALL POPREAL4(cldlw)
    END IF
    cldlwb = 0.0
  END IF
END SUBROUTINE CLDOVLP_B

!  Differentiation of getirtau1 in reverse (adjoint) mode:
!   gradient     of useful results: hydromets tcldlyr fcld enn
!                reff
!   with respect to varying inputs: hydromets tcldlyr fcld enn
!                reff
SUBROUTINE GETIRTAU1_B(ib, nlevs, dp, fcld, fcldb, reff, reffb, &
& hydromets, hydrometsb, tcldlyr, tcldlyrb, enn, ennb, aib_ir1, awb_ir1&
& , aiw_ir1, aww_ir1, aig_ir1, awg_ir1, cons_grav)
  IMPLICIT NONE
! !INPUT PARAMETERS:
!  Band number
  INTEGER, INTENT(IN) :: ib
!  Number of levels
  INTEGER, INTENT(IN) :: nlevs
!  Delta pressure in Pa
  REAL, INTENT(IN) :: dp(nlevs)
!  Cloud fraction (used sometimes)
  REAL, INTENT(IN) :: fcld(nlevs)
  REAL :: fcldb(nlevs)
!  Effective radius (microns)
  REAL, INTENT(IN) :: reff(nlevs, 4)
  REAL :: reffb(nlevs, 4)
!  Hydrometeors (kg/kg)
  REAL, INTENT(IN) :: hydromets(nlevs, 4)
  REAL :: hydrometsb(nlevs, 4)
  REAL, INTENT(IN) :: aib_ir1(3, 10), awb_ir1(4, 10), aiw_ir1(4, 10)
  REAL, INTENT(IN) :: aww_ir1(4, 10), aig_ir1(4, 10), awg_ir1(4, 10)
  REAL, INTENT(IN) :: cons_grav
! !OUTPUT PARAMETERS:
!  Flux transmissivity?
  REAL :: tcldlyr(0:nlevs)
  REAL :: tcldlyrb(0:nlevs)
!  Flux transmissivity of a cloud layer?
  REAL :: enn(0:nlevs)
  REAL :: ennb(0:nlevs)
! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for infrared wavelengths
!  Slots for reff, hydrometeors and tauall are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
!
!  In the below calculations, the constants used in the tau calculation are in 
!  m$^2$ g$^{-1}$ and m$^2$ g$^{-1}$ $\mu$m. Thus, we must convert the kg contained in the 
!  pressure (Pa = kg m$^{-1}$ s$^{-2}$) to grams.
!
! !REVISION HISTORY: 
!    2011.11.18   MAT moved to Radiation_Shared and revised arg list, units
!
!EOP
!------------------------------------------------------------------------------
!BOC
  INTEGER :: k
  REAL :: taucld1, taucld2, taucld3, taucld4
  REAL :: taucld1b, taucld2b, taucld3b, taucld4b
  REAL :: g1, g2, g3, g4, gg
  REAL :: g1b, g2b, g3b, g4b, ggb
  REAL :: w1, w2, w3, w4, ww
  REAL :: w1b, w2b, w3b, w4b, wwb
  REAL :: ff, tauc
  REAL :: ffb, taucb
  REAL :: reff_snow
  REAL :: reff_snowb
  INTRINSIC MIN
  INTRINSIC ABS
  INTRINSIC MAX
  INTRINSIC EXP
!-----Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
!     Rain optical thickness is set to 0.00307 /(gm/m**2).
!     It is for a specific drop size distribution provided by Q. Fu.
  INTEGER :: branch
  REAL :: temp3
  REAL :: temp2
  REAL :: temp1
  REAL :: temp0
  REAL :: tempb9
  REAL :: tempb8
  REAL :: tempb7
  REAL :: tempb6
  REAL :: tempb5
  REAL :: tempb4
  REAL :: tempb3
  REAL :: tempb2
  REAL :: tempb1
  REAL :: tempb0
  REAL :: tempb12
  REAL :: tempb11
  REAL :: tempb10
  REAL :: temp16
  REAL :: temp15
  REAL :: temp14
  REAL :: temp13
  REAL :: temp12
  REAL :: temp11
  REAL :: temp10
  REAL :: max1b
  REAL :: tempb
  REAL :: abs0
  REAL :: temp
  REAL :: max1
  REAL :: temp9
  REAL :: temp8
  REAL :: temp7
  REAL :: temp6
  REAL :: temp5
  REAL :: temp4
  DO k=1,nlevs
    IF (reff(k, 1) .LE. 0.0) THEN
      CALL PUSHREAL4(taucld1)
      taucld1 = 0.0
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL4(taucld1)
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*(aib_ir1(1, ib)+&
&       aib_ir1(2, ib)/reff(k, 1)**aib_ir1(3, ib))
      CALL PUSHCONTROL1B(0)
    END IF
    CALL PUSHREAL4(taucld2)
    taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_ir1(1, ib)+(&
&     awb_ir1(2, ib)+(awb_ir1(3, ib)+awb_ir1(4, ib)*reff(k, 2))*reff(k, &
&     2))*reff(k, 2))
    taucld3 = 0.00307*(dp(k)*1.0e3/cons_grav*hydromets(k, 3))
    IF (reff(k, 4) .GT. 112.0) THEN
      CALL PUSHREAL4(reff_snow)
      reff_snow = 112.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL4(reff_snow)
      reff_snow = reff(k, 4)
      CALL PUSHCONTROL1B(1)
    END IF
    IF (reff_snow .LE. 0.0) THEN
      CALL PUSHREAL4(taucld4)
      taucld4 = 0.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL4(taucld4)
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*(aib_ir1(1, ib)+&
&       aib_ir1(2, ib)/reff_snow**aib_ir1(3, ib))
      CALL PUSHCONTROL1B(1)
    END IF
!-----Compute cloud single-scattering albedo and asymmetry factor for
!     a mixture of ice particles and liquid drops following 
!     Eqs. (6.5), (6.6), (6.15) and (6.16).
!     Single-scattering albedo and asymmetry factor of rain are set
!     to 0.54 and 0.95, respectively, based on the information provided
!     by Prof. Qiang Fu.
    CALL PUSHREAL4(tauc)
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      CALL PUSHREAL4(w1)
      w1 = taucld1*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff(k, 1))*reff(k, 1))*reff(k, 1))
      CALL PUSHREAL4(w2)
      w2 = taucld2*(aww_ir1(1, ib)+(aww_ir1(2, ib)+(aww_ir1(3, ib)+&
&       aww_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reff(k, 2))
      w3 = taucld3*0.54
      CALL PUSHREAL4(w4)
      w4 = taucld4*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff_snow)*reff_snow)*reff_snow)
      ww = (w1+w2+w3+w4)/tauc
      CALL PUSHREAL4(g1)
      g1 = w1*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(4&
&       , ib)*reff(k, 1))*reff(k, 1))*reff(k, 1))
      CALL PUSHREAL4(g2)
      g2 = w2*(awg_ir1(1, ib)+(awg_ir1(2, ib)+(awg_ir1(3, ib)+awg_ir1(4&
&       , ib)*reff(k, 2))*reff(k, 2))*reff(k, 2))
      g3 = w3*0.95
      CALL PUSHREAL4(g4)
      g4 = w4*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(4&
&       , ib)*reff_snow)*reff_snow)*reff_snow)
      IF (w1 + w2 + w3 + w4 .GE. 0.) THEN
        abs0 = w1 + w2 + w3 + w4
      ELSE
        abs0 = -(w1+w2+w3+w4)
      END IF
      IF (abs0 .GT. 0.0) THEN
        CALL PUSHREAL4(gg)
        gg = (g1+g2+g3+g4)/(w1+w2+w3+w4)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHREAL4(gg)
        gg = 0.5
        CALL PUSHCONTROL1B(0)
      END IF
!-----Parameterization of LW scattering following Eqs. (6.11)
!     and (6.12). 
      ff = 0.5 + (0.3739+(0.0076+0.1185*gg)*gg)*gg
      IF (1. - ww*ff .LT. 0.0) THEN
        CALL PUSHREAL4(max1)
        max1 = 0.0
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHREAL4(max1)
        max1 = 1. - ww*ff
        CALL PUSHCONTROL1B(1)
      END IF
!ALT: temporary protection against negative cloud optical thickness
      CALL PUSHREAL4(tauc)
      tauc = max1*tauc
!-----compute cloud diffuse transmittance. It is approximated by using 
!     a diffusivity factor of 1.66.
      CALL PUSHREAL4(tcldlyr(k))
      tcldlyr(k) = EXP(-(1.66*tauc))
! N in the documentation (6.13)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO k=nlevs,1,-1
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      ennb(k) = 0.0
      tcldlyrb(k) = 0.0
      reff_snowb = 0.0
      taucld1b = 0.0
      taucld2b = 0.0
      taucld3b = 0.0
      taucld4b = 0.0
      taucb = 0.0
    ELSE
      fcldb(k) = fcldb(k) + (1.0-tcldlyr(k))*ennb(k)
      tcldlyrb(k) = tcldlyrb(k) - fcld(k)*ennb(k)
      ennb(k) = 0.0
      CALL POPREAL4(tcldlyr(k))
      taucb = -(EXP(-(1.66*tauc))*1.66*tcldlyrb(k))
      tcldlyrb(k) = 0.0
      CALL POPREAL4(tauc)
      max1b = tauc*taucb
      taucb = max1*taucb
      taucld3 = 0.00307*(dp(k)*1.0e3/cons_grav*hydromets(k, 3))
      w3 = taucld3*0.54
      ww = (w1+w2+w3+w4)/tauc
      ff = 0.5 + (0.3739+(0.0076+0.1185*gg)*gg)*gg
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4(max1)
        wwb = 0.0
        ffb = 0.0
      ELSE
        CALL POPREAL4(max1)
        wwb = -(ff*max1b)
        ffb = -(ww*max1b)
      END IF
      ggb = ((0.1185*gg+0.0076)*gg+gg*(0.1185*gg+0.0076)+gg**2*0.1185+&
&       0.3739)*ffb
      g3 = w3*0.95
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4(gg)
        w1b = 0.0
        w2b = 0.0
        w3b = 0.0
        w4b = 0.0
        g1b = 0.0
        g2b = 0.0
        g3b = 0.0
        g4b = 0.0
      ELSE
        CALL POPREAL4(gg)
        tempb11 = ggb/(w1+w2+w3+w4)
        tempb12 = -((g1+g2+g3+g4)*tempb11/(w1+w2+w3+w4))
        g1b = tempb11
        g2b = tempb11
        g3b = tempb11
        g4b = tempb11
        w1b = tempb12
        w2b = tempb12
        w3b = tempb12
        w4b = tempb12
      END IF
      tempb5 = wwb/tauc
      CALL POPREAL4(g4)
      temp16 = aig_ir1(3, ib) + aig_ir1(4, ib)*reff_snow
      temp15 = aig_ir1(2, ib) + temp16*reff_snow
      tempb4 = w4*g4b
      w4b = w4b + tempb5 + (aig_ir1(1, ib)+temp15*reff_snow)*g4b
      w3b = w3b + tempb5 + 0.95*g3b
      CALL POPREAL4(g2)
      temp14 = awg_ir1(3, ib) + awg_ir1(4, ib)*reff(k, 2)
      temp13 = awg_ir1(2, ib) + temp14*reff(k, 2)
      tempb7 = w2*reff(k, 2)*g2b
      w2b = w2b + tempb5 + (awg_ir1(1, ib)+temp13*reff(k, 2))*g2b
      reffb(k, 2) = reffb(k, 2) + w2*temp13*g2b + (temp14+reff(k, 2)*&
&       awg_ir1(4, ib))*tempb7
      CALL POPREAL4(g1)
      temp12 = aig_ir1(3, ib) + aig_ir1(4, ib)*reff(k, 1)
      temp11 = aig_ir1(2, ib) + temp12*reff(k, 1)
      tempb8 = w1*reff(k, 1)*g1b
      w1b = w1b + tempb5 + (aig_ir1(1, ib)+temp11*reff(k, 1))*g1b
      reffb(k, 1) = reffb(k, 1) + w1*temp11*g1b + (temp12+reff(k, 1)*&
&       aig_ir1(4, ib))*tempb8
      taucb = taucb - (w1+w2+w3+w4)*tempb5/tauc
      CALL POPREAL4(w4)
      temp10 = aiw_ir1(3, ib) + aiw_ir1(4, ib)*reff_snow
      temp9 = aiw_ir1(2, ib) + temp10*reff_snow
      tempb6 = taucld4*w4b
      reff_snowb = (temp9+reff_snow*temp10+reff_snow**2*aiw_ir1(4, ib))*&
&       tempb6 + (temp15+reff_snow*temp16+reff_snow**2*aig_ir1(4, ib))*&
&       tempb4
      taucld4b = (aiw_ir1(1, ib)+temp9*reff_snow)*w4b
      taucld3b = 0.54*w3b
      CALL POPREAL4(w2)
      temp8 = aww_ir1(3, ib) + aww_ir1(4, ib)*reff(k, 2)
      temp7 = aww_ir1(2, ib) + temp8*reff(k, 2)
      tempb9 = taucld2*reff(k, 2)*w2b
      taucld2b = (aww_ir1(1, ib)+temp7*reff(k, 2))*w2b
      reffb(k, 2) = reffb(k, 2) + taucld2*temp7*w2b + (temp8+reff(k, 2)*&
&       aww_ir1(4, ib))*tempb9
      CALL POPREAL4(w1)
      temp6 = aiw_ir1(3, ib) + aiw_ir1(4, ib)*reff(k, 1)
      temp5 = aiw_ir1(2, ib) + temp6*reff(k, 1)
      tempb10 = taucld1*reff(k, 1)*w1b
      taucld1b = (aiw_ir1(1, ib)+temp5*reff(k, 1))*w1b
      reffb(k, 1) = reffb(k, 1) + taucld1*temp5*w1b + (temp6+reff(k, 1)*&
&       aiw_ir1(4, ib))*tempb10
    END IF
    CALL POPREAL4(tauc)
    taucld1b = taucld1b + taucb
    taucld2b = taucld2b + taucb
    taucld3b = taucld3b + taucb
    taucld4b = taucld4b + taucb
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(taucld4)
    ELSE
      CALL POPREAL4(taucld4)
      temp4 = reff_snow**aib_ir1(3, ib)
      temp3 = aib_ir1(2, ib)/temp4
      tempb3 = dp(k)*1.0e3*taucld4b
      hydrometsb(k, 4) = hydrometsb(k, 4) + (aib_ir1(1, ib)+temp3)*&
&       tempb3/cons_grav
      IF (.NOT.(reff_snow .LE. 0.0 .AND. (aib_ir1(3, ib) .EQ. 0.0 .OR. &
&         aib_ir1(3, ib) .NE. INT(aib_ir1(3, ib))))) reff_snowb = &
&         reff_snowb - temp3*hydromets(k, 4)*aib_ir1(3, ib)*reff_snow**(&
&         aib_ir1(3, ib)-1)*tempb3/(temp4*cons_grav)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(reff_snow)
    ELSE
      CALL POPREAL4(reff_snow)
      reffb(k, 4) = reffb(k, 4) + reff_snowb
    END IF
    hydrometsb(k, 3) = hydrometsb(k, 3) + dp(k)*1.0e3*0.00307*taucld3b/&
&     cons_grav
    CALL POPREAL4(taucld2)
    temp2 = awb_ir1(3, ib) + awb_ir1(4, ib)*reff(k, 2)
    temp1 = awb_ir1(2, ib) + temp2*reff(k, 2)
    tempb0 = dp(k)*1.0e3*taucld2b
    tempb1 = hydromets(k, 2)*tempb0/cons_grav
    tempb2 = reff(k, 2)*tempb1
    hydrometsb(k, 2) = hydrometsb(k, 2) + (awb_ir1(1, ib)+temp1*reff(k, &
&     2))*tempb0/cons_grav
    reffb(k, 2) = reffb(k, 2) + temp1*tempb1 + (temp2+reff(k, 2)*awb_ir1&
&     (4, ib))*tempb2
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(taucld1)
      temp0 = reff(k, 1)**aib_ir1(3, ib)
      temp = aib_ir1(2, ib)/temp0
      tempb = dp(k)*1.0e3*taucld1b
      hydrometsb(k, 1) = hydrometsb(k, 1) + (aib_ir1(1, ib)+temp)*tempb/&
&       cons_grav
      IF (.NOT.(reff(k, 1) .LE. 0.0 .AND. (aib_ir1(3, ib) .EQ. 0.0 .OR. &
&         aib_ir1(3, ib) .NE. INT(aib_ir1(3, ib))))) reffb(k, 1) = reffb&
&         (k, 1) - temp*hydromets(k, 1)*aib_ir1(3, ib)*reff(k, 1)**(&
&         aib_ir1(3, ib)-1)*tempb/(temp0*cons_grav)
    ELSE
      CALL POPREAL4(taucld1)
    END IF
  END DO
  ennb(0) = 0.0
  tcldlyrb(0) = 0.0
END SUBROUTINE GETIRTAU1_B

!  Differentiation of sorad in reverse (adjoint) mode:
!   gradient     of useful results: wa_dev fcld_dev cwc_dev ta_dev
!                flx_dev reff_dev oa_dev
!   with respect to varying inputs: ssaa_dev wa_dev fcld_dev cwc_dev
!                ta_dev asya_dev reff_dev taua_dev oa_dev
!Number of soundings
!Number of model levels
!Number of bands
!Cosine of solar zenith angle
!Pressure (Pa)
!Temperature (K)
!Specific humidity (kgkg^-1)
!Ozone (kgkg^-1)
!CO2 (pppv)
!Cloud properties (kgkg^-1)
!Cloud fractions (1)
!Level index separating high and middle clouds
!Level index separating middle and low clouds
!Moisture refflectivity properties
!Solar UV constant
!Solar IR constant
!Aerosol optical thickness 
!Aerosol single scattering albedo
!Aerosol asymmetry factor
!Surface reflectivity in the UV+par for beam insolation
!Surface reflectivity in the UV+par for diffuse insolation
!Surface reflectivity in the near-ir region for beam insolation
!Surface reflectivity in the near-ir region for diffuse insolation
!Outputs
!
!Constants            
SUBROUTINE SORAD_B(m, np, nb, cosz_dev, pl_dev, ta_dev, ta_devb, wa_dev&
& , wa_devb, oa_dev, oa_devb, co2, cwc_dev, cwc_devb, fcld_dev, &
& fcld_devb, ict, icb, reff_dev, reff_devb, hk_uv, hk_ir, taua_dev, &
& taua_devb, ssaa_dev, ssaa_devb, asya_dev, asya_devb, rsuvbm_dev, &
& rsuvdf_dev, rsirbm_dev, rsirdf_dev, flx_dev, flx_devb, cons_grav, &
& wk_uv, zk_uv, ry_uv, xk_ir, ry_ir, cah, coa, aig_uv, awg_uv, arg_uv, &
& aib_uv, awb_uv, arb_uv, aib_nir, awb_nir, arb_nir, aia_nir, awa_nir, &
& ara_nir, aig_nir, awg_nir, arg_nir, caib, caif)
  IMPLICIT NONE
! Parameters
! ----------
  INTEGER, PARAMETER :: nu=43
  INTEGER, PARAMETER :: nw=37
  INTEGER, PARAMETER :: nx=62
  INTEGER, PARAMETER :: ny=101
  INTEGER, PARAMETER :: nband_uv=5
  INTEGER, PARAMETER :: nk_ir=10
  INTEGER, PARAMETER :: nband_ir=3
  INTEGER, PARAMETER :: nband=nband_uv+nband_ir
  REAL, PARAMETER :: dsm=0.602
!-----input values
!-----input parameters
  INTEGER :: m, np, ict, icb, nb
  REAL :: cosz_dev(m), pl_dev(m, np+1), ta_dev(m, np), wa_dev(m, np), &
& oa_dev(m, np), co2
  REAL :: ta_devb(m, np), wa_devb(m, np), oa_devb(m, np)
  REAL :: cwc_dev(m, np, 4), fcld_dev(m, np), reff_dev(m, np, 4), hk_uv(&
& 5), hk_ir(3, 10)
  REAL :: cwc_devb(m, np, 4), fcld_devb(m, np), reff_devb(m, np, 4)
  REAL :: rsuvbm_dev(m), rsuvdf_dev(m), rsirbm_dev(m), rsirdf_dev(m)
  REAL :: taua_dev(m, np, nb)
  REAL :: taua_devb(m, np, nb)
  REAL :: ssaa_dev(m, np, nb)
  REAL :: ssaa_devb(m, np, nb)
  REAL :: asya_dev(m, np, nb)
  REAL :: asya_devb(m, np, nb)
  LOGICAL :: overcast
! Constants
  REAL, INTENT(IN) :: wk_uv(5), zk_uv(5), ry_uv(5)
  REAL, INTENT(IN) :: xk_ir(10), ry_ir(3)
  REAL, INTENT(IN) :: cah(43, 37), coa(62, 101)
  REAL, INTENT(IN) :: aig_uv(3), awg_uv(3), arg_uv(3)
  REAL, INTENT(IN) :: aib_uv, awb_uv(2), arb_uv(2)
  REAL, INTENT(IN) :: aib_nir, awb_nir(3, 2), arb_nir(3, 2)
  REAL, INTENT(IN) :: aia_nir(3, 3), awa_nir(3, 3), ara_nir(3, 3)
  REAL, INTENT(IN) :: aig_nir(3, 3), awg_nir(3, 3), arg_nir(3, 3)
  REAL, INTENT(IN) :: caib(11, 9, 11), caif(9, 11)
  REAL, INTENT(IN) :: cons_grav
!-----output parameters
  REAL :: flx_dev(m, np+1), flc_dev(m, np+1)
  REAL :: flx_devb(m, np+1)
  REAL :: flxu_dev(m, np+1), flcu_dev(m, np+1)
  REAL :: fdiruv_dev(m), fdifuv_dev(m)
  REAL :: fdirpar_dev(m), fdifpar_dev(m)
  REAL :: fdirir_dev(m), fdifir_dev(m)
  REAL :: flx_sfc_band_dev(m, nband)
!-----temporary arrays
  INTEGER :: i, j, k, l, in, ntop
  REAL :: dp(np), wh(np), oh(np)
  REAL :: whb(np), ohb(np)
  REAL :: scal(np)
  REAL :: swh(np+1), so2(np+1), df(0:np+1)
  REAL :: swhb(np+1), dfb(0:np+1)
  REAL :: scal0, wvtoa, o3toa, pa
  REAL :: wvtoab, o3toab
  REAL :: snt, cnt, x, xx4, xtoa
  REAL :: xx4b
  REAL :: dp_pa(np)
!-----parameters for co2 transmission tables
  REAL :: w1, dw, u1, du
  INTEGER :: ib, rc
  REAL :: tauclb(np), tauclf(np), asycl(np)
  REAL :: tauclbb(np), tauclfb(np), asyclb(np)
  REAL :: taubeam(np, 4), taudiff(np, 4)
  REAL :: taubeamb(np, 4), taudiffb(np, 4)
  REAL :: fcld_col(np)
  REAL :: fcld_colb(np)
  REAL :: cwc_col(np, 4)
  REAL :: cwc_colb(np, 4)
  REAL :: reff_col(np, 4)
  REAL :: reff_colb(np, 4)
  REAL :: taurs, tauoz, tauwv
  REAL :: tauozb, tauwvb
  REAL :: tausto, ssatau, asysto
  REAL :: taustob, ssataub, asystob
  REAL :: tautob, ssatob, asytob
  REAL :: tautobb, ssatobb, asytobb
  REAL :: tautof, ssatof, asytof
  REAL :: tautofb, ssatofb, asytofb
  REAL :: rr(0:np+1, 2), tt(0:np+1, 2), td(0:np+1, 2)
  REAL :: rrb(0:np+1, 2), ttb(0:np+1, 2), tdb(0:np+1, 2)
  REAL :: rs(0:np+1, 2), ts(0:np+1, 2)
  REAL :: rsb(0:np+1, 2), tsb(0:np+1, 2)
  REAL :: fall(np+1), fclr(np+1), fsdir, fsdif
  REAL :: fallb(np+1)
  REAL :: fupa(np+1), fupc(np+1)
  REAL :: cc1, cc2, cc3
  REAL :: cc1b, cc2b, cc3b
  REAL :: rrt, ttt, tdt, rst, tst
  REAL :: rrtb, tttb, tdtb, rstb, tstb
  INTEGER :: iv, ik
  REAL :: ssacl(np)
  REAL :: ssaclb(np)
  INTEGER :: im
  INTEGER :: ic, iw
  REAL :: ulog, wlog, dc, dd, x0, x1, x2, y0, y1, y2, du2, dw2
  REAL :: wlogb, ddb, x2b, y2b
  INTEGER :: ih
!if (overcast == true) then
!real :: rra(0:np+1),rxa(0:np+1)
!real :: ttaold,tdaold,rsaold
!real :: ttanew,tdanew,rsanew 
!else
  REAL :: rra(0:np+1, 2, 2), tta(0:np, 2, 2)
  REAL :: rrab(0:np+1, 2, 2), ttab(0:np, 2, 2)
  REAL :: tda(0:np, 2, 2)
  REAL :: tdab(0:np, 2, 2)
  REAL :: rsa(0:np, 2, 2), rxa(0:np+1, 2, 2)
  REAL :: rsab(0:np, 2, 2), rxab(0:np+1, 2, 2)
!endif
  REAL :: flxdn
  REAL :: flxdnb
  REAL :: fdndir, fdndif, fupdif
  REAL :: fdndirb, fdndifb, fupdifb
  REAL :: denm, yy
  REAL :: denmb, yyb
  INTEGER :: is
  REAL :: ch, cm, ct
  REAL :: chb, cmb, ctb
  INTEGER :: foundtop
  REAL :: dftop
  REAL :: dftopb
!-----Variables for aerosols
  INTEGER :: ii, jj, irhp1, an
  REAL :: dum
  REAL :: dumb
  INTRINSIC MAX
  INTRINSIC EXP
  INTRINSIC MIN
  INTRINSIC SQRT
  INTRINSIC REAL
  INTRINSIC LOG10
  INTRINSIC INT
  INTRINSIC ABS
  INTRINSIC EPSILON
  REAL :: result1
  INTEGER :: branch
  REAL :: temp3
  REAL :: temp29
  REAL :: tempb52
  REAL :: temp2
  REAL :: temp28
  REAL :: tempb51
  REAL :: temp1
  REAL :: temp27
  REAL :: tempb50
  REAL :: temp0
  REAL :: temp26
  REAL :: temp25
  REAL :: temp24
  REAL :: temp23
  REAL :: temp22
  REAL :: temp21
  REAL :: temp20
  REAL :: tempb9
  REAL :: tempb8
  REAL :: tempb7
  REAL :: tempb6
  REAL :: tempb5
  REAL :: tempb4
  REAL :: tempb19
  REAL :: tempb3
  REAL :: tempb18
  REAL :: tempb2
  REAL :: tempb17
  REAL :: tempb1
  REAL :: tempb16
  REAL :: tempb0
  REAL :: tempb15
  REAL :: tempb14
  REAL :: tempb13
  REAL :: tempb12
  REAL :: tempb49
  REAL :: x6
  REAL :: tempb11
  REAL :: tempb48
  REAL :: x5
  REAL :: tempb10
  REAL :: tempb47
  REAL :: x4
  REAL :: tempb46
  REAL :: x3
  REAL :: tempb45
  REAL :: tempb44
  REAL :: tempb43
  REAL :: temp19
  REAL :: tempb42
  REAL :: temp18
  REAL :: tempb41
  REAL :: temp17
  REAL :: tempb40
  REAL :: temp16
  REAL :: temp15
  REAL :: temp14
  REAL :: temp13
  REAL :: temp12
  REAL :: temp11
  REAL :: temp10
  REAL :: tempb
  REAL :: tempb39
  REAL :: tempb38
  REAL :: tempb37
  REAL :: tempb36
  REAL :: tempb35
  REAL :: tempb34
  REAL :: tempb33
  REAL :: tempb32
  REAL :: tempb31
  REAL :: tempb30
  REAL :: x4b
  REAL :: temp38
  REAL :: temp37
  REAL :: temp36
  REAL :: temp35
  REAL :: temp34
  REAL :: temp33
  REAL :: temp32
  REAL :: temp31
  REAL :: temp30
  REAL :: abs0
  REAL :: tempb29
  REAL :: tempb28
  REAL :: tempb27
  REAL :: tempb26
  REAL :: tempb25
  REAL :: temp
  REAL :: tempb24
  REAL :: tempb23
  REAL :: tempb22
  REAL :: temp9
  REAL :: tempb21
  REAL :: temp8
  REAL :: tempb20
  REAL :: temp7
  REAL :: tempb56
  REAL :: temp6
  REAL :: tempb55
  REAL :: temp5
  REAL :: tempb54
  REAL :: temp4
  REAL :: tempb53
run_loop:DO i=1,m
    CALL PUSHINTEGER4(ntop)
    ntop = 0
!-----Beginning of sorad code
!-----wvtoa and o3toa are the water vapor and o3 amounts of the region 
!     above the pl(1) level.
!     snt is the secant of the solar zenith angle
    CALL PUSHREAL4(snt)
    snt = 1.0/cosz_dev(i)
    IF (pl_dev(i, 1) .LT. 1.e-3) THEN
      CALL PUSHREAL4(xtoa)
      xtoa = 1.e-3
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL4(xtoa)
      xtoa = pl_dev(i, 1)
      CALL PUSHCONTROL1B(1)
    END IF
    CALL PUSHREAL4(scal0)
    scal0 = xtoa*(0.5*xtoa/300.)**.8
    o3toa = 1.02*oa_dev(i, 1)*xtoa*466.7 + 1.0e-8
    wvtoa = 1.02*wa_dev(i, 1)*scal0*(1.0+0.00135*(ta_dev(i, 1)-240.)) + &
&     1.0e-9
    CALL PUSHREAL4(swh(1))
    swh(1) = wvtoa
    DO k=1,np
!-----compute layer thickness. indices for the surface level and
!     surface layer are np+1 and np, respectively.
      CALL PUSHREAL4(dp(k))
      dp(k) = pl_dev(i, k+1) - pl_dev(i, k)
! dp in pascals
      CALL PUSHREAL4(dp_pa(k))
      dp_pa(k) = dp(k)*100.
!-----compute scaled water vapor amount following Eqs. (3.3) and (3.5) 
!     unit is g/cm**2
!
      pa = 0.5*(pl_dev(i, k)+pl_dev(i, k+1))
      CALL PUSHREAL4(scal(k))
      scal(k) = dp(k)*(pa/300.)**.8
      wh(k) = 1.02*wa_dev(i, k)*scal(k)*(1.+0.00135*(ta_dev(i, k)-240.))&
&       + 1.e-9
      CALL PUSHREAL4(swh(k+1))
      swh(k+1) = swh(k) + wh(k)
!-----compute ozone amount, unit is (cm-atm)stp
!     the number 466.7 is the unit conversion factor
!     from g/cm**2 to (cm-atm)stp
      oh(k) = 1.02*oa_dev(i, k)*dp(k)*466.7 + 1.e-8
!-----Fill the reff, cwc, and fcld for the column
      CALL PUSHREAL4(fcld_col(k))
      fcld_col(k) = fcld_dev(i, k)
      DO l=1,4
        CALL PUSHREAL4(reff_col(k, l))
        reff_col(k, l) = reff_dev(i, k, l)
        CALL PUSHREAL4(cwc_col(k, l))
        cwc_col(k, l) = cwc_dev(i, k, l)
      END DO
    END DO
!-----Initialize temporary arrays to zero to avoid UMR
    CALL PUSHREAL4ARRAY(rr, (np+2)*2)
    rr = 0.0
    CALL PUSHREAL4ARRAY(tt, (np+2)*2)
    tt = 0.0
    CALL PUSHREAL4ARRAY(td, (np+2)*2)
    td = 0.0
    CALL PUSHREAL4ARRAY(rs, (np+2)*2)
    rs = 0.0
    CALL PUSHREAL4ARRAY(ts, (np+2)*2)
    ts = 0.0
    CALL PUSHREAL4ARRAY(rra, (np+2)*2**2)
    rra = 0.0
    CALL PUSHREAL4ARRAY(rxa, (np+2)*2**2)
    rxa = 0.0
!if( OVERCAST == .false. ) then
    CALL PUSHREAL4ARRAY(tta, (np+1)*2**2)
    tta = 0.0
    CALL PUSHREAL4ARRAY(tda, (np+1)*2**2)
    tda = 0.0
    CALL PUSHREAL4ARRAY(rsa, (np+1)*2**2)
    rsa = 0.0
!endif
!-----initialize fluxes for all-sky (flx), clear-sky (flc), and
!     flux reduction (df)
!
    DO k=1,np+1
      flx_dev(i, k) = 0.
    END DO
!-----Begin inline of SOLUV
!-----compute solar uv and par fluxes
!-----initialize fdiruv, fdifuv, surface reflectances and transmittances.
!     the reflectance and transmittance of the clear and cloudy portions
!     of a layer are denoted by 1 and 2, respectively.
!     cc is the maximum cloud cover in each of the high, middle, and low
!     cloud groups.
!     1/dsm=1/cos(53) = 1.66
    rr(np+1, 1) = rsuvbm_dev(i)
    rr(np+1, 2) = rsuvbm_dev(i)
    rs(np+1, 1) = rsuvdf_dev(i)
    rs(np+1, 2) = rsuvdf_dev(i)
    td(np+1, 1) = 0.0
    td(np+1, 2) = 0.0
    tt(np+1, 1) = 0.0
    tt(np+1, 2) = 0.0
    ts(np+1, 1) = 0.0
    ts(np+1, 2) = 0.0
    rr(0, 1) = 0.0
    rr(0, 2) = 0.0
    rs(0, 1) = 0.0
    rs(0, 2) = 0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
    tt(0, 1) = 1.0
    tt(0, 2) = 1.0
    ts(0, 1) = 1.0
    ts(0, 2) = 1.0
    cc1 = 0.0
    CALL PUSHREAL4(cc2)
    cc2 = 0.0
    CALL PUSHREAL4(cc3)
    cc3 = 0.0
!-----options for scaling cloud optical thickness
!if ( OVERCAST == .true. ) then
!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud asymmetry factor
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.
!         call getvistau1(np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,0,0,&
!                        taubeam,taudiff,asycl,                                  &
!                         aig_uv, awg_uv, arg_uv,                                 &
!                         aib_uv, awb_uv, arb_uv,                                 &
!                         aib_nir, awb_nir, arb_nir,                              &
!                         aia_nir, awa_nir, ara_nir,                              &
!                         aig_nir, awg_nir, arg_nir,                              &
!                         caib, caif,                                             &
!                         CONS_GRAV                                               )
!else
!-----Compute scaled cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud asymmetry factor
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.
    CALL PUSHREAL4ARRAY(asycl, np)
    CALL GETVISTAU1(np, cosz_dev(i), dp_pa, fcld_col, reff_col, cwc_col&
&             , ict, icb, taubeam, taudiff, asycl, aig_uv, awg_uv, &
&             arg_uv, aib_uv, awb_uv, arb_uv, aib_nir, awb_nir, arb_nir&
&             , aia_nir, awa_nir, ara_nir, aig_nir, awg_nir, arg_nir, &
&             caib, caif, cons_grav)
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
!     The cc1,2,3 are still needed in the flux calculations below
!MAT---DO NOT FUSE THIS LOOP
!MAT---Loop must run to completion so that cc[1,2,3] are correct.
    DO k=1,np
      IF (k .LT. ict) THEN
        IF (cc1 .LT. fcld_dev(i, k)) THEN
          cc1 = fcld_dev(i, k)
          CALL PUSHCONTROL3B(4)
        ELSE
          CALL PUSHCONTROL3B(5)
          cc1 = cc1
        END IF
      ELSE IF (k .LT. icb) THEN
        IF (cc2 .LT. fcld_dev(i, k)) THEN
          cc2 = fcld_dev(i, k)
          CALL PUSHCONTROL3B(2)
        ELSE
          CALL PUSHCONTROL3B(3)
          cc2 = cc2
        END IF
      ELSE IF (cc3 .LT. fcld_dev(i, k)) THEN
        cc3 = fcld_dev(i, k)
        CALL PUSHCONTROL3B(0)
      ELSE
        CALL PUSHCONTROL3B(1)
        cc3 = cc3
      END IF
    END DO
!MAT---DO NOT FUSE THIS LOOP
!endif !overcast
    DO k=1,np
      CALL PUSHREAL4(tauclb(k))
      tauclb(k) = taubeam(k, 1) + taubeam(k, 2) + taubeam(k, 3) + &
&       taubeam(k, 4)
      CALL PUSHREAL4(tauclf(k))
      tauclf(k) = taudiff(k, 1) + taudiff(k, 2) + taudiff(k, 3) + &
&       taudiff(k, 4)
    END DO
!-----integration over spectral bands
!-----Compute optical thickness, single-scattering albedo and asymmetry
!     factor for a mixture of "na" aerosol types. [Eqs. (4.16)-(4.18)]
    DO ib=1,nband_uv
!-----compute direct beam transmittances of the layer above pl(1)
      CALL PUSHREAL4(td(0, 1))
      td(0, 1) = EXP(-((wvtoa*wk_uv(ib)+o3toa*zk_uv(ib))/cosz_dev(i)))
      CALL PUSHREAL4(td(0, 2))
      td(0, 2) = td(0, 1)
      DO k=1,np
!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor (Eqs. 6.2-6.4)
        taurs = ry_uv(ib)*dp(k)
        tauoz = zk_uv(ib)*oh(k)
        tauwv = wk_uv(ib)*wh(k)
        tausto = taurs + tauoz + tauwv + taua_dev(i, k, ib) + 1.0e-7
        ssatau = ssaa_dev(i, k, ib) + taurs
        asysto = asya_dev(i, k, ib)
        CALL PUSHREAL4(tautob)
        tautob = tausto
        asytob = asysto/ssatau
        CALL PUSHREAL4(ssatob)
        ssatob = ssatau/tautob + 1.0e-8
        IF (ssatob .GT. 0.999999) THEN
          ssatob = 0.999999
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          ssatob = ssatob
        END IF
!-----for direct incident radiation
        CALL DELEDD(tautob, ssatob, asytob, cosz_dev(i), rrt, ttt, tdt)
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
        CALL DELEDD(tautob, ssatob, asytob, dsm, rst, tst, dum)
        CALL PUSHREAL4(rr(k, 1))
        rr(k, 1) = rrt
        CALL PUSHREAL4(tt(k, 1))
        tt(k, 1) = ttt
        CALL PUSHREAL4(td(k, 1))
        td(k, 1) = tdt
        CALL PUSHREAL4(rs(k, 1))
        rs(k, 1) = rst
        CALL PUSHREAL4(ts(k, 1))
        ts(k, 1) = tst
!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer
!-----for direct incident radiation
!     The effective layer optical properties. Eqs. (6.2)-(6.4)
        CALL PUSHREAL4(tautob)
        tautob = tausto + tauclb(k)
        CALL PUSHREAL4(ssatob)
        ssatob = (ssatau+tauclb(k))/tautob + 1.0e-8
        IF (ssatob .GT. 0.999999) THEN
          ssatob = 0.999999
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          ssatob = ssatob
        END IF
        asytob = (asysto+asycl(k)*tauclb(k))/(ssatob*tautob)
!-----for diffuse incident radiation
        CALL PUSHREAL4(tautof)
        tautof = tausto + tauclf(k)
        CALL PUSHREAL4(ssatof)
        ssatof = (ssatau+tauclf(k))/tautof + 1.0e-8
        IF (ssatof .GT. 0.999999) THEN
          ssatof = 0.999999
          CALL PUSHCONTROL1B(0)
        ELSE
          ssatof = ssatof
          CALL PUSHCONTROL1B(1)
        END IF
        asytof = (asysto+asycl(k)*tauclf(k))/(ssatof*tautof)
!-----for direct incident radiation
!     note that the cloud optical thickness is scaled differently 
!     for direct and diffuse insolation, Eqs. (7.3) and (7.4).
        CALL DELEDD(tautob, ssatob, asytob, cosz_dev(i), rrt, ttt, tdt)
!-----diffuse incident radiation is approximated by beam radiation 
!     with an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
        CALL DELEDD(tautof, ssatof, asytof, dsm, rst, tst, dum)
        CALL PUSHREAL4(rr(k, 2))
        rr(k, 2) = rrt
        CALL PUSHREAL4(tt(k, 2))
        tt(k, 2) = ttt
        CALL PUSHREAL4(td(k, 2))
        td(k, 2) = tdt
        CALL PUSHREAL4(rs(k, 2))
        rs(k, 2) = rst
        CALL PUSHREAL4(ts(k, 2))
        ts(k, 2) = tst
      END DO
!-----flux calculations
!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)
      DO k=1,np+1
        fall(k) = 0.0
      END DO
!if ( OVERCAST == .true. ) then
!-----Inline CLDFLXY
!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is either 0 or 1.
!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)
!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated by 
!         beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)
!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux
!-----ih=1 for clear sky; ih=2 for cloudy sky.
!-----First set is ih = 1
!            rra(np+1)=rr(np+1,1)
!            rxa(np+1)=rs(np+1,1)
!
!            do k=np,0,-1
!               denm=ts(k,1)/(1.-rs(k,1)*rxa(k+1))
!               rra(k)=rr(k,1)+(td(k,1)*rra(k+1)+(tt(k,1)-td(k,1))*rxa(k+1))*denm
!               rxa(k)=rs(k,1)+ts(k,1)*rxa(k+1)*denm
!            end do
!
!            do k=1,np+1
!               if (k <= np) then
!                  if (k == 1) then
!                     tdaold = td(0,1)
!                     ttaold = tt(0,1)
!                     rsaold = rs(0,1)
!
!                     tdanew = 0.0
!                     ttanew = 0.0
!                     rsanew = 0.0
!                  end if
!                  denm=ts(k,1)/(1.-rsaold*rs(k,1))
!                  tdanew=tdaold*td(k,1)
!                  ttanew=tdaold*tt(k,1)+(tdaold*rsaold*rr(k,1)+ttaold-tdaold)*denm
!                  rsanew=rs(k,1)+ts(k,1)*rsaold*denm
!               end if
!
!               denm=1./(1.-rsaold*rxa(k))
!               fdndir=tdaold
!               xx4=tdaold*rra(k)
!               yy=ttaold-tdaold
!               fdndif=(xx4*rsaold+yy)*denm
!               fupdif=(xx4+yy*rxa(k))*denm
!               flxdn=fdndir+fdndif-fupdif
!               fupc(k)=fupdif 
!               fclr(k)=flxdn
!
!               tdaold = tdanew
!               ttaold = ttanew
!               rsaold = rsanew
!
!               tdanew = 0.0
!               ttanew = 0.0
!               rsanew = 0.0
!            end do
!
!!-----Second set is ih = 2
!
!            rra(np+1)=rr(np+1,2)
!            rxa(np+1)=rs(np+1,2)
!
!            do k=np,0,-1
!               denm=ts(k,2)/(1.-rs(k,2)*rxa(k+1))
!               rra(k)=rr(k,2)+(td(k,2)*rra(k+1)+(tt(k,2)-td(k,2))*rxa(k+1))*denm
!               rxa(k)=rs(k,2)+ts(k,2)*rxa(k+1)*denm
!            end do
!
!            do k=1,np+1
!               if (k <= np) then
!                  if (k == 1) then
!                     tdaold = td(0,2)
!                     ttaold = tt(0,2)
!                     rsaold = rs(0,2)
!                     tdanew = 0.0
!                     ttanew = 0.0
!                     rsanew = 0.0
!                  end if
!                  denm=ts(k,2)/(1.-rsaold*rs(k,2))
!                  tdanew=tdaold*td(k,2)
!                  ttanew=tdaold*tt(k,2)+(tdaold*rsaold*rr(k,2)+ttaold-tdaold)*denm
!                  rsanew=rs(k,2)+ts(k,2)*rsaold*denm
!               end if
!
!               denm=1./(1.-rsaold*rxa(k))
!               fdndir=tdaold
!               xx4=tdaold*rra(k)
!               yy=ttaold-tdaold
!               fdndif=(xx4*rsaold+yy)*denm
!               fupdif=(xx4+yy*rxa(k))*denm
!               flxdn=fdndir+fdndif-fupdif
!
!               fupa(k)=fupdif
!               fall(k)=flxdn
!
!               tdaold = tdanew
!               ttaold = ttanew
!               rsaold = rsanew
!
!               tdanew = 0.0
!               ttanew = 0.0
!               rsanew = 0.0
!            end do
!
!            fsdir=fdndir
!            fsdif=fdndif
!
!!-----End CLDFLXY inline
!
!else
!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is allowed to be between 0 and 1.
!     the all-sky flux, fall is the summation inside the brackets
!     of Eq. (7.11)
!-----Inline CLDFLX
!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated 
!         by beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)
!-----To save memory space, tda, tta, and rsa are pre-computed 
!     for k<icb. The dimension of these parameters is (m,np,2,2). 
!     It would have been (m,np,2,2,2) if these parameters were 
!     computed for all k's.
!-----for high clouds
!     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition
      DO ih=1,2
        CALL PUSHREAL4(tda(0, ih, 1))
        tda(0, ih, 1) = td(0, ih)
        CALL PUSHREAL4(tta(0, ih, 1))
        tta(0, ih, 1) = tt(0, ih)
        CALL PUSHREAL4(rsa(0, ih, 1))
        rsa(0, ih, 1) = rs(0, ih)
        CALL PUSHREAL4(tda(0, ih, 2))
        tda(0, ih, 2) = td(0, ih)
        CALL PUSHREAL4(tta(0, ih, 2))
        tta(0, ih, 2) = tt(0, ih)
        CALL PUSHREAL4(rsa(0, ih, 2))
        rsa(0, ih, 2) = rs(0, ih)
        DO k=1,ict-1
          CALL PUSHREAL4(denm)
          denm = ts(k, ih)/(1.-rsa(k-1, ih, 1)*rs(k, ih))
          CALL PUSHREAL4(tda(k, ih, 1))
          tda(k, ih, 1) = tda(k-1, ih, 1)*td(k, ih)
          CALL PUSHREAL4(tta(k, ih, 1))
          tta(k, ih, 1) = tda(k-1, ih, 1)*tt(k, ih) + (tda(k-1, ih, 1)*&
&           rsa(k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)-tda(k-1, ih, 1))*&
&           denm
          CALL PUSHREAL4(rsa(k, ih, 1))
          rsa(k, ih, 1) = rs(k, ih) + ts(k, ih)*rsa(k-1, ih, 1)*denm
          CALL PUSHREAL4(tda(k, ih, 2))
          tda(k, ih, 2) = tda(k, ih, 1)
          CALL PUSHREAL4(tta(k, ih, 2))
          tta(k, ih, 2) = tta(k, ih, 1)
          CALL PUSHREAL4(rsa(k, ih, 2))
          rsa(k, ih, 2) = rsa(k, ih, 1)
        END DO
! k loop
!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition
        DO k=ict,icb-1
          DO im=1,2
            CALL PUSHREAL4(denm)
            denm = ts(k, im)/(1.-rsa(k-1, ih, im)*rs(k, im))
            CALL PUSHREAL4(tda(k, ih, im))
            tda(k, ih, im) = tda(k-1, ih, im)*td(k, im)
            CALL PUSHREAL4(tta(k, ih, im))
            tta(k, ih, im) = tda(k-1, ih, im)*tt(k, im) + (tda(k-1, ih, &
&             im)*rsa(k-1, ih, im)*rr(k, im)+tta(k-1, ih, im)-tda(k-1, &
&             ih, im))*denm
            CALL PUSHREAL4(rsa(k, ih, im))
            rsa(k, ih, im) = rs(k, im) + ts(k, im)*rsa(k-1, ih, im)*denm
          END DO
        END DO
      END DO
! im loop
! k loop
! ih loop
!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)
!-----To save memory space, rra and rxa are pre-computed for k>=icb.
!     the dimension of these parameters is (m,np,2,2). It would have
!     been (m,np,2,2,2) if these parameters were computed for all k's.
!-----for the low clouds
!     is=1 for clear-sky condition, is=2 for cloudy-sky condition
      DO is=1,2
        CALL PUSHREAL4(rra(np+1, 1, is))
        rra(np+1, 1, is) = rr(np+1, is)
        CALL PUSHREAL4(rxa(np+1, 1, is))
        rxa(np+1, 1, is) = rs(np+1, is)
        CALL PUSHREAL4(rra(np+1, 2, is))
        rra(np+1, 2, is) = rr(np+1, is)
        CALL PUSHREAL4(rxa(np+1, 2, is))
        rxa(np+1, 2, is) = rs(np+1, is)
        DO k=np,icb,-1
          CALL PUSHREAL4(denm)
          denm = ts(k, is)/(1.-rs(k, is)*rxa(k+1, 1, is))
          CALL PUSHREAL4(rra(k, 1, is))
          rra(k, 1, is) = rr(k, is) + (td(k, is)*rra(k+1, 1, is)+(tt(k, &
&           is)-td(k, is))*rxa(k+1, 1, is))*denm
          CALL PUSHREAL4(rxa(k, 1, is))
          rxa(k, 1, is) = rs(k, is) + ts(k, is)*rxa(k+1, 1, is)*denm
          CALL PUSHREAL4(rra(k, 2, is))
          rra(k, 2, is) = rra(k, 1, is)
          CALL PUSHREAL4(rxa(k, 2, is))
          rxa(k, 2, is) = rxa(k, 1, is)
        END DO
! k loop
!-----for middle clouds
        DO k=icb-1,ict,-1
          DO im=1,2
            CALL PUSHREAL4(denm)
            denm = ts(k, im)/(1.-rs(k, im)*rxa(k+1, im, is))
            CALL PUSHREAL4(rra(k, im, is))
            rra(k, im, is) = rr(k, im) + (td(k, im)*rra(k+1, im, is)+(tt&
&             (k, im)-td(k, im))*rxa(k+1, im, is))*denm
            CALL PUSHREAL4(rxa(k, im, is))
            rxa(k, im, is) = rs(k, im) + ts(k, im)*rxa(k+1, im, is)*denm
          END DO
        END DO
      END DO
! im loop
! k loop
! is loop
!-----integration over eight sky situations.
!     ih, im, is denote high, middle and low cloud groups.
      DO ih=1,2
!-----clear portion 
        IF (ih .EQ. 1) THEN
          CALL PUSHREAL4(ch)
          ch = 1.0 - cc1
!-----cloudy portion
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHREAL4(ch)
          ch = cc1
          CALL PUSHCONTROL1B(0)
        END IF
        DO im=1,2
!-----clear portion 
          IF (im .EQ. 1) THEN
            CALL PUSHREAL4(cm)
            cm = ch*(1.0-cc2)
!-----cloudy portion
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHREAL4(cm)
            cm = ch*cc2
            CALL PUSHCONTROL1B(0)
          END IF
          DO is=1,2
!-----clear portion 
            IF (is .EQ. 1) THEN
              CALL PUSHREAL4(ct)
              ct = cm*(1.0-cc3)
!-----cloudy portion
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHREAL4(ct)
              ct = cm*cc3
              CALL PUSHCONTROL1B(0)
            END IF
!-----add one layer at a time, going down.
            DO k=icb,np
              CALL PUSHREAL4(denm)
              denm = ts(k, is)/(1.-rsa(k-1, ih, im)*rs(k, is))
              CALL PUSHREAL4(tda(k, ih, im))
              tda(k, ih, im) = tda(k-1, ih, im)*td(k, is)
              CALL PUSHREAL4(tta(k, ih, im))
              tta(k, ih, im) = tda(k-1, ih, im)*tt(k, is) + (tda(k-1, ih&
&               , im)*rr(k, is)*rsa(k-1, ih, im)+tta(k-1, ih, im)-tda(k-&
&               1, ih, im))*denm
              CALL PUSHREAL4(rsa(k, ih, im))
              rsa(k, ih, im) = rs(k, is) + ts(k, is)*rsa(k-1, ih, im)*&
&               denm
            END DO
! k loop
!-----add one layer at a time, going up.
            DO k=ict-1,0,-1
              CALL PUSHREAL4(denm)
              denm = ts(k, ih)/(1.-rs(k, ih)*rxa(k+1, im, is))
              CALL PUSHREAL4(rra(k, im, is))
              rra(k, im, is) = rr(k, ih) + (td(k, ih)*rra(k+1, im, is)+(&
&               tt(k, ih)-td(k, ih))*rxa(k+1, im, is))*denm
              CALL PUSHREAL4(rxa(k, im, is))
              rxa(k, im, is) = rs(k, ih) + ts(k, ih)*rxa(k+1, im, is)*&
&               denm
            END DO
! k loop
!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux
            DO k=1,np+1
              CALL PUSHREAL4(denm)
              denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
              fdndir = tda(k-1, ih, im)
              xx4 = tda(k-1, ih, im)*rra(k, im, is)
              yy = tta(k-1, ih, im) - tda(k-1, ih, im)
              fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
              fupdif = (xx4+yy*rxa(k, im, is))*denm
              flxdn = fdndir + fdndif - fupdif
!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)
              fall(k) = fall(k) + flxdn*ct
            END DO
          END DO
        END DO
      END DO
! is loop
! im loop
! ih loop
!-----End CLDFLX inline
!endif !overcast
!-----flux integration, Eq. (6.1)
      DO k=1,np+1
        flx_dev(i, k) = flx_dev(i, k) + fall(k)*hk_uv(ib)
      END DO
    END DO
!-----Inline SOLIR
!-----compute and update solar ir fluxes
    CALL PUSHREAL4(rr(np+1, 1))
    rr(np+1, 1) = rsirbm_dev(i)
    CALL PUSHREAL4(rr(np+1, 2))
    rr(np+1, 2) = rsirbm_dev(i)
    CALL PUSHREAL4(rs(np+1, 1))
    rs(np+1, 1) = rsirdf_dev(i)
    CALL PUSHREAL4(rs(np+1, 2))
    rs(np+1, 2) = rsirdf_dev(i)
    CALL PUSHREAL4(td(np+1, 1))
    td(np+1, 1) = 0.0
    CALL PUSHREAL4(td(np+1, 2))
    td(np+1, 2) = 0.0
    CALL PUSHREAL4(tt(np+1, 1))
    tt(np+1, 1) = 0.0
    CALL PUSHREAL4(tt(np+1, 2))
    tt(np+1, 2) = 0.0
    CALL PUSHREAL4(ts(np+1, 1))
    ts(np+1, 1) = 0.0
    CALL PUSHREAL4(ts(np+1, 2))
    ts(np+1, 2) = 0.0
    CALL PUSHREAL4(rr(0, 1))
    rr(0, 1) = 0.0
    CALL PUSHREAL4(rr(0, 2))
    rr(0, 2) = 0.0
    CALL PUSHREAL4(rs(0, 1))
    rs(0, 1) = 0.0
    CALL PUSHREAL4(rs(0, 2))
    rs(0, 2) = 0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
    CALL PUSHREAL4(tt(0, 1))
    tt(0, 1) = 1.0
    CALL PUSHREAL4(tt(0, 2))
    tt(0, 2) = 1.0
    CALL PUSHREAL4(ts(0, 1))
    ts(0, 1) = 1.0
    CALL PUSHREAL4(ts(0, 2))
    ts(0, 2) = 1.0
    cc1 = 0.0
    CALL PUSHREAL4(cc2)
    cc2 = 0.0
    CALL PUSHREAL4(cc3)
    cc3 = 0.0
!-----integration over spectral bands
!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10)
!     The indices 1, 2, 3 are for ice, water, rain particles,
!     respectively.
    DO ib=1,nband_ir
      CALL PUSHINTEGER4(iv)
      iv = ib + 5
!-----options for scaling cloud optical thickness
!if ( OVERCAST == .true. ) then
!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.
!            call getnirtau1(ib,np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,0,0,&
!                           taubeam,taudiff,asycl,ssacl,                               &
!                            aig_uv, awg_uv, arg_uv,                                    &
!                            aib_uv, awb_uv, arb_uv,                                    &
!                            aib_nir, awb_nir, arb_nir,                                 &
!                            aia_nir, awa_nir, ara_nir,                                 &
!                            aig_nir, awg_nir, arg_nir,                                 &
!                            caib, caif,                                                &
!                            CONS_GRAV                                                  )
!else
!-----Compute scaled cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.
      CALL PUSHREAL4ARRAY(ssacl, np)
      CALL PUSHREAL4ARRAY(asycl, np)
      CALL GETNIRTAU1(ib, np, cosz_dev(i), dp_pa, fcld_col, reff_col, &
&               cwc_col, ict, icb, taubeam, taudiff, asycl, ssacl, &
&               aig_uv, awg_uv, arg_uv, aib_uv, awb_uv, arb_uv, aib_nir&
&               , awb_nir, arb_nir, aia_nir, awa_nir, ara_nir, aig_nir, &
&               awg_nir, arg_nir, caib, caif, cons_grav)
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
!MAT--DO NOT FUSE THIS LOOP
!MAT  Loop must run to completion so that cc[1,2,3] are correct.
      DO k=1,np
        IF (k .LT. ict) THEN
          IF (cc1 .LT. fcld_dev(i, k)) THEN
            cc1 = fcld_dev(i, k)
            CALL PUSHCONTROL3B(4)
          ELSE
            CALL PUSHCONTROL3B(5)
            cc1 = cc1
          END IF
        ELSE IF (k .LT. icb) THEN
          IF (cc2 .LT. fcld_dev(i, k)) THEN
            CALL PUSHREAL4(cc2)
            cc2 = fcld_dev(i, k)
            CALL PUSHCONTROL3B(2)
          ELSE
            CALL PUSHREAL4(cc2)
            cc2 = cc2
            CALL PUSHCONTROL3B(3)
          END IF
        ELSE IF (cc3 .LT. fcld_dev(i, k)) THEN
          CALL PUSHREAL4(cc3)
          cc3 = fcld_dev(i, k)
          CALL PUSHCONTROL3B(0)
        ELSE
          CALL PUSHREAL4(cc3)
          cc3 = cc3
          CALL PUSHCONTROL3B(1)
        END IF
      END DO
!MAT--DO NOT FUSE THIS LOOP
!endif !overcast
      DO k=1,np
        CALL PUSHREAL4(tauclb(k))
        tauclb(k) = taubeam(k, 1) + taubeam(k, 2) + taubeam(k, 3) + &
&         taubeam(k, 4)
        CALL PUSHREAL4(tauclf(k))
        tauclf(k) = taudiff(k, 1) + taudiff(k, 2) + taudiff(k, 3) + &
&         taudiff(k, 4)
      END DO
!-----integration over the k-distribution function
      DO ik=1,nk_ir
!-----compute direct beam transmittances of the layer above pl(1)
        CALL PUSHREAL4(td(0, 1))
        td(0, 1) = EXP(-(wvtoa*xk_ir(ik)/cosz_dev(i)))
        CALL PUSHREAL4(td(0, 2))
        td(0, 2) = td(0, 1)
        DO k=1,np
          taurs = ry_ir(ib)*dp(k)
          tauwv = xk_ir(ik)*wh(k)
!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor. Eqs.(6.2)-(6.4)
          tausto = taurs + tauwv + taua_dev(i, k, iv) + 1.0e-7
          ssatau = ssaa_dev(i, k, iv) + taurs + 1.0e-8
          asysto = asya_dev(i, k, iv)
          CALL PUSHREAL4(tautob)
          tautob = tausto
          asytob = asysto/ssatau
          CALL PUSHREAL4(ssatob)
          ssatob = ssatau/tautob + 1.0e-8
          IF (ssatob .GT. 0.999999) THEN
            ssatob = 0.999999
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
            ssatob = ssatob
          END IF
!-----Compute reflectance and transmittance of the clear portion 
!     of a layer
!-----for direct incident radiation
          CALL DELEDD(tautob, ssatob, asytob, cosz_dev(i), rrt, ttt, tdt&
&              )
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
          CALL DELEDD(tautob, ssatob, asytob, dsm, rst, tst, dum)
          CALL PUSHREAL4(rr(k, 1))
          rr(k, 1) = rrt
          CALL PUSHREAL4(tt(k, 1))
          tt(k, 1) = ttt
          CALL PUSHREAL4(td(k, 1))
          td(k, 1) = tdt
          CALL PUSHREAL4(rs(k, 1))
          rs(k, 1) = rst
          CALL PUSHREAL4(ts(k, 1))
          ts(k, 1) = tst
!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer
!-----for direct incident radiation. Eqs.(6.2)-(6.4)
          CALL PUSHREAL4(tautob)
          tautob = tausto + tauclb(k)
          CALL PUSHREAL4(ssatob)
          ssatob = (ssatau+ssacl(k)*tauclb(k))/tautob + 1.0e-8
          IF (ssatob .GT. 0.999999) THEN
            ssatob = 0.999999
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
            ssatob = ssatob
          END IF
          asytob = (asysto+asycl(k)*ssacl(k)*tauclb(k))/(ssatob*tautob)
!-----for diffuse incident radiation
          CALL PUSHREAL4(tautof)
          tautof = tausto + tauclf(k)
          CALL PUSHREAL4(ssatof)
          ssatof = (ssatau+ssacl(k)*tauclf(k))/tautof + 1.0e-8
          IF (ssatof .GT. 0.999999) THEN
            ssatof = 0.999999
            CALL PUSHCONTROL1B(0)
          ELSE
            ssatof = ssatof
            CALL PUSHCONTROL1B(1)
          END IF
          asytof = (asysto+asycl(k)*ssacl(k)*tauclf(k))/(ssatof*tautof)
!-----for direct incident radiation
          CALL DELEDD(tautob, ssatob, asytob, cosz_dev(i), rrt, ttt, tdt&
&              )
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs.(6.5) and (6.6)
          CALL DELEDD(tautof, ssatof, asytof, dsm, rst, tst, dum)
          CALL PUSHREAL4(rr(k, 2))
          rr(k, 2) = rrt
          CALL PUSHREAL4(tt(k, 2))
          tt(k, 2) = ttt
          CALL PUSHREAL4(td(k, 2))
          td(k, 2) = tdt
          CALL PUSHREAL4(rs(k, 2))
          rs(k, 2) = rst
          CALL PUSHREAL4(ts(k, 2))
          ts(k, 2) = tst
        END DO
!-----FLUX CALCULATIONS
!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)
        DO k=1,np+1
          fall(k) = 0.0
        END DO
!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is either 0 or 1.
!if ( OVERCAST == .true. ) then
!-----Inline CLDFLXY
!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)
!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated by 
!         beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)
!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux
!-----ih=1 for clear sky; ih=2 for cloudy sky.
!-----First set is ih = 1
!               rra(np+1)=rr(np+1,1)
!               rxa(np+1)=rs(np+1,1)
!
!               do k=np,0,-1
!                  denm=ts(k,1)/(1.-rs(k,1)*rxa(k+1))
!                  rra(k)=rr(k,1)+(td(k,1)*rra(k+1)+(tt(k,1)-td(k,1))*rxa(k+1))*denm
!                  rxa(k)=rs(k,1)+ts(k,1)*rxa(k+1)*denm
!               end do
!
!               do k=1,np+1
!                  if (k <= np) then
!                     if (k == 1) then
!                        tdaold = td(0,1)
!                        ttaold = tt(0,1)
!                        rsaold = rs(0,1)
!
!                        tdanew = 0.0
!                        ttanew = 0.0
!                        rsanew = 0.0
!                     end if
!                     denm=ts(k,1)/(1.-rsaold*rs(k,1))
!                     tdanew=tdaold*td(k,1)
!                     ttanew=tdaold*tt(k,1)+(tdaold*rsaold*rr(k,1)+ttaold-tdaold)*denm
!                     rsanew=rs(k,1)+ts(k,1)*rsaold*denm
!                  end if
!
!                  denm=1./(1.-rsaold*rxa(k))
!                  fdndir=tdaold
!                  xx4=tdaold*rra(k)
!                  yy=ttaold-tdaold
!                  fdndif=(xx4*rsaold+yy)*denm
!                  fupdif=(xx4+yy*rxa(k))*denm
!                  flxdn=fdndir+fdndif-fupdif
!
!                  fupc(k)=fupdif
!                  fclr(k)=flxdn
!
!                  tdaold = tdanew
!                  ttaold = ttanew
!                  rsaold = rsanew
!
!                  tdanew = 0.0
!                  ttanew = 0.0
!                  rsanew = 0.0
!               end do
!
!!-----Second set is ih = 2
!
!               rra(np+1)=rr(np+1,2)
!               rxa(np+1)=rs(np+1,2)
!
!               do k=np,0,-1
!                  denm=ts(k,2)/(1.-rs(k,2)*rxa(k+1))
!                  rra(k)=rr(k,2)+(td(k,2)*rra(k+1)+(tt(k,2)-td(k,2))*rxa(k+1))*denm
!                  rxa(k)=rs(k,2)+ts(k,2)*rxa(k+1)*denm
!               end do
!
!               do k=1,np+1
!                  if (k <= np) then
!                     if (k == 1) then
!                        tdaold = td(0,2)
!                        ttaold = tt(0,2)
!                        rsaold = rs(0,2)
!
!                        tdanew = 0.0
!                        ttanew = 0.0
!                        rsanew = 0.0
!                     end if
!                     denm=ts(k,2)/(1.-rsaold*rs(k,2))
!                     tdanew=tdaold*td(k,2)
!                     ttanew=tdaold*tt(k,2)+(tdaold*rsaold*rr(k,2)+ttaold-tdaold)*denm
!                     rsanew=rs(k,2)+ts(k,2)*rsaold*denm
!                  end if
!
!                  denm=1./(1.-rsaold*rxa(k))
!                  fdndir=tdaold
!                  xx4=tdaold*rra(k)
!                  yy=ttaold-tdaold
!                  fdndif=(xx4*rsaold+yy)*denm
!                  fupdif=(xx4+yy*rxa(k))*denm
!                  flxdn=fdndir+fdndif-fupdif
!
!                  fupa(k)=fupdif
!                  fall(k)=flxdn
!
!                  tdaold = tdanew
!                  ttaold = ttanew
!                  rsaold = rsanew
!
!                  tdanew = 0.0
!                  ttanew = 0.0
!                  rsanew = 0.0
!               end do
!
!               fsdir=fdndir
!               fsdif=fdndif
!
!!-----End CLDFLXY inline
!
!else
!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is allowed to be between 0 and 1.
!     the all-sky flux, fall is the summation inside the brackets
!     of Eq. (7.11)
!-----Inline CLDFLX
!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated 
!         by beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)
!-----To save memory space, tda, tta, and rsa are pre-computed 
!     for k<icb. The dimension of these parameters is (m,np,2,2). 
!     It would have been (m,np,2,2,2) if these parameters were 
!     computed for all k's.
!-----for high clouds
!     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition
        DO ih=1,2
          CALL PUSHREAL4(tda(0, ih, 1))
          tda(0, ih, 1) = td(0, ih)
          CALL PUSHREAL4(tta(0, ih, 1))
          tta(0, ih, 1) = tt(0, ih)
          CALL PUSHREAL4(rsa(0, ih, 1))
          rsa(0, ih, 1) = rs(0, ih)
          CALL PUSHREAL4(tda(0, ih, 2))
          tda(0, ih, 2) = td(0, ih)
          CALL PUSHREAL4(tta(0, ih, 2))
          tta(0, ih, 2) = tt(0, ih)
          CALL PUSHREAL4(rsa(0, ih, 2))
          rsa(0, ih, 2) = rs(0, ih)
          DO k=1,ict-1
            CALL PUSHREAL4(denm)
            denm = ts(k, ih)/(1.-rsa(k-1, ih, 1)*rs(k, ih))
            CALL PUSHREAL4(tda(k, ih, 1))
            tda(k, ih, 1) = tda(k-1, ih, 1)*td(k, ih)
            CALL PUSHREAL4(tta(k, ih, 1))
            tta(k, ih, 1) = tda(k-1, ih, 1)*tt(k, ih) + (tda(k-1, ih, 1)&
&             *rsa(k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)-tda(k-1, ih, 1)&
&             )*denm
            CALL PUSHREAL4(rsa(k, ih, 1))
            rsa(k, ih, 1) = rs(k, ih) + ts(k, ih)*rsa(k-1, ih, 1)*denm
            CALL PUSHREAL4(tda(k, ih, 2))
            tda(k, ih, 2) = tda(k, ih, 1)
            CALL PUSHREAL4(tta(k, ih, 2))
            tta(k, ih, 2) = tta(k, ih, 1)
            CALL PUSHREAL4(rsa(k, ih, 2))
            rsa(k, ih, 2) = rsa(k, ih, 1)
          END DO
! k loop
!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition
          DO k=ict,icb-1
            DO im=1,2
              CALL PUSHREAL4(denm)
              denm = ts(k, im)/(1.-rsa(k-1, ih, im)*rs(k, im))
              CALL PUSHREAL4(tda(k, ih, im))
              tda(k, ih, im) = tda(k-1, ih, im)*td(k, im)
              CALL PUSHREAL4(tta(k, ih, im))
              tta(k, ih, im) = tda(k-1, ih, im)*tt(k, im) + (tda(k-1, ih&
&               , im)*rsa(k-1, ih, im)*rr(k, im)+tta(k-1, ih, im)-tda(k-&
&               1, ih, im))*denm
              CALL PUSHREAL4(rsa(k, ih, im))
              rsa(k, ih, im) = rs(k, im) + ts(k, im)*rsa(k-1, ih, im)*&
&               denm
            END DO
          END DO
        END DO
! im loop
! k loop
! ih loop
!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)
!-----To save memory space, rra and rxa are pre-computed for k>=icb.
!     the dimension of these parameters is (m,np,2,2). It would have
!     been (m,np,2,2,2) if these parameters were computed for all k's.
!-----for the low clouds
!     is=1 for clear-sky condition, is=2 for cloudy-sky condition
        DO is=1,2
          CALL PUSHREAL4(rra(np+1, 1, is))
          rra(np+1, 1, is) = rr(np+1, is)
          CALL PUSHREAL4(rxa(np+1, 1, is))
          rxa(np+1, 1, is) = rs(np+1, is)
          CALL PUSHREAL4(rra(np+1, 2, is))
          rra(np+1, 2, is) = rr(np+1, is)
          CALL PUSHREAL4(rxa(np+1, 2, is))
          rxa(np+1, 2, is) = rs(np+1, is)
          DO k=np,icb,-1
            CALL PUSHREAL4(denm)
            denm = ts(k, is)/(1.-rs(k, is)*rxa(k+1, 1, is))
            CALL PUSHREAL4(rra(k, 1, is))
            rra(k, 1, is) = rr(k, is) + (td(k, is)*rra(k+1, 1, is)+(tt(k&
&             , is)-td(k, is))*rxa(k+1, 1, is))*denm
            CALL PUSHREAL4(rxa(k, 1, is))
            rxa(k, 1, is) = rs(k, is) + ts(k, is)*rxa(k+1, 1, is)*denm
            CALL PUSHREAL4(rra(k, 2, is))
            rra(k, 2, is) = rra(k, 1, is)
            CALL PUSHREAL4(rxa(k, 2, is))
            rxa(k, 2, is) = rxa(k, 1, is)
          END DO
! k loop
!-----for middle clouds
          DO k=icb-1,ict,-1
            DO im=1,2
              CALL PUSHREAL4(denm)
              denm = ts(k, im)/(1.-rs(k, im)*rxa(k+1, im, is))
              CALL PUSHREAL4(rra(k, im, is))
              rra(k, im, is) = rr(k, im) + (td(k, im)*rra(k+1, im, is)+(&
&               tt(k, im)-td(k, im))*rxa(k+1, im, is))*denm
              CALL PUSHREAL4(rxa(k, im, is))
              rxa(k, im, is) = rs(k, im) + ts(k, im)*rxa(k+1, im, is)*&
&               denm
            END DO
          END DO
        END DO
! im loop
! k loop
! is loop
!-----integration over eight sky situations.
!     ih, im, is denote high, middle and low cloud groups.
        DO ih=1,2
!-----clear portion 
          IF (ih .EQ. 1) THEN
            CALL PUSHREAL4(ch)
            ch = 1.0 - cc1
!-----cloudy portion
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHREAL4(ch)
            ch = cc1
            CALL PUSHCONTROL1B(0)
          END IF
          DO im=1,2
!-----clear portion 
            IF (im .EQ. 1) THEN
              CALL PUSHREAL4(cm)
              cm = ch*(1.0-cc2)
!-----cloudy portion
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHREAL4(cm)
              cm = ch*cc2
              CALL PUSHCONTROL1B(0)
            END IF
            DO is=1,2
!-----clear portion 
              IF (is .EQ. 1) THEN
                CALL PUSHREAL4(ct)
                ct = cm*(1.0-cc3)
!-----cloudy portion
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHREAL4(ct)
                ct = cm*cc3
                CALL PUSHCONTROL1B(0)
              END IF
!-----add one layer at a time, going down.
              DO k=icb,np
                CALL PUSHREAL4(denm)
                denm = ts(k, is)/(1.-rsa(k-1, ih, im)*rs(k, is))
                CALL PUSHREAL4(tda(k, ih, im))
                tda(k, ih, im) = tda(k-1, ih, im)*td(k, is)
                CALL PUSHREAL4(tta(k, ih, im))
                tta(k, ih, im) = tda(k-1, ih, im)*tt(k, is) + (tda(k-1, &
&                 ih, im)*rr(k, is)*rsa(k-1, ih, im)+tta(k-1, ih, im)-&
&                 tda(k-1, ih, im))*denm
                CALL PUSHREAL4(rsa(k, ih, im))
                rsa(k, ih, im) = rs(k, is) + ts(k, is)*rsa(k-1, ih, im)*&
&                 denm
              END DO
! k loop
!-----add one layer at a time, going up.
              DO k=ict-1,0,-1
                CALL PUSHREAL4(denm)
                denm = ts(k, ih)/(1.-rs(k, ih)*rxa(k+1, im, is))
                CALL PUSHREAL4(rra(k, im, is))
                rra(k, im, is) = rr(k, ih) + (td(k, ih)*rra(k+1, im, is)&
&                 +(tt(k, ih)-td(k, ih))*rxa(k+1, im, is))*denm
                CALL PUSHREAL4(rxa(k, im, is))
                rxa(k, im, is) = rs(k, ih) + ts(k, ih)*rxa(k+1, im, is)*&
&                 denm
              END DO
! k loop
!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux
              DO k=1,np+1
                CALL PUSHREAL4(denm)
                denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
                fdndir = tda(k-1, ih, im)
                xx4 = tda(k-1, ih, im)*rra(k, im, is)
                yy = tta(k-1, ih, im) - tda(k-1, ih, im)
                fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
                fupdif = (xx4+yy*rxa(k, im, is))*denm
                flxdn = fdndir + fdndif - fupdif
!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)
                fall(k) = fall(k) + flxdn*ct
              END DO
            END DO
          END DO
        END DO
! is loop
! im loop
! ih loop
!-----End CLDFLX inline
!endif !overcast
!-----flux integration following Eq. (6.1)
        DO k=1,np+1
          flx_dev(i, k) = flx_dev(i, k) + fall(k)*hk_ir(ib, ik)
        END DO
      END DO
    END DO
! ik loop
!-----compute pressure-scaled o2 amount following Eq. (3.5) with f=1.
!     unit is (cm-atm)stp. 165.22 = (1000/980)*23.14%*(22400/32)
!     compute flux reduction due to oxygen following Eq. (3.18). 0.0633 is the
!     fraction of insolation contained in the oxygen bands
    df(0) = 0.0
    cnt = 165.22*snt
    so2(1) = scal0*cnt
! LLT increased parameter 145 to 155 to enhance effect
    df(1) = 0.0633*(1.-EXP(-(0.000155*SQRT(so2(1)))))
    DO k=1,np
      so2(k+1) = so2(k) + scal(k)*cnt
! LLT increased parameter 145 to 155 to enhance effect
      df(k+1) = 0.0633*(1.0-EXP(-(0.000155*SQRT(so2(k+1)))))
    END DO
!-----for solar heating due to co2 scaling follows Eq(3.5) with f=1.
!     unit is (cm-atm)stp. 789 = (1000/980)*(44/28.97)*(22400/44)
    so2(1) = 789.*co2*scal0
    DO k=1,np
      so2(k+1) = so2(k) + 789.*co2*scal(k)
    END DO
!-----The updated flux reduction for co2 absorption in Band 7 where absorption due to
!     water vapor and co2 are both moderate. df is given by the second term on the
!     right-hand-side of Eq. (3.24) divided by So. so2 and swh are the co2 and
!     water vapor amounts integrated from the top of the atmosphere
    u1 = -3.0
    du = 0.15
    w1 = -4.0
    dw = 0.15
!-----Inline RFLX
    x0 = u1 + REAL(nu)*du
    y0 = w1 + REAL(nw)*dw
    x1 = u1 - 0.5*du
    y1 = w1 - 0.5*dw
    DO k=1,np+1
      x3 = LOG10(so2(k)*snt)
      IF (x3 .GT. x0) THEN
        ulog = x0
      ELSE
        ulog = x3
      END IF
      x4 = LOG10(swh(k)*snt)
      IF (x4 .GT. y0) THEN
        wlog = y0
        CALL PUSHCONTROL1B(0)
      ELSE
        wlog = x4
        CALL PUSHCONTROL1B(1)
      END IF
      CALL PUSHINTEGER4(ic)
      ic = INT((ulog-x1)/du + 1.)
      CALL PUSHINTEGER4(iw)
      iw = INT((wlog-y1)/dw + 1.)
      IF (ic .LT. 2) ic = 2
      IF (iw .LT. 2) iw = 2
      IF (ic .GT. nu) ic = nu
      IF (iw .GT. nw) iw = nw
      dc = ulog - REAL(ic-2)*du - u1
      dd = wlog - REAL(iw-2)*dw - w1
      x2 = cah(ic-1, iw-1) + (cah(ic-1, iw)-cah(ic-1, iw-1))/dw*dd
      y2 = x2 + (cah(ic, iw-1)-cah(ic-1, iw-1))/du*dc
      IF (y2 .LT. 0.0) THEN
        y2 = 0.0
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        y2 = y2
      END IF
! LLT increase CO2 effect to help reduce cold tropopause bias
      df(k) = df(k) + 1.5*y2
    END DO
!-----df is the updated flux reduction for co2 absorption
!     in Band 8 where the co2 absorption has a large impact
!     on the heating of middle atmosphere. From the table
!     given by Eq. (3.19)
    u1 = 0.000250
    du = 0.000050
    w1 = -2.0
    dw = 0.05
!-----Inline RFLX
    x0 = u1 + REAL(nx)*du
    y0 = w1 + REAL(ny)*dw
    x1 = u1 - 0.5*du
    y1 = w1 - 0.5*dw
    DO k=1,np+1
      IF (co2*snt .GT. x0) THEN
        ulog = x0
      ELSE
        ulog = co2*snt
      END IF
      x5 = LOG10(pl_dev(i, k))
      IF (x5 .GT. y0) THEN
        wlog = y0
      ELSE
        wlog = x5
      END IF
      CALL PUSHINTEGER4(ic)
      ic = INT((ulog-x1)/du + 1.)
      CALL PUSHINTEGER4(iw)
      iw = INT((wlog-y1)/dw + 1.)
      IF (ic .LT. 2) ic = 2
      IF (iw .LT. 2) iw = 2
      IF (ic .GT. nx) ic = nx
      IF (iw .GT. ny) iw = ny
      dc = ulog - REAL(ic-2)*du - u1
      dd = wlog - REAL(iw-2)*dw - w1
      x2 = coa(ic-1, iw-1) + (coa(ic-1, iw)-coa(ic-1, iw-1))/dw*dd
      y2 = x2 + (coa(ic, iw-1)-coa(ic-1, iw-1))/du*dc
      IF (y2 .LT. 0.0) THEN
        y2 = 0.0
      ELSE
        y2 = y2
      END IF
! LLT increase CO2 effect to help reduce cold tropopause bias
      df(k) = df(k) + 1.5*y2
    END DO
!-----adjust the o2-co2 reduction below cloud top following Eq. (6.18)
    foundtop = 0
    DO k=1,np
      IF (fcld_dev(i, k) .GT. 0.02 .AND. foundtop .EQ. 0) THEN
        foundtop = 1
        ntop = k
      END IF
    END DO
    IF (foundtop .EQ. 0) ntop = np + 1
    CALL PUSHREAL4(dftop)
    dftop = df(ntop)
    DO k=1,np+1
      IF (k .GT. ntop) THEN
        xx4 = flx_dev(i, k)/flx_dev(i, ntop)
        CALL PUSHREAL4(df(k))
        df(k) = dftop + xx4*(df(k)-dftop)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
!-----update the net fluxes
    DO k=1,np+1
      IF (df(k) .GT. flx_dev(i, k) - 1.0e-8) THEN
        df(k) = flx_dev(i, k) - 1.0e-8
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        df(k) = df(k)
      END IF
    END DO
  END DO run_loop
  ssaa_devb = 0.0
  asya_devb = 0.0
  taua_devb = 0.0
  dfb = 0.0
  swhb = 0.0
  reff_colb = 0.0
  ohb = 0.0
  fallb = 0.0
  asyclb = 0.0
  ssaclb = 0.0
  fcld_colb = 0.0
  cwc_colb = 0.0
  tauclbb = 0.0
  tauclfb = 0.0
  whb = 0.0
  DO i=m,1,-1
    DO k=np+1,1,-1
      dfb(k) = dfb(k) - flx_devb(i, k)
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        flx_devb(i, k) = flx_devb(i, k) + dfb(k)
        dfb(k) = 0.0
      END IF
    END DO
    dftopb = 0.0
    DO k=np+1,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        xx4 = flx_dev(i, k)/flx_dev(i, ntop)
        CALL POPREAL4(df(k))
        dftopb = dftopb + (1.0-xx4)*dfb(k)
        xx4b = (df(k)-dftop)*dfb(k)
        dfb(k) = xx4*dfb(k)
        tempb56 = xx4b/flx_dev(i, ntop)
        flx_devb(i, k) = flx_devb(i, k) + tempb56
        flx_devb(i, ntop) = flx_devb(i, ntop) - flx_dev(i, k)*tempb56/&
&         flx_dev(i, ntop)
      END IF
    END DO
    CALL POPREAL4(dftop)
    dfb(ntop) = dfb(ntop) + dftopb
    DO k=np+1,1,-1
      CALL POPINTEGER4(iw)
      CALL POPINTEGER4(ic)
    END DO
    dw = 0.15
    DO k=np+1,1,-1
      y2b = 1.5*dfb(k)
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) y2b = 0.0
      x2b = y2b
      ddb = (cah(ic-1, iw)-cah(ic-1, iw-1))*x2b/dw
      wlogb = ddb
      CALL POPINTEGER4(iw)
      CALL POPINTEGER4(ic)
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        x4b = 0.0
      ELSE
        x4b = wlogb
      END IF
      swhb(k) = swhb(k) + x4b/(swh(k)*LOG(10.0))
    END DO
    DO k=np,1,-1
      dfb(k+1) = 0.0
    END DO
    dfb(1) = 0.0
    dfb(0) = 0.0
    wvtoa = 1.02*wa_dev(i, 1)*scal0*(1.0+0.00135*(ta_dev(i, 1)-240.)) + &
&     1.0e-9
    tsb = 0.0
    ttb = 0.0
    rsab = 0.0
    rrb = 0.0
    rsb = 0.0
    wvtoab = 0.0
    cc1b = 0.0
    cc2b = 0.0
    cc3b = 0.0
    ttab = 0.0
    rxab = 0.0
    tdab = 0.0
    rrab = 0.0
    tdb = 0.0
    DO ib=nband_ir,1,-1
      DO ik=nk_ir,1,-1
        DO k=np+1,1,-1
          fallb(k) = fallb(k) + hk_ir(ib, ik)*flx_devb(i, k)
        END DO
        DO ih=2,1,-1
          chb = 0.0
          DO im=2,1,-1
            cmb = 0.0
            DO is=2,1,-1
              ctb = 0.0
              DO k=np+1,1,-1
                yy = tta(k-1, ih, im) - tda(k-1, ih, im)
                denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
                xx4 = tda(k-1, ih, im)*rra(k, im, is)
                fupdif = (xx4+yy*rxa(k, im, is))*denm
                fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
                fdndir = tda(k-1, ih, im)
                flxdn = fdndir + fdndif - fupdif
                flxdnb = ct*fallb(k)
                ctb = ctb + flxdn*fallb(k)
                fdndirb = flxdnb
                fdndifb = flxdnb
                fupdifb = -flxdnb
                tempb53 = denm*fupdifb
                denmb = (xx4*rsa(k-1, ih, im)+yy)*fdndifb + (xx4+yy*rxa(&
&                 k, im, is))*fupdifb
                tempb54 = denm*fdndifb
                xx4b = rsa(k-1, ih, im)*tempb54 + tempb53
                yyb = tempb54 + rxa(k, im, is)*tempb53
                ttab(k-1, ih, im) = ttab(k-1, ih, im) + yyb
                tdab(k-1, ih, im) = tdab(k-1, ih, im) + rra(k, im, is)*&
&                 xx4b + fdndirb - yyb
                rrab(k, im, is) = rrab(k, im, is) + tda(k-1, ih, im)*&
&                 xx4b
                CALL POPREAL4(denm)
                temp38 = -(rsa(k-1, ih, im)*rxa(k, im, is)) + 1.
                tempb55 = -(denmb/temp38**2)
                rxab(k, im, is) = rxab(k, im, is) + yy*tempb53 - rsa(k-1&
&                 , ih, im)*tempb55
                rsab(k-1, ih, im) = rsab(k-1, ih, im) + xx4*tempb54 - &
&                 rxa(k, im, is)*tempb55
              END DO
              DO k=0,ict-1,1
                CALL POPREAL4(rxa(k, im, is))
                tempb50 = rxa(k+1, im, is)*rxab(k, im, is)
                CALL POPREAL4(rra(k, im, is))
                tempb52 = denm*rrab(k, im, is)
                temp37 = rxa(k+1, im, is)
                temp36 = tt(k, ih) - td(k, ih)
                rrb(k, ih) = rrb(k, ih) + rrab(k, im, is)
                tdb(k, ih) = tdb(k, ih) + (rra(k+1, im, is)-temp37)*&
&                 tempb52
                rrab(k+1, im, is) = rrab(k+1, im, is) + td(k, ih)*&
&                 tempb52
                denmb = (td(k, ih)*rra(k+1, im, is)+temp36*temp37)*rrab(&
&                 k, im, is) + ts(k, ih)*tempb50
                ttb(k, ih) = ttb(k, ih) + temp37*tempb52
                rrab(k, im, is) = 0.0
                temp35 = -(rs(k, ih)*rxa(k+1, im, is)) + 1.
                tsb(k, ih) = tsb(k, ih) + denmb/temp35 + denm*tempb50
                tempb51 = -(ts(k, ih)*denmb/temp35**2)
                rsb(k, ih) = rsb(k, ih) + rxab(k, im, is) - rxa(k+1, im&
&                 , is)*tempb51
                rxab(k+1, im, is) = rxab(k+1, im, is) + ts(k, ih)*denm*&
&                 rxab(k, im, is)
                rxab(k, im, is) = 0.0
                rxab(k+1, im, is) = rxab(k+1, im, is) + temp36*tempb52 -&
&                 rs(k, ih)*tempb51
                CALL POPREAL4(denm)
              END DO
              DO k=np,icb,-1
                CALL POPREAL4(rsa(k, ih, im))
                tempb47 = rsa(k-1, ih, im)*rsab(k, ih, im)
                CALL POPREAL4(tta(k, ih, im))
                tempb49 = denm*ttab(k, ih, im)
                temp34 = rsa(k-1, ih, im)
                temp33 = tda(k-1, ih, im)*rr(k, is)
                tdab(k-1, ih, im) = tdab(k-1, ih, im) + td(k, is)*tdab(k&
&                 , ih, im) + (temp34*rr(k, is)-1.0)*tempb49 + tt(k, is)&
&                 *ttab(k, ih, im)
                ttb(k, is) = ttb(k, is) + tda(k-1, ih, im)*ttab(k, ih, &
&                 im)
                rrb(k, is) = rrb(k, is) + temp34*tda(k-1, ih, im)*&
&                 tempb49
                ttab(k-1, ih, im) = ttab(k-1, ih, im) + tempb49
                denmb = (temp33*temp34+tta(k-1, ih, im)-tda(k-1, ih, im)&
&                 )*ttab(k, ih, im) + ts(k, is)*tempb47
                ttab(k, ih, im) = 0.0
                CALL POPREAL4(tda(k, ih, im))
                tdb(k, is) = tdb(k, is) + tda(k-1, ih, im)*tdab(k, ih, &
&                 im)
                tdab(k, ih, im) = 0.0
                temp32 = -(rsa(k-1, ih, im)*rs(k, is)) + 1.
                tsb(k, is) = tsb(k, is) + denmb/temp32 + denm*tempb47
                tempb48 = -(ts(k, is)*denmb/temp32**2)
                rsb(k, is) = rsb(k, is) + rsab(k, ih, im) - rsa(k-1, ih&
&                 , im)*tempb48
                rsab(k-1, ih, im) = rsab(k-1, ih, im) + ts(k, is)*denm*&
&                 rsab(k, ih, im)
                rsab(k, ih, im) = 0.0
                rsab(k-1, ih, im) = rsab(k-1, ih, im) + temp33*tempb49 -&
&                 rs(k, is)*tempb48
                CALL POPREAL4(denm)
              END DO
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPREAL4(ct)
                cmb = cmb + cc3*ctb
                cc3b = cc3b + cm*ctb
              ELSE
                CALL POPREAL4(ct)
                cmb = cmb + (1.0-cc3)*ctb
                cc3b = cc3b - cm*ctb
              END IF
            END DO
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL4(cm)
              chb = chb + cc2*cmb
              cc2b = cc2b + ch*cmb
            ELSE
              CALL POPREAL4(cm)
              chb = chb + (1.0-cc2)*cmb
              cc2b = cc2b - ch*cmb
            END IF
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4(ch)
            cc1b = cc1b + chb
          ELSE
            CALL POPREAL4(ch)
            cc1b = cc1b - chb
          END IF
        END DO
        DO is=2,1,-1
          DO k=ict,icb-1,1
            DO im=2,1,-1
              CALL POPREAL4(rxa(k, im, is))
              tempb44 = rxa(k+1, im, is)*rxab(k, im, is)
              CALL POPREAL4(rra(k, im, is))
              tempb46 = denm*rrab(k, im, is)
              temp31 = rxa(k+1, im, is)
              temp30 = tt(k, im) - td(k, im)
              rrb(k, im) = rrb(k, im) + rrab(k, im, is)
              tdb(k, im) = tdb(k, im) + (rra(k+1, im, is)-temp31)*&
&               tempb46
              rrab(k+1, im, is) = rrab(k+1, im, is) + td(k, im)*tempb46
              denmb = (td(k, im)*rra(k+1, im, is)+temp30*temp31)*rrab(k&
&               , im, is) + ts(k, im)*tempb44
              ttb(k, im) = ttb(k, im) + temp31*tempb46
              rrab(k, im, is) = 0.0
              temp29 = -(rs(k, im)*rxa(k+1, im, is)) + 1.
              tsb(k, im) = tsb(k, im) + denmb/temp29 + denm*tempb44
              tempb45 = -(ts(k, im)*denmb/temp29**2)
              rsb(k, im) = rsb(k, im) + rxab(k, im, is) - rxa(k+1, im, &
&               is)*tempb45
              rxab(k+1, im, is) = rxab(k+1, im, is) + ts(k, im)*denm*&
&               rxab(k, im, is)
              rxab(k, im, is) = 0.0
              rxab(k+1, im, is) = rxab(k+1, im, is) + temp30*tempb46 - &
&               rs(k, im)*tempb45
              CALL POPREAL4(denm)
            END DO
          END DO
          DO k=icb,np,1
            CALL POPREAL4(rxa(k, 2, is))
            rxab(k, 1, is) = rxab(k, 1, is) + rxab(k, 2, is)
            rxab(k, 2, is) = 0.0
            CALL POPREAL4(rra(k, 2, is))
            rrab(k, 1, is) = rrab(k, 1, is) + rrab(k, 2, is)
            rrab(k, 2, is) = 0.0
            CALL POPREAL4(rxa(k, 1, is))
            tempb41 = rxa(k+1, 1, is)*rxab(k, 1, is)
            CALL POPREAL4(rra(k, 1, is))
            tempb43 = denm*rrab(k, 1, is)
            temp28 = rxa(k+1, 1, is)
            temp27 = tt(k, is) - td(k, is)
            rrb(k, is) = rrb(k, is) + rrab(k, 1, is)
            tdb(k, is) = tdb(k, is) + (rra(k+1, 1, is)-temp28)*tempb43
            rrab(k+1, 1, is) = rrab(k+1, 1, is) + td(k, is)*tempb43
            denmb = (td(k, is)*rra(k+1, 1, is)+temp27*temp28)*rrab(k, 1&
&             , is) + ts(k, is)*tempb41
            ttb(k, is) = ttb(k, is) + temp28*tempb43
            rrab(k, 1, is) = 0.0
            temp26 = -(rs(k, is)*rxa(k+1, 1, is)) + 1.
            tsb(k, is) = tsb(k, is) + denmb/temp26 + denm*tempb41
            tempb42 = -(ts(k, is)*denmb/temp26**2)
            rsb(k, is) = rsb(k, is) + rxab(k, 1, is) - rxa(k+1, 1, is)*&
&             tempb42
            rxab(k+1, 1, is) = rxab(k+1, 1, is) + ts(k, is)*denm*rxab(k&
&             , 1, is)
            rxab(k, 1, is) = 0.0
            rxab(k+1, 1, is) = rxab(k+1, 1, is) + temp27*tempb43 - rs(k&
&             , is)*tempb42
            CALL POPREAL4(denm)
          END DO
          CALL POPREAL4(rxa(np+1, 2, is))
          rsb(np+1, is) = rsb(np+1, is) + rxab(np+1, 2, is)
          rxab(np+1, 2, is) = 0.0
          CALL POPREAL4(rra(np+1, 2, is))
          rrb(np+1, is) = rrb(np+1, is) + rrab(np+1, 2, is)
          rrab(np+1, 2, is) = 0.0
          CALL POPREAL4(rxa(np+1, 1, is))
          rsb(np+1, is) = rsb(np+1, is) + rxab(np+1, 1, is)
          rxab(np+1, 1, is) = 0.0
          CALL POPREAL4(rra(np+1, 1, is))
          rrb(np+1, is) = rrb(np+1, is) + rrab(np+1, 1, is)
          rrab(np+1, 1, is) = 0.0
        END DO
        DO ih=2,1,-1
          DO k=icb-1,ict,-1
            DO im=2,1,-1
              CALL POPREAL4(rsa(k, ih, im))
              tempb38 = rsa(k-1, ih, im)*rsab(k, ih, im)
              CALL POPREAL4(tta(k, ih, im))
              tempb40 = denm*ttab(k, ih, im)
              temp25 = rsa(k-1, ih, im)
              temp24 = tda(k-1, ih, im)*rr(k, im)
              tdab(k-1, ih, im) = tdab(k-1, ih, im) + td(k, im)*tdab(k, &
&               ih, im) + (temp25*rr(k, im)-1.0)*tempb40 + tt(k, im)*&
&               ttab(k, ih, im)
              ttb(k, im) = ttb(k, im) + tda(k-1, ih, im)*ttab(k, ih, im)
              rrb(k, im) = rrb(k, im) + temp25*tda(k-1, ih, im)*tempb40
              ttab(k-1, ih, im) = ttab(k-1, ih, im) + tempb40
              denmb = (temp24*temp25+tta(k-1, ih, im)-tda(k-1, ih, im))*&
&               ttab(k, ih, im) + ts(k, im)*tempb38
              ttab(k, ih, im) = 0.0
              CALL POPREAL4(tda(k, ih, im))
              tdb(k, im) = tdb(k, im) + tda(k-1, ih, im)*tdab(k, ih, im)
              tdab(k, ih, im) = 0.0
              temp23 = -(rsa(k-1, ih, im)*rs(k, im)) + 1.
              tsb(k, im) = tsb(k, im) + denmb/temp23 + denm*tempb38
              tempb39 = -(ts(k, im)*denmb/temp23**2)
              rsb(k, im) = rsb(k, im) + rsab(k, ih, im) - rsa(k-1, ih, &
&               im)*tempb39
              rsab(k-1, ih, im) = rsab(k-1, ih, im) + ts(k, im)*denm*&
&               rsab(k, ih, im)
              rsab(k, ih, im) = 0.0
              rsab(k-1, ih, im) = rsab(k-1, ih, im) + temp24*tempb40 - &
&               rs(k, im)*tempb39
              CALL POPREAL4(denm)
            END DO
          END DO
          DO k=ict-1,1,-1
            CALL POPREAL4(rsa(k, ih, 2))
            rsab(k, ih, 1) = rsab(k, ih, 1) + rsab(k, ih, 2)
            rsab(k, ih, 2) = 0.0
            CALL POPREAL4(tta(k, ih, 2))
            ttab(k, ih, 1) = ttab(k, ih, 1) + ttab(k, ih, 2)
            ttab(k, ih, 2) = 0.0
            CALL POPREAL4(tda(k, ih, 2))
            tdab(k, ih, 1) = tdab(k, ih, 1) + tdab(k, ih, 2)
            tdab(k, ih, 2) = 0.0
            CALL POPREAL4(rsa(k, ih, 1))
            tempb35 = rsa(k-1, ih, 1)*rsab(k, ih, 1)
            CALL POPREAL4(tta(k, ih, 1))
            tempb37 = denm*ttab(k, ih, 1)
            temp22 = rsa(k-1, ih, 1)
            temp21 = tda(k-1, ih, 1)*rr(k, ih)
            tdab(k-1, ih, 1) = tdab(k-1, ih, 1) + td(k, ih)*tdab(k, ih, &
&             1) + (temp22*rr(k, ih)-1.0)*tempb37 + tt(k, ih)*ttab(k, ih&
&             , 1)
            ttb(k, ih) = ttb(k, ih) + tda(k-1, ih, 1)*ttab(k, ih, 1)
            rrb(k, ih) = rrb(k, ih) + temp22*tda(k-1, ih, 1)*tempb37
            ttab(k-1, ih, 1) = ttab(k-1, ih, 1) + tempb37
            denmb = (temp21*temp22+tta(k-1, ih, 1)-tda(k-1, ih, 1))*ttab&
&             (k, ih, 1) + ts(k, ih)*tempb35
            ttab(k, ih, 1) = 0.0
            CALL POPREAL4(tda(k, ih, 1))
            tdb(k, ih) = tdb(k, ih) + tda(k-1, ih, 1)*tdab(k, ih, 1)
            tdab(k, ih, 1) = 0.0
            temp20 = -(rsa(k-1, ih, 1)*rs(k, ih)) + 1.
            tsb(k, ih) = tsb(k, ih) + denmb/temp20 + denm*tempb35
            tempb36 = -(ts(k, ih)*denmb/temp20**2)
            rsb(k, ih) = rsb(k, ih) + rsab(k, ih, 1) - rsa(k-1, ih, 1)*&
&             tempb36
            rsab(k-1, ih, 1) = rsab(k-1, ih, 1) + ts(k, ih)*denm*rsab(k&
&             , ih, 1)
            rsab(k, ih, 1) = 0.0
            rsab(k-1, ih, 1) = rsab(k-1, ih, 1) + temp21*tempb37 - rs(k&
&             , ih)*tempb36
            CALL POPREAL4(denm)
          END DO
          CALL POPREAL4(rsa(0, ih, 2))
          rsb(0, ih) = rsb(0, ih) + rsab(0, ih, 2)
          rsab(0, ih, 2) = 0.0
          CALL POPREAL4(tta(0, ih, 2))
          ttb(0, ih) = ttb(0, ih) + ttab(0, ih, 2)
          ttab(0, ih, 2) = 0.0
          CALL POPREAL4(tda(0, ih, 2))
          tdb(0, ih) = tdb(0, ih) + tdab(0, ih, 2)
          tdab(0, ih, 2) = 0.0
          CALL POPREAL4(rsa(0, ih, 1))
          rsb(0, ih) = rsb(0, ih) + rsab(0, ih, 1)
          rsab(0, ih, 1) = 0.0
          CALL POPREAL4(tta(0, ih, 1))
          ttb(0, ih) = ttb(0, ih) + ttab(0, ih, 1)
          ttab(0, ih, 1) = 0.0
          CALL POPREAL4(tda(0, ih, 1))
          tdb(0, ih) = tdb(0, ih) + tdab(0, ih, 1)
          tdab(0, ih, 1) = 0.0
        END DO
        DO k=np+1,1,-1
          fallb(k) = 0.0
        END DO
        DO k=np,1,-1
          CALL POPREAL4(ts(k, 2))
          tstb = tsb(k, 2)
          tsb(k, 2) = 0.0
          CALL POPREAL4(rs(k, 2))
          rstb = rsb(k, 2)
          rsb(k, 2) = 0.0
          CALL POPREAL4(td(k, 2))
          tdtb = tdb(k, 2)
          tdb(k, 2) = 0.0
          CALL POPREAL4(tt(k, 2))
          tttb = ttb(k, 2)
          ttb(k, 2) = 0.0
          CALL POPREAL4(rr(k, 2))
          rrtb = rrb(k, 2)
          rrb(k, 2) = 0.0
          asysto = asya_dev(i, k, iv)
          asytof = (asysto+asycl(k)*ssacl(k)*tauclf(k))/(ssatof*tautof)
          tautofb = 0.0
          ssatofb = 0.0
          asytofb = 0.0
          dumb = 0.0
          CALL DELEDD_B(tautof, tautofb, ssatof, ssatofb, asytof, &
&                 asytofb, dsm, rst, rstb, tst, tstb, dum, dumb)
          asytob = (asysto+asycl(k)*ssacl(k)*tauclb(k))/(ssatob*tautob)
          tautobb = 0.0
          ssatobb = 0.0
          asytobb = 0.0
          CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, &
&                 asytobb, cosz_dev(i), rrt, rrtb, ttt, tttb, tdt, tdtb)
          tempb33 = asytofb/(ssatof*tautof)
          temp19 = asycl(k)*ssacl(k)
          tempb34 = -((asysto+temp19*tauclf(k))*tempb33/(ssatof*tautof))
          asystob = tempb33
          asyclb(k) = asyclb(k) + tauclf(k)*ssacl(k)*tempb33
          ssaclb(k) = ssaclb(k) + tauclf(k)*asycl(k)*tempb33
          tauclfb(k) = tauclfb(k) + temp19*tempb33
          ssatofb = ssatofb + tautof*tempb34
          tautofb = tautofb + ssatof*tempb34
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            taurs = ry_ir(ib)*dp(k)
            ssatau = ssaa_dev(i, k, iv) + taurs + 1.0e-8
            ssatofb = 0.0
          ELSE
            taurs = ry_ir(ib)*dp(k)
            ssatau = ssaa_dev(i, k, iv) + taurs + 1.0e-8
          END IF
          tempb31 = asytobb/(ssatob*tautob)
          CALL POPREAL4(ssatof)
          tempb30 = ssatofb/tautof
          ssataub = tempb30
          ssaclb(k) = ssaclb(k) + tauclb(k)*asycl(k)*tempb31 + tauclf(k)&
&           *tempb30
          tautofb = tautofb - (ssatau+ssacl(k)*tauclf(k))*tempb30/tautof
          tauclfb(k) = tauclfb(k) + tautofb + ssacl(k)*tempb30
          CALL POPREAL4(tautof)
          taustob = tautofb
          temp18 = asycl(k)*ssacl(k)
          tempb32 = -((asysto+temp18*tauclb(k))*tempb31/(ssatob*tautob))
          asystob = asystob + tempb31
          asyclb(k) = asyclb(k) + tauclb(k)*ssacl(k)*tempb31
          tauclbb(k) = tauclbb(k) + temp18*tempb31
          ssatobb = ssatobb + tautob*tempb32
          tautobb = tautobb + ssatob*tempb32
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) ssatobb = 0.0
          CALL POPREAL4(ssatob)
          tempb29 = ssatobb/tautob
          ssataub = ssataub + tempb29
          ssaclb(k) = ssaclb(k) + tauclb(k)*tempb29
          tautobb = tautobb - (ssatau+ssacl(k)*tauclb(k))*tempb29/tautob
          tauclbb(k) = tauclbb(k) + tautobb + ssacl(k)*tempb29
          CALL POPREAL4(tautob)
          taustob = taustob + tautobb
          CALL POPREAL4(ts(k, 1))
          tstb = tsb(k, 1)
          tsb(k, 1) = 0.0
          CALL POPREAL4(rs(k, 1))
          rstb = rsb(k, 1)
          rsb(k, 1) = 0.0
          CALL POPREAL4(td(k, 1))
          tdtb = tdb(k, 1)
          tdb(k, 1) = 0.0
          CALL POPREAL4(tt(k, 1))
          tttb = ttb(k, 1)
          ttb(k, 1) = 0.0
          CALL POPREAL4(rr(k, 1))
          rrtb = rrb(k, 1)
          rrb(k, 1) = 0.0
          asytob = asysto/ssatau
          tautobb = 0.0
          ssatobb = 0.0
          asytobb = 0.0
          dumb = 0.0
          CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, &
&                 asytobb, dsm, rst, rstb, tst, tstb, dum, dumb)
          CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, &
&                 asytobb, cosz_dev(i), rrt, rrtb, ttt, tttb, tdt, tdtb)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) ssatobb = 0.0
          CALL POPREAL4(ssatob)
          ssataub = ssataub + ssatobb/tautob - asysto*asytobb/ssatau**2
          tautobb = tautobb - ssatau*ssatobb/tautob**2
          asystob = asystob + asytobb/ssatau
          CALL POPREAL4(tautob)
          taustob = taustob + tautobb
          asya_devb(i, k, iv) = asya_devb(i, k, iv) + asystob
          ssaa_devb(i, k, iv) = ssaa_devb(i, k, iv) + ssataub
          tauwvb = taustob
          taua_devb(i, k, iv) = taua_devb(i, k, iv) + taustob
          whb(k) = whb(k) + xk_ir(ik)*tauwvb
        END DO
        CALL POPREAL4(td(0, 2))
        tdb(0, 1) = tdb(0, 1) + tdb(0, 2)
        tdb(0, 2) = 0.0
        CALL POPREAL4(td(0, 1))
        wvtoab = wvtoab - xk_ir(ik)*EXP(-(xk_ir(ik)*(wvtoa/cosz_dev(i)))&
&         )*tdb(0, 1)/cosz_dev(i)
        tdb(0, 1) = 0.0
      END DO
      taudiffb = 0.0
      taubeamb = 0.0
      DO k=np,1,-1
        CALL POPREAL4(tauclf(k))
        taudiffb(k, 1) = taudiffb(k, 1) + tauclfb(k)
        taudiffb(k, 2) = taudiffb(k, 2) + tauclfb(k)
        taudiffb(k, 3) = taudiffb(k, 3) + tauclfb(k)
        taudiffb(k, 4) = taudiffb(k, 4) + tauclfb(k)
        tauclfb(k) = 0.0
        CALL POPREAL4(tauclb(k))
        taubeamb(k, 1) = taubeamb(k, 1) + tauclbb(k)
        taubeamb(k, 2) = taubeamb(k, 2) + tauclbb(k)
        taubeamb(k, 3) = taubeamb(k, 3) + tauclbb(k)
        taubeamb(k, 4) = taubeamb(k, 4) + tauclbb(k)
        tauclbb(k) = 0.0
      END DO
      DO k=np,1,-1
        CALL POPCONTROL3B(branch)
        IF (branch .LT. 3) THEN
          IF (branch .EQ. 0) THEN
            CALL POPREAL4(cc3)
            fcld_devb(i, k) = fcld_devb(i, k) + cc3b
            cc3b = 0.0
          ELSE IF (branch .EQ. 1) THEN
            CALL POPREAL4(cc3)
          ELSE
            CALL POPREAL4(cc2)
            fcld_devb(i, k) = fcld_devb(i, k) + cc2b
            cc2b = 0.0
          END IF
        ELSE IF (branch .EQ. 3) THEN
          CALL POPREAL4(cc2)
        ELSE IF (branch .EQ. 4) THEN
          fcld_devb(i, k) = fcld_devb(i, k) + cc1b
          cc1b = 0.0
        END IF
      END DO
      CALL POPREAL4ARRAY(asycl, np)
      CALL POPREAL4ARRAY(ssacl, np)
      CALL GETNIRTAU1_B(ib, np, cosz_dev(i), dp_pa, fcld_col, fcld_colb&
&                 , reff_col, reff_colb, cwc_col, cwc_colb, ict, icb, &
&                 taubeam, taubeamb, taudiff, taudiffb, asycl, asyclb, &
&                 ssacl, ssaclb, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv&
&                 , arb_uv, aib_nir, awb_nir, arb_nir, aia_nir, awa_nir&
&                 , ara_nir, aig_nir, awg_nir, arg_nir, caib, caif, &
&                 cons_grav)
      CALL POPINTEGER4(iv)
    END DO
    CALL POPREAL4(cc3)
    CALL POPREAL4(cc2)
    CALL POPREAL4(ts(0, 2))
    tsb(0, 2) = 0.0
    CALL POPREAL4(ts(0, 1))
    tsb(0, 1) = 0.0
    CALL POPREAL4(tt(0, 2))
    ttb(0, 2) = 0.0
    CALL POPREAL4(tt(0, 1))
    ttb(0, 1) = 0.0
    CALL POPREAL4(rs(0, 2))
    rsb(0, 2) = 0.0
    CALL POPREAL4(rs(0, 1))
    rsb(0, 1) = 0.0
    CALL POPREAL4(rr(0, 2))
    rrb(0, 2) = 0.0
    CALL POPREAL4(rr(0, 1))
    rrb(0, 1) = 0.0
    CALL POPREAL4(ts(np+1, 2))
    tsb(np+1, 2) = 0.0
    CALL POPREAL4(ts(np+1, 1))
    tsb(np+1, 1) = 0.0
    CALL POPREAL4(tt(np+1, 2))
    ttb(np+1, 2) = 0.0
    CALL POPREAL4(tt(np+1, 1))
    ttb(np+1, 1) = 0.0
    CALL POPREAL4(td(np+1, 2))
    tdb(np+1, 2) = 0.0
    CALL POPREAL4(td(np+1, 1))
    tdb(np+1, 1) = 0.0
    CALL POPREAL4(rs(np+1, 2))
    rsb(np+1, 2) = 0.0
    CALL POPREAL4(rs(np+1, 1))
    rsb(np+1, 1) = 0.0
    CALL POPREAL4(rr(np+1, 2))
    rrb(np+1, 2) = 0.0
    CALL POPREAL4(rr(np+1, 1))
    rrb(np+1, 1) = 0.0
    o3toa = 1.02*oa_dev(i, 1)*xtoa*466.7 + 1.0e-8
    cc1b = 0.0
    cc2b = 0.0
    cc3b = 0.0
    o3toab = 0.0
    DO ib=nband_uv,1,-1
      DO k=np+1,1,-1
        fallb(k) = fallb(k) + hk_uv(ib)*flx_devb(i, k)
      END DO
      DO ih=2,1,-1
        chb = 0.0
        DO im=2,1,-1
          cmb = 0.0
          DO is=2,1,-1
            ctb = 0.0
            DO k=np+1,1,-1
              yy = tta(k-1, ih, im) - tda(k-1, ih, im)
              denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
              xx4 = tda(k-1, ih, im)*rra(k, im, is)
              fupdif = (xx4+yy*rxa(k, im, is))*denm
              fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
              fdndir = tda(k-1, ih, im)
              flxdn = fdndir + fdndif - fupdif
              flxdnb = ct*fallb(k)
              ctb = ctb + flxdn*fallb(k)
              fdndirb = flxdnb
              fdndifb = flxdnb
              fupdifb = -flxdnb
              tempb26 = denm*fupdifb
              denmb = (xx4*rsa(k-1, ih, im)+yy)*fdndifb + (xx4+yy*rxa(k&
&               , im, is))*fupdifb
              tempb27 = denm*fdndifb
              xx4b = rsa(k-1, ih, im)*tempb27 + tempb26
              yyb = tempb27 + rxa(k, im, is)*tempb26
              ttab(k-1, ih, im) = ttab(k-1, ih, im) + yyb
              tdab(k-1, ih, im) = tdab(k-1, ih, im) + rra(k, im, is)*&
&               xx4b + fdndirb - yyb
              rrab(k, im, is) = rrab(k, im, is) + tda(k-1, ih, im)*xx4b
              CALL POPREAL4(denm)
              temp17 = -(rsa(k-1, ih, im)*rxa(k, im, is)) + 1.
              tempb28 = -(denmb/temp17**2)
              rxab(k, im, is) = rxab(k, im, is) + yy*tempb26 - rsa(k-1, &
&               ih, im)*tempb28
              rsab(k-1, ih, im) = rsab(k-1, ih, im) + xx4*tempb27 - rxa(&
&               k, im, is)*tempb28
            END DO
            DO k=0,ict-1,1
              CALL POPREAL4(rxa(k, im, is))
              tempb23 = rxa(k+1, im, is)*rxab(k, im, is)
              CALL POPREAL4(rra(k, im, is))
              tempb25 = denm*rrab(k, im, is)
              temp16 = rxa(k+1, im, is)
              temp15 = tt(k, ih) - td(k, ih)
              rrb(k, ih) = rrb(k, ih) + rrab(k, im, is)
              tdb(k, ih) = tdb(k, ih) + (rra(k+1, im, is)-temp16)*&
&               tempb25
              rrab(k+1, im, is) = rrab(k+1, im, is) + td(k, ih)*tempb25
              denmb = (td(k, ih)*rra(k+1, im, is)+temp15*temp16)*rrab(k&
&               , im, is) + ts(k, ih)*tempb23
              ttb(k, ih) = ttb(k, ih) + temp16*tempb25
              rrab(k, im, is) = 0.0
              temp14 = -(rs(k, ih)*rxa(k+1, im, is)) + 1.
              tsb(k, ih) = tsb(k, ih) + denmb/temp14 + denm*tempb23
              tempb24 = -(ts(k, ih)*denmb/temp14**2)
              rsb(k, ih) = rsb(k, ih) + rxab(k, im, is) - rxa(k+1, im, &
&               is)*tempb24
              rxab(k+1, im, is) = rxab(k+1, im, is) + ts(k, ih)*denm*&
&               rxab(k, im, is)
              rxab(k, im, is) = 0.0
              rxab(k+1, im, is) = rxab(k+1, im, is) + temp15*tempb25 - &
&               rs(k, ih)*tempb24
              CALL POPREAL4(denm)
            END DO
            DO k=np,icb,-1
              CALL POPREAL4(rsa(k, ih, im))
              tempb20 = rsa(k-1, ih, im)*rsab(k, ih, im)
              CALL POPREAL4(tta(k, ih, im))
              tempb22 = denm*ttab(k, ih, im)
              temp13 = rsa(k-1, ih, im)
              temp12 = tda(k-1, ih, im)*rr(k, is)
              tdab(k-1, ih, im) = tdab(k-1, ih, im) + td(k, is)*tdab(k, &
&               ih, im) + (temp13*rr(k, is)-1.0)*tempb22 + tt(k, is)*&
&               ttab(k, ih, im)
              ttb(k, is) = ttb(k, is) + tda(k-1, ih, im)*ttab(k, ih, im)
              rrb(k, is) = rrb(k, is) + temp13*tda(k-1, ih, im)*tempb22
              ttab(k-1, ih, im) = ttab(k-1, ih, im) + tempb22
              denmb = (temp12*temp13+tta(k-1, ih, im)-tda(k-1, ih, im))*&
&               ttab(k, ih, im) + ts(k, is)*tempb20
              ttab(k, ih, im) = 0.0
              CALL POPREAL4(tda(k, ih, im))
              tdb(k, is) = tdb(k, is) + tda(k-1, ih, im)*tdab(k, ih, im)
              tdab(k, ih, im) = 0.0
              temp11 = -(rsa(k-1, ih, im)*rs(k, is)) + 1.
              tsb(k, is) = tsb(k, is) + denmb/temp11 + denm*tempb20
              tempb21 = -(ts(k, is)*denmb/temp11**2)
              rsb(k, is) = rsb(k, is) + rsab(k, ih, im) - rsa(k-1, ih, &
&               im)*tempb21
              rsab(k-1, ih, im) = rsab(k-1, ih, im) + ts(k, is)*denm*&
&               rsab(k, ih, im)
              rsab(k, ih, im) = 0.0
              rsab(k-1, ih, im) = rsab(k-1, ih, im) + temp12*tempb22 - &
&               rs(k, is)*tempb21
              CALL POPREAL4(denm)
            END DO
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL4(ct)
              cmb = cmb + cc3*ctb
              cc3b = cc3b + cm*ctb
            ELSE
              CALL POPREAL4(ct)
              cmb = cmb + (1.0-cc3)*ctb
              cc3b = cc3b - cm*ctb
            END IF
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4(cm)
            chb = chb + cc2*cmb
            cc2b = cc2b + ch*cmb
          ELSE
            CALL POPREAL4(cm)
            chb = chb + (1.0-cc2)*cmb
            cc2b = cc2b - ch*cmb
          END IF
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL4(ch)
          cc1b = cc1b + chb
        ELSE
          CALL POPREAL4(ch)
          cc1b = cc1b - chb
        END IF
      END DO
      DO is=2,1,-1
        DO k=ict,icb-1,1
          DO im=2,1,-1
            CALL POPREAL4(rxa(k, im, is))
            tempb17 = rxa(k+1, im, is)*rxab(k, im, is)
            CALL POPREAL4(rra(k, im, is))
            tempb19 = denm*rrab(k, im, is)
            temp10 = rxa(k+1, im, is)
            temp9 = tt(k, im) - td(k, im)
            rrb(k, im) = rrb(k, im) + rrab(k, im, is)
            tdb(k, im) = tdb(k, im) + (rra(k+1, im, is)-temp10)*tempb19
            rrab(k+1, im, is) = rrab(k+1, im, is) + td(k, im)*tempb19
            denmb = (td(k, im)*rra(k+1, im, is)+temp9*temp10)*rrab(k, im&
&             , is) + ts(k, im)*tempb17
            ttb(k, im) = ttb(k, im) + temp10*tempb19
            rrab(k, im, is) = 0.0
            temp8 = -(rs(k, im)*rxa(k+1, im, is)) + 1.
            tsb(k, im) = tsb(k, im) + denmb/temp8 + denm*tempb17
            tempb18 = -(ts(k, im)*denmb/temp8**2)
            rsb(k, im) = rsb(k, im) + rxab(k, im, is) - rxa(k+1, im, is)&
&             *tempb18
            rxab(k+1, im, is) = rxab(k+1, im, is) + ts(k, im)*denm*rxab(&
&             k, im, is)
            rxab(k, im, is) = 0.0
            rxab(k+1, im, is) = rxab(k+1, im, is) + temp9*tempb19 - rs(k&
&             , im)*tempb18
            CALL POPREAL4(denm)
          END DO
        END DO
        DO k=icb,np,1
          CALL POPREAL4(rxa(k, 2, is))
          rxab(k, 1, is) = rxab(k, 1, is) + rxab(k, 2, is)
          rxab(k, 2, is) = 0.0
          CALL POPREAL4(rra(k, 2, is))
          rrab(k, 1, is) = rrab(k, 1, is) + rrab(k, 2, is)
          rrab(k, 2, is) = 0.0
          CALL POPREAL4(rxa(k, 1, is))
          tempb14 = rxa(k+1, 1, is)*rxab(k, 1, is)
          CALL POPREAL4(rra(k, 1, is))
          tempb16 = denm*rrab(k, 1, is)
          temp7 = rxa(k+1, 1, is)
          temp6 = tt(k, is) - td(k, is)
          rrb(k, is) = rrb(k, is) + rrab(k, 1, is)
          tdb(k, is) = tdb(k, is) + (rra(k+1, 1, is)-temp7)*tempb16
          rrab(k+1, 1, is) = rrab(k+1, 1, is) + td(k, is)*tempb16
          denmb = (td(k, is)*rra(k+1, 1, is)+temp6*temp7)*rrab(k, 1, is)&
&           + ts(k, is)*tempb14
          ttb(k, is) = ttb(k, is) + temp7*tempb16
          rrab(k, 1, is) = 0.0
          temp5 = -(rs(k, is)*rxa(k+1, 1, is)) + 1.
          tsb(k, is) = tsb(k, is) + denmb/temp5 + denm*tempb14
          tempb15 = -(ts(k, is)*denmb/temp5**2)
          rsb(k, is) = rsb(k, is) + rxab(k, 1, is) - rxa(k+1, 1, is)*&
&           tempb15
          rxab(k+1, 1, is) = rxab(k+1, 1, is) + ts(k, is)*denm*rxab(k, 1&
&           , is)
          rxab(k, 1, is) = 0.0
          rxab(k+1, 1, is) = rxab(k+1, 1, is) + temp6*tempb16 - rs(k, is&
&           )*tempb15
          CALL POPREAL4(denm)
        END DO
        CALL POPREAL4(rxa(np+1, 2, is))
        rsb(np+1, is) = rsb(np+1, is) + rxab(np+1, 2, is)
        rxab(np+1, 2, is) = 0.0
        CALL POPREAL4(rra(np+1, 2, is))
        rrb(np+1, is) = rrb(np+1, is) + rrab(np+1, 2, is)
        rrab(np+1, 2, is) = 0.0
        CALL POPREAL4(rxa(np+1, 1, is))
        rsb(np+1, is) = rsb(np+1, is) + rxab(np+1, 1, is)
        rxab(np+1, 1, is) = 0.0
        CALL POPREAL4(rra(np+1, 1, is))
        rrb(np+1, is) = rrb(np+1, is) + rrab(np+1, 1, is)
        rrab(np+1, 1, is) = 0.0
      END DO
      DO ih=2,1,-1
        DO k=icb-1,ict,-1
          DO im=2,1,-1
            CALL POPREAL4(rsa(k, ih, im))
            tempb11 = rsa(k-1, ih, im)*rsab(k, ih, im)
            CALL POPREAL4(tta(k, ih, im))
            tempb13 = denm*ttab(k, ih, im)
            temp4 = rsa(k-1, ih, im)
            temp3 = tda(k-1, ih, im)*rr(k, im)
            tdab(k-1, ih, im) = tdab(k-1, ih, im) + td(k, im)*tdab(k, ih&
&             , im) + (temp4*rr(k, im)-1.0)*tempb13 + tt(k, im)*ttab(k, &
&             ih, im)
            ttb(k, im) = ttb(k, im) + tda(k-1, ih, im)*ttab(k, ih, im)
            rrb(k, im) = rrb(k, im) + temp4*tda(k-1, ih, im)*tempb13
            ttab(k-1, ih, im) = ttab(k-1, ih, im) + tempb13
            denmb = (temp3*temp4+tta(k-1, ih, im)-tda(k-1, ih, im))*ttab&
&             (k, ih, im) + ts(k, im)*tempb11
            ttab(k, ih, im) = 0.0
            CALL POPREAL4(tda(k, ih, im))
            tdb(k, im) = tdb(k, im) + tda(k-1, ih, im)*tdab(k, ih, im)
            tdab(k, ih, im) = 0.0
            temp2 = -(rsa(k-1, ih, im)*rs(k, im)) + 1.
            tsb(k, im) = tsb(k, im) + denmb/temp2 + denm*tempb11
            tempb12 = -(ts(k, im)*denmb/temp2**2)
            rsb(k, im) = rsb(k, im) + rsab(k, ih, im) - rsa(k-1, ih, im)&
&             *tempb12
            rsab(k-1, ih, im) = rsab(k-1, ih, im) + ts(k, im)*denm*rsab(&
&             k, ih, im)
            rsab(k, ih, im) = 0.0
            rsab(k-1, ih, im) = rsab(k-1, ih, im) + temp3*tempb13 - rs(k&
&             , im)*tempb12
            CALL POPREAL4(denm)
          END DO
        END DO
        DO k=ict-1,1,-1
          CALL POPREAL4(rsa(k, ih, 2))
          rsab(k, ih, 1) = rsab(k, ih, 1) + rsab(k, ih, 2)
          rsab(k, ih, 2) = 0.0
          CALL POPREAL4(tta(k, ih, 2))
          ttab(k, ih, 1) = ttab(k, ih, 1) + ttab(k, ih, 2)
          ttab(k, ih, 2) = 0.0
          CALL POPREAL4(tda(k, ih, 2))
          tdab(k, ih, 1) = tdab(k, ih, 1) + tdab(k, ih, 2)
          tdab(k, ih, 2) = 0.0
          CALL POPREAL4(rsa(k, ih, 1))
          tempb8 = rsa(k-1, ih, 1)*rsab(k, ih, 1)
          CALL POPREAL4(tta(k, ih, 1))
          tempb10 = denm*ttab(k, ih, 1)
          temp1 = rsa(k-1, ih, 1)
          temp0 = tda(k-1, ih, 1)*rr(k, ih)
          tdab(k-1, ih, 1) = tdab(k-1, ih, 1) + td(k, ih)*tdab(k, ih, 1)&
&           + (temp1*rr(k, ih)-1.0)*tempb10 + tt(k, ih)*ttab(k, ih, 1)
          ttb(k, ih) = ttb(k, ih) + tda(k-1, ih, 1)*ttab(k, ih, 1)
          rrb(k, ih) = rrb(k, ih) + temp1*tda(k-1, ih, 1)*tempb10
          ttab(k-1, ih, 1) = ttab(k-1, ih, 1) + tempb10
          denmb = (temp0*temp1+tta(k-1, ih, 1)-tda(k-1, ih, 1))*ttab(k, &
&           ih, 1) + ts(k, ih)*tempb8
          ttab(k, ih, 1) = 0.0
          CALL POPREAL4(tda(k, ih, 1))
          tdb(k, ih) = tdb(k, ih) + tda(k-1, ih, 1)*tdab(k, ih, 1)
          tdab(k, ih, 1) = 0.0
          temp = -(rsa(k-1, ih, 1)*rs(k, ih)) + 1.
          tsb(k, ih) = tsb(k, ih) + denmb/temp + denm*tempb8
          tempb9 = -(ts(k, ih)*denmb/temp**2)
          rsb(k, ih) = rsb(k, ih) + rsab(k, ih, 1) - rsa(k-1, ih, 1)*&
&           tempb9
          rsab(k-1, ih, 1) = rsab(k-1, ih, 1) + ts(k, ih)*denm*rsab(k, &
&           ih, 1)
          rsab(k, ih, 1) = 0.0
          rsab(k-1, ih, 1) = rsab(k-1, ih, 1) + temp0*tempb10 - rs(k, ih&
&           )*tempb9
          CALL POPREAL4(denm)
        END DO
        CALL POPREAL4(rsa(0, ih, 2))
        rsb(0, ih) = rsb(0, ih) + rsab(0, ih, 2)
        rsab(0, ih, 2) = 0.0
        CALL POPREAL4(tta(0, ih, 2))
        ttb(0, ih) = ttb(0, ih) + ttab(0, ih, 2)
        ttab(0, ih, 2) = 0.0
        CALL POPREAL4(tda(0, ih, 2))
        tdb(0, ih) = tdb(0, ih) + tdab(0, ih, 2)
        tdab(0, ih, 2) = 0.0
        CALL POPREAL4(rsa(0, ih, 1))
        rsb(0, ih) = rsb(0, ih) + rsab(0, ih, 1)
        rsab(0, ih, 1) = 0.0
        CALL POPREAL4(tta(0, ih, 1))
        ttb(0, ih) = ttb(0, ih) + ttab(0, ih, 1)
        ttab(0, ih, 1) = 0.0
        CALL POPREAL4(tda(0, ih, 1))
        tdb(0, ih) = tdb(0, ih) + tdab(0, ih, 1)
        tdab(0, ih, 1) = 0.0
      END DO
      DO k=np+1,1,-1
        fallb(k) = 0.0
      END DO
      DO k=np,1,-1
        CALL POPREAL4(ts(k, 2))
        tstb = tsb(k, 2)
        tsb(k, 2) = 0.0
        CALL POPREAL4(rs(k, 2))
        rstb = rsb(k, 2)
        rsb(k, 2) = 0.0
        CALL POPREAL4(td(k, 2))
        tdtb = tdb(k, 2)
        tdb(k, 2) = 0.0
        CALL POPREAL4(tt(k, 2))
        tttb = ttb(k, 2)
        ttb(k, 2) = 0.0
        CALL POPREAL4(rr(k, 2))
        rrtb = rrb(k, 2)
        rrb(k, 2) = 0.0
        asysto = asya_dev(i, k, ib)
        asytof = (asysto+asycl(k)*tauclf(k))/(ssatof*tautof)
        tautofb = 0.0
        ssatofb = 0.0
        asytofb = 0.0
        dumb = 0.0
        CALL DELEDD_B(tautof, tautofb, ssatof, ssatofb, asytof, asytofb&
&               , dsm, rst, rstb, tst, tstb, dum, dumb)
        asytob = (asysto+asycl(k)*tauclb(k))/(ssatob*tautob)
        tautobb = 0.0
        ssatobb = 0.0
        asytobb = 0.0
        CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, asytobb&
&               , cosz_dev(i), rrt, rrtb, ttt, tttb, tdt, tdtb)
        tempb6 = asytofb/(ssatof*tautof)
        tempb7 = -((asysto+asycl(k)*tauclf(k))*tempb6/(ssatof*tautof))
        asystob = tempb6
        asyclb(k) = asyclb(k) + tauclf(k)*tempb6
        tauclfb(k) = tauclfb(k) + asycl(k)*tempb6
        ssatofb = ssatofb + tautof*tempb7
        tautofb = tautofb + ssatof*tempb7
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          taurs = ry_uv(ib)*dp(k)
          ssatau = ssaa_dev(i, k, ib) + taurs
          ssatofb = 0.0
        ELSE
          taurs = ry_uv(ib)*dp(k)
          ssatau = ssaa_dev(i, k, ib) + taurs
        END IF
        CALL POPREAL4(ssatof)
        tempb3 = ssatofb/tautof
        ssataub = tempb3
        tautofb = tautofb - (ssatau+tauclf(k))*tempb3/tautof
        tauclfb(k) = tauclfb(k) + tautofb + tempb3
        CALL POPREAL4(tautof)
        taustob = tautofb
        tempb4 = asytobb/(ssatob*tautob)
        tempb5 = -((asysto+asycl(k)*tauclb(k))*tempb4/(ssatob*tautob))
        asystob = asystob + tempb4
        asyclb(k) = asyclb(k) + tauclb(k)*tempb4
        tauclbb(k) = tauclbb(k) + asycl(k)*tempb4
        ssatobb = ssatobb + tautob*tempb5
        tautobb = tautobb + ssatob*tempb5
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) ssatobb = 0.0
        CALL POPREAL4(ssatob)
        tempb2 = ssatobb/tautob
        ssataub = ssataub + tempb2
        tautobb = tautobb - (ssatau+tauclb(k))*tempb2/tautob
        tauclbb(k) = tauclbb(k) + tautobb + tempb2
        CALL POPREAL4(tautob)
        taustob = taustob + tautobb
        CALL POPREAL4(ts(k, 1))
        tstb = tsb(k, 1)
        tsb(k, 1) = 0.0
        CALL POPREAL4(rs(k, 1))
        rstb = rsb(k, 1)
        rsb(k, 1) = 0.0
        CALL POPREAL4(td(k, 1))
        tdtb = tdb(k, 1)
        tdb(k, 1) = 0.0
        CALL POPREAL4(tt(k, 1))
        tttb = ttb(k, 1)
        ttb(k, 1) = 0.0
        CALL POPREAL4(rr(k, 1))
        rrtb = rrb(k, 1)
        rrb(k, 1) = 0.0
        asytob = asysto/ssatau
        tautobb = 0.0
        ssatobb = 0.0
        asytobb = 0.0
        dumb = 0.0
        CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, asytobb&
&               , dsm, rst, rstb, tst, tstb, dum, dumb)
        CALL DELEDD_B(tautob, tautobb, ssatob, ssatobb, asytob, asytobb&
&               , cosz_dev(i), rrt, rrtb, ttt, tttb, tdt, tdtb)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) ssatobb = 0.0
        CALL POPREAL4(ssatob)
        ssataub = ssataub + ssatobb/tautob - asysto*asytobb/ssatau**2
        tautobb = tautobb - ssatau*ssatobb/tautob**2
        asystob = asystob + asytobb/ssatau
        CALL POPREAL4(tautob)
        taustob = taustob + tautobb
        asya_devb(i, k, ib) = asya_devb(i, k, ib) + asystob
        ssaa_devb(i, k, ib) = ssaa_devb(i, k, ib) + ssataub
        tauozb = taustob
        tauwvb = taustob
        taua_devb(i, k, ib) = taua_devb(i, k, ib) + taustob
        whb(k) = whb(k) + wk_uv(ib)*tauwvb
        ohb(k) = ohb(k) + zk_uv(ib)*tauozb
      END DO
      CALL POPREAL4(td(0, 2))
      tdb(0, 1) = tdb(0, 1) + tdb(0, 2)
      tdb(0, 2) = 0.0
      CALL POPREAL4(td(0, 1))
      tempb1 = -(EXP(-((wk_uv(ib)*wvtoa+zk_uv(ib)*o3toa)/cosz_dev(i)))*&
&       tdb(0, 1)/cosz_dev(i))
      wvtoab = wvtoab + wk_uv(ib)*tempb1
      o3toab = o3toab + zk_uv(ib)*tempb1
      tdb(0, 1) = 0.0
    END DO
    taudiffb = 0.0
    taubeamb = 0.0
    DO k=np,1,-1
      CALL POPREAL4(tauclf(k))
      taudiffb(k, 1) = taudiffb(k, 1) + tauclfb(k)
      taudiffb(k, 2) = taudiffb(k, 2) + tauclfb(k)
      taudiffb(k, 3) = taudiffb(k, 3) + tauclfb(k)
      taudiffb(k, 4) = taudiffb(k, 4) + tauclfb(k)
      tauclfb(k) = 0.0
      CALL POPREAL4(tauclb(k))
      taubeamb(k, 1) = taubeamb(k, 1) + tauclbb(k)
      taubeamb(k, 2) = taubeamb(k, 2) + tauclbb(k)
      taubeamb(k, 3) = taubeamb(k, 3) + tauclbb(k)
      taubeamb(k, 4) = taubeamb(k, 4) + tauclbb(k)
      tauclbb(k) = 0.0
    END DO
    DO k=np,1,-1
      CALL POPCONTROL3B(branch)
      IF (branch .LT. 3) THEN
        IF (branch .EQ. 0) THEN
          fcld_devb(i, k) = fcld_devb(i, k) + cc3b
          cc3b = 0.0
        ELSE IF (branch .NE. 1) THEN
          fcld_devb(i, k) = fcld_devb(i, k) + cc2b
          cc2b = 0.0
        END IF
      ELSE IF (branch .NE. 3) THEN
        IF (branch .EQ. 4) THEN
          fcld_devb(i, k) = fcld_devb(i, k) + cc1b
          cc1b = 0.0
        END IF
      END IF
    END DO
    CALL POPREAL4ARRAY(asycl, np)
    CALL GETVISTAU1_B(np, cosz_dev(i), dp_pa, fcld_col, fcld_colb, &
&               reff_col, reff_colb, cwc_col, cwc_colb, ict, icb, &
&               taubeam, taubeamb, taudiff, taudiffb, asycl, asyclb, &
&               aig_uv, awg_uv, arg_uv, aib_uv, awb_uv, arb_uv, aib_nir&
&               , awb_nir, arb_nir, aia_nir, awa_nir, ara_nir, aig_nir, &
&               awg_nir, arg_nir, caib, caif, cons_grav)
    CALL POPREAL4(cc3)
    CALL POPREAL4(cc2)
    DO k=np+1,1,-1
      flx_devb(i, k) = 0.0
    END DO
    CALL POPREAL4ARRAY(rsa, (np+1)*2**2)
    CALL POPREAL4ARRAY(tda, (np+1)*2**2)
    CALL POPREAL4ARRAY(tta, (np+1)*2**2)
    CALL POPREAL4ARRAY(rxa, (np+2)*2**2)
    CALL POPREAL4ARRAY(rra, (np+2)*2**2)
    CALL POPREAL4ARRAY(ts, (np+2)*2)
    CALL POPREAL4ARRAY(rs, (np+2)*2)
    CALL POPREAL4ARRAY(td, (np+2)*2)
    CALL POPREAL4ARRAY(tt, (np+2)*2)
    CALL POPREAL4ARRAY(rr, (np+2)*2)
    DO k=np,1,-1
      DO l=4,1,-1
        CALL POPREAL4(cwc_col(k, l))
        cwc_devb(i, k, l) = cwc_devb(i, k, l) + cwc_colb(k, l)
        cwc_colb(k, l) = 0.0
        CALL POPREAL4(reff_col(k, l))
        reff_devb(i, k, l) = reff_devb(i, k, l) + reff_colb(k, l)
        reff_colb(k, l) = 0.0
      END DO
      CALL POPREAL4(fcld_col(k))
      fcld_devb(i, k) = fcld_devb(i, k) + fcld_colb(k)
      fcld_colb(k) = 0.0
      oa_devb(i, k) = oa_devb(i, k) + dp(k)*466.7*1.02*ohb(k)
      ohb(k) = 0.0
      CALL POPREAL4(swh(k+1))
      swhb(k) = swhb(k) + swhb(k+1)
      whb(k) = whb(k) + swhb(k+1)
      swhb(k+1) = 0.0
      tempb0 = scal(k)*1.02*whb(k)
      wa_devb(i, k) = wa_devb(i, k) + (0.00135*(ta_dev(i, k)-240.)+1.)*&
&       tempb0
      ta_devb(i, k) = ta_devb(i, k) + wa_dev(i, k)*0.00135*tempb0
      whb(k) = 0.0
      CALL POPREAL4(scal(k))
      CALL POPREAL4(dp_pa(k))
      CALL POPREAL4(dp(k))
    END DO
    CALL POPREAL4(swh(1))
    wvtoab = wvtoab + swhb(1)
    swhb(1) = 0.0
    tempb = scal0*1.02*wvtoab
    wa_devb(i, 1) = wa_devb(i, 1) + (0.00135*(ta_dev(i, 1)-240.)+1.0)*&
&     tempb
    ta_devb(i, 1) = ta_devb(i, 1) + wa_dev(i, 1)*0.00135*tempb
    oa_devb(i, 1) = oa_devb(i, 1) + xtoa*466.7*1.02*o3toab
    CALL POPREAL4(scal0)
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(xtoa)
    ELSE
      CALL POPREAL4(xtoa)
    END IF
    CALL POPREAL4(snt)
    CALL POPINTEGER4(ntop)
  END DO
END SUBROUTINE SORAD_B

!  Differentiation of deledd in reverse (adjoint) mode:
!   gradient     of useful results: g01 tt1 td1 rr1 tau1 ssc1
!   with respect to varying inputs: g01 tau1 ssc1
!*********************************************************************
SUBROUTINE DELEDD_B(tau1, tau1b, ssc1, ssc1b, g01, g01b, cza1, rr1, rr1b&
& , tt1, tt1b, td1, td1b)
  IMPLICIT NONE
! 8 byte real
  INTEGER, PARAMETER :: real_de=8
!integer,parameter :: REAL_SP = 4 ! 4 byte real
!-----input parameters
  REAL*4, INTENT(IN) :: tau1, ssc1, g01, cza1
  REAL*4 :: tau1b, ssc1b, g01b
!-----output parameters
  REAL*4 :: rr1, tt1, td1
  REAL*4 :: rr1b, tt1b, td1b
!-----temporary parameters
  REAL*8, PARAMETER :: zero=0.0_REAL_DE
  REAL*8, PARAMETER :: one=1.0_REAL_DE
  REAL*8, PARAMETER :: two=2.0_REAL_DE
  REAL*8, PARAMETER :: three=3.0_REAL_DE
  REAL*8, PARAMETER :: four=4.0_REAL_DE
  REAL*8, PARAMETER :: fourth=0.25_REAL_DE
  REAL*8, PARAMETER :: seven=7.0_REAL_DE
  REAL*8, PARAMETER :: thresh=1.e-8_REAL_DE
  REAL*8 :: tau, ssc, g0, rr, tt, td
  REAL*8 :: taub, sscb, g0b, rrb, ttb, tdb
  REAL*8 :: zth, ff, xx, taup, sscp, gp, gm1, gm2, gm3, akk, alf1, alf2
  REAL*8 :: ffb, xxb, taupb, sscpb, gpb, gm1b, gm2b, gm3b, akkb, alf1b, &
& alf2b
  REAL*8 :: all, bll, st7, st8, cll, dll, fll, ell, st1, st2, st3, st4
  REAL*8 :: allb, bllb, st7b, st8b, cllb, dllb, fllb, ellb, st1b, st2b, &
& st3b, st4b
  INTRINSIC DBLE
  INTRINSIC SQRT
  INTRINSIC ABS
  INTRINSIC EXP
  INTRINSIC MAX
  INTRINSIC REAL
  INTEGER :: branch
  REAL*8 :: tempb9
  REAL*8 :: tempb8
  REAL*8 :: tempb7
  REAL*8 :: tempb6
  REAL*8 :: tempb5
  REAL*8 :: tempb4
  REAL*8 :: tempb3
  REAL*8 :: tempb2
  REAL*8 :: tempb1
  REAL*8 :: tempb0
  REAL*8 :: tempb10
  REAL*8 :: tempb
  REAL*8 :: abs0
  REAL*8 :: temp
!zth = real(cza1,kind=REAL_DE)
!g0  = real(g01 ,kind=REAL_DE)
!tau = real(tau1,kind=REAL_DE)
!ssc = real(ssc1,kind=REAL_DE)
  zth = DBLE(cza1)
  g0 = DBLE(g01)
  tau = DBLE(tau1)
  ssc = DBLE(ssc1)
  ff = g0*g0
  xx = one - ff*ssc
  taup = tau*xx
  sscp = ssc*(one-ff)/xx
  gp = g0/(one+g0)
  xx = three*gp
  gm1 = (seven-sscp*(four+xx))*fourth
  gm2 = -((one-sscp*(four-xx))*fourth)
  akk = SQRT((gm1+gm2)*(gm1-gm2))
  xx = akk*zth
  st7 = one - xx
  st8 = one + xx
  st3 = st7*st8
  IF (st3 .GE. 0.) THEN
    abs0 = st3
  ELSE
    abs0 = -st3
  END IF
  IF (abs0 .LT. thresh) THEN
    CALL PUSHREAL8(zth)
    zth = zth + 0.0010
    IF (zth .GT. 1.0) zth = zth - 0.0020
    xx = akk*zth
    CALL PUSHREAL8(st7)
    st7 = one - xx
    CALL PUSHREAL8(st8)
    st8 = one + xx
    st3 = st7*st8
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
  END IF
  td = EXP(-(taup/zth))
  gm3 = (two-zth*three*gp)*fourth
  xx = gm1 - gm2
  alf1 = gm1 - gm3*xx
  alf2 = gm2 + gm3*xx
  xx = akk*two
  all = (gm3-alf2*zth)*xx*td
  bll = (one-gm3+alf1*zth)*xx
  xx = akk*gm3
  cll = (alf2+xx)*st7
  dll = (alf2-xx)*st8
  xx = akk*(one-gm3)
  fll = (alf1+xx)*st8
  ell = (alf1-xx)*st7
  st2 = EXP(-(akk*taup))
  st4 = st2*st2
  st1 = sscp/((akk+gm1+(akk-gm1)*st4)*st3)
  rr = (cll-dll*st4-all*st2)*st1
  tt = -(((fll-ell*st4)*td-bll*st2)*st1)
  IF (rr .LT. zero) THEN
    rr = zero
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    rr = rr
  END IF
  IF (tt .LT. zero) THEN
    tt = zero
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
    tt = tt
  END IF
  tt = tt + td
!td1 = real(td,kind=REAL_SP)
!rr1 = real(rr,kind=REAL_SP)
!tt1 = real(tt,kind=REAL_SP)
  ttb = tt1b
  rrb = rr1b
  tdb = ttb + td1b
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) ttb = 0.0_8
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) rrb = 0.0_8
  st1b = (cll-dll*st4-all*st2)*rrb - ((fll-ell*st4)*td-bll*st2)*ttb
  tempb4 = st1*rrb
  temp = akk + gm1 + (akk-gm1)*st4
  tempb7 = st1b/(temp*st3)
  tempb8 = -(sscp*tempb7/(temp*st3))
  tempb5 = st3*tempb8
  tempb2 = -(st1*ttb)
  tempb3 = td*tempb2
  fllb = tempb3
  ellb = -(st4*tempb3)
  st4b = (akk-gm1)*tempb5 - dll*tempb4 - ell*tempb3
  bllb = -(st2*tempb2)
  st2b = 2*st2*st4b - all*tempb4 - bll*tempb2
  cllb = tempb4
  dllb = -(st4*tempb4)
  allb = -(st2*tempb4)
  sscpb = tempb7
  st3b = temp*tempb8
  tempb9 = EXP(-(akk*taup))*st2b
  xxb = st8*fllb - st7*ellb
  akkb = (one-gm3)*xxb - taup*tempb9 + (st4+1.0_8)*tempb5
  st7b = (alf1-xx)*ellb
  st8b = (alf1+xx)*fllb
  gm3b = -(akk*xxb)
  xx = akk*gm3
  xxb = st7*cllb - st8*dllb
  st8b = st8b + (alf2-xx)*dllb
  st7b = st7b + (alf2+xx)*cllb
  akkb = akkb + gm3*xxb
  xx = akk*two
  alf1b = st8*fllb + xx*zth*bllb + st7*ellb
  gm3b = gm3b + akk*xxb - xx*bllb
  tempb10 = xx*td*allb
  alf2b = st7*cllb - zth*tempb10 + st8*dllb
  tempb6 = (gm3-zth*alf2)*allb
  tdb = tdb + xx*tempb6 + (fll-ell*st4)*tempb2
  taupb = -(EXP(-(taup/zth))*tdb/zth) - akk*tempb9
  xxb = td*tempb6 + (one+zth*alf1-gm3)*bllb
  akkb = akkb + two*xxb
  xx = gm1 - gm2
  gm3b = gm3b + xx*alf2b - xx*alf1b + tempb10
  xxb = gm3*alf2b - gm3*alf1b
  gm1b = alf1b + xxb + (1.0_8-st4)*tempb5
  gm2b = alf2b - xxb
  gpb = -(fourth*zth*three*gm3b)
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) THEN
    st7b = st7b + st8*st3b
    st8b = st8b + st7*st3b
    CALL POPREAL8(st8)
    xxb = st8b - st7b
    CALL POPREAL8(st7)
    akkb = akkb + zth*xxb
    CALL POPREAL8(zth)
    st3b = 0.0_8
    st7b = 0.0_8
    st8b = 0.0_8
  END IF
  st7b = st7b + st8*st3b
  st8b = st8b + st7*st3b
  xxb = st8b - st7b
  akkb = akkb + zth*xxb
  xx = three*gp
  IF ((gm1+gm2)*(gm1-gm2) .EQ. 0.0) THEN
    tempb = 0.0
  ELSE
    tempb = akkb/(2.0*SQRT((gm1+gm2)*(gm1-gm2)))
  END IF
  gm1b = gm1b + 2*gm1*tempb
  gm2b = gm2b - 2*gm2*tempb
  sscpb = sscpb + fourth*(four-xx)*gm2b - fourth*(four+xx)*gm1b
  xxb = -(fourth*sscp*gm1b) - sscp*fourth*gm2b
  gpb = gpb + three*xxb
  xx = one - ff*ssc
  tempb0 = gpb/(one+g0)
  tempb1 = (one-ff)*sscpb/xx
  xxb = tau*taupb - ssc*tempb1/xx
  ffb = -(ssc*xxb) - ssc*sscpb/xx
  g0b = 2*g0*ffb + (1.0_8-g0/(one+g0))*tempb0
  sscb = tempb1 - ff*xxb
  taub = xx*taupb
  ssc1b = ssc1b + sscb
  tau1b = tau1b + taub
  g01b = g01b + g0b
END SUBROUTINE DELEDD_B

!  Differentiation of getvistau1 in reverse (adjoint) mode:
!   gradient     of useful results: hydromets asycl taudiff fcld
!                taubeam reff
!   with respect to varying inputs: hydromets asycl fcld reff
SUBROUTINE GETVISTAU1_B(nlevs, cosz, dp, fcld, fcldb, reff, reffb, &
& hydromets, hydrometsb, ict, icb, taubeam, taubeamb, taudiff, taudiffb&
& , asycl, asyclb, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv, arb_uv, &
& aib_nir, awb_nir, arb_nir, aia_nir, awa_nir, ara_nir, aig_nir, awg_nir&
& , arg_nir, caib, caif, cons_grav)
  IMPLICIT NONE
!EOC
! !INPUT PARAMETERS:
!  Number of levels
  INTEGER, INTENT(IN) :: nlevs
!  Cosine of solar zenith angle
  REAL, INTENT(IN) :: cosz
!  Delta pressure (Pa)
  REAL, INTENT(IN) :: dp(nlevs)
!  Cloud fraction (used sometimes)
  REAL, INTENT(IN) :: fcld(nlevs)
  REAL :: fcldb(nlevs)
!  Effective radius (microns)
  REAL, INTENT(IN) :: reff(nlevs, 4)
  REAL :: reffb(nlevs, 4)
!  Hydrometeors (kg/kg)
  REAL, INTENT(IN) :: hydromets(nlevs, 4)
  REAL :: hydrometsb(nlevs, 4)
!  Flags for various uses 
  INTEGER, INTENT(IN) :: ict, icb
!                 ict  = 0   Indicates that in-cloud values have been given
!                            and are expected
!                 ict != 0   Indicates that overlap computation is needed, and:
!                               ict is the level of the mid-high boundary
!                               icb is the level of the low-mid  boundary
!                
! !OUTPUT PARAMETERS:
!  Optical Depth for Beam Radiation
  REAL :: taubeam(nlevs, 4)
  REAL :: taubeamb(nlevs, 4)
!  Optical Depth for Diffuse Radiation
  REAL :: taudiff(nlevs, 4)
  REAL :: taudiffb(nlevs, 4)
!  Cloud Asymmetry Factor
  REAL :: asycl(nlevs)
  REAL :: asyclb(nlevs)
! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for visible wavelengths
!  In general will compute in-cloud - will do grid mean when called
!  for diagnostic use only. ict flag will indicate which to do.
!  Slots for reff, hydrometeors, taubeam, taudiff, and asycl are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
!
!  In the below calculations, the constants used in the tau calculation are in 
!  m$^2$ g$^{-1}$ and m$^2$ g$^{-1}$ $\mu$m. Thus, we must convert the kg contained in the 
!  pressure (Pa = kg m$^{-1}$ s$^{-2}$) to grams.
!
! !REVISION HISTORY: 
!    2011.10.27   Molod moved to Radiation_Shared and revised arg list, units
!    2011.11.16   MAT: Generalized to a call that is per-column
!
!EOP
!------------------------------------------------------------------------------
!BOC
  INTEGER :: k, in, im, it, ia, kk
  REAL :: fm, ft, fa, xai, tauc, asyclt
  REAL :: ftb, fab, xaib, taucb, asycltb
  REAL :: cc(3)
  REAL :: ccb(3)
  REAL :: taucld1, taucld2, taucld3, taucld4
  REAL :: taucld1b, taucld2b, taucld3b, taucld4b
  REAL :: g1, g2, g3, g4
  REAL :: g1b, g2b, g3b, g4b
  REAL :: reff_snow
  REAL :: reff_snowb
  INTEGER, PARAMETER :: nm=11, nt=9, na=11
  REAL, PARAMETER :: dm=0.1, dt=0.30103, da=0.1, t1=-0.9031
  REAL, INTENT(IN) :: aig_uv(3), awg_uv(3), arg_uv(3)
  REAL, INTENT(IN) :: aib_uv, awb_uv(2), arb_uv(2)
  REAL, INTENT(IN) :: aib_nir, awb_nir(3, 2), arb_nir(3, 2)
  REAL, INTENT(IN) :: aia_nir(3, 3), awa_nir(3, 3), ara_nir(3, 3)
  REAL, INTENT(IN) :: aig_nir(3, 3), awg_nir(3, 3), arg_nir(3, 3)
  REAL, INTENT(IN) :: caib(11, 9, 11), caif(9, 11)
  REAL, INTENT(IN) :: cons_grav
  INTRINSIC MAX
  INTRINSIC MIN
  INTRINSIC LOG10
  INTRINSIC INT
  INTRINSIC REAL
  INTEGER :: branch
  REAL :: temp3
  REAL :: temp2
  REAL :: temp1
  REAL :: temp0
  REAL :: tempb7
  REAL :: tempb6
  REAL :: tempb5
  REAL :: tempb4
  REAL :: tempb3
  REAL :: tempb2
  REAL :: tempb1
  REAL :: tempb0
  REAL :: tempb
  REAL :: temp
  IF (ict .NE. 0) THEN
!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
    cc = 0.0
    DO k=1,ict-1
      IF (cc(1) .LT. fcld(k)) THEN
        cc(1) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(1) = cc(1)
      END IF
    END DO
    DO k=ict,icb-1
      IF (cc(2) .LT. fcld(k)) THEN
        cc(2) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(2) = cc(2)
      END IF
    END DO
    DO k=icb,nlevs
      IF (cc(3) .LT. fcld(k)) THEN
        cc(3) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(3) = cc(3)
      END IF
    END DO
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
  END IF
!-----Compute cloud optical thickness.  Eqs. (4.6) and (4.10)
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     taucld1 is the optical thickness for ice particles
!     taucld2 is the optical thickness for liquid particles
!     taucld3 is the optical thickness for rain drops
!     taucld4 is the optical thickness for snow
  DO k=1,nlevs
    IF (reff(k, 1) .LE. 0.) THEN
      CALL PUSHREAL4(taucld1)
      taucld1 = 0.
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL4(taucld1)
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*aib_uv/reff(k, 1)
      CALL PUSHCONTROL1B(1)
    END IF
    IF (reff(k, 2) .LE. 0.) THEN
      CALL PUSHREAL4(taucld2)
      taucld2 = 0.
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL4(taucld2)
      taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_uv(1)+awb_uv(&
&       2)/reff(k, 2))
      CALL PUSHCONTROL1B(0)
    END IF
    CALL PUSHREAL4(taucld3)
    taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_uv(1)
    IF (reff(k, 4) .GT. 112.0) THEN
      CALL PUSHREAL4(reff_snow)
      reff_snow = 112.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL4(reff_snow)
      reff_snow = reff(k, 4)
      CALL PUSHCONTROL1B(1)
    END IF
    IF (reff_snow .LE. 0.) THEN
      CALL PUSHREAL4(taucld4)
      taucld4 = 0.
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL4(taucld4)
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*aib_uv/reff_snow
      CALL PUSHCONTROL1B(1)
    END IF
    IF (ict .NE. 0) THEN
!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
      IF (k .LT. ict) THEN
        CALL PUSHINTEGER4(kk)
        kk = 1
        CALL PUSHCONTROL2B(0)
      ELSE IF (k .GE. ict .AND. k .LT. icb) THEN
        CALL PUSHINTEGER4(kk)
        kk = 2
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHINTEGER4(kk)
        kk = 3
        CALL PUSHCONTROL2B(2)
      END IF
      tauc = taucld1 + taucld2 + taucld3 + taucld4
      IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
!-----normalize cloud cover following Eq. (7.8)
        CALL PUSHREAL4(fa)
        fa = fcld(k)/cc(kk)
        IF (tauc .GT. 32.) THEN
          tauc = 32.
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          tauc = tauc
        END IF
        fm = cosz/dm
        CALL PUSHREAL4(ft)
        ft = (LOG10(tauc)-t1)/dt
        fa = fa/da
        CALL PUSHINTEGER4(im)
        im = INT(fm + 1.5)
        CALL PUSHINTEGER4(it)
        it = INT(ft + 1.5)
        CALL PUSHINTEGER4(ia)
        ia = INT(fa + 1.5)
        IF (im .LT. 2) THEN
          im = 2
        ELSE
          im = im
        END IF
        IF (it .LT. 2) THEN
          it = 2
        ELSE
          it = it
        END IF
        IF (ia .LT. 2) THEN
          ia = 2
        ELSE
          ia = ia
        END IF
        IF (im .GT. nm - 1) THEN
          im = nm - 1
        ELSE
          im = im
        END IF
        IF (it .GT. nt - 1) THEN
          it = nt - 1
        ELSE
          it = it
        END IF
        IF (ia .GT. na - 1) THEN
          ia = na - 1
        ELSE
          ia = ia
        END IF
        fm = fm - REAL(im - 1)
        ft = ft - REAL(it - 1)
        fa = fa - REAL(ia - 1)
!-----scale cloud optical thickness for beam radiation following 
!     Eq. (7.3).
!     the scaling factor, xai, is a function of the solar zenith
!     angle, optical thickness, and cloud cover.
        CALL PUSHREAL4(xai)
        xai = (-(caib(im-1, it, ia)*(1.-fm))+caib(im+1, it, ia)*(1.+fm))&
&         *fm*.5 + caib(im, it, ia)*(1.-fm*fm)
        xai = xai + (-(caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(&
&         1.+ft))*ft*.5 + caib(im, it, ia)*(1.-ft*ft)
        xai = xai + (-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(&
&         1.+fa))*fa*.5 + caib(im, it, ia)*(1.-fa*fa)
        xai = xai - 2.*caib(im, it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
        CALL PUSHREAL4(xai)
        xai = (-(caif(it-1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ft*.5 +&
&         caif(it, ia)*(1.-ft*ft)
        xai = xai + (-(caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*&
&         fa*.5 + caif(it, ia)*(1.-fa*fa)
        xai = xai - caif(it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
!-----cloud asymmetry factor for a mixture of liquid and ice particles.
!     unit of reff is micrometers. Eqs. (4.8) and (6.4)
    CALL PUSHREAL4(tauc)
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
  END DO
  ccb = 0.0
  DO k=nlevs,1,-1
    asycltb = asyclb(k)
    asyclb(k) = 0.0
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      g1 = (aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff(k, 1))*reff(k, 1))*&
&       taucld1
      g2 = (awg_uv(1)+(awg_uv(2)+awg_uv(3)*reff(k, 2))*reff(k, 2))*&
&       taucld2
      taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_uv(1)
      g3 = arg_uv(1)*taucld3
      g4 = (aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff_snow)*reff_snow)*taucld4
      tauc = taucld1 + taucld2 + taucld3 + taucld4
      tempb7 = asycltb/tauc
      g1b = tempb7
      g2b = tempb7
      g3b = tempb7
      g4b = tempb7
      taucb = -((g1+g2+g3+g4)*tempb7/tauc)
      temp3 = aig_uv(2) + aig_uv(3)*reff_snow
      reff_snowb = (taucld4*temp3+reff_snow*taucld4*aig_uv(3))*g4b
      taucld4b = (aig_uv(1)+temp3*reff_snow)*g4b
      taucld3b = arg_uv(1)*g3b
      temp2 = awg_uv(2) + awg_uv(3)*reff(k, 2)
      reffb(k, 2) = reffb(k, 2) + (taucld2*temp2+reff(k, 2)*taucld2*&
&       awg_uv(3))*g2b
      taucld2b = (awg_uv(1)+temp2*reff(k, 2))*g2b
      temp1 = aig_uv(2) + aig_uv(3)*reff(k, 1)
      reffb(k, 1) = reffb(k, 1) + (taucld1*temp1+reff(k, 1)*taucld1*&
&       aig_uv(3))*g1b
      taucld1b = (aig_uv(1)+temp1*reff(k, 1))*g1b
    ELSE
      reff_snowb = 0.0
      taucld1b = 0.0
      taucld2b = 0.0
      taucld3b = 0.0
      taucld4b = 0.0
      taucb = 0.0
    END IF
    CALL POPREAL4(tauc)
    taucld1b = taucld1b + taucb
    taucld2b = taucld2b + taucb
    taucld3b = taucld3b + taucb
    taucld4b = taucld4b + taucb
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      taucld4b = taucld4b + xai*taudiffb(k, 4)
      xaib = taucld4*taudiffb(k, 4)
      taudiffb(k, 4) = 0.0
      taucld3b = taucld3b + xai*taudiffb(k, 3)
      xaib = xaib + taucld3*taudiffb(k, 3)
      taudiffb(k, 3) = 0.0
      taucld2b = taucld2b + xai*taudiffb(k, 2)
      xaib = xaib + taucld2*taudiffb(k, 2)
      taudiffb(k, 2) = 0.0
      taucld1b = taucld1b + xai*taudiffb(k, 1)
      xaib = xaib + taucld1*taudiffb(k, 1)
      taudiffb(k, 1) = 0.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0
      tempb5 = .5*fa*xaib
      fab = (.5*(caif(it, ia+1)*(fa+1.)-caif(it, ia-1)*(1.-fa))-caif(it&
&       , ia)*2*fa)*xaib + (caif(it, ia-1)+caif(it, ia+1))*tempb5
      CALL POPREAL4(xai)
      tempb6 = .5*ft*xaib
      ftb = (.5*(caif(it+1, ia)*(ft+1.)-caif(it-1, ia)*(1.-ft))-caif(it&
&       , ia)*2*ft)*xaib + (caif(it-1, ia)+caif(it+1, ia))*tempb6
      taucld4b = taucld4b + xai*taubeamb(k, 4)
      xaib = taucld4*taubeamb(k, 4)
      taubeamb(k, 4) = 0.0
      taucld3b = taucld3b + xai*taubeamb(k, 3)
      xaib = xaib + taucld3*taubeamb(k, 3)
      taubeamb(k, 3) = 0.0
      taucld2b = taucld2b + xai*taubeamb(k, 2)
      xaib = xaib + taucld2*taubeamb(k, 2)
      taubeamb(k, 2) = 0.0
      taucld1b = taucld1b + xai*taubeamb(k, 1)
      xaib = xaib + taucld1*taubeamb(k, 1)
      taubeamb(k, 1) = 0.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0
      tempb3 = .5*fa*xaib
      fab = fab + (.5*(caib(im, it, ia+1)*(fa+1.)-caib(im, it, ia-1)*(1.&
&       -fa))-caib(im, it, ia)*2*fa)*xaib + (caib(im, it, ia-1)+caib(im&
&       , it, ia+1))*tempb3
      tempb4 = .5*ft*xaib
      ftb = ftb + (.5*(caib(im, it+1, ia)*(ft+1.)-caib(im, it-1, ia)*(1.&
&       -ft))-caib(im, it, ia)*2*ft)*xaib + (caib(im, it-1, ia)+caib(im&
&       , it+1, ia))*tempb4
      CALL POPREAL4(xai)
      CALL POPINTEGER4(ia)
      CALL POPINTEGER4(it)
      CALL POPINTEGER4(im)
      fab = fab/da
      CALL POPREAL4(ft)
      taucb = ftb/(dt*tauc*LOG(10.0))
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) taucb = 0.0
      CALL POPREAL4(fa)
      tempb2 = fab/cc(kk)
      fcldb(k) = fcldb(k) + tempb2
      ccb(kk) = ccb(kk) - fcld(k)*tempb2/cc(kk)
    ELSE IF (branch .EQ. 1) THEN
      taucb = 0.0
    ELSE
      taucld4b = taucld4b + taubeamb(k, 4) + taudiffb(k, 4)
      taudiffb(k, 4) = 0.0
      taubeamb(k, 4) = 0.0
      taucld3b = taucld3b + taubeamb(k, 3) + taudiffb(k, 3)
      taudiffb(k, 3) = 0.0
      taubeamb(k, 3) = 0.0
      taucld2b = taucld2b + taubeamb(k, 2) + taudiffb(k, 2)
      taudiffb(k, 2) = 0.0
      taubeamb(k, 2) = 0.0
      taucld1b = taucld1b + taubeamb(k, 1) + taudiffb(k, 1)
      taudiffb(k, 1) = 0.0
      taubeamb(k, 1) = 0.0
      GOTO 100
    END IF
    taucld1b = taucld1b + taucb
    taucld2b = taucld2b + taucb
    taucld3b = taucld3b + taucb
    taucld4b = taucld4b + taucb
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER4(kk)
    ELSE IF (branch .EQ. 1) THEN
      CALL POPINTEGER4(kk)
    ELSE
      CALL POPINTEGER4(kk)
    END IF
 100 CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(taucld4)
    ELSE
      CALL POPREAL4(taucld4)
      tempb1 = dp(k)*aib_uv*1.0e3*taucld4b/(cons_grav*reff_snow)
      hydrometsb(k, 4) = hydrometsb(k, 4) + tempb1
      reff_snowb = reff_snowb - hydromets(k, 4)*tempb1/reff_snow
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(reff_snow)
    ELSE
      CALL POPREAL4(reff_snow)
      reffb(k, 4) = reffb(k, 4) + reff_snowb
    END IF
    CALL POPREAL4(taucld3)
    hydrometsb(k, 3) = hydrometsb(k, 3) + dp(k)*arb_uv(1)*1.0e3*taucld3b&
&     /cons_grav
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(taucld2)
      temp0 = awb_uv(2)/reff(k, 2)
      tempb0 = dp(k)*1.0e3*taucld2b
      hydrometsb(k, 2) = hydrometsb(k, 2) + (awb_uv(1)+temp0)*tempb0/&
&       cons_grav
      reffb(k, 2) = reffb(k, 2) - hydromets(k, 2)*temp0*tempb0/(reff(k, &
&       2)*cons_grav)
    ELSE
      CALL POPREAL4(taucld2)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(taucld1)
    ELSE
      CALL POPREAL4(taucld1)
      temp = cons_grav*reff(k, 1)
      tempb = dp(k)*aib_uv*1.0e3*taucld1b/temp
      hydrometsb(k, 1) = hydrometsb(k, 1) + tempb
      reffb(k, 1) = reffb(k, 1) - hydromets(k, 1)*cons_grav*tempb/temp
    END IF
  END DO
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) THEN
    DO k=nlevs,icb,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(3)
        ccb(3) = 0.0
      END IF
    END DO
    DO k=icb-1,ict,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(2)
        ccb(2) = 0.0
      END IF
    END DO
    DO k=ict-1,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(1)
        ccb(1) = 0.0
      END IF
    END DO
  END IF
END SUBROUTINE GETVISTAU1_B

!  Differentiation of getnirtau1 in reverse (adjoint) mode:
!   gradient     of useful results: hydromets asycl taudiff fcld
!                ssacl taubeam reff
!   with respect to varying inputs: hydromets asycl fcld ssacl
!                reff
SUBROUTINE GETNIRTAU1_B(ib, nlevs, cosz, dp, fcld, fcldb, reff, reffb, &
& hydromets, hydrometsb, ict, icb, taubeam, taubeamb, taudiff, taudiffb&
& , asycl, asyclb, ssacl, ssaclb, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv&
& , arb_uv, aib_nir, awb_nir, arb_nir, aia_nir, awa_nir, ara_nir, &
& aig_nir, awg_nir, arg_nir, caib, caif, cons_grav)
  IMPLICIT NONE
! !INPUT PARAMETERS:
!  Band number
  INTEGER, INTENT(IN) :: ib
!  Number of levels
  INTEGER, INTENT(IN) :: nlevs
!  Cosine of solar zenith angle
  REAL, INTENT(IN) :: cosz
!  Delta pressure in Pa
  REAL, INTENT(IN) :: dp(nlevs)
!  Cloud fraction (used sometimes)
  REAL, INTENT(IN) :: fcld(nlevs)
  REAL :: fcldb(nlevs)
!  Effective radius (microns)
  REAL, INTENT(IN) :: reff(nlevs, 4)
  REAL :: reffb(nlevs, 4)
!  Hydrometeors (kg/kg)
  REAL, INTENT(IN) :: hydromets(nlevs, 4)
  REAL :: hydrometsb(nlevs, 4)
!  Flags for various uses 
  INTEGER, INTENT(IN) :: ict, icb
  REAL, INTENT(IN) :: aig_uv(3), awg_uv(3), arg_uv(3)
  REAL, INTENT(IN) :: aib_uv, awb_uv(2), arb_uv(2)
  REAL, INTENT(IN) :: aib_nir, awb_nir(3, 2), arb_nir(3, 2)
  REAL, INTENT(IN) :: aia_nir(3, 3), awa_nir(3, 3), ara_nir(3, 3)
  REAL, INTENT(IN) :: aig_nir(3, 3), awg_nir(3, 3), arg_nir(3, 3)
  REAL, INTENT(IN) :: caib(11, 9, 11), caif(9, 11)
  REAL, INTENT(IN) :: cons_grav
! !OUTPUT PARAMETERS:
!  Optical depth for beam radiation
  REAL :: taubeam(nlevs, 4)
  REAL :: taubeamb(nlevs, 4)
!  Optical depth for diffuse radiation
  REAL :: taudiff(nlevs, 4)
  REAL :: taudiffb(nlevs, 4)
!  Cloud single scattering albedo
  REAL :: ssacl(nlevs)
  REAL :: ssaclb(nlevs)
!  Cloud asymmetry factor
  REAL :: asycl(nlevs)
  REAL :: asyclb(nlevs)
  INTEGER :: k, in, im, it, ia, kk
  REAL :: fm, ft, fa, xai, tauc, asyclt, ssaclt
  REAL :: ftb, fab, xaib, taucb, asycltb, ssacltb
  REAL :: cc(3)
  REAL :: ccb(3)
  REAL :: taucld1, taucld2, taucld3, taucld4
  REAL :: taucld1b, taucld2b, taucld3b, taucld4b
  REAL :: g1, g2, g3, g4
  REAL :: g1b, g2b, g3b, g4b
  REAL :: w1, w2, w3, w4
  REAL :: w1b, w2b, w3b, w4b
  REAL :: reff_snow
  REAL :: reff_snowb
  INTEGER, PARAMETER :: nm=11, nt=9, na=11
  REAL, PARAMETER :: dm=0.1, dt=0.30103, da=0.1, t1=-0.9031
  INTRINSIC MAX
  INTRINSIC MIN
  INTRINSIC LOG10
  INTRINSIC INT
  INTRINSIC REAL
  INTEGER :: branch
  REAL :: temp3
  REAL :: temp2
  REAL :: temp1
  REAL :: temp0
  REAL :: tempb9
  REAL :: tempb8
  REAL :: tempb7
  REAL :: tempb6
  REAL :: tempb5
  REAL :: tempb4
  REAL :: tempb3
  REAL :: tempb2
  REAL :: tempb1
  REAL :: tempb0
  REAL :: tempb
  REAL :: temp
  REAL :: temp6
  REAL :: temp5
  REAL :: temp4
  IF (ict .NE. 0) THEN
!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
    cc = 0.0
    DO k=1,ict-1
      IF (cc(1) .LT. fcld(k)) THEN
        cc(1) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(1) = cc(1)
      END IF
    END DO
    DO k=ict,icb-1
      IF (cc(2) .LT. fcld(k)) THEN
        cc(2) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(2) = cc(2)
      END IF
    END DO
    DO k=icb,nlevs
      IF (cc(3) .LT. fcld(k)) THEN
        cc(3) = fcld(k)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
        cc(3) = cc(3)
      END IF
    END DO
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
  END IF
!-----Compute cloud optical thickness.  Eqs. (4.6) and (4.10)
!     taucld1 is the optical thickness for ice particles
!     taucld2 is the optical thickness for liquid particles
!     taucld3 is the optical thickness for rain drops
!     taucld4 is the optical thickness for snow
  DO k=1,nlevs
    IF (reff(k, 1) .LE. 0.) THEN
      CALL PUSHREAL4(taucld1)
      taucld1 = 0.
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL4(taucld1)
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*aib_nir/reff(k, 1)
      CALL PUSHCONTROL1B(1)
    END IF
    IF (reff(k, 2) .LE. 0.) THEN
      CALL PUSHREAL4(taucld2)
      taucld2 = 0.
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL4(taucld2)
      taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_nir(ib, 1)+&
&       awb_nir(ib, 2)/reff(k, 2))
      CALL PUSHCONTROL1B(0)
    END IF
    CALL PUSHREAL4(taucld3)
    taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_nir(ib, 1)
    IF (reff(k, 4) .GT. 112.0) THEN
      CALL PUSHREAL4(reff_snow)
      reff_snow = 112.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL4(reff_snow)
      reff_snow = reff(k, 4)
      CALL PUSHCONTROL1B(1)
    END IF
    IF (reff_snow .LE. 0.) THEN
      CALL PUSHREAL4(taucld4)
      taucld4 = 0.
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHREAL4(taucld4)
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*aib_nir/reff_snow
      CALL PUSHCONTROL1B(1)
    END IF
    IF (ict .NE. 0) THEN
!-----scale cloud optical thickness in each layer from taucld (with
!     cloud amount fcld) to taubeam and taudiff (with cloud amount cc).
!     taubeam is the scaled optical thickness for beam radiation and
!     taudiff is for diffuse radiation (see section 7).
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
      IF (k .LT. ict) THEN
        CALL PUSHINTEGER4(kk)
        kk = 1
        CALL PUSHCONTROL2B(0)
      ELSE IF (k .GE. ict .AND. k .LT. icb) THEN
        CALL PUSHINTEGER4(kk)
        kk = 2
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHINTEGER4(kk)
        kk = 3
        CALL PUSHCONTROL2B(2)
      END IF
      tauc = taucld1 + taucld2 + taucld3 + taucld4
      IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
!-----normalize cloud cover following Eq. (7.8)
        IF (cc(kk) .NE. 0.0) THEN
          CALL PUSHREAL4(fa)
          fa = fcld(k)/cc(kk)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHREAL4(fa)
          fa = 0.0
          CALL PUSHCONTROL1B(0)
        END IF
        IF (tauc .GT. 32.) THEN
          tauc = 32.
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          tauc = tauc
        END IF
        fm = cosz/dm
        CALL PUSHREAL4(ft)
        ft = (LOG10(tauc)-t1)/dt
        fa = fa/da
        CALL PUSHINTEGER4(im)
        im = INT(fm + 1.5)
        CALL PUSHINTEGER4(it)
        it = INT(ft + 1.5)
        CALL PUSHINTEGER4(ia)
        ia = INT(fa + 1.5)
        IF (im .LT. 2) THEN
          im = 2
        ELSE
          im = im
        END IF
        IF (it .LT. 2) THEN
          it = 2
        ELSE
          it = it
        END IF
        IF (ia .LT. 2) THEN
          ia = 2
        ELSE
          ia = ia
        END IF
        IF (im .GT. nm - 1) THEN
          im = nm - 1
        ELSE
          im = im
        END IF
        IF (it .GT. nt - 1) THEN
          it = nt - 1
        ELSE
          it = it
        END IF
        IF (ia .GT. na - 1) THEN
          ia = na - 1
        ELSE
          ia = ia
        END IF
        fm = fm - REAL(im - 1)
        ft = ft - REAL(it - 1)
        fa = fa - REAL(ia - 1)
!-----scale cloud optical thickness for beam radiation following 
!     Eq. (7.3).
!     the scaling factor, xai, is a function of the solar zenith
!     angle, optical thickness, and cloud cover.
        CALL PUSHREAL4(xai)
        xai = (-(caib(im-1, it, ia)*(1.-fm))+caib(im+1, it, ia)*(1.+fm))&
&         *fm*.5 + caib(im, it, ia)*(1.-fm*fm)
        xai = xai + (-(caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(&
&         1.+ft))*ft*.5 + caib(im, it, ia)*(1.-ft*ft)
        xai = xai + (-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(&
&         1.+fa))*fa*.5 + caib(im, it, ia)*(1.-fa*fa)
        xai = xai - 2.*caib(im, it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
        CALL PUSHREAL4(xai)
        xai = (-(caif(it-1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ft*.5 +&
&         caif(it, ia)*(1.-ft*ft)
        xai = xai + (-(caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*&
&         fa*.5 + caif(it, ia)*(1.-fa*fa)
        xai = xai - caif(it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
          xai = xai
        END IF
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
!-----compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
    CALL PUSHREAL4(tauc)
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      CALL PUSHREAL4(w1)
      w1 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff(k, 1)&
&       )*reff(k, 1)))*taucld1
      CALL PUSHREAL4(w2)
      w2 = (1.-(awa_nir(ib, 1)+(awa_nir(ib, 2)+awa_nir(ib, 3)*reff(k, 2)&
&       )*reff(k, 2)))*taucld2
      CALL PUSHREAL4(w3)
      w3 = (1.-ara_nir(ib, 1))*taucld3
      CALL PUSHREAL4(w4)
      w4 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff_snow)&
&       *reff_snow))*taucld4
      g1 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 1))*&
&       reff(k, 1))*w1
      g2 = (awg_nir(ib, 1)+(awg_nir(ib, 2)+awg_nir(ib, 3)*reff(k, 2))*&
&       reff(k, 2))*w2
      g3 = arg_nir(ib, 1)*w3
      g4 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 4))*&
&       reff(k, 4))*w4
      IF (w1 + w2 + w3 + w4 .NE. 0.0) THEN
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
  END DO
  ccb = 0.0
  DO k=nlevs,1,-1
    asycltb = asyclb(k)
    asyclb(k) = 0.0
    ssacltb = ssaclb(k)
    ssaclb(k) = 0.0
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      w1 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff(k, 1)&
&       )*reff(k, 1)))*taucld1
      w2 = (1.-(awa_nir(ib, 1)+(awa_nir(ib, 2)+awa_nir(ib, 3)*reff(k, 2)&
&       )*reff(k, 2)))*taucld2
      taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_nir(ib, 1)
      w3 = (1.-ara_nir(ib, 1))*taucld3
      w4 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff_snow)&
&       *reff_snow))*taucld4
      g1 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 1))*&
&       reff(k, 1))*w1
      g2 = (awg_nir(ib, 1)+(awg_nir(ib, 2)+awg_nir(ib, 3)*reff(k, 2))*&
&       reff(k, 2))*w2
      g3 = arg_nir(ib, 1)*w3
      g4 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 4))*&
&       reff(k, 4))*w4
      tempb8 = asycltb/(w1+w2+w3+w4)
      tempb9 = -((g1+g2+g3+g4)*tempb8/(w1+w2+w3+w4))
      g1b = tempb8
      g2b = tempb8
      g3b = tempb8
      g4b = tempb8
      w1b = tempb9
      w2b = tempb9
      w3b = tempb9
      w4b = tempb9
    ELSE IF (branch .EQ. 1) THEN
      w1b = 0.0
      w2b = 0.0
      w3b = 0.0
      w4b = 0.0
      g1b = 0.0
      g2b = 0.0
      g3b = 0.0
      g4b = 0.0
    ELSE
      reff_snowb = 0.0
      taucld1b = 0.0
      taucld2b = 0.0
      taucld3b = 0.0
      taucld4b = 0.0
      taucb = 0.0
      GOTO 100
    END IF
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    tempb7 = ssacltb/tauc
    temp6 = aig_nir(ib, 2) + aig_nir(ib, 3)*reff(k, 4)
    reffb(k, 4) = reffb(k, 4) + (w4*temp6+reff(k, 4)*w4*aig_nir(ib, 3))*&
&     g4b
    w4b = w4b + tempb7 + (aig_nir(ib, 1)+temp6*reff(k, 4))*g4b
    w3b = w3b + tempb7 + arg_nir(ib, 1)*g3b
    temp5 = awg_nir(ib, 2) + awg_nir(ib, 3)*reff(k, 2)
    reffb(k, 2) = reffb(k, 2) + (w2*temp5+reff(k, 2)*w2*awg_nir(ib, 3))*&
&     g2b
    w2b = w2b + tempb7 + (awg_nir(ib, 1)+temp5*reff(k, 2))*g2b
    temp4 = aig_nir(ib, 2) + aig_nir(ib, 3)*reff(k, 1)
    reffb(k, 1) = reffb(k, 1) + (w1*temp4+reff(k, 1)*w1*aig_nir(ib, 3))*&
&     g1b
    w1b = w1b + tempb7 + (aig_nir(ib, 1)+temp4*reff(k, 1))*g1b
    taucb = -((w1+w2+w3+w4)*tempb7/tauc)
    CALL POPREAL4(w4)
    temp3 = aia_nir(ib, 2) + aia_nir(ib, 3)*reff_snow
    reff_snowb = (-(taucld4*temp3)-reff_snow*taucld4*aia_nir(ib, 3))*w4b
    taucld4b = (1.-temp3*reff_snow-aia_nir(ib, 1))*w4b
    CALL POPREAL4(w3)
    taucld3b = (1.-ara_nir(ib, 1))*w3b
    CALL POPREAL4(w2)
    temp2 = awa_nir(ib, 2) + awa_nir(ib, 3)*reff(k, 2)
    reffb(k, 2) = reffb(k, 2) + (-(taucld2*temp2)-reff(k, 2)*taucld2*&
&     awa_nir(ib, 3))*w2b
    taucld2b = (1.-temp2*reff(k, 2)-awa_nir(ib, 1))*w2b
    CALL POPREAL4(w1)
    temp1 = aia_nir(ib, 2) + aia_nir(ib, 3)*reff(k, 1)
    reffb(k, 1) = reffb(k, 1) + (-(taucld1*temp1)-reff(k, 1)*taucld1*&
&     aia_nir(ib, 3))*w1b
    taucld1b = (1.-temp1*reff(k, 1)-aia_nir(ib, 1))*w1b
 100 CALL POPREAL4(tauc)
    taucld1b = taucld1b + taucb
    taucld2b = taucld2b + taucb
    taucld3b = taucld3b + taucb
    taucld4b = taucld4b + taucb
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      taucld4b = taucld4b + xai*taudiffb(k, 4)
      xaib = taucld4*taudiffb(k, 4)
      taudiffb(k, 4) = 0.0
      taucld3b = taucld3b + xai*taudiffb(k, 3)
      xaib = xaib + taucld3*taudiffb(k, 3)
      taudiffb(k, 3) = 0.0
      taucld2b = taucld2b + xai*taudiffb(k, 2)
      xaib = xaib + taucld2*taudiffb(k, 2)
      taudiffb(k, 2) = 0.0
      taucld1b = taucld1b + xai*taudiffb(k, 1)
      xaib = xaib + taucld1*taudiffb(k, 1)
      taudiffb(k, 1) = 0.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0
      tempb5 = .5*fa*xaib
      fab = (.5*(caif(it, ia+1)*(fa+1.)-caif(it, ia-1)*(1.-fa))-caif(it&
&       , ia)*2*fa)*xaib + (caif(it, ia-1)+caif(it, ia+1))*tempb5
      CALL POPREAL4(xai)
      tempb6 = .5*ft*xaib
      ftb = (.5*(caif(it+1, ia)*(ft+1.)-caif(it-1, ia)*(1.-ft))-caif(it&
&       , ia)*2*ft)*xaib + (caif(it-1, ia)+caif(it+1, ia))*tempb6
      taucld4b = taucld4b + xai*taubeamb(k, 4)
      xaib = taucld4*taubeamb(k, 4)
      taubeamb(k, 4) = 0.0
      taucld3b = taucld3b + xai*taubeamb(k, 3)
      xaib = xaib + taucld3*taubeamb(k, 3)
      taubeamb(k, 3) = 0.0
      taucld2b = taucld2b + xai*taubeamb(k, 2)
      xaib = xaib + taucld2*taubeamb(k, 2)
      taubeamb(k, 2) = 0.0
      taucld1b = taucld1b + xai*taubeamb(k, 1)
      xaib = xaib + taucld1*taubeamb(k, 1)
      taubeamb(k, 1) = 0.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) xaib = 0.0
      tempb3 = .5*fa*xaib
      fab = fab + (.5*(caib(im, it, ia+1)*(fa+1.)-caib(im, it, ia-1)*(1.&
&       -fa))-caib(im, it, ia)*2*fa)*xaib + (caib(im, it, ia-1)+caib(im&
&       , it, ia+1))*tempb3
      tempb4 = .5*ft*xaib
      ftb = ftb + (.5*(caib(im, it+1, ia)*(ft+1.)-caib(im, it-1, ia)*(1.&
&       -ft))-caib(im, it, ia)*2*ft)*xaib + (caib(im, it-1, ia)+caib(im&
&       , it+1, ia))*tempb4
      CALL POPREAL4(xai)
      CALL POPINTEGER4(ia)
      CALL POPINTEGER4(it)
      CALL POPINTEGER4(im)
      fab = fab/da
      CALL POPREAL4(ft)
      taucb = ftb/(dt*tauc*LOG(10.0))
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) taucb = 0.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4(fa)
      ELSE
        CALL POPREAL4(fa)
        tempb2 = fab/cc(kk)
        fcldb(k) = fcldb(k) + tempb2
        ccb(kk) = ccb(kk) - fcld(k)*tempb2/cc(kk)
      END IF
    ELSE IF (branch .EQ. 1) THEN
      taucb = 0.0
    ELSE
      taucld4b = taucld4b + taubeamb(k, 4) + taudiffb(k, 4)
      taudiffb(k, 4) = 0.0
      taubeamb(k, 4) = 0.0
      taucld3b = taucld3b + taubeamb(k, 3) + taudiffb(k, 3)
      taudiffb(k, 3) = 0.0
      taubeamb(k, 3) = 0.0
      taucld2b = taucld2b + taubeamb(k, 2) + taudiffb(k, 2)
      taudiffb(k, 2) = 0.0
      taubeamb(k, 2) = 0.0
      taucld1b = taucld1b + taubeamb(k, 1) + taudiffb(k, 1)
      taudiffb(k, 1) = 0.0
      taubeamb(k, 1) = 0.0
      GOTO 110
    END IF
    taucld1b = taucld1b + taucb
    taucld2b = taucld2b + taucb
    taucld3b = taucld3b + taucb
    taucld4b = taucld4b + taucb
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPINTEGER4(kk)
    ELSE IF (branch .EQ. 1) THEN
      CALL POPINTEGER4(kk)
    ELSE
      CALL POPINTEGER4(kk)
    END IF
 110 CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(taucld4)
    ELSE
      CALL POPREAL4(taucld4)
      tempb1 = dp(k)*aib_nir*1.0e3*taucld4b/(cons_grav*reff_snow)
      hydrometsb(k, 4) = hydrometsb(k, 4) + tempb1
      reff_snowb = reff_snowb - hydromets(k, 4)*tempb1/reff_snow
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(reff_snow)
    ELSE
      CALL POPREAL4(reff_snow)
      reffb(k, 4) = reffb(k, 4) + reff_snowb
    END IF
    CALL POPREAL4(taucld3)
    hydrometsb(k, 3) = hydrometsb(k, 3) + dp(k)*1.0e3*arb_nir(ib, 1)*&
&     taucld3b/cons_grav
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(taucld2)
      temp0 = awb_nir(ib, 2)/reff(k, 2)
      tempb0 = dp(k)*1.0e3*taucld2b
      hydrometsb(k, 2) = hydrometsb(k, 2) + (awb_nir(ib, 1)+temp0)*&
&       tempb0/cons_grav
      reffb(k, 2) = reffb(k, 2) - hydromets(k, 2)*temp0*tempb0/(reff(k, &
&       2)*cons_grav)
    ELSE
      CALL POPREAL4(taucld2)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4(taucld1)
    ELSE
      CALL POPREAL4(taucld1)
      temp = cons_grav*reff(k, 1)
      tempb = dp(k)*aib_nir*1.0e3*taucld1b/temp
      hydrometsb(k, 1) = hydrometsb(k, 1) + tempb
      reffb(k, 1) = reffb(k, 1) - hydromets(k, 1)*cons_grav*tempb/temp
    END IF
  END DO
  CALL POPCONTROL1B(branch)
  IF (branch .NE. 0) THEN
    DO k=nlevs,icb,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(3)
        ccb(3) = 0.0
      END IF
    END DO
    DO k=icb-1,ict,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(2)
        ccb(2) = 0.0
      END IF
    END DO
    DO k=ict-1,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        fcldb(k) = fcldb(k) + ccb(1)
        ccb(1) = 0.0
      END IF
    END DO
  END IF
END SUBROUTINE GETNIRTAU1_B

