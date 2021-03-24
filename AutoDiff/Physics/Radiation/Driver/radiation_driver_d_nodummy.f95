SUBROUTINE RADIATION_DRIVER_D(im, jm, lm, runalarm, dt, ut, ptt, pttd, &
& qvt, qvtd, o3t, o3td, cflst, cflstd, cfcnt, cfcntd, qit, qitd, qlt, &
& qltd, ps, ts, emis, delt, cosz, slr, rgbuv, rgfuv, rgbir, rgfir, sc, &
& du001t, du001td, du002t, du002td, du003t, du003td, du004t, du004td, &
& du005t, du005td, taua_ir_c, ssaa_ir_c, asya_ir_c, taua_so_c, ssaa_so_c&
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
  REAL, DIMENSION(im, jm, lm), INTENT(IN) :: qvtd, o3td, cflstd, cfcntd&
& , qitd, qltd
  REAL, DIMENSION(im, jm), INTENT(IN) :: ps, ts
  REAL, DIMENSION(im, jm), INTENT(IN) :: emis, delt, cosz, slr, rgbuv, &
& rgfuv, rgbir, rgfir
  REAL, DIMENSION(im, jm, lm), INTENT(IN) :: du001t, du002t, du003t, &
& du004t, du005t
  REAL, DIMENSION(im, jm, lm), INTENT(IN) :: du001td, du002td, du003td, &
& du004td, du005td
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
  REAL :: hk_ir_tempd(3, 10)
  REAL, INTENT(IN) :: hk_uv_old(5), hk_ir_old(3, 10)
!Inouts
  REAL, DIMENSION(im, jm, lm), INTENT(INOUT) :: ptt
  REAL, DIMENSION(im, jm, lm), INTENT(INOUT) :: pttd
!Locals
  INTEGER :: i, j, k, l
  REAL :: co2
  INTEGER :: levs925, lcldmh, lcldlm
  LOGICAL, PARAMETER :: trace=.true., overcast=.false.
  INTEGER, PARAMETER :: ns=1
  REAL, DIMENSION(im, jm, lm, 10) :: taua_irt, ssaa_irt, asya_irt
  REAL, DIMENSION(im, jm, lm, 10) :: taua_irtd, ssaa_irtd, asya_irtd
  REAL, DIMENSION(im, jm, lm, 8) :: taua_sot, ssaa_sot, asya_sot
  REAL, DIMENSION(im, jm, lm, 8) :: taua_sotd, ssaa_sotd, asya_sotd
  REAL, DIMENSION(im, jm, lm, na) :: aerot, saerot
  REAL, DIMENSION(im, jm, lm, na) :: aerotd, saerotd
  CHARACTER(len=50) :: aerosols(na)
  REAL :: x
  INTEGER :: ai
  REAL, DIMENSION(im, jm, lm) :: ptt1, tt, ple, pleso, plof, pif, deltap&
& , tt_out, dtdt
  REAL, DIMENSION(im, jm, lm) :: ptt1d, ttd, pled, plesod, tt_outd, &
& dtdtd
  REAL, DIMENSION(im, jm, lm) :: o3mmt
  REAL, DIMENSION(im, jm, lm) :: o3mmtd
  REAL, DIMENSION(im, jm, lm) :: n2ot, ch4t, cfc11t, cfc12t, cfc22t
  REAL, DIMENSION(im, jm, lm) :: rad_qlt, rad_qit, rad_cft, rad_rlt, &
& rad_rit
  REAL, DIMENSION(im, jm, lm) :: rad_qltd, rad_qitd, rad_cftd, rad_rltd&
& , rad_ritd
  REAL, DIMENSION(im, jm, lm, 4) :: cwct, refft
  REAL, DIMENSION(im, jm, lm, 4) :: cwctd, refftd
  REAL, DIMENSION(im, jm, 1) :: fst, tgt, tvt
  REAL, DIMENSION(im, jm, 1) :: tgtd, tvtd
  REAL, DIMENSION(im, jm, 1, 10) :: egt, evt, rvt
  REAL, DIMENSION(im, jm, 1, 10) :: egtd
  REAL, DIMENSION(im, jm, lm) :: flwut, flwdt, flxnt, flwt_int, fswt_int
  REAL, DIMENSION(im, jm, lm) :: flwutd, flwdtd, flxntd, flwt_intd, &
& fswt_intd
  REAL, DIMENSION(im, jm, lm) :: flwavet, fswavet, dfdtst_int
  REAL, DIMENSION(im, jm, lm) :: flwavetd, fswavetd, dfdtst_intd
  REAL, DIMENSION(im, jm, lm) :: qilst, qllst, qicnt, qlcnt
  REAL, DIMENSION(im, jm, lm) :: qilstd, qllstd, qicntd, qlcntd
  REAL, DIMENSION(im, jm) :: t2mt, tempor, cosztmp
  REAL, DIMENSION(im, jm) :: t2mtd, cosztmpd
  REAL, DIMENSION(0:lm) :: ak, bk, pref
  INTRINSIC MAX
  INTRINSIC SUM
  REAL, DIMENSION(im, jm, lm) :: pwx1
  REAL :: pwy1
  REAL :: pwr1
  REAL, DIMENSION(im, jm) :: pwx10
  REAL, DIMENSION(im, jm) :: pwr10
!Compute pref from ak and bk
  ak = 1.0
  bk = 1.0
  DO l=0,lm
    pref(l) = ak(l+1) + bk(l+1)*mapl_p00
  END DO
!Pressure at the half levels from Ps
  DO l=0,lm
    pled(:, :, l) = 0.0
    ple(:, :, l) = ak(l+1) + bk(l+1)*ps(:, :)
  END DO
  deltap = ple(:, :, 1:lm) - ple(:, :, 0:lm-1)
!Pressure (hPa) and Exner pressure at the full levels
  plof(:, :, 1:lm) = 0.01*0.5*(ple(:, :, 0:lm-1)+ple(:, :, 1:lm))
  pwx1 = plof(:, :, 1:lm)/1000.
  pwy1 = mapl_rgas/mapl_cp
  pif(:, :, 1:lm) = pwx1**pwy1
!Potential temperature with p0=10000
!Place in holder so as not to overwrite
  pwr1 = mapl_p00**mapl_kappa
  ptt1d = pwr1*pttd
  ptt1 = ptt*pwr1
!Temperature
!Save input for later.
  ttd = pif*ptt1d
  tt = ptt1*pif
  IF (runalarm .EQ. 1) THEN
!2m temperature
    pwx10 = 0.5*(1.0+ple(:, :, lm-1)/ple(:, :, lm))
    pwy1 = -mapl_kappa
    pwr10 = pwx10**pwy1
    t2mtd = pwr10*ttd(:, :, lm)
    t2mt = tt(:, :, lm)*pwr10
!Pressure in mb for SORAD
    plesod(:, :, 1:lm+1) = 0.0
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
      WHERE (ut(:, :, l) .GT. 4.) tempor(:, :) = 1.
    END DO
!Convert ozone from ppmv to kgkg 
    o3mmtd = o3td/1.e6
    o3mmt = o3t/1.e6
    co2 = 1.0
!Initialize other gases, not available right now.
    n2ot = 0.0
    ch4t = 0.0
    cfc11t = 0.0
    cfc12t = 0.0
    cfc22t = 0.0
    qilstd = ilsf*qitd
    qilst = qit*ilsf
    qicntd = icnf*qitd
    qicnt = qit*icnf
    qllstd = llsf*qltd
    qllst = qlt*llsf
    qlcntd = lcnf*qltd
    qlcnt = qlt*lcnf
!Initialize RAD/CLOUD variables
    rad_cft = 0.0
    rad_qlt = 0.0
    rad_qit = 0.0
    rad_rlt = 0.0
    rad_rit = 0.0
    rad_cftd = 0.0
    rad_rltd = 0.0
    rad_ritd = 0.0
    rad_qltd = 0.0
    rad_qitd = 0.0
!Call RADcouple to produce cloud-rad variables
    DO i=1,im
      DO j=1,jm
        DO l=1,lm
          CALL RADCOUPLE_D(tt(i, j, l), ttd(i, j, l), plof(i, j, l), &
&                    cflst(i, j, l), cflstd(i, j, l), cfcnt(i, j, l), &
&                    cfcntd(i, j, l), qllst(i, j, l), qllstd(i, j, l), &
&                    qilst(i, j, l), qilstd(i, j, l), qlcnt(i, j, l), &
&                    qlcntd(i, j, l), qicnt(i, j, l), qicntd(i, j, l), &
&                    rad_qlt(i, j, l), rad_qltd(i, j, l), rad_qit(i, j, &
&                    l), rad_qitd(i, j, l), rad_cft(i, j, l), rad_cftd(i&
&                    , j, l), rad_rlt(i, j, l), rad_rltd(i, j, l), &
&                    rad_rit(i, j, l), rad_ritd(i, j, l), tempor(i, j))
        END DO
      END DO
    END DO
!QI
    cwctd = 0.0
    cwctd(:, :, :, 1) = rad_qitd
    cwct(:, :, :, 1) = rad_qit
!QL
    cwctd(:, :, :, 2) = rad_qltd
    cwct(:, :, :, 2) = rad_qlt
!QR - this is not available right now
    cwctd(:, :, :, 3) = 0.0
    cwct(:, :, :, 3) = 0.0
!QI - this is not available right now
    cwctd(:, :, :, 4) = 0.0
    cwct(:, :, :, 4) = 0.0
!RI
    refftd = 0.0
    refftd(:, :, :, 1) = 1.0e6*rad_ritd
    refft(:, :, :, 1) = rad_rit*1.0e6
!RL
    refftd(:, :, :, 2) = 1.0e6*rad_rltd
    refft(:, :, :, 2) = rad_rlt*1.0e6
!RR
    refftd(:, :, :, 3) = 0.0
    refft(:, :, :, 3) = 100.e-6*1.0e6
!RS
    refftd(:, :, :, 4) = 0.0
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
    fst = 0.0
    tgt = 0.0
    tvt = 0.0
    evt = 0.0
    rvt = 0.0
    DO l=1,10
      egtd(:, :, 1, l) = 0.0
      egt(:, :, 1, l) = emis(:, :)
    END DO
    fst = 1.0
    tgtd(:, :, 1) = 0.0
    tgt(:, :, 1) = ts
    tvtd(:, :, 1) = 0.0
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
      ssaa_irtd = 0.0
      taua_sotd = 0.0
      ssaa_sotd = 0.0
      asya_irtd = 0.0
      asya_sotd = 0.0
      taua_irtd = 0.0
    ELSE IF (na .GT. 0) THEN
!Place aerosol fields into single rank 4 array that can be looped over
      aerotd = 0.0
      aerotd(:, :, :, 1) = du001td
      aerot(:, :, :, 1) = du001t
      aerotd(:, :, :, 2) = du002td
      aerot(:, :, :, 2) = du002t
      aerotd(:, :, :, 3) = du003td
      aerot(:, :, :, 3) = du003t
      aerotd(:, :, :, 4) = du004td
      aerot(:, :, :, 4) = du004t
      aerotd(:, :, :, 5) = du005td
      aerot(:, :, :, 5) = du005t
!Names of aerosols as named in NL model
      aerosols(1) = 'du001'
      aerosols(2) = 'du002'
      aerosols(3) = 'du003'
      aerosols(4) = 'du004'
      aerosols(5) = 'du005'
      saerotd = 0.0
!Optical property calculation needs pressure as (1/g)dp/dz*q
      DO l=1,lm
        DO j=1,jm
          DO i=1,im
            x = (ple(i, j, l)-ple(i, j, l-1))*0.01*(100./mapl_grav)
            DO ai=1,na
              saerotd(i, j, l, ai) = x*aerotd(i, j, l, ai)
              saerot(i, j, l, ai) = x*aerot(i, j, l, ai)
            END DO
          END DO
        END DO
      END DO
      taua_irt = 0.0
      ssaa_irt = 0.0
      asya_irt = 0.0
      ssaa_irtd = 0.0
      asya_irtd = 0.0
      taua_irtd = 0.0
      DO j=1,nbchou_ir
        DO l=1,na
          taua_irtd(:, :, :, j) = taua_irtd(:, :, :, j) + taua_ir_c(:, :&
&           , :, j, l)*saerotd(:, :, :, l)
          taua_irt(:, :, :, j) = taua_irt(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_ir_c(:, :, :, j, l)
          ssaa_irtd(:, :, :, j) = ssaa_irtd(:, :, :, j) + taua_ir_c(:, :&
&           , :, j, l)*ssaa_ir_c(:, :, :, j, l)*saerotd(:, :, :, l)
          ssaa_irt(:, :, :, j) = ssaa_irt(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_ir_c(:, :, :, j, l)*ssaa_ir_c(:, :, :, j, l)
          asya_irtd(:, :, :, j) = asya_irtd(:, :, :, j) + taua_ir_c(:, :&
&           , :, j, l)*ssaa_ir_c(:, :, :, j, l)*asya_ir_c(:, :, :, j, l)&
&           *saerotd(:, :, :, l)
          asya_irt(:, :, :, j) = asya_irt(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_ir_c(:, :, :, j, l)*ssaa_ir_c(:, :, :, j, l)*&
&           asya_ir_c(:, :, :, j, l)
        END DO
      END DO
      taua_sot = 0.0
      ssaa_sot = 0.0
      asya_sot = 0.0
      taua_sotd = 0.0
      ssaa_sotd = 0.0
      asya_sotd = 0.0
      DO j=1,nbchou_so
        DO l=1,na
          taua_sotd(:, :, :, j) = taua_sotd(:, :, :, j) + taua_so_c(:, :&
&           , :, j, l)*saerotd(:, :, :, l)
          taua_sot(:, :, :, j) = taua_sot(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_so_c(:, :, :, j, l)
          ssaa_sotd(:, :, :, j) = ssaa_sotd(:, :, :, j) + taua_so_c(:, :&
&           , :, j, l)*ssaa_so_c(:, :, :, j, l)*saerotd(:, :, :, l)
          ssaa_sot(:, :, :, j) = ssaa_sot(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_so_c(:, :, :, j, l)*ssaa_so_c(:, :, :, j, l)
          asya_sotd(:, :, :, j) = asya_sotd(:, :, :, j) + taua_so_c(:, :&
&           , :, j, l)*ssaa_so_c(:, :, :, j, l)*asya_so_c(:, :, :, j, l)&
&           *saerotd(:, :, :, l)
          asya_sot(:, :, :, j) = asya_sot(:, :, :, j) + saerot(:, :, :, &
&           l)*taua_so_c(:, :, :, j, l)*ssaa_so_c(:, :, :, j, l)*&
&           asya_so_c(:, :, :, j, l)
        END DO
      END DO
    ELSE
      ssaa_irtd = 0.0
      taua_sotd = 0.0
      ssaa_sotd = 0.0
      asya_irtd = 0.0
      asya_sotd = 0.0
      taua_irtd = 0.0
    END IF
    IF (sc .LT. 0.0) THEN
      hk_uv_temp = hk(1:5)
      DO l=1,3
        hk_ir_tempd(l, :) = 0.0
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
          cosztmpd(i, j) = 0.0
          cosztmp(i, j) = cosz(i, j)
        ELSE
          cosztmpd(i, j) = 0.0
          cosztmp(i, j) = .0001
        END IF
      END DO
    END DO
!Initialize the fluxes
    flwut = 0.0
    flwdt = 0.0
    flxnt = 0.0
!SOLAR RADIATION (SHORT WAVE)
    CALL SORAD_D(im*jm, lm, nbchou_so, cosztmp, pleso, tt, ttd, qvt, &
&          qvtd, o3mmt, o3mmtd, co2, cwct, cwctd, rad_cft, rad_cftd, &
&          lcldmh, lcldlm, refft, refftd, hk_uv_temp, hk_ir_temp, &
&          taua_sot, taua_sotd, ssaa_sot, ssaa_sotd, asya_sot, asya_sotd&
&          , rgbuv, rgfuv, rgbir, rgfir, flxnt, flxntd, mapl_grav, wk_uv&
&          , zk_uv, ry_uv, xk_ir, ry_ir, cah, coa, aig_uv, awg_uv, &
&          arg_uv, aib_uv, awb_uv, arb_uv, aib_nir, awb_nir, arb_nir, &
&          aia_nir, awa_nir, ara_nir, aig_nir, awg_nir, arg_nir, caib, &
&          caif)
    fswt_intd = flxntd
    fswt_int = flxnt
    dfdtst_intd = 0.0
    flwdtd = 0.0
    flwutd = 0.0
    DO i=1,im
      DO j=1,jm
!INFRARED RADIATION (LONG WAVE)
        CALL IRRAD_D(1, lm, ple(i, j, :), tt(i, j, :), ttd(i, j, :), qvt&
&              (i, j, :), qvtd(i, j, :), o3mmt(i, j, :), o3mmtd(i, j, :)&
&              , t2mt(i, j), t2mtd(i, j), co2, trace, n2ot(i, j, :), &
&              ch4t(i, j, :), cfc11t(i, j, :), cfc12t(i, j, :), cfc22t(i&
&              , j, :), cwct(i, j, :, :), cwctd(i, j, :, :), rad_cft(i, &
&              j, :), rad_cftd(i, j, :), lcldmh, lcldlm, refft(i, j, :, &
&              :), refftd(i, j, :, :), ns, fst(i, j, :), tgt(i, j, :), &
&              egt(i, j, :, :), tvt(i, j, :), evt(i, j, :, :), rvt(i, j&
&              , :, :), na, nbchou_ir, taua_irt(i, j, :, :), taua_irtd(i&
&              , j, :, :), ssaa_irt(i, j, :, :), ssaa_irtd(i, j, :, :), &
&              asya_irt(i, j, :, :), asya_irtd(i, j, :, :), flwut(i, j, &
&              :), flwutd(i, j, :), flwdt(i, j, :), flwdtd(i, j, :), &
&              dfdtst_int(i, j, :), dfdtst_intd(i, j, :), overcast, &
&              aib_ir, awb_ir, aiw_ir, aww_ir, aig_ir, awg_ir, xkw, xke&
&              , mw, aw, bw, pm, fkw, gkw, cb, dcb, w11, w12, w13, p11, &
&              p12, p13, dwe, dpe, c1, c2, c3, oo1, oo2, oo3, h11, h12, &
&              h13, h21, h22, h23, h81, h82, h83)
      END DO
    END DO
    flwt_intd = flwutd + flwdtd
    flwt_int = flwut + flwdt
    flwavetd = 0.0
    fswavetd = 0.0
  ELSE
    dfdtst_intd = 0.0
    flwavetd = 0.0
    fswt_intd = 0.0
    flwt_intd = 0.0
    fswavetd = 0.0
  END IF
  DO l=0,lm
    flwavetd(:, :, l) = flwt_intd(:, :, l) + delt*dfdtst_intd(:, :, l)
    flwavet(:, :, l) = flwt_int(:, :, l) + dfdtst_int(:, :, l)*delt
    fswavetd(:, :, l) = slr*fswt_intd(:, :, l)
    fswavet(:, :, l) = fswt_int(:, :, l)*slr
  END DO
!Tangent of fluxes to DTDt (mass-weighted)
  dtdtd = mapl_grav*(flwavetd(:, :, 0:lm-1)-flwavetd(:, :, 1:lm)+&
&   fswavetd(:, :, 0:lm-1)-fswavetd(:, :, 1:lm))/mapl_cp
  dtdt = (flwavet(:, :, 0:lm-1)-flwavet(:, :, 1:lm)+(fswavet(:, :, 0:lm-&
&   1)-fswavet(:, :, 1:lm)))*(mapl_grav/mapl_cp)
!Tangent of mass-weighted DTDt to updated temperature
  tt_outd = ttd + dt*dtdtd/deltap
  tt_out = tt + dt*dtdt/deltap
!Tangent of temperature to potential temperature
  ptt1d = tt_outd/pif
  ptt1 = tt_out/pif
!Tangent of potential temperature to P0 = 1
  pwr1 = mapl_p00**mapl_kappa
  pttd = ptt1d/pwr1
  ptt = ptt1/pwr1
END SUBROUTINE RADIATION_DRIVER_D

!  Differentiation of radcouple in forward (tangent) mode:
!   variations   of useful results: rad_cf rad_qi rad_ql rad_ri
!                rad_rl
!   with respect to varying inputs: qcian af qclls qcils cf qclan
!                te
SUBROUTINE RADCOUPLE_D(te, ted, pl, cf, cfd, af, afd, qclls, qcllsd, &
& qcils, qcilsd, qclan, qcland, qcian, qciand, rad_ql, rad_qld, rad_qi, &
& rad_qid, rad_cf, rad_cfd, rad_rl, rad_rld, rad_ri, rad_rid, tempor)
  IMPLICIT NONE
!Inputs
  REAL, INTENT(IN) :: te, pl, tempor
  REAL, INTENT(IN) :: ted
  REAL, INTENT(IN) :: af, cf, qclan, qcian, qclls, qcils
  REAL, INTENT(IN) :: afd, cfd, qcland, qciand, qcllsd, qcilsd
! real, intent(in ) :: QRN_ALL, QSN_ALL
!Outputs
  REAL, INTENT(OUT) :: rad_ql, rad_qi, rad_cf, rad_rl, rad_ri
  REAL, INTENT(OUT) :: rad_qld, rad_qid, rad_cfd, rad_rld, rad_rid
! real, intent(out) :: RAD_QR,RAD_QS
!Locals
  REAL :: ss, rad_ri_an, afx, alph
  REAL :: ssd, rad_ri_and, afxd
  REAL, PARAMETER :: min_ri=20.e-6, max_ri=40.e-6, ri_anv=30.e-6
  INTRINSIC MIN
  INTRINSIC MAX
  REAL :: max2d
  REAL :: min3
  REAL :: min2
  REAL :: min1
  REAL :: min1d
  REAL :: x2
  REAL :: x2d
  REAL :: x1
  REAL :: max3d
  REAL :: max3
  REAL :: max2
  REAL :: max1
  REAL :: min2d
!Initialize outputs
  rad_ql = 0.0
  rad_qi = 0.0
  rad_cf = 0.0
  rad_rl = 0.0
  rad_ri = 0.0
!RAD_QR = 0.0
!RAD_QS = 0.0
! Adjust Anvil fractions for warm clouds
  alph = 0.1
  ssd = (-ted)/20.
  ss = (280.-te)/20.
  IF (1.0 .GT. ss) THEN
    ss = ss
  ELSE
    ss = 1.0
    ssd = 0.0
  END IF
  IF (0.0 .LT. ss) THEN
    ss = ss
  ELSE
    ss = 0.0
    ssd = 0.0
  END IF
  ssd = (1.0-alph)*3*ss**2*ssd
  ss = alph + ss**3*(1.0-alph)
  afxd = 0.5*(afd*ss+af*ssd)
  afx = af*ss*0.5
  IF (cf + afx .GT. 1.00) THEN
    rad_cf = 1.00
    rad_cfd = 0.0
  ELSE
    rad_cfd = cfd + afxd
    rad_cf = cf + afx
  END IF
!Total In-cloud liquid
  IF (rad_cf .GT. 10.0e-8) THEN
!0 -> 10e-8 FOR LINEARIZATION PROTECTION
    rad_qld = ((qcllsd+qcland)*rad_cf-(qclls+qclan)*rad_cfd)/rad_cf**2
    rad_ql = (qclls+qclan)/rad_cf
  ELSE
    rad_ql = 0.0
    rad_qld = 0.0
  END IF
  IF (rad_ql .GT. 0.01) THEN
    rad_ql = 0.01
    rad_qld = 0.0
  ELSE
    rad_ql = rad_ql
  END IF
! Total In-cloud ice
  IF (rad_cf .GT. 10.0e-8) THEN
!0 -> 10e-8 FOR LINEARIZATION PROTECTION
    rad_qid = ((qcilsd+qciand)*rad_cf-(qcils+qcian)*rad_cfd)/rad_cf**2
    rad_qi = (qcils+qcian)/rad_cf
  ELSE
    rad_qi = 0.0
    rad_qid = 0.0
  END IF
  IF (rad_qi .GT. 0.01) THEN
    rad_qi = 0.01
    rad_qid = 0.0
  ELSE
    rad_qi = rad_qi
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
      rad_ri_and = ((qcilsd+qciand)*(qcils/rad_ri+qcian/ri_anv)-(qcils+&
&       qcian)*(qcilsd/rad_ri+qciand/ri_anv))/(qcils/rad_ri+qcian/ri_anv&
&       )**2
      rad_ri_an = (qcils+qcian)/(qcils/rad_ri+qcian/ri_anv)
    ELSE
      rad_ri_and = 0.0
    END IF
  ELSE
    rad_ri_and = 0.0
  END IF
  IF (rad_ri .GT. rad_ri_an) THEN
    rad_rid = rad_ri_and
    rad_ri = rad_ri_an
  ELSE
    rad_ri = rad_ri
    rad_rid = 0.0
  END IF
  IF (rad_ri .LT. min_ri) THEN
    rad_ri = min_ri
    rad_rid = 0.0
  ELSE
    rad_ri = rad_ri
  END IF
! Implement ramps for gradual change in effective radius
  IF (pl .LT. 300.) rad_rl = 21.e-6
  IF (pl .GE. 300.) rad_rl = 21.e-6*300./pl
  IF (rad_rl .LT. 10.e-6) THEN
    rad_rl = 10.e-6
  ELSE
    rad_rl = rad_rl
  END IF
! Thicken low high lat clouds
  IF (pl .GE. 775. .AND. te .LE. 275. .AND. tempor .EQ. 1.) THEN
    IF (-(0.1*pl) + 87.5 .GT. 10.) THEN
      x1 = 10.
    ELSE
      x1 = -(0.1*pl) + 87.5
    END IF
    IF (x1 .LT. 5.) THEN
      max1 = 5.
    ELSE
      max1 = x1
    END IF
    rad_rl = max1*1.e-6
  END IF
  IF (pl .GE. 825. .AND. te .LE. 282. .AND. tempor .EQ. 1.) THEN
    IF (0.71*te - 190.25 .LT. 5.) THEN
      max2 = 5.
      max2d = 0.0
    ELSE
      max2d = 0.71*ted
      max2 = 0.71*te - 190.25
    END IF
    rad_rld = 1.e-6*max2d
    rad_rl = max2*1.e-6
  ELSE
    rad_rld = 0.0
  END IF
  IF (pl .GE. 775. .AND. pl .LT. 825. .AND. te .LE. 282. .AND. te .GT. &
&     275. .AND. tempor .EQ. 1.) THEN
    IF (-(0.1*pl) + 0.71*te - 107.75 .GT. 10.) THEN
      min1 = 10.
      min1d = 0.0
    ELSE
      min1d = 0.71*ted
      min1 = -(0.1*pl) + 0.71*te - 107.75
    END IF
    rad_rld = 1.e-6*min1d
    rad_rl = min1*1.e-6
  END IF
  IF (pl .GE. 825. .AND. te .LE. 275. .AND. tempor .EQ. 1.) THEN
    rad_rld = 0.0
    rad_rl = 5.*1.e-6
  END IF
! Thin low tropical clouds
  IF (pl .GE. 950. .AND. te .GE. 285.) THEN
    IF (2.2*te - 617. .GT. 21.) THEN
      min2 = 21.
      min2d = 0.0
    ELSE
      min2d = 2.2*ted
      min2 = 2.2*te - 617.
    END IF
    rad_rld = 1.e-6*min2d
    rad_rl = min2*1.e-6
  END IF
  IF (pl .GE. 925. .AND. te .GE. 290.) THEN
    IF (0.44*pl - 397. .GT. 21.) THEN
      min3 = 21.
    ELSE
      min3 = 0.44*pl - 397.
    END IF
    rad_rl = min3*1.e-6
    rad_rld = 0.0
  END IF
  IF (pl .GE. 925. .AND. pl .LT. 950. .AND. te .GT. 285. .AND. te .LT. &
&     290.) THEN
    IF (0.44*pl + 2.2*te - 1035. .GT. 21.) THEN
      x2 = 21.
      x2d = 0.0
    ELSE
      x2d = 2.2*ted
      x2 = 0.44*pl + 2.2*te - 1035.
    END IF
    IF (x2 .LT. 10.) THEN
      max3 = 10.
      max3d = 0.0
    ELSE
      max3d = x2d
      max3 = x2
    END IF
    rad_rld = 1.e-6*max3d
    rad_rl = max3*1.e-6
  END IF
  IF (pl .GE. 950. .AND. te .GE. 290.) THEN
    rad_rld = 0.0
    rad_rl = 21.*1.e-6
  END IF
  IF (rad_cf .LT. 1.e-5) THEN
    rad_ql = 0.
    rad_qi = 0.
    rad_cf = 0.
!RAD_QR = 0.
!RAD_QS = 0.
    rad_cfd = 0.0
    rad_qid = 0.0
    rad_qld = 0.0
  END IF
END SUBROUTINE RADCOUPLE_D

!  Differentiation of irrad in forward (tangent) mode:
!   variations   of useful results: flxu_dev dfdts_dev flxd_dev
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
SUBROUTINE IRRAD_D(m, np, ple_dev, ta_dev, ta_devd, wa_dev, wa_devd, &
& oa_dev, oa_devd, tb_dev, tb_devd, co2, trace, n2o_dev, ch4_dev, &
& cfc11_dev, cfc12_dev, cfc22_dev, cwc_dev, cwc_devd, fcld_dev, &
& fcld_devd, ict, icb, reff_dev, reff_devd, ns, fs_dev, tg_dev, eg_dev, &
& tv_dev, ev_dev, rv_dev, na, nb, taua_dev, taua_devd, ssaa_dev, &
& ssaa_devd, asya_dev, asya_devd, flxu_dev, flxu_devd, flxd_dev, &
& flxd_devd, dfdts_dev, dfdts_devd, overcastl, aib_ir, awb_ir, aiw_ir, &
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
  REAL, DIMENSION(m), INTENT(IN) :: tb_devd
!Rank 3 (Prognostic variables and tracers)
  REAL, DIMENSION(m, np), INTENT(IN) :: ta_dev, wa_dev, oa_dev, fcld_dev
  REAL, DIMENSION(m, np), INTENT(IN) :: ta_devd, wa_devd, oa_devd, &
& fcld_devd
  REAL, DIMENSION(m, np), INTENT(IN) :: n2o_dev, ch4_dev, cfc11_dev, &
& cfc12_dev, cfc22_dev
  REAL, DIMENSION(m, np + 1), INTENT(IN) :: ple_dev
!Rank 3 (surface types)
  REAL, DIMENSION(m, ns), INTENT(IN) :: fs_dev, tg_dev, tv_dev
  REAL, DIMENSION(m, ns, 10), INTENT(IN) :: eg_dev, ev_dev, rv_dev
!Rank 3 (diagnostic cloud parts)
  REAL, DIMENSION(m, np, 4), INTENT(IN) :: cwc_dev, reff_dev
  REAL, DIMENSION(m, np, 4), INTENT(IN) :: cwc_devd, reff_devd
!Rank 3 (aerosols)
  REAL, DIMENSION(m, np, nb), INTENT(IN) :: taua_dev, ssaa_dev, asya_dev
  REAL, DIMENSION(m, np, nb), INTENT(IN) :: taua_devd, ssaa_devd, &
& asya_devd
!----- OUPUTS -----
!REAL, DIMENSION(m), INTENT(OUT)         :: sfcem_dev
!, flcu_dev, flau_dev
  REAL, DIMENSION(m, np + 1), INTENT(OUT) :: flxu_dev
  REAL, DIMENSION(m, np+1), INTENT(OUT) :: flxu_devd
!, flcd_dev, flad_dev
  REAL, DIMENSION(m, np + 1), INTENT(OUT) :: flxd_dev
  REAL, DIMENSION(m, np+1), INTENT(OUT) :: flxd_devd
  REAL, DIMENSION(m, np + 1), INTENT(OUT) :: dfdts_dev
  REAL, DIMENSION(m, np+1), INTENT(OUT) :: dfdts_devd
!----- LOCALS -----
  REAL, PARAMETER :: cons_grav=9.80665
  INTEGER, PARAMETER :: nx1=26
  INTEGER, PARAMETER :: no1=21
  INTEGER, PARAMETER :: nc1=30
  INTEGER, PARAMETER :: nh1=31
!Temporary arrays
  REAL :: pa(0:np), dt(0:np)
  REAL :: pad(0:np), dtd(0:np)
  REAL :: x1, x2, x3
  REAL :: x1d, x2d, x3d
  REAL :: dh2o(0:np), dcont(0:np), dco2(0:np), do3(0:np)
  REAL :: dh2od(0:np), dcontd(0:np), dco2d(0:np), do3d(0:np)
  REAL :: dn2o(0:np), dch4(0:np)
  REAL :: dn2od(0:np), dch4d(0:np)
  REAL :: df11(0:np), df12(0:np), df22(0:np)
  REAL :: df11d(0:np), df12d(0:np), df22d(0:np)
  REAL :: th2o(6), tcon(3), tco2(6)
  REAL :: th2od(6), tcond(3), tco2d(6)
  REAL :: tn2o(4), tch4(4), tcom(6)
  REAL :: tn2od(4), tch4d(4), tcomd(6)
  REAL :: tf11, tf12, tf22
  REAL :: tf11d, tf12d, tf22d
  REAL :: blayer(0:np+1), blevel(0:np+1)
  REAL :: blayerd(0:np+1), bleveld(0:np+1)
  REAL :: cd(0:np+1), cu(0:np+1)
  REAL :: bd(0:np+1), bu(0:np+1)
  REAL :: bdd(0:np+1), bud(0:np+1)
  REAL :: ad(0:np+1), au(0:np+1)
  REAL :: bs, dbs, rflxs
  REAL :: dp(0:np)
  REAL :: dpd(0:np)
  REAL :: taant
  REAL :: trant, tranal
  REAL :: trantd, tranald
  REAL :: transfc(0:np+1), trantcr(0:np+1), trantca(0:np+1)
  REAL :: transfcd(0:np+1)
  REAL :: flau(0:np+1), flad(0:np+1)
  REAL :: flcu(0:np+1), flcd(0:np+1)
  REAL :: flxu(0:np+1), flxd(0:np+1)
  REAL :: flxud(0:np+1), flxdd(0:np+1)
  REAL :: taerlyr(0:np)
  REAL :: taerlyrd(0:np)
!OVERCAST
  INTEGER :: ncld(3)
  INTEGER :: icx(0:np)
!OVERCAST
  INTEGER :: idx, rc
  INTEGER :: i, j, k, l, ip, iw, ibn, ik, iq, isb, k1, k2, ne
  REAL :: enn(0:np)
  REAL :: ennd(0:np)
  REAL :: cldhi, cldmd, cldlw, tcldlyr(0:np), fclr, fclr_above
  REAL :: cldhid, cldmdd, cldlwd, tcldlyrd(0:np), fclrd, fclr_aboved
  REAL :: x, xx, yy, p1, a1, b1, fk1, a2, b2, fk2
  REAL :: xxd, yyd
  REAL :: w1, ff
  REAL :: ffd
  LOGICAL :: oznbnd, co2bnd, h2otable, conbnd, n2obnd
  LOGICAL :: ch4bnd, combnd, f11bnd, f12bnd, f22bnd, b10bnd
  LOGICAL :: do_aerosol
!Temp arrays and variables for consolidation of tables
  INTEGER, PARAMETER :: max_num_tables=17
  REAL :: exptbl(0:np, max_num_tables)
  REAL :: exptbld(0:np, max_num_tables)
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
  REAL :: dp_pad(np)
  REAL :: taudiaglyr(np, 4)
  REAL :: fcld_col(np)
  REAL :: fcld_cold(np)
  REAL :: reff_col(np, 4)
  REAL :: reff_cold(np, 4)
  REAL :: cwc_col(np, 4)
  REAL :: cwc_cold(np, 4)
  REAL :: h2oexp_tmp(0:np, 5), conexp_tmp(0:np), co2exp_tmp(0:np, 6), &
& n2oexp_tmp(0:np, 2)
  REAL :: h2oexp_tmpd(0:np, 5), co2exp_tmpd(0:np, 6), n2oexp_tmpd(0:np, &
& 2)
  REAL, DIMENSION(m, np, nb) :: taua_dev_tmp, ssaa_dev_tmp, asya_dev_tmp
  REAL, DIMENSION(m, np, nb) :: taua_dev_tmpd, ssaa_dev_tmpd, &
& asya_dev_tmpd
  INTRINSIC MAX
  INTRINSIC EXP
  INTRINSIC MIN
  INTRINSIC ALOG
!BEGIN CALCULATIONS ...
  taua_dev_tmpd = taua_devd
  taua_dev_tmp = taua_dev
  ssaa_dev_tmpd = ssaa_devd
  ssaa_dev_tmp = ssaa_dev
  asya_dev_tmpd = asya_devd
  asya_dev_tmp = asya_dev
  i = 1
  dh2od = 0.0
  dcontd = 0.0
  dtd = 0.0
  reff_cold = 0.0
  fcld_cold = 0.0
  cwc_cold = 0.0
  do3d = 0.0
!do i=1,m
!-----compute layer pressure (pa) and layer temperature minus 250K (dt)
  DO k=1,np
    pad(k) = 0.0
    pa(k) = 0.5*(ple_dev(i, k+1)+ple_dev(i, k))*0.01
    dpd(k) = 0.0
    dp(k) = (ple_dev(i, k+1)-ple_dev(i, k))*0.01
! dp in Pascals for getirtau
    dp_pad(k) = 0.0
    dp_pa(k) = ple_dev(i, k+1) - ple_dev(i, k)
    dtd(k) = ta_devd(i, k)
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
    dh2od(k) = 1.02*dp(k)*wa_devd(i, k)
    dh2o(k) = 1.02*wa_dev(i, k)*dp(k)
    do3d(k) = 476.*dp(k)*oa_devd(i, k)
    do3(k) = 476.*oa_dev(i, k)*dp(k)
    dco2d(k) = 0.0
    dco2(k) = 789.*co2*dp(k)
    dch4d(k) = 0.0
    dch4(k) = 789.*ch4_dev(i, k)*dp(k)
    dn2od(k) = 0.0
    dn2o(k) = 789.*n2o_dev(i, k)*dp(k)
    df11d(k) = 0.0
    df11(k) = 789.*cfc11_dev(i, k)*dp(k)
    df12d(k) = 0.0
    df12(k) = 789.*cfc12_dev(i, k)*dp(k)
    df22d(k) = 0.0
    df22(k) = 789.*cfc22_dev(i, k)*dp(k)
    IF (dh2o(k) .LT. 1.e-10) THEN
      dh2od(k) = 0.0
      dh2o(k) = 1.e-10
    ELSE
      dh2o(k) = dh2o(k)
    END IF
    IF (do3(k) .LT. 1.e-6) THEN
      do3d(k) = 0.0
      do3(k) = 1.e-6
    ELSE
      do3(k) = do3(k)
    END IF
    IF (dco2(k) .LT. 1.e-4) THEN
      dco2d(k) = 0.0
      dco2(k) = 1.e-4
    ELSE
      dco2d(k) = 0.0
      dco2(k) = dco2(k)
    END IF
!-----compute scaled water vapor amount for h2o continuum absorption
!     following eq. (4.21).
    xxd = pa(k)*0.001618*dp(k)*(wa_devd(i, k)*wa_dev(i, k)+wa_dev(i, k)*&
&     wa_devd(i, k))
    xx = pa(k)*0.001618*wa_dev(i, k)*wa_dev(i, k)*dp(k)
    dcontd(k) = xxd*EXP(1800./ta_dev(i, k)-6.081) - xx*1800.*ta_devd(i, &
&     k)*EXP(1800./ta_dev(i, k)-6.081)/ta_dev(i, k)**2
    dcont(k) = xx*EXP(1800./ta_dev(i, k)-6.081)
!-----Fill the reff, cwc, and fcld for the column
    fcld_cold(k) = fcld_devd(i, k)
    fcld_col(k) = fcld_dev(i, k)
    DO l=1,4
      reff_cold(k, l) = reff_devd(i, k, l)
      reff_col(k, l) = reff_dev(i, k, l)
      cwc_cold(k, l) = cwc_devd(i, k, l)
      cwc_col(k, l) = cwc_dev(i, k, l)
    END DO
  END DO
  IF (ple_dev(i, 1)*0.01 .LT. 0.005) THEN
    dpd(0) = 0.0
    dp(0) = 0.005
  ELSE
    dpd(0) = 0.0
    dp(0) = ple_dev(i, 1)*0.01
  END IF
  pad(0) = 0.0
  pa(0) = 0.5*dp(0)
  dtd(0) = ta_devd(i, 1)
  dt(0) = ta_dev(i, 1) - 250.0
  dh2od(0) = 1.02*dp(0)*wa_devd(i, 1)
  dh2o(0) = 1.02*wa_dev(i, 1)*dp(0)
  do3d(0) = 476.*dp(0)*oa_devd(i, 1)
  do3(0) = 476.*oa_dev(i, 1)*dp(0)
  dco2d(0) = 0.0
  dco2(0) = 789.*co2*dp(0)
  dch4d(0) = 0.0
  dch4(0) = 789.*ch4_dev(i, 1)*dp(0)
  dn2od(0) = 0.0
  dn2o(0) = 789.*n2o_dev(i, 1)*dp(0)
  df11d(0) = 0.0
  df11(0) = 789.*cfc11_dev(i, 1)*dp(0)
  df12d(0) = 0.0
  df12(0) = 789.*cfc12_dev(i, 1)*dp(0)
  df22d(0) = 0.0
  df22(0) = 789.*cfc22_dev(i, 1)*dp(0)
  IF (dh2o(0) .LT. 1.e-10) THEN
    dh2od(0) = 0.0
    dh2o(0) = 1.e-10
  ELSE
    dh2o(0) = dh2o(0)
  END IF
  IF (do3(0) .LT. 1.e-6) THEN
    do3d(0) = 0.0
    do3(0) = 1.e-6
  ELSE
    do3(0) = do3(0)
  END IF
  IF (dco2(0) .LT. 1.e-4) THEN
    dco2d(0) = 0.0
    dco2(0) = 1.e-4
  ELSE
    dco2d(0) = 0.0
    dco2(0) = dco2(0)
  END IF
  xxd = pa(0)*0.001618*dp(0)*(wa_devd(i, 1)*wa_dev(i, 1)+wa_dev(i, 1)*&
&   wa_devd(i, 1))
  xx = pa(0)*0.001618*wa_dev(i, 1)*wa_dev(i, 1)*dp(0)
  dcontd(0) = xxd*EXP(1800./ta_dev(i, 1)-6.081) - xx*1800.*ta_devd(i, 1)&
&   *EXP(1800./ta_dev(i, 1)-6.081)/ta_dev(i, 1)**2
  dcont(0) = xx*EXP(1800./ta_dev(i, 1)-6.081)
!-----the surface (np+1) is treated as a layer filled with black clouds.
!     transfc is the transmittance between the surface and a pressure
!     level.
!     trantcr is the clear-sky transmittance between the surface and a
!     pressure level.
!sfcem_dev(i) =0.0
  transfcd(np+1) = 0.0
  transfc(np+1) = 1.0
  trantcr(np+1) = 1.0
  trantca(np+1) = 1.0
!-----initialize fluxes
  DO k=1,np+1
    flxu_devd(i, k) = 0.0
    flxu_dev(i, k) = 0.0
!flcu_dev(i,k)  = 0.0
!flau_dev(i,k)  = 0.0
    flxd_devd(i, k) = 0.0
    flxd_dev(i, k) = 0.0
!flcd_dev(i,k)  = 0.0
!flad_dev(i,k)  = 0.0
    dfdts_devd(i, k) = 0.0
    dfdts_dev(i, k) = 0.0
  END DO
  n2oexp_tmpd = 0.0
  blayerd = 0.0
  tn2od = 0.0
  transfcd = 0.0
  bdd = 0.0
  h2oexp_tmpd = 0.0
  tcomd = 0.0
  tf11d = 0.0
  tf12d = 0.0
  trantd = 0.0
  bud = 0.0
  th2od = 0.0
  tcldlyrd = 0.0
  tch4d = 0.0
  tf22d = 0.0
  ennd = 0.0
  fclrd = 0.0
  tco2d = 0.0
  co2exp_tmpd = 0.0
  taerlyrd = 0.0
  bleveld = 0.0
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
      exptbl = 0.0
!-----Control packing of the new exponential tables by band
      SELECT CASE  (ibn) 
      CASE (2) 
        conexp%start = 1
        conexp%end = 1
      CASE (3) 
        h2oexp%start = 1
        h2oexp%end = 6
        conexp%start = 7
        conexp%end = 9
      CASE (4) 
        h2oexp%start = 1
        h2oexp%end = 6
        conexp%start = 7
        conexp%end = 7
        comexp%start = 8
        comexp%end = 13
        f11exp%start = 14
        f11exp%end = 14
        f12exp%start = 15
        f12exp%end = 15
        f22exp%start = 16
        f22exp%end = 16
      CASE (5) 
        h2oexp%start = 1
        h2oexp%end = 6
        conexp%start = 7
        conexp%end = 7
        comexp%start = 8
        comexp%end = 13
        f11exp%start = 14
        f11exp%end = 14
      CASE (6) 
        h2oexp%start = 1
        h2oexp%end = 6
        conexp%start = 7
        conexp%end = 7
        n2oexp%start = 8
        n2oexp%end = 11
        ch4exp%start = 12
        ch4exp%end = 15
        f12exp%start = 16
        f12exp%end = 16
        f22exp%start = 17
        f22exp%end = 17
      CASE (7) 
        h2oexp%start = 1
        h2oexp%end = 6
        conexp%start = 7
        conexp%end = 7
        n2oexp%start = 8
        n2oexp%end = 11
        ch4exp%start = 12
        ch4exp%end = 15
      CASE (9) 
        h2oexp%start = 1
        h2oexp%end = 6
      CASE (10) 
        h2oexp%start = 1
        h2oexp%end = 5
        conexp%start = 6
        conexp%end = 6
        co2exp%start = 7
        co2exp%end = 12
        n2oexp%start = 13
        n2oexp%end = 14
      END SELECT
!-----blayer is the spectrally integrated planck flux of the mean layer
!     temperature derived from eq. (3.11)
!     The fitting for the planck flux is valid for the range 160-345 K.
      DO k=1,np
        CALL PLANCK_D(ibn, cb, ta_dev(i, k), ta_devd(i, k), blayer(k), &
&               blayerd(k))
      END DO
!-----Index "0" is the layer above the top of the atmosphere.
      blayerd(0) = blayerd(1)
      blayer(0) = blayer(1)
      bleveld(0) = blayerd(1)
      blevel(0) = blayer(1)
!-----Surface emission and reflectivity. See Section 9.
!     bs and dbs include the effect of surface emissivity.
      CALL SFCFLUX(ibn, m, i, cb, dcb, ns, fs_dev, tg_dev, eg_dev, &
&            tv_dev, ev_dev, rv_dev, bs, dbs, rflxs)
      blayerd(np+1) = 0.0
      blayer(np+1) = bs
!------interpolate Planck function at model levels (linear in p)
      DO k=2,np
        bleveld(k) = (dp(k)*blayerd(k-1)+dp(k-1)*blayerd(k))/(dp(k-1)+dp&
&         (k))
        blevel(k) = (blayer(k-1)*dp(k)+blayer(k)*dp(k-1))/(dp(k-1)+dp(k)&
&         )
      END DO
!-----Extrapolate blevel(1) from blayer(2) and blayer(1)
      bleveld(1) = blayerd(1) + dp(1)*(blayerd(1)-blayerd(2))/(dp(1)+dp(&
&       2))
      blevel(1) = blayer(1) + (blayer(1)-blayer(2))*dp(1)/(dp(1)+dp(2))
      bleveld(0) = bleveld(1)
      blevel(0) = blevel(1)
!-----If the surface air temperature tb is known, compute blevel(np+1)
      CALL PLANCK_D(ibn, cb, tb_dev(i), tb_devd(i), blevel(np+1), &
&             bleveld(np+1))
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
      CALL GETIRTAU1_D(ibn, np, dp_pa, fcld_col, fcld_cold, reff_col, &
&                reff_cold, cwc_col, cwc_cold, tcldlyr, tcldlyrd, enn, &
&                ennd, aib_ir, awb_ir, aiw_ir, aww_ir, aig_ir, awg_ir, &
&                cons_grav)
!     The getirtau call returns results for a single band and column. Thus
!     we need to transfer the taudiag for that band/column to the overall array.
!MAT-- icx and ncld only used when overcast=.false.
!Overcast
      IF (overcastl .EQV. .false.) THEN
        DO k=0,np
          icx(k) = k
        END DO
        CALL MKICX(np, ict, icb, enn, icx, ncld)
      END IF
!-----Compute optical thickness, single-scattering albedo and asymmetry
!     factor for a mixture of "na" aerosol types. Eqs. (7.1)-(7.3)
      IF (do_aerosol) THEN
        taerlyrd(0) = 0.0
        taerlyr(0) = 1.0
        DO k=1,np
!-----taerlyr is the aerosol diffuse transmittance
          taerlyrd(k) = 0.0
          taerlyr(k) = 1.0
          IF (taua_dev_tmp(i, k, ibn) .GT. 0.001) THEN
            IF (ssaa_dev_tmp(i, k, ibn) .GT. 0.001) THEN
              asya_dev_tmpd(i, k, ibn) = (asya_dev_tmpd(i, k, ibn)*&
&               ssaa_dev_tmp(i, k, ibn)-asya_dev_tmp(i, k, ibn)*&
&               ssaa_dev_tmpd(i, k, ibn))/ssaa_dev_tmp(i, k, ibn)**2
              asya_dev_tmp(i, k, ibn) = asya_dev_tmp(i, k, ibn)/&
&               ssaa_dev_tmp(i, k, ibn)
              ssaa_dev_tmpd(i, k, ibn) = (ssaa_dev_tmpd(i, k, ibn)*&
&               taua_dev_tmp(i, k, ibn)-ssaa_dev_tmp(i, k, ibn)*&
&               taua_dev_tmpd(i, k, ibn))/taua_dev_tmp(i, k, ibn)**2
              ssaa_dev_tmp(i, k, ibn) = ssaa_dev_tmp(i, k, ibn)/&
&               taua_dev_tmp(i, k, ibn)
!-----Parameterization of aerosol scattering following Eqs. (6.11)
!     and (6.12). 
              ffd = (0.1185*asya_dev_tmpd(i, k, ibn)*asya_dev_tmp(i, k, &
&               ibn)+(0.0076+0.1185*asya_dev_tmp(i, k, ibn))*&
&               asya_dev_tmpd(i, k, ibn))*asya_dev_tmp(i, k, ibn) + (&
&               .3739+(0.0076+0.1185*asya_dev_tmp(i, k, ibn))*&
&               asya_dev_tmp(i, k, ibn))*asya_dev_tmpd(i, k, ibn)
              ff = .5 + (.3739+(0.0076+0.1185*asya_dev_tmp(i, k, ibn))*&
&               asya_dev_tmp(i, k, ibn))*asya_dev_tmp(i, k, ibn)
              taua_dev_tmpd(i, k, ibn) = taua_dev_tmpd(i, k, ibn)*(1.-&
&               ssaa_dev_tmp(i, k, ibn)*ff) + taua_dev_tmp(i, k, ibn)*(-&
&               (ssaa_dev_tmpd(i, k, ibn)*ff)-ssaa_dev_tmp(i, k, ibn)*&
&               ffd)
              taua_dev_tmp(i, k, ibn) = taua_dev_tmp(i, k, ibn)*(1.-&
&               ssaa_dev_tmp(i, k, ibn)*ff)
            END IF
            taerlyrd(k) = -(1.66*taua_dev_tmpd(i, k, ibn)*EXP(-(1.66*&
&             taua_dev_tmp(i, k, ibn))))
            taerlyr(k) = EXP(-(1.66*taua_dev_tmp(i, k, ibn)))
          END IF
        END DO
      END IF
!-----Compute the exponential terms (Eq. 8.21) at each layer due to
!     water vapor line absorption when k-distribution is used
      IF (.NOT.h2otable .AND. (.NOT.b10bnd)) THEN
        CALL H2OEXPS_D(ibn, np, dh2o, dh2od, pa, dt, dtd, xkw, aw, bw, &
&                pm, mw, exptbl(:, h2oexp%start:h2oexp%end), exptbld(:, &
&                h2oexp%start:h2oexp%end))
      ELSE
        exptbld = 0.0
      END IF
!-----compute the exponential terms (Eq. 4.24) at each layer due to
!     water vapor continuum absorption.
!     ne is the number of terms used in each band to compute water 
!     vapor continuum transmittance (Table 9).
      ne = 0
      IF (conbnd) THEN
        ne = 1
        IF (ibn .EQ. 3) ne = 3
        CALL CONEXPS_D(ibn, np, dcont, dcontd, xke, exptbl(:, conexp%&
&                start:conexp%end), exptbld(:, conexp%start:conexp%end))
      END IF
!----- for trace gases 
      IF (trace) THEN
!-----compute the exponential terms at each layer due to n2o absorption
        IF (n2obnd) CALL N2OEXPS_D(ibn, np, dn2o, pa, dt, dtd, exptbl(:&
&                            , n2oexp%start:n2oexp%end), exptbld(:, &
&                            n2oexp%start:n2oexp%end))
!-----compute the exponential terms at each layer due to ch4 absorption
        IF (ch4bnd) CALL CH4EXPS_D(ibn, np, dch4, pa, dt, dtd, exptbl(:&
&                            , ch4exp%start:ch4exp%end), exptbld(:, &
&                            ch4exp%start:ch4exp%end))
!-----Compute the exponential terms due to co2 minor absorption
        IF (combnd) CALL COMEXPS_D(ibn, np, dco2, dt, dtd, exptbl(:, &
&                            comexp%start:comexp%end), exptbld(:, comexp&
&                            %start:comexp%end))
!-----Compute the exponential terms due to cfc11 absorption.
!     The values of the parameters are given in Table 7.
        IF (f11bnd) THEN
          a1 = 1.26610e-3
          b1 = 3.55940e-6
          fk1 = 1.89736e+1
          a2 = 8.19370e-4
          b2 = 4.67810e-6
          fk2 = 1.01487e+1
          CALL CFCEXPS_D(ibn, np, a1, b1, fk1, a2, b2, fk2, df11, dt, &
&                  dtd, exptbl(:, f11exp%start:f11exp%end), exptbld(:, &
&                  f11exp%start:f11exp%end))
        END IF
!-----Compute the exponential terms due to cfc12 absorption.
        IF (f12bnd) THEN
          a1 = 8.77370e-4
          b1 = -5.88440e-6
          fk1 = 1.58104e+1
          a2 = 8.62000e-4
          b2 = -4.22500e-6
          fk2 = 3.70107e+1
          CALL CFCEXPS_D(ibn, np, a1, b1, fk1, a2, b2, fk2, df12, dt, &
&                  dtd, exptbl(:, f12exp%start:f12exp%end), exptbld(:, &
&                  f12exp%start:f12exp%end))
        END IF
!-----Compute the exponential terms due to cfc22 absorption.
        IF (f22bnd) THEN
          a1 = 9.65130e-4
          b1 = 1.31280e-5
          fk1 = 6.18536e+0
          a2 = -3.00010e-5
          b2 = 5.25010e-7
          fk2 = 3.27912e+1
          CALL CFCEXPS_D(ibn, np, a1, b1, fk1, a2, b2, fk2, df22, dt, &
&                  dtd, exptbl(:, f22exp%start:f22exp%end), exptbld(:, &
&                  f22exp%start:f22exp%end))
        END IF
!-----Compute the exponential terms at each layer in band 10 due to
!     h2o line and continuum, co2, and n2o absorption
        IF (b10bnd) THEN
          CALL B10EXPS_D(np, dh2o, dh2od, dcont, dcontd, dco2, dn2o, pa&
&                  , dt, dtd, h2oexp_tmp, h2oexp_tmpd, exptbl(:, conexp%&
&                  start:conexp%end), exptbld(:, conexp%start:conexp%end&
&                  ), co2exp_tmp, co2exp_tmpd, n2oexp_tmp, n2oexp_tmpd)
          exptbld(:, h2oexp%start:h2oexp%end) = h2oexp_tmpd
          exptbl(:, h2oexp%start:h2oexp%end) = h2oexp_tmp
!              exptbl(:,conexp%start:conexp%end) = conexp_tmp
          exptbld(:, co2exp%start:co2exp%end) = co2exp_tmpd
          exptbl(:, co2exp%start:co2exp%end) = co2exp_tmp
          exptbld(:, n2oexp%start:n2oexp%end) = n2oexp_tmpd
          exptbl(:, n2oexp%start:n2oexp%end) = n2oexp_tmp
        END IF
      END IF
!-----blayer(np+1) includes the effect of surface emissivity.
! ALT: this was undefined, check with Max if 0.0 is good value
      bud(0) = 0.0
      bu(0) = 0.0
      bdd(0) = blayerd(1)
      bd(0) = blayer(1)
      bud(np+1) = blayerd(np+1)
      bu(np+1) = blayer(np+1)
! ALT: this was undefined, check with Max if 0.0 is good value
      au(0) = 0.0
      ad(0) = blayer(1)
      au(np+1) = blayer(np+1)
! ALT: this was undefined, check with Max if 0.0 is good value
      cu(0) = 0.0
      cd(0) = blayer(1)
      cu(np+1) = blayer(np+1)
!-----do-loop 1500 is for computing upward (bu) and downward (bd)
!     emission of a layer following Eqs. (8.17), (8.18), (8.19).
!     Here, trant is the transmittance of the layer k2-1.
      DO k2=1,np+1
!-----for h2o line transmission
        IF (.NOT.h2otable) THEN
          th2o = 1.0
          th2od = 0.0
        END IF
!-----for h2o continuum transmission
        tcon = 1.0
        x1 = 0.0
        x2 = 0.0
        x3 = 0.0
        taant = 1.0
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
            trantd = 0.0
            x3d = 0.0
            x2d = 0.0
            x1d = 0.0
            CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2-1), pa(k2-1), &
&                   dt(k2-1), dtd(k2-1), x1, x1d, x2, x2d, x3, x3d, w11&
&                   , p11, dwe, dpe, h11, h12, h13, trant, trantd)
          ELSE
            trantd = 0.0
            x1d = 0.0
            x2d = 0.0
            x3d = 0.0
          END IF
          IF (ibn .EQ. 2) CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2-1&
&                                 ), pa(k2-1), dt(k2-1), dtd(k2-1), x1, &
&                                 x1d, x2, x2d, x3, x3d, w11, p11, dwe, &
&                                 dpe, h21, h22, h23, trant, trantd)
          IF (ibn .EQ. 8) CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2-1&
&                                 ), pa(k2-1), dt(k2-1), dtd(k2-1), x1, &
&                                 x1d, x2, x2d, x3, x3d, w11, p11, dwe, &
&                                 dpe, h81, h82, h83, trant, trantd)
!bdc
!-----for water vapor continuum absorption
          IF (conbnd) THEN
! Only the first exp
            tcond = 0.0
            tcond(1) = tcon(1)*exptbld(k2-1, conexp%start)
            tcon(1) = tcon(1)*exptbl(k2-1, conexp%start)
            trantd = trantd*tcon(1) + trant*tcond(1)
            trant = trant*tcon(1)
          END IF
        ELSE IF (.NOT.b10bnd) THEN
!-----compute water vapor transmittance using k-distribution
          tcond = 0.0
          CALL H2OKDIS_D(ibn, np, k2 - 1, fkw, gkw, ne, exptbl(:, h2oexp&
&                  %start:h2oexp%end), exptbld(:, h2oexp%start:h2oexp%&
&                  end), exptbl(:, conexp%start:conexp%end), exptbld(:, &
&                  conexp%start:conexp%end), th2o, th2od, tcon, tcond, &
&                  trant, trantd)
          x1d = 0.0
          x2d = 0.0
          x3d = 0.0
        ELSE
          trantd = 0.0
          x1d = 0.0
          x2d = 0.0
          x3d = 0.0
        END IF
        IF (co2bnd) THEN
!-----Compute co2 transmittance using table look-up method.
!     The following values are taken from Table 8.
!              w1=-4.0
!              p1=-2.0
!              dwe=0.3
!              dpe=0.2
          dco2d = 0.0
          CALL TABLUP_D(nx1, nc1, dco2(k2-1), dco2d(k2-1), pa(k2-1), dt(&
&                 k2-1), dtd(k2-1), x1, x1d, x2, x2d, x3, x3d, w12, p12&
&                 , dwe, dpe, c1, c2, c3, trant, trantd)
        END IF
!-----Always use table look-up to compute o3 transmittance.
!     The following values are taken from Table 8.
        IF (oznbnd) CALL TABLUP_D(nx1, no1, do3(k2-1), do3d(k2-1), pa(k2&
&                           -1), dt(k2-1), dtd(k2-1), x1, x1d, x2, x2d, &
&                           x3, x3d, w13, p13, dwe, dpe, oo1, oo2, oo3, &
&                           trant, trantd)
!              w1=-6.0
!              p1=-2.0
!              dwe=0.3
!              dpe=0.2
!-----include aerosol effect
        taant = trant
        IF (do_aerosol) THEN
          trantd = trantd*taerlyr(k2-1) + trant*taerlyrd(k2-1)
          trant = trant*taerlyr(k2-1)
        END IF
!-----Compute upward (bu) and downward (bd) emission of the layer k2-1
!     following Eqs.(8.17) and (8.18).
!     The effect of clouds on the transmission of a layer is taken
!     into account, following Eq. (8.19).
!     trant is the total transmittance of the layer k2-1.
        xxd = (1.-enn(k2-1))*trantd - ennd(k2-1)*trant
        xx = (1.-enn(k2-1))*trant
        IF (0.9999 .GT. xx) THEN
          yyd = xxd
          yy = xx
        ELSE
          yy = 0.9999
          yyd = 0.0
        END IF
        IF (0.00001 .LT. yy) THEN
          yy = yy
        ELSE
          yy = 0.00001
          yyd = 0.0
        END IF
        xxd = ((bleveld(k2-1)-bleveld(k2))*ALOG(yy)-(blevel(k2-1)-blevel&
&         (k2))*yyd/yy)/ALOG(yy)**2
        xx = (blevel(k2-1)-blevel(k2))/ALOG(yy)
        bdd(k2-1) = ((bleveld(k2)-bleveld(k2-1)*yy-blevel(k2-1)*yyd)*(&
&         1.0-yy)+(blevel(k2)-blevel(k2-1)*yy)*yyd)/(1.0-yy)**2 - xxd
        bd(k2-1) = (blevel(k2)-blevel(k2-1)*yy)/(1.0-yy) - xx
        bud(k2-1) = bleveld(k2-1) + bleveld(k2) - bdd(k2-1)
        bu(k2-1) = blevel(k2-1) + blevel(k2) - bd(k2-1)
        xx = trant
        IF (0.9999 .GT. xx) THEN
          yy = xx
        ELSE
          yy = 0.9999
        END IF
        IF (0.00001 .LT. yy) THEN
          yy = yy
        ELSE
          yy = 0.00001
        END IF
        xx = (blevel(k2-1)-blevel(k2))/ALOG(yy)
        cd(k2-1) = (blevel(k2)-blevel(k2-1)*yy)/(1.0-yy) - xx
        cu(k2-1) = blevel(k2-1) + blevel(k2) - cd(k2-1)
        IF (do_aerosol) THEN
          xx = taant
          IF (0.9999 .GT. xx) THEN
            yy = xx
          ELSE
            yy = 0.9999
          END IF
          IF (0.00001 .LT. yy) THEN
            yy = yy
          ELSE
            yy = 0.00001
          END IF
          xx = (blevel(k2-1)-blevel(k2))/ALOG(yy)
          ad(k2-1) = (blevel(k2)-blevel(k2-1)*yy)/(1.0-yy) - xx
          au(k2-1) = blevel(k2-1) + blevel(k2) - ad(k2-1)
        ELSE
          ad(k2-1) = cd(k2-1)
          au(k2-1) = cu(k2-1)
        END IF
      END DO
!-----initialize fluxes
      flxu = 0.0
      flxd = 0.0
      flcu = 0.0
      flcd = 0.0
      flau = 0.0
      flad = 0.0
      flxud = 0.0
      flxdd = 0.0
!-----Compute upward and downward fluxes for each spectral band, ibn.
      DO k1=0,np
!-----initialization
!
!     cldlw, cldmd, and cldhi are the equivalent black-cloud fractions
!     of low, middle, and high troposphere.
!     tranal is the aerosol transmission function
        cldlw = 0.0
        cldmd = 0.0
        cldhi = 0.0
        tranal = 1.0
!-----for h2o line transmission
        IF (.NOT.h2otable) THEN
          th2o = 1.0
          th2od = 0.0
        END IF
!-----for h2o continuum transmission
        tcon = 1.0
!----- for trace gases
        IF (trace) THEN
!-----for n2o transmission using k-distribution method.
          IF (n2obnd) THEN
            tn2o = 1.0
            tn2od = 0.0
          END IF
!-----for ch4 transmission using k-distribution method.
          IF (ch4bnd) THEN
            tch4 = 1.0
            tch4d = 0.0
          END IF
!-----for co2-minor transmission using k-distribution method.
          IF (combnd) THEN
            tcom = 1.0
            tcomd = 0.0
          END IF
!-----for cfc-11 transmission using k-distribution method.
          IF (f11bnd) THEN
            tf11 = 1.0
            tf11d = 0.0
          END IF
!-----for cfc-12 transmission using k-distribution method.
          IF (f12bnd) THEN
            tf12 = 1.0
            tf12d = 0.0
          END IF
!-----for cfc-22 transmission when using k-distribution method.
          IF (f22bnd) THEN
            tf22 = 1.0
            tf22d = 0.0
          END IF
!-----for the transmission in band 10 using k-distribution method.
          IF (b10bnd) THEN
            th2o = 1.0
            tco2 = 1.0
            tcond(1) = 0.0
            tcon(1) = 1.0
            tn2o = 1.0
            tn2od = 0.0
            th2od = 0.0
            tco2d = 0.0
          END IF
        END IF
!----- end trace gases
        x1 = 0.0
        x2 = 0.0
        x3 = 0.0
!-----do-loop 3000 are for computing (a) transmittance, trant,
!     and (b) clear line-of-sight, fclr(k2), between levels k1 and k2.
        fclr_above = 1.0
        cldhid = 0.0
        tcond = 0.0
        tranald = 0.0
        x1d = 0.0
        x2d = 0.0
        x3d = 0.0
        cldlwd = 0.0
        fclr_aboved = 0.0
        cldmdd = 0.0
!MAT--Beginning of original 3000 loop
        DO k2=k1+1,np+1
          taant = 1.0
          trant = 1.0
          fclr = 1.0
          IF (h2otable) THEN
!-----Compute water vapor transmittance using table look-up.
!     The following values are taken from Table 8.
!                 w1=-8.0
!                 p1=-2.0
!                 dwe=0.3
!                 dpe=0.2
            IF (ibn .EQ. 1) THEN
              trantd = 0.0
              CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2-1), pa(k2-1)&
&                     , dt(k2-1), dtd(k2-1), x1, x1d, x2, x2d, x3, x3d, &
&                     w11, p11, dwe, dpe, h11, h12, h13, trant, trantd)
            ELSE
              trantd = 0.0
            END IF
            IF (ibn .EQ. 2) CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2&
&                                   -1), pa(k2-1), dt(k2-1), dtd(k2-1), &
&                                   x1, x1d, x2, x2d, x3, x3d, w11, p11&
&                                   , dwe, dpe, h21, h22, h23, trant, &
&                                   trantd)
            IF (ibn .EQ. 8) CALL TABLUP_D(nx1, nh1, dh2o(k2-1), dh2od(k2&
&                                   -1), pa(k2-1), dt(k2-1), dtd(k2-1), &
&                                   x1, x1d, x2, x2d, x3, x3d, w11, p11&
&                                   , dwe, dpe, h81, h82, h83, trant, &
&                                   trantd)
            IF (conbnd) THEN
! Only the first exp
              tcond(1) = tcond(1)*exptbl(k2-1, conexp%start) + tcon(1)*&
&               exptbld(k2-1, conexp%start)
              tcon(1) = tcon(1)*exptbl(k2-1, conexp%start)
              trantd = trantd*tcon(1) + trant*tcond(1)
              trant = trant*tcon(1)
            END IF
          ELSE IF (.NOT.b10bnd) THEN
!-----compute water vapor transmittance using k-distribution
            CALL H2OKDIS_D(ibn, np, k2 - 1, fkw, gkw, ne, exptbl(:, &
&                    h2oexp%start:h2oexp%end), exptbld(:, h2oexp%start:&
&                    h2oexp%end), exptbl(:, conexp%start:conexp%end), &
&                    exptbld(:, conexp%start:conexp%end), th2o, th2od, &
&                    tcon, tcond, trant, trantd)
          ELSE
            trantd = 0.0
          END IF
          IF (co2bnd) THEN
!-----Compute co2 transmittance using table look-up method.
!     The following values are taken from Table 8.
!                 w1=-4.0
!                 p1=-2.0
!                 dwe=0.3
!                 dpe=0.2
            dco2d = 0.0
            CALL TABLUP_D(nx1, nc1, dco2(k2-1), dco2d(k2-1), pa(k2-1), &
&                   dt(k2-1), dtd(k2-1), x1, x1d, x2, x2d, x3, x3d, w12&
&                   , p12, dwe, dpe, c1, c2, c3, trant, trantd)
          END IF
!-----Always use table look-up to compute o3 transmittance.
!     The following values are taken from Table 8.
          IF (oznbnd) CALL TABLUP_D(nx1, no1, do3(k2-1), do3d(k2-1), pa(&
&                             k2-1), dt(k2-1), dtd(k2-1), x1, x1d, x2, &
&                             x2d, x3, x3d, w13, p13, dwe, dpe, oo1, oo2&
&                             , oo3, trant, trantd)
!                 w1=-6.0
!                 p1=-2.0
!                 dwe=0.3
!                 dpe=0.2
!------ for trace gases 
          IF (trace) THEN
!-----compute n2o transmittance using k-distribution method
            IF (n2obnd) CALL N2OKDIS_D(ibn, np, k2 - 1, exptbl(:, n2oexp&
&                                %start:n2oexp%end), exptbld(:, n2oexp%&
&                                start:n2oexp%end), tn2o, tn2od, trant, &
&                                trantd)
!-----compute ch4 transmittance using k-distribution method
            IF (ch4bnd) CALL CH4KDIS_D(ibn, np, k2 - 1, exptbl(:, ch4exp&
&                                %start:ch4exp%end), exptbld(:, ch4exp%&
&                                start:ch4exp%end), tch4, tch4d, trant, &
&                                trantd)
!-----compute co2-minor transmittance using k-distribution method
            IF (combnd) CALL COMKDIS_D(ibn, np, k2 - 1, exptbl(:, comexp&
&                                %start:comexp%end), exptbld(:, comexp%&
&                                start:comexp%end), tcom, tcomd, trant, &
&                                trantd)
!-----compute cfc11 transmittance using k-distribution method
            IF (f11bnd) CALL CFCKDIS_D(np, k2 - 1, exptbl(:, f11exp%&
&                                start:f11exp%end), exptbld(:, f11exp%&
&                                start:f11exp%end), tf11, tf11d, trant, &
&                                trantd)
!-----compute cfc12 transmittance using k-distribution method
            IF (f12bnd) CALL CFCKDIS_D(np, k2 - 1, exptbl(:, f12exp%&
&                                start:f12exp%end), exptbld(:, f12exp%&
&                                start:f12exp%end), tf12, tf12d, trant, &
&                                trantd)
!-----compute cfc22 transmittance using k-distribution method
            IF (f22bnd) CALL CFCKDIS_D(np, k2 - 1, exptbl(:, f22exp%&
&                                start:f22exp%end), exptbld(:, f22exp%&
&                                start:f22exp%end), tf22, tf22d, trant, &
&                                trantd)
!-----Compute transmittance in band 10 using k-distribution method.
!     For band 10, trant is the change in transmittance due to n2o 
!     absorption.
            IF (b10bnd) CALL B10KDIS_D(np, k2 - 1, exptbl(:, h2oexp%&
&                                start:h2oexp%end), exptbld(:, h2oexp%&
&                                start:h2oexp%end), exptbl(:, conexp%&
&                                start:conexp%end), exptbld(:, conexp%&
&                                start:conexp%end), exptbl(:, co2exp%&
&                                start:co2exp%end), exptbld(:, co2exp%&
&                                start:co2exp%end), exptbl(:, n2oexp%&
&                                start:n2oexp%end), exptbld(:, n2oexp%&
&                                start:n2oexp%end), th2o, th2od, tcon, &
&                                tcond, tco2, tco2d, tn2o, tn2od, trant&
&                                , trantd)
          END IF
!-----   end trace gases
!-----include aerosol effect
          taant = trant
          IF (do_aerosol) THEN
            tranald = tranald*taerlyr(k2-1) + tranal*taerlyrd(k2-1)
            tranal = tranal*taerlyr(k2-1)
            trantd = trantd*tranal + trant*tranald
            trant = trant*tranal
          END IF
!----- cloud overlapping 
!OVERCAST
          IF (overcastl .EQV. .false.) THEN
            IF (enn(k2-1) .GE. 0.001) CALL CLDOVLP_D(np, k1, k2, ict, &
&                                              icb, icx, ncld, enn, ennd&
&                                              , tcldlyr, tcldlyrd, &
&                                              cldhi, cldhid, cldmd, &
&                                              cldmdd, cldlw, cldlwd)
            fclrd = (-(cldhid*(1.0-cldmd))-(1.0-cldhi)*cldmdd)*(1.0-&
&             cldlw) - (1.0-cldhi)*(1.0-cldmd)*cldlwd
            fclr = (1.0-cldhi)*(1.0-cldmd)*(1.0-cldlw)
          ELSE
            fclrd = fclr_aboved*tcldlyr(k2-1) + fclr_above*tcldlyrd(k2-1&
&             )
            fclr = fclr_above*tcldlyr(k2-1)
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
            flau(k1) = flau(k1) - au(k1)
            flad(k2) = flad(k2) + ad(k1)
            flcu(k1) = flcu(k1) - cu(k1)
            flcd(k2) = flcd(k2) + cd(k1)
            flxud(k1) = flxud(k1) - bud(k1)
            flxu(k1) = flxu(k1) - bu(k1)
            flxdd(k2) = flxdd(k2) + bdd(k1)
            flxd(k2) = flxd(k2) + bd(k1)
          END IF
!-----The summation terms on the rhs of Eqs. (8.15) and (8.16).
!     Also see Eqs. (5.4) and (5.5) for Band 10.
          xxd = trantd*(bu(k2-1)-bu(k2)) + trant*(bud(k2-1)-bud(k2))
          xx = trant*(bu(k2-1)-bu(k2))
          flxud(k1) = flxud(k1) + xxd*fclr + xx*fclrd
          flxu(k1) = flxu(k1) + xx*fclr
          xx = trant*(cu(k2-1)-cu(k2))
          flcu(k1) = flcu(k1) + xx
          IF (do_aerosol) xx = taant*(au(k2-1)-au(k2))
          flau(k1) = flau(k1) + xx
          IF (k1 .EQ. 0) THEN
!mjs  bd(-1) is not defined
            xxd = -(trantd*bd(k1)+trant*bdd(k1))
            xx = -(trant*bd(k1))
          ELSE
            xxd = trantd*(bd(k1-1)-bd(k1)) + trant*(bdd(k1-1)-bdd(k1))
            xx = trant*(bd(k1-1)-bd(k1))
          END IF
          flxdd(k2) = flxdd(k2) + xxd*fclr + xx*fclrd
          flxd(k2) = flxd(k2) + xx*fclr
          IF (k1 .EQ. 0) THEN
!mjs  bd(-1) is not defined
            xx = -(trant*cd(k1))
          ELSE
            xx = trant*(cd(k1-1)-cd(k1))
          END IF
          flcd(k2) = flcd(k2) + xx
          IF (do_aerosol) THEN
            IF (k1 .EQ. 0) THEN
!mjs  bd(-1) is not defined
              xx = -(taant*ad(k1))
            ELSE
              xx = taant*(ad(k1-1)-ad(k1))
            END IF
          END IF
          flad(k2) = flad(k2) + xx
!MAT--End of original 4000 loop
          fclr_aboved = fclrd
          fclr_above = fclr
        END DO
!-----Here, fclr and trant are, respectively, the clear line-of-sight 
!     and the transmittance between k1 and the surface.
        trantca(k1) = taant
        trantcr(k1) = trant
        transfcd(k1) = trantd*fclr + trant*fclrd
        transfc(k1) = trant*fclr
!-----compute the partial derivative of fluxes with respect to
!     surface temperature (Eq. 3.12). 
!     Note: upward flux is negative, and so is dfdts.
        IF (k1 .GT. 0) THEN
          dfdts_devd(i, k1) = dfdts_devd(i, k1) - dbs*transfcd(k1)
          dfdts_dev(i, k1) = dfdts_dev(i, k1) - dbs*transfc(k1)
        END IF
      END DO
!-----For surface emission.
      IF (.NOT.b10bnd) THEN
!-----Note: blayer(np+1) and dbs include the surface emissivity 
!     effect. Both dfdts and sfcem are negative quantities.
        flau(np+1) = -blayer(np+1)
        flcu(np+1) = -blayer(np+1)
        flxud(np+1) = -blayerd(np+1)
        flxu(np+1) = -blayer(np+1)
!sfcem_dev(i)     = sfcem_dev(i)-blayer(np+1)
        dfdts_dev(i, np+1) = dfdts_dev(i, np+1) - dbs
!-----Add the flux reflected by the surface. (Second term on the
!     rhs of Eq. 8.16). rflxs is the surface reflectivity.
        DO k=1,np+1
          flau(k) = flau(k) - flad(np+1)*trantca(k)*rflxs
          flcu(k) = flcu(k) - flcd(np+1)*trantcr(k)*rflxs
          flxud(k) = flxud(k) - rflxs*(flxdd(np+1)*transfc(k)+flxd(np+1)&
&           *transfcd(k))
          flxu(k) = flxu(k) - flxd(np+1)*transfc(k)*rflxs
        END DO
      END IF
!-----Summation of fluxes over spectral bands
      DO k=1,np+1
!flau_dev(i,k) = flau_dev(i,k) + flau(k)
!flcu_dev(i,k) = flcu_dev(i,k) + flcu(k)
        flxu_devd(i, k) = flxu_devd(i, k) + flxud(k)
        flxu_dev(i, k) = flxu_dev(i, k) + flxu(k)
!flad_dev(i,k) = flad_dev(i,k) + flad(k)
!flcd_dev(i,k) = flcd_dev(i,k) + flcd(k)
        flxd_devd(i, k) = flxd_devd(i, k) + flxdd(k)
        flxd_dev(i, k) = flxd_dev(i, k) + flxd(k)
      END DO
    END IF
  END DO
  GOTO 110
 100 RETURN
 110 CONTINUE
END SUBROUTINE IRRAD_D

!  Differentiation of planck in forward (tangent) mode:
!   variations   of useful results: xlayer
!   with respect to varying inputs: t
!----------------------------------!
!      SUBROUTINES GO HERE         !
!----------------------------------!
SUBROUTINE PLANCK_D(ibn, cb, t, td, xlayer, xlayerd)
  IMPLICIT NONE
! spectral band index
  INTEGER :: ibn
! Planck table coefficients
  REAL :: cb(6, 10)
! temperature (K)
  REAL :: t
  REAL :: td
! planck flux (w/m2)
  REAL :: xlayer
  REAL :: xlayerd
  xlayerd = td*(t*(t*(t*(t*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+cb(3, ibn)&
&   )+cb(2, ibn)) + t*(td*(t*(t*(t*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+cb&
&   (3, ibn))+t*(td*(t*(t*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+t*(td*(t*cb&
&   (6, ibn)+cb(5, ibn))+t*cb(6, ibn)*td)))
  xlayer = t*(t*(t*(t*(t*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+cb(3, ibn))+&
&   cb(2, ibn)) + cb(1, ibn)
END SUBROUTINE PLANCK_D

!  Differentiation of h2oexps in forward (tangent) mode:
!   variations   of useful results: h2oexp
!   with respect to varying inputs: dh2o dt
SUBROUTINE H2OEXPS_D(ib, np, dh2o, dh2od, pa, dt, dtd, xkw, aw, bw, pm, &
& mw, h2oexp, h2oexpd)
  IMPLICIT NONE
  INTEGER :: ib, np, ik, k
!---- input parameters ------
  REAL :: dh2o(0:np), pa(0:np), dt(0:np)
  REAL :: dh2od(0:np), dtd(0:np)
!---- output parameters -----
  REAL :: h2oexp(0:np, 6)
  REAL :: h2oexpd(0:np, 6)
!---- static data -----
  INTEGER :: mw(9)
  REAL :: xkw(9), aw(9), bw(9), pm(9)
!---- temporary arrays -----
  REAL :: xh
  REAL :: xhd
  INTRINSIC EXP
  REAL :: pwx1
  REAL :: pwr1
  h2oexpd = 0.0
!    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
!    and bw,  therefore, h2oexp for these sub-bands are identical.
  DO k=0,np
!-----xh is the scaled water vapor amount for line absorption
!     computed from Eq. (4.4).
    pwx1 = pa(k)/500.
    pwr1 = pwx1**pm(ib)
    xhd = pwr1*(dh2od(k)*(1.+(aw(ib)+bw(ib)*dt(k))*dt(k))+dh2o(k)*(bw(ib&
&     )*dtd(k)*dt(k)+(aw(ib)+bw(ib)*dt(k))*dtd(k)))
    xh = dh2o(k)*pwr1*(1.+(aw(ib)+bw(ib)*dt(k))*dt(k))
!-----h2oexp is the water vapor transmittance of the layer k
!     due to line absorption
    h2oexpd(k, 1) = -(xkw(ib)*xhd*EXP(-(xh*xkw(ib))))
    h2oexp(k, 1) = EXP(-(xh*xkw(ib)))
!-----compute transmittances from Eq. (8.22)
    DO ik=2,6
      IF (mw(ib) .EQ. 6) THEN
        xhd = h2oexpd(k, ik-1)*h2oexp(k, ik-1) + h2oexp(k, ik-1)*h2oexpd&
&         (k, ik-1)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        h2oexpd(k, ik) = (xhd*xh+xh*xhd)*xh + xh**2*xhd
        h2oexp(k, ik) = xh*xh*xh
      ELSE IF (mw(ib) .EQ. 8) THEN
        xhd = h2oexpd(k, ik-1)*h2oexp(k, ik-1) + h2oexp(k, ik-1)*h2oexpd&
&         (k, ik-1)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        xhd = xhd*xh + xh*xhd
        xh = xh*xh
        h2oexpd(k, ik) = xhd*xh + xh*xhd
        h2oexp(k, ik) = xh*xh
      ELSE IF (mw(ib) .EQ. 9) THEN
        xhd = (h2oexpd(k, ik-1)*h2oexp(k, ik-1)+h2oexp(k, ik-1)*h2oexpd(&
&         k, ik-1))*h2oexp(k, ik-1) + h2oexp(k, ik-1)**2*h2oexpd(k, ik-1&
&         )
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)*h2oexp(k, ik-1)
        h2oexpd(k, ik) = (xhd*xh+xh*xhd)*xh + xh**2*xhd
        h2oexp(k, ik) = xh*xh*xh
      ELSE
        xhd = h2oexpd(k, ik-1)*h2oexp(k, ik-1) + h2oexp(k, ik-1)*h2oexpd&
&         (k, ik-1)
        xh = h2oexp(k, ik-1)*h2oexp(k, ik-1)
        xhd = xhd*xh + xh*xhd
        xh = xh*xh
        xhd = xhd*xh + xh*xhd
        xh = xh*xh
        h2oexpd(k, ik) = xhd*xh + xh*xhd
        h2oexp(k, ik) = xh*xh
      END IF
    END DO
  END DO
END SUBROUTINE H2OEXPS_D

!  Differentiation of conexps in forward (tangent) mode:
!   variations   of useful results: conexp
!   with respect to varying inputs: dcont conexp
SUBROUTINE CONEXPS_D(ib, np, dcont, dcontd, xke, conexp, conexpd)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters ------
  REAL :: dcont(0:np)
  REAL :: dcontd(0:np)
!---- updated parameters -----
  REAL :: conexp(0:np, 3)
  REAL :: conexpd(0:np, 3)
!---- static data -----
  REAL :: xke(9)
  INTRINSIC EXP
  DO k=0,np
    conexpd(k, 1) = -(xke(ib)*dcontd(k)*EXP(-(dcont(k)*xke(ib))))
    conexp(k, 1) = EXP(-(dcont(k)*xke(ib)))
!-----The absorption coefficients for sub-bands 3b and 3a are, respectively,
!     two and four times the absorption coefficient for sub-band 3c (Table 9).
!     Note that conexp(3) is for sub-band 3a. 
    IF (ib .EQ. 3) THEN
      conexpd(k, 2) = conexpd(k, 1)*conexp(k, 1) + conexp(k, 1)*conexpd(&
&       k, 1)
      conexp(k, 2) = conexp(k, 1)*conexp(k, 1)
      conexpd(k, 3) = conexpd(k, 2)*conexp(k, 2) + conexp(k, 2)*conexpd(&
&       k, 2)
      conexp(k, 3) = conexp(k, 2)*conexp(k, 2)
    END IF
  END DO
END SUBROUTINE CONEXPS_D

!  Differentiation of n2oexps in forward (tangent) mode:
!   variations   of useful results: n2oexp
!   with respect to varying inputs: dt n2oexp
SUBROUTINE N2OEXPS_D(ib, np, dn2o, pa, dt, dtd, n2oexp, n2oexpd)
  IMPLICIT NONE
  INTEGER :: ib, k, np
!---- input parameters -----
  REAL :: dn2o(0:np), pa(0:np), dt(0:np)
  REAL :: dtd(0:np)
!---- output parameters -----
  REAL :: n2oexp(0:np, 4)
  REAL :: n2oexpd(0:np, 4)
!---- temporary arrays -----
  REAL :: xc, xc1, xc2
  REAL :: xcd, xc1d, xc2d
  INTRINSIC EXP
!-----Scaling and absorption data are given in Table 5.
!     Transmittances are computed using Eqs. (8.21) and (8.22).
  DO k=0,np
!-----four exponential by powers of 21 for band 6.
    IF (ib .EQ. 6) THEN
      xcd = dn2o(k)*(4.3750e-6*dtd(k)*dt(k)+(1.9297e-3+4.3750e-6*dt(k))*&
&       dtd(k))
      xc = dn2o(k)*(1.+(1.9297e-3+4.3750e-6*dt(k))*dt(k))
      n2oexpd(k, 1) = -(6.31582e-2*xcd*EXP(-(xc*6.31582e-2)))
      n2oexp(k, 1) = EXP(-(xc*6.31582e-2))
      xcd = (n2oexpd(k, 1)*n2oexp(k, 1)+n2oexp(k, 1)*n2oexpd(k, 1))*&
&       n2oexp(k, 1) + n2oexp(k, 1)**2*n2oexpd(k, 1)
      xc = n2oexp(k, 1)*n2oexp(k, 1)*n2oexp(k, 1)
      xc1d = xcd*xc + xc*xcd
      xc1 = xc*xc
      xc2d = xc1d*xc1 + xc1*xc1d
      xc2 = xc1*xc1
      n2oexpd(k, 2) = (xcd*xc1+xc*xc1d)*xc2 + xc*xc1*xc2d
      n2oexp(k, 2) = xc*xc1*xc2
!-----four exponential by powers of 8 for band 7
    ELSE
      xcd = dn2o(k)*(pa(k)/500.0)**0.48*(7.4838e-6*dtd(k)*dt(k)+(&
&       1.3804e-3+7.4838e-6*dt(k))*dtd(k))
      xc = dn2o(k)*(pa(k)/500.0)**0.48*(1.+(1.3804e-3+7.4838e-6*dt(k))*&
&       dt(k))
      n2oexpd(k, 1) = -(5.35779e-2*xcd*EXP(-(xc*5.35779e-2)))
      n2oexp(k, 1) = EXP(-(xc*5.35779e-2))
      xcd = n2oexpd(k, 1)*n2oexp(k, 1) + n2oexp(k, 1)*n2oexpd(k, 1)
      xc = n2oexp(k, 1)*n2oexp(k, 1)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      n2oexpd(k, 2) = xcd*xc + xc*xcd
      n2oexp(k, 2) = xc*xc
      xcd = n2oexpd(k, 2)*n2oexp(k, 2) + n2oexp(k, 2)*n2oexpd(k, 2)
      xc = n2oexp(k, 2)*n2oexp(k, 2)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      n2oexpd(k, 3) = xcd*xc + xc*xcd
      n2oexp(k, 3) = xc*xc
      xcd = n2oexpd(k, 3)*n2oexp(k, 3) + n2oexp(k, 3)*n2oexpd(k, 3)
      xc = n2oexp(k, 3)*n2oexp(k, 3)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      n2oexpd(k, 4) = xcd*xc + xc*xcd
      n2oexp(k, 4) = xc*xc
    END IF
  END DO
END SUBROUTINE N2OEXPS_D

!  Differentiation of ch4exps in forward (tangent) mode:
!   variations   of useful results: ch4exp
!   with respect to varying inputs: dt ch4exp
SUBROUTINE CH4EXPS_D(ib, np, dch4, pa, dt, dtd, ch4exp, ch4expd)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL :: dch4(0:np), pa(0:np), dt(0:np)
  REAL :: dtd(0:np)
!---- output parameters -----
  REAL :: ch4exp(0:np, 4)
  REAL :: ch4expd(0:np, 4)
!---- temporary arrays -----
  REAL :: xc
  REAL :: xcd
  INTRINSIC EXP
!-----  Scaling and absorption data are given in Table 5 
  DO k=0,np
!-----four exponentials for band 6
    IF (ib .EQ. 6) THEN
      xcd = dch4(k)*(1.5826e-4*dtd(k)*dt(k)+(1.7007e-2+1.5826e-4*dt(k))*&
&       dtd(k))
      xc = dch4(k)*(1.+(1.7007e-2+1.5826e-4*dt(k))*dt(k))
      ch4expd(k, 1) = -(5.80708e-3*xcd*EXP(-(xc*5.80708e-3)))
      ch4exp(k, 1) = EXP(-(xc*5.80708e-3))
!-----four exponentials by powers of 12 for band 7
    ELSE
      xcd = dch4(k)*(pa(k)/500.0)**0.65*((5.9590e-4-2.2931e-6*dt(k))*dtd&
&       (k)-2.2931e-6*dtd(k)*dt(k))
      xc = dch4(k)*(pa(k)/500.0)**0.65*(1.+(5.9590e-4-2.2931e-6*dt(k))*&
&       dt(k))
      ch4expd(k, 1) = -(6.29247e-2*xcd*EXP(-(xc*6.29247e-2)))
      ch4exp(k, 1) = EXP(-(xc*6.29247e-2))
      xcd = (ch4expd(k, 1)*ch4exp(k, 1)+ch4exp(k, 1)*ch4expd(k, 1))*&
&       ch4exp(k, 1) + ch4exp(k, 1)**2*ch4expd(k, 1)
      xc = ch4exp(k, 1)*ch4exp(k, 1)*ch4exp(k, 1)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      ch4expd(k, 2) = xcd*xc + xc*xcd
      ch4exp(k, 2) = xc*xc
      xcd = (ch4expd(k, 2)*ch4exp(k, 2)+ch4exp(k, 2)*ch4expd(k, 2))*&
&       ch4exp(k, 2) + ch4exp(k, 2)**2*ch4expd(k, 2)
      xc = ch4exp(k, 2)*ch4exp(k, 2)*ch4exp(k, 2)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      ch4expd(k, 3) = xcd*xc + xc*xcd
      ch4exp(k, 3) = xc*xc
      xcd = (ch4expd(k, 3)*ch4exp(k, 3)+ch4exp(k, 3)*ch4expd(k, 3))*&
&       ch4exp(k, 3) + ch4exp(k, 3)**2*ch4expd(k, 3)
      xc = ch4exp(k, 3)*ch4exp(k, 3)*ch4exp(k, 3)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      ch4expd(k, 4) = xcd*xc + xc*xcd
      ch4exp(k, 4) = xc*xc
    END IF
  END DO
END SUBROUTINE CH4EXPS_D

!  Differentiation of comexps in forward (tangent) mode:
!   variations   of useful results: comexp
!   with respect to varying inputs: dt comexp
SUBROUTINE COMEXPS_D(ib, np, dcom, dt, dtd, comexp, comexpd)
  IMPLICIT NONE
  INTEGER :: ib, ik, np, k
!---- input parameters -----
  REAL :: dcom(0:np), dt(0:np)
  REAL :: dtd(0:np)
!---- output parameters -----
  REAL :: comexp(0:np, 6)
  REAL :: comexpd(0:np, 6)
!---- temporary arrays -----
  REAL :: xc
  REAL :: xcd
  INTRINSIC EXP
  xcd = 0.0
!-----  Scaling and absorpton data are given in Table 6
  DO k=0,np
    IF (ib .EQ. 4) THEN
      xcd = dcom(k)*(4.0447e-4*dtd(k)*dt(k)+(3.5775e-2+4.0447e-4*dt(k))*&
&       dtd(k))
      xc = dcom(k)*(1.+(3.5775e-2+4.0447e-4*dt(k))*dt(k))
    END IF
    IF (ib .EQ. 5) THEN
      xcd = dcom(k)*(3.7401e-4*dtd(k)*dt(k)+(3.4268e-2+3.7401e-4*dt(k))*&
&       dtd(k))
      xc = dcom(k)*(1.+(3.4268e-2+3.7401e-4*dt(k))*dt(k))
    END IF
    comexpd(k, 1) = -(1.922e-7*xcd*EXP(-(xc*1.922e-7)))
    comexp(k, 1) = EXP(-(xc*1.922e-7))
    DO ik=2,6
      xcd = comexpd(k, ik-1)*comexp(k, ik-1) + comexp(k, ik-1)*comexpd(k&
&       , ik-1)
      xc = comexp(k, ik-1)*comexp(k, ik-1)
      xcd = xcd*xc + xc*xcd
      xc = xc*xc
      comexpd(k, ik) = xcd*comexp(k, ik-1) + xc*comexpd(k, ik-1)
      comexp(k, ik) = xc*comexp(k, ik-1)
    END DO
  END DO
END SUBROUTINE COMEXPS_D

!  Differentiation of cfcexps in forward (tangent) mode:
!   variations   of useful results: cfcexp
!   with respect to varying inputs: dt cfcexp
SUBROUTINE CFCEXPS_D(ib, np, a1, b1, fk1, a2, b2, fk2, dcfc, dt, dtd, &
& cfcexp, cfcexpd)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL :: dcfc(0:np), dt(0:np)
  REAL :: dtd(0:np)
!---- output parameters -----
  REAL :: cfcexp(0:np)
  REAL :: cfcexpd(0:np)
!---- static data -----
  REAL :: a1, b1, fk1, a2, b2, fk2
!---- temporary arrays -----
  REAL :: xf
  REAL :: xfd
  INTRINSIC EXP
  DO k=0,np
!-----compute the scaled cfc amount (xf) and exponential (cfcexp)
    IF (ib .EQ. 4) THEN
      xfd = dcfc(k)*(b1*dtd(k)*dt(k)+(a1+b1*dt(k))*dtd(k))
      xf = dcfc(k)*(1.+(a1+b1*dt(k))*dt(k))
      cfcexpd(k) = -(fk1*xfd*EXP(-(xf*fk1)))
      cfcexp(k) = EXP(-(xf*fk1))
    ELSE
      xfd = dcfc(k)*(b2*dtd(k)*dt(k)+(a2+b2*dt(k))*dtd(k))
      xf = dcfc(k)*(1.+(a2+b2*dt(k))*dt(k))
      cfcexpd(k) = -(fk2*xfd*EXP(-(xf*fk2)))
      cfcexp(k) = EXP(-(xf*fk2))
    END IF
  END DO
END SUBROUTINE CFCEXPS_D

!  Differentiation of b10exps in forward (tangent) mode:
!   variations   of useful results: co2exp h2oexp n2oexp conexp
!   with respect to varying inputs: dh2o dcont dt co2exp h2oexp
!                n2oexp conexp
SUBROUTINE B10EXPS_D(np, dh2o, dh2od, dcont, dcontd, dco2, dn2o, pa, dt&
& , dtd, h2oexp, h2oexpd, conexp, conexpd, co2exp, co2expd, n2oexp, &
& n2oexpd)
  IMPLICIT NONE
  INTEGER :: np, k
!---- input parameters -----
  REAL :: dh2o(0:np), dcont(0:np), dn2o(0:np)
  REAL :: dh2od(0:np), dcontd(0:np)
  REAL :: dco2(0:np), pa(0:np), dt(0:np)
  REAL :: dtd(0:np)
!---- output parameters -----
  REAL :: h2oexp(0:np, 5), conexp(0:np), co2exp(0:np, 6), n2oexp(0:np, 2&
& )
  REAL :: h2oexpd(0:np, 5), conexpd(0:np), co2expd(0:np, 6), n2oexpd(0:&
& np, 2)
!---- temporary arrays -----
  REAL :: xx, xx1, xx2, xx3
  REAL :: xxd, xx1d, xx2d, xx3d
  INTRINSIC EXP
  DO k=0,np
!-----Compute scaled h2o-line amount for Band 10 (Eq. 4.4 and Table 3).
    xxd = pa(k)*(dh2od(k)*(1.+(0.0149+6.20e-5*dt(k))*dt(k))+dh2o(k)*(&
&     6.20e-5*dtd(k)*dt(k)+(0.0149+6.20e-5*dt(k))*dtd(k)))/500.0
    xx = dh2o(k)*(pa(k)/500.0)*(1.+(0.0149+6.20e-5*dt(k))*dt(k))
!-----six exponentials by powers of 8
    h2oexpd(k, 1) = -(0.10624*xxd*EXP(-(xx*0.10624)))
    h2oexp(k, 1) = EXP(-(xx*0.10624))
    xxd = h2oexpd(k, 1)*h2oexp(k, 1) + h2oexp(k, 1)*h2oexpd(k, 1)
    xx = h2oexp(k, 1)*h2oexp(k, 1)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    h2oexpd(k, 2) = xxd*xx + xx*xxd
    h2oexp(k, 2) = xx*xx
    xxd = h2oexpd(k, 2)*h2oexp(k, 2) + h2oexp(k, 2)*h2oexpd(k, 2)
    xx = h2oexp(k, 2)*h2oexp(k, 2)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    h2oexpd(k, 3) = xxd*xx + xx*xxd
    h2oexp(k, 3) = xx*xx
    xxd = h2oexpd(k, 3)*h2oexp(k, 3) + h2oexp(k, 3)*h2oexpd(k, 3)
    xx = h2oexp(k, 3)*h2oexp(k, 3)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    h2oexpd(k, 4) = xxd*xx + xx*xxd
    h2oexp(k, 4) = xx*xx
    xxd = h2oexpd(k, 4)*h2oexp(k, 4) + h2oexp(k, 4)*h2oexpd(k, 4)
    xx = h2oexp(k, 4)*h2oexp(k, 4)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    h2oexpd(k, 5) = xxd*xx + xx*xxd
    h2oexp(k, 5) = xx*xx
!-----one exponential of h2o continuum for sub-band 3a (Table 9).
    conexpd(k) = -(109.0*dcontd(k)*EXP(-(dcont(k)*109.0)))
    conexp(k) = EXP(-(dcont(k)*109.0))
!-----Scaled co2 amount for the Band 10 (Eq. 4.4, Tables 3 and 6).
    xxd = dco2(k)*(pa(k)/300.0)**0.5*(1.02e-4*dtd(k)*dt(k)+(0.0179+&
&     1.02e-4*dt(k))*dtd(k))
    xx = dco2(k)*(pa(k)/300.0)**0.5*(1.+(0.0179+1.02e-4*dt(k))*dt(k))
!-----six exponentials by powers of 8
    co2expd(k, 1) = -(2.656e-5*xxd*EXP(-(xx*2.656e-5)))
    co2exp(k, 1) = EXP(-(xx*2.656e-5))
    xxd = co2expd(k, 1)*co2exp(k, 1) + co2exp(k, 1)*co2expd(k, 1)
    xx = co2exp(k, 1)*co2exp(k, 1)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    co2expd(k, 2) = xxd*xx + xx*xxd
    co2exp(k, 2) = xx*xx
    xxd = co2expd(k, 2)*co2exp(k, 2) + co2exp(k, 2)*co2expd(k, 2)
    xx = co2exp(k, 2)*co2exp(k, 2)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    co2expd(k, 3) = xxd*xx + xx*xxd
    co2exp(k, 3) = xx*xx
    xxd = co2expd(k, 3)*co2exp(k, 3) + co2exp(k, 3)*co2expd(k, 3)
    xx = co2exp(k, 3)*co2exp(k, 3)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    co2expd(k, 4) = xxd*xx + xx*xxd
    co2exp(k, 4) = xx*xx
    xxd = co2expd(k, 4)*co2exp(k, 4) + co2exp(k, 4)*co2expd(k, 4)
    xx = co2exp(k, 4)*co2exp(k, 4)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    co2expd(k, 5) = xxd*xx + xx*xxd
    co2exp(k, 5) = xx*xx
    xxd = co2expd(k, 5)*co2exp(k, 5) + co2exp(k, 5)*co2expd(k, 5)
    xx = co2exp(k, 5)*co2exp(k, 5)
    xxd = xxd*xx + xx*xxd
    xx = xx*xx
    co2expd(k, 6) = xxd*xx + xx*xxd
    co2exp(k, 6) = xx*xx
!-----Compute the scaled n2o amount for Band 10 (Table 5).
    xxd = dn2o(k)*(3.6656e-6*dtd(k)*dt(k)+(1.4476e-3+3.6656e-6*dt(k))*&
&     dtd(k))
    xx = dn2o(k)*(1.+(1.4476e-3+3.6656e-6*dt(k))*dt(k))
!-----Two exponentials by powers of 58
    n2oexpd(k, 1) = -(0.25238*xxd*EXP(-(xx*0.25238)))
    n2oexp(k, 1) = EXP(-(xx*0.25238))
    xxd = n2oexpd(k, 1)*n2oexp(k, 1) + n2oexp(k, 1)*n2oexpd(k, 1)
    xx = n2oexp(k, 1)*n2oexp(k, 1)
    xx1d = xxd*xx + xx*xxd
    xx1 = xx*xx
    xx1d = xx1d*xx1 + xx1*xx1d
    xx1 = xx1*xx1
    xx2d = xx1d*xx1 + xx1*xx1d
    xx2 = xx1*xx1
    xx3d = xx2d*xx2 + xx2*xx2d
    xx3 = xx2*xx2
    n2oexpd(k, 2) = (xxd*xx1+xx*xx1d)*xx2*xx3 + xx*xx1*(xx2d*xx3+xx2*&
&     xx3d)
    n2oexp(k, 2) = xx*xx1*xx2*xx3
  END DO
END SUBROUTINE B10EXPS_D

!  Differentiation of tablup in forward (tangent) mode:
!   variations   of useful results: s1 s2 s3 tran
!   with respect to varying inputs: s1 s2 s3 dt dw tran
SUBROUTINE TABLUP_D(nx1, nh1, dw, dwd, p, dt, dtd, s1, s1d, s2, s2d, s3&
& , s3d, w1, p1, dwe, dpe, coef1, coef2, coef3, tran, trand)
  IMPLICIT NONE
  INTEGER :: nx1, nh1
!---- input parameters -----
  REAL :: w1, p1, dwe, dpe
  REAL :: dw, p, dt
  REAL :: dwd, dtd
  REAL :: coef1(nx1, nh1), coef2(nx1, nh1), coef3(nx1, nh1)
!---- update parameter -----
  REAL :: s1, s2, s3, tran
  REAL :: s1d, s2d, s3d, trand
!---- temporary variables -----
  REAL :: we, pe, fw, fp, pa, pb, pc, ax, ba, bb, t1, ca, cb, t2
  REAL :: wed, ped, fwd, fpd, pad, pbd, pcd, axd, bad, bbd, t1d, cad, &
& cbd, t2d
  REAL :: x1, x2, x3, xx, x1c
  REAL :: x1d, x2d, x3d, xxd, x1cd
  INTEGER :: iw, ip
  INTRINSIC LOG10
  INTRINSIC REAL
  INTRINSIC MIN
  INTRINSIC INT
  INTRINSIC MAX
  REAL :: y2
  REAL :: y1
!-----Compute effective pressure (x2) and temperature (x3) following 
!     Eqs. (8.28) and (8.29)
  s1d = s1d + dwd
  s1 = s1 + dw
  s2d = s2d + p*dwd
  s2 = s2 + p*dw
  s3d = s3d + dtd*dw + dt*dwd
  s3 = s3 + dt*dw
  x1d = s1d
  x1 = s1
  x1cd = -(s1d/s1**2)
  x1c = 1.0/s1
  x2d = s2d*x1c + s2*x1cd
  x2 = s2*x1c
  x3d = s3d*x1c + s3*x1cd
  x3 = s3*x1c
!-----normalize we and pe
!       we=(log10(x1)-w1)/dwe
!       pe=(log10(x2)-p1)/dpe
  wed = dwe*x1d/(x1*LOG(10.0))
  we = (LOG10(x1)-w1)*dwe
  ped = dpe*x2d/(x2*LOG(10.0))
  pe = (LOG10(x2)-p1)*dpe
  y1 = REAL(nh1 - 1)
  IF (we .GT. y1) THEN
    we = y1
    wed = 0.0
  ELSE
    we = we
  END IF
  y2 = REAL(nx1 - 1)
  IF (pe .GT. y2) THEN
    pe = y2
    ped = 0.0
  ELSE
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
  fwd = wed
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
  fpd = ped
  fp = pe - REAL(ip - 1)
!-----linear interpolation in pressure
!       pa = coef1(ip,iw-1)*(1.-fp)+coef1(ip+1,iw-1)*fp
!       pb = coef1(ip,  iw)*(1.-fp)+coef1(ip+1,  iw)*fp
!       pc = coef1(ip,iw+1)*(1.-fp)+coef1(ip+1,iw+1)*fp
  pad = (coef1(ip+1, iw-1)-coef1(ip, iw-1))*fpd
  pa = coef1(ip, iw-1) + (coef1(ip+1, iw-1)-coef1(ip, iw-1))*fp
  pbd = (coef1(ip+1, iw)-coef1(ip, iw))*fpd
  pb = coef1(ip, iw) + (coef1(ip+1, iw)-coef1(ip, iw))*fp
  pcd = (coef1(ip+1, iw+1)-coef1(ip, iw+1))*fpd
  pc = coef1(ip, iw+1) + (coef1(ip+1, iw+1)-coef1(ip, iw+1))*fp
!-----quadratic interpolation in absorber amount for coef1
!       ax = (-pa*(1.-fw)+pc*(1.+fw)) *fw*0.5 + pb*(1.-fw*fw)
  axd = 0.5*(((pcd+pad)*fw+(pc+pa)*fwd+pcd-pad)*fw+((pc+pa)*fw+(pc-pa))*&
&   fwd) + pbd*(1.-fw*fw) + pb*(-(fwd*fw)-fw*fwd)
  ax = ((pc+pa)*fw+(pc-pa))*fw*0.5 + pb*(1.-fw*fw)
!-----linear interpolation in absorber amount for coef2 and coef3
!       ba = coef2(ip,  iw)*(1.-fp)+coef2(ip+1,  iw)*fp
!       bb = coef2(ip,iw+1)*(1.-fp)+coef2(ip+1,iw+1)*fp
!       t1 = ba*(1.-fw) + bb*fw
  bad = (coef2(ip+1, iw)-coef2(ip, iw))*fpd
  ba = coef2(ip, iw) + (coef2(ip+1, iw)-coef2(ip, iw))*fp
  bbd = (coef2(ip+1, iw+1)-coef2(ip, iw+1))*fpd
  bb = coef2(ip, iw+1) + (coef2(ip+1, iw+1)-coef2(ip, iw+1))*fp
  t1d = bad + (bbd-bad)*fw + (bb-ba)*fwd
  t1 = ba + (bb-ba)*fw
!       ca = coef3(ip,  iw)*(1.-fp)+coef3(ip+1,  iw)*fp
!       cb = coef3(ip,iw+1)*(1.-fp)+coef3(ip+1,iw+1)*fp
!       t2 = ca*(1.-fw) + cb*fw
  cad = (coef3(ip+1, iw)-coef3(ip, iw))*fpd
  ca = coef3(ip, iw) + (coef3(ip+1, iw)-coef3(ip, iw))*fp
  cbd = (coef3(ip+1, iw+1)-coef3(ip, iw+1))*fpd
  cb = coef3(ip, iw+1) + (coef3(ip+1, iw+1)-coef3(ip, iw+1))*fp
  t2d = cad + (cbd-cad)*fw + (cb-ca)*fwd
  t2 = ca + (cb-ca)*fw
!-----update the total transmittance between levels k1 and k2
  xxd = axd + (t1d+t2d*x3+t2*x3d)*x3 + (t1+t2*x3)*x3d
  xx = ax + (t1+t2*x3)*x3
  IF (xx .GT. 0.9999999) THEN
    xx = 0.9999999
    xxd = 0.0
  ELSE
    xx = xx
  END IF
  IF (xx .LT. 0.0000001) THEN
    xx = 0.0000001
    xxd = 0.0
  ELSE
    xx = xx
  END IF
  trand = trand*xx + tran*xxd
  tran = tran*xx
END SUBROUTINE TABLUP_D

!  Differentiation of h2okdis in forward (tangent) mode:
!   variations   of useful results: tran tcon th2o
!   with respect to varying inputs: tcon h2oexp th2o conexp
SUBROUTINE H2OKDIS_D(ib, np, k, fkw, gkw, ne, h2oexp, h2oexpd, conexp, &
& conexpd, th2o, th2od, tcon, tcond, tran, trand)
  IMPLICIT NONE
!---- input parameters ------
  INTEGER :: ib, ne, np, k
  REAL :: h2oexp(0:np, 6), conexp(0:np, 3)
  REAL :: h2oexpd(0:np, 6), conexpd(0:np, 3)
  REAL :: fkw(6, 9), gkw(6, 3)
!---- updated parameters -----
  REAL :: th2o(6), tcon(3), tran
  REAL :: th2od(6), tcond(3), trand
!---- temporary arrays -----
  REAL :: trnth2o
  REAL :: trnth2od
  INTEGER :: i
!-----tco2 are the six exp factors between levels k1 and k2 
!     tran is the updated total transmittance between levels k1 and k2
!-----th2o is the 6 exp factors between levels k1 and k2 due to
!     h2o line absorption. 
!-----tcon is the 3 exp factors between levels k1 and k2 due to
!     h2o continuum absorption.
!-----trnth2o is the total transmittance between levels k1 and k2 due
!     to both line and continuum absorption.
!-----Compute th2o following Eq. (8.23).
  th2od(1) = th2od(1)*h2oexp(k, 1) + th2o(1)*h2oexpd(k, 1)
  th2o(1) = th2o(1)*h2oexp(k, 1)
  th2od(2) = th2od(2)*h2oexp(k, 2) + th2o(2)*h2oexpd(k, 2)
  th2o(2) = th2o(2)*h2oexp(k, 2)
  th2od(3) = th2od(3)*h2oexp(k, 3) + th2o(3)*h2oexpd(k, 3)
  th2o(3) = th2o(3)*h2oexp(k, 3)
  th2od(4) = th2od(4)*h2oexp(k, 4) + th2o(4)*h2oexpd(k, 4)
  th2o(4) = th2o(4)*h2oexp(k, 4)
  th2od(5) = th2od(5)*h2oexp(k, 5) + th2o(5)*h2oexpd(k, 5)
  th2o(5) = th2o(5)*h2oexp(k, 5)
  th2od(6) = th2od(6)*h2oexp(k, 6) + th2o(6)*h2oexpd(k, 6)
  th2o(6) = th2o(6)*h2oexp(k, 6)
  IF (ne .EQ. 0) THEN
!-----Compute trnh2o following Eq. (8.25). fkw is given in Table 4.
    trnth2od = fkw(1, ib)*th2od(1) + fkw(2, ib)*th2od(2) + fkw(3, ib)*&
&     th2od(3) + fkw(4, ib)*th2od(4) + fkw(5, ib)*th2od(5) + fkw(6, ib)*&
&     th2od(6)
    trnth2o = fkw(1, ib)*th2o(1) + fkw(2, ib)*th2o(2) + fkw(3, ib)*th2o(&
&     3) + fkw(4, ib)*th2o(4) + fkw(5, ib)*th2o(5) + fkw(6, ib)*th2o(6)
    trand = tran*trnth2od
    tran = tran*trnth2o
  ELSE IF (ne .EQ. 1) THEN
!-----Compute trnh2o following Eqs. (8.25) and (4.27).
    tcond(1) = tcond(1)*conexp(k, 1) + tcon(1)*conexpd(k, 1)
    tcon(1) = tcon(1)*conexp(k, 1)
    trnth2od = (fkw(1, ib)*th2od(1)+fkw(2, ib)*th2od(2)+fkw(3, ib)*th2od&
&     (3)+fkw(4, ib)*th2od(4)+fkw(5, ib)*th2od(5)+fkw(6, ib)*th2od(6))*&
&     tcon(1) + (fkw(1, ib)*th2o(1)+fkw(2, ib)*th2o(2)+fkw(3, ib)*th2o(3&
&     )+fkw(4, ib)*th2o(4)+fkw(5, ib)*th2o(5)+fkw(6, ib)*th2o(6))*tcond(&
&     1)
    trnth2o = (fkw(1, ib)*th2o(1)+fkw(2, ib)*th2o(2)+fkw(3, ib)*th2o(3)+&
&     fkw(4, ib)*th2o(4)+fkw(5, ib)*th2o(5)+fkw(6, ib)*th2o(6))*tcon(1)
    trand = tran*trnth2od
    tran = tran*trnth2o
  ELSE
!-----For band 3. This band is divided into 3 subbands.
    tcond(1) = tcond(1)*conexp(k, 1) + tcon(1)*conexpd(k, 1)
    tcon(1) = tcon(1)*conexp(k, 1)
    tcond(2) = tcond(2)*conexp(k, 2) + tcon(2)*conexpd(k, 2)
    tcon(2) = tcon(2)*conexp(k, 2)
    tcond(3) = tcond(3)*conexp(k, 3) + tcon(3)*conexpd(k, 3)
    tcon(3) = tcon(3)*conexp(k, 3)
!-----Compute trnh2o following Eqs. (4.29) and (8.25).
    trnth2od = (gkw(1, 1)*th2od(1)+gkw(2, 1)*th2od(2)+gkw(3, 1)*th2od(3)&
&     +gkw(4, 1)*th2od(4)+gkw(5, 1)*th2od(5)+gkw(6, 1)*th2od(6))*tcon(1)&
&     + (gkw(1, 1)*th2o(1)+gkw(2, 1)*th2o(2)+gkw(3, 1)*th2o(3)+gkw(4, 1)&
&     *th2o(4)+gkw(5, 1)*th2o(5)+gkw(6, 1)*th2o(6))*tcond(1) + (gkw(1, 2&
&     )*th2od(1)+gkw(2, 2)*th2od(2)+gkw(3, 2)*th2od(3)+gkw(4, 2)*th2od(4&
&     )+gkw(5, 2)*th2od(5)+gkw(6, 2)*th2od(6))*tcon(2) + (gkw(1, 2)*th2o&
&     (1)+gkw(2, 2)*th2o(2)+gkw(3, 2)*th2o(3)+gkw(4, 2)*th2o(4)+gkw(5, 2&
&     )*th2o(5)+gkw(6, 2)*th2o(6))*tcond(2) + (gkw(1, 3)*th2od(1)+gkw(2&
&     , 3)*th2od(2)+gkw(3, 3)*th2od(3)+gkw(4, 3)*th2od(4)+gkw(5, 3)*&
&     th2od(5)+gkw(6, 3)*th2od(6))*tcon(3) + (gkw(1, 3)*th2o(1)+gkw(2, 3&
&     )*th2o(2)+gkw(3, 3)*th2o(3)+gkw(4, 3)*th2o(4)+gkw(5, 3)*th2o(5)+&
&     gkw(6, 3)*th2o(6))*tcond(3)
    trnth2o = (gkw(1, 1)*th2o(1)+gkw(2, 1)*th2o(2)+gkw(3, 1)*th2o(3)+gkw&
&     (4, 1)*th2o(4)+gkw(5, 1)*th2o(5)+gkw(6, 1)*th2o(6))*tcon(1) + (gkw&
&     (1, 2)*th2o(1)+gkw(2, 2)*th2o(2)+gkw(3, 2)*th2o(3)+gkw(4, 2)*th2o(&
&     4)+gkw(5, 2)*th2o(5)+gkw(6, 2)*th2o(6))*tcon(2) + (gkw(1, 3)*th2o(&
&     1)+gkw(2, 3)*th2o(2)+gkw(3, 3)*th2o(3)+gkw(4, 3)*th2o(4)+gkw(5, 3)&
&     *th2o(5)+gkw(6, 3)*th2o(6))*tcon(3)
    trand = tran*trnth2od
    tran = tran*trnth2o
  END IF
END SUBROUTINE H2OKDIS_D

!  Differentiation of n2okdis in forward (tangent) mode:
!   variations   of useful results: tran tn2o
!   with respect to varying inputs: tran tn2o n2oexp
SUBROUTINE N2OKDIS_D(ib, np, k, n2oexp, n2oexpd, tn2o, tn2od, tran, &
& trand)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL :: n2oexp(0:np, 4)
  REAL :: n2oexpd(0:np, 4)
!---- updated parameters -----
  REAL :: tn2o(4), tran
  REAL :: tn2od(4), trand
!---- temporary arrays -----
  REAL :: xc
  REAL :: xcd
!-----tn2o is computed from Eq. (8.23). 
!     xc is the total n2o transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.
!-----band 6
  IF (ib .EQ. 6) THEN
    tn2od(1) = tn2od(1)*n2oexp(k, 1) + tn2o(1)*n2oexpd(k, 1)
    tn2o(1) = tn2o(1)*n2oexp(k, 1)
    xcd = 0.940414*tn2od(1)
    xc = 0.940414*tn2o(1)
    tn2od(2) = tn2od(2)*n2oexp(k, 2) + tn2o(2)*n2oexpd(k, 2)
    tn2o(2) = tn2o(2)*n2oexp(k, 2)
    xcd = xcd + 0.059586*tn2od(2)
    xc = xc + 0.059586*tn2o(2)
!-----band 7
  ELSE
    tn2od(1) = tn2od(1)*n2oexp(k, 1) + tn2o(1)*n2oexpd(k, 1)
    tn2o(1) = tn2o(1)*n2oexp(k, 1)
    xcd = 0.561961*tn2od(1)
    xc = 0.561961*tn2o(1)
    tn2od(2) = tn2od(2)*n2oexp(k, 2) + tn2o(2)*n2oexpd(k, 2)
    tn2o(2) = tn2o(2)*n2oexp(k, 2)
    xcd = xcd + 0.138707*tn2od(2)
    xc = xc + 0.138707*tn2o(2)
    tn2od(3) = tn2od(3)*n2oexp(k, 3) + tn2o(3)*n2oexpd(k, 3)
    tn2o(3) = tn2o(3)*n2oexp(k, 3)
    xcd = xcd + 0.240670*tn2od(3)
    xc = xc + 0.240670*tn2o(3)
    tn2od(4) = tn2od(4)*n2oexp(k, 4) + tn2o(4)*n2oexpd(k, 4)
    tn2o(4) = tn2o(4)*n2oexp(k, 4)
    xcd = xcd + 0.058662*tn2od(4)
    xc = xc + 0.058662*tn2o(4)
  END IF
  trand = trand*xc + tran*xcd
  tran = tran*xc
END SUBROUTINE N2OKDIS_D

!  Differentiation of ch4kdis in forward (tangent) mode:
!   variations   of useful results: tran tch4
!   with respect to varying inputs: tran ch4exp tch4
SUBROUTINE CH4KDIS_D(ib, np, k, ch4exp, ch4expd, tch4, tch4d, tran, &
& trand)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL :: ch4exp(0:np, 4)
  REAL :: ch4expd(0:np, 4)
!---- updated parameters -----
  REAL :: tch4(4), tran
  REAL :: tch4d(4), trand
!---- temporary arrays -----
  REAL :: xc
  REAL :: xcd
!-----tch4 is computed from Eq. (8.23). 
!     xc is the total ch4 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 5.
!-----band 6
  IF (ib .EQ. 6) THEN
    tch4d(1) = tch4d(1)*ch4exp(k, 1) + tch4(1)*ch4expd(k, 1)
    tch4(1) = tch4(1)*ch4exp(k, 1)
    xcd = tch4d(1)
    xc = tch4(1)
!-----band 7
  ELSE
    tch4d(1) = tch4d(1)*ch4exp(k, 1) + tch4(1)*ch4expd(k, 1)
    tch4(1) = tch4(1)*ch4exp(k, 1)
    xcd = 0.610650*tch4d(1)
    xc = 0.610650*tch4(1)
    tch4d(2) = tch4d(2)*ch4exp(k, 2) + tch4(2)*ch4expd(k, 2)
    tch4(2) = tch4(2)*ch4exp(k, 2)
    xcd = xcd + 0.280212*tch4d(2)
    xc = xc + 0.280212*tch4(2)
    tch4d(3) = tch4d(3)*ch4exp(k, 3) + tch4(3)*ch4expd(k, 3)
    tch4(3) = tch4(3)*ch4exp(k, 3)
    xcd = xcd + 0.107349*tch4d(3)
    xc = xc + 0.107349*tch4(3)
    tch4d(4) = tch4d(4)*ch4exp(k, 4) + tch4(4)*ch4expd(k, 4)
    tch4(4) = tch4(4)*ch4exp(k, 4)
    xcd = xcd + 0.001789*tch4d(4)
    xc = xc + 0.001789*tch4(4)
  END IF
  trand = trand*xc + tran*xcd
  tran = tran*xc
END SUBROUTINE CH4KDIS_D

!  Differentiation of comkdis in forward (tangent) mode:
!   variations   of useful results: tran tcom
!   with respect to varying inputs: tran tcom comexp
SUBROUTINE COMKDIS_D(ib, np, k, comexp, comexpd, tcom, tcomd, tran, &
& trand)
  IMPLICIT NONE
  INTEGER :: ib, np, k
!---- input parameters -----
  REAL :: comexp(0:np, 6)
  REAL :: comexpd(0:np, 6)
!---- updated parameters -----
  REAL :: tcom(6), tran
  REAL :: tcomd(6), trand
!---- temporary arrays -----
  REAL :: xc
  REAL :: xcd
!-----tcom is computed from Eq. (8.23). 
!     xc is the total co2 transmittance computed from (8.25)
!     The k-distribution functions are given in Table 6.
!-----band 4
  IF (ib .EQ. 4) THEN
    tcomd(1) = tcomd(1)*comexp(k, 1) + tcom(1)*comexpd(k, 1)
    tcom(1) = tcom(1)*comexp(k, 1)
    xcd = 0.12159*tcomd(1)
    xc = 0.12159*tcom(1)
    tcomd(2) = tcomd(2)*comexp(k, 2) + tcom(2)*comexpd(k, 2)
    tcom(2) = tcom(2)*comexp(k, 2)
    xcd = xcd + 0.24359*tcomd(2)
    xc = xc + 0.24359*tcom(2)
    tcomd(3) = tcomd(3)*comexp(k, 3) + tcom(3)*comexpd(k, 3)
    tcom(3) = tcom(3)*comexp(k, 3)
    xcd = xcd + 0.24981*tcomd(3)
    xc = xc + 0.24981*tcom(3)
    tcomd(4) = tcomd(4)*comexp(k, 4) + tcom(4)*comexpd(k, 4)
    tcom(4) = tcom(4)*comexp(k, 4)
    xcd = xcd + 0.26427*tcomd(4)
    xc = xc + 0.26427*tcom(4)
    tcomd(5) = tcomd(5)*comexp(k, 5) + tcom(5)*comexpd(k, 5)
    tcom(5) = tcom(5)*comexp(k, 5)
    xcd = xcd + 0.07807*tcomd(5)
    xc = xc + 0.07807*tcom(5)
    tcomd(6) = tcomd(6)*comexp(k, 6) + tcom(6)*comexpd(k, 6)
    tcom(6) = tcom(6)*comexp(k, 6)
    xcd = xcd + 0.04267*tcomd(6)
    xc = xc + 0.04267*tcom(6)
!-----band 5
  ELSE
    tcomd(1) = tcomd(1)*comexp(k, 1) + tcom(1)*comexpd(k, 1)
    tcom(1) = tcom(1)*comexp(k, 1)
    xcd = 0.06869*tcomd(1)
    xc = 0.06869*tcom(1)
    tcomd(2) = tcomd(2)*comexp(k, 2) + tcom(2)*comexpd(k, 2)
    tcom(2) = tcom(2)*comexp(k, 2)
    xcd = xcd + 0.14795*tcomd(2)
    xc = xc + 0.14795*tcom(2)
    tcomd(3) = tcomd(3)*comexp(k, 3) + tcom(3)*comexpd(k, 3)
    tcom(3) = tcom(3)*comexp(k, 3)
    xcd = xcd + 0.19512*tcomd(3)
    xc = xc + 0.19512*tcom(3)
    tcomd(4) = tcomd(4)*comexp(k, 4) + tcom(4)*comexpd(k, 4)
    tcom(4) = tcom(4)*comexp(k, 4)
    xcd = xcd + 0.33446*tcomd(4)
    xc = xc + 0.33446*tcom(4)
    tcomd(5) = tcomd(5)*comexp(k, 5) + tcom(5)*comexpd(k, 5)
    tcom(5) = tcom(5)*comexp(k, 5)
    xcd = xcd + 0.17199*tcomd(5)
    xc = xc + 0.17199*tcom(5)
    tcomd(6) = tcomd(6)*comexp(k, 6) + tcom(6)*comexpd(k, 6)
    tcom(6) = tcom(6)*comexp(k, 6)
    xcd = xcd + 0.08179*tcomd(6)
    xc = xc + 0.08179*tcom(6)
  END IF
  trand = trand*xc + tran*xcd
  tran = tran*xc
END SUBROUTINE COMKDIS_D

!  Differentiation of cfckdis in forward (tangent) mode:
!   variations   of useful results: tran tcfc
!   with respect to varying inputs: tran tcfc cfcexp
SUBROUTINE CFCKDIS_D(np, k, cfcexp, cfcexpd, tcfc, tcfcd, tran, trand)
  IMPLICIT NONE
!---- input parameters -----
  INTEGER :: k, np
  REAL :: cfcexp(0:np)
  REAL :: cfcexpd(0:np)
!---- updated parameters -----
  REAL :: tcfc, tran
  REAL :: tcfcd, trand
!-----tcfc is the exp factors between levels k1 and k2. 
  tcfcd = tcfcd*cfcexp(k) + tcfc*cfcexpd(k)
  tcfc = tcfc*cfcexp(k)
  trand = trand*tcfc + tran*tcfcd
  tran = tran*tcfc
END SUBROUTINE CFCKDIS_D

!  Differentiation of b10kdis in forward (tangent) mode:
!   variations   of useful results: tran tn2o tcon th2o tco2
!   with respect to varying inputs: tn2o co2exp tcon h2oexp n2oexp
!                th2o tco2 conexp
SUBROUTINE B10KDIS_D(np, k, h2oexp, h2oexpd, conexp, conexpd, co2exp, &
& co2expd, n2oexp, n2oexpd, th2o, th2od, tcon, tcond, tco2, tco2d, tn2o&
& , tn2od, tran, trand)
  IMPLICIT NONE
  INTEGER :: np, k
!---- input parameters -----
  REAL :: h2oexp(0:np, 5), conexp(0:np), co2exp(0:np, 6), n2oexp(0:np, 2&
& )
  REAL :: h2oexpd(0:np, 5), conexpd(0:np), co2expd(0:np, 6), n2oexpd(0:&
& np, 2)
!---- updated parameters -----
  REAL :: th2o(6), tcon(3), tco2(6), tn2o(4), tran
  REAL :: th2od(6), tcond(3), tco2d(6), tn2od(4), trand
!---- temporary arrays -----
  REAL :: xx
  REAL :: xxd
!-----For h2o line. The k-distribution functions are given in Table 4.
  th2od(1) = th2od(1)*h2oexp(k, 1) + th2o(1)*h2oexpd(k, 1)
  th2o(1) = th2o(1)*h2oexp(k, 1)
  xxd = 0.3153*th2od(1)
  xx = 0.3153*th2o(1)
  th2od(2) = th2od(2)*h2oexp(k, 2) + th2o(2)*h2oexpd(k, 2)
  th2o(2) = th2o(2)*h2oexp(k, 2)
  xxd = xxd + 0.4604*th2od(2)
  xx = xx + 0.4604*th2o(2)
  th2od(3) = th2od(3)*h2oexp(k, 3) + th2o(3)*h2oexpd(k, 3)
  th2o(3) = th2o(3)*h2oexp(k, 3)
  xxd = xxd + 0.1326*th2od(3)
  xx = xx + 0.1326*th2o(3)
  th2od(4) = th2od(4)*h2oexp(k, 4) + th2o(4)*h2oexpd(k, 4)
  th2o(4) = th2o(4)*h2oexp(k, 4)
  xxd = xxd + 0.0798*th2od(4)
  xx = xx + 0.0798*th2o(4)
  th2od(5) = th2od(5)*h2oexp(k, 5) + th2o(5)*h2oexpd(k, 5)
  th2o(5) = th2o(5)*h2oexp(k, 5)
  xxd = xxd + 0.0119*th2od(5)
  xx = xx + 0.0119*th2o(5)
  trand = xxd
  tran = xx
!-----For h2o continuum. Note that conexp(k,3) is for subband 3a.
  tcond(1) = tcond(1)*conexp(k) + tcon(1)*conexpd(k)
  tcon(1) = tcon(1)*conexp(k)
  trand = trand*tcon(1) + tran*tcond(1)
  tran = tran*tcon(1)
!-----For co2 (Table 6)
  tco2d(1) = tco2d(1)*co2exp(k, 1) + tco2(1)*co2expd(k, 1)
  tco2(1) = tco2(1)*co2exp(k, 1)
  xxd = 0.2673*tco2d(1)
  xx = 0.2673*tco2(1)
  tco2d(2) = tco2d(2)*co2exp(k, 2) + tco2(2)*co2expd(k, 2)
  tco2(2) = tco2(2)*co2exp(k, 2)
  xxd = xxd + 0.2201*tco2d(2)
  xx = xx + 0.2201*tco2(2)
  tco2d(3) = tco2d(3)*co2exp(k, 3) + tco2(3)*co2expd(k, 3)
  tco2(3) = tco2(3)*co2exp(k, 3)
  xxd = xxd + 0.2106*tco2d(3)
  xx = xx + 0.2106*tco2(3)
  tco2d(4) = tco2d(4)*co2exp(k, 4) + tco2(4)*co2expd(k, 4)
  tco2(4) = tco2(4)*co2exp(k, 4)
  xxd = xxd + 0.2409*tco2d(4)
  xx = xx + 0.2409*tco2(4)
  tco2d(5) = tco2d(5)*co2exp(k, 5) + tco2(5)*co2expd(k, 5)
  tco2(5) = tco2(5)*co2exp(k, 5)
  xxd = xxd + 0.0196*tco2d(5)
  xx = xx + 0.0196*tco2(5)
  tco2d(6) = tco2d(6)*co2exp(k, 6) + tco2(6)*co2expd(k, 6)
  tco2(6) = tco2(6)*co2exp(k, 6)
  xxd = xxd + 0.0415*tco2d(6)
  xx = xx + 0.0415*tco2(6)
  trand = trand*xx + tran*xxd
  tran = tran*xx
!-----For n2o (Table 5)
  tn2od(1) = tn2od(1)*n2oexp(k, 1) + tn2o(1)*n2oexpd(k, 1)
  tn2o(1) = tn2o(1)*n2oexp(k, 1)
  xxd = 0.970831*tn2od(1)
  xx = 0.970831*tn2o(1)
  tn2od(2) = tn2od(2)*n2oexp(k, 2) + tn2o(2)*n2oexpd(k, 2)
  tn2o(2) = tn2o(2)*n2oexp(k, 2)
  xxd = xxd + 0.029169*tn2od(2)
  xx = xx + 0.029169*tn2o(2)
  trand = trand*(xx-1.0) + tran*xxd
  tran = tran*(xx-1.0)
END SUBROUTINE B10KDIS_D

!  Differentiation of cldovlp in forward (tangent) mode:
!   variations   of useful results: cldhi cldlw cldmd
!   with respect to varying inputs: cldhi ett cldlw enn cldmd
!mjs
SUBROUTINE CLDOVLP_D(np, k1, k2, ict, icb, icx, ncld, enn, ennd, ett, &
& ettd, cldhi, cldhid, cldmd, cldmdd, cldlw, cldlwd)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: np, k1, k2, ict, icb, icx(0:np), ncld(3)
  REAL, INTENT(IN) :: enn(0:np), ett(0:np)
  REAL, INTENT(IN) :: ennd(0:np), ettd(0:np)
  REAL, INTENT(INOUT) :: cldhi, cldmd, cldlw
  REAL, INTENT(INOUT) :: cldhid, cldmdd, cldlwd
  INTEGER :: i, j, k, km, kx
  km = k2 - 1
  IF (km .LT. ict) THEN
! do high clouds
    kx = ncld(1)
    IF (kx .EQ. 1 .OR. cldhi .EQ. 0.) THEN
      cldhid = ennd(km)
      cldhi = enn(km)
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      cldhi = 0.0
      IF (kx .NE. 0) THEN
        cldhid = 0.0
        DO k=ict-kx,ict-1
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            cldhid = ennd(j) + ettd(j)*cldhi + ett(j)*cldhid
            cldhi = enn(j) + ett(j)*cldhi
          END IF
        END DO
      ELSE
        cldhid = 0.0
      END IF
    END IF
  ELSE IF (km .GE. ict .AND. km .LT. icb) THEN
! do middle clouds
    kx = ncld(2)
    IF (kx .EQ. 1 .OR. cldmd .EQ. 0.) THEN
      cldmdd = ennd(km)
      cldmd = enn(km)
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      cldmd = 0.0
      IF (kx .NE. 0) THEN
        cldmdd = 0.0
        DO k=icb-kx,icb-1
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            cldmdd = ennd(j) + ettd(j)*cldmd + ett(j)*cldmdd
            cldmd = enn(j) + ett(j)*cldmd
          END IF
        END DO
      ELSE
        cldmdd = 0.0
      END IF
    END IF
  ELSE
! do low clouds
    kx = ncld(3)
    IF (kx .EQ. 1 .OR. cldlw .EQ. 0.) THEN
      cldlwd = ennd(km)
      cldlw = enn(km)
    ELSE
!if ( (kx.lt.1 .or. kx.gt.1)  .and. abs(cldhi) .gt. 0.) then
      cldlw = 0.0
      IF (kx .NE. 0) THEN
        cldlwd = 0.0
        DO k=np+1-kx,np
          j = icx(k)
          IF (j .GE. k1 .AND. j .LE. km) THEN
            cldlwd = ennd(j) + ettd(j)*cldlw + ett(j)*cldlwd
            cldlw = enn(j) + ett(j)*cldlw
          END IF
        END DO
      ELSE
        cldlwd = 0.0
      END IF
    END IF
  END IF
END SUBROUTINE CLDOVLP_D

!  Differentiation of getirtau1 in forward (tangent) mode:
!   variations   of useful results: tcldlyr enn
!   with respect to varying inputs: hydromets tcldlyr fcld enn
!                reff
SUBROUTINE GETIRTAU1_D(ib, nlevs, dp, fcld, fcldd, reff, reffd, &
& hydromets, hydrometsd, tcldlyr, tcldlyrd, enn, ennd, aib_ir1, awb_ir1&
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
  REAL, INTENT(IN) :: fcldd(nlevs)
!  Effective radius (microns)
  REAL, INTENT(IN) :: reff(nlevs, 4)
  REAL, INTENT(IN) :: reffd(nlevs, 4)
!  Hydrometeors (kg/kg)
  REAL, INTENT(IN) :: hydromets(nlevs, 4)
  REAL, INTENT(IN) :: hydrometsd(nlevs, 4)
  REAL, INTENT(IN) :: aib_ir1(3, 10), awb_ir1(4, 10), aiw_ir1(4, 10)
  REAL, INTENT(IN) :: aww_ir1(4, 10), aig_ir1(4, 10), awg_ir1(4, 10)
  REAL, INTENT(IN) :: cons_grav
! !OUTPUT PARAMETERS:
!  Flux transmissivity?
  REAL, INTENT(OUT) :: tcldlyr(0:nlevs)
  REAL, INTENT(OUT) :: tcldlyrd(0:nlevs)
!  Flux transmissivity of a cloud layer?
  REAL, INTENT(OUT) :: enn(0:nlevs)
  REAL, INTENT(OUT) :: ennd(0:nlevs)
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
  REAL :: taucld1d, taucld2d, taucld3d, taucld4d
  REAL :: g1, g2, g3, g4, gg
  REAL :: g1d, g2d, g3d, g4d, ggd
  REAL :: w1, w2, w3, w4, ww
  REAL :: w1d, w2d, w3d, w4d, wwd
  REAL :: ff, tauc
  REAL :: ffd, taucd
  REAL :: reff_snow
  REAL :: reff_snowd
  INTRINSIC MIN
  INTRINSIC ABS
  INTRINSIC MAX
  INTRINSIC EXP
  REAL :: pwr1
  REAL :: pwr1d
  REAL :: max1d
  REAL :: abs0
  REAL :: max1
!-----Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
!     Rain optical thickness is set to 0.00307 /(gm/m**2).
!     It is for a specific drop size distribution provided by Q. Fu.
  tcldlyrd(0) = 0.0
  tcldlyr(0) = 1.0
  ennd(0) = 0.0
  enn(0) = 0.0
  DO k=1,nlevs
    IF (reff(k, 1) .LE. 0.0) THEN
      taucld1 = 0.0
      taucld1d = 0.0
    ELSE
      IF (reff(k, 1) .GT. 0.0 .OR. (reff(k, 1) .LT. 0.0 .AND. aib_ir1(3&
&         , ib) .EQ. INT(aib_ir1(3, ib)))) THEN
        pwr1d = aib_ir1(3, ib)*reff(k, 1)**(aib_ir1(3, ib)-1)*reffd(k, 1&
&         )
      ELSE IF (reff(k, 1) .EQ. 0.0 .AND. aib_ir1(3, ib) .EQ. 1.0) THEN
        pwr1d = reffd(k, 1)
      ELSE
        pwr1d = 0.0
      END IF
      pwr1 = reff(k, 1)**aib_ir1(3, ib)
      taucld1d = dp(k)*1.0e3*(hydrometsd(k, 1)*(aib_ir1(1, ib)+aib_ir1(2&
&       , ib)/pwr1)-hydromets(k, 1)*aib_ir1(2, ib)*pwr1d/pwr1**2)/&
&       cons_grav
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*(aib_ir1(1, ib)+&
&       aib_ir1(2, ib)/pwr1)
    END IF
    taucld2d = dp(k)*1.0e3*(hydrometsd(k, 2)*(awb_ir1(1, ib)+(awb_ir1(2&
&     , ib)+(awb_ir1(3, ib)+awb_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reff(&
&     k, 2))+hydromets(k, 2)*((awb_ir1(4, ib)*reffd(k, 2)*reff(k, 2)+(&
&     awb_ir1(3, ib)+awb_ir1(4, ib)*reff(k, 2))*reffd(k, 2))*reff(k, 2)+&
&     (awb_ir1(2, ib)+(awb_ir1(3, ib)+awb_ir1(4, ib)*reff(k, 2))*reff(k&
&     , 2))*reffd(k, 2)))/cons_grav
    taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_ir1(1, ib)+(&
&     awb_ir1(2, ib)+(awb_ir1(3, ib)+awb_ir1(4, ib)*reff(k, 2))*reff(k, &
&     2))*reff(k, 2))
    taucld3d = 0.00307*dp(k)*1.0e3*hydrometsd(k, 3)/cons_grav
    taucld3 = 0.00307*(dp(k)*1.0e3/cons_grav*hydromets(k, 3))
    IF (reff(k, 4) .GT. 112.0) THEN
      reff_snow = 112.0
      reff_snowd = 0.0
    ELSE
      reff_snowd = reffd(k, 4)
      reff_snow = reff(k, 4)
    END IF
    IF (reff_snow .LE. 0.0) THEN
      taucld4 = 0.0
      taucld4d = 0.0
    ELSE
      IF (reff_snow .GT. 0.0 .OR. (reff_snow .LT. 0.0 .AND. aib_ir1(3, &
&         ib) .EQ. INT(aib_ir1(3, ib)))) THEN
        pwr1d = aib_ir1(3, ib)*reff_snow**(aib_ir1(3, ib)-1)*reff_snowd
      ELSE IF (reff_snow .EQ. 0.0 .AND. aib_ir1(3, ib) .EQ. 1.0) THEN
        pwr1d = reff_snowd
      ELSE
        pwr1d = 0.0
      END IF
      pwr1 = reff_snow**aib_ir1(3, ib)
      taucld4d = dp(k)*1.0e3*(hydrometsd(k, 4)*(aib_ir1(1, ib)+aib_ir1(2&
&       , ib)/pwr1)-hydromets(k, 4)*aib_ir1(2, ib)*pwr1d/pwr1**2)/&
&       cons_grav
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*(aib_ir1(1, ib)+&
&       aib_ir1(2, ib)/pwr1)
    END IF
!-----Compute cloud single-scattering albedo and asymmetry factor for
!     a mixture of ice particles and liquid drops following 
!     Eqs. (6.5), (6.6), (6.15) and (6.16).
!     Single-scattering albedo and asymmetry factor of rain are set
!     to 0.54 and 0.95, respectively, based on the information provided
!     by Prof. Qiang Fu.
    taucd = taucld1d + taucld2d + taucld3d + taucld4d
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      w1d = taucld1d*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff(k, 1))*reff(k, 1))*reff(k, 1)) + taucld1*((&
&       aiw_ir1(4, ib)*reffd(k, 1)*reff(k, 1)+(aiw_ir1(3, ib)+aiw_ir1(4&
&       , ib)*reff(k, 1))*reffd(k, 1))*reff(k, 1)+(aiw_ir1(2, ib)+(&
&       aiw_ir1(3, ib)+aiw_ir1(4, ib)*reff(k, 1))*reff(k, 1))*reffd(k, 1&
&       ))
      w1 = taucld1*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff(k, 1))*reff(k, 1))*reff(k, 1))
      w2d = taucld2d*(aww_ir1(1, ib)+(aww_ir1(2, ib)+(aww_ir1(3, ib)+&
&       aww_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reff(k, 2)) + taucld2*((&
&       aww_ir1(4, ib)*reffd(k, 2)*reff(k, 2)+(aww_ir1(3, ib)+aww_ir1(4&
&       , ib)*reff(k, 2))*reffd(k, 2))*reff(k, 2)+(aww_ir1(2, ib)+(&
&       aww_ir1(3, ib)+aww_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reffd(k, 2&
&       ))
      w2 = taucld2*(aww_ir1(1, ib)+(aww_ir1(2, ib)+(aww_ir1(3, ib)+&
&       aww_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reff(k, 2))
      w3d = 0.54*taucld3d
      w3 = taucld3*0.54
      w4d = taucld4d*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff_snow)*reff_snow)*reff_snow) + taucld4*((&
&       aiw_ir1(4, ib)*reff_snowd*reff_snow+(aiw_ir1(3, ib)+aiw_ir1(4, &
&       ib)*reff_snow)*reff_snowd)*reff_snow+(aiw_ir1(2, ib)+(aiw_ir1(3&
&       , ib)+aiw_ir1(4, ib)*reff_snow)*reff_snow)*reff_snowd)
      w4 = taucld4*(aiw_ir1(1, ib)+(aiw_ir1(2, ib)+(aiw_ir1(3, ib)+&
&       aiw_ir1(4, ib)*reff_snow)*reff_snow)*reff_snow)
      wwd = ((w1d+w2d+w3d+w4d)*tauc-(w1+w2+w3+w4)*taucd)/tauc**2
      ww = (w1+w2+w3+w4)/tauc
      g1d = w1d*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(&
&       4, ib)*reff(k, 1))*reff(k, 1))*reff(k, 1)) + w1*((aig_ir1(4, ib)&
&       *reffd(k, 1)*reff(k, 1)+(aig_ir1(3, ib)+aig_ir1(4, ib)*reff(k, 1&
&       ))*reffd(k, 1))*reff(k, 1)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+&
&       aig_ir1(4, ib)*reff(k, 1))*reff(k, 1))*reffd(k, 1))
      g1 = w1*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(4&
&       , ib)*reff(k, 1))*reff(k, 1))*reff(k, 1))
      g2d = w2d*(awg_ir1(1, ib)+(awg_ir1(2, ib)+(awg_ir1(3, ib)+awg_ir1(&
&       4, ib)*reff(k, 2))*reff(k, 2))*reff(k, 2)) + w2*((awg_ir1(4, ib)&
&       *reffd(k, 2)*reff(k, 2)+(awg_ir1(3, ib)+awg_ir1(4, ib)*reff(k, 2&
&       ))*reffd(k, 2))*reff(k, 2)+(awg_ir1(2, ib)+(awg_ir1(3, ib)+&
&       awg_ir1(4, ib)*reff(k, 2))*reff(k, 2))*reffd(k, 2))
      g2 = w2*(awg_ir1(1, ib)+(awg_ir1(2, ib)+(awg_ir1(3, ib)+awg_ir1(4&
&       , ib)*reff(k, 2))*reff(k, 2))*reff(k, 2))
      g3d = 0.95*w3d
      g3 = w3*0.95
      g4d = w4d*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(&
&       4, ib)*reff_snow)*reff_snow)*reff_snow) + w4*((aig_ir1(4, ib)*&
&       reff_snowd*reff_snow+(aig_ir1(3, ib)+aig_ir1(4, ib)*reff_snow)*&
&       reff_snowd)*reff_snow+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(4&
&       , ib)*reff_snow)*reff_snow)*reff_snowd)
      g4 = w4*(aig_ir1(1, ib)+(aig_ir1(2, ib)+(aig_ir1(3, ib)+aig_ir1(4&
&       , ib)*reff_snow)*reff_snow)*reff_snow)
      IF (w1 + w2 + w3 + w4 .GE. 0.) THEN
        abs0 = w1 + w2 + w3 + w4
      ELSE
        abs0 = -(w1+w2+w3+w4)
      END IF
      IF (abs0 .GT. 0.0) THEN
        ggd = ((g1d+g2d+g3d+g4d)*(w1+w2+w3+w4)-(g1+g2+g3+g4)*(w1d+w2d+&
&         w3d+w4d))/(w1+w2+w3+w4)**2
        gg = (g1+g2+g3+g4)/(w1+w2+w3+w4)
      ELSE
        gg = 0.5
        ggd = 0.0
      END IF
!-----Parameterization of LW scattering following Eqs. (6.11)
!     and (6.12). 
      ffd = (0.1185*ggd*gg+(0.0076+0.1185*gg)*ggd)*gg + (0.3739+(0.0076+&
&       0.1185*gg)*gg)*ggd
      ff = 0.5 + (0.3739+(0.0076+0.1185*gg)*gg)*gg
      IF (1. - ww*ff .LT. 0.0) THEN
        max1 = 0.0
        max1d = 0.0
      ELSE
        max1d = -(wwd*ff) - ww*ffd
        max1 = 1. - ww*ff
      END IF
!ALT: temporary protection against negative cloud optical thickness
      taucd = max1d*tauc + max1*taucd
      tauc = max1*tauc
!-----compute cloud diffuse transmittance. It is approximated by using 
!     a diffusivity factor of 1.66.
      tcldlyrd(k) = -(1.66*taucd*EXP(-(1.66*tauc)))
      tcldlyr(k) = EXP(-(1.66*tauc))
! N in the documentation (6.13)
      ennd(k) = fcldd(k)*(1.0-tcldlyr(k)) - fcld(k)*tcldlyrd(k)
      enn(k) = fcld(k)*(1.0-tcldlyr(k))
    ELSE
      tcldlyrd(k) = 0.0
      tcldlyr(k) = 1.0
      ennd(k) = 0.0
      enn(k) = 0.0
    END IF
  END DO
END SUBROUTINE GETIRTAU1_D

!  Differentiation of sorad in forward (tangent) mode:
!   variations   of useful results: flx_dev
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
SUBROUTINE SORAD_D(m, np, nb, cosz_dev, pl_dev, ta_dev, ta_devd, wa_dev&
& , wa_devd, oa_dev, oa_devd, co2, cwc_dev, cwc_devd, fcld_dev, &
& fcld_devd, ict, icb, reff_dev, reff_devd, hk_uv, hk_ir, taua_dev, &
& taua_devd, ssaa_dev, ssaa_devd, asya_dev, asya_devd, rsuvbm_dev, &
& rsuvdf_dev, rsirbm_dev, rsirdf_dev, flx_dev, flx_devd, cons_grav, &
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
  REAL :: ta_devd(m, np), wa_devd(m, np), oa_devd(m, np)
  REAL :: cwc_dev(m, np, 4), fcld_dev(m, np), reff_dev(m, np, 4), hk_uv(&
& 5), hk_ir(3, 10)
  REAL :: cwc_devd(m, np, 4), fcld_devd(m, np), reff_devd(m, np, 4)
  REAL :: rsuvbm_dev(m), rsuvdf_dev(m), rsirbm_dev(m), rsirdf_dev(m)
  REAL :: taua_dev(m, np, nb)
  REAL :: taua_devd(m, np, nb)
  REAL :: ssaa_dev(m, np, nb)
  REAL :: ssaa_devd(m, np, nb)
  REAL :: asya_dev(m, np, nb)
  REAL :: asya_devd(m, np, nb)
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
  REAL :: flx_devd(m, np+1)
  REAL :: flxu_dev(m, np+1), flcu_dev(m, np+1)
  REAL :: fdiruv_dev(m), fdifuv_dev(m)
  REAL :: fdirpar_dev(m), fdifpar_dev(m)
  REAL :: fdirir_dev(m), fdifir_dev(m)
  REAL :: flx_sfc_band_dev(m, nband)
!-----temporary arrays
  INTEGER :: i, j, k, l, in, ntop
  REAL :: dp(np), wh(np), oh(np)
  REAL :: dpd(np), whd(np), ohd(np)
  REAL :: scal(np)
  REAL :: scald(np)
  REAL :: swh(np+1), so2(np+1), df(0:np+1)
  REAL :: swhd(np+1), so2d(np+1), dfd(0:np+1)
  REAL :: scal0, wvtoa, o3toa, pa
  REAL :: wvtoad, o3toad
  REAL :: snt, cnt, x, xx4, xtoa
  REAL :: xx4d
  REAL :: dp_pa(np)
  REAL :: dp_pad(np)
!-----parameters for co2 transmission tables
  REAL :: w1, dw, u1, du
  INTEGER :: ib, rc
  REAL :: tauclb(np), tauclf(np), asycl(np)
  REAL :: tauclbd(np), tauclfd(np), asycld(np)
  REAL :: taubeam(np, 4), taudiff(np, 4)
  REAL :: taubeamd(np, 4), taudiffd(np, 4)
  REAL :: fcld_col(np)
  REAL :: fcld_cold(np)
  REAL :: cwc_col(np, 4)
  REAL :: cwc_cold(np, 4)
  REAL :: reff_col(np, 4)
  REAL :: reff_cold(np, 4)
  REAL :: taurs, tauoz, tauwv
  REAL :: tauozd, tauwvd
  REAL :: tausto, ssatau, asysto
  REAL :: taustod, ssataud, asystod
  REAL :: tautob, ssatob, asytob
  REAL :: tautobd, ssatobd, asytobd
  REAL :: tautof, ssatof, asytof
  REAL :: tautofd, ssatofd, asytofd
  REAL :: rr(0:np+1, 2), tt(0:np+1, 2), td(0:np+1, 2)
  REAL :: rrd(0:np+1, 2), ttd(0:np+1, 2), tdd(0:np+1, 2)
  REAL :: rs(0:np+1, 2), ts(0:np+1, 2)
  REAL :: rsd(0:np+1, 2), tsd(0:np+1, 2)
  REAL :: fall(np+1), fclr(np+1), fsdir, fsdif
  REAL :: falld(np+1)
  REAL :: fupa(np+1), fupc(np+1)
  REAL :: cc1, cc2, cc3
  REAL :: cc1d, cc2d, cc3d
  REAL :: rrt, ttt, tdt, rst, tst
  REAL :: rrtd, tttd, tdtd, rstd, tstd
  INTEGER :: iv, ik
  REAL :: ssacl(np)
  REAL :: ssacld(np)
  INTEGER :: im
  INTEGER :: ic, iw
  REAL :: ulog, wlog, dc, dd, x0, x1, x2, y0, y1, y2, du2, dw2
  REAL :: wlogd, ddd, x2d, y2d
  INTEGER :: ih
!if (overcast == true) then
!real :: rra(0:np+1),rxa(0:np+1)
!real :: ttaold,tdaold,rsaold
!real :: ttanew,tdanew,rsanew 
!else
  REAL :: rra(0:np+1, 2, 2), tta(0:np, 2, 2)
  REAL :: rrad(0:np+1, 2, 2), ttad(0:np, 2, 2)
  REAL :: tda(0:np, 2, 2)
  REAL :: tdad(0:np, 2, 2)
  REAL :: rsa(0:np, 2, 2), rxa(0:np+1, 2, 2)
  REAL :: rsad(0:np, 2, 2), rxad(0:np+1, 2, 2)
!endif
  REAL :: flxdn
  REAL :: flxdnd
  REAL :: fdndir, fdndif, fupdif
  REAL :: fdndird, fdndifd, fupdifd
  REAL :: denm, yy
  REAL :: denmd, yyd
  INTEGER :: is
  REAL :: ch, cm, ct
  REAL :: chd, cmd, ctd
  INTEGER :: foundtop
  REAL :: dftop
  REAL :: dftopd
!-----Variables for aerosols
  INTEGER :: ii, jj, irhp1, an
  REAL :: dum
  REAL :: dumd
  INTRINSIC MAX
  INTRINSIC EXP
  INTRINSIC MIN
  INTRINSIC SQRT
  INTRINSIC REAL
  INTRINSIC LOG10
  INTRINSIC INT
  INTRINSIC ABS
  INTRINSIC EPSILON
  REAL :: arg1
  REAL :: arg1d
  REAL :: result1
  REAL :: x6
  REAL :: x5
  REAL :: x4
  REAL :: x3
  REAL :: x4d
  REAL :: abs0
  flx_devd = 0.0
  dfd = 0.0
  swhd = 0.0
  reff_cold = 0.0
  ohd = 0.0
  falld = 0.0
  asycld = 0.0
  ssacld = 0.0
  fcld_cold = 0.0
  cwc_cold = 0.0
  tauclbd = 0.0
  tauclfd = 0.0
  whd = 0.0
run_loop:DO i=1,m
    ntop = 0
    fdndir = 0.0
    fdndif = 0.0
!-----Beginning of sorad code
!-----wvtoa and o3toa are the water vapor and o3 amounts of the region 
!     above the pl(1) level.
!     snt is the secant of the solar zenith angle
    snt = 1.0/cosz_dev(i)
    IF (pl_dev(i, 1) .LT. 1.e-3) THEN
      xtoa = 1.e-3
    ELSE
      xtoa = pl_dev(i, 1)
    END IF
    scal0 = xtoa*(0.5*xtoa/300.)**.8
    o3toad = 1.02*xtoa*466.7*oa_devd(i, 1)
    o3toa = 1.02*oa_dev(i, 1)*xtoa*466.7 + 1.0e-8
    wvtoad = 1.02*scal0*(wa_devd(i, 1)*(1.0+0.00135*(ta_dev(i, 1)-240.))&
&     +wa_dev(i, 1)*0.00135*ta_devd(i, 1))
    wvtoa = 1.02*wa_dev(i, 1)*scal0*(1.0+0.00135*(ta_dev(i, 1)-240.)) + &
&     1.0e-9
    swhd(1) = wvtoad
    swh(1) = wvtoa
    DO k=1,np
!-----compute layer thickness. indices for the surface level and
!     surface layer are np+1 and np, respectively.
      dpd(k) = 0.0
      dp(k) = pl_dev(i, k+1) - pl_dev(i, k)
! dp in pascals
      dp_pad(k) = 0.0
      dp_pa(k) = dp(k)*100.
!-----compute scaled water vapor amount following Eqs. (3.3) and (3.5) 
!     unit is g/cm**2
!
      pa = 0.5*(pl_dev(i, k)+pl_dev(i, k+1))
      scald(k) = 0.0
      scal(k) = dp(k)*(pa/300.)**.8
      whd(k) = 1.02*scal(k)*(wa_devd(i, k)*(1.+0.00135*(ta_dev(i, k)-&
&       240.))+wa_dev(i, k)*0.00135*ta_devd(i, k))
      wh(k) = 1.02*wa_dev(i, k)*scal(k)*(1.+0.00135*(ta_dev(i, k)-240.))&
&       + 1.e-9
      swhd(k+1) = swhd(k) + whd(k)
      swh(k+1) = swh(k) + wh(k)
!-----compute ozone amount, unit is (cm-atm)stp
!     the number 466.7 is the unit conversion factor
!     from g/cm**2 to (cm-atm)stp
      ohd(k) = 1.02*dp(k)*466.7*oa_devd(i, k)
      oh(k) = 1.02*oa_dev(i, k)*dp(k)*466.7 + 1.e-8
!-----Fill the reff, cwc, and fcld for the column
      fcld_cold(k) = fcld_devd(i, k)
      fcld_col(k) = fcld_dev(i, k)
      DO l=1,4
        reff_cold(k, l) = reff_devd(i, k, l)
        reff_col(k, l) = reff_dev(i, k, l)
        cwc_cold(k, l) = cwc_devd(i, k, l)
        cwc_col(k, l) = cwc_dev(i, k, l)
      END DO
    END DO
!-----Initialize temporary arrays to zero to avoid UMR
    rr = 0.0
    tt = 0.0
    td = 0.0
    rs = 0.0
    ts = 0.0
    rra = 0.0
    rxa = 0.0
!if( OVERCAST == .false. ) then
    tta = 0.0
    tda = 0.0
    rsa = 0.0
!endif
!-----initialize fluxes for all-sky (flx), clear-sky (flc), and
!     flux reduction (df)
!
    DO k=1,np+1
      flx_devd(i, k) = 0.0
      flx_dev(i, k) = 0.
      flc_dev(i, k) = 0.
      flxu_dev(i, k) = 0.
      flcu_dev(i, k) = 0.
    END DO
!-----Initialize new per-band surface fluxes
    DO ib=1,nband
      flx_sfc_band_dev(i, ib) = 0.
    END DO
!-----Begin inline of SOLUV
!-----compute solar uv and par fluxes
!-----initialize fdiruv, fdifuv, surface reflectances and transmittances.
!     the reflectance and transmittance of the clear and cloudy portions
!     of a layer are denoted by 1 and 2, respectively.
!     cc is the maximum cloud cover in each of the high, middle, and low
!     cloud groups.
!     1/dsm=1/cos(53) = 1.66
    fdiruv_dev(i) = 0.0
    fdifuv_dev(i) = 0.0
    rrd(np+1, 1) = 0.0
    rr(np+1, 1) = rsuvbm_dev(i)
    rrd(np+1, 2) = 0.0
    rr(np+1, 2) = rsuvbm_dev(i)
    rsd(np+1, 1) = 0.0
    rs(np+1, 1) = rsuvdf_dev(i)
    rsd(np+1, 2) = 0.0
    rs(np+1, 2) = rsuvdf_dev(i)
    tdd(np+1, 1) = 0.0
    td(np+1, 1) = 0.0
    tdd(np+1, 2) = 0.0
    td(np+1, 2) = 0.0
    ttd(np+1, 1) = 0.0
    tt(np+1, 1) = 0.0
    ttd(np+1, 2) = 0.0
    tt(np+1, 2) = 0.0
    tsd(np+1, 1) = 0.0
    ts(np+1, 1) = 0.0
    tsd(np+1, 2) = 0.0
    ts(np+1, 2) = 0.0
    rrd(0, 1) = 0.0
    rr(0, 1) = 0.0
    rrd(0, 2) = 0.0
    rr(0, 2) = 0.0
    rsd(0, 1) = 0.0
    rs(0, 1) = 0.0
    rsd(0, 2) = 0.0
    rs(0, 2) = 0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
    ttd(0, 1) = 0.0
    tt(0, 1) = 1.0
    ttd(0, 2) = 0.0
    tt(0, 2) = 1.0
    tsd(0, 1) = 0.0
    ts(0, 1) = 1.0
    tsd(0, 2) = 0.0
    ts(0, 2) = 1.0
    cc1 = 0.0
    cc2 = 0.0
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
    CALL GETVISTAU1_D(np, cosz_dev(i), dp_pa, fcld_col, fcld_cold, &
&               reff_col, reff_cold, cwc_col, cwc_cold, ict, icb, &
&               taubeam, taubeamd, taudiff, taudiffd, asycl, asycld, &
&               aig_uv, awg_uv, arg_uv, aib_uv, awb_uv, arb_uv, aib_nir&
&               , awb_nir, arb_nir, aia_nir, awa_nir, ara_nir, aig_nir, &
&               awg_nir, arg_nir, caib, caif, cons_grav)
    cc1d = 0.0
    cc2d = 0.0
    cc3d = 0.0
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
          cc1d = fcld_devd(i, k)
          cc1 = fcld_dev(i, k)
        ELSE
          cc1 = cc1
        END IF
      ELSE IF (k .LT. icb) THEN
        IF (cc2 .LT. fcld_dev(i, k)) THEN
          cc2d = fcld_devd(i, k)
          cc2 = fcld_dev(i, k)
        ELSE
          cc2 = cc2
        END IF
      ELSE IF (cc3 .LT. fcld_dev(i, k)) THEN
        cc3d = fcld_devd(i, k)
        cc3 = fcld_dev(i, k)
      ELSE
        cc3 = cc3
      END IF
    END DO
!MAT---DO NOT FUSE THIS LOOP
!endif !overcast
    DO k=1,np
      tauclbd(k) = taubeamd(k, 1) + taubeamd(k, 2) + taubeamd(k, 3) + &
&       taubeamd(k, 4)
      tauclb(k) = taubeam(k, 1) + taubeam(k, 2) + taubeam(k, 3) + &
&       taubeam(k, 4)
      tauclfd(k) = taudiffd(k, 1) + taudiffd(k, 2) + taudiffd(k, 3) + &
&       taudiffd(k, 4)
      tauclf(k) = taudiff(k, 1) + taudiff(k, 2) + taudiff(k, 3) + &
&       taudiff(k, 4)
    END DO
    tsd = 0.0
    ttd = 0.0
    rsad = 0.0
    rrd = 0.0
    rsd = 0.0
    ttad = 0.0
    rxad = 0.0
    tdad = 0.0
    rrad = 0.0
    tdd = 0.0
!-----integration over spectral bands
!-----Compute optical thickness, single-scattering albedo and asymmetry
!     factor for a mixture of "na" aerosol types. [Eqs. (4.16)-(4.18)]
    DO ib=1,nband_uv
!-----compute direct beam transmittances of the layer above pl(1)
      arg1d = -((wk_uv(ib)*wvtoad+zk_uv(ib)*o3toad)/cosz_dev(i))
      arg1 = -((wvtoa*wk_uv(ib)+o3toa*zk_uv(ib))/cosz_dev(i))
      tdd(0, 1) = arg1d*EXP(arg1)
      td(0, 1) = EXP(arg1)
      tdd(0, 2) = tdd(0, 1)
      td(0, 2) = td(0, 1)
      DO k=1,np
!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor (Eqs. 6.2-6.4)
        taurs = ry_uv(ib)*dp(k)
        tauozd = zk_uv(ib)*ohd(k)
        tauoz = zk_uv(ib)*oh(k)
        tauwvd = wk_uv(ib)*whd(k)
        tauwv = wk_uv(ib)*wh(k)
        taustod = tauozd + tauwvd + taua_devd(i, k, ib)
        tausto = taurs + tauoz + tauwv + taua_dev(i, k, ib) + 1.0e-7
        ssataud = ssaa_devd(i, k, ib)
        ssatau = ssaa_dev(i, k, ib) + taurs
        asystod = asya_devd(i, k, ib)
        asysto = asya_dev(i, k, ib)
        tautobd = taustod
        tautob = tausto
        asytobd = (asystod*ssatau-asysto*ssataud)/ssatau**2
        asytob = asysto/ssatau
        ssatobd = (ssataud*tautob-ssatau*tautobd)/tautob**2
        ssatob = ssatau/tautob + 1.0e-8
        IF (ssatob .GT. 0.999999) THEN
          ssatob = 0.999999
          ssatobd = 0.0
        ELSE
          ssatob = ssatob
        END IF
!-----for direct incident radiation
        CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, asytobd&
&               , cosz_dev(i), rrt, rrtd, ttt, tttd, tdt, tdtd)
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
        CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, asytobd&
&               , dsm, rst, rstd, tst, tstd, dum, dumd)
        rrd(k, 1) = rrtd
        rr(k, 1) = rrt
        ttd(k, 1) = tttd
        tt(k, 1) = ttt
        tdd(k, 1) = tdtd
        td(k, 1) = tdt
        rsd(k, 1) = rstd
        rs(k, 1) = rst
        tsd(k, 1) = tstd
        ts(k, 1) = tst
!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer
!-----for direct incident radiation
!     The effective layer optical properties. Eqs. (6.2)-(6.4)
        tautobd = taustod + tauclbd(k)
        tautob = tausto + tauclb(k)
        ssatobd = ((ssataud+tauclbd(k))*tautob-(ssatau+tauclb(k))*&
&         tautobd)/tautob**2
        ssatob = (ssatau+tauclb(k))/tautob + 1.0e-8
        IF (ssatob .GT. 0.999999) THEN
          ssatob = 0.999999
          ssatobd = 0.0
        ELSE
          ssatob = ssatob
        END IF
        asytobd = ((asystod+asycld(k)*tauclb(k)+asycl(k)*tauclbd(k))*&
&         ssatob*tautob-(asysto+asycl(k)*tauclb(k))*(ssatobd*tautob+&
&         ssatob*tautobd))/(ssatob*tautob)**2
        asytob = (asysto+asycl(k)*tauclb(k))/(ssatob*tautob)
!-----for diffuse incident radiation
        tautofd = taustod + tauclfd(k)
        tautof = tausto + tauclf(k)
        ssatofd = ((ssataud+tauclfd(k))*tautof-(ssatau+tauclf(k))*&
&         tautofd)/tautof**2
        ssatof = (ssatau+tauclf(k))/tautof + 1.0e-8
        IF (ssatof .GT. 0.999999) THEN
          ssatof = 0.999999
          ssatofd = 0.0
        ELSE
          ssatof = ssatof
        END IF
        asytofd = ((asystod+asycld(k)*tauclf(k)+asycl(k)*tauclfd(k))*&
&         ssatof*tautof-(asysto+asycl(k)*tauclf(k))*(ssatofd*tautof+&
&         ssatof*tautofd))/(ssatof*tautof)**2
        asytof = (asysto+asycl(k)*tauclf(k))/(ssatof*tautof)
!-----for direct incident radiation
!     note that the cloud optical thickness is scaled differently 
!     for direct and diffuse insolation, Eqs. (7.3) and (7.4).
        CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, asytobd&
&               , cosz_dev(i), rrt, rrtd, ttt, tttd, tdt, tdtd)
!-----diffuse incident radiation is approximated by beam radiation 
!     with an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
        CALL DELEDD_D(tautof, tautofd, ssatof, ssatofd, asytof, asytofd&
&               , dsm, rst, rstd, tst, tstd, dum, dumd)
        rrd(k, 2) = rrtd
        rr(k, 2) = rrt
        ttd(k, 2) = tttd
        tt(k, 2) = ttt
        tdd(k, 2) = tdtd
        td(k, 2) = tdt
        rsd(k, 2) = rstd
        rs(k, 2) = rst
        tsd(k, 2) = tstd
        ts(k, 2) = tst
      END DO
!-----flux calculations
!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)
      DO k=1,np+1
        fclr(k) = 0.0
        falld(k) = 0.0
        fall(k) = 0.0
        fupa(k) = 0.0
        fupc(k) = 0.0
      END DO
      fsdir = 0.0
      fsdif = 0.0
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
        tdad(0, ih, 1) = tdd(0, ih)
        tda(0, ih, 1) = td(0, ih)
        ttad(0, ih, 1) = ttd(0, ih)
        tta(0, ih, 1) = tt(0, ih)
        rsad(0, ih, 1) = rsd(0, ih)
        rsa(0, ih, 1) = rs(0, ih)
        tdad(0, ih, 2) = tdd(0, ih)
        tda(0, ih, 2) = td(0, ih)
        ttad(0, ih, 2) = ttd(0, ih)
        tta(0, ih, 2) = tt(0, ih)
        rsad(0, ih, 2) = rsd(0, ih)
        rsa(0, ih, 2) = rs(0, ih)
        DO k=1,ict-1
          denmd = (tsd(k, ih)*(1.-rsa(k-1, ih, 1)*rs(k, ih))-ts(k, ih)*(&
&           -(rsad(k-1, ih, 1)*rs(k, ih))-rsa(k-1, ih, 1)*rsd(k, ih)))/(&
&           1.-rsa(k-1, ih, 1)*rs(k, ih))**2
          denm = ts(k, ih)/(1.-rsa(k-1, ih, 1)*rs(k, ih))
          tdad(k, ih, 1) = tdad(k-1, ih, 1)*td(k, ih) + tda(k-1, ih, 1)*&
&           tdd(k, ih)
          tda(k, ih, 1) = tda(k-1, ih, 1)*td(k, ih)
          ttad(k, ih, 1) = tdad(k-1, ih, 1)*tt(k, ih) + tda(k-1, ih, 1)*&
&           ttd(k, ih) + ((tdad(k-1, ih, 1)*rr(k, ih)+tda(k-1, ih, 1)*&
&           rrd(k, ih))*rsa(k-1, ih, 1)+tda(k-1, ih, 1)*rr(k, ih)*rsad(k&
&           -1, ih, 1)+ttad(k-1, ih, 1)-tdad(k-1, ih, 1))*denm + (tda(k-&
&           1, ih, 1)*rsa(k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)-tda(k-1&
&           , ih, 1))*denmd
          tta(k, ih, 1) = tda(k-1, ih, 1)*tt(k, ih) + (tda(k-1, ih, 1)*&
&           rsa(k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)-tda(k-1, ih, 1))*&
&           denm
          rsad(k, ih, 1) = rsd(k, ih) + (tsd(k, ih)*denm+ts(k, ih)*denmd&
&           )*rsa(k-1, ih, 1) + ts(k, ih)*denm*rsad(k-1, ih, 1)
          rsa(k, ih, 1) = rs(k, ih) + ts(k, ih)*rsa(k-1, ih, 1)*denm
          tdad(k, ih, 2) = tdad(k, ih, 1)
          tda(k, ih, 2) = tda(k, ih, 1)
          ttad(k, ih, 2) = ttad(k, ih, 1)
          tta(k, ih, 2) = tta(k, ih, 1)
          rsad(k, ih, 2) = rsad(k, ih, 1)
          rsa(k, ih, 2) = rsa(k, ih, 1)
        END DO
! k loop
!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition
        DO k=ict,icb-1
          DO im=1,2
            denmd = (tsd(k, im)*(1.-rsa(k-1, ih, im)*rs(k, im))-ts(k, im&
&             )*(-(rsad(k-1, ih, im)*rs(k, im))-rsa(k-1, ih, im)*rsd(k, &
&             im)))/(1.-rsa(k-1, ih, im)*rs(k, im))**2
            denm = ts(k, im)/(1.-rsa(k-1, ih, im)*rs(k, im))
            tdad(k, ih, im) = tdad(k-1, ih, im)*td(k, im) + tda(k-1, ih&
&             , im)*tdd(k, im)
            tda(k, ih, im) = tda(k-1, ih, im)*td(k, im)
            ttad(k, ih, im) = tdad(k-1, ih, im)*tt(k, im) + tda(k-1, ih&
&             , im)*ttd(k, im) + ((tdad(k-1, ih, im)*rr(k, im)+tda(k-1, &
&             ih, im)*rrd(k, im))*rsa(k-1, ih, im)+tda(k-1, ih, im)*rr(k&
&             , im)*rsad(k-1, ih, im)+ttad(k-1, ih, im)-tdad(k-1, ih, im&
&             ))*denm + (tda(k-1, ih, im)*rsa(k-1, ih, im)*rr(k, im)+tta&
&             (k-1, ih, im)-tda(k-1, ih, im))*denmd
            tta(k, ih, im) = tda(k-1, ih, im)*tt(k, im) + (tda(k-1, ih, &
&             im)*rsa(k-1, ih, im)*rr(k, im)+tta(k-1, ih, im)-tda(k-1, &
&             ih, im))*denm
            rsad(k, ih, im) = rsd(k, im) + (tsd(k, im)*denm+ts(k, im)*&
&             denmd)*rsa(k-1, ih, im) + ts(k, im)*denm*rsad(k-1, ih, im)
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
        rrad(np+1, 1, is) = rrd(np+1, is)
        rra(np+1, 1, is) = rr(np+1, is)
        rxad(np+1, 1, is) = rsd(np+1, is)
        rxa(np+1, 1, is) = rs(np+1, is)
        rrad(np+1, 2, is) = rrd(np+1, is)
        rra(np+1, 2, is) = rr(np+1, is)
        rxad(np+1, 2, is) = rsd(np+1, is)
        rxa(np+1, 2, is) = rs(np+1, is)
        DO k=np,icb,-1
          denmd = (tsd(k, is)*(1.-rs(k, is)*rxa(k+1, 1, is))-ts(k, is)*(&
&           -(rsd(k, is)*rxa(k+1, 1, is))-rs(k, is)*rxad(k+1, 1, is)))/(&
&           1.-rs(k, is)*rxa(k+1, 1, is))**2
          denm = ts(k, is)/(1.-rs(k, is)*rxa(k+1, 1, is))
          rrad(k, 1, is) = rrd(k, is) + (tdd(k, is)*rra(k+1, 1, is)+td(k&
&           , is)*rrad(k+1, 1, is)+(ttd(k, is)-tdd(k, is))*rxa(k+1, 1, &
&           is)+(tt(k, is)-td(k, is))*rxad(k+1, 1, is))*denm + (td(k, is&
&           )*rra(k+1, 1, is)+(tt(k, is)-td(k, is))*rxa(k+1, 1, is))*&
&           denmd
          rra(k, 1, is) = rr(k, is) + (td(k, is)*rra(k+1, 1, is)+(tt(k, &
&           is)-td(k, is))*rxa(k+1, 1, is))*denm
          rxad(k, 1, is) = rsd(k, is) + (tsd(k, is)*denm+ts(k, is)*denmd&
&           )*rxa(k+1, 1, is) + ts(k, is)*denm*rxad(k+1, 1, is)
          rxa(k, 1, is) = rs(k, is) + ts(k, is)*rxa(k+1, 1, is)*denm
          rrad(k, 2, is) = rrad(k, 1, is)
          rra(k, 2, is) = rra(k, 1, is)
          rxad(k, 2, is) = rxad(k, 1, is)
          rxa(k, 2, is) = rxa(k, 1, is)
        END DO
! k loop
!-----for middle clouds
        DO k=icb-1,ict,-1
          DO im=1,2
            denmd = (tsd(k, im)*(1.-rs(k, im)*rxa(k+1, im, is))-ts(k, im&
&             )*(-(rsd(k, im)*rxa(k+1, im, is))-rs(k, im)*rxad(k+1, im, &
&             is)))/(1.-rs(k, im)*rxa(k+1, im, is))**2
            denm = ts(k, im)/(1.-rs(k, im)*rxa(k+1, im, is))
            rrad(k, im, is) = rrd(k, im) + (tdd(k, im)*rra(k+1, im, is)+&
&             td(k, im)*rrad(k+1, im, is)+(ttd(k, im)-tdd(k, im))*rxa(k+&
&             1, im, is)+(tt(k, im)-td(k, im))*rxad(k+1, im, is))*denm +&
&             (td(k, im)*rra(k+1, im, is)+(tt(k, im)-td(k, im))*rxa(k+1&
&             , im, is))*denmd
            rra(k, im, is) = rr(k, im) + (td(k, im)*rra(k+1, im, is)+(tt&
&             (k, im)-td(k, im))*rxa(k+1, im, is))*denm
            rxad(k, im, is) = rsd(k, im) + (tsd(k, im)*denm+ts(k, im)*&
&             denmd)*rxa(k+1, im, is) + ts(k, im)*denm*rxad(k+1, im, is)
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
          chd = -cc1d
          ch = 1.0 - cc1
!-----cloudy portion
        ELSE
          chd = cc1d
          ch = cc1
        END IF
        DO im=1,2
!-----clear portion 
          IF (im .EQ. 1) THEN
            cmd = chd*(1.0-cc2) - ch*cc2d
            cm = ch*(1.0-cc2)
!-----cloudy portion
          ELSE
            cmd = chd*cc2 + ch*cc2d
            cm = ch*cc2
          END IF
          DO is=1,2
!-----clear portion 
            IF (is .EQ. 1) THEN
              ctd = cmd*(1.0-cc3) - cm*cc3d
              ct = cm*(1.0-cc3)
!-----cloudy portion
            ELSE
              ctd = cmd*cc3 + cm*cc3d
              ct = cm*cc3
            END IF
!-----add one layer at a time, going down.
            DO k=icb,np
              denmd = (tsd(k, is)*(1.-rsa(k-1, ih, im)*rs(k, is))-ts(k, &
&               is)*(-(rsad(k-1, ih, im)*rs(k, is))-rsa(k-1, ih, im)*rsd&
&               (k, is)))/(1.-rsa(k-1, ih, im)*rs(k, is))**2
              denm = ts(k, is)/(1.-rsa(k-1, ih, im)*rs(k, is))
              tdad(k, ih, im) = tdad(k-1, ih, im)*td(k, is) + tda(k-1, &
&               ih, im)*tdd(k, is)
              tda(k, ih, im) = tda(k-1, ih, im)*td(k, is)
              ttad(k, ih, im) = tdad(k-1, ih, im)*tt(k, is) + tda(k-1, &
&               ih, im)*ttd(k, is) + ((tdad(k-1, ih, im)*rr(k, is)+tda(k&
&               -1, ih, im)*rrd(k, is))*rsa(k-1, ih, im)+tda(k-1, ih, im&
&               )*rr(k, is)*rsad(k-1, ih, im)+ttad(k-1, ih, im)-tdad(k-1&
&               , ih, im))*denm + (tda(k-1, ih, im)*rr(k, is)*rsa(k-1, &
&               ih, im)+tta(k-1, ih, im)-tda(k-1, ih, im))*denmd
              tta(k, ih, im) = tda(k-1, ih, im)*tt(k, is) + (tda(k-1, ih&
&               , im)*rr(k, is)*rsa(k-1, ih, im)+tta(k-1, ih, im)-tda(k-&
&               1, ih, im))*denm
              rsad(k, ih, im) = rsd(k, is) + (tsd(k, is)*denm+ts(k, is)*&
&               denmd)*rsa(k-1, ih, im) + ts(k, is)*denm*rsad(k-1, ih, &
&               im)
              rsa(k, ih, im) = rs(k, is) + ts(k, is)*rsa(k-1, ih, im)*&
&               denm
            END DO
! k loop
!-----add one layer at a time, going up.
            DO k=ict-1,0,-1
              denmd = (tsd(k, ih)*(1.-rs(k, ih)*rxa(k+1, im, is))-ts(k, &
&               ih)*(-(rsd(k, ih)*rxa(k+1, im, is))-rs(k, ih)*rxad(k+1, &
&               im, is)))/(1.-rs(k, ih)*rxa(k+1, im, is))**2
              denm = ts(k, ih)/(1.-rs(k, ih)*rxa(k+1, im, is))
              rrad(k, im, is) = rrd(k, ih) + (tdd(k, ih)*rra(k+1, im, is&
&               )+td(k, ih)*rrad(k+1, im, is)+(ttd(k, ih)-tdd(k, ih))*&
&               rxa(k+1, im, is)+(tt(k, ih)-td(k, ih))*rxad(k+1, im, is)&
&               )*denm + (td(k, ih)*rra(k+1, im, is)+(tt(k, ih)-td(k, ih&
&               ))*rxa(k+1, im, is))*denmd
              rra(k, im, is) = rr(k, ih) + (td(k, ih)*rra(k+1, im, is)+(&
&               tt(k, ih)-td(k, ih))*rxa(k+1, im, is))*denm
              rxad(k, im, is) = rsd(k, ih) + (tsd(k, ih)*denm+ts(k, ih)*&
&               denmd)*rxa(k+1, im, is) + ts(k, ih)*denm*rxad(k+1, im, &
&               is)
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
              denmd = -((-(rsad(k-1, ih, im)*rxa(k, im, is))-rsa(k-1, ih&
&               , im)*rxad(k, im, is))/(1.-rsa(k-1, ih, im)*rxa(k, im, &
&               is))**2)
              denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
              fdndird = tdad(k-1, ih, im)
              fdndir = tda(k-1, ih, im)
              xx4d = tdad(k-1, ih, im)*rra(k, im, is) + tda(k-1, ih, im)&
&               *rrad(k, im, is)
              xx4 = tda(k-1, ih, im)*rra(k, im, is)
              yyd = ttad(k-1, ih, im) - tdad(k-1, ih, im)
              yy = tta(k-1, ih, im) - tda(k-1, ih, im)
              fdndifd = (xx4d*rsa(k-1, ih, im)+xx4*rsad(k-1, ih, im)+yyd&
&               )*denm + (xx4*rsa(k-1, ih, im)+yy)*denmd
              fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
              fupdifd = (xx4d+yyd*rxa(k, im, is)+yy*rxad(k, im, is))*&
&               denm + (xx4+yy*rxa(k, im, is))*denmd
              fupdif = (xx4+yy*rxa(k, im, is))*denm
              flxdnd = fdndird + fdndifd - fupdifd
              flxdn = fdndir + fdndif - fupdif
!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)
              IF (ih .EQ. 1 .AND. im .EQ. 1 .AND. is .EQ. 1) THEN
                fupc(k) = fupdif
                fclr(k) = flxdn
              END IF
              fupa(k) = fupa(k) + fupdif*ct
              falld(k) = falld(k) + flxdnd*ct + flxdn*ctd
              fall(k) = fall(k) + flxdn*ct
            END DO
! k loop
            fsdir = fsdir + fdndir*ct
            fsdif = fsdif + fdndif*ct
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
        flx_devd(i, k) = flx_devd(i, k) + hk_uv(ib)*falld(k)
        flx_dev(i, k) = flx_dev(i, k) + fall(k)*hk_uv(ib)
        flc_dev(i, k) = flc_dev(i, k) + fclr(k)*hk_uv(ib)
        flxu_dev(i, k) = flxu_dev(i, k) + fupa(k)*hk_uv(ib)
        flcu_dev(i, k) = flcu_dev(i, k) + fupc(k)*hk_uv(ib)
      END DO
!-----get surface flux for each band
      flx_sfc_band_dev(i, ib) = flx_sfc_band_dev(i, ib) + fall(np+1)*&
&       hk_uv(ib)
!-----compute direct and diffuse downward surface fluxes in the UV
!     and par regions
      IF (ib .LT. 5) THEN
        fdiruv_dev(i) = fdiruv_dev(i) + fsdir*hk_uv(ib)
        fdifuv_dev(i) = fdifuv_dev(i) + fsdif*hk_uv(ib)
      ELSE
        fdirpar_dev(i) = fsdir*hk_uv(ib)
        fdifpar_dev(i) = fsdif*hk_uv(ib)
      END IF
    END DO
!-----Inline SOLIR
!-----compute and update solar ir fluxes
    fdirir_dev(i) = 0.0
    fdifir_dev(i) = 0.0
    rrd(np+1, 1) = 0.0
    rr(np+1, 1) = rsirbm_dev(i)
    rrd(np+1, 2) = 0.0
    rr(np+1, 2) = rsirbm_dev(i)
    rsd(np+1, 1) = 0.0
    rs(np+1, 1) = rsirdf_dev(i)
    rsd(np+1, 2) = 0.0
    rs(np+1, 2) = rsirdf_dev(i)
    tdd(np+1, 1) = 0.0
    td(np+1, 1) = 0.0
    tdd(np+1, 2) = 0.0
    td(np+1, 2) = 0.0
    ttd(np+1, 1) = 0.0
    tt(np+1, 1) = 0.0
    ttd(np+1, 2) = 0.0
    tt(np+1, 2) = 0.0
    tsd(np+1, 1) = 0.0
    ts(np+1, 1) = 0.0
    tsd(np+1, 2) = 0.0
    ts(np+1, 2) = 0.0
    rrd(0, 1) = 0.0
    rr(0, 1) = 0.0
    rrd(0, 2) = 0.0
    rr(0, 2) = 0.0
    rsd(0, 1) = 0.0
    rs(0, 1) = 0.0
    rsd(0, 2) = 0.0
    rs(0, 2) = 0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
    ttd(0, 1) = 0.0
    tt(0, 1) = 1.0
    ttd(0, 2) = 0.0
    tt(0, 2) = 1.0
    tsd(0, 1) = 0.0
    ts(0, 1) = 1.0
    tsd(0, 2) = 0.0
    ts(0, 2) = 1.0
    cc1 = 0.0
    cc2 = 0.0
    cc3 = 0.0
    cc1d = 0.0
    cc2d = 0.0
    cc3d = 0.0
!-----integration over spectral bands
!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10)
!     The indices 1, 2, 3 are for ice, water, rain particles,
!     respectively.
    DO ib=1,nband_ir
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
      CALL GETNIRTAU1_D(ib, np, cosz_dev(i), dp_pa, fcld_col, fcld_cold&
&                 , reff_col, reff_cold, cwc_col, cwc_cold, ict, icb, &
&                 taubeam, taubeamd, taudiff, taudiffd, asycl, asycld, &
&                 ssacl, ssacld, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv&
&                 , arb_uv, aib_nir, awb_nir, arb_nir, aia_nir, awa_nir&
&                 , ara_nir, aig_nir, awg_nir, arg_nir, caib, caif, &
&                 cons_grav)
!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
!MAT--DO NOT FUSE THIS LOOP
!MAT  Loop must run to completion so that cc[1,2,3] are correct.
      DO k=1,np
        IF (k .LT. ict) THEN
          IF (cc1 .LT. fcld_dev(i, k)) THEN
            cc1d = fcld_devd(i, k)
            cc1 = fcld_dev(i, k)
          ELSE
            cc1 = cc1
          END IF
        ELSE IF (k .LT. icb) THEN
          IF (cc2 .LT. fcld_dev(i, k)) THEN
            cc2d = fcld_devd(i, k)
            cc2 = fcld_dev(i, k)
          ELSE
            cc2 = cc2
          END IF
        ELSE IF (cc3 .LT. fcld_dev(i, k)) THEN
          cc3d = fcld_devd(i, k)
          cc3 = fcld_dev(i, k)
        ELSE
          cc3 = cc3
        END IF
      END DO
!MAT--DO NOT FUSE THIS LOOP
!endif !overcast
      DO k=1,np
        tauclbd(k) = taubeamd(k, 1) + taubeamd(k, 2) + taubeamd(k, 3) + &
&         taubeamd(k, 4)
        tauclb(k) = taubeam(k, 1) + taubeam(k, 2) + taubeam(k, 3) + &
&         taubeam(k, 4)
        tauclfd(k) = taudiffd(k, 1) + taudiffd(k, 2) + taudiffd(k, 3) + &
&         taudiffd(k, 4)
        tauclf(k) = taudiff(k, 1) + taudiff(k, 2) + taudiff(k, 3) + &
&         taudiff(k, 4)
      END DO
!-----integration over the k-distribution function
      DO ik=1,nk_ir
!-----compute direct beam transmittances of the layer above pl(1)
        tdd(0, 1) = -(xk_ir(ik)*wvtoad*EXP(-(wvtoa*xk_ir(ik)/cosz_dev(i)&
&         ))/cosz_dev(i))
        td(0, 1) = EXP(-(wvtoa*xk_ir(ik)/cosz_dev(i)))
        tdd(0, 2) = tdd(0, 1)
        td(0, 2) = td(0, 1)
        DO k=1,np
          taurs = ry_ir(ib)*dp(k)
          tauwvd = xk_ir(ik)*whd(k)
          tauwv = xk_ir(ik)*wh(k)
!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor. Eqs.(6.2)-(6.4)
          taustod = tauwvd + taua_devd(i, k, iv)
          tausto = taurs + tauwv + taua_dev(i, k, iv) + 1.0e-7
          ssataud = ssaa_devd(i, k, iv)
          ssatau = ssaa_dev(i, k, iv) + taurs + 1.0e-8
          asystod = asya_devd(i, k, iv)
          asysto = asya_dev(i, k, iv)
          tautobd = taustod
          tautob = tausto
          asytobd = (asystod*ssatau-asysto*ssataud)/ssatau**2
          asytob = asysto/ssatau
          ssatobd = (ssataud*tautob-ssatau*tautobd)/tautob**2
          ssatob = ssatau/tautob + 1.0e-8
          IF (ssatob .GT. 0.999999) THEN
            ssatob = 0.999999
            ssatobd = 0.0
          ELSE
            ssatob = ssatob
          END IF
!-----Compute reflectance and transmittance of the clear portion 
!     of a layer
!-----for direct incident radiation
          CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, &
&                 asytobd, cosz_dev(i), rrt, rrtd, ttt, tttd, tdt, tdtd)
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
          CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, &
&                 asytobd, dsm, rst, rstd, tst, tstd, dum, dumd)
          rrd(k, 1) = rrtd
          rr(k, 1) = rrt
          ttd(k, 1) = tttd
          tt(k, 1) = ttt
          tdd(k, 1) = tdtd
          td(k, 1) = tdt
          rsd(k, 1) = rstd
          rs(k, 1) = rst
          tsd(k, 1) = tstd
          ts(k, 1) = tst
!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer
!-----for direct incident radiation. Eqs.(6.2)-(6.4)
          tautobd = taustod + tauclbd(k)
          tautob = tausto + tauclb(k)
          ssatobd = ((ssataud+ssacld(k)*tauclb(k)+ssacl(k)*tauclbd(k))*&
&           tautob-(ssatau+ssacl(k)*tauclb(k))*tautobd)/tautob**2
          ssatob = (ssatau+ssacl(k)*tauclb(k))/tautob + 1.0e-8
          IF (ssatob .GT. 0.999999) THEN
            ssatob = 0.999999
            ssatobd = 0.0
          ELSE
            ssatob = ssatob
          END IF
          asytobd = ((asystod+(asycld(k)*ssacl(k)+asycl(k)*ssacld(k))*&
&           tauclb(k)+asycl(k)*ssacl(k)*tauclbd(k))*ssatob*tautob-(&
&           asysto+asycl(k)*ssacl(k)*tauclb(k))*(ssatobd*tautob+ssatob*&
&           tautobd))/(ssatob*tautob)**2
          asytob = (asysto+asycl(k)*ssacl(k)*tauclb(k))/(ssatob*tautob)
!-----for diffuse incident radiation
          tautofd = taustod + tauclfd(k)
          tautof = tausto + tauclf(k)
          ssatofd = ((ssataud+ssacld(k)*tauclf(k)+ssacl(k)*tauclfd(k))*&
&           tautof-(ssatau+ssacl(k)*tauclf(k))*tautofd)/tautof**2
          ssatof = (ssatau+ssacl(k)*tauclf(k))/tautof + 1.0e-8
          IF (ssatof .GT. 0.999999) THEN
            ssatof = 0.999999
            ssatofd = 0.0
          ELSE
            ssatof = ssatof
          END IF
          asytofd = ((asystod+(asycld(k)*ssacl(k)+asycl(k)*ssacld(k))*&
&           tauclf(k)+asycl(k)*ssacl(k)*tauclfd(k))*ssatof*tautof-(&
&           asysto+asycl(k)*ssacl(k)*tauclf(k))*(ssatofd*tautof+ssatof*&
&           tautofd))/(ssatof*tautof)**2
          asytof = (asysto+asycl(k)*ssacl(k)*tauclf(k))/(ssatof*tautof)
!-----for direct incident radiation
          CALL DELEDD_D(tautob, tautobd, ssatob, ssatobd, asytob, &
&                 asytobd, cosz_dev(i), rrt, rrtd, ttt, tttd, tdt, tdtd)
!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs.(6.5) and (6.6)
          CALL DELEDD_D(tautof, tautofd, ssatof, ssatofd, asytof, &
&                 asytofd, dsm, rst, rstd, tst, tstd, dum, dumd)
          rrd(k, 2) = rrtd
          rr(k, 2) = rrt
          ttd(k, 2) = tttd
          tt(k, 2) = ttt
          tdd(k, 2) = tdtd
          td(k, 2) = tdt
          rsd(k, 2) = rstd
          rs(k, 2) = rst
          tsd(k, 2) = tstd
          ts(k, 2) = tst
        END DO
!-----FLUX CALCULATIONS
!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)
        DO k=1,np+1
          fclr(k) = 0.0
          falld(k) = 0.0
          fall(k) = 0.0
          fupc(k) = 0.0
          fupa(k) = 0.0
        END DO
        fsdir = 0.0
        fsdif = 0.0
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
          tdad(0, ih, 1) = tdd(0, ih)
          tda(0, ih, 1) = td(0, ih)
          ttad(0, ih, 1) = ttd(0, ih)
          tta(0, ih, 1) = tt(0, ih)
          rsad(0, ih, 1) = rsd(0, ih)
          rsa(0, ih, 1) = rs(0, ih)
          tdad(0, ih, 2) = tdd(0, ih)
          tda(0, ih, 2) = td(0, ih)
          ttad(0, ih, 2) = ttd(0, ih)
          tta(0, ih, 2) = tt(0, ih)
          rsad(0, ih, 2) = rsd(0, ih)
          rsa(0, ih, 2) = rs(0, ih)
          DO k=1,ict-1
            denmd = (tsd(k, ih)*(1.-rsa(k-1, ih, 1)*rs(k, ih))-ts(k, ih)&
&             *(-(rsad(k-1, ih, 1)*rs(k, ih))-rsa(k-1, ih, 1)*rsd(k, ih)&
&             ))/(1.-rsa(k-1, ih, 1)*rs(k, ih))**2
            denm = ts(k, ih)/(1.-rsa(k-1, ih, 1)*rs(k, ih))
            tdad(k, ih, 1) = tdad(k-1, ih, 1)*td(k, ih) + tda(k-1, ih, 1&
&             )*tdd(k, ih)
            tda(k, ih, 1) = tda(k-1, ih, 1)*td(k, ih)
            ttad(k, ih, 1) = tdad(k-1, ih, 1)*tt(k, ih) + tda(k-1, ih, 1&
&             )*ttd(k, ih) + ((tdad(k-1, ih, 1)*rr(k, ih)+tda(k-1, ih, 1&
&             )*rrd(k, ih))*rsa(k-1, ih, 1)+tda(k-1, ih, 1)*rr(k, ih)*&
&             rsad(k-1, ih, 1)+ttad(k-1, ih, 1)-tdad(k-1, ih, 1))*denm +&
&             (tda(k-1, ih, 1)*rsa(k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)&
&             -tda(k-1, ih, 1))*denmd
            tta(k, ih, 1) = tda(k-1, ih, 1)*tt(k, ih) + (tda(k-1, ih, 1)&
&             *rsa(k-1, ih, 1)*rr(k, ih)+tta(k-1, ih, 1)-tda(k-1, ih, 1)&
&             )*denm
            rsad(k, ih, 1) = rsd(k, ih) + (tsd(k, ih)*denm+ts(k, ih)*&
&             denmd)*rsa(k-1, ih, 1) + ts(k, ih)*denm*rsad(k-1, ih, 1)
            rsa(k, ih, 1) = rs(k, ih) + ts(k, ih)*rsa(k-1, ih, 1)*denm
            tdad(k, ih, 2) = tdad(k, ih, 1)
            tda(k, ih, 2) = tda(k, ih, 1)
            ttad(k, ih, 2) = ttad(k, ih, 1)
            tta(k, ih, 2) = tta(k, ih, 1)
            rsad(k, ih, 2) = rsad(k, ih, 1)
            rsa(k, ih, 2) = rsa(k, ih, 1)
          END DO
! k loop
!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition
          DO k=ict,icb-1
            DO im=1,2
              denmd = (tsd(k, im)*(1.-rsa(k-1, ih, im)*rs(k, im))-ts(k, &
&               im)*(-(rsad(k-1, ih, im)*rs(k, im))-rsa(k-1, ih, im)*rsd&
&               (k, im)))/(1.-rsa(k-1, ih, im)*rs(k, im))**2
              denm = ts(k, im)/(1.-rsa(k-1, ih, im)*rs(k, im))
              tdad(k, ih, im) = tdad(k-1, ih, im)*td(k, im) + tda(k-1, &
&               ih, im)*tdd(k, im)
              tda(k, ih, im) = tda(k-1, ih, im)*td(k, im)
              ttad(k, ih, im) = tdad(k-1, ih, im)*tt(k, im) + tda(k-1, &
&               ih, im)*ttd(k, im) + ((tdad(k-1, ih, im)*rr(k, im)+tda(k&
&               -1, ih, im)*rrd(k, im))*rsa(k-1, ih, im)+tda(k-1, ih, im&
&               )*rr(k, im)*rsad(k-1, ih, im)+ttad(k-1, ih, im)-tdad(k-1&
&               , ih, im))*denm + (tda(k-1, ih, im)*rsa(k-1, ih, im)*rr(&
&               k, im)+tta(k-1, ih, im)-tda(k-1, ih, im))*denmd
              tta(k, ih, im) = tda(k-1, ih, im)*tt(k, im) + (tda(k-1, ih&
&               , im)*rsa(k-1, ih, im)*rr(k, im)+tta(k-1, ih, im)-tda(k-&
&               1, ih, im))*denm
              rsad(k, ih, im) = rsd(k, im) + (tsd(k, im)*denm+ts(k, im)*&
&               denmd)*rsa(k-1, ih, im) + ts(k, im)*denm*rsad(k-1, ih, &
&               im)
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
          rrad(np+1, 1, is) = rrd(np+1, is)
          rra(np+1, 1, is) = rr(np+1, is)
          rxad(np+1, 1, is) = rsd(np+1, is)
          rxa(np+1, 1, is) = rs(np+1, is)
          rrad(np+1, 2, is) = rrd(np+1, is)
          rra(np+1, 2, is) = rr(np+1, is)
          rxad(np+1, 2, is) = rsd(np+1, is)
          rxa(np+1, 2, is) = rs(np+1, is)
          DO k=np,icb,-1
            denmd = (tsd(k, is)*(1.-rs(k, is)*rxa(k+1, 1, is))-ts(k, is)&
&             *(-(rsd(k, is)*rxa(k+1, 1, is))-rs(k, is)*rxad(k+1, 1, is)&
&             ))/(1.-rs(k, is)*rxa(k+1, 1, is))**2
            denm = ts(k, is)/(1.-rs(k, is)*rxa(k+1, 1, is))
            rrad(k, 1, is) = rrd(k, is) + (tdd(k, is)*rra(k+1, 1, is)+td&
&             (k, is)*rrad(k+1, 1, is)+(ttd(k, is)-tdd(k, is))*rxa(k+1, &
&             1, is)+(tt(k, is)-td(k, is))*rxad(k+1, 1, is))*denm + (td(&
&             k, is)*rra(k+1, 1, is)+(tt(k, is)-td(k, is))*rxa(k+1, 1, &
&             is))*denmd
            rra(k, 1, is) = rr(k, is) + (td(k, is)*rra(k+1, 1, is)+(tt(k&
&             , is)-td(k, is))*rxa(k+1, 1, is))*denm
            rxad(k, 1, is) = rsd(k, is) + (tsd(k, is)*denm+ts(k, is)*&
&             denmd)*rxa(k+1, 1, is) + ts(k, is)*denm*rxad(k+1, 1, is)
            rxa(k, 1, is) = rs(k, is) + ts(k, is)*rxa(k+1, 1, is)*denm
            rrad(k, 2, is) = rrad(k, 1, is)
            rra(k, 2, is) = rra(k, 1, is)
            rxad(k, 2, is) = rxad(k, 1, is)
            rxa(k, 2, is) = rxa(k, 1, is)
          END DO
! k loop
!-----for middle clouds
          DO k=icb-1,ict,-1
            DO im=1,2
              denmd = (tsd(k, im)*(1.-rs(k, im)*rxa(k+1, im, is))-ts(k, &
&               im)*(-(rsd(k, im)*rxa(k+1, im, is))-rs(k, im)*rxad(k+1, &
&               im, is)))/(1.-rs(k, im)*rxa(k+1, im, is))**2
              denm = ts(k, im)/(1.-rs(k, im)*rxa(k+1, im, is))
              rrad(k, im, is) = rrd(k, im) + (tdd(k, im)*rra(k+1, im, is&
&               )+td(k, im)*rrad(k+1, im, is)+(ttd(k, im)-tdd(k, im))*&
&               rxa(k+1, im, is)+(tt(k, im)-td(k, im))*rxad(k+1, im, is)&
&               )*denm + (td(k, im)*rra(k+1, im, is)+(tt(k, im)-td(k, im&
&               ))*rxa(k+1, im, is))*denmd
              rra(k, im, is) = rr(k, im) + (td(k, im)*rra(k+1, im, is)+(&
&               tt(k, im)-td(k, im))*rxa(k+1, im, is))*denm
              rxad(k, im, is) = rsd(k, im) + (tsd(k, im)*denm+ts(k, im)*&
&               denmd)*rxa(k+1, im, is) + ts(k, im)*denm*rxad(k+1, im, &
&               is)
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
            chd = -cc1d
            ch = 1.0 - cc1
!-----cloudy portion
          ELSE
            chd = cc1d
            ch = cc1
          END IF
          DO im=1,2
!-----clear portion 
            IF (im .EQ. 1) THEN
              cmd = chd*(1.0-cc2) - ch*cc2d
              cm = ch*(1.0-cc2)
!-----cloudy portion
            ELSE
              cmd = chd*cc2 + ch*cc2d
              cm = ch*cc2
            END IF
            DO is=1,2
!-----clear portion 
              IF (is .EQ. 1) THEN
                ctd = cmd*(1.0-cc3) - cm*cc3d
                ct = cm*(1.0-cc3)
!-----cloudy portion
              ELSE
                ctd = cmd*cc3 + cm*cc3d
                ct = cm*cc3
              END IF
!-----add one layer at a time, going down.
              DO k=icb,np
                denmd = (tsd(k, is)*(1.-rsa(k-1, ih, im)*rs(k, is))-ts(k&
&                 , is)*(-(rsad(k-1, ih, im)*rs(k, is))-rsa(k-1, ih, im)&
&                 *rsd(k, is)))/(1.-rsa(k-1, ih, im)*rs(k, is))**2
                denm = ts(k, is)/(1.-rsa(k-1, ih, im)*rs(k, is))
                tdad(k, ih, im) = tdad(k-1, ih, im)*td(k, is) + tda(k-1&
&                 , ih, im)*tdd(k, is)
                tda(k, ih, im) = tda(k-1, ih, im)*td(k, is)
                ttad(k, ih, im) = tdad(k-1, ih, im)*tt(k, is) + tda(k-1&
&                 , ih, im)*ttd(k, is) + ((tdad(k-1, ih, im)*rr(k, is)+&
&                 tda(k-1, ih, im)*rrd(k, is))*rsa(k-1, ih, im)+tda(k-1&
&                 , ih, im)*rr(k, is)*rsad(k-1, ih, im)+ttad(k-1, ih, im&
&                 )-tdad(k-1, ih, im))*denm + (tda(k-1, ih, im)*rr(k, is&
&                 )*rsa(k-1, ih, im)+tta(k-1, ih, im)-tda(k-1, ih, im))*&
&                 denmd
                tta(k, ih, im) = tda(k-1, ih, im)*tt(k, is) + (tda(k-1, &
&                 ih, im)*rr(k, is)*rsa(k-1, ih, im)+tta(k-1, ih, im)-&
&                 tda(k-1, ih, im))*denm
                rsad(k, ih, im) = rsd(k, is) + (tsd(k, is)*denm+ts(k, is&
&                 )*denmd)*rsa(k-1, ih, im) + ts(k, is)*denm*rsad(k-1, &
&                 ih, im)
                rsa(k, ih, im) = rs(k, is) + ts(k, is)*rsa(k-1, ih, im)*&
&                 denm
              END DO
! k loop
!-----add one layer at a time, going up.
              DO k=ict-1,0,-1
                denmd = (tsd(k, ih)*(1.-rs(k, ih)*rxa(k+1, im, is))-ts(k&
&                 , ih)*(-(rsd(k, ih)*rxa(k+1, im, is))-rs(k, ih)*rxad(k&
&                 +1, im, is)))/(1.-rs(k, ih)*rxa(k+1, im, is))**2
                denm = ts(k, ih)/(1.-rs(k, ih)*rxa(k+1, im, is))
                rrad(k, im, is) = rrd(k, ih) + (tdd(k, ih)*rra(k+1, im, &
&                 is)+td(k, ih)*rrad(k+1, im, is)+(ttd(k, ih)-tdd(k, ih)&
&                 )*rxa(k+1, im, is)+(tt(k, ih)-td(k, ih))*rxad(k+1, im&
&                 , is))*denm + (td(k, ih)*rra(k+1, im, is)+(tt(k, ih)-&
&                 td(k, ih))*rxa(k+1, im, is))*denmd
                rra(k, im, is) = rr(k, ih) + (td(k, ih)*rra(k+1, im, is)&
&                 +(tt(k, ih)-td(k, ih))*rxa(k+1, im, is))*denm
                rxad(k, im, is) = rsd(k, ih) + (tsd(k, ih)*denm+ts(k, ih&
&                 )*denmd)*rxa(k+1, im, is) + ts(k, ih)*denm*rxad(k+1, &
&                 im, is)
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
                denmd = -((-(rsad(k-1, ih, im)*rxa(k, im, is))-rsa(k-1, &
&                 ih, im)*rxad(k, im, is))/(1.-rsa(k-1, ih, im)*rxa(k, &
&                 im, is))**2)
                denm = 1./(1.-rsa(k-1, ih, im)*rxa(k, im, is))
                fdndird = tdad(k-1, ih, im)
                fdndir = tda(k-1, ih, im)
                xx4d = tdad(k-1, ih, im)*rra(k, im, is) + tda(k-1, ih, &
&                 im)*rrad(k, im, is)
                xx4 = tda(k-1, ih, im)*rra(k, im, is)
                yyd = ttad(k-1, ih, im) - tdad(k-1, ih, im)
                yy = tta(k-1, ih, im) - tda(k-1, ih, im)
                fdndifd = (xx4d*rsa(k-1, ih, im)+xx4*rsad(k-1, ih, im)+&
&                 yyd)*denm + (xx4*rsa(k-1, ih, im)+yy)*denmd
                fdndif = (xx4*rsa(k-1, ih, im)+yy)*denm
                fupdifd = (xx4d+yyd*rxa(k, im, is)+yy*rxad(k, im, is))*&
&                 denm + (xx4+yy*rxa(k, im, is))*denmd
                fupdif = (xx4+yy*rxa(k, im, is))*denm
                flxdnd = fdndird + fdndifd - fupdifd
                flxdn = fdndir + fdndif - fupdif
!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)
                IF (ih .EQ. 1 .AND. im .EQ. 1 .AND. is .EQ. 1) THEN
                  fupc(k) = fupdif
                  fclr(k) = flxdn
                END IF
                fupa(k) = fupa(k) + fupdif*ct
                falld(k) = falld(k) + flxdnd*ct + flxdn*ctd
                fall(k) = fall(k) + flxdn*ct
              END DO
! k loop
              fsdir = fsdir + fdndir*ct
              fsdif = fsdif + fdndif*ct
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
          flx_devd(i, k) = flx_devd(i, k) + hk_ir(ib, ik)*falld(k)
          flx_dev(i, k) = flx_dev(i, k) + fall(k)*hk_ir(ib, ik)
          flc_dev(i, k) = flc_dev(i, k) + fclr(k)*hk_ir(ib, ik)
          flxu_dev(i, k) = flxu_dev(i, k) + fupa(k)*hk_ir(ib, ik)
          flcu_dev(i, k) = flcu_dev(i, k) + fupc(k)*hk_ir(ib, ik)
        END DO
!-----compute downward surface fluxes in the ir region
        fdirir_dev(i) = fdirir_dev(i) + fsdir*hk_ir(ib, ik)
        fdifir_dev(i) = fdifir_dev(i) + fsdif*hk_ir(ib, ik)
!-----tabulate surface flux at ir bands
        flx_sfc_band_dev(i, iv) = flx_sfc_band_dev(i, iv) + fall(np+1)*&
&         hk_ir(ib, ik)
      END DO
    END DO
! ik loop
!-----compute pressure-scaled o2 amount following Eq. (3.5) with f=1.
!     unit is (cm-atm)stp. 165.22 = (1000/980)*23.14%*(22400/32)
!     compute flux reduction due to oxygen following Eq. (3.18). 0.0633 is the
!     fraction of insolation contained in the oxygen bands
    dfd(0) = 0.0
    df(0) = 0.0
    cnt = 165.22*snt
    so2d(1) = 0.0
    so2(1) = scal0*cnt
! LLT increased parameter 145 to 155 to enhance effect
    result1 = SQRT(so2(1))
    dfd(1) = 0.0
    df(1) = 0.0633*(1.-EXP(-(0.000155*result1)))
    DO k=1,np
      so2d(k+1) = 0.0
      so2(k+1) = so2(k) + scal(k)*cnt
! LLT increased parameter 145 to 155 to enhance effect
      result1 = SQRT(so2(k+1))
      dfd(k+1) = 0.0
      df(k+1) = 0.0633*(1.0-EXP(-(0.000155*result1)))
    END DO
!-----for solar heating due to co2 scaling follows Eq(3.5) with f=1.
!     unit is (cm-atm)stp. 789 = (1000/980)*(44/28.97)*(22400/44)
    so2d(1) = 0.0
    so2(1) = 789.*co2*scal0
    DO k=1,np
      so2d(k+1) = 0.0
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
    du2 = du*du
    dw2 = dw*dw
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
      x4d = swhd(k)/(swh(k)*LOG(10.0))
      x4 = LOG10(swh(k)*snt)
      IF (x4 .GT. y0) THEN
        wlog = y0
        wlogd = 0.0
      ELSE
        wlogd = x4d
        wlog = x4
      END IF
      ic = INT((ulog-x1)/du + 1.)
      iw = INT((wlog-y1)/dw + 1.)
      IF (ic .LT. 2) ic = 2
      IF (iw .LT. 2) iw = 2
      IF (ic .GT. nu) ic = nu
      IF (iw .GT. nw) iw = nw
      dc = ulog - REAL(ic-2)*du - u1
      ddd = wlogd
      dd = wlog - REAL(iw-2)*dw - w1
      x2d = (cah(ic-1, iw)-cah(ic-1, iw-1))*ddd/dw
      x2 = cah(ic-1, iw-1) + (cah(ic-1, iw)-cah(ic-1, iw-1))/dw*dd
      y2d = x2d
      y2 = x2 + (cah(ic, iw-1)-cah(ic-1, iw-1))/du*dc
      IF (y2 .LT. 0.0) THEN
        y2 = 0.0
        y2d = 0.0
      ELSE
        y2 = y2
      END IF
! LLT increase CO2 effect to help reduce cold tropopause bias
      dfd(k) = dfd(k) + 1.5*y2d
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
    du2 = du*du
    dw2 = dw*dw
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
      ic = INT((ulog-x1)/du + 1.)
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
    dftopd = dfd(ntop)
    dftop = df(ntop)
    DO k=1,np+1
      IF (k .GT. ntop) THEN
        xx4d = (flx_devd(i, k)*flx_dev(i, ntop)-flx_dev(i, k)*flx_devd(i&
&         , ntop))/flx_dev(i, ntop)**2
        xx4 = flx_dev(i, k)/flx_dev(i, ntop)
        dfd(k) = dftopd + xx4d*(df(k)-dftop) + xx4*(dfd(k)-dftopd)
        df(k) = dftop + xx4*(df(k)-dftop)
      END IF
    END DO
!-----update the net fluxes
    DO k=1,np+1
      IF (df(k) .GT. flx_dev(i, k) - 1.0e-8) THEN
        dfd(k) = flx_devd(i, k)
        df(k) = flx_dev(i, k) - 1.0e-8
      ELSE
        df(k) = df(k)
      END IF
!           df(k) = 0.0
      flx_devd(i, k) = flx_devd(i, k) - dfd(k)
      flx_dev(i, k) = flx_dev(i, k) - df(k)
      flc_dev(i, k) = flc_dev(i, k) - df(k)
    END DO
!-----update the downward surface fluxes 
!        xx4 = fdirir (i) + fdifir (i) +&
!              fdiruv (i) + fdifuv (i) +&
!              fdirpar(i) + fdifpar(i)
    xx4 = flx_dev(i, np+1) + df(np+1)
    IF (xx4 .GE. 0.) THEN
      abs0 = xx4
    ELSE
      abs0 = -xx4
    END IF
    result1 = EPSILON(1.0)
    IF (abs0 .GT. result1) THEN
      IF (1.0 - df(np+1)/xx4 .GT. 1.) THEN
        x6 = 1.
      ELSE
        x6 = 1.0 - df(np+1)/xx4
      END IF
      IF (x6 .LT. 0.) THEN
        xx4 = 0.
      ELSE
        xx4 = x6
      END IF
    ELSE
      xx4 = 0.0
    END IF
    fdirir_dev(i) = xx4*fdirir_dev(i)
    fdifir_dev(i) = xx4*fdifir_dev(i)
    fdiruv_dev(i) = xx4*fdiruv_dev(i)
    fdifuv_dev(i) = xx4*fdifuv_dev(i)
    fdirpar_dev(i) = xx4*fdirpar_dev(i)
    fdifpar_dev(i) = xx4*fdifpar_dev(i)
    DO ib=1,nband
      flx_sfc_band_dev(i, ib) = xx4*flx_sfc_band_dev(i, ib)
    END DO
  END DO run_loop
END SUBROUTINE SORAD_D

!  Differentiation of deledd in forward (tangent) mode:
!   variations   of useful results: tt1 td1 rr1
!   with respect to varying inputs: g01 tau1 ssc1
!*********************************************************************
SUBROUTINE DELEDD_D(tau1, tau1d, ssc1, ssc1d, g01, g01d, cza1, rr1, rr1d&
& , tt1, tt1d, td1, td1d)
  IMPLICIT NONE
! 8 byte real
  INTEGER, PARAMETER :: real_de=8
!integer,parameter :: REAL_SP = 4 ! 4 byte real
!-----input parameters
  REAL*4, INTENT(IN) :: tau1, ssc1, g01, cza1
  REAL*4, INTENT(IN) :: tau1d, ssc1d, g01d
!-----output parameters
  REAL*4, INTENT(OUT) :: rr1, tt1, td1
  REAL*4, INTENT(OUT) :: rr1d, tt1d, td1d
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
  REAL*8 :: taud, sscd, g0d, rrd, ttd, tdd
  REAL*8 :: zth, ff, xx, taup, sscp, gp, gm1, gm2, gm3, akk, alf1, alf2
  REAL*8 :: ffd, xxd, taupd, sscpd, gpd, gm1d, gm2d, gm3d, akkd, alf1d, &
& alf2d
  REAL*8 :: all, bll, st7, st8, cll, dll, fll, ell, st1, st2, st3, st4
  REAL*8 :: alld, blld, st7d, st8d, clld, dlld, flld, elld, st1d, st2d, &
& st3d, st4d
  INTRINSIC DBLE
  INTRINSIC SQRT
  INTRINSIC ABS
  INTRINSIC EXP
  INTRINSIC MAX
  INTRINSIC REAL
  REAL*8 :: arg1
  REAL*8 :: arg1d
  REAL*8 :: abs0
!zth = real(cza1,kind=REAL_DE)
!g0  = real(g01 ,kind=REAL_DE)
!tau = real(tau1,kind=REAL_DE)
!ssc = real(ssc1,kind=REAL_DE)
  zth = DBLE(cza1)
  g0d = g01d
  g0 = DBLE(g01)
  taud = tau1d
  tau = DBLE(tau1)
  sscd = ssc1d
  ssc = DBLE(ssc1)
  ffd = g0d*g0 + g0*g0d
  ff = g0*g0
  xxd = -(ffd*ssc) - ff*sscd
  xx = one - ff*ssc
  taupd = taud*xx + tau*xxd
  taup = tau*xx
  sscpd = ((sscd*(one-ff)-ssc*ffd)*xx-ssc*(one-ff)*xxd)/xx**2
  sscp = ssc*(one-ff)/xx
  gpd = (g0d*(one+g0)-g0*g0d)/(one+g0)**2
  gp = g0/(one+g0)
  xxd = three*gpd
  xx = three*gp
  gm1d = fourth*(-(sscpd*(four+xx))-sscp*xxd)
  gm1 = (seven-sscp*(four+xx))*fourth
  gm2d = -(fourth*(sscp*xxd-sscpd*(four-xx)))
  gm2 = -((one-sscp*(four-xx))*fourth)
  arg1d = (gm1d+gm2d)*(gm1-gm2) + (gm1+gm2)*(gm1d-gm2d)
  arg1 = (gm1+gm2)*(gm1-gm2)
  IF (arg1 .EQ. 0.0) THEN
    akkd = 0.0_8
  ELSE
    akkd = arg1d/(2.0*SQRT(arg1))
  END IF
  akk = SQRT(arg1)
  xxd = zth*akkd
  xx = akk*zth
  st7d = -xxd
  st7 = one - xx
  st8d = xxd
  st8 = one + xx
  st3d = st7d*st8 + st7*st8d
  st3 = st7*st8
  IF (st3 .GE. 0.) THEN
    abs0 = st3
  ELSE
    abs0 = -st3
  END IF
  IF (abs0 .LT. thresh) THEN
    zth = zth + 0.0010
    IF (zth .GT. 1.0) zth = zth - 0.0020
    xxd = zth*akkd
    xx = akk*zth
    st7d = -xxd
    st7 = one - xx
    st8d = xxd
    st8 = one + xx
    st3d = st7d*st8 + st7*st8d
    st3 = st7*st8
  END IF
  tdd = -(taupd*EXP(-(taup/zth))/zth)
  td = EXP(-(taup/zth))
  gm3d = -(fourth*zth*three*gpd)
  gm3 = (two-zth*three*gp)*fourth
  xxd = gm1d - gm2d
  xx = gm1 - gm2
  alf1d = gm1d - gm3d*xx - gm3*xxd
  alf1 = gm1 - gm3*xx
  alf2d = gm2d + gm3d*xx + gm3*xxd
  alf2 = gm2 + gm3*xx
  xxd = two*akkd
  xx = akk*two
  alld = (gm3d-zth*alf2d)*xx*td + (gm3-alf2*zth)*(xxd*td+xx*tdd)
  all = (gm3-alf2*zth)*xx*td
  blld = (zth*alf1d-gm3d)*xx + (one-gm3+alf1*zth)*xxd
  bll = (one-gm3+alf1*zth)*xx
  xxd = akkd*gm3 + akk*gm3d
  xx = akk*gm3
  clld = (alf2d+xxd)*st7 + (alf2+xx)*st7d
  cll = (alf2+xx)*st7
  dlld = (alf2d-xxd)*st8 + (alf2-xx)*st8d
  dll = (alf2-xx)*st8
  xxd = akkd*(one-gm3) - akk*gm3d
  xx = akk*(one-gm3)
  flld = (alf1d+xxd)*st8 + (alf1+xx)*st8d
  fll = (alf1+xx)*st8
  elld = (alf1d-xxd)*st7 + (alf1-xx)*st7d
  ell = (alf1-xx)*st7
  st2d = -((akkd*taup+akk*taupd)*EXP(-(akk*taup)))
  st2 = EXP(-(akk*taup))
  st4d = st2d*st2 + st2*st2d
  st4 = st2*st2
  st1d = (sscpd*(akk+gm1+(akk-gm1)*st4)*st3-sscp*((akkd+gm1d+(akkd-gm1d)&
&   *st4+(akk-gm1)*st4d)*st3+(akk+gm1+(akk-gm1)*st4)*st3d))/((akk+gm1+(&
&   akk-gm1)*st4)*st3)**2
  st1 = sscp/((akk+gm1+(akk-gm1)*st4)*st3)
  rrd = (clld-dlld*st4-dll*st4d-alld*st2-all*st2d)*st1 + (cll-dll*st4-&
&   all*st2)*st1d
  rr = (cll-dll*st4-all*st2)*st1
  ttd = -(((flld-elld*st4-ell*st4d)*td+(fll-ell*st4)*tdd-blld*st2-bll*&
&   st2d)*st1+((fll-ell*st4)*td-bll*st2)*st1d)
  tt = -(((fll-ell*st4)*td-bll*st2)*st1)
  IF (rr .LT. zero) THEN
    rr = zero
    rrd = 0.0_8
  ELSE
    rr = rr
  END IF
  IF (tt .LT. zero) THEN
    tt = zero
    ttd = 0.0_8
  ELSE
    tt = tt
  END IF
  ttd = ttd + tdd
  tt = tt + td
!td1 = real(td,kind=REAL_SP)
!rr1 = real(rr,kind=REAL_SP)
!tt1 = real(tt,kind=REAL_SP)
  td1d = tdd
  td1 = REAL(td)
  rr1d = rrd
  rr1 = REAL(rr)
  tt1d = ttd
  tt1 = REAL(tt)
END SUBROUTINE DELEDD_D

!  Differentiation of getvistau1 in forward (tangent) mode:
!   variations   of useful results: asycl taudiff taubeam
!   with respect to varying inputs: hydromets asycl fcld reff
SUBROUTINE GETVISTAU1_D(nlevs, cosz, dp, fcld, fcldd, reff, reffd, &
& hydromets, hydrometsd, ict, icb, taubeam, taubeamd, taudiff, taudiffd&
& , asycl, asycld, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv, arb_uv, &
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
  REAL, INTENT(IN) :: fcldd(nlevs)
!  Effective radius (microns)
  REAL, INTENT(IN) :: reff(nlevs, 4)
  REAL, INTENT(IN) :: reffd(nlevs, 4)
!  Hydrometeors (kg/kg)
  REAL, INTENT(IN) :: hydromets(nlevs, 4)
  REAL, INTENT(IN) :: hydrometsd(nlevs, 4)
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
  REAL, INTENT(OUT) :: taubeam(nlevs, 4)
  REAL, INTENT(OUT) :: taubeamd(nlevs, 4)
!  Optical Depth for Diffuse Radiation
  REAL, INTENT(OUT) :: taudiff(nlevs, 4)
  REAL, INTENT(OUT) :: taudiffd(nlevs, 4)
!  Cloud Asymmetry Factor
  REAL, INTENT(OUT) :: asycl(nlevs)
  REAL, INTENT(OUT) :: asycld(nlevs)
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
  REAL :: ftd, fad, xaid, taucd, asycltd
  REAL :: cc(3)
  REAL :: ccd(3)
  REAL :: taucld1, taucld2, taucld3, taucld4
  REAL :: taucld1d, taucld2d, taucld3d, taucld4d
  REAL :: g1, g2, g3, g4
  REAL :: g1d, g2d, g3d, g4d
  REAL :: reff_snow
  REAL :: reff_snowd
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
  taubeam = 0.0
  taudiff = 0.0
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
    ccd = 0.0
    DO k=1,ict-1
      IF (cc(1) .LT. fcld(k)) THEN
        ccd(1) = fcldd(k)
        cc(1) = fcld(k)
      ELSE
        cc(1) = cc(1)
      END IF
    END DO
    DO k=ict,icb-1
      IF (cc(2) .LT. fcld(k)) THEN
        ccd(2) = fcldd(k)
        cc(2) = fcld(k)
      ELSE
        cc(2) = cc(2)
      END IF
    END DO
    DO k=icb,nlevs
      IF (cc(3) .LT. fcld(k)) THEN
        ccd(3) = fcldd(k)
        cc(3) = fcld(k)
      ELSE
        cc(3) = cc(3)
      END IF
    END DO
    taudiffd = 0.0
    taubeamd = 0.0
  ELSE
    taudiffd = 0.0
    taubeamd = 0.0
    ccd = 0.0
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
      taucld1 = 0.
      taucld1d = 0.0
    ELSE
      taucld1d = (dp(k)*1.0e3*aib_uv*hydrometsd(k, 1)*reff(k, 1)/&
&       cons_grav-dp(k)*1.0e3*hydromets(k, 1)*aib_uv*reffd(k, 1)/&
&       cons_grav)/reff(k, 1)**2
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*aib_uv/reff(k, 1)
    END IF
    IF (reff(k, 2) .LE. 0.) THEN
      taucld2 = 0.
      taucld2d = 0.0
    ELSE
      taucld2d = dp(k)*1.0e3*(hydrometsd(k, 2)*(awb_uv(1)+awb_uv(2)/reff&
&       (k, 2))-hydromets(k, 2)*awb_uv(2)*reffd(k, 2)/reff(k, 2)**2)/&
&       cons_grav
      taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_uv(1)+awb_uv(&
&       2)/reff(k, 2))
    END IF
    taucld3d = dp(k)*1.0e3*arb_uv(1)*hydrometsd(k, 3)/cons_grav
    taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_uv(1)
    IF (reff(k, 4) .GT. 112.0) THEN
      reff_snow = 112.0
      reff_snowd = 0.0
    ELSE
      reff_snowd = reffd(k, 4)
      reff_snow = reff(k, 4)
    END IF
    IF (reff_snow .LE. 0.) THEN
      taucld4 = 0.
      taucld4d = 0.0
    ELSE
      taucld4d = (dp(k)*1.0e3*aib_uv*hydrometsd(k, 4)*reff_snow/&
&       cons_grav-dp(k)*1.0e3*hydromets(k, 4)*aib_uv*reff_snowd/&
&       cons_grav)/reff_snow**2
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*aib_uv/reff_snow
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
        kk = 1
      ELSE IF (k .GE. ict .AND. k .LT. icb) THEN
        kk = 2
      ELSE
        kk = 3
      END IF
      taucd = taucld1d + taucld2d + taucld3d + taucld4d
      tauc = taucld1 + taucld2 + taucld3 + taucld4
      IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
!-----normalize cloud cover following Eq. (7.8)
        fad = (fcldd(k)*cc(kk)-fcld(k)*ccd(kk))/cc(kk)**2
        fa = fcld(k)/cc(kk)
        IF (tauc .GT. 32.) THEN
          tauc = 32.
          taucd = 0.0
        ELSE
          tauc = tauc
        END IF
        fm = cosz/dm
        ftd = taucd/(tauc*LOG(10.0))/dt
        ft = (LOG10(tauc)-t1)/dt
        fad = fad/da
        fa = fa/da
        im = INT(fm + 1.5)
        it = INT(ft + 1.5)
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
        xai = (-(caib(im-1, it, ia)*(1.-fm))+caib(im+1, it, ia)*(1.+fm))&
&         *fm*.5 + caib(im, it, ia)*(1.-fm*fm)
        xaid = .5*((caib(im, it-1, ia)*ftd+caib(im, it+1, ia)*ftd)*ft+(-&
&         (caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(1.+ft))*ftd) &
&         + caib(im, it, ia)*(-(ftd*ft)-ft*ftd)
        xai = xai + (-(caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(&
&         1.+ft))*ft*.5 + caib(im, it, ia)*(1.-ft*ft)
        xaid = xaid + .5*((caib(im, it, ia-1)*fad+caib(im, it, ia+1)*fad&
&         )*fa+(-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(1.+fa)&
&         )*fad) + caib(im, it, ia)*(-(fad*fa)-fa*fad)
        xai = xai + (-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(&
&         1.+fa))*fa*.5 + caib(im, it, ia)*(1.-fa*fa)
        xai = xai - 2.*caib(im, it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          xaid = 0.0
        ELSE
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          xaid = 0.0
        ELSE
          xai = xai
        END IF
        taubeamd(k, 1) = taucld1d*xai + taucld1*xaid
        taubeam(k, 1) = taucld1*xai
        taubeamd(k, 2) = taucld2d*xai + taucld2*xaid
        taubeam(k, 2) = taucld2*xai
        taubeamd(k, 3) = taucld3d*xai + taucld3*xaid
        taubeam(k, 3) = taucld3*xai
        taubeamd(k, 4) = taucld4d*xai + taucld4*xaid
        taubeam(k, 4) = taucld4*xai
!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
        xaid = .5*((caif(it-1, ia)*ftd+caif(it+1, ia)*ftd)*ft+(-(caif(it&
&         -1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ftd) + caif(it, ia)*(&
&         -(ftd*ft)-ft*ftd)
        xai = (-(caif(it-1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ft*.5 +&
&         caif(it, ia)*(1.-ft*ft)
        xaid = xaid + .5*((caif(it, ia-1)*fad+caif(it, ia+1)*fad)*fa+(-(&
&         caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*fad) + caif(it&
&         , ia)*(-(fad*fa)-fa*fad)
        xai = xai + (-(caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*&
&         fa*.5 + caif(it, ia)*(1.-fa*fa)
        xai = xai - caif(it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          xaid = 0.0
        ELSE
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          xaid = 0.0
        ELSE
          xai = xai
        END IF
        taudiffd(k, 1) = taucld1d*xai + taucld1*xaid
        taudiff(k, 1) = taucld1*xai
        taudiffd(k, 2) = taucld2d*xai + taucld2*xaid
        taudiff(k, 2) = taucld2*xai
        taudiffd(k, 3) = taucld3d*xai + taucld3*xaid
        taudiff(k, 3) = taucld3*xai
        taudiffd(k, 4) = taucld4d*xai + taucld4*xaid
        taudiff(k, 4) = taucld4*xai
      END IF
    ELSE
! Overlap calculation scaling not needed
      taubeamd(k, 1) = taucld1d
      taubeam(k, 1) = taucld1
      taubeamd(k, 2) = taucld2d
      taubeam(k, 2) = taucld2
      taubeamd(k, 3) = taucld3d
      taubeam(k, 3) = taucld3
      taubeamd(k, 4) = taucld4d
      taubeam(k, 4) = taucld4
      taudiffd(k, 1) = taucld1d
      taudiff(k, 1) = taucld1
      taudiffd(k, 2) = taucld2d
      taudiff(k, 2) = taucld2
      taudiffd(k, 3) = taucld3d
      taudiff(k, 3) = taucld3
      taudiffd(k, 4) = taucld4d
      taudiff(k, 4) = taucld4
    END IF
!-----cloud asymmetry factor for a mixture of liquid and ice particles.
!     unit of reff is micrometers. Eqs. (4.8) and (6.4)
    asyclt = 1.0
    taucd = taucld1d + taucld2d + taucld3d + taucld4d
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      g1d = (aig_uv(3)*reffd(k, 1)*reff(k, 1)+(aig_uv(2)+aig_uv(3)*reff(&
&       k, 1))*reffd(k, 1))*taucld1 + (aig_uv(1)+(aig_uv(2)+aig_uv(3)*&
&       reff(k, 1))*reff(k, 1))*taucld1d
      g1 = (aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff(k, 1))*reff(k, 1))*&
&       taucld1
      g2d = (awg_uv(3)*reffd(k, 2)*reff(k, 2)+(awg_uv(2)+awg_uv(3)*reff(&
&       k, 2))*reffd(k, 2))*taucld2 + (awg_uv(1)+(awg_uv(2)+awg_uv(3)*&
&       reff(k, 2))*reff(k, 2))*taucld2d
      g2 = (awg_uv(1)+(awg_uv(2)+awg_uv(3)*reff(k, 2))*reff(k, 2))*&
&       taucld2
      g3d = arg_uv(1)*taucld3d
      g3 = arg_uv(1)*taucld3
      g4d = (aig_uv(3)*reff_snowd*reff_snow+(aig_uv(2)+aig_uv(3)*&
&       reff_snow)*reff_snowd)*taucld4 + (aig_uv(1)+(aig_uv(2)+aig_uv(3)&
&       *reff_snow)*reff_snow)*taucld4d
      g4 = (aig_uv(1)+(aig_uv(2)+aig_uv(3)*reff_snow)*reff_snow)*taucld4
      asycltd = ((g1d+g2d+g3d+g4d)*tauc-(g1+g2+g3+g4)*taucd)/tauc**2
      asyclt = (g1+g2+g3+g4)/tauc
    ELSE
      asycltd = 0.0
    END IF
    asycld(k) = asycltd
    asycl(k) = asyclt
  END DO
  RETURN
END SUBROUTINE GETVISTAU1_D

!  Differentiation of getnirtau1 in forward (tangent) mode:
!   variations   of useful results: asycl taudiff ssacl taubeam
!   with respect to varying inputs: hydromets asycl fcld ssacl
!                reff
SUBROUTINE GETNIRTAU1_D(ib, nlevs, cosz, dp, fcld, fcldd, reff, reffd, &
& hydromets, hydrometsd, ict, icb, taubeam, taubeamd, taudiff, taudiffd&
& , asycl, asycld, ssacl, ssacld, aig_uv, awg_uv, arg_uv, aib_uv, awb_uv&
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
  REAL, INTENT(IN) :: fcldd(nlevs)
!  Effective radius (microns)
  REAL, INTENT(IN) :: reff(nlevs, 4)
  REAL, INTENT(IN) :: reffd(nlevs, 4)
!  Hydrometeors (kg/kg)
  REAL, INTENT(IN) :: hydromets(nlevs, 4)
  REAL, INTENT(IN) :: hydrometsd(nlevs, 4)
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
  REAL, INTENT(OUT) :: taubeam(nlevs, 4)
  REAL, INTENT(OUT) :: taubeamd(nlevs, 4)
!  Optical depth for diffuse radiation
  REAL, INTENT(OUT) :: taudiff(nlevs, 4)
  REAL, INTENT(OUT) :: taudiffd(nlevs, 4)
!  Cloud single scattering albedo
  REAL, INTENT(OUT) :: ssacl(nlevs)
  REAL, INTENT(OUT) :: ssacld(nlevs)
!  Cloud asymmetry factor
  REAL, INTENT(OUT) :: asycl(nlevs)
  REAL, INTENT(OUT) :: asycld(nlevs)
  INTEGER :: k, in, im, it, ia, kk
  REAL :: fm, ft, fa, xai, tauc, asyclt, ssaclt
  REAL :: ftd, fad, xaid, taucd, asycltd, ssacltd
  REAL :: cc(3)
  REAL :: ccd(3)
  REAL :: taucld1, taucld2, taucld3, taucld4
  REAL :: taucld1d, taucld2d, taucld3d, taucld4d
  REAL :: g1, g2, g3, g4
  REAL :: g1d, g2d, g3d, g4d
  REAL :: w1, w2, w3, w4
  REAL :: w1d, w2d, w3d, w4d
  REAL :: reff_snow
  REAL :: reff_snowd
  INTEGER, PARAMETER :: nm=11, nt=9, na=11
  REAL, PARAMETER :: dm=0.1, dt=0.30103, da=0.1, t1=-0.9031
  INTRINSIC MAX
  INTRINSIC MIN
  INTRINSIC LOG10
  INTRINSIC INT
  INTRINSIC REAL
  taubeam = 0.0
  taudiff = 0.0
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
    ccd = 0.0
    DO k=1,ict-1
      IF (cc(1) .LT. fcld(k)) THEN
        ccd(1) = fcldd(k)
        cc(1) = fcld(k)
      ELSE
        cc(1) = cc(1)
      END IF
    END DO
    DO k=ict,icb-1
      IF (cc(2) .LT. fcld(k)) THEN
        ccd(2) = fcldd(k)
        cc(2) = fcld(k)
      ELSE
        cc(2) = cc(2)
      END IF
    END DO
    DO k=icb,nlevs
      IF (cc(3) .LT. fcld(k)) THEN
        ccd(3) = fcldd(k)
        cc(3) = fcld(k)
      ELSE
        cc(3) = cc(3)
      END IF
    END DO
    taudiffd = 0.0
    taubeamd = 0.0
  ELSE
    taudiffd = 0.0
    taubeamd = 0.0
    ccd = 0.0
  END IF
!-----Compute cloud optical thickness.  Eqs. (4.6) and (4.10)
!     taucld1 is the optical thickness for ice particles
!     taucld2 is the optical thickness for liquid particles
!     taucld3 is the optical thickness for rain drops
!     taucld4 is the optical thickness for snow
  DO k=1,nlevs
    IF (reff(k, 1) .LE. 0.) THEN
      taucld1 = 0.
      taucld1d = 0.0
    ELSE
      taucld1d = (dp(k)*1.0e3*aib_nir*hydrometsd(k, 1)*reff(k, 1)/&
&       cons_grav-dp(k)*1.0e3*hydromets(k, 1)*aib_nir*reffd(k, 1)/&
&       cons_grav)/reff(k, 1)**2
      taucld1 = dp(k)*1.0e3/cons_grav*hydromets(k, 1)*aib_nir/reff(k, 1)
    END IF
    IF (reff(k, 2) .LE. 0.) THEN
      taucld2 = 0.
      taucld2d = 0.0
    ELSE
      taucld2d = dp(k)*1.0e3*(hydrometsd(k, 2)*(awb_nir(ib, 1)+awb_nir(&
&       ib, 2)/reff(k, 2))-hydromets(k, 2)*awb_nir(ib, 2)*reffd(k, 2)/&
&       reff(k, 2)**2)/cons_grav
      taucld2 = dp(k)*1.0e3/cons_grav*hydromets(k, 2)*(awb_nir(ib, 1)+&
&       awb_nir(ib, 2)/reff(k, 2))
    END IF
    taucld3d = dp(k)*1.0e3*arb_nir(ib, 1)*hydrometsd(k, 3)/cons_grav
    taucld3 = dp(k)*1.0e3/cons_grav*hydromets(k, 3)*arb_nir(ib, 1)
    IF (reff(k, 4) .GT. 112.0) THEN
      reff_snow = 112.0
      reff_snowd = 0.0
    ELSE
      reff_snowd = reffd(k, 4)
      reff_snow = reff(k, 4)
    END IF
    IF (reff_snow .LE. 0.) THEN
      taucld4 = 0.
      taucld4d = 0.0
    ELSE
      taucld4d = (dp(k)*1.0e3*aib_nir*hydrometsd(k, 4)*reff_snow/&
&       cons_grav-dp(k)*1.0e3*hydromets(k, 4)*aib_nir*reff_snowd/&
&       cons_grav)/reff_snow**2
      taucld4 = dp(k)*1.0e3/cons_grav*hydromets(k, 4)*aib_nir/reff_snow
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
        kk = 1
      ELSE IF (k .GE. ict .AND. k .LT. icb) THEN
        kk = 2
      ELSE
        kk = 3
      END IF
      taucd = taucld1d + taucld2d + taucld3d + taucld4d
      tauc = taucld1 + taucld2 + taucld3 + taucld4
      IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
!-----normalize cloud cover following Eq. (7.8)
        IF (cc(kk) .NE. 0.0) THEN
          fad = (fcldd(k)*cc(kk)-fcld(k)*ccd(kk))/cc(kk)**2
          fa = fcld(k)/cc(kk)
        ELSE
          fa = 0.0
          fad = 0.0
        END IF
        IF (tauc .GT. 32.) THEN
          tauc = 32.
          taucd = 0.0
        ELSE
          tauc = tauc
        END IF
        fm = cosz/dm
        ftd = taucd/(tauc*LOG(10.0))/dt
        ft = (LOG10(tauc)-t1)/dt
        fad = fad/da
        fa = fa/da
        im = INT(fm + 1.5)
        it = INT(ft + 1.5)
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
        xai = (-(caib(im-1, it, ia)*(1.-fm))+caib(im+1, it, ia)*(1.+fm))&
&         *fm*.5 + caib(im, it, ia)*(1.-fm*fm)
        xaid = .5*((caib(im, it-1, ia)*ftd+caib(im, it+1, ia)*ftd)*ft+(-&
&         (caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(1.+ft))*ftd) &
&         + caib(im, it, ia)*(-(ftd*ft)-ft*ftd)
        xai = xai + (-(caib(im, it-1, ia)*(1.-ft))+caib(im, it+1, ia)*(&
&         1.+ft))*ft*.5 + caib(im, it, ia)*(1.-ft*ft)
        xaid = xaid + .5*((caib(im, it, ia-1)*fad+caib(im, it, ia+1)*fad&
&         )*fa+(-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(1.+fa)&
&         )*fad) + caib(im, it, ia)*(-(fad*fa)-fa*fad)
        xai = xai + (-(caib(im, it, ia-1)*(1.-fa))+caib(im, it, ia+1)*(&
&         1.+fa))*fa*.5 + caib(im, it, ia)*(1.-fa*fa)
        xai = xai - 2.*caib(im, it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          xaid = 0.0
        ELSE
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          xaid = 0.0
        ELSE
          xai = xai
        END IF
        taubeamd(k, 1) = taucld1d*xai + taucld1*xaid
        taubeam(k, 1) = taucld1*xai
        taubeamd(k, 2) = taucld2d*xai + taucld2*xaid
        taubeam(k, 2) = taucld2*xai
        taubeamd(k, 3) = taucld3d*xai + taucld3*xaid
        taubeam(k, 3) = taucld3*xai
        taubeamd(k, 4) = taucld4d*xai + taucld4*xaid
        taubeam(k, 4) = taucld4*xai
!-----scale cloud optical thickness for diffuse radiation following 
!     Eq. (7.4).
!     the scaling factor, xai, is a function of the cloud optical
!     thickness and cover but not the solar zenith angle.
        xaid = .5*((caif(it-1, ia)*ftd+caif(it+1, ia)*ftd)*ft+(-(caif(it&
&         -1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ftd) + caif(it, ia)*(&
&         -(ftd*ft)-ft*ftd)
        xai = (-(caif(it-1, ia)*(1.-ft))+caif(it+1, ia)*(1.+ft))*ft*.5 +&
&         caif(it, ia)*(1.-ft*ft)
        xaid = xaid + .5*((caif(it, ia-1)*fad+caif(it, ia+1)*fad)*fa+(-(&
&         caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*fad) + caif(it&
&         , ia)*(-(fad*fa)-fa*fad)
        xai = xai + (-(caif(it, ia-1)*(1.-fa))+caif(it, ia+1)*(1.+fa))*&
&         fa*.5 + caif(it, ia)*(1.-fa*fa)
        xai = xai - caif(it, ia)
        IF (xai .LT. 0.0) THEN
          xai = 0.0
          xaid = 0.0
        ELSE
          xai = xai
        END IF
        IF (xai .GT. 1.0) THEN
          xai = 1.0
          xaid = 0.0
        ELSE
          xai = xai
        END IF
        taudiffd(k, 1) = taucld1d*xai + taucld1*xaid
        taudiff(k, 1) = taucld1*xai
        taudiffd(k, 2) = taucld2d*xai + taucld2*xaid
        taudiff(k, 2) = taucld2*xai
        taudiffd(k, 3) = taucld3d*xai + taucld3*xaid
        taudiff(k, 3) = taucld3*xai
        taudiffd(k, 4) = taucld4d*xai + taucld4*xaid
        taudiff(k, 4) = taucld4*xai
      END IF
    ELSE
! Overlap calculation scaling not needed
      taubeamd(k, 1) = taucld1d
      taubeam(k, 1) = taucld1
      taubeamd(k, 2) = taucld2d
      taubeam(k, 2) = taucld2
      taubeamd(k, 3) = taucld3d
      taubeam(k, 3) = taucld3
      taubeamd(k, 4) = taucld4d
      taubeam(k, 4) = taucld4
      taudiffd(k, 1) = taucld1d
      taudiff(k, 1) = taucld1
      taudiffd(k, 2) = taucld2d
      taudiff(k, 2) = taucld2
      taudiffd(k, 3) = taucld3d
      taudiff(k, 3) = taucld3
      taudiffd(k, 4) = taucld4d
      taudiff(k, 4) = taucld4
    END IF
!-----compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
    ssaclt = 0.99999
    asyclt = 1.0
    taucd = taucld1d + taucld2d + taucld3d + taucld4d
    tauc = taucld1 + taucld2 + taucld3 + taucld4
    IF (tauc .GT. 0.02 .AND. fcld(k) .GT. 0.01) THEN
      w1d = (-(aia_nir(ib, 3)*reffd(k, 1)*reff(k, 1))-(aia_nir(ib, 2)+&
&       aia_nir(ib, 3)*reff(k, 1))*reffd(k, 1))*taucld1 + (1.-(aia_nir(&
&       ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff(k, 1))*reff(k, 1)))*&
&       taucld1d
      w1 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff(k, 1)&
&       )*reff(k, 1)))*taucld1
      w2d = (-(awa_nir(ib, 3)*reffd(k, 2)*reff(k, 2))-(awa_nir(ib, 2)+&
&       awa_nir(ib, 3)*reff(k, 2))*reffd(k, 2))*taucld2 + (1.-(awa_nir(&
&       ib, 1)+(awa_nir(ib, 2)+awa_nir(ib, 3)*reff(k, 2))*reff(k, 2)))*&
&       taucld2d
      w2 = (1.-(awa_nir(ib, 1)+(awa_nir(ib, 2)+awa_nir(ib, 3)*reff(k, 2)&
&       )*reff(k, 2)))*taucld2
      w3d = (1.-ara_nir(ib, 1))*taucld3d
      w3 = (1.-ara_nir(ib, 1))*taucld3
      w4d = (-(aia_nir(ib, 3)*reff_snowd*reff_snow)-(aia_nir(ib, 2)+&
&       aia_nir(ib, 3)*reff_snow)*reff_snowd)*taucld4 + (1.-(aia_nir(ib&
&       , 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff_snow)*reff_snow))*&
&       taucld4d
      w4 = (1.-(aia_nir(ib, 1)+(aia_nir(ib, 2)+aia_nir(ib, 3)*reff_snow)&
&       *reff_snow))*taucld4
      ssacltd = ((w1d+w2d+w3d+w4d)*tauc-(w1+w2+w3+w4)*taucd)/tauc**2
      ssaclt = (w1+w2+w3+w4)/tauc
      g1d = (aig_nir(ib, 3)*reffd(k, 1)*reff(k, 1)+(aig_nir(ib, 2)+&
&       aig_nir(ib, 3)*reff(k, 1))*reffd(k, 1))*w1 + (aig_nir(ib, 1)+(&
&       aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 1))*reff(k, 1))*w1d
      g1 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 1))*&
&       reff(k, 1))*w1
      g2d = (awg_nir(ib, 3)*reffd(k, 2)*reff(k, 2)+(awg_nir(ib, 2)+&
&       awg_nir(ib, 3)*reff(k, 2))*reffd(k, 2))*w2 + (awg_nir(ib, 1)+(&
&       awg_nir(ib, 2)+awg_nir(ib, 3)*reff(k, 2))*reff(k, 2))*w2d
      g2 = (awg_nir(ib, 1)+(awg_nir(ib, 2)+awg_nir(ib, 3)*reff(k, 2))*&
&       reff(k, 2))*w2
      g3d = arg_nir(ib, 1)*w3d
      g3 = arg_nir(ib, 1)*w3
      g4d = (aig_nir(ib, 3)*reffd(k, 4)*reff(k, 4)+(aig_nir(ib, 2)+&
&       aig_nir(ib, 3)*reff(k, 4))*reffd(k, 4))*w4 + (aig_nir(ib, 1)+(&
&       aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 4))*reff(k, 4))*w4d
      g4 = (aig_nir(ib, 1)+(aig_nir(ib, 2)+aig_nir(ib, 3)*reff(k, 4))*&
&       reff(k, 4))*w4
      IF (w1 + w2 + w3 + w4 .NE. 0.0) THEN
        asycltd = ((g1d+g2d+g3d+g4d)*(w1+w2+w3+w4)-(g1+g2+g3+g4)*(w1d+&
&         w2d+w3d+w4d))/(w1+w2+w3+w4)**2
        asyclt = (g1+g2+g3+g4)/(w1+w2+w3+w4)
      ELSE
        asycltd = 0.0
      END IF
    ELSE
      ssacltd = 0.0
      asycltd = 0.0
    END IF
    ssacld(k) = ssacltd
    ssacl(k) = ssaclt
    asycld(k) = asycltd
    asycl(k) = asyclt
  END DO
  RETURN
END SUBROUTINE GETNIRTAU1_D

