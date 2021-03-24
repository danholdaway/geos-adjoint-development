 subroutine sascnvn_tl(im, ix, km, jcap, delt, del, prsl, ps, sl,    &
                      ql, q1,  t1, u1, v1,  &
                      rcs, rn, jkbcon0,jkbcon1,jktcon0,jktcon1,slimsk, &
                      dot, ncloud,        &
                      t_out, q_out, u_out, v_out, ql_out,  &
                      qld_in, q1d_in, t1d_in, u1d_in, v1d_in, dotd,   &
                      t_outd, q_outd, u_outd, v_outd, ql_outd, rnd)
                    

  use kinds, only:r_single,r_kind,i_kind
  use constants, only: CP, HVAP,  RD, FV, T0C, &
                     EPS, EPSM1, G=>grav, FACT1=>FACTOR1, FACT2=>FACTOR2,&
                    ONE, ZERO, HALF, tiny_r_kind, TWO, &
                    EL2ORC,elocp, init_constants_derived, init_constants

  IMPLICIT NONE

  REAL(r_kind) :: mbdt, tem, ptem, ptem1, pgcon
  REAL(r_kind) :: delt, clam, cxlamu, xlamde, xlamdd

  REAL(r_kind), DIMENSION(im) :: cnvscl_1, cnvscl_2, cnvscl_3, cnvscl_4, cnvscl_5,  &
                cnvscl_6, cnvscl_7, cnvscl_8
  REAL(r_kind), DIMENSION(im) :: cnvscl_2d, cnvscl_5d, cnvscl_6d, cnvscl_7d, cnvscl_8d

  REAL(r_kind) :: ps(im), del(ix, km), prsl(ix, km), ql(ix, km, 2),ql_out(ix,km,2), q1(ix, km),     &
                  t1(ix, km), u1(ix, km), v1(ix, km), rcs(im), cldwrk(im), rn(im),  &
                  slimsk(im), dot(ix, km), sl(ix, km), t_out(im, km), q_out(im, km),&
                  u_out(im, km), v_out(im, km)

  REAL(r_kind),intent(in) :: t1d_in(im,km), q1d_in(im,km), u1d_in(im,km),v1d_in(im,km), &
                             qld_in(im,km,2)

  REAL(r_kind) :: qld(ix, km, 2), q1d(ix, km), t1d(ix, km), u1d(ix, km), &
                  v1d(ix, km), rnd(im), t_outd(im, km), q_outd(im, km),  &
                  u_outd(im, km), v_outd(im, km),ql_outd(ix,km,2)

  REAL(r_kind) :: adw, aup, aafac, beta, betal, betas, dellat, desdt,            &
                  deta, detad, dg, dh, dhh, dlnsig, dp, dq, dqsdp, dqsdt,        &
                  dt, dt2, dtmax, dtmin, dv1h, dv1q, dv2h, dv2q, dv1u, dv1v,     & 
                  dv2u, dv2v, dv3q, dv3h, dv3u, dv3v, dz, dz1, e1, edtmax,       &
                  edtmaxl, edtmaxs, es, es0, etah, evef, evfact, evfactl,        &
                  factor, fjcap, fkm, gamma, pprime, qlk, qrch, qs, rain, rfact, &
                  shear, tem1, tem2, val, val1, val2, w1, w1l, w1s, w2, w2l, w2s,&
                  w3, w3l, w3s, w4, w4l, w4s, xdby, xpw, xpwd, xqrch

  REAL(r_kind) :: dellatd, desdtd, dgd, dhd, dhhd, dqd, dqsdpd, dqsdtd, dtd,     &
                  dv1hd, dv1qd, dv2hd, dv2qd, dv1ud, dv1vd, dv2ud, dv2vd, dv3qd, &
                  dv3hd, dv3ud, dv3vd, dzd, dz1d, e1d, esd, es0d, etahd, evefd,  &
                  factord, gammad, pprimed, qlkd, qrchd, qsd, raind, rfactd,     &
                  sheard, tem1d, xdbyd, xpwd0, xpwdd, xqrchd, temd, ptemd, ptem1d

  REAL(r_kind) :: aa1(im), acrt(im), acrtfct(im), delhbar(im), delq(im), delq2(im),&
                  delqbar(im), delqev(im), deltbar(im), deltv(im), dtconv(im),     &
                  edt(im), edto(im), edtx(im), fld(im), hcdo(im, km), hmax(im),    &
                  hmin(im), ucdo(im, km), vcdo(im, km), aa2(im), pbcdif(im),       &
                  pdot(im),po(im, km), pwavo(im), pwevo(im), xlamud(im),           &
                  qcdo(im, km),qcond(im), qevap(im), rntot(im), vshear(im),        &
                  xaa0(im), xk(im), xlamd(im), xmb(im), xmbmax(im), xpwav(im),     &
                  xpwev(im), delubar(im), delvbar(im)

  REAL(r_kind) :: aa1d(im), delqevd(im), edtd(im), edtod(im), edtxd(im), fldd(im), &
                  hcdod(im, km), ucdod(im, km), vcdod(im, km), pod(im, km),        &
                  pwavod(im), pwevod(im), xlamudd(im), qcdod(im, km), qcondd(im),  &
                  qevapd(im), rntotd(im), vsheard(im), xaa0d(im), xkd(im),         &
                  xlamd0d(im), xmbd(im), xpwavd(im), xpwevd(im)

!  physical parameters
  REAL(r_kind), PARAMETER :: cpoel=cp/hvap
  REAL(r_kind), PARAMETER :: tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf)
  real(r_kind) cincr, cincrmax, cincrmin, dum, dthk, cthk, c0, c1, &
                   terr0, delta

!  local variables and arrays
  REAL(r_kind) :: pfld(im, km),to(im, km),qo(im, km),uo(im, km),vo(im, km),qeso(im, km)
  REAL(r_kind) :: tod(im, km), qod(im, km), uod(im, km), vod(im, km), qesod(im, km)
  real(r_kind) :: uo1(im,km),vo1(im,km),uo1d(im,km),vo1d(im,km)
!  cloud water
  REAL(r_kind) :: qlko_ktcon(im), dellal(im, km), tvo(im, km), dbyo(im, km), &
                  zo(im, km), xlamue(im, km), fent1(im, km), fent2(im, km),  &
                  frh(im, km), heo(im, km), heso(im, km), qrcd(im, km),      &
                  dellah(im, km), dellaq(im, km), dellau(im, km), dellav(im, km),&
                  hcko(im, km), ucko(im, km), vcko(im, km), qcko(im, km),    &
                  eta(im, km), etad(im, km), zi(im, km), qrcdo(im, km),      &
                  pwo(im, km), pwdo(im, km), tx1(im), sumx(im)

  REAL(r_kind) :: qlko_ktcond(im), dellald(im, km), tvod(im, km), dbyod(im, km),&
                  zod(im, km), xlamued(im, km), fent1d(im, km), fent2d(im, km), &
                  frhd(im, km), heod(im, km), hesod(im, km), qrcdd(im, km),     &
                  dellahd(im, km), dellaqd(im, km), dellaud(im, km), dellavd(im, km),&
                  hckod(im, km), uckod(im, km), vckod(im, km), qckod(im, km),   &
                  eta1d(im, km), etadd(im, km), zid(im, km), qrcdod(im, km),    &
                  pwod(im, km), pwdod(im, km), sumxd(im),&
                  dotd(im,km),pdotd(im),acrtfctd(im),dtconvd(im)

  REAL(r_kind) :: pcrit(15), acritt(15), acrit(15)

  REAL(r_kind) :: min1, min1d
  REAL(r_kind) :: max1, max2, max3, max4, max5
  REAL(r_kind) :: max1d,max2d,max3d,max4d
  REAL(r_kind) :: y1, y2
  REAL(r_kind) :: y1d,y2d
  REAL(r_kind) :: tempd

  INTEGER(i_kind) :: im, ix, km, jcap, ncloud, kbot(im), ktop(im), kcnv(im)
  INTEGER(i_kind) :: i, j, indx, jmn, k, kk, latd, lond, km1
  INTEGER(i_kind) :: kb(im), kbcon(im), kbcon1(im), ktcon(im), ktcon1(im), &
                     jmin(im), lmin(im), kbmax(im), kbm(im), kmax(im)
      integer(i_kind) :: jkbcon0(im),jkbcon1(im),jktcon0(im),jktcon1(im)

  LOGICAL :: cnvflg(im), flg(im), flg1(im), flg2(im), flg3(im)

  DATA pcrit /850.0_r_kind, 800.0_r_kind, 750.0_r_kind, 700.0_r_kind, 650.0_r_kind, 600.0_r_kind,&
              550.0_r_kind, 500.0_r_kind, 450.0_r_kind, 400.0_r_kind, &
              350.0_r_kind, 300.0_r_kind, 250.0_r_kind, 200.0_r_kind, 150.0_r_kind/
  DATA acritt /.0633_r_kind, .0445_r_kind, .0553_r_kind, .0664_r_kind, .075_r_kind, &
               .1082_r_kind, .1521_r_kind, .2216_r_kind, &
               .3151_r_kind, .3677_r_kind, .41_r_kind, .5255_r_kind, &
               .7663_r_kind, 1.1686_r_kind, 1.6851_r_kind/

!------------------------
  km1 = km - 1

 call init_constants_derived
 call init_constants(.false.)

  terr0=zero
  c0 = 0.002_r_kind
  c1= 0.002_r_kind
  delta = fv
  cthk=150.0_r_kind
  cincrmax=180.0_r_kind
  cincrmin=120.0_r_kind
  dthk=25.0_r_kind


!  initialize arrays
  DO i=1,im
    cnvflg(i) = .true.
    rnd(i)    = zero 
    rn(i)     = zero
    kbot(i)   = km + 1
    ktop(i)   = 0
    kbcon(i)  = km
    ktcon(i)  = 1
    dtconv(i) = 3600.0_r_kind
    dtconvd(i) = zero 
    cldwrk(i) = zero 
    pdot(i)   = zero 
    pdotd(i)   = zero 
    pbcdif(i) = zero 
    lmin(i)   = 1
    jmin(i)   = 1
    qlko_ktcond(i)= zero 
    qlko_ktcon(i) = zero 
    edtd(i)   = zero 
    edt(i)    = zero 
    edtod(i)  = zero 
    edto(i)   = zero
    edtxd(i)  = zero 
    edtx(i)   = zero 
    acrt(i)   = zero 
    acrtfct(i)= one
    acrtfctd(i)= zero 
    aa1d(i)   = zero 
    aa1(i)    = zero 
    aa2(i)    = zero 
    xaa0d(i)  = zero 
    xaa0(i)   = zero 
    pwavod(i) = zero 
    pwavo(i)  = zero 
    pwevod(i) = zero 
    pwevo(i)  = zero
    xpwavd(i) = zero
    xpwav(i)  = zero
    xpwevd(i) = zero
    xpwev(i)  = zero
    vsheard(i)= zero
    vshear(i) = zero
    xlamd0d(i) = zero
    fldd(i) = zero
    xkd(i) = zero
    xmbd(i) = zero
    rntotd(i) = zero
    rntot(i)  = zero
    delqevd(i)= zero
    delqev(i) = zero
    delq2(i)  = zero
    qcond(i) = zero
    qcondd(i) = zero
    qevap(i) = zero
    qevapd(i) = zero
    deltv(i) = zero
    delq(i) = zero
    xlamudd(i) = zero 
    pwavod(i) = zero
    vsheard(i) = zero
    vshear(i)  = zero
    sumxd(i) = zero
    sumx(i)  = zero
  END DO


! Initialize output
  DO k=1,km
    DO i=1,im
      t_outd(i, k) = t1d_in(i,k)
      q_outd(i, k) = q1d_in(i,k)
      u_outd(i, k) = u1d_in(i,k)
      v_outd(i, k) = v1d_in(i,k)
      ql_outd(i,k,1) = qld_in(i,k,1)
      ql_outd(i,k,2) = zero 
      t1d(i, k) = t1d_in(i,k)
      q1d(i, k) = q1d_in(i,k)
      u1d(i, k) = u1d_in(i,k)
      v1d(i, k) = v1d_in(i,k)
      qld(i,k,1) = qld_in(i,k,1)
      qld(i,k,2) = zero
      tod(i,k) = zero
      qod(i,k) = zero
      uod(i,k) = zero
      vod(i,k) = zero
      uo1d(i,k) = zero
      vo1d(i,k) = zero
    END DO
  END DO

  DO i=1,im
    DO k=1,km
      t_out(i, k)  = t1(i,k)
      q_out(i, k)  = q1(i,k)
      u_out(i, k)  = u1(i,k)
      v_out(i, k)  = v1(i,k)
      ql_out(i,k,1) = ql(i,k,1)
      ql_out(i,k,2) = zero 
    END DO
  END DO

  DO k=1,15
    acrit(k) = acritt(k)*(975.0_r_kind-pcrit(k))
  END DO

  dt2 = delt
  val = 1200.0_r_kind

  IF (dt2 .LT. val) THEN
    dtmin = val
  ELSE
    dtmin = dt2
  END IF

  val = 3600.0_r_kind

  IF (dt2 .LT. val) THEN
    dtmax = val
  ELSE
    dtmax = dt2
  END IF

!  model tunable parameters are all here
  mbdt    = 10.0_r_kind
  edtmaxl = .3_r_kind
  edtmaxs = .3_r_kind
  clam    = .1_r_kind
  aafac   = .1_r_kind
  betal   = .05_r_kind
  betas   = .05_r_kind
  evfact  = 0.3_r_kind
  evfactl = 0.3_r_kind

  cxlamu  = 1.0e-4_r_kind
  xlamde  = 1.0e-4_r_kind
  xlamdd  = 1.0e-4_r_kind

! Zhang & Wu (2003,JAS)
  pgcon = 0.55_r_kind
  fjcap = (FLOAT(jcap)/126.0_r_kind)**2
  val = one

  IF (fjcap .LT. val) THEN
    fjcap = val
  ELSE
    fjcap = fjcap
  END IF

  fkm = (FLOAT(km)/28.0_r_kind)**2

  IF (fkm .LT. val) THEN
    fkm = val
  ELSE
    fkm = fkm
  END IF

  w1l = -8.e-3_r_kind
  w2l = -4.e-2_r_kind
  w3l = -5.e-3_r_kind
  w4l = -5.e-4_r_kind
  w1s = -2.e-4_r_kind
  w2s = -2.e-3_r_kind
  w3s = -1.e-3_r_kind
  w4s = -2.e-5_r_kind

!  define top layer for search of the downdraft originating layer
!  and the maximum thetae for updraft

  DO i=1,im
    kbmax(i) = km
    kbm(i)   = km
    kmax(i)  = km
    tx1(i)   = one/ps(i)
  END DO

  DO k=1,km
    DO i=1,im
      IF (prsl(i, k)*tx1(i) .GT. 0.04_r_kind) kmax(i) = k + 1
      IF (prsl(i, k)*tx1(i) .GT. 0.45_r_kind) kbmax(i) = k + 1
      IF (prsl(i, k)*tx1(i) .GT. 0.70_r_kind) kbm(i) = k + 1
    END DO
  END DO

  DO i=1,im
    IF (kbmax(i) .GT. kmax(i)) THEN
      kbmax(i) = kmax(i)
    ELSE
      kbmax(i) = kbmax(i)
    END IF
    IF (kbm(i) .GT. kmax(i)) THEN
      kbm(i)   = kmax(i)
    ELSE
      kbm(i)   = kbm(i)
    END IF
  END DO

!   convert surface pressure to mb from cb
  DO k=1,km
    DO i=1,im
        tod(i, k)   = t1d(i, k)
        to(i, k)    = t1(i, k)
        qod(i, k)   = q1d(i, k)
        qo(i, k)    = q1(i, k)
        uod(i, k)   = rcs(i)*u1d(i, k)
        uo(i, k)    = u1(i, k)*rcs(i)
        vod(i, k)   = rcs(i)*v1d(i, k)
        vo(i, k)    = v1(i, k)*rcs(i)
        pfld(i, k)  = prsl(i, k)*10.0_r_kind
        eta1d(i, k) = zero
        eta(i, k)   = one
        fent1d(i, k)= zero
        fent1(i, k) = one
        fent2d(i, k)= zero
        fent2(i, k) = one
        frhd(i, k)  = zero
        frh(i, k)   = zero
        hckod(i, k) = zero
        hcko(i, k)  = zero
        qckod(i, k) = zero
        qcko(i, k)  = zero
        uckod(i, k) = zero
        ucko(i, k)  = zero
        vckod(i, k) = zero
        vcko(i, k)  = zero
        etadd(i, k) = zero
        etad(i, k)  = one
        hcdod(i, k) = zero
        hcdo(i, k)  = zero
        qcdod(i, k) = zero
        qcdo(i, k)  = zero
        ucdod(i, k) = zero
        ucdo(i, k)  = zero
        vcdod(i, k) = zero
        vcdo(i, k)  = zero
        qrcdd(i, k) = zero
        qrcd(i, k)  = zero
        qrcdod(i, k)= zero
        qrcdo(i, k) = zero
        dbyod(i, k) = zero
        dbyo(i, k)  = zero
        pwod(i, k)  = zero
        pwo(i, k)   = zero
        pwdod(i, k) = zero
        pwdo(i, k)  = zero
        dellald(i,k) = zero
        dellal(i,k)  = zero
        dellahd(i,k) = zero
        dellah(i,k)  = zero
        dellaqd(i,k) = zero
        dellaq(i,k)  = zero
        dellaud(i,k) = zero
        pwod(i,k)    =  zero
        dellau(i,k)  = zero
        dellavd(i,k) = zero
        dellav(i,k)  = zero
        zod(i,k) = zero
        zid(i,k)     = zero
        xlamued(i,k) = zero
        qesod(i,k) = zero
        tvod(i,k) = zero
        hesod(i,k) = zero
        heod(i,k)  = zero
        fent1d(i,k) = zero
        fent2d(i,k) = zero
        qckod(i,k) = zero
    END DO
  END DO


!  column variables
!  p is pressure of the layer (mb)
!  t is temperature at t-dt (k)..tn
!  q is mixing ratio at t-dt (kg/kg)..qn
!  to is temperature at t+dt (k)... this is after advection and turbulan
!  qo is mixing ratio at t+dt (kg/kg)..q1

  DO k=1,km
    DO i=1,im
        call fpvsx_tl(to(i, k), es0, tod(i, k),  es0d) ! fpvs is in pa
        qesod(i, k) = 10.0_r_kind*es0d
        qeso(i, k)  = es0*10.0_r_kind
        qesod(i, k) = (eps*qesod(i, k)*(pfld(i, k)+epsm1*qeso(i, k))-eps&
                      *qeso(i, k)*epsm1*qesod(i, k))/(pfld(i, k)+epsm1*qeso(i, k))**2
        qeso(i, k)  = eps*qeso(i, k)/(pfld(i, k)+epsm1*qeso(i, k))
        val1 = 1.e-8_r_kind
        IF (qeso(i, k) .LT. val1) THEN
          qesod(i, k)= zero
          qeso(i, k) = val1
        ELSE
          qeso(i, k) = qeso(i, k)
        END IF
        val2 = 1.e-10_r_kind
        IF (qo(i, k) .LT. val2) THEN
          qod(i, k)= zero
          qo(i, k) = val2
        ELSE
          qo(i, k) = qo(i, k)
        END IF
        tvod(i, k) = tod(i, k) + delta*(tod(i, k)*qo(i, k)+to(i, k)*qod(i, k))
        tvo(i, k)  = to(i, k) + delta*to(i, k)*qo(i, k)
    END DO
  END DO

!  hydrostatic height assume zero terr and initially assume
!  updraft entrainment rate as an inverse function of height

  DO i=1,im
    zod(i, 1) = -(LOG(sl(i, 1))*rd*tvod(i, 1)/g)
    zo(i, 1)  = zero - LOG(sl(i, 1))*rd/g*tvo(i, 1)
  END DO

  DO k=2,km
    DO i=1,im
      zod(i, k)= zod(i, k-1) - &
                 LOG(sl(i, k)/sl(i, k-1))*rd*half*(tvod(i, k)+tvod(i, k-1))/g
      zo(i, k) = zo(i, k-1) - &
                 LOG(sl(i, k)/sl(i, k-1))*rd/g*half*(tvo(i, k)+tvo(i, k-1))
    END DO
  END DO

  DO k=1,km1
    DO i=1,im
      zid(i, k) = half*(zod(i, k)+zod(i, k+1))
      zi(i, k)  = half*(zo(i, k)+zo(i, k+1))
      xlamued(i, k) = -(clam*zid(i, k)/zi(i, k)**2)
      xlamue(i, k)  =   clam/zi(i, k)

    END DO
  END DO

!  compute moist static energy
  DO k=1,km
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        temd = g*zod(i, k) + cp*tod(i, k)
        tem  = g*zo(i, k) + cp*to(i, k)
        heod(i, k) = temd + hvap*qod(i, k)
        heo(i, k)  = tem + hvap*qo(i, k)
        hesod(i, k)= temd + hvap*qesod(i, k)
        heso(i, k) = tem + hvap*qeso(i, k)
      END IF
    END DO
  END DO

!  determine level with largest moist static energy
!  this is the level where updraft starts

  DO i=1,im
    kb(i) = 1
  END DO

  DO k=1,km1
    DO i=1,im
      IF (k .LE. kmax(i) - 1) THEN
        dz = half*(zo(i, k+1)-zo(i, k))
        dzd = half*(zod(i, k+1)-zod(i, k))
        dp = half*(pfld(i, k+1)-pfld(i, k))
        call fpvsx_tl(to(i, k+1), es0, tod(i, k+1), es0d) ! fpvs is in pa
        esd = 10.0_r_kind*es0d
        es  = es0*10.0_r_kind
        pprimed = epsm1*esd
        pprime  = pfld(i, k+1) + epsm1*es
        qsd = (eps*esd*pprime-eps*es*pprimed)/pprime**2
        qs  = eps*es/pprime

        dqsdpd = -((qsd*pprime-qs*pprimed)/pprime**2)
        dqsdp  = -(qs/pprime)

        desdtd = esd*(fact1/to(i, k+1)+fact2/to(i, k+1)**2) + es*(-(&
                 fact1*tod(i, k+1)/to(i, k+1)**2)-fact2*2*to(i, k+1)*tod(i, k+1)/&
                 (to(i, k+1)**2)**2)
        desdt  = es*(fact1/to(i, k+1)+fact2/to(i, k+1)**2)

        dqsdtd = (pfld(i, k+1)*(qsd*desdt+qs*desdtd)*es*pprime-qs*pfld(i, k+1)*&
                  desdt*(esd*pprime+es*pprimed))/(es*pprime)**2
        dqsdt  = qs*pfld(i, k+1)*desdt/(es*pprime)


        gammad = (el2orc*qesod(i, k+1)*to(i, k+1)**2-el2orc*qeso(i, k+1)*2*&
                  to(i, k+1)*tod(i, k+1))/(to(i, k+1)**2)**2
        gamma = el2orc*qeso(i, k+1)/to(i, k+1)**2

        dtd   = ((g*dzd+hvap*dp*dqsdpd)*cp*(one+gamma)-(g*dz+hvap*dqsdp*dp)*cp*gammad)/&
                 (cp*(one+gamma))**2
        dt  = (g*dz+hvap*dqsdp*dp)/(cp*(one+gamma))

        dqd = dqsdtd*dt + dqsdt*dtd + dp*dqsdpd
        dq  = dqsdt*dt + dqsdp*dp

        tod(i, k) = tod(i, k+1) + dtd
        to(i, k)  = to(i, k+1) + dt

        qod(i, k) = qod(i, k+1) + dqd
        qo(i, k)  = qo(i, k+1) + dq

        pod(i, k) = zero
        po(i, k)  = half*(pfld(i, k)+pfld(i, k+1))

      END IF
    END DO
  END DO


  DO k=1,km1
    DO i=1,im
      IF (k .LE. kmax(i) - 1) THEN
        call fpvsx_tl(to(i, k), es0, tod(i, k), es0d)
        qesod(i, k) = 10.0_r_kind*es0d
        qeso(i, k)  = es0*10.0_r_kind
        qesod(i, k) = (eps*qesod(i, k)*(po(i, k)+epsm1*qeso(i, k))-eps*&
                       qeso(i, k)*epsm1*qesod(i, k))/(po(i, k)+epsm1*qeso(i, k))**2
        qeso(i, k)  = eps*qeso(i, k)/(po(i, k)+epsm1*qeso(i, k))
        val1 = 1.e-8_r_kind
        IF (qeso(i, k) .LT. val1) THEN
          qesod(i, k) = zero
          qeso(i, k)  = val1
        ELSE
          qeso(i, k)  = qeso(i, k)
        END IF
        val2 = 1.e-10_r_kind
        IF (qo(i, k) .LT. val2) THEN
          qod(i, k) = zero
          qo(i, k)  = val2
        ELSE
          qo(i, k) = qo(i, k)
        END IF

        IF (qo(i, k)/qeso(i, k) .GT. one) THEN
          min1  = one
          min1d = zero
        ELSE
          min1d = (qod(i, k)*qeso(i, k)-qo(i, k)*qesod(i, k))/qeso(i, k)**2
          min1 = qo(i, k)/qeso(i, k)
        END IF

        frhd(i,k) = -min1d
        frh(i,k)  = one - min1

        heod(i,k) = half*g*(zod(i, k)+zod(i, k+1)) + cp*tod(i, k) + hvap*qod(i, k)
        heo(i,k)  = half*g*(zo(i, k)+zo(i, k+1)) + cp*to(i, k) + hvap*qo(i, k)

        hesod(i,k)= half*g*(zod(i, k)+zod(i, k+1)) + cp*tod(i, k) + hvap*qesod(i, k)
        heso(i,k) = half*g*(zo(i, k)+zo(i, k+1)) + cp*to(i, k) + hvap*qeso(i, k)

        uod(i,k)  = half*(uod(i, k)+uod(i, k+1))
        uo(i,k)   = half*(uo(i, k)+uo(i, k+1))
        vod(i,k)  = half*(vod(i, k)+vod(i, k+1))
        vo(i,k)   = half*(vo(i, k)+vo(i, k+1))
      END IF
    END DO
  END DO

!  look for the level of free convection as cloud base

  DO i=1,im
    flg(i) = .true.
    cnvscl_1(i) = one
    cnvscl_2(i) = one
    cnvscl_2d(i) = zero
    cnvscl_3(i) = one
    cnvscl_4(i) = one
    cnvscl_5d(i)= zero
    cnvscl_5(i) = one
    cnvscl_6d(i)= zero
    cnvscl_6(i) = one
    cnvscl_7d(i)= zero
    cnvscl_7(i) = one
    cnvscl_8d(i)= zero
    cnvscl_8(i) = one
    kbcon(i)    = kmax(i)
  END DO

  DO k=1,km1
    DO i=1,im
      IF (flg(i) .AND. k .LE. kbmax(i)) THEN
        IF (k .GT. kb(i) .AND. heo(i, kb(i)) .GT. heso(i, k)) THEN
          kbcon(i) = k
          flg(i)   = .false.
        END IF
      END IF
    END DO
  END DO

 do i = 1, im
    kbcon(i) = jkbcon0(i)
 enddo

 DO i=1,im
    IF (kbcon(i) .EQ. kmax(i)) cnvscl_1(i) = tiny_r_kind
 END DO


!  determine critical convective inhibition

  DO i=1,im
    pdot(i) = 10.0_r_kind*dot(i, kbcon(i))
    pdotd(i) = 10.0_r_kind*dotd(i, kbcon(i))
  END DO


  DO i=1,im
    IF (slimsk(i) .EQ. one) THEN
      w1 = w1l
      w2 = w2l
      w3 = w3l
      w4 = w4l
    ELSE
      w1 = w1s
      w2 = w2s
      w3 = w3s
      w4 = w4s
    END IF

    IF (pdot(i) .LE. w4) THEN
      tem = (pdot(i)-w4)/(w3-w4)
    ELSE IF (pdot(i) .GE. -w4) THEN
      tem = -((pdot(i)+w4)/(w4-w3))
    ELSE
      tem = zero
    END IF

    val1 = -one
    IF (tem .LT. val1) THEN
      tem = val1
    ELSE
      tem = tem
    END IF

    val2 = one
    IF (tem .GT. val2) THEN
      tem = val2
    ELSE
      tem = tem
    END IF

    tem  = one - tem
    tem1 = half*(cincrmax-cincrmin)
    cincr = cincrmax - tem*tem1
    pbcdif(i) = pfld(i, kb(i)) - pfld(i, kbcon(i))
    if(pbcdif(i).gt.cincr) then
        cnvscl_2(i) = tiny_r_kind
        cnvscl_2d(i) = zero
    endif
    cnvscl_2d(i) = cnvscl_2d(i)*cnvscl_1(i)
    cnvscl_2(i) = cnvscl_2(i)*cnvscl_1(i)
  END DO


!  assume that updraft entrainment rate above cloud base is
!    same as that at cloud base

  DO k=2,km1
    DO i=1,im
      IF (k .GT. kbcon(i) .AND. k .LT. kmax(i)) THEN
        xlamued(i, k) = xlamued(i, kbcon(i))
        xlamue(i, k) = xlamue(i, kbcon(i))
      END IF
    END DO
  END DO


!  assume the detrainment rate for the updrafts to be same as
!  the entrainment rate at cloud base

  DO i=1,im
    xlamudd(i) = xlamued(i, kbcon(i))
    xlamud(i)  = xlamue(i, kbcon(i))
  END DO

!  functions rapidly decreasing with height, mimicking a cloud ensemble
!    (Bechtold et al., 2008)

  DO k=2,km1
    DO i=1,im
      IF (k .GT. kbcon(i) .AND. k .LT. kmax(i)) THEN
        temd=(qesod(i, k)*qeso(i, kbcon(i))-qeso(i, k)*qesod(i, kbcon(i)))/&
              qeso(i, kbcon(i))**2
        tem = qeso(i, k)/qeso(i, kbcon(i))
        fent1d(i, k) = 2.0_r_kind*tem*temd
        fent1(i, k)  = tem**2
        fent2d(i, k) = 3.0_r_kind*tem**2*temd
        fent2(i, k)  = tem**3
      END IF
    END DO
  END DO

!  final entrainment rate as the sum of turbulent part and organized entrainment
!    depending on the environmental relative humidity
!    (Bechtold et al., 2008)
  DO k=2,km1
    DO i=1,im
      IF (k .GE. kbcon(i) .AND. k .LT. kmax(i)) THEN
        temd = cxlamu*(frhd(i, k)*fent2(i, k)+frh(i, k)*fent2d(i, k))
        tem  = cxlamu*frh(i, k)*fent2(i, k)
        xlamued(i, k)= xlamued(i, k)*fent1(i, k)+ xlamue(i, k)*fent1d(i, k) + temd
        xlamue(i, k) = xlamue(i, k)*fent1(i, k) + tem
      END IF
    END DO
  END DO


!  determine updraft mass flux for the subcloud layers
  DO k=km1,1,-1
    DO i=1,im
      IF (k .LT. kbcon(i) .AND. k .GE. kb(i)) THEN
        dzd = zid(i, k+1) - zid(i, k)
        dz  = zi(i, k+1)  - zi(i, k)
        ptemd = half*(xlamued(i, k)+xlamued(i, k+1)) - xlamudd(i)
        ptem  = half*(xlamue(i, k)+xlamue(i, k+1)) - xlamud(i)
        eta1d(i, k) = (eta1d(i, k+1)*(one+ptem*dz)-eta(i, k+1)*(ptemd*dz+ptem*dzd))/&
                      (one+ptem*dz)**2
        eta(i, k)   = eta(i, k+1)/(one+ptem*dz)
      END IF
    END DO
  END DO

!  compute mass flux above cloud base
  DO k=2,km1
    DO i=1,im
      IF (k .GT. kbcon(i) .AND. k .LT. kmax(i)) THEN
        dzd = zid(i, k) - zid(i, k-1)
        dz  = zi(i, k)  - zi(i, k-1)
        ptemd = half*(xlamued(i, k)+xlamued(i, k-1)) - xlamudd(i)
        ptem  = half*(xlamue(i, k)+xlamue(i, k-1)) - xlamud(i)
        eta1d(i, k) = eta1d(i, k-1)*(one+ptem*dz) + eta(i, k-1)*(ptemd*dz+ptem*dzd)
        eta(i, k)   = eta(i, k-1)*(one+ptem*dz)
      END IF
    END DO
  END DO

!  compute updraft cloud properties
  DO i=1,im
    indx = kb(i)
    hckod(i, indx) = heod(i, indx)
    hcko(i, indx)  = heo(i, indx)
    uckod(i, indx) = uod(i, indx)
    ucko(i, indx)  = uo(i, indx)
    vckod(i, indx) = vod(i, indx)
    vcko(i, indx)  = vo(i, indx)
    pwavo(i)  = zero
  END DO


!  cloud property is modified by the entrainment process
  DO k=2,km1
    DO i=1,im
      IF (k .GT. kb(i) .AND. k .LT. kmax(i)) THEN
        dzd = zid(i, k) - zid(i, k-1)
        dz  = zi(i, k)  - zi(i, k-1)
        temd = half*((xlamued(i,k)+xlamued(i,k-1))*dz+(xlamue(i,k)+xlamue(i,k-1))*dzd)
        tem  = half*(xlamue(i,k)+xlamue(i,k-1))*dz
        tem1d= half*(xlamudd(i)*dz+xlamud(i)*dzd)
        tem1 = half*xlamud(i)*dz
        factord= temd - tem1d
        factor = one + tem - tem1
        ptemd  = half*temd
        ptem   = half*tem + pgcon
        ptem1d = half*temd
        ptem1  = half*tem - pgcon
        hckod(i, k)= (((one-tem1)*hckod(i, k-1)-tem1d*hcko(i, k-1)+half*(temd* &
                     (heo(i, k)+heo(i, k-1))+tem*(heod(i, k)+heod(i, k-1))))*&
                      factor-((one-tem1)*hcko(i, k-1)+tem*half*(heo(i, k)+     &
                      heo(i, k-1)))*factord)/factor**2
        hcko(i, k) = ((one-tem1)*hcko(i,k-1)+tem*half*(heo(i,k)+heo(i,k-1)))/factor

        uckod(i, k)= (((one-tem1)*uckod(i, k-1)-tem1d*ucko(i, k-1)+ptemd*&
                      uo(i, k)+ptem*uod(i, k)+ptem1d*uo(i, k-1)+ptem1*  &
                      uod(i, k-1))*factor-((one-tem1)*ucko(i, k-1)+ptem* &
                      uo(i, k)+ptem1*uo(i, k-1))*factord)/factor**2
        ucko(i, k) = ((one-tem1)*ucko(i,k-1)+ptem*uo(i,k)+ptem1*uo(i,k-1))/factor

        vckod(i, k)= (((one-tem1)*vckod(i, k-1)-tem1d*vcko(i, k-1)+ptemd*vo(i, k)+&
                      ptem*vod(i, k)+ptem1d*vo(i, k-1)+ptem1*vod(i, k-1))*factor-&
                      ((one-tem1)*vcko(i, k-1)+ptem*vo(i, k)+ptem1*vo(i, k-1))*   &
                      factord)/factor**2
        vcko(i, k) = ((one-tem1)*vcko(i,k-1)+ptem*vo(i,k)+ptem1*vo(i,k-1))/factor

        dbyod(i, k)= hckod(i, k) - hesod(i, k)
        dbyo(i, k) = hcko(i, k) - heso(i, k)
      END IF
    END DO
  END DO

!   taking account into convection inhibition due to existence of
!    dry layers below cloud base

  DO i=1,im
    kbcon1(i) = kmax(i)
    flg1(i)   = .true.
  END DO

  DO k=2,km1
    DO i=1,im
      IF (k .LT. kmax(i)) THEN
        IF (flg1(i) .AND. k .GE. kbcon(i) .AND. dbyo(i, k) .GT. zero) THEN
          kbcon1(i) = k
          flg1(i)   = .false.
        END IF
      END IF
    END DO
  END DO

 do i = 1, im
   kbcon1(i) = jkbcon1(i)
 enddo


 DO i=1,im
  IF (kbcon1(i) .EQ. kmax(i)) THEN
    cnvscl_2(i) = tiny_r_kind
    cnvscl_2d(i) = zero
  END IF
 END DO

 DO i=1,im
   tem   = pfld(i, kbcon(i)) - pfld(i, kbcon1(i))
   cnvscl_3(i) = one - (half**((dthk/tem)**200))
   cnvscl_3(i) = cnvscl_3(i)*cnvscl_2(i)
 END DO


!  determine first guess cloud top as the level of zero buoyancy

  DO i=1,im
    ktcon(i)= 1
    flg2(i) = .true.
  END DO

  DO k=2,km1
    DO i=1,im
      IF (k .LT. kmax(i)) THEN
        IF (flg2(i) .AND. k .GT. kbcon1(i) .AND. dbyo(i, k) .LT. zero) THEN
          ktcon(i)= k
          flg2(i) = .false.
        END IF
      END IF
    END DO
  END DO

  do i = 1, im
    ktcon(i) = jktcon0(i)
  enddo

if(ktcon(1) .gt. 2) then

  DO i=1,im
    tem  = pfld(i, kbcon(i)) - pfld(i, ktcon(i))
    cnvscl_4(i) = one - (half**((tem/cthk)**200))
    cnvscl_4(i) = cnvscl_4(i)*cnvscl_3(i)
  END DO


!  search for downdraft originating level above theta-e minimum

  DO i=1,im
    hmin(i) = heo(i, kbcon1(i))
    lmin(i) = kbmax(i)
    jmin(i) = kbmax(i)
  END DO

  DO k=2,km1
    DO i=1,im
      IF (k .LE. kbmax(i)) THEN
        IF (k .GT. kbcon1(i) .AND. heo(i, k) .LT. hmin(i)) THEN
          lmin(i) = k + 1
          hmin(i) = heo(i, k)
        END IF
      END IF
    END DO
  END DO

!  make sure that jmin(i) is within the cloud

  DO i=1,im
    IF (lmin(i) .GT. ktcon(i) - 1) THEN
      jmin(i) = ktcon(i) - 1
    ELSE
      jmin(i) = lmin(i)
    END IF
    IF (jmin(i) .LT. kbcon1(i) + 1) THEN
      jmin(i) = kbcon1(i) + 1
    ELSE
      jmin(i) = jmin(i)
    END IF
    IF (jmin(i) .GE. ktcon(i)) THEN
      cnvscl_4(i) = tiny_r_kind
    END IF
  END DO


!  specify upper limit of mass flux at cloud base
  DO i=1,im
    k = kbcon(i)
    dp = 1000.0_r_kind*del(i, k)*ps(i)
    xmbmax(i) = dp/(g*dt2)
  END DO


!  compute cloud moisture property and precipitation

  DO i=1,im
    qckod(i, kb(i))= qod(i, kb(i))
    qcko(i, kb(i)) = qo(i, kb(i))
  END DO


  DO k=2,km1
    DO i=1,im
      IF (k .GT. kb(i) .AND. k .LT. ktcon(i)) THEN
        dzd = zid(i, k) - zid(i, k-1)
        dz  = zi(i, k) - zi(i, k-1)
        gammad= (el2orc*qesod(i, k)*to(i, k)**2-el2orc*qeso(i, k)*2*to(i, k)*&
                 tod(i, k))/(to(i, k)**2)**2
        gamma = el2orc*qeso(i, k)/to(i, k)**2
        qrchd = qesod(i, k) + ((gammad*dbyo(i, k)+gamma*dbyod(i, k))*hvap*&
                (one+gamma)-gamma*dbyo(i, k)*hvap*gammad)/(hvap*(one+gamma))**2
        qrch  = qeso(i, k) + gamma*dbyo(i, k)/(hvap*(one+gamma))
        temd  = half*((xlamued(i,k)+xlamued(i,k-1))*dz+(xlamue(i,k)+xlamue(i,k-1))*dzd)
        tem   = half*(xlamue(i, k)+xlamue(i, k-1))*dz
        tem1d = half*(xlamudd(i)*dz+xlamud(i)*dzd)
        tem1  = half*xlamud(i)*dz
        factord= temd - tem1d
        factor = one + tem - tem1
        qckod(i, k)= (((one-tem1)*qckod(i, k-1)-tem1d*qcko(i, k-1)+half*(temd*      &
                       (qo(i, k)+qo(i, k-1))+tem*(qod(i, k)+qod(i, k-1))))*factor-&
                       ((one-tem1)*qcko(i, k-1)+tem*half*(qo(i, k)+qo(i, k-1)))*    &
                        factord)/factor**2
        qcko(i, k) = ((one-tem1)*qcko(i, k-1)+tem*half*(qo(i, k)+qo(i, k-1)))/factor
        dqd = eta1d(i, k)*(qcko(i, k)-qrch) + eta(i, k)*(qckod(i, k)-qrchd)
        dq  = eta(i, k)*(qcko(i, k)-qrch)

!  check if there is excess moisture to release latent heat
        IF (k .GE. kbcon(i) .AND. dq .GT. zero) THEN
          etahd = half*(eta1d(i, k)+eta1d(i, k-1))
          etah  = half*(eta(i, k)+eta(i, k-1))
          IF (ncloud .GT. zero .AND. k .GT. jmin(i)) THEN
            dp   = 1000.0_r_kind*del(i, k)*ps(i)
            qlkd = (dqd*(eta(i, k)+etah*(c0+c1)*dz)-dq*(eta1d(i, k)+(c0+c1)*&
                   (etahd*dz+etah*dzd)))/(eta(i, k)+etah*(c0+c1)*dz)**2
            qlk  = dq/(eta(i, k)+etah*(c0+c1)*dz)
            dellald(i, k) = c1*g*((etahd*dz+etah*dzd)*qlk+etah*dz*qlkd)/dp
            dellal(i, k)  = etah*c1*dz*qlk*g/dp
          ELSE
            qlkd= (dqd*(eta(i, k)+etah*c0*dz)-dq*(eta1d(i, k)+c0*(etahd*dz+&
                   etah*dzd)))/(eta(i, k)+etah*c0*dz)**2
            qlk = dq/(eta(i, k)+etah*c0*dz)
          END IF
          aa1d(i) = aa1d(i) - g*(dzd*qlk+dz*qlkd)
          aa1(i)  = aa1(i) - dz*g*qlk
          qckod(i, k)= qlkd + qrchd
          qcko(i, k) = qlk + qrch
          pwod(i, k) = c0*((etahd*dz+etah*dzd)*qlk+etah*dz*qlkd)
          pwo(i, k)  = etah*c0*dz*qlk
          pwavod(i)  = pwavod(i) + pwod(i, k)
          pwavo(i)   = pwavo(i) + pwo(i, k)
        ENDIF
      END IF
    END DO
  END DO

!  calculate cloud work function
  DO k=2,km1
    DO i=1,im
      IF (k .GE. kbcon(i) .AND. k .LT. ktcon(i)) THEN
        dz1d = zod(i, k+1) - zod(i, k)
        dz1  = zo(i, k+1) - zo(i, k)
        gammad = (el2orc*qesod(i, k)*to(i, k)**2-el2orc*qeso(i, k)*2*to(i, k)*&
                  tod(i, k))/(to(i, k)**2)**2
        gamma  = el2orc*qeso(i, k)/to(i, k)**2
        rfactd = delta*cp*(gammad*to(i, k)+gamma*tod(i, k))/hvap
        rfact  = one + delta*cp*gamma*to(i, k)/hvap
        aa1d(i)= aa1d(i) + (((dz1d*dbyo(i, k)+dz1*dbyod(i, k))*g/(cp*to(i, k))-&
                 dz1*dbyo(i, k)*g*tod(i, k)/(cp*to(i, k)**2))*(one+gamma)-dz1*g*&
                 dbyo(i, k)*gammad/(cp*to(i, k)))*rfact/(one+gamma)**2 + dz1*g* &
                 dbyo(i, k)*rfactd/(cp*to(i, k)*(one+gamma))
        aa1(i) = aa1(i) + dz1*(g/(cp*to(i, k)))*dbyo(i, k)/(one+gamma)*rfact
        val    = zero
        IF (val .LT. qeso(i, k) - qo(i, k)) THEN
          max1d = qesod(i, k) - qod(i, k)
          max1  = qeso(i, k) - qo(i, k)
        ELSE
          max1  = val
          max1d = zero
        END IF
        aa1d(i) = aa1d(i) + g*delta*(dz1d*max1+dz1*max1d)
        aa1(i)  = aa1(i) + dz1*g*delta*max1

      END IF
    END DO
  END DO

  DO i=1,im
    cnvscl_5d(i)= half*aa1d(i)*(one-TANH(aa1(i))**2)
    cnvscl_5(i) = half*(one+TANH(aa1(i)))
    cnvscl_5d(i)= cnvscl_4(i)*cnvscl_5d(i)
    cnvscl_5(i) = cnvscl_5(i)*cnvscl_4(i)
  END DO

!  estimate the onvective overshooting as the level
!    where the [aafac * cloud work function] becomes zero,
!    which is the final cloud top
  DO i=1,im
    ktcon1(i) =  jktcon1(i)
  END DO

!  compute cloud moisture property, detraining cloud water
!    and precipitation in overshooting layers
  DO k=2,km1
    DO i=1,im
      IF (k .GE. ktcon(i) .AND. k .LT. ktcon1(i)) THEN
        dzd = zid(i, k) - zid(i, k-1)
        dz  = zi(i, k) - zi(i, k-1)
        gammad= (el2orc*qesod(i, k)*to(i, k)**2-el2orc*qeso(i, k)*2*to(i, k)*&
                 tod(i, k))/(to(i, k)**2)**2
        gamma = el2orc*qeso(i, k)/to(i, k)**2
        qrchd = qesod(i, k) + ((gammad*dbyo(i, k)+gamma*dbyod(i, k))*hvap*&
                (one+gamma)-gamma*dbyo(i, k)*hvap*gammad)/(hvap*(one+gamma))**2
        qrch = qeso(i, k) + gamma*dbyo(i, k)/(hvap*(one+gamma))
        temd = half*((xlamued(i,k)+xlamued(i,k-1))*dz+(xlamue(i,k)+xlamue(i, k-1))*dzd)
        tem  = half*(xlamue(i, k)+xlamue(i, k-1))*dz
        tem1d= half*(xlamudd(i)*dz+xlamud(i)*dzd)
        tem1 = half*xlamud(i)*dz
        factord= temd - tem1d
        factor = one + tem - tem1
        qckod(i, k) = (((one-tem1)*qckod(i,k-1)-tem1d*qcko(i, k-1)+half*(temd*     &
                       (qo(i, k)+qo(i,k-1))+tem*(qod(i, k)+qod(i, k-1))))*factor-&
                       ((one-tem1)*qcko(i, k-1)+tem*half*(qo(i, k)+qo(i, k-1)))*    &
                        factord)/factor**2
        qcko(i, k)  = ((one-tem1)*qcko(i, k-1)+tem*half*(qo(i, k)+qo(i, k-1)))/factor
        dqd = eta1d(i, k)*(qcko(i, k)-qrch) + eta(i, k)*(qckod(i, k)-qrchd)
        dq  = eta(i, k)*(qcko(i, k)-qrch)
!  check if there is excess moisture to release latent heat

        IF (dq .GT. zero) THEN
          etahd = half*(eta1d(i, k)+eta1d(i, k-1))
          etah  = half*(eta(i, k)+eta(i, k-1))
          IF (ncloud .GT. zero) THEN
            dp   = 1000.0_r_kind*del(i, k)*ps(i)
            qlkd = (dqd*(eta(i, k)+etah*(c0+c1)*dz)-dq*(eta1d(i, k)+(c0+c1)*&
                   (etahd*dz+etah*dzd)))/(eta(i, k)+etah*(c0+c1)*dz)**2
            qlk  = dq/(eta(i, k)+etah*(c0+c1)*dz)
            dellald(i, k) = c1*g*((etahd*dz+etah*dzd)*qlk+etah*dz*qlkd)/dp
            dellal(i, k)  = etah*c1*dz*qlk*g/dp
          ELSE
            qlkd= (dqd*(eta(i, k)+etah*c0*dz)-dq*(eta1d(i, k)+c0*(etahd*dz+&
                   etah*dzd)))/(eta(i, k)+etah*c0*dz)**2
            qlk = dq/(eta(i, k)+etah*c0*dz)
          END IF
          qckod(i, k) = qlkd + qrchd
          qcko(i, k)  = qlk + qrch
          pwod(i, k)  = c0*((etahd*dz+etah*dzd)*qlk+etah*dz*qlkd)
          pwo(i, k)   = etah*c0*dz*qlk
          pwavod(i)   = pwavod(i) + pwod(i, k)
          pwavo(i)    = pwavo(i) + pwo(i, k)
        END IF
      ENDIF
    END DO
  END DO

! exchange ktcon with ktcon1
  DO i=1,im
    kk = ktcon(i)
    ktcon(i)  = ktcon1(i)
    ktcon1(i) = kk
  END DO

!  this section is ready for cloud water

  IF (ncloud .GT. 0) THEN

!  compute liquid and vapor separation at cloud top

    DO i=1,im
      k = ktcon(i) - 1

      gammad = el2orc*qesod(i, k)/to(i, k)**2 - el2orc*qeso(i,k)*2*tod(i,k)/to(i,k)**3
      gamma  = el2orc*qeso(i, k)/to(i, k)**2
      qrchd  = qesod(i, k) + ((gammad*dbyo(i, k)+gamma*dbyod(i, k))*hvap*&
               (one+gamma)-gamma*dbyo(i, k)*hvap*gammad)/(hvap*(one+gamma))**2
      qrch = qeso(i, k) + gamma*dbyo(i, k)/(hvap*(one+gamma))
      dqd  = qckod(i, k) - qrchd
      dq   = qcko(i, k) - qrch

!  check if there is excess moisture to release latent heat

      IF (dq .GT. zero) THEN
        qlko_ktcond(i) = dqd
        qlko_ktcon(i)  = dq
        qckod(i, k)    = qrchd
        qcko(i, k)     = qrch
      END IF
    END DO
  END IF

!------- downdraft calculations
!--- compute precipitation efficiency in terms of windshear
  DO k=2,km
    DO i=1,im
      IF (k .GT. kb(i) .AND. k .LE. ktcon(i)) THEN
        IF ((uo(i, k)-uo(i, k-1))**2 + (vo(i, k)-vo(i, k-1))**2 .EQ. zero) THEN
          sheard = zero
        ELSE
          sheard = (2*(uo(i, k)-uo(i, k-1))*(uod(i, k)-uod(i, k-1)) + 2*(vo(i, k)-&
                    vo(i, k-1))*(vod(i, k)-vod(i, k-1)))/(2.0*SQRT((uo(i, k)-     &
                    uo(i, k-1))**2 + (vo(i, k)-vo(i, k-1))**2))
        END IF
        shear    = SQRT((uo(i, k)-uo(i, k-1))**2 + (vo(i, k)-vo(i, k-1))**2)
        vsheard(i)= vsheard(i) + sheard
        vshear(i) = vshear(i) + shear
      END IF
    END DO
  END DO

  DO i=1,im
    vsheard(i) = (1.e3_r_kind*vsheard(i)*(zi(i, ktcon(i))-zi(i, kb(i)))-1.e3_r_kind*vshear(i)*&
                 (zid(i, ktcon(i))-zid(i, kb(i))))/(zi(i, ktcon(i))-zi(i, kb(i)))**2
    vshear(i)  = 1.e3_r_kind*vshear(i)/(zi(i, ktcon(i))-zi(i, kb(i)))
    e1d= .0953_r_kind*2*vshear(i)*vsheard(i) - .639_r_kind*vsheard(i) - .00496_r_kind*3*vshear(i)**2*&
          vsheard(i)
    e1 = 1.591_r_kind - .639_r_kind*vshear(i) + .0953_r_kind*vshear(i)**2 - .00496_r_kind*vshear(i)**3
    edtd(i) = -e1d
    edt(i)  = one - e1
    val     = .9_r_kind
    IF (edt(i) .GT. val) THEN
      edtd(i) = zero
      edt(i)  = val
    ELSE
      edt(i)  = edt(i)
    END IF
    val = zero
    IF (edt(i) .LT. val) THEN
      edtd(i) = zero
      edt(i)  = val
    ELSE
      edt(i)  = edt(i)
    END IF
    edtod(i)  = edtd(i)
    edto(i)   = edt(i)
    edtxd(i)  = edtd(i)
    edtx(i)   = edt(i)
  END DO

!  determine detrainment rate between 1 and kbcon


  DO k=1,km1
    DO i=1,im
      IF (k .GE. 1 .AND. k .LT. kbcon(i)) THEN
        dzd = zid(i, k+1) - zid(i, k)
        dz  = zi(i, k+1) - zi(i, k)
        sumxd(i) = sumxd(i) + dzd
        sumx(i)  = sumx(i) + dz
      END IF
    END DO
  END DO

  DO i=1,im
    beta = betas
    IF (slimsk(i) .EQ. one) beta = betal
    dzd = (sumxd(i)+zid(i, 1))/FLOAT(kbcon(i))
    dz  = (sumx(i)+zi(i, 1))/FLOAT(kbcon(i))
    tem = one/FLOAT(kbcon(i))
    xlamd0d(i) = -((one-beta**tem)*dzd/dz**2)
    xlamd(i)   =   (one-beta**tem)/dz
  END DO

!  determine downdraft mass flux

  DO k=km1,1,-1
    DO i=1,im
      IF (k .LE. kmax(i) - 1) THEN
        IF (k .LT. jmin(i) .AND. k .GE. kbcon(i)) THEN
          dzd = zid(i, k+1) - zid(i, k)
          dz  = zi(i, k+1) - zi(i, k)
          ptem = xlamdd - xlamde
          etadd(i, k)= etadd(i, k+1)*(one-ptem*dz) - etad(i, k+1)*ptem*dzd
          etad(i, k) = etad(i, k+1)*(one-ptem*dz)
        ELSE IF (k .LT. kbcon(i)) THEN
          dzd = zid(i, k+1) - zid(i, k)
          dz  = zi(i, k+1) - zi(i, k)
          ptemd = xlamd0d(i)
          ptem  = xlamd(i) + xlamdd - xlamde
          etadd(i, k) = etadd(i, k+1)*(one-ptem*dz) + etad(i, k+1)*(-(ptemd*dz)-&
                        ptem*dzd)
          etad(i, k)  = etad(i, k+1)*(one-ptem*dz)
        END IF
      END IF
    END DO
  END DO


!--- downdraft moisture properties
  DO i=1,im
    jmn = jmin(i)
    hcdod(i, jmn) = heod(i, jmn)
    hcdo(i, jmn)  = heo(i, jmn)
    qcdod(i, jmn) = qod(i, jmn)
    qcdo(i, jmn)  = qo(i, jmn)
    qrcdod(i, jmn)= qesod(i, jmn)
    qrcdo(i, jmn) = qeso(i, jmn)
    ucdod(i, jmn) = uod(i, jmn)
    ucdo(i, jmn)  = uo(i, jmn)
    vcdod(i, jmn) = vod(i, jmn)
    vcdo(i, jmn)  = vo(i, jmn)
  END DO

  DO k=km1,1,-1
    DO i=1,im
      IF (k .LT. jmin(i)) THEN
        dzd = zid(i, k+1) - zid(i, k)
        dz  = zi(i, k+1) - zi(i, k)
        IF (k .GE. kbcon(i)) THEN
          temd = xlamde*dzd
          tem  = xlamde*dz
          tem1d= half*xlamdd*dzd
          tem1 = half*xlamdd*dz
        ELSE
          temd = xlamde*dzd
          tem  = xlamde*dz
          tem1d= half*(xlamd0d(i)*dz+(xlamd(i)+xlamdd)*dzd)
          tem1 = half*(xlamd(i)+xlamdd)*dz
        END IF
        factord= temd - tem1d
        factor = one + tem - tem1
        ptemd  = half*temd
        ptem   = half*tem - pgcon
        ptem1d = half*temd
        ptem1  = half*tem + pgcon
        hcdod(i, k) = (((one-tem1)*hcdod(i, k+1)-tem1d*hcdo(i, k+1)+half*(temd*&
                      (heo(i, k)+heo(i, k+1))+tem*(heod(i, k)+heod(i, k+1))))*factor-&
                      ((one-tem1)*hcdo(i, k+1)+tem*half*(heo(i, k)+heo(i, k+1)))*&
                      factord)/factor**2
        hcdo(i, k)  = ((one-tem1)*hcdo(i, k+1)+tem*half*(heo(i, k)+heo(i, k+1)))/factor
        ucdod(i, k) = (((one-tem1)*ucdod(i, k+1)-tem1d*ucdo(i, k+1)+ptemd*&
                      uo(i, k+1)+ptem*uod(i, k+1)+ptem1d*uo(i, k)+ptem1*uod(i, k))*&
                      factor-((one-tem1)*ucdo(i, k+1)+ptem*uo(i, k+1)+ptem1*uo(i, k))*&
                      factord)/factor**2
        ucdo(i, k)  = ((one-tem1)*ucdo(i, k+1)+ptem*uo(i, k+1)+ptem1*uo(i, k))/factor
        vcdod(i, k) = (((one-tem1)*vcdod(i, k+1)-tem1d*vcdo(i, k+1)+ptemd*&
                      vo(i, k+1)+ptem*vod(i, k+1)+ptem1d*vo(i, k)+ptem1*vod(i, k))*&
                      factor-((one-tem1)*vcdo(i, k+1)+ptem*vo(i, k+1)+ptem1*vo(i, k))*&
                      factord)/factor**2
        vcdo(i, k)  = ((one-tem1)*vcdo(i, k+1)+ptem*vo(i, k+1)+ptem1*vo(i, k))/factor
        dbyod(i, k) = hcdod(i, k) - hesod(i, k)
        dbyo(i, k)  = hcdo(i, k) - heso(i, k)
      END IF
    END DO
  END DO

  DO k=km1,1,-1
    DO i=1,im
      IF (k .LT. jmin(i)) THEN
        gammad = (el2orc*qesod(i, k)*to(i, k)**2-el2orc*qeso(i, k)*2*to(i, k)*&
                 tod(i, k))/(to(i, k)**2)**2
        gamma  = el2orc*qeso(i, k)/to(i, k)**2


        qrcdod(i, k) = qesod(i, k) + ((gammad*(one+gamma)-gamma*gammad)*dbyo(i, k)/&
                       (one+gamma)**2+gamma*dbyod(i, k)/(one+gamma))/hvap
        qrcdo(i, k)  = qeso(i, k) + one/hvap*(gamma/(one+gamma))*dbyo(i, k)
        dzd = zid(i, k+1) - zid(i, k)
        dz  = zi(i, k+1) - zi(i, k)
        IF (k .GE. kbcon(i)) THEN
          temd  = xlamde*dzd
          tem   = xlamde*dz
          tem1d = half*xlamdd*dzd
          tem1  = half*xlamdd*dz
        ELSE
          temd  = xlamde*dzd
          tem   = xlamde*dz
          tem1d = half*(xlamd0d(i)*dz+(xlamd(i)+xlamdd)*dzd)
          tem1  = half*(xlamd(i)+xlamdd)*dz
        END IF
        factord= temd - tem1d
        factor = one + tem - tem1
        qcdod(i, k) = (((one-tem1)*qcdod(i, k+1)-tem1d*qcdo(i, k+1)+half*(temd*&
                      (qo(i, k)+qo(i, k+1))+tem*(qod(i, k)+qod(i, k+1))))*&
                      factor-((one-tem1)*qcdo(i, k+1)+tem*half*(qo(i, k)+qo(i, k+1)))*&
                      factord)/factor**2
        qcdo(i, k)  = ((one-tem1)*qcdo(i, k+1)+tem*half*(qo(i, k)+qo(i, k+1)))/factor
        pwdod(i, k) = etadd(i, k+1)*(qcdo(i, k)-qrcdo(i, k)) + etad(i, k+1)*&
                      (qcdod(i, k)-qrcdod(i, k))
        pwdo(i, k)  = etad(i, k+1)*(qcdo(i, k)-qrcdo(i, k))
        qcdod(i, k) = qrcdod(i, k)
        qcdo(i, k)  = qrcdo(i, k)
        pwevod(i)   = pwevod(i) + pwdod(i, k)
        pwevo(i)    = pwevo(i) + pwdo(i, k)
      END IF
    END DO
  END DO


!--- final downdraft strength dependent on precip
!--- efficiency (edt), normalized condensate (pwav), and
!--- evaporate (pwev)
  DO i=1,im
    edtmax = edtmaxl
    IF (slimsk(i) .EQ. zero) edtmax = edtmaxs
    IF (pwevo(i) .LT. zero) THEN
      edtod(i) = -(((edtod(i)*pwavo(i)+edto(i)*pwavod(i))*pwevo(i)-edto(i)*&
                   pwavo(i)*pwevod(i))/pwevo(i)**2)
      edto(i)  = -(edto(i)*pwavo(i)/pwevo(i))
      IF (edto(i) .GT. edtmax) THEN
        edtod(i) = zero
        edto(i)  = edtmax
      ELSE
        edto(i)  = edto(i)
      END IF
    ELSE
      edtod(i) = zero
      edto(i)  = zero
    END IF
  END DO

  DO k=km1,1,-1
    DO i=1,im
      IF (k .LT. jmin(i)) THEN
        gammad = (el2orc*qesod(i, k)*to(i, k)**2-el2orc*qeso(i, k)*2*to(i, k)*&
                  tod(i, k))/(to(i, k)**2)**2
        gamma  = el2orc*qeso(i, k)/to(i, k)**2
        dhhd   = hcdod(i, k)
        dhh    = hcdo(i, k)
        dtd    = tod(i, k)
        dt     = to(i, k)
        dgd    = gammad
        dg     = gamma
        dhd    = hesod(i, k)
        dh     = heso(i, k)
        dzd    = -(zod(i, k+1)-zod(i, k))
        dz     = -(one*(zo(i, k+1)-zo(i, k)))


        aa1d(i) = aa1d(i) + edtod(i)*dz*(g/(cp*dt))*((dhh-dh)/(one+dg))*(one+delta*cp*dg*dt/hvap) &
                          + edto(i)*dzd*(g/(cp*dt))*((dhh-dh)/(one+dg))*(one+delta*cp*dg*dt/hvap) &
                          + edto(i)*dz*(-g*dtd/(cp*dt**2))*((dhh-dh)/(one+dg))*(one+delta*cp*dg*dt/hvap) &
                          + edto(i)*dz*(g/(cp*dt))*((dhhd-dhd)/(one+dg))*(one+delta*cp*dg*dt/hvap) &
                          + edto(i)*dz*(g/(cp*dt))*(-(dhh-dh)*dgd/(one+dg)**2)*(one+delta*cp*dg*dt/hvap) &
                          + edto(i)*dz*(g/(cp*dt))*((dhh-dh)/(one+dg))*(delta*cp*dgd*dt/hvap) &
                          + edto(i)*dz*(g/(cp*dt))*((dhh-dh)/(one+dg))*(delta*cp*dg*dtd/hvap)
       aa1(i) = aa1(i) + edto(i)*dz*(g/(cp*dt))*((dhh-dh)/(one+dg))*&
                          (one+delta*cp*dg*dt/hvap)
        val    = zero
        IF (val .LT. qeso(i, k) - qo(i, k)) THEN
          max2d= qesod(i, k) - qod(i, k)
          max2 = qeso(i, k) - qo(i, k)
        ELSE
          max2 = val
          max2d= zero
        END IF
        aa1d(i)= aa1d(i) + g*delta*(edtod(i)*dz*max2+edto(i)*(dzd*max2+dz*max2d))
        aa1(i) = aa1(i) + edto(i)*dz*g*delta*max2

      END IF
    END DO
  END DO

  DO i=1,im
    cnvscl_6d(i) = half*aa1d(i)*(one-TANH(aa1(i))**2)
    cnvscl_6(i)  = half*(one+TANH(aa1(i)))
    cnvscl_6d(i) = cnvscl_6d(i)*cnvscl_5(i) + cnvscl_6(i)*cnvscl_5d(i)
    cnvscl_6(i)  = cnvscl_6(i)*cnvscl_5(i)
  END DO

!--- what would the change be, that a cloud with unit mass
!--- will do to the environment?

  DO i=1,im
    dp = 1000.0_r_kind*del(i, 1)*ps(i)
    dellahd(i, 1) = g*((edtod(i)*etad(i, 1)+edto(i)*etadd(i, 1))*(hcdo(i, 1)-&
                    heo(i, 1))+edto(i)*etad(i, 1)*(hcdod(i, 1)-heod(i, 1)))/dp
    dellah(i, 1)  = edto(i)*etad(i, 1)*(hcdo(i, 1)-heo(i, 1))*g/dp
    dellaqd(i, 1) = g*((edtod(i)*etad(i, 1)+edto(i)*etadd(i, 1))*(qcdo(i, 1)-&
                    qo(i, 1))+edto(i)*etad(i, 1)*(qcdod(i, 1)-qod(i, 1)))/dp
    dellaq(i, 1)  = edto(i)*etad(i, 1)*(qcdo(i, 1)-qo(i, 1))*g/dp
    dellaud(i, 1) = g*((edtod(i)*etad(i, 1)+edto(i)*etadd(i, 1))*(ucdo(i, 1)-&
                    uo(i, 1))+edto(i)*etad(i, 1)*(ucdod(i, 1)-uod(i, 1)))/dp
    dellau(i, 1)  = edto(i)*etad(i, 1)*(ucdo(i, 1)-uo(i, 1))*g/dp
    dellavd(i, 1) = g*((edtod(i)*etad(i, 1)+edto(i)*etadd(i, 1))*(vcdo(i, 1)-&
                    vo(i, 1))+edto(i)*etad(i, 1)*(vcdod(i, 1)-vod(i, 1)))/dp
    dellav(i, 1)  = edto(i)*etad(i, 1)*(vcdo(i, 1)-vo(i, 1))*g/dp
  END DO

!--- changed due to subsidence and entrainment
  DO k=2,km1
    DO i=1,im
      IF (k .LT. ktcon(i)) THEN
        aup = one
        IF (k .LE. kb(i))   aup = zero
        adw = one
        IF (k .GT. jmin(i)) adw = zero
        dp    = 1000.0_r_kind*del(i, k)*ps(i)
        dzd   = zid(i, k) - zid(i, k-1)
        dz    = zi(i, k) - zi(i, k-1)
        dv1hd = heod(i, k)
        dv1h  = heo(i, k)
        dv2hd = half*(heod(i, k)+heod(i, k-1))
        dv2h  = half*(heo(i, k)+heo(i, k-1))
        dv3hd = heod(i, k-1)
        dv3h  = heo(i, k-1)
        dv1qd = qod(i, k)
        dv1q  = qo(i, k)
        dv2qd = half*(qod(i, k)+qod(i, k-1))
        dv2q  = half*(qo(i, k)+qo(i, k-1))
        dv3qd = qod(i, k-1)
        dv3q  = qo(i, k-1)
        dv1ud = uod(i, k)
        dv1u  = uo(i, k)
        dv2ud = half*(uod(i, k)+uod(i, k-1))
        dv2u  = half*(uo(i, k)+uo(i, k-1))
        dv3ud = uod(i, k-1)
        dv3u  = uo(i, k-1)
        dv1vd = vod(i, k)
        dv1v  = vo(i, k)
        dv2vd = half*(vod(i, k)+vod(i, k-1))
        dv2v  = half*(vo(i, k)+vo(i, k-1))
        dv3vd = vod(i, k-1)
        dv3v  = vo(i, k-1)
        temd  = half*(xlamued(i, k)+xlamued(i, k-1))
        tem   = half*(xlamue(i, k)+xlamue(i, k-1))
        tem1d = xlamudd(i)
        tem1  = xlamud(i)

        IF (k .LE. kbcon(i)) THEN
          ptem   = xlamde
          ptem1d = xlamd0d(i)
          ptem1  = xlamd(i) + xlamdd
        ELSE
          ptem   = xlamde
          ptem1  = xlamdd
          ptem1d = zero
        END IF
        dellahd(i, k) = dellahd(i, k) + g*((aup*eta1d(i, k)-adw*(edtod(i)*etad(i, k)+&
                        edto(i)*etadd(i, k)))*dv1h+(aup*eta(i, k)-adw*edto(i)*&
                        etad(i, k))*dv1hd-(aup*eta1d(i, k-1)-adw*(edtod(i)*&
                        etad(i, k-1)+edto(i)*etadd(i, k-1)))*dv3h-(aup*eta(i, k-1)-&
                        adw*edto(i)*etad(i, k-1))*dv3hd-(aup*(temd*eta(i, k-1)+tem*&
                        eta1d(i, k-1))+adw*ptem*(edtod(i)*etad(i, k)+edto(i)*&
                        etadd(i, k)))*dv2h*dz-(aup*tem*eta(i, k-1)+adw*edto(i)*ptem*&
                        etad(i, k))*(dv2hd*dz+dv2h*dzd)+aup*half*(((tem1d*dz+tem1*dzd)*&
                        eta(i, k-1)+tem1*dz*eta1d(i, k-1))*(hcko(i, k)+hcko(i, k-1))+&
                        tem1*dz*eta(i, k-1)*(hckod(i, k)+hckod(i, k-1)))+adw*half*&
                        (((edtod(i)*etad(i, k)+edto(i)*etadd(i, k))*ptem1*dz+edto(i)*&
                        etad(i, k)*(ptem1d*dz+ptem1*dzd))*(hcdo(i, k)+hcdo(i, k-1))+&
                        edto(i)*etad(i, k)*ptem1*dz*(hcdod(i, k)+hcdod(i, k-1))))/dp
        dellah(i, k)  = dellah(i, k) + ((aup*eta(i, k)-adw*edto(i)*etad(i, k))*&
                        dv1h-(aup*eta(i, k-1)-adw*edto(i)*etad(i, k-1))*dv3h-(&
                        aup*tem*eta(i, k-1)+adw*edto(i)*ptem*etad(i, k))*dv2h*dz+aup*&
                        tem1*eta(i, k-1)*half*(hcko(i, k)+hcko(i, k-1))*dz+adw*edto(i)*&
                        ptem1*etad(i, k)*half*(hcdo(i, k)+hcdo(i, k-1))*dz)*g/dp
        dellaqd(i, k) = dellaqd(i, k) + g*((aup*eta1d(i, k)-adw*(edtod(i)*etad(i, k)+&
                        edto(i)*etadd(i, k)))*dv1q+(aup*eta(i, k)-adw*edto(i)*&
                        etad(i, k))*dv1qd-(aup*eta1d(i, k-1)-adw*(edtod(i)*&
                        etad(i, k-1)+edto(i)*etadd(i, k-1)))*dv3q-(aup*eta(i, k-1)-&
                        adw*edto(i)*etad(i, k-1))*dv3qd-(aup*(temd*eta(i, k-1)+tem*&
                        eta1d(i, k-1))+adw*ptem*(edtod(i)*etad(i, k)+edto(i)*&
                        etadd(i, k)))*dv2q*dz-(aup*tem*eta(i, k-1)+adw*edto(i)*ptem*&
                        etad(i, k))*(dv2qd*dz+dv2q*dzd)+aup*half*(((tem1d*dz+tem1*dzd)*&
                        eta(i, k-1)+tem1*dz*eta1d(i, k-1))*(qcko(i, k)+qcko(i, k-1))+&
                        tem1*dz*eta(i, k-1)*(qckod(i, k)+qckod(i, k-1)))+adw*half*&
                        (((edtod(i)*etad(i, k)+edto(i)*etadd(i, k))*ptem1*dz+edto(i)*&
                        etad(i, k)*(ptem1d*dz+ptem1*dzd))*(qrcdo(i, k)+qrcdo(i, k-1))+&
                        edto(i)*etad(i, k)*ptem1*dz*(qrcdod(i, k)+qrcdod(i, k-1))))/dp
        dellaq(i, k)  = dellaq(i, k) + ((aup*eta(i, k)-adw*edto(i)*etad(i, k))*dv1q-&
                        (aup*eta(i, k-1)-adw*edto(i)*etad(i, k-1))*dv3q-(aup*tem*&
                        eta(i, k-1)+adw*edto(i)*ptem*etad(i, k))*dv2q*dz+aup*tem1*&
                        eta(i, k-1)*half*(qcko(i, k)+qcko(i, k-1))*dz+adw*edto(i)*ptem1*&
                        etad(i, k)*half*(qrcdo(i, k)+qrcdo(i, k-1))*dz)*g/dp
        dellaud(i, k) = dellaud(i, k) + g*((aup*eta1d(i, k)-adw*(edtod(i)*etad(i, k)+&
                        edto(i)*etadd(i, k)))*dv1u+(aup*eta(i, k)-adw*edto(i)*&
                        etad(i, k))*dv1ud-(aup*eta1d(i, k-1)-adw*(edtod(i)*&
                        etad(i, k-1)+edto(i)*etadd(i, k-1)))*dv3u-(aup*eta(i, k-1)-&
                        adw*edto(i)*etad(i, k-1))*dv3ud-(aup*(temd*eta(i, k-1)+tem*&
                        eta1d(i, k-1))+adw*ptem*(edtod(i)*etad(i, k)+edto(i)*&
                        etadd(i, k)))*dv2u*dz-(aup*tem*eta(i, k-1)+adw*edto(i)*ptem*&
                        etad(i, k))*(dv2ud*dz+dv2u*dzd)+aup*half*(((tem1d*dz+tem1*dzd)*&
                        eta(i, k-1)+tem1*dz*eta1d(i, k-1))*(ucko(i, k)+ucko(i, k-1))+&
                        tem1*dz*eta(i, k-1)*(uckod(i, k)+uckod(i, k-1)))+adw*half*&
                        (((edtod(i)*etad(i, k)+edto(i)*etadd(i, k))*ptem1*dz+edto(i)*&
                        etad(i, k)*(ptem1d*dz+ptem1*dzd))*(ucdo(i, k)+ucdo(i, k-1))+&
                        edto(i)*etad(i, k)*ptem1*dz*(ucdod(i, k)+ucdod(i, k-1)))-&
                        pgcon*((aup*eta1d(i, k-1)-adw*(edtod(i)*etad(i, k)+edto(i)*&
                        etadd(i, k)))*(dv1u-dv3u)+(aup*eta(i, k-1)-adw*edto(i)*&
                        etad(i, k))*(dv1ud-dv3ud)))/dp
        dellau(i, k)  = dellau(i, k) + ((aup*eta(i, k)-adw*edto(i)*etad(i, k))*&
                        dv1u-(aup*eta(i, k-1)-adw*edto(i)*etad(i, k-1))*dv3u-(aup*&
                        tem*eta(i, k-1)+adw*edto(i)*ptem*etad(i, k))*dv2u*dz+aup*&
                        tem1*eta(i, k-1)*half*(ucko(i, k)+ucko(i, k-1))*dz+adw*edto(i)*&
                        ptem1*etad(i, k)*half*(ucdo(i, k)+ucdo(i, k-1))*dz-pgcon*(aup*&
                        eta(i, k-1)-adw*edto(i)*etad(i, k))*(dv1u-dv3u))*g/dp
        dellavd(i, k) = dellavd(i, k) + g*((aup*eta1d(i, k)-adw*(edtod(i)*etad(i, k)+&
                        edto(i)*etadd(i, k)))*dv1v+(aup*eta(i, k)-adw*edto(i)*&
                        etad(i, k))*dv1vd-(aup*eta1d(i, k-1)-adw*(edtod(i)*&
                        etad(i, k-1)+edto(i)*etadd(i, k-1)))*dv3v-(aup*eta(i, k-1)-&
                        adw*edto(i)*etad(i, k-1))*dv3vd-(aup*(temd*eta(i, k-1)+&
                        tem*eta1d(i, k-1))+adw*ptem*(edtod(i)*etad(i, k)+edto(i)*&
                        etadd(i, k)))*dv2v*dz-(aup*tem*eta(i, k-1)+adw*edto(i)*&
                        ptem*etad(i, k))*(dv2vd*dz+dv2v*dzd)+aup*half*(((tem1d*dz+&
                        tem1*dzd)*eta(i, k-1)+tem1*dz*eta1d(i, k-1))*(vcko(i, k)+&
                        vcko(i, k-1))+tem1*dz*eta(i, k-1)*(vckod(i, k)+&
                        vckod(i, k-1)))+adw*half*(((edtod(i)*etad(i, k)+edto(i)*&
                        etadd(i, k))*ptem1*dz+edto(i)*etad(i, k)*(ptem1d*dz+ptem1*&
                        dzd))*(vcdo(i, k)+vcdo(i, k-1))+edto(i)*etad(i, k)*ptem1*dz*&
                        (vcdod(i, k)+vcdod(i, k-1)))-pgcon*((aup*eta1d(i, k-1)-adw*&
                        (edtod(i)*etad(i, k)+edto(i)*etadd(i, k)))*(dv1v-dv3v)+(aup*&
                        eta(i, k-1)-adw*edto(i)*etad(i, k))*(dv1vd-dv3vd)))/dp
        dellav(i, k)  = dellav(i, k) + ((aup*eta(i, k)-adw*edto(i)*etad(i, k))*dv1v-&
                        (aup*eta(i, k-1)-adw*edto(i)*etad(i, k-1))*dv3v-(aup*tem*&
                        eta(i, k-1)+adw*edto(i)*ptem*etad(i, k))*dv2v*dz+aup*tem1*&
                        eta(i, k-1)*half*(vcko(i, k)+vcko(i, k-1))*dz+adw*edto(i)*ptem1*&
                        etad(i, k)*half*(vcdo(i, k)+vcdo(i, k-1))*dz-pgcon*(aup*&
                        eta(i, k-1)-adw*edto(i)*etad(i, k))*(dv1v-dv3v))*g/dp
      END IF
    END DO
  END DO

!------- cloud top

  DO i=1,im
    indx = ktcon(i)
    dp   = 1000.*del(i, indx)*ps(i)
    dv1hd= heod(i, indx-1)
    dv1h = heo(i, indx-1)
    dellahd(i, indx) = g*(eta1d(i, indx-1)*(hcko(i, indx-1)-dv1h)+eta(i, indx-1)*&
                      (hckod(i, indx-1)-dv1hd))/dp
    dellah(i, indx)  = eta(i, indx-1)*(hcko(i, indx-1)-dv1h)*g/dp

    dv1qd = qod(i, indx-1)
    dv1q  = qo(i, indx-1)
    dellaqd(i, indx) = g*(eta1d(i, indx-1)*(qcko(i, indx-1)-dv1q)+eta(i, indx-1)*&
                       (qckod(i, indx-1)-dv1qd))/dp
    dellaq(i, indx)  = eta(i, indx-1)*(qcko(i, indx-1)-dv1q)*g/dp

    dv1ud = uod(i, indx-1)
    dv1u  = uo(i, indx-1)
    dellaud(i, indx) = g*(eta1d(i, indx-1)*(ucko(i, indx-1)-dv1u)+eta(i, indx-1)*&
                       (uckod(i, indx-1)-dv1ud))/dp
    dellau(i, indx)  = eta(i, indx-1)*(ucko(i, indx-1)-dv1u)*g/dp

    dv1vd = vod(i, indx-1)
    dv1v  = vo(i, indx-1)
    dellavd(i, indx) = g*(eta1d(i, indx-1)*(vcko(i, indx-1)-dv1v)+eta(i, indx-1)*&
                       (vckod(i, indx-1)-dv1vd))/dp
    dellav(i, indx)  = eta(i, indx-1)*(vcko(i, indx-1)-dv1v)*g/dp


    dellald(i, indx) = g*(eta1d(i,indx-1)*qlko_ktcon(i)+eta(i,indx-1)*&
                       qlko_ktcond(i))/dp
    dellal(i, indx)  = eta(i, indx-1)*qlko_ktcon(i)*g/dp
  END DO


!------- final changed variable per unit mass flux
  DO k=1,km
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        IF (k .GT. ktcon(i)) THEN
          qod(i, k) = q1d(i, k)
          qo(i, k)  = q1(i, k)
          tod(i, k) = t1d(i, k)
          to(i, k)  = t1(i, k)
        END IF
        IF (k .LE. ktcon(i)) THEN
          qod(i, k) = mbdt*dellaqd(i, k) + q1d(i, k)
          qo(i, k)  = dellaq(i, k)*mbdt + q1(i, k)
          dellatd   = (dellahd(i, k)-hvap*dellaqd(i, k))/cp
          dellat    = (dellah(i, k)-hvap*dellaq(i, k))/cp
          tod(i, k) = mbdt*dellatd + t1d(i, k)
          to(i, k)  = dellat*mbdt + t1(i, k)
          val       = 1.e-10_r_kind
          IF (qo(i, k) .LT. val) THEN
            qod(i, k) = zero
            qo(i, k)  = val
          ELSE
            qo(i, k)  = qo(i, k)
          END IF
        END IF
      END IF
    END DO
  END DO
!--- the above changed environment is now used to calulate the
!--- effect the arbitrary cloud (with unit mass flux)
!--- would have on the stability,
!--- which then is used to calculate the real mass flux,
!--- necessary to keep this change in balance with the large-scale
!--- destabilization.

!--- environmental conditions again, first heights
  DO k=1,km
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        call fpvsx_tl(to(i, k), es0, tod(i, k), es0d)
        qesod(i, k) = 10.0_r_kind*es0d
        qeso(i, k)  = es0*10.0_r_kind
        qesod(i, k) = (eps*qesod(i, k)*(pfld(i, k)+epsm1*qeso(i, k))-eps*qeso(i, k)*&
                      epsm1*qesod(i, k))/(pfld(i, k)+epsm1*qeso(i, k))**2
        qeso(i, k)  = eps*qeso(i, k)/(pfld(i, k)+epsm1*qeso(i, k))
        val = 1.e-8_r_kind
        IF (qeso(i, k) .LT. val) THEN
          qesod(i, k) = zero
          qeso(i, k)  = val
        END IF
      END IF
    END DO
  END DO

!--- moist static energy

  DO k=1,km1
    DO i=1,im
      IF (k .LE. kmax(i) - 1) THEN
        dzd = half*(zod(i, k+1)-zod(i, k))
        dz  = half*(zo(i, k+1)-zo(i, k))
        dp  = half*(pfld(i, k+1)-pfld(i, k))
        call fpvsx_tl(to(i, k+1), es0, tod(i, k+1), es0d)
        esd = 10.0_r_kind*es0d
        es  = es0*10.0_r_kind
        pprimed = epsm1*esd
        pprime  = pfld(i, k+1) + epsm1*es
        qsd = (eps*esd*pprime-eps*es*pprimed)/pprime**2
        qs  = eps*es/pprime
        dqsdpd = -((qsd*pprime-qs*pprimed)/pprime**2)
        dqsdp  = -(qs/pprime)
        desdtd = esd*(fact1/to(i, k+1)+fact2/to(i, k+1)**2) + es*(-(fact1*&
                 tod(i, k+1)/to(i, k+1)**2)-fact2*2*to(i, k+1)*tod(i, k+1)/&
                 (to(i, k+1)**2)**2)
        desdt  = es*(fact1/to(i, k+1)+fact2/to(i, k+1)**2)
        dqsdtd = (pfld(i, k+1)*(qsd*desdt+qs*desdtd)*es*pprime-qs*pfld(i, k+1)*&
                 desdt*(esd*pprime+es*pprimed))/(es*pprime)**2
        dqsdt  = qs*pfld(i, k+1)*desdt/(es*pprime)
        gammad = (el2orc*qesod(i, k+1)*to(i, k+1)**2-el2orc*qeso(i, k+1)*2*&
                 to(i, k+1)*tod(i, k+1))/(to(i, k+1)**2)**2
        gamma  = el2orc*qeso(i, k+1)/to(i, k+1)**2
        dtd    = ((g*dzd+hvap*dp*dqsdpd)*cp*(one+gamma)-(g*dz+hvap*dqsdp*dp)*cp*&
                 gammad)/(cp*(one+gamma))**2
        dt  = (g*dz+hvap*dqsdp*dp)/(cp*(one+gamma))
        dqd = dqsdtd*dt + dqsdt*dtd + dp*dqsdpd
        dq  = dqsdt*dt + dqsdp*dp
        tod(i, k) = tod(i, k+1) + dtd
        to(i, k)  = to(i, k+1)  + dt
        qod(i, k) = qod(i, k+1) + dqd
        qo(i, k)  = qo(i, k+1)  + dq
        pod(i, k) = zero
        po(i, k)  = half*(pfld(i, k)+pfld(i, k+1))
      END IF
    END DO
  END DO

  DO k=1,km1
    DO i=1,im
      IF (k .LE. kmax(i) - 1) THEN
        call fpvsx_tl(to(i, k), es0, tod(i, k),  es0d)
        qesod(i, k) = 10.0_r_kind*es0d
        qeso(i, k)  = es0*10.0_r_kind
        qesod(i, k) = (eps*qesod(i, k)*(po(i, k)+epsm1*qeso(i, k))-eps*qeso(i, k)*&
                       epsm1*qesod(i, k))/(po(i, k)+epsm1*qeso(i, k))**2
        qeso(i, k)  = eps*qeso(i, k)/(po(i, k)+epsm1*qeso(i, k))
        val1 = 1.e-8-r_kind
        IF (qeso(i, k) .LT. val1) THEN
          qesod(i, k) = zero
          qeso(i, k)  = val1
        ELSE
          qeso(i, k)  = qeso(i, k)
        END IF
        val2 = 1.e-10_r_kind
        IF (qo(i, k) .LT. val2) THEN
          qod(i, k) = zero
          qo(i, k)  = val2
        ELSE
          qo(i, k) = qo(i, k)
        END IF
        heod(i, k) = half*g*(zod(i, k)+zod(i, k+1)) + cp*tod(i, k) + hvap*qod(i, k)
        heo(i, k)  = half*g*(zo(i, k)+zo(i, k+1)) + cp*to(i, k) + hvap*qo(i, k)
        hesod(i, k)= half*g*(zod(i, k)+zod(i, k+1)) + cp*tod(i, k) + hvap*qesod(i, k)
        heso(i, k) = half*g*(zo(i, k)+zo(i, k+1)) + cp*to(i, k) + hvap*qeso(i, k)
      END IF
    END DO
  END DO

  DO i=1,im
    k = kmax(i)
    heod(i, k) = g*zod(i, k) + hvap*qod(i, k)+ cp*tod(i,k)
    heo(i, k)  = g*zo(i, k) + cp*to(i, k) + hvap*qo(i, k)
    hesod(i, k)= g*zod(i, k) + hvap*qesod(i, k)+cp*tod(i,k)
    heso(i, k) = g*zo(i, k) + cp*to(i, k) + hvap*qeso(i, k)
  END DO


!**************************** static control
!------- moisture and cloud work functions
  DO i=1,im
    indx = kb(i)
    hckod(i, indx) = heod(i, indx)
    hcko(i, indx)  = heo(i, indx)
    qckod(i, indx) = qod(i, indx)
    qcko(i, indx)  = qo(i, indx)
  END DO

  DO k=2,km1
    DO i=1,im
      IF (k .GT. kb(i) .AND. k .LE. ktcon(i)) THEN
        dzd  = zid(i, k) - zid(i, k-1)
        dz   = zi(i, k) - zi(i, k-1)
        temd = half*((xlamued(i,k)+xlamued(i,k-1))*dz+(xlamue(i,k)+xlamue(i,k-1))*dzd)
        tem  = half*(xlamue(i, k)+xlamue(i, k-1))*dz
        tem1d= half*(xlamudd(i)*dz+xlamud(i)*dzd)
        tem1 = half*xlamud(i)*dz
        factord = temd - tem1d
        factor  = one + tem - tem1
        hckod(i, k) = (((one-tem1)*hckod(i, k-1)-tem1d*hcko(i, k-1)+half*(temd*&
                      (heo(i, k)+heo(i, k-1))+tem*(heod(i, k)+heod(i, k-1))))*&
                      factor-((one-tem1)*hcko(i, k-1)+tem*half*(heo(i, k)+heo(i, k-1)))*&
                      factord)/factor**2
        hcko(i, k)  = ((one-tem1)*hcko(i, k-1)+tem*half*(heo(i, k)+heo(i, k-1)))/factor
      END IF
    END DO
  END DO

  DO k=2,km1
    DO i=1,im
      IF (k .GT. kb(i) .AND. k .LT. ktcon(i)) THEN
        dzd = zid(i, k) - zid(i, k-1)
        dz  = zi(i, k) - zi(i, k-1)
        gammad = (el2orc*qesod(i, k)*to(i, k)**2-el2orc*qeso(i, k)*2*to(i, k)*&
                 tod(i, k))/(to(i, k)**2)**2
        gamma  = el2orc*qeso(i, k)/to(i, k)**2
        xdbyd  = hckod(i, k) - hesod(i, k)
        xdby   = hcko(i, k) - heso(i, k)
        xqrchd = qesod(i, k) + ((gammad*xdby+gamma*xdbyd)*hvap*(one+gamma)-gamma*&
                 xdby*hvap*gammad)/(hvap*(one+gamma))**2
        xqrch  = qeso(i, k) + gamma*xdby/(hvap*(one+gamma))
        temd   = half*((xlamued(i,k)+xlamued(i,k-1))*dz+(xlamue(i,k)+xlamue(i,k-1))*dzd)
        tem    = half*(xlamue(i, k)+xlamue(i, k-1))*dz
        tem1d  = half*(xlamudd(i)*dz+xlamud(i)*dzd)
        tem1   = half*xlamud(i)*dz
        factord= temd - tem1d
        factor = one + tem - tem1
        qckod(i, k) = (((one-tem1)*qckod(i, k-1)-tem1d*qcko(i, k-1)+half*(temd*&
                      (qo(i, k)+qo(i, k-1))+tem*(qod(i, k)+qod(i, k-1))))*&
                      factor-((one-tem1)*qcko(i, k-1)+tem*half*(qo(i, k)+qo(i, k-1)))*&
                      factord)/factor**2
        qcko(i, k)  = ((one-tem1)*qcko(i, k-1)+tem*half*(qo(i, k)+qo(i, k-1)))/factor
        dqd = eta1d(i, k)*(qcko(i, k)-xqrch) + eta(i, k)*(qckod(i, k)-xqrchd)
        dq  = eta(i, k)*(qcko(i, k)-xqrch)
        IF (k .GE. kbcon(i) .AND. dq .GT. zero) THEN
          etahd = half*(eta1d(i, k)+eta1d(i, k-1))
          etah  = half*(eta(i, k)+eta(i, k-1))
          IF (ncloud .GT. zero .AND. k .GT. jmin(i)) THEN
            qlkd = (dqd*(eta(i, k)+etah*(c0+c1)*dz)-dq*(eta1d(i, k)+(c0+c1)*&
                   (etahd*dz+etah*dzd)))/(eta(i, k)+etah*(c0+c1)*dz)**2
            qlk  = dq/(eta(i, k)+etah*(c0+c1)*dz)
          ELSE
            qlkd = (dqd*(eta(i, k)+etah*c0*dz)-dq*(eta1d(i, k)+c0*(etahd*dz+etah*&
                   dzd)))/(eta(i, k)+etah*c0*dz)**2
            qlk  = dq/(eta(i, k)+etah*c0*dz)
          END IF

          IF (k .LT. ktcon1(i)) THEN
            xaa0d(i) = xaa0d(i) - g*(dzd*qlk+dz*qlkd)
            xaa0(i)  = xaa0(i) - dz*g*qlk
          END IF
          qckod(i, k)= qlkd + xqrchd
          qcko(i, k) = qlk + xqrch
          xpwd0      = c0*((etahd*dz+etah*dzd)*qlk+etah*dz*qlkd)
          xpw        = etah*c0*dz*qlk
          xpwavd(i)  = xpwavd(i) + xpwd0
          xpwav(i)   = xpwav(i) + xpw
        END IF
      END IF
      IF (k .GE. kbcon(i) .AND. k .LT. ktcon1(i)) THEN
        dz1d = zod(i, k+1) - zod(i, k)
        dz1  = zo(i, k+1) - zo(i, k)
        gammad = (el2orc*qesod(i, k)*to(i, k)**2-el2orc*qeso(i, k)*2*to(i, k)*&
                 tod(i, k))/(to(i, k)**2)**2
        gamma  = el2orc*qeso(i, k)/to(i, k)**2
        rfactd = delta*cp*(gammad*to(i, k)+gamma*tod(i, k))/hvap
        rfact  = one + delta*cp*gamma*to(i, k)/hvap
        xaa0d(i) = xaa0d(i) + (((dz1d*xdby+dz1*xdbyd)*g/(cp*to(i, k))-dz1*xdby*g*&
                   tod(i, k)/(cp*to(i, k)**2))*(one+gamma)-dz1*g*xdby*gammad/(cp*&
                   to(i, k)))*rfact/(one+gamma)**2 + dz1*g*xdby*rfactd/(cp*&
                   to(i, k)*(one+gamma))
        xaa0(i)  = xaa0(i) + dz1*(g/(cp*to(i, k)))*xdby/(one+gamma)*rfact
        val = zero
        IF (val .LT. qeso(i, k) - qo(i, k)) THEN
          max3d = qesod(i, k) - qod(i, k)
          max3  = qeso(i, k) - qo(i, k)
        ELSE
          max3d = zero
          max3 = val
        END IF
        xaa0d(i) = xaa0d(i) + g*delta*(dz1d*max3+dz1*max3d)
        xaa0(i)  = xaa0(i) + dz1*g*delta*max3
      END IF
    END DO
  END DO

!------- downdraft calculations
!--- downdraft moisture properties

  DO i=1,im
    jmn = jmin(i)
    hcdod(i, jmn) = heod(i, jmn)
    hcdo(i, jmn)  = heo(i, jmn)
    qcdod(i, jmn) = qod(i, jmn)
    qcdo(i, jmn)  = qo(i, jmn)
    qrcdd(i, jmn) = qesod(i, jmn)
    qrcd(i, jmn)  = qeso(i, jmn)
  END DO

  DO k=km1,1,-1
    DO i=1,im
      IF (k .LT. jmin(i)) THEN
        dzd = zid(i, k+1) - zid(i, k)
        dz  = zi(i, k+1)  - zi(i, k)
        IF (k .GE. kbcon(i)) THEN
          temd  = xlamde*dzd
          tem   = xlamde*dz
          tem1d = half*xlamdd*dzd
          tem1  = half*xlamdd*dz
        ELSE
          temd = xlamde*dzd
          tem  = xlamde*dz
          tem1d= half*(xlamd0d(i)*dz+(xlamd(i)+xlamdd)*dzd)
          tem1 = half*(xlamd(i)+xlamdd)*dz
        END IF
        factord = temd - tem1d
        factor  = one + tem - tem1
        hcdod(i, k)= (((one-tem1)*hcdod(i, k+1)-tem1d*hcdo(i, k+1)+half*(temd*&
                     (heo(i, k)+heo(i, k+1))+tem*(heod(i, k)+heod(i, k+1))))*&
                     factor-((one-tem1)*hcdo(i, k+1)+tem*half*(heo(i, k)+heo(i, k+1)))*&
                     factord)/factor**2
        hcdo(i, k) = ((one-tem1)*hcdo(i, k+1)+tem*half*(heo(i, k)+heo(i, k+1)))/factor
      END IF
    END DO
  END DO

  DO k=km1,1,-1
    DO i=1,im
      IF (k .LT. jmin(i)) THEN
        dqd = qesod(i, k)
        dq  = qeso(i, k)
        dtd = tod(i, k)
        dt  = to(i, k)
        gammad = (el2orc*dqd*dt**2-el2orc*dq*2*dt*dtd)/(dt**2)**2
        gamma  = el2orc*dq/dt**2
        dhd = hcdod(i, k) - hesod(i, k)
        dh  = hcdo(i, k) - heso(i, k)
        qrcdd(i, k) = dqd + ((gammad*(one+gamma)-gamma*gammad)*dh/(one+gamma)**2+&
                      gamma*dhd/(one+gamma))/hvap
        qrcd(i, k)  = dq + one/hvap*(gamma/(one+gamma))*dh
        dzd = zid(i, k+1) - zid(i, k)
        dz  = zi(i, k+1) - zi(i, k)
        IF (k .GE. kbcon(i)) THEN
          temd = xlamde*dzd
          tem  = xlamde*dz
          tem1d= half*xlamdd*dzd
          tem1 = half*xlamdd*dz
        ELSE
          temd = xlamde*dzd
          tem  = xlamde*dz
          tem1d= half*(xlamd0d(i)*dz+(xlamd(i)+xlamdd)*dzd)
          tem1 = half*(xlamd(i)+xlamdd)*dz
        END IF
        factord = temd - tem1d
        factor  = one + tem - tem1
        qcdod(i, k) = (((one-tem1)*qcdod(i, k+1)-tem1d*qcdo(i, k+1)+half*(temd*(qo(i, k)&
                      +qo(i, k+1))+tem*(qod(i, k)+qod(i, k+1))))*factor-((one-tem1)*&
                      qcdo(i, k+1)+tem*half*(qo(i, k)+qo(i, k+1)))*factord)/factor**2
        qcdo(i, k)  = ((one-tem1)*qcdo(i, k+1)+tem*half*(qo(i, k)+qo(i, k+1)))/factor
        xpwdd       = etadd(i, k+1)*(qcdo(i, k)-qrcd(i, k)) + etad(i, k+1)*&
                      (qcdod(i, k)-qrcdd(i, k))
        xpwd        = etad(i, k+1)*(qcdo(i, k)-qrcd(i, k))
        qcdod(i, k) = qrcdd(i, k)
        qcdo(i, k)  = qrcd(i, k)
        xpwevd(i)   = xpwevd(i) + xpwdd
        xpwev(i)    = xpwev(i) + xpwd
      END IF
    END DO
  END DO

  DO i=1,im
    edtmax = edtmaxl
    IF (slimsk(i) .EQ. zero) edtmax = edtmaxs
    IF (xpwev(i) .GE. zero) THEN
      edtxd(i) = zero
      edtx(i)  = zero
    ELSE
      edtxd(i) = -(((edtxd(i)*xpwav(i)+edtx(i)*xpwavd(i))*xpwev(i)-edtx(i)*&
                  xpwav(i)*xpwevd(i))/xpwev(i)**2)
      edtx(i)  = -(edtx(i)*xpwav(i)/xpwev(i))
      IF (edtx(i) .GT. edtmax) THEN
        edtxd(i)= zero
        edtx(i) = edtmax
      ELSE
        edtx(i) = edtx(i)
      END IF
    END IF
  END DO

!--- downdraft cloudwork functions
  DO k=km1,1,-1
    DO i=1,im
      IF (k .LT. jmin(i)) THEN
        gammad= (el2orc*qesod(i, k)*to(i, k)**2-el2orc* &
                qeso(i, k)*2*to(i, k)*tod(i, k))/(to(i, k)**2)**2
        gamma = el2orc*qeso(i, k)/to(i, k)**2
        dhhd  = hcdod(i, k)
        dhh = hcdo(i, k)
        dtd = tod(i, k)
        dt  = to(i, k)
        dgd = gammad
        dg  = gamma
        dhd = hesod(i, k)
        dh  = heso(i, k)
        dzd = -(zod(i, k+1)-zod(i, k))
        dz  = -(one*(zo(i, k+1)-zo(i, k)))
        xaa0d(i) = xaa0d(i) + (((edtxd(i)*dz+edtx(i)*dzd)*(dhh-dh)/(one+dg)+&
                   edtx(i)*dz*((dhhd-dhd)*(one+dg)-(dhh-dh)*dgd)/(one+dg)**2)*g/&
                   (cp*dt)-edtx(i)*dz*(dhh-dh)*g*dtd/((one+dg)*cp*dt**2))*(one+delta*&
                   cp*dg*dt/hvap) + edtx(i)*dz*(dhh-dh)*g*delta*(dgd*dt+dg*dtd)/&
                   ((one+dg)*dt*hvap)
        xaa0(i)  = xaa0(i) + edtx(i)*dz*(g/(cp*dt))*((dhh-dh)/(one+dg))*&
                   (one+delta*cp*dg*dt/hvap)
        val = zero
        IF (val .LT. qeso(i, k) - qo(i, k)) THEN
          max4d = qesod(i, k) - qod(i, k)
          max4  = qeso(i, k) - qo(i, k)
        ELSE
          max4  = val
          max4d = zero
        END IF
        xaa0d(i)= xaa0d(i) + g*delta*(edtxd(i)*dz*max4+edtx(i)*(dzd*max4+dz*max4d))
        xaa0(i) = xaa0(i) + edtx(i)*dz*g*delta*max4
      END IF
    END DO
  END DO

!  calculate critical cloud work function

  DO i=1,im
    IF (pfld(i, ktcon(i)) .LT. pcrit(15)) THEN
      acrt(i) = acrit(15)*(975.0_r_kind-pfld(i, ktcon(i)))/(975.0_r_kind-pcrit(15))
    ELSE IF (pfld(i, ktcon(i)) .GT. pcrit(1)) THEN
      acrt(i) = acrit(1)
    ELSE
      k = INT((850.0_r_kind-pfld(i, ktcon(i)))/50.0_r_kind) + 2
      IF (k .GT. 15) THEN
        k = 15
      ELSE
        k = k
      END IF
      IF (k .LT. 2) THEN
        k = 2
      ELSE
        k = k
      END IF
      acrt(i)  = acrit(k) + (acrit(k-1)-acrit(k))*(pfld(i, ktcon(i))-pcrit(k))/&
                            (pcrit(k-1)-pcrit(k))
    END IF
  END DO


  DO i=1,im
    IF (slimsk(i) .EQ. one) THEN
      w1 = w1l
      w2 = w2l
      w3 = w3l
      w4 = w4l
    ELSE
      w1 = w1s
      w2 = w2s
      w3 = w3s
      w4 = w4s
    END IF

!  modify critical cloud workfunction by cloud base vertical velocity
    IF (pdot(i) .LE. w4) THEN
      acrtfct(i) = (pdot(i)-w4)/(w3-w4)
      acrtfctd(i) = pdotd(i)/(w3-w4)
    ELSE IF (pdot(i) .GE. -w4) THEN
      acrtfct(i) = -((pdot(i)+w4)/(w4-w3))
      acrtfctd(i) = -pdotd(i)/(w4-w3)
    ELSE
      acrtfct(i) = zero
      acrtfctd(i) = zero
    END IF
    val1 = -one
    IF (acrtfct(i) .LT. val1) THEN
      acrtfct(i) = val1
      acrtfctd(i) = zero
    ELSE
      acrtfct(i) = acrtfct(i)
    END IF

    val2 = one
    IF (acrtfct(i) .GT. val2) THEN
      acrtfct(i) = val2
      acrtfctd(i) = zero
    ELSE
     acrtfct(i) = acrtfct(i)
    END IF

    acrtfct(i) = one - acrtfct(i)
    acrtfctd(i) = -acrtfctd(i)

    IF (1800.0_r_kind - dt2 .LT. zero) THEN
      max5 = zero
    ELSE
      max5 = 1800.0_r_kind - dt2
    END IF


!  modify adjustment time scale by cloud base vertical velocity
    dtconv(i) = dt2 + max5*(pdot(i)-w2)/(w1-w2)
    dtconvd(i) = max5*pdotd(i)/(w1-w2)

    IF (dtconv(i) .LT. dtmin) THEN
      dtconv(i) = dtmin
      dtconvd(i) = zero
    ELSE
      dtconv(i) = dtconv(i)
    END IF
    IF (dtconv(i) .GT. dtmax) THEN
      dtconv(i) = dtmax
      dtconvd(i) = zero
    ELSE
      dtconv(i) = dtconv(i)
    END IF
 ENDDO

!--- large scale forcing

  DO i=1,im
    fldd(i) = ((aa1d(i)-acrt(i)*acrtfctd(i))*dtconv(i)- &
               dtconvd(i)*(aa1(i)-acrt(i)*acrtfct(i)))/(dtconv(i)**2)
    fld(i)  = (aa1(i)-acrt(i)*acrtfct(i))/dtconv(i)
    if(fld(i) .le. zero) then
          cnvscl_7(i) = tiny_r_kind
          cnvscl_7d(i) = zero
    endif

    cnvscl_7d(i) = cnvscl_7d(i)*cnvscl_6(i) + cnvscl_7(i)*cnvscl_6d(i)
    cnvscl_7(i)  = cnvscl_7(i)*cnvscl_6(i)

    xkd(i)  = (xaa0d(i)-aa1d(i))/mbdt
    xk(i)   = (xaa0(i)-aa1(i))/mbdt

    if(xk(i) .ge. zero) then
        cnvscl_8(i) = tiny_r_kind
        cnvscl_8d(i) =zero
    endif

    cnvscl_8d(i) = cnvscl_8d(i)*cnvscl_7(i) + cnvscl_8(i)*cnvscl_7d(i)
    cnvscl_8(i)  = cnvscl_8(i)*cnvscl_7(i)


!--- kernel, cloud base mass flux

    xmbd(i) = -((fldd(i)*xk(i)-fld(i)*xkd(i))/xk(i)**2)
    xmb(i)  = -(fld(i)/xk(i))
    IF (xmb(i) .GT. xmbmax(i)) THEN
      xmbd(i) = zero
      xmb(i)  = xmbmax(i)
    ELSE
      xmb(i)  = xmb(i)
    END IF

    xmbd(i) = xmbd(i)*cnvscl_8(i) + xmb(i)*cnvscl_8d(i)
    xmb(i)  = xmb(i)*cnvscl_8(i)

  END DO

!  restore to,qo,uo,vo to t1,q1,u1,v1 in case convection stops

  DO k=1,km
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        tod(i, k) = t1d(i, k)
        to(i, k)  = t1(i, k)
        qod(i, k) = q1d(i, k)
        qo(i, k)  = q1(i, k)
        uod(i, k) = u1d(i, k)
        uo(i, k)  = u1(i, k)
        vod(i, k) = v1d(i, k)
        vo(i, k)  = v1(i, k)
        call fpvsx_tl(t1(i, k),es0, t1d(i, k), es0d)
        esd = 10.0_r_kind*es0d
        es  = es0*10.0_r_kind
        qesod(i, k) = (eps*esd*(pfld(i, k)+epsm1*es)-eps*es*&
                      epsm1*esd)/(pfld(i, k)+epsm1*es)**2
        qeso(i, k)  = eps*es/(pfld(i, k)+epsm1*es)
        val = 1.e-8_r_kind
        IF (qeso(i, k) .LT. val) THEN
          qesod(i, k) = zero
          qeso(i, k)  = val
        ELSE
          qeso(i, k)  = qeso(i, k)
        END IF
      END IF
    END DO
  END DO

!--- feedback: simply the changes from the cloud with unit mass flux
!---           multiplied by  the mass flux necessary to keep the
!---           equilibrium with the larger-scale.

  DO k=1,km
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        IF (k .LE. ktcon(i)) THEN
          dellatd   = (dellahd(i, k)-hvap*dellaqd(i, k))/cp
          dellat    = (dellah(i, k)-hvap*dellaq(i, k))/cp
          t1d(i, k) = t1d(i, k) + dt2*(dellatd*xmb(i)+dellat*xmbd(i))
          t1(i, k)  = t1(i, k) + dellat*xmb(i)*dt2

          q1d(i, k) = q1d(i, k) + dt2*(dellaqd(i, k)*xmb(i)+dellaq(i, k)*xmbd(i))
          q1(i, k)  = q1(i, k) + dellaq(i, k)*xmb(i)*dt2
          tem = one/rcs(i)
          u1d(i, k) = u1d(i, k) + dt2*tem*(dellaud(i, k)*xmb(i)+dellau(i, k)*xmbd(i))
          u1(i, k)  = u1(i, k) + dellau(i, k)*xmb(i)*dt2*tem
          v1d(i, k) = v1d(i, k) + dt2*tem*(dellavd(i, k)*xmb(i)+dellav(i, k)*xmbd(i))
          v1(i, k)  = v1(i, k) + dellav(i, k)*xmb(i)*dt2*tem
          dp = 1000.0_r_kind*del(i, k)*ps(i)
          delhbar(i)= delhbar(i) + dellah(i, k)*xmb(i)*dp/g
          delqbar(i)= delqbar(i) + dellaq(i, k)*xmb(i)*dp/g
          deltbar(i)= deltbar(i) + dellat*xmb(i)*dp/g
          delubar(i)= delubar(i) + dellau(i, k)*xmb(i)*dp/g
          delvbar(i)= delvbar(i) + dellav(i, k)*xmb(i)*dp/g
        END IF
      END IF
    END DO
  END DO

  DO k=1,km
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        IF (k .LE. ktcon(i)) THEN
          call fpvsx_tl(t1(i, k), es0, t1d(i, k), es0d)
          qesod(i, k) = 10.0_r_kind*es0d
          qeso(i, k)  = es0*10.0_r_kind
          qesod(i, k) = (eps*qesod(i, k)*(pfld(i, k)+epsm1*qeso(i, k))-eps*&
                         qeso(i, k)*epsm1*qesod(i, k))/(pfld(i, k)+epsm1*qeso(i, k))**2
          qeso(i, k)  = eps*qeso(i, k)/(pfld(i, k)+epsm1*qeso(i, k))
          val = 1.e-8_r_kind
          IF (qeso(i, k) .LT. val) THEN
            qesod(i, k) = zero
            qeso(i, k)  = val
          ELSE
            qeso(i, k)  = qeso(i, k)
          END IF
        END IF
      END IF
    END DO
  END DO

  DO k=km,1,-1
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        IF (k .LT. ktcon(i)) THEN
          aup = one
          IF (k .LE. kb(i)) aup = zero
          adw = one
          IF (k .GE. jmin(i)) adw = zero
          raind = aup*pwod(i, k) + adw*(edtod(i)*pwdo(i, k)+edto(i)*pwdod(i, k))
          rain  = aup*pwo(i, k) + adw*edto(i)*pwdo(i, k)
          rntotd(i) = rntotd(i) + .001_r_kind*dt2*(raind*xmb(i)+rain*xmbd(i))
          rntot(i)  = rntot(i) + rain*xmb(i)*.001_r_kind*dt2
        END IF
      END IF
    END DO
  END DO
  DO k=km,1,-1
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        deltv(i) = zero
        delq(i)  = zero
        qevapd(i)= zero
        qevap(i) = zero
        dp = 1000.0_r_kind*del(i, k)*ps(i)

        IF (k .LT. ktcon(i)) THEN
          aup = one
          IF (k .LE. kb(i)) aup = zero
          adw = one
          IF (k .GE. jmin(i)) adw = zero
          raind = aup*pwod(i, k) + adw*(edtod(i)*pwdo(i, k)+edto(i)*pwdod(i, k))
          rain  = aup*pwo(i, k) + adw*edto(i)*pwdo(i, k)
          rnd(i)= rnd(i) + .001_r_kind*dt2*(raind*xmb(i)+rain*xmbd(i))
          rn(i) = rn(i) + rain*xmb(i)*.001_r_kind*dt2
        END IF

        IF (k .LT. ktcon(i)) THEN
          evefd = evfact*edtd(i)
          evef  = edt(i)*evfact
          IF (slimsk(i) .EQ. 1.) THEN
            evefd= evfactl*edtd(i)
            evef = edt(i)*evfactl
          END IF

          qcondd(i) = ((evefd*(q1(i, k)-qeso(i, k))+evef*(q1d(i, k)-&
                      qesod(i, k)))*(one+el2orc*qeso(i, k)/t1(i, k)**2)-evef*(&
                      q1(i, k)-qeso(i, k))*(el2orc*qesod(i, k)*t1(i, k)**2-el2orc*&
                      qeso(i, k)*2*t1(i, k)*t1d(i, k))/t1(i, k)**4)/(one+el2orc*&
                      qeso(i, k)/t1(i, k)**2)**2
          qcond(i)  = evef*(q1(i, k)-qeso(i, k))/(one+el2orc*qeso(i, k)/t1(i, k)**2)

          IF (rn(i) .GT. zero .AND. qcond(i) .LT. zero) THEN
            qevap(i)  = -(qcond(i)*(ONE-EXP(-(.32_r_kind*SQRT(dt2*rn(i))))))
            qevapd(i) = -qcondd(i)*(ONE-EXP(-(.32_r_kind*SQRT(dt2*rn(i))))) &
           -0.32_r_kind*qcond(i)*EXP(-(.32_r_kind*SQRT(dt2*rn(i))))*half*dt2*rnd(i)/SQRT(dt2*rn(i))
            y1d = 1000.0_r_kind*g*rnd(i)/dp
            y1  = rn(i)*1000.0_r_kind*g/dp
            IF (qevap(i) .GT. y1) THEN
              qevapd(i) = y1d
              qevap(i)  = y1
            ENDIF
          END IF

          delq2(i) = delqev(i) + 0.001_r_kind * qevap(i) * dp / g

          IF (rn(i) .GT. zero .AND. qcond(i) .LT. zero .AND. delq2(i) .GT. rntot(i)) THEN
            qevapd(i) = 1000.0_r_kind*g*(rntotd(i)-delqevd(i))/dp
            qevap(i)  = 1000.0_r_kind*g*(rntot(i)-delqev(i))/dp
          END IF

          IF (rn(i) .GT. zero .AND. qevap(i) .GT. zero) THEN
            q1d(i, k) = q1d(i, k) + qevapd(i)
            q1(i, k)  = q1(i, k) + qevap(i)
            t1d(i, k) = t1d(i, k) - elocp*qevapd(i)
            t1(i, k)  = t1(i, k) - elocp*qevap(i)
            rnd(i)    = rnd(i) - .001_r_kind*dp*qevapd(i)/g
            rn(i)     = rn(i) - .001_r_kind*qevap(i)*dp/g
            deltv(i)  = -(elocp*qevap(i)/dt2)
            delq(i)   = qevap(i)/dt2
            delqevd(i)= delqevd(i) + .001_r_kind*dp*qevapd(i)/g
            delqev(i) = delqev(i) + .001_r_kind*dp*qevap(i)/g
          END IF
        ENDIF !IF (k .LT. ktcon(i)) THEN

      END IF
    END DO
  END DO

  IF (ncloud .GT. 0) THEN
    DO k=1,km
      DO i=1,im
        IF (rn(i) .GT. zero) THEN
          IF (k .GT. kb(i) .AND. k .LE. ktcon(i)) THEN
            temd = dt2*(dellald(i, k)*xmb(i)+dellal(i, k)*xmbd(i))
            tem  = dellal(i, k)*xmb(i)*dt2
            qld(i, k, 1) = qld(i, k, 1) + temd
            ql(i, k, 1)  = ql(i, k, 1) + tem
            IF(ql(i,k,1) .lt. 1.0E-10_r_kind) then
               ql(i,k,1) = zero
               qld(i,k,1) = zero
            ENDIF
          END IF
        END IF
      END DO
    END DO
  END IF

  do k = 1, km
    do i = 1, im
       t_out(i,k) = t1(i,k)
       t_outd(i,k) = t1d(i,k)
       q_out(i,k) = q1(i,k)
       q_outd(i,k) = q1d(i,k)
       u_out(i,k) = u1(i,k)
       u_outd(i,k) = u1d(i,k)
       v_out(i,k) = v1(i,k)
       v_outd(i,k) = v1d(i,k)
       ql_out(i,k,1) = ql(i,k,1)
       ql_outd(i,k,1) = qld(i,k,1)
       ql_out(i,k,2) = ql(i,k,2)
       ql_outd(i,k,2) = zero
     enddo
  enddo

endif
  RETURN

END SUBROUTINE SASCNVN_TL
