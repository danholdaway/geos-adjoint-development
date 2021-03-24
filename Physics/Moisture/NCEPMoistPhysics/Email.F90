!Hi Dan,

!I've attached the moisture physics codes I developed so far here in NOAA.
!As I said, convection schemes need more work obviously. 
!The other two schemes (one for large scale condensation and one for precipitation) 
!seem working fine and they passed adjoint tests at least. 
!But I still need to check the jacobians for all these schemes using the methods you taught me yesterday.
!
!I will spend some time to digest what we talked about yesterday and will get back to you with more questions.
!
!Thank you and have a great weekend!
!
!Min-Jeong
!
!gfs_gscond_tl.f90




!gfs_gscond_ad.f90



!gfs_precip_tl.f90



!gfs_precip_ad.f90



!gfs_deepcnv_fw.f90
 subroutine sascnvn(im,ix,km,jcap,delt,del,prsl,ps,sl,ql, &
          q1,t1,u1,v1,rcs,rn,jkbcon0,jkbcon1,jktcon0,jktcon1,slimsk,  &
          dot,ncloud,t_out,q_out,ql_out,u_out,v_out,kmax,kbcon)

!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    sascnvn_fw    forward model for GFS deep convection scheme
!     prgmmr:   Min-Jeong Kim     org: np23                date: 2012-02-16
!
! abstract:  This subroutine contains the forward model for the
!            GFS Simplified Arakawa-Schubert (SAS) deep convection scheme
!
! program history log:
!   2012-02-25  kim - initial routine
!
!$$$

      USE kinds, only:r_single,r_kind,i_kind
      use CONSTANTS, only: CP, HVAP,  RD, FV, T0C, &
                     EPS, EPSM1, G=>grav, FACT1=>FACTOR1, FACT2=>FACTOR2,&
                    ONE, ZERO, HALF, tiny_r_kind, TWO, &
                    EL2ORC,elocp, init_constants_derived, init_constants


      implicit none

      integer            im, ix,  km, jcap, ncloud, &
                        kbot(im), ktop(im), kcnv(im) 

      real(r_kind)  :: delt
      real(r_kind) ::  ps(im),     del(ix,km),  prsl(ix,km), &
                          ql(ix,km,2), q1(ix,km),   t1(ix,km), &
                          u1(ix,km),  v1(ix,km),   rcs(im), &
                          cldwrk(im), rn(im),      slimsk(im),  &
                          dot(ix,km), sl(ix,km), &
            t_out(im,km),q_out(im,km),u_out(im,km),v_out(im,km) ,ql_out(ix,km,2)

      integer(i_kind) ::        i, j, indx, jmn, k, kk, latd, lond, km1

      real(r_kind) :: clam, cxlamu, xlamde, xlamdd
 
      real(r_kind) :: adw,     aup,     aafac, &
                          beta,    betal,   betas,  &
                                  dellat,   &
                          desdt,   deta,    detad,   dg, &
                          dh,      dhh,     dlnsig,  dp, &
                          dq,      dqsdp,   dqsdt,   dt, &
                          dt2,     dtmax,   dtmin,   dv1h, &
                          dv1q,    dv2h,    dv2q,    dv1u, &
                          dv1v,    dv2u,    dv2v,    dv3q, &
                          dv3h,    dv3u,    dv3v,   &
                          dz,      dz1,     e1,      edtmax, &
                          edtmaxl, edtmaxs,   & 
                          es,  es0,    etah,    &
                          evef,    evfact,  evfactl,  &
                             factor,  fjcap,   fkm, &
                                   gamma,   pprime, &
                          qlk,     qrch,    qs,       &
                          rain,    rfact,   shear,   tem1, &
                          tem2,        val,     val1, &
                          val2,    w1,      w1l,     w1s, &
                          w2,      w2l,     w2s,     w3, &
                          w3l,     w3s,     w4,      w4l, &
                          w4s,     xdby,    xpw,     xpwd, &
                          xqrch,   mbdt,    tem, &
                          ptem,    ptem1,   pgcon

      integer(i_kind) ::   kb(im), kbcon(im), kbcon1(im),   &
                          ktcon(im), ktcon1(im),  &
                          jmin(im), lmin(im), kbmax(im), &
                          kbm(im), kmax(im)

      integer(i_kind) :: jkbcon0(im),jkbcon1(im),jktcon0(im),jktcon1(im)
      real(r_kind) :: aa1(im),     acrt(im),   acrtfct(im),   &
                          delhbar(im), delq(im),   delq2(im), &
                          delqbar(im), delqev(im), deltbar(im), &
                          deltv(im),   dtconv(im), edt(im), &
                          edto(im),    edtx(im),   fld(im), &
                          hcdo(im,km), hmax(im),   hmin(im), &
                          ucdo(im,km), vcdo(im,km),aa2(im),  &
                          pbcdif(im),  pdot(im),   po(im,km), &
                          pwavo(im),   pwevo(im),  xlamud(im), & 
                          qcdo(im,km), qcond(im),  qevap(im), &
                          rntot(im),   vshear(im), xaa0(im), &
                          xk(im),      xlamd(im), &
                          xmb(im),     xmbmax(im), xpwav(im), &
                          xpwev(im),   delubar(im),delvbar(im)

      real(r_kind), dimension(im) :: cnvscl_1,cnvscl_2,cnvscl_3, &
                                     cnvscl_4,cnvscl_5,cnvscl_6, & 
                                     cnvscl_7,cnvscl_8 

!  physical parameters
!     real(r_kind), parameter :: cpoel=cp/hvap,elocp=hvap/cp, &
!               el2orc=hvap*hvap/(rv*cp),terr=0.,c0=.002,c1=.002,delta=fv, &
!      fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c,   &
!      cthk=150.,cincrmax=180.,cincrmin=120.,dthk=25.

!  local variables and arrays
      real(r_kind) :: pfld(im,km),    to(im,km),     qo(im,km), &
                uo(im,km),   vo(im,km),    qeso(im,km)
!  cloud water
      real(r_kind) :: qlko_ktcon(im), dellal(im,km), tvo(im,km),  &
                          dbyo(im,km),    zo(im,km),     xlamue(im,km), &
                          fent1(im,km),   fent2(im,km),  frh(im,km), &
                          heo(im,km),     heso(im,km), &
                          qrcd(im,km),    dellah(im,km), dellaq(im,km), &
                          dellau(im,km),  dellav(im,km), hcko(im,km), &
                          ucko(im,km),    vcko(im,km),   qcko(im,km), &
                          eta(im,km),     etad(im,km),   zi(im,km), &
                          qrcdo(im,km),   pwo(im,km),    pwdo(im,km), &
                          tx1(im),        sumx(im)

      logical ::  cnvflg(im), flg(im), flg1(im),flg2(im),flg3(im)

      real(r_kind) :: pcrit(15), acritt(15), acrit(15)
      DATA pcrit /850.0_r_kind, 800.0_r_kind, 750.0_r_kind, &
              700.0_r_kind, 650.0_r_kind, 600.0_r_kind,&
              550.0_r_kind, 500.0_r_kind, 450.0_r_kind, 400.0_r_kind, &
              350.0_r_kind, 300.0_r_kind, 250.0_r_kind, 200.0_r_kind, 150.0_r_kind/
      DATA acritt /.0633_r_kind, .0445_r_kind, .0553_r_kind, .0664_r_kind, .075_r_kind, &
               .1082_r_kind, .1521_r_kind, .2216_r_kind, &
               .3151_r_kind, .3677_r_kind, .41_r_kind, .5255_r_kind, &
               .7663_r_kind, 1.1686_r_kind, 1.6851_r_kind/

      real(r_kind), parameter ::tf=233.16_r_kind, tcr=263.16_r_kind, tcrf=one/(tcr-tf)

      real(r_kind), parameter :: cpoel=cp/hvap
      real(r_kind) cincr, cincrmax, cincrmin, dum, dthk, cthk, c0, c1, &
                   terr0, delta

!------------------------
 km1 = km - 1

 call init_constants_derived
 call init_constants(.false.)

      km1 = km - 1
      terr0=zero
      c0 = 0.002_r_kind
      c1= 0.002_r_kind
      delta = fv
      cthk=150.0_r_kind
      cincrmax=180.0_r_kind
      cincrmin=120.0_r_kind
      dthk=25.0_r_kind

!  initialize arrays

      do i=1,im
        cnvflg(i) = .true.
        rn(i)=zero
        kbot(i)=km+1
        ktop(i)=0
        kbcon(i)=km
        ktcon(i)=1
        dtconv(i) = 3600.0_r_kind
        cldwrk(i) = zero
        pdot(i) = zero
        pbcdif(i)= zero
        lmin(i) = 1
        jmin(i) = 1
        qlko_ktcon(i) = zero
        edt(i)  = zero
        edto(i) = zero
        edtx(i) = zero
        acrt(i) = zero
        acrtfct(i) = one
        aa1(i)  = zero
        aa2(i)  = zero
        xaa0(i) = zero
        pwavo(i)= zero
        pwevo(i)= zero
        xpwav(i)= zero
        xpwev(i)= zero 
        vshear(i) = zero
      enddo

      do k = 1, km
        do i = 1, im
          t_out(i,k) = t1(i,k)
          q_out(i,k) = q1(i,k)
          u_out(i,k) = u1(i,k)
          v_out(i,k) = v1(i,k)
          ql_out(i,k,1) = ql(i,k,1)
          ql_out(i,k,2) = ql(i,k,2)
        enddo
      enddo

      do k = 1, 15
        acrit(k) = acritt(k) * (975. - pcrit(k))
      enddo

      dt2 = delt
      val   = 1200.0_r_kind
      dtmin = max(dt2, val )
      val   = 3600.0_r_kind
      dtmax = max(dt2, val )
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

      pgcon   = 0.55_r_kind    ! Zhang & Wu (2003,JAS)
      fjcap   = (float(jcap) / 126.0_r_kind) ** 2
      val     = one 
      fjcap   = max(fjcap,val)
      fkm     = (float(km) / 28.0_r_kind) ** 2
      fkm     = max(fkm,val)
      w1l     = -8.e-3_r_kind
      w2l     = -4.e-2_r_kind
      w3l     = -5.e-3_r_kind
      w4l     = -5.e-4_r_kind
      w1s     = -2.e-4_r_kind
      w2s     = -2.e-3_r_kind
      w3s     = -1.e-3_r_kind
      w4s     = -2.e-5_r_kind

!  define top layer for search of the downdraft originating layer
!  and the maximum thetae for updraft

      do i=1,im
        kbmax(i) = km
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = one / ps(i)
      enddo
     
      do k = 1, km
        do i=1,im
          if (prsl(i,k)*tx1(i) .gt. 0.04_r_kind) kmax(i)  = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.45_r_kind) kbmax(i) = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.70_r_kind) kbm(i)   = k + 1
        enddo
      enddo

      do i=1,im
        kbmax(i) = min(kbmax(i),kmax(i))
        kbm(i)   = min(kbm(i),kmax(i))
      enddo

      do k = 1, km
        do i = 1, im
            pfld(i,k) = prsl(i,k) * 10.0_r_kind
            eta(i,k)  = one
            fent1(i,k)= one
            fent2(i,k)= one
            frh(i,k)  = zero
            hcko(i,k) = zero 
            qcko(i,k) = zero
            ucko(i,k) = zero
            vcko(i,k) = zero
            etad(i,k) = one
            hcdo(i,k) = zero
            qcdo(i,k) = zero
            ucdo(i,k) = zero
            vcdo(i,k) = zero
            qrcd(i,k) = zero
            qrcdo(i,k)= zero
            dbyo(i,k) = zero
            pwo(i,k)  = zero
            pwdo(i,k) = zero
            dellal(i,k) = zero
            to(i,k)   = t1(i,k)
            qo(i,k)   = q1(i,k)
            uo(i,k)   = u1(i,k) * rcs(i)
            vo(i,k)   = v1(i,k) * rcs(i)
        enddo
      enddo

!  column variables
!  p is pressure of the layer (mb)
!  t is temperature at t-dt (k)..tn
!  q is mixing ratio at t-dt (kg/kg)..qn
!  to is temperature at t+dt (k)... this is after advection and turbulan
!  qo is mixing ratio at t+dt (kg/kg)..q1

      do k = 1, km
        do i=1,im
            call fpvsx_ad(to(i,k),es0,dum, dum, .false.)      ! fpvs is in pa
            qeso(i,k) = es0*10.0_r_kind
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val1      = 1.e-8_r_kind
            qeso(i,k) = max(qeso(i,k), val1)
            val2      = 1.e-10_r_kind
            qo(i,k)   = max(qo(i,k), val2 )
            tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
        enddo
      enddo

!  hydrostatic height assume zero terr and initially assume
!    updraft entrainment rate as an inverse function of height 

      DO I = 1, IM   !kim
        ZO(I,1) = ZERO - log(SL(I,1))*RD/G * TVO(I,1)  !kim
      ENDDO   !kim

      do k = 2, km   !kim
        do i=1,im
        ZO(I,K) = ZO(I,K-1) -   &
         log(SL(I,K)/SL(I,K-1))*RD/G * HALF*(TVO(I,K) + TVO(I,K-1)) !kim
        enddo
      enddo

      do k = 1, km1
        do i=1,im
          zi(i,k) = half*(zo(i,k)+zo(i,k+1))
          xlamue(i,k) = clam / zi(i,k)
        enddo
      enddo

!  compute moist static energy

      do k = 1, km
        do i=1,im
          if (k .le. kmax(i)) then
            tem       = g * zo(i,k) + cp * to(i,k)
            heo(i,k)  = tem  + hvap * qo(i,k)
            heso(i,k) = tem  + hvap * qeso(i,k)
          endif
        enddo
      enddo

!  determine level with largest moist static energy
!  this is the level where updraft starts

      do i=1,im
        kb(i)   = 1
      enddo

      do k = 1, km1
        do i=1,im
          if (k .le. kmax(i)-1) then
            dz      = half * (zo(i,k+1) - zo(i,k))
            dp      = half * (pfld(i,k+1) - pfld(i,k))
            call fpvsx_ad(to(i,k+1),es0,dum, dum, .false.)      
            es = es0 * 10.0_r_kind
            pprime  = pfld(i,k+1) + epsm1 * es
            qs      = eps * es / pprime
            dqsdp   = - qs / pprime
            desdt   = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt   = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma   = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt      = (g * dz + hvap * dqsdp * dp) / (cp * (one + gamma))
            dq      = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = half * (pfld(i,k) + pfld(i,k+1))

          endif
        enddo
      enddo

      do k = 1, km1
        do i=1,im
          if (k .le. kmax(i)-1) then
            call fpvsx_ad(to(i,k),es0,dum,dum,.false.)     
            qeso(i,k) = es0 * 10.0_r_kind
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1*qeso(i,k))
            val1      = 1.e-8_r_kind
            qeso(i,k) = max(qeso(i,k), val1)
            val2      = 1.e-10_r_kind
            qo(i,k)   = max(qo(i,k), val2 )
            frh(i,k)  = one - min(qo(i,k)/qeso(i,k), one)
            heo(i,k)  = half * g * (zo(i,k) + zo(i,k+1)) +  &
                       cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = half * g * (zo(i,k) + zo(i,k+1)) +  &
                       cp * to(i,k) + hvap * qeso(i,k)
            uo(i,k)   = half * (uo(i,k) + uo(i,k+1))
            vo(i,k)   = half * (vo(i,k) + vo(i,k+1))

          endif
        enddo
      enddo

!  look for the level of free convection as cloud base

      do i=1,im
        flg(i)   = .true.
        cnvscl_1(i) = one
        cnvscl_2(i) = one
        cnvscl_3(i) = one
        cnvscl_4(i) = one
        cnvscl_5(i) = one
        cnvscl_6(i) = one
        cnvscl_7(i) = one
        cnvscl_8(i) = one
        kbcon(i) = kmax(i)
      enddo

      do k = 1, km1
        do i=1,im
          if (flg(i).and.k.le.kbmax(i)) then
            if(k.gt.kb(i).and.heo(i,kb(i)).gt.heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
    
    do i = 1, im 
      if(jkbcon0(i) .gt. 0) then
          kbcon(i) = jkbcon0(i)
      else
          jkbcon0(i)=kbcon(i)
      endif
    enddo

      do i=1,im
        if(kbcon(i).eq.kmax(i)) cnvscl_1(i) = tiny_r_kind 
      enddo



!  determine critical convective inhibition
!  as a function of vertical velocity at cloud base.

      do i=1,im
          pdot(i)  = 10.0_r_kind* dot(i,kbcon(i))
      enddo

      do i=1,im
          if(slimsk(i).eq.1.) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
          if(pdot(i).le.w4) then
            tem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            tem = - (pdot(i) + w4) / (w4 - w3)
          else
            tem = zero 
          endif
          val1    =   -one 
          tem = max(tem,val1)
          val2    =   one 
          tem = min(tem,val2)
          tem = one - tem
          tem1= half*(cincrmax-cincrmin)
          cincr = cincrmax - tem * tem1
          pbcdif(i) = pfld(i,kb(i)) - pfld(i,kbcon(i))
          if(pbcdif(i).gt.cincr) cnvscl_2(i) = tiny_r_kind 
          cnvscl_2(i) = cnvscl_2(i)*cnvscl_1(i)
      enddo


!  assume that updraft entrainment rate above cloud base is
!    same as that at cloud base

      do k = 2, km1
        do i=1,im
          if((k.gt.kbcon(i).and.k.lt.kmax(i))) then
              xlamue(i,k) = xlamue(i,kbcon(i))
          endif
        enddo
      enddo

!  assume the detrainment rate for the updrafts to be same as
!  the entrainment rate at cloud base

      do i = 1, im
          xlamud(i) = xlamue(i,kbcon(i))
      enddo

!  functions rapidly decreasing with height, mimicking a cloud ensemble
!    (Bechtold et al., 2008)

      do k = 2, km1
        do i=1,im
          if((k.gt.kbcon(i).and.k.lt.kmax(i))) then
              tem = qeso(i,k)/qeso(i,kbcon(i))
              fent1(i,k) = tem**2
              fent2(i,k) = tem**3
          endif
        enddo
      enddo

!  final entrainment rate as the sum of turbulent part and organized entrainment
!    depending on the environmental relative humidity
!    (Bechtold et al., 2008)

      do k = 2, km1
        do i=1,im
          if((k.ge.kbcon(i).and.k.lt.kmax(i))) then
              tem = cxlamu * frh(i,k) * fent2(i,k)
              xlamue(i,k) = xlamue(i,k)*fent1(i,k) + tem
          endif
        enddo
      enddo

!  determine updraft mass flux for the subcloud layers

      do k = km1, 1, -1
        do i = 1, im
            if(k.lt.kbcon(i).and.k.ge.kb(i)) then
              dz       = zi(i,k+1) - zi(i,k)
              ptem     = half*(xlamue(i,k)+xlamue(i,k+1))-xlamud(i)
              eta(i,k) = eta(i,k+1) / (one + ptem * dz)
            endif
        enddo
      enddo

!  compute mass flux above cloud base

      do k = 2, km1
        do i = 1, im
           if(k.gt.kbcon(i).and.k.lt.kmax(i)) then
              dz       = zi(i,k) - zi(i,k-1)
              ptem     = half*(xlamue(i,k)+xlamue(i,k-1))-xlamud(i)
              eta(i,k) = eta(i,k-1) * (one + ptem * dz)
           endif
        enddo
      enddo

!  compute updraft cloud properties

      do i = 1, im
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
          pwavo(i)     = zero 
      enddo

!  cloud property is modified by the entrainment process

      do k = 2, km1
        do i = 1, im
            if(k.gt.kb(i).and.k.lt.kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = half * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = half * xlamud(i) * dz
              factor = one + tem - tem1
              ptem = half * tem + pgcon
              ptem1= half * tem - pgcon
              hcko(i,k) = ((one-tem1)*hcko(i,k-1)+tem*0.5*  &
                          (heo(i,k)+heo(i,k-1)))/factor
              ucko(i,k) = ((one-tem1)*ucko(i,k-1)+ptem*uo(i,k)  &
                          +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((one-tem1)*vcko(i,k-1)+ptem*vo(i,k)  &
                          +ptem1*vo(i,k-1))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)
            endif
        enddo
      enddo

!   taking account into convection inhibition due to existence of
!   dry layers below cloud base

      do i=1,im
        kbcon1(i) = kmax(i)
        flg1(i) = .true.
      enddo

      do k = 2, km1
       do i=1,im
        if (k.lt.kmax(i)) then
          if(flg1(i) .and. k.ge.kbcon(i) .and. dbyo(i,k).gt.zero) then
            kbcon1(i) = k
            flg1(i) = .false. 
          endif
        endif
       enddo
      enddo 
     
     do i = 1, im
       if(jkbcon1(i) .gt.0 ) then
          kbcon1(i) = jkbcon1(i)
       else
          jkbcon1(i) = kbcon1(i)
       endif
     enddo 

      do i=1,im
          if(kbcon1(i).eq.kmax(i)) cnvscl_2(i) = tiny_r_kind 
         
      enddo
      do i=1,im
          tem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          cnvscl_3(i) = ONE - HALF**((dthk/tem)**200)
          cnvscl_3(i) = cnvscl_3(i)*cnvscl_2(i)
      enddo


!  determine first guess cloud top as the level of zero buoyancy

      do i = 1, im
        ktcon(i) = 1
        flg2(i) = .true.
      enddo

     do k = 2, km1
      do i = 1, im
        if (k .lt. kmax(i)) then
          if(flg2(i) .and. k.gt.kbcon1(i) .and. dbyo(i,k).lt.zero) then
             ktcon(i) = k
             flg2(i) = .false. 
          endif
        endif
      enddo
      enddo
 
     do i = 1, im
       if(jktcon0(i) .gt. 0) then
            ktcon(i) = jktcon0(i)
       else
            jktcon0(i) = ktcon(i)
       endif
     enddo 

if(ktcon(1) .gt. 2) then

      do i = 1, im
          tem = pfld(i,kbcon(i))-pfld(i,ktcon(i))
          cnvscl_4(i) = ONE - HALF**((tem/cthk)**200)
          cnvscl_4(i) = cnvscl_4(i)*cnvscl_3(i)
      enddo


!  search for downdraft originating level above theta-e minimum

      do i = 1, im
           hmin(i) = heo(i,kbcon1(i))
           lmin(i) = kbmax(i)
           jmin(i) = kbmax(i)
      enddo

      do k = 2, km1
        do i = 1, im
          if (k .le. kbmax(i)) then
            if(k.gt.kbcon1(i).and.heo(i,k).lt.hmin(i)) then
               lmin(i) = k + 1
               hmin(i) = heo(i,k)
            endif
          endif
        enddo
      enddo

!  make sure that jmin(i) is within the cloud

      do i = 1, im
          jmin(i) = min(lmin(i),ktcon(i)-1)
          jmin(i) = max(jmin(i),kbcon1(i)+1)
          if(jmin(i).ge.ktcon(i)) cnvscl_4(i) = tiny_r_kind 
      enddo

!  specify upper limit of mass flux at cloud base

      do i = 1, im
          k = kbcon(i)
          dp = 1000.0_r_kind * del(i,k)*ps(i)
          xmbmax(i) = dp / (g * dt2)
      enddo

!  compute cloud moisture property and precipitation

      do i = 1, im
          aa1(i) = zero
          qcko(i,kb(i)) = qo(i,kb(i))
      enddo

      do k = 2, km1
        do i = 1, im
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k) + gamma * dbyo(i,k) / (hvap * (one + gamma))
              tem  = half * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = half * xlamud(i) * dz
              factor = one + tem - tem1
              qcko(i,k) = ((one-tem1)*qcko(i,k-1)+tem*half*(qo(i,k)+qo(i,k-1)))/factor
              dq = eta(i,k) * (qcko(i,k) - qrch)

!  check if there is excess moisture to release latent heat

              if(k .ge. kbcon(i) .and. dq .gt. zero) then
                etah = half * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt. 0 .and. k .gt.jmin(i)) then
                  dp = 1000.0_r_kind * del(i,k)*ps(i)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                aa1(i) = aa1(i) - dz * g * qlk
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
                pwavo(i) = pwavo(i) + pwo(i,k)
              endif
            endif
        enddo
      enddo


!  calculate cloud work function

      do k = 2, km1
        do i = 1, im
            if(k.ge.kbcon(i) .and. k.lt.ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  one + delta * cp * gamma* to(i,k) / hvap
              aa1(i) = aa1(i) + dz1 * (g / (cp * to(i,k))) &
                      * dbyo(i,k) / (one + gamma) * rfact
              val = zero 
              aa1(i)=aa1(i)+ dz1 * g * delta * max(val,(qeso(i,k) - qo(i,k)))
            endif
        enddo
      enddo

      do i = 1, im
       cnvscl_5(i) = half*(one+tanh(aa1(i))) 
       cnvscl_5(i) = cnvscl_5(i)*cnvscl_4(i)
      enddo


!  estimate the onvective overshooting as the level 
!    where the [aafac * cloud work function] becomes zero,
!    which is the final cloud top

      do i = 1, im
        aa2(i) = aafac * aa1(i)
        flg3(i) = .true.
      enddo

      do i = 1, im
        ktcon1(i) = kmax(i) - 1
      enddo

      do i = 1, im
        if(jktcon1(i) .gt. 0) then
          ktcon1(i) = jktcon1(i)
          do k = 2, km1
            if(k.ge.ktcon(i).and.k.lt.kmax(i) .and. k .le. ktcon1(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  one + delta * cp * gamma * to(i,k) / hvap
              aa2(i) = aa2(i) + dz1 * (g / (cp * to(i,k)))  &
                      * dbyo(i,k) / (one + gamma) * rfact
            endif
          enddo
        else
          do k = 2, km1
            if(flg3(i) .and. k.ge.ktcon(i).and.k.lt.kmax(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  one + delta * cp * gamma * to(i,k) / hvap
              aa2(i) = aa2(i) + dz1 * (g / (cp * to(i,k)))  &
                      * dbyo(i,k) / (one + gamma) * rfact
              if(aa2(i) .lt. zero) then
                ktcon1(i) = k
                flg3(i) = .false.
              endif
            endif
          enddo
          jktcon1(i) = ktcon1(i)
        endif
      enddo 


!  compute cloud moisture property, detraining cloud water 
!    and precipitation in overshooting layers 

      do k = 2, km1
        do i = 1, im
            if(k.ge.ktcon(i).and.k.lt.ktcon1(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)+ gamma * dbyo(i,k) / (hvap * (one + gamma))

              tem  = half * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = half * xlamud(i) * dz
              factor = one + tem - tem1
              qcko(i,k) = ((one-tem1)*qcko(i,k-1)+tem*half*    &
                          (qo(i,k)+qo(i,k-1)))/factor
              dq = eta(i,k) * (qcko(i,k) - qrch)

!  check if there is excess moisture to release latent heat

              if(dq .gt. zero) then
                etah = half * (eta(i,k) + eta(i,k-1))
                if(ncloud .gt. zero) then
                  dp = 1000.0_r_kind * del(i,k)*ps(i)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
                pwavo(i) = pwavo(i) + pwo(i,k)
              endif
            endif
        enddo
      enddo


! exchange ktcon with ktcon1

      do i = 1, im
          kk = ktcon(i)
          ktcon(i) = ktcon1(i)
          ktcon1(i) = kk
      enddo


!  this section is ready for cloud water

      if(ncloud.gt.0) then
        do i = 1, im
          k = ktcon(i) - 1
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
          qrch = qeso(i,k) + gamma * dbyo(i,k) / (hvap * (one + gamma))
          dq = qcko(i,k) - qrch

!  check if there is excess moisture to release latent heat

          if(dq .gt. zero) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
          endif
        enddo
      endif

!------- downdraft calculations
!--- compute precipitation efficiency in terms of windshear

      do i = 1, im
          vshear(i) = zero 
      enddo

      do k = 2, km
        do i = 1, im
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              shear= sqrt((uo(i,k)-uo(i,k-1)) ** 2 &
                       + (vo(i,k)-vo(i,k-1)) ** 2)
              vshear(i) = vshear(i) + shear
            endif
        enddo
      enddo

      do i = 1, im
          vshear(i) = 1.e3_r_kind * vshear(i) / (zi(i,ktcon(i))-zi(i,kb(i)))
          e1=1.591_r_kind-.639_r_kind*vshear(i)  &
            +.0953_r_kind*(vshear(i)**2)-.00496*(vshear(i)**3)
          edt(i)=one-e1
          val =  0.9_r_kind
          edt(i) = min(edt(i),val)
          val =  zero 
          edt(i) = max(edt(i),val)
          edto(i)=edt(i)
          edtx(i)=edt(i)
      enddo

!  determine detrainment rate between 1 and kbcon

      do i = 1, im
          sumx(i) = zero
      enddo

      do k = 1, km1
      do i = 1, im
        if(k.ge.1.and.k.lt.kbcon(i)) then
          dz = zi(i,k+1) - zi(i,k)
          sumx(i) = sumx(i) + dz
        endif
      enddo
      enddo

      do i = 1, im
        beta = betas
        if(slimsk(i) .eq. one) beta = betal
          dz  = (sumx(i)+zi(i,1))/float(kbcon(i))
          tem = one/float(kbcon(i))
          xlamd(i) = (one-beta**tem)/dz
      enddo

!  determine downdraft mass flux

      do k = km1, 1, -1
        do i = 1, im
          if (k .le. kmax(i)-1) then
           if(k.lt.jmin(i).and.k.ge.kbcon(i)) then
              dz        = zi(i,k+1) - zi(i,k)
              ptem      = xlamdd - xlamde
              etad(i,k) = etad(i,k+1) * (one - ptem * dz)
           else if(k.lt.kbcon(i)) then
              dz        = zi(i,k+1) - zi(i,k)
              ptem      = xlamd(i) + xlamdd - xlamde
              etad(i,k) = etad(i,k+1) * (one - ptem * dz)
           endif
          endif
        enddo
      enddo

!--- downdraft moisture properties

      do i = 1, im
          jmn = jmin(i)
          hcdo(i,jmn) = heo(i,jmn)
          qcdo(i,jmn) = qo(i,jmn)
          qrcdo(i,jmn)= qeso(i,jmn)
          ucdo(i,jmn) = uo(i,jmn)
          vcdo(i,jmn) = vo(i,jmn)
          pwevo(i) = zero 
      enddo

      do k = km1, 1, -1
        do i = 1, im
          if ( k.lt.jmin(i)) then
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = half * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = half * (xlamd(i)+xlamdd) * dz
              endif
              factor = one + tem - tem1
              ptem = half * tem - pgcon
              ptem1= half * tem + pgcon
              hcdo(i,k) = ((one-tem1)*hcdo(i,k+1)+tem*half*   &
                          (heo(i,k)+heo(i,k+1)))/factor
              ucdo(i,k) = ((one-tem1)*ucdo(i,k+1)+ptem*uo(i,k+1) &
                          +ptem1*uo(i,k))/factor
              vcdo(i,k) = ((one-tem1)*vcdo(i,k+1)+ptem*vo(i,k+1)  &
                          +ptem1*vo(i,k))/factor
              dbyo(i,k) = hcdo(i,k) - heso(i,k)
          endif
        enddo
      enddo

      do k = km1, 1, -1
        do i = 1, im
          if (k.lt.jmin(i)) then
              gamma      = el2orc * qeso(i,k) / (to(i,k)**2)
              qrcdo(i,k) = qeso(i,k)+ (one/hvap)*(gamma/(one+gamma))*dbyo(i,k)

              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = half * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = half * (xlamd(i)+xlamdd) * dz
              endif
              factor = one + tem - tem1
              qcdo(i,k) = ((one-tem1)*qcdo(i,k+1)+tem*half*  &
                          (qo(i,k)+qo(i,k+1)))/factor

              pwdo(i,k)  = etad(i,k+1) * (qcdo(i,k) - qrcdo(i,k))
 
              qcdo(i,k)  = qrcdo(i,k)
              pwevo(i)   = pwevo(i) + pwdo(i,k)
          endif
        enddo
      enddo

!--- final downdraft strength dependent on precip
!--- efficiency (edt), normalized condensate (pwav), and
!--- evaporate (pwev)

      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i) .eq. zero) edtmax = edtmaxs
          if(pwevo(i) .lt. zero) then
            edto(i) = -edto(i) * pwavo(i) / pwevo(i)
            edto(i) = min(edto(i),edtmax)
          else
            edto(i) = zero 
          endif
      enddo

!--- downdraft cloudwork functions

      do k = km1, 1, -1
        do i = 1, im
          if (k .lt. jmin(i)) then
              gamma = el2orc * qeso(i,k) / to(i,k)**2
              dhh=hcdo(i,k)
              dt=to(i,k)
              dg=gamma
              dh=heso(i,k)
              dz=-one*(zo(i,k+1)-zo(i,k))
              aa1(i)=aa1(i)+edto(i)*dz*(g/(cp*dt))*((dhh-dh)/(one+dg)) &
                    *(one+delta*cp*dg*dt/hvap)
              val=zero
              aa1(i)=aa1(i)+edto(i)*  & 
             dz*g*delta*max(val,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo

      do i = 1, im
         cnvscl_6(i) = half*(one+tanh(aa1(i)))
         cnvscl_6(i) = cnvscl_6(i)*cnvscl_5(i)
      enddo


!--- what would the change be, that a cloud with unit mass
!--- will do to the environment?

      do k = 1, km
        do i = 1, im
          if(k .le. kmax(i)) then
            dellah(i,k) = zero
            dellaq(i,k) = zero 
            dellau(i,k) = zero 
            dellav(i,k) = zero 
          endif
        enddo
      enddo

      do i = 1, im
          dp = 1000. * del(i,1)*ps(i)
          dellah(i,1) = edto(i) * etad(i,1) * (hcdo(i,1)  &
                        - heo(i,1)) * g / dp
          dellaq(i,1) = edto(i) * etad(i,1) * (qcdo(i,1)  &
                        - qo(i,1)) * g / dp
          dellau(i,1) = edto(i) * etad(i,1) * (ucdo(i,1)  &
                        - uo(i,1)) * g / dp
          dellav(i,1) = edto(i) * etad(i,1) * (vcdo(i,1) &
                        - vo(i,1)) * g / dp
      enddo

!--- changed due to subsidence and entrainment

      do k = 2, km1
        do i = 1, im
          if (k.lt.ktcon(i)) then
              aup = one
              if(k.le.kb(i)) aup = zero
              adw = one
              if(k.gt.jmin(i)) adw = zero
              dp = 1000.0_r_kind * del(i,k)*ps(i)
              dz = zi(i,k) - zi(i,k-1)

              dv1h = heo(i,k)
              dv2h = half * (heo(i,k) + heo(i,k-1))
              dv3h = heo(i,k-1)
              dv1q = qo(i,k)
              dv2q = half * (qo(i,k) + qo(i,k-1))
              dv3q = qo(i,k-1)
              dv1u = uo(i,k)
              dv2u = half * (uo(i,k) + uo(i,k-1))
              dv3u = uo(i,k-1)
              dv1v = vo(i,k)
              dv2v = half * (vo(i,k) + vo(i,k-1))
              dv3v = vo(i,k-1)

              tem  = half * (xlamue(i,k)+xlamue(i,k-1))
              tem1 = xlamud(i)

              if(k.le.kbcon(i)) then
                ptem  = xlamde
                ptem1 = xlamd(i)+xlamdd
              else
                ptem  = xlamde
                ptem1 = xlamdd
              endif

              dellah(i,k) = dellah(i,k) +    &
          ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1h   &
         - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3h   &
         - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2h*dz  &
         +  aup*tem1*eta(i,k-1)*half*(hcko(i,k)+hcko(i,k-1))*dz    &
         +  adw*edto(i)*ptem1*etad(i,k)*half*(hcdo(i,k)+hcdo(i,k-1))*dz  &
              ) *g/dp

              dellaq(i,k) = dellaq(i,k) +   &
          ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1q  &
         - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3q  &
         - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2q*dz   &
         +  aup*tem1*eta(i,k-1)*half*(qcko(i,k)+qcko(i,k-1))*dz    &
         +  adw*edto(i)*ptem1*etad(i,k)*half*(qrcdo(i,k)+qrcdo(i,k-1))*dz  &
              ) *g/dp

              dellau(i,k) = dellau(i,k) +   &
          ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1u  &
         - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3u  &
         - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2u*dz  &
         +  aup*tem1*eta(i,k-1)*half*(ucko(i,k)+ucko(i,k-1))*dz   &
         +  adw*edto(i)*ptem1*etad(i,k)*half*(ucdo(i,k)+ucdo(i,k-1))*dz  &
         -  pgcon*(aup*eta(i,k-1)-adw*edto(i)*etad(i,k))*(dv1u-dv3u)  &
              ) *g/dp

              dellav(i,k) = dellav(i,k) +    &
          ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1v  &
         - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3v  &
         - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2v*dz  &
         +  aup*tem1*eta(i,k-1)*half*(vcko(i,k)+vcko(i,k-1))*dz   &
         +  adw*edto(i)*ptem1*etad(i,k)*half*(vcdo(i,k)+vcdo(i,k-1))*dz   &
         -  pgcon*(aup*eta(i,k-1)-adw*edto(i)*etad(i,k))*(dv1v-dv3v)  &
             ) *g/dp

          endif
        enddo
      enddo
!c
!c------- cloud top
!c
      do i = 1, im
          indx = ktcon(i)
          dp = 1000. * del(i,indx)*ps(i)
          dv1h = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) *   &
                          (hcko(i,indx-1) - dv1h) * g / dp
          dv1q = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *   &
                          (qcko(i,indx-1) - dv1q) * g / dp
          dv1u = uo(i,indx-1)
          dellau(i,indx) = eta(i,indx-1) *    &
                          (ucko(i,indx-1) - dv1u) * g / dp
          dv1v = vo(i,indx-1)
          dellav(i,indx) = eta(i,indx-1) *    &
                          (vcko(i,indx-1) - dv1v) * g / dp
!c
!c  cloud water
!c
          dellal(i,indx) = eta(i,indx-1) * qlko_ktcon(i) * g / dp
      enddo
!c
!c------- final changed variable per unit mass flux
!c
      do k = 1, km
        do i = 1, im
          if (k .le. kmax(i)) then
            if(k.gt.ktcon(i)) then
              qo(i,k) = q1(i,k)
              to(i,k) = t1(i,k)
            endif
            if(k.le.ktcon(i)) then
              qo(i,k) = dellaq(i,k) * mbdt + q1(i,k)
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              to(i,k) = dellat * mbdt + t1(i,k)
              val   =           1.e-10
              qo(i,k) = max(qo(i,k), val  )
            endif
          endif
        enddo
      enddo
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c
!c--- the above changed environment is now used to calulate the
!c--- effect the arbitrary cloud (with unit mass flux)
!c--- would have on the stability,
!c--- which then is used to calculate the real mass flux,
!c--- necessary to keep this change in balance with the large-scale
!c--- destabilization.
!c
!c--- environmental conditions again, first heights
!c
      do k = 1, km
        do i = 1, im
          if(k .le. kmax(i)) then
            call fpvsx_ad(to(i,k),es0,dum,dum,.false.)      ! fpvs is in pa
            qeso(i,k) = es0*10.0_r_kind
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k)+epsm1*qeso(i,k))
            val       =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
!c
!c--- moist static energy
!c
      do k = 1, km1
        do i = 1, im
          if(k .le. kmax(i)-1) then
            dz = half * (zo(i,k+1) - zo(i,k))
            dp = half * (pfld(i,k+1) - pfld(i,k))
            call fpvsx_ad(to(i,k+1),es0,dum,dum,.false.)      ! fpvs is in pa
            es = es0*10.0_r_kind
            pprime = pfld(i,k+1) + epsm1 * es
            qs = eps * es / pprime
            dqsdp = - qs / pprime
            desdt = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = half * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
      do k = 1, km1
        do i = 1, im
          if(k .le. kmax(i)-1) then
            call fpvsx_ad(to(i,k),es0,dum,dum,.false.)      ! fpvs is in pa
            qeso(i,k) = es0*10.0_r_kind
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1 * qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
            heo(i,k)   = half * g * (zo(i,k) + zo(i,k+1)) + &
                         cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = half * g * (zo(i,k) + zo(i,k+1)) +  &
                       cp * to(i,k) + hvap * qeso(i,k)
          endif
        enddo
      enddo
      do i = 1, im
          k = kmax(i)
          heo(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qo(i,k)
          heso(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qeso(i,k)
      enddo
!c
!c**************************** static control
!c
!c------- moisture and cloud work functions
!c
      do i = 1, im
          xaa0(i) = 0.
          xpwav(i) = 0.
      enddo

      do i = 1, im
          indx = kb(i)
          hcko(i,indx) = heo(i,indx)
          qcko(i,indx) = qo(i,indx)
      enddo
      do k = 2, km1
        do i = 1, im
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              tem  = half * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = half * xlamud(i) * dz
              factor = one + tem - tem1
              hcko(i,k) = ((one-tem1)*hcko(i,k-1)+tem*half*   &
                          (heo(i,k)+heo(i,k-1)))/factor
            endif
        enddo
      enddo
      do k = 2, km1
        do i = 1, im
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              xdby = hcko(i,k) - heso(i,k)
              xqrch = qeso(i,k) + gamma * xdby / (hvap * (1. + gamma))

              tem  = half * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = half * xlamud(i) * dz
              factor = one + tem - tem1
              qcko(i,k) = ((one-tem1)*qcko(i,k-1)+tem*half*(qo(i,k)+qo(i,k-1)))/factor
              dq = eta(i,k) * (qcko(i,k) - xqrch)
              if(k.ge.kbcon(i).and.dq.gt.0.) then
                etah = half * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0..and.k.gt.jmin(i)) then
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                if(k.lt.ktcon1(i)) then
                  xaa0(i) = xaa0(i) - dz * g * qlk
                endif
                qcko(i,k) = qlk + xqrch
                xpw = etah * c0 * dz * qlk
                xpwav(i) = xpwav(i) + xpw
              endif
            endif
            if(k.ge.kbcon(i).and.k.lt.ktcon1(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma  &
                      * to(i,k) / hvap
              xaa0(i) = xaa0(i)  &
                     + dz1 * (g / (cp * to(i,k)))  &
                     * xdby / (1. + gamma)  &
                     * rfact
              val=0.
              xaa0(i)=xaa0(i)+  &
                      dz1 * g * delta *  &
                      max(val,(qeso(i,k) - qo(i,k)))
            endif
        enddo
      enddo
!c
!c------- downdraft calculations
!c
!c--- downdraft moisture properties
!c
      do i = 1, im
          jmn = jmin(i)
          hcdo(i,jmn) = heo(i,jmn)
          qcdo(i,jmn) = qo(i,jmn)
          qrcd(i,jmn) = qeso(i,jmn)
          xpwev(i) = zero 
      enddo

      do k = km1, 1, -1
        do i = 1, im
          if (k.lt.jmin(i)) then
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = half * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = half * (xlamd(i)+xlamdd) * dz
              endif
              factor = one + tem - tem1
              hcdo(i,k) = ((one-tem1)*hcdo(i,k+1)+tem*half*  &
                          (heo(i,k)+heo(i,k+1)))/factor
          endif
        enddo
      enddo
!
      do k = km1, 1, -1
        do i = 1, im
          if (k .lt. jmin(i)) then
              dq = qeso(i,k)
              dt = to(i,k)
              gamma    = el2orc * dq / dt**2
              dh       = hcdo(i,k) - heso(i,k)
              qrcd(i,k)=dq+(1./hvap)*(gamma/(1.+gamma))*dh

              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = half * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = half * (xlamd(i)+xlamdd) * dz
              endif
              factor = one + tem - tem1
              qcdo(i,k) = ((1.-tem1)*qcdo(i,k+1)+tem*half*  &
                          (qo(i,k)+qo(i,k+1)))/factor
              xpwd     = etad(i,k+1) * (qcdo(i,k) - qrcd(i,k))
              qcdo(i,k)= qrcd(i,k)
              xpwev(i) = xpwev(i) + xpwd
          endif
        enddo
      enddo

      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i).eq.0.) edtmax = edtmaxs
          if(xpwev(i).ge.0.) then
            edtx(i) = 0.
          else
            edtx(i) = -edtx(i) * xpwav(i) / xpwev(i)
            edtx(i) = min(edtx(i),edtmax)
          endif
      enddo


!c--- downdraft cloudwork functions
!c
!c
      do k = km1, 1, -1
        do i = 1, im
          if ( k.lt.jmin(i)) then
              gamma = el2orc * qeso(i,k) / to(i,k)**2
              dhh=hcdo(i,k)
              dt= to(i,k)
              dg= gamma
              dh= heso(i,k)
              dz=-1.*(zo(i,k+1)-zo(i,k))
              xaa0(i)=xaa0(i)+edtx(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg))  &
                     *(one+delta*cp*dg*dt/hvap)
              val=zero
              xaa0(i)=xaa0(i)+edtx(i)*  &
             dz*g*delta*max(val,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo
!c
!c  calculate critical cloud work function
!c
      do i = 1, im
          if(pfld(i,ktcon(i)).lt.pcrit(15))then
            acrt(i)=acrit(15)*(975.0_r_kind-pfld(i,ktcon(i)))  &
                   /(975.0_r_kind-pcrit(15))
          else if(pfld(i,ktcon(i)).gt.pcrit(1))then
            acrt(i)=acrit(1)
          else
            k =  int((850.0_r_kind - pfld(i,ktcon(i)))/50.0_r_kind) + 2
            k = min(k,15)
            k = max(k,2)
            acrt(i)=acrit(k)+(acrit(k-1)-acrit(k))*  &
                (pfld(i,ktcon(i))-pcrit(k))/(pcrit(k-1)-pcrit(k))
          endif
      enddo

      do i = 1, im
          if(slimsk(i) .eq. one) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif

!  modify critical cloud workfunction by cloud base vertical velocity

          if(pdot(i).le.w4) then
            acrtfct(i) = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            acrtfct(i) = - (pdot(i) + w4) / (w4 - w3)
          else
            acrtfct(i) = zero
          endif
          val1    =   -one 
          acrtfct(i) = max(acrtfct(i),val1)
          val2    =    one
          acrtfct(i) = min(acrtfct(i),val2)
          acrtfct(i) = one - acrtfct(i)

          dtconv(i) = dt2 + max((1800.0_r_kind - dt2),zero) *  &
                     (pdot(i) - w2) / (w1 - w2)
          dtconv(i) = max(dtconv(i),dtmin)
          dtconv(i) = min(dtconv(i),dtmax)
      enddo

!--- large scale forcing

      do i= 1, im
          fld(i)=(aa1(i)-acrt(i)* acrtfct(i))/dtconv(i)
          if(fld(i) .le. zero) cnvscl_7(i) = tiny_r_kind

          cnvscl_7(i)=cnvscl_7(i)*cnvscl_6(i)
          xk(i) = (xaa0(i) - aa1(i)) / mbdt
          if(xk(i) .ge. zero) cnvscl_8(i) = tiny_r_kind 
          cnvscl_8(i)=cnvscl_8(i)*cnvscl_7(i)

          xmb(i) = -fld(i) / xk(i)
          xmb(i) = min(xmb(i),xmbmax(i))
          xmb(i) = xmb(i)*cnvscl_8(i)
      enddo


!  restore to,qo,uo,vo to t1,q1,u1,v1 in case convection stops

      do k = 1, km
        do i = 1, im
          if (k .le. kmax(i)) then
            to(i,k) = t1(i,k)
            qo(i,k) = q1(i,k)
            uo(i,k) = u1(i,k)
            vo(i,k) = v1(i,k)
            call fpvsx_ad(t1(i,k),es0,dum,dum,.false.)      ! fpvs is in pa
            qeso(i,k) = es0*10.0_r_kind
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val     =   1.e-8_r_kind
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo

!--- feedback: simply the changes from the cloud with unit mass flux
!---           multiplied by  the mass flux necessary to keep the
!---           equilibrium with the larger-scale.

      do i = 1, im
        delhbar(i) = zero
        delqbar(i) = zero
        deltbar(i) = zero
        delubar(i) = zero
        delvbar(i) = zero
        qcond(i) = zero
      enddo

      do k = 1, km
        do i = 1, im
          if (k .le. kmax(i)) then
            if(k.le.ktcon(i)) then
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
              q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
              tem = one/rcs(i)
              u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
              v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
              dp = 1000.0_r_kind * del(i,k)*ps(i)
              delhbar(i) = delhbar(i) + dellah(i,k)*xmb(i)*dp/g
              delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
              deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
              delubar(i) = delubar(i) + dellau(i,k)*xmb(i)*dp/g
              delvbar(i) = delvbar(i) + dellav(i,k)*xmb(i)*dp/g
            endif
          endif
        enddo
      enddo

      do k = 1, km
        do i = 1, im
          if ( k .le. kmax(i)) then
            if(k.le.ktcon(i)) then
              call fpvsx_ad(t1(i,k),es0,dum,dum,.false.)      ! fpvs is in pa
              qeso(i,k) = es0*10.0_r_kind
              qeso(i,k) = eps * qeso(i,k)/(pfld(i,k) + epsm1*qeso(i,k))
              val       = 1.e-8_r_kind
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo

      do i = 1, im
        rntot(i) = zero 
        delqev(i)= zero 
        delq2(i) = zero
      enddo

      do k = km, 1, -1
        do i = 1, im
          if (k .le. kmax(i)) then
            if(k.lt.ktcon(i)) then
              aup = one
              if(k.le.kb(i))   aup = zero
              adw = one
              if(k.ge.jmin(i)) adw = zero
              rain     =  aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
              rntot(i) = rntot(i) + rain * xmb(i) * .001_r_kind * dt2
            endif
          endif
        enddo
      enddo

      do k = km, 1, -1
        do i = 1, im
          if (k .le. kmax(i)) then
            deltv(i) = zero 
            delq(i) = zero
            qevap(i) = zero
            if(k.lt.ktcon(i)) then
              aup = one
              if(k.le.kb(i)) aup = zero
              adw = one
              if(k.ge.jmin(i)) adw = zero
              rain =  aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
              rn(i) = rn(i) + rain * xmb(i) * .001_r_kind * dt2
            endif
            if(k.lt.ktcon(i)) then
              evef = edt(i) * evfact
              if(slimsk(i) .eq. one) evef=edt(i) * evfactl
              qcond(i) = evef * (q1(i,k) - qeso(i,k))  &
                      / (one + el2orc * qeso(i,k) / t1(i,k)**2)
              dp = 1000.0_r_kind * del(i,k)*ps(i)
              if(rn(i).gt.zero .and. qcond(i) .lt. zero) then
                qevap(i) = -qcond(i) * (one-exp(-0.32_r_kind*sqrt(dt2*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.0_r_kind*g/dp)
                delq2(i) = delqev(i) + .001_r_kind * qevap(i) * dp / g
              endif
              if(rn(i) .gt. zero .and. qcond(i).lt.zero.and. &
                delq2(i).gt.rntot(i)) then
                qevap(i) = 1000.0_r_kind* g * (rntot(i) - delqev(i)) / dp
              endif
              if(rn(i).gt.zero.and.qevap(i).gt.zero) then
                q1(i,k) = q1(i,k) + qevap(i)
                t1(i,k) = t1(i,k) - elocp * qevap(i)
                rn(i) = rn(i) - .001_r_kind * qevap(i) * dp / g
                deltv(i) = - elocp*qevap(i)/dt2
                delq(i) =  + qevap(i)/dt2
                delqev(i) = delqev(i) + .001_r_kind*dp*qevap(i)/g
              endif
              dellaq(i,k)= dellaq(i,k) + delq(i) / xmb(i)
              delqbar(i) = delqbar(i) + delq(i)*dp/g
              deltbar(i) = deltbar(i) + deltv(i)*dp/g
            endif
          endif
        enddo
      enddo


!  precipitation rate converted to actual precip
!  in unit of m instead of kg

      do i = 1, im
          if(rn(i).lt. zero) rn(i) = zero
          if(rn(i).le. zero) then
            rn(i) = zero 
          else
            ktop(i) = ktcon(i)
            kbot(i) = kbcon(i)
            kcnv(i) = 1
            cldwrk(i) = aa1(i)
          endif
      enddo

!  cloud water

    if (ncloud.gt.0) then
      do k = 1, km
        do i = 1, im
          if (rn(i) .gt. zero) then
            if (k.gt.kb(i).and.k.le.ktcon(i)) then
              tem  = dellal(i,k) * xmb(i) * dt2
              tem1 = max(zero, min(one, (tcr-t1(i,k))*tcrf))
              if (ql(i,k,2) .gt. -999.0_r_kind) then
                ql(i,k,1) = ql(i,k,1) + tem * tem1            ! ice
                ql(i,k,2) = ql(i,k,2) + tem *(one-tem1)       ! water
              else
                ql(i,k,1) = ql(i,k,1) + tem
              endif
              if(ql(i,k,1) .lt. 1.0E-10_r_kind) ql(i,k,1) = zero 
            endif
          endif
        enddo
      enddo
    endif


      do k = 1, km
        do i = 1, im
              t_out(i,k) = t1(i,k)
              q_out(i,k) = q1(i,k)
              u_out(i,k) = u1(i,k)
              v_out(i,k) = v1(i,k)
              ql_out(i,k,1) = ql(i,k,1)
              ql_out(i,k,2) = ql(i,k,2)
        enddo
      enddo
endif  !if (ktcon .gt. 2

      return
      end



gfs_deepcnv_tl.f90

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


gfs_deepcnv_ad.f90

 subroutine sascnvn_ad(im, ix, km, jcap, delt, del, prsl, ps, sl,    &
                      ql,  q1, t1, u1, v1, &  
                      rcs, rn, jkbcon0,jkbcon1,jktcon0,jktcon1, slimsk, &
                      dot, ncloud,        &
                      t_out, q_out, u_out, v_out, ql_out, &
                      qld, q1d, t1d, u1d, v1d, dotd, &      !output
                      t_outd, q_outd, u_outd, v_outd, ql_outd, rnd)  !input

      USE kinds, only:r_single,r_kind,i_kind
      use CONSTANTS, only: CP, HVAP,  RD, FV, T0C, &
                     EPS, EPSM1, G=>grav, FACT1=>FACTOR1, FACT2=>FACTOR2,&
                    ONE, ZERO, HALF, tiny_r_kind, TWO, &
                    EL2ORC,elocp, init_constants_derived, init_constants


  IMPLICIT NONE

  REAL(r_kind) :: mbdt, tem, ptem, ptem1, pgcon
  REAL(r_kind) :: delt, clam, cxlamu, xlamde, xlamdd

  REAL(r_kind), DIMENSION(im) :: cnvscl_1, cnvscl_2, cnvscl_3, cnvscl_4, cnvscl_5,  &
                cnvscl_6,cnvscl_6_tmp1, cnvscl_7, cnvscl_7_tmp1, cnvscl_8,cnvscl_8_tmp1
  REAL(r_kind), DIMENSION(im) :: cnvscl_2d, cnvscl_5d, cnvscl_6d, cnvscl_7d, cnvscl_8d

  REAL(r_kind) :: ps(im), del(ix, km), prsl(ix, km), ql(ix, km, 2),ql_out(ix,km,2), q1(ix, km),     &
                  t1(ix, km), u1(ix, km), v1(ix, km), rcs(im), cldwrk(im), rn(im),  &
                  slimsk(im), dot(ix, km), sl(ix, km), t_out(im, km), q_out(im, km),&
                  u_out(im, km), v_out(im, km)

  REAL(r_kind),intent(out) :: qld(ix, km, 2), q1d(ix, km), t1d(ix, km), u1d(ix, km), &
                  v1d(ix, km), dotd(ix,km)
  REAL(r_kind) ::  rnd(im), t_outd(im, km), q_outd(im, km),  &
                  u_outd(im, km), v_outd(im, km),ql_outd(ix,km,2)

  REAL(r_kind) :: adw, aup, aafac, beta, betal, betas, dellat, desdt,            &
                  deta, detad, dg, dh, dhh, dlnsig, dp, dq, dqsdp, dqsdt,        &
                  dt, dt2, dtmax, dtmin, dv1h, dv1q, dv2h, dv2q, dv1u, dv1v,     & 
                  dv2u, dv2v, dv3q, dv3h, dv3u, dv3v, dz, dz1, e1, edtmax,       &
                  edtmaxl, edtmaxs, es(ix,km), es0, etah, evef, evfact, evfactl,        &
                  factor, fjcap, fkm, gamma, pprime, qlk, qrch, qs(im,km), rain, rfact, &
                  shear, tem1, tem2, val, val1, val2, val3,val4,w1, w1l, w1s, w2, w2l, w2s,&
                  w3, w3l, w3s, w4, w4l, w4s, xdby, xpw, xpwd, xqrch

  REAL(r_kind) :: dellatd, desdtd, dgd, dhd, dhhd, dqd, dqsdpd, dqsdtd, dtd,     &
                  dv1hd, dv1qd, dv2hd, dv2qd, dv1ud, dv1vd, dv2ud, dv2vd, dv3qd, &
                  dv3hd, dv3ud, dv3vd, dzd, dz1d, e1d, esd(im,km), es0d, etahd, evefd,  &
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
  REAL(r_kind) :: tod(im, km), qod(im, km), uod(im, km), vod(im, km), qesod(im, km),&
                  uo1d(im,km),vo1d(im,km),uo1(im,km),vo1(im,km)   

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
                  pwod(im, km), pwdod(im, km), sumxd(im), &
                  pdotd(im),acrtfctd(im),dtconvd(im)

  REAL(r_kind) :: pcrit(15), acritt(15), acrit(15)

  REAL(r_kind) :: min1, min1d
  REAL(r_kind) :: max1, max2, max3, max4, max5
  REAL(r_kind) :: max1d,max2d,max3d,max4d
  REAL(r_kind) :: y1, y2
  REAL(r_kind) :: y1d,y2d
  REAL(r_kind) :: tempd

  INTEGER(i_kind) :: im, ix, km, jcap, ncloud, kbot(im), ktop(im), kcnv(im)
  INTEGER(i_kind) :: i, j, indx, jmn, k, kk, latd, lond, km1
  INTEGER(i_kind) :: kb(im), kbcon(im), kbcon1(im), ktcon(im), ktcon1(im), ktcon_new(im), ktcon1_new(im),&
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


! Adjoint related variables
  REAL(r_kind), dimension(im,km) :: qeso_tmp1, qeso_tmp2, qeso_tmp3, qeso_tmp4, qeso_tmp5,qeso_tmp6,&
              qeso_tmp7, qeso_tmp8,qeso_tmp9,qeso_tmp10,qeso_tmp11,&
              qesod_tmp1, &
              qo_tmp1, qo_tmp2,qo_tmp3, qo_tmp4,qo_tmp5,qo_tmp6,qo_tmp7,&
              to_tmp1, to_tmp2, to_tmp3,to_tmp4,&
              uo_tmp1, &
              vo_tmp1, &
              es_tmp1,es_tmp2, es_tmp3, es_tmp4,es_tmp5,es_tmp6,es_tmp7,es_tmp8,&
              qs_tmp1,qs_tmp2,&
              xlamue_tmp1, xlamue_tmp2, xlamue_tmp3, xlamue_tmp4, &
              xlamue_tmp3d,&
              eta_tmp1, eta_tmp2, eta_tmp3, &
              ptem_tmp1, ptem_tmp2,ptem_tmp3, ptem_tmp4,ptem_tmp5,ptem_tmp6,&
              dz_tmp0,dz_tmp1,dz_tmp2, dz_tmp3,dz_tmp4, dz_tmp5,dz_tmp6,dz_tmp7,dz_tmp8, &
              dz_tmp9,dz_tmp10,dz_tmp11,dz_tmp12,dz_tmp13, dz_tmp14 ,dz_tmp15,dz_tmp16, &
              dz_tmp17,&
              dz1_tmp1,dz1_tmp2,&
              rfact_tmp1,rfact_tmp2,&
              tem_tmp1 , tem_tmp2, tem_tmp3,tem_tmp4,tem_tmp5,tem_tmp6, tem_tmp7, tem_tmp8, &
              tem_tmp9,tem_tmp10,tem_tmp11,&
              tem1_tmp1, tem1_tmp2, tem1_tmp3,tem1_tmp4,tem1_tmp5, tem1_tmp6,tem1_tmp7,tem1_tmp8, &
              tem1_tmp9,tem1_tmp10,tem1_tmp11,&
              ptem1_tmp1, ptem1_tmp2, ptem1_tmp3,& 
              gamma_tmp0,gamma_tmp1, gamma_tmp2,gamma_tmp3,gamma_tmp4,gamma_tmp5,gamma_tmp6,gamma_tmp7, &
              gamma_tmp8,gamma_tmp9, gamma_tmp10,gamma_tmp11,&
              qrch_tmp1, qrch_tmp2,qrch_tmp3,&
              factor_tmp1,factor_tmp2,factor_tmp3,factor_tmp4,factor_tmp5,factor_tmp6,factor_tmp7, &
              factor_tmp8, factor_tmp9,&
              hcko_tmp1, hcko_tmp2,  hcko_tmp3, &
              ucko_tmp1, ucko_tmp2, &
              vcko_tmp1, vcko_tmp2, &
              heo_tmp1, heo_tmp2, heo_tmp3, &
              heso_tmp1, heso_tmp2, &
              dq_tmp0,dq_tmp1, dq_tmp2, dq_tmp3, dq_tmp4, dq_tmp5,&
              dp_tmp0,dp_tmp1, dp_tmp2, &
              etah_tmp1, etah_tmp2,etah_tmp3,&
              qlk_tmp1, qlk_tmp2, qlk_tmp3,&
              dbyo_tmp1,dbyo_tmp2, &
              qcko_tmp1, qcko_tmp2, qcko_tmp3, qcko_tmp4,qcko_tmp5,qcko_tmp6,qcko_tmp7,&
              max1_tmp1, max3_tmp,max4_tmp,&
              etad_tmp1, &
              ucdo_tmp1, vcdo_tmp1, hcdo_tmp1, hcdo_tmp2,&
              qcdo_tmp1, qcdo_tmp2, qcdo_tmp3,qcdo_tmp4,&
              dh_tmp1,dh_tmp2, dh_tmp3, dhh_tmp1, dhh_tmp2, & 
              dt_tmp0,dt_tmp1, dg_tmp1, dg_tmp2, dt_tmp2, dt_tmp3, dt_tmp4,&
              dv1h_tmp1,dv2h_tmp1,dv3h_tmp1,&
              dv1q_tmp1,dv2q_tmp1,dv3q_tmp1,&
              dv1u_tmp1,dv2u_tmp1,dv3u_tmp1,&
              dv1v_tmp1,dv2v_tmp1,dv3v_tmp1, &
              qrcdo_tmp1,qrcdo_tmp2 , &
              dqsdt_tmp0, dqsdt_tmp1, desdt_tmp0,desdt_tmp1,dqsdp_tmp0, dqsdp_tmp1, &
              pprime_tmp0,pprime_tmp1,xdby_tmp1,xqrch_tmp1,qrcd_tmp1,dellat_tmp1, &
              t1_tmp1, t1_tmp2, &
              q1_tmp2, evef_tmp1,&
              rain_tmp1, rain_tmp2,&
              qevap_tmp1, qevap_tmp2,y1_tmp1,qcond_tmp1,rn_tmp1, &
              delq2_tmp1, rntot_tmp1

 REAL(r_kind), dimension(im) :: aa1_tmp1,aa1_tmp2, vshear_tmp1, vshear_tmp2, vshear_tmp2d, &
                edt_tmp1, edt_tmp2,pwevo_tmp1,edto_tmp1,edto_tmp2, edto_tmp3,pwavo_tmp1, &
                edtx_tmp1, edtx_tmp2, xpwav_tmp1, xpwev_tmp1,edtx_tmp3, &
                dtconv_tmp1, dtconv_tmp2, acrtfct_tmp1, acrtfct_tmp2,&
                xmb_tmp1, xmb_tmp2,rn_tmp2,rn_tmp3

 REAL(r_kind) :: amjk1 ,amjk1d,edto_tem1,edto_tem1d,pwevo_tem1d, &
                 dv1v_tmp2, dv1u_tmp2, dv1q_tmp2, dv1h_tmp2,edtx_tmpd,xmbtmp1d
 REAL(r_kind) :: ql_tmp1(im,km,2), ql_tmp2(im,km,2)
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

      do i=1,im
        cnvflg(i) = .true.
        rn(i)=zero
        kbot(i)=km+1
        ktop(i)=0
        kbcon(i)=km
        ktcon(i)=1
        dtconv(i) = 3600.0_r_kind
        cldwrk(i) = zero
        pdot(i) = zero
        pbcdif(i)= zero
        lmin(i) = 1
        jmin(i) = 1
        qlko_ktcon(i) = zero
        edt(i)  = zero
        edto(i) = zero
        edtx(i) = zero
        acrt(i) = zero
        acrtfct(i) = one
        aa1(i)  = zero
        aa2(i)  = zero
        xaa0(i) = zero
        pwavo(i)= zero
        pwevo(i)= zero
        xpwav(i)= zero
        xpwev(i)= zero
        vshear(i) = zero
        rntot(i) = zero
        delqev(i) = zero
        delq2(i) = zero
      enddo

      do k = 1, km
        do i = 1, im
          t_out(i,k) = t1(i,k)
          q_out(i,k) = q1(i,k)
          u_out(i,k) = u1(i,k)
          v_out(i,k) = v1(i,k)
          ql_out(i,k,1) = ql(i,k,1)
          ql_out(i,k,2) = ql(i,k,2)
        enddo
      enddo


      do k = 1, 15
        acrit(k) = acritt(k) * (975.0_r_kind - pcrit(k))
      enddo

      dt2 = delt
      val   =         1200.0_r_kind
      dtmin = max(dt2, val )
      val   =         3600.0_r_kind
      dtmax = max(dt2, val )
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

      pgcon   = 0.55_r_kind    ! Zhang & Wu (2003,JAS)
      fjcap   = (float(jcap) / 126.0_r_kind) ** 2
      val     = one 
      fjcap   = max(fjcap,val)
      fkm     = (float(km) / 28.0_r_kind) ** 2
      fkm     = max(fkm,val)
      w1l     = -8.e-3_r_kind
      w2l     = -4.e-2_r_kind
      w3l     = -5.e-3_r_kind
      w4l     = -5.e-4_r_kind
      w1s     = -2.e-4_r_kind
      w2s     = -2.e-3_r_kind
      w3s     = -1.e-3_r_kind
      w4s     = -2.e-5_r_kind

!  define top layer for search of the downdraft originating layer
!  and the maximum thetae for updraft
      do i=1,im
        kbmax(i) = km
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = one / ps(i)
      enddo

      do k = 1, km
        do i=1,im
          if (prsl(i,k)*tx1(i) .gt. 0.04_r_kind) kmax(i)  = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.45_r_kind) kbmax(i) = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.70_r_kind) kbm(i)   = k + 1
        enddo
      enddo
      do i=1,im
        kbmax(i) = min(kbmax(i),kmax(i))
        kbm(i)   = min(kbm(i),kmax(i))
      enddo

!   convert surface pressure to mb from cb

      do k = 1, km
        do i = 1, im
            to(i,k)   = t1(i,k)
            qo(i,k)   = q1(i,k)
            uo(i,k)   = u1(i,k) * rcs(i)
            vo(i,k)   = v1(i,k) * rcs(i)
            pfld(i,k) = prsl(i,k) * 10.0_r_kind
            eta(i,k)  = one
            fent1(i,k)= one
            fent2(i,k)= one
            frh(i,k)  = zero
            hcko(i,k) = zero
            qcko(i,k) = zero
            ucko(i,k) = zero
            vcko(i,k) = zero
            etad(i,k) = one
            hcdo(i,k) = zero
            qcdo(i,k) = zero
            ucdo(i,k) = zero
            vcdo(i,k) = zero
            qrcd(i,k) = zero
            qrcdo(i,k)= zero
            dbyo(i,k) = zero
            pwo(i,k)  = zero 
            pwdo(i,k) = zero
            dellal(i,k) = zero
        enddo
      enddo

!  column variables
!  p is pressure of the layer (mb)
!  t is temperature at t-dt (k)..tn
!  q is mixing ratio at t-dt (kg/kg)..qn
!  to is temperature at t+dt (k)... this is after advection and turbulan
!  qo is mixing ratio at t+dt (kg/kg)..q1

      do k = 1, km
        do i=1,im
          to_tmp1(i,k) = to(i,k)
            call fpvsx_ad(to(i,k),es0,dum, dum, .false.)      ! fpvs is in pa
            es(i,k) = es0*10.0_r_kind
              es_tmp1(i,k)=es(i,k)
            qeso(i,k) = eps * es(i,k) / (pfld(i,k) + epsm1*es(i,k))
               qeso_tmp1(i,k)= qeso(i,k)
            val1      =             1.e-8_r_kind
            qeso(i,k) = max(qeso(i,k), val1)
               qeso_tmp2(i,k)= qeso(i,k)
            val2      =           1.e-10_r_kind
               qo_tmp1(i,k) = qo(i,k)
            qo(i,k)   = max(qo(i,k), val2 )
               qo_tmp2(i,k)   = qo(i,k)
            tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
        enddo
      enddo

!  hydrostatic height assume zero terr and initially assume
!    updraft entrainment rate as an inverse function of height

      DO I = 1, IM   !kim
        ZO(I,1) = ZERO - log(SL(I,1))*RD/G * TVO(I,1)  !kim
      ENDDO   !kim

      do k = 2, km   !kim
        do i=1,im
        ZO(I,K) = ZO(I,K-1) -   &
         log(SL(I,K)/SL(I,K-1))*RD/G * HALF*(TVO(I,K) + TVO(I,K-1)) !kim
        enddo
      enddo

  DO k=1,km1
    DO i=1,im
      zi(i, k)  = half*(zo(i, k)+zo(i, k+1))
      xlamue(i, k)  =   clam/zi(i, k)
    END DO
  END DO

  DO k=1,km
    DO i=1,im
      xlamue_tmp1(i, k)  =   xlamue(i,k) 
    END DO
  END DO


!  compute moist static energy

  DO k=1,km
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        tem  = g*zo(i, k) + cp*to(i, k)
        heo(i, k)  = tem + hvap*qo(i, k)
        heso(i, k) = tem + hvap*qeso(i, k)
      END IF
    END DO
  END DO

  DO k=1,km
    DO i=1,im
      heo_tmp1(i,k) = heo(i,k)
      heso_tmp1(i,k) = heso(i,k)
    END DO
  END DO


!  determine level with largest moist static energy
!  this is the level where updraft starts

  DO i=1,im
    kb(i) = 1
  END DO


  do k =1, km1
    do i=1,im
       if (k .le. kmax(i)-1) then
          dz      = half * (zo(i,k+1) - zo(i,k))
             dz_tmp0(i,k) = dz
          dp      = half * (pfld(i,k+1) - pfld(i,k))
             dp_tmp0(i,k) = dp
          call fpvsx_ad(to(i,k+1),es0,dum, dum, .false.)      ! fpvs is in pa
          es(i,k) = es0 * 10.0_r_kind
          es_tmp2(i,k) = es(i,k)
          pprime  = pfld(i,k+1) + epsm1 * es(i,k)
             pprime_tmp0(i,k) = pprime
          qs(i,k)      = eps * es(i,k) / pprime
          qs_tmp1(i,k) = qs(i,k)
          dqsdp   = - qs(i,k) / pprime
              dqsdp_tmp0(i,k) = dqsdp
          desdt   = es(i,k) * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
              desdt_tmp0(i,k) = desdt
          dqsdt   = qs(i,k) * pfld(i,k+1) * desdt / (es(i,k) * pprime)
              dqsdt_tmp0(i,k) = dqsdt
          gamma   = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
              gamma_tmp0(i,k) = gamma
          dt      = (g * dz + hvap * dqsdp * dp) / (cp * (one + gamma))
              dt_tmp0(i,k) = dt
          dq      = dqsdt * dt + dqsdp * dp
              dq_tmp0(i,k) = dq
          to(i,k) = to(i,k+1) + dt
          qo(i,k) = qo(i,k+1) + dq
          po(i,k) = half  * (pfld(i,k) + pfld(i,k+1))
       endif
    enddo
  enddo

      do k = 1, km1
        do i=1,im
           to_tmp2(i,k) = to(i,k)
           qo_tmp3(i,k) = qo(i,k)
        enddo
      enddo

      do k = 1, km1
        do i=1,im
          if (k .le. kmax(i)-1) then
            call fpvsx_ad(to(i,k),es0,dum,dum,.false.)     
            es(i,k) = es0 * 10.0_r_kind
              es_tmp3(i,k) = es(i,k)
            qeso(i,k) = eps * es(i,k) / (po(i,k) + epsm1*es(i,k))
              qeso_tmp3(i,k) = qeso(i,k)
            val1      =     1.e-8_r_kind
            qeso(i,k) = max(qeso(i,k), val1)
              qeso_tmp4(i,k) = qeso(i,k)
            val2      =     1.e-10_r_kind
            qo(i,k)   = max(qo(i,k), val2 )
              qo_tmp4(i,k) = qo(i,k)

            IF (qo(i, k)/qeso(i, k) .GT. one) THEN
              min1  = one 
            ELSE
              min1 = qo(i, k)/qeso(i, k)
            END IF

            frh(i,k)  = one-min1

            heo(i,k)  = half * g * (zo(i,k) + zo(i,k+1)) +  &
                       cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = half * g * (zo(i,k) + zo(i,k+1)) +  &
                       cp * to(i,k) + hvap * qeso(i,k)
            uo1(i,k)   = half * (uo(i,k) + uo(i,k+1))
            vo1(i,k)   = half * (vo(i,k) + vo(i,k+1))
            uo(i,k)   = uo1(i,k) 
            vo(i,k)   = vo1(i,k) 
          endif
        enddo
      enddo

  DO k=1,km
    DO i=1,im
      heo_tmp2(i,k) = heo(i,k)
      heso_tmp2(i,k) = heso(i,k)
    END DO
  END DO

  DO k=1,km
    DO i=1,im
      uo_tmp1(i,k) = uo(i,k)
      vo_tmp1(i,k) = vo(i,k)
    END DO
  END DO


!  look for the level of free convection as cloud base

      do i=1,im
        flg(i)   = .true.
        cnvscl_1(i) = one
        cnvscl_2(i) = one
        cnvscl_3(i) = one
        cnvscl_4(i) = one
        cnvscl_5(i) = one
        cnvscl_6(i) = one
        cnvscl_7(i) = one
        cnvscl_8(i) = one
        kbcon(i) = kmax(i)
      enddo


      do k = 1, km1
        do i=1,im
          if (flg(i).and.k.le.kbmax(i)) then
            if(k.gt.kb(i).and.heo(i,kb(i)).gt.heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo

    do i = 1, im
      kbcon(i) = jkbcon0(i)
    enddo


      do i=1,im
        if(kbcon(i).eq.kmax(i)) cnvscl_1(i) = tiny_r_kind
      enddo

!  determine critical convective inhibition
!  as a function of vertical velocity at cloud base.

      do i=1,im
          pdot(i)  = 10.0_r_kind* dot(i,kbcon(i))
      enddo

      do i=1,im
          if(slimsk(i) .eq. one) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif

          if(pdot(i).le.w4) then
            tem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            tem = - (pdot(i) + w4) / (w4 - w3)
          else
            tem = zero
          endif
          val3    =             -one
          tem = max(tem,val3)
          val4    =             one
          tem = min(tem,val4)

          tem = one - tem
          tem1= half*(cincrmax-cincrmin)
          cincr = cincrmax - tem * tem1
          pbcdif(i) = pfld(i,kb(i)) - pfld(i,kbcon(i))
          if(pbcdif(i).gt.cincr) cnvscl_2(i) = tiny_r_kind
          cnvscl_2(i) = cnvscl_2(i)*cnvscl_1(i)
      enddo



!  assume that updraft entrainment rate above cloud base is
!    same as that at cloud base

      do k = 2, km1
        do i=1,im
          if((k.gt.kbcon(i).and.k.lt.kmax(i))) then
              xlamue(i,k) = xlamue(i,kbcon(i))
          endif
        enddo
      enddo

  DO k=1,km
    DO i=1,im
      xlamue_tmp2(i, k)  =   xlamue(i,k)
    END DO
  END DO

!  assume the detrainment rate for the updrafts to be same as
!  the entrainment rate at cloud base

      do i = 1, im
          xlamud(i) = xlamue(i,kbcon(i))
      enddo

!  functions rapidly decreasing with height, mimicking a cloud ensemble

      do k = 2, km1
        do i=1,im
          if((k.gt.kbcon(i).and.k.lt.kmax(i))) then
              tem = qeso(i,k)/qeso(i,kbcon(i))
              fent1(i,k) = tem**2
              fent2(i,k) = tem**3
          endif
        enddo
      enddo


!  final entrainment rate as the sum of turbulent part and organized entrainment
!    depending on the environmental relative humidity
!    (Bechtold et al., 2008)

      do k = 2, km1
        do i=1,im
          if((k.ge.kbcon(i).and.k.lt.kmax(i))) then
              tem = cxlamu * frh(i,k) * fent2(i,k)
              xlamue_tmp3(i,k) = xlamue(i,k)*fent1(i,k) + tem
          endif
        enddo
      enddo


      do k = 2, km1
        do i=1,im
          if((k.ge.kbcon(i).and.k.lt.kmax(i))) then
              xlamue(i,k) = xlamue_tmp3(i,k)
          endif
        enddo
      enddo

      do k = 1, km1
        do i=1,im
              xlamue_tmp4(i,k) = xlamue(i,k)
        enddo
      enddo


!  determine updraft mass flux for the subcloud layers

      do k = 1, km
        do i = 1, im
           eta_tmp1(i,k) = eta(i,k)
        enddo
      enddo

      do k = km1, 1, -1
        do i = 1, im
            if(k.lt.kbcon(i).and.k.ge.kb(i)) then
              dz       = zi(i,k+1) - zi(i,k)
                dz_tmp1(i,k) = dz
              ptem     = half*(xlamue(i,k)+xlamue(i,k+1))-xlamud(i)
                ptem_tmp1(i,k) = ptem 
              eta(i,k) = eta(i,k+1) / (one + ptem * dz)
            endif
        enddo
      enddo

      do k = 1, km
        do i = 1, im
           eta_tmp2(i,k) = eta(i,k)
        enddo
      enddo



!  compute mass flux above cloud base

      do k = 2, km1
        do i = 1, im
           if(k.gt.kbcon(i).and.k.lt.kmax(i)) then
              dz       = zi(i,k) - zi(i,k-1) 
               dz_tmp2(i,k) = dz
              ptem     = half*(xlamue(i,k)+xlamue(i,k-1))-xlamud(i)
               ptem_tmp2(i,k) = ptem
              eta(i,k) = eta(i,k-1) * (one + ptem * dz)
           endif
        enddo
      enddo
      do k = 1, km
        do i = 1, im
           eta_tmp3(i,k) = eta(i,k)
        enddo
      enddo


!  compute updraft cloud properties

      do i = 1, im
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
          pwavo(i)     = zero 
      enddo


!  cloud property is modified by the entrainment process

      do k = 2, km1
        do i = 1, im
            if(k.gt.kb(i).and.k.lt.kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
                dz_tmp3(i,k)= dz
              tem  = half * (xlamue(i,k)+xlamue(i,k-1)) * dz
                tem_tmp1(i,k) = tem
              tem1 = half * xlamud(i) * dz
                tem1_tmp1(i,k) = tem1
              factor = one + tem - tem1
                factor_tmp1(i,k) = factor
              ptem = half * tem + pgcon
                ptem_tmp3(i,k) = ptem 
              ptem1= half * tem - pgcon
                ptem1_tmp1(i,k) = ptem1 
              hcko(i,k) = ((one-tem1)*hcko(i,k-1)+tem*half*  &
                          (heo(i,k)+heo(i,k-1)))/factor
              ucko(i,k) = ((one-tem1)*ucko(i,k-1)+ptem*uo(i,k)  &
                          +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((one-tem1)*vcko(i,k-1)+ptem*vo(i,k)  &
                          +ptem1*vo(i,k-1))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)
            endif
        enddo
      enddo

      do k = 1, km
        do i = 1, im
           hcko_tmp2(i,k) = hcko(i,k)
           ucko_tmp2(i,k) = ucko(i,k)
           vcko_tmp2(i,k) = vcko(i,k)
           dbyo_tmp1(i,k) = dbyo(i,k) 
        enddo
      enddo


!   taking account into convection inhibition due to existence of
!    dry layers below cloud base

      do i=1,im
        kbcon1(i) = kmax(i)
        flg1(i) = .true.
      enddo

      do k = 2, km1
       do i=1,im
        if (k.lt.kmax(i)) then
          if(flg1(i) .and. k.ge.kbcon(i).and.dbyo(i,k).gt.zero) then
            kbcon1(i) = k
            flg1(i) = .false.
          endif
        endif
       enddo
      enddo

     do i = 1, im
        kbcon1(i) = jkbcon1(i)
     enddo

      do i=1,im
        if(kbcon1(i).eq.kmax(i)) cnvscl_2(i) = tiny_r_kind
      enddo

      do i=1,im
          tem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          cnvscl_3(i) = ONE - HALF**((dthk/tem)**200)
          cnvscl_3(i) = cnvscl_3(i)*cnvscl_2(i)
      enddo


!  determine first guess cloud top as the level of zero buoyancy

      do i = 1, im
        ktcon(i) = 1
        flg2(i) = .true.
      enddo

     do k = 2, km1
      do i = 1, im
        if (k .lt. kmax(i)) then
          if(flg2(i) .and. k.gt.kbcon1(i) .and. dbyo(i,k).lt. zero) then
             ktcon(i) = k
             flg2(i) = .false.
          endif
        endif
      enddo
      enddo

     do i = 1, im
        ktcon(i) = jktcon0(i)
     enddo

if(ktcon(1) .gt. 2) then

      do i = 1, im
          tem = pfld(i,kbcon(i))-pfld(i,ktcon(i))
          cnvscl_4(i) = ONE - HALF**((tem/cthk)**200)
          cnvscl_4(i) = cnvscl_4(i)*cnvscl_3(i)
      enddo



!  search for downdraft originating level above theta-e minimum

      do i = 1, im
           hmin(i) = heo(i,kbcon1(i))
           lmin(i) = kbmax(i)
           jmin(i) = kbmax(i)
      enddo

      do k = 2, km1
        do i = 1, im
          if (k .le. kbmax(i)) then
            if(k.gt.kbcon1(i).and.heo(i,k).lt.hmin(i)) then
               lmin(i) = k + 1
               hmin(i) = heo(i,k)
            endif
          endif
        enddo
      enddo

!  make sure that jmin(i) is within the cloud

      do i = 1, im
          jmin(i) = min(lmin(i),ktcon(i)-1)
          jmin(i) = max(jmin(i),kbcon1(i)+1)
          if(jmin(i).ge.ktcon(i)) cnvscl_4(i) = tiny_r_kind
      enddo

!  specify upper limit of mass flux at cloud base

      do i = 1, im
          k = kbcon(i)
          dp = 1000.0_r_kind * del(i,k)*ps(i)
          xmbmax(i) = dp / (g * dt2)
      enddo

!  compute cloud moisture property and precipitation

      do i = 1, im
          qcko(i,kb(i)) = qo(i,kb(i))
          qcko_tmp1(i,kb(i)) = qcko(i,kb(i))
      enddo

      do k = 2, km1
        do i = 1, im
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
                dz_tmp4(i,k) = dz
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
                gamma_tmp1(i,k) = gamma 
              qrch = qeso(i,k) + gamma * dbyo(i,k) / (hvap * (one + gamma))
                qrch_tmp1(i,k) = qrch 
              tem  = half * (xlamue(i,k)+xlamue(i,k-1)) * dz 
                tem_tmp2(i,k) = tem
              tem1 = half * xlamud(i) * dz
                tem1_tmp2(i,k) = tem1
              factor = one + tem - tem1
                factor_tmp2(i,k) = factor
              qcko(i,k) = ((one-tem1)*qcko(i,k-1)+tem*half*(qo(i,k)+qo(i,k-1)))/factor
           qcko_tmp1(i,k) = qcko(i,k)
              dq = eta(i,k) * (qcko(i,k) - qrch)
                dq_tmp1(i,k) = dq

!  check if there is excess moisture to release latent heat

              if(k.ge.kbcon(i).and.dq.gt.zero) then
                etah = half * (eta(i,k) + eta(i,k-1))
                  etah_tmp1(i,k) = etah
                if(ncloud .gt. 0  .and.k.gt.jmin(i)) then
                  dp = 1000.0_r_kind * del(i,k)*ps(i)
                    dp_tmp1(i,k) = dp
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                    qlk_tmp1(i,k) = qlk
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                    qlk_tmp1(i,k) = qlk
                endif
                aa1(i)  = aa1(i) - dz*g*qlk
                qcko(i, k) = qlk + qrch
                pwo(i, k)  = etah*c0*dz*qlk
                pwavo(i)   = pwavo(i) + pwo(i, k)
              endif
            endif
        enddo
      enddo

      do k = 1, km
        do i = 1, im
           qcko_tmp2(i,k) = qcko(i,k)
        enddo
      enddo


!  calculate cloud work function

      do k = 2, km1
        do i = 1, im
            if(k.ge.kbcon(i).and.k.lt.ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
                dz1_tmp1(i,k) = dz1
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
                gamma_tmp2(i,k) = gamma
              rfact =  one + delta * cp * gamma* to(i,k) / hvap
                rfact_tmp1(i,k) = rfact
              aa1(i) = aa1(i) + dz1 * (g / (cp * to(i,k))) &
                      * dbyo(i,k) / (one + gamma) * rfact
              val =zero 
              max1_tmp1(i,k) = max(val,(qeso(i,k) - qo(i,k)))
              aa1(i)=aa1(i)+ dz1 * g * delta * max(val,(qeso(i,k) - qo(i,k)))
            endif
        enddo
      enddo

      do i = 1, im
        aa1_tmp1(i) =  aa1(i) 
        cnvscl_5(i) = half*(one+tanh(aa1(i)))
        cnvscl_5(i) = cnvscl_5(i)*cnvscl_4(i)
      enddo


!  estimate the onvective overshooting as the level
!    where the [aafac * cloud work function] becomes zero,
!    which is the final cloud top

      do i = 1, im
        ktcon1(i) = jktcon1(i) 
      enddo


!  compute cloud moisture property, detraining cloud water
!    and precipitation in overshooting layers

      do k = 2, km1
        do i = 1, im
            if(k.ge.ktcon(i).and.k.lt.ktcon1(i)) then
              dz    = zi(i,k) - zi(i,k-1)
                dz_tmp5(i,k) = dz
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
                gamma_tmp3(i,k) = gamma
              qrch = qeso(i,k)+ gamma * dbyo(i,k) / (hvap * (one + gamma))
                qrch_tmp2(i,k) = qrch
              tem  = half * (xlamue(i,k)+xlamue(i,k-1)) * dz
                tem_tmp3(i,k) = tem
              tem1 = half * xlamud(i) * dz
                tem1_tmp3(i,k) = tem1
              factor = one + tem - tem1
                factor_tmp3(i,k) = factor
              qcko(i,k) = ((one-tem1)*qcko(i,k-1)+tem*half*    &
                          (qo(i,k)+qo(i,k-1)))/factor
                qcko_tmp3(i,k) = qcko(i,k)
              dq = eta(i,k) * (qcko(i,k) - qrch)
                dq_tmp2(i,k) = dq

!  check if there is excess moisture to release latent heat
              if(dq.gt.zero) then
                etah = half * (eta(i,k) + eta(i,k-1))
                   etah_tmp2(i,k) = etah
                if(ncloud .gt. 0) then
                  dp = 1000.0_r_kind * del(i,k)*ps(i)
                    dp_tmp2(i,k) = dp
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                    qlk_tmp2(i,k) = qlk
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                    qlk_tmp2(i,k) = qlk
                endif
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
                pwavo(i) = pwavo(i) + pwo(i,k)
              endif
            endif
        enddo
      enddo

      do i = 1, im
           pwavo_tmp1(i) = pwavo(i)
        do k = 1, km
           qcko_tmp4(i,k) = qcko(i,k)
        enddo
     enddo



! exchange ktcon with ktcon1

      do i = 1, im
          kk = ktcon(i)
          ktcon_new(i) = ktcon1(i)
          ktcon1_new(i) = kk
      enddo


!  this section is ready for cloud water

    if(ncloud .gt. 0) then
      do i = 1, im
          k = ktcon_new(i) - 1
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
               gamma_tmp4(i,k) = gamma
          qrch = qeso(i,k) + gamma * dbyo(i,k) / (hvap * (one + gamma))
               qrch_tmp3(i,k) = qrch
          dq = qcko(i,k) - qrch
               dq_tmp3(i,k) = dq
          if(dq .gt. zero) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
          endif
      enddo
    endif


      do i = 1, im
        do k = 1, km
           qcko_tmp5(i,k) = qcko(i,k)
        enddo
     enddo

!------- downdraft calculations
!--- compute precipitation efficiency in terms of windshear

      do i = 1, im
          vshear(i) = zero
      enddo

      do k = 2, km
        do i = 1, im
            if(k.gt.kb(i).and.k.le.ktcon_new(i)) then
              shear= sqrt((uo(i,k)-uo(i,k-1)) ** 2 &
                       + (vo(i,k)-vo(i,k-1)) ** 2)
              vshear(i) = vshear(i) + shear
              vshear_tmp1(i) = vshear(i)
            endif
        enddo
      enddo

      do i = 1, im
          vshear_tmp2(i) = 1.e3_r_kind * vshear(i) / (zi(i,ktcon_new(i))-zi(i,kb(i)))
          vshear(i) = vshear_tmp2(i)
          e1=1.591_r_kind-.639_r_kind*vshear(i)  &
            +.0953_r_kind*(vshear(i)**2)-.00496_r_kind*(vshear(i)**3)
          edt(i)=one-e1
            edt_tmp1(i) =edt(i)
          val =         .9_r_kind
          edt(i) = min(edt(i),val)
            edt_tmp2(i) =edt(i)
          val =         .0_r_kind
          edt(i) = max(edt(i),val)
          edto(i)=edt(i)
          edtx(i)=edt(i)
      enddo

!  determine detrainment rate between 1 and kbcon

      do i = 1, im
          sumx(i) = zero
      enddo

      do k = 1, km1
        do i = 1, im
          if(k.ge.1.and.k.lt.kbcon(i)) then
            dz = zi(i,k+1) - zi(i,k)
            sumx(i) = sumx(i) + dz
          endif
        enddo
      enddo


      do i = 1, im
        beta = betas
        if(slimsk(i).eq.1.) beta = betal
          dz  = (sumx(i)+zi(i,1))/float(kbcon(i))
            dz_tmp6(i,1) = dz
          tem = one/float(kbcon(i))
          xlamd(i) = (one-beta**tem)/dz
      enddo

!  determine downdraft mass flux

      do k = km1, 1, -1
        do i = 1, im
          if (k .le. kmax(i)-1) then
           if(k.lt.jmin(i).and.k.ge.kbcon(i)) then
              dz        = zi(i,k+1) - zi(i,k)
              dz_tmp7(i,k) = dz
              ptem      = xlamdd - xlamde
                ptem_tmp4(i,k) = ptem
              etad(i,k) = etad(i,k+1) * (one - ptem * dz)
           else if(k.lt.kbcon(i)) then
              dz        = zi(i,k+1) - zi(i,k)
              dz_tmp7(i,k) = dz
              ptem      = xlamd(i) + xlamdd - xlamde
                ptem_tmp4(i,k) = ptem
              etad(i,k) = etad(i,k+1) * (one - ptem * dz)
           endif
          endif
        enddo
      enddo

   do i = 1, im
     do k = 1, km
         etad_tmp1(i,k) =etad(i,k)
     enddo
   enddo



! downdraft moisture properties

      do i = 1, im
          jmn = jmin(i)
          hcdo(i,jmn) = heo(i,jmn)
          qcdo(i,jmn) = qo(i,jmn)
          qrcdo(i,jmn)= qeso(i,jmn)
          ucdo(i,jmn) = uo(i,jmn)
          vcdo(i,jmn) = vo(i,jmn)
      enddo

      do k = km1, 1, -1
        do i = 1, im
          if ( k.lt.jmin(i)) then
              dz = zi(i,k+1) - zi(i,k)
                  dz_tmp8(i,k) = dz
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                    tem_tmp4(i,k) = tem
                 tem1 = half * xlamdd * dz
                    tem1_tmp4(i,k) = tem1
              else
                 tem  = xlamde * dz
                    tem_tmp4(i,k) = tem
                 tem1 = half * (xlamd(i)+xlamdd) * dz
                    tem1_tmp4(i,k) = tem1
              endif
              factor = one + tem - tem1
                  factor_tmp4(i,k) = factor
              ptem = half * tem - pgcon
                  ptem_tmp5(i,k) = ptem
              ptem1= half * tem + pgcon
                  ptem1_tmp2(i,k) = ptem1
              hcdo(i,k) = ((one-tem1)*hcdo(i,k+1)+tem*half*   &
                          (heo(i,k)+heo(i,k+1)))/factor
              ucdo(i,k) = ((one-tem1)*ucdo(i,k+1)+ptem*uo(i,k+1) &
                          +ptem1*uo(i,k))/factor
              vcdo(i,k) = ((one-tem1)*vcdo(i,k+1)+ptem*vo(i,k+1)  &
                          +ptem1*vo(i,k))/factor
              dbyo(i,k) = hcdo(i,k) - heso(i,k)
          endif
        enddo
      enddo


     do i = 1, im
          do k = 1, km
            hcdo_tmp1(i,k) = hcdo(i,k)
            ucdo_tmp1(i,k) = ucdo(i,k)
            vcdo_tmp1(i,k) = vcdo(i,k)
            dbyo_tmp2(i,k) = dbyo(i,k)
          enddo
     enddo


      do k = km1, 1, -1
        do i = 1, im
          if (k.lt.jmin(i)) then
              gamma      = el2orc * qeso(i,k) / (to(i,k)**2)
              gamma_tmp5(i,k) = gamma
              qrcdo(i,k) = qeso(i,k)+ (one/hvap)*(gamma/(one+gamma))*dbyo(i,k)
                qrcdo_tmp1(i,k) = qrcdo(i,k)
              dz = zi(i,k+1) - zi(i,k)
              dz_tmp9(i,k) = dz
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = half * xlamdd * dz
                    tem_tmp5(i,k)= tem
                    tem1_tmp5(i,k)= tem1
              else
                 tem  = xlamde * dz
                 tem1 = half * (xlamd(i)+xlamdd) * dz
                    tem_tmp5(i,k)= tem
                    tem1_tmp5(i,k)= tem1
              endif
              factor = one + tem - tem1
              factor_tmp5(i,k) = factor
              qcdo(i,k) = ((one-tem1)*qcdo(i,k+1)+tem*half*  &
                          (qo(i,k)+qo(i,k+1)))/factor
              qcdo_tmp1(i,k) = qcdo(i,k)
              pwdo(i,k)  = etad(i,k+1) * (qcdo(i,k) - qrcdo(i,k))

              qcdo(i,k)  = qrcdo(i,k)
              pwevo(i)   = pwevo(i) + pwdo(i,k)
          endif
        enddo
      enddo

     do i = 1, im
       pwevo_tmp1(i) = pwevo(i)
          do k = 1, km
            qcdo_tmp2(i,k) = qcdo(i,k)
            qrcdo_tmp2(i,k) = qrcdo(i,k)
          enddo
     enddo


!--- final downdraft strength dependent on precip
!--- efficiency (edt), normalized condensate (pwav), and
!--- evaporate (pwev)

      do i = 1, im
        edtmax = edtmaxl
        edto_tmp1(i) = edto(i)
          if(pwevo(i).lt. zero) then
            edto(i) = -edto(i) * pwavo(i) / pwevo(i)
            edto_tmp2(i) = edto(i)
            edto(i) = min(edto(i),edtmax)
          else
            edto(i) = zero
          endif
            edto_tmp3(i) = edto(i)
      enddo


!--- downdraft cloudwork functions

      do k = km1, 1, -1
        do i = 1, im
          if (k .lt. jmin(i)) then
              gamma = el2orc * qeso(i,k) / to(i,k)**2
                gamma_tmp6(i,k) = gamma
              dhh=hcdo(i,k)
                dhh_tmp1(i,k) = dhh 
              dt=to(i,k)
                dt_tmp1(i,k) = dt 
              dg=gamma
                dg_tmp1(i,k) = dg 
              dh=heso(i,k)
                dh_tmp1(i,k) = dh 
              dz=-1.*(zo(i,k+1)-zo(i,k))
                dz_tmp10(i,k) = dz

              aa1(i)=aa1(i)+edto(i)*dz*(g/(cp*dt))*((dhh-dh)/(one+dg)) &
                    *(one+delta*cp*dg*dt/hvap)

              val=0.
              IF (val .LT. qeso(i, k) - qo(i, k)) THEN
                 max2 = qeso(i, k) - qo(i, k)
              ELSE
                max2 = val
              END IF

              aa1(i)=aa1(i)+edto(i)*dz*g*delta*max2
          endif
        enddo
      enddo

      do i = 1, im
         aa1_tmp2(i) = aa1(i)
         cnvscl_6(i) = half*(one+tanh(aa1(i)))
         cnvscl_6_tmp1(i) = cnvscl_6(i)
         cnvscl_6(i) = cnvscl_6(i)*cnvscl_5(i)
      enddo



!--- what would the change be, that a cloud with unit mass
!--- will do to the environment?

      do k = 1, km
        do i = 1, im
            dellah(i,k) = zero
            dellaq(i,k) = zero
            dellau(i,k) = zero
            dellav(i,k) = zero
        enddo
      enddo

      do i = 1, im
          dp = 1000.0_r_kind * del(i,1)*ps(i)
          dellah(i,1) = edto(i) * etad(i,1) * (hcdo(i,1)  &
                        - heo(i,1)) * g / dp
          dellaq(i,1) = edto(i) * etad(i,1) * (qcdo(i,1)  &
                        - qo(i,1)) * g / dp
          dellau(i,1) = edto(i) * etad(i,1) * (ucdo(i,1)  &
                        - uo(i,1)) * g / dp
          dellav(i,1) = edto(i) * etad(i,1) * (vcdo(i,1) &
                        - vo(i,1)) * g / dp
      enddo

!--- changed due to subsidence and entrainment

      do k = 2, km1
        do i = 1, im
          if (k.lt.ktcon_new(i)) then
              aup = one
              if(k.le.kb(i)) aup = zero
              adw = one
              if(k.gt.jmin(i)) adw = zero
              dp = 1000.0_r_kind * del(i,k)*ps(i)
              dz = zi(i,k) - zi(i,k-1)
              dz_tmp11(i,k) = dz

              dv1h = heo(i,k)
              dv2h = half * (heo(i,k) + heo(i,k-1))
              dv3h = heo(i,k-1)
              dv1q = qo(i,k)
              dv2q = half * (qo(i,k) + qo(i,k-1))
              dv3q = qo(i,k-1)
              dv1u = uo(i,k)
              dv2u = half * (uo(i,k) + uo(i,k-1))
              dv3u = uo(i,k-1)
              dv1v = vo(i,k)
              dv2v = half * (vo(i,k) + vo(i,k-1))
              dv3v = vo(i,k-1)

              dv1h_tmp1(i,k) = dv1h
              dv2h_tmp1(i,k) = dv2h
              dv3h_tmp1(i,k) = dv3h
              dv1q_tmp1(i,k) = dv1q
              dv2q_tmp1(i,k) = dv2q
              dv3q_tmp1(i,k) = dv3q
              dv1u_tmp1(i,k) = dv1u
              dv2u_tmp1(i,k) = dv2u
              dv3u_tmp1(i,k) = dv3u
              dv1v_tmp1(i,k) = dv1v
              dv2v_tmp1(i,k) = dv2v
              dv3v_tmp1(i,k) = dv3v

              tem  = half * (xlamue(i,k)+xlamue(i,k-1))
              tem1 = xlamud(i)

              tem_tmp6(i,k) = tem
              tem1_tmp6(i,k) = tem1

        IF (k .LE. kbcon(i)) THEN
          ptem   = xlamde
          ptem1  = xlamd(i) + xlamdd
        ELSE
          ptem   = xlamde
          ptem1  = xlamdd
        END IF

        ptem1_tmp3(i,k) = ptem1
        ptem_tmp6(i,k) = ptem

        dellah(i, k)  = dellah(i, k) + ((aup*eta(i, k)-adw*edto(i)*etad(i, k))*&
                        dv1h-(aup*eta(i, k-1)-adw*edto(i)*etad(i, k-1))*dv3h-(&
                        aup*tem*eta(i, k-1)+adw*edto(i)*ptem*etad(i, k))*dv2h*dz+aup*&
                        tem1*eta(i, k-1)*half*(hcko(i, k)+hcko(i, k-1))*dz+adw*edto(i)*&
                        ptem1*etad(i, k)*half*(hcdo(i, k)+hcdo(i, k-1))*dz)*g/dp

        dellaq(i, k)  = dellaq(i, k) + ((aup*eta(i, k)-adw*edto(i)*etad(i, k))*dv1q-&
                        (aup*eta(i, k-1)-adw*edto(i)*etad(i, k-1))*dv3q-(aup*tem*&
                        eta(i, k-1)+adw*edto(i)*ptem*etad(i, k))*dv2q*dz+aup*tem1*&
                        eta(i, k-1)*half*(qcko(i, k)+qcko(i, k-1))*dz+adw*edto(i)*ptem1*&
                        etad(i, k)*half*(qrcdo(i, k)+qrcdo(i, k-1))*dz)*g/dp
        dellau(i, k)  = dellau(i, k) + ((aup*eta(i, k)-adw*edto(i)*etad(i, k))*&
                        dv1u-(aup*eta(i, k-1)-adw*edto(i)*etad(i, k-1))*dv3u-(aup*&
                        tem*eta(i, k-1)+adw*edto(i)*ptem*etad(i, k))*dv2u*dz+aup*&
                        tem1*eta(i, k-1)*half*(ucko(i, k)+ucko(i, k-1))*dz+adw*edto(i)*&
                        ptem1*etad(i, k)*half*(ucdo(i, k)+ucdo(i, k-1))*dz-pgcon*(aup*&
                        eta(i, k-1)-adw*edto(i)*etad(i, k))*(dv1u-dv3u))*g/dp
        dellav(i, k)  = dellav(i, k) + ((aup*eta(i, k)-adw*edto(i)*etad(i, k))*dv1v-&
                        (aup*eta(i, k-1)-adw*edto(i)*etad(i, k-1))*dv3v-(aup*tem*&
                        eta(i, k-1)+adw*edto(i)*ptem*etad(i, k))*dv2v*dz+aup*tem1*&
                        eta(i, k-1)*half*(vcko(i, k)+vcko(i, k-1))*dz+adw*edto(i)*ptem1*&
                        etad(i, k)*half*(vcdo(i, k)+vcdo(i, k-1))*dz-pgcon*(aup*&
                        eta(i, k-1)-adw*edto(i)*etad(i, k))*(dv1v-dv3v))*g/dp
          endif
        enddo
      enddo


!------- cloud top

      do i = 1, im
          indx = ktcon_new(i)
          dp = 1000.0_r_kind * del(i,indx)*ps(i)
          dv1h = heo(i,indx-1)
          dv1h_tmp2 = heo(i,indx-1)

          dellah(i,indx) = eta(i,indx-1) *   &
                          (hcko(i,indx-1) - dv1h) * g / dp
          dv1q = qo(i,indx-1)
          dv1q_tmp2 = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *   &
                          (qcko(i,indx-1) - dv1q) * g / dp
          dv1u = uo(i,indx-1)
          dv1u_tmp2 = uo(i,indx-1)
          dellau(i,indx) = eta(i,indx-1) *    &
                          (ucko(i,indx-1) - dv1u) * g / dp
          dv1v = vo(i,indx-1)
          dv1v_tmp2 = vo(i,indx-1)
          dellav(i,indx) = eta(i,indx-1) *    &
                          (vcko(i,indx-1) - dv1v) * g / dp

          dellal(i,indx) = eta(i,indx-1) * qlko_ktcon(i) * g / dp
      enddo


!------- final changed variable per unit mass flux

      do k = 1, km
        do i = 1, im
          if (k .le. kmax(i)) then
            if(k.gt.ktcon_new(i)) then
              qo(i,k) = q1(i,k)
              to(i,k) = t1(i,k)
            endif
            if(k.le.ktcon_new(i)) then
              qo(i,k) = dellaq(i,k) * mbdt + q1(i,k)
                qo_tmp5(i,k) = qo(i,k)
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              to(i,k) = dellat * mbdt + t1(i,k)
              val   =    1.e-10_r_kind
              qo(i,k) = max(qo(i,k), val  )
            endif
          endif
        enddo
      enddo

   do i = 1, im
       do k = 1, km
           to_tmp3(i,k) = to(i,k)
       enddo
  enddo


!--- the above changed environment is now used to calulate the
!--- effect the arbitrary cloud (with unit mass flux)
!--- would have on the stability,
!--- which then is used to calculate the real mass flux,
!--- necessary to keep this change in balance with the large-scale
!--- destabilization.
!--- environmental conditions again, first heights

      do k = 1, km
        do i = 1, im
          if(k .le. kmax(i)) then
            call fpvsx_ad(to(i,k),es0,dum,dum,.false.) 
            es(i,k) = es0*10.0_r_kind
               es_tmp4(i,k) = es(i,k)
            qeso(i,k) = eps * es(i,k) / (pfld(i,k)+epsm1*es(i,k))
               qeso_tmp5(i,k) = qeso(i,k)
            val       =   1.e-8_r_kind
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
      do i = 1, im
        do k = 1, km
           qeso_tmp6(i,k) = qeso(i,k)
        enddo
      enddo



!--- moist static energy

      do k = 1, km1
        do i = 1, im
          if(k .le. kmax(i)-1) then
            dz  = half * (zo(i,k+1) - zo(i,k))
             dz_tmp12(i,k) = dz
            dp = half * (pfld(i,k+1) - pfld(i,k))
            call fpvsx_ad(to(i,k+1),es0,dum,dum,.false.)     
            es(i,k) = es0*10.0_r_kind
              es_tmp5(i,k) = es(i,k)
            pprime = pfld(i,k+1) + epsm1 * es(i,k)
              pprime_tmp1(i,k) = pprime
            qs(i,k) = eps * es(i,k) / pprime
              qs_tmp2(i,k) = qs(i,k)
            dqsdp = - qs(i,k) / pprime
              dqsdp_tmp1(i,k) = dqsdp
            desdt = es(i,k) * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
              desdt_tmp1(i,k) = desdt
            dqsdt = qs(i,k) * pfld(i,k+1) * desdt / (es(i,k) * pprime)
              dqsdt_tmp1(i,k) = dqsdt
            gamma = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
              gamma_tmp7(i,k) = gamma
            dt = (g * dz + hvap * dqsdp * dp) / (cp * (one + gamma))
              dt_tmp2(i,k) = dt
            dq = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = half * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo


      do i = 1, im
         do k = 1, km
            to_tmp4(i,k) = to(i,k)
            qo_tmp6(i,k) = qo(i,k)
         enddo
      enddo

      do k = 1, km1
        do i = 1, im
          if(k .le. kmax(i)-1) then
            call fpvsx_ad(to(i,k),es0,dum,dum,.false.)     
            es(i,k) = es0*10.0_r_kind
              es_tmp6(i,k) = es(i,k)
            qeso(i,k) = eps * es(i,k) / (po(i,k) + epsm1 * es(i,k))
              qeso_tmp7(i,k) = qeso(i,k)
            val1      =  1.e-8_r_kind
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =  1.e-10_r_kind
            qo(i,k)   = max(qo(i,k), val2 )
            heo(i,k)   = half * g * (zo(i,k) + zo(i,k+1)) + &
                         cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = half * g * (zo(i,k) + zo(i,k+1)) +  &
                       cp * to(i,k) + hvap * qeso(i,k)
          endif
        enddo
      enddo

      do i = 1, im
          k = kmax(i)
          heo(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qo(i,k)
          heso(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qeso(i,k)
      enddo

      do i = 1, im
         do k = 1, km
            heo_tmp3(i,k) = heo(i,k)
            qeso_tmp8(i,k) = qeso(i,k)
            qo_tmp7(i,k) = qo(i,k)
         enddo
     enddo


!**************************** static control
!------- moisture and cloud work functions

      do i = 1, im
          xaa0(i) = zero
          xpwav(i) = zero
      enddo

      do i = 1, im
          indx = kb(i)
          hcko(i,indx) = heo(i,indx)
          qcko(i,indx) = qo(i,indx)
      enddo

      do k = 2, km1
        do i = 1, im
            if(k.gt.kb(i).and.k.le.ktcon_new(i)) then
              dz = zi(i,k) - zi(i,k-1)
                 dz_tmp13(i,k) = dz
              tem  = half * (xlamue(i,k)+xlamue(i,k-1)) * dz
                 tem_tmp7(i,k) = tem
              tem1 = half * xlamud(i) * dz
                 tem1_tmp7(i,k) = tem1
              factor = one + tem - tem1
                 factor_tmp6(i,k) = factor
              hcko(i,k) = ((one-tem1)*hcko(i,k-1)+tem*half*   &
                          (heo(i,k)+heo(i,k-1)))/factor
            endif
        enddo
      enddo

      do i = 1, im
        do k = 1, km
            hcko_tmp3(i,k) = hcko(i,k)
        enddo
      enddo

      do k = 2, km1
        do i = 1, im
            if(k.gt.kb(i).and.k.lt.ktcon_new(i)) then
              dz = zi(i,k) - zi(i,k-1)
                 dz_tmp14(i,k) = dz
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
                 gamma_tmp8(i,k) = gamma
              xdby = hcko(i,k) - heso(i,k)
                 xdby_tmp1(i,k) = xdby
              xqrch = qeso(i,k) + gamma * xdby / (hvap * (one + gamma))
                 xqrch_tmp1(i,k) = xqrch

              tem  = half * (xlamue(i,k)+xlamue(i,k-1)) * dz
                tem_tmp8(i,k) = tem
              tem1 = half * xlamud(i) * dz
                tem1_tmp8(i,k) = tem1
              factor = one + tem - tem1
                factor_tmp7(i,k) = factor
              qcko(i,k) = ((one-tem1)*qcko(i,k-1)+tem*half*(qo(i,k)+qo(i,k-1)))/factor
                qcko_tmp6(i,k) = qcko(i,k)
              dq = eta(i,k) * (qcko(i,k) - xqrch)
                dq_tmp4(i,k) = dq
              if(k.ge.kbcon(i) .and. dq.gt.zero) then
                etah = half * (eta(i,k) + eta(i,k-1))
                   etah_tmp3(i,k) = etah
                if(ncloud.gt.zero.and.k.gt.jmin(i)) then
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                  qlk_tmp3(i,k) = qlk
               if(k.lt.ktcon1_new(i)) then
                  xaa0(i) = xaa0(i) - dz * g * qlk
                endif
                qcko(i,k) = qlk + xqrch
                xpw = etah * c0 * dz * qlk
                xpwav(i) = xpwav(i) + xpw
              endif
            endif

            if(k.ge.kbcon(i).and.k.lt.ktcon1_new(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
                 dz1_tmp2(i,k) = dz1
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
                 gamma_tmp9(i,k) = gamma
              rfact =  one + delta * cp * gamma  * to(i,k) / hvap
                 rfact_tmp2(i,k) = rfact
              xaa0(i) = xaa0(i)  &
                     + dz1 * (g / (cp * to(i,k)))  &
                     * xdby / (one + gamma)  &
                     * rfact

              val=zero
              IF (val .LT. qeso(i, k) - qo(i, k)) THEN
                max3  = qeso(i, k) - qo(i, k)
              ELSE
                max3 = val
              END IF
              max3_tmp(i,k) = max3
              xaa0(i)=xaa0(i)+ dz1 * g * delta *max3

            endif
        enddo
      enddo

      do i = 1, im
         do k = 1, km
           qcko_tmp7(i,k) = qcko(i,k)
         enddo
      enddo



!------- downdraft calculations
!--- downdraft moisture properties

      do i = 1, im
          jmn = jmin(i)
          hcdo(i,jmn) = heo(i,jmn)
          qcdo(i,jmn) = qo(i,jmn)
          qrcd(i,jmn) = qeso(i,jmn)
      enddo

      do k = km1, 1, -1
        do i = 1, im
          if (k.lt.jmin(i)) then
              dz = zi(i,k+1) - zi(i,k)
                dz_tmp15(i,k) = dz
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = half * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = half * (xlamd(i)+xlamdd) * dz
              endif
              tem_tmp9(i,k) = tem
              tem1_tmp9(i,k) = tem1

              factor = one + tem - tem1 
                factor_tmp8(i,k) = factor
             hcdo(i,k) = ((one-tem1)*hcdo(i,k+1)+tem*half*  &
                          (heo(i,k)+heo(i,k+1)))/factor

          endif
        enddo
      enddo

      do i = 1, im
        do k = 1, km
           hcdo_tmp2(i,k) = hcdo(i,k)
        enddo
      enddo

      do k = km1, 1, -1
        do i = 1, im
          if (k .lt. jmin(i)) then
              dq = qeso(i,k)
                dq_tmp5(i,k) = dq
              dt = to(i,k)
                dt_tmp3(i,k) = dt
              gamma    = el2orc * dq / dt**2
                gamma_tmp10(i,k) = gamma
              dh       = hcdo(i,k) - heso(i,k)
                dh_tmp2(i,k) = dh
              qrcd(i,k)=dq+(one/hvap)*(gamma/(one+gamma))*dh
                qrcd_tmp1(i,k) = qrcd(i,k)

              dz = zi(i,k+1) - zi(i,k)
                dz_tmp16(i,k) = dz

              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = half * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = half * (xlamd(i)+xlamdd) * dz
              endif
                  tem_tmp10(i,k) = tem
                  tem1_tmp10(i,k) = tem1

              factor = one + tem - tem1
                  factor_tmp9(i,k) = factor

              qcdo(i,k) = ((one-tem1)*qcdo(i,k+1)+tem*half*  &
                          (qo(i,k)+qo(i,k+1)))/factor
                  qcdo_tmp3(i,k) = qcdo(i,k)
              xpwd     = etad(i,k+1) * (qcdo(i,k) - qrcd(i,k))
              qcdo(i,k)= qrcd(i,k)
              xpwev(i) = xpwev(i) + xpwd
          endif
        enddo
      enddo

       do i = 1, im
         xpwev_tmp1(i) = xpwev(i)
         xpwav_tmp1(i) = xpwav(i)
       enddo

      do i = 1, im
        edtx_tmp1(i) = edtx(i)
        edtmax = edtmaxl
        if(slimsk(i).eq.zero) edtmax = edtmaxs
          if(xpwev(i).ge.zero) then
            edtx(i) = zero 
          else
            edtx(i) = -edtx(i) * xpwav(i) / xpwev(i)
            edtx_tmp2(i) = edtx(i)
            edtx(i) = min(edtx(i),edtmax)
          endif
      enddo

     do i = 1, im
        edtx_tmp3(i) = edtx(i)
        do k = 1, km
          qcdo_tmp4(i,k) = qcdo(i,k)
        enddo
    enddo

!--- downdraft cloudwork functions


      do k = km1, 1, -1
        do i = 1, im
          if ( k.lt.jmin(i)) then
              gamma = el2orc * qeso(i,k) / to(i,k)**2
                 gamma_tmp11(i,k) = gamma
              dhh=hcdo(i,k)
                 dhh_tmp2(i,k) = dhh
              dt= to(i,k)
                 dt_tmp4(i,k) = dt
              dg= gamma
                 dg_tmp2(i,k) = dg
              dh= heso(i,k)
                 dh_tmp3(i,k) = dh
              dz=-1.*(zo(i,k+1)-zo(i,k))
                 dz_tmp17(i,k) = dz
          xaa0(i)=xaa0(i)+edtx(i)*dz*(g/(cp*dt))*((dhh-dh)/(one+dg))  &
                     *(one+delta*cp*dg*dt/hvap)

              val=zero
              IF (val .LT. qeso(i, k) - qo(i, k)) THEN
                 max4  = qeso(i, k) - qo(i, k)
              ELSE
                 max4  = val
              END IF
              max4_tmp(i,k) = max4
              xaa0(i)=xaa0(i)+edtx(i)* dz*g*delta*max4
          endif
        enddo
      enddo

!  calculate critical cloud work function

      do i = 1, im
          if(pfld(i,ktcon_new(i)).lt.pcrit(15))then
            acrt(i)=acrit(15)*(975.0_r_kind-pfld(i,ktcon_new(i)))  &
                   /(975.0_r_kind-pcrit(15))
          else if(pfld(i,ktcon_new(i)).gt.pcrit(1))then
            acrt(i)=acrit(1)
          else
            k =  int((850.0_r_kind - pfld(i,ktcon_new(i)))/50.0_r_kind) + 2
            k = min(k,15)
            k = max(k,2)
            acrt(i)=acrit(k)+(acrit(k-1)-acrit(k))*  &
                (pfld(i,ktcon_new(i))-pcrit(k))/(pcrit(k-1)-pcrit(k))
          endif
      enddo


      do i = 1, im
          if(slimsk(i).eq.one) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif

!  modify critical cloud workfunction by cloud base vertical velocity

          if(pdot(i).le.w4) then
            acrtfct(i) = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            acrtfct(i) = - (pdot(i) + w4) / (w4 - w3)
          else
            acrtfct(i) = zero
          endif
            acrtfct_tmp1(i) = acrtfct(i)
          val1    =             -one
          acrtfct(i) = max(acrtfct(i),val1)
            acrtfct_tmp2(i) = acrtfct(i)
          val2    =             one
          acrtfct(i) = min(acrtfct(i),val2)
          acrtfct(i) = one - acrtfct(i)

          dtconv(i) = dt2 + max((1800.0_r_kind - dt2),zero) *  &
                     (pdot(i) - w2) / (w1 - w2)
            dtconv_tmp1(i) = dtconv(i)
          dtconv(i) = max(dtconv(i),dtmin)
            dtconv_tmp2(i) = dtconv(i)
          dtconv(i) = min(dtconv(i),dtmax)
      enddo


!--- large scale forcing

      do i= 1, im
          fld(i)=(aa1(i)-acrt(i)* acrtfct(i))/dtconv(i)
          if(fld(i) .le. zero) cnvscl_7(i) = tiny_r_kind
          cnvscl_7_tmp1(i) = cnvscl_7(i)
          cnvscl_7(i)=cnvscl_7(i)*cnvscl_6(i)
          xk(i) = (xaa0(i) - aa1(i)) / mbdt
          if(xk(i) .ge. zero) cnvscl_8(i) = tiny_r_kind
          cnvscl_8_tmp1(i) = cnvscl_8(i)
          cnvscl_8(i)=cnvscl_8(i)*cnvscl_7(i)

!--- kernel, cloud base mass flux

          xmb(i) = -fld(i) / xk(i)
          xmb_tmp1(i) = xmb(i)
          xmb(i) = min(xmb(i),xmbmax(i))
          xmb_tmp2(i) = xmb(i)
          xmb(i) = xmb(i)*cnvscl_8(i)
      enddo


!  restore to,qo,uo,vo to t1,q1,u1,v1 in case convection stops

  DO k=1,km
    DO i=1,im
      t1_tmp1(i,k)  = t1(i,k)
    ENDDO
  ENDDO

      do k = 1, km
        do i = 1, im
          if (k .le. kmax(i)) then
            to(i,k) = t1(i,k)
            qo(i,k) = q1(i,k)
            uo(i,k) = u1(i,k)
            vo(i,k) = v1(i,k)
            call fpvsx_ad(t1(i,k),es0,dum,dum,.false.)      ! fpvs is in pa
            es(i,k) = es0*10.0_r_kind
            es_tmp7(i,k) = es(i,k)
            qeso(i,k) = eps * es(i,k) / (pfld(i,k) + epsm1*es(i,k))
            qeso_tmp9(i,k) = qeso(i,k)
            val     =   1.e-8_r_kind
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo


!--- feedback: simply the changes from the cloud with unit mass flux
!---           multiplied by  the mass flux necessary to keep the
!---           equilibrium with the larger-scale.


      do i = 1, im
        delhbar(i) = zero
        delqbar(i) = zero
        deltbar(i) = zero
        delubar(i) = zero
        delvbar(i) = zero
        qcond(i) = zero
      enddo

      do k = 1, km
        do i = 1, im
          if (k .le. kmax(i)) then
            if(k.le.ktcon_new(i)) then
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              dellat_tmp1(i,k)  = dellat
              t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
              q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
              tem = 1./rcs(i)
              u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
              v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
              dp = 1000.0_r_kind * del(i,k)*ps(i)
              delhbar(i) = delhbar(i) + dellah(i,k)*xmb(i)*dp/g
              delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
              deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
              delubar(i) = delubar(i) + dellau(i,k)*xmb(i)*dp/g
              delvbar(i) = delvbar(i) + dellav(i,k)*xmb(i)*dp/g
            endif
          endif
        enddo
      enddo


  DO k=1,km
    DO i=1,im
      t1_tmp2(i,k)  = t1(i,k)
      q1_tmp2(i,k)  = q1(i,k)
    ENDDO
  ENDDO

      do k = 1, km
        do i = 1, im
          if ( k .le. kmax(i)) then
            if(k.le.ktcon_new(i)) then
              call fpvsx_ad(t1(i,k),es0,dum,dum,.false.)      ! fpvs is in pa
              es(i,k) = es0*10.0_r_kind
                es_tmp8(i,k) = es(i,k)
              qeso(i,k) = eps * es(i,k)/(pfld(i,k) + epsm1*es(i,k))
                qeso_tmp10(i,k) = qeso(i,k)
              val       = 1.e-8_r_kind
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo


       do i = 1, im
        do k = 1, km
            qeso_tmp11(i,k) = qeso(i,k)
        enddo
      enddo


      do k = km, 1, -1
        do i = 1, im
          if (k .le. kmax(i)) then
            if(k.lt.ktcon_new(i)) then
              aup = one
              if(k.le.kb(i))   aup = zero
              adw = one
              if(k.ge.jmin(i)) adw = zero
              rain     =  aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
              rain_tmp1(i,k) = rain
              rntot(i) = rntot(i) + rain * xmb(i) * .001_r_kind * dt2
            endif
          endif
          rntot_tmp1(i,k) = rntot(i)
        enddo
      enddo


      do k = km, 1, -1
        do i = 1, im
          if (k .le. kmax(i)) then
            deltv(i) = zero
            delq(i) = zero
            qevap(i) = zero
            dp = 1000.0_r_kind * del(i,k) *ps(i)

            if(k.lt.ktcon_new(i)) then
              aup = one 
              if(k.le.kb(i)) aup = zero
              adw = one
              if(k.ge.jmin(i)) adw = zero
              rain =  aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
                 rain_tmp2(i,k) = rain
              rn(i) = rn(i) + rain * xmb(i) * .001_r_kind * dt2
              rn_tmp1(i,k) = rn(i)
            endif
            if(k.lt.ktcon_new(i)) then
               evef = edt(i) * evfact
               if(slimsk(i) .eq. one) evef=edt(i) * evfactl
                 evef_tmp1(i,k) = evef
               qcond(i) = evef * (q1(i,k) - qeso(i,k))  &
                      / (one + el2orc * qeso(i,k) / t1(i,k)**2)
               qcond_tmp1(i,k) =qcond(i)

               if(rn(i).gt.zero .and. qcond(i).lt.zero) then
                  qevap(i) = -qcond(i) * (one-exp(-.32_r_kind*sqrt(dt2*rn(i))))
                  qevap_tmp1(i,k) = qevap(i)
                  y1  = rn(i)*1000.0_r_kind*g/dp
                  y1_tmp1(i,k) = y1
                  IF (qevap(i) .GT. y1) qevap(i)  = y1
               endif
               delq2(i) = delqev(i) + .001_r_kind * qevap(i) * dp / g
               delq2_tmp1(i,k) = delq2(i)

               if(rn(i).gt.zero .and.qcond(i).lt.zero .and. delq2(i).gt.rntot(i)) then
                  qevap(i) = 1000.0_r_kind* g * (rntot(i) - delqev(i)) / dp
               endif
               qevap_tmp2(i,k) = qevap(i)

              if(rn(i).gt.zero .and. qevap(i).gt. zero) then
                  q1(i, k)  = q1(i, k) + qevap(i)
                  t1(i, k)  = t1(i, k) - elocp*qevap(i)
                  rn(i) = rn(i) - .001_r_kind * qevap(i) * dp / g
                  deltv(i) = - elocp*qevap(i)/dt2
                  delq(i) =  + qevap(i)/dt2
                  delqev(i) = delqev(i) + .001_r_kind*dp*qevap(i)/g
               endif
            endif  !(k.lt.ktcon_new(i)) then


      END IF
    END DO
  END DO

  do i = 1, im
     rn_tmp3(i) = rn(i)
  enddo

  if (ncloud .gt. 0) then
     do k = 1, km
        do i = 1, im
          if (rn(i).gt. zero) then
            if (k.gt.kb(i) .and. k.le.ktcon_new(i)) then
               tem  = dellal(i,k) * xmb(i) * dt2
               tem_tmp11(i,k) = tem
               ql(i,k,1) = ql(i,k,1) + tem
               ql_tmp1(i,k,1)=ql(i,k,1)
               if(ql(i,k,1) .lt. 1.0E-10_r_kind) ql(i,k,1) = zero
            endif
          endif
        enddo
     enddo
  endif

      do k = 1, km
        do i = 1, im
           t_out(i,k) = t1(i,k)
           q_out(i,k) = q1(i,k)
           u_out(i,k) = u1(i,k)
           v_out(i,k) = v1(i,k)
           ql_out(i,k,1) = ql(i,k,1)
           ql_out(i,k,2) = ql(i,k,2)
        enddo
      enddo

endif !ktcno(1) .gt. 2


! Adjoint starts here

! Initialize adjoint variables 

  DO k=1,km
    DO i=1,im
      eta1d(i,k) = zero
      tod(i,k) = zero
      qod(i,k) = zero
      uod(i,k) = zero
      vod(i,k) = zero
      t1d(i,k) = zero
      q1d(i,k) = zero
      u1d(i,k) = zero
      v1d(i,k) = zero
      uo1d(i,k) = zero
      vo1d(i,k) = zero
      tvod(i,k) = zero
      qesod(i,k) = zero
      esd(i,k) = zero
      zod(i,k) = zero
      zid(i,k) = zero
      xlamued(i,k) = zero
      heod(i,k) = zero
      hesod(i,k) = zero
      hckod(i,k) = zero
      frhd(i,k) = zero
      dotd(i,k) = zero
      xlamued(i,k) = zero
      xlamue_tmp3d(i,k) = zero
      fent1d(i,k) = zero
      fent2d(i,k) = zero
      vod(i,k) = zero 
      vckod(i,k) = zero 
      uod(i,k) = zero 
      uckod(i,k) = zero 
      qckod(i,k) = zero
      dbyod(i,k) = zero
      qesod(i,k) = zero
      dellald(i,k) = zero
      pwod(i,k)= zero
      etadd(i,k) = zero
      ucdod(i,k) = zero
      vcdod(i,k) = zero
      qrcdod(i,k) = zero
      qcdod(i,k) = zero
      hcdod(i,k) = zero
      pwdod(i,k) = zero
      dellahd(i,k) = zero
      dellaqd(i,k) = zero
      dellaud(i,k) = zero
      dellavd(i,k) = zero
      eta1d(i,k) = zero
      qrcdd(i,k) = zero
      qld(i,k,1) = zero
      qld(i,k,2) = zero
    END DO
   END DO

  DO i=1,im
     pdotd(i) = zero
     cnvscl_2d(i) = zero
     cnvscl_5d(i) = zero
     cnvscl_6d(i) = zero
     cnvscl_7d(i) = zero
     cnvscl_8d(i) = zero
     pwavod(i)= zero
     pwevod(i)= zero 
     aa1d(i) = zero 
     qlko_ktcond(i) = zero
     vsheard(i) = zero
     edtd(i) = zero
     edtod(i) = zero
     edtxd(i) = zero
     xlamd0d(i) = zero
     xlamudd(i) = zero
     sumxd(i) = zero
     xpwavd(i) = zero
     xpwevd(i ) = zero
     xaa0d(i) = zero
     dtconvd(i) = zero
     acrtfctd(i) = zero
     fldd(i) =  zero
     xkd(i) = zero 
     xmbd(i) = zero
     rntotd(i) = zero
     qcondd(i) = zero
     qevapd(i) = zero
     delqevd(i) = zero
  END DO


 if(ktcon(1) .gt. 2) then

      do k = km,1,-1
        do i = 1, im
           ql_outd(i,k,2) = zero

           qld(i,k,1) = qld(i,k,1) + ql_outd(i,k,1)
           ql_outd(i,k,1) = zero

           v1d(i,k) = v1d(i,k) + v_outd(i,k)
           v_outd(i,k) = zero

           u1d(i,k) = u1d(i,k) + u_outd(i,k)
           u_outd(i,k) = zero

           q1d(i,k) = q1d(i,k) + q_outd(i,k)
           q_outd(i,k) =zero

           t1d(i,k) = t1d(i,k) + t_outd(i,k)
           t_outd(i,k) = zero
        enddo
      enddo

  IF (ncloud .GT. 0) THEN
    DO k=1,km
      DO i=1,im
         temd = zero

        IF (rn_tmp3(i) .GT. zero) THEN
          IF (k .GT. kb(i) .AND. k .LE. ktcon_new(i)) THEN

            IF(ql_tmp1(i,k,1) .lt. 1.0E-10_r_kind) then
               qld(i,k,1) = zero
            ENDIF

            temd = temd + qld(i,k,1)

            dellald(i,k) = dellald(i,k) + dt2*xmb(i)*temd
            xmbd(i) = xmbd(i) + dellal(i,k)*dt2*temd
            temd = zero

          END IF
        END IF
      END DO
    END DO
  END IF

  DO k=1,km
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        dp = 1000.0_r_kind*del(i, k)*ps(i)
        raind = zero
        evefd = zero
        y1d = zero
        tempd = zero
        dp = 1000.0_r_kind*del(i, k)*ps(i)

        IF (k .LT. ktcon_new(i)) THEN
          IF (rn_tmp1(i,k) .GT. zero .AND. qevap_tmp2(i,k) .GT. zero) THEN
            qevapd(i) = qevapd(i) + .001_r_kind*dp*delqevd(i)/g
            qevapd(i) = qevapd(i) - .001_r_kind*dp*rnd(i)/g

            qevapd(i) = qevapd(i) - elocp*t1d(i,k)
            qevapd(i) = qevapd(i) + q1d(i,k)
          END IF

          IF (rn_tmp1(i,k) .GT. zero .AND. qcond_tmp1(i,k) .LT. zero  &
          .AND. delq2_tmp1(i,k) .GT. rntot(i)) THEN
            rntotd(i) = rntotd(i) + 1000.0_r_kind*g*qevapd(i)/dp
            delqevd(i)= delqevd(i) -1000.0_r_kind*g*qevapd(i)/dp
            qevapd(i) = zero
          END IF

          IF (rn_tmp1(i,k) .GT. zero .AND. qcond_tmp1(i,k) .LT. zero) THEN
             IF (qevap_tmp1(i,k) .GT. y1_tmp1(i,k)) THEN
               y1d = y1d + qevapd(i)
               qevapd(i) = zero
             END IF
             rnd(i) = rnd(i) + 1000.0_r_kind*g*y1d/dp
             y1d = zero

             qcondd(i) = qcondd(i) - (ONE-EXP(-(.32_r_kind*SQRT(dt2*rn_tmp1(i,k)))))*qevapd(i)
             rnd(i) = rnd(i) &
             -0.32_r_kind*qcond_tmp1(i,k)*EXP(-(.32_r_kind*SQRT(dt2*rn_tmp1(i,k))))*half*dt2*qevapd(i) &
             /SQRT(dt2*rn_tmp1(i,k))
             qevapd(i) = zero
          ENDIF

          evefd = evefd + &
               qcondd(i)*(q1_tmp2(i, k)-qeso(i, k))/(one+el2orc*qeso(i, k)/t1_tmp2(i, k)**2)
          q1d(i,k) = q1d(i,k) + evef_tmp1(i,k)*qcondd(i)/(one+el2orc*qeso(i, k)/t1_tmp2(i, k)**2)
          qesod(i,k) = qesod(i,k) -  &
                evef_tmp1(i,k)*qcondd(i)/(one+el2orc*qeso(i, k)/t1_tmp2(i, k)**2)
          qesod(i,k) = qesod(i,k) -  &
                evef_tmp1(i,k)*qcondd(i)*(q1_tmp2(i, k)-qeso(i, k))*(el2orc/t1_tmp2(i,k)**2) &
                      /(one+el2orc*qeso(i, k)/t1_tmp2(i, k)**2)**2
          t1d(i,k) = t1d(i,k) -  &
                evef_tmp1(i,k)*(q1_tmp2(i, k)-qeso(i, k))*qcondd(i)*(-el2orc*qeso(i,k) &
                    *two/t1_tmp2(i,k)**3) /(one+el2orc*qeso(i, k)/t1_tmp2(i, k)**2)**2
          qcondd(i) = zero

          IF (slimsk(i) .EQ. one) THEN
             edtd(i) = edtd(i) + evfactl*evefd
             evefd = zero
          END IF
          edtd(i) = edtd(i) + evfact*evefd
          evefd=zero
        ENDIF

        IF (k .LT. ktcon_new(i)) THEN
          aup = one
          IF (k .LE. kb(i)) aup = zero
          adw =  one
          IF (k .GE. jmin(i)) adw = zero

          raind = raind + .001*dt2*xmb(i)*rnd(i)
          xmbd(i) = xmbd(i) +  .001_r_kind*dt2*rain_tmp2(i,k)*rnd(i)

          pwod(i,k)= pwod(i,k) + aup*raind
          edtod(i) = edtod(i) + adw*pwdo(i,k)*raind
          pwdod(i,k) = pwdod(i,k) + adw*edto(i)*raind
          raind = zero
        END IF

      END IF
    END DO
  END DO



  DO k=1,km
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        IF (k .LT. ktcon_new(i)) THEN
          aup = one
          IF (k .LE. kb(i)) aup = zero
          adw = one
          IF (k .GE. jmin(i)) adw = zero
             raind = zero

             raind = raind + .001_r_kind*dt2*xmb(i)*rntotd(i)
             xmbd(i) = xmbd(i) + .001_r_kind*dt2*rain_tmp1(i,k)*rntotd(i)

             pwod(i,k)=pwod(i,k) + aup*raind
             edtod(i) = edtod(i) + adw*pwdo(i,k)*raind
             pwdod(i,k) = pwdod(i,k) + adw*edto(i)*raind
             raind = zero
        END IF
      END IF
    END DO
  END DO

  DO k = km,1,-1
    DO i = 1,im
      IF (k .LE. kmax(i)) THEN
        IF (k .LE. ktcon_new(i)) THEN
          es0d = zero

          IF (qeso_tmp10(i, k) .LT. 1.e-8_r_kind) THEN
            qesod(i, k) = zero
          END IF

          esd(i,k) = esd(i,k) + eps*pfld(i,k)*qesod(i,k)/(pfld(i, k)+epsm1*es_tmp8(i, k))**2
          qesod(i,k) = zero

          es0d = es0d + 10.0_r_kind*esd(i,k)
          esd(i,k) = zero

          call fpvsx_ad(t1_tmp2(i,k), es0, t1d(i,k), es0d, .true.)      ! fpvs is in pa

        END IF
      END IF
    END DO
  END DO

  DO k=km,1,-1
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        IF (k .LE. ktcon_new(i)) THEN

          tem = one/rcs(i)
          dp = 1000.0_r_kind*del(i, k)*ps(i)
          dellatd = zero

          dellavd(i,k) = dellavd(i,k) + xmb(i)*dt2*tem*v1d(i,k)
          xmbd(i) = xmbd(i) + dellav(i,k)*dt2*tem*v1d(i,k)

          dellaud(i,k) = dellaud(i,k) + xmb(i)*dt2*tem*u1d(i,k)
          xmbd(i) = xmbd(i) + dellau(i,k)*dt2*tem*u1d(i,k)

          dellaqd(i,k) = dellaqd(i,k) + xmb(i)*dt2*q1d(i,k)
          xmbd(i) = xmbd(i) + dellaq(i,k)*dt2*q1d(i,k)

          dellatd = dellatd + xmb(i)*dt2*t1d(i,k)
          xmbd(i) = xmbd(i) + dellat_tmp1(i,k)*dt2*t1d(i,k)

          dellahd(i,k) = dellahd(i,k) + dellatd/cp
          dellaqd(i,k) = dellaqd(i,k) - hvap*dellatd/cp
          dellatd = zero

        END IF
      END IF
    END DO
  END DO

  DO k=km,1,-1
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        es0d = zero

        IF (qeso_tmp9(i, k) .LT. 1.e-8_r_kind) THEN
          qesod(i, k) = zero
        END IF

        esd(i,k) = esd(i,k)  + eps*pfld(i,k)*qesod(i,k)/(pfld(i, k)+epsm1*es_tmp7(i,k))**2
        qesod(i,k) = zero

        es0d = es0d + 10.0_r_kind*esd(i,k)
        esd(i,k)=zero

        call fpvsx_ad(t1_tmp1(i,k), es0, t1d(i,k), es0d, .true.)

        v1d(i,k) = v1d(i,k) + vod(i,k)
        vod(i,k) = zero

        u1d(i,k) = u1d(i,k) + uod(i,k)
        uod(i,k) = zero

        q1d(i,k) = q1d(i,k) + qod(i,k)
        qod(i,k) = zero

        t1d(i,k) = t1d(i,k) + tod(i,k)
        tod(i,k) = zero
      END IF
    END DO
  END DO


 do i = 1, im
    xmbtmp1d = zero

    xmbtmp1d = xmbtmp1d + xmbd(i)
    xmbd(i) = zero

    xmbd(i) = xmbd(i) + cnvscl_8(i)*xmbtmp1d
    cnvscl_8d(i) = cnvscl_8d(i) + xmb_tmp2(i)*xmbtmp1d
    xmbtmp1d  =  zero


    IF (xmb_tmp1(i) .GT. xmbmax(i)) THEN
      xmbd(i) = zero
    END IF


    fldd(i) = fldd(i) - xmbd(i)/xk(i)
    xkd(i) = xkd(i) + fld(i)*xmbd(i)/xk(i)**2
    xmbd(i) = zero

    cnvscl_7d(i) = cnvscl_7d(i) + cnvscl_8_tmp1(i)*cnvscl_8d(i)
    cnvscl_8d(i) = cnvscl_8d(i) + cnvscl_7(i)*cnvscl_8d(i)

    if(xk(i) .ge. zero) then
        cnvscl_8d(i) = zero
    endif

    xaa0d(i) = xaa0d(i) + xkd(i)/mbdt
    aa1d(i) = aa1d(i) - xkd(i)/mbdt
    xkd(i) = zero

    cnvscl_6d(i) = cnvscl_6d(i) + cnvscl_7_tmp1(i)*cnvscl_7d(i)
    cnvscl_7d(i) = cnvscl_7d(i) + cnvscl_6(i)*cnvscl_7d(i)

    if(fld(i) .le. zero) then
      cnvscl_7d(i) = zero
    endif

    aa1d(i) = aa1d(i) + fldd(i)/dtconv(i)
    acrtfctd(i) = acrtfctd(i) - acrt(i)*fldd(i)/dtconv(i)
    dtconvd(i) = dtconvd(i) - (aa1(i)-acrt(i)*acrtfct(i))*fldd(i)/dtconv(i)**2
    fldd(i) = zero
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

    IF (1800.0_r_kind - dt2 .LT. zero) THEN
      max5 = zero
    ELSE
      max5 = 1800.0_r_kind - dt2
    END IF


    IF (dtconv_tmp2(i) .GT. dtmax) THEN
      dtconvd(i) = zero
    END IF

    IF (dtconv_tmp1(i) .LT. dtmin) THEN
      dtconvd(i) = zero
    END IF

    pdotd(i) = pdotd(i) + max5*dtconvd(i)/(w1-w2)
    dtconvd(i) = zero


    acrtfctd(i) = -acrtfctd(i)

    IF (acrtfct_tmp2(i) .GT. one) THEN
      acrtfctd(i) = zero
    END IF

    IF (acrtfct_tmp1(i) .LT. -one) THEN
      acrtfctd(i) = zero
    END IF


    IF (pdot(i) .LE. w4) THEN
        pdotd(i) = pdotd(i) + acrtfctd(i)/(w3-w4)
        acrtfctd(i) = zero
    ELSE IF (pdot(i) .GE. -w4) THEN
        pdotd(i) = pdotd(i) - acrtfctd(i)/(w4-w3)
        acrtfctd(i) = zero
    ELSE
      acrtfctd(i) = zero
    END IF

 ENDDO



  DO k=1,km1
    DO i=1,im
      IF (k .LT. jmin(i)) THEN

         dzd = zero
         dhd = zero
         dhhd = zero
         dgd = zero
         dtd = zero
         gammad = zero
         max4d = zero

         edtxd(i) = edtxd(i) + g*delta*dz_tmp17(i,k)*max4_tmp(i,k)*xaa0d(i)
         dzd = dzd + edtx_tmp3(i)*g*delta*max4_tmp(i,k)*xaa0d(i)
         max4d = max4d + edtx_tmp3(i)*dz_tmp17(i,k)*g*delta*xaa0d(i)

        IF (zero .LT. qeso_tmp8(i, k) - qo_tmp7(i, k)) THEN
          qesod(i,k) = qesod(i,k) + max4d
          qod(i,k) = qod(i,k) - max4d
          max4d = zero
        ELSE
          max4d = zero
        END IF

        edtxd(i) = edtxd(i) + dz_tmp17(i,k)*(g/(cp*dt_tmp4(i,k)))*((dhh_tmp2(i,k)-dh_tmp3(i,k))/ &
           (one+dg_tmp2(i,k)))*(one+delta*cp*dg_tmp2(i,k)*dt_tmp4(i,k)/hvap)*xaa0d(i)
        dzd = dzd + edtx_tmp3(i)*(g/(cp*dt_tmp4(i,k)))*((dhh_tmp2(i,k)-dh_tmp3(i,k))/(one+dg_tmp2(i,k))) &
             * (one+delta*cp*dg_tmp2(i,k)*dt_tmp4(i,k)/hvap)*xaa0d(i)
        dtd = dtd  - edtx_tmp3(i)*dz_tmp17(i,k)*(g/(cp*dt_tmp4(i,k)*dt_tmp4(i,k)))* &
          ((dhh_tmp2(i,k)-dh_tmp3(i,k))/(one+dg_tmp2(i,k)))*(one+delta*cp*dg_tmp2(i,k) &
                *dt_tmp4(i,k)/hvap)*xaa0d(i)
        dhhd=  dhhd + edtx_tmp3(i)*dz_tmp17(i,k)*(g/(cp*dt_tmp4(i,k)))*(one/(one+dg_tmp2(i,k))) &
           *(one+delta*cp*dg_tmp2(i,k)*dt_tmp4(i,k)/hvap)*xaa0d(i)
        dhd=  dhd - edtx_tmp3(i)*dz_tmp17(i,k)*(g/(cp*dt_tmp4(i,k)))*(one/(one+dg_tmp2(i,k)))* &
            (one+delta*cp*dg_tmp2(i,k)*dt_tmp4(i,k)/hvap)*xaa0d(i)
        dgd = dgd - edtx_tmp3(i)*dz_tmp17(i,k)*(g/(cp*dt_tmp4(i,k)))*((dhh_tmp2(i,k)-dh_tmp3(i,k)) &
           /(one+dg_tmp2(i,k))**2)*(one+delta*cp*dg_tmp2(i,k)*dt_tmp4(i,k)/hvap)*xaa0d(i)
        dgd = dgd + edtx_tmp3(i)*dz_tmp17(i,k)*(g/(cp*dt_tmp4(i,k)))*((dhh_tmp2(i,k)-dh_tmp3(i,k)) &
            /(one+dg_tmp2(i,k)))*delta*cp*dt_tmp4(i,k)*xaa0d(i)/hvap
        dtd = dtd + edtx_tmp3(i)*dz_tmp17(i,k)*(g/(cp*dt_tmp4(i,k)))*((dhh_tmp2(i,k)-dh_tmp3(i,k)) &
            /(one+dg_tmp2(i,k)))*delta*cp*dg_tmp2(i,k)*xaa0d(i)/hvap

        zod(i,k+1) = zod(i,k+1) - dzd
        zod(i,k) = zod(i,k) + dzd
        dzd = zero

        hesod(i,k) = hesod(i,k) + dhd
        dhd = zero

        gammad = gammad + dgd
        dgd = zero

        tod(i,k) = tod(i,k) + dtd
        dtd = zero

        hcdod(i,k) = hcdod(i,k) + dhhd
        dhhd = zero

        qesod(i,k) = qesod(i,k) + el2orc*gammad/to_tmp4(i,k)**2
        tod(i,k) = tod(i,k) - el2orc*qeso_tmp8(i,k)*two*gammad/to_tmp4(i,k)**3
        gammad = zero

      END IF
    END DO
  END DO




  DO i = 1, im
    edtmax = edtmaxl
    IF (slimsk(i) .EQ. zero) edtmax = edtmaxs
       edtx_tmpd = zero

       IF (xpwev_tmp1(i) .GE. zero) THEN
          edtxd(i) = zero
       ELSE
          IF (edtx_tmp2(i) .GT. edtmax) THEN
           edtxd(i)= zero
          END IF

          edtx_tmpd = edtx_tmpd + edtxd(i)
          edtxd(i) = zero

          edtxd(i) = edtxd(i) - edtx_tmpd*xpwav_tmp1(i)/xpwev_tmp1(i)
          xpwavd(i) = xpwavd(i) - edtx_tmp1(i)*edtx_tmpd/xpwev_tmp1(i)
          xpwevd(i) = xpwevd(i) + edtx_tmp1(i)*xpwav_tmp1(i)*edtx_tmpd/xpwev_tmp1(i)**2
          edtx_tmpd = zero
     END IF
  ENDDO


  DO k=1,km1
    DO i=1,im
      IF (k .LT. jmin(i)) THEN
         factord = zero
         xpwdd = zero
         temd = zero
         tem1d =zero
         dzd = zero
         dqd = zero
         dtd = zero
         gammad = zero

         xpwdd = xpwdd + xpwevd(i)

         qrcdd(i,k) = qrcdd(i,k) + qcdod(i,k)
         qcdod(i,k) = zero

         etadd(i,k+1) = etadd(i,k+1) + (qcdo_tmp3(i, k)-qrcd_tmp1(i, k))*xpwdd
         qcdod(i,k) = qcdod(i,k) + etad_tmp1(i,k+1)*xpwdd
         qrcdd(i,k) = qrcdd(i,k) - etad_tmp1(i,k+1)*xpwdd
         xpwdd = zero

         qcdod(i,k+1) = qcdod(i,k+1) + (one-tem1_tmp10(i,k))*qcdod(i,k) / factor_tmp9(i,k)
         tem1d = tem1d - qcdo_tmp4(i,k+1)*qcdod(i,k)/factor_tmp9(i,k)
         temd = temd + half*(qo_tmp7(i,k) + qo_tmp7(i,k+1))*qcdod(i,k)/factor_tmp9(i,k)
         qod(i,k) = qod(i,k) + half*tem_tmp10(i,k)*qcdod(i,k)/factor_tmp9(i,k)
         qod(i,k+1) = qod(i,k+1) + half*tem_tmp10(i,k)*qcdod(i,k)/factor_tmp9(i,k)
         factord = factord -  ((one-tem1_tmp10(i,k))*qcdo_tmp4(i, k+1)+ &
           tem_tmp10(i,k)*half*(qo_tmp7(i, k)+qo_tmp7(i, k+1)))*qcdod(i,k)/factor_tmp9(i,k)**2
         qcdod(i,k) = zero

         temd = temd + factord
         tem1d = tem1d - factord
         factord = zero

        IF (k .GE. kbcon(i)) THEN
          dzd = dzd  + half*xlamdd*tem1d
          tem1d = zero
          dzd = dzd + xlamde*temd
          temd = zero
        ELSE
          xlamd0d(i) = xlamd0d(i)+half*dz_tmp16(i,k)*tem1d
          dzd = dzd + half*(xlamd(i)+xlamdd)*tem1d
          tem1d = zero
          dzd = dzd + xlamde*temd
          temd = zero
        END IF

        zid(i,k+1) = zid(i,k+1) + dzd
        zid(i,k) = zid(i,k) - dzd
        dzd = zero

        dqd = dqd + qrcdd(i,k)
        dhd = dhd + one/hvap*(gamma_tmp10(i,k)/(one+gamma_tmp10(i,k)))*qrcdd(i,k)
        gammad = gammad +  one/hvap*(qrcdd(i,k)/(one+gamma_tmp10(i,k)))*dh_tmp2(i,k)
        gammad = gammad - one/hvap*(gamma_tmp10(i,k)/(one+gamma_tmp10(i,k))**2)*dh_tmp2(i,k)*qrcdd(i,k)
        qrcdd(i,k) = zero

        hcdod(i,k) = hcdod(i,k) + dhd
        hesod(i,k) = hesod(i,k) - dhd
        dhd = zero

        dqd = dqd + el2orc*gammad/dt_tmp3(i,k)**2
        dtd = dtd - el2orc*two*dq_tmp5(i,k)*gammad/dt_tmp3(i,k)**3
        gammad = zero

        tod(i,k) = tod(i,k) + dtd
        dtd = zero

        qesod(i,k) = qesod(i,k) + dqd
        dqd = zero
      END IF
    END DO
  END DO

  DO k=1,km1
    DO i=1,im
      tem1d = zero
      temd = zero
      dzd = zero
      factord = zero


      IF (k .LT. jmin(i)) THEN
          hcdod(i,k+1) = hcdod(i,k+1) + (one-tem1_tmp9(i,k))*hcdod(i,k)/factor_tmp8(i,k)
          tem1d = tem1d - hcdo_tmp2(i,k+1)*hcdod(i,k)/factor_tmp8(i,k)
          temd = temd + half*(heo_tmp3(i,k) + heo_tmp3(i,k+1))*hcdod(i,k)/factor_tmp8(i,k)
          heod(i,k) = heod(i,k) + tem_tmp9(i,k)*half*hcdod(i,k)/factor_tmp8(i,k)
          heod(i,k+1) = heod(i,k+1) + tem_tmp9(i,k)*half*hcdod(i,k)/factor_tmp8(i,k)
          factord = factord - ((one-tem1_tmp9(i,k))*hcdo_tmp2(i, k+1)+tem_tmp9(i,k) &
             *half*(heo_tmp3(i, k)+heo_tmp3(i, k+1)))*hcdod(i,k)/factor_tmp8(i,k)**2
          hcdod(i,k) = zero

          temd = temd + factord
          tem1d = tem1d - factord
          factord = zero

          IF (k .GE. kbcon(i)) THEN
            dzd = dzd + tem1d*half*xlamdd
            tem1d = zero

            dzd = dzd + xlamde*temd
            temd = zero
          ELSE
            xlamd0d(i) = xlamd0d(i) + half*dz_tmp15(i,k)*tem1d
            dzd = dzd + half*(xlamd(i)+xlamdd)*tem1d
            tem1d = zero

            dzd = dzd + xlamde*temd
            temd = zero
          END IF

          zid(i,k+1) = zid(i,k+1) + dzd
          zid(i,k) = zid(i,k) - dzd
          dzd = zero
      END IF
    END DO
  END DO

  DO i=1,im
    jmn = jmin(i)
    heod(i,jmn) = heod(i,jmn) + hcdod(i,jmn)
    hcdod(i,jmn) = zero

    qod(i,jmn) = qod(i,jmn) + qcdod(i,jmn)
    qcdod(i,jmn) = zero

    qesod(i,jmn) = qesod(i,jmn) + qrcdd(i,jmn)
    qrcdd(i,jmn) = zero
  END DO

  do k = km1,2,-1
    do i = 1, im
       dzd = zero
       dz1d = zero
       max3d = zero
       xdbyd= zero
       gammad = zero
       xpwd0 = zero
       qlkd = zero
       xpwd0 =zero
       xqrchd = zero
       dqd = zero
       etahd = zero
       temd = zero
       tem1d= zero
       factord = zero
       rfactd= zero

       if (k .GE. kbcon(i) .AND. k .LT. ktcon1_new(i)) THEN
            dz1d = dz1d + g*delta*max3_tmp(i,k)*xaa0d(i)
            max3d = max3d + g*delta*dz1_tmp2(i,k)*xaa0d(i)

            IF (zero .LT. qeso_tmp8(i, k) - qo_tmp7(i, k)) THEN
               qesod(i,k) = qesod(i,k) + max3d
               qod(i,k) = qod(i,k) - max3d
               max3d = zero
            ELSE
               max3d = zero
            END IF

            dz1d = dz1d +  xaa0d(i)*(g/(cp*to_tmp4(i, k)))*xdby_tmp1(i,k)/(one+gamma_tmp9(i,k))*rfact_tmp2(i,k)
            tod(i,k) = tod(i,k) - xaa0d(i)*dz1_tmp2(i,k)*(g/(cp*to_tmp4(i,k)*to_tmp4(i,k))) &
                       *xdby_tmp1(i,k)/(one+gamma_tmp9(i,k))*rfact_tmp2(i,k)
            xdbyd = xdbyd +  dz1_tmp2(i,k)*(g/(cp*to_tmp4(i, k)))*xaa0d(i)/(one+gamma_tmp9(i,k))*rfact_tmp2(i,k)
            gammad = gammad - xaa0d(i)*dz1_tmp2(i,k)*(g/(cp*to_tmp4(i, k))) &
               *xdby_tmp1(i,k)/(one+gamma_tmp9(i,k))**2*rfact_tmp2(i,k)
            rfactd = rfactd + dz1_tmp2(i,k)*(g/(cp*to_tmp4(i, k)))*xdby_tmp1(i,k)/(one+gamma_tmp9(i,k))*xaa0d(i)

            gammad = gammad + delta*cp*to_tmp4(i,k)*rfactd/hvap
            tod(i,k) = tod(i,k) + delta*cp*gamma_tmp9(i,k)*rfactd/hvap
            rfactd = zero

            qesod(i,k) = qesod(i,k) + el2orc*gammad/to_tmp4(i,k)**2
            tod(i,k) = tod(i,k) - el2orc*qeso_tmp8(i,k)*two*gammad/to_tmp4(i,k)**3
            gammad = zero

            zod(i,k+1) = zod(i,k+1) + dz1d
            zod(i,k) = zod(i,k) - dz1d
            dz1d = zero
       endif

       if (k .GT. kb(i) .AND. k .LT. ktcon_new(i)) THEN
          if (k .GE. kbcon(i) .AND. dq_tmp4(i,k) .GT. zero) THEN
             xpwd0 = xpwd0 + xpwavd(i)

             etahd = etahd + c0*dz_tmp14(i,k)*qlk_tmp3(i,k)*xpwd0
             dzd = dzd + etah_tmp3(i,k)*c0*qlk_tmp3(i,k)*xpwd0
             qlkd = qlkd + etah_tmp3(i,k)*c0*dz_tmp14(i,k)*xpwd0
             xpwd0 = zero

             qlkd = qlkd + qckod(i,k)
             xqrchd = xqrchd + qckod(i,k)
             qckod(i,k) = zero

             IF (k .lt. ktcon1_new(i)) then
               dzd = dzd - g*qlk_tmp3(i,k)*xaa0d(i)
               qlkd = qlkd - dz_tmp14(i,k)*g*xaa0d(i)
             ENDIF

             IF (ncloud .GT. zero .AND. k .GT. jmin(i)) THEN
                dqd = dqd + qlkd/(eta_tmp3(i, k)+etah_tmp3(i,k)*(c0+c1)*dz_tmp14(i,k))
                eta1d(i,k) = eta1d(i,k) - dq_tmp4(i,k)*qlkd/(eta_tmp3(i, k)+etah_tmp3(i,k)*(c0+c1)*dz_tmp14(i,k))**2
                etahd = etahd - dq_tmp4(i,k)*(c0+c1)*dz_tmp14(i,k)*qlkd &
                      /(eta_tmp3(i, k)+etah_tmp3(i,k)*(c0+c1)*dz_tmp14(i,k))**2
                dzd = dzd - dq_tmp4(i,k)*(c0+c1)*etah_tmp3(i,k)*qlkd &
                      /(eta_tmp3(i, k)+etah_tmp3(i,k)*(c0+c1)*dz_tmp14(i,k))**2
                qlkd = zero
             ELSE
                dqd = dqd + qlkd/(eta_tmp3(i, k)+etah_tmp3(i,k)*c0*dz_tmp14(i,k))
                eta1d(i,k) = eta1d(i,k) - dq_tmp4(i,k)*qlkd/(eta_tmp3(i, k)+etah_tmp3(i,k)*c0*dz_tmp14(i,k))**2
                etahd = etahd - dq_tmp4(i,k)*c0*dz_tmp14(i,k)*qlkd/(eta_tmp3(i, k) &
                       +etah_tmp3(i,k)*c0*dz_tmp14(i,k))**2
                dzd = dzd - dq_tmp4(i,k)*c0*etah_tmp3(i,k)*qlkd/(eta_tmp3(i, k)+etah_tmp3(i,k)*c0*dz_tmp14(i,k))**2
                qlkd = zero
             END IF

             eta1d(i,k) =eta1d(i,k) + half*etahd
             eta1d(i,k-1) =eta1d(i,k-1) + half*etahd
             etahd = zero
           END IF


        eta1d(i,k) = eta1d(i,k) + (qcko_tmp6(i, k)-xqrch_tmp1(i,k))*dqd
        qckod(i,k) = qckod(i,k) + eta_tmp3(i,k)*dqd
        xqrchd = xqrchd - eta_tmp3(i,k)*dqd
        dqd = zero

        qckod(i,k-1) = qckod(i,k-1) + (one-tem1_tmp8(i,k))*qckod(i,k)/factor_tmp7(i,k)
        tem1d = tem1d - qcko_tmp7(i,k-1)*qckod(i,k)/factor_tmp7(i,k)
        temd = temd + half*(qo_tmp7(i, k)+qo_tmp7(i, k-1))*qckod(i,k)/factor_tmp7(i,k)
        qod(i,k) = qod(i,k) + tem_tmp8(i,k)*half*qckod(i,k)/factor_tmp7(i,k)
        qod(i,k-1) = qod(i,k-1) + tem_tmp8(i,k)*half*qckod(i,k)/factor_tmp7(i,k)
        factord = factord - ((one-tem1_tmp8(i,k))*qcko_tmp7(i, k-1)+ &
           tem_tmp8(i,k)*half*(qo_tmp7(i, k)+qo_tmp7(i, k-1)))*qckod(i,k)/factor_tmp7(i,k)**2
        qckod(i,k) = zero

        temd = temd + factord
        tem1d = tem1d - factord
        factord = zero

        xlamudd(i) = xlamudd(i) + half*tem1d*dz_tmp14(i,k)
        dzd = dzd + half*tem1d*xlamud(i)
        tem1d = zero

        xlamued(i,k) = xlamued(i,k)+ half*dz_tmp14(i,k)*temd
        xlamued(i,k-1) = xlamued(i,k-1)+ half*dz_tmp14(i,k)*temd
        dzd = dzd + half*(xlamue_tmp4(i, k)+xlamue_tmp4(i, k-1))*temd
        temd = zero

        qesod(i,k) = qesod(i,k) + xqrchd
        gammad = gammad + xdby_tmp1(i,k)*xqrchd/(hvap*(one+gamma_tmp8(i,k)))
        xdbyd = xdbyd + gamma_tmp8(i,k)*xqrchd/(hvap*(one+gamma_tmp8(i,k)))
        gammad = gammad - gamma_tmp8(i,k)*xdby_tmp1(i,k)*hvap*xqrchd/(hvap*(one+gamma_tmp8(i,k)))**2
        xqrchd = zero

        hckod(i,k) = hckod(i,k) + xdbyd
        hesod(i,k) = hesod(i,k) - xdbyd
        xdbyd = zero

        qesod(i,k) = qesod(i,k) + el2orc*gammad/to_tmp4(i,k)**2
        tod(i,k) = tod(i,k) - el2orc*two*qeso_tmp8(i,k)*gammad/to_tmp4(i,k)**3
        gammad = zero

        zid(i,k) = zid(i,k) + dzd
        zid(i,k-1) = zid(i,k-1) - dzd
        dzd = zero
       endif
    enddo
  enddo



  DO k=km1,2,-1
    DO i=1,im
      IF (k .GT. kb(i) .AND. k .LE. ktcon_new(i)) THEN
          temd = zero
          tem1d = zero
          factord = zero
          dzd = zero

          hckod(i,k-1) = hckod(i,k-1) + (one-tem1_tmp7(i,k))*hckod(i,k)/factor_tmp6(i,k)
          tem1d = tem1d  - hcko_tmp3(i,k-1)*hckod(i,k)/factor_tmp6(i,k)
          temd= temd + half*(heo_tmp3(i, k)+heo_tmp3(i, k-1))*hckod(i,k)/factor_tmp6(i,k)
          heod(i,k) = heod(i,k) + tem_tmp7(i,k)*half*hckod(i,k)/factor_tmp6(i,k)
          heod(i,k-1) = heod(i,k-1) + tem_tmp7(i,k)*half*hckod(i,k)/factor_tmp6(i,k)
          factord =  factord - hckod(i,k)*((one-tem1_tmp7(i,k))*hcko_tmp3(i, k-1)+ &
                tem_tmp7(i,k)*half*(heo_tmp3(i, k)+heo_tmp3(i, k-1)))/factor_tmp6(i,k)**2
          hckod(i,k) = zero

          temd = temd + factord
          tem1d = tem1d - factord
          factord = zero

          xlamudd(i) = xlamudd(i) + half*dz_tmp13(i,k)*tem1d
          dzd = dzd + half*xlamud(i)*tem1d
          tem1d = zero

          xlamued(i,k) = xlamued(i,k) + half*dz_tmp13(i,k)*temd
          xlamued(i,k-1) = xlamued(i,k-1) + half*dz_tmp13(i,k)*temd
          dzd = dzd +  half*(xlamue_tmp4(i, k)+xlamue_tmp4(i, k-1))*temd
          temd = zero

          zid(i,k) = zid(i,k) + dzd
          zid(i,k-1) = zid(i,k-1) - dzd
          dzd = zero
      END IF
    END DO
  END DO

  DO i=1,im
    indx = kb(i)
    heod(i,indx) = heod(i,indx) + hckod(i,indx)
    hckod(i,indx) = zero
    qod(i,indx) = qod(i,indx) + qckod(i,indx)
    qckod(i,indx) = zero
  END DO


  DO i=1,im
    k = kmax(i)
    zod(i,k) = zod(i,k) + g*hesod(i,k)
    qesod(i,k) = qesod(i,k) + hvap*hesod(i,k)
    tod(i,k) = tod(i,k) + cp*hesod(i,k)
    hesod(i,k) = zero

    zod(i,k) = zod(i,k) + g*heod(i,k)
    qod(i,k) = qod(i,k) + hvap*heod(i,k)
    tod(i,k) = tod(i,k) + cp*heod(i,k)
    heod(i,k) = zero
  END DO


  DO k=km1,1,-1
    DO i=1,im
      IF (k .LE. kmax(i) - 1) THEN
         es0d = zero
         zod(i,k) = zod(i,k) + half*g*hesod(i,k)
         zod(i,k+1) = zod(i,k+1) + half*g*hesod(i,k)
         tod(i,k) = tod(i,k) + cp*hesod(i,k)
         qesod(i,k) = qesod(i,k) + hvap*hesod(i,k)
         hesod(i,k) = zero

         zod(i,k) = zod(i,k) + half*g*heod(i,k)
         zod(i,k+1) = zod(i,k+1) + half*g*heod(i,k)
         tod(i,k) = tod(i,k) + cp*heod(i,k)
         qod(i,k) = qod(i,k) + hvap*heod(i,k)
         heod(i,k) = zero

        IF (qo_tmp6(i, k) .LT. 1.e-10_r_kind) THEN
          qod(i, k) = zero
        ENDIF
        IF (qeso_tmp7(i, k) .LT. 1.e-8_r_kind) THEN
          qesod(i, k) = zero
        ENDIF

        esd(i,k) = esd(i,k) + eps*po(i,k)*qesod(i,k)/(po(i, k)+epsm1*es_tmp6(i, k))**2
        qesod(i,k) = zero

        es0d = es0d + 10.0_r_kind*esd(i,k)
        esd(i,k) = zero

        call fpvsx_ad(to_tmp4(i,k),es0, tod(i,k),es0d,.true.)

      END IF
    END DO
  END DO

  DO k=km1,1,-1
    DO i=1,im
      IF (k .LE. kmax(i) - 1) THEN
        dp  = half*(pfld(i, k+1)-pfld(i, k))
        dtd = zero
        dqd = zero
        dzd = zero
        qsd = zero
        es0d = zero
        dqsdpd = zero
        dqsdtd = zero
        gammad = zero
        pprimed = zero
        desdtd = zero

        qod(i,k+1)  =  qod(i,k+1) + qod(i,k)
        dqd = dqd + qod(i,k)
        qod(i,k) = zero

        tod(i,k+1) = tod(i,k+1) + tod(i,k)
        dtd = dtd + tod(i,k)
        tod(i,k) = zero

        dqsdtd  = dqsdtd + dt_tmp2(i,k)*dqd
        dtd = dtd + dqsdt_tmp1(i,k)*dqd
        dqsdpd = dqsdpd + dp*dqd
        dqd = zero

        dzd = dzd + g*dtd/(cp*(one+gamma_tmp7(i,k)))
        dqsdpd = dqsdpd + hvap*dp*dtd/(cp*(one+gamma_tmp7(i,k)))
        gammad = gammad - (g*dz_tmp12(i,k)+hvap*dqsdp_tmp1(i,k)*dp)*cp*dtd/(cp*(one+gamma_tmp7(i,k)))**2
        dtd = zero

        qesod(i,k+1) = qesod(i,k+1) + el2orc*gammad/to_tmp3(i, k+1)**2
        tod(i,k+1) = tod(i,k+1) - el2orc*qeso_tmp6(i,k+1)*two*gammad/(to_tmp3(i, k+1)**3)
        gammad = zero

        qsd = qsd + pfld(i,k+1)*desdt_tmp1(i,k)*dqsdtd/(es_tmp5(i,k)*pprime_tmp1(i,k))
        desdtd = desdtd + pfld(i, k+1)*qs_tmp2(i,k)*dqsdtd/(es_tmp5(i,k)*pprime_tmp1(i,k))
        esd(i,k) = esd(i,k) - qs_tmp2(i,k)*pfld(i, k+1)*desdt_tmp1(i,k)*dqsdtd/ &
            (es_tmp5(i,k)*es_tmp5(i,k)*pprime_tmp1(i,k))
        pprimed = pprimed - qs_tmp2(i,k)*pfld(i, k+1)*desdt_tmp1(i,k)*dqsdtd/ &
            (es_tmp5(i,k)*pprime_tmp1(i,k)*pprime_tmp1(i,k))
        dqsdtd = zero

        esd(i,k) = esd(i,k)  + (fact1/to_tmp3(i, k+1)+fact2/to_tmp3(i, k+1)**2)*desdtd
        tod(i,k+1) = tod(i,k+1) - es_tmp5(i,k)*fact1*desdtd/to_tmp3(i,k+1)**2
        tod(i,k+1) = tod(i,k+1) - es_tmp5(i,k)*fact2*two*desdtd/to_tmp3(i,k+1)**3
        desdtd = zero

        qsd = qsd - dqsdpd/pprime_tmp1(i,k)
        pprimed = pprimed + qs_tmp2(i,k)*dqsdpd/pprime_tmp1(i,k)**2
        dqsdpd = zero

        esd(i,k) = esd(i,k) + eps*qsd/pprime_tmp1(i,k)
        pprimed = pprimed - eps*es_tmp5(i,k)*qsd/pprime_tmp1(i,k)**2
        qsd = zero

        esd(i,k) = esd(i,k) + epsm1*pprimed
        pprimed =  zero

        es0d = es0d + 10.0_r_kind*esd(i,k)
        esd(i,k) = zero

        call fpvsx_ad(to_tmp3(i,k+1),es0,tod(i,k+1),es0d,.true.)

        zod(i,k+1) = zod(i,k+1) + half*dzd
        zod(i,k) = zod(i,k) - half*dzd
        dzd = zero
      END IF
    END DO
  END DO



  DO k=km,1,-1
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        es0d = zero
        IF (qeso_tmp5(i, k) .LT. 1.e-8_r_kind) THEN
          qesod(i, k) = zero
        END IF

        esd(i,k) = esd(i,k) + eps*pfld(i,k)*qesod(i,k)/(pfld(i, k)+epsm1*es_tmp4(i, k))**2
        qesod(i,k) = zero

        es0d = es0d + 10.0_r_kind*esd(i,k)
        esd(i,k) = zero

        call fpvsx_ad(to_tmp3(i,k),es0,tod(i,k),es0d,.true.)

      END IF
    END DO
  END DO

  DO k=km,1,-1
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        dellatd = zero

        IF (k .LE. ktcon_new(i)) THEN
          IF (qo_tmp5(i, k) .LT. 1.0e-10) THEN
            qod(i, k) = zero
          END IF
          dellatd = dellatd + mbdt*tod(i,k)
          t1d(i,k) = t1d(i,k) + tod(i,k)
          tod(i,k) = zero

          dellahd(i,k) = dellahd(i,k) + dellatd/cp
          dellaqd(i,k) = dellaqd(i,k) - dellatd*hvap/cp
          dellatd = zero

          dellaqd(i,k) = dellaqd(i,k) + mbdt*qod(i,k)
          q1d(i,k) = q1d(i,k) + qod(i,k)
          qod(i,k) = zero
        END IF

        IF (k .GT. ktcon_new(i)) THEN
           t1d(i,k) = t1d(i,k) + tod(i,k)
           tod(i,k) = zero
           q1d(i,k) = q1d(i,k) + qod(i,k)
           qod(i,k) = zero
        END IF
      END IF
    END DO
  END DO

  DO i=1,im
    indx = ktcon_new(i)
    dp   = 1000.0_r_kind*del(i, indx)*ps(i)
     dv1hd = zero
     dv1qd = zero
     dv1ud = zero
     dv1vd = zero

     eta1d(i,indx-1) = eta1d(i,indx-1) + qlko_ktcon(i)*dellald(i,indx)*g/dp
     qlko_ktcond(i) = qlko_ktcond(i) + eta_tmp3(i,indx-1)*dellald(i,indx)*g/dp
     dellald(i,indx) = zero

     eta1d(i,indx-1) = eta1d(i,indx-1) + (vcko_tmp2(i,indx-1)-dv1v_tmp2)*dellavd(i,indx)*g/dp
     vckod(i,indx-1) = vckod(i,indx-1) + eta_tmp3(i,indx-1)*dellavd(i,indx)*g/dp
     dv1vd = dv1vd - eta_tmp3(i,indx-1)*dellavd(i,indx)*g/dp
     dellavd(i,indx) =  zero

     vod(i,indx-1) = vod(i,indx-1) + dv1vd
     dv1vd = zero

     eta1d(i,indx-1) = eta1d(i,indx-1) + (ucko_tmp2(i,indx-1)-dv1u_tmp2)*dellaud(i,indx)*g/dp
     uckod(i,indx-1) = uckod(i,indx-1) + eta_tmp3(i,indx-1)*dellaud(i,indx)*g/dp
     dv1ud = dv1ud - eta_tmp3(i,indx-1)*dellaud(i,indx)*g/dp
     dellaud(i,indx) = zero

     uod(i,indx-1) = uod(i,indx-1) + dv1ud
     dv1ud = zero

     eta1d(i,indx-1) = eta1d(i,indx-1) + (qcko_tmp5(i,indx-1)-dv1q_tmp2)*dellaqd(i,indx)*g/dp
     qckod(i,indx-1) = qckod(i,indx-1) + eta_tmp3(i,indx-1)*dellaqd(i,indx)*g/dp
     dv1qd = dv1qd - eta_tmp3(i,indx-1)*dellaqd(i,indx)*g/dp
     dellaqd(i,indx) = zero

     qod(i,indx-1) = qod(i,indx-1) + dv1qd
     dv1qd  = zero

     eta1d(i,indx-1) = eta1d(i,indx-1) + (hcko_tmp2(i, indx-1)-dv1h_tmp2)*dellahd(i,indx)*g/dp
     hckod(i,indx-1) = hckod(i,indx-1) + eta_tmp3(i,indx-1)*dellahd(i,indx)*g/dp
     dv1hd = dv1hd - eta_tmp3(i,indx-1)*dellahd(i,indx)*g/dp
     dellahd(i,indx) = zero

     heod(i,indx-1) = heod(i,indx-1) + dv1hd
     dv1hd = zero
  END DO

  DO k=km1,2,-1
    DO i=1,im
      IF (k .LT. ktcon_new(i)) THEN
        aup = one
        IF (k .LE. kb(i))   aup = zero
        adw = one
        IF (k .GT. jmin(i)) adw = zero
        dp    = 1000.0_r_kind*del(i, k)*ps(i)

        dzd = zero
        tem1d = zero
        temd = zero
        ptem1d=zero

        dv1hd = zero
        dv2hd = zero
        dv3hd = zero
        dv1qd = zero
        dv2qd = zero
        dv3qd = zero
        dv1ud = zero
        dv2ud = zero
        dv3ud = zero
        dv1vd = zero
        dv2vd = zero
        dv3vd = zero

        eta1d(i,k-1) = eta1d(i,k-1) - pgcon*aup*(dv1v_tmp1(i,k)-dv3v_tmp1(i,k))*dellavd(i,k)*g/dp
        edtod(i) = edtod(i) + pgcon*adw*etad_tmp1(i,k)*(dv1v_tmp1(i,k)-dv3v_tmp1(i,k))*dellavd(i,k)*g/dp
        etadd(i,k) = etadd(i,k) + pgcon*adw*edto_tmp3(i)*(dv1v_tmp1(i,k)-dv3v_tmp1(i,k))*dellavd(i,k)*g/dp
        dv1vd = dv1vd - pgcon*(aup*eta_tmp3(i, k-1)-adw*edto_tmp3(i)*etad_tmp1(i, k))*dellavd(i,k)*g/dp
        dv3vd = dv3vd + pgcon*(aup*eta_tmp3(i, k-1)-adw*edto_tmp3(i)*etad_tmp1(i, k))*dellavd(i,k)*g/dp

        eta1d(i,k) = eta1d(i,k) + aup*dv1v_tmp1(i,k)*dellavd(i,k)*g/dp
        edtod(i) = edtod(i) - adw*etad_tmp1(i,k)*dv1v_tmp1(i,k)*dellavd(i,k)*g/dp
        etadd(i,k) = etadd(i,k) - adw*edto_tmp3(i)*dv1v_tmp1(i,k)*dellavd(i,k)*g/dp
        dv1vd = dv1vd + (aup*eta_tmp3(i, k)-adw*edto_tmp3(i)*etad_tmp1(i, k))*dellavd(i,k)*g/dp

        eta1d(i,k-1) = eta1d(i,k-1) - aup*dv3v_tmp1(i,k)*dellavd(i,k)*g/dp
        edtod(i) = edtod(i) + adw*etad_tmp1(i,k-1)*dv3v_tmp1(i,k)*dellavd(i,k)*g/dp
        etadd(i,k-1) = etadd(i,k-1) + adw*edto_tmp3(i)*dv3v_tmp1(i,k)*dellavd(i,k)*g/dp
        dv3vd = dv3vd - (aup*eta_tmp3(i, k-1)-adw*edto_tmp3(i)*etad_tmp1(i, k-1))*dellavd(i,k)*g/dp

        temd = temd - aup*eta_tmp3(i,k-1)*dv2v_tmp1(i,k)*dz_tmp11(i,k)*dellavd(i,k)*g/dp
        eta1d(i,k-1) = eta1d(i,k-1) -aup*tem_tmp6(i,k)*dv2v_tmp1(i,k)*dz_tmp11(i,k)*dellavd(i,k)*g/dp
        edtod(i) = edtod(i) - adw*ptem_tmp6(i,k)*etad_tmp1(i,k)*dv2v_tmp1(i,k)*dz_tmp11(i,k)*dellavd(i,k)*g/dp
        etadd(i,k) = etadd(i,k) - adw*ptem_tmp6(i,k)*edto_tmp3(i)*dz_tmp11(i,k)*dv2v_tmp1(i,k)*dellavd(i,k)*g/dp
        dv2vd = dv2vd - (aup*tem_tmp6(i,k)*eta_tmp3(i,k-1)+adw*edto_tmp3(i)*ptem_tmp6(i,k)*etad_tmp1(i, k)) &
                  *dz_tmp11(i,k)*dellavd(i,k)*g/dp
        dzd = dzd - (aup*tem_tmp6(i,k)*eta_tmp3(i, k-1)+adw*edto_tmp3(i)*ptem_tmp6(i,k)*etad_tmp1(i, k)) &
                  *dv2v_tmp1(i,k)*dellavd(i,k)*g/dp

        tem1d = tem1d + aup*half*dz_tmp11(i,k)*eta_tmp3(i,k-1)*(vcko_tmp2(i, k)+vcko_tmp2(i, k-1))* dellavd(i,k)*g/dp
        dzd = dzd + aup*half*tem1_tmp6(i,k)*eta_tmp3(i,k-1) *(vcko_tmp2(i, k)+vcko_tmp2(i, k-1))* dellavd(i,k)*g/dp
        eta1d(i,k-1) = eta1d(i,k-1) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*(vcko_tmp2(i, k)+vcko_tmp2(i, k-1))* &
                       dellavd(i,k)*g/dp
        vckod(i,k) = vckod(i,k) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*eta_tmp3(i, k-1)*dellavd(i,k)*g/dp
        vckod(i,k-1) = vckod(i,k-1) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*eta_tmp3(i, k-1)*dellavd(i,k)*g/dp

        edtod(i) = edtod(i) + adw*half*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*(vcdo_tmp1(i, k) &
                   +vcdo_tmp1(i, k-1))* dellavd(i,k)*g/dp
        etadd(i,k) = etadd(i,k)  + adw*half*edto_tmp3(i)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*(vcdo_tmp1(i, k) &
                   +vcdo_tmp1(i, k-1))* dellavd(i,k)*g/dp
        ptem1d = ptem1d + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*dz_tmp11(i,k)*(vcdo_tmp1(i, k) &
                  +vcdo_tmp1(i, k-1))* dellavd(i,k)*g/dp
        dzd = dzd + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*(vcdo_tmp1(i, k)+vcdo_tmp1(i, k-1))* dellavd(i,k)*g/dp
        vcdod(i,k) = vcdod(i,k) + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*dellavd(i,k)*g/dp
        vcdod(i,k-1) = vcdod(i,k-1) + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*dellavd(i,k)*g/dp

        eta1d(i,k-1) = eta1d(i,k-1) - pgcon*aup*(dv1u_tmp1(i,k)-dv3u_tmp1(i,k))*dellaud(i,k)*g/dp
        edtod(i) = edtod(i) + pgcon*adw*etad_tmp1(i,k)*(dv1u_tmp1(i,k)-dv3u_tmp1(i,k))*dellaud(i,k)*g/dp
        etadd(i,k) = etadd(i,k) + pgcon*adw*edto_tmp3(i)*(dv1u_tmp1(i,k)-dv3u_tmp1(i,k))*dellaud(i,k)*g/dp
        dv1ud = dv1ud - pgcon*(aup*eta_tmp3(i, k-1)-adw*edto_tmp3(i)*etad_tmp1(i, k))*dellaud(i,k)*g/dp
        dv3ud = dv3ud + pgcon*(aup*eta_tmp3(i, k-1)-adw*edto_tmp3(i)*etad_tmp1(i, k))*dellaud(i,k)*g/dp

        eta1d(i,k) = eta1d(i,k) + aup*dv1u_tmp1(i,k)*dellaud(i,k)*g/dp
        edtod(i) = edtod(i) - adw*etad_tmp1(i,k)*dv1u_tmp1(i,k)*dellaud(i,k)*g/dp
        etadd(i,k) = etadd(i,k) - adw*edto_tmp3(i)*dv1u_tmp1(i,k)*dellaud(i,k)*g/dp
        dv1ud = dv1ud + (aup*eta_tmp3(i, k)-adw*edto_tmp3(i)*etad_tmp1(i, k))*dellaud(i,k)*g/dp

        eta1d(i,k-1) = eta1d(i,k-1) - aup*dv3u_tmp1(i,k)*dellaud(i,k)*g/dp
        edtod(i) = edtod(i) + adw*etad_tmp1(i,k-1)*dv3u_tmp1(i,k)*dellaud(i,k)*g/dp
        etadd(i,k-1) = etadd(i,k-1) + adw*edto_tmp3(i)*dv3u_tmp1(i,k)*dellaud(i,k)*g/dp
        dv3ud = dv3ud - (aup*eta_tmp3(i, k-1)-adw*edto_tmp3(i)*etad_tmp1(i, k-1))*dellaud(i,k)*g/dp

        temd = temd - aup*eta_tmp3(i,k-1)*dv2u_tmp1(i,k)*dz_tmp11(i,k)*dellaud(i,k)*g/dp
        eta1d(i,k-1) = eta1d(i,k-1) -aup*tem_tmp6(i,k)*dv2u_tmp1(i,k)*dz_tmp11(i,k)*dellaud(i,k)*g/dp
        edtod(i) = edtod(i) - adw*ptem_tmp6(i,k)*etad_tmp1(i,k)*dv2u_tmp1(i,k)*dz_tmp11(i,k)*dellaud(i,k)*g/dp
        etadd(i,k) = etadd(i,k) - adw*ptem_tmp6(i,k)*edto_tmp3(i)*dz_tmp11(i,k)*dv2u_tmp1(i,k)*dellaud(i,k)*g/dp
        dv2ud = dv2ud - (aup*tem_tmp6(i,k)*eta_tmp3(i, k-1)+adw*edto_tmp3(i)*ptem_tmp6(i,k)*etad_tmp1(i, k)) &
                  *dz_tmp11(i,k)*dellaud(i,k)*g/dp
        dzd = dzd - (aup*tem_tmp6(i,k)*eta_tmp3(i, k-1)+adw*edto_tmp3(i)*ptem_tmp6(i,k)*etad_tmp1(i, k)) &
                  *dv2u_tmp1(i,k)*dellaud(i,k)*g/dp

        tem1d = tem1d + aup*half*dz_tmp11(i,k)*eta_tmp3(i,k-1)*(ucko_tmp2(i, k)+ucko_tmp2(i, k-1))* dellaud(i,k)*g/dp
        dzd = dzd + aup*half*tem1_tmp6(i,k)*eta_tmp3(i,k-1) *(ucko_tmp2(i, k)+ucko_tmp2(i, k-1))* dellaud(i,k)*g/dp
        eta1d(i,k-1) = eta1d(i,k-1) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*(ucko_tmp2(i, k)+ucko_tmp2(i, k-1))* &
                       dellaud(i,k)*g/dp
        uckod(i,k) = uckod(i,k) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*eta_tmp3(i, k-1)*dellaud(i,k)*g/dp
        uckod(i,k-1) = uckod(i,k-1) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*eta_tmp3(i, k-1)*dellaud(i,k)*g/dp

        edtod(i) = edtod(i) + adw*half*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*(ucdo_tmp1(i, k) &
                   +ucdo_tmp1(i, k-1))* dellaud(i,k)*g/dp
        etadd(i,k) = etadd(i,k)  + adw*half*edto_tmp3(i)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*(ucdo_tmp1(i, k) &
                   +ucdo_tmp1(i, k-1))* dellaud(i,k)*g/dp
        ptem1d = ptem1d + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*dz_tmp11(i,k)*(ucdo_tmp1(i, k) &
                   +ucdo_tmp1(i, k-1))* dellaud(i,k)*g/dp
        dzd = dzd + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*(ucdo_tmp1(i, k)+ucdo_tmp1(i, k-1))* dellaud(i,k)*g/dp
        ucdod(i,k) = ucdod(i,k) + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*dellaud(i,k)*g/dp
        ucdod(i,k-1) = ucdod(i,k-1) + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*dellaud(i,k)*g/dp

        eta1d(i,k) = eta1d(i,k) + aup*dv1q_tmp1(i,k)*dellaqd(i,k)*g/dp
        edtod(i) = edtod(i) - adw*etad_tmp1(i,k)*dv1q_tmp1(i,k)*dellaqd(i,k)*g/dp
        etadd(i,k) = etadd(i,k) - adw*edto_tmp3(i)*dv1q_tmp1(i,k)*dellaqd(i,k)*g/dp
        dv1qd = dv1qd + (aup*eta_tmp3(i, k)-adw*edto_tmp3(i)*etad_tmp1(i, k))*dellaqd(i,k)*g/dp


        eta1d(i,k-1) = eta1d(i,k-1) - aup*dv3q_tmp1(i,k)*dellaqd(i,k)*g/dp
        edtod(i) = edtod(i) + adw*etad_tmp1(i,k-1)*dv3q_tmp1(i,k)*dellaqd(i,k)*g/dp
        etadd(i,k-1) = etadd(i,k-1) + adw*edto_tmp3(i)*dv3q_tmp1(i,k)*dellaqd(i,k)*g/dp
        dv3qd = dv3qd - (aup*eta_tmp3(i, k-1)-adw*edto_tmp3(i)*etad_tmp1(i, k-1))*dellaqd(i,k)*g/dp

        temd = temd - aup*eta_tmp3(i,k-1)*dv2q_tmp1(i,k)*dz_tmp11(i,k)*dellaqd(i,k)*g/dp
        eta1d(i,k-1) = eta1d(i,k-1) -aup*tem_tmp6(i,k)*dv2q_tmp1(i,k)*dz_tmp11(i,k)*dellaqd(i,k)*g/dp
        edtod(i) = edtod(i) - adw*ptem_tmp6(i,k)*etad_tmp1(i,k)*dv2q_tmp1(i,k)*dz_tmp11(i,k)*dellaqd(i,k)*g/dp
        etadd(i,k) = etadd(i,k) - adw*ptem_tmp6(i,k)*edto_tmp3(i)*dz_tmp11(i,k)*dv2q_tmp1(i,k)*dellaqd(i,k)*g/dp
        dv2qd = dv2qd - (aup*tem_tmp6(i,k)*eta_tmp3(i, k-1)+adw*edto_tmp3(i)*ptem_tmp6(i,k)*etad_tmp1(i, k)) &
                  *dz_tmp11(i,k)*dellaqd(i,k)*g/dp
        dzd = dzd - (aup*tem_tmp6(i,k)*eta_tmp3(i, k-1)+adw*edto_tmp3(i)*ptem_tmp6(i,k)*etad_tmp1(i, k)) &
                  *dv2q_tmp1(i,k)*dellaqd(i,k)*g/dp

        tem1d = tem1d + aup*half*dz_tmp11(i,k)*eta_tmp3(i,k-1)*(qcko_tmp5(i, k)+qcko_tmp5(i, k-1))* dellaqd(i,k)*g/dp
        dzd = dzd + aup*half*tem1_tmp6(i,k)*eta_tmp3(i,k-1) *(qcko_tmp5(i, k)+qcko_tmp5(i, k-1))* dellaqd(i,k)*g/dp
        eta1d(i,k-1) = eta1d(i,k-1) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*(qcko_tmp5(i, k)+qcko_tmp5(i, k-1))* &
                       dellaqd(i,k)*g/dp
        qckod(i,k) = qckod(i,k) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*eta_tmp3(i, k-1)*dellaqd(i,k)*g/dp
        qckod(i,k-1) = qckod(i,k-1) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*eta_tmp3(i, k-1)*dellaqd(i,k)*g/dp

        edtod(i) = edtod(i) + adw*half*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*(qrcdo_tmp2(i, k) &
                   +qrcdo_tmp2(i, k-1))* dellaqd(i,k)*g/dp
        etadd(i,k) = etadd(i,k)  + adw*half*edto_tmp3(i)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*(qrcdo_tmp2(i, k) &
                   +qrcdo_tmp2(i, k-1))* dellaqd(i,k)*g/dp
        ptem1d = ptem1d + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*dz_tmp11(i,k)*(qrcdo_tmp2(i, k) &
                    +qrcdo_tmp2(i, k-1))* dellaqd(i,k)*g/dp
        dzd = dzd + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*(qrcdo_tmp2(i, k) &
                    +qrcdo_tmp2(i, k-1))* dellaqd(i,k)*g/dp
        qrcdod(i,k) = qrcdod(i,k) + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*dellaqd(i,k)*g/dp
        qrcdod(i,k-1) = qrcdod(i,k-1) + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*dellaqd(i,k)*g/dp

        eta1d(i,k) = eta1d(i,k) + aup*dv1h_tmp1(i,k)*dellahd(i,k)*g/dp
        edtod(i) = edtod(i) - adw*etad_tmp1(i,k)*dv1h_tmp1(i,k)*dellahd(i,k)*g/dp
        etadd(i,k) = etadd(i,k) - adw*edto_tmp3(i)*dv1h_tmp1(i,k)*dellahd(i,k)*g/dp
        dv1hd = dv1hd + (aup*eta(i, k)-adw*edto_tmp3(i)*etad_tmp1(i, k))*dellahd(i,k)*g/dp

        eta1d(i,k-1) = eta1d(i,k-1) - aup*dv3h_tmp1(i,k)*dellahd(i,k)*g/dp
        edtod(i) = edtod(i) + adw*etad_tmp1(i,k-1)*dv3h_tmp1(i,k)*dellahd(i,k)*g/dp
        etadd(i,k-1) = etadd(i,k-1) + adw*edto_tmp3(i)*dv3h_tmp1(i,k)*dellahd(i,k)*g/dp
        dv3hd = dv3hd - (aup*eta(i, k-1)-adw*edto_tmp3(i)*etad_tmp1(i, k-1))*dellahd(i,k)*g/dp

        temd = temd - aup*eta(i,k-1)*dv2h_tmp1(i,k)*dz_tmp11(i,k)*dellahd(i,k)*g/dp
        eta1d(i,k-1) = eta1d(i,k-1) -aup*tem_tmp6(i,k)*dv2h_tmp1(i,k)*dz_tmp11(i,k)*dellahd(i,k)*g/dp
        edtod(i) = edtod(i) - adw*ptem_tmp6(i,k)*etad_tmp1(i,k)*dv2h_tmp1(i,k)*dz_tmp11(i,k)*dellahd(i,k)*g/dp
        etadd(i,k) = etadd(i,k) - adw*ptem_tmp6(i,k)*edto_tmp3(i)*dz_tmp11(i,k)*dv2h_tmp1(i,k)*dellahd(i,k)*g/dp
        dv2hd = dv2hd - (aup*tem_tmp6(i,k)*eta_tmp3(i, k-1)+adw*edto_tmp3(i)*ptem_tmp6(i,k)*etad_tmp1(i, k)) &
                  *dz_tmp11(i,k)*dellahd(i,k)*g/dp
        dzd = dzd - (aup*tem_tmp6(i,k)*eta_tmp3(i, k-1)+adw*edto_tmp3(i)*ptem_tmp6(i,k)*etad_tmp1(i, k)) &
                  *dv2h_tmp1(i,k)*dellahd(i,k)*g/dp
        tem1d = tem1d + aup*half*dz_tmp11(i,k)*eta_tmp3(i,k-1)*(hcko_tmp2(i, k)+hcko_tmp2(i, k-1))* dellahd(i,k)*g/dp
        dzd = dzd + aup*half*tem1_tmp6(i,k)*eta_tmp3(i,k-1) *(hcko_tmp2(i, k)+hcko_tmp2(i, k-1))* dellahd(i,k)*g/dp
        eta1d(i,k-1) = eta1d(i,k-1) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*(hcko_tmp2(i, k)+hcko_tmp2(i, k-1))* &
                       dellahd(i,k)*g/dp
        hckod(i,k) = hckod(i,k) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*eta_tmp3(i, k-1)*dellahd(i,k)*g/dp
        hckod(i,k-1) = hckod(i,k-1) + aup*half*tem1_tmp6(i,k)*dz_tmp11(i,k)*eta_tmp3(i, k-1)*dellahd(i,k)*g/dp

        edtod(i) = edtod(i) + adw*half*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*(hcdo_tmp1(i, k) &
                   +hcdo_tmp1(i, k-1))* dellahd(i,k)*g/dp
        etadd(i,k) = etadd(i,k)  + adw*half*edto_tmp3(i)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*(hcdo_tmp1(i, k) &
                   +hcdo_tmp1(i, k-1))* dellahd(i,k)*g/dp
        ptem1d = ptem1d + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*dz_tmp11(i,k)*(hcdo_tmp1(i, k) &
                   +hcdo_tmp1(i, k-1))* dellahd(i,k)*g/dp
        dzd = dzd + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*(hcdo_tmp1(i, k)+hcdo_tmp1(i, k-1))* dellahd(i,k)*g/dp
        hcdod(i,k) = hcdod(i,k) + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*dellahd(i,k)*g/dp
        hcdod(i,k-1) = hcdod(i,k-1) + adw*half*edto_tmp3(i)*etad_tmp1(i,k)*ptem1_tmp3(i,k)*dz_tmp11(i,k)*dellahd(i,k)*g/dp

        IF (k .LE. kbcon(i)) THEN
          xlamd0d(i) = xlamd0d(i) + ptem1d
          ptem1d = zero
        ELSE
          ptem1d = zero
        END IF

        xlamudd(i) = xlamudd(i) + tem1d
        tem1d = zero

        xlamued(i,k) = xlamued(i,k) + half*temd
        xlamued(i,k-1) = xlamued(i,k-1) + half*temd
        temd = zero

        vod(i,k-1) = vod(i,k-1) + dv3vd
        dv3vd = zero

        vod(i,k) = vod(i,k) + half*dv2vd
        vod(i,k-1) = vod(i,k-1) + half*dv2vd
        dv2vd=zero

        vod(i,k) = vod(i,k) + dv1vd
        dv1vd = zero

        uod(i,k-1) = uod(i,k-1) + dv3ud
        dv3ud = zero

        uod(i,k) = uod(i,k) + half*dv2ud
        uod(i,k-1) = uod(i,k-1) + half*dv2ud
        dv2ud = zero

        uod(i,k) = uod(i,k) + dv1ud
        dv1ud =zero

        qod(i,k-1) = qod(i,k-1) +  dv3qd
        dv3qd = zero

        qod(i,k) = qod(i,k) + half*dv2qd
        qod(i,k-1) = qod(i,k-1) + half*dv2qd
        dv2qd = zero

        qod(i,k) = qod(i,k) + dv1qd
        dv1qd = zero

        heod(i,k-1) = heod(i,k-1) + dv3hd
        dv3hd = zero

        heod(i,k) = heod(i,k) + half*dv2hd
        heod(i,k-1) = heod(i,k-1) + half*dv2hd
        dv2hd = zero

        heod(i,k) = heod(i,k) + dv1hd
        dv1hd = zero

        zid(i,k) = zid(i,k) + dzd
        zid(i,k-1) = zid(i,k-1) - dzd
        dzd = zero
      END IF
    END DO
  END DO


  DO i=1,im
    dp = 1000.0_r_kind*del(i, 1)*ps(i)

    edtod(i) = edtod(i) +etad_tmp1(i, 1)*(vcdo_tmp1(i, 1)-vo_tmp1(i, 1))*dellavd(i,1)*g/dp
    etadd(i,1) = etadd(i,1) + edto_tmp3(i)*dellavd(i,1)*(vcdo_tmp1(i, 1)-vo_tmp1(i, 1))*g/dp
    vcdod(i,1) = vcdod(i,1) + edto_tmp3(i)*etad_tmp1(i,1)*dellavd(i,1)*g/dp
    vod(i,1) = vod(i,1) - edto_tmp3(i)*etad_tmp1(i, 1)*dellavd(i,1)*g/dp
    dellavd(i,1) = zero

    edtod(i) = edtod(i) + etad_tmp1(i,1)*dellaud(i,1)*(ucdo_tmp1(i, 1)-uo_tmp1(i, 1))*g/dp
    etadd(i,1) = etadd(i,1) + edto_tmp3(i)*dellaud(i,1)*(ucdo_tmp1(i, 1)-uo_tmp1(i, 1))*g/dp
    ucdod(i,1) = ucdod(i,1) + edto_tmp3(i)*etad_tmp1(i,1)*dellaud(i,1)*g/dp
    uod(i,1) = uod(i,1) - edto_tmp3(i)*etad_tmp1(i, 1)*dellaud(i,1)*g/dp
    dellaud(i,1) = zero

    edtod(i) = edtod(i) + etad_tmp1(i,1)*dellaqd(i,1)*(qcdo_tmp2(i, 1)-qo_tmp4(i, 1))*g/dp
    etadd(i,1) = etadd(i,1) + edto_tmp3(i)*dellaqd(i,1)*(qcdo_tmp2(i, 1)-qo_tmp4(i, 1))*g/dp
    qcdod(i,1) = qcdod(i,1) + edto_tmp3(i)*etad_tmp1(i,1)*dellaqd(i,1)*g/dp
    qod(i,1) = qod(i,1) - edto_tmp3(i)*etad_tmp1(i, 1)*dellaqd(i,1)*g/dp
    dellaqd(i,1) = zero

    edtod(i) = edtod(i) + etad(i,1)*dellahd(i,1)*(hcdo_tmp1(i, 1)-heo_tmp2(i, 1))*g/dp
    etadd(i,1) = etadd(i,1) + edto(i)*dellahd(i,1)*(hcdo_tmp1(i, 1)-heo_tmp2(i, 1))*g/dp
    hcdod(i,1) = hcdod(i,1) + edto(i)*etad_tmp1(i,1)*dellahd(i,1)*g/dp
    heod(i,1) = heod(i,1) - edto(i)*etad_tmp1(i, 1)*dellahd(i,1)*g/dp
    dellahd(i,1) = zero

  END DO

  DO i=1,im
     cnvscl_5d(i) = cnvscl_5d(i)+cnvscl_6_tmp1(i)*cnvscl_6d(i)
     cnvscl_6d(i) = cnvscl_5(i)*cnvscl_6d(i)

     aa1d(i) = aa1d(i) + half*(one-TANH(aa1_tmp2(i))**2)*cnvscl_6d(i)
     cnvscl_6d(i) = zero
  END DO


  DO k=1,km1
    DO i=1,im
      IF (k .LT. jmin(i)) THEN

         dzd = zero
         max2d = zero
         dhhd=zero
         dhd = zero
         dgd = zero
         dtd = zero
         gammad = zero

         val=0.
         IF (val .LT. qeso_tmp4(i, k) - qo_tmp4(i, k)) THEN
            max2 = qeso_tmp4(i, k) - qo_tmp4(i, k)
         ELSE
            max2 = val
         END IF

         edtod(i) = edtod(i) + g*delta*dz_tmp10(i,k)*max2*aa1d(i)
         dzd = dzd + g*delta*edto_tmp3(i)*max2*aa1d(i)
         max2d = max2d +  g*delta*edto_tmp3(i)*dz_tmp10(i,k)*aa1d(i)

         IF (val .LT. qeso_tmp4(i, k) - qo_tmp4(i, k)) THEN
           qesod(i,k) = qesod(i,k) + max2d
           qod(i,k) = qod(i,k) - max2d
           max2d = zero
         ELSE
           max2d= zero
         END IF


         edtod(i) = edtod(i) + dz_tmp10(i,k)*(g/(cp*dt_tmp1(i,k)))*((dhh_tmp1(i,k)-dh_tmp1(i,k)) &
                   /(one+dg_tmp1(i,k)))*(one+delta*cp*dg_tmp1(i,k)*dt_tmp1(i,k)/hvap)*aa1d(i)
         dzd = dzd +  edto_tmp3(i)*aa1d(i)*(g/(cp*dt_tmp1(i,k)))*((dhh_tmp1(i,k)-dh_tmp1(i,k)) &
               /(one+dg_tmp1(i,k)))*(one+delta*cp*dg_tmp1(i,k)*dt_tmp1(i,k)/hvap)
         dtd = dtd  + edto_tmp3(i)*dz_tmp10(i,k)*(-g*aa1d(i)/(cp*dt_tmp1(i,k)**2))*((dhh_tmp1(i,k)-dh_tmp1(i,k)) &
             /(one+dg_tmp1(i,k)))*(one+delta*cp*dg_tmp1(i,k)*dt_tmp1(i,k)/hvap)
         dhhd = dhhd + edto_tmp3(i)*dz_tmp10(i,k)*(g/(cp*dt_tmp1(i,k)))*(aa1d(i)/ &
             (one+dg_tmp1(i,k)))*(one+delta*cp*dg_tmp1(i,k)*dt_tmp1(i,k)/hvap)
         dhd = dhd + edto_tmp3(i)*dz_tmp10(i,k)*(g/(cp*dt_tmp1(i,k)))*(-aa1d(i)/  &
                 (one+dg_tmp1(i,k)))*(one+delta*cp*dg_tmp1(i,k)*dt_tmp1(i,k)/hvap)
         dgd = dgd + edto_tmp3(i)*dz_tmp10(i,k)*(g/(cp*dt_tmp1(i,k)))*(-(dhh_tmp1(i,k)-dh_tmp1(i,k))*aa1d(i) &
                  /(one+dg_tmp1(i,k))**2)*(one+delta*cp*dg_tmp1(i,k)*dt_tmp1(i,k)/hvap)
         dgd = dgd + edto_tmp3(i)*dz_tmp10(i,k)*(g/(cp*dt_tmp1(i,k)))*((dhh_tmp1(i,k)-dh_tmp1(i,k)) &
             /(one+dg_tmp1(i,k)))*(delta*cp*aa1d(i)*dt_tmp1(i,k)/hvap)
         dtd = dtd + edto_tmp3(i)*dz_tmp10(i,k)*(g/(cp*dt_tmp1(i,k)))*((dhh_tmp1(i,k)-dh_tmp1(i,k))/ &
             (one+dg_tmp1(i,k)))*(delta*cp*dg_tmp1(i,k)*aa1d(i)/hvap)

        zod(i,k+1) = zod(i,k+1) - dzd
        zod(i,k) = zod(i,k) + dzd
        dzd = zero

        hesod(i,k) = hesod(i,k) + dhd
        dhd = zero

        gammad = gammad + dgd
        dgd = zero

        tod(i,k) = tod(i,k) + dtd
        dtd = zero

        hcdod(i,k) = hcdod(i,k) + dhhd
        dhhd  = zero

        qesod(i,k) = qesod(i,k) + el2orc*gammad/to_tmp2(i,k)**2
        tod(i,k) = tod(i,k) - el2orc*qeso_tmp4(i,k)*2*gammad/to_tmp2(i,k)**3
        gammad = zero

      END IF
    END DO
  END DO


  DO i=1,im
    edtmax = edtmaxl
    edto_tem1d  = zero
    IF (slimsk(i) .EQ. zero) edtmax = edtmaxs
    IF (pwevo_tmp1(i) .LT. zero) THEN
      IF (edto_tmp2(i) .GT. edtmax) THEN
        edtod(i) = zero
      END IF
       edto_tem1d = edto_tem1d + edtod(i)
       edtod(i) = zero

       edtod(i) = edtod(i) - edto_tem1d*pwavo(i)/pwevo(i)
       pwavod(i) = pwavod(i) - edto_tmp1(i)*edto_tem1d/pwevo(i)
       pwevod(i)= pwevod(i) + edto_tmp1(i)*pwavo(i)*edto_tem1d/pwevo(i)**2
       edto_tem1d = zero
    ELSE
      edtod(i) = zero
    END IF
  END DO




  DO k=1,km1
    DO i=1,im
      IF (k .LT. jmin(i)) THEN
        gammad = zero
        temd = zero
        tem1d = zero
        factord = zero
        dzd = zero

        pwdod(i,k) = pwdod(i,k) + pwevod(i)

        qrcdod(i,k) = qrcdod(i,k) + qcdod(i,k)
        qcdod(i,k) = zero

        etadd(i,k+1) = etadd(i,k+1) + (qcdo_tmp1(i, k)-qrcdo_tmp1(i, k))*pwdod(i,k)
        qcdod(i,k) = qcdod(i,k) + etad(i,k+1)*pwdod(i,k)
        qrcdod(i,k) = qrcdod(i,k) - etad(i,k+1)*pwdod(i,k)
        pwdod(i,k) = zero

        qcdod(i,k+1)=  qcdod(i,k+1) + (one-tem1_tmp5(i,k))*qcdod(i,k)/factor_tmp5(i,k)
        tem1d = tem1d - qcdo_tmp2(i, k+1)*qcdod(i,k)/factor_tmp5(i,k)
        temd = temd + half*(qo_tmp4(i, k)+qo_tmp4(i, k+1))*qcdod(i,k)/factor_tmp5(i,k)
        qod(i,k) = qod(i,k) + half*tem_tmp5(i,k)*qcdod(i,k)/factor_tmp5(i,k)
        qod(i,k+1) = qod(i,k+1) + half*tem_tmp5(i,k)*qcdod(i,k)/factor_tmp5(i,k)
        factord = factord - ((one-tem1_tmp5(i,k))*qcdo_tmp2(i, k+1) &
             +tem_tmp5(i,k)*half*(qo_tmp4(i, k)+qo_tmp4(i, k+1)))*qcdod(i,k)/factor_tmp5(i,k)**2
        qcdod(i,k) = zero

        temd = temd + factord
        tem1d = tem1d - factord
        factord = zero

        IF (k .GE. kbcon(i)) THEN
            dzd = dzd + half*xlamdd*tem1d
            tem1d = zero

            dzd = dzd + xlamde*temd
            temd = zero
        ELSE
            xlamd0d(i) = xlamd0d(i) + half*dz_tmp9(i,k)*tem1d
            dzd = dzd + half*(xlamd(i)+xlamdd)*tem1d
            tem1d = zero

            dzd = dzd + xlamde*temd
            temd = zero
        END IF

        zid(i,k+1) = zid(i,k+1) + dzd
        zid(i,k) = zid(i,k) - dzd
        dzd = zero

        qesod(i,k) = qesod(i,k) + qrcdod(i,k)
        gammad = gammad + dbyo_tmp2(i,k)*qrcdod(i,k) /((one+gamma_tmp5(i,k))*hvap)
        gammad = gammad - gamma_tmp5(i,k)*dbyo_tmp2(i,k)*qrcdod(i,k) &
                /(((one+gamma_tmp5(i,k))**2)*hvap)
        dbyod(i,k) = dbyod(i,k) + one/hvap*(gamma_tmp5(i,k)/(one+gamma_tmp5(i,k)))*qrcdod(i,k)
        qrcdod(i,k) = zero

        qesod(i,k) = qesod(i,k) + el2orc*gammad/to_tmp2(i,k)**2
        tod(i,k) = tod(i,k) - el2orc*qeso_tmp4(i,k)*2*gammad/to_tmp2(i,k)**3
        gammad = zero
      END IF
    END DO
  END DO


  DO k=1,km1
    DO i=1,im
      IF (k .LT. jmin(i)) THEN
       ptem1d = zero
       ptemd = zero
       factord = zero
       tem1d = zero
       temd = zero
       dzd = zero

       hcdod(i,k) = hcdod(i,k) + dbyod(i,k)
       hesod(i,k) = hesod(i,k) - dbyod(i,k)
       dbyod(i,k) = zero

       vcdod(i,k+1) = vcdod(i,k+1) + (one-tem1_tmp4(i,k))*vcdod(i,k)/factor_tmp4(i,k)
       tem1d = tem1d - vcdo_tmp1(i, k+1)*vcdod(i,k)/factor_tmp4(i,k)
       ptemd = ptemd + vo_tmp1(i,k+1)*vcdod(i,k)/factor_tmp4(i,k)
       vod(i,k+1) = vod(i,k+1)  + ptem_tmp5(i,k)*vcdod(i,k)/factor_tmp4(i,k)
       ptem1d = ptem1d + vo_tmp1(i,k)*vcdod(i,k)/factor_tmp4(i,k)
       vod(i,k) = vod(i,k) + ptem1_tmp2(i,k)*vcdod(i,k)/factor_tmp4(i,k)
       factord = factord -  ((one-tem1_tmp4(i,k))*vcdo_tmp1(i,k+1)+ptem_tmp5(i,k)*vo_tmp1(i,k+1) &
            +ptem1_tmp2(i,k)*vo_tmp1(i,k))*vcdod(i,k)/factor_tmp4(i,k)**2
       vcdod(i,k) = zero

       ucdod(i,k+1) = ucdod(i,k+1) + (one-tem1_tmp4(i,k))*ucdod(i,k)/factor_tmp4(i,k)
       tem1d = tem1d - ucdo_tmp1(i, k+1)*ucdod(i,k)/factor_tmp4(i,k)
       ptemd = ptemd + uo_tmp1(i,k+1)*ucdod(i,k)/factor_tmp4(i,k)
       uod(i,k+1) = uod(i,k+1)  + ptem_tmp5(i,k)*ucdod(i,k)/factor_tmp4(i,k)
       ptem1d = ptem1d + uo_tmp1(i,k)*ucdod(i,k)/factor_tmp4(i,k)
       uod(i,k) = uod(i,k) + ptem1_tmp2(i,k)*ucdod(i,k)/factor_tmp4(i,k)
       factord = factord -  ((one-tem1_tmp4(i,k))*ucdo_tmp1(i,k+1)+ptem_tmp5(i,k)*uo_tmp1(i,k+1) &
                  +ptem1_tmp2(i,k)*uo_tmp1(i,k))*ucdod(i,k)/factor_tmp4(i,k)**2
       ucdod(i,k) = zero

       hcdod(i,k+1) = hcdod(i,k+1) + (one-tem1_tmp4(i,k))*hcdod(i,k)/factor_tmp4(i,k)
       tem1d = tem1d - hcdo_tmp1(i, k+1)*hcdod(i,k)/factor_tmp4(i,k)
       temd = temd + half*(heo_tmp2(i, k)+heo_tmp2(i, k+1))*hcdod(i,k)/factor_tmp4(i,k)
       heod(i,k) = heod(i,k) + half*tem_tmp4(i,k)*hcdod(i,k)/factor_tmp4(i,k)
       heod(i,k+1) = heod(i,k+1) + half*tem_tmp4(i,k)*hcdod(i,k)/factor_tmp4(i,k)
       factord = factord - ((one-tem1_tmp4(i,k))*hcdo_tmp1(i, k+1) &
             +tem_tmp4(i,k)*half*(heo_tmp2(i, k)+heo_tmp2(i, k+1)))*hcdod(i,k) &
                /factor_tmp4(i,k)**2
       hcdod(i,k) = zero

       temd = temd + half*ptem1d
       ptem1d = zero

       temd = temd + half*ptemd
       ptemd = zero

       temd = temd + factord
       tem1d = tem1d - factord
       factord = zero


        IF (k .GE. kbcon(i)) THEN
          dzd = dzd + half*xlamdd*tem1d
          tem1d = zero

          dzd = dzd + xlamde*temd
          temd = zero
        ELSE
          xlamd0d(i) = xlamd0d(i) + half*dz_tmp8(i,k)*tem1d
          dzd = dzd + half*(xlamd(i)+xlamdd)*tem1d
          tem1d = zero

          dzd = dzd + xlamde*temd
          temd = zero
        END IF

        zid(i,k+1) = zid(i,k+1) + dzd
        zid(i,k) = zid(i,k) - dzd
        dzd = zero
      END IF
    END DO
  END DO

  DO i=1,im
    jmn = jmin(i)
    vod(i,jmn) = vod(i,jmn) + vcdod(i,jmn)
    vcdod(i,jmn) = zero
    uod(i,jmn) = uod(i,jmn) + ucdod(i,jmn)
    ucdod(i,jmn) = zero
    qesod(i,jmn) = qesod(i,jmn) + qrcdod(i,jmn)
    qrcdod(i,jmn) = zero
    qod(i,jmn) = qod(i,jmn) + qcdod(i,jmn)
    qcdod(i,jmn) = zero
    heod(i,jmn) = heod(i,jmn) + hcdod(i,jmn)
    hcdod(i,jmn) = zero
  END DO

  DO k=1,km1
    DO i=1,im
      dzd = zero
      ptemd = zero
      IF (k .LE. kmax(i) - 1) THEN
        IF (k .LT. jmin(i) .AND. k .GE. kbcon(i)) THEN
          etadd(i,k+1) = etadd(i,k+1) + (one-ptem_tmp4(i,k)*dz_tmp7(i,k)) * etadd(i,k)
          dzd = dzd - etad_tmp1(i, k+1)*ptem_tmp4(i,k)*etadd(i,k)
          etadd(i,k) = zero

          zid(i,k+1) = zid(i,k+1) + dzd
          zid(i,k) = zid(i,k) - dzd
          dzd = zero
        ELSE IF (k .LT. kbcon(i)) THEN
          etadd(i,k+1) = etadd(i,k+1) + (one-ptem_tmp4(i,k)*dz_tmp7(i,k))*etadd(i,k)
          ptemd = ptemd - etad_tmp1(i,k+1)*dz_tmp7(i,k)*etadd(i,k)
          dzd = dzd - etad_tmp1(i,k+1)*ptem_tmp4(i,k)*etadd(i,k)
          etadd(i,k) = zero

          xlamd0d(i) = xlamd0d(i) + ptemd
          ptemd = zero

          zid(i,k+1) = zid(i,k+1) + dzd
          zid(i,k) = zid(i,k) - dzd
          dzd = zero
        END IF
      END IF
    END DO
  END DO


  DO i=1,im
    beta = betas
    IF (slimsk(i) .EQ. one) beta = betal
    tem = one/FLOAT(kbcon(i))

    dzd = zero

    dzd = dzd - (one-beta**tem)*xlamd0d(i) /dz_tmp6(i,1)**2
    xlamd0d(i) = zero

    sumxd(i) = sumxd(i) + dzd/FLOAT(kbcon(i))
    zid(i,1) = zid(i,1) + dzd/FLOAT(kbcon(i))
    dzd = zero
  END DO

  DO k=km1,1,-1
    DO i=1,im
      dzd = zero
      IF (k .GE. 1 .AND. k .LT. kbcon(i)) THEN
          dzd = dzd + sumxd(i)

          zid(i,k+1) = zid(i,k+1) + dzd
          zid(i,k) = zid(i,k) - dzd
          dzd = zero
      END IF
    END DO
  END DO


  DO i=1,im
     e1d = zero

     edtd(i) =edtd(i) + edtxd(i)
     edtx(i) = zero

     edtd(i) = edtd(i) + edtod(i)
     edtod(i) = zero

     IF (edt_tmp2(i) .LT. zero) THEN
       edtd(i) = zero
     ENDIF

     IF (edt_tmp1(i) .GT. 0.9_r_kind) THEN
      edtd(i) = zero
     END IF

    e1d = e1d - edtd(i)
    edtd(i) = zero

    vsheard(i) = vsheard(i) + .0953_r_kind*2*vshear_tmp2(i)*e1d
    vsheard(i) = vsheard(i) - 0.639_r_kind*e1d
    vsheard(i) = vsheard(i) - .00496_r_kind*3*vshear_tmp2(i)**2*e1d
    e1d = zero

    vshear_tmp2d(i)= vshear_tmp2d(i) + vsheard(i)
    vsheard(i) = zero

    vsheard(i) = vsheard(i) +  1.e3_r_kind*vshear_tmp2d(i)/(zi(i, ktcon_new(i))-zi(i, kb(i)))
    zid(i,ktcon_new(i)) = zid(i,ktcon_new(i))  - 1.e3_r_kind*vshear_tmp1(i)*vshear_tmp2d(i) &
                        /(zi(i, ktcon_new(i))-zi(i, kb(i)))**2
    zid(i,kb(i)) = zid(i,kb(i)) + 1.e3_r_kind*vshear_tmp1(i)*vshear_tmp2d(i)  &
                       /(zi(i, ktcon_new(i))-zi(i, kb(i)))**2
    vshear_tmp2d(i) = zero
  END DO

  DO k=km,2,-1
    DO i=1,im
      sheard = zero

      IF (k .GT. kb(i) .AND. k .LE. ktcon_new(i)) THEN
          sheard = sheard + vsheard(i)

          uod(i,k) = uod(i,k) + (uo_tmp1(i, k)-uo_tmp1(i, k-1))*sheard/(SQRT((uo_tmp1(i, k)-     &
                    uo_tmp1(i, k-1))**2 + (vo_tmp1(i, k)-vo_tmp1(i, k-1))**2))
          uod(i,k-1) = uod(i,k-1) - (uo_tmp1(i, k)-uo_tmp1(i, k-1))*sheard/(SQRT((uo_tmp1(i, k)-     &
                    uo_tmp1(i, k-1))**2 + (vo_tmp1(i, k)-vo_tmp1(i, k-1))**2))
          vod(i,k) = vod(i,k) + (vo_tmp1(i, k)-vo_tmp1(i, k-1))*sheard/(SQRT((uo_tmp1(i, k)-     &
                    uo_tmp1(i, k-1))**2 + (vo_tmp1(i, k)-vo_tmp1(i, k-1))**2))
          vod(i,k-1) = vod(i,k-1) - (vo_tmp1(i, k)-vo_tmp1(i, k-1))*sheard/(SQRT((uo_tmp1(i, k)-     &
                    uo_tmp1(i, k-1))**2 + (vo_tmp1(i, k)-vo_tmp1(i, k-1))**2))
          sheard = zero
      END IF
    END DO
  END DO

  IF (ncloud .GT. 0) THEN
    DO i=1,im
      k = ktcon_new(i) - 1

      qrchd = zero
      dqd = zero
      gammad = zero

      IF (dq_tmp3(i,k) .GT. zero) THEN
        qrchd = qrchd + qckod(i,k)
        qckod(i,k) = zero

        dqd = dqd + qlko_ktcond(i)
        qlko_ktcond(i) = zero
      END IF

      qckod(i,k) = qckod(i,k) + dqd
      qrchd = qrchd - dqd
      dqd = zero

      qesod(i,k) = qesod(i,k) + qrchd
      gammad = gammad + dbyo_tmp1(i,k)*qrchd/(hvap*(one+gamma_tmp4(i,k)))
      dbyod(i,k) = dbyod(i,k) + gamma_tmp4(i,k)*qrchd/(hvap*(one+gamma_tmp4(i,k)))
      gammad = gammad - gamma_tmp4(i,k)*dbyo_tmp1(i, k)*hvap*qrchd/(hvap*(one+gamma_tmp4(i,k)))**2
      qrchd = zero

      qesod(i,k) = qesod(i,k) + el2orc*gammad/to_tmp2(i, k)**2
      tod(i,k) = tod(i,k) - el2orc*qeso_tmp4(i,k)*2*gammad/to_tmp2(i,k)**3
      gammad = zero
    END DO
  END IF

  DO k=km1,2,-1
    DO i=1,im
      IF (k .GE. ktcon(i) .AND. k .LT. ktcon1(i)) THEN
         dzd = zero
         gammad = zero
         qrchd = zero
         temd = zero
         tem1d = zero
         factord = zero
         etahd = zero
         qlkd = zero
         dqd = zero


         IF (dq_tmp2(i,k) .GT. zero) THEN
         pwod(i,k) = pwod(i,k) + pwavod(i)

         etahd = etahd + c0*dz_tmp5(i,k)*qlk_tmp2(i,k)*pwod(i,k)
         dzd = dzd + c0*etah_tmp2(i,k)*qlk_tmp2(i,k)*pwod(i,k)
         qlkd = qlkd + c0*etah_tmp2(i,k)*dz_tmp5(i,k)*pwod(i,k)
         pwod(i,k) = zero

         qlkd = qlkd + qckod(i,k)
         qrchd = qrchd + qckod(i,k)
         qckod(i,k) = zero

          IF (ncloud .GT. 0 ) THEN
             etahd = etahd + c1*g*dz_tmp5(i,k)*qlk_tmp2(i,k)*dellald(i,k)/dp_tmp2(i,k)
             dzd = dzd + c1*g*etah_tmp2(i,k)*qlk_tmp2(i,k)*dellald(i,k)/dp_tmp2(i,k)
             qlkd = qlkd + c1*g*etah_tmp2(i,k)*dz_tmp5(i,k)*dellald(i,k)/dp_tmp2(i,k)
             dellald(i,k) = zero

             dqd = dqd + qlkd/(eta_tmp3(i, k)+etah_tmp2(i,k)*(c0+c1)*dz_tmp5(i,k))
             eta1d(i,k) = eta1d(i,k) -dq_tmp2(i,k)*qlkd &
                      /(eta_tmp3(i, k)+etah_tmp2(i,k)*(c0+c1)*dz_tmp5(i,k))**2
             etahd  = etahd - dq_tmp2(i,k)*(c0+c1)*dz_tmp5(i,k)*qlkd &
                     /(eta_tmp3(i, k)+etah_tmp2(i,k)*(c0+c1)*dz_tmp5(i,k))**2
             dzd = dzd -dq_tmp2(i,k)*etah_tmp2(i,k)*(c0+c1)*qlkd &
                     /(eta_tmp3(i, k)+etah_tmp2(i,k)*(c0+c1)*dz_tmp5(i,k))**2
             qlkd = zero
          ELSE
             dqd = dqd + qlkd/(eta_tmp3(i,k)+etah_tmp2(i,k)*c0*dz_tmp5(i,k))
             qlkd = zero
          END IF
            eta1d(i,k) = eta1d(i,k) + half*etahd
            eta1d(i,k-1) = eta1d(i,k-1) + half*etahd
            etahd = zero
         ENDIF


         eta1d(i,k) = eta1d(i,k) + (qcko_tmp3(i,k)-qrch_tmp2(i,k))*dqd
         qckod(i,k) = qckod(i,k) + eta_tmp3(i,k)*dqd
         qrchd = qrchd - eta_tmp3(i,k)*dqd
         dqd = zero

         qckod(i,k-1) = qckod(i,k-1) + (one-tem1_tmp3(i,k))*qckod(i,k)/factor_tmp3(i,k)
         tem1d = tem1d - qcko_tmp4(i,k-1)*qckod(i,k)/factor_tmp3(i,k)
         temd = temd + half*(qo_tmp4(i, k)+qo_tmp4(i, k-1))*qckod(i,k)/factor_tmp3(i,k)
         qod(i,k) = qod(i,k) + half*tem_tmp3(i,k)*qckod(i,k)/factor_tmp3(i,k)
         qod(i,k-1) = qod(i,k-1) + half*tem_tmp3(i,k)*qckod(i,k)/factor_tmp3(i,k)
         factord = factord - ((one-tem1_tmp3(i,k))*qcko_tmp4(i,k-1) &
               + tem_tmp3(i,k)*half*(qo_tmp4(i, k)+qo_tmp4(i, k-1)))*qckod(i,k) &
                /factor_tmp3(i,k)**2
         qckod(i,k) = zero

         temd =temd + factord
         tem1d = tem1d - factord
         factord = zero

         xlamudd(i) = xlamudd(i) + half*dz_tmp5(i,k)*tem1d
         dzd = dzd + half*xlamud(i)*tem1d
         tem1d = zero

         xlamued(i,k) = xlamued(i,k) + half*dz_tmp5(i,k)*temd
         xlamued(i,k-1) = xlamued(i,k-1) + half*dz_tmp5(i,k)*temd
         dzd = dzd + half*(xlamue_tmp4(i, k)+xlamue_tmp4(i, k-1))*temd
         temd = zero


         qesod(i,k) = qesod(i,k) + qrchd
         gammad = gammad + dbyo_tmp1(i,k)*qrchd/(hvap*(one+gamma_tmp3(i,k)))
         dbyod(i,k) = dbyod(i,k) + gamma_tmp3(i,k)*qrchd/(hvap*(one+gamma_tmp3(i,k)))
         gammad = gammad - gamma_tmp3(i,k)*dbyo_tmp1(i, k)*hvap*qrchd/(hvap*(one+gamma_tmp3(i,k)))**2
         qrchd = zero

         qesod(i,k) = qesod(i,k) + el2orc*gammad/to_tmp2(i, k)**2
         tod(i,k) = tod(i,k) -el2orc*qeso_tmp4(i,k)*2*gammad/to_tmp2(i, k)**3
         gammad = zero

         zid(i,k) = zid(i,k) + dzd
         zid(i,k-1) = zid(i,k-1) -dzd
         dzd = zero
      END IF
    END DO
  END DO

  DO i=1,im
    cnvscl_5d(i) = cnvscl_4(i)*cnvscl_5d(i)
    aa1d(i) = aa1d(i) + half*(one-TANH(aa1_tmp1(i))**2)*cnvscl_5d(i)
    cnvscl_5d(i) = zero
  END DO

!  calculate cloud work function

  DO k=km1,2,-1
    DO i=1,im
      IF (k .GE. kbcon(i) .AND. k .LT. ktcon(i)) THEN
           dz1d = zero
           gammad = zero
           rfactd = zero
           max1d = zero
           val= zero

          dz1d = dz1d + g*delta*max1_tmp1(i,k)*aa1d(i)
          max1d = max1d + g*delta*dz1_tmp1(i,k)*aa1d(i)

          IF (val .LT. qeso_tmp4(i, k) - qo_tmp4(i, k)) THEN
             qesod(i,k) = qesod(i,k) + max1d
             qod(i,k) = qod(i,k) - max1d
             max1d = zero
          ELSE
             max1d = zero
          END IF

          dz1d = dz1d + (g/(cp*to_tmp2(i, k)))*dbyo_tmp1(i, k)*aa1d(i) &
                           /(one+gamma_tmp2(i,k))*rfact_tmp1(i,k)
          dbyod(i,k) = dbyod(i,k) + dz1_tmp1(i,k)*(g/(cp*to_tmp2(i, k)))*aa1d(i) &
                        /(one+gamma_tmp2(i,k))*rfact_tmp1(i,k)
          tod(i,k) = tod(i,k) -  dz1_tmp1(i,k)*(g/(cp*to_tmp2(i, k)**2))*dbyo_tmp1(i, k)*aa1d(i) &
                        /(one+gamma_tmp2(i,k))*rfact_tmp1(i,k)
          gammad = gammad - dz1_tmp1(i,k)*(g/(cp*to_tmp2(i, k)))*dbyo_tmp1(i, k)*aa1d(i) &
                      /((one+gamma_tmp2(i,k))**2)*rfact_tmp1(i,k)
          rfactd = rfactd + dz1_tmp1(i,k)*(g/(cp*to_tmp2(i, k)))*dbyo_tmp1(i, k)*aa1d(i) &
                     /(one+gamma_tmp2(i,k))

           gammad = gammad + delta*cp*to_tmp2(i,k)*rfactd/hvap
           tod(i,k) = tod(i,k) + delta*cp*gamma_tmp2(i,k)*rfactd/hvap
           rfactd = zero

           qesod(i,k) = qesod(i,k) + el2orc*gammad/to_tmp2(i, k)**2
           tod(i,k) = tod(i,k) - el2orc*qeso_tmp4(i,k)*2*gammad/to_tmp2(i, k)**3
           gammad = zero

           zod(i,k+1) = zod(i,k+1) + dz1d
           zod(i,k) = zod(i,k) - dz1d
           dz1d = zero
      END IF
    END DO
  END DO

  DO k=km1,2,-1
    DO i=1,im
      IF (k .GT. kb(i) .AND. k .LT. ktcon(i)) THEN
        etahd =  zero
        dzd = zero
        qlkd = zero
        qrchd = zero
        dqd = zero
        temd = zero
        tem1d = zero
        gammad = zero
        factord = zero
        IF (k .GE. kbcon(i) .AND. dq_tmp1(i,k) .GT. zero) THEN

           pwod(i,k) = pwod(i,k) + pwavod(i)

           etahd = etahd + c0*dz_tmp4(i,k)*qlk_tmp1(i,k)*pwod(i,k)
           dzd = dzd + c0*etah_tmp1(i,k)*qlk_tmp1(i,k)*pwod(i,k)
           qlkd = qlkd + c0*etah_tmp1(i,k)*dz_tmp4(i,k)*pwod(i,k)
           pwod(i,k) = zero


           qlkd = qlkd + qckod(i,k)
           qrchd = qrchd + qckod(i,k)
           qckod(i,k) = zero

           dzd = dzd - g*qlk_tmp1(i,k)*aa1d(i)
           qlkd = qlkd - g*dz_tmp4(i,k)*aa1d(i)

           IF (ncloud .GT. 0 .AND. k .GT. jmin(i)) THEN
             etahd = etahd + c1*g*dz_tmp4(i,k)*qlk_tmp1(i,k)*dellald(i,k)/dp_tmp1(i,k)
             dzd = dzd + c1*g*etah_tmp1(i,k)*qlk_tmp1(i,k)*dellald(i,k)/dp_tmp1(i,k)
             qlkd = qlkd + c1*g*etah_tmp1(i,k)*dz_tmp4(i,k)*dellald(i,k)/dp_tmp1(i,k)
             dellald(i,k) = zero

             dqd = dqd + qlkd /(eta_tmp3(i, k)+etah_tmp1(i,k)*(c0+c1)*dz_tmp4(i,k))
             eta1d(i,k) = eta1d(i,k)  &
              - dq_tmp1(i,k)*qlkd/(eta_tmp3(i, k)+etah_tmp1(i,k)*(c0+c1)*dz_tmp4(i,k))**2
             etahd = etahd - dq_tmp1(i,k)*(c0+c1)*dz_tmp4(i,k)*qlkd &
                   /(eta_tmp3(i, k)+etah_tmp1(i,k)*(c0+c1)*dz_tmp4(i,k))**2
             dzd = dzd - dq_tmp1(i,k)*(c0+c1)*etah_tmp1(i,k)*qlkd &
                  /(eta_tmp3(i, k)+etah_tmp1(i,k)*(c0+c1)*dz_tmp4(i,k))**2
             qlkd = zero
           ELSE
            dqd = dqd + qlkd/(eta_tmp3(i, k)+etah_tmp1(i,k)*c0*dz_tmp4(i,k))
            eta1d(i,k) = eta1d(i,k) - dq_tmp1(i,k)*qlkd &
                /(eta_tmp3(i, k)+etah_tmp1(i,k)*c0*dz_tmp4(i,k))**2
            etahd = etahd - dq_tmp1(i,k)*c0*dz_tmp4(i,k)*qlkd &
                /(eta_tmp3(i, k)+etah_tmp1(i,k)*c0*dz_tmp4(i,k))**2
            dzd = dzd - dq_tmp1(i,k)*c0*etah_tmp1(i,k)*qlkd &
                /(eta_tmp3(i, k)+etah_tmp1(i,k)*c0*dz_tmp4(i,k))**2
            qlkd = zero
           END IF

           eta1d(i,k) = eta1d(i,k) + half*etahd
           eta1d(i,k-1) = eta1d(i,k-1) + half*etahd
           etahd = zero
         ENDIF

           eta1d(i,k) = eta1d(i,k) + (qcko_tmp1(i, k)-qrch_tmp1(i,k))*dqd
           qckod(i,k) = qckod(i,k) + eta_tmp3(i,k)*dqd
           qrchd = qrchd - eta(i,k)*dqd
           dqd = zero

           qckod(i,k-1) = qckod(i,k-1) + (one-tem1_tmp2(i,k))*qckod(i,k)/factor_tmp2(i,k)
           tem1d = tem1d - qcko_tmp2(i,k-1)*qckod(i,k)/factor_tmp2(i,k)
           temd = temd + half*(qo_tmp4(i, k)+qo_tmp4(i, k-1))*qckod(i,k)/factor_tmp2(i,k)
           qod(i,k) = qod(i,k) + half*tem_tmp2(i,k)*qckod(i,k)/factor_tmp2(i,k)
           qod(i,k-1) = qod(i,k-1) + half*tem_tmp2(i,k)*qckod(i,k)/factor_tmp2(i,k)
           factord = factord - ((one-tem1_tmp2(i,k))*qcko_tmp2(i, k-1) &
                  +tem_tmp2(i,k)*half*(qo_tmp4(i, k)+qo_tmp4(i, k-1)))*qckod(i,k)/factor_tmp2(i,k)**2
           qckod(i,k) = zero

           temd = temd + factord
           tem1d = tem1d - factord
           factord = zero

           xlamudd(i) = xlamudd(i) + half*dz_tmp4(i,k)*tem1d
           dzd = dzd + half*xlamud(i)*tem1d
           tem1d = zero


           xlamued(i,k) = xlamued(i,k) + half*dz_tmp4(i,k)*temd
           xlamued(i,k-1) = xlamued(i,k-1) + half*dz_tmp4(i,k)*temd
           dzd = dzd + half*(xlamue_tmp4(i,k)+xlamue_tmp4(i,k-1))*temd
           temd = zero

           qesod(i,k) = qesod(i,k) + qrchd
           gammad = gammad + dbyo_tmp1(i,k)*qrchd/(hvap*(one+gamma_tmp1(i,k)))
           dbyod(i,k) = dbyod(i,k) + gamma_tmp1(i,k)*qrchd/(hvap*(one+gamma_tmp1(i,k)))
           gammad = gammad - gamma_tmp1(i,k)*dbyo_tmp1(i,k)*hvap*qrchd/(hvap*(one+gamma_tmp1(i,k)))**2
           qrchd = zero

           qesod(i,k) = qesod(i,k) + el2orc*gammad/to_tmp2(i, k)**2
           tod(i,k) = tod(i,k) -el2orc*qeso_tmp4(i,k)*2*gammad/to_tmp2(i,k)**3
           gammad = zero

           zid(i,k) = zid(i,k) + dzd
           zid(i,k-1) = zid(i,k-1) - dzd
           dzd = zero
      END IF
    END DO
  END DO


  DO i=1,im
    qod(i,kb(i)) = qod(i,kb(i)) + qckod(i,kb(i))
    qckod(i,kb(i)) = zero
    qod(i,kb(i)) = qod(i,kb(i)) + qckod(i,kb(i))
    qckod(i,kb(i)) = zero
  END DO


endif !ktcon(1) .gt. 2


  DO i=1,im
    IF (kbcon1(i) .EQ. kmax(i)) THEN
      cnvscl_2d(i) = zero
    END IF
  END DO



  do k = km1,2,-1
    do i = 1, im
       if(k.gt.kb(i).and.k.lt.kmax(i)) then
            temd = zero
            tem1d = zero
            ptemd = zero
            ptem1d = zero
            factord = zero
            dzd = zero

            hckod(i,k) = hckod(i,k) + dbyod(i,k)
            hesod(i,k) = hesod(i,k) - dbyod(i,k)
            dbyod(i,k) = zero

            vckod(i,k-1) = vckod(i,k-1) + (one-tem1_tmp1(i,k))*vckod(i,k)/factor_tmp1(i,k)
            tem1d = tem1d - vcko_tmp2(i,k-1)*vckod(i,k)/factor_tmp1(i,k)
            ptemd = ptemd + vo_tmp1(i,k)*vckod(i,k)/factor_tmp1(i,k)
            vod(i,k) = vod(i,k) + ptem_tmp3(i,k)*vckod(i,k)/factor_tmp1(i,k)
            ptem1d = ptem1d + vo_tmp1(i,k-1)*vckod(i,k)/factor_tmp1(i,k)
            vod(i,k-1) = vod(i,k-1) + ptem1_tmp1(i,k)*vckod(i,k)/factor_tmp1(i,k)
            factord = factord -  &
               ((one-tem1_tmp1(i,k))*vcko_tmp2(i, k-1)+ptem_tmp3(i,k)*vo_tmp1(i, k) &
                +ptem1_tmp1(i,k)*vo_tmp1(i, k-1))*vckod(i,k)/factor_tmp1(i,k)**2
            vckod(i,k) = zero

            uckod(i,k-1) = uckod(i,k-1) + (one-tem1_tmp1(i,k))*uckod(i,k)/factor_tmp1(i,k)
            tem1d = tem1d - ucko_tmp2(i,k-1)*uckod(i,k)/factor_tmp1(i,k)
            ptemd = ptemd + uo_tmp1(i,k)*uckod(i,k)/factor_tmp1(i,k)
            uod(i,k) = uod(i,k) + ptem_tmp3(i,k)*uckod(i,k)/factor_tmp1(i,k)
            ptem1d = ptem1d + uo_tmp1(i,k-1)*uckod(i,k)/factor_tmp1(i,k)
            uod(i,k-1) = uod(i,k-1) + ptem1_tmp1(i,k)*uckod(i,k)/factor_tmp1(i,k)
            factord = factord -  &
                ((one-tem1_tmp1(i,k))*ucko_tmp2(i,k-1)+ptem_tmp3(i,k)*uo_tmp1(i, k) &
                +ptem1_tmp1(i,k)*uo_tmp1(i,k-1))*uckod(i,k)/factor_tmp1(i,k)**2
            uckod(i,k) = zero

            hckod(i,k-1) = hckod(i,k-1) + (one-tem1_tmp1(i,k))*hckod(i,k)/factor_tmp1(i,k)
            tem1d = tem1d - hcko_tmp2(i,k-1)*hckod(i,k)/factor_tmp1(i,k)
            temd = temd + half*(heo_tmp2(i,k)+heo_tmp2(i,k-1))*hckod(i,k)/factor_tmp1(i,k)
            heod(i,k) = heod(i,k) + half*tem_tmp1(i,k)*hckod(i,k)/factor_tmp1(i,k)
            heod(i,k-1) = heod(i,k-1) + half*tem_tmp1(i,k)*hckod(i,k)/factor_tmp1(i,k)
            factord = factord  &
                  - ((one-tem1_tmp1(i,k))*hcko_tmp2(i, k-1)+tem_tmp1(i,k)*half*(heo_tmp2(i, k) &
                      +heo_tmp2(i,k-1)))*hckod(i,k)/factor_tmp1(i,k)**2
            hckod(i,k) = zero

            temd =  temd + half*ptem1d
            ptem1d = zero

            temd = temd + half*ptemd
            ptemd = zero

            temd = temd + factord
            tem1d = tem1d - factord
            factord = zero

            xlamudd(i) = xlamudd(i) + half*dz_tmp3(i,k)*tem1d
            dzd = dzd + half*xlamud(i)*tem1d
            tem1d = zero

            xlamued(i,k) = xlamued(i,k) + half*dz_tmp3(i,k)*temd
            xlamued(i,k-1) = xlamued(i,k-1) + half*dz_tmp3(i,k)*temd
            dzd = dzd + half*(xlamue_tmp4(i,k)+xlamue_tmp4(i,k-1))*temd
            temd = zero

            zid(i,k) = zid(i,k) + dzd
            zid(i,k-1) = zid(i,k-1) - dzd
            dzd = zero
        endif
     enddo
  enddo

  DO i=1,im
    indx = kb(i)

    vod(i,indx) = vod(i,indx) + vckod(i,indx)
    vckod(i,indx) = zero

    uod(i,indx) = uod(i,indx) + uckod(i,indx)
    uckod(i,indx) = zero

    heod(i,indx) = heod(i,indx) + hckod(i,indx)
    hckod(i,indx) = zero
  END DO


  DO k=km1,2,-1
    DO i=1,im
      IF (k .GT. kbcon(i) .AND. k .LT. kmax(i)) THEN
        dzd = zero
        ptemd = zero

        eta1d(i,k-1) = eta1d(i,k-1) + eta1d(i,k)*(one+ptem_tmp2(i,k)*dz_tmp2(i,k))
        ptemd = ptemd + eta_tmp3(i,k-1)*dz_tmp2(i,k)*eta1d(i,k)
        dzd = dzd + eta_tmp3(i,k-1)*ptem_tmp2(i,k)*eta1d(i,k)
        eta1d(i,k) = zero

        xlamued(i,k) = xlamued(i,k) + half*ptemd
        xlamued(i,k-1) = xlamued(i,k-1) + half*ptemd
        xlamudd(i) = xlamudd(i) - ptemd
        ptemd = zero

        zid(i,k) = zid(i,k) + dzd
        zid(i,k-1) = zid(i,k-1) - dzd
        dzd = zero
      END IF
    END DO
  END DO

  DO k=1,km1
    DO i=1,im
      IF (k .LT. kbcon(i) .AND. k .GE. kb(i)) THEN
        dzd = zero
        ptemd = zero

        eta1d(i,k+1) = eta1d(i,k+1) + eta1d(i,k)/(one+ptem_tmp1(i,k)*dz_tmp1(i,k))
        ptemd = ptemd - eta_tmp2(i,k+1)*dz_tmp1(i,k)*eta1d(i,k) &
                         /(one+ptem_tmp1(i,k)*dz_tmp1(i,k))**2
        dzd = dzd - eta_tmp2(i,k+1)*ptem_tmp1(i,k)*eta1d(i,k) &
                         /(one+ptem_tmp1(i,k)*dz_tmp1(i,k))**2
        eta1d(i,k) = zero

        xlamued(i,k) = xlamued(i,k) + half*ptemd
        xlamued(i,k+1) = xlamued(i,k+1) + half*ptemd
        xlamudd(i) = xlamudd(i) - ptemd
        ptemd = zero

        zid(i,k+1) = zid(i,k+1) + dzd
        zid(i,k) = zid(i,k) - dzd
        dzd = zero
      END IF
    END DO
  END DO

  do k = 2, km1
    do i=1,im
      if((k.ge.kbcon(i).and.k.lt.kmax(i))) then
         xlamue_tmp3d(i,k) = xlamue_tmp3d(i,k) + xlamued(i,k)
         xlamued(i,k) = zero
      endif
    enddo
  enddo


  DO k=km1,2,-1
    DO i=1,im
      IF (k .GE. kbcon(i) .AND. k .LT. kmax(i)) THEN
       tem  = cxlamu*frh(i, k)*fent2(i, k)
       temd = zero

        xlamued(i,k) = xlamued(i,k) + xlamue_tmp3d(i,k)*fent1(i,k)
        fent1d(i,k) = fent1d(i,k) + xlamue_tmp2(i,k)*xlamue_tmp3d(i,k)
        temd= temd + xlamue_tmp3d(i,k)
        xlamue_tmp3d(i,k) = zero

        frhd(i,k) = frhd(i,k) +  fent2(i, k)*temd*cxlamu
        fent2d(i,k) = fent2d(i,k) + frh(i,k)*temd*cxlamu
        temd = zero

      END IF
    END DO
  END DO

  do k = km1, 2, -1
    do i=1,im
      if((k.gt.kbcon(i).and.k.lt.kmax(i))) then
          temd=zero
          tem = qeso_tmp4(i,k)/qeso_tmp4(i,kbcon(i))
          temd = temd + 3*(tem**2)*fent2d(i,k)
          fent2d(i,k) = zero

          temd = temd + 2*tem*fent1d(i,k)
          fent1d(i,k) = zero

          qesod(i,k) = qesod(i,k) + temd/qeso_tmp4(i, kbcon(i))
          qesod(i,kbcon(i)) = qesod(i,kbcon(i)) - qeso_tmp4(i,k)*temd/qeso_tmp4(i, kbcon(i))**2
          temd = zero
      endif
    enddo
  enddo

  DO i=1,im
    xlamued(i,kbcon(i)) = xlamued(i,kbcon(i)) + xlamudd(i)
    xlamudd(i) = zero
  END DO

  DO k=km1,2,-1
    DO i=1,im
      IF (k .GT. kbcon(i) .AND. k .LT. kmax(i)) THEN
         xlamued(i,kbcon(i)) = xlamued(i,kbcon(i)) + xlamued(i,k)
         xlamued(i,k) = zero
      END IF
    END DO
  END DO

  DO i=1,im
    if(pbcdif(i).gt.cincr) then
        cnvscl_2d(i) = zero
    endif
    cnvscl_2d(i) = cnvscl_2d(i)*cnvscl_1(i)
  END DO

  DO i = 1, im
    dotd(i,kbcon(i)) = dotd(i,kbcon(i)) + pdotd(i)*10.0_r_kind
  ENDDO

  DO k=km1,1,-1
    DO i=1,im
      IF (k .LE. kmax(i) - 1) THEN
        min1d = zero
        es0d = zero
        uo1d(i,k) = uo1d(i,k) + uod(i,k)
        uod(i,k) = zero
        vo1d(i,k) = vo1d(i,k) + vod(i,k)
        vod(i,k) = zero

        vod(i,k+1) = vod(i,k+1) + half*vo1d(i,k)
        vod(i,k) = vod(i,k) + half*vo1d(i,k)
        vo1d(i,k) = zero

        uod(i,k+1) = uod(i,k+1) + half*uo1d(i,k)
        uod(i,k) = uod(i,k) + half*uo1d(i,k)
        uo1d(i,k) = zero

        zod(i,k) = zod(i,k)+ half*g*hesod(i,k)
        zod(i,k+1) = zod(i,k+1) + half*g*hesod(i,k)
        tod(i,k) = tod(i,k) + cp*hesod(i,k)
        qesod(i,k) = qesod(i,k) + hvap*hesod(i,k)
        hesod(i,k) = zero

        zod(i,k) = zod(i,k) + half*g*heod(i,k)
        zod(i,k+1) = zod(i,k+1) + half*g*heod(i,k)
        tod(i,k) = tod(i,k) + cp*heod(i,k)
        qod(i,k) = qod(i,k) + hvap*heod(i,k)
        heod(i,k) = zero

        min1d = min1d - frhd(i,k)
        frhd(i,k) = zero

        IF (qo_tmp4(i, k)/qeso_tmp4(i, k) .GT. 1.0_r_kind) THEN
          min1d = zero
        ELSE
          qod(i,k) = qod(i,k) + min1d/qeso_tmp4(i,k)
          qesod(i,k) = qesod(i,k) - qo_tmp4(i,k)*min1d/qeso_tmp4(i,k)**2
          min1d = zero
        END IF

        IF (qo_tmp3(i, k) .LT. 1.e-10_r_kind) THEN
          qod(i, k) = zero
        END IF

        IF (qeso_tmp3(i, k) .LT. 1.e-8_r_kind) THEN
          qesod(i, k) = zero
        END IF

        esd(i,k) = esd(i,k) + eps*po(i,k)*qesod(i,k)/(po(i, k)+epsm1*es_tmp3(i, k))**2
        qesod(i,k) = zero

        es0d = es0d + 10.0_r_kind * esd(i,k)
        esd(i,k) = zero

        call fpvsx_ad(to_tmp2(i,k),es0,tod(i,k),es0d,.true.)
      END IF
    END DO
  END DO



  DO k= km1,1,-1
    DO i=1,im
      IF (k .LE. kmax(i) - 1) THEN
        dqsdtd = zero
        dqsdpd = zero
        dtd = zero
        dqd = zero
        dzd = zero
        es0d = zero
        gammad = zero
        qsd = zero
        desdtd=zero
        pprimed=zero

        pod(i,k) = zero

        dqd = zero
        qod(i,k+1) = qod(i,k+1) + qod(i,k)
        dqd = dqd + qod(i,k)
        qod(i,k) = zero

        tod(i,k+1) = tod(i,k+1) + tod(i,k)
        dtd = dtd + tod(i,k)
        tod(i,k) = zero



        dqsdtd = dqsdtd + dt_tmp0(i,k)*dqd
        dtd = dtd + dqsdt_tmp0(i,k)*dqd
        dqsdpd = dqsdpd + dp_tmp0(i,k)*dqd
        dqd = zero

        dzd = dzd + g*dtd/ (cp*(one+gamma_tmp0(i,k)))
        dqsdpd = dqsdpd + hvap*dp_tmp0(i,k)*dtd/ (cp*(one+gamma_tmp0(i,k)))
        gammad = gammad - (g*dz_tmp0(i,k)+hvap*dqsdp_tmp0(i,k)*dp_tmp0(i,k))*cp*dtd/(cp*(one+gamma_tmp0(i,k)))**2
        dtd = zero

        qesod(i,k+1) = qesod(i,k+1) + el2orc*gammad/to_tmp1(i, k+1)**2
        tod(i,k+1) = tod(i,k+1) - el2orc*qeso_tmp2(i, k+1)*2*gammad/to_tmp1(i,k+1)**3
        gammad = zero

        qsd = qsd + pfld(i,k+1)*desdt_tmp0(i,k)*dqsdtd/(es_tmp2(i,k)*pprime_tmp0(i,k))
        desdtd = desdtd + pfld(i,k+1)*qs_tmp1(i,k)*dqsdtd/(es_tmp2(i,k)*pprime_tmp0(i,k))
        esd(i,k) = esd(i,k) - qs_tmp1(i,k)*pfld(i,k+1)* &
                 desdt_tmp0(i,k)*pprime_tmp0(i,k)*dqsdtd/(es_tmp2(i,k)*pprime_tmp0(i,k))**2
        pprimed = pprimed - qs_tmp1(i,k)*pfld(i,k+1)*desdt_tmp0(i,k)*es_tmp2(i,k)*dqsdtd &
                  /(es_tmp2(i,k)*pprime_tmp0(i,k))**2
        dqsdtd = zero

        esd(i,k) = esd(i,k) + (fact1/to_tmp1(i, k+1)+fact2/to_tmp1(i, k+1)**2)*desdtd
        tod(i,k+1) = tod(i,k+1) -  es_tmp2(i,k)*fact1*desdtd/to_tmp1(i, k+1)**2
        tod(i,k+1) = tod(i,k+1) - es_tmp2(i,k)*fact2*2*to_tmp1(i, k+1)*desdtd/(to_tmp1(i, k+1)**2)**2
        desdtd = zero

        qsd = qsd - dqsdpd/pprime_tmp0(i,k)
        pprimed = pprimed  + qs_tmp1(i,k)*dqsdpd/pprime_tmp0(i,k)**2
        dqsdpd = zero

        esd(i,k) = esd(i,k) + eps*qsd/pprime_tmp0(i,k)
        pprimed = pprimed - eps*es_tmp2(i,k)*qsd/pprime_tmp0(i,k)**2
        qsd = zero

        esd(i,k) = esd(i,k) + epsm1*pprimed
        pprimed = zero

        es0d = es0d + 10.0_r_kind*esd(i,k)
        esd(i,k) = zero

        call fpvsx_ad(to_tmp1(i,k+1), es0, tod(i,k+1), es0d, .true.)

        zod(i,k+1) = zod(i,k+1) + half*dzd
        zod(i,k) = zod(i,k) - half*dzd
        dzd = zero
      END IF
    END DO
  END DO




!  compute moist static energy

  DO k=1,km
    DO i=1,im
      IF (k .LE. kmax(i)) THEN
        temd = zero


        temd = temd + hesod(i,k)
        qesod(i,k) =qesod(i,k) + hvap*hesod(i,k)
        hesod(i,k) = zero

        temd = temd + heod(i,k)
        qod(i,k) = qod(i,k) + hvap*heod(i,k)
        heod(i,k) = zero

        zod(i,k) = zod(i,k) + g * temd
        tod(i,k) = tod(i,k) + cp * temd
        temd = zero
      END IF
    END DO
  END DO

  DO k=1,km1
    DO i=1,im
      zid(i,k) = zid(i,k) - (clam*xlamued(i, k)/zi(i, k)**2)
      xlamued(i,k) = zero
      zod(i,k) = zod(i,k) + half*zid(i,k)
      zod(i,k+1) = zod(i,k+1) + half*zid(i,k)
      zid(i,k) = zero
    END DO
  END DO

!  hydrostatic height assume zero terr and initially assume
!    updraft entrainment rate as an inverse function of height

  DO k=km,2,-1
    DO i=1,im
      zod(i,k-1) = zod(i,k-1) + zod(i,k)
      tvod(i,k) = tvod(i,k) - LOG(sl(i, k)/sl(i, k-1))*rd*half*zod(i,k)/g
      tvod(i,k-1) = tvod(i,k-1) - LOG(sl(i, k)/sl(i, k-1))*rd*half*zod(i,k)/g
      zod(i,k) = zero
    END DO
  END DO

  DO i=1,im
    tvod(i, 1) = tvod(i,1) -(LOG(sl(i, 1))*rd*zod(i, 1)/g)
    zod(i,1) = zero
  END DO

  DO k=1,km
    DO i=1,im
        tod(i,k) = tod(i,k) + (one+delta*qo_tmp2(i,k))*tvod(i,k)
        qod(i,k) = qod(i,k) + delta*to_tmp1(i,k)*tvod(i,k)
        tvod(i,k) = zero

        if(qo_tmp1(i,k) .lt. 1.e-10_r_kind) then
          qod(i,k) = zero
        endif

        if(qeso_tmp1(i,k) .lt. 1.e-8_r_kind) then
          qesod(i,k) = zero
        endif

        esd(i,k) = esd(i,k) + eps*pfld(i,k)*qesod(i,k)/(pfld(i, k)+epsm1*es_tmp1(i,k))**2
        qesod(i,k) = zero

        es0d = es0d + 10.0_r_kind*esd(i,k)
        esd(i,k) = zero

        call fpvsx_ad(to_tmp1(i,k), es0, tod(i,k), es0d, .true.)
    END DO
  END DO


  DO k=1,km
    DO i=1,im
       t1d(i,k) = t1d(i,k) + tod(i,k)
       tod(i,k) = zero
       q1d(i,k) = q1d(i,k) + qod(i,k)
       qod(i,k) = zero
       u1d(i,k) = u1d(i,k) + rcs(i)*uod(i,k)
       uod(i,k) = zero
       v1d(i,k) = v1d(i,k) + rcs(i)*vod(i,k)
       vod(i,k) = zero
    END DO
  END DO


  DO k=1,km
    DO i=1,im
      t1d(i,k) = t1d(i,k) + t_outd(i,k)
      t_outd(i,k) = zero
      q1d(i,k) = q1d(i,k) + q_outd(i,k)
      q_outd(i,k) = zero
      u1d(i,k) = u1d(i,k) + u_outd(i,k)
      u_outd(i,k) = zero
      v1d(i,k) = v1d(i,k) + v_outd(i,k)
      v_outd(i,k) = zero
      qld(i,k,1) = qld(i,k,1) + ql_outd(i,k,1)
      ql_outd(i,k,1) = zero 
    END DO
  END DO

 RETURN
END SUBROUTINE SASCNVN_AD


