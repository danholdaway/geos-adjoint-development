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

