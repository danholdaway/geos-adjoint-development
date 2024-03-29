module fv_sg_mod

!-----------------------------------------------------------------------
! FV sub-grid mixing
!-----------------------------------------------------------------------
  use fv_arrays_mod,  only: REAL4, REAL8, FVPRC
  use constants_mod,      only: rdgas, rvgas, cp_air, cp_vapor, hlv, hlf, kappa, grav
  use tracer_manager_mod, only: get_tracer_index
  use field_manager_mod,  only: MODEL_ATMOS
  use lin_cld_microphys_mod, only: wqs1, wqs2, wqsat2_moist
  use fv_mp_mod,          only: mp_reduce_min, is_master

implicit none
private

public  fv_dry_conv, qsmith, neg_adj3

  real(FVPRC), parameter:: esl = 0.621971831
  real(FVPRC), parameter:: tice = 273.16
  real(FVPRC), parameter:: c_ice = 2106.  ! Emanuel table, page 566
  real(FVPRC), parameter:: c_liq = 4190.
! real(FVPRC), parameter:: cv_vap = 1410.
! For consistency, cv_vap derived FMS constants:
   real(FVPRC), parameter:: cv_vap = cp_vapor - rvgas  ! 1384.5

  real(FVPRC), parameter:: dc_vap0 =  cp_vapor - c_liq   ! = -2368.
  real(FVPRC), parameter:: dc_vap1 =  cv_vap - c_liq
  real(FVPRC), parameter:: dc_ice =  c_liq - c_ice      ! =  2112.
  real(FVPRC), parameter:: hlv0 = 2.501e6   ! Emanual Appendix-2
  real(FVPRC), parameter:: hlf0 = 3.337e5   ! Emanual
  real(FVPRC), parameter:: t_ice = 273.15
  real(FVPRC), parameter:: Lv0 =  hlv0 - dc_vap0*t_ice   ! = 3.147782e6
  real(FVPRC), parameter:: Lv1 =  hlv0 - dc_vap1*t_ice
  real(FVPRC), parameter:: Li0 =  hlf0 - dc_ice*t_ice   ! = -2.431928e5 

  real(FVPRC), parameter:: zvir =  rvgas/rdgas - 1.     ! = 0.607789855
  real(FVPRC), allocatable:: table(:),des(:)

!---- version number -----
  character(len=128) :: version = '$Id: fv_sg.F90,v 1.1.2.1.2.1.30.1.90.1.20.1.2.1.4.2 2017/02/21 20:43:38 bmauer Exp $'
  character(len=128) :: tagname = '$Name: Heracles-UNSTABLE_ncepdyn_Feb222017 $'

contains

 subroutine fv_dry_conv( isd, ied, jsd, jed, is, ie, js, je, km, nq, dt,    &
                         tau, delp, pe, peln, pkz, ta, qa, ua, va,  &
                         hydrostatic, w, delz, u_dt, v_dt, t_dt, q_dt )
! Dry convective adjustment-mixing
!-------------------------------------------
      integer, intent(in):: is, ie, js, je, km, nq
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: tau         ! Relaxation time scale
      real(FVPRC), intent(in):: dt             ! model time step
      real(FVPRC), intent(in)::   pe(is-1:ie+1,km+1,js-1:je+1) 
      real(FVPRC), intent(in):: peln(is  :ie,  km+1,js  :je)
      real(FVPRC), intent(in):: delp(isd:ied,jsd:jed,km)      ! Delta p at each model level
      real(FVPRC), intent(in):: delz(isd:,jsd:,1:)      ! Delta z at each model level
      real(FVPRC), intent(in)::  pkz(is:ie,js:je,km)
      logical, intent(in)::  hydrostatic
! 
      real(FVPRC), intent(inout):: ua(isd:ied,jsd:jed,km)
      real(FVPRC), intent(inout):: va(isd:ied,jsd:jed,km)
      real(FVPRC), intent(inout)::  w(isd:,jsd:,1:)
      real(FVPRC), intent(inout):: ta(isd:ied,jsd:jed,km)      ! Temperature
      real(FVPRC), intent(inout):: qa(isd:ied,jsd:jed,km,nq)   ! Specific humidity & tracers
      real(FVPRC), intent(inout):: u_dt(isd:ied,jsd:jed,km) 
      real(FVPRC), intent(inout):: v_dt(isd:ied,jsd:jed,km) 
      real(FVPRC), intent(inout):: t_dt(is:ie,js:je,km) 
      real(FVPRC), intent(inout):: q_dt(is:ie,js:je,km,nq) 
!---------------------------Local variables-----------------------------
      real(FVPRC), dimension(is:ie,km):: u0, v0, w0, t0, hd, te, gz, tvm, pm, den
      real(FVPRC) q0(is:ie,km,nq), qcon(is:ie,km) 
      real(FVPRC), dimension(is:ie):: gzh, lcp2, icp2
      real(FVPRC) ri, pt1, pt2, ratio, tv, cv, tmp, q_liq, q_sol, cpm, cvm
      real(FVPRC) tv1, tv2, g2, h0, mc, fra, rk, rz, rdt
      real(FVPRC) dh, dhs, dq, qsw, dqsdt, tcp3
      integer i, j, k, kk, n, m, iq, km1
      real(FVPRC), parameter:: ustar2 = 1.E-4
      real(FVPRC):: cv_air
      integer :: sphum, rainwat, snowwat, graupel
#ifdef MAPL_MODE
      integer :: liq_wat_cn, liq_wat_ls, ice_wat_cn, ice_wat_ls, cld_amt_cn, cld_amt_ls
#else
      integer :: liq_wat, ice_wat, cld_amt
#endif

      cv_air = cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68
        rz = rvgas - rdgas          ! rz = zvir * rdgas
        rk = cp_air/rdgas + 1.
        cv = cp_air - rdgas

      g2 = 0.5*grav

      rdt = 1./ dt

#ifdef MAPL_MODE
        sphum    = 1
      liq_wat_cn = 2
      liq_wat_ls = 3
      ice_wat_cn = 4
      ice_wat_ls = 5
      cld_amt_cn = 6
      cld_amt_ls = 7
#else
        sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
      liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
      ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
      cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')

      if ( nq.ge.6 ) then
           rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
           snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
           graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
      endif
#endif

!------------------------------------------------------------------------
! The nonhydrostatic pressure changes if there is heating (under constant
! volume and mass is locally conserved).
!------------------------------------------------------------------------
   m = 4
   fra = dt/real(tau)

!$OMP parallel do default(none) shared(is,ie,js,je,nq,km,qa,ta,sphum,ua,va,delp,peln,     &
!$OMP                                  hydrostatic,pe,delz,g2,w,liq_wat,rainwat,ice_wat,  &
!$OMP                                  snowwat,cv_air,m,graupel,pkz,rk,rz,fra,cld_amt,    &
!$OMP                                  u_dt,rdt,v_dt,t_dt,q_dt)                           &
!$OMP                          private(kk,lcp2,icp2,tcp3,dh,dhs,dq,den,qsw,dqsdt,qcon,q0, &
!$OMP                                  t0,u0,v0,w0,h0,pm,gzh,tvm,tmp,cpm,cvm, q_liq,q_sol,&
!$OMP                                  tv,gz,hd,te,ratio,pt1,pt2,tv1,tv2,ri,mc,km1)
  do 1000 j=js,je  

    do iq=1, nq
       do k=1,km
          do i=is,ie
             q0(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo
    enddo

    do k=1,km
       do i=is,ie
          t0(i,k) = ta(i,j,k)
         tvm(i,k) = t0(i,k)*(1.+zvir*q0(i,k,sphum))
          u0(i,k) = ua(i,j,k)
          v0(i,k) = va(i,j,k)
          pm(i,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
       enddo
    enddo

    do i=is,ie
       gzh(i) = 0.
    enddo

    if( hydrostatic ) then
       do k=km, 1,-1
          do i=is,ie
                tv  = rdgas*tvm(i,k)
           den(i,k) = pm(i,k)/tv
            gz(i,k) = gzh(i) + tv*(1.-pe(i,k,j)/pm(i,k))
            hd(i,k) = cp_air*tvm(i,k)+gz(i,k)+0.5*(u0(i,k)**2+v0(i,k)**2)
             gzh(i) = gzh(i) + tv*(peln(i,k+1,j)-peln(i,k,j))
          enddo
       enddo
    else
       do k=km, 1, -1
          do i=is,ie
           den(i,k) = -delp(i,j,k)/(grav*delz(i,j,k))
             w0(i,k) = w(i,j,k)
             gz(i,k) = gzh(i)  - g2*delz(i,j,k)
                tmp  = gz(i,k) + 0.5*(u0(i,k)**2+v0(i,k)**2+w0(i,k)**2)
#ifdef MAPL_MODE
             q_liq = q0(i,k,liq_wat_cn) + q0(i,k,liq_wat_ls) 
             q_sol = q0(i,k,ice_wat_cn) + q0(i,k,ice_wat_ls) 
#else
             q_liq = q0(i,k,liq_wat) + q0(i,k,rainwat)
             q_sol = q0(i,k,ice_wat) + q0(i,k,snowwat)
#endif
          cpm = (1.-(q0(i,k,sphum)+q_liq+q_sol))*cp_air + q0(i,k,sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
          cvm = (1.-(q0(i,k,sphum)+q_liq+q_sol))*cv_air + q0(i,k,sphum)*cv_vap   + q_liq*c_liq + q_sol*c_ice
             hd(i,k) = cpm*t0(i,k) + tmp
             te(i,k) = cvm*t0(i,k) + tmp
              gzh(i) = gzh(i) - grav*delz(i,j,k)
          enddo
       enddo
    endif

   do n=1,m

      ratio = real(n)/real(m)

      do i=is,ie
         gzh(i) = 0.
      enddo

! Compute total condensate
#ifdef MAPL_MODE
      do k=1,km
         do i=is,ie
            qcon(i,k) = q0(i,k,liq_wat_cn) + q0(i,k,ice_wat_cn) + q0(i,k,liq_wat_ls) + q0(i,k,ice_wat_ls)
         enddo
      enddo
#else
   if ( nq .le. 3 .or. zvir .lt. 1.e-3 ) then
      do k=1,km
         do i=is,ie
            qcon(i,k) = 0.
         enddo
      enddo
   elseif ( nq .le. 5 ) then
      do k=1,km
         do i=is,ie
            qcon(i,k) = q0(i,k,liq_wat) + q0(i,k,ice_wat)
         enddo
      enddo
   else
      do k=1,km
         do i=is,ie
            qcon(i,k) = q0(i,k,liq_wat)+q0(i,k,ice_wat)+q0(i,k,snowwat)+q0(i,k,rainwat)+q0(i,k,graupel)
         enddo
      enddo
   endif
#endif

      do k=km, 2, -1
         km1 = k-1
         do i=is,ie
! Richardson number = g*delz * del_theta/theta / (del_u**2 + del_v**2)
! Use exact form for "density temperature"
            tv1 = t0(i,km1)*(1.+zvir*q0(i,km1,sphum)-qcon(i,km1))
            tv2 = t0(i,k  )*(1.+zvir*q0(i,k  ,sphum)-qcon(i,k))
            pt1 = tv1 / pkz(i,j,km1)
            pt2 = tv2 / pkz(i,j,k  )
            ri = (gz(i,k-1)-gz(i,k))*(pt1-pt2)/( 0.5*(pt1+pt2)*        &
                ((u0(i,k-1)-u0(i,k))**2+(v0(i,k-1)-v0(i,k))**2+ustar2) )

! Dry convective adjustment for K-H instability:
! Compute equivalent mass flux: mc
            mc = 0.
            if ( ri < 1. ) then                ! Dry case:
               mc = ratio*delp(i,j,km1)*delp(i,j,k)/(delp(i,j,km1)+delp(i,j,k))*(1.-max(0.0,ri))**2
#ifdef DO_MOIST_ADJ
            elseif ( pm(i,k) > 75.e2 ) then   ! Moist case:
              tcp3 = hlv + hlf * min(1., dim(tice,t0(i,km1))/44.)**2
               qsw = wqs1(t0(i,km1), den(i,km1))
               dh  = cp_air*(tv2-tv1) + (gz(i,k)-gz(i,km1))
               dhs = dh + tcp3*(q0(i,k,sphum) - qsw)
! For ri>1 generally implies dh<0
! dhs > 0 implies q0(k)>qsw (if dh<0)
! Therefore, q0(k) > qsw > q0(k-1)
               if ( dhs>0.0 .and. q0(i,km1,sphum)<qsw ) then
                    mc = delp(i,j,k)*min( ratio*dhs/(tcp3*q0(i,k,sphum)),  &
                                          0.5*delp(i,j,km1)/(delp(i,j,km1)+delp(i,j,k)) )
               endif
            endif
            if ( mc > 1.e-5 ) then
#endif
                 do iq=1,nq
                    h0 = mc*(q0(i,k,iq)-q0(i,km1,iq))
                    q0(i,km1,iq) = q0(i,km1,iq) + h0/delp(i,j,km1)
                    q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                 enddo
! Recompute qcon
#ifdef MAPL_MODE
                    qcon(i,km1) = q0(i,km1,liq_wat_cn) + q0(i,km1,ice_wat_cn) + q0(i,km1,liq_wat_ls) + q0(i,km1,ice_wat_ls)
#else
                 if ( nq .le. 3 .or. zvir .lt. 1.e-3 ) then
                    qcon(i,km1) = 0.
                 elseif ( nq .le. 5 ) then
                    qcon(i,km1) = q0(i,km1,liq_wat) + q0(i,km1,ice_wat)
                 else
                    qcon(i,km1) = q0(i,km1,liq_wat) + q0(i,km1,ice_wat) +                  &
                                  q0(i,km1,snowwat) + q0(i,km1,rainwat) + q0(i,km1,graupel)
                 endif
#endif
! u:
                 h0 = mc*(u0(i,k)-u0(i,k-1))
                 u0(i,k-1) = u0(i,k-1) + h0/delp(i,j,k-1)
                 u0(i,k  ) = u0(i,k  ) - h0/delp(i,j,k  )
! v:
                 h0 = mc*(v0(i,k)-v0(i,k-1))
                 v0(i,k-1) = v0(i,k-1) + h0/delp(i,j,k-1)
                 v0(i,k  ) = v0(i,k  ) - h0/delp(i,j,k  )

              if ( hydrostatic ) then
! Static energy
                        h0 = mc*(hd(i,k)-hd(i,k-1))
                 hd(i,k-1) = hd(i,k-1) + h0/delp(i,j,k-1)
                 hd(i,k  ) = hd(i,k  ) - h0/delp(i,j,k  )
              else
! Total energy
                        h0 = mc*(hd(i,k)-hd(i,k-1))
                 te(i,k-1) = te(i,k-1) + h0/delp(i,j,k-1)
                 te(i,k  ) = te(i,k  ) - h0/delp(i,j,k  )
! w:
                        h0 = mc*(w0(i,k)-w0(i,k-1))
                 w0(i,k-1) = w0(i,k-1) + h0/delp(i,j,k-1)
                 w0(i,k  ) = w0(i,k  ) - h0/delp(i,j,k  )
              endif
            endif
         enddo

!-------------- 
! Retrive Temp:
!--------------
       if ( hydrostatic ) then
         kk = k
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ( rk - pe(i,kk,j)/pm(i,kk) )
              gzh(i) = gzh(i) + t0(i,kk)*(peln(i,kk+1,j)-peln(i,kk,j))
            t0(i,kk) = t0(i,kk) / ( rdgas + rz*q0(i,kk,sphum) )
         enddo
         kk = k-1
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ((rk-pe(i,kk,j)/pm(i,kk))*(rdgas+rz*q0(i,kk,sphum)))
         enddo
       else
! Non-hydrostatic under constant volume heating/cooling
         do kk=k-1,k
            do i=is,ie
               tv = gz(i,kk) + 0.5*(u0(i,kk)**2+v0(i,kk)**2+w0(i,kk)**2)
#ifdef MAPL_MODE
               q_liq = q0(i,kk,liq_wat_cn) + q0(i,kk,liq_wat_ls)
               q_sol = q0(i,kk,ice_wat_cn) + q0(i,kk,ice_wat_ls)
#else
               q_liq = q0(i,kk,liq_wat) + q0(i,kk,rainwat)
               q_sol = q0(i,kk,ice_wat) + q0(i,kk,snowwat)
#endif
         cpm = (1.-(q0(i,kk,sphum)+q_liq+q_sol))*cp_air + q0(i,kk,sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
         cvm = (1.-(q0(i,kk,sphum)+q_liq+q_sol))*cv_air + q0(i,kk,sphum)*cv_vap   + q_liq*c_liq + q_sol*c_ice
               t0(i,kk) = (te(i,kk)- tv) / cvm
               hd(i,kk) = cpm*t0(i,kk) + tv
            enddo
         enddo
       endif
      enddo   ! k-loop
   enddo       ! n-loop

!--------------------
   if ( fra < 1. ) then
      do k=1, km
         do i=is,ie
            t0(i,k) = ta(i,j,k) + (t0(i,k) - ta(i,j,k))*fra
            u0(i,k) = ua(i,j,k) + (u0(i,k) - ua(i,j,k))*fra
            v0(i,k) = va(i,j,k) + (v0(i,k) - va(i,j,k))*fra
         enddo
      enddo

      if ( .not. hydrostatic ) then
         do k=1,km
            do i=is,ie
               w0(i,k) = w(i,j,k) + (w0(i,k) - w(i,j,k))*fra
            enddo
         enddo
      endif

      do iq=1,nq
         do k=1,km
            do i=is,ie
               q0(i,k,iq) = qa(i,j,k,iq) + (q0(i,k,iq) - qa(i,j,k,iq))*fra
            enddo
         enddo
      enddo
   endif

!----------------------
! Saturation adjustment
!----------------------
    do k=1, km
      if ( hydrostatic ) then
        do i=is, ie
! Compute pressure hydrostatically
           den(i,k) = pm(i,k)/(rdgas*t0(i,k)*(1.+zvir*q0(i,k,sphum)))
#ifdef MAPL_MODE
           q_liq = q0(i,k,liq_wat_cn) + q0(i,k,liq_wat_ls)
           q_sol = q0(i,k,ice_wat_cn) + q0(i,k,ice_wat_ls)
#else
           q_liq = q0(i,k,liq_wat) + q0(i,k,rainwat)
           q_sol = q0(i,k,ice_wat) + q0(i,k,snowwat)
#endif
           cpm = (1.-(q0(i,k,sphum)+q_liq+q_sol))*cp_air + q0(i,k,sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
           lcp2(i) = hlv / cpm
           icp2(i) = hlf / cpm
        enddo
      else
        do i=is, ie
           den(i,k) = -delp(i,j,k)/(grav*delz(i,j,k))
#ifdef MAPL_MODE
           q_liq = q0(i,k,liq_wat_cn) + q0(i,k,liq_wat_ls) 
           q_sol = q0(i,k,ice_wat_cn) + q0(i,k,ice_wat_ls)
#else
           q_liq = q0(i,k,liq_wat) + q0(i,k,rainwat)
           q_sol = q0(i,k,ice_wat) + q0(i,k,snowwat)
#endif
           cvm = (1.-(q0(i,k,sphum)+q_liq+q_sol))*cv_air + q0(i,k,sphum)*cv_vap + q_liq*c_liq + q_sol*c_ice
           lcp2(i) = (Lv1+dc_vap1*t0(i,k)) / cvm
           icp2(i) = (Li0+dc_ice*t0(i,k)) / cvm
        enddo
      endif

! Prevent super saturation over water:
       do i=is, ie
          qsw = wqs2(t0(i,k), den(i,k), dqsdt)
           dq = q0(i,k,sphum) - qsw
          if ( dq > 0. ) then   ! remove super-saturation
             tcp3 = lcp2(i) + icp2(i)*min(1., dim(tice,t0(i,k))/40.)
              tmp = dq/(1.+tcp3*dqsdt)
             t0(i,k) = t0(i,k) + tmp*lcp2(i)
             q0(i,k,  sphum) = q0(i,k,  sphum) - tmp
#ifndef MAPL_MODE
             q0(i,k,liq_wat) = q0(i,k,liq_wat) + tmp
! Grid box mean is saturated; 50% or higher cloud cover
             qa(i,j,k,cld_amt) = max(0.5, min(1., qa(i,j,k,cld_amt)+25.*tmp/qsw))
#else
             q0(i,k,liq_wat_ls) = q0(i,k,liq_wat_ls) + tmp
! Grid box mean is saturated; 50% or higher cloud cover
             qa(i,j,k,cld_amt_ls) = max(0.5, min(1., qa(i,j,k,cld_amt_ls)+25.*tmp/qsw))
#endif
          endif
! Freezing
          tmp = tice-40. - t0(i,k)
#ifndef MAPL_MODE
          if( tmp>0.0 .and. q0(i,k,liq_wat)>0. ) then
              dh = min( q0(i,k,liq_wat), q0(i,k,liq_wat)*tmp*0.125, tmp/icp2(i) )
              q0(i,k,liq_wat) = q0(i,k,liq_wat) - dh
              q0(i,k,ice_wat) = q0(i,k,ice_wat) + dh
              t0(i,k) =  t0(i,k) + dh*icp2(i)
          endif
#else
          if( tmp>0.0 ) then
            if( q0(i,k,liq_wat_cn)>0. ) then
              dh = min( q0(i,k,liq_wat_cn), q0(i,k,liq_wat_cn)*tmp*0.125, tmp/icp2(i) )
              q0(i,k,liq_wat_cn) = q0(i,k,liq_wat_cn) - dh
              q0(i,k,ice_wat_cn) = q0(i,k,ice_wat_cn) + dh
              t0(i,k) =  t0(i,k) + dh*icp2(i)
            endif
            if( q0(i,k,liq_wat_ls)>0. ) then
              dh = min( q0(i,k,liq_wat_ls), q0(i,k,liq_wat_ls)*tmp*0.125, tmp/icp2(i) )
              q0(i,k,liq_wat_ls) = q0(i,k,liq_wat_ls) - dh
              q0(i,k,ice_wat_ls) = q0(i,k,ice_wat_ls) + dh
              t0(i,k) =  t0(i,k) + dh*icp2(i)
            endif
          endif
#endif
       enddo
    enddo

   do k=1,km
      do i=is,ie
         u_dt(i,j,k) = rdt*(u0(i,k) - ua(i,j,k))
         v_dt(i,j,k) = rdt*(v0(i,k) - va(i,j,k))
         t_dt(i,j,k) = 0.
           ta(i,j,k) = t0(i,k)   ! *** temperature updated ***
      enddo
      do iq=1,nq
#ifndef MAPL_MODE
         if (iq .eq. cld_amt ) then
#else
         if ( (iq .eq. cld_amt_cn) .or. (iq .eq. cld_amt_ls) ) then
#endif
            do i=is,ie
               q_dt(i,j,k,iq) = 0.
            enddo
         else
            do i=is,ie
               q_dt(i,j,k,iq) = rdt*(q0(i,k,iq)-qa(i,j,k,iq))
            enddo
         endif 
      enddo
   enddo

   if ( .not. hydrostatic ) then
      do k=1,km
         do i=is,ie
            w(i,j,k) = w0(i,k)   ! w updated
         enddo
      enddo
   endif

1000 continue


 end subroutine fv_dry_conv



 real(FVPRC) function qs1d(t, p, q)
! Based on "moist" mixing ratio, p is the total (dry+vapor) pressure
  real(FVPRC), intent(in):: t, p, q
! Local:
  real(FVPRC) es, ap1
  real(FVPRC), parameter:: Tmin=tice - 160.
  integer it

       ap1 = 10.*DIM(t, Tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table(it) + (ap1-it)*des(it)
      qs1d = esl*es*(1.+zvir*q)/p

  end function qs1d


  subroutine qsmith_init
  integer, parameter:: length=2621 
  integer i

  if( .not. allocated(table) ) then
!                            Generate es table (dT = 0.1 deg. C)

       allocate ( table(length) )
       allocate (  des (length) )

       call qs_table(length, table)

       do i=1,length-1
          des(i) = table(i+1) - table(i)
       enddo
       des(length) = des(length-1)
  endif
 
  end subroutine qsmith_init


  subroutine qsmith(im, km, k1, t, p, q, qs, dqdt)
! input T in deg K; p (Pa)
  integer, intent(in):: im, km, k1
  real(FVPRC), intent(in),dimension(im,km):: t, p, q
  real(FVPRC), intent(out),dimension(im,km):: qs
  real(FVPRC), intent(out), optional:: dqdt(im,km)
! Local:
  real(FVPRC) es(im,km)
  real(FVPRC) ap1, eps10
  real(FVPRC) Tmin
  integer i, k, it

  Tmin = tice-160.
  eps10  = 10.*esl

  if( .not. allocated(table) ) call  qsmith_init
 
      do k=k1,km
         do i=1,im
            ap1 = 10.*DIM(t(i,k), Tmin) + 1.
            ap1 = min(2621., ap1)
            it = ap1
            es(i,k) = table(it) + (ap1-it)*des(it)
            qs(i,k) = esl*es(i,k)*(1.+zvir*q(i,k))/p(i,k)
         enddo
      enddo

      if ( present(dqdt) ) then
      do k=k1,km
           do i=1,im
              ap1 = 10.*DIM(t(i,k), Tmin) + 1.
              ap1 = min(2621., ap1) - 0.5
              it  = ap1
              dqdt(i,k) = eps10*(des(it)+(ap1-it)*(des(it+1)-des(it)))*(1.+zvir*q(i,k))/p(i,k)
           enddo
      enddo
      endif
 
  end subroutine qsmith
 

 subroutine qs_table(n,table)
      integer, intent(in):: n
      real(FVPRC) table (n)
      real(FVPRC):: dt=0.1
      real(FVPRC) esbasw, tbasw, esbasi, tbasi, Tmin, tem, aa, b, c, d, e, esh20 
      real(FVPRC) wice, wh2o
      integer i
! Constants
      esbasw = 1013246.0
      tbasw =   373.16
      tbasi =   273.16
      Tmin = tbasi - 160.
!  Compute es over water
!  see smithsonian meteorological tables page 350.
      do  i=1,n
          tem = Tmin+dt*real(i-1)
          aa  = -7.90298*(tbasw/tem-1)
          b   =  5.02808*log10(tbasw/tem)
          c   = -1.3816e-07*(10**((1-tem/tbasw)*11.344)-1)
          d   =  8.1328e-03*(10**((tbasw/tem-1)*(-3.49149))-1)
          e   = log10(esbasw)
          table(i)  = 0.1*10**(aa+b+c+d+e)
      enddo

 end subroutine qs_table

 subroutine qs_table_m(n,table)
! Mixed (blended) table
      integer, intent(in):: n
      real(FVPRC) table (n)
      real(FVPRC) esupc(200)
      real(FVPRC):: dt=0.1
      real(FVPRC) esbasw, tbasw, esbasi, tbasi, Tmin, tem, aa, b, c, d, e, esh20 
      real(FVPRC) wice, wh2o
      integer i

! Constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16
! ****************************************************
!  Compute es over ice between -160c and 0 c.
      Tmin = tbasi - 160.
!  see smithsonian meteorological tables page 350.
      do i=1,1600
         tem = Tmin+dt*real(i-1)
         aa  = -9.09718 *(tbasi/tem-1.0)
         b   = -3.56654 *log10(tbasi/tem)
         c   =  0.876793*(1.0-tem/tbasi)
         e   = log10(esbasi)
         table(i)=10**(aa+b+c+e)
      enddo
! *****************************************************
!  Compute es over water between -20c and 102c.
!  see smithsonian meteorological tables page 350.
      do  i=1,1221
          tem = 253.16+dt*real(i-1)
          aa  = -7.90298*(tbasw/tem-1)
          b   =  5.02808*log10(tbasw/tem)
          c   = -1.3816e-07*(10**((1-tem/tbasw)*11.344)-1)
          d   =  8.1328e-03*(10**((tbasw/tem-1)*(-3.49149))-1)
          e   = log10(esbasw)
          esh20  = 10**(aa+b+c+d+e)
          if (i <= 200) then
              esupc(i) = esh20
          else
              table(i+1400) = esh20
          endif
      enddo
!********************************************************************
!  Derive blended es over ice and supercooled water between -20c and 0c
      do i=1,200
         tem  = 253.16+dt*real(i-1)
         wice = 0.05*(273.16-tem)
         wh2o = 0.05*(tem-253.16)
         table(i+1400) = wice*table(i+1400)+wh2o*esupc(i)
      enddo

      do i=1,n
         table(i) = table(i)*0.1
      enddo

 end subroutine qs_table_m

 subroutine neg_adj3(is, ie, js, je, ng, kbot, hydrostatic,   &
                     peln, delz, pt, dp, qv, ql, qr, qi, qs, qg, qa, check_negative)

! This is designed for 6-class micro-physics schemes
 integer, intent(in):: is, ie, js, je, ng, kbot
 logical, intent(in):: hydrostatic
 real(FVPRC), intent(in):: dp(is-ng:ie+ng,js-ng:je+ng,kbot)  ! total delp-p
 real(FVPRC), intent(in):: delz(is-ng:,js-ng:,1:)
 real(FVPRC), intent(in):: peln(is:ie,kbot+1,js:je)           ! ln(pe)
 logical, intent(in):: check_negative
 real(FVPRC), intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,kbot)::    &
                                 pt, qv, ql, qr, qi, qs, qg, qa
! Local:
 logical:: sat_adj = .false.
 real(FVPRC), parameter :: t48 = tice - 48.
 real(FVPRC), dimension(is:ie,kbot):: dpk, q2
real(FVPRC), dimension(is:ie,js:je):: pt2, qv2, ql2, qi2, qs2, qr2, qg2, dp2, p2, icpk, lcpk
 real(FVPRC):: cv_air
 real(FVPRC):: dq, qsum, dq1, q_liq, q_sol, cpm, sink, qsw, dwsdt
 integer i, j, k

 cv_air = cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68

  if ( check_negative ) then
     call prt_negative('Temperature', pt, is, ie, js, je, ng, kbot, 165.0_FVPRC)
     call prt_negative('sphum',   qv, is, ie, js, je, ng, kbot, -1.e-8_FVPRC)
     call prt_negative('liq_wat', ql, is, ie, js, je, ng, kbot, -1.e-7_FVPRC)
     call prt_negative('rainwat', qr, is, ie, js, je, ng, kbot, -1.e-7_FVPRC)
     call prt_negative('ice_wat', qi, is, ie, js, je, ng, kbot, -1.e-7_FVPRC)
     call prt_negative('snowwat', qs, is, ie, js, je, ng, kbot, -1.e-7_FVPRC)
     call prt_negative('graupel', qg, is, ie, js, je, ng, kbot, -1.e-7_FVPRC)
  endif

!$OMP parallel do default(none) shared(is,ie,js,je,kbot,qv,ql,qi,qs,qr,qg,dp,pt,       &
!$OMP                                  hydrostatic,peln,delz,cv_air,sat_adj) &
!$OMP                          private(dq,dq1,qsum,dp2,p2,pt2,qv2,ql2,qi2,qs2,qg2,qr2, &
!$OMP                                  lcpk,icpk,qsw,dwsdt,sink,q_liq,q_sol,cpm)
  do k=1, kbot
     do j=js, je
        do i=is, ie
        qv2(i,j) = qv(i,j,k)
        ql2(i,j) = ql(i,j,k)
        qi2(i,j) = qi(i,j,k)
        qs2(i,j) = qs(i,j,k)
        qr2(i,j) = qr(i,j,k)
        qg2(i,j) = qg(i,j,k)
        dp2(i,j) = dp(i,j,k)
        pt2(i,j) = pt(i,j,k)
        enddo
     enddo

     if ( hydrostatic ) then
       do j=js, je
          do i=is, ie
             p2(i,j) = dp2(i,j)/(peln(i,k+1,j)-peln(i,k,j))
             q_liq = max(0., ql2(i,j) + qr2(i,j))
             q_sol = max(0., qi2(i,j) + qs2(i,j))
             cpm = (1.-(qv2(i,j)+q_liq+q_sol))*cp_air + qv2(i,j)*cp_vapor + q_liq*c_liq + q_sol*c_ice
             lcpk(i,j) = hlv / cpm
             icpk(i,j) = hlf / cpm
          enddo
       enddo
     else
       do j=js, je
          do i=is, ie
             p2(i,j) = -dp2(i,j)/(grav*delz(i,j,k))*rdgas*pt2(i,j)*(1.+zvir*qv2(i,j))
             q_liq = max(0., ql2(i,j) + qr2(i,j))
             q_sol = max(0., qi2(i,j) + qs2(i,j))
             cpm = (1.-(qv2(i,j)+q_liq+q_sol))*cv_air + qv2(i,j)*cv_vap + q_liq*c_liq + q_sol*c_ice
             lcpk(i,j) = (Lv1+dc_vap1*pt2(i,j)) / cpm
             icpk(i,j) = (Li0+dc_ice*pt2(i,j)) / cpm
          enddo
       enddo
     endif

! Fix the negatives:
!-----------
! Ice-phase:
!-----------
    do j=js, je
       do i=is, ie
        qsum = qi2(i,j) + qs2(i,j)
        if ( qsum > 0. ) then
             if ( qi2(i,j) < 0. ) then
                  qi2(i,j) = 0.
                  qs2(i,j) = qsum
             elseif ( qs2(i,j) < 0. ) then
                  qs2(i,j) = 0.
                  qi2(i,j) = qsum
             endif
        else
! borrow froom graupel
             qi2(i,j) = 0.
             qs2(i,j) = 0.
             qg2(i,j) = qg2(i,j) + qsum
        endif

! At this stage qi and qs should be positive definite
! If graupel < 0 then borrow from qs then qi
        if ( qg2(i,j) < 0. ) then
             dq = min( qs2(i,j), -qg2(i,j) )
             qs2(i,j) = qs2(i,j) - dq
             qg2(i,j) = qg2(i,j) + dq
             if ( qg2(i,j) < 0. ) then
! if qg is still negative
                  dq = min( qi2(i,j), -qg2(i,j) )
                  qi2(i,j) = qi2(i,j) - dq
                  qg2(i,j) = qg2(i,j) + dq
             endif
        endif

! If qg is still negative then borrow from rain water: phase change
        if ( qg2(i,j)<0. .and. qr2(i,j)>0. ) then
             dq = min( qr2(i,j), -qg2(i,j) )
             qg2(i,j) = qg2(i,j) + dq
             qr2(i,j) = qr2(i,j) - dq
             pt2(i,j) = pt2(i,j) + dq*icpk(i,j)  ! conserve total energy
        endif
! If qg is still negative then borrow from cloud water: phase change
        if ( qg2(i,j)<0. .and. ql2(i,j)>0. ) then
             dq = min( ql2(i,j), -qg2(i,j) )
             qg2(i,j) = qg2(i,j) + dq
             ql2(i,j) = ql2(i,j) - dq
             pt2(i,j) = pt2(i,j) + dq*icpk(i,j)
        endif
! Last resort; borrow from water vapor
        if ( qg2(i,j)<0. .and. qv2(i,j)>0. ) then
             dq = min( 0.999*qv2(i,j), -qg2(i,j) )
             qg2(i,j) = qg2(i,j) + dq
             qv2(i,j) = qv2(i,j) - dq
             pt2(i,j) = pt2(i,j) + dq*(icpk(i,j)+lcpk(i,j))
        endif

!--------------
! Liquid phase:
!--------------
        qsum = ql2(i,j) + qr2(i,j)
        if ( qsum > 0. ) then
             if ( qr2(i,j) < 0. ) then
                  qr2(i,j) = 0.
                  ql2(i,j) = qsum
             elseif ( ql2(i,j) < 0. ) then
                  ql2(i,j) = 0.
                  qr2(i,j) = qsum
             endif
        else
          ql2(i,j) = 0.
          qr2(i,j) = qsum     ! rain water is still negative
! fill negative rain with qg first
          dq = min( max(0.0, qg2(i,j)), -qr2(i,j) )
          qr2(i,j) = qr2(i,j) + dq
          qg2(i,j) = qg2(i,j) - dq
          pt2(i,j) = pt2(i,j) - dq*icpk(i,j)
          if ( qr(i,j,k) < 0. ) then
! fill negative rain with available qi & qs (cooling)
               dq = min( qi2(i,j)+qs2(i,j), -qr2(i,j) )
               qr2(i,j) = qr2(i,j) + dq
               dq1 = min( dq, qs2(i,j) )
               qs2(i,j) = qs2(i,j) - dq1
               qi2(i,j) = qi2(i,j) + dq1 - dq 
               pt2(i,j) = pt2(i,j) - dq*icpk(i,j)
          endif
! fix negative rain water with available vapor
          if ( qr2(i,j)<0. .and. qv2(i,j)>0. ) then
               dq = min( 0.999*qv2(i,j), -qr2(i,j) )
               qv2(i,j) = qv2(i,j) - dq
               qr2(i,j) = qr2(i,j) + dq
               pt2(i,j) = pt2(i,j) + dq*lcpk(i,j)
          endif
        endif
     enddo
   enddo

!******************************************
! Fast moist physics: Saturation adjustment
!******************************************
 if ( sat_adj ) then

   do j=js, je
     do i=is, ie
! Melting of cloud ice into cloud water ********
        if ( qi2(i,j)>1.e-8 .and. pt2(i,j) > tice ) then
           sink = min( qi2(i,j), (pt2(i,j)-tice)/icpk(i,j) )
           ql2(i,j) = ql2(i,j) + sink
           qi2(i,j) = qi2(i,j) - sink
           pt2(i,j) = pt2(i,j) - sink*icpk(i,j)
        endif

! vapor <---> liquid water --------------------------------
        qsw = wqsat2_moist(pt2(i,j), qv2(i,j), p2(i,j), dwsdt)
        sink = min( ql2(i,j), (qsw-qv2(i,j))/(1.+lcpk(i,j)*dwsdt) )
        qv2(i,j) = qv2(i,j) + sink
        ql2(i,j) = ql2(i,j) - sink
        pt2(i,j) = pt2(i,j) - sink*lcpk(i,j)
!-----------------------------------------------------------

! freezing of cloud water ********
        if( ql2(i,j)>1.e-8 .and. pt2(i,j) < t48 ) then
! Enforce complete freezing below t_00 (-48 C)
            sink = min( ql2(i,j), (t48-pt2(i,j))/icpk(i,j) )
            ql2(i,j) = ql2(i,j) - sink
            qi2(i,j) = qi2(i,j) + sink
            pt2(i,j) = pt2(i,j) + sink*icpk(i,j)
        endif ! significant ql existed
     enddo
   enddo
 endif

!----------------------------------------------------------------
! Update fields:
   do j=js, je
     do i=is, ie
        qv(i,j,k) = qv2(i,j)
        ql(i,j,k) = ql2(i,j)
        qi(i,j,k) = qi2(i,j)
        qs(i,j,k) = qs2(i,j)
        qr(i,j,k) = qr2(i,j)
        qg(i,j,k) = qg2(i,j)
        pt(i,j,k) = pt2(i,j)
     enddo
   enddo

 enddo

!$OMP parallel do default(none) shared(is,ie,js,je,kbot,dp,qg,qr) &
!$OMP                          private(dpk, q2)
 do j=js, je
! Graupel:
    do k=1,kbot
       do i=is,ie
          dpk(i,k) = dp(i,j,k)
           q2(i,k) = qg(i,j,k)
       enddo
    enddo
    call fillq(ie-is+1, kbot, q2, dpk)
    do k=1,kbot
       do i=is,ie
          qg(i,j,k) = q2(i,k)
       enddo
    enddo
! Rain water:
    do k=1,kbot
       do i=is,ie
          q2(i,k) = qr(i,j,k)
       enddo
    enddo
    call fillq(ie-is+1, kbot, q2, dpk)
    do k=1,kbot
       do i=is,ie
          qr(i,j,k) = q2(i,k)
       enddo
    enddo
 enddo

!-----------------------------------
! Fix water vapor
!-----------------------------------
! Top layer: borrow from below
    k = 1
!$OMP parallel do default(none) shared(is,ie,js,je,k,qv,dp)
   do j=js, je
       do i=is, ie
          if( qv(i,j,k) < 0. ) then
              qv(i,j,k+1) = qv(i,j,k+1) + qv(i,j,k)*dp(i,j,k)/dp(i,j,k+1)
              qv(i,j,k  ) = 0.
          endif
     enddo
   enddo

! this OpenMP do-loop cannot be parallelized with recursion on k/k-1
!$OMP parallel do default(none) shared(is,ie,js,je,kbot,qv,dp) &
!$OMP                          private(dq)
 do j=js, je
   do k=2,kbot-1
       do i=is, ie
          if( qv(i,j,k) < 0. .and. qv(i,j,k-1) > 0. ) then
              dq = min(-qv(i,j,k)*dp(i,j,k), qv(i,j,k-1)*dp(i,j,k-1))
              qv(i,j,k-1) = qv(i,j,k-1) - dq/dp(i,j,k-1) 
              qv(i,j,k  ) = qv(i,j,k  ) + dq/dp(i,j,k  ) 
          endif
          if( qv(i,j,k) < 0. ) then
              qv(i,j,k+1) = qv(i,j,k+1) + qv(i,j,k)*dp(i,j,k)/dp(i,j,k+1)
              qv(i,j,k  ) = 0.
          endif
       enddo
    enddo
  enddo
 
! Bottom layer; Borrow from above
!$OMP parallel do default(none) shared(is,ie,js,je,kbot,qv,dp) private(dq)
  do j=js, je
     do i=is, ie
     if( qv(i,j,kbot) < 0. ) then
         do k=kbot-1,1,-1
            if ( qv(i,j,kbot)>=0. ) goto 123
            if ( qv(i,j,k) > 0. ) then
                 dq = min(-qv(i,j,kbot)*dp(i,j,kbot), qv(i,j,k)*dp(i,j,k))
                 qv(i,j,k   ) = qv(i,j,k   ) - dq/dp(i,j,k) 
                 qv(i,j,kbot) = qv(i,j,kbot) + dq/dp(i,j,kbot) 
            endif
         enddo   ! k-loop
123      continue
     endif
     enddo ! i-loop
 enddo   ! j-loop

!-----------------------------------
! Fix negative cloud fraction
!-----------------------------------
! this OpenMP do-loop cannot be parallelized by the recursion on k/k+1
!$OMP parallel do default(none) shared(is,ie,js,je,kbot,qa,dp)
 do j=js, je
    do k=1,kbot-1
       do i=is, ie
          if( qa(i,j,k) < 0. ) then
              qa(i,j,k+1) = qa(i,j,k+1) + qa(i,j,k)*dp(i,j,k)/dp(i,j,k+1)
              qa(i,j,k  ) = 0.
          endif
     enddo
   enddo
 enddo
 
! Bottom layer; Borrow from above
!$OMP parallel do default(none) shared(is,ie,js,je,qa,kbot,dp) &
!$OMP                          private(dq)
  do j=js, je
     do i=is, ie
        if( qa(i,j,kbot) < 0. .and. qa(i,j,kbot-1)>0.) then
            dq = min(-qa(i,j,kbot)*dp(i,j,kbot), qa(i,j,kbot-1)*dp(i,j,kbot-1))
            qa(i,j,kbot-1) = qa(i,j,kbot-1) - dq/dp(i,j,kbot-1) 
            qa(i,j,kbot  ) = qa(i,j,kbot  ) + dq/dp(i,j,kbot  ) 
        endif
! if qa is still < 0
        qa(i,j,kbot) = max(0., qa(i,j,kbot))
   enddo
 enddo


 end subroutine neg_adj3

 subroutine fillq(im, km, q, dp)
! Aggresive 1D filling algorithm for qr and qg
 integer, intent(in):: im, km
 real(FVPRC), intent(inout), dimension(im,km):: q, dp
 integer:: i, k
 real(FVPRC):: sum1, sum2, dq

 do 500 i=1,im
    sum1 = 0.
    do k=1,km
       if ( q(i,k)>0. ) then
            sum1 = sum1 + q(i,k)*dp(i,k)
       endif
    enddo
    if ( sum1<1.E-12  ) goto 500
    sum2 = 0.
    do k=km,1,-1
       if ( q(i,k)<0.0 .and. sum1>0. ) then
            dq = min( sum1, -q(i,k)*dp(i,k) )
            sum1 = sum1 - dq
            sum2 = sum2 + dq
            q(i,k) = q(i,k) + dq/dp(i,k)
       endif
    enddo
    do k=km,1,-1
       if ( q(i,k)>0.0 .and. sum2>0. ) then
            dq = min( sum2, q(i,k)*dp(i,k) )
            sum2 = sum2 - dq
            q(i,k) = q(i,k) - dq/dp(i,k)
       endif
    enddo
500  continue

 end subroutine fillq

 subroutine prt_negative(qname, q, is, ie, js, je, n_g, km, threshold)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je
      integer, intent(in):: n_g, km
      real(FVPRC), intent(in)::    q(is-n_g:ie+n_g, js-n_g:je+n_g, km)
      real(FVPRC), intent(in)::    threshold
      real(FVPRC) qmin
      integer i,j,k

      qmin = q(is,js,1)
      do k=1,km
      do j=js,je
         do i=is,ie
            qmin = min(qmin, q(i,j,k))
         enddo
      enddo
      enddo
      call mp_reduce_min(qmin)
      if(is_master() .and. qmin<threshold) write(6,*) qname, ' min (negative) = ', qmin

 end subroutine prt_negative


end module fv_sg_mod
