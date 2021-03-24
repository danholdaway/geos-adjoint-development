!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
module fv_sg_mod

!-----------------------------------------------------------------------
! FV sub-grid mixing
!-----------------------------------------------------------------------
  use constants_mod,      only: rdgas, rvgas, cp_air, cp_vapor, hlv, hlf, kappa, grav
!  use tracer_manager_mod, only: get_tracer_index
!  use field_manager_mod,  only: MODEL_ATMOS
!  use lin_cld_microphys_mod, only: wqs2, wqsat2_moist
!  use fv_mp_mod,          only: mp_reduce_min, is_master

implicit none
private

public  fv_subgrid_z, neg_adj3

  real, parameter:: esl = 0.621971831
  real, parameter:: tice = 273.16
! real, parameter:: c_ice = 2106.  ! Emanuel table, page 566
  real, parameter:: c_ice = 1972.  !  -15 C
  real, parameter:: c_liq = 4.1855e+3    ! GFS
! real, parameter:: c_liq = 4218.        ! ECMWF-IFS
  real, parameter:: cv_vap = cp_vapor - rvgas  ! 1384.5
  real, parameter:: c_con = c_ice

! real, parameter:: dc_vap =  cp_vapor - c_liq   ! = -2368.
  real, parameter:: dc_vap =  cv_vap - c_liq   ! = -2368.
  real, parameter:: dc_ice =  c_liq - c_ice      ! = 2112.
! Values at 0 Deg C
  real, parameter:: hlv0 = 2.5e6
  real, parameter:: hlf0 = 3.3358e5
! real, parameter:: hlv0 = 2.501e6   ! Emanual Appendix-2
! real, parameter:: hlf0 = 3.337e5   ! Emanual
  real, parameter:: t_ice = 273.16
  real, parameter:: ri_max = 1.
  real, parameter:: ri_min = 0.25
  real, parameter:: t1_min = 160.
  real, parameter:: t2_min = 165.
  real, parameter:: t2_max = 315.
  real, parameter:: t3_max = 325.
  real, parameter:: Lv0 =  hlv0 - dc_vap*t_ice   ! = 3.147782e6
  real, parameter:: Li0 =  hlf0 - dc_ice*t_ice   ! = -2.431928e5

  real, parameter:: zvir =  rvgas/rdgas - 1.     ! = 0.607789855
  real, allocatable:: table(:),des(:)
  real:: lv00, d0_vap


contains


 subroutine fv_subgrid_z( isd, ied, jsd, jed, is, ie, js, je, km, nq, dt,    &
                         tau, nwat, delp, pe, peln, pkz, ta, qa, ua, va,  &
                         hydrostatic, w, delz, u_dt, v_dt, t_dt, k_bot )
! Dry convective adjustment-mixing
!-------------------------------------------
      integer, intent(in):: is, ie, js, je, km, nq, nwat
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: tau         ! Relaxation time scale
      real, intent(in):: dt             ! model time step
      real, intent(in)::   pe(is-1:ie+1,km+1,js-1:je+1) 
      real, intent(in):: peln(is  :ie,  km+1,js  :je)
      real, intent(in):: delp(isd:ied,jsd:jed,km)      ! Delta p at each model level
      real, intent(in):: delz(isd:,jsd:,1:)      ! Delta z at each model level
      real, intent(in)::  pkz(is:ie,js:je,km)
      logical, intent(in)::  hydrostatic
      integer, intent(in), optional:: k_bot
!
      real, intent(inout):: ua(isd:ied,jsd:jed,km)
      real, intent(inout):: va(isd:ied,jsd:jed,km)
      real, intent(inout)::  w(isd:,jsd:,1:)
      real, intent(inout):: ta(isd:ied,jsd:jed,km)      ! Temperature
      real, intent(inout):: qa(isd:ied,jsd:jed,km,nq)   ! Specific humidity & tracers
      real, intent(inout):: u_dt(isd:ied,jsd:jed,km) 
      real, intent(inout):: v_dt(isd:ied,jsd:jed,km) 
      real, intent(inout):: t_dt(is:ie,js:je,km) 
!---------------------------Local variables-----------------------------
      real, dimension(is:ie,km):: u0, v0, w0, t0, hd, te, gz, tvm, pm, den
      real q0(is:ie,km,nq), qcon(is:ie,km) 
      real, dimension(is:ie):: gzh, lcp2, icp2, cvm, cpm, qs
      real ri_ref, ri, pt1, pt2, ratio, tv, cv, tmp, q_liq, q_sol
      real tv1, tv2, g2, h0, mc, fra, rk, rz, rdt, tvd, tv_surf
      real dh, dq, qsw, dqsdt, tcp3, t_max, t_min
      integer i, j, k, kk, n, m, iq, km1, im, kbot
      real, parameter:: ustar2 = 1.E-4
      real:: cv_air, xvir
      integer :: sphum, liq_wat, rainwat, snowwat, graupel, ice_wat, cld_amt

      cv_air = cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68
        rk = cp_air/rdgas + 1.
        cv = cp_air - rdgas

      g2 = 0.5*grav

      rdt = 1./ dt
      im = ie-is+1

      if ( present(k_bot) ) then
           if ( k_bot < 3 ) return
           kbot = k_bot
      else
           kbot = km
      endif
      if ( pe(is,1,js) < 2. ) then
           t_min = t1_min
      else
           t_min = t2_min
      endif

      if ( k_bot < min(km,24)  ) then
         t_max = t2_max
      else
         t_max = t3_max
      endif


      sphum = 1
      rainwat = -1; snowwat = -1; graupel = -1
      if ( nwat == 0 ) then
         xvir = 0.
         rz = 0.
      else
         xvir = zvir
         rz = rvgas - rdgas          ! rz = zvir * rdgas
         if ( nwat == 3) then
            liq_wat = 2
            ice_wat = 3
         endif
      endif


!------------------------------------------------------------------------
! The nonhydrostatic pressure changes if there is heating (under constant
! volume and mass is locally conserved).
!------------------------------------------------------------------------
   m = 3
   fra = dt/real(tau)

!$OMP parallel do default(none) shared(im,is,ie,js,je,nq,kbot,qa,ta,sphum,ua,va,delp,peln,   &
!$OMP                                  hydrostatic,pe,delz,g2,w,liq_wat,rainwat,ice_wat,     &
!$OMP                                  snowwat,cv_air,m,graupel,pkz,rk,rz,fra, t_max, t_min, &
!$OMP                                  u_dt,rdt,v_dt,xvir,nwat)                              &
!$OMP                          private(kk,lcp2,icp2,tcp3,dh,dq,den,qs,qsw,dqsdt,qcon,q0,     &
!$OMP                                  t0,u0,v0,w0,h0,pm,gzh,tvm,tmp,cpm,cvm,q_liq,q_sol,    &
!$OMP                                  tv,gz,hd,te,ratio,pt1,pt2,tv1,tv2,ri_ref, ri,mc,km1)
  do 1000 j=js,je  

    do iq=1, nq
       do k=1,kbot
          do i=is,ie
             q0(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo
    enddo

    do k=1,kbot
       do i=is,ie
          t0(i,k) = ta(i,j,k)
         tvm(i,k) = t0(i,k)*(1.+xvir*q0(i,k,sphum))
          u0(i,k) = ua(i,j,k)
          v0(i,k) = va(i,j,k)
          pm(i,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
       enddo
    enddo

    do i=is,ie
       gzh(i) = 0.
    enddo

    if( hydrostatic ) then
       do k=kbot, 1,-1
          do i=is,ie
                tv  = rdgas*tvm(i,k)
           den(i,k) = pm(i,k)/tv
            gz(i,k) = gzh(i) + tv*(1.-pe(i,k,j)/pm(i,k))
            hd(i,k) = cp_air*tvm(i,k)+gz(i,k)+0.5*(u0(i,k)**2+v0(i,k)**2)
             gzh(i) = gzh(i) + tv*(peln(i,k+1,j)-peln(i,k,j))
          enddo
       enddo
    else
       do k=kbot, 1, -1
       if ( nwat == 0 ) then
          do i=is,ie
             cpm(i) = cp_air
             cvm(i) = cv_air
          enddo
       elseif ( nwat==1 ) then
          do i=is,ie
             cpm(i) = (1.-q0(i,k,sphum))*cp_air + q0(i,k,sphum)*cp_vapor
             cvm(i) = (1.-q0(i,k,sphum))*cv_air + q0(i,k,sphum)*cv_vap
          enddo
       elseif ( nwat==2 ) then   ! GFS
          do i=is,ie
             cpm(i) = (1.-q0(i,k,sphum))*cp_air + q0(i,k,sphum)*cp_vapor
             cvm(i) = (1.-q0(i,k,sphum))*cv_air + q0(i,k,sphum)*cv_vap
          enddo
       elseif ( nwat==3 ) then
          do i=is,ie
             q_liq = q0(i,k,liq_wat) 
             q_sol = q0(i,k,ice_wat)
             cpm(i) = (1.-(q0(i,k,sphum)+q_liq+q_sol))*cp_air + q0(i,k,sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
             cvm(i) = (1.-(q0(i,k,sphum)+q_liq+q_sol))*cv_air + q0(i,k,sphum)*cv_vap   + q_liq*c_liq + q_sol*c_ice
          enddo
       elseif ( nwat==4 ) then
          do i=is,ie
             q_liq = q0(i,k,liq_wat) + q0(i,k,rainwat)
             cpm(i) = (1.-(q0(i,k,sphum)+q_liq))*cp_air + q0(i,k,sphum)*cp_vapor + q_liq*c_liq
             cvm(i) = (1.-(q0(i,k,sphum)+q_liq))*cv_air + q0(i,k,sphum)*cv_vap   + q_liq*c_liq
          enddo
       else
          do i=is,ie
             q_liq = q0(i,k,liq_wat) + q0(i,k,rainwat)
             q_sol = q0(i,k,ice_wat) + q0(i,k,snowwat) + q0(i,k,graupel)
             cpm(i) = (1.-(q0(i,k,sphum)+q_liq+q_sol))*cp_air + q0(i,k,sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
             cvm(i) = (1.-(q0(i,k,sphum)+q_liq+q_sol))*cv_air + q0(i,k,sphum)*cv_vap   + q_liq*c_liq + q_sol*c_ice
          enddo
       endif

          do i=is,ie
           den(i,k) = -delp(i,j,k)/(grav*delz(i,j,k))
             w0(i,k) = w(i,j,k)
             gz(i,k) = gzh(i)  - g2*delz(i,j,k)
                tmp  = gz(i,k) + 0.5*(u0(i,k)**2+v0(i,k)**2+w0(i,k)**2)
             hd(i,k) = cpm(i)*t0(i,k) + tmp
             te(i,k) = cvm(i)*t0(i,k) + tmp
              gzh(i) = gzh(i) - grav*delz(i,j,k)
          enddo
       enddo
    endif

   do n=1,m

     if ( m==3 ) then
        if ( n==1) ratio = 0.25
        if ( n==2) ratio = 0.5
        if ( n==3) ratio = 0.999
     else
      ratio = real(n)/real(m)
     endif

      do i=is,ie
         gzh(i) = 0.
      enddo

! Compute total condensate
   if ( nwat<2 ) then
      do k=1,kbot
         do i=is,ie
            qcon(i,k) = 0.
         enddo
      enddo
   elseif ( nwat==2 ) then   ! GFS_2015
      do k=1,kbot
         do i=is,ie
            qcon(i,k) = q0(i,k,liq_wat)
         enddo
      enddo
   elseif ( nwat==3 ) then
      do k=1,kbot
         do i=is,ie
            qcon(i,k) = q0(i,k,liq_wat) + q0(i,k,ice_wat)
         enddo
      enddo
   elseif ( nwat==4 ) then
      do k=1,kbot
         do i=is,ie
            qcon(i,k) = q0(i,k,liq_wat) + q0(i,k,rainwat)
         enddo
      enddo
   else
      do k=1,kbot
         do i=is,ie
            qcon(i,k) = q0(i,k,liq_wat)+q0(i,k,ice_wat)+q0(i,k,snowwat)+q0(i,k,rainwat)+q0(i,k,graupel)
         enddo
      enddo
   endif

      do k=kbot, 2, -1
         km1 = k-1
         do i=is,ie
! Richardson number = g*delz * del_theta/theta / (del_u**2 + del_v**2)
! Use exact form for "density temperature"
            tv1 = t0(i,km1)*(1.+xvir*q0(i,km1,sphum)-qcon(i,km1))
            tv2 = t0(i,k  )*(1.+xvir*q0(i,k  ,sphum)-qcon(i,k))
            pt1 = tv1 / pkz(i,j,km1)
            pt2 = tv2 / pkz(i,j,k  )
!
            ri = (gz(i,km1)-gz(i,k))*(pt1-pt2)/( 0.5*(pt1+pt2)*        &
                 ((u0(i,km1)-u0(i,k))**2+(v0(i,km1)-v0(i,k))**2+ustar2) )
            if ( tv1>t_max .and. tv1>tv2 ) then
! top layer unphysically warm
               ri = 0.
            elseif ( tv2<t_min ) then
               ri = min(ri, 0.2)
            endif
! Adjustment for K-H instability:
! Compute equivalent mass flux: mc
! Add moist 2-dz instability consideration:
!!!         ri_ref = min(ri_max, ri_min + (ri_max-ri_min)*dim(500.e2,pm(i,k))/250.e2 )
            ri_ref = min(ri_max, ri_min + (ri_max-ri_min)*max(400.e2-pm(i,k),0.)/200.e2 )
! Enhancing mixing at the model top
            if ( k==2 ) then
                 ri_ref = 4.*ri_ref
            elseif ( k==3 ) then
                 ri_ref = 2.*ri_ref
            elseif ( k==4 ) then
                 ri_ref = 1.5*ri_ref
            endif

            if ( ri < ri_ref ) then
               mc = ratio*delp(i,j,km1)*delp(i,j,k)/(delp(i,j,km1)+delp(i,j,k))*(1.-max(0.0,ri/ri_ref))**2
                 do iq=1,nq
                    h0 = mc*(q0(i,k,iq)-q0(i,km1,iq))
                    q0(i,km1,iq) = q0(i,km1,iq) + h0/delp(i,j,km1)
                    q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                 enddo
! Recompute qcon
                 if ( nwat<2 ) then
                    qcon(i,km1) = 0.
                 elseif ( nwat==2 ) then  ! GFS_2015
                    qcon(i,km1) = q0(i,km1,liq_wat)
                 elseif ( nwat==3 ) then  ! AM3/AM4
                    qcon(i,km1) = q0(i,km1,liq_wat) + q0(i,km1,ice_wat)
                 elseif ( nwat==4 ) then  ! K_warm_rain scheme with fake ice
                    qcon(i,km1) = q0(i,km1,liq_wat) + q0(i,km1,rainwat)
                 else
                    qcon(i,km1) = q0(i,km1,liq_wat) + q0(i,km1,ice_wat) +                  &
                                  q0(i,km1,snowwat) + q0(i,km1,rainwat) + q0(i,km1,graupel)
                 endif
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
           if ( nwat == 0 ) then
            do i=is,ie
               cpm(i) = cp_air
               cvm(i) = cv_air
            enddo
           elseif ( nwat == 1 ) then
            do i=is,ie
               cpm(i) = (1.-q0(i,kk,sphum))*cp_air + q0(i,kk,sphum)*cp_vapor
               cvm(i) = (1.-q0(i,kk,sphum))*cv_air + q0(i,kk,sphum)*cv_vap
            enddo
           elseif ( nwat == 2 ) then
            do i=is,ie
               cpm(i) = (1.-q0(i,kk,sphum))*cp_air + q0(i,kk,sphum)*cp_vapor
               cvm(i) = (1.-q0(i,kk,sphum))*cv_air + q0(i,kk,sphum)*cv_vap
            enddo
           elseif ( nwat == 3 ) then
            do i=is,ie
               q_liq = q0(i,kk,liq_wat)
               q_sol = q0(i,kk,ice_wat)
               cpm(i) = (1.-(q0(i,kk,sphum)+q_liq+q_sol))*cp_air + q0(i,kk,sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
               cvm(i) = (1.-(q0(i,kk,sphum)+q_liq+q_sol))*cv_air + q0(i,kk,sphum)*cv_vap   + q_liq*c_liq + q_sol*c_ice
            enddo
           elseif ( nwat == 4 ) then
            do i=is,ie
               q_liq = q0(i,kk,liq_wat) + q0(i,kk,rainwat)
               cpm(i) = (1.-(q0(i,kk,sphum)+q_liq))*cp_air + q0(i,kk,sphum)*cp_vapor + q_liq*c_liq
               cvm(i) = (1.-(q0(i,kk,sphum)+q_liq))*cv_air + q0(i,kk,sphum)*cv_vap   + q_liq*c_liq
            enddo
           else
            do i=is,ie
               q_liq = q0(i,kk,liq_wat) + q0(i,kk,rainwat)
               q_sol = q0(i,kk,ice_wat) + q0(i,kk,snowwat) + q0(i,kk,graupel)
               cpm(i) = (1.-(q0(i,kk,sphum)+q_liq+q_sol))*cp_air + q0(i,kk,sphum)*cp_vapor + q_liq*c_liq + q_sol*c_ice
               cvm(i) = (1.-(q0(i,kk,sphum)+q_liq+q_sol))*cv_air + q0(i,kk,sphum)*cv_vap   + q_liq*c_liq + q_sol*c_ice
            enddo
           endif
     
            do i=is,ie
               tv = gz(i,kk) + 0.5*(u0(i,kk)**2+v0(i,kk)**2+w0(i,kk)**2)
               t0(i,kk) = (te(i,kk)- tv) / cvm(i)
               hd(i,kk) = cpm(i)*t0(i,kk) + tv
            enddo
         enddo
       endif
      enddo   ! k-loop
   enddo       ! n-loop

!--------------------
   if ( fra < 1. ) then
      do k=1, kbot
         do i=is,ie
            t0(i,k) = ta(i,j,k) + (t0(i,k) - ta(i,j,k))*fra
            u0(i,k) = ua(i,j,k) + (u0(i,k) - ua(i,j,k))*fra
            v0(i,k) = va(i,j,k) + (v0(i,k) - va(i,j,k))*fra
         enddo
      enddo

      if ( .not. hydrostatic ) then
         do k=1,kbot
            do i=is,ie
               w0(i,k) = w(i,j,k) + (w0(i,k) - w(i,j,k))*fra
            enddo
         enddo
      endif

      do iq=1,nq
         do k=1,kbot
            do i=is,ie
               q0(i,k,iq) = qa(i,j,k,iq) + (q0(i,k,iq) - qa(i,j,k,iq))*fra
            enddo
         enddo
      enddo
   endif

   do k=1,kbot
      do i=is,ie
         u_dt(i,j,k) = rdt*(u0(i,k) - ua(i,j,k))
         v_dt(i,j,k) = rdt*(v0(i,k) - va(i,j,k))
           ta(i,j,k) = t0(i,k)   ! *** temperature updated ***

           ua(i,j,k) = u0(i,k)
           va(i,j,k) = v0(i,k)

      enddo
      do iq=1,nq
         do i=is,ie
            qa(i,j,k,iq) = q0(i,k,iq)
         enddo
      enddo
   enddo

   if ( .not. hydrostatic ) then
      do k=1,kbot
         do i=is,ie
            w(i,j,k) = w0(i,k)   ! w updated
         enddo
      enddo
   endif

1000 continue

 end subroutine fv_subgrid_z



 subroutine neg_adj3(is, ie, js, je, ng, kbot, hydrostatic,   &
                     peln, delz, pt, dp, qv, ql, qr, qi, qs, qg, qa, check_negative)

! This is designed for 6-class micro-physics schemes
 integer, intent(in):: is, ie, js, je, ng, kbot
 logical, intent(in):: hydrostatic
 real, intent(in):: dp(is-ng:ie+ng,js-ng:je+ng,kbot)  ! total delp-p
 real, intent(in):: delz(is-ng:,js-ng:,1:)
 real, intent(in):: peln(is:ie,kbot+1,js:je)           ! ln(pe)
 logical, intent(in), OPTIONAL :: check_negative
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,kbot)::    &
                                 pt, qv, ql, qr, qi, qs, qg
 real, intent(inout), OPTIONAL, dimension(is-ng:ie+ng,js-ng:je+ng,kbot):: qa

!Stubbed, would require careful checking of nonlinearity

 end subroutine neg_adj3

end module fv_sg_mod
