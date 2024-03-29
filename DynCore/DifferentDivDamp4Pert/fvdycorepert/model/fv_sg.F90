module fv_sg_mod

!-----------------------------------------------------------------------
! FV sub-grid mixing
!-----------------------------------------------------------------------
use constants_mod, only: rdgas, rvgas, cp_air, hlv, hlf, kappa, grav
use fv_mp_mod,     only: gid

implicit none
private

integer:: irad = 0
public  fv_dry_conv, fv_sg_conv, qsmith, neg_adj3
public  fv_olr, fv_abs_sw, irad

      INTERFACE qsmith
        MODULE PROCEDURE qsmith_r4
        MODULE PROCEDURE qsmith_r8
      END INTERFACE

real, allocatable:: fv_olr(:,:), fv_abs_sw(:,:)

  real, parameter:: esl = 0.621971831
  real, parameter:: tice = 273.16
  real, parameter:: zvir =  rvgas/rdgas - 1.     ! = 0.607789855
  real, allocatable:: table(:),des(:)

contains

 subroutine fv_dry_conv( isd, ied, jsd, jed, is, ie, js, je, km, nq, dt,    &
                         tau, delp, pe, peln, pkz, ta, qa, ua, va,  &
                         hydrostatic, w, delz, u_dt, v_dt, t_dt, q_dt )
! Dry convective adjustment-mixing
!-------------------------------------------
      integer, intent(in):: is, ie, js, je, km, nq
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: tau         ! Relaxation time scale
      real, intent(in):: dt             ! model time step
      real, intent(in)::   pe(is-1:ie+1,km+1,js-1:je+1) 
      real, intent(in):: peln(is  :ie,  km+1,js  :je)
      real, intent(in):: delp(isd:ied,jsd:jed,km)      ! Delta p at each model level
      real, intent(in)::  pkz(is:ie,js:je,km)      ! Delta p at each model level
      real, intent(in):: delz(is:ie,js:je,km)      ! Delta p at each model level
      logical, intent(in)::  hydrostatic
! 
      real, intent(inout):: ua(isd:ied,jsd:jed,km)
      real, intent(inout):: va(isd:ied,jsd:jed,km)
      real, intent(inout)::  w(isd:ied,jsd:jed,km)      ! Delta p at each model level
      real, intent(inout):: ta(isd:ied,jsd:jed,km)      ! Temperature
      real, intent(inout):: qa(isd:ied,jsd:jed,km,nq)   ! Specific humidity & tracers
! Output:
      real, intent(out):: u_dt(isd:ied,jsd:jed,km) 
      real, intent(out):: v_dt(isd:ied,jsd:jed,km) 
      real, intent(out):: t_dt(is:ie,js:je,km) 
      real, intent(out):: q_dt(is:ie,js:je,km,nq) 
!---------------------------Local variables-----------------------------
      real, dimension(is:ie,km):: u0, v0, w0, t0, hd, te, gz, tvm, pm
      real q0(is:ie,km,nq) 
      real gzh(is:ie)
      real ri, pt1, pt2, ratio, tv, cv
      real qmix, h0, mc, fra, rk, rz, rcv, rdt
      real qs1, qs2, lf, dh, dhs
      integer mcond
      integer i, j, k, kk, n, m, iq
      real, parameter:: ustar2 = 1.E-8

        rz = rvgas - rdgas          ! rz = zvir * rdgas
        rk = cp_air/rdgas + 1.
        cv = cp_air - rdgas
       rcv = 1./cv

      rdt = 1./ dt

!------------------------------------------------------------------------
! The nonhydrostatic pressure changes if there is heating (under constant
! volume and mass is locally conserved).
!------------------------------------------------------------------------
   mcond = 1
   m = 3
   fra = dt/real(tau)

  do 1000 j=js,je       ! this main loop can be OpneMPed in j

    do iq=1,nq
       do k=mcond,km
          do i=is,ie
             q0(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo
    enddo

    do k=mcond,km
       do i=is,ie
          t0(i,k) = ta(i,j,k)
          u0(i,k) = ua(i,j,k)
          v0(i,k) = va(i,j,k)
          pm(i,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
       enddo
    enddo

    do i=is,ie
       gzh(i) = 0.
    enddo

    if( hydrostatic ) then
       do k=km, mcond,-1
          do i=is,ie
           tvm(i,k) = t0(i,k)*(1.+zvir*q0(i,k,1))
                tv  = rdgas*tvm(i,k)
            gz(i,k) = gzh(i) + tv*(1.-pe(i,k,j)/pm(i,k))
            hd(i,k) = cp_air*tvm(i,k)+gz(i,k)+0.5*(u0(i,k)**2+v0(i,k)**2)
             gzh(i) = gzh(i) + tv*(peln(i,k+1,j)-peln(i,k,j))
          enddo
       enddo
    else
       do k=km,mcond,-1
          do i=is,ie
             w0(i,k) = w(i,j,k)
             gz(i,k) = gzh(i)  - 0.5*grav*delz(i,j,k)
                 tv  = gz(i,k) + 0.5*(u0(i,k)**2+v0(i,k)**2+w0(i,k)**2)
             hd(i,k) = cp_air*t0(i,k) + tv
             te(i,k) =     cv*t0(i,k) + tv
              gzh(i) = gzh(i) - grav*delz(i,j,k)
          enddo
       enddo
    endif

   do n=1,m

      ratio = real(n)/real(m)

      do i=is,ie
         gzh(i) = 0.
      enddo

      do k=km,mcond+1,-1

         do i=is,ie
! Richardson number = g*delz * theta / ( del_theta * (del_u**2 + del_v**2) )
            pt1 = t0(i,k-1)/pkz(i,j,k-1)
            pt2 = t0(i,k  )/pkz(i,j,k  )
            ri = (gz(i,k-1)-gz(i,k))*(pt1-pt2)/( 0.5*(pt1+pt2)*        &
                ((u0(i,k-1)-u0(i,k))**2+(v0(i,k-1)-v0(i,k))**2+ustar2) )
! Dry convective adjustment for K-H instability:
! Compute equivalent mass flux: mc
            if ( ri < 0.25 ) then
                 mc = (1.-4.*max(0.0,ri)) ** 2
                 mc = ratio*mc*delp(i,j,k-1)*delp(i,j,k)/(delp(i,j,k-1)+delp(i,j,k))
                 do iq=1,nq
                    h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                    q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                    q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                 enddo
! u:
                 h0 = mc*(u0(i,k)-u0(i,k-1))
                 u0(i,k-1) = u0(i,k-1) + h0/delp(i,j,k-1)
                 u0(i,k  ) = u0(i,k  ) - h0/delp(i,j,k  )
! v:
                 h0 = mc*(v0(i,k)-v0(i,k-1))
                 v0(i,k-1) = v0(i,k-1) + h0/delp(i,j,k-1)
                 v0(i,k  ) = v0(i,k  ) - h0/delp(i,j,k  )
              if ( hydrostatic ) then
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
            t0(i,kk) = t0(i,kk) / ( rdgas + rz*q0(i,kk,1) )
         enddo
         kk = k-1
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ((rk-pe(i,kk,j)/pm(i,kk))*(rdgas+rz*q0(i,kk,1)))
         enddo
       else
! Non-hydrostatic under constant volume heating/cooling
         do kk=k-1,k
            do i=is,ie
                     tv = gz(i,kk) + 0.5*(u0(i,kk)**2+v0(i,kk)**2+w0(i,kk)**2)
               t0(i,kk) = rcv*(te(i,kk)- tv)
               hd(i,kk) = cp_air*t0(i,kk) + tv
            enddo
         enddo
       endif
      enddo   ! k-loop

   enddo       ! n-loop


!--------------------
   if ( fra < 1. ) then
      do k=mcond,km
         do i=is,ie
            t0(i,k) = ta(i,j,k) + (t0(i,k) - ta(i,j,k))*fra
            u0(i,k) = ua(i,j,k) + (u0(i,k) - ua(i,j,k))*fra
            v0(i,k) = va(i,j,k) + (v0(i,k) - va(i,j,k))*fra
         enddo
      enddo

      if ( .not. hydrostatic ) then
         do k=mcond,km
            do i=is,ie
               w0(i,k) = w(i,j,k) + (w0(i,k) - w(i,j,k))*fra
            enddo
         enddo
      endif

      do iq=1,nq
         do k=mcond,km
            do i=is,ie
               q0(i,k,iq) = qa(i,j,k,iq) + (q0(i,k,iq) - qa(i,j,k,iq))*fra
            enddo
         enddo
      enddo
   endif
!--------------------

   if ( mcond/=1 ) then
   do k=1,mcond-1
      do i=is,ie
         u_dt(i,j,k) = 0.
         v_dt(i,j,k) = 0.
         t_dt(i,j,k) = 0.
      enddo
      do iq=1,nq
         do i=is,ie
            q_dt(i,j,k,iq) = 0.
         enddo
      enddo
   enddo
   endif

   do k=mcond,km
      do i=is,ie
         u_dt(i,j,k) = rdt*(u0(i,k) - ua(i,j,k))
         v_dt(i,j,k) = rdt*(v0(i,k) - va(i,j,k))
           ta(i,j,k) = t0(i,k)   ! *** temperature updated ***
         t_dt(i,j,k) = 0.
      enddo
   enddo

   do iq=1,nq
      do k=mcond,km
         do i=is,ie
            q_dt(i,j,k,iq) = rdt*(q0(i,k,iq)-qa(i,j,k,iq))
         enddo
      enddo
   enddo

   if ( .not. hydrostatic ) then
      do k=mcond,km
         do i=is,ie
            w(i,j,k) = w0(i,k)   ! w updated
         enddo
      enddo
   endif

1000 continue


 end subroutine fv_dry_conv


 subroutine fv_sg_conv( isd, ied, jsd, jed, is, ie, js, je, km,    &
                        nq, dt, tau, delp, pe, peln, pkz, ta, qa,  &
                        ua, va, hydrostatic, w, delz, u_dt, v_dt, t_dt, q_dt )
! Non-precipitating sub-grid scale convective adjustment-mixing
!-------------------------------------------
      integer, intent(in):: is, ie, js, je, km, nq
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: tau         ! Relaxation time scale
      real, intent(in):: dt             ! model time step
      real, intent(in)::   pe(is-1:ie+1,km+1,js-1:je+1) 
      real, intent(in):: peln(is  :ie,  km+1,js  :je)
      real, intent(in):: delp(isd:ied,jsd:jed,km)      ! Delta p at each model level
      real, intent(in)::  pkz(is:ie,js:je,km)      ! Delta p at each model level
      real, intent(in):: delz(is:ie,js:je,km)      ! Delta p at each model level
      logical, intent(in)::  hydrostatic
! 
      real, intent(inout):: ua(isd:ied,jsd:jed,km)
      real, intent(inout):: va(isd:ied,jsd:jed,km)
      real, intent(inout)::  w(isd:ied,jsd:jed,km)      ! Delta p at each model level
      real, intent(inout):: ta(isd:ied,jsd:jed,km)      ! Temperature
      real, intent(inout):: qa(isd:ied,jsd:jed,km,nq)   ! Specific humidity & tracers
! Output:
      real, intent(out):: u_dt(isd:ied,jsd:jed,km) 
      real, intent(out):: v_dt(isd:ied,jsd:jed,km) 
      real, intent(out):: t_dt(is:ie,js:je,km) 
      real, intent(out):: q_dt(is:ie,js:je,km,nq) 
!---------------------------Local variables-----------------------------
      real, dimension(is:ie,km):: u0, v0, w0, t0, hd, te, gz, tvm, pm
      real q0(is:ie,km,nq) 
      real gzh(is:ie)
      real ri, pt1, pt2, ratio, tv, cv
      real qmix, h0, mc, fra, rk, rz, rcv, rdt
      real qs1, qs2, lf, dh, dhs
      integer mcond
      integer i, j, k, kk, n, m, iq
      real, parameter:: ustar2 = 1.E-8
      real, parameter:: p_crt = 100.E2
      real, parameter:: dh_min = 0.1

        rz = rvgas - rdgas          ! rz = zvir * rdgas
        rk = cp_air/rdgas + 1.
        cv = cp_air - rdgas
       rcv = 1./cv

      rdt = 1./ dt

!------------------------------------------------------------------------
! The nonhydrostatic pressure changes if there is heating (under constant
! volume and mass is locally conserved).
!------------------------------------------------------------------------
   mcond = 1
   m = 4
   fra = dt/real(tau)

  do 1000 j=js,je       ! this main loop can be OpneMPed in j

    do iq=1,nq
       do k=mcond,km
          do i=is,ie
             q0(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo
    enddo

    do k=mcond,km
       do i=is,ie
          t0(i,k) = ta(i,j,k)
          u0(i,k) = ua(i,j,k)
          v0(i,k) = va(i,j,k)
          pm(i,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
       enddo
    enddo

    do i=is,ie
       gzh(i) = 0.
    enddo

    if( hydrostatic ) then
       do k=km, mcond,-1
          do i=is,ie
           tvm(i,k) = t0(i,k)*(1.+zvir*q0(i,k,1))
                tv  = rdgas*tvm(i,k)
            gz(i,k) = gzh(i) + tv*(1.-pe(i,k,j)/pm(i,k))
            hd(i,k) = cp_air*tvm(i,k)+gz(i,k)+0.5*(u0(i,k)**2+v0(i,k)**2)
             gzh(i) = gzh(i) + tv*(peln(i,k+1,j)-peln(i,k,j))
          enddo
       enddo
    else
       do k=km,mcond,-1
          do i=is,ie
             w0(i,k) = w(i,j,k)
             gz(i,k) = gzh(i)  - 0.5*grav*delz(i,j,k)
                 tv  = gz(i,k) + 0.5*(u0(i,k)**2+v0(i,k)**2+w0(i,k)**2)
             hd(i,k) = cp_air*t0(i,k) + tv
             te(i,k) =     cv*t0(i,k) + tv
              gzh(i) = gzh(i) - grav*delz(i,j,k)
          enddo
       enddo
    endif

   do n=1,m

      ratio = real(n)/real(m)

      do i=is,ie
         gzh(i) = 0.
      enddo

      do k=km,mcond+1,-1

         do i=is,ie
! Richardson number = g*delz * theta / ( del_theta * (del_u**2 + del_v**2) )
            pt1 = t0(i,k-1)/pkz(i,j,k-1)
            pt2 = t0(i,k  )/pkz(i,j,k  )
            ri = (gz(i,k-1)-gz(i,k))*(pt1-pt2)/( 0.5*(pt1+pt2)*        &
                ((u0(i,k-1)-u0(i,k))**2+(v0(i,k-1)-v0(i,k))**2+ustar2) )
! Dry convective adjustment for K-H instability:
! Compute equivalent mass flux: mc
            if ( ri < 0.25 ) then
                 mc = (1.-4.*max(0.0,ri)) ** 2
                 mc = ratio*mc*delp(i,j,k-1)*delp(i,j,k)/(delp(i,j,k-1)+delp(i,j,k))
                 do iq=1,nq
                    h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                    q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                    q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                 enddo
! u:
                 h0 = mc*(u0(i,k)-u0(i,k-1))
                 u0(i,k-1) = u0(i,k-1) + h0/delp(i,j,k-1)
                 u0(i,k  ) = u0(i,k  ) - h0/delp(i,j,k  )
! v:
                 h0 = mc*(v0(i,k)-v0(i,k-1))
                 v0(i,k-1) = v0(i,k-1) + h0/delp(i,j,k-1)
                 v0(i,k  ) = v0(i,k  ) - h0/delp(i,j,k  )
              if ( hydrostatic ) then
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
            t0(i,kk) = t0(i,kk) / ( rdgas + rz*q0(i,kk,1) )
         enddo
         kk = k-1
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ((rk-pe(i,kk,j)/pm(i,kk))*(rdgas+rz*q0(i,kk,1)))
         enddo
       else
! Non-hydrostatic under constant volume heating/cooling
         do kk=k-1,k
            do i=is,ie
                     tv = gz(i,kk) + 0.5*(u0(i,kk)**2+v0(i,kk)**2+w0(i,kk)**2)
               t0(i,kk) = rcv*(te(i,kk)- tv)
               hd(i,kk) = cp_air*t0(i,kk) + tv
            enddo
         enddo
       endif
      enddo   ! k-loop
   enddo       ! n-loop


!--------------
! Moist mixing:
!--------------

    if( .not. allocated(table) ) then
       call  qsmith_init
    endif

    do i=is,ie
       gzh(i) = 0.
    enddo

    if( hydrostatic ) then
       do k=km, mcond,-1
          do i=is,ie
           tvm(i,k) = t0(i,k)*(1.+zvir*q0(i,k,1))
                tv  = rdgas*tvm(i,k)
            gz(i,k) = gzh(i) + tv*(1.-pe(i,k,j)/pm(i,k))
            hd(i,k) = cp_air*tvm(i,k)+gz(i,k)+0.5*(u0(i,k)**2+v0(i,k)**2)
             gzh(i) = gzh(i) + tv*(peln(i,k+1,j)-peln(i,k,j))
          enddo
       enddo
    else
       do k=km,mcond,-1
          do i=is,ie
             gz(i,k) = gzh(i)  - 0.5*grav*delz(i,j,k)
                 tv  = gz(i,k) + 0.5*(u0(i,k)**2+v0(i,k)**2+w0(i,k)**2)
             hd(i,k) = cp_air*t0(i,k) + tv
             te(i,k) =     cv*t0(i,k) + tv
              gzh(i) = gzh(i) - grav*delz(i,j,k)
          enddo
       enddo
    endif

   do n=1,m
      ratio = 0.5*real(n)/real(m)

      do i=is,ie
         gzh(i) = 0.
      enddo

      do k=km,mcond+1,-1
         do i=is,ie
            if ( pm(i,k) > p_crt ) then
               qs1 = qs1d(t0(i,k-1), pm(i,k-1), q0(i,k-1,1))
!              qs2 = qs1d(t0(i,k  ), pm(i,k  ), q0(i,k  ,1))
!           if ( q0(i,k-1,1)>qs1 .and. q0(i,k,1)>qs2 ) then
               lf = hlv + hlf*min(1., max(0., (tice-t0(i,k-1))/30.))
              dh  = hd(i,k) - hd(i,k-1)
              dhs = dh + lf*(q0(i,k,1)-qs1        )
              dh  = dh + lf*(q0(i,k,1)-q0(i,k-1,1))

              if ( dh>dh_min .and. dhs>dh_min ) then   ! layer above is also saturated
                   mc = delp(i,j,k)  *     &
                        min( ratio*dhs/dh, delp(i,j,k-1)/(delp(i,j,k-1)+delp(i,j,k)) )
! Perform local mixing of all advected tracers:
                   do iq=1,nq
                                h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                      q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                      q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                   enddo
                          h0 = mc*(u0(i,k)-u0(i,k-1))
                   u0(i,k-1) = u0(i,k-1) + h0/delp(i,j,k-1)
                   u0(i,k  ) = u0(i,k  ) - h0/delp(i,j,k  )
                          h0 = mc*(v0(i,k)-v0(i,k-1))
                   v0(i,k-1) = v0(i,k-1) + h0/delp(i,j,k-1)
                   v0(i,k  ) = v0(i,k  ) - h0/delp(i,j,k  )
              if ( hydrostatic ) then
                          h0 = mc*(hd(i,k)-hd(i,k-1))
                   hd(i,k-1) = hd(i,k-1) + h0/delp(i,j,k-1)
                   hd(i,k  ) = hd(i,k  ) - h0/delp(i,j,k  )
              else
                          h0 = mc*(hd(i,k)-hd(i,k-1))
                   te(i,k-1) = te(i,k-1) + h0/delp(i,j,k-1)
                   te(i,k  ) = te(i,k  ) - h0/delp(i,j,k  )
                          h0 = mc*(w0(i,k)-w0(i,k-1))
                   w0(i,k-1) = w0(i,k-1) + h0/delp(i,j,k-1)
                   w0(i,k  ) = w0(i,k  ) - h0/delp(i,j,k  )
              endif
              endif  ! dh check
!           endif    ! qs check
            endif    ! p_crt check
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
            t0(i,kk) = t0(i,kk) / ( rdgas + rz*q0(i,kk,1) )
         enddo
         kk = k-1
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ((rk-pe(i,kk,j)/pm(i,kk))*(rdgas+rz*q0(i,kk,1)))
         enddo
       else
! Non-hydrostatic under constant volume heating/cooling
         do kk=k-1,k
            do i=is,ie
                     tv = gz(i,kk) + 0.5*(u0(i,kk)**2+v0(i,kk)**2+w0(i,kk)**2)
               t0(i,kk) = rcv*(te(i,kk)- tv)
               hd(i,kk) = cp_air*t0(i,kk) + tv
            enddo
         enddo
       endif
      enddo    ! k-loop
   enddo       ! n-loop


   if ( fra < 1. ) then
      do k=mcond,km
         do i=is,ie
            t0(i,k) = ta(i,j,k) + (t0(i,k) - ta(i,j,k))*fra
            u0(i,k) = ua(i,j,k) + (u0(i,k) - ua(i,j,k))*fra
            v0(i,k) = va(i,j,k) + (v0(i,k) - va(i,j,k))*fra
         enddo
      enddo

      if ( .not. hydrostatic ) then
         do k=mcond,km
            do i=is,ie
               w0(i,k) = w(i,j,k) + (w0(i,k) - w(i,j,k))*fra
            enddo
         enddo
      endif

      do iq=1,nq
         do k=mcond,km
            do i=is,ie
               q0(i,k,iq) = qa(i,j,k,iq) + (q0(i,k,iq) - qa(i,j,k,iq))*fra
            enddo
         enddo
      enddo
   endif
!--------------------

   if ( mcond/=1 ) then
   do k=1,mcond-1
      do i=is,ie
         u_dt(i,j,k) = 0.
         v_dt(i,j,k) = 0.
         t_dt(i,j,k) = 0.
      enddo
      do iq=1,nq
         do i=is,ie
            q_dt(i,j,k,iq) = 0.
         enddo
      enddo
   enddo
   endif

   do k=mcond,km
      do i=is,ie
         u_dt(i,j,k) = rdt*(u0(i,k) - ua(i,j,k))
         v_dt(i,j,k) = rdt*(v0(i,k) - va(i,j,k))
           ta(i,j,k) = t0(i,k)   ! temperature updated
         t_dt(i,j,k) = 0.
      enddo
   enddo

   do iq=1,nq
      do k=mcond,km
         do i=is,ie
            q_dt(i,j,k,iq) = rdt*(q0(i,k,iq)-qa(i,j,k,iq))
         enddo
      enddo
   enddo

   if ( .not. hydrostatic ) then
      do k=mcond,km
         do i=is,ie
            w(i,j,k) = w0(i,k)   ! w updated
         enddo
      enddo
   endif

1000 continue


 end subroutine fv_sg_conv



 real function qs1d(t, p, q)
! Based on "moist" mixing ratio, p is the total (dry+vapor) pressure
  real, intent(in):: t, p, q
! Local:
  real es, ap1
  real, parameter:: Tmin=tice - 160.
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


  subroutine qsmith_r8(im, km, k1, t, p, q, qs, dqdt)
! input T in deg K; p (Pa)
  integer, intent(in):: im, km, k1
  real*8, intent(in),dimension(im,km):: t, p, q
  real*8, intent(out),dimension(im,km):: qs
  real*8, intent(out), optional:: dqdt(im,km)
! Local:
  real*8 es(im,km)
  real*8 ap1, eps10
  real*8 Tmin
  integer i, k, it

  Tmin = tice-160.
  eps10  = 10.*esl

  if( .not. allocated(table) ) then
       call  qsmith_init
  endif
 
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
 
  end subroutine qsmith_r8
 
  subroutine qsmith_r4(im, km, k1, t, p, q, qs, dqdt)
! input T in deg K; p (Pa)
  integer, intent(in):: im, km, k1
  real*4, intent(in),dimension(im,km):: t, p, q
  real*4, intent(out),dimension(im,km):: qs
  real*4, intent(out), optional:: dqdt(im,km)
! Local: 
  real*4 es(im,km)
  real*4 ap1, eps10  
  real*4 Tmin
  integer i, k, it
            
  Tmin = tice-160.
  eps10  = 10.*esl 
       
  if( .not. allocated(table) ) then
       call  qsmith_init
  endif
  
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

  end subroutine qsmith_r4


  subroutine qs_table(n,table)
      integer, intent(in):: n
      real table (n)
      real esupc(200)
      real:: dt=0.1
      real esbasw, tbasw, esbasi, tbasi, Tmin, tem, aa, b, c, d, e, esh20 
      real wice, wh2o
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
         b   = -3.56654 *alog10(tbasi/tem)
         c   =  0.876793*(1.0-tem/tbasi)
         e   = alog10(esbasi)
         table(i)=10**(aa+b+c+e)
      enddo
! *****************************************************
!  Compute es over water between -20c and 102c.
!  see smithsonian meteorological tables page 350.
      do  i=1,1221
          tem = 253.16+dt*real(i-1)
          aa  = -7.90298*(tbasw/tem-1)
          b   =  5.02808*alog10(tbasw/tem)
          c   = -1.3816e-07*(10**((1-tem/tbasw)*11.344)-1)
          d   =  8.1328e-03*(10**((tbasw/tem-1)*(-3.49149))-1)
          e   = alog10(esbasw)
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

 end subroutine qs_table

 subroutine neg_adj3(is, ie, js, je, ng, kbot,      &
                     pt, dp, qv, ql, qr, qi, qs, qg, qa)

! This is designed for 6-class micro-physics schemes
 integer, intent(in):: is, ie, js, je, ng, kbot
 real, intent(in):: dp(is-ng:ie+ng,js-ng:je+ng,kbot)
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,kbot)::    &
                                 pt, qv, ql, qr, qi, qs, qg
 real, intent(inout), optional, dimension(is-ng:ie+ng,js-ng:je+ng,kbot):: qa
! Local:
 real lcp, icp
 real dq, qsum, psum
 integer i, j, k

 lcp = hlv / cp_air
 icp = hlf / cp_air

 do k=1, kbot
    do j=js, je
       do i=is, ie
!-----------
! Ice-phase:
!-----------
! if ice<0 borrow from snow (since main source of snow is cloud ice)
          if( qi(i,j,k) < 0. ) then
              qs(i,j,k) = qs(i,j,k) + qi(i,j,k)
              qi(i,j,k) = 0.
          endif
! if snow<0 borrow from graupel (same as above)
          if( qs(i,j,k) < 0. ) then
              qg(i,j,k) = qg(i,j,k) + qs(i,j,k)
              qs(i,j,k) = 0.
          endif
! If graupel < 0 then borrow from cloud ice (loop back)
          if ( qg(i,j,k) < 0. ) then
               qi(i,j,k) = qi(i,j,k) + qg(i,j,k)
               qg(i,j,k) = 0.
          endif

! If ice < 0 then borrow from cloud water
          if ( qi(i,j,k) < 0. ) then
               ql(i,j,k) = ql(i,j,k) + qi(i,j,k)
               pt(i,j,k) = pt(i,j,k) - qi(i,j,k)*icp   ! heating
               qi(i,j,k) = 0.
          endif

! Liquid phase:
! Fix negative cloud water by borrowing from rain
          if ( ql(i,j,k) < 0. ) then
               qr(i,j,k) = qr(i,j,k) + ql(i,j,k)
               ql(i,j,k) = 0.
          endif
! fix negative rain with vapor
          if ( qr(i,j,k) < 0. ) then
               qv(i,j,k) = qv(i,j,k) + qr(i,j,k)
               pt(i,j,k) = pt(i,j,k) - qr(i,j,k)*lcp
               qr(i,j,k) = 0.
          endif
     enddo
   enddo
 enddo

!-----------------------------------
! Fix water vapor
!-----------------------------------
! Top layer: borrow from below
    k = 1
    do j=js, je
       do i=is, ie
          if( qv(i,j,k) < 0. ) then
              qv(i,j,k+1) = qv(i,j,k+1) + qv(i,j,k)*dp(i,j,k)/dp(i,j,k+1)
              qv(i,j,k  ) = 0.
          endif
     enddo
   enddo

 do k=2,kbot-1
    do j=js, je
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
  do j=js, je
     do i=is, ie
        if( qv(i,j,kbot) < 0. .and. qv(i,j,kbot-1)>0.) then
            dq = min(-qv(i,j,kbot)*dp(i,j,kbot), qv(i,j,kbot-1)*dp(i,j,kbot-1))
            qv(i,j,kbot-1) = qv(i,j,kbot-1) - dq/dp(i,j,kbot-1) 
            qv(i,j,kbot  ) = qv(i,j,kbot  ) + dq/dp(i,j,kbot  ) 
        endif
! Last attempt to fix negative qv from condensates
!       if( qv(i,j,kbot) < 0. ) then
!       endif
   enddo
 enddo

!-----------------------------------
! Fix negative cloud fraction
!-----------------------------------
 if ( present(qa) ) then
 do k=1,kbot-1
    do j=js, je
       do i=is, ie
          if( qa(i,j,k) < 0. ) then
              qa(i,j,k+1) = qa(i,j,k+1) + qa(i,j,k)*dp(i,j,k)/dp(i,j,k+1)
              qa(i,j,k  ) = 0.
          endif
     enddo
   enddo
 enddo
 
! Bottom layer; Borrow from above
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
 endif

 end subroutine neg_adj3

end module fv_sg_mod
