 subroutine Lagrangian_to_Eulerian(do_consv, consv, ps, pe, delp, pkz, pk,   &
                      pdt, km, is,ie,js,je, isd,ied,jsd,jed,       &
                      nq, sphum, u, v, w, delz, pt, q, hs, r_vir, cp, akap,  &
#ifdef MAPL_MODE
                      pi, radius, grav, &
#endif
                      kord_mt, kord_tr, kord_tm,  peln, te0_2d,              &
                      ng, ua, va, omga, te, pem, fill, reproduce_sum,        &
                      ak, bk, ks, ze0, te_method, remap_t, hydrostatic, hybrid_z, do_omega, ktop, &
                      ncnst, mfx, mfy)

 IMPLICIT NONE

  logical, intent(in):: do_consv
  real,    intent(in):: pdt                   ! phys time step
  integer, intent(in):: km
  integer, intent(in):: nq                    ! number of tracers (including h2o)
  integer, intent(in):: sphum                 ! index for water vapor (specific humidity)
  integer, intent(in):: ng
  integer, intent(in):: is,ie,isd,ied         ! starting & ending X-Dir index
  integer, intent(in):: js,je,jsd,jed         ! starting & ending Y-Dir index
  integer, intent(in):: ks, ktop
  integer, intent(in):: kord_mt               ! Mapping oder for the vector winds
  integer, intent(in):: kord_tr               ! Mapping oder for tracers
  integer, intent(in):: kord_tm               ! Mapping oder for thermodynamics
  integer, intent(in):: ncnst                 ! Added 4TAF

  real, intent(in):: consv                 ! factor for TE conservation
  real, intent(in):: r_vir
  real, intent(in):: cp
  real, intent(in):: akap
#ifdef MAPL_MODE
  real, intent(in):: pi, radius, grav
#endif
  real, intent(in):: hs(isd:ied,jsd:jed)  ! surface geopotential
  real, intent(in):: te0_2d(is:ie,js:je)

  logical, intent(in):: fill                  ! fill negative tracers
  logical, intent(in):: reproduce_sum
  logical, intent(in):: do_omega
  real, intent(in) :: ak(km+1)
  real, intent(in) :: bk(km+1)

! !INPUT/OUTPUT
  real(p_precision), intent(inout):: pk(is:ie,js:je,km+1) ! pe to the kappa
  real, intent(inout):: q(isd:ied,jsd:jed,km,*)
  real, intent(inout):: delp(isd:ied,jsd:jed,km) ! pressure thickness
  real(p_precision), intent(inout)::  pe(is-1:ie+1,km+1,js-1:je+1) ! pressure at layer edges
  real, intent(inout):: pem(is-1:ie+1,km+1,js-1:je+1)
  real(p_precision), intent(inout):: ps(isd:ied,jsd:jed)      ! surface pressure
  real, intent(inout):: ze0(is:ie,js:je,km+1)    ! Specified height at edges (m)

! u-wind will be ghosted one latitude to the north upon exit
  real, intent(inout)::  u(isd:ied  ,jsd:jed+1,km)   ! u-wind (m/s)
  real, intent(inout)::  v(isd:ied+1,jsd:jed  ,km)   ! v-wind (m/s)
  real, intent(inout)::  w(isd:ied  ,jsd:jed  ,km)   ! vertical velocity (m/s)
  real, intent(inout):: pt(isd:ied  ,jsd:jed  ,km)   ! cp*virtual potential temperature 
                                                     ! as input; output: temperature
  real, intent(inout):: delz(is:ie,js:je,km)   ! delta-height (m)
  integer, intent(in):: te_method
  logical, intent(in):: remap_t
  logical, intent(in):: hydrostatic
  logical, intent(in):: hybrid_z

  real, intent(inout)::   ua(isd:ied,jsd:jed,km)   ! u-wind (m/s) on physics grid
  real, intent(inout)::   va(isd:ied,jsd:jed,km)   ! v-wind (m/s) on physics grid
  real, intent(inout):: omga(isd:ied,jsd:jed,km)   ! vertical press. velocity (pascal/sec)
  real(p_precision), intent(inout)::   peln(is:ie,km+1,js:je)     ! log(pe)
  real(p_precision), intent(out)::    pkz(is:ie,js:je,km)       ! layer-mean pk for converting t to pt
  real, intent(out)::     te(is:ie,js:je,km)

! Mass fluxes
  real, optional, intent(inout):: mfx(is:ie+1,js:je  ,km)   ! X-dir Mass Flux
  real, optional, intent(inout):: mfy(is:ie  ,js:je+1,km)   ! Y-dir Mass Flux

! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!-----------------------------------------------------------------------
  integer :: i,j,k 
      real(p_precision) q_source(is:ie,js:je,nq)    ! numerical tracer source from surface
                                      ! in case fillz is not sufficient
      real(p_precision) te_2d(is:ie,js:je)
      real(p_precision) zsum0(is:ie,js:je)
      real(p_precision) zsum1(is:ie,js:je)
      real(p_precision)   q2(is:ie,km)
      real(p_precision)  dp2(is:ie,km)
      real(p_precision)  pe1(is:ie,km+1)
      real(p_precision)  pe2(is:ie,km+1)
      real(p_precision)  pk1(is:ie,km+1)
      real(p_precision)  pk2(is:ie,km+1)
      real(p_precision)  pn1(is:ie,km+1)
      real(p_precision)  pn2(is:ie,km+1)
      real(p_precision)  pe0(is:ie+1,km+1)
      real(p_precision)  pe3(is:ie+1,km+1)
      real(p_precision) phis(is:ie,km+1)
      real(p_precision)   gz(is:ie)
! for nonhydrostatic option with hybrid_z coordinate
      real(p_precision) ze1(is:ie,km+1), ze2(is:ie,km+1), deng(is:ie,km)
      real dz1(km), ztop, z_rat
      real(p_precision) rcp, rg, ak1, tmp, tpe, cv, rgama, rrg
      real(p_precision) bkh
      real(p_precision) dtmp
      real(p_precision) dlnp
      integer iq, n, kp, k_next
      logical te_map
      real(p_precision)  k1k, kapag

        k1k = akap / (1.-akap)    ! rg/Cv=0.4
      kapag = -akap / grav
         rg = akap * cp
         cv = cp - rg
      rgama = 1. - akap           ! cv/cp
        rcp = 1./ cp
        ak1 = (akap + 1.) / akap
        rrg = - rg/grav

      if ( kord_tm < 0 ) then
           te_map = .false.
          if ( remap_t ) then
! Note: pt at this stage is cp*Theta_v
! Transform virtual pt to virtual Temp
             if ( hydrostatic ) then
!$omp parallel do default(shared)
             do k=1,km
                do j=js,je
                   do i=is,ie
                      pt(i,j,k) = pt(i,j,k) * (pk(i,j,k+1)-pk(i,j,k)) /  &
                                 (rg*(peln(i,k+1,j)-peln(i,k,j)) )
                   enddo
                enddo
             enddo
             else
               if ( ktop>1 ) then
!$omp parallel do default(shared)
                 do k=1,ktop-1
                    do j=js,je
                    do i=is,ie
                       pt(i,j,k) = pt(i,j,k) * (pk(i,j,k+1)-pk(i,j,k)) /  &
                                  (rg*(peln(i,k+1,j)-peln(i,k,j)) )
                    enddo
                    enddo
                 enddo
               endif
!$omp parallel do default(shared)
               do k=ktop,km
                  do j=js,je
                  do i=is,ie
                     pt(i,j,k) = rcp*pt(i,j,k)*exp(k1k*log(kapag*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                  enddo
                  enddo
               enddo
             endif         ! hydro test
          endif            ! remap_t test
      else
           te_map = .true.
           call pkez(km, is, ie, js, je, pe, pk, akap, peln, pkz)

! Compute cp*T + KE
!$omp parallel do default(shared) 
           do k=1,km
              do j=js,je
                 do i=is,ie
                    te(i,j,k) = 0.25*rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                                 v(i,j,k)**2+v(i+1,j,k)**2 -  &
                               (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j))  &
                              +  pt(i,j,k)*pkz(i,j,k)
                 enddo
              enddo
           enddo
     endif

     if ( (.not.hydrostatic) .and. (.not.hybrid_z) ) then
!$omp parallel do default(shared) private(i, j, k)
           do k=1,km
              do j=js,je
                 do i=is,ie
                    delz(i,j,k) = -delz(i,j,k) / delp(i,j,k) ! ="specific volume"/grav
                 enddo
              enddo
           enddo
     endif

!$omp parallel do default(shared) private(z_rat, kp, k_next, bkh, deng, dp2, pe0, pe1, pe2, pe3, pk1, pk2, pn2, phis, q2, ze1, ze2)
  do 1000 j=js,je+1

        do k=1,km+1
           do i=is,ie
              pe1(i,k) = pe(i,k,j)
           enddo
        enddo

        do i=is,ie
           pe2(i,   1) = ptop
           pe2(i,km+1) = pe(i,km+1,j)
        enddo

  if ( j < (je+1) )  then 
! update ps
        do i=is,ie
            ps(i,j) = pe1(i,km+1)
        enddo

!
! Hybrid sigma-P coordinate:
!
        do k=2,ks+1
           do i=is,ie
              pe2(i,k) = ak(k)
           enddo
        enddo
        do k=ks+2,km
           do i=is,ie
              pe2(i,k) = ak(k) + bk(k)*pe(i,km+1,j)
           enddo
        enddo
        do k=1,km
           do i=is,ie
              dp2(i,k) = pe2(i,k+1) - pe2(i,k)
           enddo
        enddo

!------------
! update delp
!------------
      do k=1,km
         do i=is,ie
            delp(i,j,k) = dp2(i,k)
         enddo
      enddo

!----------------
! Map constituents
!----------------
       if( nq /= 0 ) then
!------------------------------------------------------------------
! Do remapping one tracer at a time; seems to be faster on the SGI
! It requires less memory than mapn_ppm
!------------------------------------------------------------------
          do iq=1,nq
             call map1_q2(km, pe1, q(isd,jsd,1,iq),     &
                          km, pe2, q2, dp2,             &
                          is, ie, 0, kord_tr, j, isd, ied, jsd, jed)
            do k=1,km
               do i=is,ie
                  q(i,j,k,iq) = q2(i,k)
               enddo
            enddo
          enddo
       endif

!------------------
! Compute p**cappa
!------------------
   do k=1,km+1
      do i=is,ie
         pk1(i,k) = pk(i,j,k)
      enddo
   enddo

   do i=is,ie
      pn1(i,   :) = peln(i,   :,j)
      pn2(i,   1) = peln(i,   1,j)
      pn2(i,km+1) = peln(i,km+1,j)
      pk2(i,   1) = pk1(i,   1)
      pk2(i,km+1) = pk1(i,km+1)
   enddo

   do k=2,km
      do i=is,ie
!        pk2(i,k) = pe2(i,k) ** akap
         pn2(i,k) = log(pe2(i,k))
         pk2(i,k) = exp(akap*pn2(i,k))
      enddo
   enddo

   if ( te_map ) then
!---------------------
! Compute Total Energy
!---------------------
        do i=is,ie
           phis(i,km+1) = hs(i,j)
        enddo
        do k=km,1,-1
           do i=is,ie
              phis(i,k) = phis(i,k+1) + pt(i,j,k)*(pk1(i,k+1)-pk1(i,k))
           enddo
        enddo
        do k=1,km+1
           do i=is,ie
              phis(i,k) = phis(i,k) * pe1(i,k)
           enddo
        enddo
        do k=1,km
           do i=is,ie
              te(i,j,k) = te(i,j,k)+(phis(i,k+1)-phis(i,k))/(pe1(i,k+1)-pe1(i,k))
           enddo
        enddo
!----------------
! Map Total Energy
!----------------
        select case ( te_method )
        case (1)
        call map1_cubic_te (km,   pe1,  te,       &
                            km,   pe2,  te,       &
                            is, ie, j, is, ie, js, je, 1, kord_tm)
        case default
        call map1_ppm (km,   pe1,  te,       &
                       km,   pe2,  te,       &
                       is, ie, j, is, ie, js, je, 1, kord_tm)
        endselect
   else
     if ( remap_t ) then
!----------------------------------
! Map t using ze1 (hybrid_z) or logp 
!----------------------------------

         call map1_ppm (km,  peln(is,1,j),  pt,    &
                        km,  pn2,           pt,    &
                        is, ie, j, isd, ied, jsd, jed, 2, abs(kord_tm))

     else
!----------------
! Map pt using pk
!----------------
       call map1_ppm (km,  pk1,  pt,           &
                      km,  pk2,  pt,           &
                      is, ie, j, isd, ied, jsd, jed, 1, abs(kord_tm))
     endif
   endif

   if ( .not. hydrostatic ) then
! Remap vertical wind:
        call map1_ppm (km,   pe1,  w,       &
                       km,   pe2,  w,       &
                       is, ie, j, isd, ied, jsd, jed, -2, kord_mt)
     if ( .not. hybrid_z ) then
! Remap delz for hybrid sigma-p coordinate
        call map1_ppm (km,   pe1, delz,    &
                       km,   pe2, delz,    &
                       is, ie, j, is,  ie,  js,  je,  1, abs(kord_tm))
        do k=1,km
           do i=is,ie
              delz(i,j,k) = -delz(i,j,k)*dp2(i,k)
           enddo
        enddo
     endif
   endif

!----------
! Update pk
!----------
   do k=2,km
      do i=is,ie
         pk(i,j,k) = pk2(i,k)
      enddo
   enddo

!----------------
   if ( do_omega ) then
! Start do_omega
! Copy omega field to pe3
      do i=is,ie
         pe3(i,1) = 0.
      enddo
      do k=2,km+1
         do i=is,ie
            pe3(i,k) = omga(i,j,k-1)
         enddo
      enddo
   endif

   do k=1,km+1
      do i=is,ie
          pe0(i,k)   = peln(i,k,j)
         peln(i,k,j) =  pn2(i,k)
      enddo
   enddo

!------------
! Compute pkz
!------------
   if ( hydrostatic ) then
      do k=1,km
         do i=is,ie
            pkz(i,j,k) = (pk2(i,k+1)-pk2(i,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
         enddo
      enddo
   else
      if ( ktop>1 ) then
         do k=1,ktop-1
         do i=is,ie
            pkz(i,j,k) = (pk2(i,k+1)-pk2(i,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
         enddo
         enddo
      endif
      if ( remap_t ) then
! Note: pt at this stage is T_v
         do k=ktop,km
         do i=is,ie
            pkz(i,j,k) = exp(akap*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
         enddo
         enddo
      else
! Note: pt at this stage is cp*Theta_v
         do k=ktop,km
         do i=is,ie
            pkz(i,j,k) = exp( k1k*log(kapag*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)) )
         enddo
         enddo
      endif
   endif

! Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
   if ( do_omega ) then
   do k=1,km
      do i=is,ie
         dp2(i,k) = 0.5*(peln(i,k,j) + peln(i,k+1,j))
      enddo
   enddo
   do i=is,ie
       k_next = 1
       do n=1,km
          kp = k_next
          do k=kp,km
             if( dp2(i,n) <= pe0(i,k+1) .and. dp2(i,n) >= pe0(i,k) ) then
                 omga(i,j,n) = pe3(i,k)  +  (pe3(i,k+1) - pe3(i,k)) *    &
                       (dp2(i,n)-pe0(i,k)) / (pe0(i,k+1)-pe0(i,k) )
                 k_next = k
                 exit
             endif
          enddo
       enddo
   enddo
   endif     ! end do_omega

  endif !(j < je+1)

 if ( .not.hybrid_z ) then
      do i=is,ie+1
         pe0(i,1) = pe(i,1,j)
      enddo
!------
! map u
!------
      do k=2,km+1
         do i=is,ie
            pe0(i,k) = 0.5*(pe(i,k,j-1)+pe1(i,k))
         enddo
      enddo


      do k=1,ks+1
         do i=is,ie+1
            pe3(i,k) = ak(k)
         enddo
      enddo

      do k=ks+2,km+1
         bkh = 0.5*bk(k)
         do i=is,ie
            pe3(i,k) = ak(k) + bkh*(pe(i,km+1,j-1)+pe1(i,km+1))
         enddo
      enddo

      call map1_ppm( km, pe0(is:ie,:),   u,       &
                     km, pe3(is:ie,:),   u,       &
                     is, ie, j, isd, ied, jsd, jed+1, -1, kord_mt)
      if (present(mfy)) then
         call map1_ppm( km, pe0(is:ie,:), mfy,       &
                        km, pe3(is:ie,:), mfy,       &
                        is, ie, j, is, ie, js, je+1, -1, kord_mt)
      endif

   if (j < je+1) then
!------
! map v
!------
       do k=2,km+1
          do i=is,ie+1
             pe0(i ,k) = 0.5*(pe(i-1,k,j)+pe(i,k,j))
          enddo
       enddo
       do k=ks+2,km+1
          bkh = 0.5*bk(k)
          do i=is,ie+1
             pe3(i,k) = ak(k) + bkh*(pe(i-1,km+1,j)+pe(i,km+1,j))
          enddo
       enddo

       call map1_ppm (km, pe0,  v,              &
                      km, pe3,  v, is, ie+1,    &
                      j, isd, ied+1, jsd, jed, -1, kord_mt)
       if (present(mfx)) then
          call map1_ppm (km, pe0, mfx,              &
                         km, pe3, mfx, is, ie+1,    &
                         j, is, ie+1, js, je, -1, kord_mt)
       endif
   endif ! (j < je+1)
 endif    ! end hybrid_z check
     do k=1,km
        do i=is,ie
           ua(i,j,k) = pe2(i,k+1)
        enddo
     enddo

1000  continue


!$omp parallel do default(shared)
  do k=2,km
     do j=js,je
        do i=is,ie
           pe(i,k,j) = ua(i,j,k-1)
        enddo
     enddo
  enddo

  if( do_consv .and. consv > 0. ) then

    if ( te_map ) then
!$omp parallel do default(shared) 
      do j=js,je
          do i=is,ie
             te_2d(i,j) = te(i,j,1)*delp(i,j,1)
          enddo
          do k=2,km
             do i=is,ie
                te_2d(i,j) = te_2d(i,j) + te(i,j,k)*delp(i,j,k)
             enddo
          enddo
      enddo
    else
!$omp parallel do default(shared) private(gz, phis)
      do j=js,je
        if ( remap_t ) then
         if ( hydrostatic ) then
         do i=is,ie
            gz(i) = hs(i,j)
            do k=1,km
               gz(i) = gz(i) + rg*pt(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))
            enddo
         enddo
         do i=is,ie
            te_2d(i,j) = pe(i,km+1,j)*hs(i,j) - pe(i,1,j)*gz(i)
         enddo

         do k=1,km
            do i=is,ie
               te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cp*pt(i,j,k) +   &
                            0.25*rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                             v(i,j,k)**2+v(i+1,j,k)**2 -  &
                           (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j)))
            enddo
         enddo
         else
! non-hydrostatic & remap_t
           do i=is,ie
              phis(i,km+1) = hs(i,j)
              do k=km,1,-1
                 phis(i,k) = phis(i,k+1) - grav*delz(i,j,k)
              enddo
           enddo
           do i=is,ie
              te_2d(i,j) = 0.
           enddo
           do k=1,km
              do i=is,ie
! KE using 3D winds:
                 te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv*pt(i,j,k) +  &
                              0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*rsin2(i,j)*( &
                              u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                             (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j))))
              enddo
           enddo
         endif
        else
         if ( hydrostatic ) then
            do i=is,ie
               gz(i) = hs(i,j)
               do k=1,km
                  gz(i) = gz(i) + pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
               enddo
            enddo

            do i=is,ie
               te_2d(i,j) = pe(i,km+1,j)*hs(i,j) - pe(i,1,j)*gz(i)
            enddo
            do k=1,km
               do i=is,ie
                  te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(pt(i,j,k)*pkz(i,j,k) +   &
                               0.25*rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                                v(i,j,k)**2+v(i+1,j,k)**2 -  &
                            (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j)))
               enddo
            enddo
         else
!-----------------
! Non-hydrostatic:
!-----------------
           do i=is,ie
              phis(i,km+1) = hs(i,j)
              do k=km,1,-1
                 phis(i,k) = phis(i,k+1) - grav*delz(i,j,k)
              enddo
           enddo
           do i=is,ie
              te_2d(i,j) = 0.
           enddo
           do k=1,km
              do i=is,ie
! KE using 3D winds:
                 te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( rgama*pt(i,j,k)*pkz(i,j,k) +  &
                              0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*rsin2(i,j)*( &
                              u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                             (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j))))
              enddo
           enddo
         endif

       endif
      enddo
    endif

!$omp parallel do default(shared)
      do j=js,je
         do i=is,ie
            zsum1(i,j) = pkz(i,j,1)*delp(i,j,1)
         enddo
         do k=2,km
            do i=is,ie
               zsum1(i,j) = zsum1(i,j) + pkz(i,j,k)*delp(i,j,k)
            enddo
         enddo

         do i=is,ie
            zsum0(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1)) + zsum1(i,j)
            te_2d(i,j) = te0_2d(i,j) - te_2d(i,j)
         enddo
      enddo

         tpe = consv*g_sum(te_2d, is, ie, js, je, ng, area, 0)
      E_Flux = tpe / (grav*pdt*4.*pi*radius**2)    ! unit: W/m**2
                                                   ! Note pdt is "phys" time step

      if ( hydrostatic ) then
           dtmp = tpe / (cp*g_sum(zsum0,  is, ie, js, je, ng, area, 0))
      else
           dtmp = tpe / (cv*g_sum(zsum1,  is, ie, js, je, ng, area, 0))
      endif
!-------------------------------------------------------------------------------
! One may use this quick fix to ensure reproducibility at the expense of a lower
! floating precision; this is fine for the TE correction
!-------------------------------------------------------------------------------
      if ( reproduce_sum ) dtmp = real(dtmp, 4) ! convert to 4-byte real
  else
      dtmp   = 0.
      E_Flux = 0.
  endif        ! end consv check

  if ( te_map ) then
!$omp parallel do default(shared) private(gz, tpe, tmp, dlnp)
      do j=js,je
         do i=is,ie
            gz(i) = hs(i,j)
         enddo
         do k=km,1,-1
            do i=is,ie
               tpe = te(i,j,k) - gz(i) - 0.25*rsin2(i,j)*(    &
                     u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                    (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j) )
               dlnp = rg*(peln(i,k+1,j) - peln(i,k,j))
#ifdef CONVERT_T
               tmp = tpe / ((cp - pe(i,k,j)*dlnp/delp(i,j,k))*(1.+r_vir*q(i,j,k,sphum)) )
               pt(i,j,k) =  tmp + dtmp*pkz(i,j,k) / (1.+r_vir*q(i,j,k,sphum))
               gz(i) = gz(i) + dlnp*tmp*(1.+r_vir*q(i,j,k,sphum))
#else
               tmp = tpe / (cp - pe(i,k,j)*dlnp/delp(i,j,k))
               pt(i,j,k) = cp*(tmp/pkz(i,j,k) + dtmp)
               gz(i) = gz(i) + dlnp*tmp
#endif
            enddo
         enddo           ! end k-loop
      enddo
  else
    if ( remap_t ) then
!$omp parallel do default(shared) 
      do k=1,km
         do j=js,je
            do i=is,ie
#ifdef CONVERT_T
               pt(i,j,k) = (pt(i,j,k) + dtmp*pkz(i,j,k))/(1.+r_vir*q(i,j,k,sphum))
#else
               pt(i,j,k) = cp*(pt(i,j,k)/pkz(i,j,k) + dtmp)
#endif
            enddo
         enddo   
      enddo
    else
!$omp parallel do default(shared) 
      do k=1,km
         do j=js,je
            do i=is,ie
#ifdef CONVERT_T
               pt(i,j,k) = (rcp*pt(i,j,k) + dtmp)*pkz(i,j,k)/(1.+r_vir*q(i,j,k,sphum))
#else
               pt(i,j,k) = pt(i,j,k) + cp*dtmp
#endif
            enddo
         enddo   
      enddo
    endif
  endif

 end subroutine Lagrangian_to_Eulerian


  subroutine pkez(km, ifirst, ilast, jfirst, jlast, &
                  pe, pk, akap, peln, pkz)

  IMPLICIT NONE

! !INPUT PARAMETERS:
   integer, intent(in):: km
   integer, intent(in):: ifirst, ilast        ! Latitude strip
   integer, intent(in):: jfirst, jlast        ! Latitude strip
   real, intent(in):: akap
   real(p_precision), intent(in):: pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1)
   real(p_precision), intent(in):: pk(ifirst:ilast,jfirst:jlast,km+1)
! !OUTPUT
   real(p_precision), intent(out):: pkz(ifirst:ilast,jfirst:jlast,km)
   real(p_precision), intent(inout):: peln(ifirst:ilast, km+1, jfirst:jlast)   ! log (pe)
! Local
   real(p_precision) pk2(ifirst:ilast, km+1)
   real(p_precision) pek
   real(p_precision) lnp
   real(p_precision) ak1
   integer i, j, k

   ak1 = (akap + 1.) / akap

!$omp parallel do default(shared) private(lnp, pek, pk2)
   do j=jfirst, jlast
        pek = pk(ifirst,j,1)
        do i=ifirst, ilast
           pk2(i,1) = pek
        enddo

        do k=2,km+1
           do i=ifirst, ilast
!             peln(i,k,j) =  log(pe(i,k,j))
              pk2(i,k) =  pk(i,j,k)
           enddo
        enddo

!---- GFDL modification
       if( ptop < ptop_min ) then
           do i=ifirst, ilast
               peln(i,1,j) = peln(i,2,j) - ak1
           enddo
       else
           lnp = log( ptop )
           do i=ifirst, ilast
              peln(i,1,j) = lnp
           enddo
       endif
!---- GFDL modification

       do k=1,km
          do i=ifirst, ilast
             pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k) )  /  &
                          (akap*(peln(i,k+1,j) - peln(i,k,j)) )
          enddo
       enddo
    enddo

 end subroutine pkez

 subroutine map1_ppm( km,   pe1,    q1,                 &
                      kn,   pe2,    q2,   i1, i2,       &
                      j,    ibeg, iend, jbeg, jend, iv,  kord)
 integer, intent(in) :: i1                ! Starting longitude
 integer, intent(in) :: i2                ! Finishing longitude
 integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == ???
                                          !       2 == remap temp with cs scheme
 integer, intent(in) :: kord              ! Method order
 integer, intent(in) :: j                 ! Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                ! Original vertical dimension
 integer, intent(in) :: kn                ! Target vertical dimension
 real(p_precision), intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the original vertical coordinate
 real(p_precision), intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the new vertical coordinate
 real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) ! Field output

! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
   real(p_precision)    dp1(i1:i2,km)
   real(p_precision)   q4(4,i1:i2,km)
   real(p_precision)    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map1_ppm

