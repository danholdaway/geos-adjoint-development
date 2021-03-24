subroutine tracer_2d_1L_dh( q, dp0, mfx, mfy, cx, cy, npx, npy, npz, nq, hord,  &
                         q_split, k, q3, dt, id_divg, is, ie, js, je, isd, ied, jsd, jed, ng, &
                         dxa, dya, sina_u, sina_v, dx, dy, area, rarea, domain)

implicit none

      integer, intent(IN) :: is, ie, js, je, isd, ied, jsd, jed, ng

      integer, intent(IN) :: npx, npy, npz
      integer, intent(IN) :: k
      integer, intent(IN) :: nq    ! number of tracers to be advected
      integer, intent(IN) :: hord
      integer, intent(IN) :: q_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt
      real   , intent(INOUT) ::   q(isd:ied,jsd:jed,nq)       ! 2D Tracers
      real   , intent(INOUT) ::  q3(isd:ied,jsd:jed,npz,nq)   ! Tracers
      real   , intent(INOUT) :: dp0(is:ie,js:je)        ! DELP before dyn_core
      real   , intent(IN)    :: mfx(is:ie+1,js:je)    ! Mass Flux X-Dir
      real   , intent(IN)    :: mfy(is:ie  ,js:je+1)    ! Mass Flux Y-Dir
      real   , intent(IN)    ::  cx(is:ie+1,jsd:jed)  ! Courant Number X-Dir
      real   , intent(IN)    ::  cy(isd:ied,js :je +1)  ! Courant Number Y-Dir

      !Unrequired inputs
      real   , intent(IN) :: dx(is:ie+1,jsd:jed), dy(isd:ied,js: je+1)
      real   , intent(IN) :: dxa(is:ie+1,jsd:jed), dya(isd:ied,js: je+1)
      real   , intent(IN) :: sina_u(is:ie+1,jsd:jed), sina_v(isd:ied,js: je+1)
      real   , intent(IN) :: area(is:ie,jsd:jed), rarea(is:ie,js:je)
      real   , intent(IN) :: domain

! Local Arrays
      real   :: mfx2(is:ie+1,js:je)
      real   :: mfy2(is:ie  ,js:je+1) 
      real   ::  cx2(is:ie+1,jsd:jed)
      real   ::  cy2(isd:ied,js :je +1)

      real   :: dp1(is:ie,js:je)
      real   :: dp2(is:ie,js:je)
      real   :: fx(is:ie+1,js:je )
      real   :: fy(is:ie , js:je+1)
      real   :: ra_x(is:ie,jsd:jed)
      real   :: ra_y(isd:ied,js:je)
      real   :: xfx(is:ie+1,jsd:jed)
      real   :: yfx(isd:ied,js: je+1)
      real   :: cmax
      real   :: frac, rdt
      integer :: nsplt
      integer :: i,j,it,iq

      do j=jsd,jed
         do i=is,ie+1
            if (cx(i,j) > 0.) then
                xfx(i,j) = cx(i,j)*dxa(i-1,j)*dy(i,j)*sina_u(i,j)
            else
                xfx(i,j) = cx(i,j)*dxa(i,j)*dy(i,j)*sina_u(i,j)
            endif
         enddo
      enddo

      do j=js,je+1
         do i=isd,ied
            if (cy(i,j) > 0.) then
                yfx(i,j) = cy(i,j)*dya(i,j-1)*dx(i,j)*sina_v(i,j)
            else
                yfx(i,j) = cy(i,j)*dya(i,j)*dx(i,j)*sina_v(i,j)
            endif
         enddo
      enddo

      if ( q_split==0 ) then
! Determine nsplt for tracer advection
         cmax = 0.
         do j=js,je
            do i=is,ie
               cmax = max(abs(cx(i,j))+(1.-sina_u(i,j)),     &
                          abs(cy(i,j))+(1.-sina_v(i,j)), cmax)
            enddo
         enddo
         call mp_reduce_max(cmax)
         nsplt = int(1.01 + cmax)
      else
         nsplt = q_split
      endif

      frac  = 1. / real(nsplt)
          do j=jsd,jed
             do i=is,ie+1
                cx2(i,j) =  cx(i,j) * frac
                xfx(i,j) = xfx(i,j) * frac
             enddo
          enddo
          do j=js,je
             do i=is,ie+1
                mfx2(i,j) = mfx(i,j) * frac
             enddo
          enddo

          do j=js,je+1
             do i=isd,ied
                cy2(i,j) =  cy(i,j) * frac
               yfx(i,j) = yfx(i,j) * frac
             enddo
          enddo

          do j=js,je+1
             do i=is,ie
                mfy2(i,j) = mfy(i,j) * frac
             enddo
          enddo

      do j=jsd,jed
         do i=is,ie
            ra_x(i,j) = area(i,j) + xfx(i,j) - xfx(i+1,j)
         enddo
      enddo
      do j=js,je
         do i=isd,ied
            ra_y(i,j) = area(i,j) + yfx(i,j) - yfx(i,j+1)
         enddo
      enddo

      do j=js,je
         do i=is,ie
            dp1(i,j) = dp0(i,j)
         enddo
      enddo

      do it=1,nsplt

         do j=js,je
            do i=is,ie
               dp2(i,j) = dp1(i,j) + (mfx2(i,j) - mfx2(i+1,j) +  &
                          mfy2(i,j) - mfy2(i,j+1)) * rarea(i,j)
            enddo
         enddo

         !call timing_on('COMM_TOTAL')
              !call timing_on('COMM_TRAC')
         call mpp_update_domains( q, domain, complete= .true. )
              !call timing_off('COMM_TRAC')
         !call timing_off('COMM_TOTAL')

         do iq=1,nq

            call fv_tp_2d_trc_dh( q(isd:,jsd:,iq), cx2, cy2, npx, npy, hord, fx, fy, &
                                  xfx, yfx, area, ra_x, ra_y, mfx2, mfy2,       &
                                  is, ie, js, je, isd, ied, jsd, jed, ng  )

            if( it==nsplt ) then
            do j=js,je
               do i=is,ie
                  q3(i,j,k,iq) = (q(i,j,iq)*dp1(i,j) + (fx(i,j)-fx(i+1,j) + &
                                  fy(i,j)-fy(i,j+1))*rarea(i,j)) / dp2(i,j)
               enddo
            enddo
            else
            do j=js,je
               do i=is,ie
                  q(i,j,iq) = (q(i,j,iq)*dp1(i,j) + (fx(i,j)-fx(i+1,j) + &
                              fy(i,j)-fy(i,j+1))*rarea(i,j)) / dp2(i,j)
               enddo
            enddo
           endif
         enddo

         if ( it/=nsplt ) then
              do j=js,je
                 do i=is,ie
                    dp1(i,j) = dp2(i,j)
                 enddo
              enddo
         endif
     enddo  ! nsplt

     if ( id_divg > 0 ) then
         rdt = 1./(frac*dt)
         do j=js,je
            do i=is,ie
               dp0(i,j) = (xfx(i+1,j)-xfx(i,j) + yfx(i,j+1)-yfx(i,j))*rarea(i,j)*rdt
            enddo
         enddo
     endif

end subroutine tracer_2d_1L_dh

subroutine mp_reduce_max(cmax)

implicit none

real :: cmax

if (cmax > 2) then
cmax = cmax**2 + cmax/2
endif

end subroutine

subroutine mpp_update_domains( q, domain, complete )

implicit none

real :: q(:,:), domain
logical :: complete

if (complete == .true. .and. domain > 1) then
q = q**2 + q/2
endif

endsubroutine



 subroutine fv_tp_2d_trc_dh( q, crx, cry, npx, npy, hord, fx, fy,  &
                             xfx, yfx, area, ra_x, ra_y, mfx, mfy, &
                             is, ie, js, je, isd, ied, jsd, jed, ng    )

   IMPLICIT NONE

   integer, intent(IN) :: is, ie, js, je, isd, ied, jsd, jed, ng

   integer, intent(in):: npx, npy
   integer, intent(in)::hord

   real, intent(in)::  crx(is:ie+1,jsd:jed)  !
   real, intent(in)::  xfx(is:ie+1,jsd:jed)  !
   real, intent(in)::  cry(isd:ied,js:je+1 )  !
   real, intent(in)::  yfx(isd:ied,js:je+1 )  !
   real, intent(in):: area(isd:ied,jsd:jed)
   real, intent(in):: ra_x(is:ie,jsd:jed)
   real, intent(in):: ra_y(isd:ied,js:je)
   real, intent(inout):: q(isd:ied,jsd:jed)  ! transported scalar
   real, intent(out)::fx(is:ie+1 ,js:je)    ! Flux in x ( E )
   real, intent(out)::fy(is:ie,   js:je+1 )    ! Flux in y ( N )

! optional Arguments:
   real, OPTIONAL, intent(in):: mfx(is:ie+1,js:je  )  ! Mass Flux X-Dir
   real, OPTIONAL, intent(in):: mfy(is:ie  ,js:je+1)  ! Mass Flux Y-Dir
! Local:
   integer ord, ord_in
   real q_i(isd:ied,js:je)
   real q_j(is:ie,jsd:jed)
   real   fx2(is:ie+1,jsd:jed)
   real   fy2(isd:ied,js:je+1)
   real   fyy(isd:ied,js:je+1)
   real   fx1(is:ie+1)
   real   ppm_limiter
   integer i, j

   if ( hord < 0 ) then
      ord_in = 2 !More dissipation
      ord    = abs(hord)
   else
      ord_in = hord
      ord    = hord
   endif

   ppm_limiter = 2.0

   !Part 1
   call copy_corners(q,2)
   call ytp_dh(fy2, q, cry, ord_in, isd, ied, js, je, npx, npy, ppm_limiter, is, ie, js, je, isd, ied, jsd, jed, ng)

   do j=js,je+1
      do i=isd,ied
         fyy(i,j) = yfx(i,j) * fy2(i,j)
      enddo
   enddo
   do j=js,je
      do i=isd,ied
         q_i(i,j) = (q(i,j)*area(i,j) + fyy(i,j)-fyy(i,j+1))/ra_y(i,j)
      enddo
  enddo
  call xtp_dh(fx , q_i, crx(is,js), ord   , is, ie, js , je , npx, npy, ppm_limiter, is, ie, js, je, isd, ied, jsd, jed, ng)

  !Part 2
  call copy_corners(q,1)
  call xtp_dh(fx2, q  , crx       , ord_in, is, ie, jsd, jed, npx, npy, ppm_limiter, is, ie, js, je, isd, ied, jsd, jed, ng)

  do j=jsd,jed
     do i=is,ie+1
        fx1(i) =  xfx(i,j) * fx2(i,j)
     enddo
     do i=is,ie
        q_j(i,j) = (q(i,j)*area(i,j) + fx1(i)-fx1(i+1))/ra_x(i,j)
     enddo
  enddo

  call ytp_dh(fy, q_j, cry, ord, is, ie, js, je, npx, npy, ppm_limiter, is, ie, js, je, isd, ied, jsd, jed, ng)


! Flux averaging:
!----------------
   if ( present(mfx) .and. present(mfy) ) then
!---------------------------------
! For transport of pt and tracers
!---------------------------------
      do j=js,je
         do i=is,ie+1
            fx(i,j) = 0.5*(fx(i,j) + fx2(i,j)) * mfx(i,j)
         enddo
      enddo
      do j=js,je+1
         do i=is,ie
            fy(i,j) = 0.5*(fy(i,j) + fy2(i,j)) * mfy(i,j)
         enddo
      enddo
   else
!---------------------------------
! For transport of delp, vorticity
!---------------------------------
      do j=js,je
         do i=is,ie+1
            fx(i,j) = 0.5*(fx(i,j) + fx2(i,j)) * xfx(i,j)
         enddo
      enddo
      do j=js,je+1
         do i=is,ie
            fy(i,j) = 0.5*(fy(i,j) + fy2(i,j)) * yfx(i,j)
         enddo
      enddo
   endif

 end subroutine fv_tp_2d_trc_dh


 subroutine copy_corners(q,indy)

 IMPLICIT NONE

 integer :: indy

 !Prognostic
 real, intent(inout):: q(:,:)

 !Locals
 integer  i,j

if (indy == 2) then
 do j=1,size(q,1)
    do i=1,size(q,2)
       q(i,j) = q(j,1-i)
    enddo
 enddo
else
 do j=1,size(q,1)
    do i=1,size(q,2)
       q(i,j) = q(j,1-2)
    enddo
 enddo
endif

 end subroutine copy_corners


 subroutine xtp_dh(fx,  q,  c, iord, ifirst, ilast, jfirst, jlast, npx, npy, ppm_limiter, is, ie, js, je, isd, ied, jsd, jed, ng )

   IMPLICIT NONE

   integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed, ng

   integer, intent(IN):: ifirst, ilast   !  X-Dir strip
   integer, intent(IN):: jfirst, jlast   !  Y-Dir strip
   integer, intent(IN):: npx, npy
   integer, intent(IN):: iord
   real   , intent(in):: c(is :ie+1, jfirst:jlast)      ! Courant numbers
   real   , intent(in):: q(isd:ied,  jfirst:jlast)
   real   , intent(in):: ppm_limiter
   real,    intent(out):: fx(ifirst:ilast+1,jfirst:jlast)           
! Local:
   real   dm(is-2:ie+2)
   real   x0, x1
   integer i, j


      do j=jfirst,jlast
         do i=ifirst,ilast+1
           if ( c(i,j)>0. ) then
                fx(i,j) = q(i-1,j)
           else
                fx(i,j) = q(i,j)
           endif
         enddo
      enddo


 end subroutine xtp_dh


 subroutine ytp_dh(fy, q, c, jord, ifirst, ilast, jfirst, jlast, npx, npy, ppm_limiter, is, ie, js, je, isd, ied, jsd, jed, ng)

 IMPLICIT NONE

 integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed, ng

 integer, intent(in) :: npx, npy
 integer, INTENT(IN) :: ifirst, ilast  !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast  !  Y-Dir strip
 integer, intent(in):: jord
 real, intent(in)::   q(ifirst:ilast,jfirst-ng:jlast+ng) 
 real, intent(in)::   c(isd:ied,js:je+1 )  ! Courant number
 real, intent(out):: fy(ifirst:ilast,jfirst:jlast+1)     !  Flux
 real, intent(in)::  ppm_limiter
! !LOCAL VARIABLES:
 real  dm(ifirst:ilast,jfirst-2:jlast+2)
 real  x0, x1
 integer i, j


      do j=jfirst,jlast+1
         do i=ifirst,ilast
            if ( c(i,j)>0. ) then
                 fy(i,j) = q(i,j-1)
            else
                 fy(i,j) = q(i,j)
            endif
         enddo
      enddo

 end subroutine ytp_dh
