subroutine tracer_2d_1L(q, dp0, mfx, mfy, cx, cy, npx, npy, npz, nq, hord,  &
                        q_split, k, q3, dt, id_divg, is, ie, js, je, isd, ied, jsd, jed, &
                        dxa, dya, sina_u, sina_v, dx, dy, area, rarea)

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

      integer, intent(IN) :: is, ie, js, je, isd, ied, jsd, jed
      real   , intent(IN) ::  dx(is:ie+1,jsd:jed), dy(isd:ied,js: je+1)
      real   , intent(IN) ::  dxa(is:ie+1,jsd:jed), dya(isd:ied,js: je+1)
      real   , intent(IN) ::  sina_u(is:ie+1,jsd:jed), sina_v(isd:ied,js: je+1)
      real   , intent(IN) ::  area(is:ie,jsd:jed), rarea(is:ie,js:je)

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
         call mpp_update_domains( q, complete= .true. )
              !call timing_off('COMM_TRAC')
         !call timing_off('COMM_TOTAL')

         do iq=1,nq

            call fv_tp_2d( q(isd:,jsd:,iq), cx2, cy2, npx, npy, hord, fx, fy, &
                           xfx, yfx, area, ra_x, ra_y, mfx=mfx2, mfy=mfy2 )

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

end subroutine tracer_2d_1L

subroutine mp_reduce_max(cmax)

implicit none

real :: cmax

cmax = cmax**2 + cmax/2

end subroutine


subroutine fv_tp_2d( q, cx2, cy2, npx, npy, hord, fx, fy, xfx, yfx, area, ra_x, ra_y, mfx, mfy )

implicit none

integer :: npx, npy, hord
real :: q(:,:), cx2(:,:), cy2(:,:), fx(:,:), fy(:,:), xfx(:,:), yfx(:,:), area(:,:), ra_x(:,:), ra_y(:,:), mfx(:,:), mfy(:,:)


q = q**2 + q/2 + q/cx2 + q/cy2 + q/ra_x + q/ra_y

fx = q/3
fy = q/4
cx2 = cx2**2
cy2 = cy2**2
fx = fx**2
fy = fy**2
xfx = xfx**2
yfx = yfx**2
ra_x = ra_x**2
ra_y = ra_y**2
mfx = mfx**2
mfy = mfy**2

end subroutine


subroutine mpp_update_domains( q, complete )

implicit none

real :: q(:,:)
logical :: complete

if (complete == .true.) then
q = q**2 + q/2
endif

endsubroutine
