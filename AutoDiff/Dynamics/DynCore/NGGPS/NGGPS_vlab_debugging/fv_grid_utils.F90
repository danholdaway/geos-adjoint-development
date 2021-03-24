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
 module fv_grid_utils_mod
 



! -*-f90-*-*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!parallel machine types




!most compilers support Cray pointers
!if you find a compiler that doesnt, #undef this inside a suitable #ifdef


!values of kind: double and long are 8-byte, float and int are 4-byte
!pointer_kind is used for storing addresses as integers

!these might be different on non-SGICRAY, I believe
! Warning: these numbers may not map to byte sizes for all compilers







!DEC$ MESSAGE:'Using 8-byte addressing'



















!DEC$ MESSAGE:'Using PURE'




!DEC$ MESSAGE:'Converting pointers to allocatable components'







! 

 use mpp_mod,         only: FATAL, mpp_error, WARNING
 use mpp_domains_mod, only: mpp_update_domains, DGRID_NE, mpp_global_sum
 use mpp_domains_mod, only: domain2d, BITWISE_EFP_SUM

 use fv_arrays_mod,   only: fv_atmos_type, fv_grid_type, fv_grid_bounds_type, &
                            R_GRID
 use fv_mp_mod,       only: ng, is_master
 use fv_mp_mod,       only: mp_reduce_sum
 use fv_mp_mod,       only: fill_corners, XDir, YDir
 use fv_timing_mod,   only: timing_on, timing_off

 implicit none
 private
 logical:: symm_grid

! Higher precision (kind=16) for grid geometrical factors:
 integer, parameter:: f_p = selected_real_kind(20)

 real, parameter::  big_number=1.e8
 real, parameter:: tiny_number=1.e-8
 real, parameter:: ptop_min=1.e-8

 public f_p 
 public ptop_min, big_number !CLEANUP: OK to keep since they are constants?
 public grid_utils_init, grid_utils_end, great_circle_dist, cubed_to_latlon
 public c2l_ord2, g_sum


 contains

   subroutine grid_utils_init(Atm, npx, npy, npz, non_ortho, grid_type, c2l_order)
! Initialize 2D memory and geometrical factors
      type(fv_atmos_type), intent(inout), target :: Atm
      logical, intent(in):: non_ortho
      integer, intent(in):: npx, npy, npz
      integer, intent(in):: grid_type, c2l_order

  end subroutine grid_utils_init

 
  subroutine grid_utils_end
 
  end subroutine grid_utils_end






 real function great_circle_dist( q1, q2, radius )
      real(kind=R_GRID), intent(IN)           :: q1(2), q2(2)
      real(kind=R_GRID), intent(IN), optional :: radius
 
      real (f_p):: p1(2), p2(2)
      real (f_p):: beta
      integer n

      do n=1,2
         p1(n) = q1(n)
         p2(n) = q2(n)
      enddo

      beta = asin( sqrt( sin((p1(2)-p2(2))/2.)**2 + cos(p1(2))*cos(p2(2))*   &
                         sin((p1(1)-p2(1))/2.)**2 ) ) * 2.

      if ( present(radius) ) then
           great_circle_dist = radius * beta
      else
           great_circle_dist = beta   ! Returns the angle
      endif

  end function great_circle_dist


 subroutine cubed_to_latlon(u, v, ua, va, gridstruct, npx, npy, km, mode, grid_type, domain, nested, c2l_ord, bd)
 type(fv_grid_bounds_type), intent(IN) :: bd 
 integer, intent(in) :: km, npx, npy, grid_type, c2l_ord
 integer, intent(in) :: mode   ! update if present
 type(fv_grid_type), intent(IN) :: gridstruct
 real, intent(inout):: u(bd%isd:bd%ied,bd%jsd:bd%jed+1,km)
 real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,km)
 real, intent(out):: ua(bd%isd:bd%ied, bd%jsd:bd%jed,km)
 real, intent(out):: va(bd%isd:bd%ied, bd%jsd:bd%jed,km)
 type(domain2d), intent(INOUT) :: domain
 logical, intent(IN) :: nested

 if ( c2l_ord == 2 ) then
      call c2l_ord2(u, v, ua, va, gridstruct, km, grid_type, bd, .false.)
 else
      call c2l_ord4(u, v, ua, va, gridstruct, npx, npy, km, grid_type, domain, nested, mode, bd)
 endif

 end subroutine cubed_to_latlon


 subroutine c2l_ord4(u, v, ua, va, gridstruct, npx, npy, km, grid_type, domain, nested, mode, bd)

 type(fv_grid_bounds_type), intent(IN) :: bd
  integer, intent(in) :: km, npx, npy, grid_type
  integer, intent(in):: mode   ! update if present
 type(fv_grid_type), intent(IN), target :: gridstruct
  real, intent(inout):: u(bd%isd:bd%ied,bd%jsd:bd%jed+1,km)
  real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,km)
  real, intent(out)::  ua(bd%isd:bd%ied, bd%jsd:bd%jed,km)
  real, intent(out)::  va(bd%isd:bd%ied, bd%jsd:bd%jed,km)
  type(domain2d), intent(INOUT) :: domain
  logical, intent(IN) :: nested
! Local
! 4-pt Lagrange interpolation
  real :: a1 =  0.5625
  real :: a2 = -0.0625
  real :: c1 =  1.125
  real :: c2 = -0.125
  real utmp(bd%is:bd%ie,  bd%js:bd%je+1)
  real vtmp(bd%is:bd%ie+1,bd%js:bd%je)
  real wu(bd%is:bd%ie,  bd%js:bd%je+1)
  real wv(bd%is:bd%ie+1,bd%js:bd%je)
  integer i, j, k

  integer :: is,  ie,  js,  je


  is  = bd%is
  ie  = bd%ie
  js  = bd%js
  je  = bd%je

  if ( mode > 0 ) then
                                   call timing_on('COMM_TOTAL')
       call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
                                  call timing_off('COMM_TOTAL')
  endif

!$OMP parallel do default(none) shared(is,ie,js,je,km,npx,npy,grid_type,nested,c2,c1, &
!$OMP                                  u,v,gridstruct,ua,va,a1,a2)         &
!$OMP                          private(utmp, vtmp, wu, wv)
 do k=1,km
   if ( grid_type < 4 ) then
    if (nested) then
     do j=max(1,js),min(npy-1,je)
        do i=max(1,is),min(npx-1,ie)
           utmp(i,j) = c2*(u(i,j-1,k)+u(i,j+2,k)) + c1*(u(i,j,k)+u(i,j+1,k))
           vtmp(i,j) = c2*(v(i-1,j,k)+v(i+2,j,k)) + c1*(v(i,j,k)+v(i+1,j,k))
        enddo
     enddo
   else
     do j=max(2,js),min(npy-2,je)
        do i=max(2,is),min(npx-2,ie)
           utmp(i,j) = c2*(u(i,j-1,k)+u(i,j+2,k)) + c1*(u(i,j,k)+u(i,j+1,k))
           vtmp(i,j) = c2*(v(i-1,j,k)+v(i+2,j,k)) + c1*(v(i,j,k)+v(i+1,j,k))
        enddo
     enddo

    if ( js==1  ) then
         do i=is,ie+1
            wv(i,1) = v(i,1,k)*gridstruct%dy(i,1)
         enddo
         do i=is,ie
            vtmp(i,1) = 2.*(wv(i,1) + wv(i+1,1)) / (gridstruct%dy(i,1)+gridstruct%dy(i+1,1))
            utmp(i,1) = 2.*(u(i,1,k)*gridstruct%dx(i,1) + u(i,2,k)*gridstruct%dx(i,2))   &
                         / (         gridstruct%dx(i,1) +          gridstruct%dx(i,2))
!!!         vtmp(i,1) = (wv(i,1) + wv(i+1,1)) * gridstruct%rdya(i,1)
!!!         utmp(i,1) = (u(i,1,k)*gridstruct%dx(i,1) + u(i,2,k)*gridstruct%dx(i,2)) * gridstruct%rdxa(i,1)
         enddo
    endif

    if ( (je+1)==npy   ) then
         j = npy-1
         do i=is,ie+1
            wv(i,j) = v(i,j,k)*gridstruct%dy(i,j)
         enddo
         do i=is,ie
            vtmp(i,j) = 2.*(wv(i,j) + wv(i+1,j)) / (gridstruct%dy(i,j)+gridstruct%dy(i+1,j))
            utmp(i,j) = 2.*(u(i,j,k)*gridstruct%dx(i,j) + u(i,j+1,k)*gridstruct%dx(i,j+1))   &
                         / (         gridstruct%dx(i,j) +            gridstruct%dx(i,j+1))
!!!         vtmp(i,j) = (wv(i,j) + wv(i+1,j)) * gridstruct%rdya(i,j)
!!!         utmp(i,j) = (u(i,j,k)*gridstruct%dx(i,j) + u(i,j+1,k)*gridstruct%dx(i,j+1)) * gridstruct%rdxa(i,j)
         enddo
    endif

    if ( is==1 ) then
      i = 1
      do j=js,je
         wv(1,j) = v(1,j,k)*gridstruct%dy(1,j)
         wv(2,j) = v(2,j,k)*gridstruct%dy(2,j)
      enddo
      do j=js,je+1
         wu(i,j) = u(i,j,k)*gridstruct%dx(i,j)
      enddo
      do j=js,je
         utmp(i,j) = 2.*(wu(i,j) + wu(i,j+1))/(gridstruct%dx(i,j)+gridstruct%dx(i,j+1))
         vtmp(i,j) = 2.*(wv(1,j) + wv(2,j  ))/(gridstruct%dy(1,j)+gridstruct%dy(2,j))
!!!      utmp(i,j) = (wu(i,j) + wu(i,  j+1)) * gridstruct%rdxa(i,j)
!!!      vtmp(i,j) = (wv(i,j) + wv(i+1,j  )) * gridstruct%rdya(i,j)
      enddo
    endif

    if ( (ie+1)==npx) then
      i = npx-1
      do j=js,je
         wv(i,  j) = v(i,  j,k)*gridstruct%dy(i,  j)
         wv(i+1,j) = v(i+1,j,k)*gridstruct%dy(i+1,j)
      enddo
      do j=js,je+1
         wu(i,j) = u(i,j,k)*gridstruct%dx(i,j)
      enddo
      do j=js,je
         utmp(i,j) = 2.*(wu(i,j) + wu(i,  j+1))/(gridstruct%dx(i,j)+gridstruct%dx(i,j+1))
         vtmp(i,j) = 2.*(wv(i,j) + wv(i+1,j  ))/(gridstruct%dy(i,j)+gridstruct%dy(i+1,j))
!!!      utmp(i,j) = (wu(i,j) + wu(i,  j+1)) * gridstruct%rdxa(i,j)
!!!      vtmp(i,j) = (wv(i,j) + wv(i+1,j  )) * gridstruct%rdya(i,j)
      enddo
    endif

 endif !nested

!Transform local a-grid winds into latitude-longitude coordinates
     do j=js,je
        do i=is,ie
           ua(i,j,k) = gridstruct%a11(i,j)*utmp(i,j) + gridstruct%a12(i,j)*vtmp(i,j)
           va(i,j,k) = gridstruct%a21(i,j)*utmp(i,j) + gridstruct%a22(i,j)*vtmp(i,j)
        enddo
     enddo
   else
! Simple Cartesian Geometry:
     do j=js,je
        do i=is,ie
           ua(i,j,k) = a2*(u(i,j-1,k)+u(i,j+2,k)) + a1*(u(i,j,k)+u(i,j+1,k))
           va(i,j,k) = a2*(v(i-1,j,k)+v(i+2,j,k)) + a1*(v(i,j,k)+v(i+1,j,k))
        enddo
     enddo
   endif
 enddo
 end subroutine c2l_ord4

 subroutine c2l_ord2(u, v, ua, va, gridstruct, km, grid_type, bd, do_halo)
 type(fv_grid_bounds_type), intent(IN) :: bd
  integer, intent(in) :: km, grid_type
  real, intent(in) ::  u(bd%isd:bd%ied,bd%jsd:bd%jed+1,km)
  real, intent(in) ::  v(bd%isd:bd%ied+1,bd%jsd:bd%jed,km)
 type(fv_grid_type), intent(IN), target :: gridstruct
 logical, intent(in) :: do_halo
!
  real, intent(out):: ua(bd%isd:bd%ied, bd%jsd:bd%jed,km)
  real, intent(out):: va(bd%isd:bd%ied, bd%jsd:bd%jed,km)
!--------------------------------------------------------------
! Local
  real wu(bd%is-1:bd%ie+1,  bd%js-1:bd%je+2)
  real wv(bd%is-1:bd%ie+2,  bd%js-1:bd%je+1)
  real u1(bd%is-1:bd%ie+1), v1(bd%is-1:bd%ie+1)
  integer i, j, k
  integer :: is,  ie,  js,  je

  real, dimension(:,:), pointer :: a11, a12, a21, a22
  real, dimension(:,:), pointer :: dx, dy, rdxa, rdya

  a11 => gridstruct%a11
  a12 => gridstruct%a12
  a21 => gridstruct%a21
  a22 => gridstruct%a22

  dx   => gridstruct%dx
  dy   => gridstruct%dy
  rdxa => gridstruct%rdxa
  rdya => gridstruct%rdya

  if (do_halo) then
     is  = bd%is-1
     ie  = bd%ie+1
     js  = bd%js-1
     je  = bd%je+1
  else
     is  = bd%is
     ie  = bd%ie
     js  = bd%js
     je  = bd%je
  endif

!$OMP parallel do default(none) shared(is,ie,js,je,km,grid_type,u,dx,v,dy,ua,va,a11,a12,a21,a22) &
!$OMP                          private(u1, v1, wu, wv)
  do k=1,km
     if ( grid_type < 4 ) then
       do j=js,je+1
          do i=is,ie
             wu(i,j) = u(i,j,k)*dx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             wv(i,j) = v(i,j,k)*dy(i,j)
          enddo
       enddo

       do j=js,je
          do i=is,ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
             u1(i) = 2.*(wu(i,j) + wu(i,j+1)) / (dx(i,j)+dx(i,j+1))
             v1(i) = 2.*(wv(i,j) + wv(i+1,j)) / (dy(i,j)+dy(i+1,j))
!!!          u1(i) = (wu(i,j) + wu(i,j+1)) * rdxa(i,j)
!!!          v1(i) = (wv(i,j) + wv(i+1,j)) * rdya(i,j)
! Cubed (cell center co-variant winds) to lat-lon:
             ua(i,j,k) = a11(i,j)*u1(i) + a12(i,j)*v1(i)
             va(i,j,k) = a21(i,j)*u1(i) + a22(i,j)*v1(i)
          enddo
       enddo
     else
! 2nd order:
       do j=js,je
          do i=is,ie
             ua(i,j,k) = 0.5*(u(i,j,k)+u(i,  j+1,k))
             va(i,j,k) = 0.5*(v(i,j,k)+v(i+1,j,  k))
          enddo
       enddo
     endif
  enddo

 end subroutine c2l_ord2

! real function g_sum(domain, p, ifirst, ilast, jfirst, jlast, ngc, area, mode, reproduce)
!! Fast version of globalsum
!      integer, intent(IN) :: ifirst, ilast
!      integer, intent(IN) :: jfirst, jlast, ngc
!      integer, intent(IN) :: mode  ! if ==1 divided by area
!      logical, intent(in), optional :: reproduce
!      real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
!      real(kind=R_GRID), intent(IN) :: area(ifirst-ngc:ilast+ngc,jfirst-ngc:jlast+ngc)
!      type(domain2d), intent(IN) :: domain
!      integer :: i,j
!      real gsum
!      logical, SAVE :: g_sum_initialized = .false.
!      real(kind=R_GRID), SAVE :: global_area
!      real :: tmp(ifirst:ilast,jfirst:jlast) 
!        
!      if ( .not. g_sum_initialized ) then
!         global_area = mpp_global_sum(domain, area, flags=BITWISE_EFP_SUM)
!         if ( is_master() ) write(*,*) 'Global Area=',global_area
!         g_sum_initialized = .true.
!      end if
! 
!!-------------------------
!! FMS global sum algorithm:
!!-------------------------
!      if ( present(reproduce) ) then
!         if (reproduce) then
!            gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast), &
!                                  flags=BITWISE_EFP_SUM)
!         else
!            gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast))
!         endif
!      else
!!-------------------------
!! Quick local sum algorithm
!!-------------------------
!         gsum = 0.
!         do j=jfirst,jlast
!            do i=ifirst,ilast
!               gsum = gsum + p(i,j)*area(i,j)
!            enddo
!         enddo
!         call mp_reduce_sum(gsum)
!      endif
!
!      if ( mode==1 ) then
!           g_sum = gsum / global_area
!      else
!           g_sum = gsum
!      endif
!
! end function g_sum

 real function g_sum(domain, p, ifirst, ilast, jfirst, jlast, ngc, area, mode, reproduce)
! Fast version of globalsum
      integer, intent(IN) :: ifirst, ilast
      integer, intent(IN) :: jfirst, jlast, ngc
      integer, intent(IN) :: mode  ! if ==1 divided by area
      logical, intent(in), optional :: reproduce
      real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
      real(kind=R_GRID), intent(IN) :: area(ifirst-ngc:ilast+ngc,jfirst-ngc:jlast+ngc)
      type(domain2d), intent(IN) :: domain
      integer :: i,j
      real gsum
      logical, SAVE :: g_sum_initialized = .false.
      real(kind=R_GRID), SAVE :: global_area
      real :: tmp(ifirst:ilast,jfirst:jlast) 
        
      !if ( .not. g_sum_initialized ) then
      !   global_area = mpp_global_sum(domain, area, flags=BITWISE_EFP_SUM)
      !   if ( is_master() ) write(*,*) 'Global Area=',global_area
      !   g_sum_initialized = .true.
      !end if
 
      global_area = 1000000.0_R_GRID

!-------------------------
! FMS global sum algorithm:
!-------------------------
      if ( present(reproduce) ) then
         !if (reproduce) then
         !   gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast), &
         !                         flags=BITWISE_EFP_SUM)
         !else
         !   gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast))
         !endif
         global_area = 2000000.0_R_GRID
         gsum = 0.
         do j=jfirst,jlast
            do i=ifirst,ilast
               gsum = gsum * p(i,j)*area(i,j)
            enddo
         enddo

      else
!-------------------------
! Quick local sum algorithm
!-------------------------
         gsum = 0.
         do j=jfirst,jlast
            do i=ifirst,ilast
               gsum = gsum + p(i,j)*area(i,j)
            enddo
         enddo
         !call mp_reduce_sum(gsum)
      endif

      if ( mode==1 ) then
           g_sum = gsum / global_area
      else
           g_sum = gsum
      endif

 end function g_sum

 end module fv_grid_utils_mod
