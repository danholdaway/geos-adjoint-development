! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Variable transforms on wind variables for fv3-jedi 
!> Daniel Holdaway, NASA/JCSDA

module wind_vt_mod

use fv3jedi_geom_mod
use mpp_domains_mod

implicit none
private
public psichi_to_udvd

contains

!----------------------------------------------------------------------------

subroutine psichi_to_udvd(geom,psi,chi,u,v)

 implicit none
 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(inout) :: psi(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Stream function
 real(kind=kind_real), intent(inout) :: chi(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(out)   ::   u(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed+1,1:geom%npz) !Dgrid winds (u)
 real(kind=kind_real), intent(out)   ::   v(geom%bd%isd:geom%bd%ied+1,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Dgrid winds (v)

 integer :: i,j,k
 real(kind=kind_real) :: chib(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)

 !       x-----------------x
 !       |                 |
 !       |                 |
 !       |                 |
 !     vd|        x        |
 !       |     psi,chi     |
 !       |                 |
 !       |                 |
 !       |                 |
 !       x-----------------x
 !               ud

 !Fill halos of psi and chi
 call mpp_update_domains(psi, geom%domain, complete=.true.)
 call mpp_update_domains(chi, geom%domain, complete=.true.)
 
 !Interpolate chi to the B grid
 call a2b_ord4(chi, chib, geom, geom%npx, geom%npy, geom%bd%isc, geom%bd%iec, geom%bd%jsc, geom%bd%jec, geom%halo)

 do k=1,geom%npz
   do j=geom%bd%jsc,geom%bd%jec
     do i=geom%bd%isc,geom%bd%iec

        u(i,j,k) =  (psi (i,j+1,k) - psi (i,j,k))/(geom%dyc(i,j)) + &
                    (chib(i+1,j,k) - chib(i,j,k))/(geom%dx (i,j))
        v(i,j,k) = -(psi (i+1,j,k) - psi (i,j,k))/(geom%dxc(i,j)) + &
                    (chib(i,j+1,k) - chib(i,j,k))/(geom%dy (i,j))

     enddo
   enddo
 enddo 

endsubroutine psichi_to_udvd

!----------------------------------------------------------------------------

  subroutine a2b_ord4(qin, qout, geom, npx, npy, is, ie, js, je, ng)

  implicit none

  integer, intent(IN):: npx, npy, is, ie, js, je, ng
  real(kind=kind_real), intent(IN)::  qin(is-ng:ie+ng,js-ng:je+ng)   ! A-grid field
  real(kind=kind_real), intent(OUT):: qout(is-ng:ie+ng,js-ng:je+ng)   ! Output  B-grid field
  type(fv3jedi_geom), intent(IN), target :: geom

  real(kind=kind_real) qx(is:ie+1,js-ng:je+ng)
  real(kind=kind_real) qy(is-ng:ie+ng,js:je+1)
  real(kind=kind_real) qxx(is-ng:ie+ng,js-ng:je+ng)
  real(kind=kind_real) qyy(is-ng:ie+ng,js-ng:je+ng)
  real(kind=kind_real) g_in, g_ou
  real(kind=kind_real):: p0(2)
  real(kind=kind_real):: q1(is-1:ie+1), q2(js-1:je+1)
  integer:: i, j, is1, js1, is2, js2, ie1, je1

  real(kind=kind_real), parameter :: c1 =  2./3.
  real(kind=kind_real), parameter :: c2 = -1./6.
  real(kind=kind_real), parameter :: r3 = 1./3.
  real(kind=kind_real), parameter :: a1 =  0.5625  !  9/16
  real(kind=kind_real), parameter :: a2 = -0.0625  ! -1/16
  real(kind=kind_real), parameter :: b1 =  7./12.     ! 0.58333333
  real(kind=kind_real), parameter :: b2 = -1./12.

    is1 = max(1,is-1)
    js1 = max(1,js-1)
    is2 = max(2,is)
    js2 = max(2,js)

    ie1 = min(npx-1,ie+1)
    je1 = min(npy-1,je+1)

! Corners:
! 3-way extrapolation

    if ( geom%sw_corner ) then
          p0(1:2) = geom%grid(1,1,1:2)
        qout(1,1) = (extrap_corner(p0, geom%agrid(1,1,1:2), geom%agrid( 2, 2,1:2), qin(1,1), qin( 2, 2)) + &
                     extrap_corner(p0, geom%agrid(0,1,1:2), geom%agrid(-1, 2,1:2), qin(0,1), qin(-1, 2)) + &
                     extrap_corner(p0, geom%agrid(1,0,1:2), geom%agrid( 2,-1,1:2), qin(1,0), qin( 2,-1)))*r3

    endif
    if ( geom%se_corner ) then
            p0(1:2) = geom%grid(npx,1,1:2)
        qout(npx,1) = (extrap_corner(p0, geom%agrid(npx-1,1,1:2), geom%agrid(npx-2, 2,1:2), &
                                     qin(npx-1,1), qin(npx-2, 2)) + &
                       extrap_corner(p0, geom%agrid(npx-1,0,1:2), geom%agrid(npx-2,-1,1:2), &
                                     qin(npx-1,0), qin(npx-2,-1)) + &
                       extrap_corner(p0, geom%agrid(npx  ,1,1:2), geom%agrid(npx+1, 2,1:2), &
                                     qin(npx  ,1), qin(npx+1, 2)))*r3
    endif
    if ( geom%ne_corner ) then
              p0(1:2) = geom%grid(npx,npy,1:2)
        qout(npx,npy) = (extrap_corner(p0, geom%agrid(npx-1,npy-1,1:2), geom%agrid(npx-2,npy-2,1:2), &
                                       qin(npx-1,npy-1), qin(npx-2,npy-2)) + &
                         extrap_corner(p0, geom%agrid(npx  ,npy-1,1:2), geom%agrid(npx+1,npy-2,1:2), &
                                       qin(npx  ,npy-1), qin(npx+1,npy-2)) + &
                         extrap_corner(p0, geom%agrid(npx-1,npy  ,1:2), geom%agrid(npx-2,npy+1,1:2), &
                                       qin(npx-1,npy  ), qin(npx-2,npy+1)))*r3
    endif
    if ( geom%nw_corner ) then
            p0(1:2) = geom%grid(1,npy,1:2)
        qout(1,npy) = (extrap_corner(p0, geom%agrid(1,npy-1,1:2), geom%agrid( 2,npy-2,1:2), qin(1,npy-1), &
                                     qin( 2,npy-2)) + &
                       extrap_corner(p0, geom%agrid(0,npy-1,1:2), geom%agrid(-1,npy-2,1:2), qin(0,npy-1), &
                                     qin(-1,npy-2)) + &
                       extrap_corner(p0, geom%agrid(1,npy,  1:2), geom%agrid( 2,npy+1,1:2), qin(1,npy  ), &
                                     qin( 2,npy+1)))*r3
    endif

!------------
! X-Interior:
!------------
    do j=max(1,js-2),min(npy-1,je+2)
       do i=max(3,is), min(npx-2,ie+1)
          qx(i,j) = b2*(qin(i-2,j)+qin(i+1,j)) + b1*(qin(i-1,j)+qin(i,j))
       enddo
    enddo

    ! *** West Edges:
    if ( is==1 ) then
       do j=js1, je1
          q2(j) = (qin(0,j)*geom%dxa(1,j) + qin(1,j)*geom%dxa(0,j))/(geom%dxa(0,j) + geom%dxa(1,j))
       enddo
       do j=js2, je1
          qout(1,j) = geom%edge_w(j)*q2(j-1) + (1.-geom%edge_w(j))*q2(j)
       enddo
!
       do j=max(1,js-2),min(npy-1,je+2)
             g_in = geom%dxa(2,j) / geom%dxa(1,j)
             g_ou = geom%dxa(-1,j) / geom%dxa(0,j)
          qx(1,j) = 0.5*( ((2.+g_in)*qin(1,j)-qin( 2,j))/(1.+g_in) +          &
                          ((2.+g_ou)*qin(0,j)-qin(-1,j))/(1.+g_ou) )
          qx(2,j) = ( 3.*(g_in*qin(1,j)+qin(2,j))-(g_in*qx(1,j)+qx(3,j)) ) / (2.+2.*g_in)
       enddo
    endif

    ! East Edges:
    if ( (ie+1)==npx ) then
       do j=js1, je1
          q2(j) = (qin(npx-1,j)*geom%dxa(npx,j) + qin(npx,j)*geom%dxa(npx-1,j))/(geom%dxa(npx-1,j) + geom%dxa(npx,j))
       enddo
       do j=js2, je1
          qout(npx,j) = geom%edge_e(j)*q2(j-1) + (1.-geom%edge_e(j))*q2(j)
       enddo
!
       do j=max(1,js-2),min(npy-1,je+2)
              g_in = geom%dxa(npx-2,j) / geom%dxa(npx-1,j)
              g_ou = geom%dxa(npx+1,j) / geom%dxa(npx,j)
          qx(npx,j) = 0.5*( ((2.+g_in)*qin(npx-1,j)-qin(npx-2,j))/(1.+g_in) +          &
                            ((2.+g_ou)*qin(npx,  j)-qin(npx+1,j))/(1.+g_ou) )
          qx(npx-1,j) = (3.*(qin(npx-2,j)+g_in*qin(npx-1,j)) - (g_in*qx(npx,j)+qx(npx-2,j)))/(2.+2.*g_in)
       enddo
    endif
    
!------------
! Y-Interior:
!------------

    do j=max(3,js),min(npy-2,je+1)
       do i=max(1,is-2), min(npx-1,ie+2)
          qy(i,j) = b2*(qin(i,j-2)+qin(i,j+1)) + b1*(qin(i,j-1) + qin(i,j))
       enddo
    enddo

    ! South Edges:
    if ( js==1 ) then
       do i=is1, ie1
          q1(i) = (qin(i,0)*geom%dya(i,1) + qin(i,1)*geom%dya(i,0))/(geom%dya(i,0) + geom%dya(i,1))
       enddo
       do i=is2, ie1
          qout(i,1) = geom%edge_s(i)*q1(i-1) + (1.-geom%edge_s(i))*q1(i)
       enddo
!
       do i=max(1,is-2),min(npx-1,ie+2)
             g_in = geom%dya(i,2) / geom%dya(i,1)
             g_ou = geom%dya(i,-1) / geom%dya(i,0)
          qy(i,1) = 0.5*( ((2.+g_in)*qin(i,1)-qin(i,2))/(1.+g_in) +          &
                          ((2.+g_ou)*qin(i,0)-qin(i,-1))/(1.+g_ou) )
          qy(i,2) = (3.*(g_in*qin(i,1)+qin(i,2)) - (g_in*qy(i,1)+qy(i,3)))/(2.+2.*g_in)
       enddo
    endif

    ! North Edges:
    if ( (je+1)==npy ) then
       do i=is1, ie1
          q1(i) = (qin(i,npy-1)*geom%dya(i,npy) + qin(i,npy)*geom%dya(i,npy-1))/(geom%dya(i,npy-1)+geom%dya(i,npy))
       enddo
       do i=is2, ie1
          qout(i,npy) = geom%edge_n(i)*q1(i-1) + (1.-geom%edge_n(i))*q1(i)
       enddo
!
       do i=max(1,is-2),min(npx-1,ie+2)
              g_in = geom%dya(i,npy-2) / geom%dya(i,npy-1)
              g_ou = geom%dya(i,npy+1) / geom%dya(i,npy)
          qy(i,npy) = 0.5*( ((2.+g_in)*qin(i,npy-1)-qin(i,npy-2))/(1.+g_in) +          &
                            ((2.+g_ou)*qin(i,npy  )-qin(i,npy+1))/(1.+g_ou) )
          qy(i,npy-1) = (3.*(qin(i,npy-2)+g_in*qin(i,npy-1)) - (g_in*qy(i,npy)+qy(i,npy-2)))/(2.+2.*g_in)
       enddo
    endif

!--------------------------------------

    do j=max(3,js),min(npy-2,je+1)
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,j) = a2*(qx(i,j-2)+qx(i,j+1)) + a1*(qx(i,j-1)+qx(i,j))
       enddo
    enddo

    if ( js==1 ) then
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,2) = c1*(qx(i,1)+qx(i,2))+c2*(qout(i,1)+qxx(i,3))
       enddo
    endif
    if ( (je+1)==npy ) then
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,npy-1) = c1*(qx(i,npy-2)+qx(i,npy-1))+c2*(qout(i,npy)+qxx(i,npy-2))
       enddo
    endif

    
    do j=max(2,js),min(npy-1,je+1)
       do i=max(3,is),min(npx-2,ie+1)
          qyy(i,j) = a2*(qy(i-2,j)+qy(i+1,j)) + a1*(qy(i-1,j)+qy(i,j))
       enddo
       if ( is==1 ) qyy(2,j) = c1*(qy(1,j)+qy(2,j))+c2*(qout(1,j)+qyy(3,j))
       if((ie+1)==npx) qyy(npx-1,j) = c1*(qy(npx-2,j)+qy(npx-1,j))+c2*(qout(npx,j)+qyy(npx-2,j))
 
       do i=max(2,is),min(npx-1,ie+1)
          qout(i,j) = 0.5*(qxx(i,j) + qyy(i,j))   ! averaging
       enddo
    enddo
    
  end subroutine a2b_ord4

!----------------------------------------------------------------------------

real(kind=kind_real) function extrap_corner ( p0, p1, p2, q1, q2 )

    implicit none
    real(kind=kind_real), intent(in ), dimension(2) :: p0, p1, p2
    real(kind=kind_real), intent(in ) :: q1, q2
    real(kind=kind_real) :: x1, x2

    x1 = great_circle_dist( p1, p0 )
    x2 = great_circle_dist( p2, p0 )

    extrap_corner = q1 + x1/(x2-x1) * (q1-q2)

  end function extrap_corner
  
!----------------------------------------------------------------------------

  real(kind=kind_real) function great_circle_dist( q1, q2 )

       implicit none
       real(kind=kind_real), intent(IN)           :: q1(2), q2(2)
 
       real (kind=kind_real):: p1(2), p2(2)
       integer n
 
       do n=1,2
          p1(n) = q1(n)
          p2(n) = q2(n)
       enddo
 
       great_circle_dist = asin( sqrt( sin((p1(2)-p2(2))/2.)**2 + cos(p1(2))*cos(p2(2))*   &
                          sin((p1(1)-p2(1))/2.)**2 ) ) * 2.
  
   end function great_circle_dist

!----------------------------------------------------------------------------

end module wind_vt_mod
