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
module fv_fill_mod

   use mpp_domains_mod,     only: mpp_update_domains, domain2D

   implicit none
   public fillz
   public fill2D


contains

 subroutine fillz(im, km, nq, q, dp)
   integer,  intent(in):: im                ! No. of longitudes
   integer,  intent(in):: km                ! No. of levels
   integer,  intent(in):: nq                ! Total number of tracers
   real , intent(in)::  dp(im,km)       ! pressure thickness
   real , intent(inout) :: q(im,km,nq)   ! tracer mixing ratio
! !LOCAL VARIABLES:
   logical:: zfix(im)
   real ::  dm(km)
   integer i, k, ic, k1
   real  qup, qly, dup, dq, sum0, sum1, fac

   !Stubbed, highly nonlinear process

 end subroutine fillz


 subroutine fill2D(is, ie, js, je, ng, km, q, delp, area, domain, nested, npx, npy)
! This is a diffusive type filling algorithm
 type(domain2D), intent(INOUT) :: domain
 integer, intent(in):: is, ie, js, je, ng, km, npx, npy
 logical, intent(IN):: nested
 real, intent(in):: area(is-ng:ie+ng, js-ng:je+ng)
 real, intent(in):: delp(is-ng:ie+ng, js-ng:je+ng, km)
 real, intent(inout):: q(is-ng:ie+ng, js-ng:je+ng, km)
! LOCAL VARIABLES:
 real, dimension(is-ng:ie+ng, js-ng:je+ng,km):: qt
 real, dimension(is:ie+1, js:je):: fx
 real, dimension(is:ie, js:je+1):: fy
 real, parameter:: dif = 0.25
 integer:: i, j, k
 integer :: is1, ie1, js1, je1

 !Stubbed, highly nonlinear process
 
 end subroutine fill2D

end module fv_fill_mod
