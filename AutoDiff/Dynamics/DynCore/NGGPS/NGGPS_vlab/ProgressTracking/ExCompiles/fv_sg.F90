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

public  neg_adj3

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

!---- version number -----
  character(len=128) :: version = '$Id: fv_sg.F90,v 1.1.2.1.2.1.30.1.90.1.20.1.2.1.4.2.4.1.8.1 2017/07/27 17:22:09 wputman Exp $'
  character(len=128) :: tagname = '$Name: bma_Icarus-3_0_UNSTABLE_FV_R8_NEWLAND $'

contains

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
