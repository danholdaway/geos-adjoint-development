!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.13 (r6666) - 27 Nov 2017 17:09
!
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
MODULE FV_GRID_UTILS_MOD_DIFF
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
  USE MPP_MOD, ONLY : fatal, mpp_error, warning
  USE MPP_DOMAINS_MOD_DIFF, ONLY : mpp_update_domains, &
& mpp_update_domains_adm, dgrid_ne, mpp_global_sum
  USE MPP_DOMAINS_MOD_DIFF, ONLY : domain2d, bitwise_efp_sum
  USE FV_ARRAYS_MOD_DIFF, ONLY : fv_atmos_type, fv_grid_type, &
& fv_grid_bounds_type, r_grid
  USE FV_MP_MOD_DIFF, ONLY : ng, is_master
  USE FV_MP_MOD_DIFF, ONLY : mp_reduce_sum
  USE FV_MP_MOD_DIFF, ONLY : fill_corners, fill_corners_adm, xdir, ydir
  USE FV_TIMING_MOD, ONLY : timing_on, timing_off
  IMPLICIT NONE
  PRIVATE 
  LOGICAL :: symm_grid
! Higher precision (kind=16) for grid geometrical factors:
  INTEGER, PARAMETER :: f_p=SELECTED_REAL_KIND(20)
  REAL, PARAMETER :: big_number=1.d8
  REAL, PARAMETER :: tiny_number=1.d-8
  REAL, PARAMETER :: ptop_min=1.d-8
  PUBLIC f_p
!CLEANUP: OK to keep since they are constants?
  PUBLIC ptop_min, big_number
  PUBLIC grid_utils_init, grid_utils_end, great_circle_dist, &
& cubed_to_latlon
  PUBLIC c2l_ord2, g_sum
  PUBLIC c2l_ord2_fwd, c2l_ord2_bwd, g_sum_adm
  EXTERNAL C2L_ORD2_ADM

CONTAINS
  SUBROUTINE GRID_UTILS_INIT(atm, npx, npy, npz, non_ortho, grid_type, &
&   c2l_order)
    IMPLICIT NONE
! Initialize 2D memory and geometrical factors
    TYPE(FV_ATMOS_TYPE), INTENT(INOUT), TARGET :: atm
    LOGICAL, INTENT(IN) :: non_ortho
    INTEGER, INTENT(IN) :: npx, npy, npz
    INTEGER, INTENT(IN) :: grid_type, c2l_order
  END SUBROUTINE GRID_UTILS_INIT
  SUBROUTINE GRID_UTILS_END()
    IMPLICIT NONE
  END SUBROUTINE GRID_UTILS_END
  REAL FUNCTION GREAT_CIRCLE_DIST(q1, q2, radius)
    IMPLICIT NONE
    REAL(kind=r_grid), INTENT(IN) :: q1(2), q2(2)
    REAL(kind=r_grid), INTENT(IN), OPTIONAL :: radius
    REAL(f_p) :: p1(2), p2(2)
    REAL(f_p) :: beta
    INTEGER :: n
    INTRINSIC SIN
    INTRINSIC COS
    INTRINSIC SQRT
    INTRINSIC ASIN
    INTRINSIC PRESENT
    DO n=1,2
      p1(n) = q1(n)
      p2(n) = q2(n)
    END DO
    beta = ASIN(SQRT(SIN((p1(2)-p2(2))/2.)**2+COS(p1(2))*COS(p2(2))*SIN(&
&     (p1(1)-p2(1))/2.)**2))*2.
    IF (PRESENT(radius)) THEN
      great_circle_dist = radius*beta
    ELSE
! Returns the angle
      great_circle_dist = beta
    END IF
  END FUNCTION GREAT_CIRCLE_DIST
  SUBROUTINE CUBED_TO_LATLON(u, v, ua, va, gridstruct, npx, npy, km, &
&   mode, grid_type, domain, nested, c2l_ord, bd)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: km, npx, npy, grid_type, c2l_ord
! update if present
    INTEGER, INTENT(IN) :: mode
    TYPE(FV_GRID_TYPE), INTENT(IN) :: gridstruct
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, km)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    LOGICAL, INTENT(IN) :: nested
    IF (c2l_ord .EQ. 2) THEN
      CALL C2L_ORD2(u, v, ua, va, gridstruct, km, grid_type, bd, .false.&
&            )
    ELSE
      CALL C2L_ORD4(u, v, ua, va, gridstruct, npx, npy, km, grid_type, &
&             domain, nested, mode, bd)
    END IF
  END SUBROUTINE CUBED_TO_LATLON
  SUBROUTINE C2L_ORD4(u, v, ua, va, gridstruct, npx, npy, km, grid_type&
&   , domain, nested, mode, bd)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: km, npx, npy, grid_type
! update if present
    INTEGER, INTENT(IN) :: mode
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    REAL, INTENT(INOUT) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, km)
    REAL, INTENT(INOUT) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    TYPE(DOMAIN2D), INTENT(INOUT) :: domain
    LOGICAL, INTENT(IN) :: nested
! Local
! 4-pt Lagrange interpolation
    REAL, SAVE :: a1=0.5625
    REAL, SAVE :: a2=-0.0625
    REAL, SAVE :: c1=1.125
    REAL, SAVE :: c2=-0.125
    REAL :: utmp(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: vtmp(bd%is:bd%ie+1, bd%js:bd%je)
    REAL :: wu(bd%is:bd%ie, bd%js:bd%je+1)
    REAL :: wv(bd%is:bd%ie+1, bd%js:bd%je)
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    IF (mode .GT. 0) THEN
      CALL TIMING_ON('COMM_TOTAL')
      CALL MPP_UPDATE_DOMAINS(u, v, domain, gridtype=dgrid_ne)
      CALL TIMING_OFF('COMM_TOTAL')
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,km,npx,npy,grid_type,nested,c2,c1, &
!$OMP                                  u,v,gridstruct,ua,va,a1,a2)         &
!$OMP                          private(utmp, vtmp, wu, wv)
    DO k=1,km
      IF (grid_type .LT. 4) THEN
!nested
        IF (nested) THEN
          IF (1 .LT. js) THEN
            max1 = js
          ELSE
            max1 = 1
          END IF
          IF (npy - 1 .GT. je) THEN
            min1 = je
          ELSE
            min1 = npy - 1
          END IF
          DO j=max1,min1
            IF (1 .LT. is) THEN
              max2 = is
            ELSE
              max2 = 1
            END IF
            IF (npx - 1 .GT. ie) THEN
              min2 = ie
            ELSE
              min2 = npx - 1
            END IF
            DO i=max2,min2
              utmp(i, j) = c2*(u(i, j-1, k)+u(i, j+2, k)) + c1*(u(i, j, &
&               k)+u(i, j+1, k))
              vtmp(i, j) = c2*(v(i-1, j, k)+v(i+2, j, k)) + c1*(v(i, j, &
&               k)+v(i+1, j, k))
            END DO
          END DO
        ELSE
          IF (2 .LT. js) THEN
            max3 = js
          ELSE
            max3 = 2
          END IF
          IF (npy - 2 .GT. je) THEN
            min3 = je
          ELSE
            min3 = npy - 2
          END IF
          DO j=max3,min3
            IF (2 .LT. is) THEN
              max4 = is
            ELSE
              max4 = 2
            END IF
            IF (npx - 2 .GT. ie) THEN
              min4 = ie
            ELSE
              min4 = npx - 2
            END IF
            DO i=max4,min4
              utmp(i, j) = c2*(u(i, j-1, k)+u(i, j+2, k)) + c1*(u(i, j, &
&               k)+u(i, j+1, k))
              vtmp(i, j) = c2*(v(i-1, j, k)+v(i+2, j, k)) + c1*(v(i, j, &
&               k)+v(i+1, j, k))
            END DO
          END DO
          IF (js .EQ. 1) THEN
            DO i=is,ie+1
              wv(i, 1) = v(i, 1, k)*gridstruct%dy(i, 1)
            END DO
            DO i=is,ie
              vtmp(i, 1) = 2.*(wv(i, 1)+wv(i+1, 1))/(gridstruct%dy(i, 1)&
&               +gridstruct%dy(i+1, 1))
              utmp(i, 1) = 2.*(u(i, 1, k)*gridstruct%dx(i, 1)+u(i, 2, k)&
&               *gridstruct%dx(i, 2))/(gridstruct%dx(i, 1)+gridstruct%dx&
&               (i, 2))
            END DO
          END IF
!!!         vtmp(i,1) = (wv(i,1) + wv(i+1,1)) * gridstruct%rdya(i,1)
!!!         utmp(i,1) = (u(i,1,k)*gridstruct%dx(i,1) + u(i,2,k)*gridstruct%dx(i,2)) * gridstruct%rdxa(i,1)
          IF (je + 1 .EQ. npy) THEN
            j = npy - 1
            DO i=is,ie+1
              wv(i, j) = v(i, j, k)*gridstruct%dy(i, j)
            END DO
            DO i=is,ie
              vtmp(i, j) = 2.*(wv(i, j)+wv(i+1, j))/(gridstruct%dy(i, j)&
&               +gridstruct%dy(i+1, j))
              utmp(i, j) = 2.*(u(i, j, k)*gridstruct%dx(i, j)+u(i, j+1, &
&               k)*gridstruct%dx(i, j+1))/(gridstruct%dx(i, j)+&
&               gridstruct%dx(i, j+1))
            END DO
          END IF
!!!         vtmp(i,j) = (wv(i,j) + wv(i+1,j)) * gridstruct%rdya(i,j)
!!!         utmp(i,j) = (u(i,j,k)*gridstruct%dx(i,j) + u(i,j+1,k)*gridstruct%dx(i,j+1)) * gridstruct%rdxa(i,j)
          IF (is .EQ. 1) THEN
            i = 1
            DO j=js,je
              wv(1, j) = v(1, j, k)*gridstruct%dy(1, j)
              wv(2, j) = v(2, j, k)*gridstruct%dy(2, j)
            END DO
            DO j=js,je+1
              wu(i, j) = u(i, j, k)*gridstruct%dx(i, j)
            END DO
            DO j=js,je
              utmp(i, j) = 2.*(wu(i, j)+wu(i, j+1))/(gridstruct%dx(i, j)&
&               +gridstruct%dx(i, j+1))
              vtmp(i, j) = 2.*(wv(1, j)+wv(2, j))/(gridstruct%dy(1, j)+&
&               gridstruct%dy(2, j))
            END DO
          END IF
!!!      utmp(i,j) = (wu(i,j) + wu(i,  j+1)) * gridstruct%rdxa(i,j)
!!!      vtmp(i,j) = (wv(i,j) + wv(i+1,j  )) * gridstruct%rdya(i,j)
          IF (ie + 1 .EQ. npx) THEN
            i = npx - 1
            DO j=js,je
              wv(i, j) = v(i, j, k)*gridstruct%dy(i, j)
              wv(i+1, j) = v(i+1, j, k)*gridstruct%dy(i+1, j)
            END DO
            DO j=js,je+1
              wu(i, j) = u(i, j, k)*gridstruct%dx(i, j)
            END DO
            DO j=js,je
              utmp(i, j) = 2.*(wu(i, j)+wu(i, j+1))/(gridstruct%dx(i, j)&
&               +gridstruct%dx(i, j+1))
              vtmp(i, j) = 2.*(wv(i, j)+wv(i+1, j))/(gridstruct%dy(i, j)&
&               +gridstruct%dy(i+1, j))
            END DO
          END IF
        END IF
!!!      utmp(i,j) = (wu(i,j) + wu(i,  j+1)) * gridstruct%rdxa(i,j)
!!!      vtmp(i,j) = (wv(i,j) + wv(i+1,j  )) * gridstruct%rdya(i,j)
!Transform local a-grid winds into latitude-longitude coordinates
        DO j=js,je
          DO i=is,ie
            ua(i, j, k) = gridstruct%a11(i, j)*utmp(i, j) + gridstruct%&
&             a12(i, j)*vtmp(i, j)
            va(i, j, k) = gridstruct%a21(i, j)*utmp(i, j) + gridstruct%&
&             a22(i, j)*vtmp(i, j)
          END DO
        END DO
      ELSE
! Simple Cartesian Geometry:
        DO j=js,je
          DO i=is,ie
            ua(i, j, k) = a2*(u(i, j-1, k)+u(i, j+2, k)) + a1*(u(i, j, k&
&             )+u(i, j+1, k))
            va(i, j, k) = a2*(v(i-1, j, k)+v(i+2, j, k)) + a1*(v(i, j, k&
&             )+v(i+1, j, k))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE C2L_ORD4
!  Differentiation of c2l_ord2 in reverse (adjoint) mode, forward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod
!.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mo
!d.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix
!_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_S
!uper fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_dynamics_mod.geos_to_fv3 fv_dynamics_mod.fv3_to_geos fv_gr
!id_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mo
!d.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_t
!racer_fb fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limi
!ters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moi
!st_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2
!d fv_tracer2d_mod.tracer_2d_nested nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_
!Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_
!solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_hal
!o_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect 
!sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_cor
!e_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_
!corners_fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.grea
!t_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: u v ua va
!   with respect to varying inputs: u v ua va
  SUBROUTINE C2L_ORD2_FWD(u, v, ua, va, gridstruct, km, grid_type, bd, &
&   do_halo)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: km, grid_type
    REAL, INTENT(IN) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, km)
    REAL, INTENT(IN) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, km)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, INTENT(IN) :: do_halo
!
    REAL :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL :: va(bd%isd:bd%ied, bd%jsd:bd%jed, km)
!--------------------------------------------------------------
! Local
    REAL :: wu(bd%is-1:bd%ie+1, bd%js-1:bd%je+2)
    REAL :: wv(bd%is-1:bd%ie+2, bd%js-1:bd%je+1)
    REAL :: u1(bd%is-1:bd%ie+1), v1(bd%is-1:bd%ie+1)
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    REAL, DIMENSION(:, :), POINTER :: a11, a12, a21, a22
    REAL, DIMENSION(:, :), POINTER :: dx, dy, rdxa, rdya
    a11 => gridstruct%a11
    a12 => gridstruct%a12
    a21 => gridstruct%a21
    a22 => gridstruct%a22
    dx => gridstruct%dx
    dy => gridstruct%dy
    IF (do_halo) THEN
      is = bd%is - 1
      ie = bd%ie + 1
      js = bd%js - 1
      je = bd%je + 1
    ELSE
      is = bd%is
      ie = bd%ie
      js = bd%js
      je = bd%je
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,km,grid_type,u,dx,v,dy,ua,va,a11,a12,a21,a22) &
!$OMP                          private(u1, v1, wu, wv)
    DO k=1,km
      IF (grid_type .LT. 4) THEN
        DO j=js,je+1
          DO i=is,ie
            wu(i, j) = u(i, j, k)*dx(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            wv(i, j) = v(i, j, k)*dy(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
            u1(i) = 2.*(wu(i, j)+wu(i, j+1))/(dx(i, j)+dx(i, j+1))
            v1(i) = 2.*(wv(i, j)+wv(i+1, j))/(dy(i, j)+dy(i+1, j))
!!!          u1(i) = (wu(i,j) + wu(i,j+1)) * rdxa(i,j)
!!!          v1(i) = (wv(i,j) + wv(i+1,j)) * rdya(i,j)
! Cubed (cell center co-variant winds) to lat-lon:
            CALL PUSHREALARRAY(ua(i, j, k))
            ua(i, j, k) = a11(i, j)*u1(i) + a12(i, j)*v1(i)
            CALL PUSHREALARRAY(va(i, j, k))
            va(i, j, k) = a21(i, j)*u1(i) + a22(i, j)*v1(i)
          END DO
        END DO
        CALL PUSHCONTROL(1,1)
      ELSE
! 2nd order:
        DO j=js,je
          DO i=is,ie
            CALL PUSHREALARRAY(ua(i, j, k))
            ua(i, j, k) = 0.5*(u(i, j, k)+u(i, j+1, k))
            CALL PUSHREALARRAY(va(i, j, k))
            va(i, j, k) = 0.5*(v(i, j, k)+v(i+1, j, k))
          END DO
        END DO
        CALL PUSHCONTROL(1,0)
      END IF
    END DO
    CALL PUSHINTEGER(je)
    !CALL PUSHPOINTER8(C_LOC(a12))
    !CALL PUSHPOINTER8(C_LOC(a11))
    CALL PUSHINTEGER(is)
    CALL PUSHINTEGER(ie)
    !CALL PUSHPOINTER8(C_LOC(dy))
    !CALL PUSHPOINTER8(C_LOC(dx))
    !CALL PUSHPOINTER8(C_LOC(a22))
    !CALL PUSHPOINTER8(C_LOC(a21))
    CALL PUSHINTEGER(js)
  END SUBROUTINE C2L_ORD2_FWD
!  Differentiation of c2l_ord2 in reverse (adjoint) mode, backward sweep (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mo
!d.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_m
!od.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mi
!x_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_
!Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_dynamics_mod.geos_to_fv3 fv_dynamics_mod.fv3_to_geos fv_g
!rid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_m
!od.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_
!tracer_fb fv_mapz_mod.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_lim
!iters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.mo
!ist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_
!2d fv_tracer2d_mod.tracer_2d_nested nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem
!_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0
!_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_ha
!lo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect
! sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_co
!re_mod.compute_divergence_damping sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy
!_corners_fb tp_core_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.gre
!at_circle_dist sw_core_mod.edge_interpolate4)):
!   gradient     of useful results: u v ua va
!   with respect to varying inputs: u v ua va
  SUBROUTINE C2L_ORD2_BWD(u, u_ad, v, v_ad, ua, ua_ad, va, va_ad, &
&   gridstruct, km, grid_type, bd, do_halo)
    !USE ISO_C_BINDING
    !USE ADMM_TAPENADE_INTERFACE
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: km, grid_type
    REAL, INTENT(IN) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, km)
    REAL :: u_ad(bd%isd:bd%ied, bd%jsd:bd%jed+1, km)
    REAL, INTENT(IN) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, km)
    REAL :: v_ad(bd%isd:bd%ied+1, bd%jsd:bd%jed, km)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, INTENT(IN) :: do_halo
    REAL :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL :: ua_ad(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL :: va(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL :: va_ad(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL :: wu(bd%is-1:bd%ie+1, bd%js-1:bd%je+2)
    REAL :: wu_ad(bd%is-1:bd%ie+1, bd%js-1:bd%je+2)
    REAL :: wv(bd%is-1:bd%ie+2, bd%js-1:bd%je+1)
    REAL :: wv_ad(bd%is-1:bd%ie+2, bd%js-1:bd%je+1)
    REAL :: u1(bd%is-1:bd%ie+1), v1(bd%is-1:bd%ie+1)
    REAL :: u1_ad(bd%is-1:bd%ie+1), v1_ad(bd%is-1:bd%ie+1)
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    REAL, DIMENSION(:, :), POINTER :: a11, a12, a21, a22
    REAL, DIMENSION(:, :), POINTER :: dx, dy, rdxa, rdya
    REAL :: temp_ad
    REAL :: temp_ad0
    INTEGER :: branch
    !TYPE(C_PTR) :: cptr
    !INTEGER :: unknown_shape_in_c2l_ord2
    CALL POPINTEGER(js)
    !CALL POPPOINTER8(cptr)
    a21 => gridstruct%a21 ! (/unknown_shape_in_c2l_ord2/))
    !CALL POPPOINTER8(cptr)
    a22 => gridstruct%a22 ! (/unknown_shape_in_c2l_ord2/))
    !CALL POPPOINTER8(cptr)
    dx => gridstruct%dx ! (/unknown_shape_in_c2l_ord2/))
    !CALL POPPOINTER8(cptr)
    dy => gridstruct%dy ! (/unknown_shape_in_c2l_ord2/))
    CALL POPINTEGER(ie)
    CALL POPINTEGER(is)
    !CALL POPPOINTER8(cptr)
    a11 => gridstruct%a11 ! (/unknown_shape_in_c2l_ord2/))
    !CALL POPPOINTER8(cptr)
    a12 => gridstruct%a12 ! (/unknown_shape_in_c2l_ord2/))
    CALL POPINTEGER(je)
    wu_ad = 0.0
    v1_ad = 0.0
    wv_ad = 0.0
    u1_ad = 0.0
    DO k=km,1,-1
      CALL POPCONTROL(1,branch)
      IF (branch .EQ. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(va(i, j, k))
            v_ad(i, j, k) = v_ad(i, j, k) + 0.5*va_ad(i, j, k)
            v_ad(i+1, j, k) = v_ad(i+1, j, k) + 0.5*va_ad(i, j, k)
            va_ad(i, j, k) = 0.0
            CALL POPREALARRAY(ua(i, j, k))
            u_ad(i, j, k) = u_ad(i, j, k) + 0.5*ua_ad(i, j, k)
            u_ad(i, j+1, k) = u_ad(i, j+1, k) + 0.5*ua_ad(i, j, k)
            ua_ad(i, j, k) = 0.0
          END DO
        END DO
      ELSE
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREALARRAY(va(i, j, k))
            u1_ad(i) = u1_ad(i) + a11(i, j)*ua_ad(i, j, k) + a21(i, j)*&
&             va_ad(i, j, k)
            v1_ad(i) = v1_ad(i) + a12(i, j)*ua_ad(i, j, k) + a22(i, j)*&
&             va_ad(i, j, k)
            va_ad(i, j, k) = 0.0
            CALL POPREALARRAY(ua(i, j, k))
            ua_ad(i, j, k) = 0.0
            temp_ad = 2.*v1_ad(i)/(dy(i, j)+dy(i+1, j))
            wv_ad(i, j) = wv_ad(i, j) + temp_ad
            wv_ad(i+1, j) = wv_ad(i+1, j) + temp_ad
            v1_ad(i) = 0.0
            temp_ad0 = 2.*u1_ad(i)/(dx(i, j)+dx(i, j+1))
            wu_ad(i, j) = wu_ad(i, j) + temp_ad0
            wu_ad(i, j+1) = wu_ad(i, j+1) + temp_ad0
            u1_ad(i) = 0.0
          END DO
        END DO
        DO j=je,js,-1
          DO i=ie+1,is,-1
            v_ad(i, j, k) = v_ad(i, j, k) + dy(i, j)*wv_ad(i, j)
            wv_ad(i, j) = 0.0
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ie,is,-1
            u_ad(i, j, k) = u_ad(i, j, k) + dx(i, j)*wu_ad(i, j)
            wu_ad(i, j) = 0.0
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE C2L_ORD2_BWD
  SUBROUTINE C2L_ORD2(u, v, ua, va, gridstruct, km, grid_type, bd, &
&   do_halo)
    IMPLICIT NONE
    TYPE(FV_GRID_BOUNDS_TYPE), INTENT(IN) :: bd
    INTEGER, INTENT(IN) :: km, grid_type
    REAL, INTENT(IN) :: u(bd%isd:bd%ied, bd%jsd:bd%jed+1, km)
    REAL, INTENT(IN) :: v(bd%isd:bd%ied+1, bd%jsd:bd%jed, km)
    TYPE(FV_GRID_TYPE), INTENT(IN), TARGET :: gridstruct
    LOGICAL, INTENT(IN) :: do_halo
!
    REAL, INTENT(OUT) :: ua(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    REAL, INTENT(OUT) :: va(bd%isd:bd%ied, bd%jsd:bd%jed, km)
!--------------------------------------------------------------
! Local
    REAL :: wu(bd%is-1:bd%ie+1, bd%js-1:bd%je+2)
    REAL :: wv(bd%is-1:bd%ie+2, bd%js-1:bd%je+1)
    REAL :: u1(bd%is-1:bd%ie+1), v1(bd%is-1:bd%ie+1)
    INTEGER :: i, j, k
    INTEGER :: is, ie, js, je
    REAL, DIMENSION(:, :), POINTER :: a11, a12, a21, a22
    REAL, DIMENSION(:, :), POINTER :: dx, dy, rdxa, rdya
    a11 => gridstruct%a11
    a12 => gridstruct%a12
    a21 => gridstruct%a21
    a22 => gridstruct%a22
    dx => gridstruct%dx
    dy => gridstruct%dy
    rdxa => gridstruct%rdxa
    rdya => gridstruct%rdya
    IF (do_halo) THEN
      is = bd%is - 1
      ie = bd%ie + 1
      js = bd%js - 1
      je = bd%je + 1
    ELSE
      is = bd%is
      ie = bd%ie
      js = bd%js
      je = bd%je
    END IF
!$OMP parallel do default(none) shared(is,ie,js,je,km,grid_type,u,dx,v,dy,ua,va,a11,a12,a21,a22) &
!$OMP                          private(u1, v1, wu, wv)
    DO k=1,km
      IF (grid_type .LT. 4) THEN
        DO j=js,je+1
          DO i=is,ie
            wu(i, j) = u(i, j, k)*dx(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            wv(i, j) = v(i, j, k)*dy(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
            u1(i) = 2.*(wu(i, j)+wu(i, j+1))/(dx(i, j)+dx(i, j+1))
            v1(i) = 2.*(wv(i, j)+wv(i+1, j))/(dy(i, j)+dy(i+1, j))
!!!          u1(i) = (wu(i,j) + wu(i,j+1)) * rdxa(i,j)
!!!          v1(i) = (wv(i,j) + wv(i+1,j)) * rdya(i,j)
! Cubed (cell center co-variant winds) to lat-lon:
            ua(i, j, k) = a11(i, j)*u1(i) + a12(i, j)*v1(i)
            va(i, j, k) = a21(i, j)*u1(i) + a22(i, j)*v1(i)
          END DO
        END DO
      ELSE
! 2nd order:
        DO j=js,je
          DO i=is,ie
            ua(i, j, k) = 0.5*(u(i, j, k)+u(i, j+1, k))
            va(i, j, k) = 0.5*(v(i, j, k)+v(i+1, j, k))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE C2L_ORD2
!  Differentiation of g_sum in reverse (adjoint) mode (with options split(a2b_edge_mod.a2b_ord4 a2b_edge_mod.a2b_ord2 dyn_core
!_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dyn_cor
!e_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_mod.g
!eopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynamics_m
!od.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_dynamics_mod.geos_to_fv3 fv_dynamics_mod.fv3_to_geos fv_grid_utils_mod.cubed
!_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_en
!ergy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mod.map_scalar fv_mapz_mod.map1_ppm fv_mapz_mod.mapn_tracer fv_mapz_m
!od.map1_q2 fv_mapz_mod.remap_2d fv_mapz_mod.scalar_profile fv_mapz_mod.cs_profile fv_mapz_mod.cs_limiters fv_mapz_mod.p
!pm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.
!moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.
!tracer_2d_nested nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver_c nh_utils_
!mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver nh_utils_mo
!d.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh sw_core_mod.c
!_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_core_mod.fill3_
!4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u sw_core_mod.ytp_v sw_core_mod.compute_dive
!rgence_damping_fb sw_core_mod.smag_corner tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d tp_core_mod.copy_corners tp_core
!_mod.xppm tp_core_mod.yppm tp_core_mod.deln_flux a2b_edge_mod.extrap_corner fv_grid_utils_mod.great_circle_dist sw_c
!ore_mod.edge_interpolate4)):
!   gradient     of useful results: g_sum
!   with respect to varying inputs: p
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
  SUBROUTINE G_SUM_ADM(domain, p, p_ad, ifirst, ilast, jfirst, jlast, &
&   ngc, area, mode, reproduce, g_sum_ad)
    IMPLICIT NONE
! Fast version of globalsum
    INTEGER, INTENT(IN) :: ifirst, ilast
    INTEGER, INTENT(IN) :: jfirst, jlast, ngc
! if ==1 divided by area
    INTEGER, INTENT(IN) :: mode
    LOGICAL, INTENT(IN), OPTIONAL :: reproduce
! field to be summed
    REAL, INTENT(IN) :: p(ifirst:ilast, jfirst:jlast)
    REAL :: p_ad(ifirst:ilast, jfirst:jlast)
    REAL(kind=r_grid), INTENT(IN) :: area(ifirst-ngc:ilast+ngc, jfirst-&
&   ngc:jlast+ngc)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER :: i, j
    REAL :: gsum
    REAL :: gsum_ad
    LOGICAL, SAVE :: g_sum_initialized=.false.
    REAL(kind=r_grid), SAVE :: global_area
    REAL :: tmp(ifirst:ilast, jfirst:jlast)
    INTRINSIC PRESENT
    REAL :: temp_ad
    INTEGER :: branch
    REAL :: g_sum
    REAL :: g_sum_ad
!if ( .not. g_sum_initialized ) then
!   global_area = mpp_global_sum(domain, area, flags=BITWISE_EFP_SUM)
!   if ( is_master() ) write(*,*) 'Global Area=',global_area
!   g_sum_initialized = .true.
!end if
    global_area = 1000000.0_R_GRID
!-------------------------
! FMS global sum algorithm:
!-------------------------
    IF (PRESENT(reproduce)) THEN
!if (reproduce) then
!   gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast), &
!                         flags=BITWISE_EFP_SUM)
!else
!   gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast))
!endif
      global_area = 2000000.0_R_GRID
      gsum = 0.
      DO j=jfirst,jlast
        DO i=ifirst,ilast
          CALL PUSHREALARRAY_ADM(gsum)
          gsum = gsum*p(i, j)*area(i, j)
        END DO
      END DO
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
!call mp_reduce_sum(gsum)
    IF (mode .EQ. 1) THEN
      gsum_ad = g_sum_ad/global_area
    ELSE
      gsum_ad = g_sum_ad
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      p_ad = 0.0
      DO j=jlast,jfirst,-1
        DO i=ilast,ifirst,-1
          CALL POPREALARRAY_ADM(gsum)
          temp_ad = area(i, j)*gsum_ad
          p_ad(i, j) = p_ad(i, j) + gsum*temp_ad
          gsum_ad = p(i, j)*temp_ad
        END DO
      END DO
    ELSE
      p_ad = 0.0
      DO j=jlast,jfirst,-1
        DO i=ilast,ifirst,-1
          p_ad(i, j) = p_ad(i, j) + area(i, j)*gsum_ad
        END DO
      END DO
    END IF
  END SUBROUTINE G_SUM_ADM
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
  REAL FUNCTION G_SUM(domain, p, ifirst, ilast, jfirst, jlast, ngc, area&
&   , mode, reproduce)
    IMPLICIT NONE
! Fast version of globalsum
    INTEGER, INTENT(IN) :: ifirst, ilast
    INTEGER, INTENT(IN) :: jfirst, jlast, ngc
! if ==1 divided by area
    INTEGER, INTENT(IN) :: mode
    LOGICAL, INTENT(IN), OPTIONAL :: reproduce
! field to be summed
    REAL, INTENT(IN) :: p(ifirst:ilast, jfirst:jlast)
    REAL(kind=r_grid), INTENT(IN) :: area(ifirst-ngc:ilast+ngc, jfirst-&
&   ngc:jlast+ngc)
    TYPE(DOMAIN2D), INTENT(IN) :: domain
    INTEGER :: i, j
    REAL :: gsum
    LOGICAL, SAVE :: g_sum_initialized=.false.
    REAL(kind=r_grid), SAVE :: global_area
    REAL :: tmp(ifirst:ilast, jfirst:jlast)
    INTRINSIC PRESENT
!if ( .not. g_sum_initialized ) then
!   global_area = mpp_global_sum(domain, area, flags=BITWISE_EFP_SUM)
!   if ( is_master() ) write(*,*) 'Global Area=',global_area
!   g_sum_initialized = .true.
!end if
    global_area = 1000000.0_R_GRID
!-------------------------
! FMS global sum algorithm:
!-------------------------
    IF (PRESENT(reproduce)) THEN
!if (reproduce) then
!   gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast), &
!                         flags=BITWISE_EFP_SUM)
!else
!   gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast))
!endif
      global_area = 2000000.0_R_GRID
      gsum = 0.
      DO j=jfirst,jlast
        DO i=ifirst,ilast
          gsum = gsum*p(i, j)*area(i, j)
        END DO
      END DO
    ELSE
!-------------------------
! Quick local sum algorithm
!-------------------------
      gsum = 0.
      DO j=jfirst,jlast
        DO i=ifirst,ilast
          gsum = gsum + p(i, j)*area(i, j)
        END DO
      END DO
    END IF
!call mp_reduce_sum(gsum)
    IF (mode .EQ. 1) THEN
      g_sum = gsum/global_area
    ELSE
      g_sum = gsum
    END IF
  END FUNCTION G_SUM
END MODULE FV_GRID_UTILS_MOD_DIFF
