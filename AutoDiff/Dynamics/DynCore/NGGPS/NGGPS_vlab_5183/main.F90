module main_mod

use fv_dynamics_mod,     only: fv_dynamics
use fv_arrays_nlm_mod,   only: fv_atmos_pert_type
use fv_arrays_mod,       only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_grid_bounds_type, fvprc
use mpp_domains_mod,     only: domain2D
use fv_sg_mod,           only: fv_subgrid_z

implicit none
private
public run

contains

subroutine run( FV_Atm, FV_AtmP, bd, npz, ncnst, DT, kappa, cp, zvir, &
                u, v, w, delz, pt, delp, q )

    type(fv_atmos_type), intent(inout) :: FV_Atm
    type(fv_atmos_pert_type), intent(inout) :: FV_AtmP

    real, intent(in) :: DT, kappa, cp, zvir
    integer, intent(in) :: npz, ncnst

    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout) :: u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)
    real, intent(inout) :: v(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst)

    real, intent(inout) :: delz(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real, intent(inout) :: w(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)

    real ::   ze0(bd%isc:bd%isc, bd%jsc: bd%jsc ,1:1)
    real ::  ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)
    real ::  pe  (bd%isc-1:bd%iec+1, npz+1,bd%jsc-1:bd%jec+1)
    real ::  pk  (bd%isc:bd%iec,bd%jsc:bd%jec, npz+1)
    real ::  peln(bd%isc:bd%iec,npz+1,bd%jsc:bd%jec)
    real ::  pkz (bd%isc:bd%iec,bd%jsc:bd%jec,npz)
    real :: q_con(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real ::  phis(bd%isd:bd%ied,bd%jsd:bd%jed)
    real ::  omga(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real ::    uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)
    real ::    vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)
    real ::    ua(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz)
    real ::    va(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz)
    real ::   mfx(bd%isc:bd%iec+1, bd%jsc:bd%jec,   npz)
    real ::   mfy(bd%isc:bd%iec  , bd%jsc:bd%jec+1, npz)
    real ::    cx(bd%isc:bd%iec+1, bd%jsd:bd%jed, npz)
    real ::    cy(bd%isd:bd%ied ,bd%jsc:bd%jec+1, npz)

    real :: u_dt(FV_Atm%bd%isd:FV_Atm%bd%ied,FV_Atm%bd%jsd:FV_Atm%bd%jed,FV_Atm%npz)
    real :: v_dt(FV_Atm%bd%isd:FV_Atm%bd%ied,FV_Atm%bd%jsd:FV_Atm%bd%jed,FV_Atm%npz)
    real :: t_dt(FV_Atm%bd%isc :FV_Atm%bd%iec ,FV_Atm%bd%jsc :FV_Atm%bd%jec ,FV_Atm%npz)

    integer :: isc,iec,jsc,jec, i,j,k

    isc = FV_Atm%bd%is
    iec = FV_Atm%bd%ie
    jsc = FV_Atm%bd%js
    jec = FV_Atm%bd%je

      ze0 = 0.0
     ps   = 0.0
     pe   = 0.0
     pk   = 0.0
     peln = 0.0
     pkz  = 0.0
    q_con = 0.0
     phis = 0.0
     omga = 0.0
       uc = 0.0
       vc = 0.0
       ua = 0.0
       va = 0.0
      mfx = 0.0
      mfy = 0.0
       cx = 0.0
       cy = 0.0

    !Compute the FV pressures
    call compute_pressures(FV_Atm%bd, FV_Atm%npz, kappa, FV_Atm%ptop, &
                           delp, pe, pk, pkz, peln)

    !Convert potential temperature to dry temperature
    do k=1,FV_Atm%npz
      do j=jsc,jec
        do i=isc,iec
           pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)
        end do
      end do
    end do

    call fv_dynamics(FV_Atm%npx, FV_Atm%npy, FV_Atm%npz, FV_Atm%ncnst, FV_Atm%ng,                      &
                               DT, FV_Atm%flagstruct%consv_te, FV_Atm%flagstruct%fill,                                    &
                               FV_Atm%flagstruct%reproduce_sum, kappa,                                                       &
                               cp, zvir, FV_Atm%ptop, FV_Atm%ks, FV_Atm%flagstruct%ncnst,                              &
                               FV_Atm%flagstruct%n_split, FV_Atm%flagstruct%q_split,                                      &
                               u, v, w,                                                           &
                               delz, FV_Atm%flagstruct%hydrostatic,                                                &
                               pt, delp,                                                                    &
                               q, ps, pe,                                                         &
                               pk, peln, pkz,                                                     &
                               phis, q_con, omga,                                                 &
                               ua, va,                                                                      &
                               uc, vc,                                                                      &
                               FV_Atm%ak, FV_Atm%bk,                                                                      &
                               mfx, mfy,                                                                    &
                               cx, cy,                                                                      &
                               ze0,                                                                                   &
                               FV_Atm%flagstruct%hybrid_z, FV_Atm%gridstruct,                                             &
                               FV_Atm%flagstruct, FV_AtmP%flagstruct, FV_Atm%neststruct,                               &
                               FV_Atm%idiag, FV_Atm%bd, FV_Atm%parent_grid, FV_Atm%domain                            )

    if ( FV_Atm%flagstruct%fv_sg_adj > 0 ) then

       u_dt(:,:,:) = 0.0
       v_dt(:,:,:) = 0.0
       t_dt(:,:,:) = 0.0

       call fv_subgrid_z(FV_Atm%bd%isd, FV_Atm%bd%ied, FV_Atm%bd%jsd, FV_Atm%bd%jed, &
                         FV_Atm%bd%isc, FV_Atm%bd%iec, FV_Atm%bd%jsc, FV_Atm%bd%jec, FV_Atm%npz, &
                         FV_Atm%ncnst, DT, FV_Atm%flagstruct%fv_sg_adj,      &
                         FV_Atm%flagstruct%nwat, delp, pe,     &
                         peln, pkz, pt, q,       &
                         ua, va, FV_Atm%flagstruct%hydrostatic,&
                         w, delz, u_dt, v_dt, t_dt, FV_Atm%flagstruct%n_zfilter)

    endif


    !Convert temperature back
    do k=1,FV_Atm%npz
      do j=jsc,jec
        do i=isc,iec
           pt(i,j,k) = pt(i,j,k)/pkz(i,j,k)
        end do
      end do
    end do

      ze0 = 0.0
     ps   = 0.0
     pe   = 0.0
     pk   = 0.0
     peln = 0.0
     pkz  = 0.0
    q_con = 0.0
     phis = 0.0
     omga = 0.0
       uc = 0.0
       vc = 0.0
       ua = 0.0
       va = 0.0
      mfx = 0.0
      mfy = 0.0
       cx = 0.0
       cy = 0.0

  end subroutine run


  subroutine compute_pressures(bd, npz, kappa, ptop, delp, pe, pk, pkz, peln)

    implicit none

    !Arguments
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: npz
    real, intent(IN) :: kappa, ptop

    real, intent(IN)    :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real, intent(OUT) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)
    real, intent(OUT) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)
    real, intent(OUT) :: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)
    real, intent(OUT) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)

    !Locals
    integer :: i, j, k, is, ie, js, je

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je

    pe(:,:,:) = 0.0
    pe(:,1,:) = ptop
    do k=2,npz+1
      do j=js,je
        do i=is,ie
          pe(i,k,j)   = pe(i, k-1, j) + delp(i,j,k-1)
        enddo
      enddo
    enddo

    do k=1,npz+1
      do j=js,je
        do i=is,ie
          peln(i,k,j) = log(pe(i,k,j))
        enddo
      enddo
    enddo

    do k=1,npz+1
      do j=js,je
        do i=is,ie
          pk(i,j,k)   = exp( kappa*peln(i,k,j) )
        enddo
      enddo
    enddo

    do k=1,npz
      do j=js,je
        do i=is,ie
          pkz(i,j,k) = (pk(i,j,k+1) - pk(i,j,k) ) / (kappa*(peln(i,k+1,j) - peln(i,k,j)) )
        end do
      end do
    end do 

  end subroutine compute_pressures



end module main_mod
