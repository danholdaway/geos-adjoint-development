  subroutine geos_to_fv3(bd, npz, kappa, ptop, delp, pe, pk, pkz, peln, pt)

    implicit none

    !Arguments
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: npz
    real, intent(IN) :: kappa, ptop

    real, intent(INOUT) :: pt  (bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real, intent(INOUT) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real, intent(INOUT) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)
    real, intent(INOUT) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)
    real, intent(INOUT) :: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)
    real, intent(INOUT) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)

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

    pt(is:ie,js:je,:) = pt(is:ie,js:je,:)*pkz(is:ie,js:je,:)

  end subroutine geos_to_fv3

  subroutine fv3_to_geos(bd, npz, pkz, pt)

    implicit none

    !Arguments
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: npz
    real, intent(INOUT) :: pt  (bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real, intent(INOUT) :: pkz (bd%is :bd%ie ,bd%js :bd%je ,npz)

    !Locals
    integer :: i, j, k, is, ie, js, je

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je

    do k=1,npz
      do j=js,je
        do i=is,ie
           pt(i,j,k) = pt(i,j,k)/pkz(i,j,k)
        end do
      end do
    end do

  end subroutine fv3_to_geos

