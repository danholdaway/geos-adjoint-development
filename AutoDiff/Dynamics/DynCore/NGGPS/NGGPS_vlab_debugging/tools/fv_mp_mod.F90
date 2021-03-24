module fv_mp_mod

use mpp_domains_mod, only : group_halo_update_type => mpp_group_update_type
use mpp_domains_mod, only : domain2d

implicit none
public

logical :: mp_gather !dummy
integer :: is, ie, js, je, isd, ied, jsd, jed, isc, iec, jsc, jec
integer :: mpp_broadcast
integer :: mp_bcst !dummy

integer, parameter :: XDir=1
integer, parameter :: YDir=2
integer, parameter:: ng    = 3

interface start_group_halo_update
  module procedure start_var_group_update_2d
  module procedure start_var_group_update_3d
  module procedure start_var_group_update_4d
  module procedure start_vector_group_update_2d
  module procedure start_vector_group_update_3d
end interface start_group_halo_update

!interface complete_group_halo_update
!  module procedure complete_var_group_update_2d
!  module procedure complete_var_group_update_3d
!  module procedure complete_var_group_update_4d
!  module procedure complete_vector_group_update_2d
!  module procedure complete_vector_group_update_3d
!end interface complete_group_halo_update

INTERFACE mp_reduce_max
  MODULE PROCEDURE mp_reduce_max_r4_1d
  MODULE PROCEDURE mp_reduce_max_r4
  MODULE PROCEDURE mp_reduce_max_r8_1d
  MODULE PROCEDURE mp_reduce_max_r8
  MODULE PROCEDURE mp_reduce_max_i4
END INTERFACE

INTERFACE mp_reduce_sum
  MODULE PROCEDURE mp_reduce_sum_r4
  MODULE PROCEDURE mp_reduce_sum_r4_1d
  MODULE PROCEDURE mp_reduce_sum_r8
  MODULE PROCEDURE mp_reduce_sum_r8_1d
END INTERFACE

INTERFACE fill_corners
  MODULE PROCEDURE fill_corners_2d_r4
  MODULE PROCEDURE fill_corners_2d_r8
  MODULE PROCEDURE fill_corners_xy_2d_r4
  MODULE PROCEDURE fill_corners_xy_2d_r8
  MODULE PROCEDURE fill_corners_xy_3d_r4
  MODULE PROCEDURE fill_corners_xy_3d_r8
END INTERFACE

INTERFACE fill_corners_agrid
  MODULE PROCEDURE fill_corners_agrid_r4
  MODULE PROCEDURE fill_corners_agrid_r8
END INTERFACE

INTERFACE fill_corners_cgrid
  MODULE PROCEDURE fill_corners_cgrid_r4
  MODULE PROCEDURE fill_corners_cgrid_r8
END INTERFACE

INTERFACE fill_corners_dgrid
  MODULE PROCEDURE fill_corners_dgrid_r4
  MODULE PROCEDURE fill_corners_dgrid_r8
END INTERFACE


contains

logical function is_master()

if (1==1) then
   is_master = .true.
endif

end function

! FILL CORNERS
! ------------

      subroutine fill_corners_2d_r4(q, npx, npy, FILL, AGRID, BGRID)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: q
         integer, intent(IN):: npx,npy
         integer, intent(IN):: FILL  ! X-Dir or Y-Dir 
         logical, OPTIONAL, intent(IN) :: AGRID, BGRID 
         integer :: i,j

         if (present(BGRID)) then
            if (BGRID) then
              select case (FILL)
              case (XDir)
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-i  ,1-j  ) = q(1-j  ,i+1    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-i  ,npy+j) = q(1-j  ,npy-i  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+i,1-j  ) = q(npx+j,i+1    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+i,npy+j) = q(npx+j,npy-i  )  !NE Corner
                    enddo
                 enddo
              case (YDir)
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-j  ,1-i  ) = q(i+1  ,1-j    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-j  ,npy+i) = q(i+1  ,npy+j  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+j,1-i  ) = q(npx-i,1-j    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+j,npy+i) = q(npx-i,npy+j  )  !NE Corner
                    enddo
                 enddo
              case default
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-i  ,1-j  ) = q(1-j  ,i+1    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-i  ,npy+j) = q(1-j  ,npy-i  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+i,1-j  ) = q(npx+j,i+1    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+i,npy+j) = q(npx+j,npy-i  )  !NE Corner
                    enddo
                 enddo
              end select
            endif
          elseif (present(AGRID)) then
            if (AGRID) then
              select case (FILL)
              case (XDir)
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) q(1-i    ,1-j    ) = q(1-j    ,i        )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-i    ,npy-1+j) = q(1-j    ,npy-1-i+1)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+i,1-j    ) = q(npx-1+j,i        )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+i,npy-1+j) = q(npx-1+j,npy-1-i+1)  !NE Corner
                    enddo
                 enddo
              case (YDir)
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) q(1-j    ,1-i    ) = q(i        ,1-j    )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-j    ,npy-1+i) = q(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+j,1-i    ) = q(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+j,npy-1+i) = q(npx-1-i+1,npy-1+j)  !NE Corner
                    enddo
                 enddo
              case default
                 do j=1,ng
                    do i=1,ng        
                       if ((is==    1) .and. (js==    1)) q(1-j    ,1-i    ) = q(i        ,1-j    )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-j    ,npy-1+i) = q(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+j,1-i    ) = q(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+j,npy-1+i) = q(npx-1-i+1,npy-1+j)  !NE Corner
                   enddo
                 enddo          
              end select
            endif
          endif

      end subroutine fill_corners_2d_r4
      subroutine fill_corners_2d_r8(q, npx, npy, FILL, AGRID, BGRID)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: q
         integer, intent(IN):: npx,npy
         integer, intent(IN):: FILL  ! X-Dir or Y-Dir 
         logical, OPTIONAL, intent(IN) :: AGRID, BGRID 
         integer :: i,j

         if (present(BGRID)) then
            if (BGRID) then
              select case (FILL)
              case (XDir)
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-i  ,1-j  ) = q(1-j  ,i+1    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-i  ,npy+j) = q(1-j  ,npy-i  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+i,1-j  ) = q(npx+j,i+1    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+i,npy+j) = q(npx+j,npy-i  )  !NE Corner
                    enddo
                 enddo
              case (YDir)
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-j  ,1-i  ) = q(i+1  ,1-j    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-j  ,npy+i) = q(i+1  ,npy+j  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+j,1-i  ) = q(npx-i,1-j    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+j,npy+i) = q(npx-i,npy+j  )  !NE Corner
                    enddo
                 enddo
              case default
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-i  ,1-j  ) = q(1-j  ,i+1    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-i  ,npy+j) = q(1-j  ,npy-i  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+i,1-j  ) = q(npx+j,i+1    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+i,npy+j) = q(npx+j,npy-i  )  !NE Corner
                    enddo
                 enddo
              end select
            endif
          elseif (present(AGRID)) then
            if (AGRID) then
              select case (FILL)
              case (XDir)
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) q(1-i    ,1-j    ) = q(1-j    ,i        )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-i    ,npy-1+j) = q(1-j    ,npy-1-i+1)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+i,1-j    ) = q(npx-1+j,i        )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+i,npy-1+j) = q(npx-1+j,npy-1-i+1)  !NE Corner
                    enddo
                 enddo
              case (YDir)
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) q(1-j    ,1-i    ) = q(i        ,1-j    )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-j    ,npy-1+i) = q(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+j,1-i    ) = q(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+j,npy-1+i) = q(npx-1-i+1,npy-1+j)  !NE Corner
                    enddo
                 enddo
              case default
                 do j=1,ng
                    do i=1,ng        
                       if ((is==    1) .and. (js==    1)) q(1-j    ,1-i    ) = q(i        ,1-j    )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-j    ,npy-1+i) = q(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+j,1-i    ) = q(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+j,npy-1+i) = q(npx-1-i+1,npy-1+j)  !NE Corner
                   enddo
                 enddo          
              end select
            endif
          endif

      end subroutine fill_corners_2d_r8
      subroutine fill_corners_xy_2d_r4(x, y, npx, npy, DGRID, AGRID, CGRID, VECTOR)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: x !(isd:ied  ,jsd:jed+1)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: y !(isd:ied+1,jsd:jed  )
         integer, intent(IN):: npx,npy
         logical, OPTIONAL, intent(IN) :: DGRID, AGRID, CGRID, VECTOR
         integer :: i,j

         real(kind=4) :: mySign

         mySign = 1.0
         if (present(VECTOR)) then
            if (VECTOR) mySign = -1.0
         endif

         if (present(DGRID)) then
            call fill_corners_dgrid(x, y, npx, npy, mySign)
         elseif (present(CGRID)) then
            call fill_corners_cgrid(x, y, npx, npy, mySign)
         elseif (present(AGRID)) then
            call fill_corners_agrid(x, y, npx, npy, mySign)
         else
            call fill_corners_agrid(x, y, npx, npy, mySign)
         endif

      end subroutine fill_corners_xy_2d_r4
      subroutine fill_corners_xy_2d_r8(x, y, npx, npy, DGRID, AGRID, CGRID, VECTOR)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: x !(isd:ied  ,jsd:jed+1)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: y !(isd:ied+1,jsd:jed  )
         integer, intent(IN):: npx,npy
         logical, OPTIONAL, intent(IN) :: DGRID, AGRID, CGRID, VECTOR
         integer :: i,j

         real(kind=8) :: mySign

         mySign = 1.0
         if (present(VECTOR)) then
            if (VECTOR) mySign = -1.0
         endif

         if (present(DGRID)) then
            call fill_corners_dgrid(x, y, npx, npy, mySign)
         elseif (present(CGRID)) then
            call fill_corners_cgrid(x, y, npx, npy, mySign)
         elseif (present(AGRID)) then
            call fill_corners_agrid(x, y, npx, npy, mySign)
         else
            call fill_corners_agrid(x, y, npx, npy, mySign)
         endif

      end subroutine fill_corners_xy_2d_r8
      subroutine fill_corners_xy_3d_r4(x, y, npx, npy, npz, DGRID, AGRID, CGRID, VECTOR)
         real(kind=4), DIMENSION(isd:,jsd:,:), intent(INOUT):: x !(isd:ied  ,jsd:jed+1)
         real(kind=4), DIMENSION(isd:,jsd:,:), intent(INOUT):: y !(isd:ied+1,jsd:jed  )
         integer, intent(IN):: npx,npy,npz
         logical, OPTIONAL, intent(IN) :: DGRID, AGRID, CGRID, VECTOR
         integer :: i,j,k

         real(kind=4) :: mySign

         mySign = 1.0
         if (present(VECTOR)) then
            if (VECTOR) mySign = -1.0
         endif

         if (present(DGRID)) then
            do k=1,npz
               call fill_corners_dgrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         elseif (present(CGRID)) then
            do k=1,npz
               call fill_corners_cgrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         elseif (present(AGRID)) then
            do k=1,npz
               call fill_corners_agrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         else
            do k=1,npz
               call fill_corners_agrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         endif

      end subroutine fill_corners_xy_3d_r4
      subroutine fill_corners_xy_3d_r8(x, y, npx, npy, npz, DGRID, AGRID, CGRID, VECTOR)
         real(kind=8), DIMENSION(isd:,jsd:,:), intent(INOUT):: x !(isd:ied  ,jsd:jed+1)
         real(kind=8), DIMENSION(isd:,jsd:,:), intent(INOUT):: y !(isd:ied+1,jsd:jed  )
         integer, intent(IN):: npx,npy,npz
         logical, OPTIONAL, intent(IN) :: DGRID, AGRID, CGRID, VECTOR
         integer :: i,j,k

         real(kind=8) :: mySign

         mySign = 1.0
         if (present(VECTOR)) then
            if (VECTOR) mySign = -1.0
         endif

         if (present(DGRID)) then
            do k=1,npz
               call fill_corners_dgrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         elseif (present(CGRID)) then
            do k=1,npz
               call fill_corners_cgrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         elseif (present(AGRID)) then
            do k=1,npz
               call fill_corners_agrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         else
            do k=1,npz
               call fill_corners_agrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         endif

      end subroutine fill_corners_xy_3d_r8


! FILL CORNERS AGRID
! ------------------
      subroutine fill_corners_agrid_r4(x, y, npx, npy, mySign)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=4), intent(IN) :: mySign
         integer :: i,j

                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) x(1-i    ,1-j    ) = mySign*y(1-j    ,i        )  !SW Corner
                       if ((is==    1) .and. (je==npy-1)) x(1-i    ,npy-1+j) =        y(1-j    ,npy-1-i+1)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) x(npx-1+i,1-j    ) =        y(npx-1+j,i        )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) x(npx-1+i,npy-1+j) = mySign*y(npx-1+j,npy-1-i+1)  !NE Corner
                    enddo
                 enddo
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) y(1-j    ,1-i    ) = mySign*x(i        ,1-j    )  !SW Corner
                       if ((is==    1) .and. (je==npy-1)) y(1-j    ,npy-1+i) =        x(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) y(npx-1+j,1-i    ) =        x(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) y(npx-1+j,npy-1+i) = mySign*x(npx-1-i+1,npy-1+j)  !NE Corner
                    enddo
                 enddo

      end subroutine fill_corners_agrid_r4
      subroutine fill_corners_agrid_r8(x, y, npx, npy, mySign)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=8), intent(IN) :: mySign
         integer :: i,j

                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) x(1-i    ,1-j    ) = mySign*y(1-j    ,i        )  !SW Corner
                       if ((is==    1) .and. (je==npy-1)) x(1-i    ,npy-1+j) =        y(1-j    ,npy-1-i+1)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) x(npx-1+i,1-j    ) =        y(npx-1+j,i        )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) x(npx-1+i,npy-1+j) = mySign*y(npx-1+j,npy-1-i+1)  !NE Corner
                    enddo
                 enddo
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) y(1-j    ,1-i    ) = mySign*x(i        ,1-j    )  !SW Corner
                       if ((is==    1) .and. (je==npy-1)) y(1-j    ,npy-1+i) =        x(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) y(npx-1+j,1-i    ) =        x(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) y(npx-1+j,npy-1+i) = mySign*x(npx-1-i+1,npy-1+j)  !NE Corner
                    enddo
                 enddo

      end subroutine fill_corners_agrid_r8

! FILL CORNERS CGRID
! ------------------
      subroutine fill_corners_cgrid_r4(x, y, npx, npy, mySign)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=4), intent(IN) :: mySign
         integer :: i,j

                  do j=1,ng
                     do i=1,ng
                        if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j    ) =        y(j      ,1-i  )  !SW Corner 
                        if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy-1+j) = mySign*y(j      ,npy+i)  !NW Corner
                        if ((ie+1==npx) .and. (js  ==  1)) x(npx+i  ,1-j    ) = mySign*y(npx-j  ,1-i  )  !SE Corner
                        if ((ie+1==npx) .and. (je+1==npy)) x(npx+i  ,npy-1+j) =        y(npx-j  ,npy+i)  !NE Corner
                     enddo
                  enddo
                  do j=1,ng
                     do i=1,ng
                        if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j  ) =        x(1-j  ,i    )  !SW Corner 
                        if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy+j) = mySign*x(1-j  ,npy-i)  !NW Corner
                        if ((ie+1==npx) .and. (js  ==  1)) y(npx-1+i,1-j  ) = mySign*x(npx+j,i    )  !SE Corner
                        if ((ie+1==npx) .and. (je+1==npy)) y(npx-1+i,npy+j) =        x(npx+j,npy-i)  !NE Corner
                     enddo
                  enddo
      
      end subroutine fill_corners_cgrid_r4
      subroutine fill_corners_cgrid_r8(x, y, npx, npy, mySign)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=8), intent(IN) :: mySign
         integer :: i,j

                  do j=1,ng
                     do i=1,ng
                        if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j    ) =        y(j      ,1-i  )  !SW Corner 
                        if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy-1+j) = mySign*y(j      ,npy+i)  !NW Corner
                        if ((ie+1==npx) .and. (js  ==  1)) x(npx+i  ,1-j    ) = mySign*y(npx-j  ,1-i  )  !SE Corner
                        if ((ie+1==npx) .and. (je+1==npy)) x(npx+i  ,npy-1+j) =        y(npx-j  ,npy+i)  !NE Corner
                     enddo
                  enddo
                  do j=1,ng
                     do i=1,ng
                        if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j  ) =        x(1-j  ,i    )  !SW Corner 
                        if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy+j) = mySign*x(1-j  ,npy-i)  !NW Corner
                        if ((ie+1==npx) .and. (js  ==  1)) y(npx-1+i,1-j  ) = mySign*x(npx+j,i    )  !SE Corner
                        if ((ie+1==npx) .and. (je+1==npy)) y(npx-1+i,npy+j) =        x(npx+j,npy-i)  !NE Corner
                     enddo
                  enddo
      
      end subroutine fill_corners_cgrid_r8

! FILL CORNERS DGRID
! ------------------
      subroutine fill_corners_dgrid_r4(x, y, npx, npy, mySign)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=4), intent(IN) :: mySign 
         integer :: i,j

               do j=1,ng
                  do i=1,ng
                   !   if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) =        y(j+1  ,1-i    )  !SW Corner 
                   !   if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) = mySign*y(j+1  ,npy-1+i)  !NW Corner
                   !   if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) = mySign*y(npx-j,1-i    )  !SE Corner
                   !   if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) =        y(npx-j,npy-1+i)  !NE Corner
                      if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) = mySign*y(1-j  ,i    )  !SW Corner 
                      if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) =        y(1-j  ,npy-i)  !NW Corner
                      if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) =        y(npx+j,i    )  !SE Corner
                      if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) = mySign*y(npx+j,npy-i)  !NE Corner
                  enddo
               enddo
               do j=1,ng
                  do i=1,ng
                   !  if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) =        x(1-j    ,i+1  )  !SW Corner 
                   !  if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) = mySign*x(1-j    ,npy-i)  !NW Corner
                   !  if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) = mySign*x(npx-1+j,i+1  )  !SE Corner
                   !  if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) =        x(npx-1+j,npy-i)  !NE Corner
                     if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) = mySign*x(j      ,1-i  )  !SW Corner 
                     if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) =        x(j      ,npy+i)  !NW Corner
                     if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) =        x(npx-j  ,1-i  )  !SE Corner
                     if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) = mySign*x(npx-j  ,npy+i)  !NE Corner
                  enddo
               enddo

      end subroutine fill_corners_dgrid_r4
      subroutine fill_corners_dgrid_r8(x, y, npx, npy, mySign)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=8), intent(IN) :: mySign 
         integer :: i,j

               do j=1,ng
                  do i=1,ng
                   !   if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) =        y(j+1  ,1-i    )  !SW Corner 
                   !   if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) = mySign*y(j+1  ,npy-1+i)  !NW Corner
                   !   if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) = mySign*y(npx-j,1-i    )  !SE Corner
                   !   if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) =        y(npx-j,npy-1+i)  !NE Corner
                      if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) = mySign*y(1-j  ,i    )  !SW Corner 
                      if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) =        y(1-j  ,npy-i)  !NW Corner
                      if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) =        y(npx+j,i    )  !SE Corner
                      if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) = mySign*y(npx+j,npy-i)  !NE Corner
                  enddo
               enddo
               do j=1,ng
                  do i=1,ng
                   !  if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) =        x(1-j    ,i+1  )  !SW Corner 
                   !  if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) = mySign*x(1-j    ,npy-i)  !NW Corner
                   !  if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) = mySign*x(npx-1+j,i+1  )  !SE Corner
                   !  if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) =        x(npx-1+j,npy-i)  !NE Corner
                     if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) = mySign*x(j      ,1-i  )  !SW Corner 
                     if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) =        x(j      ,npy+i)  !NW Corner
                     if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) =        x(npx-j  ,1-i  )  !SE Corner
                     if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) = mySign*x(npx-j  ,npy+i)  !NE Corner
                  enddo
               enddo

      end subroutine fill_corners_dgrid_r8

! MP REDUCE SUM
! -------------

      subroutine mp_reduce_sum_r4(mysum)
         real(kind=4), intent(INOUT)  :: mysum

         real(kind=4) :: gsum

         mysum = 2*mysum

      end subroutine mp_reduce_sum_r4
      subroutine mp_reduce_sum_r8(mysum)
         real(kind=8), intent(INOUT)  :: mysum

         real(kind=8) :: gsum

         mysum = 2*mysum

      end subroutine mp_reduce_sum_r8
      subroutine mp_reduce_sum_r4_1d(mysum, sum1d, npts)
         integer, intent(in)  :: npts
         real(kind=4), intent(in)     :: sum1d(npts)
         real(kind=4), intent(INOUT)  :: mysum

         real(kind=4) :: gsum
         integer :: i

         mysum = 2*mysum

      end subroutine mp_reduce_sum_r4_1d
      subroutine mp_reduce_sum_r8_1d(mysum, sum1d, npts)
         integer, intent(in)  :: npts
         real(kind=8), intent(in)     :: sum1d(npts)
         real(kind=8), intent(INOUT)  :: mysum

         real(kind=8) :: gsum
         integer :: i

         mysum = 2*mysum

      end subroutine mp_reduce_sum_r8_1d



! MP REDUCE MAX
! -------------

      subroutine mp_reduce_max_r4_1d(mymax,npts)
         integer, intent(IN)  :: npts
         real(kind=4), intent(INOUT)  :: mymax(npts)
        
         real(kind=4) :: gmax(npts)
        
         mymax = 2*mymax
        
      end subroutine mp_reduce_max_r4_1d
      subroutine mp_reduce_max_r8_1d(mymax,npts)
         integer, intent(IN)  :: npts
         real(kind=8), intent(INOUT)  :: mymax(npts)
        
         real(kind=8) :: gmax(npts)
        
         mymax = 2*mymax
        
      end subroutine mp_reduce_max_r8_1d
      subroutine mp_reduce_max_r4(mymax)
         real(kind=4), intent(INOUT)  :: mymax

         real(kind=4) :: gmax
        
         mymax = 2*mymax

      end subroutine mp_reduce_max_r4
      subroutine mp_reduce_max_r8(mymax)
         real(kind=8), intent(INOUT)  :: mymax

         real(kind=8) :: gmax
        
         mymax = 2*mymax

      end subroutine mp_reduce_max_r8
      subroutine mp_reduce_max_i4(mymax)
         integer, intent(INOUT)  :: mymax

         integer :: gmax
        
         mymax = 2*mymax

      end subroutine mp_reduce_max_i4


! start_group_halo_update
! -----------------------

      subroutine start_var_group_update_2d(group, array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
        type(group_halo_update_type), intent(inout) :: group!, group_tl
        real, dimension(:,:),         intent(inout) :: array
        type(domain2D),               intent(inout) :: domain
        integer,      optional,       intent(in)    :: flags
        integer,      optional,       intent(in)    :: position
        integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
        logical,      optional,       intent(in)    :: complete
        real                                        :: d_type
        logical                                     :: is_complete
      ! Arguments: 
      !  (inout)   group - The data type that store information for group update. 
      !                    This data will be used in do_group_pass.
      !  (inout)   array - The array which is having its halos points exchanged.
      !  (in)      domain - contains domain information.
      !  (in)      flags  - An optional integer indicating which directions the
      !                       data should be sent.  
      !  (in)      position - An optional argument indicating the position.  This is
      !                       may be CORNER, but is CENTER by default.
      !  (in)      complete - An optional argument indicating whether the halo updates
      !                       should be initiated immediately or wait for second 
      !                       pass_..._start call.  Omitting complete is the same as 
      !                       setting complete to .true.
      
        array = 2*array
      
      end subroutine start_var_group_update_2d
      
      
      subroutine start_var_group_update_3d(group, array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
        type(group_halo_update_type), intent(inout) :: group!, group_tl
        real, dimension(:,:,:),       intent(inout) :: array
        type(domain2D),               intent(inout) :: domain
        integer,           optional,  intent(in)    :: flags
        integer,           optional,  intent(in)    :: position
        integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
        logical,      optional,       intent(in)    :: complete
        real                                        :: d_type
        logical                                     :: is_complete
      
      ! Arguments: 
      !  (inout)   group - The data type that store information for group update. 
      !                    This data will be used in do_group_pass.
      !  (inout)   array - The array which is having its halos points exchanged.
      !  (in)      domain - contains domain information.
      !  (in)      flags  - An optional integer indicating which directions the
      !                       data should be sent.  
      !  (in)      position - An optional argument indicating the position.  This is
      !                       may be CORNER, but is CENTER by default.
      !  (in)      complete - An optional argument indicating whether the halo updates
      !                       should be initiated immediately or wait for second 
      !                       pass_..._start call.  Omitting complete is the same as 
      !                       setting complete to .true.
      
        array = 2*array
      
      end subroutine start_var_group_update_3d
      
      subroutine start_var_group_update_4d(group, array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
        type(group_halo_update_type), intent(inout) :: group!, group_tl
        real, dimension(:,:,:,:),     intent(inout) :: array
        type(domain2D),               intent(inout) :: domain
        integer,           optional,  intent(in)    :: flags
        integer,           optional,  intent(in)    :: position
        integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
        logical,      optional,       intent(in)    :: complete
        real                                        :: d_type
        logical                                     :: is_complete
      
      ! Arguments: 
      !  (inout)   group - The data type that store information for group update. 
      !                    This data will be used in do_group_pass.
      !  (inout)   array - The array which is having its halos points exchanged.
      !  (in)      domain - contains domain information.
      !  (in)      flags  - An optional integer indicating which directions the
      !                       data should be sent.  
      !  (in)      position - An optional argument indicating the position.  This is
      !                       may be CORNER, but is CENTER by default.
      !  (in)      complete - An optional argument indicating whether the halo updates
      !                       should be initiated immediately or wait for second 
      !                       pass_..._start call.  Omitting complete is the same as 
      !                       setting complete to .true.
      
        integer :: dirflag
      
        array = 2*array
      
      end subroutine start_var_group_update_4d
      
      
      
      subroutine start_vector_group_update_2d(group, u_cmpt, v_cmpt, domain, flags, gridtype, whalo, ehalo, &
                                               shalo, nhalo, complete)
        type(group_halo_update_type), intent(inout) :: group!, group_tl
        real,       dimension(:,:),   intent(inout) :: u_cmpt, v_cmpt
        type(domain2d),               intent(inout) :: domain
        integer,            optional, intent(in)    :: flags
        integer,            optional, intent(in)    :: gridtype
        integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
        logical,      optional,       intent(in)    :: complete
        real                                        :: d_type
        logical                                     :: is_complete
      
      ! Arguments: 
      !  (inout)   group - The data type that store information for group update. 
      !                    This data will be used in do_group_pass.
      !  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
      !                     is having its halos points exchanged.
      !  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
      !                     which is having its halos points exchanged. 
      !  (in)      domain - Contains domain decomposition information.
      !  (in)      flags - An optional integer indicating which directions the
      !                        data should be sent. 
      !  (in)      gridtype - An optional flag, which may be one of A_GRID, BGRID_NE,
      !                      CGRID_NE or DGRID_NE, indicating where the two components of the
      !                      vector are discretized. 
      !  (in)      complete - An optional argument indicating whether the halo updates
      !                       should be initiated immediately or wait for second 
      !                       pass_..._start call.  Omitting complete is the same as 
      !                       setting complete to .true.
      
        u_cmpt = 2*u_cmpt
        v_cmpt = 2*v_cmpt
      
      end subroutine start_vector_group_update_2d
      
      subroutine start_vector_group_update_3d(group, u_cmpt, v_cmpt, domain, flags, gridtype, &
                                              whalo, ehalo, shalo, nhalo, complete)
        type(group_halo_update_type), intent(inout) :: group!, group_tl
        real,       dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt
        type(domain2d),               intent(inout) :: domain
        integer,            optional, intent(in)    :: flags
        integer,            optional, intent(in)    :: gridtype
        integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
        logical,      optional,       intent(in)    :: complete
        real                                        :: d_type
        logical                                     :: is_complete
      
      ! Arguments: 
      !  (inout)   group - The data type that store information for group update. 
      !                    This data will be used in do_group_pass.
      !  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
      !                     is having its halos points exchanged.
      !  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
      !                     which is having its halos points exchanged. 
      !  (in)      domain - Contains domain decomposition information.
      !  (in)      flags - An optional integer indicating which directions the
      !                        data should be sent. 
      !  (in)      gridtype - An optional flag, which may be one of A_GRID, BGRID_NE,
      !                      CGRID_NE or DGRID_NE, indicating where the two components of the
      !                      vector are discretized. 
      !  (in)      complete - An optional argument indicating whether the halo updates
      !                       should be initiated immediately or wait for second 
      !                       pass_..._start call.  Omitting complete is the same as 
      !                       setting complete to .true.
      
        u_cmpt = 2*u_cmpt
        v_cmpt = 2*v_cmpt
      
      end subroutine start_vector_group_update_3d


! complete_group_halo_update
! --------------------------

      subroutine complete_group_halo_update(group, domain)
        type(group_halo_update_type), intent(inout) :: group!, group_tl
        !real, dimension(:,:),         intent(inout) :: array
        type(domain2D),               intent(inout) :: domain
      ! Arguments: 
      !  (inout)   group - The data type that store information for group update. 
      !                    This data will be used in do_group_pass.
      !  (inout)   array - The array which is having its halos points exchanged.
      !  (in)      domain - contains domain information.
      
        !array = 2*array
      
      end subroutine complete_group_halo_update
      
      
!      subroutine complete_var_group_update_3d(group, group_tl, domain)
!        type(group_halo_update_type), intent(inout) :: group, group_tl
!        !real, dimension(:,:,:),       intent(inout) :: array
!        type(domain2D),               intent(inout) :: domain
!      
!      ! Arguments: 
!      !  (inout)   group - The data type that store information for group update. 
!      !                    This data will be used in do_group_pass.
!      !  (inout)   array - The array which is having its halos points exchanged.
!      !  (in)      domain - contains domain information.
!      
!        !array = 2*array
!      
!      end subroutine complete_var_group_update_3d
!      
!      subroutine complete_var_group_update_4d(group, group_tl, domain)
!        type(group_halo_update_type), intent(inout) :: group, group_tl
!        !real, dimension(:,:,:,:),     intent(inout) :: array
!        type(domain2D),               intent(inout) :: domain
!      
!      ! Arguments: 
!      !  (inout)   group - The data type that store information for group update. 
!      !                    This data will be used in do_group_pass.
!      !  (inout)   array - The array which is having its halos points exchanged.
!      !  (in)      domain - contains domain information.
!      
!        integer :: dirflag
!      
!        !array = 2*array
!      
!      end subroutine complete_var_group_update_4d
!      
!      
!      
!      subroutine complete_vector_group_update_2d(group, group_tl, domain)
!        type(group_halo_update_type), intent(inout) :: group, group_tl
!        !real,       dimension(:,:),   intent(inout) :: u_cmpt, v_cmpt
!        type(domain2d),               intent(inout) :: domain
!      
!      ! Arguments: 
!      !  (inout)   group - The data type that store information for group update. 
!      !                    This data will be used in do_group_pass.
!      !  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!      !                     is having its halos points exchanged.
!      !  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!      !                     which is having its halos points exchanged. 
!      !  (in)      domain - Contains domain decomposition information.
!      
!        !u_cmpt = 2*u_cmpt
!        !v_cmpt = 2*v_cmpt
!      
!      end subroutine complete_vector_group_update_2d
!      
!      subroutine complete_vector_group_update_3d(group, group_tl, domain)
!        type(group_halo_update_type), intent(inout) :: group, group_tl
!        !real,       dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt
!        type(domain2d),               intent(inout) :: domain
!
!      
!      ! Arguments: 
!      !  (inout)   group - The data type that store information for group update. 
!      !                    This data will be used in do_group_pass.
!      !  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!      !                     is having its halos points exchanged.
!      !  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!      !                     which is having its halos points exchanged. 
!      !  (in)      domain - Contains domain decomposition information.
!      
!        !u_cmpt = 2*u_cmpt
!        !v_cmpt = 2*v_cmpt
!      
!      end subroutine complete_vector_group_update_3d



end module fv_mp_mod
