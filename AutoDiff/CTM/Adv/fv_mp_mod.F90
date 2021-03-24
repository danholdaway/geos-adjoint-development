module fv_mp_mod

use mpp_domains_mod, only: domain2D

implicit none
public

integer :: is,js,ie,je
integer :: isd,jsd,ied,jed

integer, parameter:: ng    = 3

integer :: gid

type(domain2D) :: domain

INTERFACE mp_reduce_max
  MODULE PROCEDURE mp_reduce_max_r4_1d
  MODULE PROCEDURE mp_reduce_max_r4
  MODULE PROCEDURE mp_reduce_max_r8_1d
  MODULE PROCEDURE mp_reduce_max_r8
  MODULE PROCEDURE mp_reduce_max_i4
END INTERFACE



contains


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

endmodule fv_mp_mod
