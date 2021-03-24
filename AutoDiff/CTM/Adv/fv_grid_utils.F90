module fv_grid_utils_mod

use fv_mp_mod, only: is,js,ie,je,isd,jsd,ied,jed
use fv_arrays_mod, only: FVPRC

implicit none
public

integer :: sw_corner, se_corner, ne_corner, nw_corner

real(FVPRC), parameter :: big_number = 1e10
real(FVPRC) :: da_min
real(FVPRC), allocatable, dimension(:,:) :: sina_u, sina_v
real(FVPRC), allocatable, dimension(:,:,:) :: sin_sg

endmodule fv_grid_utils_mod
