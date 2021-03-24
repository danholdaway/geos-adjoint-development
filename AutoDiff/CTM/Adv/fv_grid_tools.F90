module fv_grid_tools_mod

use fv_mp_mod, only: is,js,ie,je,isd,jsd,ied,jed
use fv_arrays_mod, only: FVPRC

implicit none
public

integer :: grid_type
real(FVPRC), allocatable, dimension(:,:) :: area, rarea, dxa, dya, dx, dy, rdxc, rdyc



end module fv_grid_tools_mod
