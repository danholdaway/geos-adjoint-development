module external_sst_mod

use fv_arrays_mod,  only: REAL4, REAL8, FVPRC

#ifdef NO_GFDL_SHARED
!----------------- Public Data -----------------------------------
integer :: i_sst = -1
integer :: j_sst = -1
logical :: forecast_mode = .false.
logical :: use_ncep_sst  = .false.
real(FVPRC), allocatable, dimension(:,:) ::  sst_ncep, sst_anom
#else
use amip_interp_mod, only: i_sst, j_sst, sst_ncep, sst_anom, &
                           forecast_mode, use_ncep_sst
#endif

public i_sst, j_sst, sst_ncep, sst_anom, forecast_mode, use_ncep_sst

!---- version number -----
character(len=128) :: version = '$Id: external_sst.F90,v 1.1.2.1.4.1 2017/02/16 03:47:48 aoloso Exp $'
character(len=128) :: tagname = '$Name: Heracles-UNSTABLE_ncepdyn_Feb222017 $'

end module external_sst_mod
