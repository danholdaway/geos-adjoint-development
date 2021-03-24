module horiz_interp_type_mod

implicit none
private

public horiz_interp_type

 type horiz_interp_type
   real,    dimension(:,:), pointer   :: faci =>NULL(), facj =>NULL()   !weights for conservative scheme
   integer, dimension(:,:), pointer   :: ilon =>NULL(), jlat =>NULL()   !indices for conservative scheme
   real,    dimension(:,:), pointer   :: area_src =>NULL()              !area of the source grid
   real,    dimension(:,:), pointer   :: area_dst =>NULL()              !area of the destination grid
   real,    dimension(:,:,:), pointer :: wti =>NULL(),wtj =>NULL()      !weights for bilinear interpolation
                                                                        !wti ist used for derivative "weights" in bicubic 
   integer, dimension(:,:,:), pointer :: i_lon =>NULL(), j_lat =>NULL() !indices for bilinear interpolation 
                                                                        !and spherical regrid
   real,    dimension(:,:,:), pointer :: src_dist =>NULL()              !distance between destination grid and 
                                                                        !neighbor source grid.
   logical, dimension(:,:), pointer   :: found_neighbors =>NULL()       !indicate whether destination grid 
                                                                        !has some source grid around it.
   real                               :: max_src_dist
   integer, dimension(:,:), pointer   :: num_found => NULL()
   integer                            :: nlon_src, nlat_src !size of source grid
   integer                            :: nlon_dst, nlat_dst !size of destination grid
   integer                            :: interp_method      !interpolation method.
                                                            !=1, conservative scheme
                                                            !=2, bilinear interpolation
                                                            !=3, spherical regrid
                                                            !=4, bicubic regrid
   real,    dimension(:,:), pointer   :: rat_x =>NULL(), rat_y =>NULL() !the ratio of coordinates of the dest grid
                                                                        ! (x_dest -x_src_r)/(x_src_l -x_src_r) and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real,    dimension(:), pointer     :: lon_in =>NULL(),  lat_in =>NULL()  !the coordinates of the source grid
   logical                            :: I_am_initialized=.false.
   integer                            :: version                            !indicate conservative interpolation version with value 1 or 2
   !--- The following are for conservative interpolation scheme version 2 ( through xgrid)
   integer                            :: nxgrid                             !number of exchange grid between src and dst grid.
   integer, dimension(:), pointer     :: i_src=>NULL(), j_src=>NULL()       !indices in source grid.
   integer, dimension(:), pointer     :: i_dst=>NULL(), j_dst=>NULL()       !indices in destination grid.
   real,    dimension(:), pointer     :: area_frac_dst=>NULL()              !area fraction in destination grid.
   real,    dimension(:,:), pointer   :: mask_in=>NULL() 
 end type

end module horiz_interp_type_mod
