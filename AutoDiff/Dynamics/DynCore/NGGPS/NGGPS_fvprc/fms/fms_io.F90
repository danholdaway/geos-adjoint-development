module fms_io_mod
 use mpp_domains_mod
 use mpp_io_mod

implicit none
private

public restart_file_type

integer, parameter, private :: max_split_file = 50
integer, parameter, private :: max_fields=400
integer, parameter, private :: max_axes=40
integer, parameter, private :: max_atts=20
integer, parameter, private :: max_domains = 10
integer, parameter, private :: MAX_TIME_LEVEL_REGISTER = 2
integer, parameter, private :: MAX_TIME_LEVEL_WRITE = 20
integer, parameter          :: max_axis_size=10000

! Index postions for axes in restart_file_type
! This is done so the user may define the axes
! in any order but a check can be performed
! to ensure no registration of duplicate axis
integer, parameter, private :: XIDX=1
integer, parameter, private :: YIDX=2
integer, parameter, private :: CIDX=3
integer, parameter, private :: ZIDX=4
integer, parameter, private :: HIDX=5
integer, parameter, private :: NIDX=5

type ax_type
   private
   character(len=128) :: name
   character(len=128) :: units
   character(len=128) :: longname
   character(len=8)   :: cartesian
   character(len=256) :: compressed
   character(len=128) :: dimlen_name
   character(len=128) :: dimlen_lname
   integer            :: sense              !Orientation of z axis definition
   integer            :: dimlen             !max dim of elements across global domain
   real               :: min             !valid min for real axis data
   integer            :: imin            !valid min for integer axis data
   integer,allocatable :: idx(:)         !compressed io-domain index vector
   integer,allocatable :: nelems(:)      !num elements for each rank in io domain
   real, pointer      :: data(:) =>NULL()    !real axis values (not used if time axis)
   type(domain2d),pointer :: domain      !domain associated with compressed axis
end type ax_type

type var_type
   private
   character(len=128)                     :: name
   character(len=128)                     :: longname
   character(len=128)                     :: units
   real, dimension(:,:,:,:), ALLOCATABLE :: buffer !_NULL
   logical                                :: domain_present
   logical                                :: write_on_this_pe
   integer                                :: domain_idx
   logical                                :: is_dimvar
   type(fieldtype)                        :: field
   type(axistype)                         :: axis
   integer                                :: position
   integer                                :: ndim
   integer                                :: siz(4)      ! X/Y/Z/T extent of fields (data domain
                                                         ! size for distributed writes;global size for reads)
   integer                                :: gsiz(4)     ! global X/Y/Z/T extent of fields
   integer                                :: csiz(4)     ! actual data size in the file
   integer                                :: id_axes(3)  ! store index for x/y/z axistype.
   logical                                :: initialized ! indicate if the field is read or not in routine save_state.
   logical                                :: mandatory   ! indicate if the field is mandatory to be when restart.
   integer                                :: is, ie, js, je  ! index of the data in compute domain
   real                                   :: default_data 
   character(len=8)                       :: compressed_axis !< If on a compressed axis, which axis
   integer, dimension(:), allocatable     :: pelist
   integer                                :: ishift, jshift ! can be used to shift indices when no_domain=T
   integer                                :: x_halo, y_halo ! can be used to indicate halo size when no_domain=T
end type var_type

type Ptr0Dr
   real,                   pointer :: p => NULL()
end type Ptr0Dr

type Ptr1Dr
   real, dimension(:),     pointer :: p => NULL()
end type Ptr1Dr

type Ptr2Dr
   real, dimension(:,:),   pointer :: p => NULL()
end type Ptr2Dr

type Ptr3Dr
   real, dimension(:,:,:), pointer :: p => NULL()
end type Ptr3Dr

type Ptr0Di
   integer,                   pointer :: p => NULL()
end type Ptr0Di

type Ptr1Di
   integer, dimension(:),     pointer :: p => NULL()
end type Ptr1Di

type Ptr2Di
   integer, dimension(:,:),   pointer :: p => NULL()
end type Ptr2Di

type Ptr3Di
   integer, dimension(:,:,:), pointer :: p => NULL()
end type Ptr3Di

type restart_file_type
   private
   integer                                  :: unit = -1 ! mpp_io unit for netcdf file
   character(len=128)                       :: name
   integer                                  :: nvar, natt, max_ntime
   logical                                  :: is_root_pe
   logical                                  :: is_compressed
   integer                                  :: tile_count
   type(ax_type),  allocatable              :: axes(:)  ! Currently define only X,Y,Compressed and maybe Z
   type(var_type), dimension(:),   pointer  :: var  => NULL()
   type(Ptr0Dr),   dimension(:,:), pointer  :: p0dr => NULL()
   type(Ptr1Dr),   dimension(:,:), pointer  :: p1dr => NULL()
   type(Ptr2Dr),   dimension(:,:), pointer  :: p2dr => NULL()
   type(Ptr3Dr),   dimension(:,:), pointer  :: p3dr => NULL()
   type(Ptr0Di),   dimension(:,:), pointer  :: p0di => NULL()
   type(Ptr1Di),   dimension(:,:), pointer  :: p1di => NULL()
   type(Ptr2Di),   dimension(:,:), pointer  :: p2di => NULL()
   type(Ptr3Di),   dimension(:,:), pointer  :: p3di => NULL()
end type restart_file_type

!contains

end module fms_io_mod

