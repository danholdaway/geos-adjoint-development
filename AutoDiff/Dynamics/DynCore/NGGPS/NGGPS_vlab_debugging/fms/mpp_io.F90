module mpp_io_mod 

 use mpp_domains_mod, only: domain1D

 implicit none
 private

 integer, parameter :: LONG_KIND = 4

 public axistype, fieldtype

  type :: atttype
     private
     integer             :: type, len
     character(len=128)  :: name
     character(len=1280) :: catt
     real, pointer       :: fatt(:) =>NULL() ! just use type conversion for integers
  end type atttype

  type :: axistype
     private
     character(len=128) :: name
     character(len=128) :: units
     character(len=256) :: longname
     character(len=8)   :: cartesian
     character(len=256) :: compressed
     character(len=24)  :: calendar
     integer            :: sense, len          !+/-1, depth or height?
     type(domain1D)     :: domain              !if pointer is associated, it is a distributed data axis
     real, pointer      :: data(:) =>NULL()    !axis values (not used if time axis)
     integer, pointer   :: idata(:) =>NULL()   !compressed axis valuesi
     integer            :: id, did, type, natt !id is the "variable ID", did is the "dimension ID": 
                                               !netCDF requires 2 IDs for axes
     integer            :: shift               !normally is 0. when domain is symmetry, its value maybe 1.
     type(atttype), pointer :: Att(:) =>NULL()
  end type axistype

  type :: fieldtype
     private
     character(len=128)      :: name
     character(len=128)      :: units
     character(len=256)      :: longname
     character(len=128)      :: standard_name   ! CF standard name
     real                    :: min, max, missing, fill, scale, add
     integer                 :: pack
     integer(LONG_KIND), dimension(3) :: checksum
     type(axistype), pointer :: axes(:) =>NULL() !axes associated with field size, time_axis_index redundantly 
                                        !hold info already contained in axes. it's clunky and inelegant, 
                                        !but required so that axes can be shared among multiple files
     integer, pointer        :: size(:) =>NULL()
     integer                 :: time_axis_index
     integer                 :: id, type, natt, ndim
     type(atttype), pointer  :: Att(:) =>NULL()
     integer                 :: position ! indicate the location of the data ( CENTER, NORTH, EAST, CORNER )
  end type fieldtype


end module mpp_io_mod

