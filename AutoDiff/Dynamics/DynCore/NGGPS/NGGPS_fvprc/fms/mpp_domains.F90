module mpp_domains_mod

implicit none
private

public DGRID_NE, CGRID_NE, &
       mpp_get_boundary, mpp_update_domains, &
       nest_domain_type, &
       domain1d, domain2d, &
       mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain, &
       CENTER, CORNER, NORTH, EAST, SOUTH, WEST, &
       mpp_get_pelist, &
       mpp_global_field, &
       mpp_get_C2F_index, mpp_get_F2C_index, &
       mpp_update_nest_fine, &
       mpp_update_nest_coarse, &
       mpp_group_update_type, &
       mpp_global_sum, BITWISE_EXACT_SUM

  !!!!!!!!!FV_ARRAYS!!!!!!!!!!!!
  integer, parameter :: FVPRC = 8
  integer, parameter :: r_grid = 8
  !!!!!!!!!FV_ARRAYS!!!!!!!!!!!!


  integer :: DGRID_NE, CGRID_NE
  integer :: mpp_get_global_domain
  integer :: CENTER, CORNER, NORTH, EAST, SOUTH, WEST
  integer :: mpp_get_pelist !dummy
  integer, parameter :: BITWISE_EXACT_SUM = 1

  integer, parameter :: MAXOVERLAP = 100
  integer, parameter :: NAME_LENGTH = 64
  integer, parameter :: MAX_DOMAIN_FIELDS=100
  integer, parameter :: MAX_REQUEST = 100
  integer, parameter :: LONG_KIND = 4

  type domain_axis_spec
   private
   integer :: dummy1
  end type domain_axis_spec

  type domain1D
   private
   integer :: dummy1
  end type domain1D

  type domain1D_spec
   private
   integer :: dummy1
  end type domain1D_spec

  type overlap_type
   private
   integer :: dummy1
  end type overlap_type

  type overlapSpec
   private
   integer :: dummy1
  end type overlapSpec

  type domain2D_spec
   private
   integer :: dummy1
  end type domain2D_spec

  type tile_type
   private
   integer :: dummy1
  end type tile_type

  type index_type
   private
   integer :: dummy1
  end type index_type

  type nestSpec
   private
   integer :: dummy1
  end type nestSpec

!  type nest_domain_type
!   private
!   integer :: dummy1
!  end type nest_domain_type

  type mpp_group_update_type
   private
   integer :: dummy1
  end type mpp_group_update_type

!  type domain_axis_spec        !type used to specify index limits along an axis of a domain
!     private
!     integer :: begin, end, size, max_size      !start, end of domain axis, size, max size in set
!     logical :: is_global       !TRUE if domain axis extent covers global domain
!  end type domain_axis_spec
!
!  type domain1D
!     private
!     type(domain_axis_spec) :: compute, data, global, memory
!     logical :: cyclic
!     type(domain1D), pointer :: list(:) =>NULL()
!     integer :: pe               !PE to which this domain is assigned
!     integer :: pos              !position of this PE within link list, i.e domain%list(pos)%pe = pe
!     integer :: goffset, loffset !needed for global sum
!  end type domain1D
!
!  type domain1D_spec
!     private
!     type(domain_axis_spec) :: compute
!     integer                :: pos
!  end type domain1D_spec
!
!  type overlap_type
!     private
!     integer                  :: count = 0                 ! number of ovrelapping
!     integer                  :: pe
!     integer                  :: start_pos                 ! start position in the buffer
!     integer                  :: totsize                   ! all message size
!     integer ,        pointer :: msgsize(:)      => NULL() ! overlapping msgsize to be sent or received
!     integer,         pointer :: tileMe(:)       => NULL() ! my tile id for this overlap
!     integer,         pointer :: tileNbr(:)      => NULL() ! neighbor tile id for this overlap
!     integer,         pointer :: is(:)           => NULL() ! starting i-index 
!     integer,         pointer :: ie(:)           => NULL() ! ending   i-index 
!     integer,         pointer :: js(:)           => NULL() ! starting j-index 
!     integer,         pointer :: je(:)           => NULL() ! ending   j-index 
!     integer,         pointer :: dir(:)          => NULL() ! direction ( value 1,2,3,4 = E,S,W,N)
!     integer,         pointer :: rotation(:)     => NULL() ! rotation angle.
!     integer,         pointer :: index(:)        => NULL() ! for refinement
!     logical,         pointer :: from_contact(:) => NULL() ! indicate if the overlap is computed from define_contact_overlap
!  end type overlap_type
!
!  type overlapSpec
!     private
!     integer                     :: whalo, ehalo, shalo, nhalo ! halo size
!     integer                     :: xbegin, xend, ybegin, yend
!     integer                     :: nsend, nrecv
!     integer                     :: sendsize, recvsize
!     type(overlap_type), pointer :: send(:) => NULL()
!     type(overlap_type), pointer :: recv(:) => NULL()
!     type(overlapSpec),  pointer :: next
!  end type overlapSpec
!
!  type domain2D_spec
!     private
!     type(domain1D_spec), pointer :: x(:)       => NULL() ! x-direction domain decomposition
!     type(domain1D_spec), pointer :: y(:)       => NULL() ! x-direction domain decomposition
!     integer,        pointer :: tile_id(:) => NULL() ! tile id of each tile
!     integer                 :: pe                   ! PE to which this domain is assigned
!     integer                 :: pos                  ! position of this PE within link list
!     integer                 :: tile_root_pe         ! root pe of tile.
!  end type domain2D_spec
!
!  type tile_type
!     integer :: xbegin, xend, ybegin, yend
!  end type tile_type

  type domain2D
     private
     character(len=NAME_LENGTH)  :: name='unnamed'          ! name of the domain, default is "unspecified"
     integer(LONG_KIND)          :: id 
     integer                     :: pe                      ! PE to which this domain is assigned
     integer                     :: fold          
     integer                     :: pos                     ! position of this PE within link list
     logical                     :: symmetry                ! indicate the domain is symmetric or non-symmetric.
     integer                     :: whalo, ehalo            ! halo size in x-direction
     integer                     :: shalo, nhalo            ! halo size in y-direction
     integer                     :: ntiles                  ! number of tiles within mosaic
     integer                     :: max_ntile_pe            ! maximum value in the pelist of number of tiles on each pe.
     integer                     :: ncontacts               ! number of contact region within mosaic.
     logical                     :: rotated_ninety          ! indicate if any contact rotate NINETY or MINUS_NINETY
     logical                     :: initialized=.FALSE.     ! indicate if the overlapping is computed or not.
     integer                     :: tile_root_pe            ! root pe of current tile.
     integer                     :: io_layout(2)            ! io_layout, will be set through mpp_define_io_domain
                                                            ! default = domain layout
     integer,            pointer :: pearray(:,:)  => NULL() ! pe of each layout position 
     integer,            pointer :: tile_id(:)    => NULL() ! tile id of each tile
     type(domain1D),     pointer :: x(:)          => NULL() ! x-direction domain decomposition
     type(domain1D),     pointer :: y(:)          => NULL() ! y-direction domain decomposition
     type(domain2D_spec),pointer :: list(:)       => NULL() ! domain decomposition on pe list
     type(tile_type),    pointer :: tileList(:)   => NULL() ! store tile information
     type(overlapSpec),  pointer :: check_C       => NULL() ! send and recv information for boundary consistency check of C-cell
     type(overlapSpec),  pointer :: check_E       => NULL() ! send and recv information for boundary consistency check of E-cell
     type(overlapSpec),  pointer :: check_N       => NULL() ! send and recv information for boundary consistency check of N-cell
     type(overlapSpec),  pointer :: bound_C       => NULL() ! send information for getting boundary value for symmetry domain.
     type(overlapSpec),  pointer :: bound_E       => NULL() ! send information for getting boundary value for symmetry domain.
     type(overlapSpec),  pointer :: bound_N       => NULL() ! send information for getting boundary value for symmetry domain.
     type(overlapSpec),  pointer :: update_T      => NULL() ! send and recv information for halo update of T-cell.
     type(overlapSpec),  pointer :: update_E      => NULL() ! send and recv information for halo update of E-cell.
     type(overlapSpec),  pointer :: update_C      => NULL() ! send and recv information for halo update of C-cell.
     type(overlapSpec),  pointer :: update_N      => NULL() ! send and recv information for halo update of N-cell.
     type(domain2d),     pointer :: io_domain     => NULL() ! domain for IO, will be set through calling mpp_set_io_domain ( this will be changed).
  end type domain2D   

  type nest_domain_type
     private
     integer                    :: tile_fine, tile_coarse
     integer                    :: istart_fine, iend_fine, jstart_fine, jend_fine
     integer                    :: istart_coarse, iend_coarse, jstart_coarse, jend_coarse
     integer                    :: x_refine, y_refine
     logical                    :: is_fine_pe, is_coarse_pe
     integer,           pointer :: pelist_fine(:) => NULL()
     integer,           pointer :: pelist_coarse(:) => NULL()
     character(len=NAME_LENGTH) :: name
     type(nestSpec), pointer :: C2F_T => NULL()
     type(nestSpec), pointer :: C2F_C => NULL()
     type(nestSpec), pointer :: C2F_E => NULL()
     type(nestSpec), pointer :: C2F_N => NULL()
     type(nestSpec), pointer :: F2C_T => NULL()
     type(nestSpec), pointer :: F2C_C => NULL()
     type(nestSpec), pointer :: F2C_E => NULL()
     type(nestSpec), pointer :: F2C_N => NULL()
     type(domain2d), pointer :: domain_fine   => NULL()
     type(domain2d), pointer :: domain_coarse => NULL()
  end type nest_domain_type

!  type index_type
!     integer :: is_me, ie_me, js_me, je_me
!     integer :: is_you, ie_you, js_you, je_you
!  end type index_type
!
!  type nestSpec
!     private
!     integer                     :: xbegin, xend, ybegin, yend
!     type(index_type)            :: west, east, south, north, center
!     integer                     :: nsend, nrecv
!     integer                     :: extra_halo
!     type(overlap_type), pointer :: send(:) => NULL()
!     type(overlap_type), pointer :: recv(:) => NULL()
!     type(nestSpec),     pointer :: next
!  end type nestSpec
!
!  type nest_domain_type
!     private
!     integer                    :: tile_fine, tile_coarse
!     integer                    :: istart_fine, iend_fine, jstart_fine, jend_fine
!     integer                    :: istart_coarse, iend_coarse, jstart_coarse, jend_coarse
!     integer                    :: x_refine, y_refine
!     logical                    :: is_fine_pe, is_coarse_pe
!     integer,           pointer :: pelist_fine(:) => NULL()
!     integer,           pointer :: pelist_coarse(:) => NULL()
!     character(len=NAME_LENGTH) :: name
!     type(nestSpec), pointer :: C2F_T => NULL()
!     type(nestSpec), pointer :: C2F_C => NULL()
!     type(nestSpec), pointer :: C2F_E => NULL()
!     type(nestSpec), pointer :: C2F_N => NULL()
!     type(nestSpec), pointer :: F2C_T => NULL()
!     type(nestSpec), pointer :: F2C_C => NULL()
!     type(nestSpec), pointer :: F2C_E => NULL()
!     type(nestSpec), pointer :: F2C_N => NULL()
!     type(domain2d), pointer :: domain_fine   => NULL()
!     type(domain2d), pointer :: domain_coarse => NULL()
!  end type nest_domain_type
!
!  type mpp_group_update_type
!     private
!     logical            :: initialized = .FALSE.
!     logical            :: k_loop_inside = .TRUE.
!     integer            :: nscalar = 0
!     integer            :: nvector = 0
!     integer            :: flags_s=0, flags_v=0
!     integer            :: whalo_s=0, ehalo_s=0, shalo_s=0, nhalo_s=0
!     integer            :: isize_s=0, jsize_s=0, ksize_s=1
!     integer            :: whalo_v=0, ehalo_v=0, shalo_v=0, nhalo_v=0
!     integer            :: isize_x=0, jsize_x=0, ksize_v=1
!     integer            :: isize_y=0, jsize_y=0
!     integer            :: position=0, gridtype=0
!     logical            :: recv_s(8), recv_v(8)
!     integer            :: is_s=0, ie_s=0, js_s=0, je_s=0
!     integer            :: is_x=0, ie_x=0, js_x=0, je_x=0
!     integer            :: is_y=0, ie_y=0, js_y=0, je_y=0
!     integer            :: nrecv=0, nsend=0
!     integer            :: npack=0, nunpack=0
!     integer            :: reset_index_s = 0
!     integer            :: reset_index_v = 0
!     integer            :: tot_msgsize = 0
!     integer            :: from_pe(MAXOVERLAP)
!     integer            :: to_pe(MAXOVERLAP)
!     integer            :: recv_size(MAXOVERLAP)
!     integer            :: send_size(MAXOVERLAP)
!     integer            :: buffer_pos_recv(MAXOVERLAP)
!     integer            :: buffer_pos_send(MAXOVERLAP)
!     integer            :: pack_type(MAXOVERLAP)
!     integer            :: pack_buffer_pos(MAXOVERLAP)
!     integer            :: pack_rotation(MAXOVERLAP)
!     integer            :: pack_size(MAXOVERLAP)
!     integer            :: pack_is(MAXOVERLAP)
!     integer            :: pack_ie(MAXOVERLAP)
!     integer            :: pack_js(MAXOVERLAP)
!     integer            :: pack_je(MAXOVERLAP)
!     integer            :: unpack_type(MAXOVERLAP)
!     integer            :: unpack_buffer_pos(MAXOVERLAP)
!     integer            :: unpack_rotation(MAXOVERLAP)
!     integer            :: unpack_size(MAXOVERLAP)
!     integer            :: unpack_is(MAXOVERLAP)
!     integer            :: unpack_ie(MAXOVERLAP)
!     integer            :: unpack_js(MAXOVERLAP)
!     integer            :: unpack_je(MAXOVERLAP)
!     integer(LONG_KIND) :: addrs_s(MAX_DOMAIN_FIELDS)
!     integer(LONG_KIND) :: addrs_x(MAX_DOMAIN_FIELDS)
!     integer(LONG_KIND) :: addrs_y(MAX_DOMAIN_FIELDS)
!     integer            :: buffer_start_pos = -1
!     integer            :: request_send(MAX_REQUEST)
!     integer            :: request_recv(MAX_REQUEST)
!     integer            :: type_recv(MAX_REQUEST)
!  end type mpp_group_update_type


  interface mpp_update_nest_fine
     module procedure mpp_update_nest_fine_2d
     module procedure mpp_update_nest_fine_3d
     module procedure mpp_update_nest_fine_4d
  end interface

  interface mpp_global_field
     module procedure mpp_global_field2D_2d
     module procedure mpp_global_field2D_3d
     module procedure mpp_global_field2D_4d
     module procedure mpp_global_field2D_5d
  end interface

  interface mpp_update_domains
     module procedure mpp_update_domain2D_2d
     module procedure mpp_update_domain2D_3d
     module procedure mpp_update_domain2D_4d
     module procedure mpp_update_domain2D_5d
     module procedure mpp_update_domain2D_2dv
     module procedure mpp_update_domain2D_3dv
     module procedure mpp_update_domain2D_4dv
     module procedure mpp_update_domain2D_5dv
  end interface

  interface mpp_global_sum
     module procedure mpp_global_sum_2d
     module procedure mpp_global_sum_3d
     module procedure mpp_global_sum_4d
     module procedure mpp_global_sum_5d
  end interface

 contains




! mpp_update_domains
! ------------------

 subroutine mpp_update_domain2D_2d(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(FVPRC), intent(inout) :: field1(:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_2d
 subroutine mpp_update_domain2D_3d(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(FVPRC), intent(inout) :: field1(:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_3d
 subroutine mpp_update_domain2D_4d(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(FVPRC), intent(inout) :: field1(:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_4d
 subroutine mpp_update_domain2D_5d(field1,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(FVPRC), intent(inout) :: field1(:,:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1

 endsubroutine mpp_update_domain2D_5d

 subroutine mpp_update_domain2D_2dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(FVPRC), intent(inout) :: field1(:,:), field2(:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2

 endsubroutine mpp_update_domain2D_2dv
 subroutine mpp_update_domain2D_3dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(FVPRC), intent(inout) :: field1(:,:,:), field2(:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2


 endsubroutine mpp_update_domain2D_3dv
 subroutine mpp_update_domain2D_4dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(FVPRC), intent(inout) :: field1(:,:,:,:), field2(:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2

 endsubroutine mpp_update_domain2D_4dv
 subroutine mpp_update_domain2D_5dv(field1,field2,domain,gridtype,complete,whalo,ehalo,shalo,nhalo)

  real(FVPRC), intent(inout) :: field1(:,:,:,:,:), field2(:,:,:,:,:)
  type(domain2d), intent(in) :: domain
  integer, optional, intent(in) :: gridtype
  logical, optional, intent(in) :: complete
  integer, optional, intent(in) :: whalo,ehalo,shalo,nhalo

  field1 = 2*field1
  field2 = 2*field2

 endsubroutine mpp_update_domain2D_5dv



! mpp_get_boundary
! ----------------

 subroutine mpp_get_boundary(field1,field2,domain,wbuffery,ebuffery,sbuffery,nbuffery,wbufferx,ebufferx,sbufferx,nbufferx,gridtype)

  real(FVPRC), intent(inout) :: field1(:,:,:),field2(:,:,:)
  type(domain2d), intent(in) :: domain
  real(FVPRC), optional, intent(inout) :: wbuffery(:,:),ebuffery(:,:),sbuffery(:,:),nbuffery(:,:)
  real(FVPRC), optional, intent(inout) :: wbufferx(:,:),ebufferx(:,:),sbufferx(:,:),nbufferx(:,:)
  integer, optional, intent(in) :: gridtype

  field1 = 2*field1
  field2 = 2*field2

  wbuffery(:,:) = field1(1,:,:)
  ebuffery(:,:) = field1(2,:,:)
  sbufferx(:,:) = field2(1,:,:)
  nbufferx(:,:) = field2(2,:,:)
  
 endsubroutine mpp_get_boundary



! mpp_get_data_domain
! -------------------

 subroutine mpp_get_data_domain(domain,isd,ied,jsd,jed)

  type(domain2d), intent(in) :: domain
  integer, intent(in) :: isd,ied,jsd,jed

 endsubroutine mpp_get_data_domain



! mpp_get_compute_domain
! ----------------------
 subroutine mpp_get_compute_domain(domain,isd,ied,jsd,jed)

  type(domain2d), intent(in) :: domain
  integer, intent(in) :: isd,ied,jsd,jed

 endsubroutine mpp_get_compute_domain



! mpp_global_field
! ----------------
 subroutine mpp_global_field2d_2d(domain,field_this_grid,field,position)

  type(domain2D), intent(in) :: domain
  real(FVPRC), intent(in) :: field_this_grid(:,:)
  real(FVPRC), intent(inout) :: field(:,:)
  integer, intent(in), optional :: position

   field = 2*field_this_grid

 endsubroutine mpp_global_field2d_2d
 subroutine mpp_global_field2d_3d(domain,field_this_grid,field,position)

  type(domain2D), intent(in) :: domain
  real(FVPRC), intent(in) :: field_this_grid(:,:,:)
  real(FVPRC), intent(inout) :: field(:,:,:)
  integer, intent(in), optional :: position

   field = 2*field_this_grid

 endsubroutine mpp_global_field2d_3d
 subroutine mpp_global_field2d_4d(domain,field_this_grid,field,position)

  type(domain2D), intent(in) :: domain
  real(FVPRC), intent(in) :: field_this_grid(:,:,:,:)
  real(FVPRC), intent(inout) :: field(:,:,:,:)
  integer, intent(in), optional :: position

   field = 2*field_this_grid

 endsubroutine mpp_global_field2d_4d
 subroutine mpp_global_field2d_5d(domain,field_this_grid,field,position)

  type(domain2D), intent(in) :: domain
  real(FVPRC), intent(in) :: field_this_grid(:,:,:,:,:)
  real(FVPRC), intent(inout) :: field(:,:,:,:,:)
  integer, intent(in), optional :: position

   field = 2*field_this_grid

 endsubroutine mpp_global_field2d_5d



! mpp_get_C2F_index
! -----------------

 subroutine mpp_get_C2F_index(nestdomain, is_f, ie_f, js_f, je_f, is_c, ie_c, js_c, je_c, SIDE, position)

  type(nest_domain_type), intent(inout) :: nestdomain
  integer, intent(out) :: is_f, ie_f, js_f, je_f, is_c, ie_c, js_c, je_c
  integer, intent(in) :: SIDE, position

   is_f = 1
   ie_f = 20
   js_f = 1
   je_f = 20
   is_c = 1
   ie_c = 20
   js_c = 1
   je_c = 20

 endsubroutine mpp_get_C2F_index



! mpp_get_F2C_index
! -----------------

 subroutine mpp_get_F2C_index(nestdomain, is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f, position)

  type(nest_domain_type), intent(inout) :: nestdomain
  integer, intent(out) :: is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f
  integer, intent(in) :: position

   is_f = 1
   ie_f = 20
   js_f = 1
   je_f = 20
   is_c = 1
   ie_c = 20
   js_c = 1
   je_c = 20

 endsubroutine mpp_get_F2C_index



! mpp_update_nest_fine
! --------------------

 subroutine mpp_update_nest_fine_2d(field, nestdomain, wbuffer, sbuffer, ebuffer, nbuffer,  position)

  real(FVPRC), intent(in) :: field(:,:)
  type(nest_domain_type), intent(inout) :: nestdomain
  real(FVPRC), intent(inout) :: wbuffer(:,:), sbuffer(:,:), ebuffer(:,:), nbuffer(:,:)
  integer, intent(in), optional :: position

  wbuffer = 2*field
  ebuffer = 2*field
  sbuffer = 2*field
  nbuffer = 2*field

 endsubroutine mpp_update_nest_fine_2d
 subroutine mpp_update_nest_fine_3d(field, nestdomain, wbuffer, sbuffer, ebuffer, nbuffer,  position)

  real(FVPRC), intent(in) :: field(:,:,:)
  type(nest_domain_type), intent(inout) :: nestdomain
  real(FVPRC), intent(inout) :: wbuffer(:,:,:), sbuffer(:,:,:), ebuffer(:,:,:), nbuffer(:,:,:)
  integer, intent(in), optional :: position

  wbuffer = 2*field
  ebuffer = 2*field
  sbuffer = 2*field
  nbuffer = 2*field

 endsubroutine mpp_update_nest_fine_3d
 subroutine mpp_update_nest_fine_4d(field, nestdomain, wbuffer, sbuffer, ebuffer, nbuffer,  position)

  real(FVPRC), intent(in) :: field(:,:,:,:)
  type(nest_domain_type), intent(inout) :: nestdomain
  real(FVPRC), intent(inout) :: wbuffer(:,:,:,:), sbuffer(:,:,:,:), ebuffer(:,:,:,:), nbuffer(:,:,:,:)
  integer, intent(in), optional :: position

  wbuffer = 2*field
  ebuffer = 2*field
  sbuffer = 2*field
  nbuffer = 2*field

 endsubroutine mpp_update_nest_fine_4d


! mpp_update_nest_coarse
! ----------------------

 subroutine mpp_update_nest_coarse(field, nestdomain, nestdat, position)

  real(FVPRC), intent(inout) :: field(:,:,:)
  type(nest_domain_type), intent(inout) :: nestdomain
  real(FVPRC), intent(in) :: nestdat(:,:,:)
  integer, intent(in) :: position

   field = 2*field

 endsubroutine mpp_update_nest_coarse


! mpp_global_sum
! --------------

 real(kind=R_GRID) function mpp_global_sum_2d( domain, field, flags, position, tile_count) 

    type(domain2D), intent(in) :: domain
    real(R_GRID), intent(in) :: field(:,:)
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: position
    integer, intent(in), optional :: tile_count

    mpp_global_sum_2d = sum(field)

 endfunction mpp_global_sum_2d

 real(kind=R_GRID) function mpp_global_sum_3d( domain, field, flags, position, tile_count) 

    type(domain2D), intent(in) :: domain
    real(R_GRID), intent(in) :: field(:,:,:)
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: position
    integer, intent(in), optional :: tile_count

    mpp_global_sum_3d = sum(field)

 endfunction mpp_global_sum_3d

 real(kind=R_GRID) function mpp_global_sum_4d( domain, field, flags, position, tile_count) 

    type(domain2D), intent(in) :: domain
    real(R_GRID), intent(in) :: field(:,:,:,:)
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: position
    integer, intent(in), optional :: tile_count

    mpp_global_sum_4d = sum(field)

 endfunction mpp_global_sum_4d

 real(kind=R_GRID) function mpp_global_sum_5d( domain, field, flags, position, tile_count) 

    type(domain2D), intent(in) :: domain
    real(R_GRID), intent(in) :: field(:,:,:,:,:)
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: position
    integer, intent(in), optional :: tile_count

    mpp_global_sum_5d = sum(field)

 endfunction mpp_global_sum_5d


end module mpp_domains_mod

