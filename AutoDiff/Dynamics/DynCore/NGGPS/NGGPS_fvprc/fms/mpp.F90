module mpp_mod

use fv_arrays_mod, only: FVPRC

implicit none
private

public mpp_error, FATAL, mpp_sum, mpp_sync, mpp_npes, mpp_broadcast, WARNING, mpp_pe, &
       mpp_send, mpp_recv, mpp_sync_self, mpp_max

integer :: FATAL !dummy
integer :: mpp_pe, mpp_npes
integer :: WARNING 
integer, parameter :: NAME_LENGTH = 64
integer :: mpp_max !dummy

contains

 subroutine mpp_error(fatal,statement)

  integer, intent(in) :: fatal
  character(len=*), intent(in) :: statement

   if (fatal > 0) then
      print*, statement
   endif

 endsubroutine mpp_error

 subroutine mpp_sum(tilelist,npes)

  integer, intent(in) :: tilelist(:), npes
  

 endsubroutine mpp_sum

 subroutine mpp_broadcast(field,fsize,n,pelist)
 
  real(FVPRC), intent(inout) :: field(:,:,:)
  integer, intent(in) :: fsize,n,pelist(:)

   field = 2*field

 endsubroutine mpp_broadcast

 subroutine mpp_send(field,fieldsize,recv_proc)

  real(FVPRC), intent(inout) :: field(:,:,:)
  integer, intent(in) :: fieldsize, recv_proc

   field = 2*field

 endsubroutine mpp_send

 subroutine mpp_recv(field,fieldsize,send_proc)

  real(FVPRC), intent(inout) :: field(:,:,:)
  integer, intent(in) :: fieldsize, send_proc

   field = 2*field

 endsubroutine mpp_recv

 subroutine mpp_sync

  print*, 'mpp_sync'

 endsubroutine mpp_sync

 subroutine mpp_sync_self

  print*, 'mpp_sync_self'

 endsubroutine mpp_sync_self

end module mpp_mod
