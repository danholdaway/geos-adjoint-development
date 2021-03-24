module group_call_mod

use group_subs_mod

implicit none
private

public group_call_test

contains

subroutine group_call_test(u2,v2,u3,v3,isd,ied,jsd,jed,npz)

 integer, intent(in) :: isd,ied,jsd,jed,npz
 real, intent(inout) :: u2(isd:ied,jsd:jed)    ,v2(isd:ied,jsd:jed)
 real, intent(inout) :: u3(isd:ied,jsd:jed,npz),v3(isd:ied,jsd:jed,npz)

 call group_update(u2,isd,ied,jsd,jed)

 call group_update(u2,v2,isd,ied,jsd,jed)

 call group_update(u3,isd,ied,jsd,jed,npz)

 call group_update(u3,v3,isd,ied,jsd,jed,npz)

end subroutine group_call_test


end module group_call_mod
