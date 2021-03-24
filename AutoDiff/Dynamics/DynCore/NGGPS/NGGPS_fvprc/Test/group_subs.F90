module group_subs_mod

implicit none
private

public group_update

interface group_update
 module procedure group_update_r2_v1
 module procedure group_update_r2_v2
 module procedure group_update_r3_v1
 module procedure group_update_r3_v2
end interface group_update

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine group_update_r2_v1(field1,is,ie,js,je)

 integer, intent(in) :: is,ie,js,je
 real, intent(inout) :: field1(is:ie,js:je)

 field1 = 2*field1

end subroutine group_update_r2_v1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine group_update_r2_v2(field1,field2,is,ie,js,je)

 integer, intent(in) :: is,ie,js,je
 real, intent(inout) :: field1(is:ie,js:je),field2(is:ie,js:je)

 field1 = 2*field1
 field2 = 2*field2

end subroutine group_update_r2_v2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine group_update_r3_v1(field1,is,ie,js,je,k)

 integer, intent(in) :: is,ie,js,je,k
 real, intent(inout) :: field1(is:ie,js:je,k)

 field1 = 2*field1

end subroutine group_update_r3_v1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine group_update_r3_v2(field1,field2,is,ie,js,je,k)

 integer, intent(in) :: is,ie,js,je,k
 real, intent(inout) :: field1(is:ie,js:je,k),field2(is:ie,js:je,k)

 field1 = 2*field1
 field2 = 2*field2

end subroutine group_update_r3_v2

end module group_subs_mod
