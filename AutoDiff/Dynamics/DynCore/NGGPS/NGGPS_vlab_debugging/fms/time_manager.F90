module time_manager_mod

implicit none
private

public time_type

type time_type
   private
   integer:: seconds
   integer:: days
   integer:: ticks
   integer:: dummy ! added as a workaround bug on IRIX64 (AP)
end type time_type

end module time_manager_mod

