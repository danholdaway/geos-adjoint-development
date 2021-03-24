module tracer_manager_mod

implicit none
private

public get_tracer_index

contains

 integer function get_tracer_index (MODEL_ATMOS, TrName) 

  integer, intent(in) :: MODEL_ATMOS
  character(len=*) :: TrName

  get_tracer_index = MODEL_ATMOS

 end function get_tracer_index

end module tracer_manager_mod

