module rave_user_events
  use iso_c_binding
  implicit none
  interface
    subroutine rave_event_and_value(event, val) bind(c, name="rave_event_and_value_f")
       use iso_c_binding
       integer(kind=c_int), value :: event, val
    end subroutine rave_event_and_value

    subroutine rave_restart_trace() bind(c, name="rave_restart_trace_f")
       use iso_c_binding
    end subroutine rave_restart_trace

    subroutine rave_enable_trace() bind(c, name="rave_enable_trace_f")
       use iso_c_binding
    end subroutine rave_enable_trace
    subroutine rave_disable_trace() bind(c, name="rave_disable_trace_f")
       use iso_c_binding
    end subroutine rave_disable_trace

    subroutine rave_enable_regions() bind(c, name="rave_enable_regions_f")
       use iso_c_binding
    end subroutine rave_enable_regions
    subroutine rave_disable_regions() bind(c, name="rave_disable_regions_f")
       use iso_c_binding
    end subroutine rave_disable_regions

    subroutine rave_begin_region_internal(name,length) &
                 bind(c, name="rave_begin_region_f")
       use iso_c_binding
       type(c_ptr), value :: name
       integer(c_int), value :: length
    end subroutine rave_begin_region_internal 

    subroutine rave_end_region_internal(name,length) &
                 bind(c, name="rave_end_region_f")
       use iso_c_binding
       type(c_ptr), value :: name
       integer(c_int), value :: length
    end subroutine rave_end_region_internal 
    
    subroutine rave_name_event_internal(event, nam, length) bind(c, name="rave_name_event_f")
       use iso_c_binding
       integer(kind=c_int), value :: event
       type(c_ptr), value :: nam
       integer(c_int), value :: length
    end subroutine rave_name_event_internal

    subroutine rave_name_value_internal(event, val, nam, length) bind(c, name="rave_name_value_f")
       use iso_c_binding
       integer(kind=c_int), value :: event,val
       type(c_ptr), value :: nam
       integer(c_int), value :: length
    end subroutine rave_name_value_internal

  end interface


contains

  !Maintaining old start/stop functions instead of enable/disable
  subroutine rave_start_trace()
    call rave_enable_trace()
  end subroutine rave_start_trace
  subroutine rave_stop_trace()
    call rave_disable_trace()
  end subroutine rave_stop_trace

  subroutine rave_enable()
    call rave_enable_trace()
    call rave_enable_regions()
  end subroutine rave_enable

  subroutine rave_disable()
    call rave_disable_trace()
    call rave_disable_regions()
  end subroutine rave_disable

  subroutine rave_begin_region(name)
    character(kind=c_char, len=*), target :: name
    type(c_ptr) :: e_cstr_ptr
    e_cstr_ptr = c_loc(name)
    call rave_begin_region_internal(e_cstr_ptr, len(name))
  end subroutine rave_begin_region

  subroutine rave_end_region(name)
    character(kind=c_char, len=*), target :: name
    type(c_ptr) :: e_cstr_ptr
    e_cstr_ptr = c_loc(name)
    call rave_end_region_internal(e_cstr_ptr, len(name))
  end subroutine rave_end_region

  subroutine rave_name_event(event, nam)
    integer(kind=c_int), value :: event
    character(kind=c_char, len=*), target :: nam
    type(c_ptr) :: cstr_ptr
    cstr_ptr = c_loc(nam)
    call rave_name_event_internal(event,cstr_ptr, len(nam))
  end subroutine rave_name_event

  subroutine rave_name_value(event,val, nam)
    integer(kind=c_int), value :: event, val
    character(kind=c_char, len=*), target :: nam
    type(c_ptr) :: cstr_ptr
    cstr_ptr = c_loc(nam)
    call rave_name_value_internal(event,val,cstr_ptr, len(nam))
  end subroutine rave_name_value

end module rave_user_events
