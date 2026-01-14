#include "rave_user_events.h"

void rave_event_and_value_f(int type, int value) {
	 rave_event_and_value(type, value)
}

void rave_name_event_f(int event, char * nam, int length){
	rave_name_event_len(event,nam,length);
}

void rave_name_value_f(int event, int value, char * nam, int length){
	rave_name_value_len(event,value,nam,length);
}

void rave_begin_region_f(char * name, int length){
	rave_begin_region_len(name,length);
}

void rave_end_region_f(char * name, int length){
	rave_end_region_len(name,length);
}

void rave_restart_trace_f(){
 	rave_restart_trace();
}
void rave_enable_trace_f(){
	rave_enable_trace();
}
void rave_disable_trace_f(){
	rave_disable_trace();
}
void rave_enable_regions_f(){
	rave_enable_regions();
}
void rave_disable_regions_f(){
	rave_disable_regions();
}
