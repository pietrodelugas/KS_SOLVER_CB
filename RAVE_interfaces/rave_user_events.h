#if 0
#define rave_write_hex(x)\
			switch (x){\
							case 0: asm volatile("lui x0, 0\n"); break;\
							case 1: asm volatile("lui x0, 1\n"); break;\
							case 2: asm volatile("lui x0, 2\n"); break;\
							case 3: asm volatile("lui x0, 3\n"); break;\
							case 4: asm volatile("lui x0, 4\n"); break;\
							case 5: asm volatile("lui x0, 5\n"); break;\
							case 6: asm volatile("lui x0, 6\n"); break;\
							case 7: asm volatile("lui x0, 7\n"); break;\
							case 8: asm volatile("lui x0, 8\n"); break;\
							case 9: asm volatile("lui x0, 9\n"); break;\
							case 10: asm volatile("lui x0, 10\n"); break;\
							case 11: asm volatile("lui x0, 11\n"); break;\
							case 12: asm volatile("lui x0, 12\n"); break;\
							case 13: asm volatile("lui x0, 13\n"); break;\
							case 14: asm volatile("lui x0, 14\n"); break;\
							case 15: asm volatile("lui x0, 15\n"); break;\
			}

#include <stdio.h>

static void rave_name(char * name){
		asm volatile("li x0, -1\n");
		for(int i=0; name[i]!='\0'; ++i){
			int y=(int)name[i];
			for(int tmp=y; tmp>0; tmp>>=4){
				rave_write_hex(tmp&0x0f);
			}
		}
		asm volatile("li x0, -1\n");
}

#define rave_name_event(x,name)\
{\
	/*asm volatile("lui x0, %0\n"::"i"(x));*/\
	asm volatile("and x0, %0, %1\n"::"r"(x), "r"(-1));\
	rave_name(name);\
}
#define rave_name_value(x,y,name)\
{\
	/*asm volatile("lui x0, %0\n"::"i"(x));\
	asm volatile("lui x0, %0\n"::"i"(y));*/\
	asm volatile("and x0, %0, %1\n"::"r"(x), "r"(y));\
	rave_name(name);\
}
#else
#define rave_name_event_len(x,name,len)\
{\
	asm volatile("and x0, %0, x0\n" :: "r"(x));\
	asm volatile("sll x0, %0, %1\n" :: "r"(&name[0]), "r"(len));\
}
#define rave_name_value_len(x,y,name,len)\
{\
	asm volatile("and x0, %0, %1\n" :: "r"(x), "r"(y));\
	asm volatile("srl x0, %0, %1\n" :: "r"(&name[0]), "r"(len));\
}

#define rave_name_event(x,name) rave_name_event_len(x,name,-1)
#define rave_name_value(x,y,name) rave_name_value_len(x,y,name,-1)

#endif

#define rave_restart_trace() asm volatile("li x0, -2\n");

//Maintaining old start/stop functions instead of enable/disable
#define rave_start_trace() asm volatile("li x0, -3\n");
#define rave_stop_trace() asm volatile("li x0, -4\n");
//New naming convention:
#define rave_enable_trace() rave_start_trace()
#define rave_disable_trace() rave_stop_trace() 
#define rave_enable_regions() asm volatile("li x0, -7\n");
#define rave_disable_regions() asm volatile("li x0, -8\n");

#define rave_enable() rave_enable_regions() rave_enable_trace()
#define rave_disable() rave_disable_regions() rave_disable_trace()

#define rave_event_and_value(x,y) asm volatile("or x0, %0, %1\n"::"r"(x),"r"(y));

#define rave_begin_region_len(name,len) asm volatile("add x0, %0, %1\n" :: "r"(&name[0]), "r"(len))
#define rave_end_region_len(name,len) asm volatile("sub x0, %0, %1\n" :: "r"(&name[0]), "r"(len))

#define rave_begin_region(name) rave_begin_region_len(name,-1)
#define rave_end_region(name) rave_end_region_len(name,-1) 
