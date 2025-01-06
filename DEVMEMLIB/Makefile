# Makefile for DAVID

include ../make.inc

# location of needed modules and included files (if any)
MODFLAGS=$(MOD_FLAG).

DEVMEMS = device_memcpy.o \
device_memcpy_mod.o 

all : libdev.a


libdev.a:  $(DEVMEMS)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend
