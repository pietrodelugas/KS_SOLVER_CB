SUBROUTINE f_print_rmss()
   USE, INTRINSIC :: iso_c_binding
   IMPLICIT NONE 
   INTERFACE 
     FUNCTION subcheckmem() RESULT (check) BIND(C, name='subcheckmem') 
                  INTEGER :: check
     END FUNCTION subcheckmem
   END INTERFACE
   integer :: retval
   retval = subcheckmem()
END SUBROUTINE f_print_rmss

