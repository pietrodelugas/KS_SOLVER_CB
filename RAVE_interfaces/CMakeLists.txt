set(src_c_raveinterface 
    rave_fortran_c_func.c) 
add_library(c_raveinterface ${src_c_raveinterface}) 
set (src_f_raveinterface rave_user_events.f90) 
qe_add_library(f_raveinterface ${src_f_raveinterface})
target_include_directories(f_raveinterface PRIVATE ${qe_BINARY_DIR}) 
target_link_libraries(f_raveinterface PRIVATE c_raveinterface)
qe_install_targets(f_raveinterface c_raveinterface) 
