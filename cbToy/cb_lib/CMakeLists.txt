
set(src_cb
  cb_module.f90 
  command_line_options.f90 
  ggen.f90  
  init_k.f90 
  init_random_wfcs.f90 
  input.f90
  int_to_char.f90  
  io_global.f90 
  mp_bands.f90  
  mp_global.f90  
  mp_world.f90 
  random_numbers.f90 
  set_cb_potential.f90 
  sort.f90
  write_bands.f90 
)
qe_add_library(cb_lib ${src_cb}) 
qe_install_targets(cb_lib)
target_link_libraries(cb_lib
	PUBLIC
	qe_fftx
        qe_utilx
        qe_lax)
