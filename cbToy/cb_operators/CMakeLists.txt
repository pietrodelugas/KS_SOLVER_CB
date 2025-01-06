
set(src_cbops
   cb_g_psi.f90  
   cb_h_psi.f90  
   cb_s_psi.f90)
qe_add_library(cb_ops ${src_cbops}) 
qe_install_targets(cb_ops)
target_link_libraries(cb_ops
	PUBLIC
	qe_fftx
        qe_utilx
        qe_lax
	PRIVATE
        cb_lib)
