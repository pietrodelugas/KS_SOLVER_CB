subroutine dxlib_errore(location, message, code) 
  implicit none
  character(*) :: location, message
  integer      :: code 
  call errore (location, message, code) 
end subroutine dxlib_errore
