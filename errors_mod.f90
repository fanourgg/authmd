module errors_mod



Contains

   subroutine handle_error(str, file, line)
   implicit none
   character(len=*) :: str, file
   integer :: line

   print*,'!! ERROR in file: ', trim(file), ' line: ', line, ' :: ', trim(str)
   stop

   end subroutine handle_error
   !
   subroutine handle_warning(str, file, line)
   implicit none
   character(len=*) :: str, file
   integer :: line

   print*,'!! WARNING in file: ', trim(file), ' line: ', line, ' :: ', trim(str)

   end subroutine handle_warning
end module errors_mod
