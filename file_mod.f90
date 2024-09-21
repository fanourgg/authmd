module file_mod
implicit none
private
!public :: iun_master, iun_out, iun_log, iun_gloc
public :: iun_out, iun_log

integer :: iun_log
!integer :: iun_master  !< the general output log-file
integer :: iun_out = 111
!integer :: iun_gloc = 1212  !< general log-file

end module file_mod
