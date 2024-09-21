#define print_error(ierr, iun_out, str1) &
if (ierr>0) then;write(iun_out, '(70("=")/"FATAL ERROR in file: ",a,2x,"line: ",i5, 2x/a,/,70("="))')__FILE__,__LINE__,str1; stop; endif

#define print_warning(ierr, iun_out, str1) &
if (ierr>0) then;write(iun_out, '(70("=")/"WARNING in file: ",a,2x,"line: ",i5, 2x/a,/,70("="))')__FILE__,__LINE__,str1; endif
