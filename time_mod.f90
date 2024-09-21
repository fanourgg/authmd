!------------------------------------------------------------------------------
!> MODULE elec_mod
!------------------------------------------------------------------------------
!>
!> @author G.S.Fanourgakis
!>
!> @brief It handles timers for measuring time spend at each part of the code
!>
!> @todo 
!> 
!>-----------------------------------------------------------------------------
module time_mod

   type t_timer
       character(len=20) :: name
       integer :: Ncalls
       double precision :: start, time
   end type t_timer

   integer, parameter :: Ntimers=100

   integer, parameter :: itm_tot=1;
   integer, parameter :: itm_tg = 2;
   integer, parameter :: itm_pot = 20;
   integer, parameter :: itm_bonds=21;
   integer, parameter :: itm_angles=22;
   integer, parameter :: itm_dihedrals=23;
   integer, parameter :: itm_impropers=24;
   integer, parameter :: itm_urey=25;
   integer, parameter :: itm_waters=26;
   integer, parameter :: itm_pairs=30;
  integer, parameter :: itm_potrec=39;
   integer, parameter :: itm_ewaldstdqq=40;
   integer, parameter :: itm_ewaldstdqd=41;
   integer, parameter :: itm_ewaldstddd=42;
   integer, parameter :: itm_ewaldsincos=43;
   integer, parameter :: itm_vlist=44;
   type (t_timer), dimension(Ntimers) :: timer

   
contains
   subroutine timer_setup
   implicit none
   integer :: i
    
   do i=1, Ntimers
      timer(i)% name=''
   enddo
   
   timer(itm_tot)% name = 'TOT'
   timer(itm_tg)% name = 'Test Gradients'
   timer(itm_pot)% name = 'POT'
   timer(itm_potrec)% name = 'POTrec'
   timer(itm_pairs)% name = 'Pairs'
   timer(itm_waters)% name = 'Waters'
   timer(itm_bonds)% name = 'Bonds'
   timer(itm_angles)% name = 'Angles'
   timer(itm_dihedrals)% name = 'Torsions'
   timer(itm_impropers)% name = 'ImpTors'
   timer(itm_urey)% name = 'Urey'
   timer(itm_urey)% name = 'Water'
   timer(itm_ewaldstdqq)% name = 'Ewald_STD_QQ'
   timer(itm_ewaldstddd)% name = 'Ewald_STD_DD'
   timer(itm_ewaldstdqd)% name = 'Ewald_STD_QD'
   timer(itm_ewaldsincos)% name = 'Ewald_SINCOS'
   timer(itm_vlist)% name = 'Verlet List'

   do i=1, Ntimers
      timer(i)% Ncalls = 0
   enddo
   end subroutine timer_setup

   subroutine timer_start(itp)
   implicit none
   integer, intent( in    ) :: itp

   timer(itp)% start = get_time()
   end subroutine timer_start

   subroutine timer_stop(itp)
   implicit none
   integer, intent( in    ) :: itp

   timer(itp)% time = timer(itp)% time + get_time() - timer(itp)% start
   timer(itp)% Ncalls = timer(itp)% Ncalls + 1
   end subroutine timer_stop

   subroutine timer_print(itp)
   implicit none
   integer :: itp

   if (timer(itp)% Ncalls>0) &
          print*,'time: ', trim(timer(itp)%name), timer(itp)%time, timer(itp)% Ncalls
   end subroutine timer_print

   subroutine timers_print
   implicit none
   integer :: itp

   write(*,'(//80("*")/10("*"),2x,"ANAYTICAL EXECUTION TIMES",2x,10("*")/80("*"))')
   write(*,99)"SUBR/TASK", "time(tot)", "#calls", "time/call", "%PERC"
   do itp=1, Ntimers
      if (timer(itp)%ncalls > 0) then
         write(*,100)trim(timer(itp)%name), timer(itp)%time, timer(itp)% Ncalls, &
                   timer(itp)%time/timer(itp)% Ncalls, timer(itp)%time/timer(itm_tot)%time*100.
      endif
   enddo
99  format(/a20,4x, a12,   4x, a8, 4x, a12,   4x, a6/)
100 format(a20,4x, f12.3, 4x, i8, 4x, f12.6, 4x, f6.2)
   end subroutine timers_print
   !--------------------------------------------------------------------------
   !>  @author  GSF
   !!  @brief get system time making a system call
   !!  @todo set a timer if no OPENMP
   !<
   !--------------------------------------------------------------------------
   double precision function get_time()
!      integer(4) ::  mclock
!      integer(4) ::  clock_count, clock_rate, clock_max
!      double precision :: aaa
#ifdef _OPENMP
   use omp_lib
!      double precision, external :: omp_get_wtime
#endif
       integer :: count_rate, count, count_max

#ifdef _OPENMP
      get_time = omp_get_wtime()
!       call cpu_time(get_time)
!       call SYSTEM_CLOCK(count, count_rate, count_max)
!       print*,'count=', count, count_rate, count_max
!       get_time=dble(count)
#else
!       get_time = 0.d0
       call cpu_time(get_time)
#endif
   end function get_time
   !--------------------------------------------------------------------------
   !>  @author  GSF
   !!  @brief It prints-out  the time in the (DD/HH/MM/SS) format
   !<
   !--------------------------------------------------------------------------
   subroutine print_curtime(time)
   use file_mod,           ONLY: iun_out
   implicit none
   double precision :: time
   double precision :: ttt
   integer :: t(4)
   
   ttt=time
   
   t(1) = int(ttt/86400.d0) ;               ! days
   ttt = ttt - dble(t(1))*86400.d0
   t(2) = int(ttt/3600.d0)   ! hours
   ttt = ttt - dble(t(2)*3600)
   t(3) = int(ttt/60.d0)      ! minutes
   ttt = ttt - dble(t(3)*60)
   t(4) = int(ttt)            ! seconds

   write(iun_out,'(5x,a,1x,i4,"d",i2,"h",i2,"m",i2,"s")')"*** Current time : ",t(1:4)
   
   end subroutine print_curtime
   !--------------------------------------------------------------------------
   !>  @author  GSF
   !!  @brief It puts to the string get_timestr  the time in the (DD/HH/MM/SS) format
   !!  @todo see why did I do that (seems to be useless)
   !<
   !--------------------------------------------------------------------------
   function get_timestr(time)
   implicit none
   double precision :: time
   double precision :: ttt
   integer :: t(4)
   character(len=13) :: get_timestr
   
   ttt = time
   get_timestr = ' '
   
   t(1) = int(ttt/86400.d0) ;               ! days
   ttt = ttt - dble(t(1))*86400.d0
   t(2) = int(ttt/3600.d0)   ! hours
   ttt = ttt - dble(t(2)*3600)
   t(3) = int(ttt/60.d0)      ! minutes
   ttt = ttt - dble(t(3)*60)
   t(4) = int(ttt)            ! seconds

   write(get_timestr,'(i3,"d",i2,"h",i2,"m",f6.3"s")')t(1:3), ttt
   
   end function get_timestr
   !--------------------------------------------------------------------------
   !>  @author  GSF
   !!  @brief print-out the ellapsed time and the percentage of all time for 
   !!       the task specified by the string "a"
   !<
   !--------------------------------------------------------------------------
   subroutine sprint_time(a, t, Nt,tt)
   use file_mod,           ONLY: iun_out
      implicit none
      character(*) :: a
      double precision :: t, tt
      double precision :: Nt

      if (Nt>0.d0) write(iun_out,200)adjustl(a),(t),int(Nt),t/(Nt), t/tt*100.d0
   200 format ( a20, 3x, f14.3, 5x, 1i8, 5x, f12.2, 5x, f6.2 )
   end subroutine sprint_time


end module time_mod
