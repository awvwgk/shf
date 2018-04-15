module timings
   use precision, only : wp => sp ! this is enough for timings
   use iso_fortran_env, only : output_unit
   implicit none

   public  :: init_timing,start_timing_run,stop_timing_run
   public  :: start_timing,stop_timing
   public  :: prdate,prtiming

   private :: wp ! avoid name clash in main
   private :: timing_wall,timing_cpu,timing_max
   private :: start_date,start_time,start_zone,start_values
   private :: stop_date,stop_time,stop_zone,stop_values
   private :: timing

   real(wp),allocatable :: timing_wall(:)
   real(wp),allocatable :: timing_cpu(:)
   integer  :: timing_max
   
   character(len=8)  :: start_date,stop_date
   character(len=10) :: start_time,stop_time
   character(len=5)  :: start_zone,stop_zone
   integer :: start_values(8),stop_values(8)

   intrinsic :: date_and_time,system_clock,cpu_time,present

contains

subroutine prtiming(i,inmsg)
   implicit none
   integer,intent(in) :: i
   character(len=*),intent(in),optional :: inmsg
   character(len=:),allocatable :: msg
   character(len=*),parameter :: cpufmt = &
   &     '(x,a,2x,''cpu-time:'',x,f9.3,x,''seconds'')'
   character(len=*),parameter :: walfmt = &
   &     '(x,a,x,''wall-time:'',x,f9.3,x,''seconds'')'
   if (present(inmsg)) then
      msg = inmsg
   else
      msg = '*'
   endif
   write(output_unit,cpufmt) msg, &
   &   timing_cpu (i)
   write(output_unit,walfmt) msg, &
   &   timing_wall(i)
end subroutine prtiming

subroutine prdate(mode)
   implicit none
   character,intent(in) :: mode
   character(len=*),parameter :: outfmt = &
   &  '('' * '',a,x,a,''/'',a,''/'',a,x,''at'',x,a,'':'',a,'':'',a)'
   character(len=8)  :: date
   character(len=10) :: time
   character(len=5)  :: zone
   integer :: values(8)
   select case(mode)
   case('S','s')
   write(output_unit,outfmt) 'started run on', &
   &   start_date(:4),start_date(5:6),start_date(7:), &
   &   start_time(:2),start_time(3:4),start_time(5:)
   case('E','e')
   write(output_unit,outfmt) 'finished run on', &
   &   stop_date(:4),stop_date(5:6),stop_date(7:), &
   &   stop_time(:2),stop_time(3:4),stop_time(5:)
   case default
   call date_and_time(date,time,zone,values)
   write(output_unit,outfmt) 'current time:', &
   &   date(:4),date(5:6),date(7:), &
   &   time(:2),time(3:4),time(5:)
   end select
end subroutine prdate

subroutine init_timing(i)
   implicit none
   integer,intent(in) :: i
   if (allocated(timing_wall)) deallocate(timing_wall)
   if (allocated(timing_cpu))  deallocate(timing_cpu)
   timing_max = i
   allocate( timing_wall(i),  &
   &         timing_cpu (i),  &
   &         source=0.0_wp )
end subroutine init_timing

subroutine start_timing_run
   call date_and_time(start_date,start_time,start_zone,start_values)
end subroutine start_timing_run

subroutine stop_timing_run
   call date_and_time(stop_date,stop_time,stop_zone,stop_values)
end subroutine stop_timing_run

subroutine start_timing(i)
   use precision, only : wp => sp ! this is enough for timings
   implicit none
   integer,intent(in) :: i
   real(wp) :: time_cpu
   real(wp) :: time_wall
   call timing(time_cpu,time_wall)
   timing_cpu (i) = timing_cpu (i) - time_cpu
   timing_wall(i) = timing_wall(i) - time_wall
end subroutine start_timing

subroutine stop_timing(i)
   use precision, only : wp => sp ! this is enough for timings
   implicit none
   integer,intent(in) :: i
   real(wp) :: time_cpu
   real(wp) :: time_wall
   call timing(time_cpu,time_wall)
   timing_cpu (i) = timing_cpu (i) + time_cpu
   timing_wall(i) = timing_wall(i) + time_wall
end subroutine stop_timing

subroutine timing(time_cpu,time_wall)
   use precision, only : wp => sp ! this is enough for timings
   implicit none
   real(wp),intent(out) :: time_cpu
   real(wp),intent(out) :: time_wall
   integer :: time_count,time_rate,time_max
   call system_clock(time_count,time_rate,time_max)
   call cpu_time(time_cpu)
   time_wall = real(time_count)/real(time_rate)
end subroutine timing

end module timings
