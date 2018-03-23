subroutine raise(mode,message)
   character,       intent(in) :: mode
   character(len=*),intent(in) :: message
   select case(mode)
   case('W','w')
      print'(''#WARNING!'',x,a)',message
   case('E','e')
      print'(''#ERROR!'',x,a)',  message
      call terminate(1)
   end select
end subroutine raise

subroutine terminate(signal)
   integer,intent(in) :: signal
   select case(signal)
   case(0)
      stop         'normal termination of s-hf'
   case(-1)
      error stop 'external termination of s-hf'
   case default
      error stop 'abnormal termination of s-hf'
   end select
end subroutine terminate

subroutine wsigint
   print'(''recieved SIGINT, terminating...'')'
   call terminate(-1)
end subroutine wsigint

subroutine wsigterm
   print'(''recieved SIGTERM, terminating...'')'
   call terminate(-1)
end subroutine wsigterm
