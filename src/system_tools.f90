module system_tools

contains

subroutine rdarg(i,arg)
   integer,intent(in) :: i
   character(len=:),allocatable,intent(out) :: arg
   integer :: l,err
   if (allocated(arg)) deallocate(arg)
   call get_command_argument(i,length=l,status=err)
   if (err.gt.0) call raise('E','Command argument corrupted')
   allocate( character(len=l) :: arg, stat=err )
   if (err.ne.0) call raise('E','could not be allocated')
   call get_command_argument(i,arg,status=err)
   if (err.gt.0) call raise('E','Command argument corrupted')
end subroutine rdarg

subroutine rdvar(name,var)
   character(len=*),intent(in) :: name
   character(len=:),allocatable,intent(out) :: var
   integer :: l,err
   if (allocated(var)) deallocate(var)
   call get_environment_variable(name,length=l,status=err)
   if (err.gt.0) call raise('E','System variable unassigned')
   allocate( character(len=l) :: var, stat=err )
   if (err.ne.0) call raise('E','could not be allocated')
   call get_environment_variable(name,var,status=err)
   if (err.gt.0) call raise('E','System variable corrupted')
end subroutine rdvar

end module system_tools
