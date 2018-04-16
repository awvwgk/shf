module readin

contains

subroutine rdargv(fname,wftlvl,extmode,nuhf,acc,maxiter, &
           &      diis,maxdiis,startdiis,direct_scf,mesh)
   use precision, only : wp => dp
   use typedef, only : tmesh
   use system_tools, only: rdarg,rdvar
   use density, only : rdmesh
   implicit none
   character(len=:), allocatable,intent(out) :: fname
   integer,intent(out) :: wftlvl,extmode,acc,maxiter,nuhf
   logical,intent(out) :: diis
   logical,intent(out) :: direct_scf
   integer,intent(out) :: maxdiis,startdiis
   type(tmesh),intent(out) :: mesh
   integer :: i,j,k,l
   integer :: err,idum
   logical :: exist
   logical :: getopts
   character(len=:), allocatable :: arg
   character(len=:), allocatable :: sec
   if (command_argument_count().eq.0) then
      call help
      call terminate(1)
   endif
   wftlvl = 0  ! hf
   extmode = 0 ! sp
   nuhf = -1   ! not set here
   direct_scf = .false.
   mesh % n = -1 ! not set here

   acc = 6
   maxiter = 25

!* DIIS
   diis = .false.
   maxdiis = 6
   startdiis = 1

   getopts = .true.
   j = 0
   do i = 1, command_argument_count()
      if (j.gt.0) then
         j = j-1
         cycle
      endif
      call rdarg(i,arg)
      if (arg.eq.'--') then
         getopts = .false.
         cycle
      endif
      if (getopts) then
         select case(arg)
         case('-v','-version','--version')
            call credits
            call terminate(0)
         case('-h','-help', '--help')
            call help
            call terminate(0)
         case('-a','-acc',  '--acc')
            j = 1
            call rdarg(i+1,sec)
            read(sec,*,iostat=err) idum
            if(err.ne.0) call raise('E','Bad input (acc)')
            acc = idum
            if (idum.lt.4) acc = 4
            if (idum.gt.9) acc = 9
         case('-d','-direct','--direct')
            direct_scf = .true.
         case(     '-diis', '--diis')
            diis = .true.
            if (i+1.le.command_argument_count()) then
            call rdarg(i+1,sec)
            read(sec,*,iostat=err) idum
            if (err.eq.0) then
            j = 1
            maxdiis = idum
            if (idum.lt.3) maxdiis = 3
            endif
            endif
         case('-i','-iter', '--iter')
            j = 1
            call rdarg(i+1,sec)
            read(sec,*,iostat=err) idum
            if(err.ne.0) call raise('E','Bad input (iter)')
            maxiter = idum
            if (idum.lt.1) call raise('E','iter<1 is not supported')
         case('-u','-uhf',  '--uhf')
            j = 1
            call rdarg(i+1,sec)
            read(sec,*,iostat=err) idum
            if(err.ne.0) call raise('E','Bad input (uhf)')
            nuhf = idum
            if (idum.lt.0) call raise('E','Bad input (uhf)')
         case('-o','-opt',  '--opt')
            extmode = 1
         case(     '-mp2',  '--mp2')
            wftlvl = 1
         case(     '-ccd', '--ccd')
            wftlvl = 2
         case default
            inquire(file=arg,exist=exist)
            if (exist) then
               fname = arg
            else
               if (arg(1:1).eq.'-') then
                  call raise('W','Flag unknown: '//arg)
               else
                  call raise('E','File not found: '//arg)
               endif
            endif
         case(     '-cis',  '--cis')
            extmode = 2
            wftlvl = 1
         case(     '-tdhf', '--tdhf')
            extmode = 2
            wftlvl = 2
         case(     '-dens', '--dens')
            j = 1
            call rdarg(i+1,sec)
            inquire(file=sec,exist=exist)
            if (.not.exist) call raise('E','File not found: '//sec)
            call rdmesh(sec,mesh)
         end select
      else
         inquire(file=arg,exist=exist)
         if (exist) then
            fname = arg
         else
            call raise('E','File not found: '//arg)
         endif
      endif
   enddo
end subroutine rdargv

subroutine rdinput(fname,nat,nel,nbf,at,xyz,bas,C)
   use iso_fortran_env, only : input_unit
   use precision, only : wp => dp
   use typedef, only : basis
   implicit none
!* INPUT
   character(len=:),allocatable,intent(in) :: fname
!* OUTPUT
   integer, intent(out) :: nat
   integer, intent(out) :: nel
   integer, intent(out) :: nbf
   integer, allocatable,intent(out) :: at(:)
   real(wp),allocatable,intent(out) :: xyz(:,:)
   type(basis),intent(out) :: bas
!  real(wp),allocatable,intent(out) :: zeta(:)
!  integer, allocatable,intent(out) :: aoc(:,:)
!  integer, allocatable,intent(out) :: ng(:)
!  integer, allocatable,intent(out) :: ityp(:)
   real(wp),allocatable,intent(out) :: C(:,:)
!* TEMPORARY
   character(len=72) :: line
   character(len=2)  :: chdum
   integer  :: i,j,k,ibf,iz,err,id,ii,idum
   logical  :: stdin
   logical  :: exist
   real(wp) :: x,y,z,dum,zeta

   id = 42
   stdin = .false.

   if (.not.allocated(fname)) then
!* in case you really want to you can provide everything by hand
      call raise('W','No input file given, reading STDIN')
      stdin = .true.
      id = input_unit
   else
!* rdargv should have catched this errors already, but to be sure
!  we check again, if an input file is provided and if it exist
      inquire(file=fname,exist=exist)
      if (.not.exist) call raise('E','File: '//fname//' not found')
      open(id,file=fname,status='old')
   endif
  
   if (stdin) print'(''<nat> <nel> <nbf>'',/,''> '',$)'
   read(id,*,iostat=err) nat,nel,nbf
   if (err.ne.0) call raise('E','Error while reading input file')
   if (nat.lt.1) call raise('E','No atoms')
   if (nel.lt.1) call raise('E','No electrons')
   if (nbf.lt.1) call raise('E','No basis functions')
!* Now make some space for your arrays
   bas % n = nbf
   allocate( xyz(3,nat), bas % zeta(nbf), &
   &         source = 0.0_wp )
   allocate( at(nat), bas % aoc(2,nat), &
   &         source = 0 )
   allocate( bas % ityp(nbf), source = 1 )
   allocate( bas % ng(nbf), source = 6 )
!* actually nbf is redundant, but lets check if the user has given 
!  the right value
   k=0
   do i = 1, nat !* loop over all atoms
      if (stdin) print'(''<x> <y> <z> <iat> <ibf>'',/,i0'' > '',$)',i
      read(id,*,iostat=err) x,y,z,iz,ibf
      if(err.ne.0) call raise('E','Error while reading input file')
      !* and save them into the array
      xyz(:,i) = (/x,y,z/)
      if(iz.lt.0) call raise('E','negative nuclear charge')
      at(i) = iz ! no fractional nuclear charge supported
      !* now check for the exponents
      if((k+ibf).gt.nbf) call raise('E','to many exponents')
      do j = 1, ibf
         if (stdin) print'(''<zeta> [<ityp> <ng>]'',/,i0,x,i0,'' > '',$)',i,j
         read(id,'(a)',iostat=err) line
         if(err.ne.0) call raise('E','Error while reading exponents')
         read(line,*,iostat=err) dum,chdum,idum
         if(err.ne.0) then
            read(line,*,iostat=err) zeta
            bas % zeta(k+j) = zeta
            if(err.ne.0) call raise('E','Error while reading exponents')
            bas % ng(j) = 6
            bas % ityp(j) = 1
         else
            bas % zeta(k+j) = dum
            if(idum.lt.1) call raise('E','negative primitive count')
            if(idum.gt.6) call raise('W','>6 primitives is not supported')
            bas % ng(j) = idum
            select case(chdum)
            case('1s'); bas % ityp(j) =  1
            case('2s'); bas % ityp(j) =  2
            case('3s'); bas % ityp(j) =  3
            case('4s'); bas % ityp(j) =  4
            case('5s'); bas % ityp(j) =  5
         !  case('2p'); bas % ityp(j) =  6
         !  case('3p'); bas % ityp(j) =  7
         !  case('4p'); bas % ityp(j) =  8
         !  case('5p'); bas % ityp(j) =  9
         !  case('3d'); bas % ityp(j) = 10
         !  case('4d'); bas % ityp(j) = 11
         !  case('5d'); bas % ityp(j) = 12
         !  case('4f'); bas % ityp(j) = 13
         !  case('5f'); bas % ityp(j) = 14
         !  case('5g'); bas % ityp(j) = 15
            case default
               call raise('W',chdum//' functions are not implemented, using 1s')
               bas % ityp(j) = 1
            end select
         endif
      enddo
      !* save some informations on the basis functions
      bas % aoc(1,i) = k+1
      bas % aoc(2,i) = k+ibf
      k=k+ibf
   enddo
!* we already cover part of this error, but lets check anyway
   if(k.ne.nbf) call raise('E','number of basis function mismatch')

!* restart section
   if (stdin) print'(''write <restart> to provide start orbitals'',/,''> '',$)'
   read(id,'(a)',iostat=err) line
   if (err.eq.0) then
      if (index(line,'restart').ne.0) then
         allocate( C(nbf,nbf), source = 0.0_wp )
         i = 0
         do
            i = i + 1
            if (i.gt.nbf) exit
            if (stdin) print'(''which MO?'',x,i0,x,''left'',/,''> '',$)', nbf-i+1
            read(id,*,iostat=err) ii
            if (err.ne.0) then
               call raise('W','could not read restart')
               deallocate(C)
               exit
            endif
            if ((ii.gt.0).and.(ii.le.nbf)) then
               if (stdin) print'(''please provide coeffients'',/,i0,x''> '',$)', ii
               read(id,*,iostat=err) (C(j,ii),j=1,nbf)
               if (err.ne.0) then
                  call raise('W','could not read restart')
                  deallocate(C)
                  exit
               endif
            endif
         enddo
      endif
      if (stdin) print'(''Thank you for your patience'')'
   endif
   close(id)

end subroutine rdinput

end module readin
