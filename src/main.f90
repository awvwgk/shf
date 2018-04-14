!!================================================================= 1803
!!     _______           _______      |   This program is free software: 
!!    |  ____ \|\     /||  ____ \     | you can redistribute it and/or 
!!    | |____\/| |___| || |__  \/     | modify it under the terms of the 
!!    |_____  \|  ___  ||  __|        | GNU General Public License as 
!!    /\____) || |   | || |           | published by the Free Software 
!!    \_______/|/     \||/            | Foundation, either version 3 of 
!!  © 2017,2018 by Sebastian Ehlert   | the License, or (at your option) 
!!------------------------------------+ any later version.
!!      This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!    GNU General Public License for more details.
!!================================================================= 1803
program shf
   use precision, only : wp => dp

!* input/output modules
   use readin
   use printout
   use timings

!* quantum mechanical calculations
   use ints
   use guess
   use scf
   use mbpt

!* some analysis
   use pop

!* optimisiation stuff
   use opt

   implicit none

!* configuration data
   character(len=:),allocatable :: fname
   character(len=3) :: chacc
   integer  :: wftlvl,extmode,acc
   real(wp) :: ethr,pthr,gthr
   logical  :: first,brsym
   integer  :: maxiter,err

!* molecule data
   integer  :: nat,nel,nbf,nuhf,nocc,nalp,nbet,ndum
   integer, allocatable :: at(:)
   real(wp),allocatable :: xyz(:,:)
   real(wp),allocatable :: zeta(:)
   integer, allocatable :: aoc(:,:)
   real(wp),allocatable :: g(:,:)
   real(wp) :: e,ecorr

!* integrals and stuff
   integer, allocatable :: ng(:),ityp(:)
   real(wp),allocatable :: eri(:)
   real(wp),allocatable :: S(:,:),V(:,:),T(:,:)
   real(wp),allocatable :: X(:,:)
   real(wp),allocatable :: F(:,:),H(:,:),P(:,:),C(:,:)
   real(wp),allocatable :: eps(:),z(:)
!  real(wp),allocatable :: eri(:,:,:,:)
   real(wp),allocatable :: Fa(:,:),Pa(:,:),Ca(:,:)
   real(wp),allocatable :: Fb(:,:),Pb(:,:),Cb(:,:)
   real(wp),allocatable :: epsa(:),epsb(:)

!* DIIS
   logical :: diis
   integer :: maxdiis,startdiis

!* OMP
   common /proc/ nproc
   integer nproc
   integer tid,omp_get_num_threads,omp_get_thread_num

!* signal handler
   external wSIGTERM,wSIGINT
   call signal(2,wSIGINT)
   call signal(15,wSIGTERM)

   call start_timing_run
   call init_timing(10)
   call start_timing(1)

!* init section for intermediate logicals
   first = .true.
   brsym = .false.

!* get command line argument, most of the stuff gets initialized here
   call rdargv(fname,wftlvl,extmode,nuhf,acc,maxiter, &
        &      diis,maxdiis,startdiis)

!* print some fancy banner
   call banner
   call prdate('S')

!$omp parallel private(tid)
      tid = omp_get_thread_num()
      if (tid.eq.0) then
         nproc = omp_get_num_threads()
         !print'(72(''-''))'
         print'('' # OMP threads:'',i25)',nproc
         !print'(72(''-''))'
      end if
!$omp end parallel 

!* set thresholds based on accuracy
   ethr = 10._wp**(-acc)
   pthr = 10._wp**(-acc+1)
   gthr = 10._wp**(-acc+1)
   write(chacc,'(i3)') acc
   chacc = trim(adjustl(chacc))

!* read the input file
   call rdinput(fname,nat,nel,nbf,at,xyz,zeta,aoc,ng,ityp,C)

!* sanity check
   if (nuhf.eq.-1) then ! not set before, set here
      nuhf = mod(nel,2)
   !* first check if unrestricted calculation is necessary
      if (nuhf.eq.0) then
         nocc = nel/2
      else
         brsym = .true. ! allow breaking of symmetry
         ndum = nel - nuhf
         nalp = ndum/2 + nuhf
         nbet = ndum/2
      endif
   else ! oh, we got user input, double check this!!
      if((mod(nel,2).eq.0).and.(mod(nuhf,2).ne.0)) &
      &   call raise('E','Even number of electrons and odd Nα-Nβ')
      if((mod(nel,2).ne.0).and.(mod(nuhf,2).eq.0)) &
      &   call raise('E','Odd number of electrons and even Nα-Nβ')
      if(nuhf.gt.nel) call raise('E','Nα-Nβ>Nel, this is impossible')
      brsym = .true. ! allow breaking of symmetry
      ndum = nel - nuhf
      nalp = ndum/2 + nuhf
      nbet = ndum/2
   endif

!* print some information on the calculation
   call prinput(fname,nat,nel,nocc,nalp,nbet,nbf,at,xyz,zeta,aoc,ng,brsym)

   select case(extmode)
!* extmode sets:
!  0 -> singel point calculation
!  1 -> geometry optimiazation
!* wftlvl sets:
!  0 -> HF-SCF
!  1 -> MP2
!  2 -> CCSD
!  3 -> CCSD(T)

   case(0) !* scf calculation

   allocate( S(nbf,nbf),V(nbf,nbf),T(nbf,nbf), &
   &         eri(nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2), &
   &         source = 0.0_wp )
   call start_timing(2)
   call integrals(nat,nbf,at,xyz,zeta,aoc,ng,ityp,S,V,T,eri)
   call stop_timing(2)
!* debug
!  call prmat(S,nbf,nbf,name='S')
!  call prmat(V,nbf,nbf,name='V')
!  call prmat(T,nbf,nbf,name='T')

   if (brsym) then !* UHF

!* no restart possible, yet
   allocate( X(nbf,nbf),H(nbf,nbf),z(nbf),g(3,nat), &
   &         Pa(nbf,nbf),Fa(nbf,nbf),Ca(nbf,nbf),epsa(nbf), &
   &         Pb(nbf,nbf),Fb(nbf,nbf),Cb(nbf,nbf),epsb(nbf), &
   &         source = 0.0_wp )

!* set up the calculation
   H=T+V
!  call prmat(H,nbf,nbf,name='H') !* debug

   call build_orthonormalizer(nbf,S,X,err)
   if (err.ne.0) call raise('E','building orthonormalizer failed')
!  call prmat(X,nbf,nbf,name='S^-1/2') !* debug

!* initial guess density (currently hcore guess)
!  if (first) then ! currently no restart possible for uhf case
      if (nalp.eq.nbet) then
         print'('' * doing hcore guess for alpha MOs'')'
      else
         print'('' * doing hcore guess'')'
      endif
      call hcore_guess(nbf,nalp,H,X,Fa,Ca)
      if (nalp.eq.nbet) then
         print'('' * doing really bad orthormalizer guess for beta MOs'')'
         Cb = X ! use really bad orthonormalizer guess to break symmetry
      else
         call hcore_guess(nbf,nbet,H,X,Fb,Cb)
      endif
!  endif

   call start_timing(3)
   call hf(nat,nbf,nalp,nbet,at,xyz,ethr,pthr,first, &
        &  chacc,maxiter,diis,maxdiis,startdiis, &
        &  S,V,T,X,Pa,Pb,H,Fa,Fb,Ca,Cb,eri,epsa,epsb,e)
   call stop_timing(3)
   call prenergy(nat,nbf,nalp,nbet,at,xyz,chacc, &
        &        S,V,T,eri,H,Fa,Fb,Pa,Pb,Ca,Cb,epsa,epsb)

   print'(a)'
   print'('' * doing Mulliken population analysis'')'
   call mulliken(nat,nbf,at,aoc,S,Pa,Pb,z)
   call prchrg(nat,at,z,chacc)

   print'(a)'
   print'('' * doing Löwdin population analysis'')'
   call loewdin(nat,nbf,at,aoc,S,Pa,Pb,X,z)
   call prchrg(nat,at,z,chacc)

!* that's it, no MP2 or CC for you in unrestricted cases
   if (wftlvl.gt.0) call raise('E','Only HF supported for unrestricted case')

   else !* RHF

!* check if we got some user input here
   if(allocated(C)) then
   first = .false.
   allocate( X(nbf,nbf), &
   &         P(nbf,nbf),H(nbf,nbf),F(nbf,nbf), &
   &         eps(nbf),z(nbf),g(3,nat), &
   &         source = 0.0_wp )
   else
   allocate( X(nbf,nbf), &
   &         P(nbf,nbf),H(nbf,nbf),F(nbf,nbf),C(nbf,nbf), &
   &         eps(nbf),z(nbf),g(3,nat), &
   &         source = 0.0_wp )
   endif

!* set up the calculation
   H=T+V
!  call prmat(H,nbf,nbf,name='H') !* debug

   call build_orthonormalizer(nbf,S,X,err)
   if (err.ne.0) call raise('E','building orthonormalizer failed')
!  call prmat(X,nbf,nbf,name='S^-1/2') !* debug

!* initial guess density (currently hcore guess)
   if (first) then
      print'('' * doing hcore guess'')'
      call hcore_guess(nbf,nocc,H,X,F,C)
   endif

   call start_timing(3)
   call hf(nat,nbf,nocc,at,xyz,ethr,pthr,first, &
        &  chacc,maxiter,diis,maxdiis,startdiis, &
        &  S,V,T,X,P,H,F,C,eri,eps,e)
   call stop_timing(3)
   call prenergy(nat,nbf,nocc,nocc,at,xyz,chacc, &
        &        S,V,T,eri,H,F,F,P,P,C,C,eps,eps)
   call onetrafo(nbf,F,C)
   call prmat(F,nbf,nbf,'Final Fock matrix')
   print'(a)'
   print'('' * dumping restart file'')'
   call prrestart('restart',nat,nel,nbf,at,xyz,zeta,aoc,ng,ityp,C,chacc)

   print'(a)'
   print'('' * doing Mulliken population analysis'')'
   call mulliken(nat,nbf,at,aoc,S,P,P,z)
   call prchrg(nat,at,z,chacc)

   print'(a)'
   print'('' * doing Löwdin population analysis'')'
   call loewdin(nat,nbf,at,aoc,S,P,P,X,z)
   call prchrg(nat,at,z,chacc)

!* experimental direct scf
!  call hf(nat,nbf,nocc,at,xyz,ethr,pthr,first, &
!       &  chacc,maxiter,diis,maxdiis,startdiis, &
!       &  zeta,aoc,ng,ityp,S,X,P,H,F,C,eps,e)
   
!* second order Møller-Plesset many-body pertubation theory
   if (wftlvl.eq.1) then !* MP2
   call mp2(nbf,nocc,C,eri,eps,ecorr,chacc)
   e = e + ecorr
   endif !* MP2

!  if (wftlvl.ge.2) then !* CCSD
!  call ccsd(nbf,nocc,nvir,F,eri,ethr,chacc,maxiter,ecorr)
!  if (wftlvl.eq.3) then !* CCSD(T)
!  call ccpt(...)
!  endif !* CCSD(T)
!  endif !* CCSD
   
   endif
   
!* Geometry optimization
   case(1)
   if(allocated(C)) then
   first = .false.
   allocate( S(nbf,nbf),V(nbf,nbf),T(nbf,nbf),X(nbf,nbf), &
   &         P(nbf,nbf),H(nbf,nbf),F(nbf,nbf), &
   &         eps(nbf),z(nbf),g(3,nat), &
   &         eri(nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2), &
   &         source = 0.0_wp )
   else
   allocate( S(nbf,nbf),V(nbf,nbf),T(nbf,nbf),X(nbf,nbf), &
   &         P(nbf,nbf),H(nbf,nbf),F(nbf,nbf),C(nbf,nbf), &
   &         eps(nbf),z(nbf),g(3,nat), &
   &         eri(nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2), &
   &         source = 0.0_wp )
   endif

   call geoopt_sd(nat,nbf,nocc,at,xyz,zeta,aoc,ng,ityp,first,ethr,pthr,gthr, &
        &         chacc,maxiter,diis,maxdiis,startdiis, &
        &         S,V,T,X,P,H,F,C,eri,eps,e,g)
   call prrestart('geoopt.in',nat,nel,nbf,at,xyz,zeta,aoc,ng,ityp,C,chacc)

   end select
   print'(a)'

!* done, print some result
   print'(72(''=''))'
   print'('' TOTAL ENERGY'',f57.'//chacc//')',e
   print'(72(''=''))'

!* give some statistics on timings and stuff
   print'(a)'
   call stop_timing_run
   call stop_timing(1)
   call prdate('E')
   call prtiming(1)
!  call prtiming(2,'int')
!  call prtiming(3,'scf')
   print'(a)'

   call terminate(0)
end program shf
