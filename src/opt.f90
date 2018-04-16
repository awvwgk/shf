module opt
   use precision, only : wp => dp
   implicit none

   real(wp),parameter :: eta=0.01_wp

contains

subroutine geoopt_sd(nat,nbf,nocc,at,xyz,zeta,aoc,ng,ityp,first, &
           &         ethr,pthr,gthr,acc,maxiter,diis,maxdiis,startdiis, &
           &         S,V,T,X,P,H,F,C,eri,eps,e,g)
   use precision, only : wp => dp
   use printout, only : prmat
   use ints, only : integrals
   use scf, only : hf
   implicit none
   integer, intent(in)    :: nat,nbf,nocc
   integer, intent(in)    :: at(nat)
   real(wp),intent(inout) :: xyz(3,nat)
   real(wp),intent(in)    :: zeta(nbf)
   integer, intent(in)    :: aoc(2,nat)
   integer, intent(in)    :: ng(nbf)
   integer, intent(in)    :: ityp(nbf)
   logical, intent(inout) :: first,diis
   real(wp),intent(in)    :: ethr,pthr,gthr
   integer, intent(in)    :: maxiter,maxdiis,startdiis
   real(wp),intent(inout) :: C(nbf,nbf)
   real(wp),intent(out)   :: S(nbf,nbf),V(nbf,nbf),T(nbf,nbf)
   real(wp),intent(in)    :: X(nbf,nbf)
   real(wp),intent(in)    :: H(nbf,nbf)
   real(wp),intent(out)   :: F(nbf,nbf)
   real(wp),intent(out)   :: P(nbf,nbf)
   real(wp),intent(out)   :: eri(nbf*(nbf+1)*(nbf*(nbf+1)/2+1)/2),eps(nbf)
   real(wp),intent(out)   :: e,g(3,nat)
   character(len=*),intent(in) :: acc

   integer  :: iter
   real(wp) :: enuc,eold,gnorm

   print'(a)'
   print'(72(''=''))'
   print'('' steepest decent geometry optimization'')'

   call integrals(nat,nbf,at,xyz,zeta,aoc,ng,ityp,S,V,T,eri)
!* debug
!  call prmat(S,nbf,nbf,name='S')
!  call prmat(V,nbf,nbf,name='V')
!  call prmat(T,nbf,nbf,name='T')

   call hf(nat,nbf,nocc,at,xyz,ethr,pthr,first, &
        &  acc,maxiter,diis,maxdiis,startdiis, &
        &  S,V,T,X,P,H,F,C,eri,eps,eold)
   iter = 0
   opt: do
      iter = iter+1
      if (iter.gt.maxiter) call raise('E','optimization did not converge')
      call rhf_numgrad(nat,nbf,nocc,at,xyz,zeta,aoc,ng,ityp,ethr,C,g)
      xyz = xyz + eta * g
      call hf(nat,nbf,nocc,at,xyz,ethr,pthr,first, &
           &  acc,maxiter,diis,maxdiis,startdiis, &
           &  S,V,T,X,P,H,F,C,eri,eps,e)
      gnorm = sqrt(sum( g**2 ))
      if (gnorm.gt.500.0_wp) call raise('W','|G|>500, something went wrong!')
      print'(a)'
      print'(''=[ITER]=========[E(SCF)]===========[|∇E|]''31(''=''))'
      print'(x,i5,x,f16.'//acc//',x,f16.'//acc//')',iter,e,gnorm
      if(abs(e-eold).lt.ethr) exit opt
      eold = e
   enddo opt

end subroutine geoopt_sd

!subroutine rhf_grad(nat,nbf,nocc,at,xyz,zeta,aoc,ethr,C,g)
!   use precision, only : wp => dp
!   implicit none
!   include 'interface_energy'
!   include 'interface_integrals'
!   integer,intent(in)  :: nat,nbf,nocc
!   integer,intent(in)  :: at(nat)
!   real(wp),intent(in)  :: xyz(3,nat)
!   real(wp),intent(in)  :: zeta(nbf)
!   integer,intent(in)  :: aoc(2,nat)
!   real(wp),intent(in)  :: ethr
!   real(wp),intent(in)  :: C(nbf,nbf)
!   real(wp),intent(out) :: g(3,nat)
!
!   integer :: i,ii
!   real(wp) :: er,el,enuc
!   real(wp),allocatable :: S(:,:),V(:,:),T(:,:)
!   real(wp),allocatable :: F(:,:),H(:,:),P(:,:)
!   real(wp),allocatable :: eri(:,:,:,:)
!   real(wp),allocatable :: xyzdup(:,:)
!
!   allocate( xyzdup(3,nat), source = xyz )
!   allocate( S(nbf,nbf),V(nbf,nbf),T(nbf,nbf), &
!   &         P(nbf,nbf),H(nbf,nbf),F(nbf,nbf), &
!   &         eri(nbf,nbf,nbf,nbf), &
!   &         source = 0.0_wp )
!
!   do i = 1, nat
!      do ii = 1, 3
!         !* move atom to the right
!         xyzdup(ii,i) = xyz(ii,i) + ethr
!         call nnrep(enuc,nat,at,xyzdup)
!         call integrals(nat,nbf,at,xyzdup,zeta,aoc,S,V,T,eri)
!         call build_density(nbf,nocc,P,C)
!         call build_fock(nbf,H,F,P,eri)
!         er = enuc+escf(H,F,P,nbf)
!         !* move atom to the left
!         xyzdup(ii,i) = xyz(ii,i) - ethr
!         call nnrep(enuc,nat,at,xyzdup)
!         call integrals(nat,nbf,at,xyzdup,zeta,aoc,S,V,T,eri)
!         call build_density(nbf,nocc,P,C)
!         call build_fock(nbf,H,F,P,eri)
!         el = enuc+escf(H,F,P,nbf)
!         !* reset atom
!         xyzdup(ii,i) = xyz(ii,i)
!         !* get numerical gradient
!         g(ii,i) = g(ii,i) + 0.5_wp*(er-el)/ethr
!      enddo
!   enddo
!
!end subroutine rhf_grad

subroutine rhf_numgrad(nat,nbf,nocc,at,xyz,zeta,aoc,ng,ityp,ethr,C,g)
   use precision, only : wp => dp
   use ints, only : integrals
   use scf, only : nnrep,build_density,build_fock,escf
   implicit none
   integer, intent(in)  :: nat,nbf,nocc
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: zeta(nbf)
   integer, intent(in)  :: aoc(2,nat)
   integer, intent(in)  :: ng(nbf)
   integer, intent(in)  :: ityp(nbf)
   real(wp),intent(in)  :: ethr
   real(wp),intent(in)  :: C(nbf,nbf)
   real(wp),intent(out) :: g(3,nat)

   integer  :: i,ii
   real(wp) :: er,el,enuc
   real(wp),allocatable :: S(:,:),V(:,:),T(:,:)
   real(wp),allocatable :: F(:,:),H(:,:),P(:,:)
   real(wp),allocatable :: eri(:)
   real(wp),allocatable :: xyzdup(:,:)

   allocate( xyzdup(3,nat), source = xyz )
   allocate( S(nbf,nbf),V(nbf,nbf),T(nbf,nbf), &
   &         P(nbf,nbf),H(nbf,nbf),F(nbf,nbf), &
   &         eri(nbf*(nbf+1)*(nbf*(nbf+1)/2+1)/2), &
   &         source = 0.0_wp )

   g = 0.0_wp

   do i = 1, nat
      do ii = 1, 3
         !* move atom to the right
         xyzdup(ii,i) = xyz(ii,i) + ethr
         call nnrep(enuc,nat,at,xyzdup)
         call integrals(nat,nbf,at,xyzdup,zeta,aoc,ng,ityp,S,V,T,eri)
         call build_density(nbf,nocc,P,C)
         call build_fock(nbf,H,F,P,eri)
         er = enuc+escf(H,F,P,nbf)
         !* move atom to the left
         xyzdup(ii,i) = xyz(ii,i) - ethr
         call nnrep(enuc,nat,at,xyzdup)
         call integrals(nat,nbf,at,xyzdup,zeta,aoc,ng,ityp,S,V,T,eri)
         call build_density(nbf,nocc,P,C)
         call build_fock(nbf,H,F,P,eri)
         el = enuc+escf(H,F,P,nbf)
         !* reset atom
         xyzdup(ii,i) = xyz(ii,i)
         !* get numerical gradient
         g(ii,i) = g(ii,i) + 0.5_wp*(er-el)/ethr
      enddo
   enddo

end subroutine rhf_numgrad

subroutine rhf_numgexp(nat,nbf,nocc,at,xyz,zeta,aoc,ng,ityp,ethr,C,g)
   use precision, only : wp => dp
   use ints, only : integrals
   use scf, only : build_density,build_fock,escf
   implicit none
   integer, intent(in)  :: nat,nbf,nocc
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: zeta(nbf)
   integer, intent(in)  :: aoc(2,nat)
   integer, intent(in)  :: ng(nbf)
   integer, intent(in)  :: ityp(nbf)
   real(wp),intent(in)  :: ethr
   real(wp),intent(in)  :: C(nbf,nbf)
   real(wp),intent(out) :: g(nbf)

   integer  :: i,ii
   real(wp) :: er,el,enuc
   real(wp),allocatable :: S(:,:),V(:,:),T(:,:)
   real(wp),allocatable :: F(:,:),H(:,:),P(:,:)
   real(wp),allocatable :: eri(:)
   real(wp),allocatable :: zetadup(:)

   allocate( zetadup(nbf), source = zeta )
   allocate( S(nbf,nbf),V(nbf,nbf),T(nbf,nbf), &
   &         P(nbf,nbf),H(nbf,nbf),F(nbf,nbf), &
   &         eri(nbf*(nbf+1)*(nbf*(nbf+1)/2+1)/2), &
   &         source = 0.0_wp )

   do i = 1, nat
      !* set higher exponent
      zetadup(i) = zeta(i) + ethr
      call integrals(nat,nbf,at,xyz,zetadup,aoc,ng,ityp,S,V,T,eri)
      call build_density(nbf,nocc,P,C)
      call build_fock(nbf,H,F,P,eri)
      er = escf(H,F,P,nbf)
      !* set lower exponent
      zetadup(i) = zeta(i) - ethr
      call integrals(nat,nbf,at,xyz,zetadup,aoc,ng,ityp,S,V,T,eri)
      call build_density(nbf,nocc,P,C)
      call build_fock(nbf,H,F,P,eri)
      el = escf(H,F,P,nbf)
      !* reset exponent
      zetadup(i) = zeta(i)
      !* get numerical gradient
      g(i) = g(i) + 0.5_wp*(er-el)/ethr
   enddo

end subroutine rhf_numgexp

end module opt
