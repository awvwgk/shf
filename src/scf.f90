module scf

!* use this nice overloading feature of FORTRAN to
!  combine all SCF variants in one name
   interface hf
   module procedure rhf_conventional
   module procedure uhf_conventional
!  module procedure rhf_direct
!  module procedure uhf_direct
   end interface hf

contains

pure subroutine nnrep(enuc,nat,at,xyz)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: enuc

   integer :: i,j

   enuc = 0.0_wp

   do i = 1, nat
      do j = 1, i-1 !* no double counting, please
         enuc = enuc + dble(at(i)*at(j)) / &
         & sqrt( sum( (xyz(1:3,j) - xyz(1:3,i) )**2 ) )
      enddo
   enddo

end subroutine nnrep

pure subroutine build_orthonormalizer(nbf,S,X,err)
   use precision, only : wp => dp
   use lapack95,  only : syev
   use blas95,    only : gemm
   use misc,      only : idx
   implicit none
   integer, intent(in)  :: nbf
   real(wp),intent(in)  :: S(nbf,nbf)
   real(wp),intent(out) :: X(nbf,nbf)
   integer, intent(out) :: err

   real(wp) :: d(nbf),tmp1(nbf,nbf),tmp2(nbf,nbf)
   logical  :: pidx(nbf,nbf)
   integer  :: i,j,l,ij

   X = 0.0_wp
   tmp1 = S

!* lapack diagonalizer
   call syev(tmp1,d,'v','u',err)
   if (err.ne.0) return
!* quick build of diagonal matrix
   forall (i=1:nbf) X(i,i) = 1.0_wp/sqrt(d(i))
!* construct orthonormalizer
   call gemm(tmp1,X,tmp2)
   call gemm(tmp2,tmp1,X,transb='t')

end subroutine build_orthonormalizer

!* scf energy for given density + fockian
pure function escf(H,F,P,nbf) result(e)
   use precision, only : wp => dp
   use blas95,    only : gemm
   implicit none
   integer, intent(in) :: nbf
   real(wp),intent(in) :: H(nbf,nbf)
   real(wp),intent(in) :: F(nbf,nbf)
   real(wp),intent(in) :: P(nbf,nbf)

   real(wp) :: e,tmp(nbf,nbf)
   integer  :: i,j

   call gemm(P,(H+F),tmp)
   e = sum((/(tmp(i,i),i=1,nbf)/))

end function escf

pure subroutine build_fock(nbf,H,F,P,eri)
   use precision, only : wp => dp
   use misc,      only : idx
   implicit none
   integer, intent(in)  :: nbf
   real(wp),intent(in)  :: H(nbf,nbf)
   real(wp),intent(in)  :: P(nbf,nbf)
   real(wp),intent(in)  :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)
   real(wp),intent(out) :: F(nbf,nbf)

   integer  :: i,j,k,l,ij,kl,ijkl,ik,jl,ikjl

   F = H
   do i = 1, nbf
      do j = 1, i
         ij = idx(i,j)
         do k = 1, nbf
            ik = idx(i,k)
            do l = 1, nbf
               kl = idx(k,l)
               jl = idx(j,l)
               ijkl = idx(ij,kl)
               ikjl = idx(ik,jl)
               F(j,i) = F(j,i) + P(l,k) * ( 2*eri(ijkl) - eri(ikjl) )
            enddo
         enddo
         F(i,j) = F(j,i)
      enddo
   enddo

end subroutine build_fock

pure subroutine build_density(nbf,nocc,P,C)
   use precision, only : wp => dp
   use blas95,    only : gemm
   implicit none
   integer, intent(in)  :: nbf,nocc
   real(wp),intent(in)  :: C(nbf,nbf)
   real(wp),intent(out) :: P(nbf,nbf)

   integer :: i,j

   call gemm(C(:,:nocc),C(:,:nocc),P,transb='t')

end subroutine build_density

pure subroutine roothaan_hall(nbf,X,F,C,eps,err)
   use precision, only : wp => dp
   use lapack95,  only : syev
   use blas95,    only : gemm
   implicit none
   integer, intent(in)  :: nbf
   real(wp),intent(in)  :: F(nbf,nbf)
   real(wp),intent(in)  :: X(nbf,nbf)
   real(wp),intent(out) :: C(nbf,nbf)
   real(wp),intent(out) :: eps(nbf)
   integer, intent(out) :: err

   real(wp) :: w(4*nbf),tmp(nbf,nbf)
   integer  :: i,j,l

   call gemm(X,F,C,transa='t')
   call gemm(C,X,tmp)
!* lapack diagonalizer
   call syev(tmp,eps,'v','u',err)
   if (err.ne.0) return
   call gemm(X,tmp,C)

end subroutine roothaan_hall

subroutine rhf_conventional &
           &  (nat,nbf,nocc,at,xyz,ethr,pthr,first, &
           &   acc,maxiter,ldiis,maxdiis,startdiis, &
           &   S,V,T,X,P,H,F,C,eri,eps,e)
   use iso_fortran_env, only : id => output_unit
   use precision, only : wp => dp
   use diis,      only : build_diis,diis_fock
   implicit none
   integer, intent(in)    :: nat
   integer, intent(in)    :: nbf
   integer, intent(in)    :: nocc
   integer, intent(in)    :: at(nat)
   real(wp),intent(in)    :: xyz(3,nat)
   real(wp),intent(in)    :: ethr
   real(wp),intent(in)    :: pthr
   logical, intent(inout) :: first
   logical, intent(inout) :: ldiis
   integer, intent(in)    :: maxiter
   integer, intent(in)    :: maxdiis
   integer, intent(in)    :: startdiis
   real(wp),intent(in)    :: X(nbf,nbf)
   real(wp),intent(in)    :: H(nbf,nbf)
   real(wp),intent(out)   :: P(nbf,nbf)
   real(wp),intent(inout) :: F(nbf,nbf)
   real(wp),intent(inout) :: C(nbf,nbf)
   real(wp),intent(out)   :: eps(nbf)
   real(wp),intent(out)   :: e
   real(wp),intent(in)    :: S(nbf,nbf)
   real(wp),intent(in)    :: V(nbf,nbf)
   real(wp),intent(in)    :: T(nbf,nbf)
   real(wp),intent(in)    :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)
   character(len=*),intent(in) :: acc

   integer  :: iter,err,i,j,ij
   real(wp) :: enuc,eold,rmsd
   real(wp) :: ekin,epot,ec,ex
   real(wp),allocatable :: P_save(:,:)

!* DIIS
   real(wp) :: emax
   real(wp),allocatable :: emat(:,:,:),F_save(:,:,:),cdiis(:),SS(:,:)

   write(id,'(a)')
   write(id,'(72(''-''))')
   write(id,'('' Restricted Hartree-Fock SCF calculation'')')

   allocate( P_save(nbf,nbf), source=0.0_wp )

   if (ldiis) then
   allocate( emat(nbf,nbf,maxdiis),F_save(nbf,nbf,maxdiis), &
   &         cdiis(maxdiis+1), source=0.0_wp )
   endif

   call nnrep(enuc,nat,at,xyz)

   call build_density(nbf,nocc,P,C)
!* debug
!  call prmat(F,nbf,nbf,name='F')
!  call prmat(P,nbf,nbf,name='P')
!  call prmat(C,nbf,nbf,name='C')

   eold = enuc+escf(H,F,P,nbf)
   iter = 0

   write(id,'(1(''-'')''[ITER]'')',   advance='no')
   write(id,'(9(''-'')''[E(SCF)]'')', advance='no')
   write(id,'(8(''-'')''[ΔE(SCF)]'')',advance='no')
   write(id,'(7(''-'')''[RMS(P)]'')', advance='no')
   if (ldiis) then
      write(id,'(5(''-'')''[max[F,P]]'')',advance='no')
      write(id,'(1(''-''))')
   else
      write(id,'(15(''-''))')
   endif
   write(id,'(x,i5,x,f16.'//acc//')'),iter,eold

   scf: do
      iter = iter+1
      if (iter.gt.maxiter) call raise('E','SCF did not converge')
      call build_fock(nbf,H,F,P,eri)

      if (ldiis) then !* DIIS
      call build_diis(nbf,iter,maxdiis,S,F,P,F_save,emat,emax,cdiis,ldiis)
      if (ldiis) then
      if (iter.ge.startdiis) then !* start DIIS
         if (iter.eq.startdiis) write(id,'('' * starting DIIS'')')
         call diis_fock(nbf,iter,maxdiis,F,F_save,cdiis)
      endif !* start DIIS
      else
         write(id,'('' * shutting down DIIS'')')
         deallocate( F_save,emat,cdiis )
      endif
      endif !* DIIS
      
      call roothaan_hall(nbf,X,F,C,eps,err)
      if(err.ne.0) call raise('E','solving Roothaan-Hall eq. failed')
      P_save = P
      call build_density(nbf,nocc,P,C)
      e = enuc+escf(H,F,P,nbf)
      rmsd = sqrt( sum( (P-P_save)**2 ) )
      write(id,'(x,i5)',advance='no') iter
      write(id,'(x,f16.'//acc//')',advance='no') e
      write(id,'(x,f16.'//acc//')',advance='no') e-eold
      write(id,'(x,f14.'//acc//')',advance='no') rmsd
      if (ldiis) write(id,'(x,f14.'//acc//')',advance='no') emax
      write(id,'(a)')
      if((abs(e-eold).lt.ethr).and.rmsd.lt.pthr) exit scf
      eold = e
   enddo scf
   first=.false.

   write(id,'(72(''-''))')
   write(id,'('' FINAL SCF ENERGY'',f53.'//acc//')') e
   write(id,'(72(''-''))')

end subroutine rhf_conventional

subroutine uhf_conventional &
           &  (nat,nbf,nalp,nbet,at,xyz,ethr,pthr,first, &
           &   acc,maxiter,ldiis,maxdiis,startdiis, &
           &   S,V,T,X,Pa,Pb,H,Fa,Fb,Ca,Cb,eri,epsa,epsb,e)
   use iso_fortran_env, only : id => output_unit
   use precision, only : wp => dp
   use diis,      only : build_diis,diis_fock
   implicit none
   integer, intent(in)    :: nat
   integer, intent(in)    :: nbf
   integer, intent(in)    :: nalp
   integer, intent(in)    :: nbet
   integer, intent(in)    :: at(nat)
   real(wp),intent(in)    :: xyz(3,nat)
   real(wp),intent(in)    :: ethr
   real(wp),intent(in)    :: pthr
   logical, intent(inout) :: first
   logical, intent(inout) :: ldiis
   integer, intent(in)    :: maxiter
   integer, intent(in)    :: maxdiis
   integer, intent(in)    :: startdiis
   real(wp),intent(in)    :: X(nbf,nbf)
   real(wp),intent(in)    :: H(nbf,nbf)
   real(wp),intent(out)   :: Pa(nbf,nbf)
   real(wp),intent(out)   :: Pb(nbf,nbf)
   real(wp),intent(inout) :: Fa(nbf,nbf)
   real(wp),intent(inout) :: Fb(nbf,nbf)
   real(wp),intent(inout) :: Ca(nbf,nbf)
   real(wp),intent(inout) :: Cb(nbf,nbf)
   real(wp),intent(out)   :: epsa(nbf)
   real(wp),intent(out)   :: epsb(nbf)
   real(wp),intent(out)   :: e
   real(wp),intent(in)    :: S(nbf,nbf)
   real(wp),intent(in)    :: V(nbf,nbf)
   real(wp),intent(in)    :: T(nbf,nbf)
   real(wp),intent(in)    :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)
   character(len=*),intent(in) :: acc

   integer  :: iter,err,i,j
   real(wp) :: enuc,eold,rmsd
   real(wp),allocatable :: Pa_save(:,:),Pb_save(:,:)
   logical  :: uidx(nbf,nbf),lidx(nbf,nbf)

!* DIIS
   logical  :: diisa,diisb
   real(wp) :: eamax,ebmax
   real(wp),allocatable :: eamat(:,:,:),Fa_save(:,:,:),cdiisa(:)
   real(wp),allocatable :: ebmat(:,:,:),Fb_save(:,:,:),cdiisb(:)

   write(id,'(a)')
   write(id,'(72(''-''))')
   write(id,'('' Unrestricted HF-SCF calculation'')')

   allocate( Pa_save(nbf,nbf),Pb_save(nbf,nbf), source=0.0_wp )

   if (ldiis) then
   allocate( eamat(nbf,nbf,maxdiis),Fa_save(nbf,nbf,maxdiis), &
   &         ebmat(nbf,nbf,maxdiis),Fb_save(nbf,nbf,maxdiis), &
   &         cdiisa(maxdiis+1),cdiisb(maxdiis+1), source=0.0_wp )
   endif

   call nnrep(enuc,nat,at,xyz)

   call build_density(nbf,nalp,Pa,Ca)
   call build_density(nbf,nbet,Pb,Cb)

   eold = enuc+uescf(H,Fa,Fb,Pa,Pb,nbf)
   iter = 0

   write(id,'(1(''-'')''[ITER]'')',   advance='no')
   write(id,'(9(''-'')''[E(SCF)]'')', advance='no')
   write(id,'(8(''-'')''[ΔE(SCF)]'')',advance='no')
   write(id,'(7(''-'')''[RMS(P)]'')', advance='no')
   if (ldiis) then
      write(id,'(5(''-'')''[max[F,P]]'')',advance='no')
      write(id,'(1(''-''))')
   else
      write(id,'(15(''-''))')
   endif
   write(id,'(x,i5,x,f16.'//acc//')') iter,eold

   scf: do
      iter = iter+1
      if (iter.gt.maxiter) call raise('E','SCF did not converge')
      call build_ufock(nbf,H,Fa,Fb,Pa,Pb,eri)
      
      if (ldiis) then !* DIIS
      call build_diis(nbf,iter,maxdiis,S,Fa,Pa,Fa_save,  &
           &          eamat,eamax,cdiisa,diisa)
      call build_diis(nbf,iter,maxdiis,S,Fb,Pb,Fb_save,  &
           &          ebmat,ebmax,cdiisb,diisb)
      ldiis = diisa.or.diisb
      if (ldiis) then
      if (iter.ge.startdiis) then !* start DIIS
         if (iter.eq.startdiis) write(id,'('' * starting DIIS'')')
         call diis_fock(nbf,iter,maxdiis,Fa,Fa_save,cdiisa)
         call diis_fock(nbf,iter,maxdiis,Fb,Fb_save,cdiisb)
      endif !* start DIIS
      else
         write(id,'('' * shutting down DIIS'')')
         deallocate( Fa_save,Fb_save,eamat,ebmat,cdiisa,cdiisb )
      endif
      endif !* DIIS

      call roothaan_hall(nbf,X,Fa,Ca,epsa,err)
      if(err.ne.0) call raise('E','solving Roothaan-Hall eq. failed')
      call roothaan_hall(nbf,X,Fb,Cb,epsb,err)
      if(err.ne.0) call raise('E','solving Roothaan-Hall eq. failed')
      Pa_save = Pa
      Pb_save = Pb
      call build_density(nbf,nalp,Pa,Ca)
      call build_density(nbf,nbet,Pb,Cb)
      e = enuc+uescf(H,Fa,Fb,Pa,Pb,nbf)
      rmsd = sqrt(sum((Pa-Pa_save)**2)) + sqrt(sum((Pb-Pb_save)**2))

      write(id,'(x,i5)',advance='no') iter
      write(id,'(x,f16.'//acc//')',advance='no') e
      write(id,'(x,f16.'//acc//')',advance='no') e-eold
      write(id,'(x,f14.'//acc//')',advance='no') rmsd
      if (ldiis) write(id,'(x,f14.'//acc//')',advance='no') max(eamax,ebmax)
      write(id,'(a)')

      if((abs(e-eold).lt.ethr).and.rmsd.lt.pthr) exit scf
      eold = e
   enddo scf
   first=.false.

   write(id,'(72(''-''))')
   write(id,'('' FINAL SCF ENERGY'',f53.'//acc//')') e
   write(id,'(72(''-''))')

end subroutine uhf_conventional

!* scf energy for given density + fockian
pure function uescf(H,Fa,Fb,Pa,Pb,nbf) result(e)
   use precision, only : wp => dp
   use blas95,    only : gemm
   implicit none
   integer, intent(in) :: nbf
   real(wp),intent(in) :: H(nbf,nbf)
   real(wp),intent(in) :: Fa(nbf,nbf)
   real(wp),intent(in) :: Fb(nbf,nbf)
   real(wp),intent(in) :: Pa(nbf,nbf)
   real(wp),intent(in) :: Pb(nbf,nbf)

   real(wp) :: e,tmp(nbf,nbf)
   integer  :: i,j

   call gemm(Pa,(H+Fa),tmp,alpha=0.5_wp)
   call gemm(Pb,(H+Fb),tmp,alpha=0.5_wp,beta=1.0_wp)
   e = sum((/(tmp(i,i),i=1,nbf)/))

end function uescf

pure subroutine build_ufock(nbf,H,Fa,Fb,Pa,Pb,eri)
   use precision, only : wp => dp
   use misc,      only : idx
   implicit none
   integer, intent(in)  :: nbf
   real(wp),intent(in)  :: H(nbf,nbf)
   real(wp),intent(in)  :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)
   real(wp),intent(in)  :: Pa(nbf,nbf)
   real(wp),intent(in)  :: Pb(nbf,nbf)
   real(wp),intent(out) :: Fa(nbf,nbf)
   real(wp),intent(out) :: Fb(nbf,nbf)

   integer  :: i,j,k,l,ij,kl,ijkl,ik,jl,ikjl
   real(wp) :: Pkl

   Fa = H
   Fb = H
   do i = 1, nbf
      do j = 1, i
         ij = idx(i,j)
         do k = 1, nbf
            ik = idx(i,k)
            do l = 1, nbf
               kl = idx(k,l)
               jl = idx(j,l)
               ijkl = idx(ij,kl)
               ikjl = idx(ik,jl)
               Pkl = Pa(l,k) + Pb(l,k)
               Fa(j,i) = Fa(j,i)+(Pkl*eri(ijkl)-Pa(l,k)*eri(ikjl))
               Fb(j,i) = Fb(j,i)+(Pkl*eri(ijkl)-Pb(l,k)*eri(ikjl))
            enddo
         enddo
         Fa(i,j) = Fa(j,i)
         Fb(i,j) = Fb(j,i)
      enddo
   enddo
end subroutine build_ufock

end module scf
