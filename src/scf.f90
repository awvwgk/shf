module scf

contains

subroutine hcore_guess(nbf,nocc,H,X,F,P,C)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf,nocc
   real(wp),intent(in)  :: H(nbf,nbf),X(nbf,nbf)
   real(wp),intent(out) :: F(nbf,nbf),P(nbf,nbf),C(nbf,nbf)
   real(wp) :: eps(nbf)
   integer  :: i,j,ii,iat,err
!* hcore_guess
   F = H
   call roothaan_hall(nbf,X,F,C,eps,err)
   if(err.ne.0) call raise('E','solving Roothaan-Hall eq. failed')
   call build_density(nbf,nocc,P,C)
end subroutine hcore_guess

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
   integer, intent(out) :: err
   real(wp),intent(out) :: X(nbf,nbf)
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
   use blas95, only : gemm
   implicit none
   integer, intent(in) :: nbf
   real(wp),intent(in) :: H(nbf,nbf),F(nbf,nbf),P(nbf,nbf)
   real(wp) :: e,tmp(nbf,nbf)
   integer  :: i,j

   call gemm(P,(H+F),tmp)
   e = 0.0_wp
   do i = 1, nbf
      e = e + tmp(i,i)
   enddo

end function escf

pure subroutine build_fock(nbf,H,F,P,eri)
   use precision, only : wp => dp
   use misc, only : idx
   implicit none
   integer, intent(in)  :: nbf
   real(wp),intent(in)  :: H(nbf,nbf),P(nbf,nbf)
   real(wp),intent(in)  :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)
   real(wp),intent(out) :: F(nbf,nbf)

   integer :: i,j,k,l,ij,kl,ijkl,ik,jl,ikjl

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
   use blas95, only : gemm
   implicit none
   integer, intent(in)  :: nbf,nocc
   real(wp),intent(in)  :: C(nbf,nbf)
   real(wp),intent(out) :: P(nbf,nbf)
   integer :: i,j

   call gemm(C(:,:nocc),C(:,:nocc),P,transb='t')

end subroutine build_density

pure subroutine roothaan_hall(nbf,X,F,C,eps,err)
   use precision, only : wp => dp
   use lapack95, only : syev
   use blas95,   only : gemm
   implicit none
   integer, intent(in)  :: nbf
   real(wp),intent(in)  :: F(nbf,nbf),X(nbf,nbf)
   real(wp),intent(out) :: C(nbf,nbf),eps(nbf)
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

subroutine rhf(nat,nbf,nocc,at,xyz,ethr,pthr,first, &
           &   acc,maxiter,ldiis,maxdiis,startdiis, &
           &   S,V,T,X,P,H,F,C,eri,eps,e)
   use precision, only : wp => dp
   use diis, only : build_diis,diis_fock
   implicit none
   integer, intent(in)  :: nat,nbf,nocc
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: ethr,pthr
   logical, intent(inout) :: first,ldiis
   integer, intent(in)  :: maxiter,maxdiis,startdiis
   real(wp),intent(out) :: X(nbf,nbf)
   real(wp),intent(out) :: F(nbf,nbf),H(nbf,nbf),P(nbf,nbf)
   real(wp),intent(inout) :: C(nbf,nbf)
   real(wp),intent(out) :: eps(nbf),e
   real(wp),intent(in)  :: S(nbf,nbf),V(nbf,nbf),T(nbf,nbf)
   real(wp),intent(in)  :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)
   character(len=*),intent(in) :: acc

   integer  :: iter,err,i,j,ij
   real(wp) :: enuc,eold,rmsd
   real(wp) :: ekin,epot,ec,ex
   real(wp),allocatable :: P_save(:,:)

!* DIIS
   real(wp),allocatable :: emat(:,:,:),F_save(:,:,:),cdiis(:),SS(:,:)

   print'(a)'
   print'(''------------------------------------------'')'
   print'('' Restricted Hartree-Fock SCF calculation'')'

   allocate( P_save(nbf,nbf), source=0.0_wp )

   if (ldiis) then
   allocate( emat(nbf,nbf,maxdiis),F_save(nbf,nbf,maxdiis), &
   &         cdiis(maxdiis+1), source=0.0_wp )
   endif

   call nnrep(enuc,nat,at,xyz)

   H=T+V
!  call prmat(H,nbf,nbf,name='H') !* debug

   call build_orthonormalizer(nbf,S,X,err)
   if (err.ne.0) call raise('E','building orthonormalizer failed')
!  call prmat(X,nbf,nbf,name='S^-1/2') !* debug

!* initial guess density (currently hcore guess)
   if (first) then
      print'('' * doing hcore guess'')'
      call hcore_guess(nbf,nocc,H,X,F,P,C)
   else ! or use provided orbitals to restart
      call build_density(nbf,nocc,P,C)
      call build_fock(nbf,H,F,P,eri)
   endif
!* debug
!  call prmat(F,nbf,nbf,name='F')
!  call prmat(P,nbf,nbf,name='P')
!  call prmat(C,nbf,nbf,name='C')
   eold = enuc+escf(H,F,P,nbf)
   iter = 0
   print'(''-[ITER]---------[E(SCF)]--------[ΔE(SCF)]-'')'
   print'(x,i5,x,f16.'//acc//')',iter,eold
   scf: do
      iter = iter+1
      if (iter.gt.maxiter) call raise('E','SCF did not converge')
      call build_fock(nbf,H,F,P,eri)

      if (ldiis) then !* DIIS
      call build_diis(nbf,iter,maxdiis,S,F,P,F_save,emat,cdiis,ldiis)
      if (ldiis) then
      if (iter.gt.startdiis) then !* start DIIS
         call diis_fock(nbf,iter,maxdiis,F,F_save,cdiis)
      elseif (iter.eq.startdiis) then
         print'('' * starting DIIS'')'
         call diis_fock(nbf,iter,maxdiis,F,F_save,cdiis)
      endif !* start DIIS
      else
         print'('' * shutting down DIIS'')'
         deallocate( F_save,emat,cdiis )
      endif
      endif !* DIIS
      
      call roothaan_hall(nbf,X,F,C,eps,err)
      if(err.ne.0) call raise('E','solving Roothaan-Hall eq. failed')
      P_save = P
      call build_density(nbf,nocc,P,C)
      e = enuc+escf(H,F,P,nbf)
      rmsd = sqrt( sum( (P-P_save)**2 ) )
      print'(x,i5,x,f16.'//acc//',x,f16.'//acc//')',iter,e,e-eold
      if((abs(e-eold).lt.ethr).and.rmsd.lt.pthr) exit scf
      eold = e
   enddo scf
   first=.false.

   print'(''------------------------------------------'')'
   print'('' FINAL SCF ENERGY'',f23.'//acc//')',e
   print'(''------------------------------------------'')'
end subroutine rhf

subroutine uhf(nat,nbf,nalp,nbet,at,xyz,ethr,pthr,first, &
           &   acc,maxiter,ldiis,maxdiis,startdiis, &
           &   S,V,T,X,Pa,Pb,H,Fa,Fb,Ca,Cb,eri,epsa,epsb,e)
   use precision, only : wp => dp
   use diis, only : build_diis,diis_fock
   implicit none
   integer, intent(in)  :: nat,nbf,nalp,nbet
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: ethr,pthr
   logical, intent(inout) :: first,ldiis
   integer, intent(in)  :: maxiter,maxdiis,startdiis
   real(wp),intent(out) :: X(nbf,nbf)
   real(wp),intent(out) :: Fa(nbf,nbf),H(nbf,nbf),Pa(nbf,nbf)
   real(wp),intent(out) :: Fb(nbf,nbf),Pb(nbf,nbf)
   real(wp),intent(inout) :: Ca(nbf,nbf),Cb(nbf,nbf)
   real(wp),intent(out) :: epsa(nbf),epsb(nbf),e
   real(wp),intent(in)  :: S(nbf,nbf),V(nbf,nbf),T(nbf,nbf)
   real(wp),intent(in)  :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)
   character(len=*),intent(in) :: acc

   integer  :: iter,err,i,j
   real(wp) :: enuc,eold,rmsd
   real(wp),allocatable :: Pa_save(:,:),Pb_save(:,:)
   logical  :: uidx(nbf,nbf),lidx(nbf,nbf)

!* DIIS
   logical :: diisa,diisb
   real(wp),allocatable :: eamat(:,:,:),Fa_save(:,:,:),cdiisa(:)
   real(wp),allocatable :: ebmat(:,:,:),Fb_save(:,:,:),cdiisb(:)

   print'(a)'
   print'(''------------------------------------------'')'
   print'('' Unrestricted HF-SCF calculation'')'

   allocate( Pa_save(nbf,nbf),Pb_save(nbf,nbf), source=0.0_wp )

   if (ldiis) then
   allocate( eamat(nbf,nbf,maxdiis),Fa_save(nbf,nbf,maxdiis), &
   &         ebmat(nbf,nbf,maxdiis),Fb_save(nbf,nbf,maxdiis), &
   &         cdiisa(maxdiis+1),cdiisb(maxdiis+1), source=0.0_wp )
   endif

   call nnrep(enuc,nat,at,xyz)

   H = T+V

   call build_orthonormalizer(nbf,S,X,err)
   if (err.ne.0) call raise('E','building orthonormalizer failed')

!* initial guess density (currently hcore guess)
   print'('' * doing hcore guess'')'
   call hcore_guess(nbf,nalp,H,X,Fa,Pa,Ca)
   call hcore_guess(nbf,nbet,H,X,Fb,Pb,Cb)
   eold = enuc+uescf(H,Fa,Fb,Pa,Pb,nbf)
   iter = 0
   print'(''-[ITER]---------[E(SCF)]--------[ΔE(SCF)]-'')'
   print'(x,i5,x,f16.'//acc//')',iter,eold
   scf: do
      iter = iter+1
      if (iter.gt.maxiter) call raise('E','SCF did not converge')
      call build_ufock(nbf,H,Fa,Fb,Pa,Pb,eri)
      
      if (ldiis) then !* DIIS
      call build_diis(nbf,iter,maxdiis,S,Fa,Pa,Fa_save,eamat,cdiisa,diisa)
      call build_diis(nbf,iter,maxdiis,S,Fb,Pb,Fb_save,ebmat,cdiisb,diisb)
      ldiis = diisa.or.diisb
      if (ldiis) then
      if (iter.gt.startdiis) then !* start DIIS
         call diis_fock(nbf,iter,maxdiis,Fa,Fa_save,cdiisa)
         call diis_fock(nbf,iter,maxdiis,Fb,Fb_save,cdiisb)
      elseif (iter.eq.startdiis) then
         print'('' * starting DIIS'')'
         call diis_fock(nbf,iter,maxdiis,Fa,Fa_save,cdiisa)
         call diis_fock(nbf,iter,maxdiis,Fb,Fb_save,cdiisb)
      endif !* start DIIS
      else
         print'('' * shutting down DIIS'')'
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
      print'(x,i5,x,f16.'//acc//',x,f16.'//acc//')',iter,e,e-eold
      if((abs(e-eold).lt.ethr).and.rmsd.lt.pthr) exit scf
      eold = e
   enddo scf
   first=.false.
   print'(''------------------------------------------'')'
   print'('' FINAL SCF ENERGY'',f23.'//acc//')',e
   print'(''------------------------------------------'')'
end subroutine uhf

!* scf energy for given density + fockian
pure function uescf(H,Fa,Fb,Pa,Pb,nbf)
   use precision, only : wp => dp
   use blas95, only : gemm
   implicit none
   integer, intent(in) :: nbf
   real(wp),intent(in) :: H(nbf,nbf),Fa(nbf,nbf),Fb(nbf,nbf)
   real(wp),intent(in) :: Pa(nbf,nbf),Pb(nbf,nbf)
   real(wp) :: uescf,tmp(nbf,nbf)
   integer  :: i,j
   call gemm(Pa,(H+Fa),tmp,alpha=0.5_wp)
   call gemm(Pb,(H+Fb),tmp,alpha=0.5_wp,beta=1.0_wp)
   uescf = 0.0_wp
   do i = 1, nbf
      uescf = uescf + tmp(i,i)
   enddo
end function uescf

pure subroutine build_ufock(nbf,H,Fa,Fb,Pa,Pb,eri)
   use precision, only : wp => dp
   use misc, only : idx
   implicit none
   integer, intent(in)  :: nbf
   real(wp),intent(in)  :: H(nbf,nbf)
   real(wp),intent(in)  :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)
   real(wp),intent(in)  :: Pa(nbf,nbf),Pb(nbf,nbf)
   real(wp),intent(out) :: Fa(nbf,nbf),Fb(nbf,nbf)

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
