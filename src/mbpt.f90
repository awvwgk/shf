module mbpt

contains

subroutine inttrafo_N8(nbf,eri,C)
   use precision, only : wp => dp
   use misc, only : idx
   implicit none
   integer, intent(in)    :: nbf
   real(wp),intent(in)    :: C(nbf,nbf)
   real(wp),intent(inout) :: eri(nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2)
   real(wp),allocatable   :: temp(:)
   real(wp) :: dum
   integer :: i,j,k,l,ii,jj,kk,ll,ij,ji,kl,lk

   allocate( temp(nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2), source = eri )

!$omp parallel private(i,j,k,l,ii,jj,kk,ll,dum) shared(eri)
!$omp do
   do i = 1,nbf; do j = 1,nbf; ij=idx(i,j)
   do k = 1,nbf; do l = 1,nbf; kl=idx(k,l)
      dum = 0.0_wp
      do ii = 1,nbf; do jj = 1,nbf; ji=idx(ii,jj)
      do kk = 1,nbf; do ll = 1,nbf; lk=idx(kk,ll)
         dum = dum+C(ii,i)*C(jj,j)*C(kk,k)*C(ll,l)*temp(idx(ji,lk))
      enddo; enddo; enddo; enddo
      eri(idx(ij,kl)) = dum
   enddo; enddo; enddo; enddo
!$omp enddo
!$omp endparallel

end subroutine inttrafo_N8

!* transform ee-repulsion integrals from ao basis to mo basis
!  use Θ(N⁵) scaling integral transformation with lots of
!  matrix multiplication.
subroutine inttrafo(nbf,eri,C)
   use precision, only : wp => dp
   use blas95, only : gemm
   use misc, only : idx
   implicit none
   integer, intent(in)    :: nbf
   real(wp),intent(in)    :: C(nbf,nbf)
   real(wp),intent(inout) :: eri(nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2)

   real(wp) :: temp(nbf*(nbf+1)/2,nbf*(nbf+1)/2),tmp1(nbf,nbf),tmp2(nbf,nbf)
   integer :: i,j,k,l,ij,kl,ijkl

   temp=0.0_wp

!$omp parallel private(i,j,k,l,ij,kl,ijkl,tmp1,tmp2) shared(temp)
!$omp do schedule(dynamic)
   do i = 1, nbf
      do j = 1, i
         ij = idx(i,j)
         do k = 1, nbf
            do l = 1, k
               kl = idx(k,l)
               ijkl = idx(ij,kl)
               tmp1(k,l) = eri(ijkl)
               tmp1(l,k) = eri(ijkl)
            enddo
         enddo
         call gemm(C,tmp1,tmp2,transa='t')
         call gemm(tmp2,C,tmp1)
         do k = 1, nbf
            do l = 1, k
               kl = idx(k,l)
               temp(kl,ij) = tmp1(k,l)
            enddo
         enddo
      enddo
   enddo
!$omp enddo
!$omp endparallel
   eri = 0.0_wp
!$omp parallel private(i,j,k,l,ij,kl,ijkl,tmp1,tmp2) shared(eri)
!$omp do schedule(dynamic)
   do k = 1, nbf
      do l = 1, k
         kl = idx(k,l)
         do i = 1, nbf
            do j = 1, i
               ij = idx(i,j)
               tmp1(i,j) = temp(kl,ij)
               tmp1(j,i) = temp(kl,ij)
            enddo
         enddo
         call gemm(C,tmp1,tmp2,transa='t')
         call gemm(tmp2,C,tmp1)
         do i = 1, nbf
            do j = 1, i
               ij = idx(i,j)
               ijkl = idx(kl,ij)
               eri(ijkl) = tmp1(i,j)
            enddo
         enddo
      enddo
   enddo
!$omp enddo
!$omp endparallel

end subroutine inttrafo

!* second order Møller-Plesset many-body pertubation theory
!  Ec = E(SS) + E(OS)
!     =  Σij Σab ( 2(ia|jb)-(ib|ja) ) · (ia|jb) / (εi+εj+εa+εb)
!  E(SS) = ½Σij εij + ½Σîĵ εîĵ
!        = Σij Σab ( (ia|jb)-(ib|ja) ) · (ia|jb) / (εi+εj+εa+εb)
!  E(OS) = Σiĵ εiĵ
!        = Σij Σab (ia|jb)² / (εi+εj+εa+εb)
!  εij = Σab (T(ij,ab) - T(ij,ba))·(ia|jb)
!  εiĵ = Σab T(iĵ,ab)·(ia|ĵb)
!  T(ij,ab) = (ia|jb)/(εi+εj+εa+εb)
subroutine mp2(nbf,nocc,C,eri,eps,e,acc)
   use precision, only : wp => dp
   use misc, only : idx
   implicit none
   integer, intent(in)    :: nbf,nocc
   real(wp),intent(in)    :: eps(nbf),C(nbf,nbf)
   real(wp),intent(inout) :: eri(nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2)
   real(wp),intent(inout) :: e
   character(len=*),intent(in) :: acc
   real(wp),allocatable   :: eridup(:)
   integer  :: i,j,a,b,ia,jb,ib,ja,iajb,ibja
   real(wp) :: dum,eridum,ess,eos,ssc,osc,dmt

   eos=0.0_wp
   ess=0.0_wp

   print'(a)'
   print'(''------------------------------------------'')'
   print'('' second order Møller-Plesset calculation'')'

!  print'('' * doing Θ(N⁸) integral transformation'')'
!  call inttrafo_N8(nbf,eri,C)
   print'('' * doing Θ(N⁵) integral transformation'')'
   call inttrafo(nbf,eri,C)
   print'(''----[i]-----[j]-------------[pair energy]-'')'
   do i = 1, nocc
      do j = 1, nocc
         osc = 0.0_wp
         ssc = 0.0_wp
         do a = nocc+1, nbf
            ia = idx(i,a)
            ja = idx(j,a)
            do b = nocc+1, nbf
               ib = idx(i,b)
               jb = idx(j,b)
               iajb = idx(ia,jb)
               ibja = idx(ib,ja)
               dmt = 1.0_wp / (eps(i)+eps(j)-eps(a)-eps(b))
               osc = osc + eri(iajb)*(eri(iajb)-eri(ibja)) * dmt
               ssc = ssc + eri(iajb)**2 * dmt
            enddo
         enddo
         print'(x,i5,3x,i5,x,f25.'//acc//')',i,j,osc+ssc
         ess=ess + ssc
         eos=eos + osc
      enddo
   enddo
   e = e + ess+eos
   print'(''------------------------------------------'')'
   print'('' same spin energy    :'',f18.'//acc//')',ess
   print'('' opposite spin energy:'',f18.'//acc//')',eos
   print'(''------------------------------------------'')'
   print'('' FINAL MP2 ENERGY'',f23.'//acc//')',e
   print'(''------------------------------------------'')'
end subroutine mp2 

end module mbpt
