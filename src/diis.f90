module diis

contains

pure subroutine diis_fock(nbf,iter,maxdiis,F,F_save,c)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)    :: nbf
   integer, intent(in)    :: iter,maxdiis
   real(wp),intent(out)   :: F(nbf,nbf)
   real(wp),intent(in)    :: F_save(nbf,nbf,maxdiis)
   real(wp),intent(in)    :: c(maxdiis+1)

   integer :: i,ii,j,jj,err,info

   F = 0.0_wp
   do i = max(1,iter-maxdiis+1), iter; ii = mod(i-1,maxdiis)+1
      F = F + c(ii)*F_save(:,:,ii)
   enddo

end subroutine diis_fock

pure subroutine build_diis(nbf,iter,maxdiis,S,F,P,F_save,emat,c,dodiis)
   use precision, only : wp => dp
   use lapack95, only : gesv
   use blas95,   only : gemm
   implicit none
   integer, intent(in)    :: nbf
   integer, intent(in)    :: iter,maxdiis
   real(wp),intent(in)    :: S(nbf,nbf),F(nbf,nbf),P(nbf,nbf)
   real(wp),intent(inout) :: F_save(nbf,nbf,maxdiis)
   real(wp),intent(inout) :: emat(nbf,nbf,maxdiis)
   real(wp),intent(out)   :: c(maxdiis+1)
   logical, intent(out)   :: dodiis

   real(wp) :: tmp1(nbf,nbf)
   real(wp),allocatable :: B(:,:)
   integer :: i,ii,j,jj,k,err,idiis
   integer :: work(maxdiis+1)

   dodiis=.true.

   idiis = min(iter,maxdiis)+1
   allocate ( B(idiis,idiis), source = 0.0_wp )

!* set up error matrix
   k = mod(iter-1,maxdiis)+1
   F_save(:,:,k) = F
   call gemm(P,F,tmp1)
   call gemm(S,tmp1,emat(:,:,k))
   call gemm(F,P,tmp1)
   call gemm(tmp1,S,emat(:,:,k),beta=-1.0_wp)

   do i = max(1,iter-maxdiis), iter; ii = mod(i-1,maxdiis)+1
      do j = max(1,iter-maxdiis), i; jj = mod(j-1,maxdiis)+1
         B(ii,jj) = sum( emat(:,:,i)*emat(:,:,j) )
         B(jj,ii) = B(ii,jj)
      enddo
   enddo

   B(:,idiis) = -1.0_wp
   B(idiis,:) = -1.0_wp
   B(idiis,idiis) = 0.0_wp
   c = 0.0_wp
   c(idiis) = -1.0_wp

   call gesv(B,c(:idiis),work,err)
   if (err.ne.0) then
      c = 0.0_wp
      c(k) = 1.0_wp
      dodiis = .false.
   endif

end subroutine build_diis

end module diis
