module pop

contains

subroutine mulliken(nat,nbf,at,aoc,S,Pa,Pb,z)
   use precision, only : wp => dp
   use blas95, only : gemm
   implicit none
   integer, intent(in)  :: nat,nbf
   integer, intent(in)  :: at(nat)
   integer, intent(in)  :: aoc(2,nat)
   real(wp),intent(in)  :: S(nbf,nbf)
   real(wp),intent(in)  :: Pa(nbf,nbf),Pb(nbf,nbf)
   real(wp),intent(out) :: z(nat)

   integer  :: i,ii,j,jj,ij
   real(wp) :: tmp(nbf,nbf)

   call gemm(Pa+Pb,S,tmp)

   do i = 1, nat
      z(i) = at(i)
      do ii = aoc(1,i), aoc(2,i)
         z(i) = z(i) - tmp(ii,ii)
      enddo
   enddo

end subroutine mulliken

subroutine loewdin(nat,nbf,at,aoc,S,Pa,Pb,X,z)
   use precision, only : wp => dp
   use blas95,   only : gemm
   implicit none
   integer, intent(in)  :: nat,nbf
   integer, intent(in)  :: at(nat)
   integer, intent(in)  :: aoc(2,nat)
   real(wp),intent(in)  :: S(nbf,nbf)
   real(wp),intent(in)  :: Pa(nbf,nbf),Pb(nbf,nbf)
   real(wp),intent(in)  :: X(nbf,nbf)
   real(wp),intent(out) :: z(nat)

   integer  :: i,ii,j,jj,ij
   real(wp) :: tmp1(nbf,nbf)
   real(wp) :: tmp2(nbf,nbf)
   real(wp) :: tmp3(nbf,nbf)

   call gemm(X,S,tmp1)
   call gemm(Pa+Pb,tmp1,tmp2)
   call gemm(tmp1,tmp2,tmp3)

   do i = 1, nat
      z(i) = at(i)
      do ii = aoc(1,i), aoc(2,i)
         z(i) = z(i) - tmp3(ii,ii)
      enddo
   enddo

end subroutine loewdin

end module pop
