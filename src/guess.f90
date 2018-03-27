module guess

contains

subroutine hcore_guess(nbf,nocc,H,X,F,C)
   use precision, only : wp => dp
   use scf,       only : roothaan_hall
   implicit none
   integer, intent(in)  :: nbf,nocc
   real(wp),intent(in)  :: H(nbf,nbf),X(nbf,nbf)
   real(wp),intent(out) :: F(nbf,nbf),C(nbf,nbf)

   real(wp) :: eps(nbf)
   integer  :: i,j,ii,iat,err

!* hcore_guess
   F = H
   call roothaan_hall(nbf,X,F,C,eps,err)
   if(err.ne.0) call raise('E','solving Roothaan-Hall eq. failed')

end subroutine hcore_guess

end module guess
