module rpa
   use precision, only : wp => dp
   implicit none

   real(wp),parameter :: autoeV = 27.21138602_wp

   interface tdhf
   module procedure solve_tdhf
   end interface tdhf

contains
subroutine solve_tdhf(nbf,nel,F,tei,acc)
   use iso_fortran_env, only : id => output_unit
   use precision, only : wp => dp
   use blas95,    only : gemm
   use lapack95,  only : geev
   use printout,  only : prmat,prvec
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   character(len=*),intent(in) :: acc

   integer  :: i,j,k,l,a,b
   integer  :: nvir,ndim
   integer  :: err
   ! solve for (A+B)(A-B)(X-Y) = E²(X-Y)
   real(wp),allocatable :: XmY(:,:)
   real(wp),allocatable :: Amat(:,:)
   real(wp),allocatable :: Bmat(:,:)
   real(wp),allocatable :: E2(:)
   real(wp),allocatable :: Ei(:)

   write(id,'(a)')
   write(id,'(72(''-''))')
   write(id,'('' excited states by time dependent HF'')')

   nvir = 2*nbf - nel
   ndim = nvir*nel

   allocate( Amat(ndim,ndim),Bmat(ndim,ndim),XmY(ndim,ndim),  &
   &         E2(ndim),Ei(ndim), source = 0.0_wp )

   k = 0
   do i = 1, nel
      do a = 1, nvir
         k = k + 1
         l = 0
         do j = 1, nel
            do b = 1, nvir
               l = l+1
               Amat(k,l) = tei(a+nel,j,i,b+nel)
               if(i.eq.j) Amat(k,l) = Amat(k,l) + F(a+nel,b+nel)
               if(a.eq.b) Amat(k,l) = Amat(k,l) - F(i,j)
               Bmat(k,l) = tei(a+nel,b+nel,i,j)
            enddo
         enddo
      enddo
   enddo
   call gemm(Amat+Bmat,Amat-Bmat,XmY)

   call prmat(XmY,ndim,ndim,'(A+B)(A-B) matrix')

   call geev(XmY,E2,Ei,info=err)
   if (err.ne.0) call raise('E','could not solve TD-HF')

   call prvec(sqrt(E2)*autoeV,ndim,'excitation energies in eV')

end subroutine solve_tdhf
end module rpa
