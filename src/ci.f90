module ci
   use precision, only : wp => dp
   implicit none

   real(wp),parameter :: autoeV = 27.21138602_wp

   interface cis
   module procedure solve_cis
   end interface cis

contains

subroutine solve_cis(nbf,nel,F,tei,acc)
   use iso_fortran_env, only : id => output_unit
   use precision, only : wp => dp
   use lapack95,  only : syev
   use printout,  only : prmat,prvec
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   character(len=*),intent(in) :: acc

   integer  :: i,j,k,l,a,b
   integer  :: nvir,ncis
   integer  :: err
   real(wp),allocatable :: H(:,:),E(:)

   write(id,'(a)')
   write(id,'(72(''-''))')
   write(id,'('' excited states by configuration interaction singles'')')

   nvir = 2*nbf - nel
   ncis = nvir*nel

   allocate( H(ncis,ncis),E(ncis), source = 0.0_wp )

   k = 0
   do i = 1, nel
      do a = 1, nvir
         k = k + 1
         l = 0
         do j = 1, nel
            do b = 1, nvir
               l = l+1
               H(k,l) = tei(a+nel,j,i,b+nel)
               if(i.eq.j) H(k,l) = H(k,l) + F(a+nel,b+nel)
               if(a.eq.b) H(k,l) = H(k,l) - F(i,j)
            enddo
         enddo
      enddo
   enddo

   call prmat(H,ncis,ncis,'CIS Hamiltonian')

   call syev(H,E,'n','u',err)
   if (err.ne.0) call raise('E','could not solve CIS Hamiltonian')

   call prvec(E*autoeV,ncis,'excitation energies in eV')

end subroutine solve_cis

end module ci
