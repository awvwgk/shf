!* calculates center of mass and shifts the geometry
!  also aligns the molecule to the eigenvectors of the
!  inertia tensor
subroutine comshift(nat,at,xyz,com)
   use precision, only : wp => dp
   use atomdata, only : ams
   use lapack95, only : syev
   use blas95,   only : gemm
   use printout, only : prmat
   implicit none
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(inout) :: xyz(3,nat)
   real(wp),intent(out) :: com(3)
   integer :: i
   real(wp) :: mass,it(3,3),moi(3)

   mass = 0.0_wp
   it = 0.0_wp
   com = 0.0_wp

   do i = 1, nat
      if(at(i).eq.0) cycle
      mass = mass + ams(at(i))
      com = com + xyz(:,i)*ams(at(i))
   enddo
   com = com/mass
   do i = 1, nat
      xyz(:,i) = xyz(:,i) - com
      it(1,1) = it(1,1) + ams(at(i)) * ( (xyz(2,i))**2 + (xyz(3,i))**2 )
      it(2,2) = it(2,2) + ams(at(i)) * ( (xyz(1,i))**2 + (xyz(3,i))**2 )
      it(3,3) = it(3,3) + ams(at(i)) * ( (xyz(1,i))**2 + (xyz(2,i))**2 )
      it(1,2) = it(1,2) + ams(at(i)) * xyz(1,i)*xyz(2,i)
      it(1,3) = it(1,3) + ams(at(i)) * xyz(1,i)*xyz(3,i)
      it(2,3) = it(2,3) + ams(at(i)) * xyz(2,i)*xyz(3,i)
   enddo
   if((it(1,2).ne.0.0_wp).or.(it(1,3).ne.0.0_wp).or.(it(2,3).ne.0.0_wp)) then
   it(2,1) = it(1,2) 
   it(3,1) = it(1,3) 
   it(3,2) = it(2,3) 
!  call prmat(it,3,3)
!  call syev(it,moi,'v')
!  call prmat(it,3,3)
   endif

end subroutine comshift
