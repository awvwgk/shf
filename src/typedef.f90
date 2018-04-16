module typedef
   use precision, only : wp => dp
   implicit none

   type :: gauss(ng)
      integer,len :: ng
      integer  :: l
      real(wp) :: point(3)
      real(wp) :: alpha(ng)
      real(wp) :: coeff(ng)
   end type gauss

   type :: tmesh
      integer :: n
      real(wp),allocatable :: w(:)     !  allocate( w(n) )
      real(wp),allocatable :: xyz(:,:) !  allocate( xyz(3,n) )
      real(wp),allocatable :: rho(:,:) !  allocate( rho(n,2) )
   end type tmesh

   type :: basis
      integer  :: n
      integer, allocatable :: ityp(:)
      integer, allocatable :: ng(:)
      integer, allocatable :: aoc(:,:)
      integer, allocatable :: sh2at(:)
      real(wp),allocatable :: zeta(:)
      real(wp),allocatable :: alpha(:,:)
      real(wp),allocatable :: coeff(:,:)
   end type basis

end module typedef
