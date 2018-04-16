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

end module typedef
