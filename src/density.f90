module density
   use precision, only : wp => dp
   use typedef,   only : tmesh

contains

subroutine rdmesh(fname,mesh)
   use precision, only : wp => dp
   implicit none
   character(len=*),intent(in) :: fname
   type(tmesh),intent(out)     :: mesh

   integer  :: i,id
   integer  :: err
   integer  :: nmesh
   real(wp) :: x,y,z,w
   logical  :: exist

   id = 42

   inquire(file=fname,exist=exist)
   if (.not.exist) call raise('E','File not found: '//fname)

   open(id,file=fname,status='old')
   read(id,*,iostat=err) nmesh
   if (err.ne.0) call raise('E','Error while reading mesh file')
   if (nmesh.lt.1) call raise('E','No mesh points')
   mesh % n = nmesh
   allocate( mesh % xyz(3,nmesh), mesh % w(nmesh), mesh % rho(nmesh,2),  &
   &         source = 0.0_wp )
   
   do i = 1, nmesh
      read(id,*,iostat=err) x,y,z,w
      if (err.ne.0) call raise('E','Error while reading mesh file')
      mesh % xyz(:,i) = (/x,y,z/)
      mesh % w(i) = w
   enddo
   
end subroutine rdmesh

subroutine prdens(mesh,acc)
   use iso_fortran_env, only : id => output_unit
   use precision, only : wp => dp
   implicit none
   type(tmesh),intent(in) :: mesh
   character(len=*),intent(in) :: acc
 
   integer  :: i

   write(id,'(''-[#]-------[x]------[y]------[z]'')',advance='no')
   write(id,'(14(''-''),''[ρ(α)]'',13(''-''),''[ρ(β)]-'')')
   do i = 1, mesh % n
      write(id,'(i3,x,3(x,f8.3),x,f19.'//acc//',f19.'//acc//')')  &
      &     i,mesh % xyz(:,i),mesh % w(i) * mesh % rho(i,1),  &
      &                       mesh % w(i) * mesh % rho(i,2)
   enddo
   write(id,'(72(''-''))')

end subroutine prdens

subroutine calc_density(mesh,nat,nbf,at,xyz,zeta,aoc,ng,ityp,Pa,Pb)
   use precision, only : wp => dp
   implicit none
   type(tmesh),intent(inout) :: mesh
   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   integer, intent(in)  :: ng(nbf)
   integer, intent(in)  :: ityp(nbf)
   integer, intent(in)  :: at(nat)
   integer, intent(in)  :: aoc(2,nat)
   real(wp),intent(in)  :: zeta(nbf)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: Pa(nbf,nbf)
   real(wp),intent(in)  :: Pb(nbf,nbf)
 
   integer  :: i
   real(wp) :: rhoa,rhob
   real(wp) :: point(3)

   do i = 1, mesh % n
      point = mesh % xyz(:,i)
      call densdriver(point,nat,nbf,xyz,zeta,aoc,ng,ityp,Pa,Pb,rhoa,rhob)
      mesh % rho(i,1) = rhoa
      mesh % rho(i,2) = rhob
   enddo

end subroutine calc_density

!* driver for calculation of density at a point
pure subroutine densdriver(point,nat,nbf,xyz,zeta,aoc,ng,ityp,  &
                &          Pa,Pb,rhoa,rhob)
   use precision, only : wp => dp
   use ints,      only : maxprim
   use stong,     only : slater
   implicit none
   real(wp),intent(in)  :: point(3)
   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   integer, intent(in)  :: ng(nbf)
   integer, intent(in)  :: ityp(nbf)
   integer, intent(in)  :: aoc(2,nat)
   real(wp),intent(in)  :: zeta(nbf)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: Pa(nbf,nbf)
   real(wp),intent(in)  :: Pb(nbf,nbf)
   real(wp),intent(out) :: rhoa
   real(wp),intent(out) :: rhob

   real(wp) :: ci(maxprim),cj(maxprim),alpi(maxprim),alpj(maxprim)
   real(wp) :: rhodum

   integer  :: i,j,ii,jj
   integer  :: ishtyp,jshtyp

   intrinsic :: dble

   rhoa = 0.0_wp
   rhob = 0.0_wp

   do i = 1, nat
      do j = 1, nat
         do ii = aoc(1,i), aoc(2,i)
            !* on-the-fly expansion
            call slater(ityp(ii),ng(ii),zeta(ii),alpi,ci) 
            do jj = aoc(1,j), aoc(2,j)
               call slater(ityp(jj),ng(jj),zeta(jj),alpj,cj)
               call onedens(ng(ii),ng(jj), & ! ,ishtyp,jshtyp, &
                    &       xyz(:,i),xyz(:,j),alpi,alpj,ci,cj, &
                    &       point,rhodum)
               rhoa = rhoa + Pa(ii,jj) * rhodum
               rhob = rhob + Pb(ii,jj) * rhodum
            enddo
         enddo
      enddo
   enddo

end subroutine densdriver

!* one electron density over s-functions
pure subroutine onedens(npa,npb,r_a,r_b,alp,bet,ci,cj,point, &
                &       rho)
   use precision, only : wp => dp
   use ints,      only : i_thr,r_thr
   implicit none

   integer, intent(in)  :: npa
   integer, intent(in)  :: npb
   real(wp),intent(in)  :: r_a(3)
   real(wp),intent(in)  :: r_b(3)
   real(wp),intent(in)  :: alp(npa)
   real(wp),intent(in)  :: bet(npb)
   real(wp),intent(in)  :: ci(npa)
   real(wp),intent(in)  :: cj(npb)
   real(wp),intent(in)  :: point(3)

   real(wp),intent(out) :: rho

   integer  :: i,j,k
   real(wp) :: rab,ab,eab,oab,xab,sqm,est
   real(wp) :: s00,fact,rcp,r_p(3),cab

   intrinsic :: sum,exp

   rab = sum( (r_a-r_b)**2 )
   if (rab.gt.r_thr) return

   rho = 0.0_wp

   do i=1,npa
      do j=1,npb
         eab = alp(i)+bet(j)
         oab = 1.0_wp/eab
         cab = ci(i)*cj(j)
         xab = alp(i)*bet(j)*oab
         est = rab*xab
         if (est.gt.i_thr) cycle ! estimate contribution, ommit if too small
         ab = exp(-est)
         r_p = (alp(i)*r_a + bet(j)*r_b)*oab
         rho = rho + cab*ab*exp(-eab*sum((point-r_p)**2))
      enddo
   enddo

end subroutine onedens

end module density
