module printout

   interface prmat
!  module procedure prgemat_e
!  module procedure prsymat_e
   module procedure prgemat_f
   module procedure prsymat_f
   end interface prmat

contains

subroutine prenergy(nat,nbf,nalp,nbet,at,xyz,acc, &
           &        S,V,T,eri,H,Fa,Fb,Pa,Pb,Ca,Cb,epsa,epsb)
   use precision, only : wp => dp
   use scf, only : nnrep,uescf
   use blas95, only : gemm
   implicit none
   integer, intent(in) :: nat,nbf
   integer, intent(in) :: nalp,nbet
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   character(len=*),intent(in) :: acc
   real(wp),intent(in) :: S(nbf,nbf)
   real(wp),intent(in) :: V(nbf,nbf)
   real(wp),intent(in) :: T(nbf,nbf)
   real(wp),intent(in) :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)
   real(wp),intent(in) :: H(nbf,nbf)
   real(wp),intent(in) :: Fa(nbf,nbf),Fb(nbf,nbf)
   real(wp),intent(in) :: Pa(nbf,nbf),Pb(nbf,nbf)
   real(wp),intent(in) :: Ca(nbf,nbf),Cb(nbf,nbf)
   real(wp),intent(in) :: epsa(nbf),epsb(nbf)

   integer  :: i,j,k,l
   real(wp) :: enuc,ekin,epot,ej,ex,s2,sc,stmp
   real(wp),allocatable :: tmp1(:,:),tmp2(:,:),tmp3(:,:)

   allocate( tmp1(nbf,nbf),tmp2(nbf,nbf), &
   &         source = 0.0_wp )

   call nnrep(enuc,nat,at,xyz)
   ekin = uescf(T,T,T,Pa,Pb,nbf)
   epot = uescf(V,Fa-T,Fb-T,Pa,Pb,nbf)

   print'('' nuclear-nuclear-rep.:'',f18.'//acc//')',enuc
   print'('' electronic energy   :'',f18.'//acc//')',epot+ekin

   print'('' virial theorem      :'',f18.'//acc//')',-(epot+enuc)/ekin
   print'(''   kinetic energy    :'',f18.'//acc//')',ekin
   print'(''   potential energy  :'',f18.'//acc//')',epot+enuc

   call gemm(Ca,S,tmp1,transa='t')
   call gemm(tmp1,Cb,tmp2)
   sc = nbet - sum( tmp2(:nalp,:nbet)**2 )
   print'('' spin contamination  :'',f18.'//acc//')',sc

!  print'('' Hartree energy      :'',f18.'//acc//')',ej
!  print'('' Fock exchange energy:'',f18.'//acc//')',ex

   print'(72(''-''))'

end subroutine prenergy

subroutine prinput(fname,nat,nel,nocc,nalp,nbet,nbf,at,xyz,zeta,aoc,ng,brsym)
   use precision, only : wp => dp
   use strings, only : capitalize
   !use geom, only : comshift
   implicit none
   character(len=:),allocatable,intent(in) :: fname
   integer, intent(in)    :: nat,nel,nbf
   integer, intent(in)    :: nalp,nbet,nocc
   logical, intent(in)    :: brsym
   integer, intent(in)    :: at(nat)
   real(wp),intent(inout) :: xyz(3,nat)
   real(wp),intent(in)    :: zeta(nbf)
   integer, intent(in)    :: aoc(2,nat)
   integer, intent(in)    :: ng(nbf)

   real(wp) :: com(3)
   integer :: i

   print'(72(''-''))'
   if (.not.allocated(fname)) then
   print'('' Geometry read from STDIN'')'
   else
   print'('' Geometry read from'',x,a)',fname
   endif
   print'(72(''-''))'
   print'('' Number of atoms     :'',i18)',nat
   print'('' Number of electrons :'',i18)',nel
   if(brsym) then
   print'('' Number of α-el.     :'',i18)',nalp
   print'('' Number of β-el.     :'',i18)',nbet
   else
   print'('' Number of el. pairs :'',i18)',nocc
   endif
   print'('' Number of AOs       :'',i18)',nbf
   print'('' Number of prim. AOs :'',i18)',sum(ng)
   print'(''-[#]---[Z]----------[x]------[y]------[z]-''30(''-''))'
   do i = 1, nat
      print'(i3,x,i5,x,a2,x,3(x,f8.3))', &
      &  i,at(i),capitalize(esym(at(i))),xyz(1:3,i)
   enddo
   call comshift(nat,at,xyz,com)
!  print'(''--------------------[x]------[y]------[z]-''30(''-''))'
   print'(4x''c.o.m.'',3x,3(x,f8.3))',com
   print'(72(''-''))'

end subroutine prinput

subroutine prchrg(nat,at,z,acc)
   use precision, only : wp => dp
   use strings, only : capitalize
   implicit none
   integer, intent(in) :: nat,at(nat)
   real(wp),intent(in) :: z(nat)
   character(len=*),intent(in) :: acc

   integer  :: i

   if (minval(at).lt.1) &
   &  call raise('W','Charges make no sense when ghost atoms are present')
   print'(''----[#]-----[Z]-----------------------[q]-''30(''-''))'
   print'(x,i5,3x,i5,x,a2,3x,f20.'//acc//')', &
   &  (i,at(i),capitalize(esym(at(i))),z(i),i=1,nat)
   print'(72(''-''))'
   
end subroutine prchrg

subroutine prrestart(fname,nat,nel,nbf,at,xyz,zeta,aoc,ng,ityp,C,acc)
   use precision, only : wp => dp
   implicit none
!* INPUT
   character(len=*),intent(in) :: fname
   integer, intent(in) :: nat,nel,nbf
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   real(wp),intent(in) :: zeta(nbf)
   integer, intent(in) :: aoc(2,nat)
   integer, intent(in) :: ng(nbf)
   integer, intent(in) :: ityp(nbf)
   real(wp),intent(in) :: C(nbf,nbf)
   character(len=*),intent(in) :: acc
!* TEMPORARY
   integer :: i,ii,j,id
   logical :: exist
   real(wp) :: x,y,z

   id = 42

!* Open the restart file and write the header
   inquire(file=fname,exist=exist)
   if (exist) then
      open(id,file=fname,status='replace',action='write')
   else
      open(id,file=fname,status='new')
   endif
  
   write(id,'(3(x,i5))') nat,nel,nbf

   do i = 1, nat !* loop over all atoms
      write(id,'(3(x,f16.'//acc//'),2(x,i5))') &
      &  xyz(:,i),at(i),aoc(2,i)-aoc(1,i)+1
      do ii = aoc(1,i), aoc(2,i)
      write(id,'(f17.'//acc//',x,a2,i5)') &
      &  zeta(ii),bsym(ityp(ii)),ng(ii)
      enddo
   enddo
   !* write orbitals to file
   write(id,'(''restart'')')
   do i = 1, nbf
      write(id,'(i5)') i
      write(id,'(4(x,e16.'//acc//'))') C(:,i)
   enddo
   close(id)

end subroutine prrestart

elemental function fsym(i) result(sym)
   integer,intent(in) :: i
   character(len=5) :: sym
   character(len=5),parameter :: par(*) = (/ &
   & 's',  'px', 'py', 'pz',  'dx²', 'dy²', 'dz²', 'dxy', 'dxz', 'dyz',   &
   & 'fx³','fy³','fz³','fx²y','fx²z','fy²x','fy²z','fxz²','fyz²','fxyz',  &
   & 'gx⁴','gy⁴','gz⁴','gx³y','gx³z','gy³x','gy³z','gz³x','gz³y','gx²y²', &
   & 'gx²z²',    'gy²z²',     'gx²yz',      'gy²xz',      'gz²xy'  /)
   if (i.eq.0) then
      sym = '{}'
   else
      sym = par(i)
   endif
end function fsym

elemental function bsym(i) result(sym)
   integer,intent(in) :: i
   character(len=2) :: sym
   character(len=2),parameter :: par(*) = (/ &
   & '1s','2s','3s','4s','5s', &
   &      '2p','3p','4p','5p', &
   &           '3d','4d','5d', &
   &                '4f','5f', &
   &                     '5g'  /)
   if (i.eq.0) then
      sym = '{}'
   else
      sym = par(i)
   endif
end function bsym

elemental function esym(i) result(sym)
   integer,intent(in) :: i
   character(len=2) :: sym
   character(len=2),parameter :: pse(*) = (/  &
   & 'h ','he',                               &
   & 'li','be','b ','c ','n ','o ','f ','ne', &
   & 'na','mg','al','si','p ','s ','cl','ar', &
   & 'ca','k ',                               &
   & 'sc','ti','v ','cr','mn','fe','co','ni','cu','zn', &
   &           'ga','ge','as','se','br','kr', &
   & 'rb','sr',                               &
   & 'y ','zr','nb','mo','tc','ru','rh','pd','ag','cd', &
   &           'in','sn','sb','te','i ','xe', &
   & 'cs','ba',                               &
   & 'la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb', &
   & 'lu','hf','ta','w ','re','os','ir','pt','au','hg', &
   &           'tl','pb','bi','po','at','rn', &
   & 'fr','ra',                               &
   & 'ac','th','pa','u ','np','pu','am','cm','bk','cf','es','fm','md','no', &
   & 'lr','rf','db','sg','bh','hs','mt','ds','rg','cn', &
   &           'nh','fl','mc','lv','ts','og'  /)
   if (i.eq.0) then
      sym = '{}'
   else
      sym = pse(i)
   endif
end function esym

subroutine prgemat_e(mat,d1,d2,name,inunit,instep)
   use iso_fortran_env, only : output_unit
   use precision, only : wp => dp
   implicit none
   integer, intent(in) :: d1
   integer, intent(in) :: d2
   real(wp),intent(in) :: mat(d1,d2)
   character(len=*),intent(in),optional :: name
   integer, intent(in),optional :: inunit
   integer, intent(in),optional :: instep
   integer :: i,j,k,l,step,unit
   if (present(inunit)) then
      unit = inunit
   else
      unit = output_unit
   endif
   if (present(instep)) then
      step = instep
   else
      step = 4
   endif
   if(present(name)) write(unit,'(/,''matrix printed:'',x,a)') name
   do i = 1, d2, step
      l = min(i+step-1,d2)
      write(unit,'(/,6x)',advance='no')
      do k = i, l
      write(unit,'(6x,i7,3x)',advance='no') k
      enddo
      write(unit,'(a)')
      do j = 1, d1
         write(unit,'(i6)',advance='no') j
         do k = i, l
         write(unit,'(x,e15.8)',advance='no') mat(j,k)
         enddo
         write(unit,'(a)')
      enddo
   enddo
   return
end subroutine prgemat_e

subroutine prsymat_e(mat,d1,name,inunit,instep)
   use iso_fortran_env, only : output_unit
   use precision, only : wp => dp
   use misc,      only : idx
   implicit none
   integer, intent(in) :: d1
   real(wp),intent(in) :: mat(d1*(d1+1))
   character(len=*),intent(in),optional :: name
   integer, intent(in),optional :: inunit
   integer, intent(in),optional :: instep
   integer :: i,j,k,l,step,unit
   if (present(inunit)) then
      unit = inunit
   else
      unit = output_unit
   endif
   if (present(instep)) then
      step = instep
   else
      step = 4
   endif
   if(present(name)) write(unit,'(/,''matrix printed:'',x,a)') name
   do i = 1, d1, step
      l = min(i+step-1,d1)
      write(unit,'(/,6x)',advance='no')
      do k = i, l
      write(unit,'(6x,i7,3x)',advance='no') k
      enddo
      write(unit,'(a)')
      do j = i, d1
         l = min(i+(step-1),j)
         write(unit,'(i6)',advance='no') j
         do k = i, l
         write(unit,'(x,e15.8)',advance='no') mat(idx(j,k))
         enddo
         write(unit,'(a)')
      enddo
   enddo
   return
end subroutine prsymat_e

subroutine prgemat_f(mat,d1,d2,name,inunit,instep)
   use iso_fortran_env, only : output_unit
   use precision, only : wp => dp
   implicit none
   integer, intent(in) :: d1
   integer, intent(in) :: d2
   real(wp),intent(in) :: mat(d1,d2)
   character(len=*),intent(in),optional :: name
   integer, intent(in),optional :: inunit
   integer, intent(in),optional :: instep
   integer :: i,j,k,l,step,unit
   if (present(inunit)) then
      unit = inunit
   else
      unit = output_unit
   endif
   if (present(instep)) then
      step = instep
   else
      step = 5
   endif
   if(present(name)) write(*,'(/,''matrix printed:'',x,a)') name
   do i = 1, d2, step
      l = min(i+step-1,d2)
      write(unit,'(/,6x)',advance='no')
      do k = i, l
      write(unit,'(3x,i7,3x)',advance='no') k
      enddo
      write(unit,'(a)')
      do j = 1, d1
         write(unit,'(i6)',advance='no') j
         do k = i, l
         write(unit,'(x,f12.7)',advance='no') mat(j,k)
         enddo
         write(unit,'(a)')
      enddo
   enddo
   return
end subroutine prgemat_f

subroutine prsymat_f(mat,d1,name,inunit,instep)
   use iso_fortran_env, only : output_unit
   use precision, only : wp => dp
   use misc, only : idx
   implicit none
   integer, intent(in) :: d1
   real(wp),intent(in) :: mat(d1*(d1+1))
   character(len=*),intent(in),optional :: name
   integer, intent(in),optional :: inunit
   integer, intent(in),optional :: instep
   integer :: i,j,k,l,step,unit
   if (present(inunit)) then
      unit = inunit
   else
      unit = output_unit
   endif
   if (present(instep)) then
      step = instep
   else
      step = 6
   endif
   if(present(name)) write(*,'(/,''matrix printed:'',x,a)') name
   do i = 1, d1, step
      l = min(i+step-1,d1)
      write(unit,'(/,6x)',advance='no')
      do k = i, l
      write(unit,'(3x,i7,3x)',advance='no') k
      enddo
      write(unit,'(a)')
      do j = i, d1
         l = min(i+(step-1),j)
         write(unit,'(i6)',advance='no') j
         do k = i, l
         write(unit,'(x,f12.7)',advance='no') mat(idx(j,k))
         enddo
         write(unit,'(a)')
      enddo
   enddo
   return
end subroutine prsymat_f

end module printout
