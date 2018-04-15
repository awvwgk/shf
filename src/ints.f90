module ints
   use precision, only : wp => dp
   implicit none

   public  :: integrals
   private :: intdriver_all,intdriver_tei,intdriver_one
   public  :: onetrafo,teitrafo,teitrafo_N8
   public  :: oneint,twoint
   private :: stvint_os2,erint_os2
   private :: factorial,binom,olap,gpt

   interface integrals
   module procedure intdriver_all
   module procedure intdriver_one
   module procedure intdriver_tei
   end interface integrals

   type :: gauss(ng)
      integer,len :: ng
      integer  :: l
      real(wp) :: point(3)
      real(wp) :: alpha(ng)
      real(wp) :: coeff(ng)
   end type gauss

   intrinsic :: shape,reshape

   integer,parameter  :: maxprim = 6

   real(wp),parameter :: r_thr = 2500._wp ! S.lt.10e-9 for r.gt.50 bohr
   real(wp),parameter :: i_thr = 25._wp ! integral cut threshold
   real(wp),parameter :: s_thr = 1e-9_wp ! Schwarz inequaltity threshold
   real(wp),parameter :: c_thr = 1.0e-8_wp ! ommit recursion coeffient
   real(wp),parameter :: t_thr = 0.18_wp ! for boys function
   real(wp),parameter :: l_thr = 19.35_wp ! for boys function
   integer, parameter :: m_max = 28 ! max iter in boys function

   real(wp),parameter :: tth = 2.0_wp/3.0_wp
   real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
   real(wp),parameter :: tpi = 2.0_wp*pi
   real(wp),parameter :: twopi25 = 34.98683665524963_wp
   real(wp),parameter :: tosqpi = 2.0_wp/sqrt(pi)
   real(wp),parameter :: pi2 = 6.36619772367582e-1_wp
   real(wp),parameter :: pi3 = pi**3
   real(wp),parameter :: oth = 1.0_wp/3.0_wp
   real(wp),parameter :: dfactorial(*) = (/  & ! see OEIS A001147
   &  1._wp,1._wp,3._wp,15._wp,105._wp,945._wp,10395._wp,135135._wp,       &
   !* I may need more elements for the evalulation of boys function...
   &  2027025._wp,34459425._wp,654729075._wp,13749310575._wp,              &
   &  316234143225._wp,7905853580625._wp,213458046676875._wp,              &
   &  6190283353629375._wp,191898783962510625._wp,6332659870762850625._wp, &
   &  221643095476699771875._wp,8200794532637891559375._wp /)
   real(wp),parameter :: ofactorial(*) = (/  & ! one over factorial
   &  1._wp, 1._wp, 1._wp/2._wp, 1._wp/6._wp, 1._wp/24._wp, 1._wp/120._wp,  &
   &  1._wp/720._wp, 1._wp/5040._wp, 1._wp/40320._wp, 1._wp/362880._wp,  &
   &  1._wp/3628800._wp, 1._wp/39916800._wp, 1._wp/479001600._wp,  &
   &  1._wp/6227020800._wp, 1._wp/87178291200._wp, 1._wp/1307674368000._wp,  &
   &  1._wp/20922789888000._wp, 1._wp/355687428096000._wp,  &
   &  1._wp/6402373705728000._wp, 1._wp/121645100408832000._wp,  &
   &  1._wp/2432902008176640000._wp, 1._wp/51090942171709440000._wp,  &
   &  1._wp/1124000727777607680000._wp /)

!* Boys function is precalculated on a grid as described in Helgaker2000
!  (see tools/boys.rb to generate anew).
   include 'boysf_grid.f90'

!* possible moments up to g functions (more then I can ever effort):
!> s    px   py   pz    dx²   dy²   dz²   dxy   dxz   dyz
!  1    2    3    4     5     6     7     8     9     10
!> fx³  fy³  fz³  fx²y  fx²z  fy²x  fy²z  fxz²  fyz²  fxyz
!  11   12   13   14    15    16    17    18    19    20
!> gx⁴  gy⁴  gz⁴  gx³y  gx³z  gy³x  gy³z  gz³x  gz³y  gx²y²
!  21   22   23   24    25    26    27    28    29    30
!> gx²z²     gy²z²      gx²yz       gy²xz       gz²xy
!  31        32         33          34          35
   integer,parameter :: lmn(3,35) = reshape((/ &
   !* moments in x-direction
   &  0,  &! s-shell
   &  1,0,0,  &! p-shell
   &  2,0,0,1,1,0,  &! d-shell
   &  3,0,0,2,2,1,0,1,0,1,  &! f-shell
   &  4,0,0,3,3,1,0,1,0,2,2,0,2,1,1,  &! g-shell
   !* moments in y-direction
   &  0,  &! s-shell
   &  0,1,0,  &! p-shell
   &  0,2,0,1,0,1,  &! d-shell
   &  0,3,0,1,0,2,2,0,1,1,  &! f-shell
   &  0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,  &! g-shell
   !* moments in z-direction
   &  0,  &! s-shell
   &  0,0,1,  &! p-shell
   &  0,0,2,0,1,1,  &! d-shell
   &  0,0,3,0,1,0,1,2,2,1,  &! f-shell
   &  0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/),&! g-shell
   !* now sort them in the right order
   &  shape(lmn), order=(/2,1/))

!* Taken from postg by Alberto Otero de la Roza
   real(wp),parameter :: s3 = sqrt(3._wp)
   real(wp),parameter :: s3_4 = sqrt(3._wp/4._wp)
   real(wp),parameter :: s3_8 = sqrt(3._wp/8._wp)
   real(wp),parameter :: s5_8 = sqrt(5._wp/8._wp)
   real(wp),parameter :: s5_16 = sqrt(5._wp/16._wp)
   real(wp),parameter :: s6 = sqrt(6._wp)
   real(wp),parameter :: s10 = sqrt(10._wp)
   real(wp),parameter :: s10_8 = sqrt(10._wp/8._wp)
   real(wp),parameter :: s15 = sqrt(15._wp)
   real(wp),parameter :: s15_4 = sqrt(15._wp/4._wp)
   real(wp),parameter :: s35_4 = sqrt(35._wp/4._wp)
   real(wp),parameter :: s35_8 = sqrt(35._wp/8._wp)
   real(wp),parameter :: s35_64 = sqrt(35._wp/64._wp)
   real(wp),parameter :: s45 = sqrt(45._wp)
   real(wp),parameter :: s45_4 = sqrt(45._wp/4._wp)
   real(wp),parameter :: s45_8 = sqrt(45._wp/8._wp)
   real(wp),parameter :: s315_8 = sqrt(315._wp/8._wp)
   real(wp),parameter :: s315_16 = sqrt(315._wp/16._wp)
   real(wp),parameter :: d32 = 3._wp/2._wp
   real(wp),parameter :: d34 = 3._wp/4._wp
   real(wp),parameter :: d38 = 3._wp/8._wp

!* dtrafo from cartesian to solid harmonics (l = 2)
!  spherical molden order: m = 0, 1, -1, 2, -2
!  Cartesian molden order: xx, yy, zz, xy, xz, yz
   real(wp),parameter :: dtrafo(5,6) = reshape((/  &
!    0      1     -1      2     -2
   -.5_wp, 0._wp, 0._wp,  s3_4, 0._wp,& ! xx
   -.5_wp, 0._wp, 0._wp, -s3_4, 0._wp,& ! yy
    1._wp, 0._wp, 0._wp, 0._wp, 0._wp,& ! zz
    0._wp, 0._wp, 0._wp, 0._wp,    s3,& ! xy
    0._wp,    s3, 0._wp, 0._wp, 0._wp,& ! xz
    0._wp, 0._wp,    s3, 0._wp, 0._wp & ! yz
   /),shape(dtrafo))

!* ftrafo from cartesian to solid harmonics (l = 3)
!  spherical molden order: m = 0, 1, -1, 2, -2, 3, -3
!  Cartesian molden order: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
   real(wp),parameter :: ftrafo(7,10) = reshape((/  &
!    0        1       -1         2       -2         3       -3 
   0._wp,   -s3_8,   0._wp,    0._wp,   0._wp,     s5_8,   0._wp,& ! xxx 
   0._wp,   0._wp,   -s3_8,    0._wp,   0._wp,    0._wp,   -s5_8,& ! yyy 
   1._wp,   0._wp,   0._wp,    0._wp,   0._wp,    0._wp,   0._wp,& ! zzz 
   0._wp,   -s3_8,   0._wp,    0._wp,   0._wp,   -s45_8,   0._wp,& ! xyy 
   0._wp,   0._wp,   -s3_8,    0._wp,   0._wp,    0._wp,   s45_8,& ! xxy 
    -d32,   0._wp,   0._wp,    s15_4,   0._wp,    0._wp,   0._wp,& ! xxz 
   0._wp,      s6,   0._wp,    0._wp,   0._wp,    0._wp,   0._wp,& ! xzz 
   0._wp,   0._wp,      s6,    0._wp,   0._wp,    0._wp,   0._wp,& ! yzz 
    -d32,   0._wp,   0._wp,   -s15_4,   0._wp,    0._wp,   0._wp,& ! yyz 
   0._wp,   0._wp,   0._wp,    0._wp,     s15,    0._wp,   0._wp & ! xyz 
   /),shape(ftrafo))

!* gtrafo from cartesian to solid harmonics (l = 4)
!  spherical molden order: m = 0, 1, -1, 2, -2, 3, -3, 4, -4
!  Cartesian molden order: xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz 
!                          yzzz xxyy xxzz yyzz xxyz xyyz xyzz
   real(wp),parameter :: gtrafo(9,15) = reshape((/  &
!    0      1     -1      2     -2       3      -3        4     -4
     d38, 0._wp, 0._wp,-s5_16, 0._wp,  0._wp,  0._wp,  s35_64, 0._wp,& ! xxxx
     d38, 0._wp, 0._wp, s5_16, 0._wp,  0._wp,  0._wp,  s35_64, 0._wp,& ! yyyy
   1._wp, 0._wp, 0._wp, 0._wp, 0._wp,  0._wp,  0._wp,   0._wp, 0._wp,& ! zzzz
   0._wp, 0._wp, 0._wp, 0._wp,-s10_8,  0._wp,  0._wp,   0._wp, s35_4,& ! xxxy
   0._wp,-s45_8, 0._wp, 0._wp, 0._wp,  s35_8,  0._wp,   0._wp, 0._wp,& ! xxxz
   0._wp, 0._wp, 0._wp, 0._wp,-s10_8,  0._wp,  0._wp,   0._wp,-s35_4,& ! xyyy
   0._wp, 0._wp,-s45_8, 0._wp, 0._wp,  0._wp, -s35_8,   0._wp, 0._wp,& ! yyyz
   0._wp,   s10, 0._wp, 0._wp, 0._wp,  0._wp,  0._wp,   0._wp, 0._wp,& ! xzzz
   0._wp, 0._wp,   s10, 0._wp, 0._wp,  0._wp,  0._wp,   0._wp, 0._wp,& ! yzzz
     d34, 0._wp, 0._wp, 0._wp, 0._wp,  0._wp,  0._wp,-s315_16, 0._wp,& ! xxyy
  -3._wp, 0._wp, 0._wp, s45_4, 0._wp,  0._wp,  0._wp,   0._wp, 0._wp,& ! xxzz
  -3._wp, 0._wp, 0._wp,-s45_4, 0._wp,  0._wp,  0._wp,   0._wp, 0._wp,& ! yyzz
   0._wp, 0._wp,-s45_8, 0._wp, 0._wp,  0._wp, s315_8,   0._wp, 0._wp,& ! xxyz
   0._wp,-s45_8, 0._wp, 0._wp, 0._wp,-s315_8,  0._wp,   0._wp, 0._wp,& ! xyyz
   0._wp, 0._wp, 0._wp, 0._wp,   s45,  0._wp,  0._wp,   0._wp, 0._wp & ! xyzz
   /),shape(gtrafo))


contains

!* expand STO-NGs to gaussians
subroutine expand_stong(nat,nbf,zeta,aoc,ng,ityp,sh2at,alpha,coeff)
   use precision, only : wp => dp
   use stong,     only : slater
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   real(wp),intent(in)  :: zeta(nbf)
   integer, intent(in)  :: aoc(2,nat)
   integer, intent(in)  :: ng(nbf)
   integer, intent(in)  :: ityp(nbf)
   integer, intent(out) :: sh2at(nbf)
   real(wp),intent(out) :: alpha(maxprim,nbf)
   real(wp),intent(out) :: coeff(maxprim,nbf)

   real(wp) :: ci(maxprim),cj(maxprim),alpi(maxprim),alpj(maxprim)
   real(wp) :: ck(maxprim),cl(maxprim),alpk(maxprim),alpl(maxprim)
   real(wp) :: sdum,tdum,vdum,eridum
   real(wp) :: chrg(nat)

   integer  :: i,j,ii,jj
   integer  :: ishtyp,jshtyp,kshtyp,lshtyp

   do i = 1, nat
      do ii = aoc(1,i), aoc(2,i)
         call slater(ityp(ii),ng(ii),zeta(ii),alpha(:,ii),coeff(:,ii)) 
         sh2at(ii) = i
      enddo
   enddo

end subroutine expand_stong

!* driver for driver for calculation of integrals
subroutine intdriver_all(nat,nbf,at,xyz,zeta,aoc,ng,ityp,S,V,T,eri)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   integer, intent(in)  :: ng(nbf)
   integer, intent(in)  :: ityp(nbf)
   integer, intent(in)  :: at(nat)
   integer, intent(in)  :: aoc(2,nat)
   real(wp),intent(in)  :: zeta(nbf)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: S(nbf,nbf)
   real(wp),intent(out) :: T(nbf,nbf)
   real(wp),intent(out) :: V(nbf,nbf)
!  real(wp),intent(out) :: eri(nbf,nbf,nbf,nbf)
!* to save some memory this packing is possible
!  real(wp),intent(out) :: S(nbf*(nbf+1)/2),V(nbf*(nbf+1)/2),T(nbf*(nbf+1)/2)
   real(wp),intent(out) :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)
   real(wp),allocatable :: qcs(:,:)

   allocate( qcs(nbf,nbf), source = 0.0_wp )

   call intdriver_one(nat,nbf,at,xyz,zeta,aoc,ng,ityp,S,V,T)

   call intdriver_qcs(nat,nbf,at,xyz,zeta,aoc,ng,ityp,qcs)
   call intdriver_tei(nat,nbf,at,xyz,zeta,aoc,ng,ityp,qcs,eri)

   deallocate(qcs)

end subroutine intdriver_all

!* driver for calculation of integrals
subroutine intdriver_qcs(nat,nbf,at,xyz,zeta,aoc,ng,ityp,qcs)
   use precision, only : wp => dp
   use stong,     only : slater
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   integer, intent(in)  :: ng(nbf)
   integer, intent(in)  :: ityp(nbf)
   integer, intent(in)  :: at(nat)
   integer, intent(in)  :: aoc(2,nat)
   real(wp),intent(in)  :: zeta(nbf)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: qcs(nbf,nbf)

   real(wp) :: ci(maxprim),cj(maxprim),alpi(maxprim),alpj(maxprim)
   real(wp) :: ck(maxprim),cl(maxprim),alpk(maxprim),alpl(maxprim)
   real(wp) :: qijij
   real(wp) :: chrg(nat)

   integer  :: i,j,ii,jj
   integer  :: ishtyp,jshtyp,kshtyp,lshtyp

!$omp parallel private(i,j,ii,jj,alpi,alpj,ci,cj,qijij) &
!$omp          &       shared(qcs)
!$omp do schedule(dynamic)
   do i = 1, nat
      do j = 1, i
         do ii = aoc(1,i), aoc(2,i)
            !* on-the-fly expansion
            call slater(ityp(ii),ng(ii),zeta(ii),alpi,ci) 
            do jj = aoc(1,j), aoc(2,j)
               call slater(ityp(jj),ng(jj),zeta(jj),alpj,cj)
               call qcsint(ng(ii),ng(jj), & ! ,ishtyp,jshtyp, &
                    &      xyz(:,i),xyz(:,j),alpi,alpj,ci,cj, &
                    &      qijij)
               qcs(ii,jj) = qijij
               qcs(jj,ii) = qijij
            enddo
         enddo
      enddo
   enddo
!$omp enddo
!$omp endparallel

end subroutine intdriver_qcs

!* driver for calculation of integrals
subroutine intdriver_one(nat,nbf,at,xyz,zeta,aoc,ng,ityp,S,V,T)
   use precision, only : wp => dp
   use stong,     only : slater
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   integer, intent(in)  :: ng(nbf)
   integer, intent(in)  :: ityp(nbf)
   integer, intent(in)  :: at(nat)
   integer, intent(in)  :: aoc(2,nat)
   real(wp),intent(in)  :: zeta(nbf)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: S(nbf,nbf)
   real(wp),intent(out) :: T(nbf,nbf)
   real(wp),intent(out) :: V(nbf,nbf)
!* to save some memory this packing is possible
!  real(wp),intent(out) :: S(nbf*(nbf+1)/2),V(nbf*(nbf+1)/2),T(nbf*(nbf+1)/2)

   real(wp) :: ci(maxprim),cj(maxprim),alpi(maxprim),alpj(maxprim)
   real(wp) :: ck(maxprim),cl(maxprim),alpk(maxprim),alpl(maxprim)
   real(wp) :: sdum,tdum,vdum,eridum
   real(wp) :: chrg(nat)

   integer  :: i,j,ii,jj
   integer  :: ishtyp,jshtyp,kshtyp,lshtyp

   chrg = real(at,wp)

!$omp parallel private(i,j,ii,jj,alpi,alpj,ci,cj,sdum,vdum,tdum) &
!$omp          &       shared(s,v,t)
!$omp do schedule(dynamic)
   do i = 1, nat
      do j = 1, i
         do ii = aoc(1,i), aoc(2,i)
            !* on-the-fly expansion
            call slater(ityp(ii),ng(ii),zeta(ii),alpi,ci) 
            do jj = aoc(1,j), aoc(2,j)
               call slater(ityp(jj),ng(jj),zeta(jj),alpj,cj)
               call oneint(ng(ii),ng(jj),nat,xyz,chrg, & ! ,ishtyp,jshtyp, &
                    &      xyz(:,i),xyz(:,j),alpi,alpj,ci,cj, &
                    &      sdum,tdum,vdum)
               s(ii,jj) = sdum
               s(jj,ii) = sdum
               t(ii,jj) = tdum
               t(jj,ii) = tdum
               v(ii,jj) = vdum
               v(jj,ii) = vdum
            enddo
         enddo
      enddo
   enddo
!$omp enddo
!$omp endparallel

end subroutine intdriver_one

!* driver for calculation of integrals
subroutine intdriver_tei(nat,nbf,at,xyz,zeta,aoc,ng,ityp,qcs,eri)
   use precision, only : wp => dp
   use misc,      only : idx
   use stong,     only : slater
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   integer, intent(in)  :: ng(nbf)
   integer, intent(in)  :: ityp(nbf)
   integer, intent(in)  :: at(nat)
   integer, intent(in)  :: aoc(2,nat)
   real(wp),intent(in)  :: zeta(nbf)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: qcs(nbf,nbf)
!  real(wp),intent(out) :: eri(nbf,nbf,nbf,nbf)
!* to save some memory this packing is possible
   real(wp),intent(out) :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)

   real(wp) :: ci(maxprim),cj(maxprim),alpi(maxprim),alpj(maxprim)
   real(wp) :: ck(maxprim),cl(maxprim),alpk(maxprim),alpl(maxprim)
   real(wp) :: sdum,tdum,vdum,eridum

   integer  :: i,j,k,l,ii,jj,kk,ll,limit
   integer  :: ishtyp,jshtyp,kshtyp,lshtyp
   integer  :: ij,kl,ijkl

!$omp parallel private(i,j,k,l,ii,jj,kk,ll,alpi,alpj,alpk,alpl, &
!$omp          &       limit,ij,kl,ijkl,ci,cj,ck,cl,eridum) &
!$omp          &       shared(eri)
!$omp do schedule(dynamic)
   do i = 1, nat
      do j = 1, i
         do k = 1, i
            if (i.eq.k) then; limit=j; else; limit=k; endif
            do l = 1, limit
               do ii = aoc(1,i), aoc(2,i)
                  !* on-the-fly expansion
                  call slater(ityp(ii),ng(ii),zeta(ii),alpi,ci)
                  do jj = aoc(1,j), aoc(2,j)
                     ij = idx(ii,jj)
                     call slater(ityp(jj),ng(jj),zeta(jj),alpj,cj)
                     do kk = aoc(1,k), aoc(2,k)
                        call slater(ityp(kk),ng(kk),zeta(kk),alpk,ck)
                        do ll = aoc(1,l), aoc(2,l)
                           kl = idx(kk,ll)
                           ijkl = idx(ij,kl)
                           call slater(ityp(ll),ng(ll),zeta(ll),alpl,cl)
                           if ((qcs(ii,jj)*qcs(kk,ll)).lt.s_thr) cycle
                           call twoint(ng(ii),ng(jj),ng(kk),ng(ll), &
                           !    &      ishtyp,jshtyp,kshtyp,lshtyp, &
                                &      xyz(:,i),xyz(:,j),xyz(:,k),xyz(:,l), &
                                &      alpi,alpj,alpk,alpl,ci,cj,ck,cl, &
                                &      eridum)
                           eri(ijkl) = eridum
                        !  eri(ii,jj,kk,ll) = eridum
                        !  eri(ii,jj,ll,kk) = eridum
                        !  eri(jj,ii,kk,ll) = eridum
                        !  eri(jj,ii,ll,kk) = eridum
                        !  eri(kk,ll,ii,jj) = eridum
                        !  eri(kk,ll,jj,ii) = eridum
                        !  eri(ll,kk,ii,jj) = eridum
                        !  eri(ll,kk,jj,ii) = eridum
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
!$omp enddo
!$omp endparallel

end subroutine intdriver_tei

!* one electron integrals over s-functions
pure subroutine oneint(npa,npb,nat,xyz,chrg,r_a,r_b,alp,bet,ci,cj, &
                &      sab,tab,vab)
   use precision, only : wp => dp
   implicit none

   integer, intent(in)  :: npa
   integer, intent(in)  :: npb
   integer, intent(in)  :: nat
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: chrg(nat)
   real(wp),intent(in)  :: r_a(3)
   real(wp),intent(in)  :: r_b(3)
   real(wp),intent(in)  :: alp(npa)
   real(wp),intent(in)  :: bet(npb)
   real(wp),intent(in)  :: ci(npa)
   real(wp),intent(in)  :: cj(npb)

   real(wp),intent(out) :: sab
   real(wp),intent(out) :: tab
   real(wp),intent(out) :: vab

   integer  :: i,j,k
   real(wp) :: rab,ab,eab,oab,xab,sqm,est
   real(wp) :: s00,fact,rcp,r_p(3),cab

   intrinsic :: sum,sqrt,exp

   sab = 0.0_wp
   tab = 0.0_wp
   vab = 0.0_wp

   rab = sum( (r_a-r_b)**2 )
   if (rab.gt.r_thr) return

   do i=1,npa
      do j=1,npb
         eab = alp(i)+bet(j)
         oab = 1.0_wp/eab
         cab = ci(i)*cj(j)
         xab = alp(i)*bet(j)*oab
         est = rab*xab
         if (est.gt.i_thr) cycle ! estimate contribution, ommit if too small
         ab = exp(-est)
         s00 = cab*ab*sqrt(pi*oab)**3

!        overlap
         sab = sab+s00

!        kinetic energy
         tab = tab + xab*(3.0_wp-2.0_wp*est)*s00

!        nuclear attraction
         fact = cab*tpi*oab*ab
         r_p = (alp(i)*r_a+bet(j)*r_b)*oab
         do k = 1, nat
            rcp = sum( (r_p-xyz(:,k))**2 )
!           vab  = vab - fact*chrg(k)*boysf0(eab*rcp)
            vab  = vab - fact*chrg(k)*boysf(0,eab*rcp)
         enddo

      enddo
   enddo

end subroutine oneint

!* Cauchy-Schwarz prescreening integrals
pure subroutine qcsint(npa,npb,r_a,r_b,alp,bet,ci,cj,qabab)
   use precision, only : wp => dp
   implicit none

   integer, intent(in)  :: npa
   integer, intent(in)  :: npb
   real(wp),intent(in)  :: r_a(3)
   real(wp),intent(in)  :: r_b(3)
   real(wp),intent(in)  :: alp(npa)
   real(wp),intent(in)  :: bet(npb)
   real(wp),intent(in)  :: ci(npa)
   real(wp),intent(in)  :: cj(npb)
   real(wp),intent(out) :: qabab

   integer  :: i,j
   real(wp) :: rab,est
   real(wp) :: eab,oab,cab
   real(wp) :: ab

   intrinsic :: sum,sqrt,exp

   qabab = 0.0_wp

!  R²(a-b)
   rab=sum( (r_a-r_b)**2 )

   do i = 1, npa
      do j = 1, npb
         cab = ci(i)*cj(j)
         eab = alp(i)+bet(j)
         oab = 1.0_wp/eab
         est = alp(i)*bet(j)*rab*oab
         if (est.gt.i_thr) then
            ab = 0.0_wp
         else
            ab = exp(-est)

            qabab = qabab + cab*ab * sqrt(twopi25/(sqrt(2*eab**5)))
         endif
      enddo
   enddo

end subroutine qcsint

!* two-electron repulsion integral [ab|cd] over s-functions
!  quantity is given in chemist's notation
pure subroutine twoint(npa,npb,npc,npd,r_a,r_b,r_c,r_d, &
                &      alp,bet,gam,del,ci,cj,ck,cl,tei)
   use precision, only : wp => dp
   implicit none

   integer, intent(in)  :: npa
   integer, intent(in)  :: npb
   integer, intent(in)  :: npc
   integer, intent(in)  :: npd
   real(wp),intent(in)  :: r_a(3)
   real(wp),intent(in)  :: r_b(3)
   real(wp),intent(in)  :: r_c(3)
   real(wp),intent(in)  :: r_d(3)
   real(wp),intent(in)  :: alp(npa)
   real(wp),intent(in)  :: bet(npb)
   real(wp),intent(in)  :: gam(npc)
   real(wp),intent(in)  :: del(npd)
   real(wp),intent(in)  :: ci(npa)
   real(wp),intent(in)  :: cj(npb)
   real(wp),intent(in)  :: ck(npc)
   real(wp),intent(in)  :: cl(npd)
   real(wp),intent(out) :: tei

   integer  :: i,j,k,l
   real(wp) :: rab,rcd,rpq,r_p(3),r_q(3),est
   real(wp) :: eab,ecd,eabcd,epq,oab,ocd,cab,ccd
   real(wp) :: ab,cd,abcd,pq,sqpq

   intrinsic :: sum,sqrt,exp

   tei = 0.0_wp

!  R²(a-b)
   rab=sum( (r_a-r_b)**2 )
!  R²(c-d)
   rcd=sum( (r_c-r_d)**2 )

   do i = 1, npa
      do j = 1, npb
         cab = ci(i)*cj(j)
         eab = alp(i)+bet(j)
         oab = 1.0_wp/eab
         est = alp(i)*bet(j)*rab*oab
         if (est.gt.i_thr) cycle
         ab = exp(-est)

!        new gaussian at r_p
         r_p = (alp(i)*r_a+bet(j)*r_b)*oab

         do k = 1, npc
            do l = 1, npd
               ccd = ck(k)*cl(l)
               ecd = gam(k)+del(l)
               ocd = 1.0_wp/ecd
               est = gam(k)*del(l)*rcd*ocd
               if (est.gt.i_thr) cycle
               cd = exp(-est)

!              new gaussian at r_q
               r_q = (gam(k)*r_c+del(l)*r_d)*ocd

!              we already have calculated the prefactors
               abcd = ab*cd

!              distance between product gaussians
               rpq = sum( (r_p-r_q)**2 )

               epq = eab*ecd
               eabcd = eab+ecd

               pq = rpq*epq/eabcd
               tei = tei + cab*ccd*abcd * twopi25/(epq*sqrt(eabcd)) &
!              &           * boysf0(pq)
               &           * boysf(0,pq)

            enddo
         enddo
      enddo
   enddo

end subroutine twoint

! ----------------------------------------------------------------------
!  integral calculation with Obara-Saika recursion relations
! ----------------------------------------------------------------------

subroutine erint_os2(tei)
   use precision, only : wp => dp
   implicit none
!* integrals
   real(wp),intent(inout) :: tei

   integer,parameter  :: max_l = 3 ! up to p-functions

!* temporary variables
   integer  :: i,j,ij,k,l,kl
   real(wp) :: tci(0:max_l,0:max_l)

! ------------------------------------------------------------------------
!  [(i+1)0|k0]⁽⁰⁾ = (P-A)[i0|k0]⁽⁰⁾ + (W-P)·[i0|k0]⁽¹⁾
!                   + 0.5·i/ζ·([(i-1)0|k0]⁽⁰⁾ - ρ/ζ[(i-1)0|k0]⁽¹⁾)
!                   + 0.5·k/(ζ+η)·[i0|(k-1)0]⁽¹⁾
! ------------------------------------------------------------------------
!  [i0|(k+1)0]⁽⁰⁾ = (Q-C)[i0|k0]⁽⁰⁾ + (W-Q)·[i0|k0]⁽¹⁾
!                   + 0.5·k/η·([i0|(k-1)0]⁽⁰⁾ - ρ/η[i0|(k-1)0]⁽¹⁾)
!                   + 0.5·i/(ζ+η)·[(i-1)0|k0]⁽¹⁾

   tci = 0.0_wp

! ------------------------------------------------------------------------
!  [s|s]⁽⁰⁾

! ------------------------------------------------------------------------
!  [p|s]⁽⁰⁾ = (P-A)·[s|s]⁽⁰⁾ + (W-P)·[s|s]⁽¹⁾

! ------------------------------------------------------------------------
!  [d|s]⁽⁰⁾ = (P-A)·[p|s]⁽⁰⁾ + (W-P)·[p|s]⁽¹⁾
!             + 0.5/ζ·([s|s]⁽⁰⁾ - ρ/ζ[s|s]⁽¹⁾)
!           = (P-A)²·[s|s]⁽⁰⁾ + 2(P-A)(W-P)·[s|s]⁽¹⁾ + (W-P)²·[s|s]⁽²⁾
!             + 0.5/ζ·([s|s]⁽⁰⁾ - ρ/ζ[s|s]⁽¹⁾)

! ------------------------------------------------------------------------
!  [p|p]⁽⁰⁾ = (Q-C)·[p|s]⁽⁰⁾ + (W-Q)·[p|s]⁽¹⁾ + 0.5/(ζ+η)·[s|s]⁽¹⁾
!           = (P-A)(Q-C)·[s|s]⁽⁰⁾ + (W-P)(Q-C)·[s|s]⁽¹⁾
!             + (P-A)(W-Q)·[s|s]⁽¹⁾ + (W-P)(W-Q)·[s|s]⁽²⁾
!             + 0.5/(ζ+η)·[s|s]⁽¹⁾

! ------------------------------------------------------------------------
!  [d|p]⁽⁰⁾ = (Q-C)·[d|s]⁽⁰⁾ + (W-Q)·[d|s]⁽¹⁾ + 1/(ζ+η)·[p|s]⁽¹⁾
!           = (P-A)²(Q-C)·[s|s]⁽⁰⁾ + (W-P)²(W-Q)·[s|s]⁽³⁾
!             + {(P-A)²(W-Q) + 2(P-A)(W-P)(Q-C)}·[s|s]⁽¹⁾
!             + {(W-P)²(Q-C) + 2(P-A)(W-P)(W-Q)}·[s|s]⁽²⁾
!             + 0.5/ζ·(Q-C)·{[s|s]⁽⁰⁾ - ρ/ζ[s|s]⁽¹⁾}
!             + 0.5/ζ·(W-Q)·{[s|s]⁽¹⁾ - ρ/ζ[s|s]⁽²⁾}
!             + 1/(ζ+η)·{(P-A)·[s|s]⁽¹⁾ + (W-P)·[s|s]⁽²⁾}

! ------------------------------------------------------------------------
!  [d|d]⁽⁰⁾ = (Q-C)[d|p]⁽⁰⁾ + (W-Q)·[d|p]⁽¹⁾ + 1/(ζ+η)·[p|p]⁽¹⁾
!             + 0.5/η·([d|s]⁽⁰⁾ - ρ/η[d|s]⁽¹⁾)
!           = (Q-C){(P-A)²(Q-C)·[s|s]⁽⁰⁾ + (W-P)²(W-Q)·[s|s]⁽³⁾
!             + {(P-A)²(W-Q) + 2(P-A)(W-P)(Q-C)}·[s|s]⁽¹⁾
!             + {(W-P)²(Q-C) + 2(P-A)(W-P)(W-Q)}·[s|s]⁽²⁾
!             + 0.5/ζ·(Q-C)·{[s|s]⁽⁰⁾ - ρ/ζ[s|s]⁽¹⁾}
!             + 0.5/ζ·(W-Q)·{[s|s]⁽¹⁾ - ρ/ζ[s|s]⁽²⁾}
!             + 1/(ζ+η)·{(P-A)·[s|s]⁽¹⁾ + (W-P)·[s|s]⁽²⁾}}
!             + (W-Q)·{(P-A)²(Q-C)·[s|s]⁽¹⁾ + (W-P)²(W-Q)·[s|s]⁽⁴⁾
!             + {(P-A)²(W-Q) + 2(P-A)(W-P)(Q-C)}·[s|s]⁽²⁾
!             + {(W-P)²(Q-C) + 2(P-A)(W-P)(W-Q)}·[s|s]⁽³⁾
!             + 0.5/ζ·(Q-C)·{[s|s]⁽¹⁾ - ρ/ζ[s|s]⁽²⁾}
!             + 0.5/ζ·(W-Q)·{[s|s]⁽²⁾ - ρ/ζ[s|s]⁽³⁾}
!             + 1/(ζ+η)·{(P-A)·[s|s]⁽²⁾ + (W-P)·[s|s]⁽³⁾}}
!             + 1/(ζ+η)·{(P-A)(Q-C)·[s|s]⁽¹⁾ + (W-P)(Q-C)·[s|s]⁽²⁾
!             + (P-A)(W-Q)·[s|s]⁽²⁾ + (W-P)(W-Q)·[s|s]⁽³⁾
!             + 0.5/(ζ+η)·[s|s]⁽²⁾}
!             + 0.5/η·({(P-A)²·[s|s]⁽⁰⁾ + 2(P-A)(W-P)·[s|s]⁽¹⁾
!             + (W-P)²·[s|s]⁽²⁾ + 0.5/ζ·([s|s]⁽⁰⁾ - ρ/ζ[s|s]⁽¹⁾)}
!             - ρ/η{(P-A)²·[s|s]⁽¹⁾ + 2(P-A)(W-P)·[s|s]⁽²⁾ + (W-P)²·[s|s]⁽³⁾
!             + 0.5/ζ·([s|s]⁽¹⁾ - ρ/ζ[s|s]⁽²⁾)}

end subroutine erint_os2

subroutine stvint_os2(s,t,v,ityp,jtyp,efact,eab,oab,xab,rab,r_p,ap,bp, &
           &         nat,chrg,xyz)
   use precision, only : wp => dp
   implicit none
!* integrals
   real(wp),intent(inout) :: s
   real(wp),intent(inout) :: t
   real(wp),intent(inout) :: v
!* data
   integer, intent(in) :: ityp ! s  px  py  pz
   integer, intent(in) :: jtyp ! 1  2   3   4
   real(wp),intent(in) :: efact ! preexponential factor
   real(wp),intent(in) :: eab ! product exponent
   real(wp),intent(in) :: oab ! one over exponent
   real(wp),intent(in) :: xab ! reduced exponent
   real(wp),intent(in) :: rab ! distance AB
   real(wp),intent(in) :: r_p(3) ! aufpunkt of product gaussian
   real(wp),intent(in) :: ap(3) ! (A-P)
   real(wp),intent(in) :: bp(3) ! (B-P)
!* molecular data for electron nuclear interaction
   integer, intent(in) :: nat
   real(wp),intent(in) :: chrg(nat)
   real(wp),intent(in) :: xyz(3,nat)

   integer,parameter  :: max_l = 3 ! up to p-functions

!* temporary variables
   integer  :: i,j,k
   real(wp) :: tmp
   real(wp) :: stmp(0:max_l,0:max_l)
   real(wp) :: ttmp(0:max_l,0:max_l)
   real(wp) :: vtmp(0:max_l,0:max_l)
   real(wp) :: rcp
   real(wp) :: cp(3)

   stmp = 0.0_wp
   ttmp = 0.0_wp
   vtmp = 0.0_wp

   i = ityp-1
   j = jtyp-1

!····················································
!· s-s
!····················································
   stmp(0,0) = sqrt(pi*oab)**3
   ttmp(0,0) = xab*(3.0_wp-2.0_wp*xab*rab)*stmp(0,0)
   tmp = 0.0_wp
   do k = 1, nat
      rcp = sum( (r_p-xyz(:,k))**2 )
      tmp = tmp + chrg(k)*boysf(0,eab*rcp)
   enddo
   vtmp(0,0) = tpi*oab*tmp
   if ((ityp.eq.0).and.(jtyp.eq.0)) goto 1

!····················································
!· p-s
!····················································
   if (i.ne.0) then
   stmp(i,0) = ap(i)*stmp(0,0)
   ttmp(i,0) = ap(i)*ttmp(0,0) + 2*xab*stmp(i,0)
   tmp = 0.0_wp
   do k = 1, nat
      cp = xyz(:,k) - r_p
      rcp = sum( cp**2 )
      tmp = tmp + chrg(k)*cp(i)*boysf(1,eab*rcp)
   enddo
   vtmp(i,0) = tpi*oab*tmp + ap(i)*vtmp(0,0)
   endif

!····················································
!· s-p
!····················································
   if (j.ne.0) then
   stmp(0,j) = bp(j)*stmp(0,0)
   ttmp(0,j) = bp(j)*ttmp(0,0) + 2*xab*stmp(0,j)
   tmp = 0.0_wp
   do k = 1, nat
      cp = xyz(:,k) - r_p
      rcp = sum( cp**2 )
      tmp = tmp + chrg(k)*cp(j)*boysf(1,eab*rcp)
   enddo
   vtmp(0,j) = tpi*oab*tmp + bp(j)*vtmp(0,0)
   endif

   if ((ityp.eq.0).or.(jtyp.eq.0)) goto 1

!····················································
!· p-p
!····················································
   stmp(i,j) = ap(i)*stmp(0,j)
   ttmp(i,j) = ap(i)*ttmp(0,j) + 2*xab*stmp(i,j)
   tmp = 0.0_wp
   do k = 1, nat
      cp = xyz(:,k) - r_p
      rcp = sum( cp**2 )
      tmp = tmp + chrg(k)*cp(i)*cp(j)*boysf(2,eab*rcp)
   enddo
   vtmp(i,j) = tpi*oab*tmp - ap(i)*bp(j)*ttmp(0,0) &
   &           + ap(i)*ttmp(0,j) + bp(j)*ttmp(i,0)

   if (ityp.ne.jtyp) goto 1

   stmp(i,j) = stmp(i,j) + 0.5_wp*oab*stmp(0,0)
   ttmp(i,j) = ttmp(i,j) + 0.5_wp*oab*xab*(5.0_wp-2.0_wp*xab*rab)*stmp(0,0)
   tmp = 0.0_wp
   do k = 1, nat
      cp = xyz(:,k) - r_p
      rcp = sum( cp**2 )
      tmp = tmp + chrg(k)*(0.5_wp*oab*boysf(0,eab*rcp)  &
      &           - 0.5_wp*oab*boysf(1,eab*rcp))
   enddo
   vtmp(i,j) = tpi*oab*tmp + vtmp(i,j)


!····················································
 1 continue

   s = s + efact*stmp(i,j)
   t = t + efact*ttmp(i,j)
   v = v + efact*vtmp(i,j)

end subroutine stvint_os2

! ----------------------------------------------------------------------
!  internal helper functions, maybe move to misc later
! ----------------------------------------------------------------------

!* uses the gaussian product theorem to calculate the center of
!  the gaussian formed by the product of two gaussians
pure subroutine gpt(r_a,alp,r_b,bet,rab,r_p,eab,efact)
   use precision, only : wp => dp
   implicit none
   real(wp),intent(in)  :: r_a(3),alp
   real(wp),intent(in)  :: r_b(3),bet
   real(wp),intent(in)  :: rab ! no need to calculate this here
   real(wp),intent(out) :: r_p(3),eab
   real(wp),intent(out) :: efact
   real(wp) :: oab,est

   intrinsic :: exp

   eab = alp+bet
   oab = 1.0_wp/eab
   r_p = (alp*r_a+bet*r_b)*oab
   est = rab*alp*bet*oab
   if (est.gt.i_thr) then
      efact = 0.0_wp
   else
      efact  = exp(-est)
   endif

end subroutine gpt

!* boys function for m=0
elemental function boysf0(arg) result(boys)
   use precision, only : wp => dp
   implicit none
   real(wp) :: boys
   real(wp),intent(in) :: arg

   intrinsic :: sqrt,erf

   if (arg.lt.t_thr) then
      boys = 1.0_wp-oth*arg
   else
      boys = 0.5_wp*sqrt(pi/arg)*erf(sqrt(arg))
   endif

end function boysf0

!* boys function driver
elemental function boysf(m,arg) result(boys)
   use precision, only : wp => dp
   implicit none
   integer, intent(in) :: m
   real(wp),intent(in) :: arg
   real(wp) :: boys
   real(wp) :: val
   integer  :: i

   intrinsic :: exp,sqrt,sum

   if (arg.gt.l_thr) then
      ! asymtotic formula for large values
      boys = sqrt(pi/arg**(2*m+1))/(2**(m+1))*dfactorial(m+2)
   else
      !boys = boysf_rec(m,arg)
      !boys = boysf_iter(m,arg)
      ! I use the grid based approach for the boys function,
      ! its systematic improvable, fast and has less numerical noise
      ! then recursive or iterative implementations
      boys = boysf_grid(m,arg)
   endif

end function boysf

!* boys function evaluated on a grid (see Helgaker2000, p. 367, eq. 9.8.12)
elemental function boysf_grid(m,arg) result(boys)
   use precision, only : wp => dp
   implicit none
   integer, intent(in) :: m
   real(wp),intent(in) :: arg
   real(wp) :: boys
   real(wp) :: delta
   integer  :: i,root

   ! tabl and fgrid are set elsewhere !!here
!  include 'boysf_grid.f90'

   intrinsic :: nint,sum

!  find the nearest gridpoint
   root = nint(arg/tabl)
!  get distance to nearest gridpoint
   delta = arg-root*tabl

!  expand with taylor series at nearest gridpoint
   boys = sum((/(fgrid(m+i,root)*(-delta)**i*ofactorial(i+1),i=5,0,-1)/))

end function boysf_grid

!* hermite coulomb integral based on boys function
!  called R-function by T. Helgaker
elemental function boysr(m,arg,eab) result(boys)
   use precision, only : wp => dp
   implicit none
   integer, intent(in) :: m
   real(wp),intent(in) :: arg
   real(wp),intent(in) :: eab
   real(wp) :: boys
   real(wp) :: val

   boys = (-2.0_wp*eab)**m * boysf(m,arg)

end function boysr

!* one center integral in one cartesian direction
elemental function olap(l,oab)
   use precision, only : wp => dp
   implicit none
   integer, intent(in) :: l
   real(wp),intent(in) :: oab
   real(wp) :: olap
   integer  :: lh

   intrinsic :: mod,sqrt
   
   if (mod(l,2).ne.0) then
   !* zero due to symmetry
      olap = 0.0_wp
   else
      lh = l/2
      olap = sqrt(pi*oab) * (0.5_wp*oab)**lh * dfactorial(lh+1)
   endif

end function olap

!* binomial coeffient (n choose k)
!  now a nice trick from wikipedia,
!  n choose l can be expressed as Π(i=1→l) (n-l+i)/i
!  this formula should garantee that there is no problem
!  with integer overflows and stuff, so I use this one.
!  It is implemented using an implicit do and an array contructor
!  feed to the intrinsic product function, that can crunch arrays
elemental function binom(n,k) result(b)
   integer,intent(in) :: n,k
   integer :: b
   integer :: i,l

   intrinsic :: product

!  n choose n-k is always equal to n choose k,
!  so if k greater n/2 we calculate n choose k instead of n choose n-k
   if ((2*k).gt.n) then
      l = n-k
   else
      l = k
   endif

   b = product((/(((n-l+i)/i),i=1,k)/))

end function binom

!* factorial (n!)
!  straight forward implementation of the factorial, the implicit
!  do expands from 1 to n inside the array constructor, which is
!  evaluated by FORTRAN's product function and returns the factorial.
!  This implementation will lead to an integer overflow, maybe for n>20.
elemental function factorial(n) result(f)
   integer,intent(in) :: n
   integer :: f
   integer :: i

   intrinsic :: product

   f = product((/1,(i,i=1,n)/))

end function factorial

! ----------------------------------------------------------------------
!  Integral transformation from AO to MO basis
! ----------------------------------------------------------------------

subroutine teitrafo_N8(nbf,eri,C)
   use precision, only : wp => dp
   use misc,      only : idx
   implicit none
   integer, intent(in)    :: nbf
   real(wp),intent(in)    :: C(nbf,nbf)
   real(wp),intent(inout) :: eri(nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2)
   real(wp),allocatable   :: temp(:)
   real(wp) :: dum
   integer :: i,j,k,l,ii,jj,kk,ll,ij,ji,kl,lk

   allocate( temp(nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2), source = eri )

!!omp parallel private(i,j,k,l,ii,jj,kk,ll,dum) shared(eri)
!!omp do
   do i = 1,nbf; do j = 1,nbf; ij=idx(i,j)
   do k = 1,nbf; do l = 1,nbf; kl=idx(k,l)
      dum = 0.0_wp
      do ii = 1,nbf; do jj = 1,nbf; ji=idx(ii,jj)
      do kk = 1,nbf; do ll = 1,nbf; lk=idx(kk,ll)
         dum = dum+C(ii,i)*C(jj,j)*C(kk,k)*C(ll,l)*temp(idx(ji,lk))
      enddo; enddo; enddo; enddo
      eri(idx(ij,kl)) = dum
   enddo; enddo; enddo; enddo
!!omp enddo
!!omp endparallel

end subroutine teitrafo_N8

!* transform ee-repulsion integrals from ao basis to mo basis
!  use Θ(N⁵) scaling integral transformation with lots of
!  matrix multiplication.
subroutine teitrafo(nbf,eri,C)
   use precision, only : wp => dp
   use blas95,    only : gemm
   use misc,      only : idx
   implicit none
   integer, intent(in)    :: nbf
   real(wp),intent(in)    :: C(nbf,nbf)
   real(wp),intent(inout) :: eri(nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2)

   real(wp) :: temp(nbf*(nbf+1)/2,nbf*(nbf+1)/2),tmp1(nbf,nbf),tmp2(nbf,nbf)
   integer :: i,j,k,l,ij,kl,ijkl

   temp=0.0_wp

!$omp parallel private(i,j,k,l,ij,kl,ijkl,tmp1,tmp2) shared(temp)
!$omp do schedule(dynamic)
   do i = 1, nbf
      do j = 1, i
         ij = idx(i,j)
         do k = 1, nbf
            do l = 1, k
               kl = idx(k,l)
               ijkl = idx(ij,kl)
               tmp1(k,l) = eri(ijkl)
               tmp1(l,k) = eri(ijkl)
            enddo
         enddo
         call gemm(C,tmp1,tmp2,transa='t')
         call gemm(tmp2,C,tmp1)
         do k = 1, nbf
            do l = 1, k
               kl = idx(k,l)
               temp(kl,ij) = tmp1(k,l)
            enddo
         enddo
      enddo
   enddo
!$omp enddo
!$omp endparallel
   eri = 0.0_wp
!$omp parallel private(i,j,k,l,ij,kl,ijkl,tmp1,tmp2) shared(eri)
!$omp do schedule(dynamic)
   do k = 1, nbf
      do l = 1, k
         kl = idx(k,l)
         do i = 1, nbf
            do j = 1, i
               ij = idx(i,j)
               tmp1(i,j) = temp(kl,ij)
               tmp1(j,i) = temp(kl,ij)
            enddo
         enddo
         call gemm(C,tmp1,tmp2,transa='t')
         call gemm(tmp2,C,tmp1)
         do i = 1, nbf
            do j = 1, i
               ij = idx(i,j)
               ijkl = idx(kl,ij)
               eri(ijkl) = tmp1(i,j)
            enddo
         enddo
      enddo
   enddo
!$omp enddo
!$omp endparallel

end subroutine teitrafo

subroutine onetrafo(nbf,M,C)
   use precision, only : wp => dp
   use blas95,    only : gemm
   integer, intent(in)    :: nbf
   real(wp),intent(in)    :: C(nbf,nbf)
   real(wp),intent(inout) :: M(nbf,nbf)

   real(wp) :: tmp(nbf,nbf)

   call gemm(C,M,tmp,transa='t')
   call gemm(tmp,C,M)

end subroutine onetrafo

end module ints
