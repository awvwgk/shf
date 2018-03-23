module ints
   use precision, only : wp => dp
   implicit none

   private :: oneint,twoint,boysf0

   real(wp),parameter :: r_thr = 2500._wp ! S.lt.10e-9 for r.gt.50 bohr
   real(wp),parameter :: i_thr = 25._wp ! integral cut threshold
   real(wp),parameter :: s_thr = 1e-9_wp ! Schwartz inequaltity threshold
   real(wp),parameter :: t_thr = 1.0e-10_wp ! switch to taylor expansion

   real(wp),parameter :: pi = 3.14159265358979_wp
   real(wp),parameter :: tpi = 6.283185307179580_wp
   real(wp),parameter :: twopi25 = 34.98683665524963_wp
   real(wp),parameter :: tosqpi = 2.0_wp/sqrt(pi)
   real(wp),parameter :: pi2 = 6.36619772367582e-1_wp
   real(wp),parameter :: pi3 = pi**3
   real(wp),parameter :: oth = 1.0_wp/3.0_wp

contains

!* driver for calculation of integrals
subroutine integrals(nat,nbf,at,xyz,zeta,aoc,ng,ityp,S,V,T,eri)
   use precision, only : wp => dp
   use misc,      only : idx
   use stong,     only : slater
   implicit none
   integer, intent(in)  :: nat,nbf,ng(nbf),ityp(nbf)
   integer, intent(in)  :: at(nat),aoc(2,nat)
   real(wp),intent(in)  :: xyz(3,nat),zeta(nbf)
   real(wp),intent(out) :: S(nbf,nbf),V(nbf,nbf),T(nbf,nbf)
!  real(wp),intent(out) :: eri(nbf,nbf,nbf,nbf)
!* to save some memory this packing is possible
!  real(wp),intent(out) :: S(nbf*(nbf+1)/2),V(nbf*(nbf+1)/2),T(nbf*(nbf+1)/2)
   real(wp),intent(out) :: eri((nbf*(nbf+1)/2)*(nbf*(nbf+1)/2+1)/2)

   integer,parameter   :: maxg=6
   real(wp) :: ci(maxg),cj(maxg),alpi(maxg),alpj(maxg)
   real(wp) :: ck(maxg),cl(maxg),alpk(maxg),alpl(maxg)
   real(wp) :: sdum,tdum,vdum,eridum
   real(wp) :: chrg(nat)

   integer :: i,j,k,l,ii,jj,kk,ll,limit
   integer :: ishtyp,jshtyp,kshtyp,lshtyp
   integer :: ij,kl,ijkl

   chrg = dble(at)

!$omp parallel private(i,j,k,l,ii,jj,kk,ll,alpi,alpj,alpk,alpl, &
!$omp          &       limit,ij,kl,ijkl,ci,cj,ck,cl,sdum,vdum,tdum, &
!$omp          &       eridum) shared(s,v,t,eri)
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
         do k = 1, i
            if (i.eq.k) then; limit=j; else; limit=i; endif
            do l = 1, limit
               do ii = aoc(1,i), aoc(2,i)
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

end subroutine integrals

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
   real(wp) :: rab,ab,eab,oab,eabm,sqm,est
   real(wp) :: s00,fact,rcp,r_p(3),cab

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
         eabm = alp(i)*bet(j)*oab
         est = rab*eabm
         if (est.gt.i_thr) cycle ! estimate contribution, ommit if too small
         ab = exp(-est)
         s00 = cab*ab*sqrt(pi*oab)**3

!        overlap
         sab = sab+s00

!        kinetic energy
         tab = tab + eabm*(3.0_wp-2.0_wp*est)*s00

!        nuclear attraction
         fact = cab*tpi*oab*ab
         r_p = (alp(i)*r_a+bet(j)*r_b)*oab
         do k = 1, nat
            rcp = sum( (r_p-xyz(:,k))**2 )
            vab  = vab - fact*chrg(k)*boysf0(eab*rcp)
         enddo

      enddo
   enddo

end subroutine oneint

!  two-electron repulsion integral [ab|cd] over s-functions
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
   real(wp) :: ab(npa,npb),cd(npc,npd),abcd,pq,sqpq
   real(wp) :: qabab,qcdcd

   tei = 0.0_wp
   qabab = 0.0_wp
   qcdcd = 0.0_wp

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
         if (est.gt.i_thr) then
            ab(i,j) = 0.0_wp
         else
            ab(i,j) = exp(-est) ! might be used later

            qabab = qabab + cab*ab(i,j) * sqrt(twopi25/(sqrt(2*eab**5)))
         endif
      enddo
   enddo
   do k = 1, npc
      do l = 1, npd
         ccd = ck(k)*cl(l)
         ecd = gam(k)+del(l)
         ocd = 1.0_wp/ecd
         est = gam(k)*del(l)*rcd*ocd
         if (est.gt.i_thr) then
            cd(k,l) = 0.0_wp
         else
            cd(k,l) = exp(-est) ! might be used later

            qcdcd = qcdcd + ccd*cd(k,l) * sqrt(twopi25/(sqrt(2*ecd**5)))
         endif
      enddo
   enddo
!* prescreening by using Schwarz inequality
!  if the upper bound of the tei falls below a threshold ommit calculation
   if ((qabab*qcdcd).lt.s_thr) return

   do i = 1, npa
      do j = 1, npb
         cab = ci(i)*cj(j)
         eab = alp(i)+bet(j)
         oab = 1.0_wp/eab

!        new gaussian at r_p
         r_p = (alp(i)*r_a+bet(j)*r_b)*oab

         do k = 1, npc
            do l = 1, npd
               ccd = ck(k)*cl(l)
               ecd = gam(k)+del(l)
               ocd = 1.0_wp/ecd

!              new gaussian at r_q
               r_q = (gam(k)*r_c+del(l)*r_d)*ocd

!              we already have calculated the prefactors
               abcd = ab(i,j)*cd(k,l)

!              distance between product gaussians
               rpq = sum( (r_p-r_q)**2 )

               epq = eab*ecd
               eabcd = eab+ecd

               pq = rpq*epq/eabcd
               tei = tei + cab*ccd*abcd * twopi25/(epq*sqrt(eabcd)) &
               &           * boysf0(pq)

            enddo
         enddo
      enddo
   enddo

end subroutine twoint

!* f0 function/boys function for m=0
elemental function boysf0(arg)
   use precision, only : wp => dp
   implicit none
   real(wp) :: boysf0
   real(wp),intent(in) :: arg

   if (arg.lt.t_thr) then
      boysf0=1.0_wp-oth*arg
   else
      boysf0=0.5_wp*sqrt(pi/arg)*erf(sqrt(arg))
   endif

end function boysf0

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

end module ints
