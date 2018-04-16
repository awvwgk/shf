module cc

   interface ccd
   module procedure solve_ccd
   end interface ccd

   interface ccsd
   module procedure solve_ccsd
   end interface ccsd

   interface build_T1
   module procedure build_T1_ccsd
   end interface build_T1

   interface build_T2
   module procedure build_T2_ccd
   module procedure build_T2_ccsd
   end interface build_T2

   interface build_Fki
   module procedure build_Fki_ccd
   module procedure build_Fki_ccsd
   end interface

   interface build_Fac
   module procedure build_Fac_ccd
   module procedure build_Fac_ccsd
   end interface

   interface build_Fkc
   module procedure build_Fkc_ccsd
   end interface

   interface build_Wklij
   module procedure build_Wklij_ccd
   module procedure build_Wklij_ccsd
   end interface

   interface build_Wabcd
   module procedure build_Wabcd_ccsd
   end interface

   interface build_Wkbcj
   module procedure build_Wkbcj_ccd
   module procedure build_Wkbcj_ccsd
   end interface

contains

!  ------------------------------------------------------------------------
!  prep stuff
!  ------------------------------------------------------------------------
subroutine chem2phys(nbf,eri,tei)
   use precision, only : wp => dp
   use misc,      only : idx
   implicit none
   integer, intent(in)  :: nbf
   real(wp),intent(in)  :: eri((nbf*(nbf+1)/2)*((nbf*(nbf+1)/2)+1)/2)
   real(wp),intent(out) :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   integer  :: i,j,k,l,ik,il,jk,jl,ikjl,iljk
   integer  :: ih,jh,kh,lh,im,jm,km,lm
   real(wp) :: eri1,eri2
   tei = 0.0_wp
   do i = 1, 2*nbf
      ih = i/2
      im = mod(i,2)
      do j = 1, 2*nbf
         jh = j/2
         jm = mod(j,2)
         do k = 1, 2*nbf
            kh = k/2
            km = mod(k,2)
            do l = 1, 2*nbf
               lh = l/2
               lm = mod(l,2)
               ik = idx(ih+im,kh+km)
               jl = idx(jh+jm,lh+lm)
               il = idx(ih+im,lh+lm)
               jk = idx(jh+jm,kh+km)
               ikjl = idx(ik,jl)
               iljk = idx(il,jk)
               if ((im.eq.km).and.(jm.eq.lm)) then
                  eri1 = eri(ikjl)
               else
                  eri1 = 0.0_wp
               endif
               if ((im.eq.lm).and.(jm.eq.km)) then
                  eri2 = eri(iljk)
               else
                  eri2 = 0.0_wp
               endif
               tei(i,j,k,l) = eri1 - eri2
            enddo
         enddo
      enddo
   enddo
end subroutine chem2phys
subroutine spinfockian(nbf,nel,F,H,tei)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   real(wp),intent(out) :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: H(nbf,nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   integer  :: i,j,k,ih,jh,im,jm
   real(wp) :: g

   F = 0.0_wp

   do i = 1, 2*nbf
      ih = i/2
      im = mod(i,2)
      do j = 1, 2*nbf
         jh = j/2
         jm = mod(j,2)
         if (im.eq.jm) F(j,i) = H(ih+im,jh+jm)
         g = 0.0_wp
         do k = 1, nel
            g = g + tei(j,k,i,k)
         enddo
         F(j,i) = F(j,i) + g
      enddo
   enddo

end subroutine spinfockian

!  ------------------------------------------------------------------------
!  CCSD driver
!  ------------------------------------------------------------------------
subroutine solve_ccsd(nbf,nel,F,tei,ethr,acc,maxiter,e)
   use iso_fortran_env, only : id => output_unit
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: ethr
   character(len=*),intent(in) :: acc
   integer, intent(in)  :: maxiter
   real(wp),intent(out) :: e

   integer  :: nvir
   integer  :: iter
   integer  :: err
   integer  :: i,j,k,l
   integer  :: a,b,c,d
   real(wp) :: eold
   real(wp) :: rmsd
   real(wp),allocatable :: D1(:,:)
   real(wp),allocatable :: T1(:,:)
   real(wp),allocatable :: T1_save(:,:)
   real(wp),allocatable :: D2(:,:,:,:)
   real(wp),allocatable :: T2(:,:,:,:)
   real(wp),allocatable :: T2_save(:,:,:,:)
   real(wp),allocatable :: Fac(:,:)
   real(wp),allocatable :: Fki(:,:)
   real(wp),allocatable :: Fkc(:,:)
   real(wp),allocatable :: Wklij(:,:,:,:)
   real(wp),allocatable :: Wabcd(:,:,:,:)
   real(wp),allocatable :: Wkbcj(:,:,:,:)

   write(id,'(a)')
   write(id,'(72(''-''))')
   write(id,'('' iterative Coupled Cluster Singles Doubles calculation'')')

   nvir = 2*nbf - nel

   allocate( D1(nvir,nel),  &
   &         T1(nvir,nel),  &
   &         T1_save(nvir,nel),  &
   &         D2(nvir,nvir,nel,nel),  &
   &         T2(nvir,nvir,nel,nel),  &
   &         T2_save(nvir,nvir,nel,nel),  &
   &         Fac(nvir,nvir),  &
   &         Fki(nel,nel),  &
   &         Fkc(nel,nvir),  &
   &         Wklij(nel,nel,nel,nel),  &
   &         Wabcd(nvir,nvir,nvir,nvir),  &
   &         Wkbcj(nel,nvir,nvir,nel),  &
   &         source = 0.0_wp, stat = err )
   if (err.ne.0) call raise('E','could not allocate for CCSD calculation')

   do i = 1, nel
      do a = 1, nvir
         D1(a,i) = 1.0_wp/( F(i,i) - F(a+nel,a+nel) )
      enddo
   enddo
 
   do i = 1, nel
      do j = 1, nel
         do a = 1, nvir
            do b = 1, nvir
               D2(a,b,i,j) = 1.0_wp/( F(i,i) + F(j,j)  &
               &   - F(a+nel,a+nel) - F(b+nel,b+nel) )
               T2(a,b,i,j) = tei(a+nel,b+nel,i,j) * D2(a,b,i,j)
            enddo
         enddo
      enddo
   enddo
   
   eold = ccsd_energy(nbf,nel,nvir,F,tei,T1,T2)
   iter = 0

   e = eold
   
!  write(id,'(72(''-''))')
!  write(id,'('' GUESS MP2 ENERGY'',f53.'//acc//')') eold
!  write(id,'(72(''-''))')
!  write(id,'(a)')

   write(id,'(1(''-'')''[ITER]'')',   advance='no')
   write(id,'(8(''-'')''[E(CCSD)]'')', advance='no')
   write(id,'(7(''-'')''[ΔE(CCSD)]'')',advance='no')
   write(id,'(7(''-'')''[RMS(T)]'')', advance='no')
   write(id,'(16(''-''))')
   write(id,'(x,i5,x,f16.'//acc//')'),iter,eold

   cciter: do
      iter = iter+1
      if (iter.gt.maxiter) call raise('E','CCSD did not converge')
      T1_save = T1
      T2_save = T2

      call build_Fac(nbf,nel,nvir,F,tei,T1,T2,Fac)
      call build_Fki(nbf,nel,nvir,F,tei,T1,T2,Fki)
      call build_Fkc(nbf,nel,nvir,F,tei,T1,T2,Fkc)
      call build_Wklij(nbf,nel,nvir,F,tei,T1,T2,Wklij)
      call build_Wabcd(nbf,nel,nvir,F,tei,T1,T2,Wabcd)
      call build_Wkbcj(nbf,nel,nvir,F,tei,T1,T2,Wkbcj)

      call build_T1(nbf,nel,nvir,F,tei,T1_save,T2_save,  &
           &        Fac,Fki,Fkc,Wklij,Wabcd,Wkbcj,D1,T1)
      call build_T2(nbf,nel,nvir,F,tei,T1_save,T2_save,  &
           &        Fac,Fki,Fkc,Wklij,Wabcd,Wkbcj,D2,T2)

      e = ccsd_energy(nbf,nel,nvir,F,tei,T1,T2)
      rmsd = sqrt( sum( (T1-T1_save)**2 )/real(nel*nvir,wp) )
      rmsd = rmsd + sqrt( sum( (T2-T2_save)**2 ) )/real(nel*nvir,wp)

      write(id,'(x,i5)',advance='no') iter
      write(id,'(x,f16.'//acc//')',advance='no') e
      write(id,'(x,f16.'//acc//')',advance='no') e-eold
      write(id,'(x,f14.'//acc//')',advance='no') rmsd
      write(id,'(a)')

      if (abs(e-eold).lt.ethr) exit cciter

      eold = e
   enddo cciter

   deallocate( D1,D2,T1,T2,T1_save,T2_save,  &
   &           Fac,Fki,Fkc,Wklij,Wabcd,Wkbcj )

   print'(72(''-''))'
   print'('' FINAL CCSD ENERGY'',f52.'//acc//')',e
   print'(72(''-''))'

end subroutine solve_ccsd

!  ------------------------------------------------------------------------
!  CCD driver
!  ------------------------------------------------------------------------
subroutine solve_ccd(nbf,nel,F,tei,ethr,acc,maxiter,e)
   use iso_fortran_env, only : id => output_unit
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: ethr
   character(len=*),intent(in) :: acc
   integer, intent(in)  :: maxiter
   real(wp),intent(out) :: e

   integer  :: nvir
   integer  :: iter
   integer  :: err
   integer  :: i,j,k,l
   integer  :: a,b,c,d
   real(wp) :: eold
   real(wp) :: rmsd
   real(wp),allocatable :: D2(:,:,:,:)
   real(wp),allocatable :: T2(:,:,:,:)
   real(wp),allocatable :: TS(:,:,:,:)
   real(wp),allocatable :: Fac(:,:)
   real(wp),allocatable :: Fki(:,:)
   real(wp),allocatable :: Wklij(:,:,:,:)
   real(wp),allocatable :: Wkbcj(:,:,:,:)

   write(id,'(a)')
   write(id,'(72(''-''))')
   write(id,'('' iterative Coupled Cluster Doubles calculation'')')

   nvir = 2*nbf - nel

   allocate( D2(nvir,nvir,nel,nel),  &
   &         T2(nvir,nvir,nel,nel),  &
   &         TS(nvir,nvir,nel,nel),  &
   &         Fac(nvir,nvir),Fki(nel,nel),  &
   &         Wklij(nel,nel,nel,nel),  &
   &         Wkbcj(nel,nvir,nvir,nel),  &
   &         source = 0.0_wp, stat = err )
   if (err.ne.0) call raise('E','could not allocate for CCD calculation')

   do i = 1, nel
      do j = 1, nel
         do a = 1, nvir
            do b = 1, nvir
               D2(a,b,i,j) = 1.0_wp/( F(i,i) + F(j,j)  &
               &   - F(a+nel,a+nel) - F(b+nel,b+nel) )
               T2(a,b,i,j) = tei(a+nel,b+nel,i,j) * D2(a,b,i,j)
            enddo
         enddo
      enddo
   enddo
   
   eold = ccd_energy(nbf,nel,nvir,tei,T2)
   iter = 0
   
!  write(id,'(72(''-''))')
!  write(id,'('' GUESS MP2 ENERGY'',f53.'//acc//')') eold
!  write(id,'(72(''-''))')
!  write(id,'(a)')

   write(id,'(1(''-'')''[ITER]'')',   advance='no')
   write(id,'(9(''-'')''[E(CCD)]'')', advance='no')
   write(id,'(8(''-'')''[ΔE(CCD)]'')',advance='no')
   write(id,'(7(''-'')''[RMS(T)]'')', advance='no')
   write(id,'(16(''-''))')
   write(id,'(x,i5,x,f16.'//acc//')'),iter,eold

   cciter: do
      iter = iter+1
      if (iter.gt.maxiter) call raise('E','CCD did not converge')
      TS = T2

      call build_Fac(nbf,nel,nvir,F,tei,T2,Fac)
      call build_Fki(nbf,nel,nvir,F,tei,T2,Fki)
      call build_Wklij(nbf,nel,nvir,F,tei,T2,Wklij)
      call build_Wkbcj(nbf,nel,nvir,F,tei,T2,Wkbcj)

      call build_T2(nbf,nel,nvir,F,tei,TS,Fac,Fki,Wklij,Wkbcj,D2,T2)

      e = ccd_energy(nbf,nel,nvir,tei,T2)
      rmsd = sqrt( sum( (T2-TS)**2 ) )/real(nel*nvir,wp)

      write(id,'(x,i5)',advance='no') iter
      write(id,'(x,f16.'//acc//')',advance='no') e
      write(id,'(x,f16.'//acc//')',advance='no') e-eold
      write(id,'(x,f14.'//acc//')',advance='no') rmsd
      write(id,'(a)')
      
      if (abs(e-eold).lt.ethr) exit cciter

      eold = e
   enddo cciter

   deallocate( D2,T2,TS,Fac,Fki,Wklij,Wkbcj )

   print'(72(''-''))'
   print'('' FINAL CCD ENERGY'',f53.'//acc//')',e
   print'(72(''-''))'

end subroutine solve_ccd

!  ------------------------------------------------------------------------
!  CCSD amplitude equations
!  ------------------------------------------------------------------------
pure subroutine build_T1_ccsd(nbf,nel,nvir,F,tei,TS,TD,  &
                &             Fac,Fki,Fkc,Wklij,Wabcd,Wkbcj,D1,T1)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: TS(nvir,nel)
   real(wp),intent(in)  :: TD(nvir,nvir,nel,nel)
   real(wp),intent(in)  :: Fac(nvir,nvir)
   real(wp),intent(in)  :: Fki(nel,nel)
   real(wp),intent(in)  :: Fkc(nel,nvir)
   real(wp),intent(in)  :: Wklij(nel,nel,nel,nel)
   real(wp),intent(in)  :: Wabcd(nvir,nvir,nvir,nvir)
   real(wp),intent(in)  :: Wkbcj(nel,nvir,nvir,nel)
   real(wp),intent(in)  :: D1(nvir,nel)
   real(wp),intent(out) :: T1(nvir,nel)

   integer  :: i,j,k,l
   integer  :: a,b,c,d

   T1 = F(nel+1:2*nbf,1:nel)

   do i = 1, nel
      do a = 1, nvir
         do c = 1, nvir
            T1(a,i) = T1(a,i) + TS(a,i) * Fac(a,c)
         enddo
         do k = 1, nel
            T1(a,i) = T1(a,i) - TS(a,i) * Fki(k,i)
         enddo
         do k = 1, nel
            do c = 1, nvir
               T1(a,i) = T1(a,i) + TD(a,c,i,k) * Fkc(k,c)
            enddo
         enddo
         do l = 1, nel
            do d = 1, nvir
               T1(a,i) = T1(a,i) - TS(d,l) * tei(l,a+nel,i,d+nel)
            enddo
         enddo
         do k = 1, nel
            do c = 1, nvir
               do d = 1, nvir
                  T1(a,i) = T1(a,i) - 0.5_wp * TD(c,d,i,k)  &
                  &                   * tei(k,a+nel,c+nel,d+nel)
               enddo
            enddo
         enddo
         do k = 1, nel
            do l = 1, nel
               do c = 1, nvir
                  T1(a,i) = T1(a,i) - 0.5_wp * TD(a,c,k,l)  &
                  &                   * tei(k,l,c+nel,i)
               enddo
            enddo
         enddo
      enddo
   enddo

   T1 = T1*D1

end subroutine build_T1_ccsd

pure subroutine build_T2_ccsd(nbf,nel,nvir,F,tei,TS,TD,  &
                &             Fac,Fki,Fkc,Wklij,Wabcd,Wkbcj,D2,T2)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: TS(nvir,nel)
   real(wp),intent(in)  :: TD(nvir,nvir,nel,nel)
   real(wp),intent(in)  :: Fac(nvir,nvir)
   real(wp),intent(in)  :: Fki(nel,nel)
   real(wp),intent(in)  :: Fkc(nel,nvir)
   real(wp),intent(in)  :: Wklij(nel,nel,nel,nel)
   real(wp),intent(in)  :: Wabcd(nvir,nvir,nvir,nvir)
   real(wp),intent(in)  :: Wkbcj(nel,nvir,nvir,nel)
   real(wp),intent(in)  :: D2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: T2(nvir,nvir,nel,nel)

   integer  :: i,j,k,l
   integer  :: a,b,c,d
   real(wp) :: tdum1,tdum2

   T2 = tei(nel+1:2*nbf,nel+1:2*nbf,1:nel,1:nel)

   do i = 1, nel
      do j = 1, nel
         do a = 1, nvir
            do b = 1, nvir
               do c = 1, nvir
                  T2(a,b,i,j) = T2(a,b,i,j) + TD(a,c,i,j) * Fac(b,c)
                  T2(a,b,i,j) = T2(a,b,i,j) - TD(b,c,i,j) * Fac(a,c)
                  tdum1 = 0.0_wp
                  tdum2 = 0.0_wp
                  do k = 1, nel
                     tdum1 = tdum1 - 0.5 * TS(b,k)*Fkc(k,c)
                     tdum2 = tdum2 - 0.5 * TS(a,k)*Fkc(k,c)
                  enddo
                  T2(a,b,i,j) = T2(a,b,i,j) + TD(a,c,i,j)*tdum1
                  T2(a,b,i,j) = T2(a,b,i,j) - TD(b,c,i,j)*tdum2
               enddo
               do k = 1, nel
                  T2(a,b,i,j) = T2(a,b,i,j) - TD(a,b,i,k) * Fki(k,j)
                  T2(a,b,i,j) = T2(a,b,i,j) + TD(a,b,j,k) * Fki(k,i)
                  tdum1 = 0.0_wp
                  tdum2 = 0.0_wp
                  do c = 1, nvir
                     tdum1 = tdum1 + 0.5 * TS(c,j)*Fkc(k,c)
                     tdum2 = tdum2 + 0.5 * TS(c,i)*Fkc(k,c)
                  enddo
                  T2(a,b,i,j) = T2(a,b,i,j) - TD(a,b,i,k)*tdum1
                  T2(a,b,i,j) = T2(a,b,i,j) + TD(a,b,j,k)*tdum2
               enddo
               tdum1 = 0.0_wp
               do k = 1, nel
                  do l = 1, nel
                     tdum1 = tdum1 + 0.5_wp * Wklij(k,l,i,j)  &
                     &       * (TD(a,b,k,l) + TS(a,k)*TS(b,l) - TS(b,k)*TS(a,l))
                  enddo
               enddo
               T2(a,b,i,j) = T2(a,b,i,j) + tdum1
               tdum2 = 0.0_wp
               do c = 1, nvir
                  do d = 1, nvir
                     tdum2 = tdum2 + 0.5_wp * Wabcd(a,b,c,d)  &
                     &       * (TD(c,d,i,j) + TS(c,i)*TS(d,j) - TS(d,i)*TS(c,j))
                  enddo
               enddo
               T2(a,b,i,j) = T2(a,b,i,j) + tdum2
               do k = 1, nel
                  do c = 1, nvir
                     T2(a,b,i,j) = T2(a,b,i,j) + (TD(a,c,i,k)*Wkbcj(k,b,c,j)  &
                     &             - TS(c,i)*TS(a,k) * tei(k,b+nel,c+nel,j))
                     T2(a,b,i,j) = T2(a,b,i,j) - (TD(b,c,i,k)*Wkbcj(k,a,c,j)  &
                     &             - TS(c,j)*TS(a,k) * tei(k,b+nel,c+nel,j))
                     T2(a,b,i,j) = T2(a,b,i,j) - (TD(a,c,j,k)*Wkbcj(k,b,c,i)  &
                     &             - TS(c,i)*TS(b,k) * tei(k,b+nel,c+nel,j))
                     T2(a,b,i,j) = T2(a,b,i,j) + (TD(b,c,j,k)*Wkbcj(k,a,c,i)  &
                     &             - TS(c,j)*TS(b,k) * tei(k,a+nel,c+nel,i))
                  enddo
               enddo
               do c = 1, nvir
                  T2(a,b,i,j) = T2(a,b,i,j) + TS(c,i) * tei(a+nel,b+nel,c+nel,j)
                  T2(a,b,i,j) = T2(a,b,i,j) - TS(c,j) * tei(a+nel,b+nel,c+nel,i)
               enddo
               do k = 1, nel
                  T2(a,b,i,j) = T2(a,b,i,j) - TS(a,i) * tei(k,b+nel,i,j)
                  T2(a,b,i,j) = T2(a,b,i,j) + TS(b,j) * tei(k,a+nel,i,j)
               enddo
            enddo
         enddo
      enddo
   enddo

   T2 = T2*D2

end subroutine build_T2_ccsd

!  ------------------------------------------------------------------------
!  CCD amplitude equations
!  ------------------------------------------------------------------------
pure subroutine build_T2_ccd(nbf,nel,nvir,F,tei,TS,Fac,Fki,Wklij,Wkbcj,D2,T2)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: TS(nvir,nvir,nel,nel)
   real(wp),intent(in)  :: Fac(nvir,nvir)
   real(wp),intent(in)  :: Fki(nel,nel)
   real(wp),intent(in)  :: Wklij(nel,nel,nel,nel)
   real(wp),intent(in)  :: Wkbcj(nel,nvir,nvir,nel)
   real(wp),intent(in)  :: D2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: T2(nvir,nvir,nel,nel)

   integer  :: i,j,k,l
   integer  :: a,b,c,d

   T2 = tei(nel+1:2*nbf,nel+1:2*nbf,1:nel,1:nel)

   do i = 1, nel
      do j = 1, nel
         do a = 1, nvir
            do b = 1, nvir
               do c = 1, nvir
                  T2(a,b,i,j) = T2(a,b,i,j) + Fac(a,c) * TS(c,b,i,j)
                  T2(a,b,i,j) = T2(a,b,i,j) - Fac(b,c) * TS(c,a,i,j)
               enddo
               do k = 1, nel
                  T2(a,b,i,j) = T2(a,b,i,j) - Fki(k,i) * TS(a,b,k,j)
                  T2(a,b,i,j) = T2(a,b,i,j) + Fki(k,j) * TS(a,b,k,i)
               enddo
               do k = 1, nel
                  do l = 1, nel
                     T2(a,b,i,j) = T2(a,b,i,j)  &
                     &   + 0.5 * Wklij(k,l,i,j) * TS(a,b,k,l)
                  enddo
               enddo
               do c = 1, nvir
                  do d = 1, nvir
                     T2(a,b,i,j) = T2(a,b,i,j)  &
                     &   + 0.5 * tei(a+nel,b+nel,c+nel,d+nel) * TS(c,d,i,j)
                  enddo
               enddo
               do k = 1, nel
                  do c = 1, nvir
                     T2(a,b,i,j) = T2(a,b,i,j)  &
                     &   + Wkbcj(k,b,c,j) * TS(a,c,i,k)  &
                     &   - Wkbcj(k,a,c,j) * TS(b,c,i,k)  &
                     &   - Wkbcj(k,b,c,i) * TS(a,c,j,k)  &
                     &   + Wkbcj(k,a,c,i) * TS(b,c,j,k)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

   T2 = T2*D2

end subroutine build_T2_ccd

!  ------------------------------------------------------------------------
!  CCSD intermediates
!  ------------------------------------------------------------------------
pure subroutine build_Fki_ccsd(nbf,nel,nvir,F,tei,T1,T2,Fki)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T1(nvir,nel)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: Fki(nel,nel)

   integer  :: i,j,k,l
   integer  :: a,b,c,d

   Fki = F(1:nel,1:nel)
   do k = 1, nel
      do i = 1, nel
         if (k.eq.i) Fki(k,i) = 0.0_wp
         do c = 1, nvir
            Fki(k,i) = Fki(k,i) + 0.5_wp * F(k,c+nel) * T1(c,i)
            do l = 1, nel
               Fki(k,i) = Fki(k,i) + T1(c,l) * tei(k,l,i,c+nel)
               do d = 1, nvir
                  Fki(k,i) = Fki(k,i) + 0.5_wp * tei(k,l,c+nel,d+nel)  &
                  &  * (T2(c,d,i,l) + 0.5_wp*(T1(c,i)*T1(d,l)-T1(c,l)*T1(d,i)))
               enddo
            enddo
         enddo
      enddo
   enddo

end subroutine build_Fki_ccsd

pure subroutine build_Fac_ccsd(nbf,nel,nvir,F,tei,T1,T2,Fac)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T1(nvir,nel)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: Fac(nvir,nvir)

   integer  :: i,j,k,l
   integer  :: a,b,c,d

   Fac = F(nel+1:2*nbf,nel+1:2*nbf)
   do a = 1, nvir
      do c = 1, nvir
         if (c.eq.a) Fac(a,c) = 0.0_wp
         do k = 1, nel
            Fac(a,c) = Fac(a,c) - 0.5_wp * F(k+nel,c+nel) * T1(a,k)
            do d = 1, nvir
               Fac(a,c) = Fac(a,c) + T1(d,k) * tei(k,a+nel,d+nel,c+nel)
               do l = 1, nel
                  Fac(a,c) = Fac(a,c) - 0.5_wp * tei(k,l,c+nel,d+nel)  &
                  &  * (T2(a,d,k,l) + 0.5_wp*(T1(a,k)*T1(d,l)-T1(d,k)*T1(a,l)))
               enddo
            enddo
         enddo
      enddo
   enddo

end subroutine build_Fac_ccsd

pure subroutine build_Fkc_ccsd(nbf,nel,nvir,F,tei,T1,T2,Fkc)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T1(nvir,nel)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: Fkc(nel,nvir)

   integer  :: i,j,k,l
   integer  :: a,b,c,d

   Fkc = F(1:nel,nel+1:2*nbf)
   do k = 1, nel
      do c = 1, nvir
         do l = 1, nvir
            do d = 1, nel
               Fkc(k,c) = Fkc(k,c) + T1(d,l) * tei(k,l,c+nel,d+nel)
            enddo
         enddo
      enddo
   enddo

end subroutine build_Fkc_ccsd

pure subroutine build_Wklij_ccsd(nbf,nel,nvir,F,tei,T1,T2,Wklij)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T1(nvir,nel)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: Wklij(nel,nel,nel,nel)

   integer  :: i,j,k,l
   integer  :: a,b,c,d
   real(wp) :: tklij

   Wklij = tei(1:nel,1:nel,1:nel,1:nel)
   do k = 1, nel
      do l = 1, nel
         do i = 1, nel
            do j = 1, nel
               tklij = 0.0_wp
               do c = 1, nvir
                  tklij = tklij + T1(c,j) * tei(k,l,i,c+nel)
                  tklij = tklij - T1(c,i) * tei(k,l,j,c+nel)
                  do d = 1, nvir
                     tklij = tklij - 0.25_wp * tei(k,l,c+nel,d+nel)  &
                     &  * ( T2(c,d,i,j) + T1(c,i)*T1(d,j) - T1(c,j)*T1(d,i) )
                  enddo
               enddo
               Wklij(k,l,i,j) = Wklij(k,l,i,j) + tklij
            enddo
         enddo
      enddo
   enddo

end subroutine build_Wklij_ccsd

pure subroutine build_Wabcd_ccsd(nbf,nel,nvir,F,tei,T1,T2,Wabcd)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T1(nvir,nel)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: Wabcd(nvir,nvir,nvir,nvir)

   integer  :: i,j,k,l
   integer  :: a,b,c,d
   real(wp) :: tabcd

   Wabcd = tei(nel+1:2*nbf,nel+1:2*nbf,nel+1:2*nbf,nel+1:2*nbf)
   do a = 1, nvir
      do b = 1, nvir
         do c = 1, nvir
            do d = 1, nvir
               tabcd = 0.0_wp
               do k = 1, nel
                  tabcd = tabcd - T1(b,k) * tei(a+nel,k,c+nel,d+nel)
                  tabcd = tabcd + T1(a,k) * tei(b+nel,k,c+nel,d+nel)
                  do l = 1, nel
                     tabcd = tabcd + 0.25_wp * tei(k,l,c+nel,d+nel)  &
                     &  * ( T2(a,b,k,l) + T1(a,k)*T1(b,l) - T1(a,l)*T1(b,k) )
                  enddo
               enddo
               Wabcd(a,b,c,d) = Wabcd(a,b,c,d) + tabcd
            enddo
         enddo
      enddo
   enddo

end subroutine build_Wabcd_ccsd

pure subroutine build_Wkbcj_ccsd(nbf,nel,nvir,F,tei,T1,T2,Wkbcj)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T1(nvir,nel)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: Wkbcj(nel,nvir,nvir,nel)

   integer  :: i,j,k,l
   integer  :: a,b,c,d
   real(wp) :: twbcj

   Wkbcj = tei(1:nel,nel+1:2*nbf,nel+1:2*nbf,1:nel)
   do k = 1, nel
      do b = 1, nvir
         do c = 1, nvir
            do j = 1, nel
               twbcj = 0.0_wp
               do d = 1, nvir
                  twbcj = twbcj + T1(d,j)*tei(k,b+nel,c+nel,d+nel)
               enddo
               do l = 1, nel
                  twbcj = twbcj - T1(b,l)*tei(k,l,c+nel,j)
               enddo
               do d = 1, nvir
                  do l = 1, nel
                     twbcj = twbcj - tei(k,l,c+nel,d+nel)  &
                     &  * ( 0.5_wp*T2(d,b,j,l) + T1(d,j)*T1(b,l) )
                  enddo
               enddo
               Wkbcj(k,b,c,j) = Wkbcj(k,b,c,j) + twbcj
            enddo
         enddo
      enddo
   enddo

end subroutine build_Wkbcj_ccsd

!  ------------------------------------------------------------------------
!  CCD intermediates
!  ------------------------------------------------------------------------
pure subroutine build_Fki_ccd(nbf,nel,nvir,F,tei,T2,Fki)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: Fki(nel,nel)

   integer  :: i,j,k,l
   integer  :: a,b,c,d

   Fki = F(1:nel,1:nel)
   do k = 1, nel
      do i = 1, nel
         if (k.eq.i) Fki(k,i) = 0.0_wp
         do l = 1, nel
            do c = 1, nvir
               do d = 1, nvir
                  Fki(k,i) = Fki(k,i) + 0.5_wp * tei(k,l,c+nel,d+nel)  &
                  &                     * T2(c,d,i,l)
               enddo
            enddo
         enddo
      enddo
   enddo

end subroutine build_Fki_ccd

pure subroutine build_Fac_ccd(nbf,nel,nvir,F,tei,T2,Fac)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: Fac(nvir,nvir)

   integer  :: i,j,k,l
   integer  :: a,b,c,d

   Fac = F(nel+1:2*nbf,nel+1:2*nbf)
   do a = 1, nvir
      do c = 1, nvir
         if (c.eq.a) Fac(a,c) = 0.0_wp
         do k = 1, nel
            do l = 1, nel
               do d = 1, nvir
                  Fac(a,c) = Fac(a,c) - 0.5_wp * tei(k,l,c+nel,d+nel)  &
                  &                     * T2(a,d,k,l)
               enddo
            enddo
         enddo
      enddo
   enddo

end subroutine build_Fac_ccd

pure subroutine build_Wklij_ccd(nbf,nel,nvir,F,tei,T2,Wklij)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: Wklij(nel,nel,nel,nel)

   integer  :: i,j,k,l
   integer  :: a,b,c,d

   Wklij = tei(1:nel,1:nel,1:nel,1:nel)
   do k = 1, nel
      do l = 1, nel
         do i = 1, nel
            do j = 1, nel
               do c = 1, nvir
                  do d = 1, nvir
                     Wklij(k,l,i,j) = Wklij(k,l,i,j)  &
                     &   + 0.5 * tei(k,l,c+nel,d+nel) * T2(c,d,i,j)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

end subroutine build_Wklij_ccd

pure subroutine build_Wkbcj_ccd(nbf,nel,nvir,F,tei,T2,Wkbcj)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp),intent(out) :: Wkbcj(nel,nvir,nvir,nel)

   integer  :: i,j,k,l
   integer  :: a,b,c,d

   Wkbcj = tei(1:nel,nel+1:2*nbf,nel+1:2*nbf,1:nel)
   do k = 1, nel
      do b = 1, nvir
         do c = 1, nvir
            do j = 1, nel
               do l = 1, nel
                  do d = 1, nvir
                     Wkbcj(k,b,c,j) = Wkbcj(k,b,c,j)  &
                     &   + 0.5 * tei(k,l,c+nel,d+nel) * T2(d,b,j,l)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

end subroutine build_Wkbcj_ccd

!  ------------------------------------------------------------------------
!  CCSD energy expression
!  ------------------------------------------------------------------------
function ccsd_energy(nbf,nel,nvir,F,tei,T1,T2) result(e)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: F(2*nbf,2*nbf)
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T1(nvir,nel)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp) :: e
   integer  :: i,j,a,b
   integer  :: ia,ja,ib,jb
   integer  :: iajb,ibja
   e = 0.0_wp
   do i = 1, nel
      do a = 1, nvir
         e = e + F(i,a+nel)*T1(a,i)
         do j = 1, nel
            do b = 1, nvir
               e = e + 0.25_wp * tei(a+nel,b+nel,i,j) * T2(a,b,i,j)
               e = e + 0.50_wp * tei(a+nel,b+nel,i,j) * T1(a,i)*T1(b,j)
            enddo
         enddo
      enddo
   enddo
end function ccsd_energy

!  ------------------------------------------------------------------------
!  CCD energy expression
!  ------------------------------------------------------------------------
function ccd_energy(nbf,nel,nvir,tei,T2) result(e)
   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nel
   integer, intent(in)  :: nvir
   real(wp),intent(in)  :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   real(wp),intent(in)  :: T2(nvir,nvir,nel,nel)
   real(wp) :: e
   integer  :: i,j,a,b
   integer  :: ia,ja,ib,jb
   integer  :: iajb,ibja
   e = 0.0_wp
   do i = 1, nel
      do j = 1, nel
         do a = 1, nvir
            do b = 1, nvir
               e = e + 0.25_wp * tei(a+nel,b+nel,i,j) * T2(a,b,i,j)
            enddo
         enddo
      enddo
   enddo
end function ccd_energy

end module cc
