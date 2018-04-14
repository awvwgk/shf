module cc

   interface ccd
   module procedure solve_ccd
   end interface ccd

contains
subroutine chem2phys(nbf,eri,tei)
   use precision, only : wp => dp
   use misc,      only : idx
   implicit none
   integer, intent(in)  :: nbf
   real(wp),intent(in)  :: eri((nbf*(nbf+1)/2)*((nbf*(nbf+1)/2)+1)/2)
   real(wp),intent(out) :: tei(2*nbf,2*nbf,2*nbf,2*nbf)
   integer  :: i,j,k,l,ik,il,jk,jl,ikjl,iljk
   real(wp) :: eri1,eri2
   tei = 0.0_wp
   do i = 1, 2*nbf
      do j = 1, 2*nbf
         do k = 1, 2*nbf
            do l = 1, 2*nbf
               ik = idx(ceiling(i/2.0),ceiling(k/2.0))
               jl = idx(ceiling(j/2.0),ceiling(l/2.0))
               il = idx(ceiling(i/2.0),ceiling(l/2.0))
               jk = idx(ceiling(j/2.0),ceiling(k/2.0))
               ikjl = idx(ik,jl)
               iljk = idx(il,jk)
               if ((mod(i,2).eq.mod(k,2)).and.(mod(j,2).eq.mod(l,2))) then
                  eri1 = eri(ikjl)
               else
                  eri1 = 0.0_wp
               endif
               if ((mod(i,2).eq.mod(l,2)).and.(mod(j,2).eq.mod(l,2))) then
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
   integer  :: i,j,k
   real(wp) :: g
   do i = 1, 2*nbf
      do j = 1, 2*nbf
         g = 0.0_wp
         do k = 1, nel
            g = g + tei(j,k,i,k)
         enddo
         F(j,i) = H(ceiling(j/2.0),ceiling(i/2.0)) + g
      enddo
   enddo
end subroutine spinfockian
   
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
               D2(a,b,i,j) = 1.0_wp/(F(i,i)+F(j,j)-F(a+nel,a+nel)-F(b+nel,b+nel))
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
   write(id,'(15(''-''))')
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
      rmsd = sqrt( sum( (T2-TS)**2 ) )

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

pure subroutine build_T2(nbf,nel,nvir,F,tei,TS,Fac,Fki,Wklij,Wkbcj,D2,T2)
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

end subroutine build_T2

pure subroutine build_Fki(nbf,nel,nvir,F,tei,T2,Fki)
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

end subroutine build_Fki

pure subroutine build_Fac(nbf,nel,nvir,F,tei,T2,Fac)
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

end subroutine build_Fac

pure subroutine build_Wklij(nbf,nel,nvir,F,tei,T2,Wklij)
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

end subroutine build_Wklij

pure subroutine build_Wkbcj(nbf,nel,nvir,F,tei,T2,Wkbcj)
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

end subroutine build_Wkbcj

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
