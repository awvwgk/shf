module mbpt

contains

!* second order Møller-Plesset many-body pertubation theory
!  Ec = E(SS) + E(OS)
!     =  Σij Σab ( 2(ia|jb)-(ib|ja) ) · (ia|jb) / (εi+εj+εa+εb)
!  E(SS) = ½Σij εij + ½Σîĵ εîĵ
!        = Σij Σab ( (ia|jb)-(ib|ja) ) · (ia|jb) / (εi+εj+εa+εb)
!  E(OS) = Σîj εîj
!        = Σîj Σâb (îâ|jb)² / (εî+εj+εâ+εb)
!  εij = Σab (T(ij,ab) - T(ij,ba))·(ia|jb)
!  εîj = Σâb T(îj,âb)·(ia|ĵb)
!  T(ij,ab) = (ia|jb)/(εi+εj+εa+εb)
subroutine mp2(nbf,nocc,C,eri,eps,e,acc)
   use precision, only : wp => dp
   use ints,      only : teitrafo
   use misc,      only : idx
   implicit none
   integer, intent(in)    :: nbf
   integer, intent(in)    :: nocc
   real(wp),intent(in)    :: eps(nbf)
   real(wp),intent(in)    :: C(nbf,nbf)
   real(wp),intent(inout) :: eri(nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2)
   real(wp),intent(inout) :: e
   character(len=*),intent(in) :: acc
   integer  :: i,j,a,b,ia,jb,ib,ja,iajb,ibja
   real(wp) :: dum,eridum,ess,eos,ssc,osc,dmt

   eos=0.0_wp
   ess=0.0_wp

   print'(a)'
   print'(72(''-''))'
   print'('' second order Møller-Plesset calculation'')'

!  print'('' * doing Θ(N⁸) integral transformation'')'
!  call teitrafo_N8(nbf,eri,C)
   print'('' * doing Θ(N⁵) integral transformation'')'
   call teitrafo(nbf,eri,C)
   print'(''----[i]-----[j]-------------[pair energy]-''30(''-''))'
   do i = 1, nocc
      do j = 1, i
         osc = 0.0_wp
         ssc = 0.0_wp
         do a = nocc+1, nbf
            ia = idx(i,a)
            ja = idx(j,a)
            do b = nocc+1, nbf
               ib = idx(i,b)
               jb = idx(j,b)
               iajb = idx(ia,jb)
               ibja = idx(ib,ja)
               dmt = 1.0_wp / (eps(i)+eps(j)-eps(a)-eps(b))
               osc = osc + eri(iajb)*(eri(iajb)-eri(ibja)) * dmt
               ssc = ssc + eri(iajb)**2 * dmt
            enddo
         enddo
         if (i.ne.j) then
            ssc = 2*ssc
            osc = 2*osc
         endif
         print'(x,i5,3x,i5,x,f25.'//acc//')',i,j,osc+ssc
         ess = ess + ssc
         eos = eos + osc
      enddo
   enddo
   e = e + ess+eos
   print'(72(''-''))'
   print'('' same spin energy    :'',f18.'//acc//')',ess
   print'('' opposite spin energy:'',f18.'//acc//')',eos
   print'(72(''-''))'
   print'('' FINAL MP2 ENERGY'',f53.'//acc//')',e
   print'(72(''-''))'
end subroutine mp2 

end module mbpt
