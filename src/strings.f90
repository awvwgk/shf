module strings

contains

elemental function findelement(elem) result(i)
   integer :: i
   character(len=2),intent(in) :: elem
   character(len=*),parameter :: pse= " "//    &
   &  "H "//                            "He"// &
   &  "LiBe"//                "B C N O F Ne"// &
   &  "NaMg"//                "AlSiP S ClAr"// &
   &  "K CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKr"// &
   &  "RbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI Xe"// &
   &  "CsBa"//                                 &
   &        "LaCePrNdPmSmEuGdTbDyHoErTmYbLu"// &
   &        "HfTaW ReOsIrPtAuHgTlPbBiPoAtRn"

   i = index(pse,elem)/2

end function findelement

elemental function selectelement(i) result(elem)
   integer,intent(in) :: i
   character(len=2) :: elem
   character(len=2),parameter :: pse(*)= (/ &
   &   "H ","He", &
   &   "Li","Be","B ","C ","N ","O ","F ","Ne", &
   &   "Na","Mg","Al","Si","P ","S ","Cl","Ar", &
   &   "K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga", &
   &   "Ge","As","Se","Br","Kr", &
   &   "Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In", &
   &   "Sn","Sb","Te","I ","Xe", &
   &   "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho", &
   &   "Er","Tm","Yb","Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg", &
   &   "Tl","Pb","Bi","Po","At","Rn" /)

   elem = pse(i)
      
end function selectelement

elemental function lowercase (str) result(final)
   implicit none
   integer :: i,il
   character(len=*),intent(in) :: str
   character(len=len(str))     :: final
   character(len=26),parameter :: cap = &
   &  'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(len=26),parameter :: low = &
   &  'abcdefghijklmnopqrstuvwxyz'
   final = str
   do i = 1, len_trim(str)
      il = INDEX(cap, str(i:i))
      if (il > 0) final(i:i) = low(il:il)
   enddo
end function lowercase

elemental function uppercase (str) result(final)
   implicit none
   integer :: i,il
   character(len=*),intent(in) :: str
   character(len=len(str))     :: final
   character(len=26),parameter :: cap = &
   &  'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(len=26),parameter :: low = &
   &  'abcdefghijklmnopqrstuvwxyz'
   final = str
   do i = 1, len_trim(str)
      il = INDEX(low, str(i:i))
      if (il.gt.0) final(i:i) = cap(il:il)
   enddo
end function uppercase

elemental function capitalize (str) result(final)
   implicit none
   integer :: i,il
   character(len=*),intent(in) :: str
   character(len=len(str))     :: final
   character(len=26),parameter :: cap = &
   &  'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(len=26),parameter :: low = &
   &  'abcdefghijklmnopqrstuvwxyz'
   final = str
   il = INDEX(low, str(1:1))
   if (il > 0) final(1:1) = cap(il:il)
   do i = 2, len_trim(str)
      il = INDEX(cap, str(i:i))
      if (il > 0) final(i:i) = low(il:il)
   enddo
end function capitalize

end module strings
