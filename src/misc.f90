module misc

contains

elemental function idx(i,j)
   implicit none
   integer,intent(in) :: i,j
   integer :: idx
   if (i.gt.j) then
      idx = i*(i-1)/2 + j
   else
      idx = j*(j-1)/2 + i
   endif
end function idx

end module misc
