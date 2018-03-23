subroutine help
   print'(''Usage: s-hf [options] <input> [options]'',/)'

   print'(''The inputfile must be provided as:'','// &
   &      '/,4x,''1> <nat> <nel> <nbf>'','// &
   &      '/,4x,''2> <x> <y> <z> <iat> <ibf>'','// &
   &      '/,4x,''3> <zeta> [<ityp> <ng>]...'','// &
   &      '/,4x,''4> ...'')'
   print'(''Start orbitals may be provided in a restart block'',/)'

   print'(''Options:'')'

   print'(3x,''-a, --acc <int>   '','// &
   &      'x,''convergence threshold'')'

   print'(3x,''    --iter <int>  '','// &
   &      'x,''maximal number of iterations'')'

   print'(3x,''-o, --opt         '','// &
   &      'x,''do a steepest decent geometry optimization'')'

   print'(3x,''    --diis [int]  '','// &
   &      'x,''use a convergence accelerator in the SCF'')'

   print'(3x,''    --uhf <na-nb> '','// &
   &      'x,''do unrestricted Hartree-Fock'')'

   print'(3x,''    --mp2         '','// &
   &      'x,''do second order Møller-Plesset MBPT'')'

!  print'(3x,''    --scs-mp2     '','// &
!  &      'x,''use the spin compound scaled MP2 variant'')'

!  print'(3x,''    --ccsd        '','// &
!  &      'x,''do coupled cluster, with single and double excitation'')'

!  print'(3x,''    --ccsd(t)     '','// &
!  &      'x,''include also triple excitation pertubatively'')'

   print'(3x,''-v, --version     '','// &
   &      'x,''prints the version number and some literature'')'

   print'(3x,''-h, --help        '','// &
   &      'x,''Show this message'',/)'
end subroutine help

subroutine banner
   print'(''=========================================='')'
   print'('' simple HF by SAW      version 0.0 (1803)'')'
   print'(''=========================================='')'
end subroutine banner

subroutine credits
   call banner
   print'(a)'

   print'(a)',&
   '------------------------------------------------------------------------'
   print'(x,''This work is mainly inspired by'','// &
   &     'x,''D. Crawford’s programming project.'')'
   print'(x,''see:'','//&
   & 'x,''http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming'')'

   print'(a)',&
   '------------------------------------------------------------------------'
   print'(x,''and of course:'','// &
   &   '/,x,''A. Szabo, N. Ostlund, Modern Quantum Chemistry:'','//&
   &     'x,''Introduction to'','//&
   &   '/,x,''Advanced Electronic Structure Theory,'','//&
   &     'x,''Dover Publications, 1989.'')'

   print'(a)',&
   '------------------------------------------------------------------------'

end subroutine credits
