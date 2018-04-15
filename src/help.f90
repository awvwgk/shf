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

   print'(3x,''    --ccd         '','// &
   &      'x,''do coupled cluster with double excitation (untested)'')'

!  print'(3x,''    --ccsd(t)     '','// &
!  &      'x,''include also triple excitation pertubatively'')'

   print'(3x,''-v, --version     '','// &
   &      'x,''prints the version number and some literature'')'

   print'(3x,''-h, --help        '','// &
   &      'x,''Show this message'',/)'
end subroutine help

subroutine banner
   use iso_fortran_env, only : output_unit
   write(output_unit,'(72(''=''))')
   write(output_unit,'(x,''simple HF by Sebastian Ehlert'')',advance='no')
   write(output_unit,'(23x,''version 1.0 (1804)'')')
   write(output_unit,'(72(''=''))')
end subroutine banner

subroutine credits
   call banner
   print'(a)'

   print'(72(''-''))'
   print'(x,''This work is mainly inspired by'','// &
   &     'x,''D. Crawford’s programming project.'')'
   print'(x,''see:'','//&
   & 'x,''http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming'')'

   print'(72(''-''))'
   print'(x,''and of course:'','// &
   &   '/,x,''A. Szabo, N. Ostlund, Modern Quantum Chemistry:'','//&
   &     'x,''Introduction to'','//&
   &   '/,x,''Advanced Electronic Structure Theory,'','//&
   &     'x,''Dover Publications, 1989.'')'

   print'(72(''-''))'
   print'(x,''integrals are drawn from this well known work:'','// &
   &   '/,x,''T. Helgaker, P. Jørgensen, J. Olsen, Molecular'','// &
   &     'x,''Electronic-Structure'','// &
   &   '/,x,''Theory, WILEY, 2000.'')'
   print'(72(''-''))'

end subroutine credits
