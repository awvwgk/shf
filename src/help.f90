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

   print'(3x,''-i, --iter <int>  '','// &
   &      'x,''maximal number of iterations'')'

   print'(3x,''-d, --direct      '','// &
   &      'x,''use a direct SCF approach'')'

   print'(3x,''    --dens <file> '','// &
   &      'x,''evaluate density on mesh specified in <file>'')'

   print'(3x,''-o, --opt         '','// &
   &      'x,''do a steepest decent geometry optimization'')'

   print'(3x,''    --diis [int]  '','// &
   &      'x,''use a convergence accelerator in the SCF'')'

   print'(3x,''-u, --uhf <na-nb> '','// &
   &      'x,''do unrestricted Hartree-Fock'')'

   print'(3x,''    --mp2         '','// &
   &      'x,''do second order Møller-Plesset MBPT'')'

!  print'(3x,''    --scs-mp2     '','// &
!  &      'x,''use the spin compound scaled MP2 variant'')'

   print'(3x,''    --ccd         '','// &
   &      'x,''do coupled cluster with double excitation (untested)'')'

!  print'(3x,''    --ccsd(t)     '','// &
!  &      'x,''include also triple excitation pertubatively'')'

   print'(3x,''    --tdhf        '','// &
   &      'x,''calculate excited states with TD-HF'')'

   print'(3x,''    --cis         '','// &
   &      'x,''calculate excited states with RPA'')'

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
   call gnu_gpl
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

subroutine gnu_gpl
   print'(a)'
   print'(" © 2017,2018 by Sebastian Ehlert")'
   print'(a)'
   print'(3x"This program is free software: you can redistribute it and/or")'
   print'(3x"modify it under the terms of the GNU General Public License as")'
   print'(3x"published by the Free Software Foundation, either version 3 of")'
   print'(3x"the License, or (at your option) any later version.")'
   print'(a)'
   print'(3x"This program is distributed in the hope that it will be useful,")'
   print'(3x"but WITHOUT ANY WARRANTY; without even the implied warranty of")'
   print'(3x"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.")'
   print'(a)'
   print'(3x"See the GNU General Public License for more details.")'
end subroutine
