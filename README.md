# simple Hartree–Fock

A dead-simple quantum chemistry code, I spend my leisure time on (in case I have some).
This code is not intended to do calculations on real systems, but to playing around
with some basic concepts well known from lectures and literature.

## Usage
```s-hf [options] <input> [options]```

The inputfile must be provided as:
```
1> <nat> <nel> <nbf>
2> <x> <y> <z> <iat> <ibf>
3> <zeta> [<ityp> <ng>]...
4> ...
```
Start orbitals may be provided in a restart block

## Options
Set the convergence thresholds with `-a`, it may take an integer 
from 4 (lowest precision) to 9 (highest precision currently available).

```
$ s-hf h2.in -a 5
$ s-hf h2.in --acc 9
```

Set the number of iterations, if the SCF is not converged in the given
number of iterations, an error is thrown.

```
$ s-hf lih.in --iter 151
```

Does a simple steepest decent geometry optimisation:

```
$ s-hf hhe.in -o
$ s-hf lih.in --opt
```

Use Pulay's DIIS scheme to accelerate SCF convergence.
The number of errormatrices used for the DIIS can be specified, default is 6.

```
$ s-hf be.in --diis
$ s-hf lih.in --diis 10
```

Enforces an unrestricted Hartree-Fock. At the moment you cannot break the
symmetry in na-nb = 0 cases, so you will get the result of the restricted HF-SCF

```
$ s-hf li.in --uhf 1
$ s-hf h2.in --uhf 2
```

Does a second order Møller-Plesset many-body pertubation theory calculation,
employing a N5 integral transformation.

```
$ s-hf be.in --mp2
```

## References

This work is mainly inspired by D. Crawford’s programming project.
See:

http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming

and of course:

A. Szabo, N. Ostlund, Modern Quantum Chemistry: Introduction to
Advanced Electronic Structure Theory, Dover Publications, 1989.
