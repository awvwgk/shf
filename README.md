# simple Hartree–Fock

A dead-simple quantum chemistry code, I spend my leisure time on.
This code is not intended to do calculations on real systems,
but to playing around with some basic concepts well known from
lectures and literature.

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

Set the number of iterations with ``--iter``, if the SCF is not
converged in the given number of iterations, an error is thrown.

```
$ s-hf lih.in --iter 151
```

Direct SCF can be requested by `--direct`

For a simple steepest decent geometry optimization use

```
$ s-hf hhe.in -o
$ s-hf lih.in --opt
```

To use Pulay's DIIS scheme to accelerate SCF convergence use the `--diis` flag.
The number of errormatrices used for the DIIS can be specified, default is 6.

```
$ s-hf be.in --diis
$ s-hf lih.in --diis 10
```

To enforce an unrestricted Hartree-Fock set `--uhf <nalp-nbet>`.
At the moment the symmetry is broken for `--uhf 0` cases by using a hcore
guess for alpha MOs and a really bad orthonormalizer guess for beta MOs
DIIS seems to fail for this bad guess, so it will be deactivated in case of
`--uhf 0`.

```
$ s-hf li.in --uhf 1
$ s-hf h2.in --uhf 2
```

### post-HF methods

To use a second order Møller-Plesset many-body pertubation theory calculation,
use the ``--mp2`` flag. Currently only a O(N5) integral transformation is
in use.

```
$ s-hf be.in --mp2
```

For Coupled Cluster Doubles in spin orbitals

```
$ s-hf h2.in --ccd
```

can be employed.

### excited states calculation

A crude implementation of CIS and TD-HF is present

For CIS use `--cis` as

```
$ s-hf h2.in --cis
```

and for TD-HF use `--tdhf`

```
$ s-hf he.in --tdhf
```

### Further analysis

For a density analysis a mesh can specified as

```
1> <nmesh>
2> <x> <y> <z> <weight>
3> <x> <y> <z> <weight>
4> ...
```

The file can be read in by `--dens <mesh.in>`

## References

This work is mainly inspired by D. Crawford’s programming project.
See:

http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming

and of course:

A. Szabo, N. Ostlund, Modern Quantum Chemistry: Introduction to
Advanced Electronic Structure Theory, Dover Publications, 1989.

more recently I also implemented integrals drawn from this well known work:

T. Helgaker, P. Jørgensen, J. Olsen, Molecular Electronic-Structure
Theory, WILEY, 2000.
