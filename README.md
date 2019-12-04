Hückel Solver
=============

Jacob Rowlands <jdr45@cam.ac.uk>

Overview
--

This program calculates and prints the Hückel energies and degeneracies of the π-system for a given molecule. Currently supported molecules are linear and cyclic poly-enes of arbitrary length, the platonic solids and buckminsterfullerene, C60. All energies are given numerically in units of β, with α = 0 and β = -1.

Requirements
--

Python3 with the networkx and numpy libraries.

Usage
--

```
./huckel.py [-l | --linear-polyene] num-sites
./huckel.py [-c | --cyclic-polyene] num-sites
./huckel.py [-p | --platonic] [4 | 6 | 8 | 12 | 20]
./huckel.py [-b | --buckyball]
```

To run unit tests:
```
./huckel-tests.py
```

Examples
--

* Calculate the Hückel MOs for benzene:
```
$ ./huckel.py -c 6
  ――     2.0

――  ――   1.0

――  ――  -1.0

  ――    -2.0

6 orbitals.
```

* Calculate the Hückel MOs for butadiene:
```
$ ./huckel.py -l 4
――   1.618

――   0.618

――  -0.618

――  -1.618

4 orbitals.
```