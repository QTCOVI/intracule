# intracule

Analytical radial intracule densities I(R) for atoms and molecules from Gaussian-basis wavefunctions.

This repository contains the Fortran program used in

> E. Francisco, J. A. Quiñonero and Á. Martín Pendás,
> *Quadrature-Free Radial Intracule Densities: A Bessel-Harmonic Formulation for Gaussian Basis Sets*,
> J. Chem. Theory Comput. (submitted).

The program reads a wavefunction in AIMPAC `.wfn` format (single determinant or CI/CAS
expansion), builds the second-order density matrix, and computes the spherically averaged
intracule density I(R) together with its Coulomb, exchange-correlation and Dirac-Fock
components. The angular average over the orientation of **R** is done analytically:
the coupling exponential exp(2μ **R**·**T**) is expanded in modified spherical
Bessel functions and real spherical harmonics, and the expansion terminates exactly at
l = p<sub>x</sub> + p<sub>y</sub> + p<sub>z</sub>. The radial moments ∫₀<sup>∞</sup> Rⁿ I(R) dR (norm N(N−1), V<sub>ee</sub>, and their
Coulomb/XC analogues) are obtained in closed form through ₁F₁, by adaptive QUADPACK
integration, or on Euler–Maclaurin/Becke radial grids. Lebedev and spherical t-design angular
quadratures are kept as numerical references.

For heavy atoms the angular-momentum coupling of the analytic norm cancels by 12–15 orders of
magnitude (Kr₂, STO-3G), and binary64 loses the norm. The program therefore carries a binary128
(`real(16)`) path for the analytic norm pass (`src/binary128.f90`, `src/computeEcoef_q.f`; ported
from the binary128 chain of J. A. Quiñonero): Hermite expansion coefficients, Gaussian prefactors,
angular tables G<sub>lm</sub>, real spherical harmonics, the Jtab coupling and the radial integrals, which
use the terminating Kummer polynomial e⁻ˣ ₁F₁(a; b; x) = ₁F₁(1−m; b; −x), m = (n−l)/2 ≥ 1
(positive terms only, no library ₁F₁). By default (`BINARY128 AUTO`) every shell quartet is
first evaluated in binary64 with a running rounding-error estimate, and only the quartets whose
estimated relative error exceeds a tolerance are recomputed in binary128; `BINARY128 ALWAYS` uses
binary128 for every quartet and `BINARY128 NEVER` the binary64 path alone. The repulsion energies
(R⁻¹ moments, non-terminating ₁F₁) are always evaluated in binary64.

## Layout

| Path | Contents |
| --- | --- |
| `src/` | Fortran 77/90 sources and `Makefile` (main program `intracule.f`, driver `otroiculedrv.f90`) |
| `src/alltdsg.f` | tabulated spherical t-design grids (63 MB of data statements; see the note below) |
| `examples/` | Kr₂ RHF/STO-3G binary128 test (`kr2/`); Ne atom, RHF/cc-pVTZ: `make_ne_wfn.py` (PySCF script that writes `ne.wfn`), two inputs, reference outputs, PySCF reference energies |
| `tools/wfn_precision.py` | checks the orthonormality of the MOs in a `.wfn` file and the sum rules they imply (NumPy) |
| `LICENSE` | GNU General Public License v3 |

## Requirements and build

- gfortran (tested with GCC 11 and 15) with OpenMP
- GNU Scientific Library (`libgsl`, `libgslcblas`; on Debian/Ubuntu `apt install libgsl-dev`)

```bash
cd src
make          # serial make; the module dependencies are not safe with make -j
```

> **Note on `alltdsg.f`.** It holds the spherical t-design angular grids as 63 MB of `data`
> statements. The `Makefile` compiles it with `-O0` (about 2 minutes and 5 GB of memory);
> with the general flags (`-g -fcheck=all -O2`) gfortran needs several hours for this single file.
> If you change `FFLAGS`, keep the separate `-O0` rule for `alltdsg.o`. The grids are used only by the
> `TDESIGN` angular quadrature (a numerical reference), but the file must be compiled to link.

`real(16)` arithmetic is done in software by gfortran (libquadmath) and is one to two orders of
magnitude slower than binary64; `BINARY128 AUTO` keeps the extra cost to the quartets that need it.

## Running

```bash
./intracule < input.inp > output.out
```

Two examples, Ne RHF/cc-pVTZ (`examples/`, 1–2 minutes each on 4 cores):

```
ne.wfn
threads 4
intracule
  bessel               ! analytical (Bessel) angular average
  quadpack             ! adaptive radial integration (omit for the default fixed grid)
  integrate analytic   ! closed-form radial moments
end
```

| | `ne-ana.out` (fixed grid) | `ne-quadpack.out` (QUADPACK) | Exact / PySCF |
| --- | --- | --- | --- |
| Norm, N(N−1) | 89.99999999 | 89.99999995 | 90 |
| Coulomb pairs, N² | 100.00000002 | 99.99999994 | 100 |
| XC pairs, abs. value (= N) | 10.00000003 | 9.99999999 | 10 |
| V<sub>ee</sub> | 54.06591036 | 54.06698053 | 54.06698054 |
| V<sub>ee</sub> Coulomb | 66.17838974 | 66.18053009 | 66.18053010 |
| V<sub>ee</sub> XC | −12.11247938 | −12.11354957 | −12.11354956 |

The XC pair count is printed with a positive sign; its magnitude is N, the −N sum rule of the
article. The fixed radial grid is good for the norms but leaves a 2×10⁻⁵ relative error in V<sub>ee</sub>;
QUADPACK reproduces the PySCF energies (`examples/reference_pyscf.txt`) to 10⁻⁸ hartree.
`make_ne_wfn.py` regenerates `ne.wfn`.

### Binary128: Kr₂ RHF/STO-3G (`examples/kr2/`)

The heavy-atom case of the article (Table 5, R = 8.06 bohr), run with the three values of the
`BINARY128` keyword. The inputs compute only the closed-form norms: `integrate numerical` with
`nrad 2` and `lebedev 6` makes the numerical I(R) step trivial, and the analytic norm pass that
follows it is the quantity of interest (ignore the 2-point numerical sums printed after it).

| `BINARY128` | N(N−1) (5112) | N² (5184) | XC (72) | quartets in binary128 | wall time, 32 cores |
| --- | --- | --- | --- | --- | --- |
| `NEVER`  | 4969.26309225 | 5041.24418188 | 71.98108964 | 0 | 8.3 min |
| `AUTO` (default, tol 10⁻¹⁰) | 5112.00006641 | 5184.00006735 | 72.00000093 | 240 453 of 1 382 976 | 9.4 min |
| `ALWAYS` | 5112.00006637 | 5184.00006730 | 72.00000093 | all | 9.4 min |

Binary64 alone loses 2.8% of the norm. The remaining 1.3×10⁻⁸ relative deviation from 5112 is the
same in every binary128 run and in an independent implementation, and it is exactly the one implied
by the 8-digit MO coefficients of the `.wfn` file (see `tools/wfn_precision.py` below). The wall times include the repulsion-energy pass,
which is always done in binary64.

### Checking the input wavefunction (`tools/wfn_precision.py`)

The small deviations from the exact sum rules (1e-8 to 1e-9 relative) come from the `.wfn` file,
whose MO coefficients have eight significant digits, not from the integration.
`python tools/wfn_precision.py file.wfn` computes, in float64 and independently of `intracule`, the
electron count N′ = tr(DS) and the pair counts N′² and −½ tr(DSDS) implied by the stored orbitals
(closed-shell single determinants). For every test system they coincide with the analytic norms of
`intracule` to all printed digits, e.g. Kr₂ 5112.00006637, 5184.00006730, −72.00000093.

Running `intracule` without arguments, or the `usage` keyword, prints all keywords
(radial mapping `RMESH`, `QUADPACK`, `LEBEDEV`/`TDESIGN` grids, tolerances `ABSQUADERR`/`RELQUADERR`, ...).

## Reproducing the paper

| Paper item | System | Input |
| --- | --- | --- |
| Test | Ne RHF/cc-pVTZ | `examples/ne-ana.inp`, `examples/ne-quadpack.inp` |
| Table 5, analytic norm (one geometry) | Kr₂ RHF/STO-3G, R = 8.06 bohr | `examples/kr2/kr2-{never,auto,always}.inp` |
<!-- TODO(authors): add wfn files and inputs/outputs for Table 1 (all molecules), Table 2 and Figs. 2, 4
     (H2 FCI scan), Table 3 (C2H6), Table 4 and Fig. 5 (CO CCSD), Table 5 and Fig. 6 (Kr₂ scan),
     Fig. 1 (F3-), Fig. 3 (Ne2 CAS) -->

## Licence

Copyright (C) 2026 E. Francisco, J. A. Quiñonero and Á. Martín Pendás,
Departamento de Química Física y Analítica, Universidad de Oviedo.

This program is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version. It is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; see `LICENSE` for details.
