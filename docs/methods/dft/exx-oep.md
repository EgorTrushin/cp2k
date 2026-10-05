# Exact-exchange optimized effective potential

The experimental `DFT / XC / HF / OEP` section enables restricted EXX-OEP for
isolated, closed-shell molecules. The exchange energy and the right-hand side of
the OEP equation use the nonlocal Fock exchange operator. The SCF Hamiltonian uses
the local, multiplicative exchange potential obtained from the OEP equation.

The implementation follows the auxiliary-basis preprocessing in Trushin and
Görling, *J. Chem. Phys.* **155**, 054109 (2021),
<https://doi.org/10.1063/5.0056431>, and the molecular `EXXOEP` implementation in
PyOEP. The auxiliary functions describe charge distributions; their Coulomb
potentials form the potential basis.

## Supported setup

- `RUN_TYPE ENERGY`, restricted spin, integer occupations, and an isolated cell
  with `PERIODIC NONE` in both `CELL` and `POISSON`.
- Standard SCF diagonalization and `ADDED_MOS -1`; no orbital level shift.
- One `HF` section with `FRACTION 1` and `POTENTIAL_TYPE COULOMB`.
- `XC_FUNCTIONAL NONE` for exchange-only calculations.
- An independent `BASIS_SET OEP` for every atomic kind, in addition to the orbital
  basis. Standard CP2K basis-file syntax is used for both.
- Full four-center exchange; ADMM and ACE are unsupported.

The initial implementation stores dense replicated auxiliary tensors and targets
small molecules. Forces, stress, periodic OEP, smearing, and unrestricted spin are
outside its current scope. GAPW permits all-electron comparisons with molecular
reference implementations. Converge the GAPW grids and the orbital and OEP bases
independently. Standard XC-potential cube output does not yet include the OEP
contribution; use the self-contained potential export described below to plot it.
Response properties have not been validated.

## Input fragment

```text
&SCF
  ADDED_MOS -1
  EPS_SCF 1.0E-9
  EPS_DIIS 10.0
  MAX_SCF 100
  &DIAGONALIZATION
    ALGORITHM STANDARD
  &END
  &MIXING
    METHOD DIRECT_P_MIXING
    ALPHA 1.0
  &END
&END
&XC
  &XC_FUNCTIONAL NONE
  &END
  &HF
    FRACTION 1.0
    &INTERACTION_POTENTIAL
      POTENTIAL_TYPE COULOMB
    &END
    &SCREENING
      EPS_SCHWARZ 1.0E-12
    &END
    &OEP
      THR_FAI_OEP 0.05
      HOMO_CONDITION FALSE
      POTENTIAL_FILE_NAME vx.oep.json
    &END
  &END
&END
```

For each `KIND`, supply the two bases, for example:

```text
BASIS_SET ORB my-orbital-basis
BASIS_SET OEP my-oep-basis
POTENTIAL ALL
```

The `OEP` section's presence activates the method. An initial guess without
orbital energies first undergoes one conventional HFX diagonalization to provide
the eigensystem needed by the OEP equation.
For an atomic initial guess, enabling DIIS early (as above) avoids undamped
fixed-point iterations before CP2K's usual DIIS threshold is reached. The
auxiliary preprocessing stabilizes the OEP linear solve; the outer SCF still
needs an appropriate mixing or DIIS setup.

For reproducible comparisons, start from a **converged HF wavefunction**, as in
PyOEP: run the same system without `OEP`, then use `SCF_GUESS RESTART` and
`DFT / WFN_RESTART_FILE_NAME` to select that HF restart file. In the constrained
method the retained auxiliary rank can depend on the starting orbitals. For
example, the HOMO-constrained water test reaches distinct 35- and 36-dimensional
solutions from an atomic and a converged-HF start, respectively; PyOEP also
retains both solutions when initialized accordingly. Preprocessing alone does
not guarantee a unique SCF solution.

## Preprocessing and diagnostics

1. Normalize the auxiliary Coulomb metric and orthogonalize it. All eigenvectors
   are retained; a nonpositive eigenvalue is an error, as this implementation
   requires a positive-definite auxiliary metric.
2. Construct a Fermi–Amaldi reference with auxiliary charge -1 and restrict the
   remainder to zero charge. Auxiliary charges are integrated analytically;
   PyOEP obtains the same quantities by numerical quadrature.
3. Optionally impose the HOMO expectation-value condition.
4. Diagonalize the **unweighted** occupied–virtual coupling Gram matrix. Retain
   eigenvalues at least `THR_FAI_OEP`, as in PyOEP, then construct and solve the
   energy-denominator-weighted response equation in that space.

Output reports the retained auxiliary dimension, the linear-equation residual,
and the charge-constraint error. The linear-equation residual measures the inner
OEP solve; SCF convergence must still be established separately. Total energies,
orbital energies, and potential shapes should be checked against basis and
threshold changes.

`DEBUG_FILE_NAME` optionally writes a binary stream for independent numerical
validation, replacing it at each OEP update. It contains five native 32-bit
integers (`nao`, `naux`, `nmo`, `nocc`, `nret`), followed by native double-precision,
Fortran-order arrays: three-center integrals `(nao,nao,naux)`, Coulomb metric,
auxiliary charges, MO coefficients, MO energies, signed nonlocal exchange matrix,
total potential coefficients, local exchange matrix, and retained transformation
`W(naux,nret)`. This is a diagnostic format, not a wavefunction restart format.

## Exporting and plotting the exchange potential

Set `DFT / XC / HF / OEP / POTENTIAL_FILE_NAME vx.oep.json` to export the local
exchange potential. The default is an empty filename, which disables export.
The exact filename is used, relative to the calculation's working directory;
existing contents are replaced at each OEP update. Thus the file contains the
last potential constructed, including if SCF subsequently fails to converge.
Check the SCF convergence message before treating it as a converged result.
Choose a different filename for each calculation and for `DEBUG_FILE_NAME`.

The versioned JSON file is self-contained: it holds the atomic positions and
explicit Gaussian expansions for the total and reference exchange potentials.
The reference is the Fermi–Amaldi reference for the default constraint, or its
HOMO-constrained variant when `HOMO_CONDITION TRUE`. The remainder is their
difference. No orbitals, response matrices, diagnostic dump, external basis files,
or CP2K installation are needed to evaluate the exported potential.

From the CP2K source directory, use `tools/oep/potential.py` with a Python
environment containing NumPy and SciPy. Matplotlib is needed only for `--plot`.
PySCF and PyOEP are not dependencies of this postprocessor.

```shell
# Sample a line, writing total/reference/remainder values and a plot.
python tools/oep/potential.py /path/to/vx.oep.json line \
  --start 0 0 -8 --end 0 0 8 --points 1601 \
  --output vx-line.csv --plot vx-line.pdf

# Write a three-dimensional total-exchange potential cube.
python tools/oep/potential.py /path/to/vx.oep.json cube \
  --spacing 0.2 --margin 5 --output vx.cube

# Evaluate any set of points, including a molecular plane.
python tools/oep/potential.py /path/to/vx.oep.json points coordinates.txt \
  --output vx-points.csv
```

All coordinates, margins, and spacings are in **bohr**, and all potential values
are in **hartree**. `coordinates.txt` has three whitespace-separated columns
`x y z`. Cube export accepts `--component total`, `reference`, or `remainder`.
The cube grid covers the molecular bounding box plus the requested margin;
its spacing is at most the requested value. Evaluation is analytic in the
Gaussian representation, including at nuclei, and does not use CP2K's real-space
grid. The grid selected here only controls sampling and can be changed without
rerunning SCF.

The Python API `OEPPotential(filename).evaluate(points)` accepts an `(npoints,3)`
array and returns three columns: total, reference, and remainder. Line CSV files
prepend `x_bohr,y_bohr,z_bohr` to these three columns. The reader checks the file
version, units, and integrated auxiliary charges. Independent evaluator and
cube-format tests can be run with:

```shell
python -m unittest discover -s tools/oep -p 'test_*.py' -v
```

### Portable format, version 1

The file identifies itself with `"format": "CP2K_OEP_POTENTIAL"` and
`"version": 1`. Each atomic center has an `element`, `atomic_number`, `position`
in bohr, and `terms`. Each term is an array
`[lx, ly, lz, exponent, total, reference]`, with column names also stored in
`term_columns`. It denotes the unnormalized charge function

$$
g(\mathbf r')=(x'-X)^{l_x}(y'-Y)^{l_y}(z'-Z)^{l_z}
\exp[-\alpha|\mathbf r'-\mathbf R|^2].
$$

The potential is the sum of **Coulomb integrals** of these functions, weighted by
the selected coefficient:

$$
v_x(\mathbf r)=\sum_t c_t\int\frac{g_t(\mathbf r')}{|\mathbf r-\mathbf r'|}\,d\mathbf r'.
$$

Evaluating `g` itself gives the auxiliary charge distribution, not the potential.
CP2K folds its primitive normalization, contractions, spherical transformations,
and final OEP coefficients into the exported Cartesian weights. This preserves
the potential for general contracted and mixed-angular-momentum auxiliary bases
without requiring another program to reproduce CP2K's basis conventions. Both
coefficient columns use the same explicit primitive basis; exactly zero terms
in both columns are omitted. The original auxiliary dimension, retained
dimension, threshold, HOMO setting, and integrated charges are also recorded.
The total and reference charges are -1; the remainder has zero net charge.

The existing `.wfn` restart does not save this export. If a calculation was run
without `POTENTIAL_FILE_NAME`, restart it with the option enabled to reconstruct
and save the OEP potential. The JSON is for postprocessing, not for restarting
SCF.

## Independent comparison with PyOEP

Two scripts in `tools/regtesting` generate matched molecular calculations and
compare a CP2K matrix dump with PySCF integrals and the original PyOEP solver.
Use a Python environment containing PySCF, NumPy, and SciPy, and pass the path to
a PyOEP checkout. These packages are used only for validation.

From the CP2K source directory:

```shell
python tools/regtesting/prepare_oep_pyoep.py \
  --pyoep /path/to/PyOEP --output /path/to/oep-tests --cases h2 h2o co
```

This exports the same Gaussian exponents and contractions and coordinates in
bohr to CP2K, sorting primitive rows by decreasing exponent as required by the
GAPW projector construction. It writes `hf.inp`, `oep.inp`, and `oep_homo.inp`,
enables both the diagnostic dump and JSON potential export, and converges
independent PyOEP references. Run `hf.inp` first in each molecular directory;
both OEP inputs explicitly read the resulting `hf-RESTART.wfn`. Then compare:

```shell
python tools/regtesting/compare_oep_pyoep.py /path/to/oep-tests/h2o \
  --pyoep /path/to/PyOEP --potential /path/to/oep-tests/h2o/oep.json
python tools/regtesting/compare_oep_pyoep.py /path/to/oep-tests/h2o \
  --pyoep /path/to/PyOEP --homo --potential /path/to/oep-tests/h2o/oep_homo.json
```

The comparison checks two- and three-center Coulomb integrals, the charge
constraint, Coulomb orthonormality, retained dimension, and the local exchange
matrix at identical orbitals. It also reports nonlocal-exchange differences and
the PySCF energy evaluated with the CP2K orbitals. Agreement of the latter with
the converged PyOEP reference tests the complete SCF solution separately from
GAPW total-energy integration errors.
The local matrix is compared both with unmodified PyOEP charges and with
independent analytic charges from PySCF's zero-frequency Fourier integrals.
The second comparison isolates the ported linear algebra from PyOEP's charge
quadrature error.
With `--potential`, it additionally checks the exported total and reference
potentials at nuclei, arbitrary off-axis points, and distant points against
PyOEP's independent PySCF-integral grid evaluation.

## Initial validation

The restricted implementation was built with GNU Fortran 13.3 and tested with
the serial/OpenMP executable. `QS/regtest-oep` passes all four energy checks with
two OpenMP threads; three existing `QS/regtest-hfx` cases retain their original
reference tolerances. Spin-polarized input and incomplete virtual spaces are
rejected with OEP diagnostics. The MPI executable also builds, but MPI execution
was unavailable in the validation environment.

Potential exports were independently checked for H2, H2O with both constraints,
CO, and a contracted H2 auxiliary basis sharing primitives across s, p, and d
functions. Total and reference potentials agree with PyOEP/PySCF grid evaluation
within 4e-12 Ha at the tested points. The molecular JSON files are approximately
6-18 kB. Five postprocessor tests cover analytic monopoles, independent quadrature
for Cartesian multipoles through degree six, component charges, cube coordinates
and values, and rejection of incompatible or damaged files.

For the default charge-only constraint, independent comparisons gave:

| Molecule | Orbital basis | Retained dimension | Maximum local-matrix difference with analytic charges (Ha) | CP2K total energy minus PyOEP (Ha) |
| --- | --- | ---: | ---: | ---: |
| H2 | cc-pVTZ | 8 | 1.7e-14 | +1.06e-6 |
| H2O | cc-pVTZ | 36 | 8.5e-13 | -2.09e-6 |
| CO | aug-cc-pwCVTZ | 64 | 3.9e-12 | -5.44e-5 |

All use aug-cc-pVDZ-RIFIT auxiliary functions and `THR_FAI_OEP 0.05`. These are
finite-grid GAPW results (`CUTOFF 600`, `REL_CUTOFF 80`, radial/Lebedev grids
200/590), with an 8-Angstrom cell for H2 and 18-Angstrom cells for H2O and CO.
The matrix comparison holds the orbitals and auxiliary charges fixed and tests
the ported OEP equation. Total-energy agreement additionally requires convergence
of CP2K's GAPW representation; ordinary HF already differs from PySCF by about
5.4e-5 Ha for the CO setup. Do not interpret the matrix precision as the accuracy
of the complete molecular calculation.

The HOMO-constrained variants also converge from the HF start and pass the same
matrix comparisons: retained dimensions 8, 36, and 63 for H2, H2O, and CO,
respectively, with maximum analytic-charge matrix differences below 3.4e-12 Ha.
