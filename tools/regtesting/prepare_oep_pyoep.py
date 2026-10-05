#!/usr/bin/env python3
"""Generate matched CP2K inputs and PyOEP references for isolated EXX-OEP tests.

Requires PySCF and the supplied PyOEP checkout; CP2K itself does not use Python.
"""

import argparse
from pathlib import Path
import json
import os
import re
import sys

os.environ["OMP_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pyoep", required=True, type=Path)
parser.add_argument("--output", required=True, type=Path)
parser.add_argument(
    "--cases", nargs="+", choices=["h2", "co", "h2o"], default=["h2", "co", "h2o"]
)
parser.add_argument(
    "--rpa", action="store_true", help="Also generate self-consistent RPA tests"
)
parser.add_argument("--ri-basis", default="aug-cc-pVQZ-RIFIT")
parser.add_argument(
    "--cell-size", type=float, default=18.0, help="CP2K cubic cell side in angstrom"
)
args = parser.parse_args()
sys.path.insert(0, str(args.pyoep.resolve()))
from pyscf import gto, scf
from methods.exxoep import EXXOEP
from methods.rpaoep import RPAOEP

root = args.output.resolve()
root.mkdir(parents=True, exist_ok=True)
cases = {
    "h2": ("H 0 0 -0.37; H 0 0 0.37", "cc-pVTZ"),
    "co": ("C 0 0 -0.646514; O 0 0 0.484886", "aug-cc-pwCVTZ"),
    "h2o": ("O 0 0 0; H 0.758602 0 0.504284; H -0.758602 0 0.504284", "cc-pVTZ"),
}

for name in args.cases:
    geom, orbital = cases[name]
    aux = "aug-cc-pVDZ-RIFIT"
    run = root / name
    run.mkdir(exist_ok=True)
    mol = gto.M(atom=geom, basis=orbital, verbose=0)
    elements = sorted(set(mol.atom_symbol(i) for i in range(mol.natm)))
    basis = []
    for el in elements:
        sources = [("ORB", orbital), ("OEP", aux)]
        if args.rpa:
            sources.append(("RI", args.ri_basis))
        for tag, source in sources:
            shells = gto.basis.load(source, el)
            basis += [f"{el} {tag}-PYOEP", str(len(shells))]
            for shell in shells:
                l = shell[0]
                # GAPW uses the last primitive as the most diffuse one when
                # constructing projectors. PySCF basis files need not be sorted.
                rows = sorted(shell[1:], key=lambda row: row[0], reverse=True)
                basis.append(f"{l + 1} {l} {l} {len(rows)} {len(rows[0]) - 1}")
                basis += [" ".join(f"{v:.17g}" for v in row) for row in rows]
    (run / "BASIS").write_text("\n".join(basis) + "\n")
    kinds = "\n".join(
        f"""    &KIND {el}
      BASIS_SET ORB ORB-PYOEP
      BASIS_SET OEP OEP-PYOEP
      {"BASIS_SET RI_AUX RI-PYOEP" if args.rpa else ""}
      POTENTIAL ALL
      RADIAL_GRID 200
      LEBEDEV_GRID 590
    &END KIND"""
        for el in elements
    )
    # Avoid differences between CP2K's and PySCF's Angstrom-to-bohr constants.
    coords = "\n".join(
        f"      {mol.atom_symbol(i)} " + " ".join(f"{x:.17g}" for x in xyz)
        for i, xyz in enumerate(mol.atom_coords())
    )
    inp = f"""&GLOBAL
  PROJECT {name}
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS
    POTENTIAL_FILE_NAME ALL_POTENTIALS
    &QS
      METHOD GAPW
      EPS_DEFAULT 1e-12
      EPS_PGF_ORB 1e-16
      EPSFIT 1e-10
      EPSISO 1e-12
      EPSRHO0 1e-12
    &END QS
    &MGRID
      CUTOFF 600
      REL_CUTOFF 80
    &END MGRID
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER MT
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      ADDED_MOS -1
      EPS_SCF 1e-9
      EPS_DIIS 10.0
      MAX_SCF 100
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD DIRECT_P_MIXING
        ALPHA 1.0
      &END MIXING
    &END SCF
    &XC
      &XC_FUNCTIONAL NONE
      &END XC_FUNCTIONAL
      &HF
        FRACTION 1.0
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE COULOMB
        &END INTERACTION_POTENTIAL
        &SCREENING
          EPS_SCHWARZ 1e-12
          SCREEN_ON_INITIAL_P FALSE
        &END SCREENING
        &MEMORY
          MAX_MEMORY 1000
        &END MEMORY
        &OEP
          THR_FAI_OEP 0.05
          DEBUG_FILE_NAME oep.bin
          POTENTIAL_FILE_NAME oep.json
        &END OEP
      &END HF
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC {args.cell_size} {args.cell_size} {args.cell_size}
      PERIODIC NONE
    &END CELL
    &COORD
      UNIT BOHR
{coords}
    &END COORD
{kinds}
  &END SUBSYS
&END FORCE_EVAL
"""
    # PyOEP starts from a converged HF state. Use the same starting point in CP2K;
    # the rank selected by the OEP filter can otherwise lead to another solution.
    hf_input = re.sub(r"        &OEP\n.*?        &END OEP\n", "", inp, flags=re.S)
    hf_input = hf_input.replace(f"PROJECT {name}", "PROJECT hf").replace(
        "ALPHA 1.0", "ALPHA 0.4"
    )
    (run / "hf.inp").write_text(hf_input)
    inp = inp.replace("SCF_GUESS ATOMIC", "SCF_GUESS RESTART").replace(
        "    POTENTIAL_FILE_NAME ALL_POTENTIALS",
        "    POTENTIAL_FILE_NAME ALL_POTENTIALS\n    WFN_RESTART_FILE_NAME hf-RESTART.wfn",
    )
    (run / "oep.inp").write_text(inp)
    homo_inp = inp.replace(
        "THR_FAI_OEP 0.05", "THR_FAI_OEP 0.05\n          HOMO_CONDITION TRUE"
    )
    homo_inp = homo_inp.replace(
        "DEBUG_FILE_NAME oep.bin", "DEBUG_FILE_NAME oep_homo.bin"
    )
    homo_inp = homo_inp.replace(
        "POTENTIAL_FILE_NAME oep.json", "POTENTIAL_FILE_NAME oep_homo.json"
    )
    (run / "oep_homo.inp").write_text(homo_inp)
    mf = scf.RHF(mol).run(conv_tol=1e-12)
    references = {
        "geometry": geom,
        "orbital_basis": orbital,
        "oep_basis": aux,
        "hf": mf.e_tot,
    }
    for homo in [False, True]:
        oep = EXXOEP(mf, aux, use_HOMO_condition=homo)
        oep.run(maxit=100, thr_fai_oep=0.05, e_conv_thr=1e-12)
        references[str(homo)] = {
            "energy": oep.e_tot,
            "eigenvalues": oep.mf.mo_energy.tolist(),
            "converged": oep.converged,
        }
    assert all(references[str(h)]["converged"] for h in [False, True])
    if args.rpa:
        references["ri_basis"] = args.ri_basis
        references["rpa"] = {}
        for homo in [False, True]:
            tag = "rpa_homo" if homo else "rpa"
            source = homo_inp if homo else inp
            rpa_input = source.replace(f"PROJECT {name}", f"PROJECT {tag}")
            rpa_input = rpa_input.replace(
                "          THR_FAI_OEP",
                """          &RPA
            QUADRATURE_POINTS 20
            FREQUENCY_SCALE 2.5
            THR_FAI_RI 1e-8
            DEBUG_FILE_NAME rpa.bin
          &END RPA
          THR_FAI_OEP""",
            )
            rpa_input = rpa_input.replace(
                "DEBUG_FILE_NAME rpa.bin", f"DEBUG_FILE_NAME {tag}.bin"
            )
            rpa_input = rpa_input.replace(
                "POTENTIAL_FILE_NAME oep.json", f"POTENTIAL_FILE_NAME {tag}.json"
            )
            rpa_input = rpa_input.replace(
                "POTENTIAL_FILE_NAME oep_homo.json", f"POTENTIAL_FILE_NAME {tag}.json"
            )
            rpa_input = rpa_input.replace("hf-RESTART.wfn", f"{name}-RESTART.wfn")
            (run / f"{tag}.inp").write_text(rpa_input)
            start = EXXOEP(mf, aux, use_HOMO_condition=homo)
            start.run(maxit=100, e_conv_thr=1e-12)
            rpa = RPAOEP(start.mf, aux, args.ri_basis, use_HOMO_condition=homo)
            rpa.run(maxit=100, e_conv_thr=1e-12)
            assert rpa.converged
            references["rpa"][str(homo)] = {
                "energy": rpa.e_tot,
                "correlation_energy": rpa.E_corr,
                "eigenvalues": rpa.mf.mo_energy.tolist(),
                "converged": rpa.converged,
            }
    (run / "reference.json").write_text(json.dumps(references, indent=2))
    print(
        name,
        references["hf"],
        references["False"]["energy"],
        references["True"]["energy"],
        flush=True,
    )
