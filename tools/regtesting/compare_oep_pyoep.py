#!/usr/bin/env python3
"""Compare CP2K's OEP dump with independent PySCF integrals and the original PyOEP solver."""

import os

os.environ["OMP_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import scipy.linalg
from pyscf import gto, scf
from pyscf.gto import ft_ao

parser = argparse.ArgumentParser()
parser.add_argument("case", type=Path)
parser.add_argument("--homo", action="store_true")
parser.add_argument("--pyoep", required=True, type=Path)
parser.add_argument("--dump", type=Path)
parser.add_argument("--threshold", type=float, default=0.05)
args = parser.parse_args()
sys.path.insert(0, str(args.pyoep.resolve()))
from methods.exxoep import EXXOEP

ref = json.loads((args.case / "reference.json").read_text())
dump = args.dump or args.case / ("oep_homo.bin" if args.homo else "oep.bin")
with dump.open("rb") as f:
    nao, naux, nmo, nocc, rank = np.fromfile(f, dtype=np.int32, count=5)

    def read(shape):
        return np.fromfile(f, dtype=np.float64, count=int(np.prod(shape))).reshape(
            shape, order="F"
        )

    b = read((nao, nao, naux))
    metric = read((naux, naux))
    charge = read((naux,))
    c = read((nao, nmo))
    eps = read((nmo,))
    k = read((nao, nao))
    coefficients = read((naux,))
    v = read((nao, nao))
    w = read((naux, rank))
    assert f.read() == b""


def cp2k_order(mol):
    result = []
    loc = mol.ao_loc_nr()
    for s in range(mol.nbas):
        l = mol.bas_angular(s)
        for ctr in range(mol.bas_nctr(s)):
            order = [1, 2, 0] if l == 1 else list(range(2 * l + 1))
            result.extend(loc[s] + ctr * (2 * l + 1) + np.array(order))
    return np.array(result)


mol = gto.M(atom=ref["geometry"], basis=ref["orbital_basis"], verbose=0)
p = cp2k_order(mol)
mf = scf.RHF(mol)
mf.mo_coeff = np.empty_like(c)
mf.mo_coeff[p, :] = c
mf.mo_energy = eps
mf.mo_occ = np.r_[np.full(nocc, 2.0), np.zeros(nmo - nocc)]
oep = EXXOEP(mf, ref["oep_basis"], use_HOMO_condition=args.homo)
q = cp2k_order(oep.auxmol)


def err(a, b):
    return float(np.max(np.abs(a - b)))


results = {
    "nao": int(nao),
    "naux": int(naux),
    "retained": int(rank),
    "metric_error": err(metric, oep.auxmol.intor("int2c2e")[np.ix_(q, q)]),
    "three_center_error": err(b.transpose(2, 0, 1), oep.ints_3c_ao[np.ix_(q, p, p)]),
    "charge_error_vs_grid": err(charge, oep.y[q]),
    "exchange_matrix_error": err(k, oep.vxnl_ao[np.ix_(p, p)]),
    "charge_constraint": float(charge @ coefficients + 1),
    "coulomb_orthonormality": err(w.T @ metric @ w, np.eye(rank)),
    "pyscf_energy_at_cp2k_orbitals": float(oep.e_tot),
    "energy_vs_pyoep": float(oep.e_tot - ref[str(args.homo)]["energy"]),
    "occupied_eigenvalue_error": err(
        eps[:nocc], np.array(ref[str(args.homo)]["eigenvalues"])[:nocc]
    ),
}
ints = mf.mo_coeff.T @ oep.ints_3c_ao @ mf.mo_coeff
# Use CP2K's exchange matrix, so density-mixing errors can be distinguished
# from errors in the ported OEP solver.
kp = np.empty_like(k)
kp[np.ix_(p, p)] = k


def pyoep_matrix():
    if args.homo:
        z, z2 = oep.get_z_and_zII(ints, eps, nocc)
        vr = oep.get_v_ref_w_homo(z2, ints, mf.mo_coeff, nocc, kp)
        w3 = oep.get_W3_charge_and_homo(z2)
    else:
        vr = oep.get_v_ref(ints, nocc)
        w3 = oep.get_W3_charge()
    wp = oep.get_W(w3, ints, nocc, args.threshold)
    x = oep.get_X0(ints, eps, nocc, wp)
    vra = np.einsum("Pij,P->ij", oep.ints_3c_ao, vr)
    rhs = oep.get_rhs(vra, kp, ints, mf.mo_coeff, eps, nocc, wp)
    coeff = vr + wp @ scipy.linalg.solve(x, rhs)
    return np.einsum("Pij,P->ij", oep.ints_3c_ao, coeff), wp


vp, wp = pyoep_matrix()
results["pyoep_retained"] = int(wp.shape[1])
results["local_matrix_vs_pyoep"] = err(v, vp[np.ix_(p, p)])
fmo = mf.mo_coeff.T @ (mf.get_hcore() + oep.vj_ao + vp) @ mf.mo_coeff
results["occupied_virtual_residual_pyscf"] = float(np.max(np.abs(fmo[nocc:, :nocc])))
# A second comparison removes the quadrature error in PyOEP's auxiliary charges.
# The zero-frequency Fourier integral is an independent analytic evaluation.
analytic_charge = ft_ao.ft_ao(oep.auxmol, np.zeros((1, 3)))[0].real
results["charge_error_vs_analytic"] = err(charge, analytic_charge[q])
oep.y = analytic_charge
oep.yII = oep.WII.T @ oep.y
oep.charge_norm = oep.yII @ oep.yII
vp_analytic, _ = pyoep_matrix()
results["local_matrix_vs_pyoep_analytic_charges"] = err(v, vp_analytic[np.ix_(p, p)])
print(json.dumps(results, indent=2))
assert results["metric_error"] < 1e-8
assert results["three_center_error"] < 1e-8
assert results["local_matrix_vs_pyoep"] < 2e-7
assert results["pyoep_retained"] == rank
assert abs(results["charge_constraint"]) < 1e-9
assert results["coulomb_orthonormality"] < 1e-8
assert results["charge_error_vs_analytic"] < 1e-9
assert results["local_matrix_vs_pyoep_analytic_charges"] < 2e-9
