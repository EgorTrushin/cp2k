#!/usr/bin/env python3
"""Check native scRPA against PyOEP and independent potential finite differences.

Run prepare_oep_pyoep.py --rpa, then CP2K hf.inp, oep.inp and rpa.inp.
The exchange dump provides the orbitals; the RPA dump provides the new kernel.
"""

import argparse
import json
import os
import sys
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"

import numpy as np
import scipy.linalg
from pyscf import gto, scf
from pyscf.gto import ft_ao


def cp2k_order(mol):
    """CP2K's p functions are y,z,x; higher shells use ascending m."""
    result = []
    loc = mol.ao_loc_nr()
    for shell in range(mol.nbas):
        angular = mol.bas_angular(shell)
        for contraction in range(mol.bas_nctr(shell)):
            order = [1, 2, 0] if angular == 1 else list(range(2 * angular + 1))
            result.extend(
                loc[shell] + contraction * (2 * angular + 1) + np.array(order)
            )
    return np.array(result)


def read_array(stream, shape):
    return np.fromfile(stream, dtype=np.float64, count=int(np.prod(shape))).reshape(
        shape, order="F"
    )


def error(a, b):
    return float(np.max(np.abs(a - b)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--pyoep", type=Path, required=True)
    parser.add_argument("--homo", action="store_true")
    parser.add_argument("--rpa-dump", type=Path)
    parser.add_argument("--exchange-dump", type=Path)
    parser.add_argument("--potential", type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.pyoep.resolve()))
    from methods.rpaoep import RPAOEP

    ref = json.loads((args.case / "reference.json").read_text())
    tag = "rpa_homo" if args.homo else "rpa"
    exchange = args.exchange_dump or args.case / (
        "oep_homo.bin" if args.homo else "oep.bin"
    )
    with exchange.open("rb") as stream:
        nao, naux, nmo, nocc, nx = np.fromfile(stream, dtype=np.int32, count=5)
        b = read_array(stream, (nao, nao, naux))
        metric = read_array(stream, (naux, naux))
        charge = read_array(stream, (naux,))
        c = read_array(stream, (nao, nmo))
        eps = read_array(stream, (nmo,))
        read_array(stream, (nao, nao))  # nonlocal exchange
        cx = read_array(stream, (naux,))
        vx = read_array(stream, (nao, nao))
        read_array(stream, (naux, nx))
        assert not stream.read()
    with (args.rpa_dump or args.case / f"{tag}.bin").open("rb") as stream:
        version, na, nb, nr, nm, no, nc, nrr, nw = np.fromfile(
            stream, dtype=np.int32, count=9
        )
        assert (version, na, nb, nm, no) == (1, nao, naux, nmo, nocc)
        x0, threshold, threshold_ri, ec = read_array(stream, (4,))
        br = read_array(stream, (nao, nao, nr))
        mr = read_array(stream, (nr, nr))
        qr = read_array(stream, (nr,))
        wc = read_array(stream, (naux, nc))
        wr = read_array(stream, (nr, nrr))
        rhs = read_array(stream, (nc,))
        cc = read_array(stream, (naux,))
        vc = read_array(stream, (nao, nao))
        assert not stream.read()
    mol = gto.M(atom=ref["geometry"], basis=ref["orbital_basis"], verbose=0)
    p = cp2k_order(mol)
    mf = scf.RHF(mol)
    mf.mo_coeff = np.empty_like(c)
    mf.mo_coeff[p] = c
    mf.mo_energy = eps
    mf.mo_occ = np.r_[np.full(nocc, 2.0), np.zeros(nmo - nocc)]
    oep = RPAOEP(mf, ref["oep_basis"], ref["ri_basis"], use_HOMO_condition=args.homo)
    q = cp2k_order(oep.auxmol)
    r = cp2k_order(oep.auxmol_ri)
    results = {"ranks": [int(nc), int(nrr)], "correlation_energy": float(ec)}
    results["ri_metric_error"] = error(mr, oep.auxmol_ri.intor("int2c2e")[np.ix_(r, r)])
    results["ri_integrals_error"] = error(
        br.transpose(2, 0, 1), oep.ints_3c_ao_ri[np.ix_(r, p, p)]
    )
    results["oep_integrals_error"] = error(
        b.transpose(2, 0, 1), oep.ints_3c_ao[np.ix_(q, p, p)]
    )
    results["correlation_charge"] = float(charge @ cc)
    results["oep_orthonormality"] = error(wc.T @ metric @ wc, np.eye(nc))
    results["ri_orthonormality"] = error(wr.T @ mr @ wr, np.eye(nrr))
    wc_py = np.empty_like(wc)
    wr_py = np.empty_like(wr)
    wc_py[q] = wc
    wr_py[r] = wr
    ints = mf.mo_coeff.T @ oep.ints_3c_ao @ mf.mo_coeff
    ints_ri = mf.mo_coeff.T @ oep.ints_3c_ao_ri @ mf.mo_coeff
    n = np.einsum("pk,pij->kij", wc_py, ints)
    m = np.einsum("pk,pij->kij", wr_py, ints_ri)
    rhs_py, ec_py = oep.get_scrpa_rhs(n, m, eps, nocc, nw, x0)
    results["kernel_rhs_error"] = error(rhs, rhs_py)
    results["kernel_energy_error"] = abs(float(ec - ec_py))
    # Recompute both retained spaces independently, with analytic auxiliary charges.
    oep.y = ft_ao.ft_ao(oep.auxmol, np.zeros((1, 3)))[0].real
    oep.yII = oep.WII.T @ oep.y
    oep.charge_norm = oep.yII @ oep.yII
    yri = ft_ao.ft_ao(oep.auxmol_ri, np.zeros((1, 3)))[0].real
    results["ri_charge_error"] = error(qr, yri[r])
    y2 = oep.WII_ri.T @ yri
    projector = np.eye(nr) - np.outer(y2, y2) / (y2 @ y2)
    oep.W3_ri = scipy.linalg.eigh(projector)[1][:, 1:]
    wcp = oep.get_W(oep.get_W3_charge(), ints, nocc, threshold)
    wrp = oep._get_W_ri(ints_ri, nocc, threshold_ri)
    assert (wcp.shape[1], wrp.shape[1]) == (nc, nrr)
    npy = np.einsum("pk,pij->kij", wcp, ints)
    mpy = np.einsum("pk,pij->kij", wrp, ints_ri)
    rhs_independent, ec_independent = oep.get_scrpa_rhs(npy, mpy, eps, nocc, nw, x0)
    cp = wcp @ scipy.linalg.solve(
        4 * oep.get_X0(ints, eps, nocc, wcp), -rhs_independent
    )
    vp = np.einsum("p,pij->ij", cp, oep.ints_3c_ao)
    results["potential_error"] = error(vc, vp[np.ix_(p, p)])
    results["energy_error"] = abs(float(ec - ec_independent))
    # An independent central finite difference varies a local potential, rediagonalizes
    # the Hamiltonian, and recomputes only the RPA energy in the FIXED response space.
    direction = np.random.default_rng(531).normal(size=nc)
    direction /= np.linalg.norm(direction)
    perturbation = np.einsum("k,kij->ij", direction, n)
    nodes, weights = np.polynomial.legendre.leggauss(nw)
    frequencies = x0 * (1 + nodes) / (1 - nodes)
    weights *= 2 * x0 / (1 - nodes) ** 2 / (2 * np.pi)

    def perturbed_energy(step):
        energies, rotation = scipy.linalg.eigh(np.diag(eps) + step * perturbation)
        mm = rotation.T @ m @ rotation
        ov = mm[:, nocc:, :nocc].reshape(nrr, -1)
        gaps = (energies[nocc:, None] - energies[None, :nocc]).reshape(-1)
        value = 0.0
        for omega, weight in zip(frequencies, weights):
            lam = 4 * gaps / (gaps**2 + omega**2)
            eigenvalues = scipy.linalg.eigvalsh((ov * lam) @ ov.T)
            value += weight * np.sum(np.log1p(eigenvalues) - eigenvalues)
        return value

    step = 2e-5
    derivative = (perturbed_energy(step) - perturbed_energy(-step)) / (2 * step)
    results["finite_difference_error"] = abs(float(derivative + direction @ rhs))
    results["pyscf_total_at_cp2k_orbitals"] = float(oep.e_tot + ec)
    results["total_vs_pyoep"] = float(
        oep.e_tot + ec - ref["rpa"][str(args.homo)]["energy"]
    )
    results["occupied_eigenvalues_vs_pyoep"] = error(
        eps[:nocc], ref["rpa"][str(args.homo)]["eigenvalues"][:nocc]
    )
    vxc = np.empty_like(vc)
    vxc[np.ix_(p, p)] = vx + vc
    fmo = mf.mo_coeff.T @ (mf.get_hcore() + oep.vj_ao + vxc) @ mf.mo_coeff
    results["occupied_virtual_residual_pyscf"] = error(fmo[nocc:, :nocc], 0)
    potential_file = args.potential or args.case / f"{tag}.json"
    if potential_file.exists():
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "oep"))
        from potential import OEPPotential
        from utils.coulomb_potential_on_grid import coulomb_potential_on_grid

        potential = OEPPotential(potential_file)
        grid = np.vstack(
            (
                mol.atom_coords(),
                np.random.default_rng(21).uniform(-6, 6, (100, 3)),
                500 * np.eye(3),
            )
        )
        integrals = coulomb_potential_on_grid(oep.auxmol, grid)[:, q]
        values = potential.evaluate(grid)
        results["export_exchange_error"] = error(values[:, 3], integrals @ cx)
        results["export_correlation_error"] = error(values[:, 4], integrals @ cc)
        results["export_total_error"] = error(values[:, 0], integrals @ (cx + cc))
        np.testing.assert_allclose(potential.charges, [-1, -1, -1, 0], atol=1e-9)
    print(json.dumps(results, indent=2))
    # Total SCF energies retain finite GAPW/grid errors; the fixed-orbital kernel
    # and integral checks below isolate correctness of the RPA port itself.
    for name, value in results.items():
        if name.endswith(("error", "orthonormality")) or name == "correlation_charge":
            assert abs(value) < 2e-8, (name, value)


if __name__ == "__main__":
    main()
