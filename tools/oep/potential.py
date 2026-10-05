#!/usr/bin/env python3
"""Evaluate a self-contained CP2K OEP potential export without running SCF.

Requires NumPy and SciPy; Matplotlib is optional for line plots. All coordinates
are in bohr and potentials in hartree. No CP2K basis files or PySCF are needed.
"""

import argparse
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.special import gamma, hyp1f1


def hermite_coefficients(power, exponent):
    """Expand x**power exp(-a*x*x) in derivatives with respect to its center."""
    coefficients = np.array([1.0])
    for _ in range(power):
        next_coefficients = np.zeros(len(coefficients) + 1)
        for order, coefficient in enumerate(coefficients):
            next_coefficients[order + 1] += coefficient / (2 * exponent)
            if order:
                next_coefficients[order - 1] += order * coefficient
        coefficients = next_coefficients
    return coefficients


def gaussian_potentials(displacement, exponent, powers):
    """Coulomb integrals of unnormalized Cartesian Gaussians at given points.

    displacement is Gaussian center minus observation point. Hermite Coulomb
    recurrences act on (2*pi/a)*F_0(a*R**2), where
    F_n(t) = hyp1f1(n+1/2, n+3/2, -t)/(2*n+1).
    This expression also handles points exactly at a Gaussian center.
    """
    maximum = max(sum(power) for power in powers)
    orders = np.arange(maximum + 1)[:, None]
    argument = exponent * np.einsum("ij,ij->i", displacement, displacement)
    boys = hyp1f1(orders + 0.5, orders + 1.5, -argument) / (2 * orders + 1)
    base = (-2 * exponent) ** orders * boys

    @lru_cache(None)
    def coulomb(tx, ty, tz, order):
        indices = [tx, ty, tz]
        for axis, power in enumerate(indices):
            if power:
                indices[axis] -= 1
                value = displacement[:, axis] * coulomb(*indices, order + 1)
                if power > 1:
                    indices[axis] -= 1
                    value = value + (power - 1) * coulomb(*indices, order + 1)
                return value
        return base[order]

    expansions = {
        power: hermite_coefficients(power, exponent)
        for power in {p for triple in powers for p in triple}
    }
    result = np.empty((len(displacement), len(powers)))
    for column, (lx, ly, lz) in enumerate(powers):
        value = np.zeros(len(displacement))
        for tx in range(lx % 2, lx + 1, 2):
            for ty in range(ly % 2, ly + 1, 2):
                for tz in range(lz % 2, lz + 1, 2):
                    weight = (
                        expansions[lx][tx] * expansions[ly][ty] * expansions[lz][tz]
                    )
                    value += weight * coulomb(tx, ty, tz, 0)
        result[:, column] = (2 * np.pi / exponent) * value
    coulomb.cache_clear()
    return result


class OEPPotential:
    """Portable Gaussian expansion; evaluate() returns total, reference, remainder."""

    def __init__(self, filename):
        self.metadata = json.loads(Path(filename).read_text())
        data = self.metadata
        if data.get("format") != "CP2K_OEP_POTENTIAL" or data.get("version") != 1:
            raise ValueError("Unsupported CP2K OEP potential format or version")
        if data.get("length_unit") != "bohr" or data.get("potential_unit") != "hartree":
            raise ValueError("Expected bohr/hartree units")
        if data.get("representation") != "unnormalized_cartesian_gaussians":
            raise ValueError("Unsupported Gaussian representation")
        if data.get("term_columns") != [
            "lx",
            "ly",
            "lz",
            "exponent",
            "total",
            "reference",
        ]:
            raise ValueError("Unexpected Gaussian term columns")
        self.atoms = data["atoms"]
        if not self.atoms:
            raise ValueError("The potential has no centers")
        self.groups = []
        charges = [[], []]
        for atom in self.atoms:
            center = np.array(atom["position"], dtype=float)
            if center.shape != (3,) or not np.all(np.isfinite(center)):
                raise ValueError("Invalid atomic coordinates")
            groups = {}
            for term in atom["terms"]:
                if len(term) != 6 or not np.all(np.isfinite(term)):
                    raise ValueError("Invalid Gaussian term")
                powers = tuple(int(p) for p in term[:3])
                if any(p < 0 or p != raw for p, raw in zip(powers, term[:3])):
                    raise ValueError("Cartesian powers must be nonnegative integers")
                exponent = float(term[3])
                if exponent <= 0:
                    raise ValueError("Gaussian exponents must be positive")
                weights = np.array(term[4:], dtype=float)
                group = groups.setdefault(exponent, {})
                group[powers] = group.get(powers, np.zeros(2)) + weights
            for exponent, terms in groups.items():
                powers = list(terms)
                weights = np.array(list(terms.values()))
                self.groups.append((center, exponent, powers, weights))
                for triple, pair in terms.items():
                    moment = (
                        0.0
                        if any(p % 2 for p in triple)
                        else math.prod(
                            gamma((p + 1) / 2) / exponent ** ((p + 1) / 2)
                            for p in triple
                        )
                    )
                    for component in range(2):
                        charges[component].append(moment * pair[component])
        self.charges = np.array([math.fsum(values) for values in charges])
        expected = np.array([data["total_charge"], data["reference_charge"]])
        if not np.all(np.isfinite(expected)) or not np.allclose(
            self.charges, expected, rtol=1e-9, atol=1e-8
        ):
            raise ValueError(
                "Gaussian expansion disagrees with the exported charge checks"
            )

    def evaluate(self, points, batch_size=1024):
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
            raise ValueError("Expected finite coordinates with shape (npoints, 3)")
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        result = np.zeros((len(points), 3))
        for start in range(0, len(points), batch_size):
            stop = min(start + batch_size, len(points))
            for center, exponent, powers, weights in self.groups:
                integrals = gaussian_potentials(
                    center - points[start:stop], exponent, powers
                )
                result[start:stop, :2] += integrals @ weights
        result[:, 2] = result[:, 0] - result[:, 1]
        return result

    def write_cube(self, filename, spacing=0.2, margin=5.0, component="total"):
        """Write a molecular cube in bounded batches, including its atomic centers."""
        if spacing <= 0 or margin <= 0:
            raise ValueError("Cube spacing and margin must be positive")
        column = ["total", "reference", "remainder"].index(component)
        positions = np.array([atom["position"] for atom in self.atoms])
        origin = positions.min(axis=0) - margin
        extent = positions.max(axis=0) - positions.min(axis=0) + 2 * margin
        shape = np.ceil(extent / spacing).astype(int) + 1
        steps = extent / (shape - 1)
        with Path(filename).open("w") as stream:
            stream.write(f"CP2K EXX-OEP {component} exchange potential\n")
            stream.write("Coordinates in bohr; potential in hartree\n")
            stream.write(
                f"{len(self.atoms):5d}" + "".join(f" {v:.16e}" for v in origin) + "\n"
            )
            stream.writelines(
                f"{size:5d}" + "".join(f" {v:.16e}" for v in vector) + "\n"
                for size, vector in zip(shape, np.diag(steps))
            )
            stream.writelines(
                f"{atom['atomic_number']:5d} 0.0"
                + "".join(f" {v:.16e}" for v in atom["position"])
                + "\n"
                for atom in self.atoms
            )
            size = int(np.prod(shape))
            for start in range(0, size, 1020):
                indices = np.array(
                    np.unravel_index(np.arange(start, min(start + 1020, size)), shape)
                ).T
                values = self.evaluate(origin + indices * steps)[:, column]
                full = len(values) // 6 * 6
                np.savetxt(stream, values[:full].reshape(-1, 6), fmt="%.10e")
                if full < len(values):
                    np.savetxt(stream, values[full:][None, :], fmt="%.10e")
        return tuple(shape)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("potential", type=Path)
    commands = parser.add_subparsers(dest="command", required=True)
    line = commands.add_parser(
        "line", help="Sample a line and optionally plot all three components"
    )
    line.add_argument(
        "--start", type=float, nargs=3, required=True, metavar=("X", "Y", "Z")
    )
    line.add_argument(
        "--end", type=float, nargs=3, required=True, metavar=("X", "Y", "Z")
    )
    line.add_argument("--points", type=int, default=1601)
    line.add_argument("--output", type=Path, required=True, help="CSV filename")
    line.add_argument("--plot", type=Path, help="Optional PNG, PDF, or SVG filename")
    points = commands.add_parser(
        "points", help="Evaluate arbitrary points from a three-column text file"
    )
    points.add_argument(
        "coordinates", type=Path, help="Whitespace-separated x y z in bohr"
    )
    points.add_argument("--output", type=Path, required=True)
    cube = commands.add_parser("cube", help="Write a Gaussian cube file")
    cube.add_argument("--output", type=Path, required=True)
    cube.add_argument(
        "--spacing", type=float, default=0.2, help="Maximum spacing in bohr"
    )
    cube.add_argument(
        "--margin", type=float, default=5, help="Margin around atoms in bohr"
    )
    cube.add_argument(
        "--component", choices=["total", "reference", "remainder"], default="total"
    )
    args = parser.parse_args()
    potential = OEPPotential(args.potential)
    if args.command == "cube":
        shape = potential.write_cube(
            args.output, args.spacing, args.margin, args.component
        )
        print(f"Wrote {args.output}: {shape[0]} x {shape[1]} x {shape[2]} grid")
        return
    if args.command == "line":
        if args.points < 2 or np.array_equal(args.start, args.end):
            parser.error("A line needs distinct endpoints and at least two points")
        coordinates = np.linspace(args.start, args.end, args.points)
    else:
        coordinates = np.loadtxt(args.coordinates, ndmin=2)
    values = potential.evaluate(coordinates)
    np.savetxt(
        args.output,
        np.column_stack((coordinates, values)),
        delimiter=",",
        header="x_bohr,y_bohr,z_bohr,vx_total_hartree,vx_reference_hartree,vx_remainder_hartree",
        comments="",
        fmt="%.16e",
    )
    print(f"Wrote {args.output}: {len(coordinates)} points")
    if args.command == "line" and args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        distances = np.linalg.norm(coordinates - coordinates[0], axis=1)
        direction = np.array(args.end) - args.start
        axes = np.flatnonzero(direction)
        abscissa = coordinates[:, axes[0]] if len(axes) == 1 else distances
        xlabel = (
            f"{'xyz'[axes[0]]} (bohr)"
            if len(axes) == 1
            else "Distance along line (bohr)"
        )
        fig, ax = plt.subplots(figsize=(8, 4.8), layout="constrained")
        for column, label, color, style in [
            (0, "Total", "#126e82", "-"),
            (1, "Reference", "#c86420", "--"),
            (2, "Remainder", "#7570b3", ":"),
        ]:
            ax.plot(
                abscissa, values[:, column], label=label, color=color, linestyle=style
            )
        ax.set(
            xlabel=xlabel,
            ylabel=r"Exchange potential $v_x$ (hartree)",
            title="CP2K EXX-OEP",
        )
        ax.legend(frameon=False)
        ax.grid(alpha=0.15)
        fig.savefig(args.plot, dpi=180)
        plt.close(fig)
        print(f"Wrote {args.plot}")


if __name__ == "__main__":
    main()
