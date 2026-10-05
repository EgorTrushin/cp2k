#!/usr/bin/env python3
"""Independent Gaussian-integral and portable-file tests; NumPy/SciPy only."""

import copy
import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
from potential import OEPPotential, gaussian_potentials
from scipy.integrate import quad
from scipy.special import erf


def quadrature_potential(point, exponent, powers):
    """Independent Laplace-transform quadrature with analytic Gaussian moments."""

    def integrand(u):
        product = math.exp(-exponent * u * u * np.dot(point, point))
        for coordinate, power in zip(point, powers):
            moment = 0.0
            for k in range(power // 2 + 1):
                coefficient = math.factorial(power) / (
                    math.factorial(power - 2 * k) * math.factorial(k)
                )
                moment += (
                    coefficient
                    * ((1 - u * u) / (4 * exponent)) ** k
                    * (u * u * coordinate) ** (power - 2 * k)
                )
            product *= moment
        return product

    return (2 * np.pi / exponent) * quad(integrand, 0, 1, epsabs=1e-12, epsrel=1e-12)[0]


class GaussianTests(unittest.TestCase):
    def test_monopole_at_center_near_and_far(self):
        for exponent in (1e-4, 0.7, 1e4):
            radii = np.array([0, 1e-10, 0.1, 1, 10, 1000]) / np.sqrt(exponent)
            coordinates = np.zeros((len(radii), 3))
            coordinates[:, 2] = radii
            values = gaussian_potentials(-coordinates, exponent, [(0, 0, 0)])[:, 0]
            charge = (np.pi / exponent) ** 1.5
            expected = np.empty_like(radii)
            expected[0] = 2 * np.pi / exponent
            expected[1:] = charge * erf(np.sqrt(exponent) * radii[1:]) / radii[1:]
            np.testing.assert_allclose(values, expected, rtol=3e-13)

    def test_cartesian_multipoles_against_quadrature(self):
        powers = [
            (1, 0, 0),
            (0, 1, 0),
            (2, 0, 0),
            (1, 1, 0),
            (3, 0, 0),
            (1, 1, 1),
            (3, 1, 0),
            (1, 2, 2),
            (6, 0, 0),
            (2, 2, 2),
        ]
        points = np.array([[0, 0, 0], [0.3, -0.4, 0.2], [-2.1, 0.4, 1.3]])
        for exponent in (0.03, 0.7, 10):
            values = gaussian_potentials(-points, exponent, powers)
            expected = [
                [quadrature_potential(point, exponent, power) for power in powers]
                for point in points
            ]
            np.testing.assert_allclose(values, expected, rtol=2e-11, atol=1e-10)


class FileTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "potential.json"
        self.data = {
            "format": "CP2K_OEP_POTENTIAL",
            "version": 1,
            "length_unit": "bohr",
            "potential_unit": "hartree",
            "representation": "unnormalized_cartesian_gaussians",
            "term_columns": ["lx", "ly", "lz", "exponent", "total", "reference"],
            "total_charge": -1.0,
            "reference_charge": -1.0,
            "atoms": [
                {
                    "element": "H",
                    "atomic_number": 1,
                    "position": [0.3, -0.2, 0.5],
                    "terms": [
                        [
                            0,
                            0,
                            0,
                            0.7,
                            -0.5 * (0.7 / np.pi) ** 1.5,
                            -0.5 * (0.7 / np.pi) ** 1.5,
                        ],
                        [1, 0, 0, 0.7, 0.2, 0.1],
                        [1, 0, 0, 0.7, 0.1, -0.1],
                    ],
                },
                {
                    "element": "H",
                    "atomic_number": 1,
                    "position": [-0.4, 0.1, -0.6],
                    "terms": [
                        [
                            0,
                            0,
                            0,
                            0.3,
                            -0.5 * (0.3 / np.pi) ** 1.5,
                            -0.5 * (0.3 / np.pi) ** 1.5,
                        ]
                    ],
                },
            ],
        }
        self.path.write_text(json.dumps(self.data))

    def test_components_batching_and_tail(self):
        potential = OEPPotential(self.path)
        points = np.random.default_rng(1).normal(size=(17, 3))
        values = potential.evaluate(points)
        np.testing.assert_allclose(
            values, potential.evaluate(points, batch_size=3), rtol=1e-13, atol=1e-14
        )
        np.testing.assert_allclose(
            values[:, 0], values[:, 1] + values[:, 2], atol=1e-15
        )
        np.testing.assert_allclose(potential.charges, [-1, -1], atol=1e-14)
        far = np.r_[np.eye(3), -np.eye(3)] * 1e4
        np.testing.assert_allclose(
            potential.evaluate(far).mean(axis=0) * 1e4, [-1, -1, 0], atol=1e-10
        )
        self.assertEqual(potential.evaluate(np.empty((0, 3))).shape, (0, 3))

    def test_cube_coordinates_ordering_and_values(self):
        potential = OEPPotential(self.path)
        cube = Path(self.directory.name) / "potential.cube"
        shape = potential.write_cube(
            cube, spacing=0.45, margin=0.4, component="remainder"
        )
        lines = cube.read_text().splitlines()
        natoms, *origin = lines[2].split()
        vectors = np.array(
            [[float(v) for v in line.split()[1:]] for line in lines[3:6]]
        )
        parsed_shape = tuple(int(line.split()[0]) for line in lines[3:6])
        self.assertEqual(parsed_shape, shape)
        indices = np.indices(shape).reshape(3, -1).T
        coordinates = np.array(origin, dtype=float) + indices @ vectors
        values = np.fromstring(" ".join(lines[6 + int(natoms) :]), sep=" ")
        np.testing.assert_allclose(
            values, potential.evaluate(coordinates)[:, 2], atol=1e-11, rtol=1e-10
        )

    def test_reject_incompatible_or_damaged_export(self):
        changes = [
            ("version", 2),
            ("length_unit", "angstrom"),
            ("total_charge", -2),
            ("term_columns", ["lx", "ly", "lz", "exponent", "reference", "total"]),
        ]
        for key, value in changes:
            data = copy.deepcopy(self.data)
            data[key] = value
            self.path.write_text(json.dumps(data))
            with self.assertRaises(ValueError):
                OEPPotential(self.path)
        data = copy.deepcopy(self.data)
        data["atoms"][0]["terms"][0][3] = -0.7
        self.path.write_text(json.dumps(data))
        with self.assertRaises(ValueError):
            OEPPotential(self.path)


if __name__ == "__main__":
    unittest.main()
