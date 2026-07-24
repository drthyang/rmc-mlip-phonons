"""Unit tests for milestone1_bands.symmetrize_lattice.

The routine removes ~fmax-level numerical noise from a relaxed cell's metric
tensor so that seekpath — which phonopy calls at a hard-coded symprec of 1e-5
— classifies the Bravais lattice correctly. Units: Å throughout.
"""

from __future__ import annotations

import numpy as np
import pytest

from ase import Atoms

import milestone1_bands as m1

CU_A = 3.61


def fcc_primitive(a=CU_A, noise=0.0, seed=0):
    """One-atom fcc primitive cell, optionally with noisy lattice vectors."""
    cell = 0.5 * a * np.array([[0.0, 1.0, 1.0],
                               [1.0, 0.0, 1.0],
                               [1.0, 1.0, 0.0]])
    if noise:
        cell = cell + np.random.default_rng(seed).normal(0.0, noise, (3, 3))
    return Atoms("Cu", scaled_positions=[[0.0, 0.0, 0.0]], cell=cell, pbc=True)


def test_exact_cell_is_a_fixed_point():
    """An already-symmetric cell must come back unchanged."""
    atoms = fcc_primitive()
    out, rep = m1.symmetrize_lattice(atoms, symprec=1e-3)
    assert rep["max_lattice_shift_A"] < 1e-12
    assert np.allclose(out.cell.array, atoms.cell.array, atol=1e-12)


def test_noise_is_removed_and_metric_becomes_exact():
    """Noise at 1e-5 Å is removed; the metric acquires exact cubic form."""
    atoms = fcc_primitive(noise=2e-5, seed=1)
    g_before = atoms.cell.array @ atoms.cell.array.T
    # off-diagonals of the fcc primitive metric are all a²/4, diagonals a²/2
    assert np.ptp(np.diag(g_before)) > 1e-6      # noise present

    out, rep = m1.symmetrize_lattice(atoms, symprec=1e-3)
    g = out.cell.array @ out.cell.array.T
    assert np.ptp(np.diag(g)) < 1e-10
    off = g[~np.eye(3, dtype=bool)]
    assert np.ptp(off) < 1e-10
    assert np.diag(g).mean() == pytest.approx(2.0 * off.mean(), rel=1e-9)
    assert 0.0 < rep["max_lattice_shift_A"] < 1e-3
    assert rep["n_operations"] == 48             # m-3m


def test_fractional_coordinates_and_order_are_preserved():
    """Multi-site cell: scaled positions, order and species must not move."""
    cell = np.diag([4.0, 4.0 + 3e-5, 4.0 - 2e-5])
    frac = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    atoms = Atoms("NaCl", scaled_positions=frac, cell=cell, pbc=True)

    out, _ = m1.symmetrize_lattice(atoms, symprec=1e-3)
    assert out.get_chemical_symbols() == ["Na", "Cl"]
    assert np.allclose(out.get_scaled_positions(), frac, atol=1e-12)
    assert len(out) == len(atoms)


def test_orientation_is_preserved():
    """A' = S·A with S -> I, so the frame must not rotate or reflect."""
    atoms = fcc_primitive(noise=2e-5, seed=3)
    out, _ = m1.symmetrize_lattice(atoms, symprec=1e-3)
    # each lattice vector stays essentially parallel to its original
    for v0, v1 in zip(atoms.cell.array, out.cell.array):
        cos = v0 @ v1 / (np.linalg.norm(v0) * np.linalg.norm(v1))
        assert cos > 1.0 - 1e-8
    assert np.linalg.det(out.cell.array) > 0     # handedness kept


def test_no_symmetry_returns_input_unchanged():
    """P1 with a tight tolerance: nothing to average, cell passes through."""
    rng = np.random.default_rng(7)
    cell = np.eye(3) * 5.0 + rng.normal(0, 0.3, (3, 3))
    frac = rng.random((4, 3))
    atoms = Atoms("HHeLiBe", scaled_positions=frac, cell=cell, pbc=True)

    out, rep = m1.symmetrize_lattice(atoms, symprec=1e-6)
    assert rep["n_operations"] <= 1
    assert np.allclose(out.cell.array, cell, atol=1e-12)


def test_seekpath_recovers_the_cubic_path():
    """The regression this fixes: seekpath at its default 1e-5 must see cF."""
    seekpath = pytest.importorskip("seekpath")
    noisy = fcc_primitive(noise=1e-5, seed=11)

    def bravais(a):
        return seekpath.get_path(
            (a.cell.array, a.get_scaled_positions(),
             a.get_atomic_numbers()))["bravais_lattice_extended"]

    assert bravais(noisy) != "cF2"               # the bug
    out, _ = m1.symmetrize_lattice(noisy, symprec=1e-3)
    assert bravais(out) == "cF2"                 # the fix
