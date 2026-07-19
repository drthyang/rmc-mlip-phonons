"""Unit tests for milestone1_bands.auto_dim (phonopy supercell heuristic).

auto_dim picks the smallest supercell whose every axis reaches ~`target`
Angstrom, never below 1, and warns (but still returns) when the resulting
supercell would exceed `max_atoms`.
"""

from __future__ import annotations

import numpy as np
import pytest


def _atoms(cell_lengths, n=1):
    from ase import Atoms

    cell = np.diag(np.asarray(cell_lengths, dtype=float))
    return Atoms("Cu" * n, scaled_positions=np.zeros((n, 3)), cell=cell, pbc=True)


@pytest.mark.parametrize(
    "lengths, expected",
    [
        ([3.6, 3.6, 3.6], [4, 4, 4]),  # ceil(12/3.6)=4
        ([7.2, 7.2, 7.2], [2, 2, 2]),  # ceil(12/7.2)=2
        ([3.0, 6.0, 13.0], [4, 2, 1]),  # anisotropic, one axis already >12
        ([20.0, 20.0, 20.0], [1, 1, 1]),  # never smaller than 1
    ],
)
def test_auto_dim_reaches_target(m1, lengths, expected):
    np.testing.assert_array_equal(m1.auto_dim(_atoms(lengths)), expected)


def test_auto_dim_respects_custom_target(m1):
    np.testing.assert_array_equal(
        m1.auto_dim(_atoms([4.0, 4.0, 4.0]), target=6.0), [2, 2, 2]
    )


def test_auto_dim_warns_when_over_max_atoms(m1, capsys):
    dim = m1.auto_dim(_atoms([1.2, 1.2, 1.2]))  # ceil(12/1.2)=10 -> 1000 atoms
    np.testing.assert_array_equal(dim, [10, 10, 10])
    out = capsys.readouterr().out
    assert "consider --dim" in out
    assert "1000 atoms" in out
