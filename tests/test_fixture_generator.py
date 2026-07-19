"""Tests for the synthetic Cu fcc ensemble generator.

The generator's contract: writing N noisy configs and then running the
milestone-1 fold/circular-average over them must recover the ideal fcc unit
cell (4 sites, all Cu, at the conventional fcc basis positions). Output must
be deterministic in the seed.
"""

from __future__ import annotations

import numpy as np

from fixtures.make_synthetic_ensemble import (
    CU_A,
    FCC_BASIS,
    build_supercell,
    make_fcc_cu_ensemble,
)


def _wrap_dist(a, b):
    d = (np.asarray(a) - np.asarray(b) + 0.5) % 1.0 - 0.5
    return np.abs(d)


def test_build_supercell_atom_count_and_ids():
    dims = (2, 2, 2)
    cell, frac, site_ids, offsets = build_supercell(dims)
    n = 4 * 2 * 2 * 2
    assert frac.shape == (n, 3)
    np.testing.assert_allclose(cell, np.diag([2 * CU_A] * 3))
    # Each of the 4 basis sites appears once per cell (8 cells).
    counts = np.bincount(site_ids)[1:]
    np.testing.assert_array_equal(counts, [8, 8, 8, 8])
    assert offsets.shape == (n, 3)
    assert set(map(tuple, offsets)) == {
        (i, j, k) for i in range(2) for j in range(2) for k in range(2)
    }


def test_generator_writes_n_parseable_configs(tmp_path, m1):
    paths = make_fcc_cu_ensemble(tmp_path, n_configs=5, dims=(2, 2, 2), seed=0)
    assert len(paths) == 5
    for p in paths:
        cfg = m1.parse_rmc6f(p)
        assert cfg["elements"] == ["Cu"] * 32
        assert cfg["site_ids"] is not None
        np.testing.assert_array_equal(cfg["dims"], [2, 2, 2])
        np.testing.assert_allclose(cfg["cell"], np.diag([2 * CU_A] * 3), atol=1e-4)


def test_ensemble_average_recovers_fcc_cell(tmp_path, m1):
    """Folding + circular averaging the noisy ensemble returns the ideal
    4-site fcc basis."""
    paths = make_fcc_cu_ensemble(
        tmp_path, n_configs=32, dims=(2, 2, 2), sigma=0.08, seed=1
    )
    configs = [m1.parse_rmc6f(p) for p in paths]
    unit_cell, elements, site_frac, report = m1.fold_and_average(configs, None)

    np.testing.assert_allclose(unit_cell, np.diag([CU_A] * 3), atol=1e-4)
    assert elements == ["Cu", "Cu", "Cu", "Cu"]
    assert report["mixed_occupancy_sites"] == []
    # Sites are sorted by site id -> aligned with FCC_BASIS order.
    for recovered, ideal in zip(site_frac, FCC_BASIS):
        assert np.all(_wrap_dist(recovered, ideal) < 0.02)


def test_generator_is_deterministic(tmp_path):
    a = make_fcc_cu_ensemble(tmp_path / "a", n_configs=3, seed=7)
    b = make_fcc_cu_ensemble(tmp_path / "b", n_configs=3, seed=7)
    c = make_fcc_cu_ensemble(tmp_path / "c", n_configs=3, seed=8)
    a_txt = [p.read_text() for p in a]
    b_txt = [p.read_text() for p in b]
    c_txt = [p.read_text() for p in c]
    assert a_txt == b_txt  # same seed -> identical bytes
    assert a_txt != c_txt  # different seed -> different clouds


def test_noise_free_ensemble_is_exactly_ideal(tmp_path, m1):
    """With sigma=0 the recovered sites are exact (no averaging error)."""
    paths = make_fcc_cu_ensemble(tmp_path, n_configs=2, sigma=0.0, seed=0)
    configs = [m1.parse_rmc6f(p) for p in paths]
    _, _, site_frac, _ = m1.fold_and_average(configs, None)
    for recovered, ideal in zip(site_frac, FCC_BASIS):
        assert np.all(_wrap_dist(recovered, ideal) < 1e-6)
