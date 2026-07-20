"""Tests for hiphive_fit.py — geometry mapping (unit) and the EMT
end-to-end fit on the Cu fixture (slow).

The physics guard: effective FCs fitted from noisy fixture snapshots +
EMT forces must reproduce the EMT harmonic bands (small displacements probe
the harmonic regime), and the outputs must exist per the contract.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

import hiphive_fit as hf
import milestone1_bands as m1
from fixtures.make_synthetic_ensemble import make_fcc_cu_ensemble


def test_config_site_cell_roundtrip(tmp_path):
    """(unit_frac, sid, ijk) derived from frac*dims reconstructs positions."""
    make_fcc_cu_ensemble(tmp_path, n_configs=1, dims=(2, 3, 2), seed=0)
    cfg = m1.parse_rmc6f(next(tmp_path.glob("*.rmc6f")))
    unit, sid, ijk = hf.config_site_cell(cfg)
    assert unit.shape == cfg["frac"].shape
    assert np.all((unit >= 0) & (unit < 1))
    assert np.all(ijk >= 0) and np.all(ijk < np.array([2, 3, 2]))
    np.testing.assert_allclose((unit + ijk) / np.array([2, 3, 2]),
                               cfg["frac"], atol=1e-9)


def test_box_order_is_a_permutation():
    rng = np.random.default_rng(0)
    dims = (2, 2, 2)
    n_sites = 4
    sid = np.tile(np.arange(1, 5), 8)
    ijk = np.repeat(np.array([(i, j, k) for i in range(2) for j in range(2)
                              for k in range(2)]), 4, axis=0)
    p = rng.permutation(len(sid))
    order = hf.box_order(sid[p], ijk[p], dims, n_sites)
    assert sorted(order) == list(range(32))


def test_ideal_box_matches_config_order(tmp_path):
    """A noise-free fixture box, mapped through build_box, must coincide
    with ideal_box built from the folded reference — same order, same
    positions."""
    make_fcc_cu_ensemble(tmp_path, n_configs=1, dims=(2, 2, 2), sigma=0.0,
                         seed=0)
    cfg = m1.parse_rmc6f(next(tmp_path.glob("*.rmc6f")))
    unit_cell, elems, frac, _ = m1.fold_and_average([cfg], None)
    unit, sid, ijk = hf.config_site_cell(cfg)
    box = hf.build_box(unit, sid, ijk, cfg["dims"], list(elems), unit_cell)
    ideal = hf.ideal_box(frac, list(elems), cfg["dims"], unit_cell)
    assert box.get_chemical_symbols() == ideal.get_chemical_symbols()
    d = box.get_scaled_positions() - ideal.get_scaled_positions()
    d -= np.round(d)
    assert np.abs(d).max() < 1e-6


@pytest.mark.slow
def test_emt_end_to_end_fit(tmp_path):
    pytest.importorskip("hiphive")
    ens = tmp_path / "ens"
    out = tmp_path / "out"
    # 3x3x3 box (10.8 Å) so the cutoff can cover EMT's full ~4.1 Å range —
    # in a 2x2x2 box the second shell (3.61 Å) exceeds the half-box and the
    # fit is legitimately incomplete (rmse ~0.07 eV/Å).
    make_fcc_cu_ensemble(ens, n_configs=6, dims=(3, 3, 3), sigma=0.05, seed=0)
    hf.main([str(ens), "--calc", "emt", "--nconfigs", "6",
             "--cutoff", "4.2", "--band-dim", "3", "3", "3",
             "--npoints", "11", "-o", str(out)])
    assert (out / "band_rmc.yaml").is_file()
    rep = json.loads((out / "fit_report.json").read_text())
    assert rep["n_snapshots"] == 6
    # the rmse floor is the ANHARMONIC force component (~Phi3*sigma^2),
    # which a 2nd-order model cannot represent — not a fit-quality failure
    assert rep["rmse_test_eV_A"] < 0.10
    # small-amplitude snapshots probe the harmonic regime: effective bands
    # track the EMT harmonic bands
    assert rep["max_abs_dw_vs_harmonic_THz"] < 1.0
    # force cache written and reused
    assert len(list(out.glob("forces_*.npy"))) == 6
