"""Unit tests for milestone1_bands.fold_and_average.

Covers the circular (wrap-around) site mean, --ref CIF site assignment when
rmc6f files lack a site-id column, and mixed-occupancy majority reduction.
"""

from __future__ import annotations

import numpy as np
import pytest

from conftest import make_config

CUBIC = np.diag([4.0, 4.0, 4.0])


def _wrap_dist(a, b):
    """Minimum-image distance between two fractional coordinate arrays."""
    d = (np.asarray(a) - np.asarray(b) + 0.5) % 1.0 - 0.5
    return np.abs(d)


def test_circular_mean_wraps_around_zero(m1):
    """Coordinates straddling the 0/1 seam average to 0, not to their
    arithmetic mean (~1/3)."""
    dims = [1, 1, 1]
    configs = [
        make_config([[x, 0.0, 0.0]], ["Cu"], dims, CUBIC, site_ids=[1])
        for x in (0.9, 0.0, 0.1)
    ]
    _, elements, site_frac, report = m1.fold_and_average(configs, None)

    assert elements == ["Cu"]
    assert report["n_sites"] == 1
    # Circular mean sits at the seam (~0.0 ≡ 1.0), far from the linear mean 1/3.
    assert np.all(_wrap_dist(site_frac[0], [0.0, 0.0, 0.0]) < 1e-9)
    assert _wrap_dist(site_frac[0, 0], 1.0 / 3.0)[()] > 0.3


def test_circular_mean_matches_linear_away_from_seam(m1):
    """Away from the wrap seam the circular mean equals the ordinary mean."""
    dims = [1, 1, 1]
    xs = (0.30, 0.40, 0.50)
    configs = [
        make_config([[x, 0.25, 0.75]], ["Cu"], dims, CUBIC, site_ids=[1])
        for x in xs
    ]
    _, _, site_frac, _ = m1.fold_and_average(configs, None)
    np.testing.assert_allclose(site_frac[0], [np.mean(xs), 0.25, 0.75], atol=1e-9)


def test_fold_collapses_supercell_images(m1):
    """Two atoms one unit cell apart in a 2×1×1 supercell fold to one site."""
    cell = np.diag([8.0, 4.0, 4.0])  # 2×1×1 of a 4 Å cube
    dims = [2, 1, 1]
    # Same unit-cell site (0.25, 0.5, 0.5): image 0 at x=0.125, image 1 at 0.625.
    cfg = make_config(
        [[0.125, 0.5, 0.5], [0.625, 0.5, 0.5]], ["Cu", "Cu"], dims, cell,
        site_ids=[1, 1],
    )
    unit_cell, _, site_frac, report = m1.fold_and_average([cfg], None)

    np.testing.assert_allclose(unit_cell, np.diag([4.0, 4.0, 4.0]))
    assert report["n_sites"] == 1
    assert np.all(_wrap_dist(site_frac[0], [0.25, 0.5, 0.5]) < 1e-9)


def test_ref_cif_site_assignment(tmp_path, m1):
    """Without site ids, atoms are grouped by nearest reference-CIF site,
    independent of their order within a configuration."""
    from ase import Atoms
    from ase.io import write as ase_write

    ref = Atoms(
        symbols=["Cu", "O"],
        scaled_positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        cell=CUBIC,
        pbc=True,
    )
    ref_cif = tmp_path / "parent.cif"
    ase_write(str(ref_cif), ref)

    # Config A: [near Cu-site, near O-site]; config B swaps the atom order.
    cfg_a = make_config(
        [[0.02, 0.0, 0.98], [0.48, 0.52, 0.50]], ["Cu", "O"], [1, 1, 1], CUBIC
    )
    cfg_b = make_config(
        [[0.49, 0.51, 0.50], [0.98, 0.01, 0.02]], ["O", "Cu"], [1, 1, 1], CUBIC
    )
    _, elements, site_frac, report = m1.fold_and_average([cfg_a, cfg_b], ref_cif)

    assert report["n_sites"] == 2
    assert report["mixed_occupancy_sites"] == []
    # Identify sites by position (robust to CIF round-trip atom ordering):
    # the atoms near (0,0,0) are Cu, those near (½,½,½) are O — regardless of
    # their order within each configuration.
    at_origin = [el for el, f in zip(elements, site_frac)
                 if np.all(_wrap_dist(f, [0.0, 0.0, 0.0]) < 0.05)]
    at_body_center = [el for el, f in zip(elements, site_frac)
                      if np.all(_wrap_dist(f, [0.5, 0.5, 0.5]) < 0.05)]
    assert at_origin == ["Cu"]
    assert at_body_center == ["O"]


def test_missing_ref_cif_raises(m1):
    """A missing --ref for id-less rmc6f is a hard, actionable error."""
    cfg = make_config([[0.0, 0.0, 0.0]], ["Cu"], [1, 1, 1], CUBIC)
    with pytest.raises(SystemExit, match="site-id column"):
        m1.fold_and_average([cfg], None)


def test_mixed_occupancy_reduces_to_majority(m1):
    """A site occupied by different elements across configs collapses to the
    majority element and is reported for the caller's warning."""
    dims = [1, 1, 1]
    occupants = ["Cu", "Cu", "Cu", "Au", "Au"]  # majority Cu (3 vs 2)
    configs = [
        make_config([[0.0, 0.0, 0.0]], [el], dims, CUBIC, site_ids=[1])
        for el in occupants
    ]
    _, elements, _, report = m1.fold_and_average(configs, None)

    assert elements == ["Cu"]  # majority wins
    assert len(report["mixed_occupancy_sites"]) == 1
    sid, votes = report["mixed_occupancy_sites"][0]
    assert sid == 1
    assert votes == {"Cu": 3, "Au": 2}


def test_disagreeing_supercell_dims_raise(m1):
    cfg_a = make_config([[0.0, 0.0, 0.0]], ["Cu"], [2, 2, 2], CUBIC, site_ids=[1])
    cfg_b = make_config([[0.0, 0.0, 0.0]], ["Cu"], [2, 2, 1], CUBIC, site_ids=[1])
    with pytest.raises(ValueError, match="supercell dimensions"):
        m1.fold_and_average([cfg_a, cfg_b], None)
