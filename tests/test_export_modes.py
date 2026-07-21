"""Tests for export_modes.py — the phonopy-format mode emission.

Checks the band.yaml-schema structure parses with a standard YAML reader,
eigenvectors are unit-norm complex per phonopy's convention, the q-points
match the declared star arms, and the animations exist. Requires the GTS
data + reference (skipped otherwise); no MLIP involved.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import export_modes as em

REPO = Path(__file__).resolve().parent.parent
NEED = [REPO / "data/ensemble_20A_5K/GTS_5K_1.rmc6f",
        REPO / "reference/gts_mode_patterns.json"]
pytestmark = pytest.mark.skipif(not all(p.is_file() for p in NEED),
                                reason="GTS data/reference not present")


@pytest.fixture(scope="module")
def emitted(tmp_path_factory):
    out = tmp_path_factory.mktemp("modes")
    em.main(["--outdir", str(out), "--nframes", "8",
             "--verdicts", str(REPO / "results/verdicts.json")])
    return out


def test_yaml_parses_with_standard_reader(emitted):
    import yaml

    d = yaml.safe_load((emitted / "modes_irrep.yaml").read_text())
    assert d["natom"] == 52
    assert d["nqpoint"] == len(d["phonon"]) == 7   # 6 irreps + TOTAL
    assert len(d["points"]) == 52
    np.testing.assert_allclose(np.array(d["lattice"]),
                               np.eye(3) * em.A_CUB)


def test_eigenvectors_unit_norm_and_shape(emitted):
    import yaml

    d = yaml.safe_load((emitted / "modes_irrep.yaml").read_text())
    for p in d["phonon"]:
        for b in p["band"]:
            ev = np.array(b["eigenvector"])       # (52, 3, 2)
            assert ev.shape == (52, 3, 2)
            z = ev[..., 0] + 1j * ev[..., 1]
            assert np.linalg.norm(z) == pytest.approx(1.0, abs=1e-6)


def test_qpoints_match_declared_stars(emitted):
    import yaml

    d = yaml.safe_load((emitted / "modes_irrep.yaml").read_text())
    for p in d["phonon"]:
        key = p["label"].split()[0]
        arms = [list(a) for a in em.STAR_ARMS[key]]
        assert [round(x, 4) for x in p["q-position"]] in arms, p["label"]


def test_gamma_pattern_is_real_and_x_pattern_captured(emitted):
    import yaml

    d = yaml.safe_load((emitted / "modes_irrep.yaml").read_text())
    by = {p["label"].split()[0]: p for p in d["phonon"]}
    # Gamma modes: q = 0 -> eigenvector purely real
    for key in ("G1", "G3"):
        ev = np.array(by[key]["band"][0]["eigenvector"])
        assert np.abs(ev[..., 1]).max() < 1e-10, key
    # X5: a single arm captures the (parity-even) pattern fully
    w = float(by["X5"]["label"].split("weight")[1].strip(" )"))
    assert w > 0.8


def test_animations_written(emitted):
    xyz = sorted(emitted.glob("*.xyz"))
    assert len(xyz) == 7
    from ase.io import read
    frames = read(str(emitted / "X5.xyz"), index=":")
    assert len(frames) == 8
    assert len(frames[0]) == 104
