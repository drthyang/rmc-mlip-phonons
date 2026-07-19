"""Real-dependency integration test: MACE-MP-0 (small) on the fcc-Cu fixture.

This is the milestone-1 acceptance run — a universal MLIP (not a toy potential)
on the synthetic ensemble. fcc Cu is dynamically stable, so MACE must produce
no imaginary modes (min band frequency > -0.05 THz) and acoustic branches that
vanish at Γ.

Marked ``slow`` and ``mace``. Skipped automatically when mace-torch is not
installed. On first run MACE downloads the ~31 MB MACE-MP-0 small weights and
caches them (~/.cache/mace); afterwards it is offline.

Regression guard: with phonopy's DFT-tuned 0.01 Å displacement this run showed
a spurious ~-0.33 THz mode just off Γ (MLIP force noise). The pipeline default
is 0.03 Å; see milestone1_bands.py --displacement.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("mace", reason="mace-torch not installed")

from fixtures.make_synthetic_ensemble import make_fcc_cu_ensemble  # noqa: E402

pytestmark = [pytest.mark.slow, pytest.mark.mace]


@pytest.fixture(scope="module")
def mace_outputs(tmp_path_factory, m1):
    import json

    import yaml

    root = tmp_path_factory.mktemp("mace")
    ens = root / "ens"
    out = root / "out"
    make_fcc_cu_ensemble(ens, n_configs=4, dims=(2, 2, 2), sigma=0.06, seed=0)

    # Default supercell (auto) and default 0.03 Å displacement — this exercises
    # the shipped defaults end to end.
    m1.main([
        str(ens), "--calc", "mace", "--model", "small", "-o", str(out),
        "--npoints", "21", "--no-eigenvectors",
    ])

    band = yaml.safe_load((out / "band.yaml").read_text())
    summary = json.loads((out / "summary.json").read_text())
    return out, band, summary


def test_mace_writes_all_outputs(mace_outputs):
    out, _, summary = mace_outputs
    for name in ("band.yaml", "relaxed.cif", "summary.json"):
        assert (out / name).is_file(), f"missing {name}"
    assert summary["calculator"]["name"] == "mace"
    assert summary["calculator"]["dtype"] == "float64"  # CLAUDE.md hard rule


def test_mace_three_branches_vanish_at_gamma(mace_outputs):
    _, band, _ = mace_outputs
    assert band["natom"] == 1  # fcc primitive
    assert len(band["phonon"][0]["band"]) == 3
    for p in band["phonon"]:
        if np.allclose(p["q-position"], 0.0, atol=1e-9):
            freqs = [b["frequency"] for b in p["band"]]
            assert np.allclose(freqs, 0.0, atol=2e-2), freqs
            break
    else:
        pytest.fail("band path contained no Γ point")


def test_mace_dynamically_stable(mace_outputs):
    _, band, summary = mace_outputs
    all_freqs = [b["frequency"] for p in band["phonon"] for b in p["band"]]
    assert min(all_freqs) > -0.05  # THz; no imaginary modes for stable fcc Cu
    assert summary["min_band_frequency_THz"] > -0.05
    assert summary["dynamically_stable_at_0K"] is True
