"""End-to-end integration test of milestone1_bands with the EMT calculator.

EMT is ASE's dependency-free metals-only toy potential (Cu, Al, Ag, Au, Ni,
Pd, Pt), so the smoke-test fixture is fcc Cu. This exercises the whole
pipeline against the *installed* ase/spglib/phonopy — the leg CLAUDE.md flags
as previously unrun — and asserts the fcc known answer:

    3 acoustic branches (monatomic primitive -> 3 bands), all -> 0 at Γ,
    and no imaginary modes (dynamically stable).

Marked ``slow``: it runs a real relaxation + phonopy finite-displacement
calculation. Deselect with ``pytest -m "not slow"``.
"""

from __future__ import annotations

import numpy as np
import pytest

from fixtures.make_synthetic_ensemble import make_fcc_cu_ensemble

pytestmark = pytest.mark.slow


def _gamma_frequencies(phonon_points):
    """Frequencies (THz) at every q = Γ found in a band.yaml phonon list."""
    gammas = []
    for p in phonon_points:
        if np.allclose(p["q-position"], 0.0, atol=1e-9):
            gammas.append([b["frequency"] for b in p["band"]])
    return gammas


@pytest.fixture(scope="module")
def emt_outputs(tmp_path_factory, m1):
    """Run the full pipeline once (module-scoped) and return the output dir."""
    import yaml

    root = tmp_path_factory.mktemp("emt")
    ens = root / "ens"
    out = root / "out"
    make_fcc_cu_ensemble(ens, n_configs=4, dims=(2, 2, 2), sigma=0.06, seed=0)

    # --dim 2 2 2 keeps the phonopy supercell small (fast) while still
    # resolving the acoustic branches; --no-eigenvectors shrinks band.yaml.
    m1.main([
        str(ens), "--calc", "emt", "-o", str(out),
        "--dim", "2", "2", "2", "--npoints", "21", "--no-eigenvectors",
    ])

    band = yaml.safe_load((out / "band.yaml").read_text())
    import json
    summary = json.loads((out / "summary.json").read_text())
    return out, band, summary


def test_emt_writes_all_outputs(emt_outputs):
    out, _, _ = emt_outputs
    for name in ("band.yaml", "relaxed.cif", "summary.json"):
        assert (out / name).is_file(), f"missing {name}"


def test_emt_relaxed_cif_is_fcc_copper(emt_outputs):
    from ase.io import read as ase_read

    out, _, _ = emt_outputs
    atoms = ase_read(str(out / "relaxed.cif"))
    assert set(atoms.get_chemical_symbols()) == {"Cu"}
    # Standardized fcc conventional cell: cubic, ~3.6 Å edges.
    angles = atoms.cell.angles()
    np.testing.assert_allclose(angles, [90, 90, 90], atol=1.0)
    np.testing.assert_allclose(atoms.cell.lengths(), atoms.cell.lengths()[0],
                               atol=1e-2)


def test_emt_three_acoustic_branches_vanish_at_gamma(emt_outputs):
    _, band, _ = emt_outputs
    # Monatomic fcc primitive -> exactly 3 phonon branches.
    assert band["natom"] == 1
    assert len(band["phonon"][0]["band"]) == 3

    gammas = _gamma_frequencies(band["phonon"])
    assert gammas, "band path contained no Γ point"
    for freqs in gammas:
        # All three acoustic branches touch zero at Γ.
        assert np.allclose(freqs, 0.0, atol=2e-2), freqs


def test_emt_dynamically_stable(emt_outputs):
    _, band, summary = emt_outputs
    all_freqs = [b["frequency"] for p in band["phonon"] for b in p["band"]]
    assert min(all_freqs) > -0.05  # THz; no imaginary modes
    assert summary["dynamically_stable_at_0K"] is True
    assert summary["min_band_frequency_THz"] > -0.05
