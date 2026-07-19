"""Verification of the ingested GaTa4Se8 supplemental mode patterns.

The strong test: Table IV (total normalized displacement at 10 K) must be a
linear combination of the six irrep patterns (Tables V–X), and the fitted
coefficients must be proportional to the published mode amplitudes
(main-text Table II of PRR 4, 033123). This validates the PDF extraction,
the internal consistency of the SM, and pins the empirical normalization
convention — without assuming it.

Skipped when reference/gts_mode_patterns.json has not been generated
(requires the git-ignored data/GaTa4Se8_SM.pdf).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

REF = Path(__file__).resolve().parent.parent / "reference/gts_mode_patterns.json"
pytestmark = pytest.mark.skipif(not REF.is_file(),
                                reason="gts_mode_patterns.json not generated")


@pytest.fixture(scope="module")
def ref():
    return json.loads(REF.read_text())


def _cart(frac, cell):
    """Fractional (x/a, y/b, z/c) triple -> Cartesian Å for the tetragonal cell."""
    a, c = cell["a_A"], cell["c_A"]
    f = np.asarray(frac, dtype=float)
    return f * np.array([a, a, c])


MULT = {"4e": 4, "8f": 8}


def test_counts_and_multiplicities(ref):
    for name in ("table_II_refined", "table_III_pdf_refined", "table_IV_total"):
        assert len(ref[name]) == 20, name
    assert sum(MULT[r["wp"]] for r in ref["table_IV_total"].values()) == 104
    for key, mode in ref["modes"].items():
        assert len(mode["pattern_frac"]) == 20, key
        # identical labeling family as Table IV
        assert set(mode["pattern_frac"]) == set(ref["table_IV_total"]), key


def test_table_iv_magnitudes_consistent(ref):
    """|u| column must equal the fractional components in Å (rounding tol)."""
    for lab, row in ref["table_IV_total"].items():
        u = np.linalg.norm(_cart(row["u_frac"], ref["cell"]))
        assert u == pytest.approx(row["u_abs_A"], abs=0.004), (lab, u)


def test_table_iv_decomposes_onto_patterns(ref):
    """Multiplicity-weighted LSQ: IV = sum_m c_m * pattern_m, small residual,
    and c_m ratios proportional to the published amplitudes."""
    labels = sorted(ref["table_IV_total"])
    w = np.sqrt([MULT[ref["table_IV_total"][l]["wp"]] for l in labels])
    target = np.concatenate([
        w[i] * _cart(ref["table_IV_total"][l]["u_frac"], ref["cell"])
        for i, l in enumerate(labels)])
    keys = sorted(ref["modes"])
    A = np.column_stack([
        np.concatenate([w[i] * _cart(ref["modes"][k]["pattern_frac"][l]["d_frac"],
                                     ref["cell"])
                        for i, l in enumerate(labels)])
        for k in keys])
    c, res, *_ = np.linalg.lstsq(A, target, rcond=None)
    fit = A @ c
    rel_resid = np.linalg.norm(target - fit) / np.linalg.norm(target)
    assert rel_resid < 0.10, rel_resid  # rounding-level misfit only
    # every site at rounding level (the Se3 sign correction is validated here)
    site_resid = np.linalg.norm((target - fit).reshape(len(labels), 3), axis=1)
    assert site_resid.max() < 0.015, dict(zip(labels, site_resid.round(4)))

    # Amplitude convention (empirically established): the published
    # amplitudes are the mode-field norms over the PARENT PRIMITIVE cell,
    # i.e. A_m = c_m * ||P_m||_cell / sqrt(8)  (104-atom cell = 8 primitives).
    pub = ref["published_amplitudes_A"]
    implied = {}
    for j, k in enumerate(keys):
        implied[k] = float(c[j] * np.linalg.norm(A[:, j]) / np.sqrt(8))
    for k in ("X5", "X3", "W4", "D"):   # major modes; Gammas are at the
        assert implied[k] == pytest.approx(pub[k], rel=0.20), (  # rounding floor
            k, implied[k], pub[k])
    print("\nimplied amplitudes (A, parent-primitive convention):",
          {k: round(v, 4) for k, v in implied.items()},
          "vs published", pub)


def test_patterns_mutually_orthogonal(ref):
    """Different irreps are orthogonal under the multiplicity-weighted
    Cartesian inner product (up to table rounding)."""
    labels = sorted(ref["table_IV_total"])
    mult = np.array([MULT[ref["table_IV_total"][l]["wp"]] for l in labels])
    vecs = {}
    for k, mode in ref["modes"].items():
        v = np.concatenate([
            np.sqrt(mult[i]) * _cart(mode["pattern_frac"][l]["d_frac"],
                                     ref["cell"])
            for i, l in enumerate(labels)])
        vecs[k] = v / np.linalg.norm(v)
    keys = sorted(vecs)
    overlaps = {}
    for i, k1 in enumerate(keys):
        for k2 in keys[i + 1:]:
            dot = abs(float(vecs[k1] @ vecs[k2]))
            overlaps[(k1, k2)] = round(dot, 3)
            # Exact irrep fields are orthogonal; the printed tables carry only
            # ~2 significant figures, which limits achievable orthogonality to
            # the ~0.1-0.2 level (worst: G1-X3). mode_project.py therefore
            # re-orthogonalises the basis before projecting; this test only
            # guards against gross (wrong-table/sign-family) errors.
            assert dot < 0.20, (k1, k2, dot)
    print("\npattern overlaps (rounding-limited):", overlaps)
