"""Tests for mode_project.py — the M3 projection engine.

End-to-end validations against the published GaTa4Se8 numbers:
  1. the label->orbit mapping + pattern expansion reproduces the refined
     total distortion (implied amplitudes == published Table II);
  2. a synthetic box carrying the full refined distortion field projects to
     the published amplitudes per irrep (validates star pooling + Gram);
  3. an injected single mode along a DIFFERENT cubic axis is recovered at
     the right amplitude in the right channel (powder-degeneracy pooling).

Requires reference/gts_mode_patterns.json and data/GTS_5K.cif.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import mode_project as mp

REPO = Path(__file__).resolve().parent.parent
CIF = REPO / "data/GTS_5K.cif"
pytestmark = pytest.mark.skipif(
    not (mp.REF_JSON.is_file() and CIF.is_file()),
    reason="needs gts_mode_patterns.json and data/GTS_5K.cif")

A_CUB = 10.3563


def _cubic_sites():
    from ase.io import read

    at = read(str(CIF))
    return at.get_scaled_positions() % 1.0, at.get_chemical_symbols()


@pytest.fixture(scope="module")
def setup():
    ref = mp.load_reference()
    a = ref["cell"]["a_A"]
    c = ref["cell"]["c_A"]
    frac52_raw, elem52 = _cubic_sites()
    frac52_ideal = mp.idealize_parent(frac52_raw, elem52)
    refined, refined_elem = mp.expand_refined(ref, a, c)
    frac52, _op, worst = mp.align_parent_frame(frac52_ideal, elem52, refined,
                                               refined_elem, a, c)
    assert worst < 0.15, worst          # physical distortion scale only
    slab = mp.build_slab(frac52)
    slab_elem = list(elem52) * 2
    D = mp.displacement_field(slab, slab_elem, refined, refined_elem, a, c)
    mapping, orbits = mp.map_labels_to_orbits(ref, slab, slab_elem, D, a, c)
    fields = mp.expand_patterns(ref, slab, slab_elem, mapping, orbits, a, c)
    return ref, frac52, elem52, slab, slab_elem, D, fields


def _iv_field(setup):
    """The printed Table-IV total distortion expanded to 104 atoms (Å).

    Used as the validation target instead of the geometric D field: D carries
    a Γ-space offset from the difference between OUR cubic parent (RMC
    average) and the AMPLIMODES cubic reference (their refinement), which the
    Γ₁/Γ₃ patterns alone cannot absorb. The staggered physics channels
    (X₅/X₃/W₄/Δ) are reference-independent, so this distinction only matters
    for Γ-channel verdicts (documented in the design doc).
    """
    ref, frac52, elem52, slab, slab_elem, D, fields = setup
    a, c = ref["cell"]["a_A"], ref["cell"]["c_A"]
    mapping, orbits = mp.map_labels_to_orbits(ref, slab, slab_elem, D, a, c)
    F = np.zeros((len(slab), 3))
    for lab, (oi, anchor) in mapping.items():
        u = np.array(ref["table_IV_total"][lab]["u_frac"]) * \
            np.array([a, a, c])
        orb = orbits[oi]
        F[orb["members"]] = mp.expand_vec_over_orbit(orb, anchor, u)
    return F


def test_expansion_reproduces_published_amplitudes(setup):
    """LSQ of the expanded Table-IV field onto the six expanded mode fields
    must give the published amplitudes (parent-primitive convention)."""
    ref, *_, fields = setup
    T = _iv_field(setup)
    keys = sorted(fields)
    A = np.column_stack([fields[k].ravel() for k in keys])
    coef, *_ = np.linalg.lstsq(A, T.ravel(), rcond=None)
    resid = np.linalg.norm(T.ravel() - A @ coef) / np.linalg.norm(T)
    assert resid < 0.06, resid
    pub = ref["published_amplitudes_A"]
    for j, k in enumerate(keys):
        implied = float(coef[j] * np.linalg.norm(A[:, j]) / np.sqrt(8))
        if pub[k] >= 0.02:                       # major modes
            assert implied == pytest.approx(pub[k], rel=0.10), (k, implied)


def _box_S_from_compact(compacts, coeffs):
    """Parity sums S for a synthetic box whose displacement field is
    sum_i coeffs[i] * compact_i (each compact = (axis, G[52,2,3]))."""
    S = np.zeros((3, 2, 52, 3))
    for (axF, G), cf in zip(compacts, coeffs):
        for px in range(2):
            for py in range(2):
                for pz in range(2):
                    par = (px, py, pz)
                    for ax in range(3):
                        S[ax, par[ax]] += 64.0 * cf * G[:, par[axF]]
    return S


@pytest.fixture(scope="module")
def projector(setup):
    ref, frac52, elem52, *_ , fields = setup
    return mp.build_projector(fields, frac52, elem52)


def test_full_distortion_box_projects_to_published(setup, projector):
    """A box tiled with the complete published distortion field reads the
    published amplitude in every major channel."""
    ref, *_ = setup
    T = _iv_field(setup)
    S = _box_S_from_compact([mp.field_to_compact(T)], [1.0])
    got = mp.project_all(S, projector)
    pub = ref["published_amplitudes_A"]
    for key in ("X5", "X3", "W4", "D"):
        assert got[key] == pytest.approx(pub[key], rel=0.15), (key, got[key])


def test_injected_mode_recovered_across_axes(setup, projector):
    """Inject X5 at a known amplitude with the doubling axis rotated to x:
    star pooling must recover the amplitude in the X5 channel only."""
    i0, i1 = projector["blocks"]["X5"]
    var_x = next(v for v in projector["variants"][i0:i1] if v[0] == 0)
    norm2 = mp.compact_inner(var_x, var_x)          # box field norm^2
    target = 0.10                                    # Å, published convention
    lam = target * np.sqrt(2048.0 / norm2)
    S = _box_S_from_compact([var_x], [lam])
    got = mp.project_all(S, projector)
    assert got["X5"] == pytest.approx(target, rel=0.05)
    for key in ("X3", "W4", "D"):
        assert got[key] < 0.15 * target, (key, got[key])


def test_variant_ranks_positive(projector):
    assert projector["rank"] >= 6
    assert len(projector["variants"]) >= projector["rank"]
