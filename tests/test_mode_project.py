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


# ---------------------------------------------------------------- windows

def _synthetic_box(setup, field_fn):
    """Explicit (unit_frac, sid, ijk) arrays for an 8×8×8 box whose
    displacement field is field_fn(site_index, i, j, k) -> (3,) Å."""
    ref, frac52, *_ = setup
    a = A_CUB
    sid, ijk, unit = [], [], []
    for i in range(8):
        for j in range(8):
            for k in range(8):
                for s in range(52):
                    sid.append(s + 1)
                    ijk.append((i, j, k))
                    unit.append(frac52[s] + field_fn(s, i, j, k) / a)
    return (np.array(unit) % 1.0, np.array(sid), np.array(ijk, dtype=int),
            frac52)


def _x5_variant_z(projector, target=0.10):
    """An X5 variant with doubling axis z, scaled to `target` Å (published
    convention)."""
    i0, i1 = projector["blocks"]["X5"]
    var = next(v for v in projector["variants"][i0:i1] if v[0] == 2)
    lam = target * np.sqrt(2048.0 / mp.compact_inner(var, var))
    return var, lam


def test_windowed_matches_global_for_coherent_field(setup, projector):
    """A fully coherent X5 box reads the same amplitude at every window
    scale (normalization is scale-invariant)."""
    var, lam = _x5_variant_z(projector, target=0.10)
    ax, G = var

    unit, sid, ijk, frac52 = _synthetic_box(
        setup, lambda s, i, j, k: lam * G[s, (i, j, k)[ax] & 1])
    for w in (2, 4, 8):
        S_w = mp.window_parity_sums(unit, sid, ijk, frac52, A_CUB, w)
        out = mp.project_all_windows(S_w, projector, w)
        assert out["X5"].mean() == pytest.approx(0.10, rel=0.05), (w,
                                                                   out["X5"])


def test_domain_structure_detected_by_windows(setup, projector):
    """Two anti-phase X5 domains (z-halves of the box): the global
    projection cancels, matched windows read the full amplitude."""
    var, lam = _x5_variant_z(projector, target=0.10)
    ax, G = var

    def field(s, i, j, k):
        sign = 1.0 if k < 4 else -1.0
        return sign * lam * G[s, (i, j, k)[ax] & 1]

    unit, sid, ijk, frac52 = _synthetic_box(setup, field)
    S = mp.parity_sums(unit, sid, ijk, frac52, A_CUB)
    global_amp = mp.project_all(S, projector)["X5"]
    assert global_amp < 0.02          # anti-phase halves cancel

    S_w = mp.window_parity_sums(unit, sid, ijk, frac52, A_CUB, 4)
    out = mp.project_all_windows(S_w, projector, 4)
    # every 4-cell window lies inside one domain -> full local amplitude
    assert out["X5"].mean() == pytest.approx(0.10, rel=0.06), out["X5"]


def test_noise_nulls_have_the_right_selectivity(setup, projector):
    """The two nulls behave as designed on a coherent X5 box:

    - random SIGN flip kills every coherent channel (the universal noise
      floor): X5 collapses to the pedestal;
    - cell SHUFFLE preserves X5 exactly (an X-point pattern is constant
      across cells — it is a coherence *diagnostic*, not a noise floor).
    """
    var, lam = _x5_variant_z(projector, target=0.10)
    ax, G = var
    unit, sid, ijk, frac52 = _synthetic_box(
        setup, lambda s, i, j, k: lam * G[s, (i, j, k)[ax] & 1])
    rng = np.random.default_rng(0)

    # universal null: random signs -> pedestal (~target/sqrt(512), pooled)
    signs = mp.random_signs(len(sid), rng)
    S = mp.parity_sums(unit, sid, ijk, frac52, A_CUB, signs=signs)
    killed = mp.project_all(S, projector)["X5"]
    assert killed < 0.04, killed

    # diagnostic: cell shuffle preserves the X-point (intra-cell) pattern
    ijk_sh = mp.shuffle_cells(unit, sid, ijk, rng)
    S2 = mp.parity_sums(unit, sid, ijk_sh, frac52, A_CUB)
    kept = mp.project_all(S2, projector)["X5"]
    assert kept == pytest.approx(0.10, rel=0.05), kept
