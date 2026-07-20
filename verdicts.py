#!/usr/bin/env python3
"""verdicts.py — milestone 3: assemble the verdicts.json sidecar.

Combines three same-pipeline measurements into per-irrep static/dynamic
verdicts (schema v0.1, docs/verdicts-schema.md):

    measured   per-config windowed irrep amplitudes of the RMC ensemble
    quantum    the same projections on quantum-thermal snapshots of the
               MLIP null model (zero-point correct)
    noise      the random-sign null of the measured configs

Verdict ratio (three-component expectation; see the design doc):

    r_m(w) = A²_meas(m,w) / [ A²_qh(m,w) + f_noise · A²_null(m,w) ]
    f_noise = (σ²_tot − σ²_qh) / σ²_tot        (measured, not assumed)

dynamic: r < 1.5;  static: r > 3;  else mixed. Confidence = bootstrap over
configs and snapshots. Optional E(Q) mapping fills omega_THz by projecting
each irrep pattern onto the null-model phonon eigenvectors over its full
k-star. Units: Å, THz, meV.

The heavy inputs (windowed npz, quantum npz, σ's) are produced by the run
scripts; this module holds the pure assembly logic so it is testable.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

R_DYNAMIC, R_STATIC = 1.5, 3.0
THZ_TO_MEV = 4.135667

K_INFO = {
    "X5": {"star": [[0, 1, 0], [1, 0, 0], [0, 0, 1]], "irrep": "X5"},
    "X3": {"star": [[0, 1, 0], [1, 0, 0], [0, 0, 1]], "irrep": "X3"},
    "W4": {"star": [[0.5, 1, 0]], "irrep": "W4"},
    "D": {"star": [[0, 0, 0.5]], "irrep": "Delta (SM: Delta3)"},
    "G1": {"star": [[0, 0, 0]], "irrep": "Gamma1"},
    "G3": {"star": [[0, 0, 0]], "irrep": "Gamma3"},
}


def classify(r: float) -> str:
    """Threshold verdict for an amplitude ratio."""
    if r < R_DYNAMIC:
        return "dynamic"
    if r > R_STATIC:
        return "static"
    return "mixed"


def noise_fraction(sigma_tot: float, sigma_qh: float) -> float:
    """RMC-noise share of the total displacement variance (per component)."""
    return max(sigma_tot**2 - sigma_qh**2, 0.0) / sigma_tot**2


def ratio(a2_meas: float, a2_qh: float, a2_null: float,
          f_noise: float) -> float:
    """Three-component amplitude ratio r (dimensionless)."""
    return a2_meas / max(a2_qh + f_noise * a2_null, 1e-12)


def assemble(win, qh, sigma_tot, sigma_qh, published, extra_provenance=None,
             w_verdict=4, scales=(2, 4, 8), n_boot=400, seed=0):
    """Build the verdicts dict from the measurement npz mappings.

    win / qh: mappings with keys 'keys' and 'rms_w{w}' (+ win: 'shuffle_w{w}')
    as produced by the run scripts. Returns the verdicts dict (schema v0.1).
    """
    keys = [str(k) for k in win["keys"]]
    assert keys == [str(k) for k in qh["keys"]]
    f_noise = noise_fraction(sigma_tot, sigma_qh)
    rng = np.random.default_rng(seed)

    def a2(arr, j):
        return float((arr[:, j]**2).mean())

    modes = []
    for j, k in enumerate(keys):
        a2m = a2(win[f"rms_w{w_verdict}"], j)
        a2n = a2(win[f"shuffle_w{w_verdict}"], j)
        a2q = a2(qh[f"rms_w{w_verdict}"], j)
        r = ratio(a2m, a2q, a2n, f_noise)
        verdict = classify(r)

        n, nq = len(win[f"rms_w{w_verdict}"]), len(qh[f"rms_w{w_verdict}"])
        agree = 0
        for _ in range(n_boot):
            i = rng.integers(0, n, n)
            iq = rng.integers(0, nq, nq)
            rb = ratio(float((win[f"rms_w{w_verdict}"][i, j]**2).mean()),
                       float((qh[f"rms_w{w_verdict}"][iq, j]**2).mean()),
                       a2n, f_noise)
            agree += classify(rb) == verdict
        modes.append({
            "q": K_INFO[k]["star"][0],
            "branch": None,
            "omega_THz": None,
            "amplitude_ratio": round(r, 3),
            "well": None,
            "barrier_meV": None,
            "verdict": verdict,
            "confidence": round(agree / n_boot, 3),
            "irrep": K_INFO[k]["irrep"],
            "star": K_INFO[k]["star"],
            "amplitude_measured_A": round(float(np.sqrt(a2m)), 4),
            "amplitude_quantum_A": round(float(np.sqrt(a2q)), 4),
            "amplitude_noise_A": round(float(np.sqrt(a2n)), 4),
            "static_amplitude_A": round(float(np.sqrt(max(
                a2m - a2q - f_noise * a2n, 0.0))), 4),
            "target_amplitude_A": published.get(k),
            "window_cells": w_verdict,
            "ratio_vs_quantum_by_scale": {
                f"w{w}": round(ratio(a2(win[f"rms_w{w}"], j),
                                     a2(qh[f"rms_w{w}"], j),
                                     a2(win[f"shuffle_w{w}"], j),
                                     f_noise), 3) for w in scales},
            "notes": ("Gamma channels are reference-dependent (measured "
                      "against the ensemble mean); interpret staggered "
                      "channels only." if k in ("G1", "G3") else ""),
        })

    prov = {
        "expectation": "A2_qh + f_noise*A2_null (three-component)",
        "f_noise": round(f_noise, 3),
        "sigma_tot_A": round(sigma_tot, 4),
        "sigma_quantum_A": round(sigma_qh, 4),
        "thresholds": {"dynamic_r_below": R_DYNAMIC,
                       "static_r_above": R_STATIC},
        "generated": time.strftime("%Y-%m-%d"),
    }
    prov.update(extra_provenance or {})
    return {"schema_version": "0.1", "material": "GaTa4Se8",
            "temperature_K": 5.0, "provenance": prov, "modes": modes}


def eq_map_pattern(phonon, pattern104, slab_positions_cart, masses,
                   arms_conv, b_conv):
    """Overlap of a 104-atom pattern with null-model phonons over a k-star.

    phonon: Phonopy object with force constants (2×2×2 suffices — X, W and
    Δ are all commensurate). Returns (freqs_THz, overlaps) summed over arms
    and normalised.
    """
    prim = phonon.primitive
    b_prim = 2 * np.pi * np.linalg.inv(prim.cell).T
    conv2prim = np.linalg.inv(b_prim) @ b_conv
    inv_prim = np.linalg.inv(prim.cell.T)
    total = None
    for arm in arms_conv:
        q_prim = conv2prim @ np.asarray(arm, dtype=float)
        phonon.run_qpoints([q_prim], with_eigenvectors=True)
        freqs = phonon.qpoints.frequencies[0]
        eigs = phonon.qpoints.eigenvectors[0]
        V = np.zeros((len(prim.positions), 3), dtype=complex)
        for i in range(len(pattern104)):
            r = slab_positions_cart[i]
            dfrac = inv_prim @ (r[:, None] - prim.positions.T)
            drem = dfrac - np.rint(dfrac)
            jat = int(np.argmin((drem**2).sum(axis=0)))
            R = r - prim.positions[jat]
            phase = np.exp(-2j * np.pi * (q_prim @ (inv_prim @ R)))
            V[jat] += np.sqrt(masses[i]) * pattern104[i] * phase
        v = V.ravel()
        n2 = float(np.vdot(v, v).real)
        if n2 < 1e-12:
            continue
        o = np.abs(v.conj() @ eigs)**2 / n2
        total = o if total is None else total + o
    total = total / total.sum()
    return freqs, total
