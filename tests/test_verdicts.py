"""Unit tests for verdicts.py — the assembly logic on synthetic inputs.

Constructs windowed/quantum npz-like mappings with known amplitude
relationships and checks the three-component ratios, verdict thresholds,
confidence behaviour, and schema shape.
"""

from __future__ import annotations

import numpy as np
import pytest

import verdicts as vd


def _mk(keys, amp, n=200, jitter=0.02, seed=0):
    """Mapping like the run-script npz: rms_w{2,4,8} (+ shuffle) columns."""
    rng = np.random.default_rng(seed)
    out = {"keys": np.array(keys)}
    for w in (2, 4, 8):
        base = np.array([amp[k][f"w{w}"] for k in keys])
        out[f"rms_w{w}"] = base * (1 + jitter * rng.standard_normal(
            (n, len(keys))))
        if f"shuffle_w{w}" not in out:
            noise = np.array([amp[k].get("noise", 0.01) for k in keys])
            out[f"shuffle_w{w}"] = noise * np.ones((60, len(keys)))
    return out


def test_classify_thresholds():
    assert vd.classify(0.5) == "dynamic"
    assert vd.classify(1.49) == "dynamic"
    assert vd.classify(2.0) == "mixed"
    assert vd.classify(3.01) == "static"


def test_noise_fraction():
    assert vd.noise_fraction(0.1, 0.1) == 0.0
    assert vd.noise_fraction(0.1, 0.0) == 1.0
    assert vd.noise_fraction(0.0954, 0.0400) == pytest.approx(0.824, abs=0.001)


def test_assemble_verdicts_and_schema():
    keys = ["X5", "X3"]
    # X5: measured far above quantum+noise -> static; X3: equal -> dynamic
    quantum = {k: {"w2": 0.05, "w4": 0.03, "w8": 0.02} for k in keys}
    measured = {
        "X5": {"w2": 0.20, "w4": 0.15, "w8": 0.10, "noise": 0.02},
        "X3": {"w2": 0.055, "w4": 0.033, "w8": 0.021, "noise": 0.02},
    }
    win = _mk(keys, measured)
    qh = _mk(keys, quantum, n=64)
    # sigma_tot == sigma_qh -> f_noise = 0 -> pure measured/quantum ratios
    v = vd.assemble(win, qh, sigma_tot=0.05, sigma_qh=0.05,
                    published={"X5": 0.12, "X3": 0.07})
    assert v["schema_version"] == "0.1"
    m = {x["irrep"]: x for x in v["modes"]}
    assert m["X5"]["verdict"] == "static"
    assert m["X5"]["confidence"] > 0.9
    assert m["X3"]["verdict"] == "dynamic"
    # required contract fields present on every mode
    for x in v["modes"]:
        for f in ("q", "branch", "omega_THz", "amplitude_ratio", "well",
                  "barrier_meV", "verdict", "confidence"):
            assert f in x
    # ratio-by-scale carried and ordered
    assert set(m["X5"]["ratio_vs_quantum_by_scale"]) == {"w2", "w4", "w8"}
    # amplitude ratio ~ (0.15/0.03)^2 = 25 at w4 with f_noise = 0
    assert m["X5"]["amplitude_ratio"] == pytest.approx(25.0, rel=0.1)


def test_noise_component_lowers_ratio():
    keys = ["X5"]
    measured = {"X5": {"w2": 0.2, "w4": 0.15, "w8": 0.1, "noise": 0.12}}
    quantum = {"X5": {"w2": 0.05, "w4": 0.03, "w8": 0.02}}
    win = _mk(keys, measured)
    qh = _mk(keys, quantum, n=64)
    # large noise fraction: expectation gains f*A2_null and r drops
    v0 = vd.assemble(win, qh, 0.05, 0.05, {})     # f = 0
    v1 = vd.assemble(win, qh, 0.10, 0.04, {})     # f = 0.84
    r0 = v0["modes"][0]["amplitude_ratio"]
    r1 = v1["modes"][0]["amplitude_ratio"]
    assert r1 < r0
    f = vd.noise_fraction(0.10, 0.04)
    expect = 0.15**2 / (0.03**2 + f * 0.12**2)
    assert r1 == pytest.approx(expect, rel=0.1)
