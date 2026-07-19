"""Unit tests for md_run's pair-distribution / closure machinery.

Analytic known answers: fcc first-shell geometry, Faber–Ziman weight closure,
sine-FT round trip, and the scale+offset Rw fit. All pure numpy, no MLIP.
"""

from __future__ import annotations

import numpy as np
import pytest

import md_run as m2


# ---------------------------------------------------------------- weights

def test_weights_monatomic_sum_to_one():
    w, species = m2.neutron_weights(["Cu"] * 8)
    assert species == ["Cu"]
    assert w[("Cu", "Cu")] == pytest.approx(1.0)


def test_weights_binary_sum_to_one():
    symbols = ["Ga"] * 4 + ["Ta"] * 16 + ["Se"] * 32
    w, species = m2.neutron_weights(symbols)
    assert species == ["Ga", "Se", "Ta"]
    assert sum(w.values()) == pytest.approx(1.0)
    # cross pairs carry the factor 2: w_ab / (c_a c_b b_a b_b) is 2x diagonal
    c = {s: symbols.count(s) / len(symbols) for s in species}
    b = m2.BCOH_FM
    ratio_diag = w[("Se", "Se")] / (c["Se"]**2 * b["Se"]**2)
    ratio_cross = w[("Ga", "Ta")] / (c["Ga"] * c["Ta"] * b["Ga"] * b["Ta"])
    assert ratio_cross == pytest.approx(2 * ratio_diag)


def test_weights_unknown_element_raises():
    with pytest.raises(SystemExit, match="b_coh"):
        m2.neutron_weights(["Xx"])


# ---------------------------------------------------------------- g(r)

def _fcc_positions(a, n):
    """Ideal fcc conventional supercell positions (Cartesian, Å)."""
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    pos = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for b in basis:
                    pos.append((b + [i, j, k]) * a)
    return np.array(pos)


def test_fcc_first_shell_position_and_coordination():
    a, n = 3.61, 3
    pos = _fcc_positions(a, n)
    cell = np.diag([n * a] * 3)
    r, g = m2.pair_histograms([pos], ["Cu"] * len(pos), cell,
                              r_max=5.0, dr=0.02)
    gcc = g[("Cu", "Cu")]
    # first peak at a/sqrt(2)
    first = r[np.argmax(gcc)]
    assert first == pytest.approx(a / np.sqrt(2), abs=0.02)
    # coordination number: rho * int 4 pi r^2 g dr over the first shell = 12
    rho = len(pos) / cell.diagonal().prod()
    mask = (r > first - 0.1) & (r < first + 0.1)
    cn = rho * (4 * np.pi * r[mask]**2 * gcc[mask] * 0.02).sum()
    assert cn == pytest.approx(12.0, rel=0.02)


def test_gr_normalisation_identity():
    """rho * int 4 pi r^2 g(r) dr over [0, R] equals the exact average
    neighbour count within R — the ideal-gas normalisation identity.
    (A crystal's g(r) does NOT flatten to 1 at accessible r; the integral is
    the correct invariant.)"""
    rng = np.random.default_rng(0)
    a, n = 3.61, 4
    pos = _fcc_positions(a, n)
    cell = np.diag([n * a] * 3)
    L = np.diag(cell)
    # exact per-atom neighbour count within R on the ideal lattice (MIC)
    R = 7.0
    d = pos[:, None, :] - pos[None, :, :]
    d -= np.round(d / L) * L
    rr = np.sqrt((d**2).sum(-1))
    exact = ((rr > 1e-8) & (rr < R)).sum() / len(pos)
    # histogram from smeared frames (smearing conserves counts)
    frames = [pos + rng.normal(0, 0.1, pos.shape) for _ in range(4)]
    r, g = m2.pair_histograms(frames, ["Cu"] * len(pos), cell,
                              r_max=R, dr=0.05)
    rho = len(pos) / L.prod()
    integral = rho * (4 * np.pi * r**2 * g[("Cu", "Cu")] * 0.05).sum()
    assert integral == pytest.approx(exact, rel=0.02)


def test_rmax_beyond_half_box_raises():
    pos = _fcc_positions(3.61, 2)
    with pytest.raises(SystemExit, match="half the box"):
        m2.pair_histograms([pos], ["Cu"] * len(pos), np.diag([7.22] * 3),
                           r_max=4.0, dr=0.02)


# ---------------------------------------------------------------- FT

def test_ft_roundtrip_gaussian():
    """G -> F -> G reproduces a smooth Gaussian well inside the windows."""
    rho0 = 0.05
    r = np.arange(0.005, 30.0, 0.01)
    G = 0.4 * np.exp(-((r - 5.0) ** 2) / (2 * 0.4**2))
    Q = np.arange(0.02, 40.0, 0.01)
    F = m2.gr_to_fq(r, G, Q, rho0)
    G2 = m2.fq_to_gr(Q, F, r, rho0)
    core = (r > 3.0) & (r < 7.0)
    assert np.max(np.abs(G2[core] - G[core])) < 0.02


def test_fq_static_lattice_peaks_at_bragg_positions():
    """A static fcc lattice's F(Q) peaks at the (111) Bragg position."""
    a, n = 3.61, 4
    pos = _fcc_positions(a, n)
    cell = np.diag([n * a] * 3)
    r, g = m2.pair_histograms([pos], ["Cu"] * len(pos), cell,
                              r_max=7.0, dr=0.02)
    G = m2.total_G(r, g, ["Cu"] * len(pos))
    Q = np.arange(2.0, 6.0, 0.01)
    F = m2.gr_to_fq(r, G, Q, len(pos) / cell.diagonal().prod())
    q111 = 2 * np.pi * np.sqrt(3) / a  # 3.014 A^-1 for a=3.61
    qpk = Q[np.argmax(F)]
    assert qpk == pytest.approx(q111, abs=0.15)  # truncation-broadened


# ---------------------------------------------------------------- closure fit

def test_fit_scale_offset_recovers_parameters():
    rng = np.random.default_rng(1)
    sim = np.sin(np.linspace(0, 20, 500)) * np.exp(-np.linspace(0, 2, 500))
    data = 2.0 * sim + 0.3 + rng.normal(0, 1e-4, sim.size)
    s, o, rw = m2.fit_scale_offset(data, sim)
    assert s == pytest.approx(2.0, abs=1e-3)
    assert o == pytest.approx(0.3, abs=1e-3)
    assert rw < 1e-3


def test_parse_fq(tmp_path):
    p = tmp_path / "t.fq"
    p.write_text("        3\nrmc S(Q)#L\n  0.80  -0.87\n  0.81  -0.86\n"
                 "  0.82  -0.85\n")
    Q, F = m2.parse_fq(p)
    np.testing.assert_allclose(Q, [0.80, 0.81, 0.82])
    np.testing.assert_allclose(F, [-0.87, -0.86, -0.85])
