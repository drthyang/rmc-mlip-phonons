#!/usr/bin/env python3
"""export_modes.py — milestone 3: emit the irrep modes in phonopy's format.

Writes `modes_irrep.yaml` in the band.yaml schema (nqpoint / natom /
lattice / points / phonon → band → frequency + eigenvector), so viewer/
and any phonopy-aware tool can animate the atomic motion of each
distortion mode directly.

Each published irrep pattern (X5, X3, W4, Δ, Γ1, Γ3 — plus TOTAL, the full
refined distortion) is converted to a proper Bloch eigenvector at the arm
of its k-star that carries the pattern:

    e_κ(q) ∝ √m_κ · Σ_cells u(cell, κ) · exp(−i q·R)        (unit norm)

which is exactly phonopy's eigenvector convention: a viewer reconstructs
the real-space motion as Re[e_κ/√m_κ · exp(i(q·r − ωt))] and recovers the
staggered pattern. Frequencies are the E(Q)-weighted null-model values
from results/verdicts.json when present (omega_source mace-null), else 0.

Also writes per-mode 40-frame extended-XYZ animations of the 104-atom cell
(`<mode>.xyz`, displacement exaggerated for visibility) for OVITO/ASE.
Units: Å, THz; masses amu.

Usage:
    python export_modes.py -o results/modes    # needs data/ + reference/
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import milestone1_bands as m1
import mode_project as mp

REPO = Path(__file__).resolve().parent
A_CUB = 10.3563
STAR_ARMS = {
    "X5": [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    "X3": [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    "W4": [(0.5, 1, 0), (1, 0.5, 0), (0, 1, 0.5), (1, 0, 0.5),
           (0.5, 0, 1), (0, 0.5, 1)],
    "D": [(0, 0, 0.5), (0, 0.5, 0), (0.5, 0, 0)],
    "G1": [(0, 0, 0)],
    "G3": [(0, 0, 0)],
    "TOTAL": [(0, 0, 0.5), (0, 1, 0), (0.5, 1, 0), (0, 0, 0)],
}


def build_patterns(ensemble_dir: Path, ref_config: Path):
    """Reconstruct the 104-atom pattern fields in the RMC frame.

    Returns (fields {key: (104,3) Å}, slab_cart (104,3) Å, elements,
    primitive-cell info via a phonopy-free description: the cubic cell and
    site list). Uses the ensemble circular mean as the parent.
    """
    cfgs = [m1.parse_rmc6f(ref_config)]
    _, elems, frac, _ = m1.fold_and_average(cfgs, None)
    mean_site = frac                      # single-config circular mean
    elem52 = list(elems)

    ref = mp.load_reference()
    a_t, c_t = ref["cell"]["a_A"], ref["cell"]["c_A"]
    parent_ideal = mp.idealize_parent(mean_site, elem52)
    refined, refined_elem = mp.expand_refined(ref, a_t, c_t)
    aligned, _, _ = mp.align_parent_frame(parent_ideal, elem52, refined,
                                          refined_elem, a_t, c_t)
    slab = mp.build_slab(aligned)
    slab_elem = list(elem52) * 2
    D = mp.displacement_field(slab, slab_elem, refined, refined_elem,
                              a_t, c_t)
    mapping, orbits = mp.map_labels_to_orbits(ref, slab, slab_elem, D,
                                              a_t, c_t)
    fields = mp.expand_patterns(ref, slab, slab_elem, mapping, orbits,
                                a_t, c_t)
    fields["TOTAL"] = D
    # cartesian in the CUBIC frame: (x, y, z + p) * a
    slab_cub = np.array([[*aligned[i % 52][:2], aligned[i % 52][2] + i // 52]
                         for i in range(104)])
    return fields, slab_cub * A_CUB, slab_elem, aligned, elem52


def bloch_eigenvector(F104, pos_cart, masses, q_conv, a=A_CUB):
    """Pattern → unit-norm complex eigenvector on the 52-site cubic cell.

    q_conv in conventional-cubic reciprocal units; the phase uses the
    conventional-cell decomposition r = r_site + R (R the cubic lattice
    vector, here only (0,0,p·a) within the slab). Returns (52, 3) complex
    and the projected weight (fraction of the pattern captured at this q).
    """
    V = np.zeros((52, 3), dtype=complex)
    q = np.asarray(q_conv, dtype=float)
    for i in range(len(F104)):
        s = i % 52
        R = pos_cart[i] - pos_cart[s]          # (0,0,p·a)
        phase = np.exp(-2j * np.pi * (q @ (R / a)))
        V[s] += np.sqrt(masses[s]) * F104[i] * phase
    n2 = float(np.vdot(V.ravel(), V.ravel()).real)
    tot = sum(masses[i % 52] * float(F104[i] @ F104[i])
              for i in range(len(F104)))
    if n2 < 1e-14:
        return None, 0.0
    return V / np.sqrt(n2), n2 / (2.0 * tot)   # slab counts each site twice


def emit_yaml(entries, elem52, frac52, a, path: Path):
    """Write the band.yaml-schema file (plain text, no yaml lib needed).

    `a` is either a scalar (cubic lattice a·I) or a full (3,3) lattice."""
    from ase.data import atomic_masses, atomic_numbers

    lattice = np.eye(3) * a if np.isscalar(a) else np.asarray(a)
    L = []
    L.append(f"nqpoint: {len(entries)}")
    L.append(f"natom: {len(elem52)}")
    L.append("lattice:")
    for row in lattice:
        L.append(f"- [ {row[0]:20.15f}, {row[1]:20.15f}, {row[2]:20.15f} ]")
    L.append("points:")
    for el, x in zip(elem52, frac52):
        L.append(f"- symbol: {el}")
        L.append(f"  coordinates: [ {x[0]:18.15f}, {x[1]:18.15f}, "
                 f"{x[2]:18.15f} ]")
        L.append(f"  mass: {atomic_masses[atomic_numbers[el]]:.6f}")
    L.append("phonon:")
    for e in entries:
        q = e["q"]
        L.append(f"- q-position: [ {q[0]:12.7f}, {q[1]:12.7f}, "
                 f"{q[2]:12.7f} ]")
        L.append(f"  distance: {e['distance']:.7f}")
        L.append(f"  label: '{e['label']}'")
        L.append("  band:")
        L.append("  - # 1")
        L.append(f"    frequency: {e['frequency']:.10f}")
        L.append("    eigenvector:")
        for at in e["eigenvector"]:
            L.append("    - # atom")
            for comp in at:
                L.append(f"      - [ {comp.real:17.14f}, "
                         f"{comp.imag:17.14f} ]")
    path.write_text("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ensemble", type=Path,
                    default=REPO / "data/ensemble_20A_5K")
    ap.add_argument("--config", type=Path, default=None,
                    help="rmc6f used for the parent average "
                    "(default: GTS_5K_1.rmc6f in the ensemble dir)")
    ap.add_argument("--verdicts", type=Path,
                    default=REPO / "results/verdicts.json")
    ap.add_argument("--exaggeration", type=float, default=12.0)
    ap.add_argument("--nframes", type=int, default=40)
    ap.add_argument("-o", "--outdir", type=Path,
                    default=REPO / "results/modes")
    args = ap.parse_args(argv)

    cfg = args.config or (args.ensemble / "GTS_5K_1.rmc6f")
    fields, pos_cart, slab_elem, aligned52, elem52 = build_patterns(
        args.ensemble, cfg)

    from ase.data import atomic_masses, atomic_numbers  # noqa: F401 (main scope)
    masses52 = np.array([atomic_masses[atomic_numbers[e]] for e in elem52])

    freqs = {}
    if args.verdicts.is_file():
        v = json.loads(args.verdicts.read_text())
        irr2key = {"X5": "X5", "X3": "X3", "W4": "W4", "Delta": "D",
                   "Gamma1": "G1", "Gamma3": "G3"}
        for m in v["modes"]:
            for irr, key in irr2key.items():
                if m["irrep"].startswith(irr) and m.get("omega_THz"):
                    freqs[key] = float(m["omega_THz"])

    args.outdir.mkdir(parents=True, exist_ok=True)
    entries = []
    dist = 0.0
    for key, F in fields.items():
        # pick the star arm that carries the pattern best
        best = None
        for arm in STAR_ARMS[key]:
            V, w = bloch_eigenvector(F, pos_cart, masses52, arm)
            if V is not None and (best is None or w > best[2]):
                best = (arm, V, w)
        arm, V, w = best
        entries.append({"q": list(arm), "distance": dist,
                        "label": f"{key} (arm weight {w:.2f})",
                        "frequency": freqs.get(key, 0.0),
                        "eigenvector": V})
        dist += 0.1
        print(f"  {key:>6}: q_conv = {list(arm)}, captured weight "
              f"{w:.2f}, freq {freqs.get(key, 0.0):.3f} THz")

    yaml_path = args.outdir / "modes_irrep.yaml"
    emit_yaml(entries, elem52, aligned52, A_CUB, yaml_path)
    print(f"wrote {yaml_path} (compact: cubic cell, doubling in q)")

    # ---- Γ-folded 1x1x2 supercell representation ------------------------
    # Viewers that draw only the base cell and ignore Bloch phases need the
    # doubling explicit (viewer/'s CrystalViewer does apply the Bloch phase,
    # so for it either representation works):
    # the 104-atom tetragonal cell with every mode at q = 0 and the
    # inter-cell alternation baked into a purely real eigenvector.
    slab_frac = np.array([[*aligned52[i % 52][:2],
                           (aligned52[i % 52][2] + i // 52) / 2.0]
                          for i in range(104)])
    masses104 = np.array([atomic_masses[atomic_numbers[e]]
                          for e in slab_elem])
    entries112 = []
    for e_c, (key, F) in zip(entries, fields.items()):
        V = np.sqrt(masses104)[:, None] * F
        V = (V / np.linalg.norm(V)).astype(complex)
        entries112.append({"q": [0.0, 0.0, 0.0], "distance": e_c["distance"],
                           "label": f"{key} (1x1x2 Gamma-folded; "
                                    f"star arm {e_c['q']})",
                           "frequency": e_c["frequency"],
                           "eigenvector": V})
    yaml112 = args.outdir / "modes_irrep_112.yaml"
    emit_yaml(entries112, slab_elem, slab_frac,
              np.diag([A_CUB, A_CUB, 2 * A_CUB]), yaml112)
    print(f"wrote {yaml112} (explicit 1x1x2 cell, all modes at Gamma)")

    # extended-XYZ animations of the 104-atom cell
    from ase import Atoms
    from ase.io import write as ase_write
    ref = mp.load_reference()
    cell = np.diag([ref["cell"]["a_A"], ref["cell"]["a_A"],
                    ref["cell"]["c_A"]])
    for key, F in fields.items():
        frames = []
        for i in range(args.nframes):
            s = np.sin(2 * np.pi * i / args.nframes)
            frames.append(Atoms(symbols=slab_elem,
                                positions=pos_cart + args.exaggeration * s * F,
                                cell=cell, pbc=True))
        ase_write(str(args.outdir / f"{key}.xyz"), frames, format="extxyz")
    print(f"wrote {len(fields)} animated .xyz files "
          f"(exaggeration ×{args.exaggeration})")


if __name__ == "__main__":
    sys.exit(main())
