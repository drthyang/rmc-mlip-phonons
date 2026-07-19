#!/usr/bin/env python3
"""make_synthetic_ensemble.py — synthetic Cu fcc rmc6f ensemble for tests.

Writes N noisy ``.rmc6f`` configurations of a face-centred-cubic copper
supercell. Each configuration is the ideal fcc lattice plus independent
Gaussian displacements (a crude thermal cloud), so that folding + circular
averaging across the ensemble recovers the ideal fcc unit cell.

Known answer (used by the EMT / MACE integration tests): fcc Cu is a simple
monatomic metal, so its phonon dispersion has exactly **3 acoustic branches**
and all frequencies vanish at Γ (ω → 0). There are no optical branches.

Units: Å for lengths, fractional coordinates for positions. Deterministic for
a given ``seed``.

Standalone use::

    python tests/fixtures/make_synthetic_ensemble.py OUTDIR --n 24 \
        --dims 2 2 2 --a 3.61 --sigma 0.08 --seed 0
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

# Conventional fcc basis (4 atoms), fractional coordinates of the unit cell.
FCC_BASIS = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ]
)
# Experimental fcc Cu conventional lattice parameter (Å).
CU_A = 3.61


def build_supercell(dims, a=CU_A):
    """Ideal fcc supercell.

    Parameters
    ----------
    dims : (3,) int
        Supercell repetitions (Nx, Ny, Nz).
    a : float
        Conventional cubic lattice parameter, Å.

    Returns
    -------
    cell : (3, 3) float
        Supercell lattice vectors (rows), Å.
    frac : (N, 3) float
        Ideal fractional coordinates in the *supercell* frame.
    site_ids : (N,) int
        Unit-cell basis-site id (1..4) each atom folds onto.
    offsets : (N, 3) int
        The (i, j, k) unit-cell offset each atom sits in (for the rmc6f
        cell-index columns; the parser ignores them).
    """
    dims = np.asarray(dims, dtype=int)
    cell = np.diag(dims.astype(float) * a)
    frac, site_ids, offsets = [], [], []
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                offset = np.array([i, j, k])
                for b, basis in enumerate(FCC_BASIS, start=1):
                    # supercell frac = (unit_frac + cell_offset) / dims
                    frac.append((basis + offset) / dims)
                    site_ids.append(b)
                    offsets.append(offset)
    return (
        cell,
        np.array(frac),
        np.array(site_ids, dtype=int),
        np.array(offsets, dtype=int),
    )


def _format_config(cell, frac, site_ids, offsets, dims, config_index):
    """Render one rmc6f configuration as text (phonopy-agnostic, tolerant
    format matching milestone1_bands.parse_rmc6f)."""
    n = len(frac)
    lengths = np.linalg.norm(cell, axis=1)
    lines = [
        "(Version 6 format configuration file)",
        "Metadata generated_by: tests/fixtures/make_synthetic_ensemble.py",
        f"Metadata config_index: {config_index}",
        f"Number of atoms: {n}",
        f"Supercell dimensions: {dims[0]} {dims[1]} {dims[2]}",
        "Cell (Ang/deg): "
        f"{lengths[0]:.6f} {lengths[1]:.6f} {lengths[2]:.6f} 90.0 90.0 90.0",
        "Lattice vectors (Ang):",
    ]
    for row in cell:
        lines.append(f"    {row[0]:12.6f} {row[1]:12.6f} {row[2]:12.6f}")
    lines.append("Atoms:")
    for idx in range(n):
        x, y, z = frac[idx]
        sid = site_ids[idx]
        oi, oj, ok = offsets[idx]
        lines.append(
            f"{idx + 1:6d}   Cu   [{sid}]   "
            f"{x:11.6f} {y:11.6f} {z:11.6f}   "
            f"{sid:3d} {oi:3d} {oj:3d} {ok:3d}"
        )
    return "\n".join(lines) + "\n"


def make_fcc_cu_ensemble(
    outdir, n_configs=24, dims=(2, 2, 2), a=CU_A, sigma=0.08, seed=0
):
    """Write *n_configs* noisy fcc-Cu rmc6f files into *outdir*.

    Parameters
    ----------
    outdir : path-like
        Destination directory (created if absent).
    n_configs : int
        Number of configurations to write.
    dims : (3,) int
        Supercell repetitions.
    a : float
        fcc conventional lattice parameter, Å.
    sigma : float
        Standard deviation of the per-atom, per-axis Cartesian Gaussian
        displacement, Å (a crude thermal cloud; averages out over the
        ensemble).
    seed : int
        Seed for the (single) numpy Generator — output is bit-for-bit
        reproducible for a given (seed, n_configs, dims, a, sigma).

    Returns
    -------
    list[pathlib.Path]
        The written files, sorted, named ``config_000.rmc6f`` ...
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    dims = np.asarray(dims, dtype=int)
    cell, ideal_frac, site_ids, offsets = build_supercell(dims, a)
    lengths = np.linalg.norm(cell, axis=1)  # supercell edge lengths, Å

    rng = np.random.default_rng(seed)
    paths = []
    for c in range(n_configs):
        disp_cart = rng.normal(0.0, sigma, size=ideal_frac.shape)  # Å
        noisy_frac = (ideal_frac + disp_cart / lengths) % 1.0
        text = _format_config(cell, noisy_frac, site_ids, offsets, dims, c)
        path = outdir / f"config_{c:03d}.rmc6f"
        path.write_text(text)
        paths.append(path)
    return sorted(paths)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("outdir", type=Path, help="destination directory")
    ap.add_argument("--n", type=int, default=24, help="number of configs")
    ap.add_argument("--dims", type=int, nargs=3, default=(2, 2, 2),
                    metavar=("NX", "NY", "NZ"))
    ap.add_argument("--a", type=float, default=CU_A,
                    help="fcc lattice parameter, Å")
    ap.add_argument("--sigma", type=float, default=0.08,
                    help="Cartesian displacement std, Å")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    paths = make_fcc_cu_ensemble(
        args.outdir, n_configs=args.n, dims=tuple(args.dims),
        a=args.a, sigma=args.sigma, seed=args.seed,
    )
    print(f"wrote {len(paths)} configs to {args.outdir} "
          f"({4 * int(np.prod(args.dims))} atoms each, fcc Cu, "
          f"sigma={args.sigma} Å, seed={args.seed})")


if __name__ == "__main__":
    main()
