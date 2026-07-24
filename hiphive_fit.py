#!/usr/bin/env python3
"""hiphive_fit.py — milestone 3: experiment-constrained effective force constants.

RMC configurations serve as displacement snapshots; the MLIP evaluates
forces on them; hiPhive fits 2nd-order effective force constants; phonopy
renders the resulting bands. Output is a NEW sidecar (`band_rmc.yaml`) —
`band.yaml` stays frozen per the interchange contract.

    rmc6f configs ──> displacements about the (fixed-cell MLIP-relaxed)
    experimental-average reference ──> MLIP forces per config ──> hiPhive
    2nd-order fit ──> band_rmc.yaml + fit_report.json

Interpretation note (docs/verdicts-schema.md): the RMC displacement
distribution contains genuine quantum+static content AND RMC fitting noise
(~80 % of the variance for the GTS 5 K set), so these effective FCs sample
the MLIP surface over a broader distribution than physical motion — compare
against `band.yaml` (harmonic) and `band_T.yaml` (quantum-sampled) to
bracket the effect. Units: Å, eV, THz; MLIP forces in float64.

Usage:
    python hiphive_fit.py run_dir/ --skip-nonconverged --nconfigs 8 -o m3_out
    python hiphive_fit.py ens/ --calc emt --cutoff 3.0 -o m3_out   # smoke
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

import milestone1_bands as m1
import mode_project as mpj


# ----------------------------------------------------------------------------
# geometry mapping (pure functions)
# ----------------------------------------------------------------------------

def config_site_cell(cfg):
    """Per-atom (unit_frac, sid, ijk) from a parsed rmc6f config.

    Derived from the fractional supercell coordinates and dims — no reliance
    on the file's trailing cell-index columns, so it works for any diagonal
    supercell (synthetic fixtures included). unit_frac is in the unit-cell
    frame [0, 1).
    """
    dims = cfg["dims"]
    x = cfg["frac"] * dims
    ijk = np.floor(x + 1e-9).astype(int) % dims
    unit = x - np.floor(x + 1e-9)
    return unit, cfg["site_ids"].astype(int), ijk


def box_order(sid, ijk, dims, n_sites):
    """Deterministic atom ordering for the supercell: cell-major, then site.

    index = ((i*ny + j)*nz + k) * n_sites + (sid − 1). Guarantees config and
    ideal supercells share atom order (required by hiPhive).
    """
    nx, ny, nz = (int(d) for d in dims)
    cell = (ijk[:, 0] * ny + ijk[:, 1]) * nz + ijk[:, 2]
    return cell * n_sites + (sid - 1)


def build_box(unit_frac, sid, ijk, dims, elems_by_sid, unit_cell):
    """ASE Atoms of the full box in `box_order` ordering (positions Å)."""
    from ase import Atoms

    n_sites = len(elems_by_sid)
    order = box_order(sid, ijk, dims, n_sites)
    inv = np.argsort(order)
    frac_box = (unit_frac[inv] + ijk[inv]) / np.asarray(dims)
    symbols = [elems_by_sid[s - 1] for s in sid[inv]]
    cell = np.asarray(dims)[:, None] * np.asarray(unit_cell)
    return Atoms(symbols=symbols, scaled_positions=frac_box % 1.0,
                 cell=cell, pbc=True)


def ideal_box(ref_frac, elems, dims, unit_cell):
    """The reference structure tiled over the box, in `box_order` order."""
    from ase import Atoms

    n_sites = len(elems)
    dims = np.asarray(dims, dtype=int)
    sid, ijk, unit = [], [], []
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                for s in range(n_sites):
                    sid.append(s + 1)
                    ijk.append((i, j, k))
                    unit.append(ref_frac[s])
    return build_box(np.array(unit), np.array(sid),
                     np.array(ijk, dtype=int), dims, elems, unit_cell)


# ----------------------------------------------------------------------------
# reference structure
# ----------------------------------------------------------------------------

def build_reference(configs, calc, fmax):
    """Experimental-average cell, idealized, fixed-cell MLIP-relaxed.

    Site order = rmc6f site-id order throughout (fold_and_average returns
    sites sorted by id). Returns (atoms, elems_by_sid).
    """
    from ase import Atoms
    from ase.optimize import FIRE

    unit_cell, elems, frac, _ = m1.fold_and_average(configs, None)
    ideal = mpj.idealize_parent(frac, elems)
    atoms = Atoms(symbols=elems, scaled_positions=ideal,
                  cell=unit_cell, pbc=True)
    atoms.calc = calc
    FIRE(atoms, logfile=None).run(fmax=fmax, steps=1000)
    shift = np.abs(((atoms.get_scaled_positions() - ideal + 0.5) % 1.0)
                   - 0.5).max()
    print(f"  reference: {len(atoms)} sites, fixed-cell relax shift "
          f"{shift:.4f} frac")
    return atoms, list(elems)


# ----------------------------------------------------------------------------
# fit and bands
# ----------------------------------------------------------------------------

def fit_effective_fcs(structures, ideal, unit_atoms, cutoff):
    """hiPhive 2nd-order fit. Returns (fcp, report dict)."""
    from hiphive import ClusterSpace, StructureContainer
    from hiphive.utilities import prepare_structures
    from trainstation import Optimizer

    cs = ClusterSpace(unit_atoms, [cutoff])
    sc = StructureContainer(cs)
    for s in prepare_structures(structures, ideal):
        sc.add_structure(s)
    opt = Optimizer(sc.get_fit_data())
    opt.train()
    from hiphive import ForceConstantPotential
    fcp = ForceConstantPotential(cs, opt.parameters)
    report = {"n_parameters": int(cs.n_dofs),
              "rmse_train_eV_A": float(opt.rmse_train),
              "rmse_test_eV_A": float(opt.rmse_test)}
    print(f"  fit: {report['n_parameters']} params, rmse train/test = "
          f"{report['rmse_train_eV_A']:.4f}/{report['rmse_test_eV_A']:.4f}"
          " eV/Å")
    return fcp, report


def bands_from_fcp(fcp, unit_atoms, band_dim, npoints, symprec, outpath,
                   eigenvectors=True):
    """Render band_rmc.yaml from the fitted FCP on a small supercell.

    eigenvectors=True (default) writes the complex eigenvectors along the
    standard path — the same convention as milestone-1 band.yaml."""
    from ase import Atoms
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    # See m1.symmetrize_lattice: the band path comes from seekpath at its
    # hard-coded symprec=1e-5, so the cell metric is idealized first. Lattice
    # only — fractional coordinates and atom order are untouched, which the
    # FCP evaluation on ph.supercell below depends on.
    unit_atoms, _ = m1.symmetrize_lattice(unit_atoms, symprec)
    unit = PhonopyAtoms(symbols=unit_atoms.get_chemical_symbols(),
                        cell=unit_atoms.cell.array,
                        scaled_positions=unit_atoms.get_scaled_positions())
    ph = Phonopy(unit, supercell_matrix=np.diag(band_dim),
                 primitive_matrix="auto", symprec=symprec)
    sc = ph.supercell
    sc_ase = Atoms(symbols=sc.symbols, scaled_positions=sc.scaled_positions,
                   cell=sc.cell, pbc=True)
    ph.force_constants = fcp.get_force_constants(sc_ase).get_fc_array(order=2)
    ph.auto_band_structure(npoints=npoints, with_eigenvectors=eigenvectors,
                           write_yaml=True, filename=str(outpath))
    freqs = np.concatenate([f.ravel() for f in ph.band_structure.frequencies])
    return ph, freqs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--exclude", default=None,
                    help="drop input files whose NAME matches this regex "
                    "(e.g. AVERAGE to keep only instantaneous configs)")
    ap.add_argument("--skip-nonconverged", action="store_true")
    ap.add_argument("--nconfigs", type=int, default=8,
                    help="snapshots used for the fit (seeded random subset)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--calc", default="mace", choices=["mace", "chgnet", "emt"])
    ap.add_argument("--model", default="small",
                    choices=["small", "medium", "large"])
    ap.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    ap.add_argument("--cutoff", type=float, default=6.0,
                    help="hiPhive 2nd-order cutoff, Å (< half the band-dim "
                    "supercell edge)")
    ap.add_argument("--fmax", type=float, default=1e-3)
    ap.add_argument("--band-dim", type=int, nargs=3, default=(2, 2, 2))
    ap.add_argument("--npoints", type=int, default=51)
    ap.add_argument("--no-eigenvectors", action="store_true",
                    help="omit eigenvectors (much smaller band_rmc.yaml)")
    ap.add_argument("--symprec", type=float, default=0.1)
    ap.add_argument("-o", "--outdir", type=Path, default=Path("m3_out"))
    args = ap.parse_args(argv)

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    import re as _re
    files = m1.collect_inputs(args.inputs)
    if args.exclude:
        n0 = len(files)
        files = [f for f in files if not _re.search(args.exclude, f.name)]
        print(f"  --exclude {args.exclude!r}: {n0} -> {len(files)} files")
    if args.skip_nonconverged:
        files, dropped = m1.drop_nonconverged(files)
        if dropped:
            print(f"  dropped {len(dropped)} non-converged config(s)")
    files = m1.select_configs(files, 1, args.nconfigs, args.seed)
    print(f"[1/4] {len(files)} snapshot configs; parsing")
    configs = [m1.parse_rmc6f(f) for f in files]

    print(f"[2/4] reference structure ({args.calc})")
    calc = m1.get_calculator(args.calc, args.device, args.model)
    ref_atoms, elems = build_reference(configs, calc, args.fmax)
    ref_frac = ref_atoms.get_scaled_positions() % 1.0
    dims = configs[0]["dims"]
    ideal = ideal_box(ref_frac, elems, dims, ref_atoms.cell.array)

    print(f"[3/4] MLIP forces on {len(configs)} boxes "
          f"({len(ideal)} atoms each)")
    from ase.calculators.singlepoint import SinglePointCalculator
    structures = []
    for f, cfg in zip(files, configs):
        cache = outdir / f"forces_{f.stem}.npy"
        unit, sid, ijk = config_site_cell(cfg)
        atoms = build_box(unit, sid, ijk, dims, elems, ref_atoms.cell.array)
        if cache.is_file():
            forces = np.load(cache)
            print(f"    {f.name}: cached")
        else:
            atoms.calc = calc
            t0 = time.time()
            forces = atoms.get_forces()
            np.save(cache, forces)
            print(f"    {f.name}: forces in {time.time()-t0:.0f}s")
        atoms.calc = SinglePointCalculator(atoms, energy=0.0, forces=forces)
        structures.append(atoms)

    print("[4/4] hiPhive fit -> band_rmc.yaml")
    fcp, report = fit_effective_fcs(structures, ideal, ref_atoms, args.cutoff)
    band_path = outdir / "band_rmc.yaml"
    ph_eff, f_eff = bands_from_fcp(fcp, ref_atoms, args.band_dim,
                                   args.npoints, args.symprec, band_path,
                                   eigenvectors=not args.no_eigenvectors)

    # harmonic comparison with the same MLIP/reference
    import md_run as m2
    ph_h = m2.harmonic_model(ref_atoms, calc, np.array(args.band_dim), 0.03,
                             args.symprec)
    ph_h.auto_band_structure(npoints=args.npoints)
    f_h = np.concatenate([x.ravel() for x in ph_h.band_structure.frequencies])
    report.update({
        "n_snapshots": len(structures),
        "cutoff_A": args.cutoff,
        "reference": "experimental-average, idealized, fixed-cell relaxed",
        "band_rmc_freq_range_THz": [float(f_eff.min()), float(f_eff.max())],
        "harmonic_freq_range_THz": [float(f_h.min()), float(f_h.max())],
        "max_abs_dw_vs_harmonic_THz": float(np.abs(f_eff - f_h).max()),
        "mean_dw_vs_harmonic_THz": float((f_eff - f_h).mean()),
        "inputs": [str(f) for f in files],
    })
    (outdir / "fit_report.json").write_text(json.dumps(report, indent=2))
    print(f"  band_rmc vs harmonic: max |Δω| = "
          f"{report['max_abs_dw_vs_harmonic_THz']:.3f} THz, mean Δω = "
          f"{report['mean_dw_vs_harmonic_THz']:+.3f} THz")
    print(f"  wrote {band_path} and fit_report.json")


if __name__ == "__main__":
    sys.exit(main())
