#!/usr/bin/env python3
"""md_run.py — milestone 2: finite-temperature closure and effective force constants.

Two sampling modes produce configuration ensembles from the MLIP model of the
averaged structure; both feed the same G(r)/F(Q) closure against measured
total-scattering data and (optionally) a hiPhive effective-force-constant fit:

    sample : harmonic QUANTUM sampling — random supercell displacements drawn
             from the phonopy model with quantum statistics (zero-point
             included) at temperature T. The physically correct sampler at
             low T (design decision D1, docs/milestone2-plan.md).
    md     : classical ASE Langevin NVT — kept as an explicitly classical
             cross-check. At low T classical equipartition amplitudes are
             badly wrong (no zero-point); a loud warning is printed.

Pipeline:

    rmc6f ensemble (or --cif)  ->  averaged/symmetrized cell at the
    EXPERIMENTAL lattice (fixed-cell MLIP relax; design D2)
      ->  phonopy finite-displacement force constants on the sampling
          supercell  ->  snapshots  ->  partial g_ij(r)  ->  neutron-weighted
          G(r)  ->  F(Q) = S(Q) - 1 on the measured Q grid
      ->  scale+offset fit vs measured data (mirrors the RMC treatment)
      ->  closure.json + gr_sim.dat + sq_sim.dat (+ band_T.yaml via hiPhive)

Units: Å, eV, THz throughout (phonopy defaults); neutron b_coh in fm
(Sears 1992). All MLIP evaluations use default_dtype="float64".

Usage examples:
    python md_run.py sample run_dir/ --skip-nonconverged \\
        --data data/5K_ini/scale_ft_rmc.fq -T 5 -o m2_out
    python md_run.py sample --cif m2_out/relaxed_expt.cif --data F.fq -T 5
    python md_run.py md run_dir/ --data F.fq -T 300 --md-steps 20000
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

import milestone1_bands as m1

# Coherent neutron scattering lengths, fm (V. F. Sears, Neutron News 3, 26
# (1992)). Extend as needed for new materials.
BCOH_FM = {
    "H": -3.739, "C": 6.646, "N": 9.36, "O": 5.803, "Al": 3.449,
    "Si": 4.1491, "S": 2.847, "Cu": 7.718, "Ga": 7.288, "Ge": 8.185,
    "Se": 7.970, "Nb": 7.054, "Mo": 6.715, "Ta": 6.91, "W": 4.86,
    "Pt": 9.60, "Au": 7.63, "Pb": 9.405, "V": -0.3824, "Ti": -3.438,
    "Ni": 10.3, "Fe": 9.45, "Zn": 5.680, "Zr": 7.16, "Ag": 5.922,
}


# ----------------------------------------------------------------------------
# measured-data I/O and closure metrics
# ----------------------------------------------------------------------------

def parse_fq(path: Path):
    """Read an RMCProfile/STOG-style .fq file.

    Format: first line = number of points, second line = title, then two
    columns Q (Å⁻¹) and F(Q) = S(Q) − 1 (dimensionless).

    Returns
    -------
    (Q, F) : two (n,) float arrays.
    """
    lines = Path(path).read_text().strip().splitlines()
    n = int(lines[0].split()[0])
    rows = [l.split() for l in lines[2:2 + n]]
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1]


def fit_scale_offset(f_data: np.ndarray, f_sim: np.ndarray):
    """Least-squares scale s and offset o minimising ‖F_data − s·F_sim − o‖.

    Mirrors the RMC fit convention (scale and offset free). Returns
    (s, o, Rw) with Rw = sqrt(Σ residual² / Σ F_data²), dimensionless.
    """
    A = np.column_stack([f_sim, np.ones_like(f_sim)])
    (s, o), *_ = np.linalg.lstsq(A, f_data, rcond=None)
    resid = f_data - (s * f_sim + o)
    rw = float(np.sqrt((resid**2).sum() / (f_data**2).sum()))
    return float(s), float(o), rw


# ----------------------------------------------------------------------------
# pair distribution machinery (pure numpy, unit-tested)
# ----------------------------------------------------------------------------

def neutron_weights(symbols):
    """Faber–Ziman pair weights for a composition.

    Parameters
    ----------
    symbols : list[str]
        Per-atom element symbols of one configuration.

    Returns
    -------
    weights : dict[(el_a, el_b)] -> float, for a <= b (cross terms carry the
        factor 2); Σ weights = 1 by construction.
    species : sorted list of the distinct elements.
    """
    symbols = list(symbols)
    species = sorted(set(symbols))
    missing = [s for s in species if s not in BCOH_FM]
    if missing:
        raise SystemExit(f"no b_coh tabulated for {missing}; extend BCOH_FM")
    n = len(symbols)
    c = {s: symbols.count(s) / n for s in species}
    b = {s: BCOH_FM[s] for s in species}
    bbar = sum(c[s] * b[s] for s in species)
    weights = {}
    for i, a in enumerate(species):
        for bb in species[i:]:
            mult = 1.0 if a == bb else 2.0
            weights[(a, bb)] = mult * c[a] * c[bb] * b[a] * b[bb] / bbar**2
    return weights, species


def pair_histograms(frames, symbols, cell, r_max, dr):
    """Partial pair-distribution functions g_ab(r) from snapshots.

    Parameters
    ----------
    frames : list of (N, 3) arrays
        Cartesian positions per snapshot, Å. Orthorhombic cell assumed.
    symbols : list[str]
        Per-atom elements (identical ordering across frames).
    cell : (3, 3) array
        Supercell lattice vectors (rows), Å; must be diagonal.
    r_max, dr : float
        Histogram range and bin width, Å. r_max must be ≤ half the shortest
        box edge (minimum-image validity).

    Returns
    -------
    r : (nbins,) bin centres, Å.
    g : dict[(a, b)] -> (nbins,) partial g_ab(r), a <= b, normalised to 1 at
        large r.
    """
    cell = np.asarray(cell, dtype=float)
    if not np.allclose(cell, np.diag(np.diag(cell)), atol=1e-8):
        raise SystemExit("pair_histograms: only diagonal (orthorhombic) cells")
    L = np.diag(cell)
    if r_max > L.min() / 2 + 1e-9:
        raise SystemExit(f"r_max={r_max} exceeds half the box ({L.min()/2:.2f} Å)")
    vol = float(np.prod(L))
    symbols = np.asarray(symbols)
    species = sorted(set(symbols.tolist()))
    idx = {s: np.where(symbols == s)[0] for s in species}
    nbins = int(round(r_max / dr))
    edges = np.linspace(0.0, r_max, nbins + 1)
    r = 0.5 * (edges[1:] + edges[:-1])
    counts = {}
    for i, a in enumerate(species):
        for b in species[i:]:
            counts[(a, b)] = np.zeros(nbins)

    chunk = max(1, int(2**22 // max(len(symbols), 1)))
    for pos in frames:
        pos = np.asarray(pos, dtype=float)
        for i, a in enumerate(species):
            pa = pos[idx[a]]
            for b in species[i:]:
                pb = pos[idx[b]]
                h = np.zeros(nbins)
                for i0 in range(0, len(pa), chunk):
                    d = pa[i0:i0 + chunk, None, :] - pb[None, :, :]
                    d -= np.round(d / L) * L
                    rr = np.sqrt((d * d).sum(-1)).ravel()
                    if a == b:
                        rr = rr[rr > 1e-8]  # drop self-pairs
                    h += np.histogram(rr, bins=edges)[0]
                counts[(a, b)] += h

    nf = len(frames)
    g = {}
    for (a, b), h in counts.items():
        if a == b:
            n_ord = len(idx[a]) * (len(idx[a]) - 1)
        else:
            n_ord = len(idx[a]) * len(idx[b])  # one direction only
        shell = 4.0 * np.pi * r**2 * dr
        g[(a, b)] = vol * h / (nf * n_ord * shell)
    return r, g


def total_G(r, g_partials, symbols):
    """Neutron-weighted total G(r) = Σ w_ab (g_ab(r) − 1)  (Faber–Ziman)."""
    weights, _ = neutron_weights(symbols)
    G = np.zeros_like(r)
    for pair, w in weights.items():
        G += w * (g_partials[pair] - 1.0)
    return G


def gr_to_fq(r, G, Q, rho0):
    """F(Q) = S(Q) − 1 = 4πρ₀ ∫ r² G(r) sinc(Qr) dr  (trapezoid on r grid).

    rho0 is the atomic number density, Å⁻³.
    """
    r = np.asarray(r)
    dr = r[1] - r[0]
    Qr = np.outer(Q, r)
    with np.errstate(invalid="ignore", divide="ignore"):
        sinc = np.where(Qr > 1e-12, np.sin(Qr) / Qr, 1.0)
    return 4.0 * np.pi * rho0 * (sinc * (r**2 * G)[None, :]).sum(axis=1) * dr


def fq_to_gr(Q, F, r, rho0):
    """Inverse transform: G(r) = 1/(2π²ρ₀) ∫ Q² F(Q) sinc(Qr) dQ."""
    Q = np.asarray(Q)
    dQ = Q[1] - Q[0]
    Qr = np.outer(r, Q)
    with np.errstate(invalid="ignore", divide="ignore"):
        sinc = np.where(Qr > 1e-12, np.sin(Qr) / Qr, 1.0)
    return (sinc * (Q**2 * F)[None, :]).sum(axis=1) * dQ / (2.0 * np.pi**2 * rho0)


# ----------------------------------------------------------------------------
# structure preparation and harmonic model
# ----------------------------------------------------------------------------

def build_structure(args, calc):
    """Averaged, symmetrized cell relaxed at the closure lattice.

    Default (design D2): the folded cell keeps the *experimental* lattice
    (the RMC box is the experimental lattice); only internal coordinates are
    relaxed. --free-lattice relaxes the cell too (MLIP lattice).
    """
    from ase import Atoms
    from ase.io import read as ase_read
    from ase.optimize import FIRE

    if args.cif is not None:
        atoms = ase_read(str(args.cif))
        print(f"  structure from {args.cif}: {len(atoms)} atoms")
    else:
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
        files = m1.select_configs(files, args.stride, args.max_configs,
                                  args.seed)
        print(f"  folding {len(files)} configuration(s)")
        configs = [m1.parse_rmc6f(f) for f in files]
        unit_cell, elems, frac, _ = m1.fold_and_average(configs, args.ref)
        atoms = Atoms(symbols=elems, scaled_positions=frac,
                      cell=unit_cell, pbc=True)
        atoms = m1.symmetrize(atoms, args.symprec)

    atoms.calc = calc
    if args.free_lattice:
        try:
            from ase.filters import FrechetCellFilter as CellFilter
        except ImportError:
            from ase.constraints import ExpCellFilter as CellFilter
        FIRE(CellFilter(atoms), logfile=None).run(fmax=args.fmax, steps=1500)
    FIRE(atoms, logfile=None).run(fmax=args.fmax, steps=1000)
    fres = float(np.abs(atoms.get_forces()).max())
    print(f"  relaxed ({'free' if args.free_lattice else 'fixed'} lattice) "
          f"a = {atoms.cell.lengths().round(4)}  max|F| = {fres:.1e} eV/Å")
    return atoms


def sampling_dim(atoms, r_max, fc_target=12.0):
    """Supercell that supports both the FC range (≥ fc_target Å) and the
    pair histogram (box ≥ 2·r_max)."""
    lengths = atoms.cell.lengths()
    dim = np.maximum(np.ceil(fc_target / lengths),
                     np.ceil(2.0 * r_max / lengths)).astype(int)
    return np.maximum(dim, 1)


def harmonic_model(atoms, calc, dim, displacement, symprec):
    """Phonopy object with finite-displacement force constants on `dim`."""
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from ase import Atoms

    unit = PhonopyAtoms(symbols=atoms.get_chemical_symbols(),
                        cell=atoms.cell.array,
                        scaled_positions=atoms.get_scaled_positions())
    phonon = Phonopy(unit, supercell_matrix=np.diag(dim),
                     primitive_matrix="auto", symprec=symprec)
    phonon.generate_displacements(distance=displacement)
    scells = phonon.supercells_with_displacements
    print(f"  {len(scells)} FC displacements on supercell {dim.tolist()} "
          f"({len(scells[0])} atoms)")
    forces = []
    for n, sc in enumerate(scells, 1):
        a = Atoms(symbols=sc.symbols, scaled_positions=sc.scaled_positions,
                  cell=sc.cell, pbc=True)
        a.calc = calc
        t0 = time.time()
        forces.append(a.get_forces())
        print(f"    FC forces {n}/{len(scells)}  ({time.time()-t0:.1f}s)")
    phonon.forces = np.array(forces)
    phonon.produce_force_constants()
    return phonon


def quantum_snapshots(phonon, n_snapshots, temperature, seed):
    """Random supercell displacements with quantum statistics at T.

    Uses phonopy's random-displacement sampler (harmonic oscillator
    distribution including zero-point occupation). Returns (snapshots,
    ideal) as ASE Atoms; displacement statistics follow
    <u²> ∝ (ħ/2ω) coth(ħω/2k_BT) per mode.
    """
    from ase import Atoms

    # phonopy >=4: assigning the random-displacement dataset invalidates the
    # derived force constants / dynamical matrix — save and restore them so
    # the phonon object stays usable (band structure, thermal displacements).
    fc = phonon.force_constants
    phonon.generate_displacements(number_of_snapshots=n_snapshots,
                                  temperature=temperature,
                                  random_seed=seed)
    scells = phonon.supercells_with_displacements
    phonon.force_constants = fc
    snaps = [Atoms(symbols=sc.symbols, scaled_positions=sc.scaled_positions,
                   cell=sc.cell, pbc=True) for sc in scells]
    sc0 = phonon.supercell
    ideal = Atoms(symbols=sc0.symbols, scaled_positions=sc0.scaled_positions,
                  cell=sc0.cell, pbc=True)
    return snaps, ideal


def md_snapshots(atoms, calc, dim, temperature, timestep_fs, n_steps,
                 n_equil, interval, seed):
    """Classical Langevin NVT snapshots (cross-check sampler).

    WARNING printed at low T: classical equipartition misses zero-point
    motion entirely — amplitudes of modes with ħω >> k_BT are badly
    underestimated (docs/milestone2-plan.md, D1).
    """
    from ase import units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

    if temperature < 100.0:
        print(f"  ! classical MD at T = {temperature} K: amplitudes are "
              "equipartition (no zero-point) — use `sample` for physics; "
              "this mode is a classical cross-check only.")
    sup = atoms.repeat(tuple(int(x) for x in dim))
    sup.calc = calc
    rng = np.random.RandomState(seed)
    MaxwellBoltzmannDistribution(sup, temperature_K=temperature, rng=rng)
    dyn = Langevin(sup, timestep=timestep_fs * units.fs,
                   temperature_K=temperature, friction=0.02, rng=rng)
    frames = []
    dyn.run(n_equil)
    for _ in range(n_steps // interval):
        dyn.run(interval)
        frames.append(sup.get_positions().copy())
    ideal = atoms.repeat(tuple(int(x) for x in dim))
    print(f"  {len(frames)} MD snapshots ({len(sup)} atoms, "
          f"{timestep_fs} fs, friction 0.02)")
    return frames, ideal, sup.get_chemical_symbols()


# ----------------------------------------------------------------------------
# effective force constants (band_T.yaml) via hiPhive
# ----------------------------------------------------------------------------

def effective_bandT(phonon, atoms_unit, snapshots, forces, ideal, cutoff,
                    npoints, outpath, symprec, eigenvectors=True):
    """Fit 2nd-order effective FCs to (displacement, force) pairs; write
    band_T.yaml (with eigenvectors along the standard path by default).
    Returns max |ω_T − ω_harmonic| over the band path, THz."""
    from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential
    from hiphive.utilities import prepare_structures
    from trainstation import Optimizer
    from ase.calculators.singlepoint import SinglePointCalculator
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    for s, f in zip(snapshots, forces):
        s.calc = SinglePointCalculator(s, energy=0.0, forces=f)
    cs = ClusterSpace(atoms_unit, [cutoff])
    structures = prepare_structures(snapshots, ideal)
    sc = StructureContainer(cs)
    for s in structures:
        sc.add_structure(s)
    opt = Optimizer(sc.get_fit_data(), train_size=1.0)
    opt.train()
    print(f"  hiPhive fit: rmse train = {opt.rmse_train:.2e} eV/Å, "
          f"{cs.n_dofs} parameters, cutoff {cutoff} Å")
    fcp = ForceConstantPotential(cs, opt.parameters)
    fc2 = fcp.get_force_constants(ideal).get_fc_array(order=2)

    unit = PhonopyAtoms(symbols=atoms_unit.get_chemical_symbols(),
                        cell=atoms_unit.cell.array,
                        scaled_positions=atoms_unit.get_scaled_positions())
    dim = np.rint(np.diag(ideal.cell.array) /
                  np.diag(atoms_unit.cell.array)).astype(int)
    ph_T = Phonopy(unit, supercell_matrix=np.diag(dim),
                   primitive_matrix="auto", symprec=symprec)
    ph_T.force_constants = fc2
    ph_T.auto_band_structure(npoints=npoints, with_eigenvectors=eigenvectors,
                             write_yaml=True, filename=str(outpath))
    phonon.auto_band_structure(npoints=npoints)
    f_h = np.concatenate([f.ravel() for f in phonon.band_structure.frequencies])
    f_t = np.concatenate([f.ravel() for f in ph_T.band_structure.frequencies])
    dmax = float(np.abs(f_t - f_h).max())
    print(f"  band_T vs harmonic: max |Δω| = {dmax:.4f} THz")
    return dmax


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=["sample", "md"],
                    help="sample = harmonic quantum sampling (physical at low "
                    "T); md = classical Langevin cross-check")
    ap.add_argument("inputs", nargs="*",
                    help=".rmc6f files/dirs (or use --cif)")
    ap.add_argument("--cif", type=Path, default=None,
                    help="skip folding; use this structure as the unit cell")
    ap.add_argument("--ref", type=Path, default=None)
    ap.add_argument("--exclude", default=None,
                    help="drop input files whose NAME matches this regex "
                    "(e.g. AVERAGE to keep only instantaneous configs)")
    ap.add_argument("--skip-nonconverged", action="store_true")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-configs", type=int, default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--calc", default="mace", choices=["mace", "chgnet", "emt"])
    ap.add_argument("--model", default="small",
                    choices=["small", "medium", "large"])
    ap.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    ap.add_argument("--symprec", type=float, default=1e-3)
    ap.add_argument("--fmax", type=float, default=1e-3)
    ap.add_argument("--free-lattice", action="store_true",
                    help="relax the cell too (default: keep the experimental "
                    "lattice of the folded RMC box; design D2)")
    ap.add_argument("--displacement", type=float, default=0.03)
    ap.add_argument("-T", "--temperature", type=float, default=5.0,
                    help="sampling temperature, K")
    ap.add_argument("--nsnapshots", type=int, default=32)
    ap.add_argument("--rmax", type=float, default=20.0,
                    help="G(r) range, Å (sets the sampling supercell)")
    ap.add_argument("--dr", type=float, default=0.02, help="G(r) bin, Å")
    ap.add_argument("--data", type=Path, default=None,
                    help="measured .fq file (F(Q)=S(Q)−1) for closure")
    ap.add_argument("--no-band-t", action="store_true",
                    help="skip the hiPhive effective-FC fit / band_T.yaml")
    ap.add_argument("--cutoff", type=float, default=6.0,
                    help="hiPhive 2nd-order cutoff, Å")
    ap.add_argument("--npoints", type=int, default=51)
    ap.add_argument("--md-steps", type=int, default=4000)
    ap.add_argument("--md-equil", type=int, default=1000)
    ap.add_argument("--md-interval", type=int, default=100)
    ap.add_argument("--md-timestep", type=float, default=2.0, help="fs")
    ap.add_argument("-o", "--outdir", type=Path, default=Path("m2_out"))
    args = ap.parse_args(argv)

    if args.cif is None and not args.inputs:
        ap.error("give rmc6f inputs or --cif")
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[1/5] structure ({args.calc})")
    calc = m1.get_calculator(args.calc, args.device, args.model)
    atoms = build_structure(args, calc)
    from ase.io import write as ase_write
    ase_write(str(outdir / "relaxed_expt.cif"), atoms)
    ladder = m1.spacegroup_ladder(atoms)
    print("  spacegroup ladder:", ladder)

    dim = sampling_dim(atoms, args.rmax)
    print(f"[2/5] sampling supercell {dim.tolist()} "
          f"({int(np.prod(dim))*len(atoms)} atoms), mode = {args.mode}")

    msd = None
    if args.mode == "sample":
        phonon = harmonic_model(atoms, calc, dim, args.displacement,
                                args.symprec)
        snaps, ideal = quantum_snapshots(phonon, args.nsnapshots,
                                         args.temperature, args.seed)
        frames = [s.get_positions() for s in snaps]
        symbols = ideal.get_chemical_symbols()
        disp = np.array([s.get_positions() - ideal.get_positions()
                         for s in snaps])
        msd = float((disp**2).sum(axis=2).mean())
        print(f"  {len(frames)} quantum snapshots at T = {args.temperature} K"
              f";  MSD = {msd:.5f} Å²  (u_rms = {np.sqrt(msd):.3f} Å)")
    else:
        phonon = None
        frames, ideal, symbols = md_snapshots(
            atoms, calc, dim, args.temperature, args.md_timestep,
            args.md_steps, args.md_equil, args.md_interval, args.seed)

    print(f"[3/5] G(r) / F(Q)  (rmax = {args.rmax} Å, dr = {args.dr} Å)")
    r, gpart = pair_histograms(frames, symbols, ideal.cell.array,
                               args.rmax, args.dr)
    G = total_G(r, gpart, symbols)
    rho0 = len(symbols) / ideal.get_volume()
    print(f"  rho0 = {rho0:.6f} Å⁻³")

    fit = None
    if args.data is not None:
        Qd, Fd = parse_fq(args.data)
        Fs = gr_to_fq(r, G, Qd, rho0)
        s, o, rw_q = fit_scale_offset(Fd, Fs)
        Gd = fq_to_gr(Qd, Fd, r, rho0)
        mask = r > 1.5
        s_r, o_r, rw_r = fit_scale_offset(Gd[mask], G[mask])
        fit = {"scale": s, "offset": o, "Rw_Q": rw_q, "Rw_r": rw_r,
               "n_data": len(Qd), "Q_range": [float(Qd[0]), float(Qd[-1])]}
        print(f"  closure vs {args.data.name}: scale = {s:.3f}, "
              f"offset = {o:+.3f}, Rw(Q) = {rw_q:.3f}, Rw(r) = {rw_r:.3f}")
        np.savetxt(outdir / "sq_sim.dat",
                   np.column_stack([Qd, s * Fs + o, Fd]),
                   header="Q(1/A)  F_sim_scaled  F_data   [F=S(Q)-1]")
        np.savetxt(outdir / "gr_sim.dat", np.column_stack([r, G, Gd]),
                   header="r(A)  G_sim  G_data_FT   [Faber-Ziman G(r)]")
    else:
        Qs = np.arange(0.8, 27.0, 0.01)
        Fs = gr_to_fq(r, G, Qs, rho0)
        np.savetxt(outdir / "sq_sim.dat", np.column_stack([Qs, Fs]),
                   header="Q(1/A)  F_sim   [F=S(Q)-1]")
        np.savetxt(outdir / "gr_sim.dat", np.column_stack([r, G]),
                   header="r(A)  G_sim   [Faber-Ziman G(r)]")

    band_t_delta = None
    if args.mode == "sample" and not args.no_band_t:
        print("[4/5] effective FCs -> band_T.yaml (hiPhive)")
        forces = []
        for n, s in enumerate(snaps, 1):
            s.calc = calc
            forces.append(s.get_forces())
            if n % 8 == 0 or n == len(snaps):
                print(f"    snapshot forces {n}/{len(snaps)}")
        band_t_delta = effective_bandT(
            phonon, atoms, snaps, forces, ideal, args.cutoff,
            args.npoints, outdir / "band_T.yaml", args.symprec)
    else:
        print("[4/5] band_T skipped")

    print("[5/5] closure.json")
    closure = {
        "mode": args.mode,
        "structure": {"natom": len(atoms),
                      "a_lattice": [float(x) for x in atoms.cell.lengths()],
                      "lattice": "free" if args.free_lattice else
                                 "experimental (fixed)",
                      "spacegroup_ladder": ladder},
        "sampling": {"supercell": dim.tolist(),
                     "n_snapshots": len(frames),
                     "temperature_K": args.temperature,
                     "seed": args.seed,
                     "msd_A2": msd},
        "gr": {"rmax_A": args.rmax, "dr_A": args.dr,
               "rho0_A-3": float(rho0)},
        "closure_fit": fit,
        "band_T": {"written": band_t_delta is not None,
                   "max_abs_dw_THz_vs_harmonic": band_t_delta,
                   "cutoff_A": args.cutoff},
        "calculator": {"name": args.calc, "model": args.model,
                       "device": args.device, "dtype": "float64"},
    }
    (outdir / "closure.json").write_text(json.dumps(closure, indent=2))
    print(f"  wrote {outdir}/closure.json, gr_sim.dat, sq_sim.dat"
          + (", band_T.yaml" if band_t_delta is not None else ""))


if __name__ == "__main__":
    sys.exit(main())
