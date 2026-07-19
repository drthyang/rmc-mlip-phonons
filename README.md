# rmc-mlip-phonons

Correct phonon bands and static-vs-dynamic mode classification for
[RMCProfile](https://rmcprofile.ornl.gov) ensembles, using machine-learned
interatomic potentials (MLIPs).

This is the **native compute pipeline** companion to
[`rmc-phonon-dynamics`](https://github.com/drthyang/rmc-phonon-dynamics),
which stays a pure client-side viewer. The contract between the two repos is
the file bundle this pipeline emits:

| File           | Consumed by the viewer as                                  |
| -------------- | ---------------------------------------------------------- |
| `band.yaml`    | a phonopy band structure, overlaid on the covariance bands |
| `relaxed.cif`  | the equilibrium displacement reference                     |
| `summary.json` | provenance and stability report                            |
| `verdicts.json`| per-mode static/dynamic badges *(milestone 3)*             |

## Why

Covariance-derived phonons invert amplitudes into frequencies, and RMC
amplitudes are noise-limited and under-constrained (transverse correlations,
near-Γ acoustics in a fixed box, static-disorder inflation). This pipeline
derives frequencies from **MLIP forces** instead; the RMC ensemble supplies
the experiment-constrained structures, sampling geometry, and amplitudes.

## Milestones

1. **Harmonic bands** *(this repo, now)* — fold + circular-average the
   ensemble to one unit cell, MLIP-relax, phonopy finite displacements,
   export `band.yaml` + `relaxed.cif`. Script: `milestone1_bands.py`.
2. **Finite-temperature bands** — MLIP-MD at the experimental temperature;
   renormalized dispersions and linewidths (TDEP / DynaPhoPy-style), plus
   G(r)/S(Q) closure against the measured data.
3. **Experiment-constrained force constants** — hiPhive fit with RMC
   configurations as displacement snapshots and MLIP forces evaluated on
   them; per-mode amplitude vs (ħ/2ω)coth(ħω/2k_BT), E(Q) mode mapping,
   `verdicts.json`.
4. **Application + methods paper** — real datasets, active-learning
   fine-tuning where the foundation model falls short.

## Quickstart (milestone 1)

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# all *.rmc6f in a run folder, MACE-MP-0 medium on CPU:
python milestone1_bands.py path/to/run_dir/ -o m1_out

# GPU, explicit phonopy supercell, smaller band.yaml:
python milestone1_bands.py run_dir/ --device cuda --dim 3 3 2 --no-eigenvectors

# large ensemble: keep every 2nd config, then a deterministic 100-config sample:
python milestone1_bands.py run_dir/ --stride 2 --max-configs 100 --seed 0
```

Then open `m1_out/band.yaml` in the `rmc-phonon-dynamics` viewer next to the
covariance bands, and use `m1_out/relaxed.cif` as the displacement reference.

`--calc emt` runs a dependency-free smoke test of the plumbing (metals-only
toy potential — not for science).

## Caveats

- Universal MLIPs (MACE-MP-0, CHGNet) inherit PBE-level systematics and a
  few-percent global softening: judge band **shape** first; fine-tune on DFT
  for absolute frequencies (milestone 4).
- Finite displacements default to **0.03 Å**, not phonopy's DFT-tuned 0.01 Å:
  at 0.01 Å a universal MLIP's forces sit near its noise floor and spurious
  imaginary modes appear just off Γ even for a stable crystal. Override with
  `--displacement`.
- Phonons need an **ordered** cell: mixed-occupancy sites are reduced to the
  majority element with a warning.
- A diagonal RMCProfile supercell (Nx Ny Nz) is assumed.
- Imaginary modes on the averaged/symmetrized cell are a *result*, not an
  error: they flag a candidate frozen distortion. Re-run on the distorted
  subgroup cell, or handle at temperature in milestones 2–3.

## License

MIT © 2026 Tsung-Han Yang
