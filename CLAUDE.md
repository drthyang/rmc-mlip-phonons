# CLAUDE.md — rmc-mlip-phonons

## What this project is

Native Python pipeline computing **correct phonon bands** and
**static-vs-dynamic mode classification** for RMCProfile ensembles using
machine-learned interatomic potentials. It is the compute companion to
`drthyang/rmc-phonon-dynamics` (browser viewer, 100% client-side — that repo
must stay that way). Rationale in `README.md`, plan in `ROADMAP.md`.

Core idea: the viewer's covariance route inverts noisy **amplitudes** into
frequencies; this pipeline derives frequencies from MLIP **forces**, using the
RMC ensemble for experiment-constrained structures, sampling geometry, and
amplitudes.

## Hard boundaries — do not violate

- **No web UI in this repo.** Visualization lives in `rmc-phonon-dynamics`;
  this repo only emits files the viewer loads.
- **The interchange contract is frozen:** `band.yaml` stays phonopy-standard
  (`auto_band_structure` output). Extensions go in **new sidecar files**
  (`verdicts.json`), never by mutating `band.yaml`.
- Always emit `relaxed.cif` — it is the viewer's displacement reference.
- `default_dtype="float64"` for every MLIP force/phonon evaluation.
- `data/`, `results/`, `m1_out*/` are git-ignored; never commit ensembles,
  trajectories, or model weights.
- Units are Å, eV, THz (phonopy defaults); state units in every docstring.

## Current state

`milestone1_bands.py` is a complete end-to-end script
(rmc6f → fold/circular-average → spglib symmetrize → MLIP relax → phonopy →
`band.yaml` + `relaxed.cif` + `summary.json`). The parser and circular
averaging are unit-verified; **the ASE/phonopy/MACE leg has never run against
installed dependencies** — expect minor API friction and fix forward.

## Environment / commands

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python milestone1_bands.py <run_dir>/ -o m1_out     # real run (MACE-MP-0)
python milestone1_bands.py <run_dir>/ --calc emt    # plumbing smoke test
pytest -q                                           # once tests exist
```

Note: ASE's EMT is a metals-only toy potential (Cu, Al, Ag, Au, Ni, Pd, Pt).
Smoke-test fixtures must therefore be a **Cu fcc** synthetic ensemble, not an
ionic compound.

## Conventions

- Python ≥ 3.11, type hints on public functions, numpy-style docstrings.
- Pure functions + a thin argparse CLI. Milestone scripts stay runnable
  standalone; refactor into a `rmc_mlip_phonons/` package only at milestone 3.
- Every physics routine gets a synthetic-data unit test (follow the
  wrap-around circular-mean pattern).
- Keep `ROADMAP.md` checkboxes and `CHANGELOG.md` current — update both at the
  end of every completed task, in the same commit.

## Immediate tasks — Milestone-1 hardening (do in order, stop after each)

1. **Test suite.** `tests/` with pytest: rmc6f parser fixtures (with/without
   site-id columns, bracketed labels, `Cell (Ang/deg)` vs `Lattice vectors`
   headers), circular-mean wrap-around, `--ref` CIF site assignment,
   `auto_dim`, mixed-occupancy majority reduction + warning.
2. **Fixture generator.** `tests/fixtures/make_synthetic_ensemble.py`:
   writes N noisy rmc6f configs of a Cu fcc supercell (known answer:
   3 acoustic branches, ω → 0 at Γ).
3. **EMT end-to-end.** Run the full script on the Cu fixture with
   `--calc emt`; fix any ASE/phonopy API friction until it produces
   `band.yaml`/`relaxed.cif`/`summary.json`. Add this as an integration test
   (marked `slow`).
4. **Real-deps run.** Install requirements; run MACE-MP-0 `--model small` on
   the same fixture. Assert: dynamically stable (min band frequency
   > −0.05 THz), acoustic branches vanish at Γ.
5. **Ensemble ergonomics.** `--max-configs` and `--stride` subsampling for
   500+ config ensembles; deterministic with `--seed`.
6. Update `ROADMAP.md` + `CHANGELOG.md`.

## Milestone 2 (only after M1 is green)

`md_run.py`: MLIP-MD at the experimental temperature (ASE Langevin, then
NPT), trajectory → G(r) and S(Q) for closure against measured data, and
temperature-effective force constants → `band_T.yaml`. Design notes first in
`docs/milestone2-plan.md`; get the plan reviewed before implementing.

## Milestone 3 (design doc before code)

`docs/verdicts-schema.md` — the sidecar contract consumed by the viewer:
per-mode `{q, branch, omega_THz, amplitude_ratio, well: single|double,
barrier_meV, verdict: dynamic|static|mixed, confidence}`. Then
`hiphive_fit.py` (RMC displacements + MLIP forces → effective FCs → bands),
mode-projection amplitudes vs (ħ/2ω)coth(ħω/2k_BT), and E(Q) mode mapping.
