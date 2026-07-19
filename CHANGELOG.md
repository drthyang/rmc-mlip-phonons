# Changelog

## [Unreleased]
### Added
- `md_run.py` (milestone 2): `sample` mode draws random supercell snapshots
  from the phonopy model with **quantum statistics** (zero-point included;
  phonopy's finite-temperature random displacements — the force-constant
  cache is saved/restored around the dataset swap, phonopy ≥4 invalidates
  it). `md` mode is a classical Langevin NVT cross-check that warns loudly
  below 100 K. Native closure machinery: species-partial g_ij(r) (chunked
  minimum-image histograms), Faber–Ziman neutron weighting (Sears 1992
  b_coh), sine transforms G(r) ↔ F(Q)=S(Q)−1 on the measured STOG grid, and
  a scale+offset Rw fit mirroring the RMC convention. Optional hiPhive
  2nd-order effective-FC fit → `band_T.yaml` sidecar (band.yaml untouched).
  Outputs: `closure.json`, `gr_sim.dat`, `sq_sim.dat`, `relaxed_expt.cif`,
  `band_T.yaml`.
- Tests: 10 analytic unit tests (fcc first-shell + coordination, g(r)
  normalisation identity, FT round-trip, Bragg-position check, weight
  closure, scale/offset fit, .fq parser) and 6 slow integration tests on
  Cu/EMT — including the D1 physics guard: sampled MSD matches phonopy's
  quantum ⟨u²⟩ at 5 K and 300 K, and is zero-point-dominated at 5 K.
  Suite: 66 passed.
- `requirements.txt`: + `hiphive>=1.5`.
- First M2 acceptance run on the GaTa₄Se₈ 5 K dataset (results git-ignored):
  quantum u_rms = 0.040 Å/component vs the RMC ensemble's 0.095 — the
  displacement-variance gap M3 will decompose (with the caveat that RMC
  single-atom moves decorrelate neighbours, inflating site variance while
  still fitting pair widths). Baseline closure Rw(Q) = 0.74, Rw(r) = 0.58;
  measured G(r) peaks systematically broader than the quantum-harmonic null
  model. band_T tracks band to 0.18 THz at 5 K.
- `docs/milestone2-plan.md`: milestone-2 design draft (awaiting review).
  Decisions put forward: harmonic **quantum** sampling instead of classical
  Langevin MD at 5 K (zero-point-dominated regime), closure computed at the
  experimental lattice (MACE +1.6 % would shift G(r) peaks), hiPhive as the
  shared M2/M3 effective-FC engine, native unit-tested G(r)/S(Q) with Rw
  against the measured `scale_ft_rmc.fq` grid.
- `--skip-nonconverged`: drop input configs whose rmc6f header reports
  `Number of moves generated: 0` (RMC runs that never started). Cheap
  header-only scan (`read_moves_generated` / `drop_nonconverged`); files
  without the header line are kept, so synthetic ensembles pass through.
  Drop count recorded in `summary.json["sampling"]`. 6 unit tests.
- First real-dataset runs (GaTa₄Se₈ 5 K RMC ensemble, 26,624-atom configs):
  single converged config and all 493 converged configs both produce a
  dynamically stable F-4̄3m band structure at the same MACE-MP-0 minimum
  (max band Δω 0.0009 THz); the ensemble average recovers F-4̄3m at 10×
  tighter symprec than a single config. Not committed (data/ and results/
  are git-ignored); see ROADMAP for caveats (lattice +1.6 %, no NAC yet).

## [0.2.0] — 2026-07-19

Milestone-1 hardening. The pipeline now runs end-to-end against the installed
ase / spglib / phonopy / MACE stack, with a full test suite and two physics
fixes that the first real runs surfaced.

### Added
- `.gitignore` covering ensembles, results, `m1_out*/`, venv, and model weights.
- `tests/` pytest suite for `milestone1_bands`: rmc6f parser fixtures
  (with/without site-id columns, bracketed labels, `Cell (Ang/deg)` vs
  `Lattice vectors` headers, ionic-label element cleaning, missing-section
  errors), circular-mean wrap-around + supercell folding, `--ref` CIF site
  assignment, `auto_dim` heuristic, and mixed-occupancy majority reduction.
  `tests/conftest.py` makes the standalone script importable; `pytest.ini`.
- `tests/fixtures/make_synthetic_ensemble.py`: deterministic generator of N
  noisy fcc-Cu `.rmc6f` configs (known answer: 3 acoustic branches, ω → 0 at
  Γ). Runnable standalone or imported. Tests verify parseability, that folding
  recovers the ideal 4-site fcc basis, seed-determinism, and the noise-free
  exact limit.
- `tests/test_emt_end_to_end.py` (marked `slow`): EMT full-pipeline integration
  test — asserts all three outputs are written, `relaxed.cif` is cubic Cu,
  exactly 3 branches vanish at Γ, and no imaginary modes. First run of the
  ase/spglib/phonopy leg against installed deps (phonopy 4.4, ase 3.29,
  spglib 2.7).
- `tests/test_mace_end_to_end.py` (marked `slow` + `mace`): real-dependency
  integration test running MACE-MP-0 (small, float64) on the fcc-Cu fixture —
  no imaginary modes (min > −0.05 THz), 3 branches vanishing at Γ, float64.
  Auto-skips when mace-torch is absent (first run downloads/caches the ~31 MB
  MACE-MP-0 weights).
- Ensemble subsampling for large (500+ config) runs: `--stride` (ordered
  decimation), `--max-configs` (seeded random cap), and `--seed`, via the pure
  `select_configs` helper. Deterministic; recorded in `summary.json["sampling"]`.
- `pytest.ini`: `slow` / `mace` markers and filters for benign spglib / phonopy
  / torch deprecation warnings.

### Changed
- Default `--displacement` raised from 0.01 to **0.03 Å**. phonopy's 0.01 Å is
  DFT-tuned; at that distance a universal MLIP's forces sit near the model
  noise floor. In the first real MACE-MP-0 (small) run this produced a spurious
  ~−0.33 THz mode just off Γ for dynamically-stable fcc Cu; 0.03 Å removes it
  (verified: ASR / `symmetrize_force_constants` does not help, displacement does).

### Fixed
- **Phonopy used the wrong symmetry tolerance.** `phonopy_bands` now threads
  `--symprec` into `Phonopy(...)`; phonopy's default (1e-5) is tighter than the
  ~fmax residual asymmetry of a relaxed cell, so `primitive_matrix="auto"` read
  fcc as P1 and produced folded bands (12 modes for Cu instead of 3). With the
  fix the fcc primitive is found: 3 acoustic branches, ω → 0 at Γ.
- Migrated off phonopy's deprecated `get_band_structure_dict()` to the
  `band_structure` property (phonopy ≥4).

## [0.1.0] — 2026-07-19
### Added
- Repo scaffold: README, ROADMAP, CLAUDE.md, requirements.
- `milestone1_bands.py`: rmc6f ensemble → folded/circular-averaged unit cell →
  spglib symmetrization → MLIP relaxation (MACE-MP-0, float64) → phonopy
  finite displacements → `band.yaml` + `relaxed.cif` + `summary.json`.
- Parser and circular averaging verified on synthetic wrap-around data.
