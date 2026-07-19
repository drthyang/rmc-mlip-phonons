# Changelog

## [Unreleased]
### Fixed
- **Phonopy used the wrong symmetry tolerance.** `phonopy_bands` now threads
  `--symprec` into `Phonopy(...)`; phonopy's default (1e-5) is tighter than the
  ~fmax residual asymmetry of a relaxed cell, so `primitive_matrix="auto"` read
  fcc as P1 and produced folded bands (12 modes for Cu instead of 3). With the
  fix the fcc primitive is found: 3 acoustic branches, ω → 0 at Γ.
- Migrated off phonopy's deprecated `get_band_structure_dict()` to the
  `band_structure` property (phonopy ≥4).

### Added
- `tests/test_emt_end_to_end.py`: EMT full-pipeline integration test (marked
  `slow`) on the synthetic fcc-Cu ensemble — asserts all three outputs are
  written, `relaxed.cif` is cubic Cu, exactly 3 branches vanish at Γ, and no
  imaginary modes. First run of the ase/spglib/phonopy leg against installed
  deps (phonopy 4.4, ase 3.29, spglib 2.7). `pytest.ini` gains the `slow`
  marker and filters for two benign library deprecations.
- `tests/fixtures/make_synthetic_ensemble.py`: deterministic generator of N
  noisy fcc-Cu `.rmc6f` configs (known answer: 3 acoustic branches, ω → 0 at
  Γ). Runnable standalone or imported; `build_supercell` + `make_fcc_cu_ensemble`.
  Tests verify parseability, that folding recovers the ideal 4-site fcc basis,
  seed-determinism, and the noise-free exact limit.
- `tests/` pytest suite (21 tests) for `milestone1_bands`: rmc6f parser
  fixtures (with/without site-id columns, bracketed labels, `Cell (Ang/deg)`
  vs `Lattice vectors` headers, ionic-label element cleaning, missing-section
  errors), circular-mean wrap-around + supercell folding, `--ref` CIF site
  assignment, `auto_dim` supercell heuristic, and mixed-occupancy majority
  reduction with reporting. `tests/conftest.py` makes the standalone script
  importable and provides rmc6f/config builders; `pytest.ini`.
- `.gitignore` covering ensembles, results, `m1_out*/`, venv, and model weights.

## [0.1.0] — 2026-07-19
### Added
- Repo scaffold: README, ROADMAP, CLAUDE.md, requirements.
- `milestone1_bands.py`: rmc6f ensemble → folded/circular-averaged unit cell →
  spglib symmetrization → MLIP relaxation (MACE-MP-0, float64) → phonopy
  finite displacements → `band.yaml` + `relaxed.cif` + `summary.json`.
- Parser and circular averaging verified on synthetic wrap-around data.
