# Changelog

## [Unreleased]
### Added
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
