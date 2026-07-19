# Changelog

## [0.1.0] — 2026-07-19
### Added
- Repo scaffold: README, ROADMAP, CLAUDE.md, requirements.
- `milestone1_bands.py`: rmc6f ensemble → folded/circular-averaged unit cell →
  spglib symmetrization → MLIP relaxation (MACE-MP-0, float64) → phonopy
  finite displacements → `band.yaml` + `relaxed.cif` + `summary.json`.
- Parser and circular averaging verified on synthetic wrap-around data.
