# rmc-mlip-phonons — Development Roadmap

*Last updated: 2026-07-19. Compute pipeline companion to
[`rmc-phonon-dynamics`](https://github.com/drthyang/rmc-phonon-dynamics)
(the browser viewer). Contract: this repo emits `band.yaml`, `relaxed.cif`,
`summary.json`, and later `verdicts.json`; the viewer loads them.*

---

## Done

- [x] **Repo scaffold** — README, requirements, `CLAUDE.md`, and
      `milestone1_bands.py` (rmc6f → averaged unit cell → MLIP relax →
      phonopy `band.yaml` + `relaxed.cif` + `summary.json`). Parser and
      circular averaging unit-verified on synthetic wrap-around data.

## In progress — Milestone 1 hardening

- [x] pytest suite: parser fixtures, circular mean, `--ref` CIF assignment,
      `auto_dim`, mixed-occupancy handling (`tests/`, 21 tests)
- [x] synthetic Cu fcc ensemble generator (`tests/fixtures/`)
- [x] EMT end-to-end integration test; fix ASE/phonopy API friction
      (threaded `--symprec` into phonopy so `primitive_matrix="auto"` finds the
      fcc primitive; migrated off deprecated `get_band_structure_dict()`)
- [x] first real MACE-MP-0 run on the fixture (stability + Γ acoustics
      asserts) — surfaced + fixed the 0.01 Å MLIP-noise imaginary-mode artifact
      (default displacement now 0.03 Å). Target-dataset run pending real data.
- [ ] `--max-configs` / `--stride` / `--seed` for 500+ config ensembles

## Milestone 2 — finite temperature

- [ ] `docs/milestone2-plan.md` (design first)
- [ ] `md_run.py`: MLIP-MD at experimental T (Langevin → NPT)
- [ ] closure: G(r) / S(Q) from MD vs measured data
- [ ] temperature-effective force constants → `band_T.yaml`

## Milestone 3 — experiment-constrained FCs + verdicts

- [ ] `docs/verdicts-schema.md` — sidecar contract with the viewer
- [ ] `hiphive_fit.py`: RMC displacement snapshots + MLIP forces →
      effective FCs → bands
- [ ] mode-projection amplitudes vs (ħ/2ω)coth(ħω/2k_BT); E(Q) mode mapping
- [ ] emit `verdicts.json` (per-mode static/dynamic badges)
- [ ] viewer-side overlay + badge panels (tracked in `rmc-phonon-dynamics`)

## Milestone 4 — application & methods paper

- [ ] real datasets end-to-end; error budget of covariance vs MLIP bands
- [ ] DFT fine-tuning / active learning where the foundation model falls short
- [ ] manuscript draft
