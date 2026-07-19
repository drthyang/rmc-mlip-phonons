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

## Done — Milestone 1 hardening

*Green end-to-end on synthetic fcc Cu (EMT + MACE-MP-0) **and on the real
target dataset**: GaTa₄Se₈ 5 K RMC ensemble, single config and all 493
converged configs (F-4̄3m, dynamically stable, 39 branches). `pytest -q` is
50 tests (unit + `slow`/`mace` integration).*

- [x] pytest suite: parser fixtures, circular mean, `--ref` CIF assignment,
      `auto_dim`, mixed-occupancy handling (`tests/`, 21 tests)
- [x] synthetic Cu fcc ensemble generator (`tests/fixtures/`)
- [x] EMT end-to-end integration test; fix ASE/phonopy API friction
      (threaded `--symprec` into phonopy so `primitive_matrix="auto"` finds the
      fcc primitive; migrated off deprecated `get_band_structure_dict()`)
- [x] first real MACE-MP-0 run on the fixture (stability + Γ acoustics
      asserts) — surfaced + fixed the 0.01 Å MLIP-noise imaginary-mode artifact
      (default displacement now 0.03 Å)
- [x] target-dataset run: GaTa₄Se₈ 5 K (`data/`, 8×8×8 = 26,624-atom configs).
      Config 1 and the full 493-converged-config ensemble both → F-4̄3m,
      dynamically stable, same MACE minimum (max band Δω 0.0009 THz); ensemble
      average recovers F-4̄3m at 10× tighter symprec (0.01 vs 0.1). Caveats to
      carry: +1.6 % lattice overexpansion (MACE-MP-0 small), no NAC/LO–TO yet.
      Added `--skip-nonconverged` (drops 0-move configs; found exactly the 8
      bad ones: 0, 141, 161, 170, 413, 419, 442, 492).
- [x] `--max-configs` / `--stride` / `--seed` for 500+ config ensembles
      (`select_configs`, recorded in `summary.json["sampling"]`)

## Milestone 2 — finite temperature

- [x] `docs/milestone2-plan.md` — reviewed over three rounds (quantum sampling
      at 5 K, closure at the experimental lattice, hiPhive engine, powder
      q-star degeneracy, MACE-as-null-model confirmed by probes)
- [x] `md_run.py`: `sample` mode = harmonic **quantum** sampling (zero-point
      correct at 5 K; MSD verified against phonopy's analytic ⟨u²⟩);
      `md` mode = classical Langevin cross-check with low-T warning.
      NPT deferred until a higher-T dataset exists.
- [x] closure: native partial-g_ij(r) → Faber–Ziman G(r) → F(Q)=S(Q)−1 on the
      measured STOG grid, scale+offset fitted (mirrors RMC), Rw in Q and r
- [x] temperature-effective force constants → `band_T.yaml` (hiPhive 2nd
      order; validated band_T ≈ band on Cu/EMT)
- [x] GTS acceptance run (5 K dataset, cubic null model): u_rms(quantum) =
      0.040 Å/comp vs RMC ensemble 0.095 Å/comp (the M3 target gap);
      baseline closure Rw(Q) = 0.74 / Rw(r) = 0.58 (Q-space dominated by
      box-truncation vs instrument peak-shape mismatch — needs resolution
      convolution before the residual is fully interpretable); measured G(r)
      peaks systematically broader than the quantum null model = static/
      anharmonic excess. band_T ≈ band (max |Δω| = 0.18 THz at 5 K) —
      machinery validated on the real system.
- [ ] closure refinement: instrument-resolution / box-size-matched F(Q)
      comparison so the null-model residual isolates the distortion

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
