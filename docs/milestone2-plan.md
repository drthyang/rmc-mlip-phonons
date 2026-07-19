# Milestone 2 design — finite-temperature closure and effective force constants

*Status: **draft for review** (CLAUDE.md gates M2 implementation on this
review). Author: pipeline; target dataset: GaTa₄Se₈ 5 K RMC ensemble.*

## Goal

Close the loop between the MLIP model and the measured data, and produce
temperature-effective force constants:

1. **Closure**: simulate G(r) and S(Q) from MLIP-sampled configurations and
   compare quantitatively against the measured `scale_ft_rmc.fq`
   (neutron S(Q), Q = 0.8–27 Å⁻¹, ΔQ = 0.01 Å⁻¹, 2619 points).
2. **`band_T.yaml`**: effective (temperature-renormalized) phonon bands at the
   experimental temperature, as a **new sidecar** — `band.yaml` stays frozen.

Non-goals for M2: verdicts/static-dynamic classification (M3), NAC/Born
charges (optional M3+), viewer changes (other repo), DFT fine-tuning (M4).

## Physical design decisions (the ones review should challenge)

### D1 — 5 K is a quantum problem: sampling method

At 5 K, k_BT ≈ 0.43 meV ≈ 0.1 THz — essentially every phonon in this material
(lowest optical 1.3 THz, acoustic zone-boundary ~1–2 THz) is in its ground
state. Mode amplitudes follow ⟨u²⟩ ∝ (ħ/2ω)·coth(ħω/2k_BT) → ħ/2ω (zero-point
motion), while **classical MD gives equipartition k_BT/ω²** — wrong by orders
of magnitude for the high-ω modes at this temperature. A classical Langevin
trajectory at 5 K would produce far-too-sharp G(r) peaks and a bad closure.

**Decision:** primary sampler = **harmonic quantum sampling**: random
supercell displacements drawn from the M1 phonopy model with quantum
statistics at T = 5 K (phonopy's random-displacement machinery includes
zero-point occupation). MLIP-MD (ASE Langevin) is kept in `md_run.py` as an
explicitly *classical* cross-check and for future higher-T datasets, with its
5 K caveat printed loudly.

### D2 — lattice for closure: experimental, not MLIP-relaxed

MACE-MP-0 (small) overexpands a by +1.6 % (10.526 vs 10.356 Å). G(r) peak
*positions* scale with a, so closure must be computed **at the experimental
lattice** (fixed cell, a = 10.3563 Å; internal coordinates re-relaxed under
that constraint). The pure MLIP-relaxed structure remains the M1 theory-side
reference. Both provenances recorded in the output JSON.

### D3 — effective-FC engine: hiPhive, introduced in M2, reused in M3

`band_T.yaml` needs a displacement/force → force-constant fit. Proposal: use
**hiPhive** as the single fitting engine for both M2 (fit on our own
quantum-sampled snapshots + MLIP forces — a self-consistency loop with known
provenance) and M3 (same machinery, but RMC snapshots as the displacement
source — the experiment-constrained fit this project exists for). At 5 K the
M2 renormalization should be tiny (band_T ≈ band); *that agreement is itself a
validation* of the sampling + fitting machinery before M3 trusts it.

### D4 — S(Q)/G(r) computation: small native implementation

Partial pair histograms g_ij(r) → coherent-neutron-weighted total G(r) →
windowed Fourier transform to S(Q) on exactly the measured Q-grid. Native
(numpy) rather than an external PDF package: keeps deps light and every step
unit-testable against analytic cases, matching repo conventions. Neutron
b_coh from standard tables (Ga 7.288, Ta 6.91, Se 7.97 fm — confirm source).
Closure metric: Rw over the measured Q-window (same definition RMCProfile
uses, so numbers are comparable to the RMC fit's own residual).

## Deliverables

```
md_run.py            # CLI: sample (quantum) | md (classical Langevin/NPT)
                     #      → snapshots → forces → G(r), S(Q), band_T.yaml
m2_out/
  gr_sim.dat         # r (Å), G(r) — simulated, neutron-weighted
  sq_sim.dat         # Q (Å⁻¹), S(Q) — on the measured grid
  closure.json       # Rw, grids, sampling provenance, lattice used
  band_T.yaml        # hiPhive effective FCs → phonopy band structure at T
```

Units: Å, eV, THz throughout (phonopy defaults); stated in every docstring.
All MLIP evaluations `default_dtype="float64"`.

## Validation / tests (same discipline as M1)

- Unit: g(r) histogram on an ideal fcc lattice (analytic peak positions);
  neutron weighting closure (weights sum to ⟨b⟩²-normalized 1); FT round-trip
  on a synthetic Gaussian; quantum ⟨u²⟩ of generated snapshots vs the analytic
  coth expression per mode; Rw metric on known inputs.
- Integration (slow): Cu/EMT end-to-end sample → G(r)/S(Q) → band_T; assert
  band_T ≈ band at low T and G(r) first peak at a/√2.
- Real-data acceptance: GaTa₄Se₈ closure Rw reported and compared against the
  RMC fit's own residual; not a pass/fail gate (MLIP systematics expected),
  but recorded in `closure.json`.

## Resolved in review (2026-07-19)

1. **Measured data**: `scale_ft_rmc.fq` is STOG-generated. Verified high-Q
   level ≈ 0 (mean +0.0004 over the last 500 points, oscillations decaying)
   and low-Q → −0.87, i.e. the convention is **F(Q) = S(Q) − 1**
   (Faber–Ziman). The RMC fit ran with **scale and offset free**, so the
   closure comparison fits (and reports) the same two parameters —
   shape-plus-scale Rw, mirroring the RMC treatment.
2. **G(r) and S(Q) are the same measurement** (r-space = FT of this F(Q)).
   Primary Rw is quoted in Q-space; G(r) shown for interpretation only.
3. **The 5 K phase is NOT cubic** (user-confirmed lower symmetry). This is
   *the physics target of the whole project*, and it reframes M2–M3:
   - The cubic F-4̄3m MLIP model (M1) is the **null model**, not the answer.
   - The cubic RMC box suppresses global strain, so the distortion lives in
     **internal displacement patterns** — consistent with the residual P1
     signal at symprec 10⁻³ surviving the 493-config average.
   - **Closure becomes a discriminator, not a validation**: the residual
     between quantum-sampled-cubic F(Q) and the measured F(Q) *contains the
     distortion signature*. A perfect fit is not expected and not the goal.
   - M3 verdicts then decompose it: which modes carry static (frozen)
     displacement excess vs quantum-dynamic amplitudes.

## Still open for review

4. Neutron b_coh values/source to standardize on
   (proposal: Sears 1992 tables — Ga 7.288, Ta 6.91, Se 7.97 fm).
5. Is the classical-MD leg worth running at 5 K at all beyond the caveat demo,
   or reserve MD for a future higher-T dataset?
6. hiPhive cutoffs/order for the effective fit (start: 2nd order only,
   cutoff ~6 Å, 2×2×2 supercell of the conventional cell)?
7. Candidate low-T space group(s) from literature / your diffraction, and
   whether the RMC fit co-fitted a Bragg profile — both shape the M3
   symmetry-mode analysis (see questions posed in review).
