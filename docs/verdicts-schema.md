# Milestone 3 design — verdicts.json contract and mode-projection method

*Status: **draft for review**; see also the 2026-07-23 scope pivot in README —
the three-component ratio below is the total-scattering-only route to a
question S(Q,E) answers directly, and its per-mode noise floor is still
unvalidated (see the control experiment in ROADMAP). The sidecar is consumed
by `viewer/`; per CLAUDE.md it is a new file — `band.yaml` is not mutated.
Ground truth and mode patterns: Yang et al., PRR 4, 033123 (2022) — the
displacement patterns are taken **directly from the paper's supplemental
Tables IV–X** (user decision, 2026-07-19), not regenerated from group theory.
A Bilbao/spglib re-derivation is kept only as a *sanity check* that the
ingested patterns reproduce Table II amplitudes when applied to
`data/GTS_P-421m_LT.cif`.*

## What a verdict is

For each phonon mode considered, compare the **measured** ensemble amplitude
of that mode (projected from the RMC configs — model-free, no relaxation,
experimental reference positions) against the **quantum-harmonic
expectation** at the measurement temperature:

    ⟨Q²⟩_qh(ω, T) = (ħ / 2ω) · coth(ħω / 2k_B T)      [normal-coordinate units]

At 5 K every mode in GaTa₄Se₈ is zero-point dominated (coth → 1), so the
yardstick is ħ/2ω — weakly sensitive to the ω source. The amplitude ratio

    r = ⟨Q²⟩_measured / ⟨Q²⟩_qh

drives the verdict: r ≈ 1 → **dynamic** (pure quantum motion), r ≫ 1 with the
excess localized in a symmetry channel → **static** (frozen distortion),
intermediate or channel-ambiguous → **mixed**.

## Method decisions (accumulated through M1–M2 review; binding)

1. **Model-free projection.** Displacements u = (RMC positions − experimental
   cubic reference); reference = the ensemble's own circular-mean structure
   symmetrized to F-4̄3m (never MLIP-relaxed positions). Patterns = the
   supplemental Tables IV–X displacement vectors (X5, X3, W4, Δ1, Γ1, Γ3).
2. **Star pooling (powder degeneracy).** All arms of each k-star and all
   domain orientations are pooled; verdicts never assign an arm. Static order
   in a multi-domain box appears as excess *variance* (and per-config
   coherence), not as an ensemble-mean pattern — the grand mean cancels by
   construction.
3. **Per-config statistics.** The observable is the distribution of
   per-config projected amplitudes {Q(c)}, c = 1..493 converged configs
   (duplicate chains 105/206, 439/404, 459/208 down-weighted; partial runs
   excluded). Plus a local/windowed projection to detect intra-box domains
   (coherence length), since global projections undercount mixed-arm order.
4. **RMC noise calibration — the stiff-mode ruler.** RMC single-atom moves
   inflate site variance. The highest-ω optical modes (≳ 30 meV) must be
   pure zero-point at 5 K; any measured excess there IS the RMC noise floor.
   All amplitude ratios are reported both raw and noise-corrected using this
   internal ruler.
5. **ω sources, in order of preference:** measured GDOS constraints → the
   paper's QE-PBE frequencies → cubic MACE null model. The source is recorded
   per mode (`omega_source`); at 5 K the verdict threshold moves little
   between sources (coth ≈ 1).
6. **Well/barrier fields are Tier-2 statements** (SOC caveat, see
   ROADMAP M4): E(λ) scans along the mode pattern on the declared surface
   (currently MACE null → single wells by construction; meaningful
   double-well claims only after SOC-level fine-tuning). Fields are nullable
   and carry their `surface` provenance.

## verdicts.json schema (v0.1)

Core per-mode fields exactly as CLAUDE.md prescribes; extensions optional.

```jsonc
{
  "schema_version": "0.1",
  "material": "GaTa4Se8",
  "temperature_K": 5.0,
  "reference": {                       // displacement reference, experimental
    "spacegroup": "F-43m",
    "a_A": 10.3563,
    "source": "ensemble circular mean, symmetrized (no relaxation)"
  },
  "provenance": {
    "ensemble": {"n_configs": 493, "dedup": [...], "excluded": [...]},
    "patterns": "PRR 4.033123 supplemental Tables IV-X",
    "code": {"repo": "...", "commit": "...", "script": "mode_project.py"},
    "noise_ruler": {"modes_used": [...], "excess_ratio": 1.7}
  },
  "modes": [
    {
      // ---- required (CLAUDE.md contract) ----
      "q": [0.5, 0.5, 0.0],            // band.yaml (primitive) convention
      "branch": 7,
      "omega_THz": 3.15,
      "amplitude_ratio": 6.8,          // noise-corrected r
      "well": "single",                // "single" | "double" | null
      "barrier_meV": null,
      "verdict": "static",             // "static" | "dynamic" | "mixed"
      "confidence": 0.91,              // bootstrap x omega-source spread, 0-1
      // ---- optional extensions ----
      "irrep": "X5",
      "star": [[0,1,0],[1,0,0],[0,0,1]],   // conventional-cell arms pooled
      "omega_source": "gdos",              // "gdos" | "qe-dft" | "mace-null"
      "amplitude_ratio_raw": 11.2,         // before noise correction
      "static_amplitude_A": 0.096,         // per-atom-scale, for intuition
      "target_amplitude_A": 0.1196,        // published refined value (Table II)
      "coherence_cells": 6.5,              // windowed-projection length scale
      "well_surface": "mace-mp0-small",    // provenance of well/barrier
      "notes": "..."
    }
  ]
}
```

Conventions: q in the **primitive reciprocal basis of band.yaml** (so the
viewer can badge the dispersion directly); star arms additionally given in
conventional-cell notation for human reading. Amplitudes in the normal-mode
convention are recorded internally; user-facing fields are per-atom Å scales.
Units: THz, meV, Å, K (stated per field in the JSON schema file).

## Deliverables and order of work

1. `reference/gts_mode_patterns.json` — ingested supplemental Tables IV–X
   (small, published-derived, committed; new top-level `reference/` dir for
   such files). Ingest script + the Table-II sanity check against
   `data/GTS_P-421m_LT.cif`.
2. `mode_project.py` — model-free projection engine: reference build, star
   pooling, per-config distributions, windowed/local projections, stiff-mode
   noise ruler. Pure functions + CLI, unit-tested on synthetic ensembles.
3. `verdicts.py` — assembles `verdicts.json` (thresholds, confidence,
   optional E(λ) wells on a declared surface).
4. `hiphive_fit.py` (CLAUDE.md M3 item) — RMC displacements + MLIP forces →
   effective FCs → bands; now explicitly Tier-1/Tier-2-scoped per the SOC
   tiering in ROADMAP.
5. E(Q) mapping: projected-mode energies vs the measured GDOS features
   (the 13-meV story), using the SEQUOIA resolution polynomial.

## Validation (synthetic, before real data)

- **Dynamic-only ensemble**: M2 quantum sampling at 5 K → all verdicts must
  read *dynamic* with r ≈ 1 (closes the loop with the M2 machinery).
- **Static injection**: add a frozen X5-pattern of known amplitude (with
  random arms/domains per config) on top of quantum sampling → verdict
  *static* in the X5 channel only, recovered amplitude within errors.
- **Mixed**: partial freezing → *mixed*, with confidence responding to the
  injected fraction.
- **Noise ruler**: inflate per-atom white noise → raw ratios rise everywhere,
  corrected ratios stay put.

## Open questions for review

1. Supplemental Tables IV–X: waiting on the local file (APS blocks fetching).
   Preferred drop location: `data/` (gitignored) or `~/Downloads`.
2. Verdict thresholds: propose dynamic r < 1.5, static r > 3 (after noise
   correction) — tune on the synthetic suite?
3. Should `verdicts.json` also carry the *irrep-level* summary (X5/X3/W4/Δ1
   aggregate rows) alongside per-(q, branch) rows, for the viewer's badge
   panel? (Proposed: yes, as `"summary"` block.)
4. Down-weighting duplicates vs dropping: propose dropping 206/404/208
   (keeping 105/439/459) — simpler provenance.
5. The windowed/local projection window size (propose 2–4 conventional cells,
   scanned) — sensitivity study on synthetics first.
