# Concept note — model-space RMC / "dynamic EPSR"

*Status: concept (not a scheduled milestone). Captured 2026-07-19 from
discussion; owner: Tsung-Han Yang. This is the long-range vision the
milestone pipeline builds toward; nothing here changes M2/M3 scope.*

## The idea

Invert the RMC philosophy. Standard RMC moves **atoms** and scores *static*
snapshots; the accepted configuration ensemble is the model, and static vs
dynamic disorder stay conflated. Here instead the Monte Carlo variable is a
**model** — an equilibrium (P1) structure whose *dynamics* the MLIP supplies
— and every observable is an **ensemble average over that model's thermal
motion**:

1. Perturb the model (equilibrium sites / distortion pattern).
2. Generate its thermal ensemble with the MLIP — phonons + **quantum**
   sampling (zero-point correct; exactly the `md_run.py` machinery).
3. Compute ensemble-averaged G(r), S(Q), **Bragg**, and — with the measured
   INS data — **S(Q, ω)**.
4. Joint χ² over all channels → Metropolis accept/reject (or gradient /
   Bayesian updates where differentiable).

Converged result: an equilibrium structure **plus** a phonon spectrum that
jointly reproduce diffraction and spectroscopy. Static vs dynamic disorder
is separated **by construction**: the refined equilibrium structure carries
the static part; the sampled dynamics carries the rest. This is the clean
resolution of exactly the ambiguity the GaTa₄Se₈ 5 K problem exhibits.

## Two flavors, one loop

- **(a) Structure refinement, fixed potential** — MC over mean positions /
  order-parameter amplitudes; the MLIP is trusted for dynamics. Cheap(er),
  and sufficient where the potential is adequate.
- **(b) Potential refinement** — the force constants / MLIP parameters join
  the refined variables. INS data forces this flavor: if the structure is
  right but phonon energies are wrong, only the model can move. For GTS
  this is ultimately required (MACE-MP-0 hosts no correlation/SOC-driven
  distortion well; see docs/milestone2-plan.md probe results), via
  fine-tuning or a Δ-correction refined in the loop.

## Literature map (searched 2026-07-19; pieces exist, package does not)

| Method | Shares | Lacks |
|---|---|---|
| [EPSR (Soper)](https://www.tandfonline.com/doi/abs/10.1080/00268970110056889) / [Dissolve](https://github.com/disorderedmaterials/dissolve) | The philosophy: refine the *model*, let simulation generate the ensemble that fits S(Q) | Liquids/glasses; classical; no phonons/Bragg/INS |
| [DiffTRe (Thaler & Zavadlav 2021)](https://www.nature.com/articles/s41467-021-27241-4) / [chemtrain](https://www.sciencedirect.com/science/article/pii/S0010465525000153) | Gradient version of the loop; the reweighting cost trick | No crystal-scattering channels |
| [Force coefficients from X-ray TDS (2026)](https://arxiv.org/abs/2603.28683) (Wehinger/Bosak lineage) | Inverse problem: dynamics from scattering via model refinement | Single-crystal TDS only; harmonic; no structure/PDF/Bragg joint fit |
| [Phonons-from-powder-diffraction RMC (Goodwin/Tucker/Dove/Keen)](https://iopscience.iop.org/article/10.1088/0953-8984/19/33/335218/meta) | Dynamics from this data class | Forward-free covariance route; noisy; see [limits paper](https://arxiv.org/pdf/cond-mat/0209540) |
| [Hybrid RMC + MLIP in RMCProfile (2024)](https://onlinelibrary.wiley.com/iucr/doi/10.1107/S1600576724009282) | MLIP inside an RMC loop | MLIP is only an energy prior; observables still static |
| [uMLIP–INS benchmark (Han & Cheng 2025)](https://arxiv.org/abs/2506.01860), [ML unfolding of powder INS](https://arxiv.org/pdf/2404.13507) | MLIP phonons vs INS | One-shot forward comparison, no refinement loop |

**Gap = this idea**: EPSR-style refinement for crystals, ensemble generated
by MLIP lattice dynamics with quantum statistics, χ² jointly over total
scattering + Bragg + S(Q, ω). "Dynamic EPSR" / "model-space RMC".

## Why it is tractable (and where it bites)

- **Cost.** Naive per-atom moves × full phonon recompute per move is
  unaffordable. Rescues (compatible, cumulative):
  1. **Collective move space** — for GTS: the P-4̄2₁m order-parameter
     amplitude(s), a handful of coordinates instead of 26k positions.
  2. **Reweighting** — most trial moves rescore the previous ensemble with
     perturbation weights (DiffTRe trick); fresh sampling only on accepted
     drift beyond a trust radius.
  3. Gradient/Bayesian optimization where the chain is differentiable;
     surrogate caching of χ²(model).
- **Identifiability.** Diffraction constrains dynamics only through
  amplitudes; **the INS data is what makes the inverse problem well-posed**
  (ω's pinned directly; PDF/Bragg pin structure). Without INS, flavor (b)
  is under-determined.
- **Quantum statistics are mandatory at 5 K** (zero-point dominates —
  established in milestone 2); classical MD inside the loop would corrupt
  every channel.

## Relation to the milestone pipeline

- M2 built the **forward evaluator** (structure → quantum ensemble →
  G(r)/F(Q) → χ² with scale+offset). This idea wraps an outer refinement
  loop around it. Missing forward pieces: Bragg profile and S(Q, ω)
  calculators.
- M3's mode-projected verdicts are the **linearized first pass** of the
  same question; this note is the full nonlinear inverse. Do M3 first —
  its outputs (which modes carry static excess) initialize the move space.
- M4's fine-tuning is flavor (b)'s enabler.

## First demonstration target

GaTa₄Se₈ 5 K: move space = the **six irrep amplitudes of Yang et al.,
PRR 4, 033123 (2022), Table II** — X5 (0.120 Å), X3 (0.072 Å), W4
(0.026 Å), Δ1 (0.021 Å), Γ1, Γ3 — i.e. the published static mode
decomposition becomes the refinement coordinates, with MLIP dynamics on
top and χ² = F(Q) + Bragg (superlattice intensities) + INS (SEQUOIA
Ei = 60 meV; instrument resolution polynomial in the paper's Appendix E).
Few-dimensional, all data in hand, and the paper provides the static
answer to initialize and benchmark against — the ideal demonstration. The
paper's single-crystal refinement is effectively the *static* version of
this loop; the addition is joint dynamics (quantum ensemble) + total
scattering + INS self-consistency.

## Open questions

1. ~~INS data: instrument, form?~~ **Answered (2026-07-19): powder S(Q,E)
   from SEQUOIA and ARCS (SNS direct-geometry TOF), with reduced phonon
   DOS — for the whole lacunar-spinel series, not just GaTa₄Se₈.**
   Consequences:
   - The cheap first INS channel is the **neutron-weighted GDOS**
     (incoherent approximation, σ_s/M-weighted partial DOS) — our forward
     model produces it in minutes; series-wide MACE benchmarking becomes a
     methods-paper pillar (V/Nb/Ta × S/Se chemical trends).
   - The full-χ² channel later is coherent one-phonon powder S(Q,E)
     (Euphonic / OCLIMAX-style), resolution-convolved per instrument.
   - Remaining detail needed per dataset: incident energies, T points, and
     whether the DOS reductions are multiphonon/multiple-scattering
     corrected.
2. Bragg forward model: profile convolution parameters to mirror the
   diffractometer (reuse the Rietveld TOF profile from the user's fits?).
3. Acceptance temperature / error model for the joint χ² (channel weights).
4. Where flavor (b) parameters live: hiPhive FC corrections vs MLIP
   fine-tune vs Δ-model on top of MACE.
