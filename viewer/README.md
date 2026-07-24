# viewer/ — browser front end

Static, client-side viewer for the band and mode files this pipeline emits.
Nothing is uploaded; nothing is fetched from a network at runtime.

```bash
npm install
npm run dev        # http://localhost:5178
npm run validate   # numerical contract tests (units, eigenvector mapping, BZ, DOS)
npm run build      # -> dist/, static
```

Load any phonopy-standard band yaml — `band.yaml`, `band_T.yaml`,
`band_rmc.yaml`, `modes_irrep.yaml` — by drag-and-drop, the file picker, or
the `?load=<url>` deep link. Click a point on the band plot to animate that
mode in 3D.

## What it shows

- **Band structure** — SVG plot, drag to zoom, click to select a (q, branch).
- **3D mode** — Three.js animator. Displacements follow phonopy's convention:
  eigenvectors are unit-norm for the *mass-weighted* dynamical matrix, so the
  real-space motion is `Re[(e/√m)·exp(i(q·(n+τ) − ωt))]`. The Bloch phase uses
  the full fractional position `n + τ`, which is why a staggered X-point
  pattern renders correctly when the cell is tiled.
- **Mode details** — per-element participation in the selected mode.
- **Simulated INS** — powder S(|Q|, E) and phonon DOS from the loaded bands,
  computed in a worker.

Frequencies are converted THz → meV on load (`src/io/viewermodel.js`); a
`frequency_unit` key in the file overrides the default.

## Provenance

Vendored 2026-07-23 from
[`drthyang/rmc-phonon-dynamics`](https://github.com/drthyang/rmc-phonon-dynamics)
(MIT, © 2026 Tsung-Han Yang) — `web/src` at the time of the split.

Deliberately **not** carried over: the RMC covariance route
(`math/symmetrize.js`, `math/symmetry.js`, `math/cells.js`, `compute/engine.js`,
`compute/pipeline.js`, `compute/Sk_kernel.wgsl`) and the RMC dataset shell
(`pages/`, `io/readers.js`, `io/sqgr.js`). This pipeline derives frequencies
from MLIP forces rather than by inverting displacement covariance, so that
code has no role here. `src/App.jsx` is new.

## Boundary

The Python pipeline must never depend on this directory — every milestone runs
without node installed. Data flows one way: pipeline → files → viewer.
