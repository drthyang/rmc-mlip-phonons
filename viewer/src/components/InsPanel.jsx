import React, { useEffect, useMemo, useRef, useState } from 'react';
import { buildInsData } from '../compute/ins';
import { downloadString } from '../io/writers';
import { CMAP_NAMES, colormap, cmapGradient, intensityTransform } from './sqeColormaps';

/* ── Cobalt theme tokens ───────────────────────────────────────────────── */
const INK = 'var(--ink)', DIM = 'var(--dim)', FAINT = 'var(--faint)';
const ACCENT = 'var(--accent)', BORDER = 'var(--border)', INSET = 'var(--inset)', INSET2 = 'var(--inset2)';
const cardStyle = { background: 'var(--card)', border: `1px solid ${BORDER}`, borderRadius: 10 };
const insetInput = { width: '100%', boxSizing: 'border-box', background: INSET, border: `1px solid ${BORDER}`, borderRadius: 7, padding: '8px 10px', font: "13px 'Space Mono'", color: INK };

const HBAR2_2MN = 2.0723;   // meV·Å² — direct-geometry kinematic constant
const BG = [11, 14, 22];    // masked-pixel colour (matches the dark canvas #0b0e16)
const DOS_STROKE = '#2f6df0', DOS_FILL = 'rgba(47,109,240,0.12)';

// Kinematic |Q| window accessible at energy transfer E for incident energy Ei
// (direct geometry). Shared by the heatmap and the integrated-intensity curve so
// the curve is the exact energy-marginal of the pixels drawn. Semantics match the
// original inline test: [-1, ∞] ⇒ all Q (Ei≤0, direct); [1, -1] ⇒ none (E>Ei).
function qRange(E, Ei, ki) {
  if (Ei <= 0) return [-1, Infinity];
  if (E > Ei) return [1, -1];
  const kf = Math.sqrt(Math.max(0, Ei - E) / HBAR2_2MN);
  return [Math.abs(ki - kf), ki + kf];
}

/**
 * INS (simulated inelastic neutron scattering) panel: powder S(|Q|,E) heatmap +
 * the Q-integrated intensity curve ∫S(|Q|,E)d|Q|, drawn on ONE shared energy axis.
 * S(Q,E)/DOS math runs in io/sqeworker.js (unchanged); the integrated curve is a
 * cheap main-thread energy-marginal of the heatmap the worker returns.
 */
export default function InsPanel({ results, temperature }) {
  const workerRef = useRef(null);
  const canvasRef = useRef(null);
  const [running, setRunning] = useState(false);
  const [out, setOut] = useState(null);
  const [error, setError] = useState(null);

  // Default energy window from the BULK of the spectrum (90th percentile of
  // positive energies) so a few near-zero-eigenvalue outlier modes don't cram
  // the real bands into the bottom.
  const maxE = useMemo(() => {
    const vals = [];
    for (const row of results.bands) for (const v of row) if (isFinite(v) && v > 0) vals.push(v);
    if (!vals.length) return 50;
    vals.sort((a, b) => a - b);
    const p90 = vals[Math.min(vals.length - 1, Math.floor(vals.length * 0.90))];
    return Math.max(5, Math.ceil(p90 * 1.1));
  }, [results]);

  // Eᵢ (incident energy) defaults to just above the band top so the kinematic
  // (energy-conservation) cutoff frames the spectrum.
  const [params, setParams] = useState(() => ({
    T: temperature ?? 5, Emin: 0, Emax: maxE, sigma: Math.max(0.3, maxE / 100),
    nE: 160, nQbins: 140, Ei: Math.max(5, Math.ceil(maxE * 1.25)),
  }));
  // Display controls (heatmap-only — the intensity curve stays strictly linear).
  const [cmap, setCmap] = useState('viridis');
  const [scale, setScale] = useState('log');          // 'linear' | 'sqrt' | 'log'
  const [contrast, setContrast] = useState(1);         // gamma knob, 0.4..3
  const [showModeDos, setShowModeDos] = useState(false);

  useEffect(() => {
    setParams(p => ({ ...p, Emax: maxE, sigma: Math.max(0.3, maxE / 100), Ei: Math.max(5, Math.ceil(maxE * 1.25)) }));
  }, [maxE]);

  useEffect(() => {
    workerRef.current = new Worker(new URL('../io/sqeworker.js', import.meta.url), { type: 'module' });
    workerRef.current.onmessage = (e) => {
      setRunning(false);
      if (e.data.success) { setOut(e.data); setError(null); }
      else setError(e.data.error || 'INS computation failed');
    };
    return () => workerRef.current?.terminate();
  }, []);

  const run = () => {
    // Validate so an emptied/invalid field (parseFloat → NaN) can't drive a NaN
    // computation (NaN energies → NaN S → downstream render errors).
    const bad = ['T', 'Emin', 'Emax', 'sigma', 'Ei', 'nE', 'nQbins'].find(k => !Number.isFinite(params[k]));
    if (bad) { setError('Please enter a valid number for every field.'); return; }
    if (!(params.Emax > params.Emin)) { setError('E max must be greater than E min.'); return; }
    if (!(params.sigma > 0)) { setError('σ must be greater than 0.'); return; }
    if (params.nE < 2 || params.nQbins < 1) { setError('nE must be ≥ 2 and nQ ≥ 1.'); return; }
    setRunning(true);
    setError(null);
    try {
      const { data, transfer } = buildInsData(results);
      workerRef.current.postMessage({ data, params }, transfer);
    } catch (err) {
      setRunning(false);
      setError(err.message);
    }
  };

  // Integrated intensity g(E) = ∫ S(|Q|,E) d|Q|, summed over the SAME |Q| that
  // the heatmap actually shows (kinematic mask honoured via qRange) so the curve
  // is the literal energy-marginal of the map. Strictly LINEAR real data — the
  // scale/contrast controls never touch it (they'd distort the peak-area ratios).
  const spec = useMemo(() => {
    if (!out?.powResult) return null;
    // Ei from the COMPUTED result (not the live field): S's Q-axis (xMax/dQ) was
    // frozen at compute time to this Ei, so the mask/integration must use it too.
    const { S, nX, nE, xMax, Ei, Eaxis } = out.powResult;
    const dQ = xMax / nX;
    const ki = Ei > 0 ? Math.sqrt(Ei / HBAR2_2MN) : 0;
    const g = new Float64Array(nE);
    for (let ei = 0; ei < nE; ei++) {
      const [qlo, qhi] = qRange(Eaxis[ei], Ei, ki);
      if (qlo > qhi) continue;                         // fully inaccessible row
      let s = 0;
      for (let qi = 0; qi < nX; qi++) {
        const Q = (qi + 0.5) * dQ;
        if (Q < qlo || Q > qhi) continue;
        s += S[qi * nE + ei];
      }
      g[ei] = s * dQ;
    }
    let gmax = 0;
    for (let i = 0; i < nE; i++) if (g[i] > gmax) gmax = g[i];
    return { g, gmax, nE };
  }, [out]);

  // Draw S(Q,E) heatmap, masking the kinematically inaccessible region for a
  // direct-geometry spectrometer with incident energy Eᵢ (see qRange). Intensity
  // → colour goes through the shared scale + contrast transform.
  useEffect(() => {
    if (!out?.powResult || !canvasRef.current) return;
    const { S, nX, nE, Smax, Eaxis, xMax, Ei } = out.powResult;   // Ei from the run, not the live field
    const cv = canvasRef.current;
    cv.width = nX; cv.height = nE;
    const ctx = cv.getContext('2d');
    const img = ctx.createImageData(nX, nE);
    const inv = Smax > 0 ? 1 / Smax : 0;
    const dQ = xMax / nX;
    const ki = Ei > 0 ? Math.sqrt(Ei / HBAR2_2MN) : 0;

    for (let ei = 0; ei < nE; ei++) {
      const [qlo, qhi] = qRange(Eaxis[ei], Ei, ki);
      for (let qi = 0; qi < nX; qi++) {
        const px = ((nE - 1 - ei) * nX + qi) * 4;      // cell-centred: row 0 = Emax
        const Q = (qi + 0.5) * dQ;
        let r, g, b;
        if (Q < qlo || Q > qhi) { [r, g, b] = BG; }
        else {
          const raw = Math.max(0, S[qi * nE + ei] * inv);
          [r, g, b] = colormap(intensityTransform(raw, scale, contrast), cmap);
        }
        img.data[px] = r; img.data[px + 1] = g; img.data[px + 2] = b; img.data[px + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);
  }, [out, cmap, scale, contrast]);

  const exportCsv = () => {
    if (!out?.powResult) return;
    const { S, nX, nE, Eaxis, xMax } = out.powResult;
    const dQ = xMax / nX;
    let csv = 'Q_invA,E_meV,S\n';
    for (let qi = 0; qi < nX; qi++) for (let ei = 0; ei < nE; ei++) {
      csv += `${((qi + 0.5) * dQ).toFixed(5)},${Eaxis[ei].toFixed(5)},${S[qi * nE + ei].toExponential(6)}\n`;
    }
    downloadString(csv, 'sqe.csv');
  };

  // Two-column (E, ∫S dQ) export protects the quantitative curve.
  const exportDos = () => {
    if (!spec || !out?.powResult) return;
    const { Eaxis } = out.powResult;
    let csv = 'E_meV,integrated_intensity\n';
    for (let i = 0; i < spec.nE; i++) csv += `${Eaxis[i].toFixed(5)},${spec.g[i].toExponential(6)}\n`;
    downloadString(csv, 'integrated_intensity.csv');
  };

  const fields = [
    ['T (K)', 'T', 1], ['E min', 'Emin', 1], ['E max', 'Emax', 1], ['σ (meV)', 'sigma', 0.1],
    ['Eᵢ (meV)', 'Ei', 1], ['nE', 'nE', 1], ['nQ', 'nQbins', 1],
  ];
  const scaleLabel = scale === 'log' ? 'log' : scale === 'sqrt' ? '√' : 'linear';

  // Energy range for ALL display mapping comes from the COMPUTED axis, not the
  // live fields, so ticks/heatmap/curve stay consistent even before a re-run.
  const pr = out?.powResult;
  const Emin = pr ? pr.Eaxis[0] : 0;
  const Emax = pr ? pr.Eaxis[pr.nE - 1] : 1;
  const eTicks = pr ? energyTicks(Emin, Emax, pr.nE, 5) : [];

  return (
    <div style={{ ...cardStyle, padding: 20 }}>
      {/* header */}
      <div style={{ display: 'flex', alignItems: 'center', gap: 10, marginBottom: 16, flexWrap: 'wrap' }}>
        <span style={{ width: 24, height: 24, borderRadius: 6, background: ACCENT, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="#fff" strokeWidth="2"><path d="M2 12c2-4 4-4 6 0s4 4 6 0 4-4 6 0" /></svg>
        </span>
        <span style={{ font: "600 15px 'Space Grotesk'", color: INK }}>Simulated INS · S(|Q|,E) &amp; integrated intensity</span>
        <span style={{ font: "11px 'Space Mono'", color: FAINT }}>powder-averaged from eigenvectors</span>
      </div>

      {/* params */}
      <div style={{ display: 'flex', flexWrap: 'wrap', gap: 10, marginBottom: 16 }}>
        {fields.map(([label, key, step]) => (
          <div key={key} style={{ flex: 1, minWidth: 84 }}>
            <div style={{ font: "10px 'Space Mono'", color: FAINT, marginBottom: 5 }}>{label}</div>
            <input type="number" step={step} value={Number.isFinite(params[key]) ? params[key] : ''}
              onChange={e => setParams(p => ({ ...p, [key]: parseFloat(e.target.value) }))} style={insetInput} />
          </div>
        ))}
      </div>

      {/* actions + display controls */}
      <div style={{ display: 'flex', flexWrap: 'wrap', alignItems: 'center', gap: 12, marginBottom: 18 }}>
        <button onClick={run} disabled={running} className="rnr-btn"
          style={{ background: running ? INSET2 : ACCENT, color: running ? DIM : '#fff', border: 'none', borderRadius: 8, padding: '9px 20px', font: "700 13px 'Space Grotesk'", cursor: running ? 'default' : 'pointer' }}>
          {running ? 'Computing…' : 'Run INS'}
        </button>

        <label style={{ display: 'flex', alignItems: 'center', gap: 7, font: "11px 'Space Mono'", color: DIM }}>colormap
          <select value={cmap} onChange={e => setCmap(e.target.value)}
            style={{ background: INSET, border: `1px solid ${BORDER}`, borderRadius: 6, padding: '6px 9px', font: "12px 'Space Mono'", color: INK, cursor: 'pointer' }}>
            {CMAP_NAMES.map(n => <option key={n} value={n}>{n}</option>)}
          </select>
        </label>

        {/* intensity scale (heatmap only) */}
        <div style={{ display: 'flex', alignItems: 'center', gap: 7, font: "11px 'Space Mono'", color: DIM }}>scale
          <div style={{ display: 'flex', gap: 2, background: INSET2, borderRadius: 7, padding: 2 }}>
            {[['linear', 'Lin'], ['sqrt', '√'], ['log', 'Log']].map(([val, lab]) => (
              <button key={val} onClick={() => setScale(val)} className="rnr-btn"
                style={{ border: 'none', borderRadius: 5, padding: '5px 10px', cursor: 'pointer', font: "600 11px 'Space Mono'",
                  background: scale === val ? 'var(--card)' : 'transparent', color: scale === val ? 'var(--accentInk)' : DIM,
                  boxShadow: scale === val ? '0 1px 2px rgba(16,24,38,0.12)' : 'none' }}>{lab}</button>
            ))}
          </div>
        </div>

        {/* contrast (gamma) */}
        <label style={{ display: 'flex', alignItems: 'center', gap: 8, font: "11px 'Space Mono'", color: DIM }} title="Higher lifts faint features (heatmap only)">
          contrast
          <input type="range" min={0.4} max={3} step={0.1} value={contrast} onChange={e => setContrast(parseFloat(e.target.value))}
            style={{ width: 92, accentColor: ACCENT }} />
          <span style={{ color: 'var(--accentInk)', minWidth: 26 }}>×{contrast.toFixed(1)}</span>
        </label>

        <label style={{ display: 'flex', alignItems: 'center', gap: 6, font: "11px 'Space Mono'", color: DIM, cursor: 'pointer' }} title="Overlay the true unweighted phonon DOS (shape only)">
          <input type="checkbox" checked={showModeDos} onChange={e => setShowModeDos(e.target.checked)} /> mode DOS
        </label>

        <span style={{ font: "10.5px 'Space Mono'", color: FAINT }}>Eᵢ = 0 ⇒ direct (full Q range)</span>

        {pr && (
          <div style={{ marginLeft: 'auto', display: 'flex', gap: 8 }}>
            <button onClick={exportDos} className="rnr-btn"
              style={{ background: INSET, border: `1px solid ${BORDER}`, borderRadius: 7, padding: '8px 12px', font: "600 12px 'Space Grotesk'", color: INK, cursor: 'pointer' }}>
              Export g(E) CSV
            </button>
            <button onClick={exportCsv} className="rnr-btn"
              style={{ background: INSET, border: `1px solid ${BORDER}`, borderRadius: 7, padding: '8px 12px', font: "600 12px 'Space Grotesk'", color: INK, cursor: 'pointer' }}>
              Export S(Q,E) CSV
            </button>
          </div>
        )}
      </div>

      {error && <div style={{ color: 'var(--warnInk)', font: "13px 'Space Mono'", marginBottom: 12 }}>{error}</div>}

      {pr ? (
        <>
          <div style={{ display: 'flex', justifyContent: 'space-between', font: "11px 'Space Mono'", color: DIM, marginBottom: 10, flexWrap: 'wrap', gap: 8 }}>
            <span>S(|Q|,E) — powder-averaged simulated INS <span style={{ color: FAINT }}>· colour = intensity ({scaleLabel})</span></span>
            <span style={{ color: FAINT }}>∫S(|Q|,E) d|Q| → <span style={{ color: DIM }}>neutron-weighted, Bose-populated — not a DOS</span></span>
          </div>

          {/* One aligned figure: a single grid row makes the heatmap and the
              integrated-intensity curve share an identical energy axis + height. */}
          <div className="sqe-fig">
            <div className="sqe-grid">
              {/* shared energy tick column (the ONLY E-axis renderer) */}
              <div style={{ gridArea: 'ticks', position: 'relative', minWidth: 30, height: '100%' }}>
                {eTicks.map((t, i) => (
                  <span key={i} style={{ position: 'absolute', right: 4, top: `${t.frac * 100}%`, transform: 'translateY(-50%)', font: "10px 'Space Mono'", color: FAINT, whiteSpace: 'nowrap' }}>{t.label}</span>
                ))}
                <span style={{ position: 'absolute', left: -2, top: '50%', transform: 'translateY(-50%) rotate(-90deg)', transformOrigin: 'left center', font: "10px 'Space Mono'", color: DIM }}>E (meV)</span>
              </div>

              {/* heatmap */}
              <canvas ref={canvasRef} style={{ gridArea: 'heat', width: '100%', height: '100%', minWidth: 0, minHeight: 0, display: 'block', border: `1px solid ${BORDER}`, borderRadius: 8, imageRendering: 'pixelated', background: '#0b0e16' }} />

              {/* colorbar (same gradient stops as the pixel LUT) */}
              <div className="sqe-bar" style={{ gridArea: 'bar', display: 'flex', flexDirection: 'column', alignItems: 'center', minHeight: 0 }}>
                <div style={{ flex: 1, width: 13, borderRadius: 4, border: `1px solid ${BORDER}`, background: cmapGradient(cmap) }} />
                <span style={{ font: "10px 'Space Mono'", color: FAINT, marginTop: 5 }}>S↑</span>
              </div>

              {/* integrated-intensity curve — same energy axis as the heatmap */}
              <div className="sqe-curve" style={{ gridArea: 'curve', minWidth: 0, minHeight: 0 }}>
                <DosPlot spec={spec} Emin={Emin} Emax={Emax}
                  modeDos={showModeDos ? out.dosResult : null} />
              </div>

              {/* row 2 — axis labels under each panel */}
              <div style={{ gridArea: 'qlab' }}>
                <div style={{ display: 'flex', justifyContent: 'space-between', marginTop: 5, font: "10px 'Space Mono'", color: FAINT }}>
                  {ticks(0, pr.xMax, 5).map((t, i) => <span key={i}>{t.v.toFixed(1)}</span>)}
                </div>
                <div style={{ textAlign: 'center', font: "10px 'Space Mono'", color: DIM, marginTop: 2 }}>|Q| (Å⁻¹)</div>
              </div>
              <div style={{ gridArea: 'glab', textAlign: 'center', font: "10px 'Space Mono'", color: DIM, marginTop: 5, alignSelf: 'start' }}>
                intensity →
              </div>
            </div>
          </div>
        </>
      ) : (
        <div style={{ font: "12px 'Spline Sans'", color: FAINT }}>Press <b style={{ color: DIM }}>Run INS</b> to compute the powder S(|Q|,E) map and the integrated-intensity curve.</div>
      )}
    </div>
  );
}

function ticks(min, max, n) {
  const arr = [];
  for (let i = 0; i <= n; i++) { const frac = i / n; arr.push({ frac, v: min + frac * (max - min) }); }
  return arr;
}

// Energy ticks positioned by the SAME cell-centred fraction the heatmap uses:
// frac 0 = top edge (Emax + dE/2), frac 1 = bottom edge (Emin − dE/2).
function energyTicks(Emin, Emax, nE, n) {
  const dE = (Emax - Emin) / (nE - 1);
  const Etop = Emax + dE / 2, span = nE * dE || 1;
  const out = [];
  for (let i = 0; i <= n; i++) {
    const v = Emin + (i / n) * (Emax - Emin);
    out.push({ label: v.toFixed(0), frac: (Etop - v) / span });
  }
  return out;
}

/**
 * Integrated-intensity curve: g(E) = ∫S d|Q| with energy on the VERTICAL axis,
 * mapped by the exact same cell-centred fraction as the heatmap so a given energy
 * sits at the same pixel row in both. Rendered at devicePixelRatio for crispness.
 * `modeDos` (optional) overlays the true unweighted phonon DOS (shape reference).
 */
function DosPlot({ spec, Emin, Emax, modeDos }) {
  const ref = useRef(null);
  useEffect(() => {
    const cv = ref.current;
    // Degenerate/empty energy window: blank the canvas so a stale curve from a
    // prior valid run doesn't linger next to a freshly redrawn heatmap.
    if (!cv || !spec || !(isFinite(Emin) && isFinite(Emax) && Emax > Emin)) {
      if (cv) cv.getContext('2d')?.clearRect(0, 0, cv.width, cv.height);
      return;
    }
    const { g, gmax, nE } = spec;
    const draw = () => {
      const dpr = window.devicePixelRatio || 1;
      const rect = cv.getBoundingClientRect();          // includes .rnr zoom (dpr does not)
      const W = Math.max(1, Math.round(rect.width * dpr));
      const H = Math.max(1, Math.round(rect.height * dpr));
      if (cv.width !== W) cv.width = W;
      if (cv.height !== H) cv.height = H;
      const ctx = cv.getContext('2d');
      ctx.clearRect(0, 0, W, H);
      const dE = (Emax - Emin) / (nE - 1);
      const Etop = Emax + dE / 2, span = nE * dE || 1;
      const yOf = (E) => ((Etop - E) / span) * H;       // matches heatmap row nE-1-ei
      const gm = gmax || 1;

      // optional true-DOS overlay (faint dashed, own normalization, shape only)
      if (modeDos && modeDos.dos) {
        const mm = modeDos.dosMax || 1;
        ctx.strokeStyle = 'rgba(107,116,136,0.75)'; ctx.setLineDash([4 * dpr, 3 * dpr]); ctx.lineWidth = 1 * dpr;
        ctx.beginPath();
        for (let i = 0; i < nE; i++) { const x = (modeDos.dos[i] / mm) * W, y = yOf(Emin + i * dE); if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y); }
        ctx.stroke(); ctx.setLineDash([]);
      }

      // integrated intensity — filled area from the left baseline + solid line
      ctx.beginPath();
      ctx.moveTo(0, yOf(Emin));
      for (let i = 0; i < nE; i++) ctx.lineTo((g[i] / gm) * W, yOf(Emin + i * dE));
      ctx.lineTo(0, yOf(Emax));
      ctx.closePath();
      ctx.fillStyle = DOS_FILL; ctx.fill();
      ctx.strokeStyle = DOS_STROKE; ctx.lineWidth = 1.6 * dpr; ctx.beginPath();
      for (let i = 0; i < nE; i++) { const x = (g[i] / gm) * W, y = yOf(Emin + i * dE); if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y); }
      ctx.stroke();
    };
    draw();
    const ro = new ResizeObserver(draw);
    ro.observe(cv);
    return () => ro.disconnect();
  }, [spec, Emin, Emax, modeDos]);
  return <canvas ref={ref} style={{ width: '100%', height: '100%', display: 'block', border: '1px solid var(--border)', borderRadius: 8, background: 'var(--inset)' }} />;
}
