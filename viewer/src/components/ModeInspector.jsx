import React, { useMemo } from 'react';
import { DEFAULT_COLORS } from '../constants';

const DIM = 'var(--dim)', INK = 'var(--ink)', FAINT = 'var(--faint)', INSET2 = 'var(--inset2)';
const ACCENTINK = 'var(--accentInk)', WARNINK = 'var(--warnInk)';

/**
 * Selected-mode readout (Cobalt light footer): mode/energy/k on the left and
 * per-element band character on the right. Weights are the eigenvector content
 * per species (Σ_sites∈el √Σ_xyz|e|², normalised) — which species dominate the
 * mode. `unit`/`scale` convert the displayed energy; `colors` overrides the bar
 * colours to match the viewer's editable element colours.
 */
export default function ModeInspector({ results, selectedK, selectedMode, unit = 'meV', scale = 1, colors = {} }) {
  const info = useMemo(() => {
    if (!results) return null;
    const { bands, qPoints, eigvecs, baseStructure } = results;
    const e = bands[selectedK]?.[selectedMode];
    const q = qPoints[selectedK];
    const ev = eigvecs?.[selectedK]?.[selectedMode];
    if (e === undefined || !ev) return null;

    const { atomDic, uniqueRN } = baseStructure;
    const rnToRow = new Map((uniqueRN || []).map((rn, r) => [rn, r]));
    const weights = [];
    let total = 0;
    for (const el of Object.keys(atomDic)) {
      let w = 0;
      for (const rn of atomDic[el]) {
        const r = rnToRow.get(rn);
        if (r === undefined) continue;
        let s = 0;
        for (let c = 0; c < 3; c++) {
          const re = ev.real[r * 3 + c], im = ev.imag[r * 3 + c];
          s += re * re + im * im;
        }
        w += Math.sqrt(s);
      }
      weights.push({ el, w });
      total += w;
    }
    weights.forEach(x => { x.frac = total > 0 ? x.w / total : 0; });
    weights.sort((a, b) => b.frac - a.frac);
    return { e, q, weights };
  }, [results, selectedK, selectedMode]);

  if (!info) return null;
  const colorOf = (el) => colors[el] || DEFAULT_COLORS[el] || '#9aa1b2';
  const row = { display: 'flex', justifyContent: 'space-between', gap: 14, font: "11px 'Space Mono'", color: DIM };

  return (
    <div style={{ borderTop: '1px solid var(--border)', padding: '12px 16px', display: 'flex', alignItems: 'center', gap: 22, flexWrap: 'wrap' }}>
      <div style={{ display: 'flex', flexDirection: 'column', gap: 4, minWidth: 148 }}>
        <div style={row}><span>mode</span><span style={{ color: INK }}>band {selectedMode + 1} · k{selectedK + 1}</span></div>
        <div style={row}><span>energy</span><span style={{ color: info.e < 0 ? WARNINK : ACCENTINK, fontWeight: 700 }}>{(info.e * scale).toFixed(3)} {unit}</span></div>
        <div style={row}><span>k</span><span style={{ color: INK, fontSize: 10 }}>[{info.q.map(v => v.toFixed(3)).join(', ')}]</span></div>
      </div>
      <div style={{ flex: 1, minWidth: 180, display: 'flex', flexDirection: 'column', gap: 5 }}>
        <div style={{ font: "10px 'Space Mono'", letterSpacing: '.12em', color: FAINT }}>BAND CHARACTER</div>
        {info.weights.map(({ el, frac }) => (
          <div key={el} style={{ display: 'flex', alignItems: 'center', gap: 8, font: "11px 'Space Mono'" }}>
            <span style={{ width: 18, color: DIM }}>{el}</span>
            <div style={{ flex: 1, height: 7, borderRadius: 4, background: INSET2, overflow: 'hidden' }}>
              <div style={{ height: '100%', borderRadius: 4, width: `${frac * 100}%`, background: colorOf(el) }} />
            </div>
            <span style={{ width: 30, textAlign: 'right', color: DIM }}>{(frac * 100).toFixed(0)}%</span>
          </div>
        ))}
      </div>
    </div>
  );
}
