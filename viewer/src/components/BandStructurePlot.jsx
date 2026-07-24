import React, { useMemo, useRef, useState, useEffect } from 'react';
import { conventionalLattice, reciprocalLattice, pathDistances } from '../math/reciprocal';
import { clientToSvg } from './svgCoords';

/**
 * Physical band-structure plot (custom SVG, Cobalt light theme).
 *
 * x-axis: cumulative reciprocal-path distance (Å⁻¹) with high-symmetry tick
 * labels (Noto Sans) and vertical segment dividers. y-axis: energy (`unit`)
 * with adaptive gridlines and a dashed zero line. Branches are drawn in a red
 * family; sub-zero portions are overdrawn dashed (soft-mode / imaginary freq).
 * A blue marker flags the selected (k, mode).
 *
 * Interactions: drag a box *inside the plot frame* to zoom — the boxed data
 * region is remapped to fill the whole frame and the axis tick values recompute
 * for the zoomed range (a true data-domain zoom, not a pixel crop). A single
 * click (no drag) selects the nearest (k, mode). `resetSignal` (a changing
 * nonce) clears the zoom back to the full range.
 *
 * Axis math (reciprocalLattice / pathDistances / kpathMeta.hsymIndex) is the
 * same as before — only the styling/zoom changed.
 */
const W = 760, H = 880;
const XL = 64, XR = 740, YT = 30, YB = 810;   // plot frame (pixel/user space)
const PW = XR - XL, PH = YB - YT;
const CLIP = 'bsp-frame-clip';
// Red family — distinct hue/lightness per branch while reading as "red".
const branchColor = (m) => `hsl(${(m * 6) % 16}, 72%, ${40 + (m % 7) * 3.2}%)`;

// A "nice" gridline step (1/2/5 × 10ⁿ) for ~`target` divisions over `range`.
function niceStep(range, target = 6) {
  const raw = Math.max(range, 1e-9) / target;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const n = raw / mag;
  return (n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10) * mag;
}

export default function BandStructurePlot({ bands, qPoints, baseStructure, kpathMeta, selected, onPick, eMin, eMax, unit = 'meV', resetSignal = 0 }) {
  const svgRef = useRef(null);
  const [domain, setDomain] = useState(null);    // zoom window {xMin,xMax,yMin,yMax} in DATA coords, or null = full
  const [sel, setSel] = useState(null);          // live drag box {x0,y0,x1,y1} in pixel space
  const [hover, setHover] = useState(null);      // nearest point under the cursor {k,m,x,y}
  const dragRef = useRef(null);

  // Reset zoom whenever the parent bumps resetSignal.
  useEffect(() => { setDomain(null); }, [resetSignal]);

  const model = useMemo(() => {
    if (!bands?.length || !qPoints?.length) return null;
    // q-points are in conventional coords; prefer the conventional reciprocal the
    // runner attaches (bandRecip) so path distances are right even when the 3D
    // viewer cell is a custom/primitive computation cell.
    const recip = baseStructure.bandRecip || reciprocalLattice(conventionalLattice(baseStructure.v1, baseStructure.v2, baseStructure.v3, baseStructure.dim));
    const dist = pathDistances(qPoints, recip, kpathMeta?.segSizes);
    const xMax = dist[dist.length - 1] || 1;
    let yMin = 0, yMax = 1;
    for (const row of bands) for (const v of row) { if (!isFinite(v)) continue; if (v < yMin) yMin = v; if (v > yMax) yMax = v; }
    yMax *= 1.05; if (yMin < 0) yMin *= 1.05;
    if (Number.isFinite(eMin)) yMin = eMin;
    if (Number.isFinite(eMax)) yMax = eMax;
    const nModes = bands[0].length;

    // Effective (possibly zoomed) data window mapped onto the full plot frame.
    const xLo = domain ? domain.xMin : 0, xHi = domain ? domain.xMax : xMax;
    const yLo = domain ? domain.yMin : yMin, yHi = domain ? domain.yMax : yMax;
    const xSpan = (xHi - xLo) || 1, ySpan = (yHi - yLo) || 1;

    const xOf = (d) => XL + ((d - xLo) / xSpan) * PW;
    const yOf = (e) => YB - ((e - yLo) / ySpan) * PH;
    // Inverse (pixel → data) for converting the drag box into a data window.
    const invX = (px) => xLo + ((px - XL) / PW) * xSpan;
    const invY = (py) => yLo + ((YB - py) / PH) * ySpan;

    // One polyline per branch, plus a separate path of only its sub-zero spans.
    const branches = [];
    for (let m = 0; m < nModes; m++) {
      let d = '', neg = '', penNeg = false;
      for (let k = 0; k < bands.length; k++) {
        const v = bands[k][m];
        if (!isFinite(v)) { penNeg = false; continue; }
        const X = xOf(dist[k]).toFixed(1), Y = yOf(v).toFixed(1);
        d += (d ? 'L' : 'M') + X + ' ' + Y + ' ';
        if (v < 0) { neg += (penNeg ? 'L' : 'M') + X + ' ' + Y + ' '; penNeg = true; }
        else penNeg = false;
      }
      branches.push({ d, neg, color: branchColor(m) });
    }

    // High-symmetry ticks — only those whose distance falls in the x-window.
    // Ticks that land at the same x (a break/discontinuity, where two labels
    // would overlap) are merged into one, e.g. "X,L".
    const hsym = kpathMeta?.hsymIndex || {};
    const eps = xSpan * 1e-6;
    const xmap = new Map();
    for (const k of Object.keys(hsym).map(Number).sort((a, b) => a - b)) {
      const dx = dist[k];
      if (!(dx >= xLo - eps && dx <= xHi + eps)) continue;
      const x = xOf(dx), key = Math.round(x);
      const e = xmap.get(key) || { x, labels: [] };
      if (!e.labels.includes(hsym[k])) e.labels.push(hsym[k]);
      xmap.set(key, e);
    }
    const xticks = [...xmap.values()].map(e => ({ x: e.x, label: e.labels.join(',') }));

    // Adaptive y gridlines across the (zoomed) y-window, with matching decimals.
    const step = niceStep(yHi - yLo);
    const dec = step >= 1 ? 0 : step >= 0.1 ? 1 : 2;
    const ygrid = [];
    for (let e = Math.ceil(yLo / step) * step; e <= yHi + step * 1e-6; e += step) {
      ygrid.push({ y: yOf(e), label: e.toFixed(dec) });
    }

    const zeroIn = 0 >= yLo && 0 <= yHi;
    return { dist, xOf, yOf, invX, invY, branches, xticks, ygrid, nModes, zeroY: yOf(0), zeroIn };
  }, [bands, qPoints, baseStructure, kpathMeta, eMin, eMax, domain]);

  if (!model) return null;

  // Map client coords → viewBox user space (== pixel space here, viewBox is the
  // full frame). Uses getBoundingClientRect, not getScreenCTM, so it stays
  // aligned with the cursor under the app's ancestor `zoom` (see svgCoords.js).
  const toUser = (evt) => clientToSvg(svgRef.current, evt.clientX, evt.clientY);
  const clampX = (x) => Math.max(XL, Math.min(XR, x));
  const clampY = (y) => Math.max(YT, Math.min(YB, y));

  // Nearest selectable (k, mode) to a point, with its pixel position.
  const pickNearest = (loc) => {
    let best = null, bestD = Infinity;
    for (let k = 0; k < bands.length; k++) {
      const x = model.xOf(model.dist[k]);
      for (let m = 0; m < model.nModes; m++) {
        const v = bands[k][m];
        if (!isFinite(v)) continue;
        const y = model.yOf(v);
        const dd = (x - loc.x) ** 2 + (y - loc.y) ** 2;
        if (dd < bestD) { bestD = dd; best = { k, m, x, y }; }
      }
    }
    return best;
  };

  const onDown = (e) => { const p = toUser(e); if (p) dragRef.current = { x0: clampX(p.x), y0: clampY(p.y), moved: false }; };
  const onMove = (e) => {
    const p = toUser(e); if (!p) return;
    if (dragRef.current) {
      const d = dragRef.current;
      const x1 = clampX(p.x), y1 = clampY(p.y);
      if (Math.abs(x1 - d.x0) > 4 || Math.abs(y1 - d.y0) > 4) d.moved = true;
      if (d.moved) { setSel({ x0: d.x0, y0: d.y0, x1, y1 }); setHover(null); }
      return;
    }
    // Not dragging: preview the nearest selectable point under the cursor.
    setHover(pickNearest(p));
  };
  const onUp = (e) => {
    const d = dragRef.current; dragRef.current = null; setSel(null);
    const p = toUser(e); if (!p || !d) return;
    if (d.moved) {
      // Box-zoom: pixel box (clamped to frame) → data window. Require a usable size.
      const xa = clampX(d.x0), xb = clampX(p.x), ya = clampY(d.y0), yb = clampY(p.y);
      if (Math.abs(xb - xa) > 12 && Math.abs(yb - ya) > 12) {
        setDomain({
          xMin: model.invX(Math.min(xa, xb)), xMax: model.invX(Math.max(xa, xb)),
          yMin: model.invY(Math.max(ya, yb)), yMax: model.invY(Math.min(ya, yb)),   // pixel-y is inverted
        });
      }
    } else {
      const hit = pickNearest(p);
      if (hit && onPick) onPick(hit.k, hit.m);
    }
  };

  const selPt = selected && isFinite(bands[selected.k]?.[selected.m])
    ? { x: model.xOf(model.dist[selected.k]), y: model.yOf(bands[selected.k][selected.m]) }
    : null;
  const inFrame = (p) => p && p.x >= XL - 1 && p.x <= XR + 1 && p.y >= YT - 1 && p.y <= YB + 1;

  return (
    <svg ref={svgRef} viewBox={`0 0 ${W} ${H}`} style={{ display: 'block', width: '100%', height: 'auto', cursor: 'crosshair', userSelect: 'none' }}
      onMouseDown={onDown} onMouseMove={onMove} onMouseUp={onUp} onMouseLeave={() => { dragRef.current = null; setSel(null); setHover(null); }}>

      <defs><clipPath id={CLIP}><rect x={XL} y={YT} width={PW} height={PH} /></clipPath></defs>

      {/* y grid + labels */}
      {model.ygrid.map((g, i) => (
        <g key={`y${i}`}>
          <line x1={XL} x2={XR} y1={g.y} y2={g.y} stroke="var(--border)" strokeWidth="1" />
          <text x={56} y={g.y + 4} textAnchor="end" fontFamily="STIX Two Text, serif" fontSize="11" fill="var(--faint)">{g.label}</text>
        </g>
      ))}

      {/* zero line (only when 0 is in range) */}
      {model.zeroIn && <line x1={XL} x2={XR} y1={model.zeroY} y2={model.zeroY} stroke="var(--dim)" strokeWidth="0.8" strokeDasharray="3 3" opacity="0.6" />}

      {/* high-symmetry dividers + labels */}
      {model.xticks.map((t, i) => (
        <g key={`x${i}`}>
          <line x1={t.x} x2={t.x} y1={YT} y2={YB} stroke="var(--border)" strokeWidth="1" />
          <text x={t.x} y={832} textAnchor="middle" fontFamily="Noto Sans, sans-serif" fontSize="16" fontWeight="700" fill="var(--ink)">{t.label}</text>
        </g>
      ))}

      {/* axis frame */}
      <line x1={XL} x2={XL} y1={YT} y2={YB} stroke="var(--dim)" strokeWidth="1.1" />
      <line x1={XL} x2={XR} y1={YB} y2={YB} stroke="var(--dim)" strokeWidth="1.1" />

      {/* bands (clipped to the frame so zoomed-out-of-range parts don't overflow) */}
      <g clipPath={`url(#${CLIP})`}>
        {model.branches.map((b, m) => (
          <g key={m}>
            <path d={b.d} fill="none" stroke={b.color} strokeWidth="1.6" strokeLinejoin="round" opacity="0.9" />
            {b.neg && <path d={b.neg} fill="none" stroke="#e0564b" strokeWidth="2" strokeLinejoin="round" strokeDasharray="3 2" />}
          </g>
        ))}
        {/* hover marker — translucent spot on the nearest selectable point */}
        {hover && !sel && <circle cx={hover.x} cy={hover.y} r="7" fill="rgba(47,109,240,0.22)" stroke="rgba(47,109,240,0.5)" strokeWidth="1" style={{ pointerEvents: 'none' }} />}

        {/* selection marker (blue, for contrast against the red branches) */}
        {inFrame(selPt) && <circle cx={selPt.x} cy={selPt.y} r="6" fill="#2f6df0" stroke="#fff" strokeWidth="2" />}
      </g>

      {/* drag-to-zoom selection box */}
      {sel && (
        <rect x={Math.min(sel.x0, sel.x1)} y={Math.min(sel.y0, sel.y1)} width={Math.abs(sel.x1 - sel.x0)} height={Math.abs(sel.y1 - sel.y0)}
          fill="rgba(47,109,240,0.10)" stroke="#2f6df0" strokeWidth="1" strokeDasharray="4 3" />
      )}

      {/* axis title */}
      <text x={20} y={(YT + YB) / 2} fontFamily="STIX Two Text, serif" fontSize="14" fill="var(--dim)" textAnchor="middle" transform={`rotate(-90 20 ${(YT + YB) / 2})`}>Energy ({unit})</text>

      {/* zoom hint */}
      <text x={XR} y={YT + 4} textAnchor="end" fontFamily="Space Mono, monospace" fontSize="11" fill="var(--faint)">{domain ? 'zoomed · Reset to clear' : 'drag a box to zoom'}</text>
    </svg>
  );
}
