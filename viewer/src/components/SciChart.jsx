import React, { useMemo, useRef, useState, useEffect, useId } from 'react';
import { clientToSvg } from './svgCoords';

/**
 * Journal-style scientific line chart (Cobalt redesign): real axes with tick
 * marks, numeric labels and unit titles set in Open Sans, light gridlines, an
 * optional dashed zero line, and an in-plot data/fit legend.
 *
 * series = [{ points:[[x,y],…], color, width?, dashed?, marker?, name?, vy? }].
 * A series with `marker` is drawn as open circles (observed data); others as
 * solid/dashed lines. `vy` (parallel to points) supplies the value reported on
 * hover when it differs from the plotted y (e.g. an offset difference curve).
 * `baselines` draws dashed grey horizontal reference lines. Auto-ranges axes.
 *
 * Interactive:
 *   • hover snaps a crosshair to the nearest x and shows a value tooltip,
 *   • drag a box to zoom into a region; double-click (or the reset button)
 *     restores the full view.
 *
 * Responsive: the viewBox tracks the rendered pixel width (1 unit = 1px), so
 * the plot fills its container at any size with no letterboxing.
 */
const M = { l: 48, r: 14, t: 16, b: 38 };
const clamp = (v, a, b) => Math.max(a, Math.min(b, v));

// "Nice" round tick values (1·2·5 × 10ⁿ) inside [min,max] — avoids the random
// decimals you get from evenly slicing the raw data range.
function axisTicks(min, max, n = 5) {
  if (!(max > min)) return { ticks: [min], decimals: 0 };
  const raw = (max - min) / n;
  const exp = Math.floor(Math.log10(raw));
  const base = Math.pow(10, exp);
  const f = raw / base;
  const step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * base;
  const ticks = [];
  for (let v = Math.ceil(min / step - 1e-9) * step; v <= max + step * 1e-6; v += step) {
    ticks.push(Math.abs(v) < step * 1e-6 ? 0 : v);
  }
  return { ticks, decimals: Math.max(0, -Math.floor(Math.log10(step) + 1e-9)) };
}
const fmtTick = (v, d) => (Object.is(v, -0) || Math.abs(v) < 1e-12 ? 0 : v).toFixed(d);

export default function SciChart({ series = [], xLabel, yLabel, zeroLine = false, baselines = [], height = 220, resetKey }) {
  const wrapRef = useRef(null);
  const svgRef = useRef(null);
  const clip = useId().replace(/:/g, '');
  const [w, setW] = useState(620);
  const [hover, setHover] = useState(null);
  const [zoom, setZoom] = useState(null);   // {xmin,xmax,ymin,ymax} | null (auto)
  const [drag, setDrag] = useState(null);    // {x0,y0,x1,y1} pixels while selecting
  const H = height;

  useEffect(() => {
    const el = wrapRef.current; if (!el || typeof ResizeObserver === 'undefined') return;
    const ro = new ResizeObserver(([e]) => setW(Math.max(240, Math.round(e.contentRect.width))));
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  // Reset the zoom when the caller signals new data (e.g. a different config).
  // Keyed on `resetKey` (not the series array, which is a fresh reference each
  // render) so an unrelated parent re-render doesn't wipe the user's zoom.
  useEffect(() => { setZoom(null); }, [resetKey]);

  const model = useMemo(() => {
    let axmin = Infinity, axmax = -Infinity, aymin = Infinity, aymax = -Infinity;
    for (const s of series) for (const [x, y] of s.points) {
      if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
      if (x < axmin) axmin = x; if (x > axmax) axmax = x;
      if (y < aymin) aymin = y; if (y > aymax) aymax = y;
    }
    if (!Number.isFinite(axmin)) return null;
    if (axmax === axmin) axmax = axmin + 1;
    let xmin, xmax, ymin, ymax;
    if (zoom) { xmin = zoom.xmin; xmax = zoom.xmax; ymin = zoom.ymin; ymax = zoom.ymax; }
    else { xmin = axmin; xmax = axmax; const pad = (aymax - aymin) * 0.06 || 1; ymin = aymin - pad; ymax = aymax + pad; }
    const pw = w - M.l - M.r, ph = H - M.t - M.b;
    const xOf = (x) => M.l + ((x - xmin) / (xmax - xmin)) * pw;
    const yOf = (y) => M.t + ph - ((y - ymin) / (ymax - ymin)) * ph;
    const xInv = (px) => xmin + ((px - M.l) / pw) * (xmax - xmin);
    const yInv = (py) => ymin + ((M.t + ph - py) / ph) * (ymax - ymin);
    const paths = series.map(s => {
      let d = '';
      for (const [x, y] of s.points) {
        if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
        d += (d ? 'L' : 'M') + xOf(x).toFixed(1) + ' ' + yOf(y).toFixed(1) + ' ';
      }
      return d.trim();
    });
    const xa = axisTicks(xmin, xmax, 6), ya = axisTicks(ymin, ymax, 5);
    const xt = xa.ticks.map(v => ({ p: xOf(v), v }));
    const yt = ya.ticks.map(v => ({ p: yOf(v), v }));
    return { xOf, yOf, xInv, yInv, paths, xt, yt, xdec: xa.decimals, ydec: ya.decimals,
      x0: M.l, x1: w - M.r, y0: H - M.b, y1: M.t, zeroY: yOf(0), hasZero: ymin < 0 && ymax > 0 };
  }, [series, w, H, zoom]);

  if (!model) return <div style={{ height: H, display: 'flex', alignItems: 'center', justifyContent: 'center', color: 'var(--faint)', font: "12px 'Spline Sans'" }}>no data</div>;

  const fmt = (v) => Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 1 ? v.toFixed(2) : Math.abs(v) < 1e-9 ? '0' : v.toFixed(3);
  const font = 'Open Sans, sans-serif';
  const xsym = (xLabel || 'x').split('(')[0].trim() || 'x';
  const inPlot = (p) => p.x >= model.x0 - 1 && p.x <= model.x1 + 1 && p.y >= model.y1 - 1 && p.y <= model.y0 + 1;

  // getBoundingClientRect-based (not getScreenCTM) so hover/zoom-box picking
  // tracks the cursor under the app's ancestor `zoom` (see svgCoords.js).
  const toUser = (e) => clientToSvg(svgRef.current, e.clientX, e.clientY);

  const onDown = (e) => {
    const u = toUser(e); if (!u || !inPlot(u)) return;
    setHover(null); setDrag({ x0: u.x, y0: u.y, x1: u.x, y1: u.y });
  };
  const onMove = (e) => {
    const u = toUser(e); if (!u) return;
    if (drag) { setDrag(d => ({ ...d, x1: clamp(u.x, model.x0, model.x1), y1: clamp(u.y, model.y1, model.y0) })); return; }
    if (u.x < model.x0 - 2 || u.x > model.x1 + 2) { if (hover) setHover(null); return; }
    const items = [];
    for (const s of series) {
      let best = -1, bd = Infinity;
      for (let k = 0; k < s.points.length; k++) {
        const xx = s.points[k][0]; if (!Number.isFinite(xx)) continue;
        const d = Math.abs(model.xOf(xx) - u.x); if (d < bd) { bd = d; best = k; }
      }
      if (best < 0) continue;
      const [x, y] = s.points[best];
      items.push({ name: s.name, color: s.color, cx: model.xOf(x), cy: model.yOf(y), x, vy: s.vy ? s.vy[best] : y });
    }
    if (!items.length) { if (hover) setHover(null); return; }
    setHover({ mx: items[0].cx, x: items[0].x, items });
  };
  const onUp = () => {
    if (!drag) return;
    const { x0, y0, x1, y1 } = drag;
    if (Math.abs(x1 - x0) > 6 && Math.abs(y1 - y0) > 6) {
      setZoom({
        xmin: model.xInv(Math.min(x0, x1)), xmax: model.xInv(Math.max(x0, x1)),
        ymin: model.yInv(Math.max(y0, y1)), ymax: model.yInv(Math.min(y0, y1)),
      });
    }
    setDrag(null);
  };
  const onLeave = () => { setHover(null); setDrag(null); };

  return (
    <div ref={wrapRef} style={{ position: 'relative', width: '100%' }}>
      {zoom && (
        <button onClick={() => setZoom(null)} className="rnr-btn"
          style={{ position: 'absolute', top: 6, right: 8, zIndex: 2, font: "600 10px 'Space Mono'", color: 'var(--accentInk)', background: 'var(--soft)', border: '1px solid var(--border)', borderRadius: 6, padding: '3px 8px', cursor: 'pointer' }}>
          ⤢ reset zoom
        </button>
      )}
      <svg ref={svgRef} width="100%" height={H} viewBox={`0 0 ${w} ${H}`}
        style={{ display: 'block', cursor: drag ? 'crosshair' : 'crosshair', userSelect: 'none', touchAction: 'none' }}
        onMouseDown={onDown} onMouseMove={onMove} onMouseUp={onUp} onMouseLeave={onLeave} onDoubleClick={() => setZoom(null)}>
        <defs><clipPath id={clip}><rect x={model.x0} y={model.y1} width={model.x1 - model.x0} height={model.y0 - model.y1} /></clipPath></defs>

        {/* gridlines */}
        <g stroke="var(--border)" strokeWidth="0.7">
          {model.yt.map((t, i) => <line key={`yg${i}`} x1={model.x0} y1={t.p} x2={model.x1} y2={t.p} />)}
          {model.xt.map((t, i) => <line key={`xg${i}`} x1={t.p} y1={model.y1} x2={t.p} y2={model.y0} />)}
        </g>
        {zeroLine && model.hasZero && <line x1={model.x0} y1={model.zeroY} x2={model.x1} y2={model.zeroY} stroke="var(--dim)" strokeWidth="0.8" strokeDasharray="2 2" />}
        <g clipPath={`url(#${clip})`}>
          {baselines.map((b, i) => <line key={`bl${i}`} x1={model.x0} y1={model.yOf(b)} x2={model.x1} y2={model.yOf(b)} stroke="var(--bar)" strokeWidth="0.9" strokeDasharray="3 2.5" />)}
        </g>

        {/* axes + ticks */}
        <g stroke="var(--dim)" strokeWidth="1.1">
          <line x1={model.x0} y1={model.y1} x2={model.x0} y2={model.y0} />
          <line x1={model.x0} y1={model.y0} x2={model.x1} y2={model.y0} />
          {model.yt.map((t, i) => <line key={`yt${i}`} x1={model.x0 - 4} y1={t.p} x2={model.x0} y2={t.p} />)}
          {model.xt.map((t, i) => <line key={`xt${i}`} x1={t.p} y1={model.y0} x2={t.p} y2={model.y0 + 4} />)}
        </g>
        <g fill="var(--faint)" fontFamily={font} fontSize="10">
          {model.yt.map((t, i) => <text key={`yl${i}`} x={model.x0 - 7} y={t.p + 3} textAnchor="end">{fmtTick(t.v, model.ydec)}</text>)}
          {model.xt.map((t, i) => <text key={`xl${i}`} x={t.p} y={model.y0 + 15} textAnchor="middle">{fmtTick(t.v, model.xdec)}</text>)}
        </g>
        {xLabel && <text x={(model.x0 + model.x1) / 2} y={H - 5} textAnchor="middle" fill="var(--dim)" fontFamily={font} fontSize="11">{xLabel}</text>}
        {yLabel && <text x="13" y={(model.y0 + model.y1) / 2} textAnchor="middle" fill="var(--dim)" fontFamily={font} fontSize="11" transform={`rotate(-90 13 ${(model.y0 + model.y1) / 2})`}>{yLabel}</text>}

        {/* series (clipped to the plot so zoomed lines don't spill over axes) */}
        <g clipPath={`url(#${clip})`}>
          {series.map((s, i) => {
            if (s.marker) {
              // Every measured point (no downsampling); only those within the
              // current x-domain are emitted so a deep zoom stays light.
              const c = [];
              for (let k = 0; k < s.points.length; k++) {
                const [x, y] = s.points[k];
                if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
                const px = model.xOf(x);
                if (px < model.x0 - 4 || px > model.x1 + 4) continue;
                c.push(<circle key={k} cx={px.toFixed(1)} cy={model.yOf(y).toFixed(1)} r={s.r || 1.9} fill="none" stroke={s.color} strokeWidth="1" vectorEffect="non-scaling-stroke" />);
              }
              return <g key={i}>{c}</g>;
            }
            return <path key={i} d={model.paths[i]} fill="none" stroke={s.color} strokeWidth={s.width || 1.4}
              strokeDasharray={s.dashed ? '4 2.5' : ''} vectorEffect="non-scaling-stroke" />;
          })}
        </g>

        {/* drag-to-zoom selection rectangle */}
        {drag && (
          <rect x={Math.min(drag.x0, drag.x1)} y={Math.min(drag.y0, drag.y1)}
            width={Math.abs(drag.x1 - drag.x0)} height={Math.abs(drag.y1 - drag.y0)}
            fill="rgba(47,109,240,0.10)" stroke="var(--accent)" strokeWidth="0.8" strokeDasharray="3 2" pointerEvents="none" />
        )}

        {/* hover crosshair + value tooltip (replaces the legend while hovering) */}
        {hover && !drag ? (
          <g fontFamily={font} pointerEvents="none">
            <line x1={hover.mx} y1={model.y1} x2={hover.mx} y2={model.y0} stroke="var(--faint)" strokeWidth="0.8" strokeDasharray="3 2" />
            <g clipPath={`url(#${clip})`}>
              {hover.items.map((it, i) => <circle key={i} cx={it.cx} cy={it.cy} r="3" fill={it.color} stroke="var(--card)" strokeWidth="0.8" />)}
            </g>
            {(() => {
              const lines = [`${xsym} = ${fmt(hover.x)}`, ...hover.items.map(it => `${it.name}  ${fmt(it.vy)}`)];
              const tw = Math.max(64, 12 + lines.reduce((m, l) => Math.max(m, l.length), 0) * 5.2);
              const th = 8 + lines.length * 12;
              let tx = hover.mx + 10; if (tx + tw > model.x1) tx = hover.mx - tw - 10; if (tx < model.x0) tx = model.x0 + 2;
              const ty = model.y1 + 4;
              return (
                <g>
                  <rect x={tx} y={ty} width={tw} height={th} rx="4" fill="var(--card)" stroke="var(--border)" strokeWidth="0.8" opacity="0.97" />
                  {lines.map((ln, i) => (
                    <text key={i} x={tx + 7} y={ty + 13 + i * 12} fontSize="10" fontWeight={i === 0 ? 600 : 400} fill={i === 0 ? 'var(--ink)' : hover.items[i - 1].color}>{ln}</text>
                  ))}
                </g>
              );
            })()}
          </g>
        ) : (!drag && series.length > 0 && (
          <g fontFamily={font} fontSize="10">
            {series.map((s, i) => (
              <g key={`lg${i}`} transform={`translate(${model.x1 - 78}, ${model.y1 + 7 + i * 12})`}>
                {s.marker
                  ? <circle cx={7} cy={3.5} r={2.2} fill="none" stroke={s.color} strokeWidth="1" />
                  : <line x1={0} y1={3.5} x2={14} y2={3.5} stroke={s.color} strokeWidth={s.width || 1.4} strokeDasharray={s.dashed ? '4 2.5' : ''} />}
                <text x={19} y={6.5} fill="var(--dim)">{s.name}</text>
              </g>
            ))}
          </g>
        ))}
      </svg>
    </div>
  );
}
