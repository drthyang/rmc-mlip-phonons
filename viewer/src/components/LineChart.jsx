import React, { useMemo, useRef } from 'react';
import { clientToSvg } from './svgCoords';

/**
 * Lightweight multi-series SVG line chart (axes, gridlines, legend, optional
 * click-to-select). series = [{ name, points:[[x,y],...], color, dashed }].
 */
const W = 720, H = 320, M = { l: 52, r: 14, t: 14, b: 36 };

export default function LineChart({ series, xLabel, yLabel, onPick, height = 300, markerSeries = -1 }) {
  const ref = useRef(null);
  const model = useMemo(() => {
    let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    for (const s of series) for (const [x, y] of s.points) {
      if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
      if (x < xmin) xmin = x; if (x > xmax) xmax = x; if (y < ymin) ymin = y; if (y > ymax) ymax = y;
    }
    if (!Number.isFinite(xmin)) return null;
    if (xmax === xmin) xmax = xmin + 1;
    const pad = (ymax - ymin) * 0.05 || 1; ymin -= pad; ymax += pad;
    const pw = W - M.l - M.r, ph = H - M.t - M.b;
    const xOf = (x) => M.l + ((x - xmin) / (xmax - xmin)) * pw;
    const yOf = (y) => M.t + ph - ((y - ymin) / (ymax - ymin)) * ph;
    const paths = series.map(s => {
      let d = '';
      for (const [x, y] of s.points) { if (!Number.isFinite(x) || !Number.isFinite(y)) continue; d += (d ? 'L' : 'M') + xOf(x).toFixed(1) + ' ' + yOf(y).toFixed(1) + ' '; }
      return d;
    });
    const xticks = [], yticks = [];
    for (let i = 0; i <= 5; i++) { const x = xmin + (i / 5) * (xmax - xmin); xticks.push({ x: xOf(x), v: x }); }
    for (let i = 0; i <= 5; i++) { const y = ymin + (i / 5) * (ymax - ymin); yticks.push({ y: yOf(y), v: y }); }
    return { xOf, yOf, paths, xticks, yticks, xmin, xmax, includesZero: ymin < 0 && ymax > 0 };
  }, [series]);

  if (!model) return <div className="h-full flex items-center justify-center text-gray-500 text-xs">no data</div>;

  const fmt = (v) => Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 1 ? v.toFixed(1) : v.toFixed(2);
  const pick = (e) => {
    if (!onPick) return;
    // Map to viewBox user space (accurate under letterboxing AND the app's
    // ancestor `zoom`) — getBoundingClientRect, not getScreenCTM (see svgCoords.js).
    const u = clientToSvg(ref.current, e.clientX, e.clientY);
    if (!u) return;
    const px = u.x;
    const s = series[markerSeries >= 0 ? markerSeries : 0];
    let best = -1, bd = Infinity;
    for (let i = 0; i < s.points.length; i++) {
      const x = model.xOf(s.points[i][0]);
      const d = Math.abs(x - px);
      if (d < bd) { bd = d; best = i; }
    }
    if (best >= 0) onPick(best, s.points[best][0]);
  };

  return (
    <svg ref={ref} viewBox={`0 0 ${W} ${H}`} style={{ width: '100%', height }} onClick={pick} className={onPick ? 'cursor-pointer' : ''}>
      {model.yticks.map((t, i) => (
        <g key={`y${i}`}>
          <line x1={M.l} x2={W - M.r} y1={t.y} y2={t.y} stroke="rgba(255,255,255,0.06)" />
          <text x={M.l - 6} y={t.y + 3} textAnchor="end" fontSize="10" fill="#9ca3af">{fmt(t.v)}</text>
        </g>
      ))}
      {model.includesZero && <line x1={M.l} x2={W - M.r} y1={model.yOf(0)} y2={model.yOf(0)} stroke="rgba(255,255,255,0.2)" strokeDasharray="3 3" />}
      {model.xticks.map((t, i) => (
        <text key={`x${i}`} x={t.x} y={H - M.b + 15} textAnchor="middle" fontSize="10" fill="#9ca3af">{fmt(t.v)}</text>
      ))}
      {model.paths.map((d, i) => (
        <path key={i} d={d} fill="none" stroke={series[i].color} strokeWidth={series[i].width || 1.4} strokeDasharray={series[i].dashed ? '4 3' : ''} opacity="0.95" />
      ))}
      {/* legend */}
      {series.map((s, i) => (
        <g key={`L${i}`} transform={`translate(${M.l + 8 + i * 96}, ${M.t + 10})`}>
          <line x1={0} x2={16} y1={0} y2={0} stroke={s.color} strokeWidth="2" strokeDasharray={s.dashed ? '4 3' : ''} />
          <text x={20} y={3} fontSize="10" fill="#cbd5e1">{s.name}</text>
        </g>
      ))}
      <text x={14} y={H / 2} fontSize="11" fill="#9ca3af" transform={`rotate(-90 14 ${H / 2})`} textAnchor="middle">{yLabel}</text>
      <text x={(M.l + W - M.r) / 2} y={H - 4} fontSize="10" fill="#6b7280" textAnchor="middle">{xLabel}</text>
    </svg>
  );
}
