// Perceptual colormaps + intensity transform for the S(Q,E) heatmap and the
// Q-integrated DOS. Shared so both panels render with the identical colormap and
// contrast mapping (keeping them visually consistent). RGB stops are sampled
// from the matplotlib maps (0→1, low→high).

// turbo is deliberately omitted: its non-monotonic luminance can't be honestly
// reproduced with a few linear stops (it would misorder perceived intensity).
export const CMAPS = {
  viridis: [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]],
  inferno: [[0, 0, 4], [87, 16, 110], [188, 55, 84], [249, 142, 10], [252, 255, 164]],
  magma:   [[0, 0, 4], [81, 18, 124], [183, 55, 121], [252, 137, 97], [252, 253, 191]],
  plasma:  [[13, 8, 135], [126, 3, 168], [204, 71, 120], [248, 149, 64], [240, 249, 33]],
  cividis: [[0, 32, 76], [0, 67, 116], [87, 118, 131], [188, 186, 142], [255, 234, 70]],
  gray:    [[10, 10, 12], [128, 128, 132], [245, 245, 247]],
};

export const CMAP_NAMES = Object.keys(CMAPS);

// Linear-interpolate a colormap at t∈[0,1] → [r,g,b] (0–255, rounded).
// NaN-safe: a non-finite t clamps to 0 (never indexes past the stop array).
export function colormap(t, name = 'viridis') {
  const stops = CMAPS[name] || CMAPS.viridis;
  const tc = t > 0 ? (t < 1 ? t : 1) : 0;         // clamp to [0,1]; NaN → 0
  const x = tc * (stops.length - 1);
  const i = Math.max(0, Math.min(stops.length - 2, Math.floor(x)));
  const f = x - i;
  const a = stops[i], b = stops[i + 1];
  return [
    Math.round(a[0] + f * (b[0] - a[0])),
    Math.round(a[1] + f * (b[1] - a[1])),
    Math.round(a[2] + f * (b[2] - a[2])),
  ];
}

// CSS bottom→top gradient for a colorbar, from the same stops.
export function cmapGradient(name = 'viridis', dir = 'to top') {
  const s = CMAPS[name] || CMAPS.viridis;
  return `linear-gradient(${dir}, ${s.map(c => `rgb(${c[0]},${c[1]},${c[2]})`).join(', ')})`;
}

/**
 * Map a normalized intensity v∈[0,1] to a displayed value in [0,1], via a
 * scale mode and a contrast (gamma) knob. Shared by the heatmap and the DOS so
 * both emphasize the same features.
 *   scale: 'linear' | 'sqrt' | 'log'
 *   contrast > 1 brightens faint signal, < 1 suppresses it.
 */
export function intensityTransform(v, scale = 'log', contrast = 1) {
  let t = v > 0 ? (v < 1 ? v : 1) : 0;            // clamp to [0,1]; NaN → 0
  if (scale === 'log') t = Math.log1p(t * 1000) / Math.log1p(1000);
  else if (scale === 'sqrt') t = Math.sqrt(t);
  // contrast as gamma: exponent 1/contrast (contrast>1 ⇒ exponent<1 ⇒ brighter).
  const g = contrast > 0 ? 1 / contrast : 1;
  const tc = t > 0 ? (t < 1 ? t : 1) : 0;
  return Math.pow(tc, g);
}
