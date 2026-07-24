/**
 * Map a pointer event's client coordinates into an <svg>'s viewBox user space.
 *
 * Why not `svg.getScreenCTM()`? The whole app renders inside `.rnr` which
 * carries a CSS `zoom: 1.15` (see index.css). Chromium's `getScreenCTM()` does
 * NOT fold an ancestor's `zoom` into its matrix, but `event.clientX/clientY`
 * DO reflect it — so a CTM-based mapping disagrees with the cursor by the zoom
 * factor and the picked point visibly drifts (worse the farther you click from
 * the origin).
 *
 * `getBoundingClientRect()` lives in the SAME (zoom-consistent) space as the
 * event coords, so deriving the transform from it cancels the zoom exactly —
 * the same robust pattern BrillouinZoneViewer already uses for raycasting.
 *
 * Assumes the default preserveAspectRatio ("xMidYMid meet"): a uniform scale
 * that fits the viewBox inside the element, centred, letterboxed on the longer
 * axis. That makes this correct whether or not the element's box matches the
 * viewBox aspect ratio (e.g. LineChart's fixed-height, letterboxed case).
 *
 * @returns {{x:number,y:number}|null} point in viewBox user units, or null.
 */
export function clientToSvg(svg, clientX, clientY) {
  if (!svg) return null;
  const rect = svg.getBoundingClientRect();
  if (!rect.width || !rect.height) return null;   // hidden / not laid out yet

  const vb = svg.viewBox?.baseVal;
  const vbX = vb ? vb.x : 0, vbY = vb ? vb.y : 0;
  const vbW = vb && vb.width ? vb.width : rect.width;
  const vbH = vb && vb.height ? vb.height : rect.height;

  // xMidYMid meet: uniform fit-inside scale, then centre the letterbox slack.
  const scale = Math.min(rect.width / vbW, rect.height / vbH);
  const offX = (rect.width - vbW * scale) / 2;
  const offY = (rect.height - vbH * scale) / 2;

  return {
    x: vbX + (clientX - rect.left - offX) / scale,
    y: vbY + (clientY - rect.top - offY) / scale,
  };
}
