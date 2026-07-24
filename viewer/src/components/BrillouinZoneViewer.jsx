import React, { useEffect, useRef, useState } from 'react';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls';

/**
 * Primitive Brillouin-zone k-path picker (seekpath-style).
 *
 * Renders the true primitive Wigner-Seitz polyhedron (faces + edges) with the
 * standard high-symmetry points at their cartesian positions and the suggested
 * path. Clicking a point extends the path from the current tip. Reports the path
 * as label segments {from,to} plus the label→conventional-fractional map (which
 * the calculation consumes).
 *
 * Props: bzModel = { points:{label:{cart,fracConv,display}}, path:[[a,b]], bz, code }
 */
export default function BrillouinZoneViewer({ bzModel, system, onPathChange }) {
  const mountRef = useRef(null);
  const sceneApi = useRef(null);   // { drawPath } set by the build effect
  const tipRef = useRef(null);
  const [segments, setSegments] = useState([]);

  const reportRef = useRef(onPathChange);
  reportRef.current = onPathChange;

  const emit = (segs) => {
    if (!bzModel || !reportRef.current) return;
    // Emit each point in the CELL's own reciprocal-fractional coords (conventional,
    // primitive, or supercell) — the pipeline feeds these straight to the phase.
    const conv = {};
    for (const [l, p] of Object.entries(bzModel.points)) conv[l] = p.frac;
    reportRef.current(segs, conv);
  };

  const resetToDefault = () => {
    const segs = (bzModel?.path || []).map(([from, to]) => ({ from, to }));
    tipRef.current = segs.length ? segs[segs.length - 1].to : null;
    setSegments(segs); emit(segs);
  };
  const clearPath = () => { tipRef.current = null; setSegments([]); emit([]); };
  // Break the path: the next clicked point starts a NEW segment not joined to
  // the current tip (a discontinuity, rendered as "|" in the path string).
  const breakPath = () => { tipRef.current = null; sceneApi.current?.drawPath(segments); };

  // Build scene once per model.
  useEffect(() => {
    if (!mountRef.current || !bzModel) return;
    const { points, bz } = bzModel;
    const w = mountRef.current.clientWidth, h = mountRef.current.clientHeight;
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(45, w / h, 0.001, 1000);
    const maxR = Math.max(0.1, ...bz.vertices.map(v => Math.hypot(...v)));
    camera.position.set(maxR * 2.2, maxR * 1.6, maxR * 2.6);
    camera.lookAt(0, 0, 0);
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(w, h); renderer.setPixelRatio(window.devicePixelRatio);
    mountRef.current.innerHTML = '';
    mountRef.current.appendChild(renderer.domElement);
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;

    // depthWrite:false so the translucent faces blend as one tint instead of
    // showing their fan-triangulation seams as faint triangles.
    const faceMat = new THREE.MeshBasicMaterial({ color: 0x2f6df0, transparent: true, opacity: 0.08, side: THREE.DoubleSide, depthWrite: false });
    for (const face of bz.faces) {
      const pos = [];
      for (let i = 1; i < face.length - 1; i++) pos.push(...face[0], ...face[i], ...face[i + 1]);
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
      scene.add(new THREE.Mesh(g, faceMat));
    }
    const edgePos = [];
    for (const [a, b] of bz.edges) edgePos.push(...a, ...b);
    const eg = new THREE.BufferGeometry();
    eg.setAttribute('position', new THREE.Float32BufferAttribute(edgePos, 3));
    scene.add(new THREE.LineSegments(eg, new THREE.LineBasicMaterial({ color: 0x9bb0d6 })));

    const sphereGeo = new THREE.SphereGeometry(maxR * 0.03, 16, 16);
    const baseMat = new THREE.MeshBasicMaterial({ color: 0xe06a3b });
    // One colour per disconnected sub-path (a "Break" starts a new group).
    // Group 0 = the original path (red); later groups cycle through these.
    const GROUP_COLORS = [0xdc2626, 0x13a07f, 0x8b5cf6, 0xec4899, 0xd97706];
    const groupMats = GROUP_COLORS.map(c => new THREE.MeshBasicMaterial({ color: c }));
    const groupLineMats = GROUP_COLORS.map(c => new THREE.LineBasicMaterial({ color: c, linewidth: 2 }));
    const pointsGroup = new THREE.Group();
    scene.add(pointsGroup);
    for (const [label, p] of Object.entries(points)) {
      const mesh = new THREE.Mesh(sphereGeo, baseMat);
      mesh.position.set(...p.cart);
      mesh.userData = { label };
      pointsGroup.add(mesh);
      const spr = makeLabel(p.display, maxR * 0.16);
      spr.position.set(p.cart[0] * 1.14, p.cart[1] * 1.14 + maxR * 0.04, p.cart[2] * 1.14);
      scene.add(spr);
    }

    const pathGroup = new THREE.Group();
    scene.add(pathGroup);
    // Shared arrowhead (cone) marking each segment's direction.
    const arrowGeo = new THREE.ConeGeometry(maxR * 0.035, maxR * 0.11, 14);
    const YUP = new THREE.Vector3(0, 1, 0);
    const drawPath = (segs) => {
      pathGroup.clear();
      const ptGroup = new Map();   // label -> group index (for point colour)
      let group = 0;
      for (let si = 0; si < segs.length; si++) {
        const s = segs[si];
        if (si > 0 && segs[si - 1].to !== s.from) group++;   // a break starts a new colour group
        const gi = group % GROUP_COLORS.length;
        const a = points[s.from]?.cart, b = points[s.to]?.cart;
        if (!a || !b) continue;
        ptGroup.set(s.from, gi); ptGroup.set(s.to, gi);
        const va = new THREE.Vector3(...a), vb = new THREE.Vector3(...b);
        const g = new THREE.BufferGeometry().setFromPoints([va, vb]);
        pathGroup.add(new THREE.Line(g, groupLineMats[gi]));
        // Arrowhead ~55% along from→to, pointing toward `to`.
        const dir = new THREE.Vector3().subVectors(vb, va);
        const len = dir.length();
        if (len > 1e-6) {
          dir.divideScalar(len);
          const cone = new THREE.Mesh(arrowGeo, groupMats[gi]);
          cone.position.copy(va).addScaledVector(dir, len * 0.55);
          cone.quaternion.setFromUnitVectors(YUP, dir);
          pathGroup.add(cone);
        }
      }
      // Current tip: if it begins a fresh (post-break) sub-path, highlight it in the next colour.
      const tip = tipRef.current;
      if (tip && !ptGroup.has(tip)) ptGroup.set(tip, (segs.length ? group + 1 : 0) % GROUP_COLORS.length);
      pointsGroup.children.forEach(m => {
        const gi = ptGroup.get(m.userData.label);
        m.material = gi == null ? baseMat : groupMats[gi];
      });
    };
    sceneApi.current = { drawPath };
    drawPath(segments);

    const raycaster = new THREE.Raycaster();
    const mouse = new THREE.Vector2();
    const onClick = (ev) => {
      // Normalise by the canvas's actual on-screen rect (not the captured layout
      // w/h) so picking is correct under CSS zoom and after resizes.
      const rect = renderer.domElement.getBoundingClientRect();
      mouse.x = ((ev.clientX - rect.left) / rect.width) * 2 - 1;
      mouse.y = -((ev.clientY - rect.top) / rect.height) * 2 + 1;
      raycaster.setFromCamera(mouse, camera);
      const hit = raycaster.intersectObjects(pointsGroup.children)[0];
      if (!hit) return;
      const label = hit.object.userData.label;
      if (tipRef.current == null) { tipRef.current = label; drawPath(segments); return; }
      if (label === tipRef.current) return;
      const from = tipRef.current;
      tipRef.current = label;
      setSegments(prev => { const next = [...prev, { from, to: label }]; emit(next); return next; });
    };
    renderer.domElement.addEventListener('click', onClick);

    let id;
    const loop = () => { id = requestAnimationFrame(loop); controls.update(); if (!renderer.getContext().isContextLost()) renderer.render(scene, camera); };
    loop();
    const onResize = () => {
      if (!mountRef.current) return;
      const W = mountRef.current.clientWidth, H = mountRef.current.clientHeight;
      camera.aspect = W / H; camera.updateProjectionMatrix(); renderer.setSize(W, H);
    };
    window.addEventListener('resize', onResize);
    return () => {
      window.removeEventListener('resize', onResize);
      renderer.domElement.removeEventListener('click', onClick);
      cancelAnimationFrame(id);
      renderer.dispose();
      try { renderer.forceContextLoss(); } catch { /* free the WebGL context so contexts don't pile up */ }
      sceneApi.current = null;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [bzModel]);

  // Initialise path from the suggested path when the model loads.
  useEffect(() => { if (bzModel) resetToDefault(); /* eslint-disable-next-line */ }, [bzModel]);

  // Redraw path on segment change.
  useEffect(() => { sceneApi.current?.drawPath(segments); }, [segments]);

  if (!bzModel) return (
    <div style={{ flex: 1, minHeight: 268, display: 'flex', alignItems: 'center', justifyContent: 'center', color: 'var(--faint)', font: "12px 'Spline Sans'" }}>
      Load a dataset to build the Brillouin zone.
    </div>
  );

  const seq = pathLabelSequence(segments);
  let pathStr = '— none —';
  if (seq.length) {
    pathStr = '';
    for (let i = 0; i < seq.length; i++) {
      const p = seq[i];
      if (p === '|') { pathStr += ' | '; continue; }
      if (i > 0 && seq[i - 1] !== '|') pathStr += '→';
      pathStr += (bzModel.points[p]?.display || p);
    }
  }

  return (
    <div style={{ flex: 1, minWidth: 0, display: 'flex', flexDirection: 'column', background: 'var(--card)', border: '1px solid var(--border)', borderRadius: 10, overflow: 'hidden', position: 'relative' }}>
      <div style={{ position: 'absolute', top: 14, left: 16, zIndex: 2, font: "600 13px 'Space Grotesk'", letterSpacing: '.01em', color: 'var(--ink)' }}>Brillouin zone</div>
      <span style={{ position: 'absolute', top: 17, right: 16, zIndex: 2, font: "10px 'Space Mono'", color: 'var(--faint)' }}>
        {bzModel.code}{system ? ` · ${system}` : ''}
      </span>
      <div style={{ flex: 1, minHeight: 268, background: 'var(--inset)', position: 'relative' }}>
        <div ref={mountRef} style={{ position: 'absolute', inset: 0, cursor: 'crosshair' }} />
      </div>
      <div style={{ display: 'flex', gap: 9, padding: '9px 12px 9px 16px', borderTop: '1px solid var(--border)', font: "11px 'Space Mono'", color: 'var(--dim)', alignItems: 'center' }}>
        <span style={{ color: 'var(--faint)' }}>path</span>
        <span style={{ color: 'var(--ink)', fontFamily: "'Noto Sans', sans-serif", fontWeight: 600, letterSpacing: '.02em', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>{pathStr}</span>
        <div style={{ marginLeft: 'auto', display: 'flex', gap: 6 }}>
          <button onClick={resetToDefault} className="rnr-btn" style={{ background: 'var(--soft)', color: 'var(--accentInk)', border: 'none', borderRadius: 6, padding: '5px 11px', font: "600 11px 'Space Grotesk'", cursor: 'pointer' }}>Default path</button>
          <button onClick={breakPath} title="Start a disconnected segment — the next point you click won't join the current one" className="rnr-btn" style={{ background: 'transparent', color: 'var(--dim)', border: '1px solid var(--border)', borderRadius: 6, padding: '5px 11px', font: "600 11px 'Space Grotesk'", cursor: 'pointer' }}>Break</button>
          <button onClick={clearPath} className="rnr-btn" style={{ background: 'transparent', color: 'var(--dim)', border: '1px solid var(--border)', borderRadius: 6, padding: '5px 11px', font: "600 11px 'Space Grotesk'", cursor: 'pointer' }}>Clear</button>
        </div>
      </div>
    </div>
  );
}

function pathLabelSequence(segments) {
  const out = [];
  for (let i = 0; i < segments.length; i++) {
    const s = segments[i];
    if (i === 0) out.push(s.from);
    else if (segments[i - 1].to !== s.from) out.push('|', s.from);
    out.push(s.to);
  }
  return out;
}

function makeLabel(text, size) {
  const cv = document.createElement('canvas');
  cv.width = 128; cv.height = 128;
  const ctx = cv.getContext('2d');
  ctx.fillStyle = '#2257cf'; ctx.font = "bold 80px 'Noto Sans', sans-serif"; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
  ctx.fillText(text, 64, 64);
  const tex = new THREE.CanvasTexture(cv);
  const spr = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true }));
  spr.scale.set(size, size, size);
  return spr;
}
