import { ATOMIC_MASS, THZ_TO_MEV } from '../constants.js';
import { conventionalLattice, reciprocalLattice } from '../math/reciprocal.js';

/**
 * Build the phonopy band-data object (the common structure that both band.yaml
 * and the fast band.json serialize). Atoms are basis sites (one per RMC
 * reference number), element-grouped; eigenvector rows are ordered by sorted
 * reference number (uniqueRN) — so each site's 3-vector is row `rank[rn]`.
 */
export function buildBandData(baseStructure, qPoints, bands, eigvecs, kpathMeta) {
  const { atomDic, dim, v1, v2, v3, hsym_xyz, atomType, uniqueRN } = baseStructure;
  const A = conventionalLattice(v1, v2, v3, dim);
  const recip = reciprocalLattice(A);

  const rnToRow = new Map((uniqueRN || []).map((rn, r) => [rn, r]));
  const rnToFirstAtom = new Map();
  if (atomType) for (let i = 0; i < atomType.length; i++) if (!rnToFirstAtom.has(atomType[i])) rnToFirstAtom.set(atomType[i], i);

  const atomList = [];
  for (const el of Object.keys(atomDic)) for (const rn of atomDic[el]) {
    let pos = [0, 0, 0];
    const ai = rnToFirstAtom.get(rn);
    if (ai !== undefined && hsym_xyz) pos = [hsym_xyz[ai * 3] % 1, hsym_xyz[ai * 3 + 1] % 1, hsym_xyz[ai * 3 + 2] % 1];
    atomList.push({ el, rn, pos });
  }

  const nQ = qPoints.length;
  const nModes = bands[0]?.length || 0;
  const segSizes = kpathMeta?.segSizes || [nQ];
  const hsymIndex = kpathMeta?.hsymIndex || {};

  const segStarts = new Set();
  { let off = 0; for (const sz of segSizes) { segStarts.add(off); off += sz; } }
  const distances = [0];
  for (let i = 1; i < nQ; i++) {
    if (segStarts.has(i)) { distances.push(distances[i - 1]); continue; }
    const dq = [qPoints[i][0] - qPoints[i - 1][0], qPoints[i][1] - qPoints[i - 1][1], qPoints[i][2] - qPoints[i - 1][2]];
    const cx = dq[0] * recip[0][0] + dq[1] * recip[1][0] + dq[2] * recip[2][0];
    const cy = dq[0] * recip[0][1] + dq[1] * recip[1][1] + dq[2] * recip[2][1];
    const cz = dq[0] * recip[0][2] + dq[1] * recip[1][2] + dq[2] * recip[2][2];
    distances.push(distances[i - 1] + Math.sqrt(cx * cx + cy * cy + cz * cz));
  }

  const points = atomList.map(a => ({ symbol: a.el, coordinates: a.pos, mass: ATOMIC_MASS[a.el] || 0 }));
  const phonon = [];
  for (let qi = 0; qi < nQ; qi++) {
    const band = [];
    for (let m = 0; m < nModes; m++) {
      // Model energies are meV; the phonopy schema stores THz. The explicit
      // frequency_unit key below makes the file self-describing either way.
      const rec = { frequency: bands[qi][m] / THZ_TO_MEV };
      const ev = eigvecs?.[qi]?.[m];
      if (ev) rec.eigenvector = atomList.map(a => {
        const r = rnToRow.get(a.rn);
        return [[ev.real[r * 3], ev.imag[r * 3]], [ev.real[r * 3 + 1], ev.imag[r * 3 + 1]], [ev.real[r * 3 + 2], ev.imag[r * 3 + 2]]];
      });
      band.push(rec);
    }
    const p = { 'q-position': qPoints[qi], distance: distances[qi], band };
    if (hsymIndex[qi] !== undefined) p.label = hsymIndex[qi];
    phonon.push(p);
  }

  return { nqpoint: nQ, npath: segSizes.length, segment_nqpoint: segSizes, frequency_unit: 'THz', reciprocal_lattice: recip, natom: atomList.length, lattice: A, points, phonon };
}

/** Serialize band-data to a phonopy-style band.yaml string. */
export function generatePhonopyBandYaml(baseStructure, qPoints, bands, eigvecs, kpathMeta) {
  const d = buildBandData(baseStructure, qPoints, bands, eigvecs, kpathMeta);
  const f = (x) => Number(x).toFixed(10);
  const L = [];
  L.push(`nqpoint: ${d.nqpoint}`);
  L.push(`npath: ${d.npath}`);
  L.push(`frequency_unit: ${d.frequency_unit}`);
  L.push('segment_nqpoint:');
  for (const sz of d.segment_nqpoint) L.push(`- ${sz}`);
  L.push('reciprocal_lattice:');
  for (let r = 0; r < 3; r++) L.push(`- [ ${f(d.reciprocal_lattice[r][0])}, ${f(d.reciprocal_lattice[r][1])}, ${f(d.reciprocal_lattice[r][2])} ] # ${['a*', 'b*', 'c*'][r]}`);
  L.push(`natom: ${d.natom}`);
  L.push('lattice:');
  for (let r = 0; r < 3; r++) L.push(`- [ ${f(d.lattice[r][0])}, ${f(d.lattice[r][1])}, ${f(d.lattice[r][2])} ] # ${['a', 'b', 'c'][r]}`);
  L.push('points:');
  for (const p of d.points) {
    L.push(`- symbol: ${p.symbol}`);
    L.push(`  coordinates: [ ${f(p.coordinates[0])}, ${f(p.coordinates[1])}, ${f(p.coordinates[2])} ]`);
    L.push(`  mass: ${p.mass.toFixed(5)}`);
  }
  L.push('phonon:');
  for (const q of d.phonon) {
    L.push(`- q-position: [ ${q['q-position'][0].toFixed(8)}, ${q['q-position'][1].toFixed(8)}, ${q['q-position'][2].toFixed(8)} ]`);
    L.push(`  distance: ${q.distance.toFixed(8)}`);
    if (q.label !== undefined) L.push(`  label: '${q.label}'`);
    L.push('  band:');
    for (let m = 0; m < q.band.length; m++) {
      L.push(`  - # ${m + 1}`);
      L.push(`    frequency: ${Number(q.band[m].frequency).toFixed(10)}`);
      const ev = q.band[m].eigenvector;
      if (ev) {
        L.push('    eigenvector:');
        for (let a = 0; a < ev.length; a++) {
          L.push(`    - # atom ${a + 1}`);
          for (let c = 0; c < 3; c++) L.push(`      - [ ${f(ev[a][c][0])}, ${f(ev[a][c][1])} ]`);
        }
      }
    }
  }
  return L.join('\n') + '\n';
}

/** Serialize band-data to a compact JSON (≈10× faster to reload than YAML). */
export function generateBandJson(baseStructure, qPoints, bands, eigvecs, kpathMeta) {
  return JSON.stringify(buildBandData(baseStructure, qPoints, bands, eigvecs, kpathMeta));
}

export function downloadString(text, filename, type = 'text/plain') {
  const blob = new Blob([text], { type });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(url);
}
