import { analyzeBravais, primToConv, fracToCart } from '../src/math/bravais.js';
import { buildBZModel } from '../src/math/highsym.js';
import { conventionalLattice } from '../src/math/reciprocal.js';

let fail = 0;
function test(name, A_conv, basis, expectCode) {
  const br = analyzeBravais(A_conv, basis);
  const m = buildBZModel(br);
  const facesBySides = m.bz.faces.reduce((a, f) => { a[f.length] = (a[f.length] || 0) + 1; return a; }, {});
  const codeOk = br.code === expectCode;
  console.log(`${codeOk ? 'ok  ' : 'FAIL'} ${name}: ${br.code} (${br.system}/${br.centering}) BZ ${m.bz.vertices.length}v ${m.bz.faces.length}f ${JSON.stringify(facesBySides)}`);
  if (!codeOk) fail++;
  // Cartesian k-vector preserved under prim->conv transform:
  //   frac_prim via B_prim  ==  frac_conv via B_conv  (Å⁻¹).
  let maxErr = 0;
  for (const [, p] of Object.entries(m.points)) {
    const c1 = p.cart;                                  // = fracToCart(frac, B_prim)
    const c2 = fracToCart(primToConv(p.frac, br.T), br.B_conv);
    maxErr = Math.max(maxErr, Math.hypot(c1[0] - c2[0], c1[1] - c2[1], c1[2] - c2[2]));
  }
  const ok = maxErr < 1e-9;
  console.log(`     ${ok ? 'ok  ' : 'FAIL'} cartesian k-vector preserved (maxErr=${maxErr.toExponential(2)})`);
  if (!ok) fail++;
}
// conventionalLattice must divide each supercell VECTOR by its own repeat count
// (A_super = diag(dim)·A_conv), not each component by dim[j] — the two differ
// for non-orthogonal cells with anisotropic dim (regression: monoclinic 4×3×2).
{
  const am = 5, bm = 6, cm = 7, beta = 100 * Math.PI / 180;
  const Am = [[am, 0, 0], [0, bm, 0], [cm * Math.cos(beta), 0, cm * Math.sin(beta)]];
  const dims = [4, 3, 2];
  const V = [0, 1, 2].map(i => Am[i].map(x => x * dims[i]));
  const got = conventionalLattice(V[0], V[1], V[2], dims);
  let err = 0;
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) err = Math.max(err, Math.abs(got[i][j] - Am[i][j]));
  const ok = err < 1e-9;
  console.log(`${ok ? 'ok  ' : 'FAIL'} conventionalLattice: monoclinic 4×3×2 supercell recovers A_conv (maxErr=${err.toExponential(2)})`);
  if (!ok) fail++;
}

const a = 4.0, D = [[a, 0, 0], [0, a, 0], [0, 0, a]];
test('Simple cubic', D, [{ el: 'A', frac: [0, 0, 0] }], 'CUB');
test('FCC', D, [{ el: 'A', frac: [0, 0, 0] }, { el: 'A', frac: [0, .5, .5] }, { el: 'A', frac: [.5, 0, .5] }, { el: 'A', frac: [.5, .5, 0] }], 'FCC');
test('BCC', D, [{ el: 'A', frac: [0, 0, 0] }, { el: 'A', frac: [.5, .5, .5] }], 'BCC');
test('Hexagonal', [[3, 0, 0], [-1.5, 2.598, 0], [0, 0, 5]], [{ el: 'A', frac: [0, 0, 0] }], 'HEX');
console.log(fail === 0 ? '\n✅ BZ checks passed' : `\n❌ ${fail} failed`);
process.exit(fail ? 1 : 0);
