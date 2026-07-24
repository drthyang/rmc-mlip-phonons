import { phononDOS } from '../src/math/dos.js';
let fail = 0;
const ok = (c, m) => { if (!c) { console.error('  FAIL ' + m); fail++; } else console.log('  ok   ' + m); };
// 3 modes at 10 meV, 1 at 20 meV -> peak at 10 should be ~3x peak at 20.
const { E, dos, dosMax, count } = phononDOS([10, 10, 10, 20], { sigma: 0.5, Emin: 0, Emax: 30, nE: 600 });
ok(count === 4, 'counts 4 finite modes');
const at = (e) => dos[Math.round((e - 0) / (30 / (600 - 1)))];
ok(Math.abs(at(10) / at(20) - 3) < 0.1, `peak ratio g(10)/g(20) ≈ 3 (got ${(at(10) / at(20)).toFixed(2)})`);
ok(at(15) < 0.01 * dosMax, 'gap region ~0 between peaks');
// near-zero modes skipped (acoustic-at-Γ artifacts)
ok(phononDOS([0, 1e-7, 12], { sigma: 1, Emin: 0, Emax: 20, nE: 100 }).count === 1, 'sub-1e-6 energies skipped');
void E;
console.log(fail === 0 ? '✅ DOS ok' : `❌ ${fail} failed`);
process.exit(fail ? 1 : 0);
