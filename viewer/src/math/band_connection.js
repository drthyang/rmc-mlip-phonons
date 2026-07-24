import munkres from 'munkres-js';
import { eigh } from './diagonalize';

/**
 * Calculates the absolute overlap |v1.conj().T @ v2| between two complex vectors.
 */
function getOverlap(v1Real, v1Imag, v2Real, v2Imag) {
  let sumReal = 0;
  let sumImag = 0;
  for (let i = 0; i < v1Real.length; i++) {
    // v1.conj() * v2 = (r1 - i * i1) * (r2 + i * i2)
    // = (r1*r2 + i1*i2) + i * (r1*i2 - i1*r2)
    sumReal += v1Real[i] * v2Real[i] + v1Imag[i] * v2Imag[i];
    sumImag += v1Real[i] * v2Imag[i] - v1Imag[i] * v2Real[i];
  }
  return Math.sqrt(sumReal * sumReal + sumImag * sumImag);
}

/** Groups of mode indices whose frequencies are within `tol` (legacy _degenerate_groups). */
function degenerateGroups(freqs, tol) {
  const n = freqs.length;
  const visited = new Array(n).fill(false);
  const groups = [];
  for (let i = 0; i < n; i++) {
    if (visited[i]) continue;
    const grp = [i];
    visited[i] = true;
    for (let j = i + 1; j < n; j++) {
      if (!visited[j] && isFinite(freqs[i]) && isFinite(freqs[j]) && Math.abs(freqs[i] - freqs[j]) <= tol) {
        grp.push(j);
        visited[j] = true;
      }
    }
    groups.push(grp);
  }
  return groups;
}

/**
 * SVD-rotate the current eigenvectors within one degenerate subspace to best
 * align with the previous q-point's basis (legacy Writers.connect_bands step 2).
 *   M = P_dag @ C  (g x g),  M = U S V_dag,  R = V U_dag,  C_new = C @ R.
 * Complex SVD is done via the eigendecomposition of M_dag M (Hermitian) using
 * the existing complex `eigh`. Returns new {real,imag} arrays for the group, or
 * null on any numerical issue (caller then keeps the un-rotated vectors).
 */
function rotateDegenerate(prevG, currG) {
  const g = currG.length;
  const N = currG[0].real.length;

  // M[i][j] = <P_i | C_j> = sum_n conj(P_i[n]) C_j[n]
  const Mre = new Float64Array(g * g), Mim = new Float64Array(g * g);
  for (let i = 0; i < g; i++) {
    const pr = prevG[i].real, pi = prevG[i].imag;
    for (let j = 0; j < g; j++) {
      const cr = currG[j].real, ci = currG[j].imag;
      let re = 0, im = 0;
      for (let n = 0; n < N; n++) {
        re += pr[n] * cr[n] + pi[n] * ci[n];
        im += pr[n] * ci[n] - pi[n] * cr[n];
      }
      Mre[i * g + j] = re; Mim[i * g + j] = im;
    }
  }

  // H = M_dag @ M  (Hermitian g x g): H[a][b] = sum_k conj(M[k][a]) M[k][b]
  const Hre = new Float64Array(g * g), Him = new Float64Array(g * g);
  for (let a = 0; a < g; a++) {
    for (let b = 0; b < g; b++) {
      let re = 0, im = 0;
      for (let k = 0; k < g; k++) {
        const mar = Mre[k * g + a], mai = Mim[k * g + a]; // M[k][a]
        const mbr = Mre[k * g + b], mbi = Mim[k * g + b]; // M[k][b]
        // conj(M[k][a]) * M[k][b] = (mar - i mai)(mbr + i mbi)
        re += mar * mbr + mai * mbi;
        im += mar * mbi - mai * mbr;
      }
      Hre[a * g + b] = re; Him[a * g + b] = im;
    }
  }

  let evd;
  try { evd = eigh(Hre, Him, g); } catch { return null; }
  const { eigenvalues, eigenvectors } = evd; // V columns = eigenvectors[k]
  const sigma = eigenvalues.map(v => Math.sqrt(Math.max(v, 0)));
  for (const s of sigma) if (!(s > 1e-9)) return null; // near-singular: skip

  // U[:,k] = (M @ V[:,k]) / sigma_k
  const Ure = [], Uim = [];
  for (let k = 0; k < g; k++) {
    const vr = eigenvectors[k].real, vi = eigenvectors[k].imag;
    const ur = new Float64Array(g), ui = new Float64Array(g);
    for (let r = 0; r < g; r++) {
      let re = 0, im = 0;
      for (let c = 0; c < g; c++) {
        const mr = Mre[r * g + c], mi = Mim[r * g + c];
        re += mr * vr[c] - mi * vi[c];
        im += mr * vi[c] + mi * vr[c];
      }
      ur[r] = re / sigma[k]; ui[r] = im / sigma[k];
    }
    Ure.push(ur); Uim.push(ui);
  }

  // R[i][j] = sum_k V[i][k] * conj(U[j][k])
  const Rre = new Float64Array(g * g), Rim = new Float64Array(g * g);
  for (let i = 0; i < g; i++) {
    for (let j = 0; j < g; j++) {
      let re = 0, im = 0;
      for (let k = 0; k < g; k++) {
        const vir = eigenvectors[k].real[i], vii = eigenvectors[k].imag[i]; // V[i][k]
        const ujr = Ure[k][j], uji = Uim[k][j];                            // U[j][k]
        // V[i][k] * conj(U[j][k]) = (vir + i vii)(ujr - i uji)
        re += vir * ujr + vii * uji;
        im += vii * ujr - vir * uji;
      }
      Rre[i * g + j] = re; Rim[i * g + j] = im;
    }
  }

  // C_new[:,j] = sum_i C[:,i] * R[i][j]
  const out = [];
  for (let j = 0; j < g; j++) {
    const nr = new Float64Array(N), ni = new Float64Array(N);
    for (let i = 0; i < g; i++) {
      const cr = currG[i].real, ci = currG[i].imag;
      const rr = Rre[i * g + j], ri = Rim[i * g + j];
      for (let n = 0; n < N; n++) {
        nr[n] += cr[n] * rr - ci[n] * ri;
        ni[n] += cr[n] * ri + ci[n] * rr;
      }
    }
    out.push({ real: nr, imag: ni });
  }
  return out;
}

/**
 * Connects bands across k-points using Hungarian assignment on eigenvector overlaps.
 * Matches `connect_bands` in Python.
 * 
 * @param {Array<Float64Array>} phBands - Energies for each k-point
 * @param {Array<Array<{real: Float64Array, imag: Float64Array}>>} eigvecsAll - Eigenvectors for each k-point
 * @param {number} nPasses - Number of smoothing passes
 * @param {number} freqWeight - Weight for frequency penalty (optional)
 * @returns {Object} { connectedBands, connectedEigvecs }
 */
export function connectBands(phBands, eigvecsAll, nPasses = 2, freqWeight = 0.0, degenerateTol = 5e-3) {
  let currentBands = phBands.map(arr => new Float64Array(arr));
  let currentEigvecs = eigvecsAll.map(modes => modes.map(m => ({
    real: new Float64Array(m.real),
    imag: new Float64Array(m.imag)
  })));

  const nModes = phBands[0].length;

  for (let pass = 0; pass < nPasses; pass++) {
    let newBands = [new Float64Array(currentBands[0])];
    let newEigvecs = [currentEigvecs[0].map(m => ({
      real: new Float64Array(m.real),
      imag: new Float64Array(m.imag)
    }))];

    for (let qi = 1; qi < currentBands.length; qi++) {
      const evPrev = newEigvecs[qi - 1];
      const evCurr = currentEigvecs[qi];
      const freqsPrev = newBands[qi - 1];
      const freqsCurr = currentBands[qi];

      // Build cost matrix: 1 - overlap
      const costMatrix = [];
      let freqScale = 1e-12;
      
      if (freqWeight > 0) {
        for (let j = 0; j < nModes; j++) {
          if (Math.abs(freqsCurr[j]) > freqScale) {
            freqScale = Math.abs(freqsCurr[j]);
          }
        }
      }

      for (let i = 0; i < nModes; i++) { // i is prev mode
        const row = [];
        for (let j = 0; j < nModes; j++) { // j is curr mode
          let overlap = getOverlap(evPrev[i].real, evPrev[i].imag, evCurr[j].real, evCurr[j].imag);
          
          if (freqWeight > 0) {
            const deltaFreq = Math.abs(freqsPrev[i] - freqsCurr[j]);
            overlap *= Math.exp(-freqWeight * deltaFreq / freqScale);
          }
          
          // Munkres minimizes cost, so cost = 1.0 - overlap
          row.push(1.0 - overlap);
        }
        costMatrix.push(row);
      }

      const indices = munkres(costMatrix); // returns array of [row, col]

      // Extract reordered current bands and eigenvectors
      const reorderedBands = new Float64Array(nModes);
      const reorderedEigvecs = [];

      for (let i = 0; i < nModes; i++) {
        const matchingCol = indices[i][1];
        reorderedBands[i] = freqsCurr[matchingCol];
        reorderedEigvecs[i] = {
          real: new Float64Array(evCurr[matchingCol].real),
          imag: new Float64Array(evCurr[matchingCol].imag)
        };
      }

      // Step 2: SVD-rotate eigenvectors within near-degenerate subspaces to
      // improve continuity (does NOT change frequencies). Falls back to the
      // un-rotated vectors if the rotation is numerically unstable.
      if (degenerateTol > 0) {
        let freqScale = 1e-12;
        for (let i = 0; i < nModes; i++) {
          if (isFinite(reorderedBands[i])) freqScale = Math.max(freqScale, Math.abs(reorderedBands[i]));
        }
        const tol = degenerateTol * freqScale;
        for (const grp of degenerateGroups(reorderedBands, tol)) {
          if (grp.length < 2) continue;
          const prevG = grp.map(idx => evPrev[idx]);
          const currG = grp.map(idx => reorderedEigvecs[idx]);
          const rotated = rotateDegenerate(prevG, currG);
          if (rotated) grp.forEach((idx, gi) => { reorderedEigvecs[idx] = rotated[gi]; });
        }
      }

      newBands.push(reorderedBands);
      newEigvecs.push(reorderedEigvecs);
    }

    currentBands = newBands;
    currentEigvecs = newEigvecs;
  }

  return { connectedBands: currentBands, connectedEigvecs: currentEigvecs };
}
