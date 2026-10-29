// Stackwise secondary-structure engine: nearest-neighbor DP for dimers and hairpins
(function (root) {
  const P = root.STACKWISE_PARAMS;
  const WCP = { AT: 1, TA: 1, GC: 1, CG: 1 };
  const wc = (x, y) => WCP[x + y] === 1;
  const INIT = [200, -5.7], ATP = [2200, 6.9];
  const MAXLOOP = 16;
  const isAT = (x, y) => (x + y === 'AT' || x + y === 'TA');
  const g = (h, s, T) => h - T * s;
  const lk = (t, k) => t[k];

  function ctx(c) {
    const T = c.temp + 273.15;
    const free = Math.max(0, c.mg - c.dntp);
    const na = (c.na + 120 * Math.sqrt(free)) / 1000;
    return { T, salt: na > 0 ? 0.368 * Math.log(na) : 0 };
  }
  function add(a, b) { return [a[0] + b[0], a[1] + b[1]]; }
  function atp(x, y) { return isAT(x, y) ? ATP : [0, 0]; }

  // generic interior transition from inner/left pair (p,q) to (i,j) on strands X (top) and Y (bottom, 3'->5')
  function link(X, Y, p, q, i, j, K) {
    const a = i - p - 1, b = j - q - 1;
    let e;
    if (a === 0 && b === 0) e = lk(P.st, X[p] + X[i] + '_' + Y[q] + Y[j]);
    else if (a === 1 && b === 1) {
      if (wc(X[p + 1], Y[q + 1])) return null;
      const e1 = lk(P.st, X[p] + X[p + 1] + '_' + Y[q] + Y[q + 1]), e2 = lk(P.st, X[p + 1] + X[i] + '_' + Y[q + 1] + Y[j]);
      e = e1 && e2 ? add(e1, e2) : null;
    } else if (a === 0 || b === 0) {
      const k = a + b, ls = P.loops.bul[k - 1];
      if (ls == null) return null;
      e = [0, ls];
      if (k === 1) { const s = lk(P.st, X[p] + X[i] + '_' + Y[q] + Y[j]); if (!s) return null; e = add(e, s); }
      else e = add(add(e, atp(X[p], Y[q])), atp(X[i], Y[j]));
    } else {
      const ls = P.loops.int[a + b - 1];
      if (ls == null) return null;
      e = [0, ls];
      if (wc(X[p + 1], Y[q + 1]) || wc(X[i - 1], Y[j - 1])) return null; // loop must open and close on mismatches
      const t1 = lk(P.tmh, Y[q + 1] + Y[q] + '_' + X[p + 1] + X[p]);
      const t2 = lk(P.tmh, X[i - 1] + X[i] + '_' + Y[j - 1] + Y[j]);
      if (t1) e = add(e, t1); if (t2) e = add(e, t2);
    }
    return e ? [e[0], e[1] + K.salt] : null;
  }

  function leftEnd(X, Y, i, j) {
    const hx = i > 0, hy = j > 0;
    if (hx && hy) return lk(P.tm, X[i - 1] + X[i] + '_' + Y[j - 1] + Y[j]) || [0, 0];
    if (hx) return lk(P.dg, Y[j] + '_' + X[i] + X[i - 1]) || [0, 0];
    if (hy) return lk(P.dg, Y[j] + Y[j - 1] + '_' + X[i]) || [0, 0];
    return [0, 0];
  }
  function rightEnd(X, Y, i, j) {
    const hx = i < X.length - 1, hy = j < Y.length - 1;
    if (hx && hy) return lk(P.tm, Y[j + 1] + Y[j] + '_' + X[i + 1] + X[i]) || [0, 0];
    if (hx) return lk(P.dg, X[i] + X[i + 1] + '_' + Y[j]) || [0, 0];
    if (hy) return lk(P.dg, X[i] + '_' + Y[j] + Y[j + 1]) || [0, 0];
    return [0, 0];
  }

  // Bimolecular duplex between A (5'→3') and B (5'→3'); returns best structures
  function dimer(A, B, c, opts = {}) {
    const K = ctx(c), Y = B.split('').reverse().join(''), n = A.length, m = Y.length;
    const best = Array.from({ length: n }, () => new Array(m).fill(null));
    for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) {
      if (!wc(A[i], Y[j])) continue;
      let cur = add(add(INIT, atp(A[i], Y[j])), leftEnd(A, Y, i, j)), from = null;
      for (let p = Math.max(0, i - MAXLOOP - 1); p < i; p++) for (let q = Math.max(0, j - MAXLOOP - 1); q < j; q++) {
        const prev = best[p][q]; if (!prev) continue;
        if ((i - p - 1) + (j - q - 1) > MAXLOOP) continue;
        const e = link(A, Y, p, q, i, j, K); if (!e) continue;
        const cand = add(prev.e, e);
        if (g(cand[0], cand[1], K.T) < g(cur[0], cur[1], K.T)) { cur = cand; from = [p, q]; }
      }
      best[i][j] = { e: cur, from };
    }
    const finals = [];
    for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) {
      const b = best[i][j]; if (!b) continue;
      const e = add(add(b.e, rightEnd(A, Y, i, j)), atp(A[i], Y[j]));
      const pairs = []; let cur = [i, j];
      while (cur) { pairs.unshift(cur); cur = best[cur[0]][cur[1]].from; }
      if (pairs.length < 2) continue;
      finals.push({ dH: e[0] / 1000, dS: e[1], dG: g(e[0], e[1], K.T) / 1000, pairs });
    }
    finals.sort((x, y) => x.dG - y.dG);
    // keep distinct structures (different first pair)
    const out = [], seen = new Set();
    for (const f of finals) {
      const k = f.pairs[0].join() ; const k2 = f.pairs[f.pairs.length - 1].join();
      if (seen.has(k) || seen.has('e' + k2)) continue;
      seen.add(k); seen.add('e' + k2);
      const ct = c.oligo * 1e-9;
      f.tm = f.dH * 1000 / (f.dS + 1.98720425864083 * Math.log(ct / 4)) - 273.15;
      f.end3A = f.pairs.some(([i]) => i === n - 1);
      f.end3B = f.pairs.some(([, j]) => j === 0);
      f.A = A; f.Y = Y;
      out.push(f); if (out.length >= (opts.limit || 6)) break;
    }
    return out;
  }

  // Unimolecular hairpins on S
  function hairpin(S, c, opts = {}) {
    const K = ctx(c), n = S.length;
    const H = Array.from({ length: n }, () => new Array(n).fill(null));
    // Y view for link(): treat pair (i,j) as X=S at i, Y=S at j with Y read 3'->5' going inward => indices reversed
    const R = S.split('').reverse().join('');
    for (let span = 4; span < n; span++) for (let i = 0; i + span < n; i++) {
      const j = i + span; if (!wc(S[i], S[j])) continue;
      const L = j - i - 1;
      let cur = null, inner = null;
      if (L >= 3 && P.loops.hp[L - 1] != null) {
        let e = [0, P.loops.hp[L - 1]];
        if (L === 3) e = add(e, atp(S[i], S[j]));
        else { const t = lk(P.tm, S[j - 1] + S[j] + '_' + S[i + 1] + S[i]); if (t) e = add(e, t); }
        cur = e;
      }
      // inner closing pairs (p,q), i<p<q<j
      for (let p = i + 1; p <= Math.min(i + MAXLOOP + 1, j - 4); p++) for (let q = j - 1; q >= Math.max(p + 4, j - MAXLOOP - 1); q--) {
        const prev = H[p][q]; if (!prev) continue;
        if ((p - i - 1) + (j - q - 1) > MAXLOOP) continue;
        // map to dimer coordinates: outer pair (i, n-1-j) precedes inner pair (p, n-1-q)
        const e = link(S, R, i, n - 1 - j, p, n - 1 - q, K); if (!e) continue;
        const cand = add(prev.e, e);
        if (!cur || g(cand[0], cand[1], K.T) < g(cur[0], cur[1], K.T)) { cur = cand; inner = [p, q]; }
      }
      if (cur) H[i][j] = { e: cur, inner };
    }
    const finals = [];
    for (let i = 0; i < n; i++) for (let j = i + 4; j < n; j++) {
      const h = H[i][j]; if (!h) continue;
      let e = add(h.e, atp(S[i], S[j]));
      const hl = i > 0, hr = j < n - 1;
      let t = null;
      if (hl && hr) t = lk(P.tm, S[i - 1] + S[i] + '_' + S[j + 1] + S[j]);
      else if (hl) t = lk(P.dg, S[j] + '_' + S[i] + S[i - 1]);
      else if (hr) t = lk(P.dg, S[j] + S[j + 1] + '_' + S[i]);
      if (t) e = add(e, t);
      const pairs = []; let cur = [i, j];
      while (cur) { pairs.push(cur); cur = H[cur[0]][cur[1]].inner; }
      if (pairs.length < 2) continue;
      const dG = g(e[0], e[1], K.T) / 1000;
      finals.push({ dH: e[0] / 1000, dS: e[1], dG, tm: e[1] < 0 ? e[0] / e[1] - 273.15 : NaN, pairs, end3: pairs[0][1] >= n - 3 });
    }
    finals.sort((x, y) => x.dG - y.dG);
    const out = [], seen = new Set();
    for (const f of finals) {
      const last = f.pairs[f.pairs.length - 1].join();
      if (seen.has(last)) continue; seen.add(last);
      out.push(f); if (out.length >= (opts.limit || 6)) break;
    }
    return out;
  }

  root.StackwiseStruct = { dimer, hairpin };
  if (typeof module !== 'undefined') module.exports = root.StackwiseStruct;
})(typeof window !== 'undefined' ? window : globalThis);
