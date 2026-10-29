// Stackwise thermodynamics engine — SantaLucia 1998 unified nearest-neighbor model
(function (root) {
  const R = 1.98720425864083; // cal/(mol·K)
  const K0 = 273.15;
  const COMP = { A: 'T', T: 'A', G: 'C', C: 'G' };
  // ΔH kcal/mol, ΔS cal/(mol·K), keyed by top strand 5'XY3'
  const NN = {
    AA: [-7.9, -22.2], TT: [-7.9, -22.2], AT: [-7.2, -20.4], TA: [-7.2, -21.3],
    CA: [-8.5, -22.7], TG: [-8.5, -22.7], GT: [-8.4, -22.4], AC: [-8.4, -22.4],
    CT: [-7.8, -21.0], AG: [-7.8, -21.0], GA: [-8.2, -22.2], TC: [-8.2, -22.2],
    CG: [-10.6, -27.2], GC: [-9.8, -24.4], GG: [-8.0, -19.9], CC: [-8.0, -19.9]
  };
  const INIT_GC = [0.1, -2.8], INIT_AT = [2.3, 4.1], SYM_S = -1.4;

  const clean = s => (s || '').toUpperCase().replace(/U/g, 'T').replace(/[^ACGT]/g, '');
  const revcomp = s => s.split('').reverse().map(b => COMP[b] || 'N').join('');
  const isComp = (a, b) => COMP[a] === b;
  const gcFrac = s => s.length ? (s.match(/[GC]/g) || []).length / s.length : 0;

  function naEqM(c) { // von Ahsen 2001 monovalent equivalent, mol/L
    const free = Math.max(0, c.mg - c.dntp);
    return (c.na + 120 * Math.sqrt(free)) / 1000;
  }
  function saltDS(c, phosphates) { // SantaLucia 1998 entropy salt correction
    const na = naEqM(c);
    return na > 0 ? 0.368 * phosphates * Math.log(na) : 0;
  }

  // Raw duplex (1 M NaCl) for a fully complementary top strand
  function duplexRaw(s) {
    let dH = 0, dS = 0;
    for (let i = 0; i < s.length - 1; i++) { const p = NN[s.substr(i, 2)]; dH += p[0]; dS += p[1]; }
    for (const end of [s[0], s[s.length - 1]]) {
      const p = (end === 'G' || end === 'C') ? INIT_GC : INIT_AT; dH += p[0]; dS += p[1];
    }
    const sym = s === revcomp(s);
    if (sym) dS += SYM_S;
    return { dH, dS, sym };
  }

  function owczarzyTm(tm1, s, c) {
    const na = c.na / 1000, mg = Math.max(0, c.mg - c.dntp) / 1000;
    const f = gcFrac(s), n = s.length, inv1 = 1 / (tm1 + K0);
    const mono = () => {
      if (na <= 0) return tm1;
      const ln = Math.log(na);
      return 1 / (inv1 + (4.29 * f - 3.95) * 1e-5 * ln + 9.40e-6 * ln * ln) - K0;
    };
    if (mg <= 0) return mono();
    const ratio = na > 0 ? Math.sqrt(mg) / na : Infinity;
    if (ratio < 0.22) return mono();
    let a = 3.92e-5, b = -9.11e-6, cc = 6.26e-5, d = 1.42e-5, e = -4.82e-4, ff = 5.25e-4, g = 8.31e-5;
    if (ratio < 6) {
      const ln = Math.log(na);
      a = 3.92e-5 * (0.843 - 0.352 * Math.sqrt(na) * ln);
      d = 1.42e-5 * (1.279 - 4.03e-3 * ln - 8.03e-3 * ln * ln);
      g = 8.31e-5 * (0.486 - 0.258 * ln + 5.25e-3 * ln * ln * ln);
    }
    const lm = Math.log(mg);
    const inv = inv1 + a + b * lm + f * (cc + d * lm) + (1 / (2 * (n - 1))) * (e + ff * lm + g * lm * lm);
    return 1 / inv - K0;
  }

  // Full oligo thermodynamics against its perfect complement
  function thermo(s, c) {
    if (s.length < 2) return null;
    const raw = duplexRaw(s);
    const dS = raw.dS + saltDS(c, s.length - 1);
    const ct = c.oligo * 1e-9;
    const x = raw.sym ? 1 : 4;
    const tmSL = (raw.dH * 1000) / (dS + R * Math.log(ct / x)) - K0;
    let tm = tmSL;
    if (c.method === 'owczarzy') {
      const tm1 = (raw.dH * 1000) / (raw.dS + R * Math.log(ct / x)) - K0;
      tm = owczarzyTm(tm1, s, c);
    }
    const T = c.temp + K0;
    return {
      dH: raw.dH, dS, dS1M: raw.dS, sym: raw.sym, tm, tmSL,
      dG: raw.dH - T * dS / 1000, dG37: raw.dH - 310.15 * dS / 1000
    };
  }

  function basicTm(s) {
    const n = s.length, gc = (s.match(/[GC]/g) || []).length, at = n - gc;
    return n < 14 ? 2 * at + 4 * gc : 64.9 + 41 * (gc - 16.4) / n;
  }

  // Molecular weight (anhydrous, 5'-OH, ssDNA) and nearest-neighbor ε260
  const MWB = { A: 313.21, T: 304.2, C: 289.18, G: 329.21 };
  const EPS1 = { A: 15.4, C: 7.4, G: 11.5, T: 8.7 };
  const EPS2 = { AA: 27.4, AC: 21.2, AG: 25.0, AT: 22.8, CA: 21.2, CC: 14.6, CG: 18.0, CT: 15.2,
    GA: 25.2, GC: 17.6, GG: 21.6, GT: 20.0, TA: 23.4, TC: 16.2, TG: 19.0, TT: 16.8 };
  function physical(s) {
    const n = s.length; let mw = -61.96; for (const b of s) mw += MWB[b];
    let eps = 0;
    if (n === 1) eps = EPS1[s];
    else { for (let i = 0; i < n - 1; i++) eps += EPS2[s.substr(i, 2)]; for (let i = 1; i < n - 1; i++) eps -= EPS1[s[i]]; }
    eps *= 1000;
    const counts = { A: 0, C: 0, G: 0, T: 0 }; for (const b of s) counts[b]++;
    return { n, mw, eps, counts, gc: gcFrac(s) * 100, nmolPerOD: 1e6 / eps, ugPerOD: mw / eps * 1000 };
  }

  // ΔG of every 5-nt window (stacks only, 1 M + salt), for the stability profile
  function windows(s, c, w = 5) {
    const out = [], T = c.temp + K0;
    for (let i = 0; i + w <= s.length; i++) {
      let dH = 0, dS = 0; const sub = s.substr(i, w);
      for (let k = 0; k < w - 1; k++) { const p = NN[sub.substr(k, 2)]; dH += p[0]; dS += p[1]; }
      dS += saltDS(c, w - 1);
      out.push(dH - T * dS / 1000);
    }
    return out;
  }

  function meltCurve(s, c, t0, t1, step = 0.5) {
    const raw = duplexRaw(s);
    const dS = raw.dS + saltDS(c, s.length - 1);
    const ct = c.oligo * 1e-9, pts = [];
    for (let t = t0; t <= t1 + 1e-9; t += step) {
      const T = t + K0, K = Math.exp(-(raw.dH * 1000 - T * dS) / (R * T));
      let th;
      if (raw.sym) { const kc = K * ct; th = (4 * kc + 1 - Math.sqrt(8 * kc + 1)) / (4 * kc); }
      else { const kc = K * ct / 2; th = (2 * kc + 1 - Math.sqrt(4 * kc + 1)) / (2 * kc); }
      pts.push([t, th]);
    }
    return pts;
  }

  function repeats(s) {
    const runs = [], di = [];
    for (const m of s.matchAll(/([ACGT])\1{3,}/g)) runs.push({ at: m.index, seq: m[0] });
    for (const m of s.matchAll(/([ACGT]{2})\1{3,}/g)) if (m[1][0] !== m[1][1]) di.push({ at: m.index, seq: m[0] });
    return { runs, di };
  }

  root.Stackwise = { clean, revcomp, thermo, basicTm, physical, windows, meltCurve, repeats, gcFrac, naEqM, NN, R };
  if (typeof module !== 'undefined') module.exports = root.Stackwise;
})(typeof window !== 'undefined' ? window : globalThis);
