(function () {
const S = window.Stackwise, St = window.StackwiseStruct;
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fx = (v, d = 1) => (v == null || !isFinite(v)) ? '—' : (Math.abs(v) < 0.5 * Math.pow(10, -d) ? 0 : v).toFixed(d).replace('-', '−');
const store = {
  get(k, d) { try { const v = localStorage.getItem('stackwise:' + k); return v ? JSON.parse(v) : d; } catch (e) { return d; } },
  set(k, v) { try { localStorage.setItem('stackwise:' + k, JSON.stringify(v)); } catch (e) {} }
};
const uid = () => Math.random().toString(36).slice(2, 9);
const LONG = 80, STRUCT_MAX = 160;
const COMP = { A: 'T', T: 'A', G: 'C', C: 'G' };
const reduceMotion = () => { try { return matchMedia('(prefers-reduced-motion: reduce)').matches; } catch (e) { return false; } };

/* ---------- examples & state ---------- */
function mulberry(a) { return function () { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
function makeExamples() {
  const rnd = mulberry(20261029); let s = '';
  while (s.length < 720) s += 'ACGT'[Math.floor(rnd() * 4)];
  const F = 'GAAGGTGAAGGTCGGAGTCA', R = 'TTGAGGTCAATGAAGGGGTC';
  s = s.slice(0, 40) + F + s.slice(60, 640) + S.revcomp(R) + s.slice(660);
  return [
    { id: 'exT', kind: 'seq', name: 'Example template (synthetic)', seq: s, ex: true },
    { id: 'exF', kind: 'oligo', name: 'Example F 41–60', seq: F, role: 'F', ex: true, parent: 'exT', start: 41, end: 60, strand: '+' },
    { id: 'exR', kind: 'oligo', name: 'Example R 641–660', seq: R, role: 'R', ex: true, parent: 'exT', start: 641, end: 660, strand: '-' }
  ];
}
const OLIGO_EX = [
  { id: 'ex1', kind: 'oligo', name: 'M13 forward (−20)', seq: 'GTAAAACGACGGCCAGT', role: 'F', ex: true },
  { id: 'ex2', kind: 'oligo', name: 'M13 reverse', seq: 'CAGGAAACAGCTATGAC', role: 'R', ex: true },
  { id: 'ex3', kind: 'oligo', name: 'T7 promoter', seq: 'TAATACGACTCACTATAGGG', role: 'F', ex: true },
  { id: 'ex4', kind: 'oligo', name: 'Self-complementary test', seq: 'GCTAGCTTAGGCCTAAGC', role: '', ex: true }
];
let lib = store.get('lib', null);
let sel = store.get('sel', null);
if (!Array.isArray(lib) || !lib.length) { lib = makeExamples().concat(OLIGO_EX); sel = 'exT'; }
else {
  lib.forEach(o => { if (!o.kind) o.kind = (o.seq || '').length > LONG ? 'seq' : 'oligo'; });
  if (store.get('v', 1) < 2 && !lib.some(o => o.kind === 'seq')) { lib = makeExamples().concat(lib); sel = 'exT'; }
}
store.set('v', 2);
let partner = store.get('partner', null);
let cond = Object.assign({ na: 50, mg: 1.5, dntp: 0.2, oligo: 250, temp: 37, method: 'santalucia' }, store.get('cond', {}));
const PANES = [
  { id: 'melt', title: 'Melting & thermodynamics' },
  { id: 'checks', title: 'Primer checks' },
  { id: 'hairpin', title: 'Hairpins' },
  { id: 'self', title: 'Self-dimers' },
  { id: 'cross', title: 'Cross-dimer & pair' },
  { id: 'stab', title: '3′ stability profile' },
  { id: 'blast', title: 'NCBI BLAST' }
];
let order = store.get('order', PANES.map(p => p.id)).filter(id => PANES.some(p => p.id === id));
PANES.forEach(p => { if (!order.includes(p.id)) order.push(p.id); });
let closed = new Set(store.get('closed', []));
let wide = new Set(store.get('wide', ['blast']));
if (!store.get('stabNarrow', false)) { wide.delete('stab'); store.set('stabNarrow', true); }
if (!store.get('seenBlast', false)) { wide.add('blast'); store.set('seenBlast', true); }
const DRAFT = '__draft';
let draft = null;
const save = () => { store.set('lib', lib); store.set('sel', sel === DRAFT ? ((draft && draft.parent) || (lib[0] && lib[0].id)) : sel); store.set('partner', partner); store.set('cond', cond); store.set('order', order); store.set('closed', [...closed]); store.set('wide', [...wide]); };
const cur = () => (sel === DRAFT && draft) || lib.find(o => o.id === sel) || lib[0];
const byId = id => id === DRAFT ? draft : lib.find(o => o.id === id);

/* ---------- severity rules ---------- */
// Conservative limits (kcal/mol at the ΔG temperature). Four grades, each step needing a more negative ΔG:
// ok > first limit ≥ warn > second ≥ bad > third ≥ crit. A 3′ end inside the structure uses the stricter set,
// because a polymerase can extend it.
const STRUCT_LIMITS = {
  hairpin: { end3: [-1, -2, -3], inner: [-2, -3, -4.5] },
  dimer: { end3: [-3, -5, -7], inner: [-5, -7, -9] }
};
const GRADES = ['ok', 'warn', 'bad', 'crit'];
const GRADE_LABEL = { ok: 'Negligible', warn: 'Moderate', bad: 'Strong', crit: 'Severe' };
// Share of strands held in the structure at the ΔG temperature and oligo concentration.
// Hairpin: f = 1 / (1 + e^(ΔG/RT)). Dimer: solves the binding equilibrium (self: 2A ⇌ A₂, cross: A + B ⇌ AB, equal amounts).
function boundFraction(dG, kind) {
  const RT = 1.98720425864083e-3 * (cond.temp + 273.15), Ct = Math.max(cond.oligo, 1e-6) * 1e-9;
  const x = Math.min(dG / RT, 700);
  if (kind === 'hairpin') return 1 / (1 + Math.exp(x));
  const K = Math.exp(-x);
  if (kind === 'self') { const free = 2 * Ct / (1 + Math.sqrt(1 + 8 * K * Ct)); return 1 - free / Ct; }
  const free = 2 * Ct / (1 + Math.sqrt(1 + 4 * K * Ct)); return 1 - free / Ct;
}
const rateStruct = (dG, end3, kind) => {
  const lim = STRUCT_LIMITS[kind === 'hairpin' ? 'hairpin' : 'dimer'][end3 ? 'end3' : 'inner'];
  let g = dG > lim[0] ? 0 : dG > lim[1] ? 1 : dG > lim[2] ? 2 : 3;
  const f = boundFraction(dG, kind);
  if (kind === 'hairpin') { if (f >= 0.5) g = Math.max(g, 1); if (f >= 0.9 && end3) g = Math.max(g, 2); }
  else { if (f >= 0.01) g = Math.max(g, 1); if (f >= 0.1) g = Math.max(g, 2); if (f >= 0.3) g = Math.max(g, 3); }
  return GRADES[g];
};
const fracText = (dG, kind) => {
  const f = boundFraction(dG, kind) * 100;
  const v = f >= 99.95 ? '>99.9' : f < 0.1 ? '<0.1' : f < 10 ? f.toFixed(1) : f.toFixed(0);
  return kind === 'hairpin' ? `≈${v}% folded at ${cond.temp} °C` : `≈${v}% of strands paired at ${cond.temp} °C, ${cond.oligo} nM`;
};
const tmClass = tm => tm >= 52 && tm <= 65 ? 'ok' : tm >= 45 && tm <= 72 ? 'warn' : 'bad';
const pill = (st, t) => `<span class="pill ${st}">${t}</span>`;
const hpPill = h => h === undefined ? pill('neutral', '—') : h ? pill(rateStruct(h.dG, h.end3, 'hairpin'), fx(h.dG, 2)) : pill('ok', 'none');
const sdPill = d => d === undefined ? pill('neutral', '—') : d ? pill(rateStruct(d.dG, d.end3A || d.end3B, 'self'), fx(d.dG, 2)) : pill('ok', 'none');
const colorBase = b => `<span class="b${b}">${b}</span>`;
const colored = s => s.split('').map(colorBase).join('');

/* quick analysis with cache (used by selection bar and fragment table) */
const qcache = new Map();
function quick(seq, withStruct = true) {
  if (seq.length < 2) return null;
  const key = seq + '|' + cond.na + cond.mg + cond.dntp + cond.oligo + cond.temp + cond.method + '|' + withStruct;
  if (qcache.has(key)) return qcache.get(key);
  const th = S.thermo(seq, cond), w = S.windows(seq, cond, 5);
  const r = { th, gc: S.gcFrac(seq) * 100, w3: w.length ? w[w.length - 1] : null };
  if (withStruct && seq.length <= STRUCT_MAX) {
    r.hp = St.hairpin(seq, cond, { limit: 1 })[0] || false;
    r.sd = St.dimer(seq, seq, cond, { limit: 1 })[0] || false;
  }
  if (qcache.size > 400) qcache.clear();
  qcache.set(key, r);
  return r;
}

/* ---------- where oligos sit on a sequence ---------- */
function findOn(x, p) {
  if (!p.seq || x.seq.length < 2) return null;
  if (x.parent === p.id && x.start) {
    const a = x.start - 1, b = x.end - 1, want = x.strand === '-' ? S.revcomp(x.seq) : x.seq;
    if (p.seq.slice(a, b + 1) === want) return { a, b, strand: x.strand || '+' };
  }
  if (x.parent !== p.id && x.seq.length < 12) return null;
  let i = p.seq.indexOf(x.seq);
  if (i >= 0) return { a: i, b: i + x.seq.length - 1, strand: '+' };
  const rc = S.revcomp(x.seq); i = p.seq.indexOf(rc);
  if (i >= 0) return { a: i, b: i + rc.length - 1, strand: '-' };
  return null;
}
function seqAnnots(p) {
  const out = [];
  if (!p || !p.seq) return out;
  for (const x of lib) {
    if (x.kind !== 'oligo') continue;
    const h = findOn(x, p);
    if (h) {
      if (x.parent === p.id) { x.start = h.a + 1; x.end = h.b + 1; x.strand = h.strand; }
      out.push(Object.assign({ o: x }, h));
    }
  }
  return out.sort((u, v) => u.a - v.a || u.b - v.b);
}

/* ---------- double-stranded map ---------- */
function chWidth(host) {
  const p = document.createElement('span'); p.className = 'probe'; p.textContent = '0'.repeat(40);
  host.appendChild(p); const w = p.getBoundingClientRect().width / 40; p.remove();
  return w > 2 ? w : 8.1;
}
function fit(host, n, maxBlocks) {
  const chPx = chWidth(host), pw = String(Math.max(n, 1)).length + 1;
  const width = host.clientWidth || 600;
  const avail = width / chPx - (2 * pw + 6) - 1;
  return { chPx, pw, blocks: Math.max(1, Math.min(maxBlocks || 10, Math.floor((avail + 1) / 11))) };
}
function laneHTML(items, r0, r1, col, cls) {
  const lanes = [];
  let html = '';
  for (const it of items) {
    const a = Math.max(it.a, r0), b = Math.min(it.b, r1);
    let L = 0; while (lanes[L] != null && lanes[L] >= a) L++;
    lanes[L] = b + 1;
    const fwd = it.strand === '+';
    const h3 = fwd ? b === it.b : a === it.a;
    const left = col(a), w = col(b) + 1 - left;
    html += `<button class="bar ${fwd ? 'fwd' : 'rev'}${h3 ? ' h3' : ''}${it.o.draft ? ' dr' : ''}" data-oid="${it.o.id}" style="left:calc(${left} * var(--chpx));width:calc(${w} * var(--chpx));top:${L * 16 + 1}px" title="${esc(it.o.name)} · ${it.a + 1}–${it.b + 1} ${fwd ? '→' : '←'}">${esc(it.o.name)}</button>`;
  }
  return `<div class="${cls}" style="height:${lanes.length ? lanes.length * 16 + 2 : 0}px">${html}</div>`;
}
function renderMap(host, seq, opts = {}) {
  const n = seq.length;
  if (!n) { host.innerHTML = '<p class="empty">No sequence yet.</p>'; return { top: [], bot: [], blocks: 0 }; }
  const f = opts.fit || fit(host, n, opts.maxBlocks);
  const bpr = f.blocks * 10;
  host.style.setProperty('--pw', f.pw + 'ch');
  host.style.setProperty('--chpx', f.chPx + 'px');
  const ann = opts.annots || [];
  const parts = [];
  for (let r0 = 0; r0 < n; r0 += bpr) {
    const r1 = Math.min(n, r0 + bpr) - 1;
    const col = k => (k - r0) + Math.floor((k - r0) / 10);
    let top = '', bot = '', ru = '';
    for (let k = r0; k <= r1; k++) {
      if (k > r0 && (k - r0) % 10 === 0) { top += ' '; bot += ' '; }
      const b = seq[k], c = COMP[b];
      top += `<i class="b${b}" data-i="${k}">${b}</i>`;
      bot += `<i class="b${c}" data-i="${k}">${c}</i>`;
      if ((k + 1) % 10 === 0 || k === r1) ru += `<span style="left:calc(${col(k) + 1} * var(--chpx))">${k + 1}</span>`;
    }
    const inRow = ann.filter(x => x.b >= r0 && x.a <= r1);
    parts.push(`<div class="row">${laneHTML(inRow.filter(x => x.strand === '+'), r0, r1, col, 'af')}<div class="ru">${ru}</div>` +
      `<span class="pl">${r0 + 1}</span><span class="e1">5′</span><div class="st up">${top}</div><span class="e2">3′</span><span class="pr">${r1 + 1}</span>` +
      `<span class="e3">3′</span><div class="st dn">${bot}</div><span class="e4">5′</span>` +
      `${laneHTML(inRow.filter(x => x.strand === '-'), r0, r1, col, 'ar')}</div>`);
  }
  host.innerHTML = parts.join('');
  return { top: host.querySelectorAll('.st.up i'), bot: host.querySelectorAll('.st.dn i'), blocks: f.blocks };
}

/* ---------- duplex text rendering (structures) ---------- */
function duplexHTML(top, bot, pairs, opts = {}) {
  const T = [], B = [], M = [];
  const push = (t, b, m, ti, bi) => { T.push([t, ti]); B.push([b, bi]); M.push(m); };
  const [i0, j0] = pairs[0];
  const lead = Math.max(i0, j0);
  for (let k = 0; k < lead; k++) {
    const ti = k - (lead - i0), bi = k - (lead - j0);
    push(ti >= 0 ? top[ti] : ' ', bi >= 0 ? bot[bi] : ' ', ' ', ti >= 0 ? ti : -1, bi >= 0 ? bi : -1);
  }
  for (let k = 0; k < pairs.length; k++) {
    const [i, j] = pairs[k];
    if (k > 0) {
      const [p, q] = pairs[k - 1], a = i - p - 1, b = j - q - 1, w = Math.max(a, b);
      for (let x = 0; x < w; x++) push(x < a ? top[p + 1 + x] : '-', x < b ? bot[q + 1 + x] : '-', ' ', x < a ? p + 1 + x : -1, x < b ? q + 1 + x : -1);
    }
    push(top[i], bot[j], '|', i, j);
  }
  const [iN, jN] = pairs[pairs.length - 1];
  if (opts.loop) {
    const lp = opts.loop, half = Math.ceil(lp.length / 2);
    for (let x = 0; x < half; x++) push(lp[x] || ' ', lp[lp.length - 1 - x] && (lp.length - 1 - x) >= half ? lp[lp.length - 1 - x] : ' ', ' ', -2, -2);
  } else {
    const ta = top.length - 1 - iN, tb = bot.length - 1 - jN, w = Math.max(ta, tb);
    for (let x = 0; x < w; x++) push(x < ta ? top[iN + 1 + x] : ' ', x < tb ? bot[jN + 1 + x] : ' ', ' ', x < ta ? iN + 1 + x : -1, x < tb ? jN + 1 + x : -1);
  }
  const cell = ([c, idx], hot) => c === ' ' || c === '-' ? `<span class="lab">${c}</span>` : `<span class="b${c}${hot && hot.has(idx) ? ' hot' : ''}">${c}</span>`;
  const L = opts.labels || ['5′', '3′', '3′', '5′'];
  const tail = opts.loop ? [' ╮', ' │', ' ╯'] : [' ' + L[1], '', ' ' + L[3]];
  return `<pre class="duplex"><span class="lab">${L[0]} </span>${T.map(c => cell(c, opts.hotTop)).join('')}<span class="lab">${tail[0]}</span>\n<span class="lab">   </span><span class="pr">${M.join('')}</span><span class="lab">${tail[1]}</span>\n<span class="lab">${L[2]} </span>${B.map(c => cell(c, opts.hotBot)).join('')}<span class="lab">${tail[2]}</span></pre>`;
}
function dimerCard(d, A, B, kind) {
  const st = rateStruct(d.dG, d.end3A || d.end3B, kind);
  const Y = B.split('').reverse().join('');
  const hotTop = new Set([A.length - 1, A.length - 2, A.length - 3]);
  const hotBot = new Set([0, 1, 2]);
  const flags = [];
  if (d.end3A) flags.push(kind === 'self' ? '3′ end' : '3′ of this oligo');
  if (d.end3B && kind !== 'self') flags.push('3′ of partner');
  return `<div class="struct ${st}"><div class="struct-h"><span>ΔG <b>${fx(d.dG, 2)}</b> kcal/mol</span><span>ΔH <b>${fx(d.dH, 1)}</b></span><span>ΔS <b>${fx(d.dS, 1)}</b></span><span>${d.pairs.length} bp</span>${flags.map(f => `<span class="pill bad">${f}</span>`).join('')}<span class="hint">${fracText(d.dG, kind)}</span><span class="pill ${st}" style="margin-left:auto">${GRADE_LABEL[st]}</span></div>${duplexHTML(A, Y, d.pairs, { hotTop, hotBot })}</div>`;
}
function hairpinCard(h, s) {
  const st = rateStruct(h.dG, h.end3, 'hairpin');
  const n = s.length, Y = s.split('').reverse().join('');
  const pairs = h.pairs.map(([i, j]) => [i, n - 1 - j]);
  const [iK, jK] = h.pairs[h.pairs.length - 1];
  const loop = s.slice(iK + 1, jK);
  return `<div class="struct ${st}"><div class="struct-h"><span>ΔG <b>${fx(h.dG, 2)}</b> kcal/mol</span><span>Tm <b>${fx(h.tm, 1)}</b> °C</span><span>${h.pairs.length} bp stem · ${loop.length} nt loop</span>${h.end3 ? '<span class="pill bad">3′ end in stem</span>' : ''}<span class="hint">${fracText(h.dG, 'hairpin')}</span><span class="pill ${st}" style="margin-left:auto">${GRADE_LABEL[st]}</span></div>${duplexHTML(s.slice(0, iK + 1), Y.slice(0, n - jK), pairs, { loop, hotBot: new Set([0, 1, 2]), labels: ['5′', '', '3′', ''] })}</div>`;
}

/* ---------- charts ---------- */
function meltChart(seq, th) {
  const W = 520, H = 230, l = 40, r = 14, t = 26, b = 30;
  const t0 = Math.floor((th.tmSL - 30) / 5) * 5, t1 = t0 + 60;
  const pts = S.meltCurve(seq, cond, t0, t1, 0.5);
  const X = v => l + (v - t0) / (t1 - t0) * (W - l - r), Yv = v => t + (1 - v) * (H - t - b);
  let grid = '';
  for (let v = 0; v <= 1.0001; v += 0.25) grid += `<line class="grid" x1="${l}" x2="${W - r}" y1="${Yv(v)}" y2="${Yv(v)}"/><text x="${l - 6}" y="${Yv(v) + 3}" text-anchor="end">${Math.round(v * 100)}</text>`;
  for (let v = t0; v <= t1; v += 10) grid += `<line class="grid" x1="${X(v)}" x2="${X(v)}" y1="${t}" y2="${H - b}"/><text x="${X(v)}" y="${H - b + 14}" text-anchor="middle">${v}</text>`;
  const path = pts.map(([x, y], k) => (k ? 'L' : 'M') + X(x).toFixed(1) + ' ' + Yv(y).toFixed(1)).join('');
  const area = path + `L${X(t1)} ${Yv(0)}L${X(t0)} ${Yv(0)}Z`;
  const tmx = X(th.tm), showT = th.tm >= t0 && th.tm <= t1;
  return `<svg class="chart" viewBox="0 0 ${W} ${H}" style="max-width:720px" role="img" aria-label="Fraction of oligo in duplex versus temperature">
    ${grid}<line class="ax" x1="${l}" x2="${W - r}" y1="${H - b}" y2="${H - b}"/>
    <path d="${area}" fill="var(--accent)" fill-opacity="0.10"/>
    <path d="${path}" fill="none" stroke="var(--accent)" stroke-width="2"/>
    ${showT ? `<line x1="${tmx}" x2="${tmx}" y1="${t}" y2="${H - b}" stroke="var(--fg)" stroke-dasharray="3 3" stroke-width="1"/><circle cx="${tmx}" cy="${Yv(0.5)}" r="4" fill="var(--panel)" stroke="var(--fg)" stroke-width="1.5"/><text x="${tmx + 7}" y="${t - 8}" style="fill:var(--fg)">Tm ${fx(th.tm, 1)} °C</text>` : ''}
    <text x="${l}" y="${H - 4}">°C</text><text x="${l - 6}" y="12" text-anchor="end">%</text><text x="${l + 4}" y="12">fraction in duplex</text></svg>`;
}
function stabChart(seq) {
  const w = S.windows(seq, cond, 5);
  if (!w.length) return '<p class="empty">Needs at least 5 nt.</p>';
  const W = Math.max(300, w.length * 18 + 54), H = 170, l = 40, r = 10, t = 12, b = 26;
  const minG = Math.min(-10, ...w), Yv = v => t + (v / minG) * (H - t - b);
  const bw = (W - l - r) / w.length;
  let g = '';
  for (let v = 0; v >= minG - 0.01; v -= 2) g += `<line class="grid" x1="${l}" x2="${W - r}" y1="${Yv(v)}" y2="${Yv(v)}"/><text x="${l - 6}" y="${Yv(v) + 3}" text-anchor="end">${fx(v, 0)}</text>`;
  const bars = w.map((v, k) => {
    const last = k === w.length - 1, col = v < -9 ? 'var(--bad)' : v < -7 ? 'var(--warn)' : 'var(--accent)';
    return `<rect x="${(l + k * bw + 1.5).toFixed(1)}" y="${Yv(0)}" width="${Math.max(2, bw - 3).toFixed(1)}" height="${(Yv(v) - Yv(0)).toFixed(1)}" rx="2" fill="${col}" fill-opacity="${last ? 1 : 0.55}"><title>${k + 1}–${k + 5}: ${seq.substr(k, 5)} ΔG ${fx(v, 2)} kcal/mol</title></rect>`;
  }).join('');
  const ticks = w.map((v, k) => (k % 5 === 0 || k === w.length - 1) ? `<text x="${l + k * bw + bw / 2}" y="${H - b + 14}" text-anchor="middle">${k + 1}</text>` : '').join('');
  const lim = `<line x1="${l}" x2="${W - r}" y1="${Yv(-9)}" y2="${Yv(-9)}" stroke="var(--bad)" stroke-dasharray="4 3"/><text x="${W - r}" y="${Yv(-9) - 4}" text-anchor="end" style="fill:var(--bad)">−9 limit</text>`;
  return `<div style="overflow-x:auto"><svg class="chart" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" style="width:${W}px;max-width:none;height:${H}px" role="img" aria-label="Pentamer duplex stability along the sequence">${g}${bars}${lim}<line class="ax" x1="${l}" x2="${W - r}" y1="${Yv(0)}" y2="${Yv(0)}"/>${ticks}</svg></div>`;
}

/* ---------- oligo analysis ---------- */
function analyze(o) {
  const seq = o.seq, res = { seq };
  if (seq.length < 2) return res;
  res.th = S.thermo(seq, cond);
  res.alt = S.thermo(seq, Object.assign({}, cond, { method: cond.method === 'santalucia' ? 'owczarzy' : 'santalucia' }));
  res.ph = S.physical(seq);
  res.basic = S.basicTm(seq);
  res.win = S.windows(seq, cond, 5);
  res.rep = S.repeats(seq);
  if (seq.length <= STRUCT_MAX) { res.hp = St.hairpin(seq, cond, { limit: 4 }); res.sd = St.dimer(seq, seq, cond, { limit: 4 }); }
  return res;
}
function checks(o, r) {
  const s = o.seq, out = [], n = s.length;
  const add = (st, name, detail, val) => out.push({ st, name, detail, val });
  add(n >= 18 && n <= 30 ? 'ok' : n >= 15 && n <= 35 ? 'warn' : 'bad', 'Length', '18–30 nt for PCR primers', n + ' nt');
  const gc = r.ph.gc;
  add(gc >= 40 && gc <= 60 ? 'ok' : gc >= 30 && gc <= 70 ? 'warn' : 'bad', 'GC content', '40–60 %', fx(gc, 1) + ' %');
  add(tmClass(r.th.tm), 'Melting temperature', '52–65 °C (nearest-neighbor)', fx(r.th.tm, 1) + ' °C');
  const gc5 = (s.slice(-5).match(/[GC]/g) || []).length, endGC = /[GC]$/.test(s);
  add(gc5 > 3 ? 'warn' : endGC ? 'ok' : 'warn', 'GC clamp', 'G or C at the 3′ end, at most 3 G/C in the last 5', `${gc5}/5, ends ${s.slice(-1)}`);
  const w3 = r.win.length ? r.win[r.win.length - 1] : 0;
  add(w3 > -9 ? 'ok' : 'warn', '3′ end stability', 'ΔG of the 3′ pentamer above −9 kcal/mol', fx(w3, 2));
  add(r.rep.runs.length ? 'warn' : 'ok', 'Single-base runs', 'No run of 4 or more of one base', r.rep.runs.length ? r.rep.runs.map(x => x.seq).join(', ') : 'none');
  add(r.rep.di.length ? 'warn' : 'ok', 'Dinucleotide repeats', 'No dinucleotide repeated 4 or more times', r.rep.di.length ? r.rep.di.map(x => x.seq).join(', ') : 'none');
  if (r.hp) { const h = r.hp[0]; add(h ? rateStruct(h.dG, h.end3, 'hairpin') : 'ok', 'Hairpins', 'ΔG above −1 (3′ end in stem) / −2 kcal/mol (internal); more negative is graded harder', h ? fx(h.dG, 2) : 'none'); }
  if (r.sd) { const d = r.sd[0]; add(d ? rateStruct(d.dG, d.end3A || d.end3B, 'self') : 'ok', 'Self-dimers', 'ΔG above −3 (3′ end paired) / −5 kcal/mol (internal); more negative is graded harder', d ? fx(d.dG, 2) : 'none'); }
  return out;
}

/* ---------- library ---------- */
function renderLib() {
  const seqs = lib.filter(o => o.kind === 'seq'), oligos = lib.filter(o => o.kind !== 'seq');
  const seqItem = o => {
    const nOl = seqAnnots(o).length;
    return `<button class="item" aria-current="${o.id === sel}" data-id="${o.id}">
      <span class="role ds">ds</span><span class="nm">${esc(o.name || 'Untitled')}</span>
      <span class="tmchip num pill neutral">${o.seq.length ? fx(S.gcFrac(o.seq) * 100, 0) + '% GC' : '—'}</span>
      <span class="sub">${o.ex ? 'ex · ' : ''}${o.seq.length.toLocaleString('en-US')} bp · ${nOl} oligo${nOl === 1 ? '' : 's'}</span></button>`;
  };
  const oligoItem = o => {
    const th = o.seq.length >= 2 ? S.thermo(o.seq, cond) : null;
    const p = o.parent && byId(o.parent);
    const sub = p && o.start ? `from ${esc(p.name)} ${o.start}–${o.end} ${o.strand === '-' ? '←' : '→'}` : `${o.seq.length} nt · ${esc(o.seq.slice(0, 22))}${o.seq.length > 22 ? '…' : ''}`;
    return `<button class="item" aria-current="${o.id === sel}" data-id="${o.id}">
      <span class="role ${o.role || ''}">${o.role || '·'}</span><span class="nm">${esc(o.name || 'Untitled')}</span>
      <span class="tmchip num pill ${th ? tmClass(th.tm) : 'neutral'}">${th ? fx(th.tm, 1) + '°' : '—'}</span>
      <span class="sub">${o.ex ? 'ex · ' : ''}${sub}</span></button>`;
  };
  const draftItem = draft ? (() => {
    const p = draft.parent && byId(draft.parent);
    const sub = p && draft.start ? `from ${esc(p.name)} ${draft.start}–${draft.end} ${draft.strand === '-' ? '←' : '→'}` : `${draft.seq.length} nt`;
    return `<button class="item draft" aria-current="${sel === DRAFT}" data-id="${DRAFT}"><span class="role ${draft.role || ''}">${draft.role || '·'}</span><span class="nm">${esc(draft.name || 'Untitled')}</span><span class="tmchip pill warn">unsaved</span><span class="sub">${sub}</span></button>`;
  })() : '';
  $('lib').innerHTML = (seqs.length ? '<div class="grp">Sequences</div>' + seqs.map(seqItem).join('') : '') +
    (oligos.length || draft ? '<div class="grp">Oligos</div>' + draftItem + oligos.map(oligoItem).join('') : '');
  $('lib-count').textContent = `${seqs.length} seq · ${oligos.length} oligo${oligos.length === 1 ? '' : 's'}`;
}

/* ================= SEQUENCE VIEW ================= */
const SV = { sel: null, anchor: null, dragging: false, touchPending: false, hits: [], hitIdx: -1, editing: false, nameDirty: false, pick: '+', pending: null };
let mapEls = { top: [], bot: [], blocks: 0 }, painted = null, paintedHits = [], mapKey = '';

function hashSeq(s) { let h = 2166136261; for (let i = 0; i < s.length; i++) h = Math.imul(h ^ s.charCodeAt(i), 16777619); return h >>> 0; }
function renderSeqMap(force) {
  const o = cur(), host = $('s-map');
  const ann = seqAnnots(o);
  if (draft && draft.parent === o.id) { const h = findOn(draft, o); if (h) ann.push(Object.assign({ o: draft }, h)); ann.sort((u, v) => u.a - v.a); }
  const f = o.seq.length ? fit(host, o.seq.length, 10) : null;
  const key = [o.id, o.seq.length, hashSeq(o.seq), f && f.blocks, f && f.chPx.toFixed(2), ann.map(x => x.o.id + ':' + x.a + ':' + x.b + x.strand + x.o.name).join()].join('|');
  if (!force && key === mapKey) return ann;
  mapKey = key;
  mapEls = renderMap(host, o.seq, { annots: ann, fit: f });
  painted = null; paintedHits = [];
  paintSel(); paintHits();
  return ann;
}
function paintSel() {
  const { top, bot } = mapEls;
  if (painted) for (let k = painted.a; k <= painted.b; k++) { if (top[k]) top[k].classList.remove('sel', 'sels'); if (bot[k]) bot[k].classList.remove('sel', 'sels'); }
  painted = SV.sel ? Object.assign({}, SV.sel) : null;
  if (painted) for (let k = painted.a; k <= painted.b; k++) {
    if (top[k]) top[k].classList.add(painted.strand === '+' ? 'sels' : 'sel');
    if (bot[k]) bot[k].classList.add(painted.strand === '-' ? 'sels' : 'sel');
  }
}
function paintHits() {
  paintedHits.forEach(el => el.classList.remove('hit'));
  paintedHits = [];
  for (const h of SV.hits) for (let k = h.a; k <= h.b; k++) {
    const el = (h.strand === '+' ? mapEls.top : mapEls.bot)[k];
    if (el) { el.classList.add('hit'); paintedHits.push(el); }
  }
}
function fragName(o, s) { const short = (o.name || 'Seq').split(/[\s(,]/)[0].slice(0, 16) || 'Seq'; return `${short} ${s.strand === '+' ? 'F' : 'R'} ${s.a + 1}–${s.b + 1}`; }
function selSeq() { const o = cur(), s = SV.sel, f = o.seq.slice(s.a, s.b + 1); return s.strand === '+' ? f : S.revcomp(f); }
function setPick(st) { SV.pick = st; $('s-plus').setAttribute('aria-pressed', st === '+'); $('s-minus').setAttribute('aria-pressed', st === '-'); }
function setSel(a, b, strand) {
  const n = cur().seq.length; if (!n) return;
  const lo = Math.max(0, Math.min(a, b)), hi = Math.min(n - 1, Math.max(a, b));
  const prev = SV.sel;
  if (!prev || prev.a !== lo || prev.b !== hi || prev.strand !== strand) SV.nameDirty = false;
  SV.sel = { a: lo, b: hi, strand };
  $('s-from').value = lo + 1; $('s-to').value = hi + 1; setPick(strand);
  paintSel(); renderSelBar(false);
}
function clearSel() { SV.sel = null; SV.touchPending = false; paintSel(); renderSelBar(); }
function scrollToSel() {
  const s = SV.sel; if (!s) return;
  const el = mapEls.top[s.a]; if (el) el.scrollIntoView({ block: 'center', behavior: reduceMotion() ? 'auto' : 'smooth' });
}

function renderSelBar(full) {
  const bar = $('selbar'), s = SV.sel;
  if (!s) { bar.hidden = true; bar.innerHTML = ''; return; }
  const o = cur(), q = selSeq(), n = q.length;
  const r = n >= 2 ? quick(q, !!full) : null;
  let struct = '';
  if (r && full && n <= STRUCT_MAX) struct = `<span class="hint">Hairpin</span>${hpPill(r.hp)}<span class="hint">Self-dimer</span>${sdPill(r.sd)}`;
  else if (n > STRUCT_MAX) struct = `<span class="hint">Structure search runs up to ${STRUCT_MAX} nt</span>`;
  const defName = fragName(o, s);
  const prevInput = $('sel-name');
  const nameVal = SV.nameDirty && prevInput ? prevInput.value : defName;
  const shown = n > 70 ? colored(q.slice(0, 32)) + '<span class="lab"> … </span>' + colored(q.slice(-32)) : colored(q);
  const fwd = s.strand === '+';
  bar.hidden = false;
  bar.innerHTML = `<div class="row2"><span class="dir ${fwd ? 'fwd' : 'rev'}">${fwd ? '→ Top strand' : '← Bottom strand'}</span>
      <b class="mono num">${s.a + 1}–${s.b + 1}</b><span class="hint">${n} nt</span>
      ${r ? pill(tmClass(r.th.tm), 'Tm ' + fx(r.th.tm, 1) + ' °C') + pill('neutral', 'GC ' + fx(r.gc, 0) + ' %') + pill('neutral', 'ΔG ' + fx(r.th.dG, 2)) : '<span class="hint">Select at least 2 bases</span>'}
      ${struct}<span class="sp"></span>${SV.touchPending ? '<span class="hint">Tap the last base</span>' : ''}<button class="btn ghost" data-sb="clear">Clear</button></div>
    <div class="sq"><span class="lab">5′-</span>${shown}<span class="lab">-3′</span></div>
    <div class="row2"><div class="field grow"><label for="sel-name">Name</label><input id="sel-name" value="${esc(nameVal)}" spellcheck="false"></div>
      <button class="btn" data-sb="flip">Other strand</button>
      <button class="btn" data-sb="copy">Copy</button>
      <button class="btn" data-sb="save">Save oligo</button>
      <button class="btn primary" data-sb="analyze">Analyze</button></div>`;
}

function saveSel(open) {
  const o = cur(), s = SV.sel; if (!s) return;
  const q = selSeq();
  if (q.length < 2) { toast('Select at least 2 bases'); return; }
  const name = ($('sel-name') && $('sel-name').value.trim()) || fragName(o, s);
  let item = lib.find(x => x.kind === 'oligo' && x.parent === o.id && x.start === s.a + 1 && x.end === s.b + 1 && x.strand === s.strand);
  if (!item) { item = { id: uid(), kind: 'oligo', name, seq: q, role: s.strand === '+' ? 'F' : 'R', parent: o.id, start: s.a + 1, end: s.b + 1, strand: s.strand }; lib.push(item); toast('Saved ' + name); }
  else toast('Already saved as ' + item.name);
  save();
  if (open) { sel = item.id; renderAll(true); window.scrollTo(0, 0); return; }
  SV.sel = null; renderSelBar(); renderSeqMap(); renderFragments(); renderLib();
}

function openDraftFromSel() {
  const o = cur(), s = SV.sel; if (!s) return;
  const q = selSeq(); if (q.length < 2) { toast('Select at least 2 bases'); return; }
  const name = ($('sel-name') && $('sel-name').value.trim()) || fragName(o, s);
  draft = { id: DRAFT, kind: 'oligo', draft: true, name, seq: q, role: s.strand === '+' ? 'F' : 'R', parent: o.id, start: s.a + 1, end: s.b + 1, strand: s.strand };
  sel = DRAFT; renderAll(true); window.scrollTo(0, 0);
}
function saveDraft() {
  if (!draft) return;
  const p = draft.parent && byId(draft.parent);
  const dup = p && lib.find(x => x.kind === 'oligo' && x.parent === p.id && x.start === draft.start && x.end === draft.end && x.strand === draft.strand && x.seq === draft.seq);
  if (dup) { draft = null; sel = dup.id; save(); renderAll(true); toast('Already saved as ' + dup.name); return; }
  const item = Object.assign({}, draft, { id: uid() }); delete item.draft;
  lib.push(item); draft = null; sel = item.id; save(); renderAll(true); toast('Saved ' + item.name);
}
function nudge(end, d) {
  const o = draft, p = o && byId(o.parent); if (!p) return;
  let a = o.start - 1, b = o.end - 1; const plus = o.strand === '+';
  if ((end === 5) === plus) a -= d; else b += d;
  if (a < 0 || b >= p.seq.length) { toast('Reached the end of the sequence'); return; }
  if (b - a + 1 < 5) { toast('Oligo must stay at least 5 nt'); return; }
  const f = p.seq.slice(a, b + 1);
  if (/ [FR] \d+–\d+$/.test(o.name)) o.name = fragName(p, { a, b, strand: o.strand });
  o.seq = plus ? f : S.revcomp(f); o.start = a + 1; o.end = b + 1;
  renderAll(true);
}
function renderFragments() {
  const o = cur(), list = seqAnnots(o);
  $('frag-sum').innerHTML = list.length ? pill('neutral', list.length + ' oligo' + (list.length === 1 ? '' : 's')) : '';
  if (!list.length) {
    $('frag-body').innerHTML = '<p class="empty">No oligos here yet. Select a stretch on the map and save it. Primers already in your library show up here automatically when they match either strand (12 nt or longer).</p>';
    return;
  }
  const rows = list.map(x => {
    const r = quick(x.o.seq), fwd = x.strand === '+';
    return `<tr><td class="nm" title="${esc(x.o.name)}">${esc(x.o.name)}</td>
      <td class="n">${x.a + 1}–${x.b + 1}</td><td><span class="dir ${fwd ? 'fwd' : 'rev'}">${fwd ? '→ top' : '← bottom'}</span></td>
      <td class="n">${x.o.seq.length}</td><td>${pill(tmClass(r.th.tm), fx(r.th.tm, 1))}</td><td class="n">${fx(r.gc, 1)}</td>
      <td class="n">${fx(r.th.dG, 2)}</td><td>${hpPill(r.hp)}</td><td>${sdPill(r.sd)}</td>
      <td>${pill(r.w3 > -9 ? 'ok' : 'warn', fx(r.w3, 2))}</td>
      <td><button class="btn ghost" data-fa="show" data-oid="${x.o.id}">Show</button><button class="btn ghost" data-fa="open" data-oid="${x.o.id}">Analyze</button></td></tr>`;
  }).join('');
  const F = list.filter(x => x.strand === '+'), R = list.filter(x => x.strand === '-'), amps = [];
  for (const f of F) for (const r of R) if (r.b > f.b && r.a > f.a) amps.push({ f, r, size: r.b - f.a + 1 });
  amps.sort((u, v) => u.size - v.size);
  const ampHTML = amps.slice(0, 8).map(({ f, r, size }) => {
    const tf = quick(f.o.seq).th.tm, tr = quick(r.o.seq).th.tm, d = Math.abs(tf - tr);
    let cd = '';
    if (f.o.seq.length <= STRUCT_MAX && r.o.seq.length <= STRUCT_MAX) { const x = St.dimer(f.o.seq, r.o.seq, cond, { limit: 1 })[0]; cd = 'Cross-dimer ' + (x ? pill(rateStruct(x.dG, x.end3A || x.end3B, 'cross'), fx(x.dG, 2)) : pill('ok', 'none')); }
    return `<div><b>${esc(f.o.name)}</b><span class="hint">×</span><b>${esc(r.o.name)}</b><span class="mono num">${size.toLocaleString('en-US')} bp</span>${pill(d <= 2 ? 'ok' : d <= 5 ? 'warn' : 'bad', 'ΔTm ' + fx(d, 1) + ' °C')}${cd}<button class="btn ghost" data-blpair="${f.o.id},${r.o.id}">BLAST pair</button></div>`;
  }).join('');
  $('frag-body').innerHTML = `<div class="tablewrap"><table class="frag"><thead><tr><th>Name</th><th class="n">Position</th><th>Strand</th><th class="n">nt</th><th>Tm °C</th><th class="n">GC %</th><th class="n">ΔG ${cond.temp}°</th><th>Hairpin</th><th>Self-dimer</th><th>3′ ΔG</th><th></th></tr></thead><tbody>${rows}</tbody></table></div>` +
    (amps.length ? `<div class="amp"><span class="eyebrow">Products (forward × reverse)</span>${ampHTML}</div>` : '');
}

function seqStats(o) {
  const s = o.seq, n = s.length;
  if (!n) return '<span>Paste a strand to begin.</span>';
  const gc = S.gcFrac(s) * 100, c = { A: 0, C: 0, G: 0, T: 0 }; for (const b of s) c[b]++;
  const tmLong = n >= 2 ? 81.5 + 16.6 * Math.log10(S.naEqM(cond)) + 0.41 * gc - 600 / n : NaN;
  const mw = n >= 2 ? (S.physical(s).mw + S.physical(S.revcomp(s)).mw) : NaN;
  return `<span>Length <b>${n.toLocaleString('en-US')}</b> bp</span><span>GC <b>${fx(gc, 1)}</b> %</span>
    <span>Tm (long duplex, salt-adjusted) <b>${fx(tmLong, 1)}</b> °C</span><span>MW dsDNA <b>${(mw / 1000).toFixed(1)}</b> kDa</span>
    <span><b class="bA">A ${c.A}</b> · <b class="bC">C ${c.C}</b> · <b class="bG">G ${c.G}</b> · <b class="bT">T ${c.T}</b></span>`;
}

function renderSeqView(sync) {
  const o = cur();
  if (sync) {
    $('s-name').value = o.name; $('s-seq').value = o.seq;
    SV.sel = null; SV.hits = []; SV.hitIdx = -1; SV.touchPending = false; SV.editing = !o.seq;
    $('s-find').value = ''; $('s-hits').textContent = ''; $('s-next').hidden = true; $('s-msg').hidden = true;
    $('s-from').value = ''; $('s-to').value = '';
    $('s-del').textContent = 'Delete'; $('s-del').dataset.armed = '';
    mapKey = '';
  }
  if (SV.sel && SV.sel.b >= o.seq.length) SV.sel = null;
  $('s-seq').hidden = !SV.editing; $('s-editbtn').textContent = SV.editing ? 'Done editing' : 'Edit sequence';
  $('s-stats').innerHTML = seqStats(o);
  renderSeqMap(false);
  renderFragments();
  if (SV.pending) { const p = SV.pending; SV.pending = null; setSel(p.a, p.b, p.strand); renderSelBar(true); requestAnimationFrame(scrollToSel); }
  else renderSelBar(true);
}

function computeHits() {
  const q = S.clean($('s-find').value), p = cur().seq;
  SV.hits = []; SV.hitIdx = -1;
  if (q.length >= 3) {
    const rc = S.revcomp(q);
    for (const [pat, strand] of [[q, '+'], [rc, '-']]) {
      if (strand === '-' && rc === q) continue;
      let i = p.indexOf(pat);
      while (i >= 0 && SV.hits.length < 1000) { SV.hits.push({ a: i, b: i + pat.length - 1, strand }); i = p.indexOf(pat, i + 1); }
    }
    SV.hits.sort((u, v) => u.a - v.a);
  }
  paintHits();
  const nt = SV.hits.filter(h => h.strand === '+').length, nb = SV.hits.length - nt;
  $('s-hits').textContent = !q.length ? '' : q.length < 3 ? 'Type at least 3 bases' : SV.hits.length ? `${SV.hits.length} match${SV.hits.length > 1 ? 'es' : ''}: ${nt} top, ${nb} bottom` : 'No match on either strand';
  $('s-next').hidden = !SV.hits.length;
}

/* ================= OLIGO BENCH ================= */
function renderMetrics(r) {
  if (!r.th) { $('metrics').innerHTML = '<div class="m"><div class="k">Enter at least 2 nt</div></div>'; return; }
  const ph = r.ph, th = r.th, c = ph.counts;
  const m = (k, v, u = '', cls = '') => `<div class="m ${cls}"><div class="k">${k}</div><div class="v num">${v}${u ? `<small>${u}</small>` : ''}</div></div>`;
  $('metrics').innerHTML =
    m(`Tm · ${cond.method === 'owczarzy' ? 'Owczarzy' : 'SantaLucia'}`, fx(th.tm, 1), '°C', 'hero') +
    m('Length', ph.n, 'nt') + m('GC', fx(ph.gc, 1), '%') +
    m(`ΔG ${cond.temp}°`, fx(th.dG, 2), 'kcal/mol') + m('ΔH', fx(th.dH, 1), 'kcal/mol') + m('ΔS', fx(th.dS, 1), 'cal/K·mol') +
    m('MW', ph.mw.toFixed(1), 'g/mol') + m('ε260', Math.round(ph.eps).toLocaleString('en-US'), 'L/mol·cm') +
    m('nmol / OD', fx(ph.nmolPerOD, 2)) + m('µg / OD', fx(ph.ugPerOD, 2)) +
    `<div class="m"><div class="k">Composition</div><div class="v num counts" style="font-size:13px"><span class="bA">A${c.A}</span><span class="bC">C${c.C}</span><span class="bG">G${c.G}</span><span class="bT">T${c.T}</span></div></div>`;
}
function paneBody(id, o, r) {
  if (!r.th) return '<p class="empty">Enter a sequence to analyze.</p>';
  const tooLong = o.seq.length > STRUCT_MAX ? `<p class="empty">Structure search runs on sequences up to ${STRUCT_MAX} nt.</p>` : '';
  if (id === 'melt') {
    const th = r.th, alt = r.alt, other = cond.method === 'santalucia' ? 'Owczarzy' : 'SantaLucia';
    return meltChart(o.seq, th) + `<table class="kv">
      <tr><td>Tm, ${cond.method === 'santalucia' ? 'SantaLucia 1998 salt correction' : 'Owczarzy 2004/2008 salt correction'}</td><td>${fx(th.tm, 2)} °C</td></tr>
      <tr><td>Tm, ${other} salt correction</td><td>${fx(alt.tm, 2)} °C</td></tr>
      <tr><td>Tm, basic (Wallace / GC formula)</td><td>${fx(r.basic, 1)} °C</td></tr>
      <tr><td>ΔH° (nearest-neighbor sum)</td><td>${fx(th.dH, 1)} kcal/mol</td></tr>
      <tr><td>ΔS° (1 M NaCl → salt-corrected)</td><td>${fx(th.dS1M, 1)} → ${fx(th.dS, 1)} cal/K·mol</td></tr>
      <tr><td>ΔG° at ${cond.temp} °C</td><td>${fx(th.dG, 2)} kcal/mol</td></tr>
      <tr><td>ΔG° at 37 °C</td><td>${fx(th.dG37, 2)} kcal/mol</td></tr>
      <tr><td>Equivalent Na⁺ (Na⁺ + 120·√free Mg²⁺)</td><td>${fx(S.naEqM(cond) * 1000, 1)} mM</td></tr>
      <tr><td>Self-complementary</td><td>${th.sym ? 'yes (Ct used, not Ct/4)' : 'no'}</td></tr></table>`;
  }
  if (id === 'checks') {
    return `<div class="checks">${checks(o, r).map(c => `<div class="chk"><span class="pill ${c.st}">${c.st === 'ok' ? 'Pass' : c.st === 'warn' ? 'Check' : c.st === 'crit' ? 'Fail ✕' : 'Fail'}</span><div><div>${c.name}</div><div class="d">${c.detail}</div></div><div class="val">${esc(c.val)}</div></div>`).join('')}</div>`;
  }
  if (id === 'hairpin') {
    if (!r.hp) return tooLong;
    return r.hp.length ? r.hp.map(h => hairpinCard(h, o.seq)).join('') : '<p class="empty">No hairpin with a stem of 2 bp or more and a loop of 3 nt or more.</p>';
  }
  if (id === 'self') {
    if (!r.sd) return tooLong;
    return r.sd.length ? r.sd.map(d => dimerCard(d, o.seq, o.seq, 'self')).join('') : '<p class="empty">No self-dimer found.</p>';
  }
  if (id === 'cross') {
    const others = lib.filter(x => x.kind === 'oligo' && x.id !== o.id && x.seq.length >= 2);
    if (!others.length) return '<p class="empty">Add a second oligo to the library to check a primer pair.</p>';
    let p = others.find(x => x.id === partner);
    if (!p) p = others.find(x => o.role === 'F' ? x.role === 'R' : o.role === 'R' ? x.role === 'F' : false) || others[0];
    const opts = others.map(x => `<option value="${x.id}" ${x.id === p.id ? 'selected' : ''}>${esc(x.name || 'Untitled')} (${x.seq.length} nt)</option>`).join('');
    const pt = S.thermo(p.seq, cond), dTm = Math.abs(r.th.tm - pt.tm);
    let body = `<div class="partner"><label class="eyebrow" for="partner">Partner</label><select id="partner">${opts}</select>
      <span class="pill ${dTm <= 2 ? 'ok' : dTm <= 5 ? 'warn' : 'bad'}">ΔTm ${fx(dTm, 1)} °C</span>
      <span class="note">${fx(r.th.tm, 1)} vs ${fx(pt.tm, 1)} °C</span></div>`;
    if (o.seq.length > STRUCT_MAX || p.seq.length > STRUCT_MAX) return body + tooLong;
    const cd = St.dimer(o.seq, p.seq, cond, { limit: 4 });
    return body + (cd.length ? cd.map(d => dimerCard(d, o.seq, p.seq, 'cross')).join('') : '<p class="empty">No cross-dimer found.</p>');
  }
  if (id === 'blast') return blastPaneBody(o);
  if (id === 'stab') return `<p class="note" style="margin:0">ΔG of each 5-nt window paired with its complement, at ${cond.temp} °C. The solid bar is the 3′ pentamer; a very stable 3′ end raises the risk of mispriming.</p>` + stabChart(o.seq);
  return '';
}
function paneSummary(id, r) {
  if (!r.th) return '';
  if (id === 'melt') return pill(tmClass(r.th.tm), 'Tm ' + fx(r.th.tm, 1) + ' °C');
  if (id === 'checks') { const cs = checks(cur(), r), w = cs.filter(c => c.st !== 'ok').length; return pill(w ? (cs.some(c => c.st === 'crit') ? 'crit' : cs.some(c => c.st === 'bad') ? 'bad' : 'warn') : 'ok', w ? w + ' to review' : 'All pass'); }
  if (id === 'hairpin' && r.hp) return hpPill(r.hp[0] || false);
  if (id === 'self' && r.sd) return sdPill(r.sd[0] || false);
  if (id === 'blast') return blastSummary(cur());
  if (id === 'stab' && r.win.length) { const v = r.win[r.win.length - 1]; return pill(v > -9 ? 'ok' : 'warn', '3′ ' + fx(v, 2)); }
  return '';
}
const ICON = {
  grip: '<svg width="12" height="16" viewBox="0 0 12 16" aria-hidden="true"><g fill="currentColor"><circle cx="3" cy="3" r="1.4"/><circle cx="9" cy="3" r="1.4"/><circle cx="3" cy="8" r="1.4"/><circle cx="9" cy="8" r="1.4"/><circle cx="3" cy="13" r="1.4"/><circle cx="9" cy="13" r="1.4"/></g></svg>',
  wide: '<svg width="16" height="16" viewBox="0 0 16 16" aria-hidden="true"><path d="M2 8h12M2 8l3-3M2 8l3 3M14 8l-3-3M14 8l-3 3" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round"/></svg>',
  fold: '<svg width="16" height="16" viewBox="0 0 16 16" aria-hidden="true"><path d="M4 6l4 4 4-4" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round"/></svg>',
  unfold: '<svg width="16" height="16" viewBox="0 0 16 16" aria-hidden="true"><path d="M6 4l4 4-4 4" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round"/></svg>'
};
let lastRes = null;
function renderPanes(o, r) {
  const scroll = window.scrollY;
  $('panes').innerHTML = order.map(id => {
    const p = PANES.find(x => x.id === id), c = closed.has(id), w = wide.has(id);
    return `<section class="card pane ${c ? 'closed' : ''} ${w ? 'wide' : ''}" data-pane="${id}">
      <div class="pane-h"><span class="grip" draggable="true" title="Drag to reorder">${ICON.grip}</span><h2>${p.title}</h2>
        <span class="sum">${paneSummary(id, r)}</span>
        <button data-act="wide" aria-pressed="${w}" title="${w ? 'Half width' : 'Full width'}" aria-label="Toggle width">${ICON.wide}</button>
        <button data-act="fold" aria-expanded="${!c}" aria-label="${c ? 'Expand' : 'Collapse'} ${p.title}">${c ? ICON.unfold : ICON.fold}</button></div>
      <div class="pane-b">${c ? '' : paneBody(id, o, r)}</div></section>`;
  }).join('');
  try { window.scrollTo(0, scroll); } catch (e) {}
}
function renderOrigin(o) {
  const spots = [];
  if (o.draft && o.parent) { $('o-origin').hidden = true; $('o-origin').innerHTML = ''; return; }
  for (const p of lib) if (p.kind === 'seq') { const h = findOn(o, p); if (h) spots.push({ p, h }); }
  const el = $('o-origin');
  if (!spots.length) { el.hidden = true; el.innerHTML = ''; return; }
  el.hidden = false;
  el.innerHTML = spots.slice(0, 4).map(({ p, h }) => `<span>On <b>${esc(p.name)}</b> · <span class="mono num">${h.a + 1}–${h.b + 1}</span> · <span class="dir ${h.strand === '+' ? 'fwd' : 'rev'}">${h.strand === '+' ? '→ top strand' : '← bottom strand'}</span></span><button class="btn ghost" data-goto="${p.id}" data-a="${h.a}" data-b="${h.b}" data-s="${h.strand}">Show in sequence</button>`).join('');
}
function renderOligoView(sync) {
  const o = cur();
  if (sync) { $('name').value = o.name; $('seq').value = o.seq; $('role').value = o.role || ''; $('msg').hidden = true; }
  $('del').textContent = 'Delete'; $('del').dataset.armed = '';
  $('del').hidden = !!o.draft;
  const db = $('o-draft');
  if (o.draft) {
    const p = o.parent && byId(o.parent), at = p && findOn(o, p);
    const exact = at && o.start && at.a === o.start - 1 && at.b === o.end - 1 && at.strand === o.strand;
    db.hidden = false;
    db.innerHTML = `<div class="dtext"><b>Not saved yet.</b> ${p ? `Selection from ${esc(p.name)}, <span class="mono num">${o.start}–${o.end}</span> ${o.strand === '+' ? '→ top strand' : '← bottom strand'}.` : 'Scratch oligo.'} Save it once it looks good.</div>
      ${exact ? `<div class="nudge"><span class="hint">5′ end</span><button class="btn" data-nudge="5,-1" title="Trim one base at the 5′ end">−1</button><button class="btn" data-nudge="5,1" title="Add one base at the 5′ end">+1</button>
        <span class="hint">3′ end</span><button class="btn" data-nudge="3,-1" title="Trim one base at the 3′ end">−1</button><button class="btn" data-nudge="3,1" title="Add one base at the 3′ end">+1</button></div>` : ''}
      <div class="dact">${p ? '<button class="btn ghost" data-dr="back">Back to sequence</button>' : ''}<button class="btn ghost" data-dr="discard">Discard</button><button class="btn primary" data-dr="save">Save to library</button></div>`;
  } else { db.hidden = true; db.innerHTML = ''; }
  renderOrigin(o);
  const ban = $('o-banner');
  if (o.seq.length > LONG) { ban.hidden = false; ban.innerHTML = `<span>${o.seq.length.toLocaleString('en-US')} nt looks like a template rather than a primer. Open it as a sequence to see both strands and save primers from it.</span><button class="btn primary" id="to-seq">Open as sequence</button>`; }
  else { ban.hidden = true; ban.innerHTML = ''; }
  renderMap($('o-map'), o.seq, { maxBlocks: 10 });
  const r = analyze(o);
  lastRes = r;
  renderMetrics(r);
  renderPanes(o, r);
}

/* ================= shared ================= */
function renderAll(sync) {
  const o = cur(); if (!o) return;
  sel = o.id;
  const isSeq = o.kind === 'seq';
  $('seqview').hidden = !isSeq; $('oligoview').hidden = isSeq;
  if (isSeq) renderSeqView(sync); else renderOligoView(sync);
  renderLib();
}
let timer = null;
const schedule = (ms = 90) => { clearTimeout(timer); timer = setTimeout(() => { save(); renderAll(false); }, ms); };
function toast(t) { const el = $('toast'); el.textContent = t; el.hidden = false; clearTimeout(toast.t); toast.t = setTimeout(() => el.hidden = true, 1900); }
async function copy(text, done) {
  try { await navigator.clipboard.writeText(text); toast(done); }
  catch (e) { $('imp').hidden = false; $('imp-t').textContent = 'Copy this text'; $('imp-note').hidden = true; $('imp-go').hidden = true; const ta = $('imp-txt'); ta.value = text; ta.select(); }
}
function readSeqInput(raw, o, msgEl, nameEl) {
  let body = raw, header = null;
  if (/^LOCUS/m.test(raw) && /^ORIGIN/m.test(raw)) { header = (raw.match(/^LOCUS\s+(\S+)/m) || [])[1]; body = (raw.split(/^ORIGIN.*$/m)[1] || '').split('//')[0].replace(/\d/g, ''); }
  else { const fa = raw.match(/^>(.+)$/m); if (fa) header = fa[1].trim(); body = raw.split('\n').filter(l => !l.startsWith('>')).join(''); }
  const bad = body.toUpperCase().replace(/[\sACGTU0-9]/g, '');
  o.seq = S.clean(body); delete o.ex;
  if (header && (!o.name || /^New (oligo|sequence)$/.test(o.name))) { o.name = header; nameEl.value = header; }
  if (bad) { msgEl.hidden = false; msgEl.textContent = `Ignored ${bad.length} character${bad.length > 1 ? 's' : ''} (${[...new Set(bad)].slice(0, 8).join(' ')}). Only A, C, G, T and U are kept; N and other IUPAC codes are removed, which shifts positions after them.`; }
  else msgEl.hidden = true;
}
function parseText(txt, fallback) {
  const out = [];
  txt = txt.replace(/\r/g, '');
  if (/^LOCUS/m.test(txt) && /^ORIGIN/m.test(txt)) {
    txt.split(/^\/\/\s*$/m).forEach(rec => {
      if (!/^ORIGIN/m.test(rec)) return;
      const name = (rec.match(/^LOCUS\s+(\S+)/m) || [])[1] || fallback;
      const s = S.clean((rec.split(/^ORIGIN.*$/m)[1] || '').replace(/\d/g, ''));
      if (s) out.push({ name, seq: s });
    });
    return out;
  }
  if (/^>/m.test(txt)) {
    txt.split(/^>/m).filter(x => x.trim()).forEach(rec => { const [h, ...rest] = rec.split('\n'); const s = S.clean(rest.join('')); if (s) out.push({ name: h.trim() || fallback, seq: s }); });
    return out;
  }
  const lines = txt.split('\n').map(l => l.trim()).filter(Boolean);
  const listy = lines.some(l => /[\t,;]/.test(l));
  const wrapped = !listy && lines.length > 1 && lines.slice(0, -1).every(l => l.length === lines[0].length) && lines[0].length >= 50;
  if (lines.length === 1 || wrapped) { const s = S.clean(lines.join('')); if (s) out.push({ name: fallback, seq: s }); return out; }
  lines.forEach((l, k) => { const parts = l.split(/[\t,;]+/).map(x => x.trim()).filter(Boolean); const s = S.clean(parts[parts.length - 1]); if (s) out.push({ name: parts.length > 1 ? parts[0] : `${fallback} ${k + 1}`, seq: s }); });
  return out;
}
function addParsed(recs) {
  const added = recs.map(r => {
    const kind = r.seq.length > LONG ? 'seq' : 'oligo';
    const o = { id: uid(), kind, name: r.name, seq: r.seq };
    if (kind === 'oligo') o.role = /(^|[\s_-])(f|fw|fwd|forward)$/i.test(r.name) ? 'F' : /(^|[\s_-])(r|rv|rev|reverse)$/i.test(r.name) ? 'R' : '';
    return o;
  });
  if (!added.length) { toast('No DNA sequence found'); return false; }
  lib = added.filter(o => o.kind === 'seq').concat(lib, added.filter(o => o.kind === 'oligo'));
  sel = added[0].id; save(); renderAll(true);
  toast(`Added ${added.length} item${added.length > 1 ? 's' : ''}`);
  return true;
}

/*BLAST*/

/* ---------- events: sequence view ---------- */
$('s-name').addEventListener('input', e => { cur().name = e.target.value; clearTimeout(timer); timer = setTimeout(() => { save(); renderLib(); renderSelBar(true); mapKey = ''; renderSeqMap(); }, 250); });
$('s-seq').addEventListener('input', e => { readSeqInput(e.target.value, cur(), $('s-msg'), $('s-name')); if (SV.hits.length || $('s-find').value) SV.hits = []; schedule(cur().seq.length > 20000 ? 400 : 160); });
$('s-editbtn').onclick = () => { SV.editing = !SV.editing; renderSeqView(false); if (SV.editing) $('s-seq').focus(); };
$('s-find').addEventListener('input', computeHits);
$('s-find').addEventListener('keydown', e => { if (e.key === 'Enter' && SV.hits.length) { e.preventDefault(); $('s-next').click(); } });
$('s-next').onclick = () => { if (!SV.hits.length) return; SV.hitIdx = (SV.hitIdx + 1) % SV.hits.length; const h = SV.hits[SV.hitIdx]; setSel(h.a, h.b, h.strand); renderSelBar(true); scrollToSel(); };
$('s-plus').onclick = () => { setPick('+'); if (SV.sel) { setSel(SV.sel.a, SV.sel.b, '+'); renderSelBar(true); } };
$('s-minus').onclick = () => { setPick('-'); if (SV.sel) { setSel(SV.sel.a, SV.sel.b, '-'); renderSelBar(true); } };
$('s-selrange').onclick = () => {
  const n = cur().seq.length, a = parseInt($('s-from').value, 10), b = parseInt($('s-to').value, 10);
  if (!n || !isFinite(a) || !isFinite(b)) { toast('Enter both positions'); return; }
  if (a < 1 || b < 1 || a > n || b > n) { toast(`Positions must be between 1 and ${n}`); return; }
  setSel(a - 1, b - 1, SV.pick); renderSelBar(true); scrollToSel();
};
['s-from', 's-to'].forEach(id => $(id).addEventListener('keydown', e => { if (e.key === 'Enter') $('s-selrange').click(); }));
$('s-copy').onclick = () => copy(cur().seq, 'Top strand copied (5′→3′)');
$('s-copyrc').onclick = () => copy(S.revcomp(cur().seq), 'Bottom strand copied (5′→3′)');
$('s-flip').onclick = () => { const o = cur(); o.seq = S.revcomp(o.seq); delete o.ex; $('s-seq').value = o.seq; SV.sel = null; SV.hits = []; save(); renderAll(false); computeHits(); toast('Strands flipped; positions now count from the other end'); };
$('s-del').onclick = e => {
  const b = e.currentTarget;
  if (!b.dataset.armed) { b.dataset.armed = '1'; b.textContent = 'Confirm delete'; return; }
  const id = sel; lib = lib.filter(o => o.id !== id);
  lib.forEach(o => { if (o.parent === id) delete o.parent; });
  if (!lib.length) lib.push({ id: uid(), kind: 'seq', name: 'New sequence', seq: '' });
  sel = lib[0].id; save(); renderAll(true);
};

// map selection: drag with mouse/pen, two taps on touch
const map = $('s-map');
let lastPointer = 'mouse';
const baseAt = (x, y) => { const el = document.elementFromPoint(x, y); return el && el.closest ? el.closest('#s-map .st i') : null; };
map.addEventListener('pointerdown', e => {
  lastPointer = e.pointerType;
  const el = e.target.closest('.st i'); if (!el || e.pointerType === 'touch' || e.button !== 0) return;
  const i = +el.dataset.i, strand = el.parentElement.classList.contains('dn') ? '-' : '+';
  if (e.shiftKey && SV.sel) { setSel(SV.anchor != null ? SV.anchor : SV.sel.a, i, SV.sel.strand); renderSelBar(true); return; }
  e.preventDefault();
  SV.anchor = i; SV.dragging = true; SV.touchPending = false;
  setSel(i, i, strand);
  try { map.setPointerCapture(e.pointerId); } catch (x) {}
});
map.addEventListener('pointermove', e => {
  if (!SV.dragging) return;
  const el = baseAt(e.clientX, e.clientY);
  if (el) setSel(SV.anchor, +el.dataset.i, SV.sel.strand);
});
const endDrag = () => { if (!SV.dragging) return; SV.dragging = false; renderSelBar(true); };
map.addEventListener('pointerup', endDrag);
map.addEventListener('pointercancel', endDrag);
map.addEventListener('click', e => {
  const bar = e.target.closest('.bar');
  if (bar) { sel = bar.dataset.oid; save(); renderAll(true); window.scrollTo(0, 0); return; }
  if (lastPointer !== 'touch') return;
  const el = e.target.closest('.st i'); if (!el) return;
  const i = +el.dataset.i, strand = el.parentElement.classList.contains('dn') ? '-' : '+';
  if (SV.touchPending && SV.sel) { SV.touchPending = false; setSel(SV.anchor, i, SV.sel.strand); renderSelBar(true); }
  else { SV.anchor = i; SV.touchPending = true; setSel(i, i, strand); }
});
$('selbar').addEventListener('input', e => { if (e.target.id === 'sel-name') SV.nameDirty = true; });
$('selbar').addEventListener('keydown', e => { if (e.target.id === 'sel-name' && e.key === 'Enter') saveSel(false); });
$('selbar').addEventListener('click', e => {
  const b = e.target.closest('[data-sb]'); if (!b) return;
  const act = b.dataset.sb;
  if (act === 'clear') clearSel();
  else if (act === 'flip') { setSel(SV.sel.a, SV.sel.b, SV.sel.strand === '+' ? '-' : '+'); renderSelBar(true); }
  else if (act === 'copy') copy(selSeq(), 'Copied 5′→3′');
  else if (act === 'save') saveSel(false);
  else if (act === 'analyze') openDraftFromSel();
});
$('frag-body').addEventListener('click', e => {
  const bp = e.target.closest('[data-blpair]');
  if (bp) {
    const [f, r] = bp.dataset.blpair.split(','); sel = f; partner = r; closed.delete('blast'); save(); renderAll(true);
    const pane = document.querySelector('[data-pane=blast]'); if (pane) pane.scrollIntoView({ block: 'start', behavior: reduceMotion() ? 'auto' : 'smooth' });
    return;
  }
  const b = e.target.closest('[data-fa]'); if (!b) return;
  const x = byId(b.dataset.oid); if (!x) return;
  if (b.dataset.fa === 'open') { sel = x.id; save(); renderAll(true); window.scrollTo(0, 0); return; }
  const h = findOn(x, cur()); if (h) { setSel(h.a, h.b, h.strand); renderSelBar(true); scrollToSel(); }
});
let roT = null;
if (window.ResizeObserver) new ResizeObserver(() => { clearTimeout(roT); roT = setTimeout(() => { if (!$('seqview').hidden) renderSeqMap(false); else if (cur()) renderMap($('o-map'), cur().seq, { maxBlocks: 10 }); }, 120); }).observe($('s-map'));
window.addEventListener('resize', () => { clearTimeout(roT); roT = setTimeout(() => { if (!$('seqview').hidden) renderSeqMap(false); else if (cur()) renderMap($('o-map'), cur().seq, { maxBlocks: 10 }); }, 150); });

/* ---------- events: oligo bench ---------- */
$('seq').addEventListener('input', e => { readSeqInput(e.target.value, cur(), $('msg'), $('name')); schedule(); });
$('name').addEventListener('input', e => { cur().name = e.target.value; schedule(); });
$('role').addEventListener('change', e => { cur().role = e.target.value; schedule(); });
function setSeq(fn) { const o = cur(); o.seq = fn(o.seq); delete o.ex; $('seq').value = o.seq; save(); renderAll(false); }
$('rc').onclick = () => setSeq(S.revcomp);
$('comp').onclick = () => setSeq(s => s.split('').map(b => COMP[b]).join(''));
$('rev').onclick = () => setSeq(s => s.split('').reverse().join(''));
$('copyseq').onclick = () => copy(cur().seq, 'Sequence copied');
$('report').onclick = () => {
  const o = cur(), r = lastRes; if (!r || !r.th) return;
  const L = [`${o.name} (${o.seq.length} nt)`, `5'-${o.seq}-3'`, '',
    `Conditions: Na+ ${cond.na} mM, Mg2+ ${cond.mg} mM, dNTP ${cond.dntp} mM, oligo ${cond.oligo} nM; ${cond.method === 'owczarzy' ? 'Owczarzy' : 'SantaLucia'} salt correction`,
    `Tm ${fx(r.th.tm, 1)} °C | GC ${fx(r.ph.gc, 1)} % | ΔG(${cond.temp} °C) ${fx(r.th.dG, 2)} kcal/mol | ΔH ${fx(r.th.dH, 1)} kcal/mol | ΔS ${fx(r.th.dS, 1)} cal/K·mol`,
    `MW ${r.ph.mw.toFixed(1)} g/mol | ε260 ${Math.round(r.ph.eps)} L/mol·cm | ${fx(r.ph.nmolPerOD, 2)} nmol/OD | ${fx(r.ph.ugPerOD, 2)} µg/OD`];
  if (r.hp) L.push(`Strongest hairpin: ${r.hp[0] ? fx(r.hp[0].dG, 2) + ' kcal/mol' : 'none'}`);
  if (r.sd) L.push(`Strongest self-dimer: ${r.sd[0] ? fx(r.sd[0].dG, 2) + ' kcal/mol' : 'none'}`);
  L.push('', 'Checks:'); checks(o, r).forEach(c => L.push(`  [${c.st === 'ok' ? 'pass' : c.st === 'warn' ? 'check' : 'fail'}] ${c.name}: ${c.val}`));
  copy(L.join('\n'), 'Report copied');
};
$('del').onclick = e => {
  const b = e.currentTarget;
  if (!b.dataset.armed) { b.dataset.armed = '1'; b.textContent = 'Confirm delete'; return; }
  const gone = cur(); lib = lib.filter(o => o.id !== gone.id);
  if (!lib.length) lib.push({ id: uid(), kind: 'oligo', name: 'New oligo', seq: '', role: 'F' });
  sel = (gone.parent && byId(gone.parent)) ? gone.parent : lib[0].id; save(); renderAll(true);
};
$('oligoview').addEventListener('click', e => {
  if (e.target.id === 'to-seq') {
    let o = cur();
    if (o.draft) { o = Object.assign({}, o, { id: uid() }); delete o.draft; delete o.parent; delete o.start; delete o.end; delete o.strand; lib.push(o); draft = null; }
    o.kind = 'seq'; delete o.role; sel = o.id; save(); renderAll(true); return;
  }
  const nb = e.target.closest('[data-nudge]');
  if (nb) { const [end, d] = nb.dataset.nudge.split(',').map(Number); nudge(end, d); return; }
  const dr = e.target.closest('[data-dr]');
  if (dr) {
    const act = dr.dataset.dr;
    if (act === 'save') saveDraft();
    else if (act === 'discard') { const back = draft && draft.parent; draft = null; sel = back && byId(back) ? back : lib[0].id; save(); renderAll(true); }
    else if (act === 'back') { const p = byId(draft.parent), h = p && findOn(draft, p); if (h) SV.pending = { a: h.a, b: h.b, strand: h.strand }; sel = draft.parent; save(); renderAll(true); }
    return;
  }
  const g = e.target.closest('[data-goto]');
  if (g) { SV.pending = { a: +g.dataset.a, b: +g.dataset.b, strand: g.dataset.s }; sel = g.dataset.goto; save(); renderAll(true); }
});
$('panes').addEventListener('click', e => {
  const b = e.target.closest('button[data-act]'); if (!b) return;
  const id = b.closest('.pane').dataset.pane;
  if (b.dataset.act === 'fold') closed.has(id) ? closed.delete(id) : closed.add(id);
  else wide.has(id) ? wide.delete(id) : wide.add(id);
  save(); renderPanes(cur(), lastRes);
});
$('panes').addEventListener('change', e => { if (e.target.id === 'partner') { partner = e.target.value; save(); renderPanes(cur(), lastRes); } });
let dragId = null;
$('panes').addEventListener('dragstart', e => { const p = e.target.closest('.pane'); if (!p) return; dragId = p.dataset.pane; e.dataTransfer.effectAllowed = 'move'; try { e.dataTransfer.setData('text/plain', dragId); e.dataTransfer.setDragImage(p, 20, 20); } catch (x) {} });
$('panes').addEventListener('dragover', e => { const p = e.target.closest('.pane'); if (!p || !dragId) return; e.preventDefault(); document.querySelectorAll('.pane.drag-over').forEach(x => x !== p && x.classList.remove('drag-over')); p.classList.add('drag-over'); });
$('panes').addEventListener('drop', e => {
  const p = e.target.closest('.pane'); if (!p || !dragId) return; e.preventDefault();
  const to = p.dataset.pane; if (to !== dragId) { order = order.filter(x => x !== dragId); order.splice(order.indexOf(to), 0, dragId); save(); }
  dragId = null; renderPanes(cur(), lastRes);
});
$('panes').addEventListener('dragend', () => { dragId = null; document.querySelectorAll('.pane.drag-over').forEach(x => x.classList.remove('drag-over')); });

/* ---------- events: library & import ---------- */
$('lib').addEventListener('click', e => { const it = e.target.closest('.item'); if (!it) return; sel = it.dataset.id; save(); renderAll(true); });
$('new').onclick = () => { draft = { id: DRAFT, kind: 'oligo', draft: true, name: 'New oligo', seq: '', role: 'F' }; sel = DRAFT; renderAll(true); $('seq').focus(); };
$('new-seq').onclick = () => { const o = { id: uid(), kind: 'seq', name: 'New sequence', seq: '' }; lib.unshift(o); sel = o.id; save(); renderAll(true); $('s-seq').focus(); };
$('export').onclick = () => copy(lib.map(o => `>${o.name || 'Untitled'}\n${(o.seq.match(/.{1,70}/g) || []).join('\n')}`).join('\n'), 'Library copied as FASTA');
$('import').onclick = () => { $('imp-t').textContent = 'Paste sequences'; $('imp-note').hidden = false; $('imp-go').hidden = false; $('imp').hidden = false; $('imp-txt').value = ''; $('imp-txt').focus(); };
$('imp-x').onclick = () => $('imp').hidden = true;
$('imp').addEventListener('click', e => { if (e.target.id === 'imp') $('imp').hidden = true; });
$('imp-go').onclick = () => { const t = $('imp-txt').value.trim(); if (!t) return; if (addParsed(parseText(t, 'Pasted'))) $('imp').hidden = true; };
$('open').onclick = () => $('file').click();
$('file').addEventListener('change', async e => {
  const recs = [];
  for (const f of e.target.files) { try { const t = await f.text(); recs.push(...parseText(t, f.name.replace(/\.[^.]+$/, ''))); } catch (x) {} }
  e.target.value = '';
  addParsed(recs);
});
document.addEventListener('dragover', e => { if (e.dataTransfer && [...e.dataTransfer.types].includes('Files')) e.preventDefault(); });
document.addEventListener('drop', async e => {
  if (!e.dataTransfer || !e.dataTransfer.files.length) return;
  e.preventDefault(); const recs = [];
  for (const f of e.dataTransfer.files) { try { recs.push(...parseText(await f.text(), f.name.replace(/\.[^.]+$/, ''))); } catch (x) {} }
  addParsed(recs);
});
document.addEventListener('keydown', e => {
  if (e.key !== 'Escape') return;
  if (!$('imp').hidden) { $('imp').hidden = true; return; }
  if (SV.sel && !$('seqview').hidden) clearSel();
});

/* ---------- conditions ---------- */
const CF = { na: 'c-na', mg: 'c-mg', dntp: 'c-dntp', oligo: 'c-oligo', temp: 'c-temp' };
Object.entries(CF).forEach(([k, id]) => {
  $(id).value = cond[k];
  $(id).addEventListener('input', e => { const v = parseFloat(e.target.value); if (isFinite(v) && (k === 'temp' || v >= 0)) { cond[k] = k === 'oligo' ? Math.max(v, 0.001) : v; schedule(); } });
});
const setMethod = m => { cond.method = m; $('m-sl').setAttribute('aria-pressed', m === 'santalucia'); $('m-ow').setAttribute('aria-pressed', m === 'owczarzy'); };
setMethod(cond.method);
$('m-sl').onclick = () => { setMethod('santalucia'); schedule(); };
$('m-ow').onclick = () => { setMethod('owczarzy'); schedule(); };

wireBlast();
if (!lib.some(o => o.id === sel)) sel = lib[0].id;
renderAll(true);
save();
if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { mapKey = ''; renderAll(false); });
})();
