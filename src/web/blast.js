/* ================= NCBI BLAST =================
   Three ways to compare a pair of oligos with NCBI:
   - bind:  BLAST both oligos against an NCBI nucleotide database. In the desktop app the search runs
            through NCBI's Common URL API and the hits are shown here; a page in a browser cannot reach
            NCBI directly, so there it opens NCBI's BLAST page with the query filled in.
   - pair:  Primer-BLAST specificity check for the pair. Primer-BLAST has no public API, so this opens
            NCBI's Primer-BLAST form with both primers filled in.
   - align: align the two oligos to each other on both strands (computed here with blastn's default
            scoring), with a button that opens NCBI's "Align two sequences" page. */
const IS_DESKTOP = !!(window.__TAURI__ && window.__TAURI__.core && window.__TAURI__.core.invoke);
const BL = Object.assign({ mode: 'bind', db: 'core_nt', org: '', pbdb: 'refseq_mrna', email: '' }, store.get('blast', {}));
const saveBL = () => store.set('blast', { mode: BL.mode, db: BL.db, org: BL.org, pbdb: BL.pbdb, email: BL.email });
const BJ = new Map();      // job key -> job
const openHits = new Set(); // expanded alignment rows
const BLAST_DBS = [['core_nt', 'Core nucleotide (core_nt)'], ['refseq_rna', 'RefSeq RNA (refseq_rna)'], ['refseq_genomic', 'RefSeq genomic (refseq_genomic)'], ['nt', 'Nucleotide collection (nt)']];
const PB_DBS = [['refseq_mrna', 'RefSeq mRNA'], ['refseq_representative_genomes', 'RefSeq representative genomes'], ['nt', 'Nucleotide collection (nt)']];
const enc = encodeURIComponent;
const fasta = qs => qs.map(q => `>${q.label}\n${q.seq}`).join('\n');
const extLink = (url, label, cls = 'btn') => `<a class="${cls}" href="${esc(url)}" target="_blank" rel="noopener" data-ext>${label}</a>`;

function pickPartner(o) {
  const others = lib.filter(x => x.kind === 'oligo' && x.id !== o.id && x.seq.length >= 2);
  let p = others.find(x => x.id === partner);
  if (!p) p = others.find(x => o.role === 'F' ? x.role === 'R' : o.role === 'R' ? x.role === 'F' : false) || others[0];
  return { others, p };
}
function partnerSelect(id, others, p) {
  return `<select id="${id}">${others.map(x => `<option value="${x.id}" ${x.id === p.id ? 'selected' : ''}>${esc(x.name || 'Untitled')} (${x.seq.length} nt)</option>`).join('')}</select>`;
}
// forward primer first, if the roles say so
function orderPair(o, p) { return (o.role === 'R' && p && p.role === 'F') ? [p, o] : [o, p]; }

/* ---------- NCBI URLs (open in the browser) ---------- */
function blastWebURL(qs) {
  return `https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome&DATABASE=${enc(BL.db)}&QUERY=${enc(fasta(qs))}` +
    (BL.org ? `&EQ_MENU=${enc(BL.org)}` : '');
}
function primerBlastURL(f, r) {
  return `https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?PRIMER_LEFT_INPUT=${enc(f.seq)}&PRIMER_RIGHT_INPUT=${enc(r.seq)}&PRIMER_SPECIFICITY_DATABASE=${enc(BL.pbdb)}` +
    (BL.org ? `&ORGANISM=${enc(BL.org)}` : '');
}
function align2URL(a, b) {
  return `https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&BLAST_SPEC=blast2seq&LINK_LOC=align2seq&QUERY=${enc(fasta([a]))}&SUBJECTS=${enc(fasta([b]))}`;
}
function openExternal(url) {
  if (IS_DESKTOP) return window.__TAURI__.core.invoke('plugin:opener|open_url', { url }).catch(() => toast('Could not open the browser'));
}

/* ---------- in-app BLAST (desktop) ---------- */
const sleep = ms => new Promise(r => setTimeout(r, ms));
async function ncbi(params) {
  const p = params.concat([['TOOL', 'stackwise']]);
  if (BL.email && /.+@.+\..+/.test(BL.email)) p.push(['EMAIL', BL.email.trim()]);
  return window.__TAURI__.core.invoke('ncbi_blast', { params: p });
}
function ncbiError(html) {
  const m = html.match(/Message ID#\d+ Error:\s*([^<\n]+)/i) || html.match(/class="error"[^>]*>\s*([^<]+)/i);
  return m ? m[1].trim() : 'NCBI did not accept the search.';
}
function jobKey(qs) { return [BL.db, BL.org.trim().toLowerCase(), ...qs.map(q => q.seq)].join('|'); }

async function runBlast(qs) {
  const key = jobKey(qs);
  const job = { key, qs: qs.map(q => ({ label: q.label, name: q.name, seq: q.seq })), db: BL.db, org: BL.org.trim(), status: 'submitting', started: Date.now(), rid: null, rtoe: null, results: null, error: null, cancelled: false };
  BJ.set(key, job); refreshBlast(true);
  try {
    const put = [['CMD', 'Put'], ['PROGRAM', 'blastn'], ['DATABASE', job.db], ['QUERY', fasta(job.qs)],
      ['WORD_SIZE', '7'], ['EXPECT', '1000'], ['FILTER', 'F'], ['NUCL_REWARD', '1'], ['NUCL_PENALTY', '-3'], ['GAPCOSTS', '5 2'],
      ['SHORT_QUERY_ADJUST', 'true'], ['HITLIST_SIZE', '50'], ['FORMAT_TYPE', 'JSON2_S']];
    if (job.org) put.push(['ENTREZ_QUERY', `"${job.org.replace(/"/g, '')}"[organism]`]);
    const html = await ncbi(put);
    const rid = (html.match(/RID\s*=\s*([A-Z0-9-]+)/) || [])[1];
    if (!rid) throw new Error(ncbiError(html));
    job.rid = rid; job.rtoe = +((html.match(/RTOE\s*=\s*(\d+)/) || [])[1] || 20);
    job.status = 'waiting'; refreshBlast(true);
    await sleep(Math.max(10, Math.min(job.rtoe, 60)) * 1000);
    for (;;) {
      if (job.cancelled) return;
      if (Date.now() - job.started > 20 * 60 * 1000) throw new Error('NCBI took longer than 20 minutes. Try again later or open the search on NCBI.');
      const info = await ncbi([['CMD', 'Get'], ['FORMAT_OBJECT', 'SearchInfo'], ['RID', rid]]);
      const st = ((info.match(/Status=(\w+)/) || [])[1] || 'UNKNOWN').toUpperCase();
      if (st === 'READY') {
        if (/ThereAreHits=no/i.test(info)) { job.results = job.qs.map(q => ({ q, len: q.seq.length, hits: [] })); break; }
        const out = await ncbi([['CMD', 'Get'], ['RID', rid], ['FORMAT_TYPE', 'JSON2_S'], ['DESCRIPTIONS', '50'], ['ALIGNMENTS', '50']]);
        job.results = parseBlast(out, job.qs); break;
      }
      if (st === 'FAILED') throw new Error('NCBI reported that the search failed.');
      if (st === 'UNKNOWN') throw new Error('NCBI no longer knows this search (it may have expired).');
      await sleep(15000);
    }
    if (job.cancelled) return;
    job.status = 'done'; job.finished = Date.now();
  } catch (e) {
    if (job.cancelled) return;
    job.status = 'error'; job.error = String(e && e.message || e).replace(/^Error:\s*/, '');
  }
  refreshBlast(true);
}

function parseBlast(text, qs) {
  const i = text.indexOf('{'); if (i < 0) throw new Error('NCBI sent results in an unexpected format.');
  const j = JSON.parse(text.slice(i));
  let reps = j.BlastOutput2; if (!Array.isArray(reps)) reps = [reps];
  return reps.filter(Boolean).map((rep, k) => {
    const s = (rep.report && rep.report.results && rep.report.results.search) || {};
    const q = qs.find(x => s.query_title && s.query_title.split(/\s/)[0] === x.label) || qs[k];
    const len = s.query_len || q.seq.length;
    const hits = (s.hits || []).map(h => {
      const d = (h.description || [])[0] || {};
      const hsps = (h.hsps || []).map(x => {
        const lo = Math.min(x.hit_from, x.hit_to), hi = Math.max(x.hit_from, x.hit_to);
        const minus = /minus/i.test(x.hit_strand || '') || x.hit_from > x.hit_to;
        const mid = x.midline || '';
        const end3 = x.query_to === len && /\|\|\|$/.test(mid);
        return Object.assign({}, x, { lo, hi, minus, end3, mism: (x.align_len || 0) - (x.identity || 0) - (x.gaps || 0) });
      });
      return { acc: d.accession || d.id || '', title: d.title || '', sci: d.sciname || '', len: h.len, hsps };
    });
    return { q, len, hits, msg: s.message };
  });
}

// targets hit by both oligos, facing each other, close enough to make a product
function sharedProducts(results) {
  if (!results || results.length < 2) return [];
  const [A, B] = results, out = [];
  const byAcc = new Map(B.hits.map(h => [h.acc, h]));
  for (const ha of A.hits) {
    const hb = byAcc.get(ha.acc); if (!hb) continue;
    let best = null;
    for (const x of ha.hsps) for (const y of hb.hsps) {
      if (x.minus === y.minus) continue;
      const f = x.minus ? y : x, r = x.minus ? x : y;
      const size = r.hi - f.lo + 1;
      if (f.lo < r.hi && size > 0 && size <= 10000) {
        const score = (f.end3 && r.end3 ? 2 : 0) + (f.mism + r.mism === 0 ? 1 : 0);
        if (!best || score > best.score || (score === best.score && size < best.size)) best = { size, f, r, score, fIsA: f === x };
      }
    }
    if (best) out.push(Object.assign({ hit: ha }, best));
  }
  return out.sort((u, v) => v.score - u.score || u.size - v.size);
}

/* ---------- local two-sequence alignment (blastn scoring) ---------- */
function localAlign(a, b) {
  const MA = 2, MI = -3, GO = 5, GE = 2, n = a.length, m = b.length, W = m + 1, NEG = -1e9;
  const H = new Float64Array((n + 1) * W), E = new Float64Array((n + 1) * W).fill(NEG), F = new Float64Array((n + 1) * W).fill(NEG);
  const tH = new Uint8Array((n + 1) * W), tE = new Uint8Array((n + 1) * W), tF = new Uint8Array((n + 1) * W);
  let best = 0, bi = 0, bj = 0;
  for (let i = 1; i <= n; i++) for (let j = 1; j <= m; j++) {
    const k = i * W + j;
    const eo = H[k - 1] - GO - GE, ee = E[k - 1] - GE; E[k] = Math.max(eo, ee); tE[k] = eo >= ee ? 1 : 2;
    const fo = H[k - W] - GO - GE, fe = F[k - W] - GE; F[k] = Math.max(fo, fe); tF[k] = fo >= fe ? 1 : 2;
    const d = H[k - W - 1] + (a[i - 1] === b[j - 1] ? MA : MI);
    let h = 0, t = 0;
    if (d > h) { h = d; t = 1; } if (E[k] > h) { h = E[k]; t = 2; } if (F[k] > h) { h = F[k]; t = 3; }
    H[k] = h; tH[k] = t;
    if (h > best) { best = h; bi = i; bj = j; }
  }
  if (best <= 0) return null;
  let i = bi, j = bj, st = 'H', top = '', bot = '';
  for (;;) {
    const k = i * W + j;
    if (st === 'H') { const t = tH[k]; if (t === 0) break; if (t === 1) { top = a[i - 1] + top; bot = b[j - 1] + bot; i--; j--; } else st = t === 2 ? 'E' : 'F'; }
    else if (st === 'E') { top = '-' + top; bot = b[j - 1] + bot; const t = tE[k]; j--; if (t === 1) st = 'H'; }
    else { top = a[i - 1] + top; bot = '-' + bot; const t = tF[k]; i--; if (t === 1) st = 'H'; }
    if (i <= 0 && j <= 0) break;
  }
  let mid = '', ident = 0, gaps = 0;
  for (let x = 0; x < top.length; x++) { const s = top[x] === bot[x]; mid += s ? '|' : ' '; ident += s; gaps += (top[x] === '-' || bot[x] === '-'); }
  return { score: best, aFrom: i + 1, aTo: bi, bFrom: j + 1, bTo: bj, top, mid, bot, ident, gaps, len: top.length };
}

/* ---------- rendering ---------- */
function alnPre(l1, s1, from1, to1, mid, l2, s2, from2, to2) {
  const w = Math.max(l1.length, l2.length), pw = Math.max(String(from1).length, String(from2).length);
  const pad = (s, n) => s + ' '.repeat(Math.max(0, n - s.length)), lp = (s, n) => ' '.repeat(Math.max(0, n - String(s).length)) + s;
  const col = (s, other) => s.split('').map((c, k) => c === '-' ? '<span class="lab">-</span>' : `<span class="b${c}${other && other[k] !== c && other[k] !== '-' ? ' hot' : ''}">${c}</span>`).join('');
  return `<pre class="duplex"><span class="lab">${pad(l1, w)} ${lp(from1, pw)} </span>${col(s1, s2)}<span class="lab"> ${to1}</span>\n<span class="lab">${' '.repeat(w + pw + 2)}</span><span class="pr">${mid}</span>\n<span class="lab">${pad(l2, w)} ${lp(from2, pw)} </span>${col(s2, s1)}<span class="lab"> ${to2}</span></pre>`;
}
const fmtE = e => e == null ? '—' : e === 0 ? '0' : e < 0.01 ? e.toExponential(1).replace('e-', 'e−') : e < 100 ? e.toFixed(2) : Math.round(e).toLocaleString('en-US');
const mmss = ms => { const s = Math.round(ms / 1000); return `${Math.floor(s / 60)}:${String(s % 60).padStart(2, '0')}`; };
const hitName = h => (h.sci && !h.title.startsWith(h.sci) ? h.sci + ' · ' : '') + h.title;
const nuccore = acc => `https://www.ncbi.nlm.nih.gov/nuccore/${enc(acc)}`;

function hitsTable(res, jk) {
  if (!res.hits.length) return `<p class="empty">No hits${res.msg ? ' (' + esc(res.msg) + ')' : ''}. With an organism filter, check its spelling.</p>`;
  const perfect = res.hits.filter(h => h.hsps[0] && h.hsps[0].identity === res.len && h.hsps[0].align_len === res.len).length;
  const rows = res.hits.map((h, n) => {
    const x = h.hsps[0]; if (!x) return '';
    const id = `${jk}#${res.q.label}#${n}`, open = openHits.has(id);
    const cover = (x.query_to - x.query_from + 1) / res.len * 100;
    const full = x.identity === res.len && x.align_len === res.len;
    return `<tr data-hit="${esc(id)}" class="hitrow${open ? ' open' : ''}"><td class="nm" title="${esc(h.title)}">${esc(hitName(h))}</td>
      <td>${extLink(nuccore(h.acc), esc(h.acc), 'acc')}</td>
      <td class="n">${x.identity}/${res.len}</td><td class="n">${fx(cover, 0)}%</td>
      <td>${full ? pill('neutral', 'exact match') : x.end3 ? pill('warn', '3′ matches') : pill('ok', '3′ free')}</td>
      <td><span class="dir ${x.minus ? 'rev' : 'fwd'}">${x.minus ? '← minus' : '→ plus'}</span></td>
      <td class="n">${x.lo.toLocaleString('en-US')}–${x.hi.toLocaleString('en-US')}</td><td class="n">${fmtE(x.evalue)}</td></tr>` +
      (open ? `<tr class="alnrow"><td colspan="8">${h.hsps.slice(0, 3).map(y => alnPre(res.q.name, y.qseq, y.query_from, y.query_to, y.midline, h.acc, y.hseq, y.hit_from, y.hit_to)).join('')}${h.hsps.length > 3 ? `<p class="hint">${h.hsps.length - 3} more alignments on this sequence.</p>` : ''}</td></tr>` : '');
  }).join('');
  return `<p class="hint" style="margin:0">${res.hits.length} sequence${res.hits.length === 1 ? '' : 's'} hit · ${perfect} with a full-length perfect match. Click a row to see the alignment.</p>
    <div class="tablewrap"><table class="frag"><thead><tr><th>Sequence</th><th>Accession</th><th class="n">Identity</th><th class="n">Cover</th><th>3′ end</th><th>Strand</th><th class="n">Position</th><th class="n">E-value</th></tr></thead><tbody>${rows}</tbody></table></div>`;
}

function blastOutHTML(qs) {
  const job = BJ.get(jobKey(qs));
  if (!job) return '';
  if (job.status === 'submitting') return `<div class="bl-status">${pill('neutral', 'Sending')}<span>Sending the search to NCBI…</span></div>`;
  if (job.status === 'waiting') return `<div class="bl-status"><span class="spin" aria-hidden="true"></span><span>Running at NCBI · request <b class="mono">${esc(job.rid)}</b> · ${mmss(Date.now() - job.started)} elapsed (NCBI estimated ${job.rtoe} s). Checked every 15 s.</span><button class="btn ghost" data-bl="cancel">Stop waiting</button></div>`;
  if (job.status === 'error') return `<div class="bl-status err">${pill('bad', 'Failed')}<span>${esc(job.error)}</span></div>`;
  const res = job.results || [];
  const prods = sharedProducts(res);
  const jk = job.rid || job.key;
  let html = `<div class="bl-status">${pill('ok', 'Done')}<span>Searched ${esc(BLAST_DBS.find(d => d[0] === job.db)?.[1] || job.db)}${job.org ? ' · ' + esc(job.org) : ''} in ${mmss(job.finished - job.started)} · request <b class="mono">${esc(job.rid)}</b></span>
    ${extLink(`https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=${enc(job.rid)}`, 'Open results on NCBI', 'btn ghost')}</div>`;
  if (res.length >= 2) {
    html += `<div class="bl-q"><h3>Targets both oligos hit, facing each other</h3>` + (prods.length
      ? `<div class="tablewrap"><table class="frag"><thead><tr><th>Sequence</th><th>Accession</th><th class="n">Product</th><th>Forward</th><th>Reverse</th><th>Both 3′ ends match</th></tr></thead><tbody>${prods.slice(0, 25).map(p => `<tr><td class="nm" title="${esc(p.hit.title)}">${esc(hitName(p.hit))}</td><td>${extLink(nuccore(p.hit.acc), esc(p.hit.acc), 'acc')}</td><td class="n">${p.size.toLocaleString('en-US')} bp</td><td>${esc(p.fIsA ? res[0].q.name : res[1].q.name)} <span class="hint">${p.f.mism + p.f.gaps} mm</span></td><td>${esc(p.fIsA ? res[1].q.name : res[0].q.name)} <span class="hint">${p.r.mism + p.r.gaps} mm</span></td><td>${p.f.end3 && p.r.end3 ? pill('warn', 'yes') : pill('ok', 'no')}</td></tr>`).join('')}</tbody></table></div>
         <p class="hint" style="margin:0">Estimated from the top ${50} hits per oligo, products up to 10 kb. Primer-BLAST gives the full specificity check.</p>`
      : '<p class="empty">No sequence among the top hits is bound by both oligos in a product-forming orientation.</p>') + '</div>';
  }
  for (const r of res) html += `<div class="bl-q"><h3>${esc(r.q.name)} <span class="hint mono">${r.len} nt</span></h3>${hitsTable(r, jk)}</div>`;
  return html;
}

function alignHTML(a, b) {
  const A = a.seq.slice(0, 2000), B = b.seq.slice(0, 2000), m = B.length;
  const plus = localAlign(A, B), minus = localAlign(A, S.revcomp(B));
  const block = (r, strand) => {
    if (!r || r.score < 10) return `<div class="struct ok"><div class="struct-h"><span><b>${strand === '+' ? 'Plus / Plus' : 'Plus / Minus'}</b></span><span>No similarity worth reporting${r ? ` (best local score ${r.score})` : ''}</span></div></div>`;
    const bf = strand === '+' ? r.bFrom : m - r.bFrom + 1, bt = strand === '+' ? r.bTo : m - r.bTo + 1;
    const st = r.ident >= 12 && r.ident / r.len >= 0.9 ? 'bad' : r.ident >= 8 ? 'warn' : 'ok';
    return `<div class="struct ${st}"><div class="struct-h"><span><b>${strand === '+' ? 'Plus / Plus' : 'Plus / Minus'}</b></span><span>Score <b>${r.score}</b></span><span>Identities <b>${r.ident}/${r.len}</b> (${fx(r.ident / r.len * 100, 0)}%)</span><span>Gaps <b>${r.gaps}</b></span>
      <span class="hint">${strand === '+' ? 'same strand' : 'second oligo’s reverse complement'}</span></div>${alnPre(a.name, r.top, r.aFrom, r.aTo, r.mid, b.name, r.bot, bf, bt)}</div>`;
  };
  return `<p class="hint" style="margin:0">Both orientations, aligned here with blastn’s default scoring (match +2, mismatch −3, gap 5 + 2 per base). Plus / Minus similarity means the two oligos are complementary there; see the cross-dimer pane for its ΔG.</p>
    ${block(plus, '+')}${block(minus, '-')}
    <div class="bench-row">${extLink(align2URL({ label: 'oligo_1', seq: a.seq }, { label: 'oligo_2', seq: b.seq }), 'Align on NCBI', 'btn primary')}<button class="btn ghost" data-bl="copyfasta">Copy both as FASTA</button>
    <span class="hint">Opens NCBI’s “Align two sequences” page. If its boxes come up empty, paste the copied FASTA.</span></div>`;
}

function blastPaneBody(o) {
  if (o.seq.length < 2) return '<p class="empty">Enter a sequence to search.</p>';
  const { others, p } = pickPartner(o);
  const mode = BL.mode;
  const tabs = [['bind', 'Where each binds'], ['pair', 'Primer pair specificity'], ['align', 'Align the two']];
  let h = `<div class="bl-ctl"><div class="seg" role="group" aria-label="BLAST mode">${tabs.map(([k, t]) => `<button data-bl-mode="${k}" aria-pressed="${mode === k}">${t}</button>`).join('')}</div>
    ${others.length ? `<div class="partner"><label class="eyebrow" for="bl-partner">With</label>${partnerSelect('bl-partner', others, p)}</div>` : ''}</div>`;
  const pair = p ? orderPair(o, p) : [o];
  const qs = pair.map((x, k) => ({ label: `oligo_${k + 1}`, name: x.name || `Oligo ${k + 1}`, seq: x.seq }));
  if (mode === 'align') {
    if (!p) return h + '<p class="empty">Add a second oligo to the library to align the two.</p>';
    return h + alignHTML({ name: o.name || 'Oligo 1', seq: o.seq }, { name: p.name || 'Oligo 2', seq: p.seq });
  }
  const orgField = `<div class="field grow"><label for="bl-org">Organism</label><input id="bl-org" value="${esc(BL.org)}" placeholder="${mode === 'pair' ? 'Homo sapiens (NCBI default)' : 'all organisms'}" spellcheck="false" autocomplete="off"></div>`;
  if (mode === 'pair') {
    if (!p) return h + '<p class="empty">Add a second oligo to the library to check the pair.</p>';
    const [f, r] = pair;
    return h + `<div class="bl-form"><div class="field"><label for="bl-pbdb">Database</label><select id="bl-pbdb">${PB_DBS.map(([v, t]) => `<option value="${v}" ${BL.pbdb === v ? 'selected' : ''}>${t}</option>`).join('')}</select></div>${orgField}</div>
      <div class="bl-q"><div class="sq"><span class="eyebrow">Forward</span> <b>${esc(f.name)}</b> <span class="mono"><span class="lab">5′-</span>${colored(f.seq)}<span class="lab">-3′</span></span></div>
      <div class="sq"><span class="eyebrow">Reverse</span> <b>${esc(r.name)}</b> <span class="mono"><span class="lab">5′-</span>${colored(r.seq)}<span class="lab">-3′</span></span></div></div>
      <div class="bench-row">${extLink(primerBlastURL(f, r), 'Check pair in Primer-BLAST', 'btn primary')}
      <span class="hint">Opens NCBI Primer-BLAST with both primers filled in. Click “Get Primers” there; it lists every product the pair could make in that database. Primer-BLAST has no public API, so its results stay on NCBI’s site.</span></div>`;
  }
  // bind
  const form = `<div class="bl-form"><div class="field"><label for="bl-db">Database</label><select id="bl-db">${BLAST_DBS.map(([v, t]) => `<option value="${v}" ${BL.db === v ? 'selected' : ''}>${t}</option>`).join('')}</select></div>${orgField}</div>
    <div class="bl-q">${qs.map(q => `<div class="sq"><b>${esc(q.name)}</b> <span class="mono"><span class="lab">5′-</span>${colored(q.seq)}<span class="lab">-3′</span></span></div>`).join('')}</div>`;
  if (!IS_DESKTOP) {
    return h + form + `<div class="bench-row">${extLink(blastWebURL(qs), qs.length > 1 ? 'BLAST both on NCBI' : 'BLAST on NCBI', 'btn primary')}
      <span class="hint">Opens NCBI BLAST with ${qs.length > 1 ? 'both oligos' : 'the oligo'} filled in. In the Stackwise desktop app the hits are shown right here.</span></div>`;
  }
  const job = BJ.get(jobKey(qs)), busy = job && (job.status === 'submitting' || job.status === 'waiting');
  return h + form + `<div class="bench-row"><button class="btn primary" data-bl="run" ${busy ? 'disabled' : ''}>${job && job.status === 'done' ? 'Search again' : qs.length > 1 ? 'BLAST both oligos' : 'BLAST this oligo'}</button>
      <div class="field"><label for="bl-email">Email for NCBI</label><input id="bl-email" class="wide" type="email" value="${esc(BL.email)}" placeholder="optional" autocomplete="off"></div>
      <span class="hint">Sent to NCBI with the search (blastn tuned for short sequences, word size 7). Usually takes 30 s to a few minutes.</span></div>
    <div id="bl-out">${blastOutHTML(qs)}</div>`;
}
function blastSummary(o) {
  if (BL.mode !== 'bind' || !IS_DESKTOP) return '';
  const { p } = pickPartner(o), pair = p ? orderPair(o, p) : [o];
  const job = BJ.get(jobKey(pair.map((x, k) => ({ label: `oligo_${k + 1}`, seq: x.seq }))));
  if (!job) return '';
  if (job.status === 'done') { const prods = sharedProducts(job.results); return job.results.length > 1 ? pill(prods.some(x => x.f.end3 && x.r.end3) ? 'warn' : 'ok', prods.length + ' shared target' + (prods.length === 1 ? '' : 's')) : pill('neutral', job.results[0].hits.length + ' hits'); }
  if (job.status === 'error') return pill('bad', 'failed');
  return pill('neutral', 'running…');
}
function refreshBlast(full) {
  if ($('oligoview').hidden) return;
  const sec = document.querySelector('[data-pane=blast]'); if (!sec) return;
  const sum = sec.querySelector('.sum'); if (sum) sum.innerHTML = blastSummary(cur());
  if (closed.has('blast')) return;
  if (full) { const body = sec.querySelector('.pane-b'); if (body) body.innerHTML = blastPaneBody(cur()); return; }
  const out = $('bl-out'); if (!out) return;
  const o = cur(), { p } = pickPartner(o), pair = p ? orderPair(o, p) : [o];
  out.innerHTML = blastOutHTML(pair.map((x, k) => ({ label: `oligo_${k + 1}`, name: x.name, seq: x.seq })));
}
setInterval(() => { for (const j of BJ.values()) if (j.status === 'waiting') { refreshBlast(false); break; } }, 1000);

function wireBlast() {
  const panes = $('panes');
  panes.addEventListener('click', e => {
    const m = e.target.closest('[data-bl-mode]');
    if (m) { BL.mode = m.dataset.blMode; saveBL(); refreshBlast(true); return; }
    const b = e.target.closest('[data-bl]');
    if (b) {
      const o = cur(), { p } = pickPartner(o), pair = p ? orderPair(o, p) : [o];
      const qs = pair.map((x, k) => ({ label: `oligo_${k + 1}`, name: x.name || `Oligo ${k + 1}`, seq: x.seq }));
      if (b.dataset.bl === 'run') { const em = $('bl-email'); if (em) { BL.email = em.value.trim(); saveBL(); } runBlast(qs); }
      else if (b.dataset.bl === 'cancel') { const j = BJ.get(jobKey(qs)); if (j) { j.cancelled = true; BJ.delete(j.key); } refreshBlast(true); }
      else if (b.dataset.bl === 'copyfasta') copy(fasta(qs.map(q => ({ label: q.name.replace(/\s+/g, '_'), seq: q.seq }))), 'Both oligos copied as FASTA');
      return;
    }
    const row = e.target.closest('tr.hitrow');
    if (row && !e.target.closest('a')) { const id = row.dataset.hit; openHits.has(id) ? openHits.delete(id) : openHits.add(id); refreshBlast(false); }
  });
  panes.addEventListener('change', e => {
    const t = e.target;
    if (t.id === 'bl-partner') { partner = t.value; save(); renderPanes(cur(), lastRes); }
    else if (t.id === 'bl-db') { BL.db = t.value; saveBL(); refreshBlast(true); }
    else if (t.id === 'bl-pbdb') { BL.pbdb = t.value; saveBL(); refreshBlast(true); }
    else if (t.id === 'bl-org') { BL.org = t.value.trim(); saveBL(); refreshBlast(true); }
    else if (t.id === 'bl-email') { BL.email = t.value.trim(); saveBL(); }
  });
  panes.addEventListener('input', e => { if (e.target.id === 'bl-org') { BL.org = e.target.value.trim(); saveBL(); } });
  panes.addEventListener('keydown', e => { if (e.target.id === 'bl-org' && e.key === 'Enter') { e.preventDefault(); e.target.blur(); } });
  document.addEventListener('click', e => {
    const a = e.target.closest('a[data-ext]');
    if (a && IS_DESKTOP) { e.preventDefault(); openExternal(a.href); }
  });
}
