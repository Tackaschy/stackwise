import requests, re, time, json, os, traceback, threading
OUT = 'out2'; os.makedirs(OUT, exist_ok=True)
F, R = 'GAAGGTGAAGGTCGGAGTCA', 'TTGAGGTCAATGAAGGGGTC'
BL = 'https://blast.ncbi.nlm.nih.gov/'
log = []; lock = threading.Lock()
def note(s):
    with lock: print(s, flush=True); log.append(time.strftime('%H:%M:%S ') + s)
def save(name, text): open(f'{OUT}/{name}', 'w').write(text); note(f'saved {name} ({len(text)})')
def S():
    s = requests.Session(); s.headers['User-Agent'] = 'Stackwise-probe/0.2'; return s

def poll(s, tag, rid, rtoe, maxmin=28):
    t0 = time.time(); time.sleep(max(10, rtoe)); k = 0
    while time.time() - t0 < maxmin * 60:
        i = s.get(BL + 'Blast.cgi', params=dict(CMD='Get', FORMAT_OBJECT='SearchInfo', RID=rid), timeout=60).text
        st = (re.search(r'Status=(\w+)', i) or [None, None])[1]; k += 1
        if st != 'WAITING':
            note(f'{tag}: {st} after {int(time.time()-t0)} s ({k} polls); hits: {re.search(r"ThereAreHits=(\w+)", i) and re.search(r"ThereAreHits=(\w+)", i)[1]}')
            for fmt in ['JSON2_S']:
                g = s.get(BL + 'Blast.cgi', params=dict(CMD='Get', RID=rid, FORMAT_TYPE=fmt, DESCRIPTIONS=50, ALIGNMENTS=50), timeout=120)
                save(f'{tag}_result_{fmt}.txt', g.text)
            return
        time.sleep(15)
    note(f'{tag}: still WAITING after {maxmin} min')

short = dict(PROGRAM='blastn', WORD_SIZE='7', EXPECT='1000', FILTER='F', NUCL_REWARD='1', NUCL_PENALTY='-3', GAPCOSTS='5 2', HITLIST_SIZE='50', FORMAT_TYPE='JSON2_S')
def bind(tag, extra):
    try:
        s = S(); r = s.post(BL + 'Blast.cgi', data=dict(short, CMD='Put', TOOL='stackwise-probe', QUERY=f'>oligo_1\n{F}\n>oligo_2\n{R}', **extra), timeout=120)
        rid = re.search(r'RID = (\S+)', r.text); rtoe = re.search(r'RTOE = (\d+)', r.text)
        note(f'{tag}: rid={rid and rid[1]} rtoe={rtoe and rtoe[1]}')
        if rid: poll(s, tag, rid[1], int(rtoe[1]) if rtoe else 15)
        else: save(f'{tag}_put.html', r.text)
    except Exception: note(f'{tag}: EXC ' + traceback.format_exc()[-500:])

def align(tag, data):
    try:
        s = S(); r = s.post(BL + 'BlastAlign.cgi', data=data, timeout=120)
        save(f'{tag}_put.html', r.text)
        rid = re.search(r'RID = (\S+)', r.text) or re.search(r'RID=([A-Z0-9]{8,})', r.text); rtoe = re.search(r'RTOE = (\d+)', r.text)
        note(f'{tag}: status {r.status_code} url {r.url} rid={rid and rid[1]}')
        if rid: poll(s, tag, rid[1], int(rtoe[1]) if rtoe else 5, maxmin=10)
    except Exception: note(f'{tag}: EXC ' + traceback.format_exc()[-500:])

A, B = 'GAAGGTGAAGGTCGGAGTCAACGGATTTGGTCGTATTGG', 'CCAATACGACCAAATCCGTTGACTCCGACCTTCACCTTC'
jobs = [
  ('bind_refseqrna_human', bind, dict(DATABASE='refseq_rna', ENTREZ_QUERY='"Homo sapiens"[organism]')),
  ('bind_refseqrna_all', bind, dict(DATABASE='refseq_rna')),
  ('bind_corent_all', bind, dict(DATABASE='core_nt')),
  ('bind_megaoff_refseqrna', bind, dict(DATABASE='refseq_rna', MEGABLAST='off', ENTREZ_QUERY='"Homo sapiens"[organism]')),
  ('align_put', align, dict(CMD='Put', PROGRAM='blastn', QUERY=f'>a\n{A}', SUBJECTS=f'>b\n{B}', WORD_SIZE='7', EXPECT='1000', FILTER='F', FORMAT_TYPE='JSON2_S')),
  ('align_put_short', align, dict(CMD='Put', PROGRAM='blastn', QUERY=f'>a\n{F}', SUBJECTS=f'>b\n{R}', WORD_SIZE='7', EXPECT='1000', FILTER='F', FORMAT_TYPE='JSON2_S')),
  ('align_request', align, dict(CMD='request', PROGRAM='blastn', BLAST_SPEC='blast2seq', QUERY=f'>a\n{A}', SUBJECTS=f'>b\n{B}', WORD_SIZE='7', EXPECT='1000', FILTER='F', PAGE_TYPE='BlastSearch')),
]

def primerblast(tag, left, right, extra):
    try:
        s = S(); PB = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/'
        data = dict(PRIMER_LEFT_INPUT=left, PRIMER_RIGHT_INPUT=right, ORGANISM='Homo sapiens', PRIMER_SPECIFICITY_DATABASE='refseq_mrna',
                    SEARCH_SPECIFIC_PRIMER='on', CMD='request', TOTAL_PRIMER_SPECIFICITY_MISMATCH='1', PRIMER_3END_SPECIFICITY_MISMATCH='1',
                    MISMATCH_REGION_LENGTH='5', TOTAL_MISMATCH_IGNORE='6', MAX_TARGET_SIZE='4000', HITSIZE='50000', EVALUE='30000', WORD_SIZE='7',
                    NUM_TARGETS_WITH_PRIMERS='1000', MAX_TARGET_PER_TEMPLATE='100', PAGE_TYPE='PrimerDesign',
                    MONO_CATIONS='50.0', DIVA_CATIONS='1.5', CON_ANEAL_OLIGO='250.0', CON_DNTPS='0.2', SALT_FORMULAR='1', TM_METHOD='1')
        data.update(extra)
        t0 = time.time(); r = s.post(PB + 'primertool.cgi', files={k: (None, v) for k, v in data.items()}, timeout=180)
        jk = re.search(r'job_key=([A-Za-z0-9_\-]+)', r.text); jk = jk and jk[1]
        note(f'{tag}: submit {r.status_code} in {time.time()-t0:.1f}s job_key={jk} done_now={"prPairInfo" in r.text}')
        seen = set(); page = r.text; k = 0
        while jk and 'prPairInfo' not in page and time.time() - t0 < 900:
            time.sleep(5); k += 1
            page = s.get(PB + 'primertool.cgi', params=dict(job_key=jk), timeout=120).text
            stat = re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', ' ', page))
            m = re.search(r'Status (.{0,120})', stat); key = m and m[1][:60]
            if key not in seen: seen.add(key); save(f'{tag}_poll_{k}.html', page); note(f'{tag}: poll {k} status "{key}"')
            if re.search(r'class="error"|No primer|not found', page) and 'prPairInfo' not in page and k > 2: break
        save(f'{tag}_final.html', page); note(f'{tag}: final after {time.time()-t0:.0f}s pairs={page.count("prPairInfo")}')
    except Exception: note(f'{tag}: EXC ' + traceback.format_exc()[-500:])

jobs += [
  ('pb_gapdh_nt', primerblast, (F, R, dict(PRIMER_SPECIFICITY_DATABASE='nt'))),
  ('pb_actb', primerblast, ('CATGTACGTTGCTATCCAGGC', 'CTCCTTAATGTCACGCACGAT', {})),
  ('pb_random', primerblast, ('ACGTTGCATCGATCGGATCA', 'TTGACCGATGCATGCTAGCA', {})),
  ('pb_bad', primerblast, ('ACGT', 'TTGA', {})),
]
th = []
for tag, fn, arg in jobs:
    t = threading.Thread(target=(lambda f=fn, a=arg, g=tag: f(g, *a) if isinstance(a, tuple) else f(g, a))); t.start(); th.append(t); time.sleep(4)
for t in th: t.join()
open(f'{OUT}/LOG.txt', 'w').write('\n'.join(log))
