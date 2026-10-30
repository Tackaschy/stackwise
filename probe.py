# Runs real requests against NCBI and saves every reply, so Stackwise can be built against them.
import requests, re, time, json, os, sys, traceback
OUT = 'out'; os.makedirs(OUT, exist_ok=True)
S = requests.Session(); S.headers['User-Agent'] = 'Stackwise-probe/0.1'
F, R = 'GAAGGTGAAGGTCGGAGTCA', 'TTGAGGTCAATGAAGGGGTC'   # human GAPDH primer pair
BL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
log = []
def save(name, text): open(f'{OUT}/{name}', 'w').write(text); log.append(f'{name}: {len(text)} chars')
def note(s): print(s, flush=True); log.append(s)

def blast_run(tag, put):
    try:
        r = S.post(BL, data=dict(put, CMD='Put', TOOL='stackwise-probe'), timeout=120)
        save(f'{tag}_put.html', r.text)
        rid = re.search(r'RID = (\S+)', r.text); rtoe = re.search(r'RTOE = (\d+)', r.text)
        note(f'{tag}: put status {r.status_code} rid={rid and rid.group(1)} rtoe={rtoe and rtoe.group(1)}')
        if not rid: return
        rid = rid.group(1); time.sleep(max(10, int(rtoe.group(1)) if rtoe else 15))
        for k in range(60):
            i = S.get(BL, params=dict(CMD='Get', FORMAT_OBJECT='SearchInfo', RID=rid), timeout=60).text
            st = re.search(r'Status=(\w+)', i); st = st and st.group(1)
            if k == 0: save(f'{tag}_info_first.txt', i)
            if st != 'WAITING': save(f'{tag}_info_last.txt', i); note(f'{tag}: status {st} after {k} polls'); break
            time.sleep(15)
        for fmt in ['JSON2_S', 'Text']:
            g = S.get(BL, params=dict(CMD='Get', RID=rid, FORMAT_TYPE=fmt, DESCRIPTIONS=50, ALIGNMENTS=50), timeout=120)
            save(f'{tag}_result_{fmt}.txt', g.text)
    except Exception: note(f'{tag}: EXCEPTION ' + traceback.format_exc()[-600:])

short = dict(PROGRAM='blastn', WORD_SIZE='7', EXPECT='1000', FILTER='F', NUCL_REWARD='1', NUCL_PENALTY='-3', GAPCOSTS='5 2', SHORT_QUERY_ADJUST='true', HITLIST_SIZE='50', FORMAT_TYPE='JSON2_S')
# 1. both primers against core_nt, human only
blast_run('bind', dict(short, DATABASE='core_nt', QUERY=f'>oligo_1\n{F}\n>oligo_2\n{R}', ENTREZ_QUERY='"Homo sapiens"[organism]'))
# 2. two-sequence alignment variants
A, B = 'GAAGGTGAAGGTCGGAGTCAACGGATTTGGTCGTATTGG', 'CCAATACGACCAAATCCGTTGACTCCGACCTTCACCTTC'  # B = revcomp(A)
blast_run('bl2seq_subjects', dict(short, QUERY=f'>a\n{A}', SUBJECTS=f'>b\n{B}'))
blast_run('bl2seq_spec', dict(short, QUERY=f'>a\n{A}', SUBJECTS=f'>b\n{B}', BLAST_SPEC='blast2seq', DATABASE='n/a'))
blast_run('bl2seq_short', dict(short, QUERY=f'>a\n{F}', SUBJECTS=f'>b\n{R}', BLAST_SPEC='blast2seq'))

# 3. Primer-BLAST, the way its own form submits
PB = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/'
try:
    form0 = S.get(PB + 'index.cgi', timeout=60).text; save('pb_form.html', form0)
    fields = dict(re.findall(r'<input[^>]*name="([^"]+)"[^>]*value="([^"]*)"', form0))
    note('pb form inputs: ' + json.dumps(fields)[:3000])
    sel = re.findall(r'<select[^>]*name="([^"]+)"', form0); note('pb selects: ' + ', '.join(sel))
    data = dict(PRIMER_LEFT_INPUT=F, PRIMER_RIGHT_INPUT=R, ORGANISM='Homo sapiens', PRIMER_SPECIFICITY_DATABASE='refseq_mrna',
                SEARCH_SPECIFIC_PRIMER='on', CMD='request', PRIMER_PRODUCT_MIN='70', PRIMER_PRODUCT_MAX='1000', PRIMER_NUM_RETURN='10',
                TOTAL_PRIMER_SPECIFICITY_MISMATCH='1', PRIMER_3END_SPECIFICITY_MISMATCH='1', MISMATCH_REGION_LENGTH='5',
                TOTAL_MISMATCH_IGNORE='6', MAX_TARGET_SIZE='4000', HITSIZE='50000', EVALUE='30000', WORD_SIZE='7',
                NUM_TARGETS_WITH_PRIMERS='1000', NUM_TARGETS='20', MAX_TARGET_PER_TEMPLATE='100', PAGE_TYPE='PrimerDesign',
                MONO_CATIONS='50.0', DIVA_CATIONS='1.5', CON_ANEAL_OLIGO='50.0', CON_DNTPS='0.6', SALT_FORMULAR='1', TM_METHOD='1',
                PRIMER_MIN_TM='57.0', PRIMER_OPT_TM='60.0', PRIMER_MAX_TM='63.0', PRIMER_MAX_DIFF_TM='3')
    r = S.post(PB + 'primertool.cgi', files={k: (None, v) for k, v in data.items()}, timeout=180)
    save('pb_submit.html', r.text); note(f'pb submit status {r.status_code} url {r.url}')
    hdr = {k: v for k, v in r.headers.items() if k.lower().startswith('access-control') or k.lower() in ('content-type', 'refresh')}
    note('pb headers ' + json.dumps(hdr))
    jk = re.search(r'job_key=([A-Za-z0-9_\-]+)', r.text); jk = jk and jk.group(1); note(f'pb job_key={jk}')
    page = r.text
    for k in range(40):
        if not jk: break
        if not re.search(r'(?i)(being processed|please wait|refresh|still running)', page) and k > 0: break
        time.sleep(15)
        page = S.get(PB + 'primertool.cgi', params=dict(job_key=jk), timeout=120).text
        if k == 0: save('pb_poll_first.html', page)
    save('pb_result.html', page); note(f'pb polls done, result has "Primer pair 1": {"Primer pair 1" in page}, products: {page.count("product length")}')
except Exception: note('pb: EXCEPTION ' + traceback.format_exc()[-800:])

open(f'{OUT}/LOG.txt', 'w').write('\n'.join(log))
