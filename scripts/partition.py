# --------------------------------------------------------------------------
# partition.py
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
import csv, collections, re, sys
WRITE = '--write' in sys.argv
rows = list(csv.DictReader(open("metadata/FINAL_sample_set_M1.tsv", newline=''), delimiter='\t'))
def a(r,k):
    for kv in (r.get('all_attributes') or '').split(' | '):
        x,_,y = kv.partition('=')
        if x.strip().lower()==k: return y.strip()
    return ''
def low(r,k): return a(r,k).strip().lower()
def al(r): return r.get('sample_alias','')

def classify(r):
    c = r['cohort_id']
    if c=='CABRAL17_US': return 'exclude','mixed amplicon'
    if c=='ZHU20_CN' and not a(r,'isolate').startswith('MSWS'):
        return 'exclude','repeated collection method'
    if c=='AGP_ORAL':
        return ('train','') if all(low(r,k)=='true' for k in
                ('subset_ibd','subset_diabetes','subset_antibiotic_history')) \
            else ('exclude','AGP health flags fail or absent')
    if c=='LEE17_CN':
        v=a(r,'host_disease')
        if 'Healthy' in v: return 'train',''
        return ('validate','OSCC') if 'OSCC' in v else ('validate','precancerous lesion')
    if c=='KUMPITSCH20_EU':
        if 'ITS' in a(r,'target gene'): return 'exclude','fungal ITS amplicon'
        return ('train','') if low(r,'nose/mouth/teeth/throat disorder')=='healthy' else ('validate','OSCC')
    if c=='USA_30897':
        if not al(r).endswith('_Sa'): return 'exclude','cheek mucosa, wrong body site'
        return ('train','') if low(r,'study_hivstatus')=='hivneg' else ('validate','HIV')
    if c=='USA_74379':
        return ('validate','GERD') if low(r,'host_disease')=='gerd' else ('train','')
    if c=='CHEN18_US':
        return ('train','') if low(r,'diseased')=='n' else ('validate','oral disease, unspecified')
    if c=='TAIWAN_39064':
        v=a(r,'host disease status')
        if v=='Health': return 'train',''
        return ('validate','OSCC') if v=='OSCC' else ('validate','oral verrucous hyperplasia')
    if c=='RELVAS21_EU':
        v=a(r,'grades_global_oral_health')
        return ('train','') if v=='Global_Grade_0' else ('validate','oral health '+v)
    if c=='TAIWAN_37501':
        s=a(r,'host subject id')
        if s.startswith('SCC-'): return 'validate','OSCC'
        if s.startswith('OSF-'): return 'validate','oral submucous fibrosis'
        return 'exclude','status unverified'
    if c=='SOHAIL19_QA':
        return ('train','') if low(r,'experimental_factor')=='control' else ('validate','type 2 diabetes')
    if c in ('SPAIN_74796','PORTUGAL_74247'):
        return ('train','') if low(r,'health_condition')=='healthy' else ('validate','periodontitis')
    if c=='WANG16_CN':
        return ('train','') if low(r,'host_disease')=='healthy' else ('validate','oral lichen planus')
    if c=='HIJAZI20_EU':
        return ('train','') if 'healthy' in low(r,'host_disease') else ('validate',a(r,'host_disease'))
    if c=='CHINA_43872':
        return ('train','') if a(r,'isolate').startswith('RC-HD') else ('exclude','status unverified')
    if c=='WOLF17_EU':
        return ('train','') if low(r,'nose/mouth/teeth/throat disorder')=='none' else ('validate','OSCC')
    return 'train',''

def pid(r):
    c=r['cohort_id']
    if c=='USA_27304':  return a(r,'subjectid')
    if c=='USA_30897':  return a(r,'study_subjectid')
    if c=='CHEN18_US':  return a(r,'subject')
    if c=='ZHU20_CN_2': return re.sub(r'^DAY\d(UWS|MSWS)','',a(r,'isolate'))
    if c=='INDIA_36346':
        m=re.match(r'^[A-Za-z]+(\d+)_',al(r)); return m.group(1) if m else al(r)
    return r['run_accession']
def rank(r):
    c=r['cohort_id']
    if c=='USA_27304':
        m=re.search(r'_T(\d+)_',al(r)); return int(m.group(1)) if m else 999
    if c=='CHEN18_US':
        m=re.match(r'D(\d+)S',al(r)); return int(m.group(1)) if m else 0
    if c=='ZHU20_CN_2':
        v=a(r,'isolate'); return 0 if v.startswith('DAY1') and 'UWS' in v else 1
    if c=='INDIA_36346': return 0 if al(r).startswith('A') else 1
    return 0

out = [(r,)+classify(r) for r in rows]
groups = collections.defaultdict(list)
for t in out:
    if t[1]=='exclude': continue
    groups[(t[0]['cohort_id'], t[1])].append(t)
keep=[]
for grp in groups.values():
    byp=collections.defaultdict(list)
    for t in grp: byp[pid(t[0])].append(t)
    for v in byp.values(): keep.append(sorted(v,key=lambda t:rank(t[0]))[0])

tr=[t for t in keep if t[1]=='train']; va=[t for t in keep if t[1]=='validate']
print(f"{'cohort':18s} {'rows':>5s} {'TRAIN':>6s} {'VALID':>6s}  validation labels")
for c in sorted({t[0]['cohort_id'] for t in out}):
    n=sum(1 for t in out if t[0]['cohort_id']==c)
    a1=sum(1 for t in tr if t[0]['cohort_id']==c); a2=sum(1 for t in va if t[0]['cohort_id']==c)
    labs=collections.Counter(t[2] for t in va if t[0]['cohort_id']==c)
    mark = '' if a1==n and not a2 else '  *'
    print(f"{c:18s} {n:5d} {a1:6d} {a2:6d}  {', '.join(f'{k}×{v}' for k,v in labs.most_common())}{mark}")
print(f"\nTRAIN {len(tr)} samples / {len({t[0]['cohort_id'] for t in tr})} cohorts")
print(f"VALID {len(va)} samples / {len({t[0]['cohort_id'] for t in va})} cohorts")
print("conditions:", dict(collections.Counter(t[2] for t in va)))
small=[(c,sum(1 for t in tr if t[0]['cohort_id']==c)) for c in {t[0]['cohort_id'] for t in tr}]
print("training cohorts under 30 samples:", sorted([x for x in small if x[1]<30], key=lambda x:x[1]))
if WRITE:
    with open("metadata/partition.tsv","w") as f:
        f.write("run_accession\tcohort_id\tclass\tlabel\n")
        for t in sorted(keep, key=lambda t:(t[0]['cohort_id'], t[0]['run_accession'])):
            f.write(f"{t[0]['run_accession']}\t{t[0]['cohort_id']}\t{t[1]}\t{t[2]}\n")
    print("\nwritten: metadata/partition.tsv")
