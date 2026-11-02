from pathlib import Path
import hashlib,json,subprocess,tempfile
root=Path(__file__).resolve().parents[3]; work=Path(tempfile.mkdtemp(prefix='iwgc-telomere-final-')); observations=[]

def run(label,argv,cwd=work):
 p=subprocess.run(list(map(str,argv)),cwd=cwd,capture_output=True,text=True,timeout=30)
 record=dict(case=label,argv=list(map(str,argv)),rc=p.returncode,stdout=p.stdout,stderr=p.stderr)
 observations.append(record); return record

# >50 exact motif copies at each end; 30 kb allows separate 10 kb end bins.
sequence='TTTAGGG'*100+'ACGT'*7150+'CCCTAAA'*100
assert len(sequence)==30000 and sequence[:10000].count('TTTAGGG')==100 and sequence[-10000:].count('CCCTAAA')==100
one=work/'single.fa'; one.write_text('>Chr1\n'+sequence+'\n')
two=work/'double.fa'; two.write_text('>Chr1\n'+sequence+'\n>Chr2\n'+sequence+'\n')
small=work/'minimal.fa'; small.write_text('>Chr1\n'+'TTTAGGG'*60+'ACGT'*5+'CCCTAAA'*60+'\n')
for name,fasta in [('minimal',small),('single',one),('double',two)]:
 out=work/(name+'-explore'); out.mkdir()
 r=run(name+' tidk explore',['tidk','explore','-d',out,'-f',fasta,'-o','probe','-e','tsv','-l','7','-t','50'])
 r['candidate_table']=(out/'probe.txt').read_text(); r['candidate_locations']=(out/'probe_telomeric_locations.tsv').read_text()
 assert r['rc']==0
 if name in ('minimal','single'): assert len(r['candidate_table'].splitlines())==1
 else: assert len(r['candidate_table'].splitlines())>=2
search=work/'search'; search.mkdir()
r=run('single explicit known motif search',['tidk','search','-d',search,'-f',one,'-o','known','-e','csv','-s','TTTAGGG','-w','10000'])
r['window_counts']=(search/'known_telomeric_repeat_windows.csv').read_text(); assert r['rc']==0
for name,fasta in [('single',one),('double',two)]:
 r=run(name+' iwgc-circos wrapper',['bash',root/'iwgc_circos_tracks.sh',fasta,'-telomere','-keep-temp','-out',work/(name+'-wrapper')])
 kary=work/(name+'-wrapper')/(fasta.name+'_karyotype.circos')
 r['band_rows']=[line for line in kary.read_text().splitlines() if line.startswith('band ')]
 if name=='single': assert r['rc']!=0 and 'no telomere report' in r['stderr'] and not r['band_rows']
 else: assert r['rc']==0 and len(r['band_rows'])==4
report=dict(expected_biology='The one-record 30 kb sequence contains 100 exact plant motif copies at each end. If detected, two marker bands are expected; missing detection is not biological absence.',
 expected_candidate_contract='Retain the current detector. A missing report must cause nonzero exit and must be reported as unavailable; no fallback or duplicated biological input is introduced.',
 observations=observations,input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (one,two,small)})
text=json.dumps(report,indent=2).replace(str(root),'<repo>').replace(str(work),'<scratch>')
(work/'report.json').write_text(text+'\n'); print(work)
for r in observations:
 print(r['case'],r['rc'],r.get('candidate_table','').strip(),r.get('window_counts','').strip())
