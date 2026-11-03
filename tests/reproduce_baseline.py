#!/usr/bin/env python3
"""Run selected defect probes against frozen original Git objects, not candidate code.

Needs the dependency image/native environment. All outputs go to a new scratch
folder. Reports observed values; this is not a candidate test or release gate.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
BASE='fe224c81c41e4aee5905e587ed01e0c1829ca0a3'


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output',type=Path); args=p.parse_args()
    work=args.output or Path(tempfile.mkdtemp(prefix='iwgc-baseline-repro-'))
    work.mkdir(parents=True,exist_ok=True)
    if any(work.iterdir()): p.error('output must be empty')
    script=work/'original_tracks.sh'
    script.write_text(subprocess.run(['git','-C',str(ROOT),'show',BASE+':iwgc_circos_tracks.sh'],check=True,capture_output=True,text=True).stdout)
    for name in ('genome.fa','genes.gff3','mummer.coords','links.coords'):
        shutil.copy(ROOT/'tests/fixtures'/name,work/name)
    # Coordinate-only probe excludes the attribute-selection trap.
    genes=work/'coordinate.gff3'; genes.write_text('\n'.join((work/'genes.gff3').read_text().splitlines()[:4])+'\n')
    report=[]
    def run(name,*opts,fasta=None,timeout=10):
        out=work/name
        argv=['bash',str(script),str(fasta or work/'genome.fa'),'-out',str(out),*map(str,opts)]
        try:
            result=subprocess.run(argv,capture_output=True,text=True,timeout=timeout)
            row=dict(case=name,argv=argv,rc=result.returncode,stdout=result.stdout,stderr=result.stderr,
                     tracks={f.name:f.read_text() for f in out.glob('*.circos')})
        except subprocess.TimeoutExpired:
            row=dict(case=name,argv=argv,timeout_seconds=timeout)
        report.append(row)
    run('gff-coordinate','-gene',genes,'-window',4)
    run('gff-attribute-selection','-gene',work/'genes.gff3','-window',4)
    run('mummer-reverse','-links',work/'mummer.coords')
    run('ascending-top-n','-links',work/'links.coords','-top-n',2)
    run('zero-window','-gc','-window',0)
    run('missing-value','-gene',timeout=.25)
    (work/'empty.list').touch()
    run('empty-ltr','-ltr-dating',work/'empty.list','-window',4)
    stale=work/'stale-intact'; stale.mkdir()
    (stale/'genome.fa_intactTE_Gypsy.gff').write_text('Chr1\ts\trepeat_region\t1\t4\t.\t+\t.\tClassification=LTR/Gypsy\n')
    run('stale-intact','-intact',work/'empty.list','-window',4)

    (work/'merge.coords').write_text('Chr1 0 4 Chr2 7 11\nChr1 3 6 Chr2 0 3\n')
    run('disjoint-second-merge','-links',work/'merge.coords')
    seq=work/'precision.fa'; seq.write_text('>Chr1\n'+'A'*20000000+'\n')
    g=work/'precision.gff3'; g.write_text('Chr1\ts\tgene\t1\t9999999\t.\t+\t.\tID=a\nChr1\ts\tgene\t10000001\t20000000\t.\t+\t.\tID=b\n')
    run('precision','-gene',g,'-window',10000000,fasta=seq)
    # Harmless command-injection marker stays within this run's private subtree.
    unsafe=work/'genome$(touch${IFS}IWGC_MARKER).fa'; unsafe.write_text('>Chr1\n'+'ACGT'*100+'\n')
    if shutil.which('quartet.py'):
        run('quartet-shell-filename','-telomere','-keep-temp',fasta=unsafe)
        report[-1]['marker_created']=bool(list(work.rglob('IWGC_MARKER')))
    config=work/'original_configs.sh'
    config.write_text(subprocess.run(['git','-C',str(ROOT),'show',BASE+':create_configs.sh'],check=True,capture_output=True,text=True).stdout)
    shutil.copytree(ROOT/'iwgc_circos',work/'iwgc_circos',ignore=shutil.ignore_patterns('tmp'))
    shutil.copy(ROOT/'iwgc_circos_template.config',work)
    data=work/'iwgc_circos_data'; data.mkdir()
    kary=data/'g_karyotype.circos'
    kary.write_text('chr - Chr1 1 0 12 chr2\nchr - scaffold_T1 sc1 0 8 chr3\n')
    result=subprocess.run(['bash',str(config)],cwd=work,capture_output=True,text=True)
    generated=work/'iwgc_circos/iwgc_circos.conf'
    report.append(dict(case='unusual-chromosome-id',rc=result.returncode,
                       chromosome_lines=[l for l in generated.read_text().splitlines() if l.startswith('chromosomes ')],stdout=result.stdout,stderr=result.stderr))
    kary.write_text('chr - Chr1 1 0 12 chr2\nchr - Chr2 2 0 12 chr3\n')
    suffixes=['gene_coverage','repeat_coverage','LTRdating_coverage',
              *[c+'_intactTE_coverage' for c in ('Gypsy','Copia','unknown','DTA','DTC','DTH')]]
    for suffix in suffixes: (data/f'g_{suffix}.circos').write_text('Chr1 0 4 0.5\n')
    (data/'g_links.circos').write_text('Chr1 0 4 Chr2 0 4\n')
    result=subprocess.run(['bash',str(config)],cwd=work,capture_output=True,text=True)
    report.append(dict(case='nine-track-link-radius',rc=result.returncode,
                       radius_lines=[l for l in generated.read_text().splitlines() if l.strip().startswith('radius ')],stdout=result.stdout,stderr=result.stderr))
    (data/'z_karyotype.circos').write_text('chr - Other other 0 10 chr4\n')
    result=subprocess.run(['bash',str(config)],cwd=work,capture_output=True,text=True)
    report.append(dict(case='ambiguous-karyotype',rc=result.returncode,
                       chosen_lines=[l for l in generated.read_text().splitlines() if l.startswith('karyotype ')],stdout=result.stdout,stderr=result.stderr))
    legend=work/'original_legend.sh'
    legend.write_text(subprocess.run(['git','-C',str(ROOT),'show',BASE+':add_legend.sh'],check=True,capture_output=True,text=True).stdout)
    svg=ROOT/'iwgc_circos/tmp/iwgc_circos.svg'; out=work/'original-angle.svg'
    result=subprocess.run(['bash',str(legend),'-svg',str(svg),'-out',str(out),'--theta-deg','0'],cwd=work,capture_output=True,text=True)
    labels=[x for x in ET.parse(out).getroot().iter('{http://www.w3.org/2000/svg}text') if x.get('class')=='iwgc-gap-label']
    report.append(dict(case='legend-angle-zero',rc=result.returncode,
                       xs=[float(x.get('x')) for x in labels],ys=[float(x.get('y')) for x in labels],stdout=result.stdout,stderr=result.stderr))
    text=json.dumps(dict(baseline_commit=BASE,observations=report),indent=2)
    text=text.replace(str(ROOT),'<repo>').replace(str(work),'<scratch>')
    (work/'report.json').write_text(text+'\n')
    print(work/'report.json')
    for row in report: print(row['case'],row.get('rc','timeout'),row.get('marker_created',''))


if __name__=='__main__': main()
