#!/usr/bin/env python3
"""Render real Circos and verify SVG geometry; requires the supported dependencies.

Use --output DIR to preserve generated artifacts and commands. Nothing in the
repository's example output directory is modified.
"""
import argparse
import json
import math
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
NS='{http://www.w3.org/2000/svg}'


def run(work, args, log):
    p=subprocess.run(list(map(str,args)),cwd=work,text=True,capture_output=True,timeout=120)
    log.append(dict(argv=list(map(str,args)),rc=p.returncode,stdout=p.stdout,stderr=p.stderr))
    if p.returncode:
        raise RuntimeError(f'{args[0]} exited {p.returncode}:\n{p.stdout}\n{p.stderr}')
    return p


def check_svg(path, count):
    root=ET.parse(path).getroot()
    assert root.tag==NS+'svg'
    groups={g.get('id'):g for g in root.iter(NS+'g') if g.get('id')}
    plots=[g for key,g in groups.items() if re.fullmatch('plot[0-9]+',key)]
    assert len(plots)==count, groups.keys()
    assert len(groups['ideograms'].findall(NS+'path'))>=2
    assert len(groups['track_0'].findall(NS+'path'))>=1, 'no rendered link ribbon'
    for plot in plots:
        assert plot.findall('.//'+NS+'path'), 'empty plotted ring'
    for el in root.iter():
        for key in ('d','x','y','width','height'):
            assert not re.search(r'\b(nan|inf)\b',el.get(key,''),re.I)
    return root,groups


def smoke(work, legend=True):
    work=Path(work); work.mkdir(parents=True,exist_ok=True)
    # Refuse existing inputs/outputs: a smoke run must not overwrite user results.
    for name in ('genome.fa','iwgc_circos','iwgc_circos_data','fonts'):
        if (work/name).exists(): raise ValueError(f'Smoke directory is not fresh: {work/name}')
    shutil.copytree(ROOT/'iwgc_circos',work/'iwgc_circos',ignore=shutil.ignore_patterns('tmp'))
    shutil.copytree(ROOT/'fonts',work/'fonts',ignore=shutil.ignore_patterns('ArialBold.ttf'))
    shutil.copy(ROOT/'iwgc_circos_template.config',work)
    (work/'iwgc_circos/tmp').mkdir()
    (work/'genome.fa').write_text('>Chr1\n'+('ACGT'*30000)+'\n>scaffold_T1\n'+('GCGC'*20000)+('ATAT'*10000)+'\n')
    def gff(kind,a,b,attrs): return f'Chr1\tsynthetic\t{kind}\t{a}\t{b}\t.\t+\t.\t{attrs}\n'
    (work/'genes.gff3').write_text(gff('gene',1,30000,'ID=a')+gff('gene',40001,80000,'ID=b'))
    (work/'repeats.gff3').write_text(gff('repeat_region',1,40000,'Classification=LTR/Gypsy')+gff('repeat_region',80001,120000,'Classification=LTR/Copia'))
    (work/'ages.list').write_text('Chr1:1..40000\t'+'\t'.join(['.']*10)+'\t0\nChr1:40001..80000\t'+'\t'.join(['.']*10)+'\t20\n')
    (work/'links.coords').write_text('Chr1 0 40000 scaffold_T1 60000 100000\n')
    log=[]
    try:
        run(work,['bash',ROOT/'iwgc_circos_tracks.sh','genome.fa','-gene','genes.gff3','-repeat','repeats.gff3',
                  '-intact','repeats.gff3','-ltr-dating','ages.list','-links','links.coords','-gc','-window',40000,'-out','iwgc_circos_data'],log)
        run(work,['bash',ROOT/'create_configs.sh','-gap'],log)
        run(work,['circos','-conf','iwgc_circos/iwgc_circos.conf','-outputdir','iwgc_circos/tmp'],log)
        svg=work/'iwgc_circos/tmp/iwgc_circos.svg'
        root,groups=check_svg(svg,6)
        png=work/'iwgc_circos/tmp/iwgc_circos.png'
        data=png.read_bytes(); assert data[:8]==b'\x89PNG\r\n\x1a\n'
        assert struct.unpack('>II',data[16:24])==(3000,3000)
        if legend:
            out=work/'iwgc_circos/tmp/labeled.svg'
            run(work,['bash',ROOT/'add_legend.sh','-svg',svg,'-out',out],log)
            labeled,lg=check_svg(out,6)
            labels=lg['iwgc-gap-labels'].findall(NS+'text')
            assert [x.text for x in labels]==['a.','b.','c.','d.','e.','f.','g.']
            ys=[float(x.get('y')) for x in labels]
            assert ys==sorted(ys) and len(set(ys))==7, ys
            for label in labels:
                assert 0<float(label.get('x'))<3000 and 0<float(label.get('y'))<3000
                assert 10<=float(label.get('font-size'))<=48
            assert len(list(labeled.iter()))>len(list(root.iter()))
            # 0 degrees must put every label on the right-hand radial line.
            angle=work/'iwgc_circos/tmp/angle.svg'
            run(work,['bash',ROOT/'add_legend.sh','-svg',svg,'-out',angle,'--theta-deg','0'],log)
            ar=ET.parse(angle).getroot()
            rotated=[x for x in ar.iter(NS+'text') if x.get('class')=='iwgc-gap-label']
            assert max(float(x.get('y')) for x in rotated)-min(float(x.get('y')) for x in rotated)<.01
            assert all(float(x.get('x'))>1500 for x in rotated)
        return {'status':'passed','tracks':6,'chromosomes':2,'links_rendered':True,'legend':legend}
    finally:
        (work/'commands.json').write_text(json.dumps(log,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output',type=Path); p.add_argument('--no-legend',action='store_true'); args=p.parse_args()
    if args.output: print(json.dumps(smoke(args.output,not args.no_legend)))
    else:
        with tempfile.TemporaryDirectory(prefix='iwgc-smoke-') as work:
            print(json.dumps(smoke(work,not args.no_legend)))


if __name__=='__main__': main()
