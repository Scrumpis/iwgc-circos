#!/usr/bin/env python3
"""Bounded release check against base sets and high-precision arithmetic.

Does not import iwgc-circos code or use bedtools to derive expected values.
Run in an existing dependency environment; --output must be a fresh directory.
"""
import argparse
from decimal import Decimal, getcontext
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

getcontext().prec = 45
ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    work = args.output or Path(tempfile.mkdtemp(prefix='iwgc-independent-'))
    work.mkdir(parents=True, exist_ok=True)
    if any(work.iterdir()):
        parser.error('output directory must be empty')
    observations, commands = [], []

    def run(argv, cwd=work, env=None):
        result = subprocess.run(list(map(str, argv)), cwd=cwd, env=env,
                                capture_output=True, text=True, timeout=60)
        commands.append(dict(argv=list(map(str, argv)), rc=result.returncode,
                             stdout=result.stdout, stderr=result.stderr))
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)

    def tracks(fasta, out, *flags):
        run(['bash', ROOT/'iwgc_circos_tracks.sh', fasta, '-out', out, *flags])

    def read(path):
        return [line.split() for line in path.read_text().splitlines()]

    def compare(name, path, expected, tolerance='1e-12'):
        got = read(path)
        assert [r[:3] for r in got] == [[str(x) for x in row[:3]] for row in expected], name
        errors = [abs(Decimal(row[3])-Decimal(str(want[3]))) for row, want in zip(got, expected)]
        assert all(error <= Decimal(tolerance) for error in errors), (name, got, expected)
        observations.append(dict(check=name, expected=[[str(x) for x in r] for r in expected],
                                 actual=got, tolerance=tolerance,
                                 maximum_absolute_error=str(max(errors, default=Decimal(0)))))

    def scale(values, root=False):
        transformed = [v.sqrt() if root else v for v in values]
        lo, hi = min(transformed), max(transformed)
        return [Decimal(0) if hi == lo else (v-lo)/(hi-lo) for v in transformed]

    sequences = {'ChrA': 'gCcNATGCNNATATGCGCA', 'ChrB': 'NNGCGCATATA'}
    fasta = work/'genome.fa'
    fasta.write_text(''.join(f'>{chrom}\n{seq}\n' for chrom, seq in sequences.items()))
    genes = [('ChrA',1,1), ('ChrA',4,8), ('ChrA',7,10), ('ChrA',18,19), ('ChrB',2,4)]
    repeats = [('ChrA',2,5,'Gypsy'), ('ChrA',4,9,'Gypsy'),
               ('ChrA',12,14,'Copia'), ('ChrA',14,19,'Copia'), ('ChrB',1,6,'Gypsy')]
    def gff(chrom, kind, start, end, attrs):
        return f'{chrom}\tdefined\t{kind}\t{start}\t{end}\t.\t+\t.\t{attrs}\n'
    gene_file = work/'genes.gff3'
    gene_file.write_text(''.join(gff(c,'gene',s,e,f'ID=g{i}') for i,(c,s,e) in enumerate(genes))
                         + gff('ChrA','mRNA',1,19,'ID=mrna;Parent=gene-parent'))
    repeat_file = work/'repeats.gff3'
    repeat_file.write_text(''.join(gff(c,'repeat_region',s,e,f'Classification=LTR/{family}')
                                  for c,s,e,family in repeats))
    ages = [('ChrA',1,2,0), ('ChrA',6,9,100), ('ChrA',8,14,300), ('ChrB',2,6,150)]
    age_file = work/'ages.list'
    age_file.write_text(''.join(f'{c}:{e if c=="ChrB" else s}..{s if c=="ChrB" else e}\t'
                               + '\t'.join(['.']*10) + f'\t{age}\n' for c,s,e,age in ages))
    windows = {
        'fixed': [('ChrA',0,6),('ChrA',6,12),('ChrA',12,18),('ChrA',18,19),('ChrB',0,6),('ChrB',6,11)],
        'sliding': [('ChrA',0,6),('ChrA',3,9),('ChrA',6,12),('ChrA',9,15),('ChrA',12,18),('ChrA',15,19),
                    ('ChrB',0,6),('ChrB',3,9),('ChrB',6,11)]}
    for mode, intervals in windows.items():
        out = work/mode
        flags = ['-gene',gene_file,'-repeat',repeat_file,'-intact',repeat_file,
                 '-ltr-dating',age_file,'-gc','-window',6,'-keep-temp']
        if mode == 'sliding': flags += ['-sliding','-step','0.5']
        tracks(fasta,out,*flags)
        actual_windows = read(out/'genome.fa_windows.bed')
        assert actual_windows == [[str(x) for x in r] for r in intervals]
        observations.append(dict(check=mode+' literal window coordinates', actual=actual_windows))
        # Enumerate actual one-based genomic positions. Overlaps are unions of sets.
        def coverage(features):
            covered = {c:set() for c in sequences}
            for c,start,end in features: covered[c].update(range(start,end+1))
            return [Decimal(len(covered[c] & set(range(start+1,end+1))))/Decimal(end-start)
                    for c,start,end in intervals]
        for label, features, square_root, suffix in [
            ('genes',genes,True,'gene_coverage'),
            ('repeats',[(c,s,e) for c,s,e,_ in repeats],True,'repeat_coverage'),
            *[(family,[(c,s,e) for c,s,e,f in repeats if f==family],False,family+'_intactTE_coverage')
              for family in ('Gypsy','Copia')]]:
            expected = [(*r,v) for r,v in zip(intervals,scale(coverage(features),square_root))]
            compare(mode+' '+label,out/f'genome.fa_{suffix}.circos',expected)
        gc = [Decimal(sum(base.upper() in ('G','C') for base in sequences[c][start:end]))/Decimal(end-start)
              for c,start,end in intervals]
        compare(mode+' GC',out/'genome.fa_gc.circos',[(*r,v) for r,v in zip(intervals,scale(gc))])
        means = []
        for c,start,end in intervals:
            relevant = [Decimal(age) for chrom,s,e,age in ages
                        if chrom==c and set(range(s,e+1)) & set(range(start+1,end+1))]
            means.append(sum(relevant)/Decimal(len(relevant)) if relevant else None)
        present = [x for x in means if x is not None]; lo,hi = min(present),max(present)
        expected = [(*r,0 if m is None else Decimal('1.5') if hi==lo else 2-(m-lo)/(hi-lo))
                    for r,m in zip(intervals,means)]
        compare(mode+' LTR means/inversion',out/'genome.fa_LTRdating_coverage.circos',expected,'1e-6')
        bed_expected = [['ChrA','0','1'],['ChrA','3','8'],['ChrA','6','10'],['ChrA','17','19'],['ChrB','1','4']]
        assert read(out/'genome.fa_genes_coords.bed') == bed_expected
    constant = work/'constant.fa'; constant.write_text('>ChrC\nGCGCGCGC\n')
    full = work/'full.gff3'; full.write_text(gff('ChrC','gene',1,8,'ID=full'))
    tracks(constant,work/'constant','-gene',full,'-gc','-window',4)
    expected = [('ChrC',0,4,0),('ChrC',4,8,0)]
    for suffix in ('gene_coverage','gc'):
        compare('constant '+suffix,work/'constant'/f'constant.fa_{suffix}.circos',expected)
    precise = work/'precision.fa'; precise.write_text('>ChrP\n'+'A'*3000000+'\n')
    close = work/'precision.gff3'
    close.write_text(gff('ChrP','gene',1,999998,'ID=a')+gff('ChrP','gene',1000001,1999999,'ID=b')
                     +gff('ChrP','gene',2000001,3000000,'ID=c'))
    tracks(precise,work/'precision','-gene',close,'-window',1000000)
    # Cancellation amplifies double precision by about 10^6; use a bounded 1e-9 tolerance.
    wanted = scale([Decimal(999998)/1000000,Decimal(999999)/1000000,Decimal(1)],True)
    compare('near-constant high precision',work/'precision/precision.fa_gene_coverage.circos',
            [('ChrP',0,1000000,wanted[0]),('ChrP',1000000,2000000,wanted[1]),('ChrP',2000000,3000000,wanted[2])],'1e-9')
    links = work/'links.coords'
    links.write_text('ChrA 12 16 ChrB 0 3\nChrA 3 7 ChrB 3 7\nChrB 0 2 ChrA 15 17\n'
                     'ChrA 8 11 ChrB 8 11\nChrA 0 4 ChrB 0 4\nChrA 1 2 ChrB 1 2\nChrA 16 18 ChrA 0 2\n')
    # Hand-constructed unions, plus deliberately separate spans and unequal side lengths.
    long = ['ChrA','0','7','ChrB','0','7']; next_long = ['ChrA','12','16','ChrB','0','3']
    for name,flags,want in [
        ('top2 ascending',['-top-n',2],[next_long,long]),
        ('top2 descending',['-top-n',2,'-link-order','dsc'],[long,next_long]),
        ('inclusive minimum',['-min-tlen',7],[next_long,long]),
        ('retained IDs before top-N',['-top-n',1,'-filter-chrs','^ChrA$'],[['ChrA','16','18','ChrA','0','2']])]:
        out=work/('link-'+str(len(observations))); tracks(fasta,out,'-links',links,*flags)
        got=read(out/'genome.fa_links.circos'); assert got==want,(name,got,want)
        observations.append(dict(check=name,expected=want,actual=got))
    for columns in (11,13):
        mummer=work/f'mummer{columns}.coords'
        tail=' 21.0 36.0' if columns==13 else ''
        mummer.write_text(f'1 4 4 1 4 4 99 19 11{tail} ChrA ChrB\n4 7 7 4 4 4 99 19 11{tail} ChrA ChrB\n')
        out=work/f'mummer{columns}'; tracks(fasta,out,'-links',mummer)
        got=read(out/'genome.fa_links.circos'); assert got==[long],got
        observations.append(dict(check=f'MUMmer {columns} closed/reverse/merge',expected=[long],actual=got))
    paf=work/'align.paf'; paf.write_text('ChrA\t19\t0\t19\t+\tChrA\t19\t0\t19\t19\t19\t60\n'
                                      'ChrA\t19\t0\t1\t-\tChrB\t11\t10\t11\t1\t1\t60\n')
    fake=work/'bin'; fake.mkdir(); exe=fake/'minimap2'
    exe.write_text('#!/bin/sh\ncat "'+str(paf)+'"\n'); exe.chmod(0o755)
    run(['bash',ROOT/'preprocess/minimap2_self_alignment.sh','-f',fasta,'-p',work/'alignment'],
        env=dict(os.environ,PATH=str(fake)+os.pathsep+os.environ['PATH']))
    got=read(work/'alignment.coords'); want=[['ChrA','0','1','ChrB','10','11']]; assert got==want
    observations.append(dict(check='PAF single-base and self-diagonal',expected=want,actual=got))
    report=dict(status='passed',checks=len(observations),oracle='literal windows; enumerated one-based base sets; Decimal precision 45; hand-constructed link unions',
                reviewer='same implementation-author Codex session; independent numerical construction, not independent peer review',
                observations=observations,commands=commands,
                software_version=(ROOT/'VERSION').read_text().strip(),
                source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
                               ['iwgc_circos_tracks.sh','tools/prepare_inputs.py','preprocess/minimap2_self_alignment.sh',
                                'tests/independent_check.py']},
                input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in work.iterdir() if p.is_file()})
    text=json.dumps(report,indent=2).replace(str(ROOT),'<repo>').replace(str(work),'<scratch>')
    (work/'report.json').write_text(text+'\n')
    print(json.dumps({'status':'passed','checks':len(observations),'report':str(work/'report.json')}))


if __name__=='__main__': main()
