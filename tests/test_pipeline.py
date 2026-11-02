"""Known-answer tests. Expectations are interval/base arithmetic, not snapshots."""
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / 'tests' / 'fixtures'


class Pipeline(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='iwgc-test-')
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.fa = self.work / 'genome.fa'
        shutil.copyfile(FIX / 'genome.fa', self.fa)
        self.out = self.work / 'tracks'
        self.out.mkdir()

    def run_tracks(self, *args, ok=True, env=None):
        p = subprocess.run(['bash', str(ROOT / 'iwgc_circos_tracks.sh'), str(self.fa),
                            '-out', str(self.out), *map(str, args)],
                           text=True, capture_output=True, timeout=15, env=env)
        if ok:
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        else:
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        return p

    def file(self, name, text):
        p = self.work / name
        p.write_text(text)
        return p

    def rows(self, suffix):
        return [x.split() for x in (self.out / ('genome.fa_' + suffix)).read_text().splitlines()]

    def values(self, suffix):
        return [float(x[3]) for x in self.rows(suffix)]

    def close(self, actual, expected):
        self.assertEqual(len(actual), len(expected))
        for a, e in zip(actual, expected):
            self.assertAlmostEqual(a, e, delta=3e-6)

    def gff(self, start, end, chrom='Chr1', kind='gene', attrs='ID=a'):
        return f'{chrom}\ts\t{kind}\t{start}\t{end}\t.\t+\t.\t{attrs}\n'

    def test_karyotype_only(self):
        self.run_tracks()
        rows = self.rows('karyotype.circos')
        self.assertEqual([r[2:6] for r in rows], [['Chr1', '1', '0', '12'], ['Chr2', '2', '0', '12']])
        self.assertEqual(len(list(self.out.iterdir())), 1)

    def test_gene_exact_feature_coordinates_and_sqrt(self):
        self.run_tracks('-gene', FIX/'genes.gff3', '-window', 4, '-keep-temp')
        self.close(self.values('gene_coverage.circos'), [math.sqrt(.75), 1, math.sqrt(.5), 0, 0, 0])
        self.assertEqual((self.out/'genome.fa_genes_coords.bed').read_text(), 'Chr1\t1\t4\nChr1\t4\t8\nChr1\t9\t11\n')

    def test_overlapping_repeat_union_not_sum(self):
        self.run_tracks('-repeat', FIX/'repeats.gff3', '-window', 4)
        self.close(self.values('repeat_coverage.circos'), [math.sqrt(.75), 1, math.sqrt(.5), 0, 0, 0])

    def test_independent_bed_coverage(self):
        p = subprocess.run(['bedtools','coverage','-a',str(self.file('windows.bed','Chr1\t0\t4\nChr1\t4\t8\nChr1\t8\t12\n')),
                            '-b',str(FIX/'features.bed')],capture_output=True,text=True,check=True)
        self.close([float(r.split()[6]) for r in p.stdout.splitlines()], [.75,1,.5])

    def test_minmax_uses_all_windows(self):
        g = self.file('all.gff', self.gff(1,1)+self.gff(5,6)+self.gff(9,12)+self.gff(1,1,'Chr2')+self.gff(5,5,'Chr2')+self.gff(9,9,'Chr2'))
        self.run_tracks('-gene',g,'-window',4)
        self.close(self.values('gene_coverage.circos'),[0, (math.sqrt(.5)-.5)/.5, 1, 0,0,0])

    def test_near_identical_coverage_preserves_variation(self):
        self.fa.write_text('>Chr1\n'+'A'*20000000+'\n')
        g=self.file('precise.gff',self.gff(1,9999999)+self.gff(10000001,20000000))
        self.run_tracks('-gene',g,'-window',10000000)
        self.close(self.values('gene_coverage.circos'),[0,1])

    def test_gc_counts_preserve_small_differences(self):
        self.fa.write_text('>Chr1\n'+'G'*999999+'A'+'G'*1000000+'\n')
        self.run_tracks('-gc','-window',1000000)
        self.close(self.values('gc.circos'),[0,1])

    def test_constant_gene_and_gc_zero(self):
        self.fa.write_text('>Chr1\nACGTACGT\n')
        g = self.file('full.gff',self.gff(1,8))
        self.run_tracks('-gene',g,'-gc','-window',4)
        self.close(self.values('gene_coverage.circos'),[0,0])
        self.close(self.values('gc.circos'),[0,0])

    def test_single_base_and_window_boundary(self):
        self.run_tracks('-gene',self.file('single.gff',self.gff(4,4)+self.gff(5,5)),'-window',4)
        self.close(self.values('gene_coverage.circos'),[1,1,0,0,0,0])

    def test_gff_comments_hash_and_embedded_fasta(self):
        g=self.file('tail.gff','##gff-version 3\r\n'+self.gff(1,4,attrs='ID=a;Note=#gene')+'##FASTA\n>Chr1\nACGT\n')
        self.run_tracks('-gene',g,'-window',4)
        self.close(self.values('gene_coverage.circos'),[1,0,0,0,0,0])

    def test_no_gene_features_is_zero_coverage(self):
        self.run_tracks('-gene',self.file('empty.gff','# no features\n'),'-window',4)
        self.close(self.values('gene_coverage.circos'),[0]*6)

    def test_gc_ambiguities_denominator_and_case(self):
        self.fa.write_text('>Chr1\ngcNNATNN\n>Chr2\nGCGCAAAA\n')
        self.run_tracks('-gc','-window',4,'-keep-temp')
        # N contributes to window length: raw .5,0,1,0; min-max is unchanged.
        self.close(self.values('gc.circos'),[.5,0,1,0])
        raw=(self.out/'gc_content.bed').read_text().splitlines()[1:]
        self.close([float(x.split()[4]) for x in raw],[.5,0,1,0])

    def test_fixed_partial_tail(self):
        self.fa.write_text('>Chr1\nACGTACGTAA\n')
        self.run_tracks('-gc','-window',4,'-keep-temp')
        self.assertEqual((self.out/'genome.fa_windows.bed').read_text(),'Chr1\t0\t4\nChr1\t4\t8\nChr1\t8\t10\n')
        self.close(self.values('gc.circos'),[1,1,0])

    def test_sliding_tail_cutoff(self):
        self.fa.write_text('>Chr1\nACGTACGTA\n')
        self.run_tracks('-gc','-window',4,'-sliding','-step','0.5','-keep-temp')
        self.assertEqual([r[:3] for r in self.rows('gc.circos')],[['Chr1','0','4'],['Chr1','2','6'],['Chr1','4','8'],['Chr1','6','9']])

    def test_sliding_zero_step_rejected(self):
        for step in ('0','0.01','1.1','junk'):
            with self.subTest(step=step): self.run_tracks('-gc','-window',4,'-sliding','-step',step,ok=False)

    def test_too_short_sliding_has_no_windows_error(self):
        self.run_tracks('-gc','-window',100,'-sliding',ok=False)

    def test_intact_classes_and_overlap(self):
        self.run_tracks('-intact',FIX/'repeats.gff3','-window',4)
        self.close(self.values('Gypsy_intactTE_coverage.circos'),[1,0,0,0,0,0])
        self.close(self.values('Copia_intactTE_coverage.circos'),[0,1,0,0,0,0])
        self.close(self.values('DTA_intactTE_coverage.circos'),[0,0,1,0,0,0])

    def test_stale_intact_intermediate_not_read(self):
        (self.out/'genome.fa_intactTE_Gypsy.gff').write_text(self.gff(1,12))
        self.run_tracks('-intact',self.file('empty.gff','#empty\n'),'-window',4)
        self.assertFalse((self.out/'genome.fa_Gypsy_intactTE_coverage.circos').exists())

    def test_intact_path_traversal_rejected(self):
        g=self.file('bad.gff',self.gff(1,4,attrs='Classification=LTR/../../escape'))
        self.run_tracks('-intact',g,'-window',4,ok=False)

    def test_intact_missing_class_rejected(self):
        self.run_tracks('-intact',self.file('bad.gff',self.gff(1,4)),'-window',4,ok=False)

    def test_ltr_zero_age_missing_and_reverse(self):
        self.run_tracks('-ltr-dating',FIX/'ltr.pass.list','-window',4,'-keep-temp')
        self.close(self.values('LTRdating_coverage.circos'),[2,1.5,1,0,0,0])
        self.assertEqual([r.split()[3] for r in (self.out/'genome.fa_LTR_age.bed').read_text().splitlines()],['0','10','20','.','.','.'])

    def test_ltr_mean_not_length_weighted(self):
        def row(loc,age): return loc+'\t'+('\t'.join(['.']*10))+f'\t{age}\n'
        f=self.file('ages',row('Chr1:1..4',0)+row('Chr1:4..4',20)+row('Chr1:5..8',20))
        self.run_tracks('-ltr-dating',f,'-window',4,'-keep-temp')
        raw=(self.out/'genome.fa_LTR_age.bed').read_text().splitlines()
        self.close([float(x.split()[3]) for x in raw[:2]],[10,20])
        self.close(self.values('LTRdating_coverage.circos'),[2,1,0,0,0,0])

    def test_ltr_identical_zero_ages(self):
        f=self.file('one', 'Chr1:1..4\t'+('\t'.join(['.']*10))+'\t0\n')
        self.run_tracks('-ltr-dating',f,'-window',4)
        self.close(self.values('LTRdating_coverage.circos'),[1.5,0,0,0,0,0])

    def test_empty_ltr_all_missing(self):
        self.run_tracks('-ltr-dating',self.file('empty',''),'-window',4)
        self.close(self.values('LTRdating_coverage.circos'),[0]*6)

    def test_ltr_bad_age_and_schema(self):
        for text in ('bad\n',*[f'Chr1:1..4\t'+('\t'.join(['.']*10))+f'\t{x}\n' for x in ('nan','-1','bad')]):
            with self.subTest(text=text): self.run_tracks('-ltr-dating',self.file('bad',text),'-window',4,ok=False)

    def test_mummer_reverse_closed_to_bed(self):
        self.run_tracks('-links',FIX/'mummer.coords')
        self.assertEqual(self.rows('links.circos'),[['Chr1','0','4','Chr2','6','10']])

    def test_mummer_13_columns(self):
        f=self.file('13','1 4 10 7 4 4 99 12 12 33 33 Chr1 Chr2\n')
        self.run_tracks('-links',f)
        self.assertEqual(self.rows('links.circos'),[['Chr1','0','4','Chr2','6','10']])

    def test_sixcol_reverse_normalized(self):
        self.run_tracks('-links',self.file('rev','Chr1 4 0 Chr2 10 6\n'))
        self.assertEqual(self.rows('links.circos'),[['Chr1','0','4','Chr2','6','10']])

    def test_top_longest_then_order(self):
        self.run_tracks('-links',FIX/'links.coords','-top-n',2)
        self.assertEqual([r[1] for r in self.rows('links.circos')],['3','7'])
        self.run_tracks('-links',FIX/'links.coords','-top-n',2,'-link-order','dsc')
        self.assertEqual([r[1] for r in self.rows('links.circos')],['7','3'])

    def test_min_total_length_inclusive(self):
        self.run_tracks('-links',FIX/'links.coords','-min-tlen',4)
        self.assertEqual([r[1] for r in self.rows('links.circos')],['3','7'])

    def test_merge_both_spans_and_second_lower_bound(self):
        f=self.file('merge','Chr1 0 4 Chr2 7 11\nChr1 3 6 Chr2 6 9\n')
        self.run_tracks('-links',f)
        self.assertEqual(self.rows('links.circos'),[['Chr1','0','6','Chr2','6','11']])
        f.write_text('Chr1 0 4 Chr2 7 11\nChr1 3 6 Chr2 0 3\n')
        self.run_tracks('-links',f)
        self.assertEqual(len(self.rows('links.circos')),2)

    def test_adjacent_links_merge(self):
        self.run_tracks('-links',self.file('adj','Chr1 0 4 Chr2 0 4\nChr1 4 8 Chr2 4 8\n'))
        self.assertEqual(self.rows('links.circos'),[['Chr1','0','8','Chr2','0','8']])

    def test_filter_before_link_topn(self):
        f=self.file('filter','Chr1 0 1 Chr1 2 3\nChr1 0 12 Chr2 0 12\n')
        self.run_tracks('-links',f,'-top-n',1,'-filter-chrs','^chr1$')
        self.assertEqual(self.rows('links.circos'),[['Chr1','0','1','Chr1','2','3']])

    def test_filter_case_insensitive_and_no_match(self):
        self.run_tracks('-gc','-window',4,'-filter-chrs','^cHr2$')
        self.assertEqual([r[0] for r in self.rows('gc.circos')],['Chr2']*3)
        self.run_tracks('-filter-chrs','absent',ok=False)

    def test_empty_links_is_valid(self):
        self.run_tracks('-links',self.file('empty','# no alignments\n'))
        self.assertEqual(self.rows('links.circos'),[])

    def test_malformed_links_not_silently_skipped(self):
        for text in ('garbage\n','Chr1 0 1 Chr2 0 1\nbad\n','Chr1 -1 4 Chr2 0 4\n','Chr1 0 99 Chr2 0 1\n','missing 0 1 Chr2 0 1\n'):
            with self.subTest(text=text): self.run_tracks('-links',self.file('bad',text),ok=False)

    def test_malformed_gff_rejected(self):
        for text in ('not a GFF\n',self.gff(0,3),self.gff(4,2),self.gff(1,13),self.gff('x',4),self.gff(1,4,'absent')):
            with self.subTest(text=text): self.run_tracks('-gene',self.file('bad.gff',text),'-window',4,ok=False)

    def test_missing_and_invalid_option_values(self):
        for args in (('-gene',),('-window',),('-bogus',),('-window','0'),('-window','x'),('-top-n','1;exit'),('-link-order','desc')):
            with self.subTest(args=args): self.run_tracks(*args,ok=False)

    def test_missing_empty_and_stale_fasta(self):
        self.fa.unlink(); self.run_tracks(ok=False)
        self.fa.write_text(''); self.run_tracks(ok=False)
        self.fa.write_text('>Chr1\nACGT\n')
        self.run_tracks()
        os.utime(self.fa,(2000000000,2000000000))
        self.run_tracks(ok=False)

    def test_duplicate_index_ids(self):
        p=Path(str(self.fa)+'.fai'); p.write_text('Chr1\t12\t6\t12\t13\nChr1\t12\t6\t12\t13\n')
        self.run_tracks(ok=False)

    def test_duplicate_fasta_records_rejected(self):
        self.fa.write_text('>Chr1\nACGT\n>Chr1\nGGGG\n')
        self.run_tracks(ok=False)

    @unittest.skipUnless(shutil.which('quartet.py') and shutil.which('tidk'), 'optional quarTeT/tidk not installed')
    def test_real_plant_telomere_detector(self):
        # Two 30 kb chromosomes, >50 plant motifs at each end.
        sequence='TTTAGGG'*100+'ACGT'*7150+'CCCTAAA'*100
        self.fa.write_text('>Chr1\n'+sequence+'\n>Chr2\n'+sequence+'\n')
        self.run_tracks('-telomere','-keep-temp')
        bands=[r for r in self.rows('karyotype.circos') if r[0]=='band']
        self.assertEqual(len(bands),4)
        self.assertEqual([r[4:6] for r in bands],[['0','300'],['29700','30000']]*2)

    @unittest.skipUnless(shutil.which('quartet.py') and shutil.which('tidk'), 'optional quarTeT/tidk not installed')
    def test_single_chromosome_detector_limitation_is_failure(self):
        # The tested tidk 0.2.0 fails to nominate a motif for this one-record input.
        # This checks failure propagation, not a correct biological absence call.
        version=subprocess.run(['tidk','--version'],capture_output=True,text=True).stdout.strip()
        if version != 'tidk 0.2.0': self.skipTest('limitation is specific to the audited tidk version')
        self.fa.write_text('>Chr1\n'+'TTTAGGG'*100+'ACGT'*7150+'CCCTAAA'*100+'\n')
        p=self.run_tracks('-telomere',ok=False)
        self.assertIn('no telomere report',p.stderr)

    def test_paths_with_spaces(self):
        self.fa=self.fa.rename(self.work/'genome spaced.fa')
        p=self.run_tracks('-gc','-window',4)
        self.assertTrue((self.out/'genome spaced.fa_gc.circos').is_file(),p.stderr)

    def test_external_failure_is_nonzero(self):
        fake=self.work/'bin'; fake.mkdir()
        b=fake/'bedtools'; b.write_text('#!/bin/sh\necho deliberate failure >&2\nexit 23\n'); b.chmod(0o755)
        env=dict(os.environ,PATH=str(fake)+os.pathsep+os.environ['PATH'])
        p=self.run_tracks('-gc',env=env,ok=False)
        self.assertIn('deliberate failure',p.stderr)

    def test_telomere_uses_fixed_dependency_arguments(self):
        self.fa=self.fa.rename(self.work/'genome$(touch${IFS}MARKER).fa')
        fake=self.work/'bin'; fake.mkdir()
        b=fake/'quartet.py'; b.write_text('#!/bin/sh\n[ "$3" = input.fa ] && [ "$9" = iwgc ] || exit 90\nprintf "Chr1 12 both . . .\\nChr2 12 no . . .\\n" > iwgc.telo.info\n'); b.chmod(0o755)
        self.run_tracks('-telomere',env=dict(os.environ,PATH=str(fake)+os.pathsep+os.environ['PATH']))
        self.assertFalse(list(self.work.rglob('MARKER')))

    def test_dependency_zero_status_without_data_is_failure(self):
        real=shutil.which('bedtools')
        fake=self.work/'bin'; fake.mkdir()
        b=fake/'bedtools'; b.write_text(f'#!/bin/sh\nif [ "$1" = nuc ]; then echo "# header only"; exit 0; fi\nexec "{real}" "$@"\n'); b.chmod(0o755)
        p=self.run_tracks('-gc','-window',4,env=dict(os.environ,PATH=str(fake)+os.pathsep+os.environ['PATH']),ok=False)
        self.assertIn('Incomplete',p.stderr)

    def test_palette_cycles_after_sixty_contigs(self):
        self.fa.write_text(''.join(f'>contig{i:03}\nACGT\n' for i in range(61)))
        self.run_tracks()
        rows=self.rows('karyotype.circos')
        self.assertEqual(len(rows),61)
        self.assertTrue(all(len(r)==7 and r[6] for r in rows))
        self.assertEqual(rows[0][6],rows[60][6])

    def test_help_and_no_args(self):
        for script in ('iwgc_circos_tracks.sh','create_configs.sh','add_legend.sh'):
            p=subprocess.run(['bash',str(ROOT/script),'--help'],capture_output=True,text=True,timeout=5)
            self.assertEqual(p.returncode,0,p.stderr)
            self.assertIn('Usage:',p.stdout)
        p=subprocess.run(['bash',str(ROOT/'iwgc_circos_tracks.sh')],capture_output=True,text=True,timeout=5)
        self.assertNotEqual(p.returncode,0)

    def test_paf_helper_with_fixed_aligner_output(self):
        fake=self.work/'bin'; fake.mkdir()
        b=fake/'minimap2'; b.write_text(f'#!/bin/sh\ncat "{FIX / "alignment.paf"}"\n'); b.chmod(0o755)
        prefix=self.work/'alignment'
        p=subprocess.run(['bash',str(ROOT/'preprocess/minimap2_self_alignment.sh'),'-f',str(self.fa),'-p',str(prefix)],
                         env=dict(os.environ,PATH=str(fake)+os.pathsep+os.environ['PATH']),capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertEqual(prefix.with_suffix('.coords').read_text(),'Chr1\t0\t4\tChr2\t6\t10\n')

    def test_telomere_marker_bounds_with_quartet_contract(self):
        fake=self.work/'bin'; fake.mkdir()
        b=fake/'quartet.py'; b.write_text('#!/bin/sh\nwhile [ "$#" -gt 0 ]; do if [ "$1" = -p ]; then shift; prefix=$1; fi; shift; done\nprintf "Chr1 12 both . . .\\nChr2 12 left . . .\\n" > "$prefix.telo.info"\n'); b.chmod(0o755)
        self.run_tracks('-telomere','-ts','1',env=dict(os.environ,PATH=str(fake)+os.pathsep+os.environ['PATH']))
        bands=[r for r in self.rows('karyotype.circos') if r[0]=='band']
        self.assertEqual(len(bands),3)
        for r in bands: self.assertEqual(r[4:6],['0','12'])


if __name__ == '__main__': unittest.main()
