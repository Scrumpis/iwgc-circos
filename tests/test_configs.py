"""Check generated layout invariants instead of matching a whole config snapshot."""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
SUFFIXES=['gene_coverage','repeat_coverage','LTRdating_coverage',
          *[c+'_intactTE_coverage' for c in ['Gypsy','Copia','unknown','DTA','DTC','DTH','DTM','DTT','Helitron']], 'gc']


class Configs(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='iwgc-config-'); self.addCleanup(self.tmp.cleanup)
        self.w=Path(self.tmp.name)
        shutil.copyfile(ROOT/'iwgc_circos_template.config', self.w/'iwgc_circos_template.config')
        self.data=self.w/'iwgc_circos_data'; self.data.mkdir()
        (self.data/'g_karyotype.circos').write_text('chr - Chr1 1 0 120000 chr2\nchr - Chr2 2 0 120000 chr3\n')

    def run_config(self,*args,ok=True):
        p=subprocess.run(['bash',str(ROOT/'create_configs.sh'),*args],cwd=self.w,text=True,capture_output=True,timeout=15)
        self.assertEqual(p.returncode==0,ok,p.stdout+p.stderr)
        return p

    def text(self): return (self.w/'iwgc_circos/iwgc_circos.conf').read_text()

    def test_all_track_counts_geometry(self):
        for n in range(14):
            with self.subTest(n=n):
                for p in self.data.glob('*_coverage.circos'): p.unlink()
                for p in self.data.glob('*_gc.circos'): p.unlink()
                for suffix in SUFFIXES[:n]:
                    (self.data/f'g_{suffix}.circos').write_text('Chr1 0 40000 0.5\n')
                (self.data/'g_links.circos').write_text('Chr1 0 40000 Chr2 0 40000\n')
                self.run_config('-gap')
                s=self.text(); self.assertNotIn('${chrs}',s); self.assertNotIn('*_coverage',s)
                blocks=re.findall(r'<plot>\s*(.*?)</plot>',s,re.S)
                self.assertEqual(len(blocks),n)
                previous=1
                for block in blocks:
                    r0=float(re.search(r'\br0\s*=\s*([\d.]+)r',block)[1])
                    r1=float(re.search(r'\br1\s*=\s*([\d.]+)r',block)[1])
                    self.assertTrue(0<r0<r1<=previous)
                    previous=r0
                if 1<=n<=8: self.assertAlmostEqual(previous,.5,delta=.0002)
                if n:
                    link=float(re.search(r'\bradius\s*=\s*([\d.]+)r',s)[1])
                    self.assertAlmostEqual(link,previous-.005,delta=.0002)
                ideo=(self.w/'iwgc_circos/ideogram.conf').read_text()
                self.assertEqual(ideo.count('<pairwise '),1)
                self.assertIn('<pairwise Chr2 Chr1>',ideo)
                self.assertIn('spacing        = 1u',(self.w/'iwgc_circos/ticks.conf').read_text())

    def test_ambiguous_glob_is_error(self):
        (self.data/'second_karyotype.circos').write_text('chr - Chr9 9 0 100 chr2\n')
        p=self.run_config(ok=False); self.assertIn('Ambiguous',p.stderr)

    def test_empty_tracks_and_links_omitted(self):
        (self.data/'g_gene_coverage.circos').touch(); (self.data/'g_links.circos').touch()
        self.run_config(); s=self.text()
        self.assertNotIn('<plot>',s); self.assertNotIn('<links>',s)

    def test_ambiguous_track_is_error(self):
        for base in ('a','b'): (self.data/f'{base}_gc.circos').write_text('Chr1 0 40000 0\n')
        self.run_config(ok=False)

    def test_malformed_karyotype_rejected(self):
        for text in ('chr - Chr1 1 0 0 chr2\n','chr - Chr1 1 0 12 chr2\nchr - Chr1 1 0 12 chr3\n',
                     'chr - Chr1 1 0 12 chr2\nband absent a a 0 12 grey\n'):
            with self.subTest(text=text):
                (self.data/'g_karyotype.circos').write_text(text)
                self.run_config(ok=False)

    def test_valid_chr_id_resembling_band_is_retained(self):
        (self.data/'g_karyotype.circos').write_text('chr - Chr1 1 0 12 chr2\nchr - scaffold_T1 sc1 0 8 chr3\nband Chr1 Chr1_T1 Chr1_T1 0 1 grey\n')
        p=self.run_config('-gap')
        self.assertIn('chromosomes                 = Chr1;scaffold_T1',self.text())
        self.assertIn('2 chr',p.stdout)
        self.assertIn('<pairwise scaffold_T1 Chr1>',(self.w/'iwgc_circos/ideogram.conf').read_text())

    def test_custom_outdir_include_files(self):
        self.run_config('-outdir','custom')
        text=(self.w/'custom/iwgc_circos.conf').read_text()
        self.assertIn('<<include custom/ideogram.conf>>',text)
        self.assertIn('<<include custom/ticks.conf>>',text)
        for f in ('ideogram.conf','ideogram.position.conf','ideogram.label.conf','bands.conf','ticks.conf','iwgc_circos.conf'):
            self.assertTrue((self.w/'custom'/f).is_file())

    def test_missing_values_unknown_and_missing_template(self):
        for args in (('-outdir',),('-unknown',),('-template','absent')):
            with self.subTest(args=args): self.run_config(*args,ok=False)


if __name__=='__main__': unittest.main()
