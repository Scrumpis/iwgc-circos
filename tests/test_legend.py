from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]

class LegendErrors(unittest.TestCase):
    def run_bad(self,*args):
        p=subprocess.run(['bash',str(ROOT/'add_legend.sh'),*map(str,args)],capture_output=True,text=True,timeout=10)
        self.assertNotEqual(p.returncode,0,p.stdout+p.stderr)
        return p

    def test_missing_value_and_unknown(self):
        for args in (('-out',),('--dx-frac','junk'),('-unknown',)):
            with self.subTest(args=args): self.run_bad(*args)

    def test_missing_or_malformed_svg(self):
        with tempfile.TemporaryDirectory(prefix='iwgc-legend-') as w:
            f=Path(w)/'bad.svg'
            self.run_bad('-svg',f)
            f.write_text('<svg not closed')
            self.run_bad('-svg',f)
            f.write_text('<svg xmlns="http://www.w3.org/2000/svg"></svg>')
            self.run_bad('-svg',f)

    def test_refuses_already_labeled_or_same_output(self):
        with tempfile.TemporaryDirectory(prefix='iwgc-legend-') as w:
            f=Path(w)/'labeled.svg'; f.write_text('<svg xmlns="http://www.w3.org/2000/svg"><g id="iwgc-gap-labels"/></svg>')
            p=self.run_bad('-svg',f,'-out',Path(w)/'out.svg'); self.assertIn('already',p.stderr)
            self.run_bad('-svg',f,'-out',f)

if __name__=='__main__': unittest.main()
