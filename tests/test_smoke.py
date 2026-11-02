"""Regression for Circos's observed zero-status/missing-module failure."""
from pathlib import Path
import sys
import tempfile
import unittest

from smoke import run


class SmokeDependencyFailure(unittest.TestCase):
    def test_zero_status_without_outputs_preserves_dependency_diagnostic(self):
        with tempfile.TemporaryDirectory(prefix='iwgc-smoke-error-') as directory:
            log=[]
            # Actual subprocess: reproduce the CI dependency's misleading status.
            command=[sys.executable,'-c',"print('missing List::MoreUtils')"]
            with self.assertRaisesRegex(RuntimeError,'(?s)exited 0 without required output:.*missing List::MoreUtils'):
                run(directory,command,log,required_outputs=(Path(directory)/'figure.svg',))
            self.assertEqual(log[0]['rc'],0)
            self.assertIn('missing List::MoreUtils',log[0]['stdout'])


if __name__=='__main__': unittest.main()
