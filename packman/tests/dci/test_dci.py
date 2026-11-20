from ... import molecule
from ...apps import DCI
import unittest
import numpy
from pathlib import Path
from tempfile import TemporaryDirectory


class Test_DCI(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif')

    def test_DCI(self):
        obj = DCI(self.mol, chain='A', n_com=4)
        crosscorrelation = obj.get_crosscorrelation()
        self.assertEqual(crosscorrelation.shape, (99, 99))
        self.assertTrue(numpy.isfinite(crosscorrelation).all())
        self.assertTrue(numpy.allclose(crosscorrelation, crosscorrelation.T, atol=1e-12))
        self.assertTrue(numpy.allclose(numpy.diag(crosscorrelation), numpy.ones(99), atol=1e-12))

        labels = obj.get_communities()[4]
        self.assertEqual(labels.shape, (99,))
        self.assertEqual(len(numpy.unique(labels)), 4)
        self.assertEqual(sorted(numpy.bincount(labels)[1:].tolist()), [15, 22, 31, 31])

        with TemporaryDirectory() as output_dir:
            output_path = Path(output_dir) / 'dci_pymol_output.txt'
            with output_path.open('w') as output_file:
                self.assertTrue(obj.calcualte_pymol_commands(output_file))
            output = output_path.read_text()

        self.assertIn('Number of communities - 4', output)
        self.assertIn('CH Score - user-key', output)
        for label in numpy.unique(labels):
            self.assertIn('cluster {}'.format(label), output)

if(__name__=='__main__'):
    unittest.main()