from ... import molecule
from ...entropy import PackingEntropy
import unittest
import numpy

import logging

class Test_Entropy(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif')

    def test_PackingEntropy(self):
        obj = PackingEntropy(self.mol[0].get_atoms())

        self.assertIsInstance(obj, PackingEntropy)
        surface_points = numpy.asarray(obj.get_surafacepoints())
        self.assertEqual(surface_points.ndim, 2)
        self.assertEqual(surface_points.shape[1], 3)
        self.assertGreater(len(surface_points), 0)
        self.assertTrue(numpy.isfinite(surface_points).all())
        self.assertIsNotNone( obj.get_total_chain_entropy('A') )
        self.assertIsNotNone( obj.get_total_entropy() )

        self.assertAlmostEqual(obj.get_total_entropy(), 3087.1469554213118, places=7)
        self.assertAlmostEqual(obj.get_total_chain_entropy('A'), 1546.7759202113048, places=7)
        self.assertAlmostEqual(obj.get_total_chain_entropy('B'), 1540.3710352100072, places=7)
        self.assertAlmostEqual(
            obj.get_total_chain_entropy('A') + obj.get_total_chain_entropy('B'),
            obj.get_total_entropy(),
            places=7,
        )

        self.assertIsInstance( self.mol[0].get_entropy('PackingEntropy') , float )
        self.assertIsInstance( self.mol[0]['A'].get_entropy('PackingEntropy') , float )
        self.assertIsInstance( [i for i in self.mol[0]['A'].get_residues()][0].get_entropy('PackingEntropy') , float )
        self.assertNotEqual( self.mol[0]['A'].get_entropy('PackingEntropy') , self.mol[0]['B'].get_entropy('PackingEntropy') )

    def tearDown(self):
        logging.info('Entropy Test Done.')

if(__name__=='__main__'):
    unittest.main()