from ... import molecule
from ...entropy import PackingEntropy
import unittest

import logging
from os import remove as rm

class Test_Entropy(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif')

    def test_PackingEntropy(self):
        obj = PackingEntropy(self.mol[0].get_atoms())

        self.assertIsInstance(obj, PackingEntropy)
        self.assertIsNotNone( obj.get_surafacepoints() )
        self.assertIsNotNone( obj.get_total_chain_entropy('A') )
        self.assertIsNotNone( obj.get_total_entropy() )

        # If total entropy is constant to the strictly set threshold, everything else most likely will be will be okay.
        self.assertTrue( (obj.get_total_entropy() - 3087.1469554213118) < 1 )

        self.assertIsInstance( self.mol[0].get_entropy('PackingEntropy') , float )
        self.assertIsInstance( self.mol[0]['A'].get_entropy('PackingEntropy') , float )
        self.assertIsInstance( [i for i in self.mol[0]['A'].get_residues()][0].get_entropy('PackingEntropy') , float )
        self.assertNotEqual( self.mol[0]['A'].get_entropy('PackingEntropy') , self.mol[0]['B'].get_entropy('PackingEntropy') )

    def tearDown(self):
        logging.info('Entropy Test Done.')

if(__name__=='__main__'):
    unittest.main()