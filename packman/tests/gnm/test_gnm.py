from ... import molecule
from ...gnm import GNM
import unittest

import logging
from os import remove as rm

class TestMolecule(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif')
        self.calpha = [i for i in self.mol[0]['A'].get_calpha() if i is not None]

    def test_GNM(self):
        self.Model=GNM(self.calpha,gamma=1,dr=7.3,power=1)
        self.assertTrue( self.Model.calculate_kirchhoff() )
        self.assertTrue( self.Model.calculate_decomposition() )

        self.assertIsNotNone( self.Model.get_eigenvalues() )
        self.assertIsNotNone( self.Model.get_eigenvectors() )
        
        self.assertTrue( self.Model.calculate_crosscorrelation() )
        self.assertTrue( self.Model.calculate_fluctuations() )

        self.assertIsNotNone( self.Model.get_crosscorrelation() )
        self.assertIsNotNone( self.Model.get_fluctuations )
        self.assertIsNotNone( self.Model.get_pseudoinverse )


    def tearDown(self):
        logging.info('GNM test Done.')

if(__name__=='__main__'):
    unittest.main()