from ... import molecule
from ...anm import ANM, hdANM
import unittest

import logging
from os import remove as rm

class TestMolecule(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif')
        self.calpha = [i for i in self.mol[0]['A'].get_calpha() if i is not None]

    def test_hdANM(self):
        self.Model=hdANM(self.calpha,dr=15,power=0,hng_file='packman/tests/data/1prw.hng')
        self.assertTrue( self.Model.calculate_hessian() )
        self.assertTrue( self.Model.calculate_decomposition() )

        self.assertIsNotNone( self.Model.get_eigenvalues() )
        self.assertIsNotNone( self.Model.get_eigenvectors() )

        self.assertTrue( self.Model.calculate_movie(6,scale=2,n=10) )
        self.assertTrue( self.Model.calculate_movie(6,scale=2,n=10, ftype='pdb') )
        rm('6.cif')
        rm('6.pdb')

        self.assertIsNotNone( self.Model.get_hessian_pseudoinverse() )
        self.assertIsNotNone( self.Model.get_RT_eigen_vectors() )

        #Add corsscorrelation after the publication
        #self.assertIsNotNone( [i for i in self.Model.get_crosscorrelation_matrix()] )
    
    def test_ANM_Compliance(self):
        self.ANM_MODEL = ANM(self.calpha)

        self.assertTrue( self.ANM_MODEL.calculate_hessian() )
        self.assertTrue( self.ANM_MODEL.calculate_decomposition() )
        self.assertTrue( self.ANM_MODEL.calculate_stiffness_compliance() )
        self.assertTrue( self.ANM_MODEL.calculate_fluctuations() )

        self.assertIsNotNone( self.ANM_MODEL.get_stiffness_map() )
        self.assertIsNotNone( self.ANM_MODEL.get_compliance_map() )

        self.assertIsNotNone( self.ANM_MODEL.get_fluctuations() )
        self.assertIsNotNone( self.ANM_MODEL.get_stiffness_profile() )
        self.assertIsNotNone( self.ANM_MODEL.get_compliance_profile() )

    def tearDown(self):
        logging.info('ANM, hdANM and Compliance Test Done.')

if(__name__=='__main__'):
    unittest.main()