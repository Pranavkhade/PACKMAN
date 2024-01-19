from ... import molecule, utilities, geometry
import unittest
import numpy

import logging
from os import remove as rm

class Test_Utilities(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif')
    
    def test_superimpose(self):
        # Double check later.
        expected_a1 = numpy.array( [[-4.99437660e-01, -8.66349433e-01, -8.26228195e-04], [-8.66349792e-01,  4.99437142e-01, 7.60871525e-04], [-2.46531567e-04,  1.09581052e-03, -9.99999369e-01]] )
        expected_b1 = numpy.array( [[-0.02974615, -0.02141248,  0.65529778]] )

        expected_a2 = numpy.array( [[ 9.99999994e-01,  6.44035795e-05,  8.82404191e-05], [-6.44835190e-05,  9.99999587e-01,  9.06225051e-04], [-8.81820186e-05, -9.06230736e-04,  9.99999585e-01]] )
        expected_b2 = numpy.array( [[ 6.31828525e-03,  9.14265397e-05, -4.94164731e-03]] )

        a1, b1 = utilities.superimporse(self.mol[0]['A'], self.mol[0]['B'], use='calpha')
        self.assertTrue( numpy.allclose( a1, expected_a1, rtol=1e-14 ) )
        self.assertTrue( numpy.allclose( b1, expected_b1, rtol=1e-14 ) )

        # since superimpose change the coordinates, keep it mind for future debugging in case that is changed.
        a2, b2 = utilities.superimporse(self.mol[0]['A'], self.mol[0]['B'], use='backbone')
        self.assertTrue( numpy.allclose( a2, expected_a2, rtol=1e-14 ) )
        self.assertTrue( numpy.allclose( b2, expected_b2, rtol=1e-14 ) )
    
    def test_RMSD(self):
        self.assertTrue( (utilities.RMSD(self.mol[0]['A'], self.mol[0]['B'], use='calpha') - 0.10400695339185048) < 0.001 )
        self.assertTrue( (utilities.RMSD(self.mol[0]['A'], self.mol[0]['B'], use='backbone')- 0.11195505579911495) < 0.001  )

    def test_WriteOBJ(self):
        AS = geometry.AlphaShape([i for i in self.mol[0].get_atoms()], float('inf'))
        self.assertTrue( utilities.WriteOBJ( self.mol[0].get_atoms(), AS[0], open('test.obj','wb') ) )
        rm('test.obj')

    def tearDown(self):
        logging.info('Utilities Test Done.')

if(__name__=='__main__'):
    unittest.main()

