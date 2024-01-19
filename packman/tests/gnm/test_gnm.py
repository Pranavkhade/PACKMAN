from ... import molecule
from ...gnm import GNM
from os import remove as rm

import unittest
import logging
import numpy

class Test_GNM(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif')
        self.calpha = [i for i in self.mol[0]['A'].get_calpha() if i is not None]

    def test_GNM(self):
        self.Model=GNM(self.calpha,gamma=1,dr=7.3,power=1)
        self.assertTrue( self.Model.calculate_kirchhoff() )
        self.assertTrue( self.Model.calculate_decomposition() )

        self.assertIsNotNone( self.Model.get_eigenvalues() )
        self.assertIsNotNone( self.Model.get_eigenvectors() )

        expected_gnm_ev = [1.12895519e-15,2.49413543e-01,2.84146316e-01,4.32358642e-01,8.09232620e-01,8.96564677e-01,1.19722332e+00,1.30674130e+00,1.50117123e+00,1.87918334e+00,2.14839686e+00,2.36525651e+00,2.56994595e+00,2.87622645e+00,3.00401321e+00,3.34689286e+00,3.57275469e+00,3.93783245e+00,4.13784251e+00,4.28656258e+00,4.62316006e+00,4.71178990e+00,4.72274835e+00,4.77412491e+00,4.97971268e+00,5.12367322e+00,5.27519690e+00,5.33975077e+00,5.52526809e+00,5.67956033e+00,5.94606030e+00,5.98902622e+00,6.05394435e+00,6.20546701e+00,6.26936239e+00,6.50364633e+00,6.56224734e+00,6.70144860e+00,6.83899518e+00,6.93572063e+00,7.09919839e+00,7.23951211e+00,7.39065055e+00,7.41404445e+00,7.51028761e+00,7.59095332e+00,7.76515669e+00,7.80007298e+00,7.93583219e+00,8.14110715e+00,8.25953973e+00,8.31984590e+00,8.48294596e+00,8.50431418e+00,8.59466891e+00,8.71431965e+00,8.76821347e+00,8.78712171e+00,8.95597310e+00,9.10975255e+00,9.22031020e+00,9.33254577e+00,9.34079860e+00,9.42583082e+00,9.57883523e+00,9.76040483e+00,9.81839847e+00,9.86826133e+00,1.00317865e+01,1.02142485e+01,1.03721970e+01,1.05352368e+01,1.05955304e+01,1.07620216e+01,1.09612822e+01,1.12621794e+01,1.13802414e+01,1.15446499e+01,1.16513029e+01,1.17125632e+01,1.19510066e+01,1.20444117e+01,1.21134754e+01,1.23445128e+01,1.24573379e+01,1.25274086e+01,1.27247934e+01,1.27458249e+01,1.28568826e+01,1.31053796e+01,1.32193123e+01,1.33110903e+01,1.35277914e+01,1.38394477e+01,1.40065788e+01,1.43322515e+01,1.46453305e+01,1.53969481e+01,1.55354177e+01]
        expected_gnm_ev = numpy.array(expected_gnm_ev)

        # Are the eigenvalues almost same? everything else tends to be equal if so.
        self.assertTrue( numpy.allclose(self.Model.get_eigenvalues(), expected_gnm_ev, rtol=1e-7, atol=1e-7) )
        
        self.assertTrue( self.Model.calculate_crosscorrelation() )
        self.assertTrue( self.Model.calculate_fluctuations() )

        self.assertIsNotNone( self.Model.get_crosscorrelation() )
        self.assertIsNotNone( self.Model.get_fluctuations )
        self.assertIsNotNone( self.Model.get_pseudoinverse )


    def tearDown(self):
        logging.info('GNM test Done.')

if(__name__=='__main__'):
    unittest.main()