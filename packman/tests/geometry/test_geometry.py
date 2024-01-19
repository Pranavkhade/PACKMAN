from networkx.utils.misc import flatten
from ... import molecule
from ... import geometry
import unittest

import logging
from os import remove as rm

class Test_Geometry(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif')
    
    def test_Circumsphere(self):
        self.assertIsInstance( geometry.Circumsphere( [i for i in self.mol[0].get_atoms()][:4] )[0][0], float )
        self.assertIsInstance( geometry.Circumsphere( [i for i in self.mol[0].get_atoms()][:4] )[0][1], float )
        self.assertIsInstance( geometry.Circumsphere( [i for i in self.mol[0].get_atoms()][:4] )[0][2], float )
        self.assertIsInstance( geometry.Circumsphere( [i for i in self.mol[0].get_atoms()][:4] )[1], float )

    def test_AlphaShape(self):
        #Checked only one instance of the atom
        try:
            self.assertIsInstance( geometry.AlphaShape( [j for i in self.mol[0].get_backbone() for j in i], 4 )[0][0], molecule.Atom )
        except:
            self.assertIsInstance( geometry.AlphaShape( [j for i in self.mol[0].get_backbone() for j in i], 4 )[0][0][0], molecule.Atom )
        
    def tearDown(self):
        logging.info('Geometry Test Done.')

if(__name__=='__main__'):
    unittest.main()