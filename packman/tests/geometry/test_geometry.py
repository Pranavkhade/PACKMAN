from networkx.utils.misc import flatten
from ... import molecule
from ... import geometry
import unittest

import logging
import numpy

class Test_Geometry(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif')
    
    def test_Circumsphere(self):
        atoms = list(self.mol[0].get_atoms())[:4]
        center, radius = geometry.Circumsphere(atoms)
        expected_center = numpy.array([19.851322237851143, -22.995767842220495, -1.779101667088373])
        self.assertEqual(center.shape, (3,))
        self.assertTrue(numpy.allclose(center, expected_center, rtol=0, atol=1e-10))
        self.assertAlmostEqual(radius, 1.4243694807576892, places=10)

    def test_AlphaShape(self):
        backbone = [atom for residue_atoms in self.mol[0].get_backbone() for atom in residue_atoms]
        simplices, graph = geometry.AlphaShape(backbone, 4)
        self.assertEqual(len(backbone), 792)
        self.assertEqual(len(simplices), 2833)
        self.assertEqual(graph.number_of_nodes(), len(backbone))
        self.assertEqual(graph.number_of_edges(), 4494)
        self.assertTrue(all(len(simplex) == 4 for simplex in simplices))
        self.assertIsInstance(simplices[0][0], molecule.Atom)
        
    def tearDown(self):
        logging.info('Geometry Test Done.')

if(__name__=='__main__'):
    unittest.main()