from ... import molecule
from ...apps import DCI
import unittest

import logging
from os import remove as rm

class TestMolecule(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif')

    def test_DCI(self):
        obj = DCI(self.mol)
        
    def tearDown(self):
        rm('DCI_pymol_output.txt')
        logging.info('DCI Test Done.')

if(__name__=='__main__'):
    unittest.main()