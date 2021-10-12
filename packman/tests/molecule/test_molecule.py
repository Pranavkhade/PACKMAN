from ... import molecule
import unittest

import logging
from os import remove as rm

class TestMolecule(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif')

    def test_load_pdb(self):
        self.assertTrue( molecule.load_structure('packman/tests/data/1prw.pdb',ftype='pdb') )
    
    def test_load_cif(self):
        self.assertTrue( molecule.load_structure('packman/tests/data/4hla.cif',ftype='cif') )
    
    def test_Protein(self):
        #Basic
        self.assertIsInstance( self.mol, molecule.Protein )

        #Write Methods
        self.assertTrue( self.mol.write_structure( 'test.cif',ftype='cif') )
        self.assertTrue( self.mol.write_structure( 'test.pdb',ftype='pdb') )
        with self.assertRaises(Exception) as file_format_problem:
            self.assertTrue( self.mol.write_structure( 'test',ftype='random') )
        self.assertTrue( 'Please provide appropriate "ftype" argument. (cif/pdb).' in str(file_format_problem.exception) )
        rm( 'test.cif' )
        rm( 'test.pdb' )

        self.assertNotEqual( len([i for i in self.mol]) , 0 )
    
    def test_Model(self):
        #Basic
        self.assertIsInstance( self.mol[0], molecule.Model )
        self.assertIsInstance( self.mol[0].get_parent(), molecule.Protein )
        self.assertIsNotNone( self.mol[0].get_id() )

        #Get Methods
        self.assertNotEqual( len( [i for i in self.mol[0].get_chains()] )  , 0 )
        self.assertNotEqual( len( [i for i in self.mol[0].get_residues()] ), 0 )
        self.assertNotEqual( len( [i for i in self.mol[0].get_atoms()] )   , 0 )
        self.assertNotEqual( len( [i for i in self.mol[0].get_calpha()] )  , 0 )
        self.assertNotEqual( len( [i for i in self.mol[0].get_backbone()] ), 0 )
    
    def test_Chain(self):
        #Basic
        self.assertIsInstance( self.mol[0]['A'], molecule.Chain )
        self.assertIsInstance( self.mol[0]['A'].get_parent(), molecule.Model )
        self.assertIsNotNone( self.mol[0].get_id() )

        #Get Methods
        self.assertNotEqual( len( [i for i in self.mol[0].get_residues()] ), 0 )
        self.assertNotEqual( len( [i for i in self.mol[0].get_atoms()] )   , 0 )
        self.assertNotEqual( len( [i for i in self.mol[0].get_calpha()] )  , 0 )
        self.assertNotEqual( len( [i for i in self.mol[0].get_backbone()] ), 0 )
    
    def test_Residue(self):
        #Basic
        Residues = [ i for i in self.mol[0]['A'].get_residues() ]
        self.assertIsInstance( Residues[0], molecule.Residue )
        self.assertIsInstance( Residues[0].get_parent(), molecule.Chain )
        self.assertIsNotNone( Residues[0].get_id() )

        #Get Methods
        self.assertNotEqual( len( [i for i in Residues[0].get_atoms()] )   , 0 )
        self.assertNotEqual( len( [i for i in Residues[0].get_backbone()] ), 0 )

        self.assertIsNotNone( Residues[0].get_calpha() )
        self.assertIsNotNone( Residues[0].get_tip() )
        self.assertIsNotNone( Residues[0].get_centerofgravity() )
    
    def test_Atom(self):
        #Basic
        Residues = [ i for i in self.mol[0]['A'].get_residues() ]
        self.assertIsInstance( Residues[0], molecule.Residue )
        self.assertIsInstance( Residues[0].get_parent(), molecule.Chain )
        self.assertIsNotNone( Residues[0].get_id() )

        #Get Methods
        self.assertNotEqual( len( [i for i in Residues[0].get_atoms()] )   , 0 )
        self.assertNotEqual( len( [i for i in Residues[0].get_backbone()] ), 0 )

        self.assertIsNotNone( Residues[0].get_calpha() )
        self.assertIsNotNone( Residues[0].get_tip() )
        self.assertIsNotNone( Residues[0].get_centerofgravity() )
    
    def test_Bond(self):
        self.assertEqual( len([i.get_id() for i in self.mol[0].get_bonds()]) , 1542 )
    
    def tearDown(self):
        logging.info('Molecule Test Done.')

if(__name__=='__main__'):
    unittest.main()