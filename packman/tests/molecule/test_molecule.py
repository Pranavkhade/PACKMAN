from ... import molecule
import unittest
import numpy

import logging
from os import remove as rm

class Test_Molecule(unittest.TestCase):

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
        self.assertNotEqual( len( [i for i in self.mol.get_sequence()] ), 0 )
    
    def test_Model(self):
        #Basic
        self.assertIsInstance( self.mol[0], molecule.Model )
        self.assertIsInstance( self.mol[0].get_parent(), molecule.Protein )
        self.assertIsNotNone( self.mol[0].get_id() )

        #Get Methods
        self.assertEqual( len( [i for i in self.mol[0].get_chains()] )  , 2 )
        self.assertEqual( len( [i for i in self.mol[0].get_residues()] ), 198 )
        self.assertEqual( len( [i for i in self.mol[0].get_atoms()] )   , 1516 )
        self.assertEqual( len( [i for i in self.mol[0].get_calpha()] )  , 198 )
        self.assertEqual( len([j for i in self.mol[0].get_backbone() for j in i]), 792 )
        expected_seq = '>packman/tests/data/4hla.cif_1_A\nPQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF\n>packman/tests/data/4hla.cif_1_B\nPQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF\n'
        self.assertNotEqual( self.mol[0].get_sequence(), expected_seq )
    
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
        self.assertNotEqual( len( [i for i in self.mol[0]['A'].get_sequence()] ), 0 )
    
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
        Atoms = [ i for i in self.mol[0]['A'].get_atoms() ]
        self.assertIsInstance( Atoms[0], molecule.Atom )
        self.assertIsInstance( Atoms[0].get_parent(), molecule.Residue )
        self.assertIsNotNone( Atoms[0].get_id() )

        #Get Methods
        self.assertIsInstance( Atoms[0].get_location(), numpy.ndarray )


    def test_Bond(self):
        Bonds = [i for i in self.mol[0].get_bonds()]
        self.assertEqual( len(Bonds), 1542 )
        self.assertIsInstance( Bonds[0].get_id(), int)
        
        self.assertIsInstance(Bonds[0].get_atoms()[0], molecule.Atom)
        self.assertIsInstance(Bonds[0].get_atoms()[1], molecule.Atom)
    
    def tearDown(self):
        logging.info('Molecule Test Done.')

if(__name__=='__main__'):
    unittest.main()