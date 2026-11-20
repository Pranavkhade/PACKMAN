from ... import molecule
import unittest
import numpy
from pathlib import Path
from collections import Counter

import logging
from os import remove as rm

DATA_DIR = Path(__file__).resolve().parents[1] / 'data'

# Fixture expectations are cross-checked against source records and Biopython.
PDB_EXPECTED_OUTPUTS = {
    '1prw.pdb': {
        'model_count': 1,
        'model_id': 0,
        'chain_ids': ('A',),
        'residue_count': 147,
        'atom_count': 1157,
        'calpha_count': 147,
        'backbone_count': 588,
        'sequence': 'ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGELTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK',
        'first_atom': (4, 'N', 'N', (57.976, 24.604, 1.622)),
    },
}

PDB_EDGE_CASE_OUTPUTS = {
    '1crn.pdb': {
        'model_ids': (0,),
        'chain_ids': ('A',),
        'residue_count': 46,
        'atom_count': 327,
    },
    '1ubq.pdb': {
        'model_ids': (0,),
        'chain_ids': ('A',),
        'residue_count': 76,
        'atom_count': 602,
        'hetmol_counts': {'HOH': 58},
    },
    '4hhb.pdb': {
        'model_ids': (0,),
        'chain_ids': ('A', 'B', 'C', 'D'),
        'residue_count': 574,
        'atom_count': 4384,
        'hetmol_counts': {'HEM': 4, 'HOH': 221, 'PO4': 2},
    },
    '1a3n.pdb': {
        'model_ids': (0,),
        'chain_ids': ('A', 'B', 'C', 'D'),
        'residue_count': 572,
        'atom_count': 4370,
        'het_atom_count': 627,
        'altloc_record_count': 8,
    },
    '1g03.pdb': {
        'model_ids': tuple(range(20)),
        'chain_ids': ('A',),
        'residue_count': 134,
        'atom_count': 2081,
    },
    '1bna.pdb': {
        'model_ids': (0,),
        'chain_ids': ('A', 'B'),
        'residue_count': 24,
        'atom_count': 486,
        'hetmol_counts': {'HOH': 80},
    },
}

CIF_EXPECTED_OUTPUTS = {
    '4hla.cif': {
        'model_count': 1,
        'model_id': 1,
        'chain_sequences': {
            'A': 'PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF',
            'B': 'PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF',
        },
        'residue_count': 198,
        'atom_count': 1516,
        'calpha_count': 198,
        'backbone_count': 792,
        'bond_count': 1542,
        'first_atom': (1, 'N', 'N', (27.247, -30.31, -3.159)),
    },
}

CIF_EDGE_CASE_OUTPUTS = {
    '1crn.cif': {
        'model_ids': (1,),
        'chain_ids': ('A',),
        'residue_count': 46,
        'atom_count': 327,
        'het_atom_count': 0,
        'hetmol_counts': {},
    },
    '1ubq.cif': {
        'model_ids': (1,),
        'chain_ids': ('A',),
        'residue_count': 76,
        'atom_count': 602,
        'het_atom_count': 58,
        'hetmol_counts': {'HOH': 58},
    },
    '4hhb.cif': {
        'model_ids': (1,),
        'chain_ids': ('A', 'B', 'C', 'D'),
        'residue_count': 574,
        'atom_count': 4384,
        'het_atom_count': 395,
        'hetmol_counts': {'HEM': 4, 'HOH': 221, 'PO4': 2},
    },
    '1a3n.cif': {
        'model_ids': (1,),
        'chain_ids': ('A', 'B', 'C', 'D'),
        'residue_count': 572,
        'atom_count': 4370,
        'het_atom_count': 623,
        'hetmol_counts': {'HEM': 4, 'HOH': 451},
    },
    '1g03.cif': {
        'model_ids': tuple(range(1, 21)),
        'chain_ids': ('A',),
        'residue_count': 134,
        'atom_count': 2081,
        'het_atom_count': 0,
        'hetmol_counts': {},
    },
    '1bna.cif': {
        'model_ids': (1,),
        'chain_ids': ('A', 'B'),
        'residue_count': 24,
        'atom_count': 486,
        'het_atom_count': 80,
        'hetmol_counts': {'HOH': 80},
    },
}

class Test_Molecule(unittest.TestCase):

    def setUp(self):
        self.mol = molecule.load_structure('packman/tests/data/4hla.cif', ftype='cif')

    def test_load_pdb(self):
        self.assertTrue( molecule.load_structure('packman/tests/data/1prw.pdb', ftype='pdb') )
    
    def test_pdb_fixture_matches_expected_output(self):
        for filename, expected in PDB_EXPECTED_OUTPUTS.items():
            with self.subTest(filename=filename):
                structure = molecule.load_structure(str(DATA_DIR / filename), ftype='pdb')
                self.assertIsInstance(structure, molecule.Protein)

                models = list(structure)
                self.assertEqual(len(models), expected['model_count'])
                model = models[0]
                self.assertEqual(model.get_id(), expected['model_id'])

                chains = list(model.get_chains())
                self.assertEqual(tuple(chain.get_id() for chain in chains), expected['chain_ids'])
                self.assertEqual(len(list(model.get_residues())), expected['residue_count'])
                self.assertEqual(len(list(model.get_atoms())), expected['atom_count'])
                self.assertEqual(len(list(model.get_calpha())), expected['calpha_count'])
                self.assertEqual(
                    sum(len(list(residue.get_backbone())) for residue in model.get_residues()),
                    expected['backbone_count'],
                )
                self.assertEqual(
                    ''.join(chains[0].get_sequence().splitlines()[1:]),
                    expected['sequence'],
                )

                first_atom = next(model.get_atoms())
                self.assertEqual(
                    (
                        first_atom.get_id(),
                        first_atom.get_name(),
                        first_atom.get_element(),
                        tuple(first_atom.get_location()),
                    ),
                    expected['first_atom'],
                )
                self.assertIs(first_atom.get_parent().get_parent(), chains[0])
                self.assertIs(chains[0].get_parent(), model)
                self.assertIs(model.get_parent(), structure)

    def test_pdb_edge_case_fixtures(self):
        for filename, expected in PDB_EDGE_CASE_OUTPUTS.items():
            with self.subTest(filename=filename):
                path = DATA_DIR / filename
                self.assertTrue(path.is_file(), 'Expected PDB edge-case fixture is missing.')
                structure = molecule.load_structure(str(path), ftype='pdb')
                models = list(structure)

                self.assertEqual(tuple(model.get_id() for model in models), expected['model_ids'])
                model = models[0]
                for model in models:
                    with self.subTest(model_id=model.get_id()):
                        chains = list(model.get_chains())
                        self.assertEqual(tuple(chain.get_id() for chain in chains), expected['chain_ids'])
                        self.assertEqual(len(list(model.get_residues())), expected['residue_count'])
                        self.assertEqual(len(list(model.get_atoms())), expected['atom_count'])

                if 'hetmol_counts' in expected:
                    self.assertEqual(
                        Counter(hetmol.get_name() for hetmol in model.get_hetmols()),
                        Counter(expected['hetmol_counts']),
                    )
                if 'het_atom_count' in expected:
                    self.assertEqual(len(list(model.get_hetatoms())), expected['het_atom_count'])
                if 'altloc_record_count' in expected:
                    # Atom does not expose altloc labels; ensure alternate records are present and retained.
                    lines = path.read_text(encoding='ascii').splitlines()
                    alternate_records = [
                        line for line in lines
                        if line.startswith(('ATOM  ', 'HETATM')) and line[16:17].strip()
                    ]
                    self.assertEqual(len(alternate_records), expected['altloc_record_count'])

    def test_all_pdb_files_in_data_folder_parse(self):
        pdb_files = sorted(DATA_DIR.glob('*.pdb'))
        self.assertTrue(pdb_files, 'No PDB files found in the test data folder.')
        self.assertTrue(
            set(PDB_EXPECTED_OUTPUTS).issubset({path.name for path in pdb_files}),
            'A PDB fixture listed in PDB_EXPECTED_OUTPUTS is missing from the test data folder.',
        )

        for pdb_file in pdb_files:
            with self.subTest(filename=pdb_file.name):
                structure = molecule.load_structure(str(pdb_file), ftype='pdb')
                self.assertIsInstance(structure, molecule.Protein)
                models = list(structure)
                self.assertTrue(models, 'Parsed structure contains no models.')

                for model in models:
                    self.assertIs(model.get_parent(), structure)
                    chains = list(model.get_chains())
                    self.assertTrue(chains, 'Parsed model contains no chains.')
                    atoms = list(model.get_atoms())
                    self.assertTrue(atoms, 'Parsed model contains no atoms.')

                    for chain in chains:
                        self.assertIs(chain.get_parent(), model)
                        self.assertTrue(list(chain.get_residues()))
                    for atom in atoms:
                        self.assertIs(atom.get_parent().get_parent().get_parent(), model)
                        self.assertEqual(atom.get_location().shape, (3,))
                        self.assertTrue(numpy.isfinite(atom.get_location()).all())

    def test_load_cif(self):
        self.assertTrue( molecule.load_structure('packman/tests/data/4hla.cif', ftype='cif') )

    def test_cif_fixture_matches_expected_output(self):
        for filename, expected in CIF_EXPECTED_OUTPUTS.items():
            with self.subTest(filename=filename):
                structure = molecule.load_structure(str(DATA_DIR / filename), ftype='cif')
                self.assertIsInstance(structure, molecule.Protein)

                models = list(structure)
                self.assertEqual(len(models), expected['model_count'])
                model = models[0]
                self.assertEqual(model.get_id(), expected['model_id'])

                chains = list(model.get_chains())
                self.assertEqual(
                    {chain.get_id(): ''.join(chain.get_sequence().splitlines()[1:]) for chain in chains},
                    expected['chain_sequences'],
                )
                self.assertEqual(len(list(model.get_residues())), expected['residue_count'])
                self.assertEqual(len(list(model.get_atoms())), expected['atom_count'])
                self.assertEqual(len(list(model.get_calpha())), expected['calpha_count'])
                self.assertEqual(
                    sum(len(list(residue.get_backbone())) for residue in model.get_residues()),
                    expected['backbone_count'],
                )
                self.assertEqual(len(list(model.get_bonds())), expected['bond_count'])

                first_atom = next(model.get_atoms())
                self.assertEqual(
                    (
                        first_atom.get_id(),
                        first_atom.get_name(),
                        first_atom.get_element(),
                        tuple(first_atom.get_location()),
                    ),
                    expected['first_atom'],
                )
                self.assertIs(first_atom.get_parent().get_parent(), chains[0])
                self.assertIs(chains[0].get_parent(), model)
                self.assertIs(model.get_parent(), structure)

    def test_cif_edge_case_fixtures(self):
        for filename, expected in CIF_EDGE_CASE_OUTPUTS.items():
            with self.subTest(filename=filename):
                path = DATA_DIR / filename
                self.assertTrue(path.is_file(), 'Expected CIF edge-case fixture is missing.')
                structure = molecule.load_structure(str(path), ftype='cif')
                models = list(structure)

                self.assertEqual(tuple(model.get_id() for model in models), expected['model_ids'])
                for model in models:
                    with self.subTest(model_id=model.get_id()):
                        chains = list(model.get_chains())
                        self.assertEqual(tuple(chain.get_id() for chain in chains), expected['chain_ids'])
                        self.assertEqual(len(list(model.get_residues())), expected['residue_count'])
                        self.assertEqual(len(list(model.get_atoms())), expected['atom_count'])
                        self.assertEqual(len(list(model.get_hetatoms())), expected['het_atom_count'])
                        self.assertEqual(
                            Counter(hetmol.get_name() for hetmol in model.get_hetmols()),
                            Counter(expected['hetmol_counts']),
                        )

    def test_all_cif_files_in_data_folder_parse(self):
        cif_files = sorted(DATA_DIR.glob('*.cif'))
        self.assertTrue(cif_files, 'No CIF files found in the test data folder.')
        expected_fixtures = set(CIF_EXPECTED_OUTPUTS) | set(CIF_EDGE_CASE_OUTPUTS)
        self.assertTrue(
            expected_fixtures.issubset({path.name for path in cif_files}),
            'A CIF fixture listed in CIF_EXPECTED_OUTPUTS is missing from the test data folder.',
        )

        for cif_file in cif_files:
            with self.subTest(filename=cif_file.name):
                structure = molecule.load_structure(str(cif_file), ftype='cif')
                self.assertIsInstance(structure, molecule.Protein)
                models = list(structure)
                self.assertTrue(models, 'Parsed structure contains no models.')

                for model in models:
                    self.assertIs(model.get_parent(), structure)
                    chains = list(model.get_chains())
                    self.assertTrue(chains, 'Parsed model contains no chains.')
                    atoms = list(model.get_atoms())
                    self.assertTrue(atoms, 'Parsed model contains no atoms.')

                    for chain in chains:
                        self.assertIs(chain.get_parent(), model)
                        self.assertTrue(list(chain.get_residues()))
                    for atom in atoms:
                        self.assertIs(atom.get_parent().get_parent().get_parent(), model)
                        self.assertEqual(atom.get_location().shape, (3,))
                        self.assertTrue(numpy.isfinite(atom.get_location()).all())

    def test_Protein(self):
        #Basic
        self.assertIsInstance( self.mol, molecule.Protein )

        #Write Methods
        self.assertTrue( self.mol.write_structure( 'test.cif', ftype='cif') )
        self.assertTrue( self.mol.write_structure( 'test.pdb', ftype='pdb') )
        with self.assertRaises(Exception) as file_format_problem:
            self.assertTrue( self.mol.write_structure( 'test', ftype='random') )
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