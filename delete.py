


from typing import Counter
from packman.molecule import Bond
from packman import molecule

import logging

#This connectivity information is collected from the following file: http://ftp.wwpdb.org/pub/pdb/data/monomers/aa-variants-v1.cif.gz
aa_connectivity= {
'ALA': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'HB1'), ('CB', 'HB2'), ('CB', 'HB3'), ('OXT', 'HXT')], 
'ARG': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'CD'), ('CG', 'HG2'), ('CG', 'HG3'), ('CD', 'NE'), ('CD', 'HD2'), ('CD', 'HD3'), ('NE', 'CZ'), ('NE', 'HE'), ('CZ', 'NH1'), ('CZ', 'NH2'), ('NH1', 'HH11'), ('NH1', 'HH12'), ('NH2', 'HH21'), ('NH2', 'HH22'), ('OXT', 'HXT')], 
'ASN': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'OD1'), ('CG', 'ND2'), ('ND2', 'HD21'), ('ND2', 'HD22'), ('OXT', 'HXT')], 
'ASP': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'OD1'), ('CG', 'OD2'), ('OD2', 'HD2'), ('OXT', 'HXT')], 
'CYS': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'SG'), ('CB', 'HB2'), ('CB', 'HB3'), ('SG', 'HG'), ('OXT', 'HXT')], 
'GLN': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'CD'), ('CG', 'HG2'), ('CG', 'HG3'), ('CD', 'OE1'), ('CD', 'NE2'), ('NE2', 'HE21'), ('NE2', 'HE22'), ('OXT', 'HXT')], 
'GLU': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'CD'), ('CG', 'HG2'), ('CG', 'HG3'), ('CD', 'OE1'), ('CD', 'OE2'), ('OE2', 'HE2'), ('OXT', 'HXT')], 
'GLY': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'HA2'), ('CA', 'HA3'), ('C', 'O'), ('C', 'OXT'), ('OXT', 'HXT')],
'HIS': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'ND1'), ('CG', 'CD2'), ('ND1', 'CE1'), ('ND1', 'HD1'), ('CD2', 'NE2'), ('CD2', 'HD2'), ('CE1', 'NE2'), ('CE1', 'HE1'), ('NE2', 'HE2'), ('OXT', 'HXT')],
'ILE': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG1'), ('CB', 'CG2'), ('CB', 'HB'), ('CG1', 'CD1'), ('CG1', 'HG12'), ('CG1', 'HG13'), ('CG2', 'HG21'), ('CG2', 'HG22'), ('CG2', 'HG23'), ('CD1', 'HD11'), ('CD1', 'HD12'), ('CD1', 'HD13'), ('OXT', 'HXT')],
'LEU': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'CD1'), ('CG', 'CD2'), ('CG', 'HG'), ('CD1', 'HD11'), ('CD1', 'HD12'), ('CD1', 'HD13'), ('CD2', 'HD21'), ('CD2', 'HD22'), ('CD2', 'HD23'), ('OXT', 'HXT')], 
'LYS': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'CD'), ('CG', 'HG2'), ('CG', 'HG3'), ('CD', 'CE'), ('CD', 'HD2'), ('CD', 'HD3'), ('CE', 'NZ'), ('CE', 'HE2'), ('CE', 'HE3'), ('NZ', 'HZ1'), ('NZ', 'HZ2'), ('NZ', 'HZ3'), ('OXT', 'HXT')], 
'MET': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'SD'), ('CG', 'HG2'), ('CG', 'HG3'), ('SD', 'CE'), ('CE', 'HE1'), ('CE', 'HE2'), ('CE', 'HE3'), ('OXT', 'HXT')],
'PHE': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'CD1'), ('CG', 'CD2'), ('CD1', 'CE1'), ('CD1', 'HD1'), ('CD2', 'CE2'), ('CD2', 'HD2'), ('CE1', 'CZ'), ('CE1', 'HE1'), ('CE2', 'CZ'), ('CE2', 'HE2'), ('CZ', 'HZ'), ('OXT', 'HXT')],
'PRO': [('N', 'CA'), ('N', 'CD'), ('N', 'H'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'CD'), ('CG', 'HG2'), ('CG', 'HG3'), ('CD', 'HD2'), ('CD', 'HD3'), ('OXT', 'HXT')],
'SER': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'OG'), ('CB', 'HB2'), ('CB', 'HB3'), ('OG', 'HG'), ('OXT', 'HXT')],
'THR': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'OG1'), ('CB', 'CG2'), ('CB', 'HB'), ('OG1', 'HG1'), ('CG2', 'HG21'), ('CG2', 'HG22'), ('CG2', 'HG23'), ('OXT', 'HXT')],
'TRP': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'CD1'), ('CG', 'CD2'), ('CD1', 'NE1'), ('CD1', 'HD1'), ('CD2', 'CE2'), ('CD2', 'CE3'), ('NE1', 'CE2'), ('NE1', 'HE1'), ('CE2', 'CZ2'), ('CE3', 'CZ3'), ('CE3', 'HE3'), ('CZ2', 'CH2'), ('CZ2', 'HZ2'), ('CZ3', 'CH2'), ('CZ3', 'HZ3'), ('CH2', 'HH2'), ('OXT', 'HXT')],
'TYR': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG'), ('CB', 'HB2'), ('CB', 'HB3'), ('CG', 'CD1'), ('CG', 'CD2'), ('CD1', 'CE1'), ('CD1', 'HD1'), ('CD2', 'CE2'), ('CD2', 'HD2'), ('CE1', 'CZ'), ('CE1', 'HE1'), ('CE2', 'CZ'), ('CE2', 'HE2'), ('CZ', 'OH'), ('OH', 'HH'), ('OXT', 'HXT')],
'VAL': [('N', 'CA'), ('N', 'H'), ('N', 'H2'), ('CA', 'C'), ('CA', 'CB'), ('CA', 'HA'), ('C', 'O'), ('C', 'OXT'), ('CB', 'CG1'), ('CB', 'CG2'), ('CB', 'HB'), ('CG1', 'HG11'), ('CG1', 'HG12'), ('CG1', 'HG13'), ('CG2', 'HG21'), ('CG2', 'HG22'), ('CG2', 'HG23'), ('OXT', 'HXT')]
}


mol = molecule.load_structure('1exr.cif')
resi = [i for i in mol[0].get_residues()]
counter = 0
AllBonds = []
for i in resi:
    try:
        for j in aa_connectivity[i.get_name()]:
            try:
                atom1, atom2 =  i.get_atom(j[0]), i.get_atom(j[1])
                if(atom1!=None and atom2!=None and atom1.get_parent().get_parent().get_id() == atom2.get_parent().get_parent().get_id() ):
                    counter += 1
                    bond = Bond(counter,atom1,atom2,source='RCSB/aa-variants-v1.cif')
                    AllBonds.append(bond)
                    atom1.set_bond(bond)
                    atom2.set_bond(bond)
            except Exception as e:
                #Bond pair not found
                print(e)
                None
    except:
        logging.warn('Residue Number|Name|Chain '+str(i.get_id())+'|'+str(i.get_name())+'|'+str(i.get_parent().get_id())+' is not a standard amino acid.')

print( [i.get_bonds() for i in resi[0].get_atoms()] )