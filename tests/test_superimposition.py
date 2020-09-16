

from packman.utilities import superimporse


from packman.molecule import load_structure



mol1=load_structure('1K20/1k20.pdb')
mol2=load_structure('1K20/1k23.pdb')

superimporse(mol1[0]['A'],mol2[0]['A'],use='backbone')