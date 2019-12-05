'''
import molecule
mol=molecule.load_structure('1prw.pdb')
print([i.get_id() for i in mol[0].get_chains()])
'''

from packman.bin import PACKMAN


PACKMAN.main()