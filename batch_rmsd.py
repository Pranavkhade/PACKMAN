from packman.utilities import RMSD
from packman.molecule import load_structure

import sys

def main():
    #Input format: abc.pdb_A (filename_chain)
    filename1, chain1 = sys.argv[1].split('_')
    filename2, chain2 = sys.argv[2].split('_')

    #File Change
    mol1=load_structure(filename1)
    mol2=load_structure(filename2)

    print( RMSD(mol1[0][chain1],mol2[0][chain2]) )


    return True

if(__name__ == "__main__"):
    main()