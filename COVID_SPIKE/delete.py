

files=['6X6P','7BYR','6LZG','6M0J','6VSB','6VXX','6YBB','6YOR','7BZ5','6M17','6W41','6YM0','6YLA','6LVN','6WPS','6WPT','6X29','6X2A','6X2B','6X2C','6Z43']

from packman import molecule


all_conformers=[]

for i in files:
    try:
        molecule.download_structure(i,'delete.pdb')
        mol=molecule.load_structure('delete.pdb')
        chains=[j for j in mol[0].get_chains()]
        if(len(chains)>1):
            all_conformers.extend([i+'_'+j.get_id() for j in chains])
    except:
        print('Problem>>',i)
    

for i in all_conformers:
    print(i)

