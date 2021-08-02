import time
start_time = time.time()

from packman.molecule import load_structure
print("--- %s seconds ---" % (time.time() - start_time))

mol = load_structure('1exr.cif')
print("--- %s seconds ---" % (time.time() - start_time))

mol[0]['A'].calculate_bonds()
print( mol[0]['A'].get_bonds() )
print("--- %s seconds ---" % (time.time() - start_time))