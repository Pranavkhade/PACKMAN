#import cProfile

#import time
#start_time = time.time()

from packman.molecule import load_structure

def main():
    #print("--- %s seconds ---" % (time.time() - start_time))

    mol = load_structure('1exr.cif')
    #print("--- %s seconds ---" % (time.time() - start_time))


    #Although the bonds are in the part of the chain, the model object should have the function as the inter chain bonds can exist
    mol[0].calculate_bonds()
    #print( mol[0]['A'].get_bonds() )
    #print("--- %s seconds ---" % (time.time() - start_time))


main()

#cProfile.run('main()')