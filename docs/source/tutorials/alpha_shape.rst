.. _tutorials_alpha_shape:

Calculate Alpha Shape of the Protein
====================================

This tutorial is to show how an API user can get the alpha shape of a protein.
    

Step 1
------
1. Import the submodules::
   
    from packman.molecule import  load_structure, download_structure
    from packman.geometry import AlphaShape, Circumsphere

Step 2
------
1. Download the PDB file
2. Load the PDB file into the 'mol' object::
   
    download_structure('1exr')
    mol = load_structure('1exr.cif')

Step 3
------
Get the atoms you wish to calculate the alpha shape from.

Here, we are selecting all the atoms of the protein. However, check the 'molecule' tutorials to select various other subgroups such as c-alpha atoms.::

    atoms = [i for i in mol[0].get_atoms()]

Step 4
------
Calculate the alpha shape with infinite cutoff parameter (Delaunay Tessellations)

It is important to note that the alpha shapes are just the subsets of Delaunay tessellations. Users can also calculate the radius of each tessellation in the Delaunay Tessellation by running Step 4.2::

    #4.1
    Tessellations = AlphaShape( atoms,float('Inf') )

    #4.2
    Circumsphere( Tessellations[0] )

Step 5
------
We can also calculate the Alpha shape with specific Alpha value as a second positional argument as following::
    
    #Here, Alpha = 3
    Tessellations = AlphaShape( atoms, 3 )

The output of the AlphaShape will be a 2D array consisting of 1D arrays of the set of 4 atoms as components of the tessellations. ie.. [ [atom1 of tessellation_1,...,atom4 of tessellation_1], ..... , [atom1 of tessellation_n,...,atom4 of tessellation_n] ]