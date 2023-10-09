.. _tutorials_gnm_api:

Using Gaussian Network Model (GNM)
==================================

This tutorial aims to get the user to familiarize with the concept of GNM API. GNM is available only in API.

How to cite::

    - For original GNM: Tirion, M.M. (1996). "Large amplitude elastic motions in proteins from a single-parameter, atomic analysis". Phys. Rev. Lett. 77 (9): 1905-1908.
    - For using GNM with this package: Pranav M Khade, Robert L Jernigan, PACKMAN-Molecule: Python Toolbox for Structural Bioinformatics, Bioinformatics Advances, 2022;, vbac007, https://doi.org/10.1093/bioadv/vbac007
    - For 7.8 Angstrom cutoff: Yang, L., Song, G. & Jernigan, R. L. Protein elastic network models and the ranges of cooperativity. Proc Natl Acad Sci U S A 106, 12347-52 (2009).
    - For PACKMAN molecule API: Pranav M Khade, Robert L Jernigan, PACKMAN-Molecule: Python toolbox for structural bioinformatics, Bioinformatics Advances, Volume 2, Issue 1, 2022, vbac007, https://doi.org/10.1093/bioadv/vbac007

Note: If PACKMAN is not installed, please follow the link: https://github.com/Pranavkhade/PACKMAN


Step 1: Loading Structure
-------------------------

In this section, we load the file containing the molecular structure. To understand in dept about the :mod:`packman.molecule` object, pelase visit :ref:`tutorials_molecule`

Code Example::

    from packman import molecule

    #File loading 
    mol = molecule.load_structure('1prw.cif')
    #OR
    #mol = molecule.load_structure('1prw.pdb')

    #If the c-alpha atoms of the ALL the chains need to be studied,
    #calpha = [i for i in mol[0].get_calpha()  if i is not None]

    #If the c-alpha atoms of the specific chain 'A' needs to be studied,
    calpha = [i for i in mol[0]['A'].get_calpha() if i is not None]


Step 2: Loading the GNM
------------------------

The GNM object has one required argument (set of atoms) and several optional arguments that can be changed according to the user's needs. Please check :mod:`packman.gnm.GNM` for more details.

Code Example::

    from packman.gnm import GNM
    model = GNM(calpha)


Step 3: Calculating Kirchhoff's Matrix
--------------------------------------

Please refer to the publication mentioned at the top of this page for more details.

In future, we will allow users to provide amino acid residue specific molecular masses, atom specific atomic masses, or any other type of weights that needs to be assigned to the atoms/residues. To give custom residue-specific values as mass to the method, simply create a standard python dictionary with 20 standard amino acids as keys and assign corresponding values to it (mass_type parameter).

Code Example::

    model.calculate_kirchhoff()
    kirchhoff = model.get_kirchhoff()
    print(kirchhoff)


Step 4: Eigenvalue decomposition
--------------------------------

The Eigenvalue decomposition is carried out using the Kirchhoff's matrix and mass matrix to calculate eigenvalues and eigenvectors.

Code Example::

    model.calculate_decomposition()

Eigenvalues and Eigenvectors. can be obtained by::

	model.get_eigenvalues()
	model.get_eigenvectors()

Respectively. Other calculated attributes and properties of the hd-ANM built can be obtained by its 'get' methods. Please refer to :mod:`packman.gnm.GNM` documentation for more details.


Step 6: Getting other data from GNM
-----------------------------------

This step demonstrates how one can obtain other data from GNM.


Code Example::

    model.calculate_fluctuations()
    print(model.get_fluctuations())

    model.calculate_crosscorrelation()
    print(model.get_pseudoinverse())

    from matplotlib import pyplot as plt
    plt.imshow( model.get_crosscorrelation(), cmap='hot')
    plt.show()