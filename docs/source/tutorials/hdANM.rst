.. _tutorials_hdANM:

Using hd-ANM to study the global dynamics of the molecular structures
=====================================================================

This tutorial aims to get the user to familiarize with the concept of hdANM. 

How to cite::

    Paper Under Review.


Note: If PACKMAN is not installed, please follow the link: https://github.com/Pranavkhade/PACKMAN


Step 1: Loading Structure
-------------------------

In this section, we load the file containing the molecular structure. To understand in dept about the :mod:`packman.molecule` object, pelase visit :ref:`tutorials_molecule`

Code Example::

    from packman import molecule

    #File loading 
    mol=molecule.load_structure('1exr.pdb')

    #If the c-alpha atoms of the ALL the chains need to be studied,
    #calpha=[i for i in mol[0].get_calpha()  if i is not None]

    #If the c-alpha atoms of the specific chain 'A' needs to be studied,
    calpha=[i for i in mol[0]['A'].get_calpha() if i is not None]


Step 2: Loading the hdANM Object
--------------------------------

The hd-ANM object has one required argument (set of atoms) and several optional arguments that can be changed according to the user's needs. Please check :mod:`packman.anm.hd_anm` for more details.

Important Note: hdANM also requires the user to provide the .hng file containing the information about which parts of the protein are hinges and which parts are rigid domains. This information can be obtained from the PACKMAN's hinge prediction algorithm. Please read the following for more details: :ref:`tutorials_hinge_prediction`

Hinge Information file format (.hng) Example::

    1exr.pdb_A  D1  1:70
    1exr.pdb_A  H1  71:90
    1exr.pdb_A  D2  91:148

The first column in the .hng file is Filename_ChainID, the second column is Domain/Hinge ID, and the third column is the residues in the particular domain/hinge. The .hng file is tab-separated. (Tabs separate columns)

Code Example::

    from packman.anm import hdANM
    Model=hdANM(calpha,dr=15,power=0,hng_file='1exr.pdb.hng')

Step 3: Calculating Hessian Matrix
----------------------------------

Unlike the Hessian matrix of the regular ANM, the hd-ANM Hessian matrix has reduced dimensions. Please refer to the publication mentioned at the top of this page for more details.

Along with the difference mentioned above, we allow users to provide amino acid residue specific molecular masses, atom specific atomic masses, or any other type of weights that needs to be assigned to the atoms/residues. To give custom residue-specific values as mass to the method, simply create a standard python dictionary with 20 standard amino acids as keys and assign corresponding values to it (mass_type parameter).

Code Example::

    #Unit mass for all the components
    #Model.calculate_hessian(mass_type='unit')

    #Atomic Weight Assignment
    #Model.calculate_hessian(mass_type='atom')

    #Residue Molecular Weight Assignment
    Model.calculate_hessian(mass_type='residue')


Step 4: Eigenvalue decomposition
--------------------------------

The Eigenvalue decomposition is carried out using the hessian matrix and mass matrix to calculate eigenvalues and eigenvectors.

Code Example::

    Model.calculate_decomposition()

Eigenvalues and Eigenvectors. can be obtained by::

	Model.get_eigenvalues()
	Model.get_eigenvectors()

Respectively. Other calculated attributes and properties of the hd-ANM built can be obtained by its 'get' methods. Please refer to :mod:`packman.anm.hd_anm` documentation for more details.


Step 5: Eigenvector Motion Extrapolation
-----------------------------------------

The modes obtained from Step 4 can be visualized on the molecular structure by extrapolating them linearly on curvilinearly by adjusting the parameters of :func:`packman.anm.hd_anm.calculate_movie`. By default, the program gives curvilinear extrapolation of the Eigenvector motions.

Important Note: 7th Mode is the first non-rigid mode (0 to 6 indices are not excluded)

Code Example::
    
    Model.calculate_movie(6,scale=2,n=10)

The '6.pdb' file will be saved on the present working directory containing the motion for the 7th (First Non Rigid) Mode.


Step 6: Getting hdANM output matrices (Hessian Pseudoinverse / Cross-Correlation Matrix)
----------------------------------------------------------------------------------------

This step can be done before generating movies as well. In order to get the hdANM output matrices such as Hessian Pseudoinverse / Cross-Correlation Matrix/ Reverse Transformed Eigenvectors. Please read the paper for more details about the theory and importance of these matrices.

Note: Reverse Transformed Eigenvectors has dimension: 3N x 6D+3H (N= Number of atoms, D= Number of domains & H= Number of hinge Atoms)

Code Example::

    #Here, 'n_modes' variable is number of first non-rigid modes to get the result matrices

    #For the Hessian Pseudoinverse,
    Model.get_hessian_pseudoinverse(n_modes)

    #For the Correlation Matrix,
    Model.get_crosscorrelation_matrix(n_modes=10)

    #For the Reverse Transformed Eigenvectors,
    Model.get_RT_eigen_vectors()