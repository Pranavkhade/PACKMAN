.. _quickstart:

Quickstart
==========

This short walkthrough downloads a small example structure, explores its object
hierarchy, and runs a basic Gaussian Network Model (GNM) calculation.

Install PACKMAN in the Python environment you plan to use:

.. code-block:: console

   python -m pip install py-packman

Download and load a PDB structure:

.. code-block:: python

   from packman import molecule

   molecule.download_structure("1prw", save_name="1prw", ftype="pdb")
   protein = molecule.load_structure("1prw.pdb")

Inspect the first model and its chains:

.. code-block:: python

   model = protein[0]
   chains = list(model.get_chains())
   print("Model:", model.get_id())
   print("Chains:", [chain.get_id() for chain in chains])
   print("Residues:", len(list(model.get_residues())))
   print("Atoms:", len(list(model.get_atoms())))

The object hierarchy and the functions that load structures are shown in the
:doc:`api-map`.

Run a GNM calculation on chain A's C-alpha atoms:

.. code-block:: python

   from packman.gnm import GNM

   calpha = [atom for atom in model["A"].get_calpha() if atom is not None]
   gnm = GNM(calpha)
   gnm.calculate_kirchhoff()
   gnm.calculate_decomposition()

   print("First eigenvalues:", gnm.get_eigenvalues()[:7])

For a connected GNM network, the lowest eigenvalue is the translational zero
mode; subsequent modes describe internal fluctuations. See
:doc:`tutorials/GNM` for further analysis.

Next steps
----------

* :doc:`tutorials/installation` for installation options and prerequisites.
* :doc:`tutorials/main` for GUI, CLI, and API workflows.
* :doc:`api-map` for the object/function tree and full API reference.
