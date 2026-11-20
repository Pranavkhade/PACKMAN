.. _troubleshooting:

Troubleshooting
===============

Installation and imports
------------------------

* **``No module named packman``** — install PACKMAN into the same environment
  used to run your script: ``python -m pip install py-packman``. If using a
  virtual environment, activate it before installing and running PACKMAN.
* **GUI does not start on Linux** — the GUI requires Tkinter. Install the
  operating system's Tkinter package (for example, ``python3-tk`` on Debian or
  Ubuntu) and retry ``python -m packman gui``.

Structure files
---------------

* **The structure does not load** — check that the path exists and that the
  file is a supported PDB or mmCIF file. If the file extension is missing or
  nonstandard, pass ``ftype="pdb"`` or ``ftype="cif"`` explicitly.
* **A chain or residue cannot be found** — inspect the parsed identifiers
  first. Chain IDs and residue numbering are taken from the input structure;
  they may differ from a viewer's display labels.
* **A hinge/domain analysis does not match the structure** — verify that the
  ``.hng`` file uses the correct chain IDs and residue-number ranges for the
  structure being analyzed.

Analysis and performance
------------------------

Packing entropy and alpha-shape calculations operate on atomic coordinates and
can require substantial memory and time for large structures or assemblies.
Start with a single chain or a smaller structure to confirm the workflow before
processing a large assembly.

Reporting a problem
-------------------

When opening an issue, include the PACKMAN version, Python version, operating
system, command or minimal code used, and the complete error message. Avoid
sharing private or unpublished structure data; a public PDB identifier or a
small reproducible example is usually sufficient.

See the `GitHub issue tracker <https://github.com/Pranavkhade/PACKMAN/issues>`_
for existing reports and support.
