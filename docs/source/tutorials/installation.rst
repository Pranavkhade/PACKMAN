.. _tutorials_installation:

Install PACKMAN
===============

Requirements
------------

Use a supported Python version; the project CI currently tests Python 3.8
through 3.12. A working C compiler is not normally required when installing
the published package because pip installs compatible wheels for its scientific
Python dependencies.

Install from PyPI
-----------------

Create and activate a virtual environment if you want to keep PACKMAN's
dependencies separate, then install the package:

.. code-block:: console

   python -m pip install --upgrade pip
   python -m pip install py-packman

Check that the command-line interface is available:

.. code-block:: console

   python -m packman --help

Install from source
-------------------

To work from a local checkout:

.. code-block:: console

   git clone https://github.com/Pranavkhade/PACKMAN.git
   cd PACKMAN
   python -m pip install -e .

The editable install makes changes in the checkout available to the active
Python environment.

Graphical interface
-------------------

The GUI is launched with:

.. code-block:: console

   python -m packman gui

On some Linux distributions, Tkinter is a separate operating-system package.
Install the matching Tkinter package if the GUI reports that ``tkinter`` is
unavailable.

Next steps
----------

* :doc:`../quickstart` for a first structure and analysis.
* :doc:`gui` for graphical workflows.
* :doc:`main` for all tutorials.
