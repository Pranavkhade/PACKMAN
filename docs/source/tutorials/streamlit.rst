.. _tutorials_streamlit:

PACKMAN Streamlit web app
=========================

The repository includes a browser-based Streamlit interface for three
workflows: hinge prediction, hdANM analysis, and Voronoi packing entropy. Each
workflow accepts an uploaded PDB/mmCIF file or a PDB ID.

What the app provides
---------------------

* **Hinge prediction:** scan alpha parameters, review predicted hinges, view
  their locations in an interactive structure, and download a ``.hng`` file
  for hdANM.
* **hdANM:** analyze selected chains using an uploaded ``.hng`` file; inspect
  eigenvalue and fluctuation plots; download eigenvalue/eigenvector CSVs and
  generated motion CIF files.
* **Packing entropy:** calculate per-residue and total entropy, with controls
  for chain selection, probe size, and surface-point sampling.
* **3D visualization:** inspect structures and generated hdANM motion in the
  browser.

For packing-entropy calculations, remove hydrogen atoms from the input
structure as recommended by the app.

Run locally
-----------

Clone the repository, install the app dependencies in a Python environment,
and launch Streamlit from the repository root:

.. code-block:: console

   git clone https://github.com/Pranavkhade/PACKMAN.git
   cd PACKMAN
   python -m pip install py-packman streamlit stmol py3Dmol ipython_genutils
   python -m streamlit run "additional codes/streamlit_app.py"

Streamlit prints a local URL (normally ``http://localhost:8501``) when the app
starts.

Run with Docker
---------------

Docker avoids installing the app's Python dependencies locally. Build the
included image from the repository root and start it:

.. code-block:: console

   docker build -f "additional codes/streamlit_app_dockerfile" -t packman-streamlit .
   docker run --rm -p 8501:8501 packman-streamlit

Open ``http://localhost:8501`` in a browser. The project also publishes a
`PACKMAN Streamlit image on Docker Hub
<https://hub.docker.com/r/pranavkhade/packman-streamlit>`_.
