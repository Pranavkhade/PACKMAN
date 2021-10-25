.. PACKing and Motion ANalysis (PACKMAN) documentation master file, created by
   sphinx-quickstart on Fri Mar 27 02:56:53 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

py-PACKMAN Documentation
========================

.. toctree::
   :maxdepth: 2

This package focuses on studying molecular structures and their dynamics using a simple yet informative property known as Protein Packing. Over the last few years, we have worked on several techniques to capture and quantify the protein packing, resulting in a few publications. This package has all the code to repeat and further develop these techniques.

:ref:`Tutorials <tutorials_main>`
---------------------------------
The guides to using the PACKMAN GUI, CLI, and API for different submodules are included in the :ref:`tutorials <tutorials_main>`.

Documentation
-------------
- :ref:`genindex`
- :ref:`modindex`

Web Servers
-----------
- `Hinge Prediction <https://packman.bb.iastate.edu/>`_
- `hdANM <https://hdanm.bb.iastate.edu/>`_
- `Packing Entropy <https://packing-entropy.bb.iastate.edu/>`_

Publications
------------
- Hinge Prediction (https://doi.org/10.1016/j.jmb.2019.11.018)
- Structural Compliance (https://doi.org/10.1002/prot.25968)
- hdANM (https://doi.org/10.1016/j.bpj.2021.10.017)

Project Goal
------------
There are only limited methods available to study the global motions of the protein such as their hinge motions and shear motions. These motions take place over a broad range of time scales, from microseconds to seconds; however, molecular dynamics methods can only model easily the motions occurring on the time scale from picoseconds to microseconds, and in addition, such simulations require that replicas be run. Thus, extracting the meaningful slow motions is difficult. Hence, there is a need to explore other ways to utilize the protein structures to model/ predict the global/large scale motions of the protein. The important motions depend on the protein packing as a multiscale phenomenon that can influence the global or local motions in the proteins. To model the protein packing, an efficient, robust and simple mathematical method is needed. In this study, we have explored alpha shapes (a subset of Delaunay tessellations) for the protein backbone coordinates as a model of protein packing. We demonstrate that the method can predict the protein hinges which are responsible for the global motions of the proteins. The method is named PACKMAN (PACking and Motion ANalyses). From a literature survey for randomly selected protein structures and another 367 protein structure pairs having known open and closed conformations, PACKMAN can predict hinge locations accurately on the proteins for both the open and closed conformations outperforming existing hinge predicting methods that usually require either the open form or both conformations to predict the protein hinges. The successful implementation of the mathematical method to model a multiscale phenomenon such as protein packing to predict the hotspots of the global/large scale motions in proteins shows promise for the further exploration of other types of protein and supramolecular dynamics.

Acknowledgments
---------------
The development of the PACKMAN package is supported by `NSF <https://www.nsf.gov/>`_ grant DBI 1661391. The authors also thank `ResearchIT@Iowa State University <https://researchit.las.iastate.edu/about-research-it-iowa-state-university>`_ for helping with many aspects of computing.