[![PyPI version](https://badge.fury.io/py/py-packman.svg)](https://badge.fury.io/py/py-packman) [![Continuous Integration](https://github.com/Pranavkhade/PACKMAN/actions/workflows/python-package.yml/badge.svg)](https://github.com/Pranavkhade/PACKMAN/actions/workflows/python-package.yml) [![Documentation Status](https://readthedocs.org/projects/py-packman/badge/?version=latest)](https://py-packman.readthedocs.io/en/latest/?badge=latest) [![Downloads](https://pepy.tech/badge/py-packman)](https://pepy.tech/project/py-packman)


PACKMAN: PACKing and Motion ANalysis
------------------------------------
<img src="https://github.com/Pranavkhade/PACKMAN/blob/master/docs/_static/gallary/logo.gif" width="250">

This package focuses on studying molecular structures and their dynamics using a simple yet informative property known as Protein Packing. Over the last few years, we have worked on several techniques to capture and quantify the protein packing, resulting in a few publications. This package has all the code to repeat and further develop these techniques.

✨ What's new?
* New Streamlit based Web-GUI is launched which allows users to scan multiple alpha values in range and combine hinges along with option to download the `.hng` file for downstream hdANM analysis. Please read the instructions below for more details.
* Improved support for hinge prediction across multiple alpha parameters. You can now run and combine results for multiple structures using a single script: additional codes/ScanAlpha.py. This streamlines batch processing and saves time. Run `python3 ScanAlpha.py -h` for more information (Copy the script to PWD).

DOCUMENTATION
-------------
* [IMPORTANT : Documentation, Tutorials & More](https://py-packman.readthedocs.io)

WEB-GUI
-------
This option allows user to use the hinge prediction, hdANM and Entropy calculation all in one place with few clicks. If users do not wish to use any commandline and code, this is the best option. Please install [docker](https://www.docker.com/) software and run the following commands.
```
docker pull pranavkhade/packman-streamlit
docker run -p 8501:8501 pranavkhade/packman-streamlit:latest
```

If a webpage doesn't automatically open use either of the following address on your web-browser (Google Chrome/ Firefox etc.)
```
http://localhost:8501
OR
http://0.0.0.0:8501
```

INSTALLATION
------------

1. Installing with pip (Recommended)
```
pip install py-packman
```

2. Installing from source
```
git clone https://github.com/Pranavkhade/PACKMAN
cd PACKMAN
python setup.py install
```

HOW TO USE
----------

PACKMAN, along with its components, can be accessed via Graphical User Interface (GUI), Command-line Interface (CLI), and Application programming interface (API).

For the GUI, please run the following command.
```
python -m packman gui
```
OR
```
python3 -m packman gui
```

For the CLI and API, please read the tutorials & documentation.

REPORTING ISSUES
----------------
[Click here](https://github.com/Pranavkhade/PACKMAN/issues/new) to report issues with any aspect of the package. Please search for the resolved issues first.

WEB SERVERS
-----------
* [Hinge Prediction Web Server](https://packman.bb.iastate.edu/)
* [hd-ANM Web Server](https://hdanm.bb.iastate.edu/)
* [Packing Entropy Web Server](https://packing-entropy.bb.iastate.edu/)
* [DCI web server](https://dci.bb.iastate.edu/)

RELATED RESEARCH ARTICLES
-------------------------
* [PACKMAN-molecule](https://doi.org/10.1093/bioadv/vbac007)
* [Hinge Prediction](https://doi.org/10.1016/j.jmb.2019.11.018)
* [Compliance](https://doi.org/10.1002/prot.25968)
* [hd-ANM](https://doi.org/10.1016/j.bpj.2021.10.017)
* [Packing Entropy](https://doi.org/10.1021/acsomega.2c00999)
* [DCI](https://doi.org/10.1093/bioinformatics/btac159)

ACKNOWLEDGEMENTS
----------------

The development of the PACKMAN package was in past from 2017 to 2022 and webservers still is supported by [NSF](https://www.nsf.gov/) grant DBI 1661391. The authors also thank [ResearchIT@Iowa State University](https://researchit.las.iastate.edu/about-research-it-iowa-state-university) for helping with server deployment @Iowa State University.

Development past April 2022 are done without any financial or any other support of any kind.
