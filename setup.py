import packman
from setuptools import setup

PACKAGES=['packman',
          'packman.molecule',
          'packman.anm',
          'packman.apps',
          'packman.bin',
          'packman.constants',
          'packman.utilities'
          ]


SCRIPTS=['packman=packman.bin.PACKMAN:main']


setup(name='py-packman',
      version=packman.__version__,
      description='A software package for molecular PACKing and Motion ANalysis (PACKMAN)',
      url='https://github.com/Pranavkhade/PACKMAN',
      author='Pranav Khade',
      author_email='pranavk@iastate.edu',
      license='MIT',
      packages=PACKAGES,
      long_description = 'This package focuses on studying molecular structures and their dynamics using a simple yet informative property known as Protein Packing. Over the last few years, we have worked on several techniques to capture and quantify the protein packing, resulting in a few publications. This package has all the code to repeat and further develop these techniques.\n\nPlease visit for more details: https://github.com/Pranavkhade/PACKMAN',
      keywords=('protein, dynamics, protein packing, protein domain, protein hinge'),
      classifiers=[
              'Intended Audience :: Education',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: MIT License',
              'Operating System :: MacOS',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Programming Language :: Python',
              'Programming Language :: Python :: 2',
              'Programming Language :: Python :: 3',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Topic :: Scientific/Engineering :: Chemistry',
              'Topic :: Scientific/Engineering :: Mathematics'
                ],
      entry_points = {
              'console_scripts': SCRIPTS,
                },
    install_requires=['numpy', 'scipy', 'networkx', 'mlxtend', 'scikit-learn'],
      )