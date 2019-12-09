from setuptools import setup

PACKAGES=['packman',
          'packman.molecule',
          'packman.apps',
          'packman.bin'
          ]


SCRIPTS=['packman=packman.bin.PACKMAN:main']


setup(name='py-packman',
      version='1.0',
      description='A software package for molecular PACKing and Motion ANalysis (PACKMAN)',
      url='https://github.com/Pranavkhade/PACKMAN',
      author='Pranav Khade',
      author_email='pranavk@iastate.edu',
      license='MIT',
      packages=PACKAGES,
      keywords=('protein, dynamics, protein packing'),
      classifiers=[
              'Intended Audience :: Education',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: MIT License',
              'Operating System :: MacOS',
              'Operating System :: Microsoft :: Windows',
              'Programming Language :: Python',
              'Programming Language :: Python :: 2',
              'Programming Language :: Python :: 3',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Topic :: Scientific/Engineering :: Structural Biology',
              'Topic :: Scientific/Engineering :: Chemistry'
                ],
      entry_points = {
              'console_scripts': SCRIPTS,
                },
    install_requires=['numpy', 'scipy', 'networkx', 'mlxtend', 'sklearn'],
      )