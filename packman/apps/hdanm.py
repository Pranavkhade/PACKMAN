# -*- coding: utf-8 -*-
"""The 'hdANM' object host file.
This is file information, not the class information. This information is only for the API developers.
Please read the 'DCI' object documentation for details.

Example::
    from packman.apps import hdanm_cli
    help( hdanm_cli )

Note:
    * Tutorial: https://py-packman.readthedocs.io/en/1.4.10/tutorials/hdANM_cli.html#tutorials-hdanm-cli
    * For more details about the parameters for compliance, or to site this, read the following paper: 

Authors:
    * Ambuj Kumar  (ambuj@iastate.edu)
    * Pranav Khade (https://github.com/Pranavkhade)
"""
import numpy
import tarfile
import tempfile
import os

from ..anm import hdANM

def hdanm_cli(args,mol):

    #Load the Hinge Information
    chain    = args.chain
    dr       = float(args.dr)
    power    = float(args.power)

    if(chain is not None):
        calpha=[i for i in mol[0][chain].get_calpha() if i is not None]
    else:
        calpha=[i for i in mol[0].get_calpha()  if i is not None]

    Model = hdANM(calpha,dr=dr,power=power,hng_file=args.hngfile)
    Model.calculate_hessian(mass_type=args.mass)
    Model.calculate_decomposition()

    def generate_output():
        numpy.savetxt("eigenvalues.csv", Model.get_eigenvalues(), delimiter=",")
        numpy.savetxt("eigenvectors.csv", Model.get_eigenvectors(), delimiter=",")

        for i in range(6,6+args.modes,1):
            Model.calculate_movie(i,scale=args.scale,n=args.frames,ca_to_aa=args.ca_to_aa)

    if args.make_tar:

        with tempfile.TemporaryDirectory() as tmpdir:
            with open(args.outputfile.name, 'wb') as output_wb:

                os.chdir(tmpdir)
                
                generate_output()
            
                tar = tarfile.open(mode="w:gz", fileobj=output_wb)
                _, _, filenames = next(os.walk(tmpdir))
                for name in filenames:
                    tar.add(name)
                tar.close()

    else:
        generate_output()