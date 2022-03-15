#!/usr/bin/env python
"""The PACKMAN Command-line Interface (CLI)

How to use::
    python -m packman gui #(For GUI)
    python -m packman     #(For CLI)

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""
import numpy

import traceback
import logging

from .. import molecule
from ..apps import hinge_cli, hdanm_cli, entropy_cli, dci_cli

import operator
import argparse
import os
import sys

'''
##################################################################################################
#                                          Interface                                             #
##################################################################################################
'''

def IO():
    """User interface for the user to provide the parameters

    Todo:
        * Make sure iterative alpha shape search over the values are accomodated in this interface
    
    Returns:
        Namespace: Various arguments in various formats
    """
    parser=argparse.ArgumentParser(description='PACKMAN: PACKing and Motion ANalysis. (https://github.com/Pranavkhade/PACKMAN)\n\nFollowing Apps Available: \n1. hinge \n2. hdanm\n3. entropy\n4. dci\n\nHow to run an app: python -m packman <app name>\nExample: python -m packman hinge', formatter_class=argparse.RawTextHelpFormatter )
    subparsers = parser.add_subparsers(dest='command')

    #Hinge Prediction
    hinge_app_io = subparsers.add_parser('hinge')
    hinge_app_io.add_argument('-pdbid','--pdbid', metavar='PDB_ID', type=str, help='If provided, the PBD with this ID will be downloaded and saved to FILENAME.')
    hinge_app_io.add_argument('alpha', metavar='AlphaValue', help='Recommended: 2.8 for closed; 4.5 for open form, Please refer to the paper for more details')
    hinge_app_io.add_argument('filename', metavar='FILENAME', help='Path and filename of the PDB file.')
    hinge_app_io.add_argument('--e_clusters',metavar='NumberOfEccentricityClusters',type=int,default=4,help='Recommended: 4, Please refer to the paper for more details')
    hinge_app_io.add_argument('--minhnglen',metavar='MinimumHingeLength',type=int,default=5,help='Recommended: 5, Please refer to the paper for more details')
    hinge_app_io.add_argument("--chain", help='Enter The Chain ID')
    hinge_app_io.add_argument('--generateobj', type=argparse.FileType('wb', 0), help='Path and filename to save the .obj file at. Ignored unless --chain is provided.')

    #hdanm
    hd_anm_io = subparsers.add_parser('hdanm')
    hd_anm_io.add_argument('-pdbid','--pdbid', metavar='PDB_ID', type=str, help='If provided, the PBD with this ID will be downloaded and saved to FILENAME.')
    hd_anm_io.add_argument('filename', metavar='FILENAME', help='Path and filename of the PDB file.')
    hd_anm_io.add_argument('hngfile', metavar='HNG', help='Path and filename of the corresponding HNG file.')
    hd_anm_io.add_argument("--chain", help='Enter The Chain ID')
    hd_anm_io.add_argument("--dr", type=float, default=15, help='Distance cutoff for the ANM.')
    hd_anm_io.add_argument("--power", type=float, default=0, help='Power of the distance in non-parametric ANM.')
    hd_anm_io.add_argument("--mass", default='residue', help='Mass of the residue; unit or molecular weight')
    hd_anm_io.add_argument("--scale", type=int, default=2, help='movie scale')
    hd_anm_io.add_argument("--frames", type=int, default=10, help='number of frames')
    hd_anm_io.add_argument("--modes", type=int, default=10, help='how many modes')
    hd_anm_io.add_argument("--ca_to_aa", action=argparse.BooleanOptionalAction, type=bool, default=False, help='Project CA motion on all atoms.')
    hd_anm_io.add_argument("--make_tar", action='store_true', help='package output files into a tar.gz file')

    #Entropy
    entropy_app_io = subparsers.add_parser('entropy')
    entropy_app_io.add_argument('-type','--type', metavar='entropy_type', type=str, help='Provide the Entropy type (Options: 1. PackingEntropy)')
    entropy_app_io.add_argument('-pdbid','--pdbid', metavar='PDB_ID', type=str, help='If provided, the PBD with this ID will be downloaded and saved to FILENAME.')
    entropy_app_io.add_argument('filename', metavar='FILENAME', help='Path and filename of the PDB file.')
    entropy_app_io.add_argument('--chains',metavar='Chains to be used for the entropy calculation',type=str,default=None, help='Recommended: None. Chain IDs for the Entropy calculation (None means all the chains are included; single string means only one chain ID; multiple chains should be comma separated).')
    entropy_app_io.add_argument('--probe_size',metavar='Size surface probe radius',type=float,default=1.4, help='Recommended: 1.4 (radius of a water molecule), Please refer to the paper for more details')
    entropy_app_io.add_argument('--onspherepoints',metavar='Number of points on a sphere',type=int,default=30, help='Recommended: 30. Number of points to be generated around each point for the surface (Read the Publication for more details)')

    #DCI
    dci_app_io = subparsers.add_parser('dci')
    dci_app_io.add_argument('-pdbid','--pdbid', metavar='PDB_ID', type=str, help='If provided, the PBD with this ID will be downloaded and saved to FILENAME.')
    dci_app_io.add_argument('filename', metavar='FILENAME', help='Path and filename of the PDB file.')
    dci_app_io.add_argument('-chain','--chain', help='Enter The Chain ID', default=None)
    dci_app_io.add_argument('-cutoff','--cutoff', type=float, help='Enter the cutoff for DCI. (Read the Publication for more details)', default=7.0)
    dci_app_io.add_argument('-n_com','--n_com', type=int, help='Enter the number of communities. (Read the Publication for more details)', default=None)

    # web server parameters
    web_server_group = parser.add_argument_group('Web server parameters', 'Used by the web form')
    web_server_group.add_argument('--outputfile', type=argparse.FileType('w', 1), default=sys.stdout, help='Path and filename write output to')
    web_server_group.add_argument('--logfile', type=argparse.FileType('w', 1), default=sys.stderr, help='Path and filename write log messages to')
    
    args=parser.parse_args()
    return args

'''
##################################################################################################
#                                             CLI                                                #
##################################################################################################
'''

def load_cli():
    """Main Function for PACKMAN.

    Kernel of the PACKMAN interface.

    Todo:
        * Change the main() to accomodate iterative alpha shape search
    
    """
    args=IO()

    if(args.command==None):
        logging.error('Please provide the appropriate input. Enter "python -m packman -h" for more details.')
        exit()

    logging.basicConfig(stream=args.logfile)

    if(args.pdbid is not None):
        molecule.download_structure(args.pdbid, save_name=args.filename.split('.')[0], ftype=args.filename.split('.')[1])

    try:
        extension = args.filename.split('.')[-1]
        mol = molecule.load_structure(args.filename,ftype=extension)
    except:
        logging.warning("The filename provided does not appear to have a format extension.")
        mol = molecule.load_structure(args.filename)
    
    if(args.command == 'hinge'):
        hinge_cli(args,mol)
    elif(args.command == 'hdanm'):
        hdanm_cli(args,mol)
    elif(args.command == 'entropy'):
        entropy_cli(args,mol)
    elif(args.command == 'dci'):
        dci_cli(args,mol)

    return True


'''
##################################################################################################
#                                            Main                                                #
##################################################################################################
'''


def main():
    """Gatewayto CLI and GUI
    """
    if(len(sys.argv) >  1 and sys.argv[1] == "gui"):
        from .GUI import load_gui
        load_gui()
    else:
        load_cli()


if(__name__ == '__main__'):
    main()
