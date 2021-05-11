#!/usr/bin/env python
"""The PACKMAN Command-line Interface (CLI) and Graphical User Interface (GUI) host file.

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
from ..apps import hinge_cli, entropy_cli

import operator
import argparse
import os
import sys

try:
    from urllib.request import urlopen
except Exception:
    from urllib2 import urlopen


'''
##################################################################################################
#                                    Non Algorithm Functions                                     #
##################################################################################################
'''

def WriteOBJ(atoms,faces, fh):
    """Write the .obj file to visualize the obtain alpha shape tesselations.
    
    Args:
        atoms (packman.molecule.Atom): Atoms (Just for the node records)
        faces ([float])              : SelectedTesselations (See the packman.apps.predict_hinge)
        fh (file)                    : Output file with .obj extension
    
    """
    NewIDs={i.get_id():numi+1 for numi,i in enumerate(atoms)}
    fh.write('mtllib master.mtl\ng\n'.encode())
    fh.write('usemtl atoms\n'.encode())
    for i in atoms:
        x,y,z=i.get_location()
        fh.write("v %f %f %f\n".encode()%(x,y,z))
    
    line='usemtl bonds\nl'
    for i in atoms:
        line=line+" "+str(NewIDs[i.get_id()])
    line=line+'\n'
    fh.write(line.encode())
    
    fh.write('usemtl faces\n'.encode())
    for i in faces:
        faces=[NewIDs[j.get_id()] for j in i]
        fh.write("f %i %i %i %i\n".encode()%(faces[0],faces[1],faces[2],faces[3]))
        #fh.write("l %i %i %i %i\n"%(faces[0],faces[1],faces[2],faces[3]))
    return True

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
    parser=argparse.ArgumentParser(description='PACKMAN: PACKing and Motion ANalysis. (https://github.com/Pranavkhade/PACKMAN)\n\nFollowing Apps Available: \n1. hinge \n2. entropy\n\nHow to run an app: python -m packman <app name>\nExample: python -m packman hinge', formatter_class=argparse.RawTextHelpFormatter )
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

    #Entropy
    entropy_app_io = subparsers.add_parser('entropy')
    entropy_app_io.add_argument('-type','--type', metavar='entropy_type', type=str, help='Provide the Entropy type (Options: 1. PackingEntropy)')
    entropy_app_io.add_argument('-pdbid','--pdbid', metavar='PDB_ID', type=str, help='If provided, the PBD with this ID will be downloaded and saved to FILENAME.')
    entropy_app_io.add_argument('filename', metavar='FILENAME', help='Path and filename of the PDB file.')
    entropy_app_io.add_argument('--chains',metavar='Chains to be used for the entropy calculation',type=str,default=None, help='Recommended: None. Chain IDs for the Entropy calculation (None means all the chains are included; single string means only one chain ID; multiple chains should be comma separated).')
    entropy_app_io.add_argument('--probe_size',metavar='Size surface probe radius',type=float,default=1.4, help='Recommended: 1.4 (radius of a water molecule), Please refer to the paper for more details')
    entropy_app_io.add_argument('--onspherepoints',metavar='Number of points on a sphere',type=int,default=30, help='Recommended: 30. Number of points to be generated around each point for the surface (Read the Publication for more details)')

    # web server parameters
    web_server_group = parser.add_argument_group('Web server parameters', 'Used by the web form')
    web_server_group.add_argument('--outputfile', type=argparse.FileType('w', 1), default=sys.stdout, help='Path and filename write output to')
    web_server_group.add_argument('--logfile', type=argparse.FileType('w', 1), default=sys.stderr, help='Path and filename write log messages to')
    web_server_group.add_argument('--callbackurl', type=str, help='Optional callback url if this script was called from Drupal.')
    web_server_group.add_argument('--nodeid', type=int, help='Optional node id if this script was called from Drupal.')
    
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

    if(args.pdbid is not None):
        molecule.download_structure(args.pdbid, args.filename)

    try:
        extension = args.filename.split('.')[-1]
        mol = molecule.load_structure(args.filename,ftype=extension)
    except:
        logging.warning("The filename provided does not appear to have a format extension.")
        mol = molecule.load_structure(args.filename)
    
    if(args.command == 'hinge'):
        hinge_cli(args,mol)
    elif(args.command == 'entropy'):
        entropy_cli(args,mol)

    return True


'''
##################################################################################################
#                                            Main                                                #
##################################################################################################
'''


def main():
    """Gatewayto CLI and GUI
    """
    try:
        if(len(sys.argv) >  1 and sys.argv[1] == "gui"):
            from .GUI import load_gui
            load_gui()
        else:
            load_cli()
    except Exception as e:
        print(e)
        logging.error("Please provide a valid option. Enter 'python -m packman gui' for the GUI. Otherwise, please check the documentation for the CLI options. This function is for the CLI and not an integral function for the API.")



if(__name__ == '__main__'):
    main()