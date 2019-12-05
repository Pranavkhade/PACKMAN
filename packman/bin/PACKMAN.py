#!/usr/bin/env python
'''
Author: Pranav Khade(pranavk@iastate.edu)
'''

import traceback

from .. import molecule
from ..apps import predict_hinge

import operator
import argparse
import os
import sys

import urllib


'''
##################################################################################################
#                                    Non Algorithm Functions                                     #
##################################################################################################
'''

def WriteOBJ(atoms,faces, fh):
    """
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
    """
    INFO: Argument parser to the program for now.
    """
    parser=argparse.ArgumentParser(description='PACKMAN: PACKing and Motion ANalysis. (https://github.com/Pranavkhade/PACKMAN)')

    parser.add_argument('-pdbid','--pdbid', metavar='PDB_ID', type=str, help='If provided, the PBD with this ID will be downloaded and saved to FILENAME.')

    parser.add_argument('alpha', metavar='AlphaValue', help='Recommended: 2.8 for closed; 4.5 for open form, Please refer to the paper for more details')
    parser.add_argument('filename', metavar='FILENAME', help='Path and filename of the PDB file.')

    parser.add_argument('--e_clusters',metavar='NumberOfEccentricityClusters',type=int,default=4,help='Recommended: 4, Please refer to the paper for more details')
    parser.add_argument('--minhnglen',metavar='MinimumHingeLength',type=int,default=5,help='Recommended: 5, Please refer to the paper for more details')
    parser.add_argument("--chain", help='Enter The Chain ID')
    parser.add_argument('--generateobj', type=argparse.FileType('wb', 0), help='Path and filename to save the .obj file at. Ignored unless --chain is provided.')

    # web server parameters
    web_server_group = parser.add_argument_group('Web server parameters', 'Used by the web form')
    web_server_group.add_argument('--outputfile', type=argparse.FileType('wb', 0), default=sys.stdout, help='Path and filename write output to')
    web_server_group.add_argument('--logfile', type=argparse.FileType('wb', 0), default=sys.stderr, help='Path and filename write log messages to')
    web_server_group.add_argument('--callbackurl', type=str, help='Optional callback url if this script was called from Drupal.')
    web_server_group.add_argument('--nodeid', type=int, help='Optional node id if this script was called from Drupal.')
    
    args=parser.parse_args()
    return args

'''
##################################################################################################
#                                             Main                                               #
##################################################################################################
'''

def main():
    """
    """
    args=IO()

    if(args.pdbid is not None):
        molecule.download_structure(args.pdbid, args.filename)

    mol = molecule.load_structure(args.filename)

    try:
        if(args.chain):
            Backbone = [item for sublist in mol[0][args.chain].get_backbone() for item in sublist]
            SelectedTesselations= predict_hinge(Backbone, args.outputfile, Alpha=float(args.alpha), filename=args.filename,nclusters=args.e_clusters,MinimumHingeLength=args.minhnglen)
            if(args.generateobj is not None):
                WriteOBJ(Backbone, SelectedTesselations, args.generateobj)
        else:
            for i in mol[0].get_chains():
                Backbone = [item for sublist in mol[0][i.get_id()].get_backbone() for item in sublist]
                SelectedTesselations = predict_hinge(Backbone, args.outputfile, Alpha=float(args.alpha), filename=args.filename,nclusters=args.e_clusters,MinimumHingeLength=args.minhnglen)
        
        if args.nodeid is not None:
            urllib.urlopen(args.callbackurl + '/' + str(args.nodeid) + "/0")

    except Exception:
        print(traceback.print_exc())
        if args.nodeid is not None:
            urllib.urlopen(args.callbackurl + '/' + str(args.nodeid) + "/1")
    
    finally:
        print_footnotes(args.outputfile)
        args.outputfile.flush()
        args.logfile.flush()

        if(args.generateobj is not None):
            args.generateobj.flush()
    return True

'''
##################################################################################################
#                                             End                                                #
##################################################################################################
'''
def print_footnotes(outputfile):
    outputfile.write('Footnotes:\n\nSTATISTICS Section Legend:\nN: Number of residues\nMin: Minimum B-factor value\nMax: Maximum B-factor value\nMean: Mean B-factor value\nMode: Mode B-factor value\nMedian: Median B-factor value\nSTDDev: Standard Deviation of the B-factor values\n')
    outputfile.write('\nRMSF Plane Visualization:\nDownload plane.py file from https://pymolwiki.org/index.php/Plane_Wizard and place it in the PyMol working directory if the RMSF plane needs to be visualized.\n')
    return True



if(__name__=='__main__'):
    main()
