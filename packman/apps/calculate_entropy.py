# -*- coding: utf-8 -*-
"""The 'predict_hinge' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'predict_hinge' object documentation for details.

Example::

    from packman.apps import calculate_entropy
    help( calculate_entropy )

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""
from ..molecule import download_structure
from ..molecule import load_structure
from ..entropy import PackingEntropy

import argparse
import logging



def entropy_cli(args,mol):
    """Command-line Interface for the 'entropy' command. Please check the packman.bin.PACKMAN file for more details.

    This function is for the CLI and not an integral function for the API.

    Args:
        args (parser.parse_args())     : The arguments that were passed by the user to the PACKMAN-hinge app.
        mol (packman.molecule.Protein) : The 'Protein' object for the anaylsis.
    """
    try:
        input_chains = args.chains.split(',')
    except:
        input_chains = None
    
    if(args.type=='PackingEntropy'):
        args.outputfile.write('Frame\tChain\tResidueID\tResidueName\tPackingEntropy\n')
        for i in mol:
            try:
                result = PackingEntropy(i.get_atoms(),chains=input_chains,probe_size=args.probe_size,onspherepoints=args.onspherepoints)
            except:
                logging.error('Please check the parameters.')
                exit()
            
            for j in i.get_residues():
                try:
                    args.outputfile.write( str(i.get_id())+'\t'+str(j.get_parent().get_id())+'\t'+str(j.get_id())+'\t'+str(j.get_name())+'\t'+str(j.get_entropy('PackingEntropy'))+'\n' )
                except:
                    None
        try:
            args.outputfile.write('Total (Frame '+str(i.get_id())+')\t'+str(result.get_total_entropy())+'\n' )
            args.outputfile.write('Total Normalized Entropy (3N-6) (Frame '+str(i.get_id())+')\t'+str( result.get_total_entropy() / ( 3*len([j for j in i.get_atoms()]) -6 ) ) +'\n' )

        except:
            None
        args.outputfile.flush()
        args.outputfile.close()
    else:
        logging.error('Please provide a valid entropy type (-type option). (Options: 1. PackingEntropy )')
    return True
