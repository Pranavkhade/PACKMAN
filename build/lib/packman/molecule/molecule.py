# -*- coding: utf-8 -*-
"""The 'Molecule' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Model' object documentation for details. [ help(packman.molecule) ]

Example:
    >>>from packman.molecule import molecule
    >>>help( Model )
    OR

    >>>from packman import molecule
    >>>help( molecule )
    
Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""

import numpy
import warnings

from .protein import Protein
from .model import Model

from .chain import Chain
from .residue import Residue
from .atom import Atom

from .hetmol import HetMol
from .hetatom import HetAtom


'''
##################################################################################################
#                                          File Load                                             #
##################################################################################################
'''

def load_structure(filename):
    """This class helps user to load the 3D structure of the protein onto a packman.molecule.Protein object.

    This file is basis for all the PACKMAN API. Loading file is the first thing user does.

    Example:

        To Load the file    
        >>>from packman import molecule
        >>>molecule.download_structure('1prw','1prw.pdb')
        >>>molecule.load_structure('1prw.pdb')
        <packman.molecule.protein.Protein object at Memory>

    Args:
        filename (str)          : Name of the input file
    
    Returns:
        packman.molecule.Protein: Protein object containing all the information about the Protein
    """
    Models=[]
    start=0

    fh=open(filename).read()

    frames=fh.split('\nMODEL')

    if(len(frames)>1):
        start=1
    else:
        start=0
    for FrameNumber,frame in enumerate(frames[start:]):
        #Map
        AllAtoms={}
        AllResidues={}
        AllChains={}

        AllHetAtoms={}
        AllHetMols={}

        lines=frame.split('\n')   
        for _ in lines:
            if(_[0:4]=='ATOM'):
                #NOTE: MAPS CAN BE REMOVED SAFELY
                #Chain Defined
                ChainID=_[21]
                if(ChainID not in AllChains.keys()):AllChains[ChainID]=Chain(ChainID)

                #Residue Defined
                ResidueNumber=int(_[22:26].strip())
                ResidueName=_[17:20]
                if(str(ResidueNumber)+ChainID not in AllResidues.keys()):AllResidues[str(ResidueNumber)+ChainID]=Residue(ResidueNumber,ResidueName,AllChains[ChainID])

                #Residue Added to the chain
                AllChains[ChainID].__setitem__(ResidueNumber,AllResidues[str(ResidueNumber)+ChainID],Type='Residue')

                #Atom Defined
                id=int(_[6:11])
                AtomName=_[12:16].strip()
                Coordinates=numpy.array([float(_[30:38]),float(_[38:46]),float(_[46:54])])
                Occupancy=float(_[54:60])
                bfactor=float(_[60:66])
                Element=_[76:78].strip()
                Charge=_[78:80]
                AllAtoms[id]=Atom(id,AtomName,Coordinates,Occupancy,bfactor,Element,Charge,AllResidues[str(ResidueNumber)+ChainID])
                
                #Atom added to the residue
                AllResidues[str(ResidueNumber)+ChainID].__setitem__(id,AllAtoms[id])

                #What to do with these?
                AlternateLocationIndicator=_[16]
                CodeForInsertions=_[26]
                SegmentIdentifier=_[72:76]
            elif(_[0:6]=='HETATM'):
                ChainID=_[21]
                if(ChainID not in AllChains.keys()):AllChains[ChainID]=Chain(ChainID)

                #HetMol Defined
                HetMolNumber=int(_[22:26].strip())
                HetMolName=_[17:20]
                if(str(HetMolNumber)+ChainID not in AllHetMols.keys()):AllHetMols[str(HetMolNumber)+ChainID]=HetMol(HetMolNumber,HetMolName,AllChains[ChainID])

                #HetMol Added to the chain
                AllChains[ChainID].__setitem__(HetMolNumber,AllHetMols[str(HetMolNumber)+ChainID],Type='HetMol')

                #HetAtom Defined
                id=int(_[6:11])
                AtomName=_[12:16].strip()
                Coordinates=numpy.array([float(_[30:38]),float(_[38:46]),float(_[46:54])])
                Occupancy=float(_[54:60])
                bfactor=float(_[60:66])
                Element=_[76:78].strip()
                Charge=_[78:80]
                AllHetAtoms[id]=HetAtom(id,AtomName,Coordinates,Occupancy,bfactor,Element,Charge,AllHetMols[str(HetMolNumber)+ChainID])
                
                #HetAtom added to the residue
                AllHetMols[str(HetMolNumber)+ChainID].__setitem__(id,AllHetAtoms[id])

                #What to do with these?
                AlternateLocationIndicator=_[16]
                CodeForInsertions=_[26]
                SegmentIdentifier=_[72:76]

        #print
        Models.append(Model(FrameNumber,AllAtoms,AllResidues,AllChains,AllHetAtoms,AllHetMols))
        for i in AllChains:AllChains[i].set_parent(Models[FrameNumber])

    if(len(Models)>1):
        #NMR
        warnings.warn('Multiple models/frames are detected (B-Factor field is turned to a calculated parameter)',UserWarning)
        All_Coords=[]
        for i in Models:
            All_Coords.append(numpy.array([j.get_location() for j in i.get_atoms()]))
        All_Coords=numpy.array(All_Coords)
        
        flattened_std=[]
        for i in range(0,All_Coords.shape[1]):
            xyz_var=0
            for j in All_Coords[:,i].T:
                xyz_var=xyz_var+numpy.var(j)
            flattened_std.append(numpy.sqrt(xyz_var))
        
        for i in Models:
            for numj,j in enumerate(i.get_atoms()):
                j.set_bfactor(flattened_std[numj])
        
    return Protein(filename,None,Models)


'''
##################################################################################################
#                                        Download                                                #
##################################################################################################
'''

def download_structure(pdbid,save_name):
    """This class helps user to download the 3D structure of the protein and save it on a disk.

    Example:

        To Download the file    
        >>>from packman import molecule
        >>>molecule.download_structure('1prw','1prw.pdb')

    Args:
        pdbid (str)    : A Unique 4 Letter PDB ID (eg.. 1PRW) 
        save_name (str): Save name of the downloaded file (Don't forget the extension). (eg... 1PRW.pdb)
    """
    import urllib.request as ur
    response=ur.urlopen('https://files.rcsb.org/view/'+pdbid+'.pdb')
    try:
        open(save_name,'wb').write(response.read())
    except(IOError):
        None
    return True