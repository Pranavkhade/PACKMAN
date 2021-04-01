# -*- coding: utf-8 -*-
"""The 'Molecule' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Model' object documentation for details.

Example::

    from packman.molecule import Model
    help( Model )

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

def load_pdb(filename):
    """
    Load the PDB (.pdb) file into the 'Protein' Object.
    """

    Models=[]
    AllAnnotations = []
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
                AtomID=int(_[6:11])
                AtomName=_[12:16].strip()
                Coordinates=numpy.array([float(_[30:38]),float(_[38:46]),float(_[46:54])])
                Occupancy=float(_[54:60])
                bfactor=float(_[60:66])
                Element=_[76:78].strip()
                Charge=_[78:80]
                AllHetAtoms[AtomID]=Atom(AtomID,AtomName,Coordinates,Occupancy,bfactor,Element,Charge,AllHetMols[str(HetMolNumber)+ChainID])
                
                #HetAtom added to the residue
                AllHetMols[str(HetMolNumber)+ChainID].__setitem__(AtomID,AllHetAtoms[AtomID])

                #What to do with these?
                AlternateLocationIndicator=_[16]
                CodeForInsertions=_[26]
                SegmentIdentifier=_[72:76]
            else:
                AllAnnotations.append(_)

        Models.append(Model(FrameNumber,AllAtoms,AllResidues,AllChains,AllHetAtoms,AllHetMols))
        for i in AllChains:AllChains[i].set_parent(Models[FrameNumber])

    if(len(Models)>2):
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
        
    prot = Protein(filename,Models)
    prot.set_data(AllAnnotations)
    #Setting parent to the model object
    for i in prot: i.set_parent(prot)
    return prot


def load_cif(filename):
    """
    Load the CIF (.cif) file into the 'Protein' Object.
    """

    #Global Variables
    AllChains = []
    AllResidues = []
    AllAtoms = []

    AllHetMols = []
    AllHetAtoms =[]

    AllAnnotations = []

    sections = open(filename,'r').read().split('loop_\n')
    for section in sections:
        for subsection in section.split('#'):
            column_indices = {}
            column_names = {}
            for n_line, line in enumerate(subsection.split('\n')):
                try:
                    _ = line.split()
                    
                    # When the data is next to the columm description
                    if(line[0] == '_' and len(_) > 1):
                        AllAnnotations.append(line.strip())

                    # When the data columns are below the description section
                    elif(line[0] == '_' and len(_) == 1):
                        AllAnnotations.append(line.strip())
                        column_indices[n_line] = line.strip()
                        column_names[line.strip()] = n_line

                    else:
                        #All the molecule objects are built here
                        if(_[0] == 'ATOM'):
                            #Initiate Model Number
                            FrameNumber = int(_[ column_names['_atom_site.pdbx_PDB_model_num'] ])
                            
                            #Initiate Chain (Frame number, Chain)
                            ChainID = _[ column_names['_atom_site.label_asym_id'] ]
                            try:
                                if( ChainID not in AllChains[FrameNumber-1].keys() ): AllChains[FrameNumber-1][ChainID] = Chain(ChainID)
                            except:
                                AllChains.append( {} )
                                AllChains[FrameNumber-1][ChainID] = Chain(ChainID)
                            
                            #Initiate Residue
                            ResidueNumber = int( _[ column_names['_atom_site.label_seq_id'] ] )
                            ResidueName   = _[ column_names['_atom_site.label_comp_id'] ]
                            try:
                                if( str(ResidueNumber)+ChainID not in AllResidues[FrameNumber-1].keys() ): AllResidues[FrameNumber-1][str(ResidueNumber)+ChainID] = Residue( ResidueNumber, ResidueName, AllChains[FrameNumber-1][ChainID] )
                            except:
                                AllResidues.append( {} )
                                AllResidues[FrameNumber-1][str(ResidueNumber)+ChainID] = Residue( ResidueNumber, ResidueName, AllChains[FrameNumber-1][ChainID] )
                            
                            #Initiate Atom
                            AtomID = int(_[ column_names['_atom_site.id'] ])
                            AtomName = _[ column_names['_atom_site.label_atom_id'] ]
                            Coordinates = numpy.array( [ float(_[ column_names['_atom_site.Cartn_x'] ]), float(_[ column_names['_atom_site.Cartn_y'] ]), float(_[ column_names['_atom_site.Cartn_z'] ]) ] )
                            Occupancy = float(_[ column_names['_atom_site.occupancy'] ])
                            bfactor = float(_[ column_names['_atom_site.B_iso_or_equiv'] ])
                            Element = _[ column_names['_atom_site.type_symbol'] ]
                            Charge = _[ column_names['_atom_site.pdbx_formal_charge'] ]
                            try:
                                AllAtoms[FrameNumber-1][AtomID] = Atom(AtomID,AtomName,Coordinates,Occupancy,bfactor,Element,Charge, AllResidues[FrameNumber-1][str(ResidueNumber)+ChainID] )
                            except:
                                AllAtoms.append( {} )
                                AllAtoms[FrameNumber-1][AtomID] = Atom(AtomID,AtomName,Coordinates,Occupancy,bfactor,Element,Charge, AllResidues[FrameNumber-1][str(ResidueNumber)+ChainID] )
                            
                            #Connect objects to each other to create a tree (Ideally should not happen every iteration)
                            #Residue Added to the chain
                            AllChains[FrameNumber-1][ChainID].__setitem__( ResidueNumber, AllResidues[FrameNumber-1][str(ResidueNumber)+ChainID], Type='Residue' )
                            #Atom added to the residue
                            AllResidues[FrameNumber-1][str(ResidueNumber)+ChainID].__setitem__( AtomID, AllAtoms[FrameNumber-1][AtomID] )

                        #Heteroatoms
                        elif(_[0] == 'HETATM'):
                            #Initiate Model Number
                            FrameNumber = int(_[ column_names['_atom_site.pdbx_PDB_model_num'] ])
                            
                            #Initiate Chain (Frame number, Chain)
                            ChainID = _[ column_names['_atom_site.label_asym_id'] ]
                            try:
                                if( ChainID not in AllChains[FrameNumber-1].keys() ): AllChains[FrameNumber-1][ChainID] = Chain(ChainID)
                            except:
                                AllChains.append( {} )
                                AllChains[FrameNumber-1][ChainID] = Chain(ChainID)
                            
                            
                            #Initiate HetMol
                            #print('All ok till here', _[ column_names['_atom_site.label_seq_id'] ])
                            try:
                                HetMolNumber = int( _[ column_names['_atom_site.label_seq_id'] ] )
                            except:
                                HetMolNumber = 'HOH'+str(n_line)

                            HetMolName   = _[ column_names['_atom_site.label_comp_id'] ]
                            try:
                                if( str(HetMolNumber)+ChainID not in AllHetMols[FrameNumber-1].keys() ): AllHetMols[FrameNumber-1][str(HetMolNumber)+ChainID] = HetMol( HetMolNumber, HetMolName, AllChains[FrameNumber-1][ChainID] )
                            except:
                                AllHetMols.append( {} )
                                AllHetMols[FrameNumber-1][str(HetMolNumber)+ChainID] = Residue( HetMolNumber, HetMolName, AllChains[FrameNumber-1][ChainID] )
                            
                            #Initiate Atom
                            AtomID = int(_[ column_names['_atom_site.id'] ])
                            AtomName = _[ column_names['_atom_site.label_atom_id'] ]
                            Coordinates = numpy.array( [ float(_[ column_names['_atom_site.Cartn_x'] ]), float(_[ column_names['_atom_site.Cartn_y'] ]), float(_[ column_names['_atom_site.Cartn_z'] ]) ] )
                            Occupancy = float(_[ column_names['_atom_site.occupancy'] ])
                            bfactor = float(_[ column_names['_atom_site.B_iso_or_equiv'] ])
                            Element = _[ column_names['_atom_site.type_symbol'] ]
                            Charge = _[ column_names['_atom_site.pdbx_formal_charge'] ]

                            try:
                                AllHetAtoms[FrameNumber-1][AtomID] = Atom(AtomID,AtomName,Coordinates,Occupancy,bfactor,Element,Charge, AllHetMols[FrameNumber-1][str(HetMolNumber)+ChainID] )
                            except:
                                AllHetAtoms.append( {} )
                                AllHetAtoms[FrameNumber-1][AtomID] = Atom(AtomID,AtomName,Coordinates,Occupancy,bfactor,Element,Charge, AllHetMols[FrameNumber-1][str(HetMolNumber)+ChainID] )

                            #Connect objects to each other to create a tree (Ideally should not happen every iteration)
                            #HetMol Added to the chain
                            AllChains[FrameNumber-1][ChainID].__setitem__( HetMolNumber, AllHetAtoms[FrameNumber-1][str(HetMolNumber)+ChainID], Type='HetMol' )
                            #Atom added to the HetMol
                            AllHetAtoms[FrameNumber-1][str(HetMolNumber)+ChainID].__setitem__( AtomID, AllAtoms[FrameNumber-1][AtomID] )

                        #Annotations additions
                        else:
                            AllAnnotations.append( line.strip() )
                            #This is when we will develop more annotation objects in future
                            #for n_cols, cols in enumerate( line.split() ):
                            #    try:
                            #        subsection_data[ column_indices[n_cols] ].append(cols)
                            #    except:
                            #        subsection_data[ column_indices[n_cols] ] = []
                            #        subsection_data[ column_indices[n_cols] ].append(cols)
                except:
                    None
            if(AllAnnotations[-1]!='#'): AllAnnotations.append('#')
        AllAnnotations.append('loop_')
            
    
        
    AllModels =[]
    for i in range(0,len(AllAtoms)):
        #In case protein part is missing
        try:
            current_atoms    = AllAtoms[i]
            current_residues = AllResidues[i]
            current_chains   = AllChains[i]
        except:
            current_atoms    = [None]
            current_residues = [None]
            current_chains   = [None]

        #In case the hetatoms are not present at all
        try:
            current_hetatms = AllHetAtoms[i]
            current_hetmols = AllHetMols[i]
        except:
            current_hetatms = [None]
            current_hetmols = [None]

        AllModels.append( Model(i+1, current_atoms, current_residues, current_chains, current_hetatms, current_hetmols) )
        for j in AllModels[i].get_chains(): j.set_parent(AllModels[i])
    

    if(len(AllModels)>2):
        #NMR
        warnings.warn('Multiple models/frames are detected (B-Factor field is turned to a calculated parameter)',UserWarning)
        All_Coords=[]
        for i in AllModels:
            All_Coords.append(numpy.array([j.get_location() for j in i.get_atoms()]))
        All_Coords=numpy.array(All_Coords)
        
        flattened_std=[]
        for i in range(0,All_Coords.shape[1]):
            xyz_var=0
            for j in All_Coords[:,i].T:
                xyz_var=xyz_var+numpy.var(j)
            flattened_std.append(numpy.sqrt(xyz_var))
        
        for i in AllModels:
            for numj,j in enumerate(i.get_atoms()):
                j.set_bfactor(flattened_std[numj])

    prot = Protein( filename, AllModels )
    prot.set_data(AllAnnotations)
    #Setting parent to the model object
    for i in prot: i.set_parent(prot)
    return prot


'''
##################################################################################################
#                                           Entry                                                #
##################################################################################################
'''


def load_structure(filename, ftype = 'cif'):
    """Load a Molecule from a file.

    This class helps user to load the 3D structure of the protein onto a packman.molecule.Protein object.

    Example::

        from packman import molecule
        molecule.download_structure('1prw')
        molecule.load_structure('1prw.cif')

    Args:
        filename (str)          : Name of the input file
        ftype    (str)          : Format name ('cif' or 'pdb'); Default: cif
    
    Returns:
        packman.molecule.Protein: Protein object containing all the information about the Protein
    """
    try:
        possible_ftype = filename.split('.')[1]
        if(possible_ftype == 'cif' or possible_ftype == 'pdb'):
            ftype = possible_ftype
    except:
        None

    if(ftype == 'cif'):
        return load_cif(filename)
    elif(ftype == 'pdb'):
        return load_pdb(filename)
    else:
        print('Please provide appropriate "ftype" argument. (cif/pdb).')


'''
##################################################################################################
#                                        Download                                                #
##################################################################################################
'''

def download_structure(pdbid,save_name=None,ftype='cif'):
    """This class helps user to download the 3D structure of the protein and save it on a disk.

    Example::

        from packman import molecule
        molecule.download_structure('1prw')

    Args:
        pdbid     (str) : A Unique 4 Letter PDB ID (eg.. 1PRW) 
        save_name (str) : Save name of the downloaded file (extension will be added automatically depending on the ftype argument).
        ftype     (str) : Format name ('.cif' or '.pdb')
    """
    import urllib.request as ur

    if(ftype=='cif'):
        response=ur.urlopen('https://files.rcsb.org/view/'+pdbid+'.cif')
    elif(ftype=='pdb'):
        response=ur.urlopen('https://files.rcsb.org/view/'+pdbid+'.pdb')
    else:
        print('Please provide appropriate "ftype" argument. (cif/pdb).')

    if(save_name==None):
        try:
            open(pdbid+'.'+ftype,'wb').write(response.read())
        except(IOError):
            None
    else:
        try:
            open(save_name+'.'+ftype,'wb').write(response.read())
        except(IOError):
            None
    return True