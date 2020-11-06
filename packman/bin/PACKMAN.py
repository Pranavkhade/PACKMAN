#!/usr/bin/env python
"""The PACKMAN Command-line Interface (CLI) and Graphical User Interface (GUI) host file.

How to use:
python -m packman gui (For GUI)
python -m packman     (For CLI)

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""
import numpy

import traceback

from .. import molecule
from ..apps import predict_hinge
from ..anm import hdANM

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
    web_server_group.add_argument('--outputfile', type=argparse.FileType('w', 1), default=sys.stdout, help='Path and filename write output to')
    web_server_group.add_argument('--logfile', type=argparse.FileType('w', 1), default=sys.stderr, help='Path and filename write log messages to')
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
    """Main Function for PACKMAN.

    Kernel of the PACKMAN interface.

    Todo:
        * Change the main() to accomodate iterative alpha shape search
    
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
            urlopen(args.callbackurl + '/' + str(args.nodeid) + "/0")

    except Exception:
        print(traceback.print_exc())
        if args.nodeid is not None:
            urlopen(args.callbackurl + '/' + str(args.nodeid) + "/1")
    
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
    """Add footnotes to the output file.
    
    Args:
        outputfile (file): The file to which the footnotes will be written.
    """
    outputfile.write('Footnotes:\n\nSTATISTICS Section Legend:\nN: Number of residues\nMin: Minimum B-factor value\nMax: Maximum B-factor value\nMean: Mean B-factor value\nMode: Mode B-factor value\nMedian: Median B-factor value\nSTDDev: Standard Deviation of the B-factor values\n')
    outputfile.write('\nRMSF Plane Visualization:\nDownload plane.py file from https://pymolwiki.org/index.php/Plane_Wizard and place it in the PyMol working directory if the RMSF plane needs to be visualized.\n')
    return True



'''
##################################################################################################
#                                             GUI                                                #
##################################################################################################
'''

try:
    import tkinter as tk
    from tkinter.messagebox import showinfo, showerror
    from tkinter import N, S, E, W, Grid, END, BooleanVar, StringVar, DISABLED
except:
    print("The Python version doesn't the (built in) tkinter package. Please make sure you are using right Python interpreter. ( try: python3 -m packman gui )")


class Skeleton(tk.Tk):
    def __init__(self,frame_name):
        tk.Tk.__init__(self)

        #Global Variable
        numberofcolumns = 3
        self.title_content = ''

        for i in range(0,numberofcolumns):
            #Grid.rowconfigure(self, i+1, weight=1)
            Grid.columnconfigure(self, i, weight=1)
        #All the frames
        self.frames = {}
        self.frames['HomePage'] = HomePage(self)
        self.frames['HingePrediction'] = HingePrediction(self)
        self.frames['hdANM_GUI'] = hdANM_GUI(self)
        self.frames['TopMenu'] = top_menu(self)

        self.geometry("800x600")

        #Icon for windows
        try:
            self.iconbitmap(False,'logo.ico')
        except:
            None
        
        self.show_frame(frame_name)

        try:
            self.title_content = "GUI for PACKMAN Version: "+packman.__version__
            self.title( self.title_content )
        except:
            self.title_content = "PACKMAN GUI"
            self.title( self.title_content )
            showinfo('Notification', 'PACKMAN was not loaded properly or has older version than 1.1.3\n\nTry:\n\npython -n pip install packman\npython -n pip install packman --upgrade')
            self.focus_force()
        
    
    def show_frame(self,frame_name):
        for i in self.frames:
            if(i!='TopMenu'):self.frames[i].hide()
        self.frames[frame_name].show()

class top_menu(tk.Frame):
    def __init__(self,parent):
        tk.Frame.__init__(self)
        tk.Button(parent, text = 'Home',             command = lambda: parent.show_frame('HomePage')        ).grid(row=0,column=0,sticky=E+W)
        tk.Button(parent, text = 'Hinge Prediction', command = lambda: parent.show_frame('HingePrediction') ).grid(row=0,column=1,sticky=E+W)
        tk.Button(parent, text = 'hdANM',            command = lambda: parent.show_frame('hdANM_GUI')       ).grid(row=0,column=2,sticky=E+W)



'''
##################################################################################################
#                                               Home                                             #
##################################################################################################
'''



class HomePage(tk.Frame):
    def __init__(self,parent):
        tk.Frame.__init__(self)
        self.Label1 = tk.Label(parent, text="PACKMAN")
        self.Label2 = tk.Label(parent, text="PACKMAN is a multi-utility tool to study protein packing and its effect on protein dynamics.\n\nImportant Links:")
        
        #Input
        self.Text1 = tk.Text(parent, height=10)
        self.Text1.insert(END,"1. Tutorials: https://py-packman.readthedocs.io/en/latest/index.html\n2. GitHub: https://github.com/Pranavkhade/PACKMAN\n3. Online Hinge Prediction Server: https://packman.bb.iastate.edu/\n4. Online hdANM Server: https://hdanm.bb.iastate.edu/\n5. Reference for Hinge Prediction: https://doi.org/10.1016/j.jmb.2019.11.018\n6. Reference for Compliance: https://doi.org/10.1002/prot.25968\n7. Reference for hd-ANM: Coming Soon")
        self.Text1.config(state=DISABLED)

        self.all_objects = self.__dict__

    def show(self):
        self.Label1.grid(row=1,columnspan=3)
        self.Label2.grid(row=2,columnspan=3)
        self.Text1.grid(row=3,columnspan=3)
    
    def hide(self):
        for i in self.all_objects:
            try:
                self.all_objects[i].grid_forget()
            except:
                continue

'''
##################################################################################################
#                                      Hinge Prediction                                          #
##################################################################################################
'''

class HingePrediction(tk.Frame):
    def __init__(self,parent):
        tk.Frame.__init__(self)
        self.grid(padx=50)
        self.Text1 =  tk.Text(parent, height=10)
        self.Label1 = tk.Label(parent, text="Hinge Prediction")
        self.Label2 = tk.Label(parent, text="Note: Keep increasing the Alpha Value if hinges do not show up in the results.", borderwidth=2, relief="groove")
        self.Label3 = tk.Label(parent, text="PDB ID:")
        self.Label4 = tk.Label(parent, text="Chain ID:")
        self.Label5 = tk.Label(parent, text="Alpha Value:")
        self.Label6 = tk.Label(parent, text="Number of Eccentricity Clusters:")
        self.Label7 = tk.Label(parent, text="Minimum Hinge Length:")

        #Input
        self.Text1.insert(END,"Interface to the functionality to identify the protein hinges (separating the domains) using PACKMAN. It can also be used to read, write, manipulate and analyze protein molecules and it's properties through its API.\n\nCitation:\n\nKhade PM, Kumar A, Jernigan RL. Characterizing and Predicting Protein Hinges for Mechanistic\nInsight. J Mol Biol. November 2019. doi:10.1016/j.jmb.2019.11.018" )
        self.Text1.config(state=DISABLED)
        self.Box1 = tk.Entry(parent)
        self.Box2 = tk.Entry(parent)
        self.Box3 = tk.Entry(parent)
        self.Box4 = tk.Entry(parent)
        self.Box5 = tk.Entry(parent)

        self.Box3.insert(END, '2.8')
        self.Box4.insert(END, '4')
        self.Box5.insert(END, '5')

        self.Button1 = tk.Button(parent, text = "Run", command = lambda: self.run_PACKMAN() )

        self.all_objects = self.__dict__

    def show(self):
        self.Text1.grid(row=2,columnspan=3,sticky=E+W, padx=10, pady=10 )
        self.Label1.grid(row=1,columnspan=3,sticky=E+W, padx=10, pady=10 )
        self.Label2.grid(row=3,columnspan=3,sticky=E+W, padx=10, pady=10 )
        self.Label3.grid(row=4,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label4.grid(row=5,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label5.grid(row=6,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label6.grid(row=7,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label7.grid(row=8,column=0,sticky=N+S+E+W, padx=10, pady=10 )

        self.Box1.grid(row=4,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        self.Box2.grid(row=5,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        self.Box3.grid(row=6,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        self.Box4.grid(row=7,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        self.Box5.grid(row=8,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        
        self.Button1.grid(row=9, column=2)
    
    def hide(self):
        for i in self.all_objects:
            try:
                self.all_objects[i].grid_forget()
            except:
                continue
    
    def run_PACKMAN(self):
        """
        PACKMAN Hinge Prediction Algorithm
        """
        try:
            mol = molecule.load_structure(self.Box1.get()+'.pdb')
        except:
            showinfo('Notification','File on found in the folder, downloading it from the PDB.')
            
            #For file availability
            try:
                molecule.download_structure(self.Box1.get(), self.Box1.get()+'.pdb')
            except:
                showerror('Error','Please check the internet connection or the validity of the PDB ID')
                exit()
        mol = molecule.load_structure(self.Box1.get()+'.pdb')

        #Other parameter validity
        try:
            alpha = float(self.Box3.get())
            ecc_clust = int(self.Box4.get())
            min_length = int(self.Box5.get())
        except:
            showerror('Error','Check the value of the parameters:\n\nAlpha Value(float)\nNumber of Eccentricity Clusters (int)\nMinimum Hinge Length (int)')

        #For checking if the format is right and chain id is right
        try:
            backbone = [j for i in mol[0][self.Box2.get()].get_backbone() for j in i if j is not None]
            predict_hinge(backbone, Alpha=alpha, nclusters=ecc_clust, MinimumHingeLength=min_length, outputfile=open('output.txt', 'w'))
        except:
            showinfo('Notification','Chain ID is either not specified or not valid. Running PACKMAN on all the chains.')
            try:
                for chain in mol[0].get_chains():
                    backbone = [j for i in mol[0][chain.get_id()].get_backbone() for j in i if j is not None]
                    predict_hinge(backbone, Alpha=alpha, nclusters=ecc_clust, MinimumHingeLength=min_length, outputfile=open('output.txt', 'w'))

            except:    
                showerror('Error','Possible Solutions:\n1. Provide valid chain ID or check the PDB file format.\n2. Try changing Alpha Value')
                exit()
        
        #At this point, everything is okay
        
        
        #Results
        pop_up1 = tk.Toplevel()
        pop_up1.title('PACKMAN Results')

        tk.Label(pop_up1, text = 'NOTE: Please check in front of the hinges you want to keep and click on the "Save HNG file" button to generate input file for hdANM. Please provide a valid filename\n Tip: p-value ideally should be <0.05 ').grid(row=0,columnspan=5)
        tk.Label(pop_up1, text = 'Hinge ID').grid(row=1,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(pop_up1, text = 'PDB ID').grid(row=1,column=1,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(pop_up1, text = 'Chain ID').grid(row=1,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(pop_up1, text = 'Residue IDs').grid(row=1,column=3,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(pop_up1, text = 'p-value').grid(row=1,column=4,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(pop_up1, text = 'Select').grid(row=1,column=5,sticky=N+S+E+W, padx=10, pady=10 )

        All_Hinges = []
        if(self.Box2.get()!=''):
            All_Hinges = mol[0][self.Box2.get()].get_hinges()
        else:
            for i in mol[0].get_chains():
                All_Hinges.extend( mol[0][i.get_id()].get_hinges() )


        Selected_Hinges = [False]*len(All_Hinges)
        for numi, i in enumerate(All_Hinges):
            Residue_IDs = sorted([j.get_id() for j in i.get_elements()])
            ChainOfHinge = i.get_elements()[0].get_parent().get_id()

            tk.Label(pop_up1, text = numi+1).grid(row=numi+2,column=0,sticky=N+S+E+W, padx=10, pady=10 )
            tk.Label(pop_up1, text = self.Box1.get() ).grid(row=numi+2,column=1,sticky=N+S+E+W, padx=10, pady=10 )
            tk.Label(pop_up1, text = ChainOfHinge    ).grid(row=numi+2,column=2,sticky=N+S+E+W, padx=10, pady=10 )
            tk.Label(pop_up1, text =  str(Residue_IDs[0])+':'+str(Residue_IDs[-1]) ).grid(row=numi+2,column=3,sticky=N+S+E+W, padx=10, pady=10 )
            tk.Label(pop_up1, text = i.get_pvalue()).grid(row=numi+2,column=4,sticky=N+S+E+W, padx=10, pady=10 )

            Selected_Hinges[numi] = BooleanVar(pop_up1)
            Check = tk.Checkbutton(pop_up1, variable=Selected_Hinges[numi], onvalue=True, offvalue=False)
            Check.grid(row=numi+2,column=5,sticky=N+S+E+W, padx=10, pady=10 )

            #Details Button
            tk.Button(pop_up1,text = 'Details', command= lambda i=i:show_details( i.get_stats() ) ).grid(row=numi+2,column=6,sticky=N+S+E+W, padx=10, pady=10 )
            
        
        tk.Label(pop_up1, text = 'output filename (.hng)').grid(row=len(All_Hinges)+2,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        self.output_filename = tk.Entry(pop_up1)
        tk.Button(pop_up1,text = 'Save HNG file', command=lambda: write_hng_file() ).grid(row=len(All_Hinges)+2,column=4,sticky=N+S+E+W, padx=10, pady=10 )

        self.output_filename.insert(END,self.Box1.get()+'.hng')
        self.output_filename.grid(row=len(All_Hinges)+2,column=1,columnspan=2,sticky=N+S+E+W, padx=10, pady=10 )

        def write_hng_file():
            ALL_RESIDUES = {}
            for i in mol[0].get_chains():
                try:
                    ALL_RESIDUES[i.get_id()] = sorted([i.get_id() for i in mol[0][i.get_id()].get_residues() if i!=None])
                except:
                    None

            select_count = 0
            last_hinge_end = 0
            fh = open(self.output_filename.get(), 'w')
            for numi, i in enumerate(Selected_Hinges):
                current_hinge = All_Hinges[numi]
                ChainOfHinge = current_hinge.get_elements()[0].get_parent().get_id()

                if( i.get() and select_count==0 ):
                    hinge_res_ids = sorted([j.get_id() for j in current_hinge.get_elements()])
                    
                    select_count += 1
                    if(ALL_RESIDUES[ChainOfHinge][0]!=hinge_res_ids[0]):
                        fh.write(self.Box1.get()+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(ALL_RESIDUES[ChainOfHinge][0])+':'+str(hinge_res_ids[0]-1)+'\n' )
                        fh.write(self.Box1.get()+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
                    else:
                        fh.write(self.Box1.get()+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
                    last_hinge_end = hinge_res_ids[-1]
                elif( i.get() ):
                    hinge_res_ids = sorted([j.get_id() for j in current_hinge.get_elements()])
                    select_count += 1
                    fh.write(self.Box1.get()+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(last_hinge_end+1)+':'+str(hinge_res_ids[0]-1)+'\n' )
                    fh.write(self.Box1.get()+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
                    last_hinge_end = hinge_res_ids[-1]
                try:
                    if(ChainOfHinge != All_Hinges[numi+1].get_elements()[0].get_parent().get_id() ):
                        select_count += 1
                        fh.write(self.Box1.get()+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(last_hinge_end)+':'+str(ALL_RESIDUES[ChainOfHinge][-1])+'\n' )
                        last_hinge_end = 0
                except:
                    None

            if(last_hinge_end != ALL_RESIDUES[ChainOfHinge][-1]):
                select_count += 1
                fh.write(self.Box1.get()+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(last_hinge_end)+':'+str(ALL_RESIDUES[ChainOfHinge][-1])+'\n' )
            fh.flush()
            fh.close()
            showinfo('Notification',self.output_filename.get()+' file saved!')
        
        def show_details(content):
            pop_up1_2 = tk.Toplevel()
            pop_up1_2.title('Hinge Details')

            for numi, i in enumerate(content):
                for numj,j in  enumerate(i):
                    tk.Label(pop_up1_2, text = j).grid(row=numi+2,column=numj,sticky=N+S+E+W, padx=10, pady=10 )

            return True

        #on_click() Ends here
        return True
    
    


'''
##################################################################################################
#                                            hd-ANM                                              #
##################################################################################################
'''


class hdANM_GUI(tk.Frame):
    def __init__(self,parent):
        #Varible
        self.mass_type = StringVar(parent)
        self.mass_type.set("Molecular Weight") # default value

        #
        tk.Frame.__init__(self)
        self.Text1 =  tk.Text(parent, height=10)
        self.Label1 = tk.Label(parent, text="hdANM")
        self.Label2 = tk.Label(parent, text="Note: .hng file can be obtained from the 'Hinge Prediction' section. ", borderwidth=2, relief="groove")

        self.Label3 = tk.Label(parent, text="PDB ID:")
        self.Label4 = tk.Label(parent, text="Chain ID:")
        
        self.Label5 = tk.Label(parent, text="Cutoff Value")
        self.Label6 = tk.Label(parent, text="Power of Distance")
        self.Label7 = tk.Label(parent, text="Mass of the residue")
        self.Label7 = tk.Label(parent, text="Input (.hng) Filename")
        self.Label8 = tk.Label(parent, text="Mass of the residue")

        #Input
        self.Text1.insert(END,"hdANM is a comprehensive Elastic Network Model. Please read the paper for more details.\n\nCitation:\n\nPaper Submitted." )
        self.Text1.config(state=DISABLED)
        
        self.Box1 = tk.Entry(parent)
        self.Box2 = tk.Entry(parent)
        self.Box3 = tk.Entry(parent)
        self.Box4 = tk.Entry(parent)
        self.Box5 = tk.Entry(parent)

        self.Box3.insert(END, '15')
        self.Box4.insert(END, '0')
        self.Box5.insert(END, 'output.hng')
        self.dropdown1 = tk.OptionMenu(parent, self.mass_type, "Molecular Weight", "Unit", "Atomics Mass")

        self.Button1 = tk.Button(parent, text = "Run", command = lambda: self.run_hdANM() )

        self.all_objects = self.__dict__
    
    def show(self):
        self.Text1.grid(row=2,columnspan=3,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label1.grid(row=1,columnspan=3,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label2.grid(row=3,columnspan=3,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label3.grid(row=4,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label4.grid(row=5,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label5.grid(row=6,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label6.grid(row=7,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label7.grid(row=8,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label8.grid(row=9,column=0,sticky=N+S+E+W, padx=10, pady=10 )


        self.Box1.grid(row=4,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        self.Box2.grid(row=5,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        self.Box3.grid(row=6,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        self.Box4.grid(row=7,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        self.Box5.grid(row=8,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        self.dropdown1.grid(row=9,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        
        self.Button1.grid(row=10, column=2)
    
    
    def hide(self):
        for i in self.all_objects:
            try:
                self.all_objects[i].grid_forget()
            except:
                continue
    
    def run_hdANM(self):
        try:
            mol = molecule.load_structure(self.Box1.get()+'.pdb')
        except:
            showinfo('Notification','File on found in the folder, downloading it from the PDB.')
            
            #For file availability
            try:
                molecule.download_structure(self.Box1.get(), self.Box1.get()+'.pdb')
            except:
                showerror('Error','Please check the internet connection or the validity of the PDB ID')
                exit()
        mol = molecule.load_structure(self.Box1.get()+'.pdb')

        #Other parameter validity
        try:
            cutoff_value = float(self.Box3.get())
            power_of_distance = int(self.Box4.get())
        except:
            showerror('Error','Check the value of the parameters:\n\nCutoff Value (float)\nPower of Distance (float)')

        #For checking if the format is right and chain id is right
        try:
            backbone = [j for i in mol[0][self.Box2.get()].get_backbone() for j in i if j is not None]
        except:
            showinfo('Notification','Chain ID is either not specified or not valid. Running hdANM on all the chains.')
            try:
                backbone = []
                for chain in mol[0].get_chains():
                    backbone.extend( [j for i in mol[0][chain.get_id()].get_backbone() for j in i if j is not None] )

            except:    
                showerror('Error','Provide valid chain ID or check the PDB file format.')
                exit()
        
        calpha=[i for i in mol[0]['A'].get_calpha() if i is not None]

        Model=hdANM( calpha,dr=cutoff_value,power=power_of_distance,hng_file=self.Box5.get() )

        if(self.mass_type.get() == "Molecular Weight"):
            Model.calculate_hessian(mass_type='residue')
        elif(self.mass_type.get() == "Unit"):
            Model.calculate_hessian(mass_type='unit')
        elif(self.mass_type.get() == "Atomics Mass"):
            Model.calculate_hessian(mass_type='atom')
        
        Model.calculate_decomposition()
        
        #Results
        pop_up1 = tk.Toplevel()
        pop_up1.title('hdANM Results')
        tk.Label(pop_up1,text="Save the hdANM output matrices", borderwidth=2, relief="groove").grid(row=0,columnspan=3,sticky=N+S+E+W, padx=10, pady=10 )

        tk.Label(pop_up1,text="Select the matrix ").grid(row=2,column=0,sticky=N+S+E+W, padx=10, pady=10 )

        #Section 1
        user_matrix_choise = StringVar(pop_up1)
        user_matrix_choise.set("hdANM Eigen Vectors")
        tk.OptionMenu(pop_up1, user_matrix_choise, "hdANM Eigen Vectors", "hdANM Eigen Values ").grid(row=2,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Button(pop_up1,text='Save .CSV file',command=lambda user_matrix_choise=user_matrix_choise : save_matrix(user_matrix_choise.get()) ).grid(row=3,column=2,sticky=N+S+E+W, padx=10, pady=10 )


        #Section 2
        tk.Label(pop_up1,text="Generate Movies from Eigenvectors",borderwidth=2, relief="groove").grid(row=5,columnspan=3,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(pop_up1, text="Note: Non rigid modes start from 6").grid(row=6,columnspan=3,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(pop_up1, text="Projection Method").grid(row=7,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(pop_up1, text="Mode Number").grid(row=8,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(pop_up1, text="Number of Frames").grid(row=9,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(pop_up1, text="Movie Scale").grid(row=10,column=0,sticky=N+S+E+W, padx=10, pady=10 )

        user_projection_choise = StringVar(pop_up1)
        user_projection_choise.set("Curvilinear")
        tk.OptionMenu(pop_up1, user_projection_choise, "Curvilinear", "Linear").grid(row=7,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        n_mode   = tk.Entry(pop_up1)
        n_mode.grid(row=8 ,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        n_mode.insert(END,6)
        n_frames = tk.Entry(pop_up1)
        n_frames.grid(row=9 ,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        n_frames.insert(END,10)
        n_scale  = tk.Entry(pop_up1)
        n_scale.grid(row=10,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        n_scale.insert(END,2)

        
        self.save_movie_button = tk.Button(pop_up1,text='Save Movie',command=lambda : save_movie( n_mode, n_frames, n_scale, user_projection_choise.get() ) )
        self.save_movie_button.grid(row=11,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        

        def save_matrix(choise):
            if(choise  == "hdANM Eigen Vectors"):
                numpy.savetxt("hdANM Eigen Vectors.csv",Model.get_eigenvectors(),delimiter=",")
                showinfo('Notification','hdANM Eigen Vectors.csv file saved!')
            elif(choise== "hdANM Eigen Values" ):
                numpy.savetxt("hdANM Eigen Values.csv",Model.get_eigenvalues(),delimiter=",")
                showinfo('Notification','hdANM Eigen Values.csv file saved!')
            return True
        
        def save_movie(n_mode, n_frames, n_scale, projection):
            try:
                n_mode = int(n_mode.get())
                n_frames = int(n_frames.get())
                n_scale = float(n_scale.get())
            except:
                showinfo('Notification','Please check the parameters.\n\nMode Number (int)\nNumber of Frames (int)\nMovie Scale (float)')

            if(projection == "Curvilinear"):
                Model.calculate_movie(n_mode, n=n_frames,scale=n_scale,extrapolation="curvilinear")
            elif(projection =="Linear"):
                Model.calculate_movie(n_mode, n=n_frames,scale=n_scale,extrapolation="linear")
            showinfo('Notification',str(n_mode)+'.pdb movie generated & saved.')
            return True
        
        return True

def load_gui():
    app = Skeleton('HomePage')
    app.mainloop()


if(__name__ == '__main__'):
    GUI_OPTION = False
    try:
        if(sys.argv[1] == "gui"):GUI_OPTION = True
    except:
        main()
    
    if(GUI_OPTION):load_gui()