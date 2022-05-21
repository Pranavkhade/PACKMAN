#!/usr/bin/env python
"""The PACKMAN Graphical User Interface (GUI) host file.

How to use::
    python -m packman gui #(For GUI)
    python -m packman     #(For CLI)

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""

'''
##################################################################################################
#                                             GUI                                                #
##################################################################################################
'''
from xmlrpc.client import boolean
import numpy
import logging

from .. import molecule
from ..apps import predict_hinge
from ..anm import hdANM
from ..entropy import PackingEntropy

try:
    import tkinter as tk
    from tkinter.messagebox import showinfo, showerror
    from tkinter import N, S, E, W, Grid, END, BooleanVar, StringVar, DISABLED
    from tkinter.filedialog import askopenfilename
    import codecs
    import os.path
except:
    print("The Python version doesn't the (built in) tkinter package. Please make sure you are using right Python interpreter. ( try: python3 -m packman gui )")


#read and get_version is to print the PACKMAN version on the GUI
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    '''To find out the version of the 
    '''
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

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
        self.frames['Voronoi_Packing_Entropy_GUI'] = Voronoi_Packing_Entropy_GUI(self)
        self.frames['TopMenu'] = top_menu(self)

        self.geometry("800x600")

        #Icon for windows
        try:
            self.iconbitmap(False, os.path.abspath(os.path.dirname(__file__)+'/logo.ico') )
        except:
            None
        
        self.show_frame(frame_name)

        try:
            self.title_content = "PACKMAN "+get_version( os.path.abspath(os.path.dirname(__file__)).replace('bin','__init__.py') )+" GUI"
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
        number_of_buttons = 4
        Canvas_top_menu = tk.Canvas(parent)
        Canvas_top_menu.grid(row=0,column=0,columnspan=4, sticky=E+W )

        for i in range(0,number_of_buttons):
            Grid.columnconfigure(Canvas_top_menu, i, weight=1)
        
        #tk.Frame.__init__(self)
        tk.Button(Canvas_top_menu, text = 'Home',             command = lambda: parent.show_frame('HomePage')        ).grid( row=0,column=0,sticky=E+W )
        tk.Button(Canvas_top_menu, text = 'Hinge Prediction', command = lambda: parent.show_frame('HingePrediction') ).grid( row=0,column=1,sticky=E+W )
        tk.Button(Canvas_top_menu, text = 'hdANM',            command = lambda: parent.show_frame('hdANM_GUI')       ).grid( row=0,column=2,sticky=E+W )
        tk.Button(Canvas_top_menu, text = 'Packing Entropy (Voronoi)',          command = lambda: parent.show_frame('Voronoi_Packing_Entropy_GUI')     ).grid( row=0,column=3,sticky=E+W )

def load_gui():
    app = Skeleton('HomePage')
    app.mainloop()

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
        self.Text1.insert(END,"1. Tutorials: https://py-packman.readthedocs.io/en/latest/index.html\n2. GitHub: https://github.com/Pranavkhade/PACKMAN\n3. Online Hinge Prediction Server: https://packman.bb.iastate.edu/\n4. Online hdANM Server: https://hdanm.bb.iastate.edu/\n5. Online Voronoi Entropy Server: https://packing-entropy.bb.iastate.edu/")
        self.Text1.config(state=DISABLED)

        self.all_objects = self.__dict__

    def show(self):
        self.Label1.grid(row=1,columnspan=3, sticky=E+W, padx=10, pady=10 )
        self.Label2.grid(row=2,columnspan=3, sticky=E+W, padx=10, pady=10 )
        self.Text1.grid(row=3,columnspan=3, sticky=E+W, padx=10, pady=10 )
    
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
        self.Label3 = tk.Label(parent, text="Filename/PDB ID:")
        self.Label4 = tk.Label(parent, text="Chain ID:")
        self.Label5 = tk.Label(parent, text="Alpha Value:")
        self.Label6 = tk.Label(parent, text="Number of Eccentricity Clusters:")
        self.Label7 = tk.Label(parent, text="Minimum Hinge Length:")

        #Input
        self.Text1.insert(END,"Interface to the functionality to identify the protein hinges (separating the domains) using PACKMAN. It can also be used to read, write, manipulate and analyze protein molecules and it's properties through its API.\n\nCitation:\nKhade PM, Kumar A, Jernigan RL. Characterizing and Predicting Protein Hinges for Mechanistic\nInsight. J Mol Biol. November 2019. doi:10.1016/j.jmb.2019.11.018\n\nNOTE: Please specify filename with appropriate extension (.pdb/.cif); if the file specified is not present, it will be downloaded automatically.(if the PDB ID is valid)" )
        self.Text1.config(state=DISABLED)
        self.Box1 = tk.Entry(parent)
        self.Box2 = tk.Entry(parent)
        self.Box3 = tk.Entry(parent)
        self.Box4 = tk.Entry(parent)
        self.Box5 = tk.Entry(parent)

        self.Box3.insert(END, '2.8')
        self.Box4.insert(END, '4')
        self.Box5.insert(END, '5')

        self.Browse1 = tk.Button(parent, text = "Browse", command = lambda: self.browseFiles() )
        self.Button1 = tk.Button(parent, text = "Run", command = lambda: self.run_PACKMAN() )

        self.all_objects = self.__dict__

    def show(self):
        self.Text1.grid(row=2,columnspan=4,sticky=E+W, padx=10, pady=10 )
        self.Label1.grid(row=1,columnspan=4,sticky=E+W, padx=10, pady=10 )
        self.Label2.grid(row=3,columnspan=4,sticky=E+W, padx=10, pady=10 )
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
        
        self.Browse1.grid(row=4,column=3,sticky=N+S+E+W, padx=10, pady=10 )
        self.Button1.grid(row=9, column=2)
    
    def hide(self):
        for i in self.all_objects:
            try:
                self.all_objects[i].grid_forget()
            except:
                continue
    
    def browseFiles(self):
        filename = askopenfilename( title = "Select a File", filetypes = [ ('mmCIF','*.cif'),('PDB','*.pdb'),("All files", "*") ] )
        self.Box1.delete(0,END)
        self.Box1.insert(END,filename)
    
    def run_PACKMAN(self):
        """
        PACKMAN Hinge Prediction Algorithm
        """
        try:
            mol = molecule.load_structure(self.Box1.get())
        except:
            showinfo('Notification','File not found in the folder, downloading it from the PDB.')
            
            #For file availability
            try:
                try:
                    molecule.download_structure(self.Box1.get().split('.')[0], ftype= self.Box1.get().split('.')[1])
                except:
                    molecule.download_structure( self.Box1.get() )
            except:
                showerror('Error','Please check the internet connection or the validity of the PDB ID')
                exit()        
        try:
            mol = molecule.load_structure( self.Box1.get() )
        except:
            mol = molecule.load_structure(self.Box1.get()+'.cif')

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
                    try:
                        backbone = [j for i in mol[0][chain.get_id()].get_backbone() for j in i if j is not None]
                        predict_hinge(backbone, Alpha=alpha, nclusters=ecc_clust, MinimumHingeLength=min_length, outputfile=open('output.txt', 'w'))
                    except:
                        None

            except:    
                showerror('Error','Possible Solutions:\n1. Provide valid chain ID or check the PDB file format.\n2. Try changing Alpha Value')
                exit()
        
        #At this point, everything is okay
        
        
        #Hinge Prediction Results
        pop_up1 = tk.Toplevel()
        pop_up1.title('PACKMAN Results')

        total_columns_in_canvas1 = 8
        Canvas1 = tk.Canvas(pop_up1)
        Canvas1.grid(row=0,column=0,columnspan=total_columns_in_canvas1,sticky=N+S+E+W, padx=10, pady=10 )

        #Scrollbar
        vsb = tk.Scrollbar(pop_up1, orient="vertical", command=Canvas1.yview)
        vsb.grid(row=0, column=total_columns_in_canvas1, rowspan=10 , sticky=N+S)
        Canvas1.configure(yscrollcommand=vsb.set)

        #Frame
        canvas1_frame1 = tk.Frame(Canvas1)
        Canvas1.create_window((0, 0), window=canvas1_frame1, anchor='nw')

        tk.Label(canvas1_frame1, text = 'NOTE: Please check in front of the hinges you want to keep and click on the "Save HNG file" button to generate input file for hdANM. Please provide a valid filename\n Tip: p-value ideally should be <0.05 ').grid(row=0,columnspan=total_columns_in_canvas1)
        tk.Label(canvas1_frame1, text = 'Hinge ID').grid(row=1,column=0,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(canvas1_frame1, text = 'PDB ID/Filename').grid(row=1,column=1,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(canvas1_frame1, text = 'Chain ID').grid(row=1,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(canvas1_frame1, text = 'Residue IDs').grid(row=1,column=3,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(canvas1_frame1, text = 'p-value').grid(row=1,column=4,sticky=N+S+E+W, padx=10, pady=10 )
        tk.Label(canvas1_frame1, text = 'Select').grid(row=1,column=5,sticky=N+S+E+W, padx=10, pady=10 )

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

            tk.Label(canvas1_frame1, text = numi+1).grid(row=numi+2,column=0,sticky=N+S+E+W, padx=10, pady=10 )
            tk.Label(canvas1_frame1, text = self.Box1.get() ).grid(row=numi+2,column=1,sticky=N+S+E+W, padx=10, pady=10 )
            tk.Label(canvas1_frame1, text = ChainOfHinge    ).grid(row=numi+2,column=2,sticky=N+S+E+W, padx=10, pady=10 )
            tk.Label(canvas1_frame1, text =  str(Residue_IDs[0])+':'+str(Residue_IDs[-1]) ).grid(row=numi+2,column=3,sticky=N+S+E+W, padx=10, pady=10 )
            tk.Label(canvas1_frame1, text = i.get_pvalue()).grid(row=numi+2,column=4,sticky=N+S+E+W, padx=10, pady=10 )

            Selected_Hinges[numi] = BooleanVar(canvas1_frame1)
            Check = tk.Checkbutton(canvas1_frame1, variable=Selected_Hinges[numi], onvalue=True, offvalue=False)
            Check.grid(row=numi+2,column=5,sticky=N+S+E+W, padx=10, pady=10 )

            #Details Button
            tk.Button(canvas1_frame1,text = 'Details', command= lambda i=i:show_details( i.get_stats() ) ).grid(row=numi+2,column=6,sticky=N+S+E+W, padx=10, pady=10 )

        #Fix the length of the result
        canvas1_frame1.update_idletasks()
        Canvas1.config(width=1000 ,height = 400)
        Canvas1.config(scrollregion = Canvas1.bbox("all"))
        
        tk.Label(pop_up1, text = 'output filename (.hng)').grid(row=2,column=1,sticky=N+S+E+W, padx=10, pady=10 )
        self.output_filename = tk.Entry(pop_up1)
        tk.Button(pop_up1,text = 'Save HNG file', command=lambda: write_hng_file() ).grid(row=2,column=3,sticky=N+S+E+W, padx=10, pady=10 )

        self.output_filename.insert(END,self.Box1.get()+'.hng')
        self.output_filename.grid(row=2,column=2,sticky=N+S+E+W, padx=10, pady=10 )

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
                        fh.write(self.Box1.get().split('/')[-1].replace(' ','+')+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(ALL_RESIDUES[ChainOfHinge][0])+':'+str(hinge_res_ids[0]-1)+'\n' )
                        fh.write(self.Box1.get().split('/')[-1].replace(' ','+')+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
                    else:
                        fh.write(self.Box1.get().split('/')[-1].replace(' ','+')+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
                    last_hinge_end = hinge_res_ids[-1]
                elif( i.get() ):
                    hinge_res_ids = sorted([j.get_id() for j in current_hinge.get_elements()])
                    select_count += 1
                    fh.write(self.Box1.get().split('/')[-1].replace(' ','+')+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(last_hinge_end+1)+':'+str(hinge_res_ids[0]-1)+'\n' )
                    fh.write(self.Box1.get().split('/')[-1].replace(' ','+')+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
                    last_hinge_end = hinge_res_ids[-1]
                try:
                    if(ChainOfHinge != All_Hinges[numi+1].get_elements()[0].get_parent().get_id() ):
                        select_count += 1
                        fh.write(self.Box1.get().split('/')[-1].replace(' ','+')+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(last_hinge_end+1)+':'+str(ALL_RESIDUES[ChainOfHinge][-1])+'\n' )
                        last_hinge_end = 0
                except:
                    None

            if(last_hinge_end != ALL_RESIDUES[ChainOfHinge][-1]):
                select_count += 1
                fh.write(self.Box1.get().split('/')[-1].replace(' ','+')+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(last_hinge_end+1)+':'+str(ALL_RESIDUES[ChainOfHinge][-1])+'\n' )
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

        self.Label3 = tk.Label(parent, text="Filename/PDB ID:")
        self.Label4 = tk.Label(parent, text="Chain ID:")
        
        self.Label5 = tk.Label(parent, text="Cutoff Value")
        self.Label6 = tk.Label(parent, text="Power of Distance")
        self.Label7 = tk.Label(parent, text="Mass of the residue")
        self.Label7 = tk.Label(parent, text="Input (.hng) Filename")
        self.Label8 = tk.Label(parent, text="Mass of the residue")

        #Input
        self.Text1.insert(END,"hdANM is a comprehensive Elastic Network Model. Please read the paper for more details.\n\nCitation:\n\nPranav M. Khade, Domenico Scaramozzino, Ambuj Kumar, Giuseppe Lacidogna, Alberto Carpinteri, Robert L. Jernigan, \nhdANM: a new comprehensive dynamics model for protein hinges,\nBiophysical Journal, 2021, https://doi.org/10.1016/j.bpj.2021.10.017.\n\nNOTE: Please specify filename with appropriate extension (.pdb/.cif); if the file specified is not present, it will be downloaded automatically.(if the PDB ID is valid)" )
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

        self.Browse1 = tk.Button(parent, text = "Browse", command = lambda: self.browseFiles(1, filetypes = [ ('mmCIF','*.cif'),('PDB','*.pdb'),("All files", "*") ] ) )
        self.Browse2 = tk.Button(parent, text = "Browse", command = lambda: self.browseFiles(5, filetypes = [ ('HNG','*.hng'),("All files", "*") ]) )
        self.Button1 = tk.Button(parent, text = "Run", command = lambda: self.run_hdANM() )

        self.all_objects = self.__dict__
    
    def show(self):
        self.Text1.grid(row=2,columnspan=4,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label1.grid(row=1,columnspan=4,sticky=N+S+E+W, padx=10, pady=10 )
        self.Label2.grid(row=3,columnspan=4,sticky=N+S+E+W, padx=10, pady=10 )
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
        
        self.Browse1.grid(row=4,column=3,sticky=N+S+E+W, padx=10, pady=10 )
        self.Browse2.grid(row=8,column=3,sticky=N+S+E+W, padx=10, pady=10 )
        self.Button1.grid(row=10, column=2)
    
    
    def hide(self):
        for i in self.all_objects:
            try:
                self.all_objects[i].grid_forget()
            except:
                continue
    
    def browseFiles(self,fill_box_num,filetypes):
        filename = askopenfilename( title = "Select a File", filetypes=filetypes )
        if(fill_box_num==1):
            self.Box1.delete(0,END)
            self.Box1.insert(END,filename)
        elif(fill_box_num==5):
            self.Box5.delete(0,END)
            self.Box5.insert(END,filename)


    def run_hdANM(self):
        try:
            mol = molecule.load_structure(self.Box1.get())
        except:
            showinfo('Notification','File not found in the folder, downloading it from the PDB.')
            
            #For file availability
            try:
                try:
                    molecule.download_structure(self.Box1.get().split('.')[0], ftype= self.Box1.get().split('.')[1])
                except:
                    molecule.download_structure( self.Box1.get() )
            except:
                showerror('Error','Please check the internet connection or the validity of the PDB ID')
                exit()
        try:
            mol = molecule.load_structure( self.Box1.get() )
        except:
            mol = molecule.load_structure(self.Box1.get()+'.cif')

        #Other parameter validity
        try:
            cutoff_value = float(self.Box3.get())
            power_of_distance = int(self.Box4.get())
        except:
            showerror('Error','Check the value of the parameters:\n\nCutoff Value (float)\nPower of Distance (float)')

        #For checking if the format is right and chain id is right
        try:
            calpha=[i for i in mol[0][self.Box2.get()].get_calpha() if i is not None]
        except:
            showinfo('Notification','Chain ID is either not specified or not valid. Running hdANM on all the chains.')
            try:
                calpha = []
                for chain in mol[0].get_chains():
                    try:
                        calpha.extend( [i for i in mol[0][chain.get_id()].get_calpha() if i is not None] )
                    except:
                        None

            except:    
                showerror('Error','Provide valid chain ID or check the PDB file format.')
                exit()
        

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
        tk.Label(pop_up1, text="Project C-Alpha motions to all atoms").grid(row=11,column=0,sticky=N+S+E+W, padx=10, pady=10 )

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
        n_scale.insert(END,10)

        #C-Alpha to all atom projection
        ca_to_aa = BooleanVar()
        n_ca_to_aa = tk.Checkbutton(pop_up1, variable=ca_to_aa, onvalue=True, offvalue= False)
        n_ca_to_aa.grid(row=11,column=2,sticky=N+S+E+W, padx=10, pady=10 )
      
        self.save_movie_button = tk.Button(pop_up1,text='Save Movie',command=lambda : save_movie( n_mode, n_frames, n_scale, user_projection_choise.get(), ca_to_aa ) )
        self.save_movie_button.grid(row=12,column=2,sticky=N+S+E+W, padx=10, pady=10 )
        
        def save_matrix(choise):
            if(choise  == "hdANM Eigen Vectors"):
                numpy.savetxt("hdANM Eigen Vectors.csv",Model.get_eigenvectors(),delimiter=",")
                showinfo('Notification','hdANM Eigen Vectors.csv file saved!')
            elif(choise== "hdANM Eigen Values" ):
                numpy.savetxt("hdANM Eigen Values.csv",Model.get_eigenvalues(),delimiter=",")
                showinfo('Notification','hdANM Eigen Values.csv file saved!')
            return True
        
        def save_movie(n_mode, n_frames, n_scale, projection, ca_to_aa):
            try:
                n_mode = int(n_mode.get())
                n_frames = int(n_frames.get())
                n_scale = float(n_scale.get())
            except:
                showinfo('Notification','Please check the parameters.\n\nMode Number (int)\nNumber of Frames (int)\nMovie Scale (float)')

            if(projection == "Curvilinear"):
                Model.calculate_movie(n_mode, n=n_frames,scale=n_scale,extrapolation="curvilinear",ca_to_aa=ca_to_aa.get())
            elif(projection =="Linear"):
                Model.calculate_movie(n_mode, n=n_frames,scale=n_scale,extrapolation="linear",ca_to_aa=ca_to_aa.get())
            showinfo('Notification',str(n_mode)+'.pdb /.cif movie generated & saved.')
            return True
        
        return True

'''
##################################################################################################
#                                           Entropy                                              #
##################################################################################################
'''

class Voronoi_Packing_Entropy_GUI(tk.Frame):
    def __init__(self,parent):
        tk.Frame.__init__(self)
        self.grid(padx=50)
        self.Text1 =  tk.Text(parent, height=10)
        self.Label1 = tk.Label(parent, text="Entropy")
        self.Label2 = tk.Label(parent, text="Note: Important information in the section above. (Use scrolling)", borderwidth=2, relief="groove")
        self.Label3 = tk.Label(parent, text="Filename/PDB ID:")
        self.Label4 = tk.Label(parent, text="Chain IDs (comma separated):")
        self.Label5 = tk.Label(parent, text="Probe size:")
        self.Label6 = tk.Label(parent, text="On Sphere Points:")
        self.Label7 = tk.Label(parent, text="Output Filename:")

        #Input
        self.Text1.insert(END,"Interface to the functionality to calculate Packing Entropy using PACKMAN. The same program is available via API and CLI.\n\nCitation:\nPaper Under Review\n\nNOTE: Please specify filename with appropriate extension (.pdb/.cif); if the file specified is not present, it will be downloaded automatically.(if the PDB ID is valid).\n\nPlease note that selecting Chain IDs option is used to calculate the entropy with or without the other chains. For example, if a protein has three chains A, B & C, and the user wants to calculate the entropy of chains A & B in the absence of chain C, the user can use Chain IDs A,B parameter-option do so. However, the presence of the chain C is not an issue, but the user wants to calculate the entropy of chain A & B (with C present), Chain IDs option can be ignored, and chain column in the output should be used to select only A & B chains. The use of API is recommended to control these types of situations efficiently that require more control. Also, we recommend reading the publication for more details on the other parameters." )
        self.Text1.config(state=DISABLED)
        self.Box1 = tk.Entry(parent)
        self.Box2 = tk.Entry(parent)
        self.Box3 = tk.Entry(parent)
        self.Box4 = tk.Entry(parent)
        self.Box5 = tk.Entry(parent)

        self.Box3.insert(END, '1.4')
        self.Box4.insert(END, '30')
        self.Box5.insert(END, 'voronoi_packing_entropy_output.txt')

        self.Browse1 = tk.Button(parent, text = "Browse", command = lambda: self.browseFiles() )
        self.Button1 = tk.Button(parent, text = "Run", command = lambda: self.run_Entropy() )

        self.all_objects = self.__dict__

    def show(self):
        self.Text1.grid(row=2,columnspan=4,sticky=E+W, padx=10, pady=10 )
        self.Label1.grid(row=1,columnspan=4,sticky=E+W, padx=10, pady=10 )
        self.Label2.grid(row=3,columnspan=4,sticky=E+W, padx=10, pady=10 )
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
        
        self.Browse1.grid(row=4,column=3,sticky=N+S+E+W, padx=10, pady=10 )
        self.Button1.grid(row=9, column=2)
    
    def hide(self):
        for i in self.all_objects:
            try:
                self.all_objects[i].grid_forget()
            except:
                continue
    
    def browseFiles(self):
        filename = askopenfilename( title = "Select a File", filetypes = [ ('mmCIF','*.cif'),('PDB','*.pdb'),("All files", "*") ] )
        self.Box1.delete(0,END)
        self.Box1.insert(END,filename)

    
    def run_Entropy(self):
        """
        PACKMAN Hinge Prediction Algorithm
        """
        try:
            mol = molecule.load_structure(self.Box1.get())
        except:
            showinfo('Notification','File not found in the folder, downloading it from the PDB.')
            
            #For file availability
            try:
                try:
                    molecule.download_structure(self.Box1.get().split('.')[0], ftype= self.Box1.get().split('.')[1])
                except:
                    molecule.download_structure( self.Box1.get() )
            except:
                showerror('Error','Please check the internet connection or the validity of the PDB ID')
                exit()
        try:
            mol = molecule.load_structure( self.Box1.get() )
        except:
            mol = molecule.load_structure(self.Box1.get()+'.cif')

        #Other parameter validity
        try:
            input_probe_size = float(self.Box3.get())
            input_onspherepoints = int(self.Box4.get())
            output_file = open(self.Box5.get(),'w')
        except:
            showerror('Error','Check the value of the parameters:\n\nAlpha Value(float)\nNumber of Eccentricity Clusters (int)\nMinimum Hinge Length (int)')
        
        try:
            input_chains = self.Box2.get().split(',')
        except:
            input_chains = None
        
        if(input_chains==['']):
            input_chains=None
    
        output_file.write('Frame\tChain\tResidueID\tResidueName\tPackingEntropy\n')
        for i in mol:
            try:
                result = PackingEntropy(i.get_atoms(), chains=input_chains, probe_size=input_probe_size, onspherepoints=input_onspherepoints)
            except:
                logging.error('Please check the parameters.')
                exit()
            
            for j in i.get_residues():
                try:
                    output_file.write( str(i.get_id())+'\t'+str(j.get_parent().get_id())+'\t'+str(j.get_id())+'\t'+str(j.get_name())+'\t'+str(j.get_entropy('PackingEntropy'))+'\n' )
                except:
                    None
            try:
                output_file.write('Total (Frame '+str(i.get_id())+')\t'+str(result.get_total_entropy())+'\n' )
                output_file.write('Total Normalized Entropy (3N-6) (Frame '+str(i.get_id())+')\t'+str( result.get_total_entropy() / ( 3*len([j for j in i.get_atoms()]) -6 ) ) +'\n' )
            except:
                None
        output_file.flush()
        output_file.close()
        showinfo('Notification','The Voronoi Packing Entropy calculation is completed successfully. Please check the output file.')
        #on_click() Ends here
        return True
