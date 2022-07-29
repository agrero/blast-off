import tkinter as tk
from tkinter import ttk
from tkinter import *
from tkinter import messagebox
import os

import time

from Bio import SeqIO, Entrez, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datetime import date


#put somewhere that actually makes sense later
def dir_check(parent_directory:str, directory_name:str): 
    if os.path.isdir(os.path.join(parent_directory, directory_name)):
        pass
    else:
        os.makedirs(os.path.join(parent_directory, directory_name))



class Parent(tk.Tk):
    def __init__(self):
        super().__init__()      

        #notebook used as the apps directory

        self.directory_bar  = ttk.Notebook(self)
        self.directory_bar.pack(expand=True)

        #initialization procedures

        self.parent_directory = os.getcwd() #alter this to be taken from settings!
        self.accession_codes = self.get_current_codes()
        dir_check(self.parent_directory, 'genomes')
        self.genome_directory = os.path.join(self.parent_directory, 'genomes')
        self.genome_names = self.getGenomeNames()
        self.records_to_parse = []

        #landing page variables and widgets
        
        self.landing_frame = ttk.Frame(
            self.directory_bar,
            width=400,
            height=280
        )
        self.landing_frame.pack()

        #editing accession code variables and widgets

        self.code_editor_frame = ttk.Frame(
            self.directory_bar,
            width=400,
            height=280
        )
        self.code_editor_frame.pack(fill='both',expand=True)

        self.code_box = Listbox(
            self.code_editor_frame,
            width=250,
            height=10, 
            bd=1,
            highlightthickness=2,
            activestyle = "none"
        )
        self.code_box.pack(side='left', fill='both', expand=True)

        self.code_scrollbar = Scrollbar(self.code_box)
        self.code_scrollbar.pack(side='right', fill='y')

        self.code_box.config(yscrollcommand=self.code_scrollbar.set)
        self.code_scrollbar.config(command=self.code_box.yview)

        self.broken_codes = []

        for item in self.accession_codes:
            self.code_box.insert(END, item)

        self.code_entry_label = Label(self.code_editor_frame, text="Code Entry")
        self.code_entry_label.pack(side="top")

        self.code_entry = Entry(self.code_editor_frame)
        self.code_entry.pack(side='top', fill="x")

        self.add_code_btn = Button(
            self.code_editor_frame,
            text = "Add Code",
            command = self.newCode
        )
        self.add_code_btn.pack(expand=False, side='top')

        self.save_btn = Button(
            self.code_editor_frame,
            text = "Save",
            command = self.save_current_codes
        )
        self.save_btn.pack(expand=False, side='top')

        self.update_genomes_btn = Button(
            self.code_editor_frame,
            text = "Update",
            command = self.updateGenomes
        )
        self.update_genomes_btn.pack(expand=False, side='top')

        #name search variables and widgets

        self.name_search_frame = ttk.Frame(
            self.directory_bar,
            width=400,
            height=280
        )
        self.name_search_frame.pack(fill='both', expand=True)

        self.input_search_label = Label(
            self.name_search_frame,
            text = "Search Term"
        )
        self.input_search_label.pack(fill='both', expand=True, side='top', anchor='s')

        self.protein_entry = Entry(self.name_search_frame)
        self.protein_entry.pack(fill='x', expand=True, side='top', anchor='n')
        
        self.add_protein_botton = Button(
            self.name_search_frame,
            text = "Add Protein",
            command = self.addProtein
        )
        self.add_protein_botton.pack(fill='none',side='top',anchor='n')

        self.protein_names_box = Listbox(
            self.name_search_frame,
            width=50,
            height=8
        )
        self.protein_names_box.pack(fill='both',expand=True)

        self.add_protein_label = Label(
            self.name_search_frame,
            text = "Proteins to Search For"
        )
        self.add_protein_label.pack(expand=False, anchor='n', side='top')

        self.name_search_button = Button(
            self.name_search_frame,
            text = 'search',
            command = self.name_search
        )
        self.name_search_button.pack(expand=True, anchor='n', side='right')

        self.proteins_to_search_box = Listbox(
            self.name_search_frame,
            width=50,
            heigh=8,
            selectmode=MULTIPLE
        )
        self.proteins_to_search_box.pack(fill='both', side='bottom',anchor='center')

        self.protein_names = self.getAllGeneNames()
        self.update(self.protein_names)

        self.protein_names_box.bind('<<ListboxSelect>>', self.fillout)

        self.protein_entry.bind('<KeyRelease>', self.check)

        #sequence search variables and widgets

        self.sequence_search_frame = ttk.Frame(
            self.directory_bar,
            width=400,
            height=280
        )

        self.sequence_search_label = Label(
            self.sequence_search_frame,
            text = "Sequence Search Parameters"
        )
        self.sequence_search_label.grid(row=0, column=2, sticky='w')


        ##make scrolledtext
        self.sequence_input = Text(
            self.sequence_search_frame,
            height=18,
            width=50
        )
        self.sequence_input.grid(row=2, column= 2, sticky="e", rowspan=15)
        self.sequence_input_label = Label(
            self.sequence_search_frame,
            text = "Sequence Input"
        )
        self.sequence_input_label.grid(row=1, column=2, sticky="n")
        
        #can likely make a function to do all of this

        self.score_threshold = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.score_threshold.grid(row=2, column = 1, sticky="W")
        self.score_threshold_label = Label(
            self.sequence_search_frame,
            text = "Score Threshold",
        )
        self.score_threshold_label.grid(row=2, column=0, sticky="E")

        self.match_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.match_score.grid(row=5, column = 1, sticky="W")
        self.match_score_label = Label(
            self.sequence_search_frame,
            text = "Match Score",
        )
        self.match_score_label.grid(row=5, column=0, sticky="E")

        self.mismatch_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.mismatch_score.grid(row=6, column = 1, sticky="W")        
        self.mismatch_score_label = Label(
            self.sequence_search_frame,
            text = "Mismatch Score",
        )
        self.mismatch_score_label.grid(row=6, column=0, sticky="E")

        self.target_internal_open_gap_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.target_internal_open_gap_score.grid(row=7, column = 1, sticky="W")
        self.target_internal_open_gap_score_label = Label(
            self.sequence_search_frame,
            text = "T. Int. Open Gap"
        )
        self.target_internal_open_gap_score_label.grid(row=7, column=0, sticky="E")

        self.target_left_open_gap_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.target_left_open_gap_score.grid(row=8, column = 1, sticky="W")
        self.target_left_open_gap_score_label = Label(
            self.sequence_search_frame,
            text = "T. Left Open Gap"
        )
        self.target_left_open_gap_score_label.grid(row=8, column=0, sticky="E")

        self.target_left_extend_gap_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.target_left_extend_gap_score.grid(row=9, column = 1, sticky="W")
        self.target_left_extend_gap_score_label = Label(
            self.sequence_search_frame,
            text = "T. Left Ext. Gap"
        )
        self.target_left_extend_gap_score_label.grid(row=9, column=0, sticky="E")

        self.target_right_open_gap_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.target_right_open_gap_score.grid(row=10, column=1, sticky="W")
        self.target_right_open_gap_score_label = Label(
            self.sequence_search_frame,
            text = "T. Right Open Gap"
        )
        self.target_right_open_gap_score_label.grid(row=10, column=0, sticky="E")

        self.target_right_extend_gap_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.target_right_extend_gap_score.grid(row=11, column = 1, sticky="W")
        self.target_right_extend_gap_score_label = Label(
            self.sequence_search_frame,
            text = "T. Right Ext. Gap"
        )
        self.target_right_extend_gap_score_label.grid(row=11, column=0, sticky="E")

        self.query_left_extend_gap_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.query_left_extend_gap_score.grid(row=12, column = 1, sticky="W")
        self.query_left_extend_gap_score_label = Label(
            self.sequence_search_frame,
            text = "Q. Left Ext. Gap"
        )
        self.query_left_extend_gap_score_label.grid(row=12, column=0, sticky="E")

        self.query_left_open_gap_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.query_left_open_gap_score.grid(row=13, column = 1, sticky="W")
        self.query_left_open_gap_score_label = Label(
            self.sequence_search_frame,
            text = "Q. Left Open Gap"
        )
        self.query_left_open_gap_score_label.grid(row=13, column=0, sticky="E")

        self.query_left_extend_gap_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.query_left_extend_gap_score.grid(row=14, column = 1, sticky="W")
        self.query_left_extend_gap_score_label = Label(
            self.sequence_search_frame,
            text = "Q. Left Ext. Gap"
        )
        self.query_left_extend_gap_score_label.grid(row=14, column=0, sticky="E")

        self.query_right_open_gap_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.query_right_open_gap_score.grid(row=15, column = 1, sticky="W")
        self.query_right_open_gap_score_label = Label(
            self.sequence_search_frame,
            text = "Q. Right Open Gap"
        )
        self.query_right_open_gap_score_label.grid(row=15, column=0, sticky="E")

        self.query_right_extend_gap_score = Entry(
            self.sequence_search_frame,
            width=7
        )
        self.query_right_extend_gap_score.grid(row=16, column = 1, sticky="W")
        self.query_right_extend_gap_score_label = Label(
            self.sequence_search_frame,
            text = "Q. Right Ext. Gap"
        )
        self.query_right_extend_gap_score_label.grid(row=16, column=0, sticky="E")
        
        self.sequence_search_btn = Button(
            self.sequence_search_frame,
            text="Search",
            command = self.sequence_search
        )
        self.sequence_search_btn.grid(row=1, column=2, sticky="e")

        self.name_sequence_search_entry = Entry(
            self.sequence_search_frame,
            width=80
        )
        self.name_sequence_search_entry.grid(row=17, column=1, columnspan=2)
        self.name_sequence_search_label = Label(
            self.sequence_search_frame,
            text = "Save Name"
        ).grid(row=17, column=0)
        #Settings variables and widgets

        self.settings_frame = ttk.Frame(
            self.directory_bar,
            width=400,
            height=280
        )
        self.settings_frame.pack(expand=True, fill='both')

        self.settings_label = Label(
            self.settings_frame,
            text = 'settings'
        )
        self.settings_label.grid(row=0, column=0)

        self.save_directory_label = Label(
            self.settings_frame,
            text='Save Folder Path'
        )
        self.save_directory_label.grid(column=0, row=1)

        self.save_directory_entry = Entry(
            self.settings_frame,
            width=80,
        )
        self.save_directory_entry.grid(column=1, row=1)

        self.save_settings_button = Button(
            self.settings_frame,
            text= "Save Settings",
            command=None
        )
        self.save_settings_button.grid(row=2, column=1)


        #notebook directory items added here
        self.directory_bar.add(self.landing_frame, text='Landing')
        self.directory_bar.add(self.name_search_frame, text='Name Search')
        self.directory_bar.add(self.sequence_search_frame, text="Sequence Search")
        self.directory_bar.add(self.code_editor_frame, text="Code Entry")
        self.directory_bar.add(self.settings_frame, text="Settings")

    #initialization methods

    def get_current_codes(self):
        """Retrieves most up to date codes from the specified directory"""
        codes_path = os.path.join(self.parent_directory, "genome accession codes") ## can make these first three lines another function
        code_lists = os.listdir(codes_path)##
        codes_to_use = os.path.join(codes_path, code_lists[0])##

        with open(os.path.join(codes_path, codes_to_use)) as file: #this might be able to be done as a lambda function
                    text = file.read()
                    text_2 = text.strip("current genomes:")
                    text_3 = text_2.strip("()")
                    text_4 = text_3.replace("'","")
                    file.close()
                    return text_4.split(",")

    def save_setttings(self):
        setting_names = ("Save Path")
        settingls = []
        settings = {}
        for setting in settingls:
            pass

    #editing accession code methods

    def save_current_codes(self):
        from datetime import date
        current_codes = self.code_box.get(0,END)
        file_name = "current_genomes_{}.txt".format(date.today())
        writing_directory = "genome accession codes"
        complete_name = os.path.join(self.parent_directory, writing_directory, file_name)

        file = open(complete_name, "w")
        file.write("current genomes:" + str(current_codes))
        file.close()

    def newCode(self):
        code = self.code_entry.get()
        if code != "":
            self.code_box.insert(END, code)
            self.code_entry.delete(0, "end")
        else:
            messagebox.showwarning("Warning", "Please enter some code")

    def deleteCode(self):
        self.code_box.delete(ANCHOR)

    def updateGenomes(self):
        Entrez.email = "guerale2@isu.edu" #make this entry on landing page
        dir_check(self.parent_directory, "genomes")
        genome_dir = os.path.join(self.parent_directory, "genomes")

        for record in self.accession_codes:
            try:
                handle = Entrez.efetch(db="nucleotide", id=record, rettype="gb", retmode="text")
            except:
                time.sleep(2)
                pass
            try:
                handle = Entrez.efetch(db="nucleotide", id=record, rettype="gb", retmode="text")
            except:
                self.broken_codes.append(record)
         
            handle = Entrez.efetch(db="nucleotide", id=record, rettype="gb", retmode="text")
            species_read = SeqIO.read(handle, "genbank")
            species_name = species_read.description
            filename = species_name + ".gbk"
            filepath = os.path.join(genome_dir, filename)
            SeqIO.write(species_read, filepath, "genbank")
            handle.close()

    #name search methods

    def getAllFromListbox(self, listbox:Listbox) -> list:
        listbox_contents = []
        listbox.select_set(0, END)
        selected = listbox.get(0, END)
        for i in selected:
            listbox_contents.append(i)
        return listbox_contents

    def addProtein(self):
        protein = self.protein_entry.get()
        if protein != "":
            self.proteins_to_search_box.insert(END, protein)
            self.protein_entry.delete(0, "end")
        else:
            messagebox.showwarning("Warning", "Please enter a protein name!")
    
    def getGenomeNames(self): 
        species = []
        for genome in os.listdir(self.genome_directory):
            genome_name = genome.strip('.gbk')
            species.append(genome_name)
        return species

    def isCDS(self, record_feature) -> bool:
        if record_feature.type == "CDS":
            return True
        else:
            return False
    
    def getAllGeneNames(self):
        all_gene_names = []
        for genome in os.listdir(self.genome_directory):
            genome_path = os.path.join(self.genome_directory, genome)
            record = SeqIO.read(genome_path, "genbank")
            for feature in record.features:
                if self.isCDS(feature):
                    try:
                        feature.qualifiers["product"]
                    except:
                        continue
                    all_gene_names.append(feature.qualifiers["product"][0])
        return set(all_gene_names)

    def update(self, data):
        #Clear the listbox
        self.protein_names_box.delete(0,END)

        #Add data to listbox
        for item in data:
            self.protein_names_box.insert(END, item)
    
    def fillout(self, event):
        #delete whatever is in the entry box
        self.protein_entry.delete(0, END)

        #add clicked list item to entry box
        self.protein_entry.insert(0, self.protein_names_box.get(ANCHOR))
    
    def check(self, event):
        #grab what was typed
        typed = self.protein_entry.get()

        if typed == '':
            data = self.protein_names
        else:
            data = []
            for item in self.protein_names:
                if typed.lower() in item.lower():
                    data.append(item)
        
        self.update(data)

    def score_check(self, score_input): #change to have it pop up a messagebox
        try:
            int(score_input)
        except:
            return 0
        return int(score_input)
    
    def seq_input_check(self, sequence):
        if type(sequence) == str:
            return sequence
        else:
            messagebox.showwarning("Warning", "Invalid Sequence Entered")

    def name_search(self):
        feat_count = 1
        rec_list = self.getAllFromListbox(self.proteins_to_search_box)
        prot_to_update = {}
        for genome in os.listdir(self.genome_directory):
            path = os.path.join(self.genome_directory, genome)
            record = SeqIO.read(path, "genbank")
            for feature in record.features:
                try:
                    feature.qualifiers["product"][0]
                    feature.qualifiers["translation"][0]
                except:
                    continue
                if self.isCDS(feature) and feature.qualifiers["product"][0] in rec_list:
                    name = feat_count
                    rec = SeqRecord(
                        Seq(feature.qualifiers["translation"][0],),
                        id = feature.qualifiers["product"][0],
                        description = genome
                    )
                    prot_to_update[name] = rec
                feat_count += 1
        dir_check(self.parent_directory, "saved searches")
        filename = "name_search-{}".format(date.today())
        filepath = os.path.join(self.parent_directory, "saved searches", filename)
        with open(filepath, "w") as handle:
            SeqIO.write(prot_to_update.values(), handle, "fasta")


    
    #sequence search methods

    def sequence_search(self):
        count = 1
        prot_to_update = {} 
        filename = self.name_sequence_search_entry.get()
        if filename.strip(" ") == "":
            messagebox.showwarning("Warning", "Please enter a savename")
            return None
        query_sequence = self.seq_input_check(self.sequence_input.get("1.0", END))
        score_threshold = self.score_check(self.score_threshold.get())
        match_score = self.score_check(self.match_score.get())
        mismatch_score = self.score_check(self.mismatch_score.get())
        target_internal_open_gap_score = self.score_check(self.target_internal_open_gap_score.get())
        target_left_open_gap_score = self.score_check(self.target_left_open_gap_score.get())
        target_left_extend_gap_score = self.score_check(self.target_left_extend_gap_score.get())
        target_right_open_gap_score = self.score_check(self.target_right_open_gap_score.get())
        target_right_extend_gap_score = self.score_check(self.target_right_extend_gap_score.get())
        query_left_open_gap_score = self.score_check(self.query_left_open_gap_score.get())
        query_left_extend_gap_score = self.score_check(self.query_left_extend_gap_score.get())
        query_right_open_gap_score = self.score_check(self.query_right_open_gap_score.get())
        query_right_extend_gap_score = self.score_check(self.query_right_extend_gap_score.get())
        wildcard = "?"

        aligner = Align.PairwiseAligner()
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        aligner.target_internal_open_gap_score = target_internal_open_gap_score
        aligner.target_left_open_gap_score = target_left_open_gap_score
        aligner.target_left_extend_gap_score = target_left_extend_gap_score
        aligner.query_right_open_gap_score = target_right_open_gap_score
        aligner.target_right_extend_gap_score = target_right_extend_gap_score
        aligner.query_left_open_gap_score = target_left_open_gap_score
        aligner.query_left_open_gap_score = query_left_open_gap_score
        aligner.query_left_extend_gap_score = query_left_extend_gap_score
        aligner.query_right_open_gap_score = query_right_open_gap_score
        aligner.query_right_extend_gap_score = query_right_extend_gap_score
        aligner.wildcard = wildcard

        for genome in os.listdir(self.genome_directory):
            path = os.path.join(self.genome_directory, genome)
            record = SeqIO.read(path, "genbank")
            for feature in record.features:
                try:
                    feature.qualifiers["product"][0]
                    feature.qualifiers["translation"][0]
                except:
                    continue
                if self.isCDS(feature) and aligner.score(query_sequence, feature.qualifiers["translation"][0]) >= score_threshold:
                    name = count
                    rec = SeqRecord(
                        Seq(feature.qualifiers["translation"][0],),
                        id = feature.qualifiers["product"][0],
                        description = genome
                    )
                    prot_to_update[name] = rec
                count += 1
        dir_check(self.parent_directory, "saved searches")
        filepath = os.path.join(self.parent_directory, "saved searches", filename)
        with open(filepath, "w") as handle:
            SeqIO.write(prot_to_update.values(), handle, "fasta")


if __name__ == "__main__":
    app = Parent()
    app.mainloop()


#things to add for a settigns
#place to pick your save directory
#place to pick your genome load directory
#button to make them the same
#system to write this to a file
#make the intitialization run the reading of the file
#----------------------------------------------------
#start application
#application checks for save file
#if none exists creates default directories
#store the paths for the settings page
#if one exists read in the settings
#-----------------------------------------------------
#settings 
#search save path
#---------
