# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:34:47 2020

@author: pedro
"""


import pandas as pd
import os
from subprocess import run, PIPE
from Bio.Align.Applications import ClustalOmegaCommandline
import sys
from weblogo import *

def vfp_seqs_all(filename="Seqs_prots_fp.csv"):
    dataset = pd.read_csv(filename)
    
    dataset['Family'] = dataset['Family'].str.strip()
    
    file = open("all_proteins.fasta", "w")
        
    index_family = {}
    for i in dataset.Family.unique():
        index_family[i] = 1
    
    for index, row in dataset.iterrows():
        file.write(">" + row['Family'] + " " + str(index_family[row['Family']]) + "\n")
        index_family[row['Family']] += 1
        file.write(row['Sequence_fusogenic'] + "\n")
        
    clustalomega_cline = ClustalOmegaCommandline(infile="all_proteins.fasta", 
                                                         outfile="all_proteins.clustal_num", 
                                                         verbose=True, 
                                                         auto=False)        
        
    os.system('cmd /c crmapp\clustal-omega-1.2.2-win64\\' 
                  + str(clustalomega_cline) + ' --force')
    
    file_out = open("all_proteins.clustal_num", "r")
    seqs = read_seq_data(file_out)
    logodata = LogoData.from_seqs(seqs)
    logooptions = LogoOptions()
    logooptions.title = "VFP WEBSERVER"
    logoformat = LogoFormat(logodata, logooptions)
    eps = txt_formatter(logodata, logoformat)

    weblogo_file = "all_family_weblogo.txt"
    weblogo = open(weblogo_file, "w")
    data_weblogo = str(eps)[2:len(str(eps))-1].replace('\\n', '\n').replace('\\t', '\t')
    weblogo.write(data_weblogo)
    weblogo.close()
    
    
vfp_seqs_all()
        