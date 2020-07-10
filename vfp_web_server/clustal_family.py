# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 16:24:53 2020

@author: pedro
"""

import pandas as pd
import os
from subprocess import run, PIPE
from Bio.Align.Applications import ClustalOmegaCommandline
import sys




def clustal(in_file = "unaligned.fasta", out_file = "out_filename.fasta"):
    # file = open("unaligned.fasta", "w")
    # file.write(seqs)
    # file.close()
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=False)
    os.system('cmd /c crmapp\clustal-omega-1.2.2-win64\\' + 
              str(clustalomega_cline) + ' --outfmt clustal --force')


def vfp_seqs(filename="Seqs_prots.csv"):
    dataset = pd.read_csv(filename)
    
    dataset['Family'] = dataset['Family'].str.strip()
    
    dataset_3_in_family =  pd.DataFrame({'Family' : [], 
                                           'Sequence_fusogenic': [],
                                           'Sequence': []})
    
    """
    path = "/tmp/fastas/"

    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
        
    """
    
    families_multiple_align = []
    
    for i in dataset.Family.unique():
        dataset_family = dataset.loc[dataset['Family'] == i]
        if dataset_family.shape[0] > 3:
            dataset_3_in_family.append(dataset_family, ignore_index=True)
            file = open(i + ".fasta", "w")
            index = 1
            for j in dataset_family.values:
                file.write(">" + i + " " + str(index) + "\n")
                index += 1
                file.write(j[1] + "\n")
            families_multiple_align.append(i)
             
            file.close()
            
    for i in families_multiple_align:
        clustalomega_cline = ClustalOmegaCommandline(infile=i + ".fasta", 
                                                         outfile=i + ".clustal_num", 
                                                         verbose=True, 
                                                         auto=False)        
        
        os.system('cmd /c crmapp\clustal-omega-1.2.2-win64\\' 
                  + str(clustalomega_cline) + ' --force')
         
        

vfp_seqs()
    
    
    