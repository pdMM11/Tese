# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 16:50:48 2020

@author: pedro
"""

from Bio import AlignIO
from Bio import SeqIO
"""
with open("clustalw.aln", "r") as handle:
    for record in SeqIO.parse(handle, "clustal") :
        print(record)

filename = "out_filename.clustal_num"
format = "clustal"
alignment = AlignIO.read(filename, format)
# print("Alignment of length %i" % alignment.get_alignment_length())        
print(dir(alignment))
print(alignment)

"""


from collections import Counter



def seqs_peptides_conserv(data, n_seq, window_size):
    
    AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", '-']
    
    protein_gaps = ''
    conservation_gaps = ''
    
    #alter this to function to every sequence of the alignment
    for i in range(len(data)-n_seq):
        if data[i+n_seq][0] not in AALetter: 
            protein_gaps += data[i]
            conservation_gaps += data[i+n_seq]
    
    peptides = {}
    pos = 0
    pos_gap = []
       
    for i in range(len(protein_gaps)):
        if protein_gaps[i] == '-': pos_gap.append(i)

    protein = '' 
    conservation = ''
    
    for i in range(len(protein_gaps)):
        if i not in pos_gap:
            protein += protein_gaps[i]
            conservation += conservation_gaps[i]
    
    max_conserv = 0
    for i in range(len(protein)-window_size):
        peptide_aux = protein[i:i+window_size]
        conser_aux = conservation[i:i+window_size]
            
        conver_value = 1
        for j in conser_aux:
            if j == '*': conver_value *= 1
            elif j == ':': conver_value *= 0.5
            elif j == '.': conver_value *= 0.25
            elif j == ' ': conver_value *= 0.05 
            
        peptides[str(pos)+'-'+str(pos+window_size)] = [peptide_aux,conser_aux,conver_value]
        pos += 1
            
        if max_conserv<conver_value:
            max_conserv = conver_value
            
    for i in peptides.keys():
        peptides[i][2] = (peptides[i][2] / max_conserv) * window_size
        
    return peptides
    


def conservation_sequence(filename="out_filename.clustal_num", window_size = 15):
    """
    Parameters
    ----------
    filename : FILE, optional
        DESCRIPTION. The default is "out_filename.clustal_num".
    window_size : INTEGER, optional
        DESCRIPTION. The default is 15.

    Returns a dict containing all the peptides and their conservative legend.
    -------
    None.

    """
    
    f = open(filename, "r")  
    
    filename = "out_filename.clustal_num"
    format = "clustal"
    alignment = AlignIO.read(filename, format)

        
    ids_seqs = []
    
    for i in alignment._records:
        ids_seqs.append(i.id)
    

    AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", '-']

    clustal_content = f.readlines()

    data_raw = []

    for i in clustal_content[1:len(clustal_content)]:# it removes the first line, that has Clustal version
        if sum(Counter(list(i)).values()) > 1: 
            data_raw.append(i) #clears all the blank lines
        
    not_seq = True
    count = 0
    while not_seq: #finds the lenght of the string until the start of the sequence
        if data_raw[0][count] not in AALetter:
            count += 1
        else: 
            not_seq = False

    data = []

    for i in data_raw: data.append(i[count:len(i)-1]) #only keeps the sequence / conservation legend    
    
    legend = 0
    legend_found = False
    seqs = {}     

    while not legend_found: 
        if data[legend][0] not in AALetter: 
            legend_found = True
        else: legend += 1
               
    for i in range(0, legend):
        seqs[ids_seqs[legend-i-1]] = seqs_peptides_conserv(data, i+1, window_size)
        

if __name__ == '__main__':
    conservation_sequence()
        




    






