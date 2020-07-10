# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 12:01:31 2020

@author: pedro
"""

import pandas as pd
import numpy as np

def weblogo_sequence(sequence, filename='weblogo.txt', n_seqs=5, window_size = 15):
    
    weblogoDf = pd.read_csv(filename, skiprows=7, sep='\t')
    
    weblogoDf = weblogoDf[:-1]
    
    columns = []
    for i in weblogoDf.columns:
        j = i.replace(' ','')
        columns.append(j)
    weblogoDf.columns = columns
    
    weblogo_entropyes = weblogoDf.loc[:, weblogoDf.columns[1:len(weblogoDf.columns)-4]]
    
    entropies = list((np.log2(20) - weblogoDf.loc[:,'Entropy']) / n_seqs)
    
    weblogo_entropyes = weblogo_entropyes.mul(entropies, axis=0)
    
    
    """
    for i in range(weblogo_entropyes.shape[0]):
        for j in range(1, weblogo_entropyes.shape[1]): #excludes label line
            weblogo_entropyes.iloc[i,j] =  (np.log2(20) - weblogoDf.loc[i,'Entropy']) * weblogo_entropyes.iloc[i,j] / n_seqs
        
    """   
        
    
    return  weblogo_entropyes
    
    
    
    
    
    """
    # max_align = weblogoDf[columns[1:len(columns)-4]].max().max()
    
    results = {}
    for i in range(len(sequence)-window_size + 1):

        results[str(i)+'-'+str(i+window_size-1)] = read_weblogo(weblogoDf, sequence[i:i+window_size], n_seqs) /window_size        
    
    return results
    """
    
print(weblogo_sequence('', 'arenaviridae_weblogo.txt'))