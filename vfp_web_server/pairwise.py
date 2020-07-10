# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 16:26:19 2020

@author: pedro
"""

import pandas as pd
import numpy as np


# Import pairwise2 module
from Bio import pairwise2

# Import format_alignment method
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

import pickle



def vfp_seqs(filename="dataset1_usar.csv", families = "Seqs_prots.csv"):
    colnames=['meta', 'fp'] 
    dataset = pd.read_csv(filename, names=colnames, header=None)
    
    dataset_family = pd.read_csv(families)
    
    dataset_family['Family'] = dataset_family['Family'].str.strip()
    
    families = [x.lower() for x in list(dataset_family['Family'].unique())]
    
    for i in range(len(dataset['meta'])):
        for j in families:
            if j in str.lower(dataset['meta'][i]):
                dataset.loc[i, 'meta'] = j
                break
        if '_' in dataset['meta'][i]:
            temp = dataset['meta'][i].split('_')
            dataset['meta'][i] = str.lower(temp[len(temp)-1])
            
    dataset = dataset.sort_values(by=['meta'])

    
    return dataset



def fps_per_family(dataset):
    families = {}
    for i in dataset['meta'].unique():
        families[i] = dataset.loc[dataset['meta'] == i, 'fp']
    matrix = matlist.blosum62   
    
    results = {}
    scores = {}
    
    for i in families.keys():
        results[i] = {}
        scores[i] = {}
        
        #print(results)  
        
        current_family = families[i]
        indexes = current_family.index 

        for j in range(len(current_family)-1):
            for k in range(j+1,len(current_family)):
                
                if indexes[j] < indexes[k]:
                    
                    results[i][str(indexes[j])+'-'+str(indexes[k])]=pairwise2.align.globaldx(
                        current_family.iloc[j],current_family.iloc[k], matrix)
                    
                    scores[i][str(indexes[j])+'-'+str(indexes[k])] = (
                        results[i][str(indexes[j])+'-'+str(indexes[k])][0][2] /
                        results[i][str(indexes[j])+'-'+str(indexes[k])][0][4] ,
                        results[i][str(indexes[j])+'-'+str(indexes[k])][0][0],
                        results[i][str(indexes[j])+'-'+str(indexes[k])][0][1]
                        )
                    
                else: 
                    results[i][str(indexes[k])+'-'+str(indexes[j])]=pairwise2.align.globaldx(
                        current_family.iloc[k],current_family.iloc[j], matrix)
                    
                    scores[i][str(indexes[k])+'-'+str(indexes[j])] = (
                        results[i][str(indexes[k])+'-'+str(indexes[j])][0][2] /
                        results[i][str(indexes[k])+'-'+str(indexes[j])][0][4] ,
                        results[i][str(indexes[k])+'-'+str(indexes[j])][0][0],
                        results[i][str(indexes[k])+'-'+str(indexes[j])][0][1]
                        )
                    
    return results, scores
                
if __name__ == '__main__':
    results, scores  = fps_per_family(vfp_seqs())
    f = open("pairwise_results.pkl","wb")
    pickle.dump(scores,f, protocol=pickle.HIGHEST_PROTOCOL)
    f.close()

#    results = pickle.load(open('pairwise_results.pkl', 'rb'))
    
    
    

#python emboss_needle.py --email pedrodmmoreira@gmail.com --stype protein --asequence ASHTTLGVQDISTTAMSW --bsequence TSTPIVVDCSTYVCNGNVRCV