# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 15:19:37 2020

@author: pedro
"""

import pandas as pd
import numpy as np
import os


def size_families_4_min(filename="Seqs_prots.csv"):
    dataset = pd.read_csv(filename)
    
    dataset['Family'] = dataset['Family'].str.strip()
    
    dataset_3_in_family =  pd.DataFrame({'Family' : [], 
                                           'Sequence_fusogenic': [],
                                           'Sequence': []})    
    families_multiple_align = []
    
    for i in dataset.Family.unique():
        dataset_family = dataset.loc[dataset['Family'] == i]
        if dataset_family.shape[0] > 3:
            families_multiple_align.append((i, dataset_family.shape[0]))
    return families_multiple_align


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


def read_weblogo(weblogoDf, fusionpeptide, max_align):
    
    best_score = 0.0
    pos = 0

    
    for i in range(weblogoDf.shape[0]-len(fusionpeptide)):
        current_score = 0.0
        
        for j in range(len(fusionpeptide)):
            #current_score += (weblogoDf.loc[i+j,str(fusionpeptide[j])]/max_align
            #                  )*weblogoDf.loc[i+j,'Entropy']
            current_score += ((np.log2(20) - weblogoDf.loc[i+j,'Entropy']) 
                              * weblogoDf.loc[i+j,str(fusionpeptide[j])]/max_align)
            
            if current_score + np.log2(20) * (len(fusionpeptide) - j) < best_score:
                break
            
        if current_score > best_score:
            best_score = current_score
            pos = i

    return best_score


def weblogo_sequence(sequence, filename='weblogo.txt', n_seqs=5, window_size = 15):
    
    weblogoDf = pd.read_csv(filename, skiprows=7, sep='\t')
    
    weblogoDf = weblogoDf[:-1]
    
    columns = []
    for i in weblogoDf.columns:
        j = i.replace(' ','')
        columns.append(j)
    weblogoDf.columns = columns
    
    # max_align = weblogoDf[columns[1:len(columns)-4]].max().max()
    
    results = {}
    for i in range(len(sequence)-window_size + 1):

        results[str(i)+'-'+str(i+window_size-1)] = read_weblogo(weblogoDf, sequence[i:i+window_size], n_seqs) /window_size        
    
    return results


def weblogo_family(sequence, window_size=15):
    families = size_families_4_min()
    
    seq_dict = {}
    
    print(sequence)
    
    for i in families:
        seq_dict[i[0]] = weblogo_sequence(sequence, i[0] + '_weblogo.txt', i[1], window_size)
        
    return seq_dict


def scores_dataset():
    dataset = vfp_seqs()
    results_output = {}
    for index, row in dataset.iterrows():     
        results_output[row['meta']+ str(index)]=weblogo_family(row['fp'], len(row['fp']))
    return results_output


if __name__ == '__main__':
    results_output=scores_dataset()


