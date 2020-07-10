# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 11:18:54 2020

@author: pedro
"""

import pandas as pd
import numpy as np
import os

import time


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

def weblogo_family():
    families = size_families_4_min()
    
    # print(sorted(families))
    
    family_weblogo = {}
    
    for i in families:
        family = i[0]
        filename = family + '_weblogo.txt'
        n_seqs = i[1]
        
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
        
        family_weblogo[family] = weblogo_entropyes
        
    return family_weblogo


def read_weblogo(weblogoDf, fusionpeptide, positions):
    
    best_score = 0.0
    pos = 0
    
    tam_check = len(fusionpeptide) / 5

    for i in positions:
        
        if i + len(fusionpeptide) + 1 >= weblogoDf.shape[0]: pass
        else:
            current_score = 0.0
            
            for j in range(len(fusionpeptide)):
                
                current_score += weblogoDf.loc[i+j,str(fusionpeptide[j])]
                if (j-i) > tam_check and current_score == 0.0: break
                else:
                    if current_score + np.log2(20) * (len(fusionpeptide) - j) < best_score:
                        break
            
            if current_score > best_score:
                best_score = current_score
                pos = i

    return best_score / len(fusionpeptide)


def max_coluns(dataset):
    indexes_top = {}
    for i in dataset.columns:
        indexes_top[i] = list(dataset.nlargest(20,i)[i].index)
    return indexes_top
    

def scores_dataset():
    
    # start_time = time.time()
    
    
    dataset = vfp_seqs()
    results_output = {}
    family_weblogo = weblogo_family()
    family_column_max_index = {}
    for i in family_weblogo.keys():
        family_column_max_index[i] = max_coluns(family_weblogo[i])
    # print(family_column_max_index['Rhabdoviridae'])
    
    for index, row in dataset.iterrows():
        dict_seq = {}
        sequence=row['fp']
        # print(sequence)
    
        # print("--- %s seconds ---" % (time.time() - start_time))
        for i in family_weblogo.keys():
            current_family = family_column_max_index[i]
    
            pos_analysis = []
    
            for j in range(0,len(sequence)):
                if j == 0: pos_analysis = current_family[sequence[j]]
                else:
                    pos_temp = list(np.array(current_family[sequence[j]]) - j)
                    pos_temp = np.array(pos_temp)
                    pos_temp = pos_temp[pos_temp >=0]
                    pos_analysis=list(sorted(pos_analysis + list(set(pos_temp) - set(pos_analysis))))

            dict_seq[i]=read_weblogo(family_weblogo[i], sequence, pos_analysis)
        results_output[row['meta']+ str(index)]= dict_seq
    
    # print("FINAL: --- %s seconds ---" % (time.time() - start_time))
    return results_output



def scores_sequence(protein, window_size = 15):

    # dataset = vfp_seqs()
    results_output = {}
    family_weblogo = weblogo_family()
    family_column_max_index = {}
    for i in family_weblogo.keys():
        family_column_max_index[i] = max_coluns(family_weblogo[i])

    for i in range(len(protein) - window_size):
        pos_end = i+window_size
        sequence = protein[i:pos_end]
        dict_seq = {}
        #sequence=row['fp']

        for k in family_weblogo.keys():
            current_family = family_column_max_index[k]
    
            pos_analysis = []
    
            for j in range(0,len(sequence)):
                if j == 0: pos_analysis = current_family[sequence[j]]
                else:
                    pos_temp = list(np.array(current_family[sequence[j]]) - j)
                    pos_temp = np.array(pos_temp)
                    pos_temp = pos_temp[pos_temp >=0]
                    pos_analysis=list(sorted(pos_analysis + list(set(pos_temp) - set(pos_analysis))))

            dict_seq[k]=read_weblogo(family_weblogo[k], sequence, pos_analysis)
        results_output[str(i) + '-' + str(pos_end)]= dict_seq
    
    # print("FINAL: --- %s seconds ---" % (time.time() - start_time))
    return results_output



sequence = 'EPVSLTLALLLGGLTMGGIAAGVGTGTTALVATQQFQQLQAAMHDDLKEVEKSITNLEKSLTSLSEVVLQNRRGLDLLFLKEGGLCAALKEECCFYADHTGLVRDSMAKLRERLSQRQKLFESQQGWFEGLFNKSPWFTTLISTIMGPLIILLLILLFGPCILNRLVQFIKDRISVVQAL'
results_output = scores_sequence(sequence)