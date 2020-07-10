# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 10:45:20 2020

@author: pedro
"""

# from datasets_diana import *
# from readWeblogo import *
import pandas as pd
import numpy as np

def vfp_seqs(filename="dataset1_usar.csv", families = "Seqs_prots_fp.csv"):
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


def read_weblogo(weblogoDf, max_align, fusionpeptide = 'VSLTLALLLGGLTMGGIAAGV'):
    
    best_score = 0.0
    pos = 0
    
    for i in range(weblogoDf.shape[0]-len(fusionpeptide)):
        current_score = 0.0
        
        for j in range(len(fusionpeptide)):
            #current_score += (weblogoDf.loc[i+j,str(fusionpeptide[j])]/max_align
            #                  )*weblogoDf.loc[i+j,'Entropy']
            current_score += np.log2(20) - weblogoDf.loc[i+j,'Entropy']
            
        
        
        
        if current_score > best_score:
            best_score = current_score
            pos = i

    return best_score
    
    
def weblogo_sequence(sequence, filename='weblogo.txt', window_size = 15):
    
    weblogoDf = pd.read_csv(filename, skiprows=7, sep='\t')
    
    weblogoDf = weblogoDf[:-1]
    
    # print(weblogoDf)
    
    columns = []
    for i in weblogoDf.columns:
        j = i.replace(' ','')
        columns.append(j)
    weblogoDf.columns = columns
    
    max_align = weblogoDf[columns[1:len(columns)-4]].max().max()
    
    results = 0
    for i in range(len(sequence)-window_size + 1):
        #results[str(i)+'-'+str(i+window_size-1)] = read_weblogo(weblogoDf, max_align,
        #                                            sequence[i:i+window_size])
        results = read_weblogo(weblogoDf, max_align, sequence[i:i+window_size]) /window_size
    
    #max_key = max(results.items(), key=operator.itemgetter(1))[0]
    #max_score = results[max_key]

    #for i in results.keys():
    #    results[i] = results[i] / max_score
    
    #for i in results.keys():
    #    results[i] = results[i] / window_size
        
    
    return results
    

def weblogo_family(family, sequence, window_size=15):
    if family.lower() in ['arenaviridae','bornaviridae','bunyaviridae',
                          'coronaviridae','filoviridae',
                          'flaviviridae','herpesviridae','orthomyxoviridae',
                          'paramyxoviridae','peribunyaviridae',
                          'pneumoviridae','retroviridae','togaviridae']:
        
        
        print(sequence)
        return weblogo_sequence(sequence, family + '_weblogo.txt', len(sequence))

    else: return np.nan
    

def weblogo_all(sequence):
    return weblogo_sequence(sequence, 'all_family_weblogo.txt', len(sequence))



def scores_dataset():
    dataset = vfp_seqs()
    results_output = {}
    for index, row in dataset.iterrows():        
        # results_output[row['meta']+ str(index)]=weblogo_family(row['meta'][0].upper() + row['meta'][1:], row['fp'])
        results_output[row['meta']+ str(index)]=weblogo_all(row['fp'])
    return results_output

if __name__ == '__main__':
    results_output = scores_dataset()
