# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 10:52:28 2020

@author: pedro
"""


import pandas as pd
import operator
import os

def read_weblogo(weblogoDf, max_align, fusionpeptide = 'VSLTLALLLGGLTMGGIAAGV'):
    
    best_score = 0.0
    pos = 0
    
    for i in range(weblogoDf.shape[0]-len(fusionpeptide)):
        current_score = 0.0
        
        for j in range(len(fusionpeptide)):
            current_score += (weblogoDf.loc[i+j,str(fusionpeptide[j])]/max_align
                              )*weblogoDf.loc[i+j,'Entropy']
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
    
    results = {}
    for i in range(len(sequence)-window_size):
        results[str(i)+'-'+str(i+window_size-1)] = read_weblogo(weblogoDf, max_align,
                                                    sequence[i:i+window_size])
    
    # max_key = max(results.items(), key=operator.itemgetter(1))[0]
    # max_score = results[max_key]

    for i in results.keys():
        results[i] = results[i] / window_size
    
    return results
    

def weblogo_family(family, sequence, window_size=15):
    if family.lower() in ['arenaviridae','bornaviridae','bunyaviridae',
                          'coronaviridae','filoviridae',
                          'flaviviridae','herpesviridae','orthomyxoviridae',
                          'paramyxoviridae','peribunyaviridae',
                          'pneumoviridae','retroviridae','togaviridae']:
      return weblogo_sequence(sequence, family + '_weblogo.txt')  
        


if __name__ == '__main__':
    seq = 'FVLGAIALGVATAAAVTAGVAIAKTIRLEGEVAAIKGALRKTNEAVSTLGNGVRVLATAVNDLKDFISKKLTPAINKNKCDISDLKMAVSFGQYNRRFLNVVRQFSDNAGITPAISLDLMTDAELVRAVSNMPTSSGQINLMLENRAMVRRKGFGILIGVYGSSVVYMVQLPIFGVIDTPCWKVKAAPLCSGKDGSYACLLREDQGWYCQNAGSTVYYPNEEDCEVRSDHVFCDTAAGINVAKESEECNRNISTTKYPCKVSTGRHPISMVALSPLGALVACYDGVSCSIGSNKVGIIRPLGKGCSYISNQDADTVTIDNPVYQLSKVEGEQHTIKGKPVSSNFDPIEFPEDQFNIALDQVFESVEKSKNLIDQSNKILDSTEKGNAGFVMVIVLIVLLMLAAVGVGIFFVVKKRKAAPKFPMEMNGVNNKGFIP'
    results = weblogo_family('Pneumoviridae', seq)