# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 10:34:07 2020

@author: pedro
"""

import pickle



def import_pairwise():
    with open('pairwise_results.pkl', 'rb') as handle:
        results_pair = pickle.load(handle)
        return results_pair
    
    
def best_alignements(results):
    
    best_align = []
    
    for i in results.keys():
        for j in results[i].keys():
            if results[i][j][0] >= 2.0:
                best_align.append((results[i][j][0], 
                                   results[i][j][1], 
                                   results[i][j][2],
                                   j))
                
    best_align_sorted = sorted(best_align,reverse=True)
    for i in best_align_sorted: 
        #if i[1]==i[2]: print(i)
        # print(i)
        print(i)
        
        
def differences_between_seqs(seq1, seq2):
    counter=0
    for i in range(len(seq1)):
        if seq1==seq2: counter += 1
    return counter
    


if __name__ == '__main__':
    results_pair = import_pairwise()
    best_alignements(results_pair)
    