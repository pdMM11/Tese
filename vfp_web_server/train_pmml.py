# -*- coding: utf-8 -*-
"""
Created on Wed May 27 14:50:30 2020

@author: pedro
"""

from Propythia_ML_PMML import MachineLearning_PMML
from sklearn2pmml import sklearn2pmml
import pandas as pd
from propythia.machine_learning import MachineLearning
import pickle

def training():
    filename = 'crmapp/ml_models/dataset3_all_'
    dataset_in = r'dataset3_all_svc.csv'
    dataset = pd.read_csv(dataset_in, delimiter=',')
    x_original = dataset.loc[:, dataset.columns != 'labels']
    labels = dataset['labels']
    # create Machine learning object
    # ml = MachineLearning(x_original, labels, classes=['non_vfp', 'vfp'])
    ml = MachineLearning(x_original, labels, classes=['non_vfp', 'vfp'])
    best_model = ml.train_best_model('svm')

    # save the model to disk
    filename2 = filename + 'svm_model_13062020.sav'
    # pickle.dump(best_model, open(filename2, 'wb'))
    # sklearn2pmml(best_model, "Model_pmml_svm_14062020.pmml", with_repr=True)
    
    print(ml)
    
training()
