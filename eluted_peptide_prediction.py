#!/usr/bin/python

import pandas as pd
from sklearn import ensemble
from sklearn import linear_model
from sklearn.cross_validation import cross_val_score
import data
import reduced_alphabet

"""
Build models off Dana-Farber Repository for Machine Learning in Immunology

For details of the repository, please refer to the papers below:
Zhang GL, Lin HH, Keskin DB, Reinherz EL, Brusic V. (2011) Dana-Farber repository for machine learning in immunology. J Immunol Methods. 2011; 374(1-2):18-25. 

2nd Machine Learning Competition in Immunology 2012

"""
ELUTED = "eluted"
BINDING = "binding"
NON_BINDING = "nonbinding" 

def get_url(base_url, group, allele):
  return base_url + group + "_" + allele + ".htm"

def get_data(base_url='http://bio.dfci.harvard.edu/DFRMLI/datasets/', allele='HLA-A0201'):
  eluted = pd.read_html(get_url(base_url, ELUTED, allele), infer_types=False, header=0)[0]
  binding = pd.read_html(get_url(base_url, BINDING, allele), infer_types=False, header=0)[0]
  nonbinding = pd.read_html(get_url(base_url, BINDING, allele), infer_types=False, header=0)[0]

  return eluted, binding, nonbinding

if __name__ == '__main__':

  E, B, N = get_data()

  model = ensemble.RandomForestClassifier(n_estimators = 50)

  EB = pd.concat([E, B])
  print EB.count()
  print N.count()

  X, Y = data.make_ngram_dataset(EB.Peptide, N.Peptide, max_ngram=2, normalize_row=True, rebalance=True)

  print cross_val_score(model, X, Y, scoring='roc_auc')
