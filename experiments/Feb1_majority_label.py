import numpy as np


import sklearn
import sklearn.cross_validation 
import sklearn.ensemble
import sklearn.linear_model

import eval_dataset 
from ..data import iedb

print 
print "---"
print "Human MHC1"
X_human_mhc1_filter, Y_human_mhc1_filter = iedb.load_dataset(
                 noisy_labels = 'majority',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False)
eval_dataset.eval_cv(X_human_mhc1_filter, Y_human_mhc1_filter)



print 
print "---"
print "No HLA-A2"
X_no_hla_a2, Y_no_hla_a2 = iedb.load_dataset(
                 noisy_labels = 'majority',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = True, 
                 only_hla_a2 = False)
eval_dataset.eval_cv(X_no_hla_a2, Y_no_hla_a2)


print 
print "---"
print "Cross-accuracy for HLA-A2 data"
X_hla_a2, Y_hla_a2 = iedb.load_dataset(
                 noisy_labels = 'majority',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = True)
eval_dataset.eval_split(X_no_hla_a2, Y_no_hla_a2, X_hla_a2, Y_hla_a2)





print 
print "---"
print "Cross-accuracy for HLA-A2 data filtered (assay_group = cytotoxity)"
X_no_hla_a2_cytotoxicity, Y_no_hla_a2_cytotoxicity = iedb.load_dataset(
                 noisy_labels = 'majority',
                 assay_group = 'cytotoxicity', 
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = True, 
                 only_hla_a2 = False)
X_hla_a2_cytotoxicity, Y_hla_a2_cytotoxicity = iedb.load_dataset(
                 noisy_labels = 'majority',
                 assay_group = 'cytotoxicity', 
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = True)
eval_dataset.eval_split(X_no_hla_a2_cytotoxicity, Y_no_hla_a2_cytotoxicity, X_hla_a2_cytotoxicity, Y_hla_a2_cytotoxicity)

