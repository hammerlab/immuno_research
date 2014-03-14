import numpy as np


import sklearn
import sklearn.cross_validation
import sklearn.ensemble
import sklearn.linear_model

from epitopes import iedb
import eval_dataset

"""
Instead of dropping or keeping the noisy labels, started
trying to just the majority vote. This is saner and became the default
"""


print
print "---"
print "Human MHC1"
X_human_mhc1_filter, Y_human_mhc1_filter = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority',
                 human = True,
                 mhc_class = 1)
eval_dataset.eval_cv(X_human_mhc1_filter, Y_human_mhc1_filter)



print
print "---"
print "No HLA-A2"
X_no_hla_a2, Y_no_hla_a2 = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority',
                 human = True,
                 mhc_class = 1,
                 exclude_hla_type = 'HLA-A2$|A-\*02')
eval_dataset.eval_cv(X_no_hla_a2, Y_no_hla_a2)


print
print "---"
print "Cross-accuracy for HLA-A2 data"
X_hla_a2, Y_hla_a2 = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority',
                 human = True,
                 mhc_class = 1,
                 hla_type =  'HLA-A2$|A-\*02')
eval_dataset.eval_split(X_no_hla_a2, Y_no_hla_a2, X_hla_a2, Y_hla_a2)





print
print "---"
print "Cross-accuracy for HLA-A2 data filtered (assay_group = cytotoxity)"
X_no_hla_a2_cytotoxicity, Y_no_hla_a2_cytotoxicity = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority',
                 assay_group = 'cytotoxicity',
                 human = True,
                 mhc_class = 1,
                 exclude_hla_type = 'HLA-A2$|A-\*02')
X_hla_a2_cytotoxicity, Y_hla_a2_cytotoxicity = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority',
                 assay_group = 'cytotoxicity',
                 human = True,
                 mhc_class = 1,

                 hla_type = 'HLA-A2$|A-\*02')
eval_dataset.eval_split(X_no_hla_a2_cytotoxicity, Y_no_hla_a2_cytotoxicity, X_hla_a2_cytotoxicity, Y_hla_a2_cytotoxicity)

