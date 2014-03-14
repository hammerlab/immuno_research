import numpy as np
import sklearn
import sklearn.cross_validation
import sklearn.ensemble
import sklearn.linear_model

from epitopes import iedb, amino_acid, features, reduced_alphabet

import eval_dataset

"""
Do results from a restrict HLA sample (only A2) generalize to all the other HLA types?
"""
A2 = 'A2$|A\*02'

print
print "---"
print "Human MHC1 (keep)"
X_human_mhc1, Y_human_mhc1 = iedb.load_tcell_ngrams(
                 noisy_labels = 'keep',
                 human = True,
                 mhc_class = 1)
eval_dataset.eval_cv(X_human_mhc1, Y_human_mhc1)


print
print "---"
print "Human MHC1 (drop)"
X_human_mhc1_filter, Y_human_mhc1_filter = iedb.load_tcell_ngrams(
                 noisy_labels = 'drop',
                 human = True,
                 mhc_class = 1)
eval_dataset.eval_cv(X_human_mhc1_filter, Y_human_mhc1_filter)


print
print "---"
print "Human MHC1 noisy = positive"
X_human_mhc1_positive, Y_human_mhc1_positive = iedb.load_tcell_ngrams(
                 noisy_labels = 'positive',
                 human = True,
                 mhc_class = 1)
eval_dataset.eval_cv(X_human_mhc1_positive, Y_human_mhc1_positive)

print
print "---"
print "Human MHC1 noisy = negative"
X_human_mhc1_negative, Y_human_mhc1_negative = iedb.load_tcell_ngrams(
                 noisy_labels = 'negative',
                 human = True,
                 mhc_class = 1)
eval_dataset.eval_cv(X_human_mhc1_positive, Y_human_mhc1_positive)

print
print "---"
print "No HLA-A2"
X_no_hla_a2, Y_no_hla_a2 = iedb.load_tcell_ngrams(
                 noisy_labels = 'keep',
                 human = True,
                 mhc_class = 1,
                 exclude_hla_type = A2)
eval_dataset.eval_cv(X_no_hla_a2, Y_no_hla_a2)


print
print "---"
print "No HLA-A2 filtered"
X_no_hla_a2_filter, Y_no_hla_a2_filter = iedb.load_tcell_ngrams(
                 noisy_labels = 'drop',
                 human = True,
                 mhc_class = 1,
                 exclude_hla_type = A2)
eval_dataset.eval_cv(X_no_hla_a2_filter, Y_no_hla_a2_filter)


print
print "---"
print "No HLA-A2 noisy = positive"
X_no_hla_a2_positive, Y_no_hla_a2_positive = iedb.load_tcell_ngrams(
                 noisy_labels = 'positive',
                 human = True,
                 mhc_class = 1,
                 exclude_hla_type = A2)
eval_dataset.eval_cv(X_no_hla_a2_positive, Y_no_hla_a2_positive)



print
print "---"
print "No HLA-A2 noisy = negative"
X_no_hla_a2_negtive, Y_no_hla_a2_negative = iedb.load_tcell_ngrams(
                 noisy_labels = 'negative',
                 human = True,
                 mhc_class = 1,
                 exclude_hla_type = A2)
eval_dataset.eval_cv(X_no_hla_a2_positive, Y_no_hla_a2_positive)


print
print "---"
print "Cross-accuracy for HLA-A2 data"
X_hla_a2, Y_hla_a2 = iedb.load_tcell_ngrams(
                 noisy_labels = 'keep',
                 human = True,
                 mhc_class = 1,
                 hla_type = A2)
eval_dataset.eval_split(X_no_hla_a2, Y_no_hla_a2, X_hla_a2, Y_hla_a2)


print
print "---"
print "Cross-accuracy for HLA-A2 data filtered"
X_hla_a2_filtered, Y_hla_a2_filtered = iedb.load_tcell_ngrams(
                 noisy_labels = 'drop',
                 human = True,
                 mhc_class = 1,
                 hla_type = A2)
eval_dataset.eval_split(X_no_hla_a2_filter, Y_no_hla_a2_filter, X_hla_a2_filtered, Y_hla_a2_filtered)



print
print "---"
print "Cross-accuracy for HLA-A2 data noisy = positive"
X_hla_a2_positive, Y_hla_a2_positive = iedb.load_tcell_ngrams(
                 noisy_labels = 'positive',
                 human = True,
                 mhc_class = 1,
                 hla_type = A2)
eval_dataset.eval_split(X_no_hla_a2_positive, Y_no_hla_a2_positive, X_hla_a2_positive, Y_hla_a2_positive)



print
print "---"
print "Cross-accuracy for HLA-A2 data filtered (assay_group = cytotoxity)"
X_no_hla_a2_cytotoxicity, Y_no_hla_a2_cytotoxicity = iedb.load_tcell_ngrams(
                 noisy_labels = 'drop',
                 assay_group = 'cytotoxicity',
                 human = True,
                 mhc_class = 1,
                 exclude_hla_type = A2)
X_hla_a2_cytotoxicity, Y_hla_a2_cytotoxicity = iedb.load_tcell_ngrams(
                 noisy_labels = 'drop',
                 assay_group = 'cytotoxicity',
                 human = True,
                 mhc_class = 1,
                 hla_type = A2)
eval_dataset.eval_split(X_no_hla_a2_cytotoxicity, Y_no_hla_a2_cytotoxicity, X_hla_a2_cytotoxicity, Y_hla_a2_cytotoxicity)



print
print "---"
print "Cross-accuracy for HLA-A2 data (noisy = positive, assay_group = cytotoxity)"
X_no_hla_a2_positive_cytotoxicity, Y_no_hla_a2_positive_cytotoxicity = iedb.load_tcell_ngrams(
                 noisy_labels = 'positive',
                 assay_group = 'cytotoxicity',
                 human = True,
                 mhc_class = 1,
                 exclude_hla_type = A2)

X_hla_a2_positive_cytotoxicity, Y_hla_a2_positive_cytotoxicity = iedb.load_tcell_ngrams(
                 noisy_labels = 'positive',
                 assay_group = 'cytotoxicity',
                 human = True,
                 mhc_class = 1,
                 hla_type = A2)

eval_dataset.eval_split(X_no_hla_a2_positive_cytotoxicity, Y_no_hla_a2_positive_cytotoxicity, X_hla_a2_positive_cytotoxicity, Y_hla_a2_positive_cytotoxicity)