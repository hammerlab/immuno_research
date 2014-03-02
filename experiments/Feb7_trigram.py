import numpy as np


import sklearn
import sklearn.cross_validation 
import sklearn.ensemble
import sklearn.linear_model

from ..data import iedb, amino_acid, peptide_to_indices, reduced_alphabet

import eval_dataset

print 
print "---"
print "6-letter unigram"
X6, Y6 = iedb.load_dataset(
                 noisy_labels = 'majority',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 1,
                 reduced_alphabet= reduced_alphabet.alex6)
eval_dataset.eval_cv(X6, Y6)

print 
print "---"
print "6-letter bigram"
X6, Y6 = iedb.load_dataset(
                 noisy_labels = 'majority',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 2,
                 reduced_alphabet= reduced_alphabet.alex6)
eval_dataset.eval_cv(X6, Y6)


print 
print "---"
print "6-letter trigram"
X6, Y6 = iedb.load_dataset(
                 noisy_labels = 'majority',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 3,
                 reduced_alphabet= reduced_alphabet.alex6)
eval_dataset.eval_cv(X6, Y6)


print 
print "---"
print "2-letter trigram"
X2, Y2 = iedb.load_dataset(
                 noisy_labels = 'majority',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 3,
                 reduced_alphabet = reduced_alphabet.hp2)
eval_dataset.eval_cv(X2, Y2)

print 
print "---"
print "3-letter trigram"
X3, Y3 = iedb.load_dataset(
                 noisy_labels = 'majority',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 3,
                 reduced_alphabet= reduced_alphabet.gbmr4)
eval_dataset.eval_cv(X3, Y3)



