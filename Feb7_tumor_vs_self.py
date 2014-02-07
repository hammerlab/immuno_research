import numpy as np


import sklearn
import sklearn.cross_validation 
import sklearn.ensemble
import sklearn.linear_model

import data 
import amino_acid
import iedb
import reduced_alphabet 
import eval_dataset


"""
Conclusion: we can't use an immunogenicity predictor to distinguish potentially immunogenic self epitopes from novel mutant epitopes. This makes sense, since all the source data is just t-cell response assays, without any label for whether it's a useful t-cell response. 
"""

with open("Tumor_Mutant_Antigens_HLA_I.txt") as f:
  cancer_peptides = f.read().splitlines()

with open("Tumor_Self_Antigens_HLA_I.txt") as f:
  self_peptides = f.read().splitlines()


def run(x,y, f):
  x_test_true = f(cancer_peptides)
  x_test_false = f(self_peptides)
  x_test = np.vstack([x_test_true, x_test_false])
  y_test = np.ones(x_test.shape[0], dtype='bool')
  y_test[len(x_test_true):] = 0
  eval_dataset.eval_split(x,y,x_test,y_test)
  
ASSAY = 'cytotoxicity'


print 
print "---"
print "aromatic unigram"
X, Y, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 1,
                 reduced_alphabet= reduced_alphabet.aromatic2,
                 return_transformer = True)

eval_dataset.eval_cv(X, Y)
print "Tumor-specific antigens"
run(X,Y,f)

print 
print "---"
print "aromatic bigram"
X, Y, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 2,
                 reduced_alphabet= reduced_alphabet.aromatic2,
                 return_transformer = True)

eval_dataset.eval_cv(X, Y)
print "Tumor-specific antigens"
run(X, Y, f)

print 
print "---"
print "aromatic trigram"
X, Y, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 3,
                 reduced_alphabet= reduced_alphabet.aromatic2,
                 return_transformer = True)

eval_dataset.eval_cv(X, Y)
print "Tumor-specific antigens"
run(X, Y, f)

print 
print "---"
print "6-letter unigram"
X6, Y6, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 1,
                 reduced_alphabet= reduced_alphabet.alex6,
                 return_transformer = True)

eval_dataset.eval_cv(X6, Y6)
print "Tumor-specific antigens"
run(X6,Y6,f)

print 
print "---"
print "6-letter bigram"
X6, Y6, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 2,
                 reduced_alphabet= reduced_alphabet.alex6,
                 return_transformer = True)

eval_dataset.eval_cv(X6, Y6)
print "Tumor-specific antigens"
run(X6, Y6, f)


print 
print "---"
print "6-letter trigram"
X6, Y6, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 3,
                 reduced_alphabet= reduced_alphabet.alex6,
                 return_transformer = True)

eval_dataset.eval_cv(X6, Y6)
print "Tumor-specific antigens"
run(X6, Y6, f)


print 
print "---"
print "2-letter unigram"
X2, Y2, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 1,
                 reduced_alphabet = reduced_alphabet.hp2,
                 return_transformer = True)

eval_dataset.eval_cv(X2, Y2)
print "Tumor-specific antigens"
run(X2, Y2, f)

print 
print "---"
print "2-letter bigram"
X2, Y2, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 2,
                 reduced_alphabet = reduced_alphabet.hp2,
                 return_transformer = True)

eval_dataset.eval_cv(X2, Y2)
print "Tumor-specific antigens"
run(X2, Y2, f)


print 
print "---"
print "2-letter trigram"
X2, Y2, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 3,
                 reduced_alphabet = reduced_alphabet.hp2,
                 return_transformer = True)

eval_dataset.eval_cv(X2, Y2)
print "Tumor-specific antigens"
run(X2, Y2, f)

print 
print "---"
print "2-letter 4-gram"
X2, Y2, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 4,
                 reduced_alphabet = reduced_alphabet.hp2,
                 return_transformer = True)

eval_dataset.eval_cv(X2, Y2)
print "Tumor-specific antigens"
run(X2, Y2, f)

print 
print "---"
print "3-letter unigram"
X3, Y3, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 1,
                 reduced_alphabet= reduced_alphabet.gbmr4,
                 return_transformer = True)

eval_dataset.eval_cv(X3, Y3)
print "Tumor-specific antigens"
run(X3, Y3, f)

print 
print "---"
print "3-letter bigram"
X3, Y3, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 2,
                 reduced_alphabet= reduced_alphabet.gbmr4,
                 return_transformer = True)

eval_dataset.eval_cv(X3, Y3)
print "Tumor-specific antigens"
run(X3, Y3, f)

print 
print "---"
print "3-letter trigram"
X3, Y3, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 3,
                 reduced_alphabet= reduced_alphabet.gbmr4,
                 return_transformer = True)

eval_dataset.eval_cv(X3, Y3)
print "Tumor-specific antigens"
run(X3, Y3, f)


print 
print "---"
print "3-letter 4-gram"
X3, Y3, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 4,
                 reduced_alphabet= reduced_alphabet.gbmr4,
                 return_transformer = True)

eval_dataset.eval_cv(X3, Y3)
print "Tumor-specific antigens"
run(X3, Y3, f)


print 
print "---"
print "12-letter unigram"
X12, Y12, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 1,
                 reduced_alphabet= reduced_alphabet.sdm12,
                 return_transformer = True)

eval_dataset.eval_cv(X12, Y12)
print "Tumor-specific antigens"
run(X12, Y12, f)

print 
print "---"
print "12-letter bigram"
X12, Y12, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 2,
                 reduced_alphabet= reduced_alphabet.sdm12,
                 return_transformer = True)

eval_dataset.eval_cv(X12, Y12)
print "Tumor-specific antigens"
run(X12, Y12, f)

print 
print "---"
print "17-letter unigram"
X17, Y17, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 1,
                 reduced_alphabet= reduced_alphabet.hsdm17,
                 return_transformer = True)

eval_dataset.eval_cv(X17, Y17)
print "Tumor-specific antigens"
run(X17, Y17, f)

print 
print "---"
print "17-letter bigram"
X17, Y17, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 2,
                 reduced_alphabet= reduced_alphabet.hsdm17,
                 return_transformer = True)

eval_dataset.eval_cv(X17, Y17)
print "Tumor-specific antigens"
run(X17, Y17, f)


print 
print "---"
print "AA unigram"
X, Y, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 1,
                 reduced_alphabet= None,
                 return_transformer = True)

eval_dataset.eval_cv(X, Y)
print "Tumor-specific antigens"
run(X,Y,f)

print 
print "---"
print "AA bigram"
X, Y, f = iedb.load_dataset(
                 noisy_labels = 'majority', assay_group = ASSAY, rebalance = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = 2,
                 reduced_alphabet= None,
                 return_transformer = True)

eval_dataset.eval_cv(X, Y)
print "Tumor-specific antigens"
run(X, Y, f)


