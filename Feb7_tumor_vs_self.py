
"""
Conclusion: we can't use an immunogenicity predictor to distinguish potentially immunogenic self epitopes from novel mutant epitopes. This makes sense, since all the source data is just t-cell response assays, without any label for whether it's a useful t-cell response.
"""


from os.path import join

import numpy as np
import sklearn
import sklearn.cross_validation
import sklearn.ensemble
import sklearn.linear_model

from epitopes import (
    cri_tumor_antigens, iedb, features, reduced_alphabet, static_data
)
import eval_dataset

cancer_peptides = cri_tumor_antigens.load_peptides(mhc_class = 1)

self_peptides_file = join(static_data.DATA_DIR, 'Tumor_Self_Antigens_HLA_I.txt')
with open(self_peptides_file) as f:
  self_peptides = f.read().splitlines()

def run(x,y, f):
  x_test_true = f.transform(cancer_peptides)
  x_test_false = f.transform(self_peptides)
  x_test = np.vstack([x_test_true, x_test_false])
  y_test = np.ones(x_test.shape[0], dtype='bool')
  y_test[len(x_test_true):] = 0
  eval_dataset.eval_split(x,y,x_test,y_test)

ASSAY = 'cytotoxicity'

print
print "---"
print "aromatic unigram"
X, Y, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 1,
                 reduced_alphabet= reduced_alphabet.aromatic2,
                 return_transformer = True)

eval_dataset.eval_cv(X, Y)
print "Tumor-specific antigens"
run(X,Y,f)

print
print "---"
print "aromatic bigram"
X, Y, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 2,
                 reduced_alphabet= reduced_alphabet.aromatic2,
                 return_transformer = True)

eval_dataset.eval_cv(X, Y)
print "Tumor-specific antigens"
run(X, Y, f)

print
print "---"
print "aromatic trigram"
X, Y, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 3,
                 reduced_alphabet= reduced_alphabet.aromatic2,
                 return_transformer = True)

eval_dataset.eval_cv(X, Y)
print "Tumor-specific antigens"
run(X, Y, f)

print
print "---"
print "6-letter unigram"
X6, Y6, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 1,
                 reduced_alphabet= reduced_alphabet.alex6,
                 return_transformer = True)

eval_dataset.eval_cv(X6, Y6)
print "Tumor-specific antigens"
run(X6,Y6,f)

print
print "---"
print "6-letter bigram"
X6, Y6, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 2,
                 reduced_alphabet= reduced_alphabet.alex6,
                 return_transformer = True)

eval_dataset.eval_cv(X6, Y6)
print "Tumor-specific antigens"
run(X6, Y6, f)


print
print "---"
print "6-letter trigram"
X6, Y6, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 3,
                 reduced_alphabet= reduced_alphabet.alex6,
                 return_transformer = True)

eval_dataset.eval_cv(X6, Y6)
print "Tumor-specific antigens"
run(X6, Y6, f)


print
print "---"
print "2-letter unigram"
X2, Y2, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 1,
                 reduced_alphabet = reduced_alphabet.hp2,
                 return_transformer = True)

eval_dataset.eval_cv(X2, Y2)
print "Tumor-specific antigens"
run(X2, Y2, f)

print
print "---"
print "2-letter bigram"
X2, Y2, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 2,
                 reduced_alphabet = reduced_alphabet.hp2,
                 return_transformer = True)

eval_dataset.eval_cv(X2, Y2)
print "Tumor-specific antigens"
run(X2, Y2, f)


print
print "---"
print "2-letter trigram"
X2, Y2, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 3,
                 reduced_alphabet = reduced_alphabet.hp2,
                 return_transformer = True)

eval_dataset.eval_cv(X2, Y2)
print "Tumor-specific antigens"
run(X2, Y2, f)

print
print "---"
print "2-letter 4-gram"
X2, Y2, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 4,
                 reduced_alphabet = reduced_alphabet.hp2,
                 return_transformer = True)

eval_dataset.eval_cv(X2, Y2)
print "Tumor-specific antigens"
run(X2, Y2, f)

print
print "---"
print "3-letter unigram"
X3, Y3, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 1,
                 reduced_alphabet= reduced_alphabet.gbmr4,
                 return_transformer = True)

eval_dataset.eval_cv(X3, Y3)
print "Tumor-specific antigens"
run(X3, Y3, f)

print
print "---"
print "3-letter bigram"
X3, Y3, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 2,
                 reduced_alphabet= reduced_alphabet.gbmr4,
                 return_transformer = True)

eval_dataset.eval_cv(X3, Y3)
print "Tumor-specific antigens"
run(X3, Y3, f)

print
print "---"
print "3-letter trigram"
X3, Y3, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 3,
                 reduced_alphabet= reduced_alphabet.gbmr4,
                 return_transformer = True)

eval_dataset.eval_cv(X3, Y3)
print "Tumor-specific antigens"
run(X3, Y3, f)


print
print "---"
print "3-letter 4-gram"
X3, Y3, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 4,
                 reduced_alphabet= reduced_alphabet.gbmr4,
                 return_transformer = True)

eval_dataset.eval_cv(X3, Y3)
print "Tumor-specific antigens"
run(X3, Y3, f)


print
print "---"
print "12-letter unigram"
X12, Y12, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 1,
                 reduced_alphabet= reduced_alphabet.sdm12,
                 return_transformer = True)

eval_dataset.eval_cv(X12, Y12)
print "Tumor-specific antigens"
run(X12, Y12, f)

print
print "---"
print "12-letter bigram"
X12, Y12, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 2,
                 reduced_alphabet= reduced_alphabet.sdm12,
                 return_transformer = True)

eval_dataset.eval_cv(X12, Y12)
print "Tumor-specific antigens"
run(X12, Y12, f)

print
print "---"
print "17-letter unigram"
X17, Y17, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 1,
                 reduced_alphabet= reduced_alphabet.hsdm17,
                 return_transformer = True)

eval_dataset.eval_cv(X17, Y17)
print "Tumor-specific antigens"
run(X17, Y17, f)

print
print "---"
print "17-letter bigram"
X17, Y17, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 2,
                 reduced_alphabet= reduced_alphabet.hsdm17,
                 return_transformer = True)

eval_dataset.eval_cv(X17, Y17)
print "Tumor-specific antigens"
run(X17, Y17, f)


print
print "---"
print "AA unigram"
X, Y, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 1,
                 reduced_alphabet= None,
                 return_transformer = True)

eval_dataset.eval_cv(X, Y)
print "Tumor-specific antigens"
run(X,Y,f)

print
print "---"
print "AA bigram"
X, Y, f = iedb.load_tcell_ngrams(
                 noisy_labels = 'majority', assay_group = ASSAY, subsample_bigger_class = True,
                 human = True,
                 mhc_class = 1,
                 max_ngram = 2,
                 reduced_alphabet= None,
                 return_transformer = True)

eval_dataset.eval_cv(X, Y)
print "Tumor-specific antigens"
run(X, Y, f)




