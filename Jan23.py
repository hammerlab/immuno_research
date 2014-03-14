import numpy as np
import sklearn
import sklearn.cross_validation
import sklearn.ensemble

from epitopes import iedb, features, amino_acid

imm, non = iedb.load_tcell_classes(peptide_length = 9)
X, Y = features.make_kmer_dataset(imm, non)
X_1gram, Y_1gram = features.make_ngram_dataset(imm, non, max_ngram = 1)
X_2gram, Y_2gram = features.make_ngram_dataset(imm, non, max_ngram = 2)

n_classifiers = 100

clf = sklearn.ensemble.RandomForestClassifier(n_classifiers)

print "Amino acid 9mers"
print np.mean(sklearn.cross_validation.cross_val_score(clf, X, Y, cv = 10, scoring = 'roc_auc'))

print "Amino acid unigram frequency"
print np.mean(sklearn.cross_validation.cross_val_score(clf, X_1gram, Y_1gram, cv = 10, scoring = 'roc_auc'))

print "Amino acid bigram frequency"
print np.mean(sklearn.cross_validation.cross_val_score(clf, X_2gram, Y_2gram, cv = 10, scoring = 'roc_auc'))
