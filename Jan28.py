import numpy as np
import sklearn
import sklearn.cross_validation
import sklearn.ensemble
import sklearn.linear_model

from epitopes import iedb, amino_acid, features, reduced_alphabet

"""
Maybe we're overfitting when using larger n-grams of the 20-letter AA alphabet,
can we generalize better using other simpler alphabets?
"""

def run_classifiers(X,Y):
  lr = sklearn.linear_model.LogisticRegression()
  print "LR Accuracy", np.mean(sklearn.cross_validation.cross_val_score(lr, X, Y, cv = 10))

  lr.fit(X,Y)
  print "LR coefs", lr.coef_


  n_classifiers = 200

  rf = sklearn.ensemble.RandomForestClassifier(n_classifiers)
  print "RF Accuracy", np.mean(sklearn.cross_validation.cross_val_score(rf, X, Y, cv = 10))

  rf.fit(X,Y)
  print "RF Features", rf.feature_importances_


print "4 letter alphabet:"
X4,Y4 = iedb.load_tcell_ngrams(
  assay_group = 'cytotoxicity',
  reduced_alphabet = reduced_alphabet.gbmr4,
)
run_classifiers(X4, Y4)

print "---"
print
print "12 letter alphabet:"
X12,Y12 = iedb.load_tcell_ngrams(
  assay_group = 'cytotoxicity',
  reduced_alphabet = reduced_alphabet.sdm12,
)

run_classifiers(X12, Y12)

print "---"
print
print "17 letter alphabet:"
X17,Y17 = iedb.load_tcell_ngrams(
  assay_group = 'cytotoxicity',
  reduced_alphabet = reduced_alphabet.hsdm17,
)

run_classifiers(X17, Y17)


print "---"
print
print "full alphabet:"
X20,Y20 = iedb.load_tcell_ngrams(
  assay_group = 'cytotoxicity',
  reduced_alphabet = None,
)
run_classifiers(X20, Y20)



print "4 letter alphabet pairs:"
X4,Y4 = iedb.load_tcell_ngrams(
  assay_group = 'cytotoxicity',
  reduced_alphabet = reduced_alphabet.gbmr4,
  max_ngram = 2,
)
run_classifiers(X4, Y4)

print "---"
print
print "12 letter alphabet bigrams:"
X12,Y12 = iedb.load_tcell_ngrams(
  assay_group = 'cytotoxicity',
  reduced_alphabet = reduced_alphabet.sdm12,
  max_ngram = 2,
)

run_classifiers(X12, Y12)

print "---"
print
print "17 letter alphabet bigrams:"
X17,Y17 = iedb.load_tcell_ngrams(
  assay_group = 'cytotoxicity',
  reduced_alphabet = reduced_alphabet.hsdm17,
  max_ngram = 2,
)

run_classifiers(X17, Y17)


print "---"
print
print "full alphabet bigrams:"
X20,Y20 = iedb.load_tcell_ngrams(
  assay_group = 'cytotoxicity',
  reduced_alphabet = None,
  max_ngram = 2,
)
run_classifiers(X20, Y20)

