import numpy as np
import sklearn
import sklearn.cross_validation
import sklearn.ensemble
import sklearn.svm

from epitopes import features, imma2

print "Loading data and transforming to toxin features"
imm, non = imma2.load_classes()
X, Y = features.toxin_features(imm, non, substring_length = 3, positional=True)

def run_classifiers(X,Y):
  print "Data shape", X.shape
  for c in [0.0001, 0.001, 0.01, 0.1, 1]:
    svm = sklearn.svm.LinearSVC(C=c)
    print "SVM C =", c
    print np.mean(sklearn.cross_validation.cross_val_score(svm, X, Y, cv = 10))

  n_classifiers = 1000
  rf = sklearn.ensemble.RandomForestClassifier(n_classifiers)
  print "Random Forest"
  print np.mean(sklearn.cross_validation.cross_val_score(rf, X, Y, cv = 10))

"""
Toxin features alone seem to do terribly
"""
run_classifiers(X,Y)

