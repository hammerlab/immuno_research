import numpy as np
import sklearn
import sklearn.cross_validation
import sklearn.ensemble
import sklearn.svm

from epitopes import imma2, features

print "Loading data and transforming to toxin features"
imm, non = imma2.load_classes()
X, Y = features.toxin_features(imm, non, substring_length = 2, positional=False)

def run_classifiers(X,Y):
  print "Data shape", X.shape
  for c in [0.0001, 0.001, 0.01]:#, 0.1, 1, 10]:
    svm = sklearn.svm.LinearSVC(C=c)
    print "SVM C =", c
    print np.mean(sklearn.cross_validation.cross_val_score(svm, X, Y, cv = 10))

  n_classifiers = 2000
  rf = sklearn.ensemble.RandomForestClassifier(n_classifiers)
  print "Random Forest"
  print np.mean(sklearn.cross_validation.cross_val_score(rf, X, Y, cv = 10))

"""
Toxin features alone seem to do terribly
"""
run_classifiers(X,Y)

"""
Try reducing the dimensionality
"""
import sklearn.decomposition

print
print "PCA n = 25"
pca = sklearn.decomposition.PCA(n_components=25)
X_small = pca.fit_transform(X)
run_classifiers(X_small,Y)


print
print "PCA n = 50"
pca = sklearn.decomposition.PCA(n_components=50)
X_small = pca.fit_transform(X)
run_classifiers(X_small,Y)


print
print "PCA n = 100"
pca = sklearn.decomposition.PCA(n_components=100)
X_small = pca.fit_transform(X)
run_classifiers(X_small,Y)

print
print "PCA n = 100"
pca = sklearn.decomposition.PCA(n_components=200)
X_small = pca.fit_transform(X)
run_classifiers(X_small,Y)
