import numpy as np
import sklearn
import sklearn.cross_validation
import sklearn.ensemble
import sklearn.linear_model

from epitopes import iedb, amino_acid, features

"""
Better performance when filtering the assay group? cytotoxicity looks cleanest
"""

imm, non = iedb.load_tcell_classes(assay_group = 'cytotoxicity')
X_1gram, Y_1gram = features.make_ngram_dataset(imm, non, max_ngram = 1)
X_2gram, Y_2gram = features.make_ngram_dataset(imm, non, max_ngram = 2)

lr = sklearn.linear_model.LogisticRegression()


print "Amino acid unigrams w/ Logistic Regression"
print "LR Accuracy", np.mean(sklearn.cross_validation.cross_val_score(lr, X_1gram, Y_1gram, cv = 10))
lr.fit(X_1gram,Y_1gram)
#print "LR coefs", lr.coef_

print "Amino acid bigrams w/ Logistic Regression"
print "LR Accuracy", np.mean(sklearn.cross_validation.cross_val_score(lr, X_2gram, Y_2gram, cv = 10))
lr.fit(X_2gram,Y_2gram)
#print "LR coefs", lr.coef_


n_classifiers = 200

rf = sklearn.ensemble.RandomForestClassifier(n_classifiers)


print "Amino acid unigram frequency"
aucs = sklearn.cross_validation.cross_val_score(rf, X_1gram, Y_1gram, cv = 10, scoring = 'roc_auc')
print "AUCs", aucs
print "mean AUC", np.mean(aucs)
#rf.fit(X_1gram,Y_1gram)
#print "Features", rf.feature_importances_

print "Amino acid bigram frequency"
aucs = sklearn.cross_validation.cross_val_score(rf, X_2gram, Y_2gram, cv = 10, scoring = 'roc_auc')
print "AUCs", aucs
print "mean AUC", np.mean(aucs)
#rf.fit(X_2gram,Y_2gram)
#print "Features", rf.feature_importances_





