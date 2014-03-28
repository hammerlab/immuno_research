# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
Are ensembles of logistic regression classifiers as good as RF?
"""
import numpy as np


import pandas as pd
import sklearn
import sklearn.cross_validation


from epitopes import \
    (cri_tumor_antigens, iedb, features,
    reduced_alphabet, reference,
    hiv_frahm)

import eval_dataset
from balanced_ensemble import BalancedEnsembleClassifier

cancer_peptides = cri_tumor_antigens.load_peptides(mhc_class = None)

print "Positive validation set size: %d" % len(cancer_peptides)


# sweet over filtering criteria and n-gram transformation
# parameters and for all parameter combos evaluate
# - cross-validation accuracy
# - cross-validation area under ROC curve
# - accuracy on validation set of cancer peptides



best_model = None
best_vectorizer = None
best_params = None

param_count = 0
for assay in  ('cytotoxicity',):
    for alphabet in ('hp2', 'sdm12', None):
        for max_ngram in (2,3,):
            for mhc_class in (1,):
                if alphabet is None:
                    alphabet_dict = None
                    n_letters = 20
                else:
                    alphabet_dict = getattr(reduced_alphabet, alphabet)
                    n_letters = len(set(alphabet_dict.values()))
                n_features = 0
                for i in xrange(max_ngram):
                    n_features += n_letters ** (i+1)

                if n_features > 500:
                    continue
                else:
                    param_count += 1
                param_str =  \
                    "%d: Assay = '%s', ngram %s, alphabet %s, mhc_class %s" % \
                    (param_count, assay, max_ngram, alphabet, mhc_class)
                print param_str

                X, Y, vectorizer = iedb.load_tcell_ngrams(
                    assay_group = assay,
                    human = True,
                    mhc_class = 1,
                    max_ngram = max_ngram,
                    reduced_alphabet = alphabet_dict,
                    min_count = None,
                    return_transformer = True)
                print "Data shape", X.shape, "n_true", np.sum(Y)
                ensemble = BalancedEnsembleClassifier(logistic_regression=True)

                accs = sklearn.cross_validation.cross_val_score(
                    ensemble, X, Y, cv = 3)
                acc = np.mean(accs)
                print "LR CV accuracy %0.4f (std %0.4f)" % \
                    (acc, np.std(accs))

                aucs = sklearn.cross_validation.cross_val_score(
                    ensemble, X, Y, cv = 5, scoring='roc_auc')
                auc = np.mean(aucs)
                print "LR CV AUC %0.4f (std %0.4f)" % \
                    (auc, np.std(aucs))

                ensemble.fit(X, Y)

                X_pos_test = vectorizer.transform(cancer_peptides)
                Y_pos_pred = ensemble.predict(X_pos_test)
                pos_acc = np.mean(Y_pos_pred)
                print "LR Tumor antigen accuracy %0.4f" % (pos_acc,)

                ensemble = \
                    BalancedEnsembleClassifier(logistic_regression=False)

                accs = sklearn.cross_validation.cross_val_score(
                    ensemble, X, Y, cv = 3)
                acc = np.mean(accs)
                print "RF CV accuracy %0.4f (std %0.4f)" % \
                    (acc, np.std(accs))

                aucs = sklearn.cross_validation.cross_val_score(
                    ensemble, X, Y, cv = 5, scoring='roc_auc')
                auc = np.mean(aucs)
                print "RF CV AUC %0.4f (std %0.4f)" % \
                    (auc, np.std(aucs))

                ensemble.fit(X, Y)

                Y_pos_pred = ensemble.predict(X_pos_test)
                pos_acc = np.mean(Y_pos_pred)
                print "RF Tumor antigen accuracy %0.4f" % (pos_acc,)

                print "---"
                print
