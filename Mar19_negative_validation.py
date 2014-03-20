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
The CRI cancer epitopes are a good positive validation set
but need to also make sure we have a low false-positive rate.
Using non-reactive HIV epitopes as a non-immunogenic peptide set
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

cancer_peptides = cri_tumor_antigens.load_peptides(mhc_class = 1)
non_immunogenic_hiv_peptides = \
    hiv_frahm.load_set(max_count = 0)

# self_peptides = reference.load_peptide_set(nrows = 1000)

# sweet over filtering criteria and n-gram transformation
# parameters and for all parameter combos evaluate
# - cross-validation accuracy
# - cross-validation area under ROC curve
# - accuracy on validation set of cancer peptides

params = []
aucs = []
accs = []
recalls = []
d = {

    'assay': [],
    'alphabet' : [],
    'ngram' : [],
    'mhc': [],
    'cv_auc':[],
    'cv_acc':[],
    'pos_acc': [],
    'neg_acc' : [],
    'f1_score':[],
    'f_score': [],
}


best_model = None
best_vectorizer = None
best_params = None

param_count = 0
for assay in ('cytotoxicity', None, 'cytokine release IFNg',  ):
    for alphabet in \
        ('hp2', 'aromatic2', 'gbmr4', 'hp_vs_aromatic', 'sdm12', None):
        for max_ngram in (1, 2, 3):
            for mhc_class in (1,2,None):
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
                d['assay'].append(assay)
                d['alphabet'].append(alphabet)
                d['ngram'].append(max_ngram)
                d['mhc'].append(mhc_class)

                X, Y, vectorizer = iedb.load_tcell_ngrams(
                    assay_group = assay,
                    human = True,
                    mhc_class = 1,
                    max_ngram = max_ngram,
                    reduced_alphabet = alphabet_dict,
                    min_count = None,
                    return_transformer = True)
                print "Data shape", X.shape, "n_true", np.sum(Y)
                ensemble = BalancedEnsembleClassifier()

                accs = sklearn.cross_validation.cross_val_score(
                    ensemble, X, Y, cv = 3)
                acc = np.mean(accs)
                print "CV accuracy %0.4f (std %0.4f)" % \
                    (acc, np.std(accs))
                d['cv_acc'].append(acc)

                aucs = sklearn.cross_validation.cross_val_score(
                    ensemble, X, Y, cv = 5, scoring='roc_auc')
                auc = np.mean(aucs)
                print "CV AUC %0.4f (std %0.4f)" % \
                    (auc, np.std(aucs))
                d['cv_auc'].append(auc)

                ensemble.fit(X, Y)

                X_pos_test = vectorizer.transform(cancer_peptides)
                Y_pos_pred = ensemble.predict(X_pos_test)
                pos_acc = np.mean(Y_pos_pred)
                print "Tumor antigen accuracy %0.4f" % (pos_acc,)
                d['pos_acc'].append(pos_acc)

                X_neg_test = vectorizer.transform(
                    non_immunogenic_hiv_peptides)
                Y_neg_pred = ensemble.predict(X_neg_test)
                neg_acc = 1.0 - np.mean(Y_neg_pred)
                print "Non-immunogenic accuracy %0.4f" % (neg_acc,)
                d['neg_acc'].append(neg_acc)

                n_pos_pred = np.sum(Y_pos_pred)
                n_neg_pred = np.sum(Y_neg_pred)
                precision = n_pos_pred / float(n_pos_pred + n_neg_pred)
                recall = pos_acc
                f1_score = 2 * (precision * recall) / (precision + recall)
                print "F-1 score: %0.4f" % f1_score
                d['f1_score'] = f1_score

                f_half_score = 1.25 * \
                    (precision * recall) / ((0.25 * precision) + recall)

                print "F-0.5 score: %0.4f" % f_half_score
                d['f_score'].append(f_half_score)
                print "---"
                print

df = pd.DataFrame(d)
df = df.sort('f_score', ascending=False)
print df.to_string()
with open('Mar19_negative_validation.csv', 'w') as f:
    df.to_csv(f)
with open('Mar19_negative_validation.html', 'w') as f:
    df.to_html(f)
