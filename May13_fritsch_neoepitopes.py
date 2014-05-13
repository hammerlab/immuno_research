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


import numpy as np


import pandas as pd
import sklearn
import sklearn.cross_validation


from epitopes import \
    (fritsch_neoepitopes, iedb, features,
    reduced_alphabet, reference,
    hiv_frahm)

import eval_dataset
from balanced_ensemble import BalancedEnsembleClassifier

df = fritsch_neoepitopes.load_dataframe()
pos_peptides = df['Mutated Epitope']
neg_peptides = df['Native Epitope']


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
    'min_count': [],
    'pos_acc': [],
    'neg_acc' : [],
    'f1_score':[],
    'f_score': [],
}


best_model = None
best_vectorizer = None
best_params = None

param_count = 0
for assay in ('cytotoxicity', None):
    for alphabet in ('gbmr4', 'sdm12', 'hsdm17', None):
        for max_ngram in (2, 3):
            for mhc_class in (1, None):
                for min_count in (None,3,7):
                    if alphabet is None:
                        alphabet_dict = None
                        n_letters = 20
                    else:
                        alphabet_dict = getattr(reduced_alphabet, alphabet)
                        n_letters = len(set(alphabet_dict.values()))
                    n_features = 0
                    for i in xrange(max_ngram):
                        n_features += n_letters ** (i+1)

                    if n_features > 200:
                        continue
                    else:
                        param_count += 1
                    param_str =  \
                        "%d: Assay = '%s', ngram %s, min_count %s, alphabet %s, mhc_class %s" % \
                        (param_count, assay, max_ngram, min_count, alphabet, mhc_class)
                    print param_str
                    d['assay'].append(assay)
                    d['alphabet'].append(alphabet)
                    d['ngram'].append(max_ngram)
                    d['mhc'].append(mhc_class)
                    d['min_count'].append(min_count)

                    X, Y, vectorizer = iedb.load_tcell_ngrams(
                        assay_group = assay,
                        human = True,
                        mhc_class = mhc_class,
                        max_ngram = max_ngram,
                        reduced_alphabet = alphabet_dict,
                        min_count = min_count,
                        return_transformer = True)

                    assert len(X) > 0, "No samples!"
                    assert len(Y) > 0, "No labels!"
                    assert np.sum(Y) > 0, "No positive samples!"
                    assert np.sum(Y) < len(Y), "No negative samples!"

                    print "Data shape", X.shape, "n_true", np.sum(Y)
                    ensemble = BalancedEnsembleClassifier()
                    ensemble.fit(X, Y)

                    X_pos_test = vectorizer.transform(pos_peptides)
                    Y_pos_pred = ensemble.predict(X_pos_test)
                    pos_acc = np.mean(Y_pos_pred)
                    print "Tumor antigen accuracy %0.4f" % (pos_acc,)
                    d['pos_acc'].append(pos_acc)

                    X_neg_test = vectorizer.transform(neg_peptides)
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
                    d['f1_score'].append(f1_score)

                    f_half_score = 1.25 * \
                        (precision * recall) / ((0.25 * precision) + recall)

                    print "F-0.5 score: %0.4f" % f_half_score
                    d['f_score'].append(f_half_score)
                    print "---"
                    print

df = pd.DataFrame(d)
df = df.sort('f_score', ascending=False)
print df.to_string()
with open('May13_validation.csv', 'w') as f:
    df.to_csv(f)
with open('May13_validation.html', 'w') as f:
    df.to_html(f)
