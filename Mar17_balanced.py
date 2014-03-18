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
    (cri_tumor_antigens, iedb, features, reduced_alphabet, reference)

import eval_dataset
from balanced_ensemble import BalancedEnsembleClassifier

cancer_peptides = cri_tumor_antigens.load_peptides(mhc_class = 1)

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
d = {'param' : [], 'auc':[], 'acc':[], 'recall': []}

for assay in ('cytotoxicity', None ):
    for alphabet in (None, 'hp2', 'aromatic2', 'gbmr4'):
        for max_ngram in (1, 2, 3):
            for min_count in (None, 3):
                if alphabet is None and max_ngram > 1:
                    continue
                param_str =  \
                    "Assay %s, ngram %s, min_count %s, alphabet %s" % \
                    (assay, max_ngram, min_count,
                     alphabet)
                d['param'].append(param_str)
                print param_str
                if alphabet == 'hp2':
                    alphabet_dict =  reduced_alphabet.hp2
                elif alphabet == 'aromatic2':
                    alphabet_dict = reduced_alphabet.aromatic2
                elif alphabet == 'gbmr4':
                    alphabet_dict = reduced_alphabet.gbmr4
                else:
                    assert alphabet is None, alphabet
                    alphabet_dict = None
                X, Y, vectorizer = iedb.load_tcell_ngrams(
                    assay_group = assay,
                    human = True,
                    mhc_class = 1,
                    max_ngram = max_ngram,
                    reduced_alphabet = alphabet_dict,
                    min_count = min_count,
                    return_transformer = True)
                print "Data shape", X.shape, "n_true", np.sum(Y)
                ensemble = BalancedEnsembleClassifier()

                accs = sklearn.cross_validation.cross_val_score(
                    ensemble, X, Y, cv = 5)
                print "CV accuracy %0.4f (std %0.4f)" % \
                    (np.mean(accs), np.std(accs))
                d['acc'].append(np.mean(accs))

                aucs = sklearn.cross_validation.cross_val_score(
                    ensemble, X, Y, cv = 5, scoring='roc_auc')
                print "CV AUC %0.4f (std %0.4f)" % \
                    (np.mean(aucs), np.std(aucs))
                d['auc'].append(np.mean(aucs))

                ensemble.fit(X, Y)

                #X_self = vectorizer.transform(self_peptides)
                #Y_pred = ensemble.predict(X_self)
                #print "Self epitope accuracy %0.4f" % \
                #    (1.0 - np.mean(Y_pred))
                X_test = vectorizer.transform(cancer_peptides)
                Y_pred = ensemble.predict(X_test)
                recall = np.mean(Y_pred)
                print "Tumor antigen accuracy %0.4f" % (recall,)
                d['recall'].append(recall)
                print "---"
                print

df = pd.DataFrame(d)
print df.sort('auc')



