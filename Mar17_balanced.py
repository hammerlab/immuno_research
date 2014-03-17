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


import sklearn
import sklearn.cross_validation


from epitopes import cri_tumor_antigens, iedb, features, reduced_alphabet

import eval_dataset
from balanced_ensemble import BalancedEnsembleClassifier

cancer_peptides = cri_tumor_antigens.load_peptides(mhc_class = 1)


# sweet over filtering criteria and n-gram transformation
# parameters and for all parameter combos evaluate
# - cross-validation accuracy
# - cross-validation area under ROC curve
# - accuracy on validation set of cancer peptides

for assay in ('cytotoxicity', None ):
    for mhc_class in (1, None):
        for alphabet in (reduced_alphabet.hp2, None):
            for max_ngram in (1, 2, 3):
                param_str =  \
                    "Assay %s, MHC Class %s, n-gram %s, HP2 %s" % \
                    (assay, mhc_class, max_ngram, alphabet is not None)
                print param_str
                X, Y, vectorizer = iedb.load_tcell_ngrams(
                    assay_group = assay,
                    human = True,
                    mhc_class = mhc_class,
                    max_ngram = max_ngram,
                    reduced_alphabet= alphabet,
                    return_transformer = True)


                ensemble = BalancedEnsembleClassifier()


                accs = sklearn.cross_validation.cross_val_score(
                    ensemble, X, Y, cv = 5)
                print "CV accuracy %s (std %s)" % (np.mean(accs), np.std(accs))

                aucs = sklearn.cross_validation.cross_val_score(
                    ensemble, X, Y, cv = 5, scoring='roc_auc')
                print "CV AUC %s (std %s)" % (np.mean(aucs), np.std(aucs))

                X_test = vectorizer.transform(cancer_peptides)
                Y_test = np.array([True] * len(cancer_peptides))
                ensemble.fit(X, Y)
                Y_pred = ensemble.predict(X_test)
                print "Validation accuracy", np.mean(Y_pred == Y_test)
                print "---"
                print
