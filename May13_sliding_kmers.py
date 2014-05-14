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
Instead of looking at positionally invariant n-grams (n=1,2,3),
let's split each string into 8-mers and treat each one as a 
distinct sample (each weighted inversely with number of 8-mers
originating from single peptide).
"""

import numpy as np


import pandas as pd
import sklearn
import sklearn.cross_validation
from sklearn.ensemble import RandomForestClassifier

from epitopes import \
    (fritsch_neoepitopes, iedb, features,
    reduced_alphabet, reference,
    hiv_frahm)

import eval_dataset
from balanced_ensemble import BalancedEnsembleClassifier

df = fritsch_neoepitopes.load_dataframe()
pos_peptides = df['Mutated Epitope']
neg_peptides = df['Native Epitope']

kmer_length = 9 

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
    'mhc': [],
    'min_count': [],
    'pos_acc': [],
    'neg_acc' : [],
    'f1_score':[],
    'f_score': [],
    'acc' : [],
    'cv_auc': [], 
}


best_model = None
best_vectorizer = None
best_params = None

param_count = 0
for assay in ('cytotoxicity', None, ):
    for mhc_class in (1, None):
        for min_count in (3, 5, 7,  None):

            imm, non = iedb.load_tcell_classes(
                assay_group = assay,
                human = True,
                mhc_class = mhc_class,
                min_count = min_count)

            for alphabet in \
                    ('hp2', 'gbmr4', 'hp_vs_aromatic', 'sdm12', 'hsdm17'):
                
                transformer = reduced_alphabet.make_alphabet_transformer(alphabet)
                param_str =  \
                    "%d: Assay = '%s', min_count %s, alphabet %s, mhc_class %s" % \
                    (param_count, assay, min_count, alphabet, mhc_class)
                print param_str
                
                d['assay'].append(assay)
                d['alphabet'].append(alphabet)
                d['mhc'].append(mhc_class)
                d['min_count'].append(min_count)
                param_count += 1
                
                X_combined = []
                Y_combined = []
                W_combined = []
                dtype = 'S%d' % kmer_length
                    
                def expand(peptide_set):
                    reduced = [transformer.transform(p) for p in peptide_set]
                    X = []
                    Counts = []
                    Indices = []
                    for peptide_idx, p in enumerate(reduced):
                        n = len(p)
                        if n < kmer_length:
                            continue
                        n_substrings = n - kmer_length + 1
                        for i in xrange(0, n_substrings):
                            substr = p[i:i+kmer_length]
                            X.append(substr)
                            Counts.append(n_substrings)
                            Indices.append(peptide_idx)
                    Counts = np.array(Counts)
                    Indices = np.array(Indices)
                    return X, Counts, Indices


                X_imm, Counts_imm, _ = expand(imm)
                
                print "Min/Max/Median counts_imm", \
                    Counts_imm.min(), Counts_imm.max(), np.median(Counts_imm)
                W_imm = 1.0 / np.array(Counts_imm)
                X_combined.extend(X_imm)
                W_combined.extend(W_imm)
                Y_combined.extend([1] * len(X_imm))
                
                X_non, Counts_non, _ = expand(non)

                print "Min/Max/Median counts_non", \
                    Counts_non.min(), Counts_non.max(), np.median(Counts_non)

                W_non = 1.0 / np.array(Counts_non)
                X_combined.extend(X_non)
                W_combined.extend(W_non)
                Y_combined.extend([0] * len(X_non))

                def strings_to_array(strings):
                    all_strings = ''.join(strings)
                    X = np.fromstring(all_strings, dtype='uint8')
                    m = len(X) / kmer_length
                    X = X.reshape((m, kmer_length))
                    X -= ord('0')
                    return X

                X = strings_to_array(X_combined)
                Y = np.array(Y_combined)
                W = np.array(W_combined)
                print "# imm = %d, # non = %d" % (len(imm), len(non))
                print "Data shape", X.shape, "n_true", np.sum(Y)
                
                rf = BalancedEnsembleClassifier(n_estimators = 200)
                aucs = sklearn.cross_validation.cross_val_score(
	          rf, X, Y, cv = 5, scoring='roc_auc')
		print "CV AUC %0.4f (std %0.4f)" % (np.mean(aucs), np.std(aucs))
                d['cv_auc'].append(np.mean(aucs))
                #rf = RandomForestClassifier(n_estimators = 100)
                rf.fit(X, Y, W)
                def predict(peptides):
                    Y_pred = np.zeros(len(peptides), dtype=float)
                    counts = np.zeros(len(peptides), dtype=int)
                    X_test, _, Indices = expand(peptides)
                    X_test = strings_to_array(X_test)
                    #Y_pred_raw = rf.predict(X_test)

                    Y_pred_prob = rf.predict_proba(X_test)[:, 1]
                    Y_pred_rescaled = (2 * (Y_pred_prob - 0.5))
                    Y_pred_weight = np.sign(Y_pred_rescaled) * Y_pred_rescaled ** 2
                    # group outputs by the sample they came from, 
                    # at the end we'll have the majority vote 
                    #Y_pred = rf.predict(X_test)
                    for (y,i) in zip(Y_pred_weight, Indices):
                        Y_pred[i] += y
                        counts[i] += 1
                    Y_pred /= counts 
                    return Y_pred >= 0
                Y_pos_pred = predict(pos_peptides)
                pos_acc = np.mean(Y_pos_pred)
                print "Tumor antigen accuracy %0.4f" % (pos_acc,)
                d['pos_acc'].append(pos_acc)

                Y_neg_pred = predict(neg_peptides)
                neg_acc = 1.0 - np.mean(Y_neg_pred)
                print "Non-immunogenic accuracy %0.4f" % (neg_acc,)
                d['neg_acc'].append(neg_acc)

                tp = np.sum(Y_pos_pred)
                fp = np.sum(Y_neg_pred)
                fn = np.sum(~Y_pos_pred)
                tn = np.sum(~Y_neg_pred)

                precision = tp / float(tp + fp)
                recall = tp / float(tp + fn)
                print "tp = %d, fp = %d, fn = %d, tn = %d" % \
                    (tp, fp, fn, tn)
                f1_score = 2 * (precision * recall) / (precision + recall)
                print "F-1 score: %0.4f" % f1_score
                d['f1_score'].append(f1_score)

                f_half_score = 1.25 * \
                    (precision * recall) / ((0.25 * precision) + recall)

                print "F-0.5 score: %0.4f" % f_half_score
                d['f_score'].append(f_half_score)
                acc = np.sqrt(pos_acc * neg_acc)
                print "sqrt(pos_acc * neg_acc) =", acc
                d['acc'].append(acc)
                print "---"
                print

df = pd.DataFrame(d)
df = df.sort('acc', ascending=False)
print df.to_string()
with open('May13_sliding_kmers.csv', 'w') as f:
    df.to_csv(f)
with open('May13_sliding_kmers.html', 'w') as f:
    df.to_html(f)
