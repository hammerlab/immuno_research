import numpy as np 
from ..data import amino_acid, load_imma2, transform_features

X,Y = load_imma2()

import sklearn
import sklearn.cross_validation 
import sklearn.ensemble

n_classifiers = 5000

clf = sklearn.ensemble.RandomForestClassifier(n_classifiers)

fns = [amino_acid.hydropathy, 
       amino_acid.volume, 
       amino_acid.pK_side_chain,
       amino_acid.polarity, 
       amino_acid.prct_exposed_residues,
       amino_acid.hydrophilicity, 
       amino_acid.accessible_surface_area,
       amino_acid.local_flexibility,
       amino_acid.accessible_surface_area_folded,
       amino_acid.refractivity
       ]

print "Pairwise ratios"
X2 = transform_features(X, fns, pairwise_ratios = True)
print X2.shape
print np.mean(sklearn.cross_validation.cross_val_score(clf, X2, Y, cv = 10))

print "Mean per feature"
X3 = transform_features(X, fns, mean = True)
print np.mean(sklearn.cross_validation.cross_val_score(clf, X3, Y, cv = 10))

print "Positions 4,6,8,9"
X4 = transform_features(X, fns, positions = (4,6,8,9))
print np.mean(sklearn.cross_validation.cross_val_score(clf, X4, Y, cv = 10))


