import numpy as np
import sklearn
import sklearn.cross_validation
import sklearn.ensemble

from epitopes import iedb, features, amino_acid

imm, non = iedb.load_tcell_classes(peptide_length = 9)
X, Y = features.make_kmer_dataset(imm, non)

n_classifiers = 50
clf = sklearn.ensemble.RandomForestClassifier(n_classifiers)

print "Amino acids"
print np.mean(sklearn.cross_validation.cross_val_score(clf, X, Y, cv = 10))

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

print "All features, all positions"
X2 = features.transform_rows(X, fns)
print np.mean(sklearn.cross_validation.cross_val_score(clf, X2, Y, cv = 10))

print "Pairwise ratios"
X3 = features.transform_rows(X, fns, pairwise_ratios = True)
print X3.shape
print np.mean(sklearn.cross_validation.cross_val_score(clf, X3, Y, cv = 10))

print "Mean per feature"
X3 = features.transform_rows(X, fns, mean = True)
print np.mean(sklearn.cross_validation.cross_val_score(clf, X3, Y, cv = 10))

print "Positions 4,6,8,9"
X4 = features.transform_rows(X, fns, positions = (4,6,8,9))
print np.mean(sklearn.cross_validation.cross_val_score(clf, X4, Y, cv = 10))


