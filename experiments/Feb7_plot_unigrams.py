import scipy 
import scipy.sparse
import numpy as np 

from ..data import iedb, amino_acid, peptide_to_indices, reduced_alphabet


X, Y = iedb.load_dataset(human=True, max_ngram=3, reduced_alphabet = reduced_alphabet.gbmr4)
n_imm = np.sum(Y)

print "PCA"
import sklearn.decomposition

pca = sklearn.decomposition.RandomizedPCA(n_components=2, whiten=True)

X_2d = pca.fit_transform(X)
print "New dims", X_2d.shape

X_imm_2d = X_2d[:n_imm, :]
X_non_2d = X_2d[n_imm:, :]

import pylab

pylab.scatter(X_imm_2d[:, 0], X_imm_2d[:, 1], c='r')
pylab.scatter(X_non_2d[:, 0], X_non_2d[:, 1], c='g')
pylab.legend(('IMM', 'NON'), loc='best')
pylab.show()


print "MDS"
import sklearn.manifold
mds = sklearn.manifold.MDS(max_iter=100)

X_2d = mds.fit_transform(X)
print "New dims", X_2d.shape

X_imm_2d = X_2d[:n_imm, :]
X_non_2d = X_2d[n_imm:, :]

import pylab

pylab.scatter(X_imm_2d[:, 0], X_imm_2d[:, 1], c='r')
pylab.scatter(X_non_2d[:, 0], X_non_2d[:, 1], c='g')
pylab.legend(('IMM', 'NON'), loc='best')
pylab.show()