# Authors: G. Varoquaux, A. Gramfort, A. Passos, O. Grisel
# License: BSD
# Taken from https://github.com/scikit-learn/scikit-learn/blob/babe4a5d0637ca172d47e1dfdd2f6f3c3ecb28db/scikits/learn/utils/

import math
import numpy as np
from scipy import linalg


def check_random_state(seed):
    """Turn seed into a np.random.RandomState instance
    If seed is None, return the RandomState singleton used by np.random.
    If seed is an int, return a new RandomState instance seeded with seed.
    If seed is already a RandomState instance, return it.
    Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, int):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
' instance' % seed)


def safe_sparse_dot(a, b, dense_output=False):
    """Dot product that handle the sparse matrix case correctly"""
    from scipy import sparse
    if sparse.issparse(a) or sparse.issparse(b):
        ret = a * b
        if dense_output and hasattr(ret, "toarray"):
            ret = ret.toarray()
        return ret
    else:
		return np.dot(a,b)


def fast_svd(M, k, p=None, q=0, transpose='auto', random_state=0):
    """Computes the k-truncated randomized SVD
    Parameters
    ===========
    M: ndarray or sparse matrix
        Matrix to decompose
    k: int
        Number of singular values and vectors to extract.
    p: int (default is k)
        Additional number of samples of the range of M to ensure proper
        conditioning. See the notes below.
    q: int (default is 0)
        Number of power iterations (can be used to deal with very noisy
        problems).
    transpose: True, False or 'auto' (default)
        Whether the algorithm should be applied to M.T instead of M. The
        result should approximately be the same. The 'auto' mode will
        trigger the transposition if M.shape[1] > M.shape[0] since this
        implementation of randomized SVD tend to be a little faster in that
        case).
    random_state: RandomState or an int seed (0 by default)
        A random number generator instance to make behavior
    Notes
    =====
    This algorithm finds the exact truncated singular values decomposition
    using randomization to speed up the computations. It is particularly
    fast on large matrices on which you whish to extract only a small
    number of components.
    (k + p) should be strictly higher than the rank of M. This can be
    checked by ensuring that the lowest extracted singular value is on
    the order of the machine precision of floating points.
    References
    ==========
    Finding structure with randomness: Stochastic algorithms for constructing
    approximate matrix decompositions
    Halko, et al., 2009 (arXiv:909)
    A randomized algorithm for the decomposition of matrices
    Per-Gunnar Martinsson, Vladimir Rokhlin and Mark Tygert
    """
    if p == None:
        p = k

    random_state = check_random_state(random_state)
    n_samples, n_features = M.shape

    if transpose == 'auto' and n_samples > n_features:
        transpose = True
    if transpose:
        # this implementation is a bit faster with smaller shape[1]
        M = M.T

   # generating random gaussian vectors r with shape: (M.shape[1], k + p)
    r = random_state.normal(size=(M.shape[1], k + p))

    # sampling the range of M using by linear projection of r
    Y = safe_sparse_dot(M, r)
    del r

    # apply q power iterations on Y to make to further 'imprint' the top
    # singular values of M in Y
    for i in xrange(q):
        Y = safe_sparse_dot(M, safe_sparse_dot(M.T, Y))

    # extracting an orthonormal basis of the M range samples
    from .fixes import qr_economic
    Q, R = qr_economic(Y)
    del R

    # project M to the (k + p) dimensional space using the basis vectors
    B = safe_sparse_dot(Q.T, M)

    # compute the SVD on the thin matrix: (k + p) wide
    from scipy import linalg
    Uhat, s, V = linalg.svd(B, full_matrices=False)
    del B
    U = np.dot(Q, Uhat)

    if transpose:
        # transpose back the results according to the input convention
        return V[:k, :].T, s[:k], U[:, :k].T
    else:
return U[:, :k], s[:k], V[:k, :]