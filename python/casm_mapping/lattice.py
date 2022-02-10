import numpy as np

class LatticeMapping(object):
    """
    Stores lattice mapping solutions, L1 @ T @ N = F @ L2

    L1: (2d array) reference lattice vectors, with lattice vectors stored as columns
    L2: (2d array) unmapped lattice vectors, with lattice vectors stored as columns

    Attributes:

        T: (2d array), integer transformation matrix to supercell lattice vectors
        N: (2d array), integer reorientation matrix
        F: (2d array), deformation matrix

    """
    def __init__(self, T, N, F):
        self.T = np.array(T)
        self.N = np.array(N)
        self.F = np.array(F)
