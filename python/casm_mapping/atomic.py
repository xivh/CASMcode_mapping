import numpy as np

class AtomicMapping(object):
    """
    Stores atomic mapping solutions, r1[:,l] + disp[:,l] = F @ r2[:,perm[l]] + trans

    r1: (2d array) reference atomic positions, with the Cartesian coordinates of site l stored as
        column l
    r2: (2d array) unmapped atomic positions, with the Cartesian coordinates of site l stored as
        column l
    F: (2d array), lattice mapping solution deformation

    Attributes:

        disp: (2d array), array with Cartesian atomic displacement of site l stored as column l
        perm: (1d array), integer array, the assignment solution order permutation
        trans: (1d array) array describing a rigid translation of all sites

    """
    def __init__(self, disp, perm, trans):
        self.disp = np.array(disp)
        self.perm = np.array(perm)
        self.trans = np.array(trans)
