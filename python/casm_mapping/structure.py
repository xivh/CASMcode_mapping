import numpy as np

class Structure(object):
    """
    A crystal structure

    Attributes:
        L: (2d array) Lattice vectors, stored as columns
        atom_coords_frac: (2d array) Atom fractional coordinates, stored as columns
        atom_type: (list of str) Strings with the names of the atomic species on each site
        atom_properties: (dict of 2d array) The vector `atom_properties[propname][:,l]` is the
            value of local continuous property `propname` on the `l`-th site.
        global_properties: (dict of 1d array) The vector `global_properties[propname]` is the
            value of global continuous property `propname`.

    """
    def __init__(self, L, atom_coords_frac, atom_type, atom_properties={}, global_properties={}):
        self.L = L
        self.atom_coords_frac = atom_coords_frac
        self.atom_type = atom_type
        self.atom_properties = atom_properties
        self.global_properties = global_properties

    @property
    def atom_coords_cart(self):
        return self.L @ self.atom_coords_frac

    @staticmethod
    def from_dict(struct):
        """Read from Structure JSON formatted dictionary"""

        L = np.array(struct["lattice_vectors"]).transpose()

        N_sites = len(struct["atom_coords"])
        B = np.array(struct["atom_coords"]).transpose()
        coordinate_mode = struct["coordinate_mode"]
        c = coordinate_mode[0].lower()
        if c == "c":
            # If B = B_cart = L @ B_frac, convert to B_frac
            B = np.linalg.pinv(L) @ B

        atom_type = copy.copy(struct["atom_type"])

        atom_properties = {}
        if "atom_properties" in struct:
            for key, obj in struct["atom_properties"].items():
                atom_properties[key] = np.array(obj["value"]).transpose()

        global_properties = {}
        if "atom_properties" in struct:
            for key, obj in struct["atom_properties"].items():
                atom_properties[key] = np.array(obj["value"]).transpose()


    def to_dict(self):
        struct = dict()
        struct["lattice_vectors"] = self.L.transpose().tolist()
        struct["coordinate_mode"] = "Fractional"
        struct["atom_type"] = self.atom_type

        if len(self.atom_properties):
            d = dict()
            for key, val in self.local_dof_values.items():
                d[key] = {"value": val.transpose().tolist()}
            struct["atom_properties"] = d


class StructureMapping(obj):
    """
    A structure mapping solution, consisting of lattice and atomic mapping

        Lattice mapping: L1 @ T1 @ N = F @ L2
        Atomic mapping: r1[:,l] + disp[:,l] = F @ r2[:,perm[l]] + trans

    L1: (2d array) reference lattice vectors, with lattice vectors stored as columns
    L2: (2d array) unmapped lattice vectors, with lattice vectors stored as columns
    r1: (2d array) reference atomic positions, with the Cartesian coordinates of site l stored as
        column l
    r2: (2d array) unmapped atomic positions, with the Cartesian coordinates of site l stored as
        column l

    Attributes:

        lattice_mapping: (LatticeMapping) The lattice mapping portion of the solution
        atomic_mapping: (AtomicMapping) The atomic mapping portion of the solution

    """
    def __init__(self):
        self.lattice_mapping = None
        self.atomic_mapping = None

    def apply_to(self, structure):
        return None


def make_factor_group(structure):
    return None

def transform_structure(symop, structure):
    return None

def find_structure_mappings(structure1, structure2):
    return None

def is_equivalent(structure1, structure2):
    return None

def make_variants(structure1, structure2):
    return None
