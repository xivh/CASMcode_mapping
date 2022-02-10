import numpy as np

class DoFSet(object):
    """
        `x^{standard} = B @ x^{prim}`
    """
    def __init__(self, axis_names, B):
        self.axis_names = axis_names
        self.B = B

    @staticmethod
    def from_dict(val):
        return DoFSet(val["axis_names"], np.array(val["basis"]).transpose())

    def to_dict(self):
        return {
            "axis_names": self.axis_names,
            "basis": self.B.transpose().tolist()
        }

class Prim(object):
    """Primitive crystal structure and allowed DoF (`prim`)

    Attributes:
        L: (array) lattice vectors, stored as columns
        B_frac: (array) basis site fractional coordinates, stored as columns
        B_cart: (array) basis site Cartesian coordinates, stored as columns
        occ_dof: (list of list) occupation DoF, occ_dof[b] = [species1, species2, ...], where
            `species1`, `species2`, etc. are strings with the names of the species allowed on
            sublattice `b`.
        local_dof: (list of dict of DoFSet) continuous local DoFSet, `local_dof[b][dofname]` is the
            DoFSet for `dofname` on sublattice `b`.
        global_dof: (dict of DoFSet) continuous global DoF, `global_dof[dofname]` is the DoFSet for
            `dofname`.
        species: (dict) Dictionary describing non-atomic species. See `Prim` JSON documentation.
        title: (str) A title for the project. Must consist of alphanumeric characters and
            underscores only. The first character may not be a number.
        description: (str or None) Project description.

    """
    def __init__(self, title, L, B_frac, occ_dof, local_dof=None, global_dof={}, species={},
                 description=None):
        self.title = title
        self.L = L
        self.B_frac = B_frac
        self.occ_dof = occ_dof
        if local_dof is None:
            local_dof = [{} for b in self.B_frac.shape[1]]
        self.local_dof = local_dof
        self.global_dof = global_dof
        self.species = species
        self.description = description

    @property
    def B_cart(self):
        return self.L @ self.B_frac

    @property
    def B_size(self):
        return self.B_frac.shape[1]

    @staticmethod
    def from_dict(prim):
        """Read from Prim JSON formatted dictionary"""

        title = copy.deepcopy(prim["title"])

        L = np.array(prim["lattice_vectors"]).transpose()

        N_sites = len(prim["basis"])
        B = np.zeros((3, N_sites))
        for l in range(N_sites):
            B[:,l] = np.array(prim["basis"][l]["coordinate"])
        coordinate_mode = prim["coordinate_mode"]
        c = coordinate_mode[0].lower()
        if c == "c":
            # If B = B_cart = L @ B_frac, convert to B_frac
            B = np.linalg.pinv(L) @ B

        occ_dof = []
        for l in range(N_sites):
            occ_dof.append(prim["basis"][l]["occupants"])

        local_dof = []
        for l in range(N_sites):
            dofsets = {}
            if "dofs" in prim["basis"][l]:
                dofs = prim["basis"][l]["dofs"]
                for key, val in dofs.items():
                    dofsets[key] = DoFSet.from_dict(val)
            local_dof.append(dofsets)

        global_dof = {}
        if "dofs" in prim:
            dofs = prim["dofs"]
            for key, val in dofs.items():
                global_dof[key] = DoFSet.from_dict(val)

        return Prim(title, L, B, occ_dof, local_dof, global_dof, species=prim.get("species",{}),
            description=prim.get("description", None))


    def to_dict(self):
        """Convert to Prim JSON formatted dictionary"""

        prim = dict()
        prim["title"] = self.title
        prim["lattice_vectors"] = self.L.transpose().tolist()
        prim["coordinate_mode"] = "Fractional"
        prim["basis"] = []
        for l in range(self.B_size):
            site = dict()
            site["coordinate"] = self.B_frac[:,l].tolist()
            site["occupants"] = self.occ_dof[l]
            if len(self.local_dof[l]):
                site["dofs"] = {
                    dofname: dofset.to_dict() for dofname, dofset in self.local_dof[l].items()}
            prim["basis"].append(site)
        if len(self.global_dof):
            prim["dofs"] = {
                dofname: dofset.to_dict() for dofname, dofset in self.global_dof.items()}
        if len(self.species):
            prim["species"] = self.species
        if self.description is not None:
            prim["description"] = self.description
        return prim
