import numpy as np

class Configuration(object):
    """Configuration

    Attributes:
        prim: (Prim) Primitive crystal structure and allowed DoF
        T: (array) Integer array giving the supercell lattice vectors in terms of the prim lattice
            vectors, $L^{super} = L^{prim} @ T$.
        occ: (array of int) Integer 1d array specifying the occupants on each site in
            the configuration. The species `prim.occ_dof[occ[l]]` is occupying the `l`-th site of
            the configuration.
        local_dof_values: (dict of 2d np.array) The DoF value on the `l`-th site is $x^{standard}$ =
            `prim.local_dof[b][dofname] @ local_dof_values[dofname][:,l]`, where `b` is the
            sublattice containing site `l`.
        global_dof_values: (dict of 1d np.array) The DoF value is $x^{standard}$ =
            `prim.global_dof[dofname] @ global_dof_values[dofname]`.

    """
    def __init__(self, prim, T, occ=None, local_dof_values={}, global_dof_values={}):
        self.prim = prim
        self.T = np.array(T, dtype=int)
        if occ is None:
            occ = np.zeros(self.N_sites, dtype=int)
        else:
            self.occ = np.array(occ, dtype=int)
        self.local_dof_values = local_dof_values
        self.global_dof_values = global_dof_values

    @property
    def N_sites(self):
        return np.linalg.det(T) * self.prim.B_frac.shape[1]

    @staticmethod
    def from_dict(prim, config):
        """Read from Configuration JSON formatted dictionary"""

        T = np.array(config["transformation_matrix_to_supercell"])

        configdof = config["dof"]

        occ = copy.deepcopy(configdof["occ"])

        local_dof_values = {}
        if "local_dofs" in configdof:
            for key, obj in configdof["local_dofs"].items():
                local_dof_values[key] = np.array(obj["values"]).transpose()

        global_dof_values = {}
        if "global_dofs" in configdof:
            for key, obj in configdof["global_dofs"].items():
                global_dof_values[key] = np.array(obj["values"]).transpose()

        return Configuration(prim, T, occ=occ, local_dof_values=local_dof_values,
            global_dof_values=global_dof_values)



    def to_dict(self):
        """Convert to Configuration JSON formatted dictionary"""

        config = dict()

        config["transformation_matrix_to_supercell"] = self.T.tolist()

        config["dof"] = dict()

        config["dof"]["occ"] = self.occ

        if len(self.local_dof_values):
            d = dict()
            for key, val in self.local_dof_values.items():
                d[key] = {"values": val.transpose().tolist()}
            config["dof"]["local_dofs"] = d

        if len(self.global_dof_values):
            d = dict()
            for key, val in self.global_dof_values.items():
                d[key] = {"values": val.transpose().tolist()}
            config["dof"]["global_dofs"] = d

        return config
