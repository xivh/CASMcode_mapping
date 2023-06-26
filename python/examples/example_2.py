import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims
from scipy.spatial.transform import Rotation

prim = xtal_prims.cubic(a=1.0, occ_dof=["A"])
prim_factor_group = xtal.make_factor_group(prim)
print("prim:", prim.to_dict())

Qi = Rotation.from_euler("z", 30, degrees=True).as_matrix()
structure_lattice = xtal.Lattice(Qi @ prim.lattice().column_vector_matrix())

structure = xtal.Structure(
    lattice=structure_lattice,
    atom_coordinate_frac=prim.coordinate_frac(),
    atom_type=["A"],
)
print("structure:", structure.to_dict())

structure_mappings = mapmethods.map_structures(
    prim,
    structure,
    prim_factor_group=prim_factor_group,
    max_vol=1,
    max_cost=0.0,
    min_cost=0.0,
)
print(len(structure_mappings))
for i, m in enumerate(structure_mappings):
    print("mapping:", i, "total_cost:", m.total_cost())