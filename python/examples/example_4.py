import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims
import libcasm.xtal.structures as xtal_structures

prim = xtal_prims.BCC(r=1.0, occ_dof=["A", "B", "Va"])
prim_factor_group = xtal.make_factor_group(prim)

hcp_structure = xtal_structures.HCP(r=1.0, atom_type="A")

structure_mappings = mapmethods.map_structures(
    prim,
    hcp_structure,
    prim_factor_group=prim_factor_group,
    max_vol=4,
    k_best=1,
)

m = structure_mappings[0]

# f=0.0, end point == parent
structure = mapmethods.make_mapped_structure(m.interpolated(0.0), hcp_structure)
print("2-atom bcc:")
print(structure.to_poscar_str())

# f=1.0, end point == child
print("hcp:")
structure = mapmethods.make_mapped_structure(m.interpolated(1.0), hcp_structure)
print(structure.to_poscar_str())

# f=0.5, intermediate point
print("2-atom intermediate:")
structure = mapmethods.make_mapped_structure(m.interpolated(0.5), hcp_structure)
print(structure.to_poscar_str())
