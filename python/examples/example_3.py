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

print(len(structure_mappings))
for i, m in enumerate(structure_mappings):
    print("mapping:", i, "total_cost:", m.total_cost())
    data = m.to_dict()
    for key in data:
        print(key, data[key])
    print()
