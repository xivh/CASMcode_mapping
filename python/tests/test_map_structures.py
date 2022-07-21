import math
import numpy as np
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal


def test_map_structures_0():
    """Map to self: max_cost==0.0 -> # of mappings == 48"""
    lattice = xtal.Lattice(
        np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]).transpose())
    coordinate_frac = np.array([
        [0., 0., 0.],
    ]).transpose()
    occ_dof = [
        ["A"],
    ]
    prim = xtal.Prim(lattice=lattice,
                     coordinate_frac=coordinate_frac,
                     occ_dof=occ_dof)

    structure = xtal.Structure(lattice=lattice,
                               atom_coordinate_frac=coordinate_frac,
                               atom_type=["A"])

    structure_mappings = mapmethods.map_structures(prim,
                                                   structure,
                                                   max_vol=1,
                                                   max_cost=0.,
                                                   min_cost=0.)
    assert len(structure_mappings) == 48


def test_map_structures_1():
    """Map to self: Use prim factor group -> # of tied-for-lowest cost mappings == 1"""
    lattice = xtal.Lattice(
        np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]).transpose())
    coordinate_frac = np.array([
        [0., 0., 0.],
    ]).transpose()
    occ_dof = [
        ["A"],
    ]
    prim = xtal.Prim(lattice=lattice,
                     coordinate_frac=coordinate_frac,
                     occ_dof=occ_dof)
    prim_factor_group = xtal.make_factor_group(prim)

    structure = xtal.Structure(lattice=lattice,
                               atom_coordinate_frac=coordinate_frac,
                               atom_type=["A"])

    structure_mappings = mapmethods.map_structures(
        prim,
        structure,
        prim_factor_group=prim_factor_group,
        max_vol=1,
        max_cost=0.,
        min_cost=0.)
    assert len(structure_mappings) == 1


def test_map_structures_2():
    """Map to superstructure of self: 2 equally perfect permutations"""
    lattice = xtal.Lattice(
        np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]).transpose())
    coordinate_frac = np.array([
        [0., 0., 0.],
    ]).transpose()
    occ_dof = [
        ["A"],
    ]
    prim = xtal.Prim(lattice=lattice,
                     coordinate_frac=coordinate_frac,
                     occ_dof=occ_dof)
    prim_factor_group = xtal.make_factor_group(prim)

    structure_lattice = xtal.Lattice(
        np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 2.0],
        ]).transpose())
    atom_coordinate_frac = np.array([
        [0., 0., 0.],
        [0., 0., 0.5],
    ]).transpose()
    structure = xtal.Structure(lattice=structure_lattice,
                               atom_coordinate_frac=atom_coordinate_frac,
                               atom_type=["A", "A"])

    structure_mappings = mapmethods.map_structures(
        prim,
        structure,
        prim_factor_group=prim_factor_group,
        max_vol=2,
        max_cost=0.,
        min_cost=0.)
    assert len(structure_mappings) == 2


def test_map_structures_2():
    """Map to ordered structure to self: 1 perfect mapping"""
    prim_lattice = xtal.Lattice(
        np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 2.0],
        ]).transpose())
    coordinate_frac = np.array([
        [0., 0., 0.],
        [0., 0., 0.5],
    ]).transpose()
    occ_dof = [
        ["A"],
        ["B"],
    ]
    prim = xtal.Prim(lattice=prim_lattice,
                     coordinate_frac=coordinate_frac,
                     occ_dof=occ_dof)
    prim_factor_group = xtal.make_factor_group(prim)

    structure_lattice = xtal.Lattice(
        np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 2.0],
        ]).transpose())
    atom_coordinate_frac = np.array([
        [0., 0., 0.],
        [0., 0., 0.5],
    ]).transpose()
    structure = xtal.Structure(lattice=structure_lattice,
                               atom_coordinate_frac=atom_coordinate_frac,
                               atom_type=["A", "B"])

    structure_mappings = mapmethods.map_structures(
        prim,
        structure,
        prim_factor_group=prim_factor_group,
        max_vol=2,
        max_cost=0.,
        min_cost=0.)
    assert len(structure_mappings) == 1


def test_map_structures_3():
    """Map to Ezz && Eyz"""
    lattice = xtal.Lattice(
        np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]).transpose())
    coordinate_frac = np.array([
        [0., 0., 0.],
    ]).transpose()
    occ_dof = [
        ["A"],
    ]
    prim = xtal.Prim(lattice=lattice,
                     coordinate_frac=coordinate_frac,
                     occ_dof=occ_dof)
    prim_factor_group = xtal.make_factor_group(prim)

    structure_lattice = xtal.Lattice(
        np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.1],
            [0.0, 0.0, 1.1],
        ]).transpose())
    atom_coordinate_frac = np.array([
        [0., 0., 0.],
    ]).transpose()
    structure = xtal.Structure(lattice=structure_lattice,
                               atom_coordinate_frac=atom_coordinate_frac,
                               atom_type=["A"])

    structure_mappings = mapmethods.map_structures(
        prim,
        structure,
        prim_factor_group=prim_factor_group,
        max_vol=1,
        max_cost=0.1,
        min_cost=0.)
    assert len(structure_mappings) == 1
    strucscore, structure_mapping = structure_mappings[0]
    U = structure_mapping.lattice_mapping().right_stretch()
    assert np.allclose(
        U,
        np.array([
            [1., 0., 0.],
            [0., 1.00362465, 0.05232166],
            [0., 0.05232166, 1.09875495],
        ]))
