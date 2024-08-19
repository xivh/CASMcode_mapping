import numpy as np

import libcasm.xtal as xtal
from libcasm.mapping.info import (
    isotropic_strain_cost,
    symmetry_breaking_strain_cost,
)


def _deformation_gradient(lattice1, lattice2):
    # F * L1 = L2
    L1 = lattice1.column_vector_matrix()
    L2 = lattice2.column_vector_matrix()
    return L2 @ np.linalg.inv(L1)


def test_strain_cost_1():
    lattice1 = xtal.Lattice(np.eye(3))
    lattice2 = xtal.Lattice(np.eye(3))
    lattice1_point_group = xtal.make_point_group(lattice1)
    F = _deformation_gradient(lattice1, lattice2)

    iso_cost = isotropic_strain_cost(F)
    symbreak_cost = symmetry_breaking_strain_cost(F, lattice1_point_group)
    print(iso_cost, symbreak_cost)
    assert np.isclose(iso_cost, 0.0)
    assert np.isclose(symbreak_cost, 0.0)


def test_strain_cost_2():
    L1 = np.eye(3)
    L2 = L1 * 1.1

    lattice1 = xtal.Lattice(L1)
    lattice2 = xtal.Lattice(L2)

    lattice1_point_group = xtal.make_point_group(lattice1)
    F = _deformation_gradient(lattice1, lattice2)

    iso_cost = isotropic_strain_cost(F)
    symbreak_cost = symmetry_breaking_strain_cost(F, lattice1_point_group)
    print(iso_cost, symbreak_cost)
    assert np.isclose(iso_cost, 0.0)
    assert np.isclose(symbreak_cost, 0.0)


def test_strain_cost_2b():
    L1 = np.eye(3)
    L2 = np.array(
        [
            [1.0 / 1.1, 0.0, 0.0],  # a
            [0.0, 1.0, 0.0],  # a
            [0.0, 0.0, 1.1],  # c
        ]
    ).transpose()

    lattice1 = xtal.Lattice(L1)
    lattice2 = xtal.Lattice(L2)

    lattice1_point_group = xtal.make_point_group(lattice1)
    F = _deformation_gradient(lattice1, lattice2)

    iso_cost = isotropic_strain_cost(F)
    symbreak_cost = symmetry_breaking_strain_cost(F, lattice1_point_group)
    print(iso_cost, symbreak_cost)
    assert not np.isclose(iso_cost, 0.0)
    assert not np.isclose(symbreak_cost, 0.0)


def test_strain_cost_3():
    L1 = np.eye(3)
    L2 = np.array(
        [
            [1.0, 0.0, 0.0],  # a
            [0.0, 1.0, 0.0],  # a
            [0.0, 0.0, 1.1],  # c
        ]
    ).transpose()

    lattice1 = xtal.Lattice(L1)
    lattice2 = xtal.Lattice(L2)

    lattice1_point_group = xtal.make_point_group(lattice1)
    F = _deformation_gradient(lattice1, lattice2)

    iso_cost = isotropic_strain_cost(F)
    symbreak_cost = symmetry_breaking_strain_cost(F, lattice1_point_group)
    print(iso_cost, symbreak_cost)
    assert not np.isclose(iso_cost, 0.0)
    assert not np.isclose(symbreak_cost, 0.0)


def test_strain_cost_4():
    L1 = np.array(
        [
            [1.0, 0.0, 0.0],  # a
            [0.0, 1.0, 0.0],  # a
            [0.0, 0.0, 1.5],  # c
        ]
    ).transpose()
    L2 = np.array(
        [
            [0.9, 0.0, 0.0],  # a
            [0.0, 0.9, 0.0],  # a
            [0.0, 0.0, 1.6],  # c
        ]
    ).transpose()

    lattice1 = xtal.Lattice(L1)
    lattice2 = xtal.Lattice(L2)

    lattice1_point_group = xtal.make_point_group(lattice1)
    F = _deformation_gradient(lattice1, lattice2)

    iso_cost = isotropic_strain_cost(F)
    symbreak_cost = symmetry_breaking_strain_cost(F, lattice1_point_group)
    print(iso_cost, symbreak_cost)
    assert not np.isclose(iso_cost, 0.0)
    assert np.isclose(symbreak_cost, 0.0)


def test_strain_cost_5():
    L1 = np.array(
        [
            [1.0, 0.0, 0.0],  # a
            [0.0, 1.5, 0.0],  # b
            [0.0, 0.0, 1.8],  # c
        ]
    ).transpose()
    L2 = np.array(
        [
            [0.9, 0.0, 0.0],  # a
            [0.0, 1.4, 0.0],  # b
            [0.0, 0.0, 2.0],  # c
        ]
    ).transpose()

    lattice1 = xtal.Lattice(L1)
    lattice2 = xtal.Lattice(L2)

    lattice1_point_group = xtal.make_point_group(lattice1)
    F = _deformation_gradient(lattice1, lattice2)

    iso_cost = isotropic_strain_cost(F)
    symbreak_cost = symmetry_breaking_strain_cost(F, lattice1_point_group)
    print(iso_cost, symbreak_cost)
    assert not np.isclose(iso_cost, 0.0)
    assert np.isclose(symbreak_cost, 0.0)

    L3 = np.array(
        [
            [0.9, 0.0, 0.0],  # a
            [0.0, 1.4, 0.1],  # b
            [0.0, 0.0, 2.0],  # c
        ]
    ).transpose()
    lattice3 = xtal.Lattice(L3)
    F = _deformation_gradient(lattice1, lattice3)
    symbreak_cost = symmetry_breaking_strain_cost(F, lattice1_point_group)
    print(iso_cost, symbreak_cost)
    assert not np.isclose(symbreak_cost, 0.0)
