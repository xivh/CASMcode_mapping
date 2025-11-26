import typing

import numpy as np

import libcasm.mapping.info as mapinfo
import libcasm.xtal as xtal


def map_lattices_without_reorientation(
    lattice1: xtal.Lattice,
    lattice2: xtal.Lattice,
    transformation_matrix_to_super: typing.Optional[np.ndarray] = None,
) -> mapinfo.LatticeMapping:
    """Map lattices without reorienting the child.

    This function may be used to find the lattice mapping from an ideal lattice or
    superlattice to a deformed lattice. The lattice mapping without reorientation has
    the form

    .. math::

        F L_1 T = L_2,

    where:

    - :math:`L_1` is a shape=(3,3) matrix with columns containing the
      reference "parent" lattice vectors
    - :math:`L_2` is a shape=(3,3) matrix with columns containing the
      "child" lattice vectors
    - :math:`F` is the parent-to-child deformation gradient tensor,
      a shape=(3,3) matrix.
    - :math:`T` is an integer transformation matrix that generates a
      superlattice of :math:`L_1`.

    This is equivalent to a lattice mapping with reorientation matrix, :math:`N`, equal
    to the identity matrix, :math:`I` (see
    :class:`~libcasm.mapping.info.LatticeMapping`).


    Parameters
    ----------
    lattice1: libcasm.xtal.Lattice
        The parent lattice, :math:`L_1`.
    lattice2: libcasm.xtal.Lattice
        The child lattice, :math:`L_2`.
    transformation_matrix_to_super: Optional[np.ndarray] = None
        A shape=(3,3) integer transformation matrix, :math:`T` that generates a
        superlattice of `lattice1`. If None, :math:`T` is set to the identity matrix.

    Returns
    -------
    lattice_mapping: libcasm.mapping.info.LatticeMapping
        The lattice mapping from the parent to the child, with :math:`N = I`.
    """

    T = transformation_matrix_to_super
    if T is None:
        T = np.eye(3, dtype=int)

    # calculate deformation gradient
    L1 = lattice1.column_vector_matrix()
    L2 = lattice2.column_vector_matrix()

    # F @ (L1 @ T) = L_2
    F = L2 @ np.linalg.inv(L1 @ T)

    return mapinfo.LatticeMapping(
        deformation_gradient=F,
        transformation_matrix_to_super=T,
        reorientation=np.eye(3, dtype=float),
    )


def direct_structure_mapping(
    structure1: xtal.Structure,
    structure2: xtal.Structure,
    remove_mean_displacement: bool = True,
):
    """Map lattices and atoms without reorienting the lattice vectors or permuting atom
    indices.

    Parameters
    ----------
    structure1 : libcasm.xtal.Structure
        The reference structure (parent).
    structure2 : libcasm.xtal.Structure
        The target structure (child).
    remove_mean_displacement: bool = True
       Displacements are first calculated under periodic boundary conditions. If True,
       the mean displacement is removed and a corresponding translation added to the
       atom mapping. If False, displacements are included as calculated and a
       zero-valued translation is used.


    Returns
    -------
    lmap: libcasm.mapping.info.LatticeMapping
        The lattice mapping from parent to child.
    amap: libcasm.mapping.info.AtomMapping
        The atom mapping from parent to child.
    """
    # Lattice mapping:
    lmap = map_lattices_without_reorientation(
        lattice1=structure1.lattice(),
        lattice2=structure2.lattice(),
    )
    F = lmap.deformation_gradient()

    # The mapped structure is constructed as:
    # F @ (r1 + d) = r2 + 0
    # r1 + d = F_inv @ r2
    # d = F_inv @ r2 - r1
    #
    # If removing mean displacement:
    # d' = F_inv @ r2 - r1 - d_mean
    # =>
    # F @ (r1 + d') = r2 + translation
    # F @ r1 + r2 - F @ r1 - F @ d_mean = r2 + translation
    # translation = -F @ d_mean
    #

    # Atom mapping:
    r1 = structure1.atom_coordinate_cart()
    r2 = structure2.atom_coordinate_cart()
    F_inv = np.linalg.inv(F)
    r2_ref = F_inv @ r2
    d_direct = F_inv @ r2 - r1
    d_pbc = np.zeros_like(d_direct)
    for i in range(r1.shape[1]):
        d_pbc[:, i] = xtal.min_periodic_displacement(
            lattice=structure1.lattice(),
            r1=r1[:, i],
            r2=r2_ref[:, i],
        )
    if remove_mean_displacement:
        mean_d = np.mean(d_pbc, axis=1)  # shape (3,)
        d_pbc -= mean_d[:, np.newaxis]  # shape (3,1) for broadcasting
        translation = -F @ mean_d
    else:
        translation = np.zeros(3)
    amap = mapinfo.AtomMapping(
        displacement=d_pbc,
        permutation=[i for i in range(r1.shape[1])],
        translation=translation,
    )

    return (lmap, amap)
