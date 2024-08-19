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
