import typing

import libcasm.mapping.info as mapinfo
import libcasm.xtal as xtal
import numpy as np


def map_lattices_without_reorientation(
    lattice1: xtal.Lattice,
    lattice2: xtal.Lattice,
    transformation_matrix_to_super: typing.Optional[np.ndarray] = None,
) -> mapinfo.LatticeMapping:
    """Map lattices without reorienting the child.

    This function may be used to find the lattice mapping from
    an ideal structure to a deformed structure during geometric
    relaxation via DFT. The lattice mapping has the form

    .. math::

        FL_1T = L_2

    where :math:`F` is the deformation tensor, :math:`L_1` is
    the ideal parent lattice, :math:`L_2` is the deformed child
    lattice, and :math:`T` is an optional transformation matrix
    from the parent to the child which may be provided if :math:`L2`
    is a superlattice of :math:`L1`.

    Parameters
    ----------
    lattice1: xtal.Lattice
        The parent crystal lattice.
    lattice2: xtal.Lattice
        The child crystal lattice.
    transformation_matrix_to_super: Optional[np.ndarray] = None
        An integer shape=(3,3) transformation matrix that generates a superlattice of `lattice1`.
        The transformation matrix that generates a superlattice of `lattice1`.
        If None, :math:`T` is set to the identity matrix.

    Returns
    -------
    lattice_mapping: mapinfo.LatticeMapping
        The lattice mapping from the parent to the child.
    """
    if transformation_matrix_to_super is None:
        transformation_matrix_to_super = np.eye(3, dtype=int)

    # calculate deformation gradient
    l1 = lattice1.column_vector_matrix()
    l2 = lattice2.column_vector_matrix()
    deformation_gradient = (
        l2 @ np.linalg.inv(transformation_matrix_to_super) @ np.linalg.inv(l1)
    )
    lattice_mapping = mapinfo.LatticeMapping(
        deformation_gradient=deformation_gradient,
        transformation_matrix_to_super=transformation_matrix_to_super,
        reorientation=np.eye(3, dtype=float),
    )
    return lattice_mapping
