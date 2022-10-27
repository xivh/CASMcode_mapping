#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/mapping/AtomMapping.hh"
#include "casm/mapping/LatticeMapping.hh"
#include "casm/mapping/StructureMapping.hh"
#include "casm/mapping/map_atoms.hh"
#include "casm/mapping/map_lattices.hh"
#include "casm/mapping/map_structures.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;
using namespace CASM::mapping;

}  // namespace CASMpy

PYBIND11_MODULE(_methods, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Lattice, atom, and structure mapping methods

        libcasm.mapping.methods
        -----------------------

        The libcasm.mapping.methods module contains lattice, atom, and
        structure mapping methods.

    )pbdoc";
  py::module::import("libcasm.xtal");
  py::module::import("libcasm.mapping.info");

  m.def("map_lattices", &map_lattices, R"pbdoc(
      Find mappings between two lattices

      This method finds mappings from a superlattice of a reference "parent"
      lattice to a "child" lattice. The lattice mappings have the form:

      .. math::

          F L_1 T N = L_2

      where:

      - :math:`L_1` is a shape=(3,3) matrix with columns containing the
        reference "parent" lattice vectors.
      - :math:`L_2` is a shape=(3,3) matrix with columns containing the
        "child" lattice vectors.
      - :math:`F` is the parent-to-child deformation gradient tensor,
        a shape=(3,3) matrix. :math:`F` can be decomposed as:

        .. math::

            F = Q U = V Q

        into a pure rigid transformation component, :math:`Q`, (called
        an isometry) and either a right stretch tensor, :math:`U`, or
        a left stretch tensor, :math:`V`. :math:`U` and :math:`V` are
        both positive definite and symmetric and include all of the
        non-rigid lattice deformation. Isometries, :math:`Q`, have the
        property that :math:`Q^{-1} = Q^{\mathsf{T}}`, where :math:`^{\mathsf{T}}`
        is used to indicate matrix transpose.
      - :math:`T` is an integer transformation matrix that generates a
        superlattice of :math:`L_1`.
      - :math:`N` is a unimodular reorientation matrix that generates a
        lattice equivalent to :math:`L_1 T` with reoriented lattice
        vectors.

      Lattice mapping costs are calculated and ranked according to one of:

      - "isotropic_strain_cost": a strain cost, calculated to be
        volume-normalized and invariant to which structure is the
        parent/child. Calculated as:

        .. math::

            \mathrm{cost} &= \frac{1}{2}*\left( \frac{1}{3}*\mathrm{tr}\left(\tilde{B}_{c \to p}^{2} \right) + \frac{1}{3}*\mathrm{tr}\left(\tilde{B}_{p \to c}^{2} \right) \right)

            \tilde{F} &= \frac{1}{\det{F}^{1/3}} F = \tilde{Q} \tilde{U}

            \tilde{B} &= \tilde{U} - I,

        where :math:`F` is the deformation gradient tensor, which can be defined
        in either the parent-to-child (:math:`p \to c`) sense as
        :math:`F_{p \to c} L_1 T N = L_2`, or in the child-to-parent sense
        (:math:`c \to p`) according to :math:`F_{c \to p} = F_{p \to c}^{-1}`,
        :math:`B` is the Biot strain tensor, and the use of (:math:`\tilde{X}`)
        indicates that the quantity has been normalized to be volume
        invariant.

      - "symmetry_breaking_strain_cost": a strain cost, including only
        the components of the strain that break the symmetry of the parent
        point group. Calculated as:

        .. math::

            \mathrm{cost} &= \frac{1}{2}*\left( \frac{1}{3}*\mathrm{tr}\left( B_{sym-break, p \to c}^2 \right) + \frac{1}{3}*\mathrm{tr}\left( B_{sym-break, c \to p}^2 \right) \right)

            B_{sym-break} &= B - B_{sym}

            B_{sym} &= \frac{1}{N_{G_1}} * \sum_i ( G_1(i) * B * G_1^{\mathsf{T}}(i) )

        where :math:`B_{sym-break}`, is the symmetry-breaking Biot strain,
        :math:`G_1(i)` are parent point group operations, and
        :math:`N_{G_1}` is the total number of operations. Similar relations hold
        for calculating :math:`B_{sym-break, p \to c}` and
        :math:`B_{sym-break, c \to p}` from :math:`B`.

      For more details, see :cite:t:`THOMAS2021a`.

      Notes
      -----
      Deformations can be validly defined as parent-to-child or
      child-to-parent. Be careful as to which convention is being used.
      Lattice mappings use the parent-to-child definition.


      Parameters
      ----------
      lattice1 : libcasm.xtal.Lattice
          The reference "parent" lattice, :math:`L_1`.
      lattice2 : libcasm.xtal.Lattice
          The "child" lattice, :math:`L_2`.
      transformation_matrix_to_super : Optional[array_like, shape=(3,3)], optional
          An approximately integer transformation matrix that generates a
          superlattice of :math:`L_1`. The default value is the identity
          matrix.
      lattice1_point_group : List[libcasm.xtal.SymOp], optional
          Used to skip reorientation matrices that result in symmetrically
          equivalent mappings. The default (empty), is equivalent to only
          including the identity operation.
      lattice2_point_group : List[libcasm.xtal.SymOp], optional
          Used to skip reorientation matrices that result in symmetrically
          equivalent mappings. The default (empty), is equivalent to just
          including the identity operation.
      min_cost : float, default=0.
          Keep lattice mappings with cost >= min_cost
      max_cost : float, default=1e20
          Keep results with cost <= max_cost
      cost_method : str, default="isotropic_strain_cost"
          Selects the method used to calculate lattice mapping costs. One of
          "isotropic_strain_cost" or "symmetry_breaking_strain_cost".
      k_best : Optional[int], default=None
          If not None, then only keep the k-best results (i.e. k lattice mappings
          with minimum cost) satisfying the min_cost and max_cost constraints.
          If there are approximate ties, those will also be kept.
      reorientation_range : int, default=1
          The absolute value of the maximum element in the lattice mapping
          reorientation matrix, :math:`N`. This determines how many
          equivalent lattice vector reorientations are checked. Increasing
          the value results in more checks. The value 1 is expected to be
          sufficient because reduced cell lattices are compared internally.
      cost_tol : float, default=1e-5
          Tolerance for checking if lattice mapping costs are approximately
          equal.

      Returns
      -------
      lattice_mappings : List[Tuple[float, libcasm.mapping.LatticeMapping]]
          A list of tuple of lattice mapping cost (float) and
          `libcasm.mapping.LatticeMapping`, giving possible lattice
          mappings, sorted by lattice mapping cost.
      )pbdoc",
        py::arg("lattice1"), py::arg("lattice2"),
        py::arg("transformation_matrix_to_super") = std::nullopt,
        py::arg("reorientation_range") = 1,
        py::arg("lattice1_point_group") = std::vector<xtal::SymOp>{},
        py::arg("lattice2_point_group") = std::vector<xtal::SymOp>{},
        py::arg("min_cost") = 0.0, py::arg("max_cost") = 1e20,
        py::arg("cost_method") = std::string("isotropic_strain_cost"),
        py::arg("k_best") = std::nullopt, py::arg("cost_tol") = 1e-5);

  m.def("map_structures", &map_structures, R"pbdoc(
      Find mappings between two structures

      This method finds mappings from a superstructure of a reference "parent"
      structure to a "child" structure. It works by finding lattice mappings
      (:class:`~cast.mapping.LatticeMapping`) for symmetrically unique
      superlattices of the "parent" lattice for a range of supercell volumes,
      and for each potential lattice mapping finding atom mappings
      (:class:`~cast.mapping.AtomMapping`).

      The total structure mapping cost, total_cost, is a weighted mixture of
      the lattice mapping cost, lattice_cost, and the atom mapping cost,
      atom_cost:

      .. code-block:: Python

          total_cost = lattice_cost_weight*lattice_cost + (1.0 - lattice_cost_weight)*atom_cost

      where lattice_cost_weight is an input parameter.

      For more details, see :cite:t:`THOMAS2021a`.

      Parameters
      ----------
      prim : libcasm.xtal.Prim
          The reference "parent" structure, with lattice :math:`L_1`,
          represented as Prim, with occupation DoF indicating which atom
          types are allowed to map to each basis site.
      structure : libcasm.xtal.Lattice
          The "child" structure, with lattice :math:`L_2`.
      max_vol : int
          The maximum parent superstructure volume to consider, as a
          multiple of the parent structure volume.
      prim_factor_group : List[libcasm.xtal.SymOp], optional
          Used to skip symmetrically equivalent mappings. The default
          (empty), is equivalent to only including the identity operation.
      structure_factor_group : List[libcasm.xtal.SymOp], optional
          Used to skip symmetrically equivalent mappings. The default
          (empty), is equivalent to just including the identity operation.
      min_vol : int, default=1
          The minimum parent superstructure volume to consider, as a
          multiple of the parent structure volume.
      min_cost : float, default=0.
          Keep structure mappings with cost >= min_cost
      max_cost : float, default=1e20
          Keep structure mappings with cost <= max_cost
      lattice_cost_weight : float, default=0.5
          The fraction of the total cost due to the lattice strain cost.
          The remaining fraction (1.-lattice_cost_weight) is due to the
          atom cost.
      lattice_cost_method : str, default="isotropic_strain_cost"
          Selects the method used to calculate lattice mapping costs. One of
          "isotropic_strain_cost" or "symmetry_breaking_strain_cost".
      atom_cost_method : str, default="isotropic_disp_cost"
          Selects the method used to calculate atom mapping costs. One of
          "isotropic_disp_cost" or "symmetry_breaking_disp_cost".
      k_best : int, default=1
          Only keep the k-best results (i.e. k mappings with minimum cost)
          satisfying the min_cost and max_cost constraints. If there are
          approximate ties, those will also be kept.
      cost_tol : float, default=1e-5
          Tolerance for checking if structure mappings costs are approximately
          equal.

      Returns
      -------
      structure_mappings : List[Tuple[libcasm.mapping.StructureMappingCost, libcasm.mapping.StructureMapping]]
          A list of tuple of :class:`~libcasm.mapping.StructureMappingCost`
          and `libcasm.mapping.StructureMapping`, giving possible structure
          mappings, sorted by total cost.
      )pbdoc",
        py::arg("prim"), py::arg("structure"), py::arg("max_vol"),
        py::arg("prim_factor_group") = std::vector<xtal::SymOp>{},
        py::arg("structure_factor_group") = std::vector<xtal::SymOp>{},
        py::arg("min_vol") = 1, py::arg("min_cost") = 0.0,
        py::arg("max_cost") = 1e20, py::arg("lattice_cost_weight") = 0.5,
        py::arg("lattice_cost_method") = std::string("isotropic_strain_cost"),
        py::arg("atom_cost_method") = std::string("isotropic_disp_cost"),
        py::arg("k_best") = 1, py::arg("cost_tol") = 1e-5);

  m.def("map_atoms", &map_atoms, R"pbdoc(
      Find atom mappings between two structures, given a particular lattice mapping

      This method finds atom mappings from a superstructure of a reference
      "parent" structure to a "child" structure. It works by checking
      atom mappings (:class:`~cast.mapping.AtomMapping`) given a particular
      lattice mapping (:class:`~cast.mapping.LatticeMapping`).

      Atom mapping costs are calculated and ranked according to one of:

      - "isotropic_disp_cost": an atom mapping cost, calculated to be
        volume-normalized and invariant to which structure is the
        parent/child. Calculated as:

        .. math::

            \mathrm{cost} &= \frac{1}{2}*\left( \frac{1}{3}*\mathrm{tr}\left(\tilde{B}_{c \to p}^{2} \right) + \frac{1}{3}*\mathrm{tr}\left(\tilde{B}_{p \to c}^{2} \right) \right)

            \tilde{F} &= \frac{1}{\det{F}^{1/3}} F = \tilde{Q} \tilde{U}

            \tilde{B} &= \tilde{U} - I,

        where :math:`F` is the deformation gradient tensor, which can be defined
        in either the parent-to-child (:math:`p \to c`) sense as
        :math:`F_{p \to c} L_1 T N = L_2`, or in the child-to-parent sense
        (:math:`c \to p`) according to :math:`F_{c \to p} = F_{p \to c}^{-1}`,
        :math:`B` is the Biot strain tensor, and the use of (:math:`\tilde{X}`)
        indicates that the quantity has been normalized to be volume
        invariant.

      - "symmetry_breaking_strain_cost": a strain cost, including only
        the components of the strain that break the symmetry of the parent
        point group. Calculated as:

        .. math::

            \mathrm{cost} &= \frac{1}{2}*\left( \frac{1}{3}*\mathrm{tr}\left( B_{sym-break, p \to c}^2 \right) + \frac{1}{3}*\mathrm{tr}\left( B_{sym-break, c \to p}^2 \right) \right)

            B_{sym-break} &= B - B_{sym}

            B_{sym} &= \frac{1}{N_{G_1}} * \sum_i ( G_1(i) * B * G_1^{\mathsf{T}}(i) )

        where :math:`B_{sym-break}`, is the symmetry-breaking Biot strain,
        :math:`G_1(i)` are parent point group operations, and
        :math:`N_{G_1}` is the total number of operations. Similar relations hold
        for calculating :math:`B_{sym-break, p \to c}` and
        :math:`B_{sym-break, c \to p}` from :math:`B`.

      For more details, see :cite:t:`THOMAS2021a`.

      Parameters
      ----------
      prim : libcasm.xtal.Prim
          The reference "parent" structure, with lattice :math:`L_1`,
          represented as Prim, with occupation DoF indicating which atom
          types are allowed to map to each basis site.
      structure : libcasm.xtal.Lattice
          The "child" structure, with lattice :math:`L_2`.
      lattice_mapping : libcasm.mapping.LatticeMapping
          Defines the lattice mapping from the "parent" structure to
          the "child" structure.
      prim_factor_group : List[libcasm.xtal.SymOp], optional
          Used to skip symmetrically equivalent mappings. The default
          (empty), is equivalent to only including the identity operation.
      min_cost : float, default=0.
          Keep atom mappings with cost >= min_cost
      max_cost : float, default=1e20
          Keep atom mappings with cost <= max_cost
      atom_cost_method : str, default="isotropic_disp_cost"
          Selects the method used to calculate atom mapping costs. One of
          "isotropic_disp_cost" or "symmetry_breaking_disp_cost".
      k_best : int, default=1
          Only keep the k-best results (i.e. k mappings with minimum cost)
          satisfying the min_cost and max_cost constraints. If there are
          approximate ties, those will also be kept.
      cost_tol : float, default=1e-5
          Tolerance for checking if atom mapping costs are approximately
          equal.

      Returns
      -------
      atom_mappings : List[Tuple[float, libcasm.mapping.AtomMapping]]
          A list of tuple of atom_cost and `libcasm.mapping.AtomMapping`,
          giving possible atom mappings, sorted by atom mapping cost.
      )pbdoc",
        py::arg("prim"), py::arg("structure"), py::arg("lattice_mapping"),
        py::arg("prim_factor_group") = std::vector<xtal::SymOp>{},
        py::arg("min_cost") = 0.0, py::arg("max_cost") = 1e20,
        py::arg("atom_cost_method") = std::string("isotropic_disp_cost"),
        py::arg("k_best") = 1, py::arg("cost_tol") = 1e-5);

  // Apply structure mapping
  m.def("make_mapped_structure", &mapping::make_mapped_structure,
        py::arg("structure_mapping"), py::arg("unmapped_structure"),
        "Returns a structure equivalent to `unmapped_structure`, but "
        "translated and rotated into alignment with the reference `prim`. "
        "Strain, measured relative to the reference `prim`, is included as a "
        "global property `Ustrain`, and displacements, measured relative to "
        "the ideal site coordinates are included as the atom property `disp`. "
        "All other properties of the unmapped structure are also transformed "
        "into alignment.");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
