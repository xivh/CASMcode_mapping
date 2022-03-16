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

StructureMapping make_structure_mapping(xtal::BasicStructure const &prim,
                                        LatticeMapping const &lattice_mapping,
                                        AtomMapping const &atom_mapping) {
  return StructureMapping(std::make_shared<xtal::BasicStructure const>(prim),
                          lattice_mapping, atom_mapping);
}

}  // namespace CASMpy

PYBIND11_MODULE(mapping, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        casm.mapping
        ------------

        The casm.mapping module is a Python interface to the mapping
        classes and methods in the CASM::mapping namespace of the CASM C++ libraries.
        This includes:

        - Methods for finding the mapping transformations that relate one structure to another

    )pbdoc";
  py::module::import("casm.xtal");
  m.attr("TOL") = TOL;

  py::class_<LatticeMapping>(m, "LatticeMapping", R"pbdoc(
      A mapping between two lattices

      A lattice mapping has the form:

      .. math::

          F L_1 T N = L_2,

      where:

      - :math:`L_1` is a shape=(3,3) matrix with columns containing the
        reference "parent" lattice vectors
      - :math:`L_2` is a shape=(3,3) matrix with columns containing the
        "child" lattice vectors
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
        superlattice of :math:`L_1`
      - :math:`N` is a unimodular reorientation matrix that generates a
        lattice equivalent to :math:`L_1 T` with reoriented lattice
        vectors

      Notes
      -----
      Deformations can be validly defined as parent-to-child or
      child-to-parent. Be careful as to which convention is being used.
      )pbdoc")
      .def(py::init<Eigen::Matrix3d const &, Eigen::Matrix3d const &,
                    Eigen::Matrix3d const &>(),
           py::arg("deformation_gradient"),
           py::arg("transformation_matrix_to_super"), py::arg("reorientation"),
           R"pbdoc(
          Construct a lattice mapping

          Parameters
          ----------

          deformation_gradient : array_like, shape=(3,3)
              The parent-to-child deformation gradient tensor,:math:`F`, a shape=(3,3)
              matrix.
          transformation_matrix_to_super : array_like, shape=(3,3), dtype=int, optional
              The transformation matrix, :math:`T`, that generates a
              superlattice of the parent lattice, :math:`L_1`. The default
              value is np.eye(3).astype(int).
          reorientation : array_like, shape=(3,3), dtype=int, optional
              The unimodular matrix, :math:`N`, that generates a lattice
              equivalent to the parent superlattice, :math:`L_1 T`, with
              different lattice vectors. The default value is np.eye(3).astype(int).
          )pbdoc")
      .def(
          "deformation_gradient",
          [](LatticeMapping const &m) { return m.deformation_gradient; },
          "Return the shape=(3,3) parent-to-child deformation gradient tensor.")
      .def(
          "transformation_matrix_to_super",
          [](LatticeMapping const &m) {
            return m.transformation_matrix_to_super;
          },
          "Return the shape=(3,3) parent supercell transformation matrix, "
          ":math:`T`.")
      .def(
          "reorientation",
          [](LatticeMapping const &m) { return m.reorientation; },
          "Return the shape=(3,3) unimodular matrix, :math:`N`.")
      .def(
          "isometry", [](LatticeMapping const &m) { return m.isometry; },
          "Return the shape=(3,3) isometry matrix, :math:`Q`, of the "
          "parent-to-child deformation gradient tensor.")
      .def(
          "left_stretch",
          [](LatticeMapping const &m) { return m.left_stretch; },
          "Return the shape=(3,3) left symmetric stretch tensor, :math:`V`, "
          "of the parent-to-child deformation gradient tensor.")
      .def(
          "right_stretch",
          [](LatticeMapping const &m) { return m.right_stretch; },
          "Return the shape=(3,3) right symmetric stretch tensor, :math:`U`, "
          "of the parent-to-child deformation gradient tensor.");

  py::class_<AtomMapping>(m, "AtomMapping", R"pbdoc(
     A mapping of atoms between two structures

     An atom mapping is defined in the context of an existing
     :class:`~casm.mapping.LatticeMapping`. An atom mapping has the form:

     .. math::

         F \left(\vec{r_1}(i) + \vec{d}(i) \right) = \vec{r_2}(p_i) + \vec{t}

     where:

     - :math:`\vec{r_1}(i)` is the Cartesian coordinates of the i-th atom in
       the parent superstructure. The parent superstructure can be constructed
       using `casm.xtal.make_superstructure`:

       .. code-block:: Python

           import casm.xtal as xtal
           parent_superstructure = xtal.make_superstructure(
               T * N, parent_structure)

       where :math:`T`, and :math:`N` are from a
       :class:`~casm.mapping.LatticeMapping`. Then the i-th atom coordinate,
       :math:`\vec{r_1}(i)`, is equal to:

       .. code-block:: Python

           parent_superstructure.atom_coordinate_cart()[:,i]

     - :math:`\vec{r_2}(i)` is the Cartesian coordinates of i-th atom in
       the child structure.
     - :math:`F` is the parent-to-child deformation gradient tensor from a
       :class:`~casm.mapping.LatticeMapping`.
     - :math:`p_i` is a permutation vector, specifying which atom in the
       child structure (:math:`p_i`) is mapped to the i-th site of the parent
       superstructure. Values of :math:`p_i` greater than the number of atoms
       in the child structure indicate inferred vacancies in mappings to
       :class:`~casm.xtal.Prim` with vacancies allowed.
     - :math:`\vec{t}` is a translation vector, in Cartesian coordinates, usually
       chosen so the average displacement is zero.
     - :math:`\vec{d}(i)`: The Cartesian displacement associated with the atom
        at the i-th site in the parent superstructure.

     Notes
     -----
     Displacements can be validly defined with different references. This
     method returns displacements that are consistent with the convention
     used by CASM for displacement DoFs: displacements are added to the
     ideal basis coordinates, then strain is applied.

     )pbdoc")
      .def(py::init<Eigen::MatrixXd const &, std::vector<Index> const &,
                    Eigen::Vector3d const &>(),
           py::arg("displacement"), py::arg("permutation"),
           py::arg("translation"), R"pbdoc(
          Construct an atom mapping

          Parameters
          ----------

          displacement : array_like, shape=(3,n)
             The shape=(3,n) matrix whose columns are the atom displacements
             :math:`\vec{d}(i)`.
          permutation : List[int], size=n
             The permutation vector, :math:`p_i`.
          translation : array_like, shape=(3,1),
             The translation vector, :math:`\vec{t}`.
          )pbdoc")
      .def(
          "displacement", [](AtomMapping const &m) { return m.displacement; },
          R"pbdoc(Return the shape=(3,n) matrix whose columns are the "
           "Cartesian atom displacements :math:`\vec{d}(i)`.)pbdoc")
      .def(
          "permutation", [](AtomMapping const &m) { return m.permutation; },
          R"pbdoc(Return the permutation vector, :math:`p_i`.)pbdoc")
      .def(
          "translation", [](AtomMapping const &m) { return m.translation; },
          R"pbdoc(Return the translation vector, :math:`\vec{t}`.)pbdoc");

  py::class_<StructureMapping>(m, "StructureMapping", R"pbdoc(
    A mapping between two structures

    A structure mapping is a combination of:

    - A :class:`~casm.xtal.Prim`
    - A :class:`~casm.mapping.LatticeMapping`
    - An :class:`~casm.mapping.AtomMapping`

    See those class descriptions for details of the mappings.

    )pbdoc")
      .def(py::init(&make_structure_mapping), py::arg("prim"),
           py::arg("lattice_mapping"), py::arg("atom_mapping"), R"pbdoc(
          Construct a structure mapping

          Parameters
          ----------

          prim : casm.xtal.Prim
              A :class:`~casm.xtal.Prim`
          lattice_mapping : casm.mapping.LatticeMapping
              A :class:`~casm.mapping.LatticeMapping`
          atom_mapping : casm.mapping.AtomMapping
              An :class:`~casm.mapping.AtomMapping`
          )pbdoc")
      .def(
          "prim", [](StructureMapping const &m) { return *m.shared_prim; },
          "Return the :class:`~casm.xtal.Prim`.")
      .def(
          "lattice_mapping",
          [](StructureMapping const &m) { return m.lattice_mapping; },
          "Return the :class:`~casm.mapping.LatticeMapping`.")
      .def(
          "atom_mapping",
          [](StructureMapping const &m) { return m.atom_mapping; },
          "Return the :class:`~casm.mapping.AtomMapping`.");

  m.def(
      "has_same_prim",
      [](StructureMapping const &first, StructureMapping const &second) {
        return first.shared_prim == second.shared_prim;
      },
      R"pbdoc(
      Check if StructureMapping share the same prim

      Notes
      -----
      Currently this check is not valid for StructureMapping manually constructed
      in Python, only StructureMapping output from a mapping method.


      Parameters
      ----------
      first : casm.mapping.StructureMapping
          The first StructureMapping
      second : casm.mapping.StructureMapping
          The second StructureMapping

      Returns
      -------
      bool : True if first and second share the same prim
      )pbdoc",
      py::arg("first"), py::arg("second"));

  py::class_<StructureMappingCost>(m, "StructureMappingCost", R"pbdoc(
   Holds mapping cost between two structures

   The structure mapping cost is a combination of:

   - A lattice mapping cost, lattice_cost
   - An atom mapping cost, atom_cost

   The total structure mapping cost, total_cost, is lattice_cost_weight*lattice_cost + (1.0 - lattice_cost_weight)*atom_cost, where lattice_cost_weight is an input to the :func:`~casm.mapping.map_structures` method.

   See those :class:`~casm.mapping.LatticeMapping` and :class:`~casm.mapping.AtomMapping` for details of the lattice and atom mappings.

   )pbdoc")
      .def(py::init<double, double, double>(), py::arg("lattice_cost"),
           py::arg("atom_cost"), py::arg("total_cost"), R"pbdoc(
         Construct a structure mapping

         Parameters
         ----------

         lattice_cost : float
             The lattice mapping cost
         atom_cost : float
             The atom mapping cost
         total_cost : float
             The total mapping cost
         )pbdoc")
      .def(
          "lattice_cost",
          [](StructureMappingCost const &s) { return s.lattice_cost; },
          "Return the lattice mapping cost.")
      .def(
          "atom_cost",
          [](StructureMappingCost const &s) { return s.atom_cost; },
          "Return the atom mapping cost.")
      .def(
          "total_cost",
          [](StructureMappingCost const &s) { return s.total_cost; },
          "Return the total mapping cost.");

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
      lattice1 : casm.xtal.Lattice
          The reference "parent" lattice, :math:`L_1`.
      lattice2 : casm.xtal.Lattice
          The "child" lattice, :math:`L_2`.
      transformation_matrix_to_super : Optional[array_like, shape=(3,3)], optional
          An approximately integer transformation matrix that generates a
          superlattice of :math:`L_1`. The default value is the identity
          matrix.
      lattice1_point_group : List[casm.xtal.SymOp], optional
          Used to skip reorientation matrices that result in symmetrically
          equivalent mappings. The default (empty), is equivalent to only
          including the identity operation.
      lattice2_point_group : List[casm.xtal.SymOp], optional
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
      lattice_mappings : List[Tuple[float, casm.mapping.LatticeMapping]]
          A list of tuple of lattice mapping cost (float) and
          `casm.mapping.LatticeMapping`, giving possible lattice
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
      prim : casm.xtal.Prim
          The reference "parent" structure, with lattice :math:`L_1`,
          represented as Prim, with occupation DoF indicating which atom
          types are allowed to map to each basis site.
      structure : casm.xtal.Lattice
          The "child" structure, with lattice :math:`L_2`.
      max_vol : int
          The maximum parent superstructure volume to consider, as a
          multiple of the parent structure volume.
      prim_factor_group : List[casm.xtal.SymOp], optional
          Used to skip symmetrically equivalent mappings. The default
          (empty), is equivalent to only including the identity operation.
      structure_factor_group : List[casm.xtal.SymOp], optional
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
      strain_cost_method : str, default="isotropic_strain_cost"
          Selects the method used to calculate lattice mapping costs. One of
          "isotropic_strain_cost" or "symmetry_breaking_strain_cost".
      atom_cost_method : str, default="isotropic_strain_cost"
          Selects the method used to calculate atom mapping costs. One of
          "isotropic_atom_cost" or "symmetry_breaking_atom_cost".
      k_best : int, default=1
          Only keep the k-best results (i.e. k mappings with minimum cost)
          satisfying the min_cost and max_cost constraints. If there are
          approximate ties, those will also be kept.
      cost_tol : float, default=1e-5
          Tolerance for checking if structure mappings costs are approximately
          equal.

      Returns
      -------
      structure_mappings : List[Tuple[casm.mapping.StructureMappingCost, casm.mapping.StructureMapping]]
          A list of tuple of :class:`~casm.mapping.StructureMappingCost`
          and `casm.mapping.StructureMapping`, giving possible structure
          mappings, sorted by total cost.
      )pbdoc",
        py::arg("prim"), py::arg("structure"), py::arg("max_vol"),
        py::arg("prim_factor_group") = std::vector<xtal::SymOp>{},
        py::arg("structure_factor_group") = std::vector<xtal::SymOp>{},
        py::arg("min_vol") = 1, py::arg("min_cost") = 0.0,
        py::arg("max_cost") = 1e20, py::arg("lattice_cost_weight") = 0.5,
        py::arg("strain_cost_method") = std::string("isotropic_strain_cost"),
        py::arg("atom_cost_method") = std::string("isotropic_atom_cost"),
        py::arg("k_best") = 1, py::arg("cost_tol") = 1e-5);

  m.def("map_atoms", &map_atoms, R"pbdoc(
      Find atom mappings between two structures, given a particular lattice mapping

      This method finds atom mappings from a superstructure of a reference
      "parent" structure to a "child" structure. It works by checking
      atom mappings (:class:`~cast.mapping.AtomMapping`) given a particular
      lattice mapping (:class:`~cast.mapping.LatticeMapping`).

      Atom mapping costs are calculated and ranked according to one of:

      - "isotropic_atom_cost": an atom mapping cost, calculated to be
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
      prim : casm.xtal.Prim
          The reference "parent" structure, with lattice :math:`L_1`,
          represented as Prim, with occupation DoF indicating which atom
          types are allowed to map to each basis site.
      structure : casm.xtal.Lattice
          The "child" structure, with lattice :math:`L_2`.
      lattice_mapping : casm.mapping.LatticeMapping
          Defines the lattice mapping from the "parent" structure to
          the "child" structure.
      prim_factor_group : List[casm.xtal.SymOp], optional
          Used to skip symmetrically equivalent mappings. The default
          (empty), is equivalent to only including the identity operation.
      min_cost : float, default=0.
          Keep atom mappings with cost >= min_cost
      max_cost : float, default=1e20
          Keep atom mappings with cost <= max_cost
      atom_cost_method : str, default="isotropic_strain_cost"
          Selects the method used to calculate atom mapping costs. One of
          "isotropic_atom_cost" or "symmetry_breaking_atom_cost".
      k_best : int, default=1
          Only keep the k-best results (i.e. k mappings with minimum cost)
          satisfying the min_cost and max_cost constraints. If there are
          approximate ties, those will also be kept.
      cost_tol : float, default=1e-5
          Tolerance for checking if atom mapping costs are approximately
          equal.

      Returns
      -------
      atom_mappings : List[Tuple[float, casm.mapping.AtomMapping]]
          A list of tuple of atom_cost and `casm.mapping.AtomMapping`,
          giving possible atom mappings, sorted by atom mapping cost.
      )pbdoc",
        py::arg("prim"), py::arg("structure"), py::arg("lattice_mapping"),
        py::arg("prim_factor_group") = std::vector<xtal::SymOp>{},
        py::arg("min_cost") = 0.0, py::arg("max_cost") = 1e20,
        py::arg("atom_cost_method") = std::string("isotropic_atom_cost"),
        py::arg("k_best") = 1, py::arg("cost_tol") = 1e-5);

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
