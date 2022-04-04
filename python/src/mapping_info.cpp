#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casm/crystallography/BasicStructure.hh"
#include "casm/mapping/AtomMapping.hh"
#include "casm/mapping/LatticeMapping.hh"
#include "casm/mapping/StructureMapping.hh"

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

PYBIND11_MODULE(info, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Data structures for representing mapping transformations

        casm.mapping.info
        -----------------

        The casm.mapping.info module contains data structures representing mapping
        transformations.

    )pbdoc";
  py::module::import("casm.xtal");

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
          "Returns the shape=(3,3) parent-to-child deformation gradient "
          "tensor.")
      .def(
          "transformation_matrix_to_super",
          [](LatticeMapping const &m) {
            return m.transformation_matrix_to_super;
          },
          "Returns the shape=(3,3) parent supercell transformation matrix, "
          ":math:`T`.")
      .def(
          "reorientation",
          [](LatticeMapping const &m) { return m.reorientation; },
          "Returns the shape=(3,3) unimodular matrix, :math:`N`.")
      .def(
          "isometry", [](LatticeMapping const &m) { return m.isometry; },
          "Returns the shape=(3,3) isometry matrix, :math:`Q`, of the "
          "parent-to-child deformation gradient tensor.")
      .def(
          "left_stretch",
          [](LatticeMapping const &m) { return m.left_stretch; },
          "Returns the shape=(3,3) left symmetric stretch tensor, :math:`V`, "
          "of the parent-to-child deformation gradient tensor.")
      .def(
          "right_stretch",
          [](LatticeMapping const &m) { return m.right_stretch; },
          "Returns the shape=(3,3) right symmetric stretch tensor, :math:`U`, "
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
          R"pbdoc(Returns the shape=(3,n) matrix whose columns are the "
           "Cartesian atom displacements :math:`\vec{d}(i)`.)pbdoc")
      .def(
          "permutation", [](AtomMapping const &m) { return m.permutation; },
          R"pbdoc(Returns the permutation vector, :math:`p_i`.)pbdoc")
      .def(
          "translation", [](AtomMapping const &m) { return m.translation; },
          R"pbdoc(Returns the translation vector, :math:`\vec{t}`.)pbdoc");

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
          "Returns the :class:`~casm.xtal.Prim`.")
      .def(
          "lattice_mapping",
          [](StructureMapping const &m) { return m.lattice_mapping; },
          "Returns the :class:`~casm.mapping.LatticeMapping`.")
      .def(
          "atom_mapping",
          [](StructureMapping const &m) { return m.atom_mapping; },
          "Returns the :class:`~casm.mapping.AtomMapping`.");

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
          "Returns the lattice mapping cost.")
      .def(
          "atom_cost",
          [](StructureMappingCost const &s) { return s.atom_cost; },
          "Returns the atom mapping cost.")
      .def(
          "total_cost",
          [](StructureMappingCost const &s) { return s.total_cost; },
          "Returns the total mapping cost.");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
