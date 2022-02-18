#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

/// \brief Lattice mapping transformation
///
/// The mapping transformation has the form:
///
/// \f[
///     L_1 T N = V Q L_2.
/// \f]
/// Variables:
/// - \f$L_1\f$: a reference parent lattice, represented as a 3x3
///   matrix whose columns are the lattice vectors.
/// - \f$(L_1 * T)\f$: a superlattice of the reference parent lattice, to be
///   mapped to, represented as a 3x3 matrix whose columns are the lattice
///   vectors.
/// - \f$L_2\f$: the child lattice to be mapped, represented
///   as a 3x3 matrix whose columns are the lattice vectors
/// - \f$N\f$: a unimodular matrix (integer valued, with \f$\det{N}=1\f$,
///   though represented with a floating point matrix for multiplication
///   purposes) that generates lattice vectors \f$(L_1 T N)\f$ of
///   lattices that are equivalent to a superlattice of the reference parent
///   lattice \f$(L_1 T)\f$
///
struct LatticeMapping {
  Eigen::Matrix3d Q;
  Eigen::Matrix3d V;
  Eigen::Matrix3d N;
  Eigen::Matrix3d T;
};

Eigen::Matrix3d get_latticemapping_Q(LatticeMapping const &lattice_mapping) {
  return lattice_mapping.Q;
}

Eigen::Matrix3d get_latticemapping_V(LatticeMapping const &lattice_mapping) {
  return lattice_mapping.V;
}

Eigen::Matrix3d get_latticemapping_N(LatticeMapping const &lattice_mapping) {
  return lattice_mapping.N;
}

Eigen::Matrix3d get_latticemapping_T(LatticeMapping const &lattice_mapping) {
  return lattice_mapping.T;
}

/// \brief Deformation gradient, child to parent, F = V * Q
Eigen::Matrix3d get_latticemapping_F(LatticeMapping const &lattice_mapping) {
  return lattice_mapping.V * lattice_mapping.Q;
}

/// \brief Atom mapping transformation
///
/// The assignment portion of the structure mapping algorithm finds
/// solutions \f$(p_i, \vec{t}, \vec{d}(i))\f$ of:
/// \f[
///     \vec{r_1}(i) + \vec{d}(i) = V * Q * \vec{r_2}(p_i) + \vec{t}
/// \f]
///
/// where:
/// - \f$\vec{r_1}(i)\f$: Vector of coordinates of atoms in the parent
///   superstructure. The value \f$\vec{r_1}(i)\f$ represents the Cartesian
///   coordinate of the \f$i\f$-th atom in the parent superstructure. The parent
///   superstructure is not returned directly as part of the mapping results,
///   but it can be constructed using:
///
///       xtal::SimpleStructure parent_superstructure = make_superstructure(
///           T * N, parent_structure);
///
///   where \f$T\f$, and\f$N\f$ come from the lattice mapping solution. Then
///   the \f$i\f$-th atom coordinate, \f$\vec{r_1}(i)\f$, is equal to:
///
///       parent_superstructure.atom_info.coords.col(i)
///
/// - \f$\vec{r_2}(i)\f$: Vector of coordinates of atoms in the unmapped child
///   structure. The value \f$\vec{r_2}(i)\f$ represents the Cartesian
///   coordinate of the \f$i\f$-th atom in the unmapped child structure.
/// - \f$V * Q\f$: Lattice transformation, from the unmapped child superlattice
///   to the parent superlattice, as determined by a solution to the lattice
///   mapping problem.
/// - \f$p_i\f$: A permutation vector, describes which atom in the unmapped
///   child structure (\f$p_i\f$) is mapped to the i-th site of the mapped
///   structure. Values of \f$p_i\f$ greater than the number of atoms in the
///   unmapped structure indicate inferred vacancies.
/// - \f$\vec{t}\f$: A translation vector, in Cartesian coordinates, of the de-
///   rotated and undeformed (mapped) child superstructure that minimizes the
///   atomic displacement cost.
/// - \f$\vec{d}(i)\f$: The displacement associated with the atom at the i-th
///   site in parent superstructure.
///
/// Additionally, structures with magnetic spin may have time reversal symmetry
/// which may relate the child structure to the lattice structure.
struct AtomMapping {
  Eigen::MatrixXd displacement;
  std::vector<Index> permutation;
  Eigen::Vector3d translation;
  bool time_reversal;
};

Eigen::MatrixXd get_atom_mapping_displacement(AtomMapping const &atom_mapping) {
  return atom_mapping.displacement;
}

std::vector<Index> get_atom_mapping_permutation(
    AtomMapping const &atom_mapping) {
  return atom_mapping.permutation;
}

Eigen::Vector3d get_atom_mapping_translation(AtomMapping const &atom_mapping) {
  return atom_mapping.translation;
}

bool get_atom_mapping_time_reversal(AtomMapping const &atom_mapping) {
  return atom_mapping.time_reversal;
}

/// \brief Structure mapping transformation
///
/// Combines lattice mapping and atom mapping to form a complete structure
/// mapping transformation
struct StructureMapping {
  LatticeMapping lattice_mapping;
  AtomMapping atom_mapping;
};

LatticeMapping get_lattice_mapping(StructureMapping const &structure_mapping) {
  return structure_mapping.lattice_mapping;
}

AtomMapping get_atom_mapping(StructureMapping const &structure_mapping) {
  return structure_mapping.atom_mapping;
}

struct StructureMappingScore {
  double lattice_cost;
  double atom_cost;
  double total_cost;
};

}  // namespace CASMpy

PYBIND11_MODULE(mapping, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        casm.mapping
        ---------

        The casm.mapping module is a Python interface to the mapping
        classes and methods in the CASM::mapping namespace of the CASM C++ libraries.
        This includes:

        - Methods for finding the mapping transformations that relate one structure to another

    )pbdoc";

  m.attr("TOL") = TOL;

  py::class_<LatticeMapping>(m, "LatticeMapping", R"pbdoc(
      A lattice mapping transformation
      )pbdoc")
      .def(py::init<Eigen::Matrix3d const &, Eigen::Matrix3d const &>(),
           "Construct a LatticeMapping transformation", py::arg("isometry"),
           py::arg("left_stretch"), py::arg("transformation_matrix_to_super"),
           py::arg("reorientation"), R"pbdoc(
      Construct a lattice mapping transformation

      The lattice mapping transformation has the form:

      .. math::

          L_1 T N = V Q L_2,

      where :math:`L_1` are the reference "parent" lattice vectors, and
      :math:`L_2` are the "child" lattice vectors being mapped to the
      parent vectors, as columns of shape=(3,3) matrices. The other
      shape=(3,3) matrices are:

      - :math:`V`, a symmetric stretch tensor, representing
        the symmetric strain
      - :math:`Q` matrix, a rigid transformation of :math:`L_2`
      - :math:`T`, an integer transformation matrix that generates a
        superlattice of :math:`L_1`
      - :math:`N`, a unimodular reorientation matrix that generates a
        lattice equivalent to :math:`L_1 T` with reoriented lattice
        vectors


      Parameters
      ----------

      left_stretch : array_like, shape=(3,3)
          The :math:`V` matrix, a symmetric stretch tensor, representing
          the symmetric strain
      isometry : array_like, shape=(3,3), optional
          The :math:`Q` matrix, representing a rigid transformation of
          :math:`L_2`. The default value is np.eye(3).
      transformation_matrix_to_super : array_like, shape=(3,3), dtype=int, optional
          The transformation matrix, :math:`T`, that generates a
          superlattice of the parent lattice, :math:`L_1`. The default
          value is np.eye(3).astype(int).
      reorientation : array_like, shape=(3,3), dtype=int, optional
          The unimodular matrix, :math:`N`, that generates a lattice
          equivalent to the parent superlattice, :math:`L_1 T`. The
          default value is np.eye(3).astype(int).
      )pbdoc")
      .def("isometry", &get_latticemapping_Q,
           "Return the shape=(3,3) isometry matrix, :math:`Q`.")
      .def("left_stretch", &get_latticemapping_V,
           "Return the shape=(3,3) left symmetric stretch tensor, :math:`V`.")
      .def("transformation_matrix_to_super", &get_latticemapping_T,
           "Return the shape=(3,3) parent supercell transformation matrix, "
           ":math:`T`.")
      .def("reorienation", &get_latticemapping_N,
           "Return the shape=(3,3) unimodular matrix, :math:`N`.")
      .def("deformation_gradient", &get_latticemapping_F,
           "Return the shape=(3,3) deformation gradient tensor, :math:`F = V "
           "Q`.");

  py::class_<AtomMapping>(m, "AtomMapping", R"pbdoc(
     An atom mapping transformation
     )pbdoc")
      .def(py::init<Eigen::Matrix3d const &, Eigen::Matrix3d const &>(),
           "Construct an AtomMapping transformation", py::arg("displacement"),
           py::arg("permutation"), py::arg("translation"),
           py::arg("time_reversal"), R"pbdoc(
     Construct an atom mapping transformation

     The atom mapping transformation has the form:

     .. math::

         \vec{r_1}(i) + \vec{d}(i) = V * Q * \vec{r_2}(p_i) + \vec{t}


      where:
         / - \f$\vec{r_1}(i)\f$: Vector of coordinates of atoms in the parent
         /   superstructure. The value \f$\vec{r_1}(i)\f$ represents the Cartesian
         /   coordinate of the \f$i\f$-th atom in the parent superstructure. The parent
         /   superstructure is not returned directly as part of the mapping results,
         /   but it can be constructed using:
         /
         /       xtal::SimpleStructure parent_superstructure = make_superstructure(
         /           T * N, parent_structure);
         /
         /   where \f$T\f$, and\f$N\f$ come from the lattice mapping solution. Then
         /   the \f$i\f$-th atom coordinate, \f$\vec{r_1}(i)\f$, is equal to:
         /
         /       parent_superstructure.atom_info.coords.col(i)
         /
         / - \f$\vec{r_2}(i)\f$: Vector of coordinates of atoms in the unmapped child
         /   structure. The value \f$\vec{r_2}(i)\f$ represents the Cartesian
         /   coordinate of the \f$i\f$-th atom in the unmapped child structure.
         / - \f$V * Q\f$: Lattice transformation, from the unmapped child superlattice
         /   to the parent superlattice, as determined by a solution to the lattice
         /   mapping problem.
         / - \f$p_i\f$: A permutation vector, describes which atom in the unmapped
         /   child structure (\f$p_i\f$) is mapped to the i-th site of the mapped
         /   structure. Values of \f$p_i\f$ greater than the number of atoms in the
         /   unmapped structure indicate inferred vacancies.
         / - \f$\vec{t}\f$: A translation vector, in Cartesian coordinates, of the de-
         /   rotated and undeformed (mapped) child superstructure that minimizes the
         /   atomic displacement cost.
         / - \f$\vec{d}(i)\f$: The displacement associated with the atom at the i-th
         /   site in parent superstructure.
         /
         / Additionally, structures with magnetic spin may have time reversal symmetry
         / which may relate the child structure to the lattice structure.



     Parameters
     ----------

     left_stretch : array_like, shape=(3,3)
         The :math:`V` matrix, a symmetric stretch tensor, representing
         the symmetric strain
     isometry : array_like, shape=(3,3), optional
         The :math:`Q` matrix, representing a rigid transformation of
         :math:`L_2`. The default value is np.eye(3).
     transformation_matrix_to_super : array_like, shape=(3,3), dtype=int, optional
         The transformation matrix, :math:`T`, that generates a
         superlattice of the parent lattice, :math:`L_1`. The default
         value is np.eye(3).astype(int).
     reorientation : array_like, shape=(3,3), dtype=int, optional
         The unimodular matrix, :math:`N`, that generates a lattice
         equivalent to the parent superlattice, :math:`L_1 T`. The
         default value is np.eye(3).astype(int).
     )pbdoc")
      .def("isometry", &get_latticemapping_Q,
           "Return the shape=(3,3) isometry matrix, :math:`Q`.")
      .def("left_stretch", &get_latticemapping_V,
           "Return the shape=(3,3) left symmetric stretch tensor, :math:`V`.")
      .def("transformation_matrix_to_super", &get_latticemapping_T,
           "Return the shape=(3,3) parent supercell transformation matrix, "
           ":math:`T`.")
      .def("reorienation", &get_latticemapping_N,
           "Return the shape=(3,3) unimodular matrix, :math:`N`.")
      .def("deformation_gradient", &get_latticemapping_F,
           "Return the shape=(3,3) deformation gradient tensor, :math:`F = V "
           "Q`.");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
