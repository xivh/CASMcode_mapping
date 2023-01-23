#include "casm/mapping/AtomMapping.hh"

#include "casm/global/eigen.hh"

namespace CASM {
namespace mapping {

/// \struct AtomMapping
/// \brief Atom mapping transformation
///
/// The assignment portion of the structure mapping algorithm finds
/// solutions \f$(p_i, \vec{t}, \vec{d}(i))\f$ of:
/// \f[
///     F \left( \vec{r_1}(i) + \vec{d}(i) \right) = \vec{r_2}(p_i) + \vec{t}
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
///   structure. Values of \f$p_i\f$ greater than or equal to the number of
///   atoms in the unmapped structure indicate inferred vacancies.
/// - \f$\vec{t}\f$: A translation vector, in Cartesian coordinates, of the de-
///   rotated and undeformed (mapped) child superstructure that minimizes the
///   atomic displacement cost.
/// - \f$\vec{d}(i)\f$: The displacement associated with the atom at the i-th
///   site in parent superstructure.
///
/// Additionally, structures with magnetic spin may have time reversal symmetry
/// which may relate the child structure to the parent structure.
///
/// There are in general many potential choices of atom mapping, with different
/// permutations and displacements. An assignment algorithm such as the
/// Hungarian method, with a cost function that depends on the displacements,
/// \f$\vec{d}(i)\f$, is used to score mappings.
///

/// \brief Constructor
///
/// Note: See `AtomMapping` class documentation for the atom
/// mapping transformation definition.
///
/// \param _displacement A (3,N) matrix of displacements, in Cartesian
///     coordinates, where N is the number of sites in the parent
///     superstructure
/// \param _permutation A size N permutation vector
/// \param _translation A translation vector, in Cartesian coordiantes
AtomMapping::AtomMapping(Eigen::MatrixXd const &_displacement,
                         std::vector<Index> const &_permutation,
                         Eigen::Vector3d const &_translation)
    : displacement(_displacement),
      permutation(_permutation),
      translation(_translation) {}

/// \brief Return mappings that result in atom positions along the
///     transformation pathway from the parent to the aligned child
///     structure
///
/// \param atom_mapping Original atom mapping
/// \param interpolation_factor Interpolation factor. The value 0.0
///     corresponds to the child sites mapped to the ideal parent
///     sites; and the value 1.0 corresponds to the child sites,
///     transformed by isometry and translation to align
///     with the ideal parent structure.
AtomMapping interpolated_mapping(AtomMapping const &atom_mapping,
                                 double interpolation_factor) {
  double f = interpolation_factor;
  return AtomMapping(atom_mapping.displacement * f, atom_mapping.permutation,
                     atom_mapping.translation);
}

}  // namespace mapping
}  // namespace CASM
