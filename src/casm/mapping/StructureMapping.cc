#include "casm/mapping/StructureMapping.hh"

#include "casm/casm_io/Log.hh"
#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/mapping/LatticeMapping.hh"
#include "casm/misc/CASM_Eigen_math.hh"

// debug
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {
namespace mapping {

/// \brief Return the mapped structure, with implied vacancies,
///     strain, and atomic displacement
///
/// The "mapped structure" lattice and site coordinates are constructed
/// from the parent structure by applying the lattice and atom mappings
/// without isometry. Atom names and atom properties are
/// determined from the unmapped structure, permutation, and
/// inverse isometry. Global properties are determined
/// from the unmapped global properties and inverse isometry. Strain,
/// using the right stretch tensor as the strain metric, is stored as
/// the global property "Ustrain". Displacement is stored as
/// atom properties.
///
/// Notes:
/// - This method is only implemented for atomic structures, not
///   molecular structures
/// - Implicit vacancies are added as "Va"
/// - Throws if unmapped_structure already has a strain or disp
///   property
///
/// \returns mapped_structure, including vacancies
///
xtal::SimpleStructure make_mapped_structure(
    StructureMapping const &structure_mapping,
    xtal::SimpleStructure const &unmapped_structure) {
  if (unmapped_structure.mol_info.size()) {
    throw std::runtime_error(
        "Error: CASM::mapping::make_mapped_structure is only implemented for "
        "atomic structures");
  }
  if (unmapped_structure.atom_info.coords.cols() !=
      unmapped_structure.atom_info.names.size()) {
    throw std::runtime_error(
        "Error in CASM::mapping::make_mapped_structure: unmapped_structure has "
        "inconsistent atom coords and names");
  }
  for (auto const &value : unmapped_structure.properties) {
    if (value.first.find("strain") != std::string::npos) {
      throw std::runtime_error(
          "Error in CASM::mapping::make_mapped_structure: unmapped_structure "
          "already has a strain property");
    }
  }
  for (auto const &value : unmapped_structure.atom_info.properties) {
    if (AnisoValTraits::name_suffix(value.first) == "disp") {
      throw std::runtime_error(
          "Error in CASM::mapping::make_mapped_structure: unmapped_structure "
          "already has a disp property");
    }
  }

  xtal::BasicStructure const &prim = *structure_mapping.shared_prim;

  // make the parent structure from the prim for the ideal lattice
  // and site coordinates - occupants are set to "Va"
  xtal::SimpleStructure parent_structure;
  parent_structure.lat_column_mat = prim.lattice().lat_column_mat();
  parent_structure.atom_info.resize(prim.basis().size());

  xtal::SimpleStructure mapped_structure;

  // make mapped lattice
  auto const &lattice_mapping = structure_mapping.lattice_mapping;
  auto const &Q = lattice_mapping.isometry;
  auto const &U = lattice_mapping.right_stretch;
  auto const &L1 = parent_structure.lat_column_mat;
  auto const &T = lattice_mapping.transformation_matrix_to_super;
  auto const &N = lattice_mapping.reorientation;

  mapped_structure.lat_column_mat = U * L1 * T * N;

  // atomic coordinate relations:
  //
  ///    F \left(\vec{r_1}(i) + \vec{d}(i) \right) = \vec{r_2}(p_i) + \vec{t}
  //     U \left(\vec{r_1}(i) + \vec{d}(i) \right) = \vec{r_3}(i)
  //
  //     F: lattice_mapping.deformation_gradient
  //     U: lattice_mapping.right_stretch
  //     \vec{r_1}(i): parent_superstructure.atom_info.coords.col(i)
  //     \vec{d}(i): atom_mapping.displacement.col(i)
  //     \vec{r_2}(p_i): unmapped_structure.atom_info.coords.col(p_i)
  //     \vec{r_3}(i): mapped_struture.atom_info.coords.col(i)
  //     p_i: atom_mapping.permutation[i]
  //     \vec{t}: atom_mapping.translation
  //
  //     remember: due to mapping within periodic boundaries, coordinates
  //     comparisons must be made checking for equality up to a lattice
  //     translation.

  xtal::SimpleStructure parent_superstructure =
      make_superstructure(lround(T * N).cast<int>(), parent_structure);

  // Note, n_atoms includes implied vacancies
  auto const &atom_mapping = structure_mapping.atom_mapping;
  Index n_sites = parent_superstructure.atom_info.size();
  Index n_atoms = unmapped_structure.atom_info.size();
  mapped_structure.atom_info.resize(n_sites);
  mapped_structure.atom_info.coords =
      U * (parent_superstructure.atom_info.coords + atom_mapping.displacement);

  // atomic property relations (includes atom type/name):
  //
  //     q_3(i) = M(Q^{-1}) q_2(p_i)
  //
  //     q_2(p_i): Property at site p_i of the unmapped structure
  //     q_3(i): Property at site i of the mapped structure
  //     Q: lattice_mapping.isometry
  //     M(Q^{-1}): representation matrix for transforming property q
  //                according to inverse of lattice_mapping.isometry
  //

  // Q.transpose() == Q^{-1}
  xtal::SymOp Q_inv = xtal::SymOp::point_operation(Q.transpose());

  // Create transformed atom properties, with columns of zeros for implied
  // vacancies
  std::map<std::string, Eigen::MatrixXd> transformed_atom_properties;
  for (auto const &pair : unmapped_structure.atom_info.properties) {
    std::string key = pair.first;
    Eigen::MatrixXd const &q2 = pair.second;
    Eigen::MatrixXd q2_with_Va = Eigen::MatrixXd::Zero(q2.rows(), n_sites);
    q2_with_Va.block(0, q2.rows(), 0, q2.cols()) = q2;
    try {
      AnisoValTraits traits(key);
      Eigen::MatrixXd M = traits.symop_to_matrix(
          get_matrix(Q_inv), get_translation(Q_inv), get_time_reversal(Q_inv));
      transformed_atom_properties.emplace(key, M * q2);
    } catch (std::exception &e) {
      std::stringstream msg;
      msg << "CASM does not know how to transform the local property '" << key
          << "'. The property name suffix must be the name of a local property "
             "that CASM can transform. Local properties include 'coordinate', "
             "'disp', 'force', 'Cmagspin', 'Cunitmagspin', 'NCmagspin', "
             "'NCunitmagspin', 'SOmagspin', 'SOunitmagspin', "
             "'dorbitaloccupation', and 'dorbitaloccupationspinpolarized'.";
      CASM::err_log().paragraph(msg.str());
      throw std::runtime_error(
          std::string("Cannot transform local property '") + key + "'");
    }
  }

  // Set zeroed atom properties
  for (auto const &pair : transformed_atom_properties) {
    mapped_structure.atom_info.properties.emplace(pair.first,
                                                  0.0 * pair.second);
  }

  // Set mapped atom names and properties
  for (Index i = 0; i < n_sites; ++i) {
    Index p_i = atom_mapping.permutation[i];
    if (p_i >= n_atoms) {
      // implied vacancy - has name "Va"
      mapped_structure.atom_info.names[i] = "Va";
    } else {
      // atoms or explicit vacancies
      mapped_structure.atom_info.names[i] =
          unmapped_structure.atom_info.names[p_i];
    }

    for (auto const &pair : transformed_atom_properties) {
      mapped_structure.atom_info.properties.at(pair.first).col(i) =
          pair.second.col(p_i);
    }
  }

  // Add displacement
  mapped_structure.atom_info.properties.emplace("disp",
                                                atom_mapping.displacement);

  // global property relations:
  //
  //     q_3 = M(Q^{-1}) q_2
  //
  //     q_2: Global property of the unmapped structure
  //     q_3: Global property of the mapped structure
  //     Q: lattice_mapping.isometry
  //     M(Q^{-1}): representation matrix for transforming property q
  //                according to inverse of lattice_mapping.isometry
  //

  // Set mapped global properties
  for (auto const &pair : unmapped_structure.properties) {
    std::string key = pair.first;
    Eigen::MatrixXd const &q2 = pair.second;
    try {
      AnisoValTraits traits(key);
      Eigen::MatrixXd M = traits.symop_to_matrix(
          get_matrix(Q_inv), get_translation(Q_inv), get_time_reversal(Q_inv));
      mapped_structure.properties.emplace(key, M * q2);
    } catch (std::exception &e) {
      std::stringstream msg;
      msg << "CASM does not know how to transform the global property '" << key
          << "'. The property name suffix must be the name of a global "
             "property that CASM can transform. Global properties include "
             "'energy', 'Ustrain', 'Bstrain', 'GLstrain', 'AEstrain', "
             "'Hstrain', and 'cost'.";
      CASM::err_log().paragraph(msg.str());
      throw std::runtime_error(
          std::string("Cannot transform global property '") + key + "'");
    }
  }

  // Add Ustrain (vector in standard basis)
  Eigen::Matrix3d const &e = U;
  Eigen::VectorXd E_vector = Eigen::VectorXd::Zero(6);
  double w = std::sqrt(2.);
  E_vector << e(0, 0), e(1, 1), e(2, 2), w * e(1, 2), w * e(0, 2), w * e(0, 1);
  mapped_structure.properties.emplace("Ustrain", E_vector);

  return mapped_structure;
}

}  // namespace mapping
}  // namespace CASM
