#ifndef CASM_mapping_SearchData
#define CASM_mapping_SearchData

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/global/eigen.hh"
#include "casm/mapping/LatticeMapping.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
}

namespace mapping {

/// Mapping an atomic structure to another structure is a
/// heirarchical search problem:
/// - Try a prim
/// - Try lattice mappings to supercells of the prim
/// - Try translations of the lattice-mapped structures
/// - Try assignments of the translated, lattice-mapped structures
///
/// The data structures here hold information specific to
/// the various levels of this search:
/// - StructureSearchData: Data related to the structure
///   being mapped
/// - PrimSearchData: Data related to a prim the structure
///   might be mapped to
/// - LatticeMappingSearchData: Data related to a particular
///   lattice mapping of a structure to a superstructure prim
/// - TrialTranslationSearchData: Data related to a
///   particular trial translation of a lattice mapping

/// \brief Holds data on the structure being mapped
struct StructureSearchData {
  /// \brief Constructor
  StructureSearchData(xtal::Lattice const &_lattice,
                      Eigen::MatrixXd const &_atom_coordinate_cart,
                      std::vector<std::string> const &_atom_type,
                      std::vector<xtal::SymOp> _structure_factor_group =
                          std::vector<xtal::SymOp>{});

  /// \brief The structure's lattice
  xtal::Lattice lattice;

  /// \brief The number of atoms (includes explicitly included vacancies)
  Index N_atom;

  /// \brief Shape=(3,N_atom), with Cartesian coordinates of atoms as columns
  Eigen::MatrixXd atom_coordinate_cart;

  /// \brief Size=N_atom, with the name of the atom at each site
  std::vector<std::string> atom_type;

  /// \brief Symmetry operations that may be used to skip symmetrically
  ///     equivalent mappings
  std::vector<xtal::SymOp> structure_factor_group;
};

/// \brief Holds prim-related data used for mapping searches
struct PrimSearchData {
  /// \brief Constructor
  PrimSearchData(
      std::shared_ptr<xtal::BasicStructure const> const &_shared_prim,
      std::vector<xtal::SymOp> _prim_factor_group = std::vector<xtal::SymOp>{},
      bool make_prim_sym_invariant_displacement_modes = true);

  /// \brief The prim
  std::shared_ptr<xtal::BasicStructure const> shared_prim;

  /// \brief The prim lattice
  xtal::Lattice prim_lattice;

  /// \brief The number of prim sites
  Index N_prim_site;

  /// \brief Shape=(3,N_prim_site), with prim site Cartesian
  ///     coordinates as columns
  Eigen::MatrixXd prim_site_coordinate_cart;

  /// \brief Size=N_prim_site, with names of atoms allowed on
  ///     each prim site
  std::vector<std::vector<std::string>> prim_allowed_atom_types;

  /// \brief True if vacancies are allowed on any prim site
  bool vacancies_allowed;

  /// \brief Symmetry operations that may be used to skip symmetrically
  ///     equivalent mappings
  std::vector<xtal::SymOp> prim_factor_group;

  /// \brief Size=N_mode, with shape=(3,N_prim_site) matrices, giving the
  ///     symmetry invariant displacement modes, with columns
  ///     giving the displacement associated with each site
  ///     for a given mode.
  ///
  /// For example:
  ///
  ///     Eigen::Vector3d site_displacement =
  ///         sym_invariant_displacement_modes[mode_index].col(site_index)
  ///
  /// This member may be expensive so it is optional to construct, though
  /// it is required for calculating the symmetry_breaking_atom_cost.
  ///
  std::optional<std::vector<Eigen::MatrixXd>>
      prim_sym_invariant_displacement_modes;
};

/// \brief Holds prim and lattice mapping-specific data used
///     for mapping searches
struct LatticeMappingSearchData {
  /// \brief Constructor
  LatticeMappingSearchData(std::shared_ptr<PrimSearchData> _prim_data,
                           std::shared_ptr<StructureSearchData> _structure_data,
                           LatticeMapping const &_lattice_mapping);

  /// \brief Holds prim-related data used for mapping searches
  std::shared_ptr<PrimSearchData> prim_data;

  /// \brief Holds structure-related data
  std::shared_ptr<StructureSearchData> structure_data;

  /// \brief A specific lattice mapping
  LatticeMapping lattice_mapping;

  /// \brief The transformation matrix that gives the ideal
  ///     superstructure lattice, for this lattice mapping,
  ///     from the prim lattice (i.e. T*N of the lattice
  ///     mapping)
  ///
  /// This is equivalent to:
  ///
  ///     lround(lattice_mapping.transformation_matrix_to_super *
  ///         lattice_mapping.reorientation)
  ///
  Eigen::Matrix3l transformation_matrix_to_super;

  /// \brief The lattice of the ideal supercell.
  xtal::Lattice supercell_lattice;

  /// \brief Generates ideal superstructure sites
  xtal::UnitCellCoordIndexConverter unitcellcoord_index_converter;

  /// \brief Number of sites in the ideal superstructure specified
  ///     by the lattice_mapping
  Index N_supercell_site;

  /// \brief Matrix of shape=(3,N_supercell_site) containing
  ///     the Cartesian coordinates of the child structure atoms,
  ///     in the state after the inverse lattice mapping deformation
  ///     is applied (\f$F^{-1}\vec{r_2}\f$). The supercell refers
  ///     to the ideal superlattice, S1 = L1 * T * N.
  Eigen::MatrixXd atom_coordinate_cart_in_supercell;

  /// \brief A shape=(3,N_supercell_site) matrix of Cartesian
  ///     coordinates of the sites in the ideal supercell.
  Eigen::MatrixXd supercell_site_coordinate_cart;

  /// \brief Size=N_supercell_site, with names of atoms allowed on
  ///     each supercell site
  std::vector<std::vector<std::string>> supercell_allowed_atom_types;
};

/// \brief Make possible atom -> site translations to bring atoms into
///     registry with the sites.
std::vector<Eigen::Vector3d> make_trial_translations(
    LatticeMappingSearchData const &lattice_mapping_data);

/// \brief A function, such as `mapping_impl::make_atom_mapping_cost`,
///     which calculates the atom-to-site mapping cost given the
///     site-to-atom displacement, atom_type, allowed_atom_types,
///     and value to use unallowed mappings (infinity).
using AtomMappingCostFunction = std::function<double(
    Eigen::Vector3d const &displacement, std::string const &atom_type,
    std::vector<std::string> const &allowed_atom_types, double infinity)>;

/// \brief Make the atom mapping cost for a particular atom
///     to a particular structure site
double make_atom_mapping_cost(
    Eigen::Vector3d const &displacement, std::string const &atom_type,
    std::vector<std::string> const &allowed_atom_types, double infinity);

/// \brief Holds data shared amongst all potential atom-to-site
///     assignment problems making use of the same trial
///     translation
struct AtomMappingSearchData {
  /// \brief Constructor
  AtomMappingSearchData(
      std::shared_ptr<LatticeMappingSearchData> _lattice_mapping_data,
      Eigen::Vector3d const &_trial_translation_cart,
      AtomMappingCostFunction _atom_mapping_cost_f = make_atom_mapping_cost,
      double _infinity = 1e20);

  /// \brief Holds lattice mapping-specific data used
  ///     for mapping searches
  std::shared_ptr<LatticeMappingSearchData> lattice_mapping_data;

  /// \brief A Cartesian translation applied to atom
  ///     coordinates in the ideal superstructure setting
  ///     (i.e. atom_coordinate_cart_in_supercell) to
  ///     bring the atoms into alignment with ideal
  ///     superstructure sites.
  Eigen::Vector3d trial_translation_cart;

  /// \brief Shape=(N_supercell_site, N_atom), the
  ///     site-to-atom displacements, of minimum length under
  ///     periodic boundary conditions of the ideal
  ///     superstructure.
  std::vector<std::vector<Eigen::Vector3d>> site_displacements;

  /// \brief Shape=(N_supercell_site, N_supercell_site) cost
  ///     matrix used in atom to site assignment problem. The
  ///     element `cost_matrix(i, j)` is set to the cost of
  ///     mapping the j-th atom to the i-th site. If there
  ///     are more sites than atoms, vacancies are added.
  Eigen::MatrixXd cost_matrix;
};

}  // namespace mapping
}  // namespace CASM

#endif
