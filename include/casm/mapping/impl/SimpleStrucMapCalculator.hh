#ifndef CASM_mapping_SimpleStrucMapCalculator
#define CASM_mapping_SimpleStrucMapCalculator

#include "casm/crystallography/Adapter.hh"
#include "casm/mapping/impl/StrucMapCalculatorInterface.hh"

namespace CASM {

namespace xtal {
class SimpleStructure;
struct SymOp;
typedef std::vector<SymOp> SymOpVector;
}  // namespace xtal

namespace mapping_impl {

// In this file:
struct MappingNode;
class StrucMapCalculatorInterface;
class SimpleStrucMapCalculator;

// TODO:
// - Explain _factor_group default / non-default behavior
// - Does species_mode SpeciesMode::MOL work? how does it behave?

class SimpleStrucMapCalculator : public StrucMapCalculatorInterface {
 public:
  /// StrucMapCalculatorInterface constructor
  ///
  /// \param _parent Reference structure to be mapped to
  /// \param _factor_group Factor group of the parent structure.
  /// \param _species_mode Specifies whether to map to parent atoms or
  ///     molecules. Use `SimpleStructure::SpeciesMode::ATOM`.
  /// \param allowed_species Names of allowed species on each parent structure
  ///     site. If empty, `allowed_species[site_index]` is set to the 1-element
  ///     vector with the name of the current species on the parent structure
  ///     site. Use `allowed_molecule_names` to use the names of
  ///     `xtal::Molecule` allowed on `xtal::BasicStructure` sites.
  ///
  SimpleStrucMapCalculator(
      xtal::SimpleStructure _parent,
      xtal::SymOpVector const &_factor_group = {xtal::SymOp::identity()},
      xtal::SimpleStructure::SpeciesMode species_mode =
          xtal::SimpleStructure::SpeciesMode::ATOM,
      AllowedSpecies allowed_species = {})
      : StrucMapCalculatorInterface(std::move(_parent), _factor_group,
                                    species_mode, std::move(allowed_species)) {}

  /// StrucMapCalculatorInterface constructor
  ///
  /// \param _parent Reference structure to be mapped to
  /// \param _factor_group Factor group of the parent structure.
  /// \param _species_mode Specifies whether to map to parent atoms or
  ///     molecules. Use `SimpleStructure::SpeciesMode::ATOM`.
  /// \param allowed_species Names of allowed species on each parent structure
  ///     site. If empty, `allowed_species[site_index]` is set to the 1-element
  ///     vector with the name of the current species on the parent structure
  ///     site. Use `allowed_molecule_names` to use the names of
  ///     `xtal::Molecule` allowed on `xtal::BasicStructure` sites.
  template <typename ExternSymOpVector>
  SimpleStrucMapCalculator(
      xtal::SimpleStructure _parent,
      ExternSymOpVector const &_factor_group = {xtal::SymOp::identity()},
      xtal::SimpleStructure::SpeciesMode species_mode =
          xtal::SimpleStructure::SpeciesMode::ATOM,
      AllowedSpecies allowed_species = {})
      : SimpleStrucMapCalculator(
            _parent,
            adapter::Adapter<xtal::SymOpVector, ExternSymOpVector>()(
                _factor_group),
            species_mode, allowed_species) {}

  virtual ~SimpleStrucMapCalculator() {}

  /// Constructs a list of prospective mapping translations
  std::vector<Eigen::Vector3d> translations(
      MappingNode const &_node,
      xtal::SimpleStructure const &child_struc) const override;

  /// Creates a copy of the child structure and applies mapping
  virtual xtal::SimpleStructure resolve_setting(
      MappingNode const &_node,
      xtal::SimpleStructure const &_child_struc) const override;

  /// Sets MappingNode data based on lattice and atomic mapping results
  void finalize(MappingNode &_node, xtal::SimpleStructure const &child_struc,
                bool const &symmetrize_atomic_cost = false) const override;

  /// Populates the cost matrix for the atomic assignment problem
  bool populate_cost_mat(
      MappingNode &_node,
      xtal::SimpleStructure const &child_struc) const override;

 private:
  /// \brief Make an exact copy of the calculator (including any initialized
  /// members)
  virtual StrucMapCalculatorInterface *_clone() const override {
    return new SimpleStrucMapCalculator(*this);
  }

  /// \brief Make an exact copy of the calculator (including any initialized
  /// members)
  virtual StrucMapCalculatorInterface *_quasi_clone(
      xtal::SimpleStructure _parent,
      xtal::SymOpVector const &_factor_group = {xtal::SymOp::identity()},
      xtal::SimpleStructure::SpeciesMode _species_mode =
          xtal::SimpleStructure::SpeciesMode::ATOM,
      AllowedSpecies _allowed_species = {}) const override {
    return new SimpleStrucMapCalculator(std::move(_parent), _factor_group,
                                        _species_mode,
                                        std::move(_allowed_species));
  }
};
}  // namespace mapping_impl
}  // namespace CASM
#endif
