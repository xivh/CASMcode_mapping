#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/definitions.hh"
#include "casm/mapping/LatticeMap.hh"
#include "casm/mapping/SimpleStrucMapCalculator.hh"
#include "casm/mapping/StrucMapping.hh"
#include "casm/misc/CASM_Eigen_math.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

/// \brief Maintain k-best results
///
/// \param k_best The optional number of solutions, T, to keep. Solutions
/// approximately equal to
///     the k-th best score are also kept, in `overflow`.
/// \param cost_tol Tolerance used for comparing costs
/// \param results A map of cost -> solution
/// \param overflow A map of cost -> solution for solutions approximately equal
///     to results.rbegin()
/// \param max_cost The max cost of solution to keep in results (kept if cost
///     < max_cost + cost_tol). If results.size() > *k_best, the max_cost is
///     shrunk to results.rbegin()->first.
///
/// This does the following:
/// - Nothing if !k_best.has_value()
/// - While results.size() > *k_best:
///   - Move results approximately equal to results.rbegin() to overflow, until
///     results.size() == *k_best
///   - Shrink max_cost to equal results.rbegin()->first
template <typename K, typename T>
void maintain_k_best_results(std::optional<int> const &k_best, double cost_tol,
                             std::multimap<K, T> &results,
                             std::multimap<K, T> &overflow, double &max_cost) {
  if (!k_best.has_value()) {
    return;
  }
  // only keep k-best results, but include those ~equal to k-th best result
  while (static_cast<int>(results.size()) > *k_best) {
    // check if last two results are approximately equal
    // reverse iterator, points at extra node
    auto last_it = --results.end();

    // reverse iterator, points at k-th best result
    auto nexttolast_it = --(--results.end());

    // if extra is not approximately equal cost to k-th best, erase extra and
    // overflow
    if ((last_it->first.cost - nexttolast_it->first.cost) > cost_tol) {
      results.erase(last_it);
      overflow.clear();
    }
    // if extra is approximately equal to last, move extra to overflow
    else {
      overflow.insert(results.extract(last_it));
    }
    // shrink max_cost
    max_cost = results.rbegin()->first.cost;
  }
}

/// \brief Lattice mapping transformation
///
/// The mapping transformation has the form:
///
/// \f[
///     F L_1 T N = L_2.
/// \f]
///
/// where \f$L_1\f$ are the reference "parent" lattice vectors, and
/// \f$L_2\f$ are the "child" lattice vectors being mapped to the
/// parent vectors, as columns of shape=(3,3) matrices. The other
/// shape=(3,3) matrices are:
///
/// - \f$F\f$, the shape=(3,3) parent-to-child deformation tensor
/// - \f$T\f$, an integer transformation matrix that generates a
///   superlattice of \f$L_1\f$
/// - \f$N\f$, a unimodular reorientation matrix that generates a
///   lattice equivalent to \f$L_1 T\f$ with reoriented lattice
///   vectors
///
/// There are infinitely many potential choices of N, which leads to
/// infinitely many potential lattice mapping transformations between
/// two lattices. A cost function of the stretch, \f$V\f$, is used
/// to score mappings.
///
struct LatticeMapping {
  /// See class documentation for definitions
  LatticeMapping(Eigen::Matrix3d const &_deformation_tensor,
                 Eigen::Matrix3d const &_transformation_matrix_to_super,
                 Eigen::Matrix3d const &_reorientation)
      : deformation_tensor(_deformation_tensor),
        transformation_matrix_to_super(_transformation_matrix_to_super),
        reorientation(_reorientation),
        right_stretch(polar_decomposition(deformation_tensor)),
        isometry(deformation_tensor * right_stretch.inverse()),
        left_stretch(deformation_tensor * isometry.transpose()) {}

  /// \brief Parent-to-child deformation tensor, F
  Eigen::Matrix3d deformation_tensor;

  /// \brief Parent lattice transformation_matrix_to_super, T
  Eigen::Matrix3d transformation_matrix_to_super;

  /// \brief Parent superlattice reorienation matrix, N
  Eigen::Matrix3d reorientation;

  /// \brief U, of F = Q * U = V * Q
  Eigen::Matrix3d right_stretch;

  /// \brief Rigid transformation (isometry) matrix Q, of F = Q * U = V * Q
  ///
  /// Note: Q.inverse() == Q.transpose()
  Eigen::Matrix3d isometry;

  /// \brief V, of F = Q * U = V * Q
  Eigen::Matrix3d left_stretch;
};

/// \brief Parent-to-child deformation tensor, F
Eigen::Matrix3d get_latticemapping_deformation_tensor(
    LatticeMapping const &lattice_mapping) {
  return lattice_mapping.deformation_tensor;
}

Eigen::Matrix3d get_latticemapping_transformation_matrix_to_super(
    LatticeMapping const &lattice_mapping) {
  return lattice_mapping.transformation_matrix_to_super;
}

Eigen::Matrix3d get_latticemapping_reorientation(
    LatticeMapping const &lattice_mapping) {
  return lattice_mapping.reorientation;
}

/// \brief Isometry of parent-to-child deformation tensor, F = Q U = V Q
Eigen::Matrix3d get_latticemapping_isometry(
    LatticeMapping const &lattice_mapping) {
  return lattice_mapping.isometry;
}

/// \brief Left stretch of parent-to-child deformation tensor, F = V Q
Eigen::Matrix3d get_latticemapping_left_stretch(
    LatticeMapping const &lattice_mapping) {
  return lattice_mapping.left_stretch;
}

/// \brief Right stretch of parent-to-child deformation tensor, F = Q U
Eigen::Matrix3d get_latticemapping_right_stretch(
    LatticeMapping const &lattice_mapping) {
  return lattice_mapping.right_stretch;
}

struct LatticeMappingKey {
  LatticeMappingKey(double _cost, double _cost_tol, xtal::Lattice _lattice)
      : cost(_cost), cost_tol(_cost_tol), lattice(_lattice) {}

  double cost;
  double cost_tol;
  xtal::Lattice lattice;

  bool operator<(LatticeMappingKey const &RHS) const {
    if (this->cost < RHS.cost - cost_tol) {
      return true;
    }
    if (this->cost > RHS.cost + cost_tol) {
      return false;
    }
    return this->lattice > RHS.lattice;
  }
};

/// \param lattice1 The referece "parent" lattice
/// \param lattice2 The "child" lattice
/// \param T Mapping is performed for L1 * T * N <-> L2, where T is an
///     approximately integer transformation matrix to a supercell of L1.
///     The default value is the identity matrix.
/// \param reorientation_range The absolute value of the maximum element in
///     reorientation matrix, N. This determines how many equivalent lattice
///     vector reorientations are checked. Usually 1 is sufficient.
/// \param L1_point_group Used to skip reorientation matrices that result in
///     symmetrically equivalent mappings. The default (empty), is
///     equivalent to only including the identity operation.
/// \param L2_point_group Used to skip reorientation matrices that
///     result in symmetrically equivalent mappings. The default (empty), is
///     equivalent to just including the identity operation.
/// \param min_cost Keep results with cost >= min_cost
/// \param max_cost Keep results with cost <= max_cost
/// \param cost_method One of "isotropic_strain_cost" or
///     "symmetry_breaking_strain_cost"
/// \param k_best If k_best.has_value(), then only keep the k_best results
///     satisfying the min_cost and max_cost constraints. If there are
///     approximate ties, those will also be kept.
/// \param cost_tol Tolerance for checking if lattice mapping costs are
///     approximately equal
std::vector<std::pair<double, LatticeMapping>> map_lattices(
    xtal::Lattice const &lattice1, xtal::Lattice const &lattice2,
    std::optional<Eigen::Matrix3d> T = std::nullopt,
    int reorientation_range = 1,
    std::vector<xtal::SymOp> lattice1_point_group = std::vector<xtal::SymOp>{},
    std::vector<xtal::SymOp> lattice2_point_group = std::vector<xtal::SymOp>{},
    double min_cost = 0.0, double max_cost = 1e20,
    std::string cost_method = std::string("isotropic_strain_cost"),
    std::optional<int> k_best = std::nullopt, double cost_tol = 1e-5) {
  double init_better_than = 1e20;
  bool symmetrize_strain_cost;
  if (cost_method == "isotropic_strain_cost") {
    symmetrize_strain_cost = false;
  } else if (cost_method == "symmetry_breaking_strain_cost") {
    symmetrize_strain_cost = true;
  } else {
    throw std::runtime_error(
        "Error in map_lattices: cost_method not recognized");
  }
  if (k_best.has_value() && *k_best < 1) {
    throw std::runtime_error(
        "Error in map_lattices: k_best < 1 is not allowed");
  }
  if (lattice1_point_group.empty()) {
    lattice1_point_group.push_back(xtal::SymOp::identity());
  }
  if (lattice2_point_group.empty()) {
    lattice2_point_group.push_back(xtal::SymOp::identity());
  }
  if (!T.has_value()) {
    T = Eigen::Matrix3d::Identity();
  }
  Eigen::Matrix3d L1 = lattice1.lat_column_mat();
  xtal::Lattice parent_superlattice(L1 * T.value(), lattice1.tol());
  using mapping_v1::LatticeMap;
  LatticeMap latmap(parent_superlattice, lattice2, reorientation_range,
                    lattice1_point_group, lattice2_point_group,
                    init_better_than, symmetrize_strain_cost, cost_tol);

  // the k-best results
  std::multimap<LatticeMappingKey, LatticeMapping> results;

  // results that are approximately equal to last-place result
  std::multimap<LatticeMappingKey, LatticeMapping> overflow;

  while (latmap) {
    double cost = latmap.strain_cost();
    if (cost > (min_cost - cost_tol) && cost < (max_cost + cost_tol)) {
      results.emplace(
          LatticeMappingKey(cost, cost_tol,
                            xtal::Lattice(L1 * T.value() * latmap.matrixN())),
          LatticeMapping(latmap.deformation_gradient(), T.value(),
                         latmap.matrixN()));

      // maintain results.size() <= *k_best, keep approximately equal results,
      // shrinks max_cost to results.rbegin() if results.size() == *k_best
      maintain_k_best_results(k_best, cost_tol, results, overflow, max_cost);
    }

    // iterate reorientation matrices until finding one with cost < max_cost +
    // cost_tol
    latmap.next_mapping_better_than(max_cost);
  }

  while (overflow.size()) {
    results.insert(overflow.extract(overflow.begin()));
  }

  std::vector<std::pair<double, LatticeMapping>> final;
  auto it = results.begin();
  while (it != results.end()) {
    final.emplace_back(it->first.cost, it->second);
    ++it;
  }
  return final;
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
///
/// There are in general many potential choices of atom mapping, with different
/// permutations and displacements. An assignment algorithm such as the
/// Hungarian method, with a cost function that depends on the displacements,
/// \f$\vec{d}(i)\f$, is used to score mappings.
///
struct AtomMapping {
  AtomMapping(Eigen::MatrixXd const &_displacement,
              std::vector<Index> const &_permutation,
              Eigen::Vector3d const &_translation, bool _time_reversal)
      : displacement(_displacement),
        permutation(_permutation),
        translation(_translation),
        time_reversal(_time_reversal) {}

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
  StructureMapping(LatticeMapping const &_lattice_mapping,
                   AtomMapping const &_atom_mapping)
      : lattice_mapping(_lattice_mapping), atom_mapping(_atom_mapping) {}

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
  StructureMappingScore(double _lattice_cost, double _atom_cost,
                        double _total_cost)
      : lattice_cost(_lattice_cost),
        atom_cost(_atom_cost),
        total_cost(_total_cost) {}
  double lattice_cost;
  double atom_cost;
  double total_cost;
};

double get_lattice_cost(StructureMappingScore const &score) {
  return score.lattice_cost;
}

double get_atom_cost(StructureMappingScore const &score) {
  return score.atom_cost;
}

double get_total_cost(StructureMappingScore const &score) {
  return score.total_cost;
}

/// \param prim The reference "parent" structure
/// \param structure2 The "child" structure
/// \param max_vol The maximum parent superstructure volume to consider, as
///     a multiple of the parent structure volume
/// \param prim_factor_group Used to skip mappings that are
///     symmetrically equivalent mappings. The default (empty), is
///     equivalent to only including the identity operation.
/// \param structure2_factor_group Used to skip mappings that are
///     symmetrically equivalent mappings. The default (empty), is
///     equivalent to only including the identity operation.
/// \param min_vol The minimum parent superstructure volume to consider, as
///     a multiple of the parent structure volume. Default=1.
/// \param min_cost Keep results with total cost >= min_cost
/// \param max_cost Keep results with total cost <= max_cost
/// \param lattice_cost_weight The fraction of the total cost due to
///     the lattice strain cost. The remaining fraction
///     (1.-lattice_cost_weight) is due to the atom cost. Default=0.5.
/// \param strain_cost_method One of "isotropic_strain_cost" or
///     "symmetry_breaking_strain_cost"
/// \param atom_cost_method One of "isotropic_atom_cost" or
///     "symmetry_breaking_atom_cost"
/// \param k_best Keep the k_best results satisfying the min_cost and
///     max_cost constraints. If there are approximate ties, those
///     will also be kept.
/// \param cost_tol Tolerance for checking if lattice mapping costs are
///     approximately equal
std::vector<std::pair<StructureMappingScore, StructureMapping>> map_structures(
    xtal::BasicStructure const &prim, xtal::SimpleStructure const &structure2,
    Index max_vol,
    std::vector<xtal::SymOp> prim_factor_group = std::vector<xtal::SymOp>{},
    std::vector<xtal::SymOp> structure2_factor_group =
        std::vector<xtal::SymOp>{},
    Index min_vol = 1, double min_cost = 0.0, double max_cost = 1e20,
    double lattice_cost_weight = 0.5,
    std::string strain_cost_method = std::string("isotropic_strain_cost"),
    std::string atom_cost_method = std::string("isotropic_atom_cost"),
    int k_best = 1, double cost_tol = 1e-5) {
  bool symmetrize_strain_cost;
  if (strain_cost_method == "isotropic_strain_cost") {
    symmetrize_strain_cost = false;
  } else if (strain_cost_method == "symmetry_breaking_strain_cost") {
    symmetrize_strain_cost = true;
  } else {
    throw std::runtime_error(
        "Error in map_lattices: strain_cost_method not recognized");
  }

  bool symmetrize_atom_cost;
  if (atom_cost_method == "isotropic_atom_cost") {
    symmetrize_atom_cost = false;
  } else if (atom_cost_method == "symmetry_breaking_atom_cost") {
    symmetrize_atom_cost = true;
  } else {
    throw std::runtime_error(
        "Error in map_lattices: atom_cost_method not recognized");
  }

  if (k_best < 1) {
    throw std::runtime_error(
        "Error in map_lattices: k_best < 1 is not allowed");
  }
  if (prim_factor_group.empty()) {
    prim_factor_group.push_back(xtal::SymOp::identity());
  }
  if (structure2_factor_group.empty()) {
    structure2_factor_group.push_back(xtal::SymOp::identity());
  }
  if (min_vol < 1) {
    throw std::runtime_error("Error in map_structures: min_vol < 1");
  }
  if (max_vol < min_vol) {
    throw std::runtime_error("Error in map_structures: max_vol < min_vol");
  }

  using mapping_v1::MappingNode;
  using mapping_v1::SimpleStrucMapCalculator;
  using mapping_v1::StrucMapper;

  /// For the StrucMapper::map_deformed_struc_impose_lattice_vols method:
  /// - If invalid values of `min_vol` or `max_vol` are provided (negative
  /// values
  ///   or max_vol < min_vol), then this method throws.
  /// - Parameters `min_va_frac`, `max_va_frac`, `max_volume_change`,
  ///   and `soft_va_limit` have no effect.
  /// - `robust` search is used anyway if k_best > 1

  SimpleStrucMapCalculator calculator(
      xtal::make_simple_structure(prim), prim_factor_group,
      CASM::xtal::SimpleStructure::SpeciesMode::ATOM,
      xtal::allowed_molecule_names(prim));
  double _max_volume_change = 0.5;  // no effect
  bool _robust = true;              // no effect if k_best > 1
  bool _soft_va_limit = false;      // no effect
  double _min_va_frac = 0.;         // no effect
  double _max_va_frac = 1.;         // no effect
  std::cout << "lattice_cost_weight: " << lattice_cost_weight << std::endl;
  std::cout << "cost_tol: " << cost_tol << std::endl;
  StrucMapper strucmap(calculator, lattice_cost_weight, _max_volume_change,
                       _robust, _soft_va_limit, cost_tol, _min_va_frac,
                       _max_va_frac);

  if (symmetrize_strain_cost) {
    std::cout << "symmetrize lattice cost" << std::endl;
    strucmap.set_symmetrize_lattice_cost(true);
  }
  if (symmetrize_atom_cost) {
    std::cout << "symmetrize atom cost" << std::endl;
    auto prim_permute_group =
        xtal::make_permutation_representation(prim, prim_factor_group);
    strucmap.set_symmetrize_atomic_cost(true, prim_factor_group,
                                        prim_permute_group);
  }

  bool keep_invalid = false;
  std::cout << "min_vol: " << min_vol << std::endl;
  std::cout << "max_vol: " << max_vol << std::endl;
  std::cout << "k_best: " << k_best << std::endl;
  std::cout << "min_cost: " << min_cost << std::endl;
  std::cout << "max_cost: " << max_cost << std::endl;
  std::set<MappingNode> mappings =
      strucmap.map_deformed_struc_impose_lattice_vols(
          structure2, min_vol, max_vol, k_best, max_cost, min_cost,
          keep_invalid, structure2_factor_group);

  std::cout << "mappings.size(): " << mappings.size() << std::endl;
  std::vector<std::pair<StructureMappingScore, StructureMapping>> results;
  for (auto const &mapping : mappings) {
    if (!(mapping.cost > (min_cost - cost_tol) &&
          mapping.cost < (max_cost + cost_tol))) {
      continue;
    }
    // std::cout << "%%%%%%%%%%%%" << std::endl;
    // Get LatticeMapping data
    //
    // LatticeMapping convention is F * L1 * T * N = L2
    // LatticeNode convention is L1 * T * N = stretch * isometry * L2
    // LatticeNode convention is T * N =
    // lattice_node.parent.transformation_matrix_to_super()
    auto const &lattice_node = mapping.lattice_node;
    Eigen::Matrix3d F =
        (lattice_node.stretch * lattice_node.isometry).inverse();
    // std::cout << "lattice_node.cost: " << lattice_node.cost << std::endl;
    // std::cout << "F:\n" << F << std::endl;
    // std::cout << "L:\n" <<
    // lattice_node.parent.prim_lattice().lat_column_mat() << std::endl;
    // std::cout << "S:\n" <<
    // lattice_node.parent.superlattice().lat_column_mat() << std::endl;
    Eigen::Matrix3d T =
        lattice_node.parent.transformation_matrix_to_super().cast<double>();
    // std::cout << "T:\n" << T << std::endl;
    Eigen::Matrix3d N = Eigen::Matrix3d::Identity();
    // std::cout << "vol:" << T.determinant() << std::endl;
    LatticeMapping lattice_mapping(F, T, N);

    // Get AtomMapping data
    auto const &atomic_node = mapping.atomic_node;
    Eigen::MatrixXd disp = mapping.atom_displacement;
    std::vector<Index> perm = mapping.atom_permutation;
    Eigen::Vector3d trans = atomic_node.translation;
    bool time_reversal = atomic_node.time_reversal;
    AtomMapping atom_mapping(disp, perm, trans, time_reversal);

    results.emplace_back(StructureMappingScore(lattice_node.cost,
                                               atomic_node.cost, mapping.cost),
                         StructureMapping(lattice_mapping, atom_mapping));

    // std::cout << std::endl;
  }

  return results;
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
      - :math:`F` is the parent-to-child deformation tensor,
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
           py::arg("deformation_tensor"),
           py::arg("transformation_matrix_to_super"), py::arg("reorientation"),
           R"pbdoc(
          Construct a lattice mapping

          Parameters
          ----------

          deformation_tensor : array_like, shape=(3,3)
              The parent-to-child deformation tensor,:math:`F`, a shape=(3,3)
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
      .def("deformation_tensor", &get_latticemapping_deformation_tensor,
           "Return the shape=(3,3) parent-to-child deformation tensor.")
      .def("transformation_matrix_to_super",
           &get_latticemapping_transformation_matrix_to_super,
           "Return the shape=(3,3) parent supercell transformation matrix, "
           ":math:`T`.")
      .def("reorientation", &get_latticemapping_reorientation,
           "Return the shape=(3,3) unimodular matrix, :math:`N`.")
      .def("isometry", &get_latticemapping_isometry,
           "Return the shape=(3,3) isometry matrix, :math:`Q`, of the "
           "parent-to-child deformation tensor.")
      .def("left_stretch", &get_latticemapping_left_stretch,
           "Return the shape=(3,3) left symmetric stretch tensor, :math:`V`, "
           "of the parent-to-child deformation tensor.")
      .def("right_stretch", &get_latticemapping_right_stretch,
           "Return the shape=(3,3) right symmetric stretch tensor, :math:`U`, "
           "of the parent-to-child deformation tensor.");

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
     - :math:`F` is the parent-to-child deformation tensor from a
       :class:`~casm.mapping.LatticeMapping`.
     - :math:`p_i` is a permutation vector, specifying which atom in the
       child structure (:math:`p_i`) is mapped to the i-th site of the parent
       superstructure. Values of :math:`p_i` greater than the number of atoms
       in the child structure indicate inferred vacancies in mappings to
       :class:`~casm.xtal.Prim` with vacancies allowed.
     - :math:`\vec{t}` is a translation vector, in Cartesian coordinates, usually
       chosen so the average displacement is zero.
     - :math:`\vec{d}(i)`: The displacement associated with the atom at the
       i-th site in parent superstructure.

     Additionally, structures with magnetic spin may be related by time
     reversal symmetry.
     )pbdoc")
      .def(py::init<Eigen::MatrixXd const &, std::vector<Index> const &,
                    Eigen::Vector3d const &, bool>(),
           py::arg("displacement"), py::arg("permutation"),
           py::arg("translation"), py::arg("time_reversal"), R"pbdoc(
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
          time_reversal : bool, default=False
             If True, the parent and child structures are related by time
             reversal symmetry.
           )pbdoc")
      .def("displacement", &get_atom_mapping_displacement,
           "Return the shape=(3,n) matrix whose columns are the atom "
           "displacements :math:`\vec{d}(i)`.")
      .def("permutation", &get_atom_mapping_permutation,
           "Return the permutation vector, :math:`p_i`.")
      .def("translation", &get_atom_mapping_translation,
           "Return the translation vector, :math:`\vec{t}`.")
      .def("time_reversal", &get_atom_mapping_time_reversal,
           "Return True if the parent and child structures are related by time "
           "reversal symmetry.");

  py::class_<StructureMapping>(m, "StructureMapping", R"pbdoc(
    A mapping between two structures

    A structure mapping is a combination of:

    - A :class:`~casm.mapping.LatticeMapping`
    - An :class:`~casm.mapping.AtomMapping`

    See those class descriptions for details of the mappings.

    )pbdoc")
      .def(py::init<LatticeMapping const &, AtomMapping const &>(),
           py::arg("lattice_mapping"), py::arg("atom_mapping"), R"pbdoc(
          Construct a structure mapping

          Parameters
          ----------

          lattice_mapping : casm.mapping.LatticeMapping
              A :class:`~casm.mapping.LatticeMapping`
          atom_mapping : casm.mapping.AtomMapping
              An :class:`~casm.mapping.AtomMapping`
          )pbdoc")
      .def("lattice_mapping", &get_lattice_mapping,
           "Return the :class:`~casm.mapping.LatticeMapping`.")
      .def("atom_mapping", &get_atom_mapping,
           "Return the :class:`~casm.mapping.AtomMapping`.");

  py::class_<StructureMappingScore>(m, "StructureMappingScore", R"pbdoc(
   Holds mapping scores between two structures

   The structure mapping score is a combination of:

   - A lattice mapping score, lattice_cost
   - An atom mapping score, atom_cost

   The total structure mapping score, total_cost, is lattice_cost_weight*lattice_cost + (1.0 - lattice_cost_weight)*atom_cost, where lattice_cost_weight is an input to the :func:`~casm.mapping.map_structures` method.

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
      .def("lattice_cost", &get_lattice_cost,
           "Return the lattice mapping cost.")
      .def("atom_cost", &get_atom_cost, "Return the atom mapping cost.")
      .def("total_cost", &get_total_cost, "Return the total mapping cost.");

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
      - :math:`F` is the parent-to-child deformation tensor,
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

      Lattice mappings are scored and ranked according to one of:

      - "isotropic_strain_cost": a strain cost, calculated to be
        volume-normalized and invariant to which structure is the
        parent/child. Calculated as:

        .. math::

            \mathrm{cost} &= \frac{1}{2}*\left( \frac{1}{3}*\mathrm{tr}\left(\tilde{B}_{c \to p}^{2} \right) + \frac{1}{3}*\mathrm{tr}\left(\tilde{B}_{p \to c}^{2} \right) \right)

            \tilde{F} &= \frac{1}{\det{F}^{1/3}} F = \tilde{Q} \tilde{U}

            \tilde{B} &= \tilde{U} - I,

        where :math:`F` is the deformation tensor, which can be defined
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
          Selects the method used to score lattice mappings. One of
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
      (:class:`~cast.mapping.LatticeMapping`) and for each potential lattice
      mapping finding atom mappings (:class:`~cast.mapping.AtomMapping`).

      The total structure mapping score, total_cost, is a weighted mixture of
      the lattice mapping score, lattice_cost, and the atom mapping score,
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
          Keep lattice mappings with cost >= min_cost
      max_cost : float, default=1e20
          Keep results with cost <= max_cost
      lattice_cost_weight : float, default=0.5
          The fraction of the total cost due to the lattice strain cost.
          The remaining fraction (1.-lattice_cost_weight) is due to the
          atom cost.
      strain_cost_method : str, default="isotropic_strain_cost"
          Selects the method used to score lattice mappings. One of
          "isotropic_strain_cost" or "symmetry_breaking_strain_cost".
      atom_cost_method : str, default="isotropic_strain_cost"
          Selects the method used to score atom mappings. One of
          "isotropic_atom_cost" or "symmetry_breaking_atom_cost".
      k_best : int, default=1
          Only keep the k-best results (i.e. k mappings with minimum cost)
          satisfying the min_cost and max_cost constraints. If there are
          approximate ties, those will also be kept.
      cost_tol : float, default=1e-5
          Tolerance for checking if lattice mapping costs are approximately
          equal.

      Returns
      -------
      structure_mappings : List[Tuple[casm.mapping.StructureMappingScore, casm.mapping.StructureMapping]]
          A list of tuple of :class:`~casm.mapping.StructureMappingScore`
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

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
