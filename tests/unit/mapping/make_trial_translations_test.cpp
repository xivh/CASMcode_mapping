#include "SearchTestData.hh"
#include "casm/mapping/SearchData.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;
using namespace CASM::xtal;
using namespace CASM::mapping;

namespace test {}  // namespace test

// Test make_trial_translations a structure
// generated from the conventional BCC unit cell
TEST(MakeTrialTranslationsTest, Test1) {
  double latparam_a = 4.0;

  // F * L1 * T * N = L2
  Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d N = Eigen::Matrix3d::Identity();

  // F (r1_supercell[i] + disp) = r2[perm[i]] + trans
  // structure1_supercell_atom_type[i] = structure2_atom_type[perm[i]]
  Eigen::MatrixXd disp;
  disp.resize(3, 1);
  disp.col(0) << 0.1, 0., 0.;
  std::vector<std::string> structure1_supercell_atom_type({"A"});
  std::vector<Index> perm({0});
  Eigen::Vector3d trans(0., 0., 0.);

  test::SearchTestData d(test::make_search_prim_binary_BCC(latparam_a), F, T, N,
                         disp, structure1_supercell_atom_type, perm, trans);

  auto lattice_mapping_data = std::make_shared<LatticeMappingSearchData const>(
      d.prim_data, d.structure_data, d.lattice_mapping);

  EXPECT_EQ(d.prim_data->prim_factor_group.size(), 48);
  EXPECT_EQ(d.structure_data->structure_factor_group.size(), 48);

  // only makes symmetrically unique trial translations
  std::vector<Eigen::Vector3d> trial_translations =
      make_trial_translations(*lattice_mapping_data);

  // there is 1 atom -> site translation, (0., 0., 0.)
  EXPECT_EQ(trial_translations.size(), 1);
  EXPECT_TRUE(almost_equal(trial_translations[0], -disp.col(0)))
      << trial_translations[0] << " != " << -disp.col(0);
}

// Test make_trial_translations a structure
// generated from the conventional BCC unit cell
TEST(MakeTrialTranslationsTest, Test2) {
  double latparam_a = 4.0;

  // F * L1 * T * N = L2
  Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d N = Eigen::Matrix3d::Identity();

  // F (r1_supercell[i] + disp) = r2[perm[i]] + trans
  // structure1_supercell_atom_type[i] = structure2_atom_type[perm[i]]
  Eigen::MatrixXd disp;
  disp.resize(3, 2);
  disp.col(0) << 0., 0., 0.;
  disp.col(1) << 0.01, 0.0, 0.0;
  std::vector<std::string> structure1_supercell_atom_type({"A", "A"});
  std::vector<Index> perm({0, 1});
  Eigen::Vector3d trans(0., 0., 0.);

  test::SearchTestData d(
      test::make_search_prim_binary_conventional_BCC(latparam_a), F, T, N, disp,
      structure1_supercell_atom_type, perm, trans);

  auto lattice_mapping_data = std::make_shared<LatticeMappingSearchData const>(
      d.prim_data, d.structure_data, d.lattice_mapping);

  EXPECT_EQ(d.prim_data->prim_factor_group.size(), 48 * 2);
  EXPECT_EQ(d.structure_data->structure_factor_group.size(), 16);

  // only makes symmetrically unique trial translations
  std::vector<Eigen::Vector3d> trial_translations =
      make_trial_translations(*lattice_mapping_data);

  // there are 2 atom -> site translations (0., 0., 0.) and (2., 2., 2.),
  // but only 1 is unique
  EXPECT_EQ(trial_translations.size(), 1);
  EXPECT_TRUE(almost_equal(trial_translations[0], Eigen::Vector3d(0., 0., 0.)));
}
