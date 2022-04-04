#include "SearchTestData.hh"
#include "casm/mapping/SearchData.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;
using namespace CASM::xtal;
using namespace CASM::mapping;

namespace test {}  // namespace test

// Test LatticeMappingSearchData construction using a structure
// generated from the conventional BCC unit cell
TEST(AtomMappingSearchDataTest, Test1) {
  double latparam_a = 4.0;

  // F * L1 * T * N = L2
  Eigen::Matrix3d F;
  F << 1.0, 0.0, 0.0,  //
      0.0, 1.0, 0.1,   //
      0.0, 0.0, 1.1;   //
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

  // make AtomMappingSearchData, with trial_translation=(0., 0., 0.)
  auto atom_mapping_data = std::make_shared<AtomMappingSearchData const>(
      lattice_mapping_data, Eigen::Vector3d(0., 0., 0.));

  EXPECT_TRUE(almost_equal(atom_mapping_data->trial_translation_cart,
                           Eigen::Vector3d(0., 0., 0.)));
  Index N_supercell_site = lattice_mapping_data->N_supercell_site;
  Index N_atom = d.structure_data->N_atom;
  EXPECT_EQ(atom_mapping_data->site_displacements.size(), N_supercell_site);
  for (Index site_index = 0; site_index < N_supercell_site; ++site_index) {
    EXPECT_EQ(atom_mapping_data->site_displacements[site_index].size(), N_atom);
  }
  EXPECT_EQ(atom_mapping_data->cost_matrix.rows(), N_supercell_site);
  EXPECT_EQ(atom_mapping_data->cost_matrix.cols(), N_supercell_site);
}

// Test AtomMappingSearchData construction using a structure
// generated from a supercell of the conventional BCC unit cell
// (T = [[1, 0, 0], [0, 1, 0], [0, 0, 2]])
TEST(AtomMappingSearchDataTest, Test2) {
  double latparam_a = 4.0;

  // F * L1 * T * N = L2
  Eigen::Matrix3d F;
  F << 1.0, 0.0, 0.0,  //
      0.0, 1.0, 0.1,   //
      0.0, 0.0, 1.1;   //
  Eigen::Matrix3d T;
  T << 1.0, 0.0, 0.0,  //
      0.0, 1.0, 0.0,   //
      0.0, 0.0, 2.0;   //
  Eigen::Matrix3d N = Eigen::Matrix3d::Identity();

  // F (r1_supercell[i] + disp) = r2[perm[i]] + trans
  // structure1_supercell_atom_type[i] = structure2_atom_type[perm[i]]
  Eigen::MatrixXd disp;
  disp.resize(3, 4);
  disp.col(0) << 0., 0., 0.;
  disp.col(1) << 0., 0., 0.;
  disp.col(2) << 0.01, 0.0, 0.0;
  disp.col(3) << 0., 0., 0.;
  std::vector<std::string> structure1_supercell_atom_type({"A", "A", "A", "A"});
  std::vector<Index> perm({0, 1, 2, 3});
  Eigen::Vector3d trans(0., 0., 0.);

  test::SearchTestData d(
      test::make_search_prim_binary_conventional_BCC(latparam_a), F, T, N, disp,
      structure1_supercell_atom_type, perm, trans);

  auto lattice_mapping_data = std::make_shared<LatticeMappingSearchData const>(
      d.prim_data, d.structure_data, d.lattice_mapping);

  EXPECT_EQ(d.prim_data->prim_factor_group.size(), 48 * 2);

  // make AtomMappingSearchData, with trial_translation=(0., 0., 0.)
  auto atom_mapping_data = std::make_shared<AtomMappingSearchData const>(
      lattice_mapping_data, Eigen::Vector3d(0., 0., 0.));

  EXPECT_TRUE(almost_equal(atom_mapping_data->trial_translation_cart,
                           Eigen::Vector3d(0., 0., 0.)));
  Index N_supercell_site = lattice_mapping_data->N_supercell_site;
  Index N_atom = d.structure_data->N_atom;
  EXPECT_EQ(atom_mapping_data->site_displacements.size(), N_supercell_site);
  for (Index site_index = 0; site_index < N_supercell_site; ++site_index) {
    EXPECT_EQ(atom_mapping_data->site_displacements[site_index].size(), N_atom);
  }
  EXPECT_EQ(atom_mapping_data->cost_matrix.rows(), N_supercell_site);
  EXPECT_EQ(atom_mapping_data->cost_matrix.cols(), N_supercell_site);
}
