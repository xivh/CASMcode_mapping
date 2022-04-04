#include "SearchTestData.hh"
#include "casm/mapping/SearchData.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;
using namespace CASM::xtal;
using namespace CASM::mapping;

// Test LatticeMappingSearchData construction using a structure
// generated from the conventional BCC unit cell
TEST(LatticeMappingSearchDataTest, Test1) {
  double latparam_a = 4.0;

  // F * L1 * T * N = L2
  Eigen::Matrix3d F;
  F << 1.0, 0.0, 0.0,  //
      0.0, 1.0, 0.1,   //
      0.0, 0.0, 1.1;   //
  Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d N = Eigen::Matrix3d::Identity();

  // F (r1_supercell[i] + disp) = r2[perm[i]] + trans
  Eigen::MatrixXd disp;
  disp.resize(3, 2);
  disp.col(0) << 0., 0., 0.;
  disp.col(1) << 0.01, 0.0, 0.0;
  std::vector<std::string> atom_type({"A", "A"});
  std::vector<Index> perm({0, 1});
  Eigen::Vector3d trans(0., 0., 0.);

  test::SearchTestData d(
      test::make_search_prim_binary_conventional_BCC(latparam_a), F, T, N, disp,
      atom_type, perm, trans);

  auto lattice_mapping_data = std::make_shared<LatticeMappingSearchData const>(
      d.prim_data, d.structure_data, d.lattice_mapping);

  EXPECT_TRUE(almost_equal(lattice_mapping_data->transformation_matrix_to_super,
                           lround(d.T * d.N)));
  EXPECT_TRUE(
      almost_equal(lattice_mapping_data->supercell_lattice.lat_column_mat(),
                   d.L1 * d.T * d.N));
  EXPECT_EQ(lattice_mapping_data->unitcellcoord_index_converter.total_sites(),
            2);
  EXPECT_EQ(lattice_mapping_data->N_supercell_site, 2);
  EXPECT_TRUE(
      almost_equal(lattice_mapping_data->atom_coordinate_cart_in_supercell,
                   d.F.inverse() * d.r2));
  EXPECT_TRUE(almost_equal(lattice_mapping_data->supercell_site_coordinate_cart,
                           d.r1_supercell));

  {
    std::vector<std::vector<std::string>> expected({{"A", "B"}, {"A", "B"}});
    EXPECT_EQ(lattice_mapping_data->supercell_allowed_atom_types, expected);
  }
}

// Test LatticeMappingSearchData construction using a structure
// generated from a supercell of the conventional BCC unit cell
// (T = [[1, 0, 0], [0, 1, 0], [0, 0, 2]])
TEST(LatticeMappingSearchDataTest, Test2) {
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
  Eigen::MatrixXd disp;
  disp.resize(3, 4);
  disp.col(0) << 0., 0., 0.;
  disp.col(1) << 0., 0., 0.;
  disp.col(2) << 0.01, 0.0, 0.0;
  disp.col(3) << 0., 0., 0.;
  std::vector<std::string> atom_type({"A", "A", "A", "A"});
  std::vector<Index> perm({0, 1, 2, 3});
  Eigen::Vector3d trans(0., 0., 0.);

  test::SearchTestData d(
      test::make_search_prim_binary_conventional_BCC(latparam_a), F, T, N, disp,
      atom_type, perm, trans);

  auto lattice_mapping_data = std::make_shared<LatticeMappingSearchData const>(
      d.prim_data, d.structure_data, d.lattice_mapping);

  EXPECT_TRUE(almost_equal(lattice_mapping_data->transformation_matrix_to_super,
                           lround(d.T * d.N)));
  EXPECT_TRUE(
      almost_equal(lattice_mapping_data->supercell_lattice.lat_column_mat(),
                   d.L1 * d.T * d.N));
  EXPECT_EQ(lattice_mapping_data->unitcellcoord_index_converter.total_sites(),
            4);
  EXPECT_EQ(lattice_mapping_data->N_supercell_site, 4);
  EXPECT_TRUE(
      almost_equal(lattice_mapping_data->atom_coordinate_cart_in_supercell,
                   d.F.inverse() * d.r2));
  EXPECT_TRUE(almost_equal(lattice_mapping_data->supercell_site_coordinate_cart,
                           d.r1_supercell));

  {
    std::vector<std::vector<std::string>> expected(
        {{"A", "B"}, {"A", "B"}, {"A", "B"}, {"A", "B"}});
    EXPECT_EQ(lattice_mapping_data->supercell_allowed_atom_types, expected);
  }
}
