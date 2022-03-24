#include "casm/crystallography/SimpleStructure.hh"
#include "casm/mapping/SearchData.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"

using namespace CASM;
using namespace CASM::mapping;

xtal::SimpleStructure make_conventional_fcc_ideal() {
  xtal::SimpleStructure simple;

  // clang-format off
  simple.lat_column_mat <<
    4., 0., 0.,
    0., 4., 0.,
    0., 0., 4.;

  simple.atom_info.coords = Eigen::MatrixXd(3, 4);
  simple.atom_info.coords <<
    0., 2., 0., 0.,
    0., 2., 2., 2.,
    0., 0., 2., 2.;

  simple.atom_info.names = {"A", "A", "A", "B"};
  // clang-format on

  return simple;
}

TEST(StructureSearchDataTest, Test1) {
  xtal::SimpleStructure s = make_conventional_fcc_ideal();
  auto structure_data = std::make_shared<StructureSearchData const>(
      xtal::Lattice(s.lat_column_mat), s.atom_info.coords, s.atom_info.names);
  EXPECT_TRUE(
      almost_equal(structure_data->lattice.lat_column_mat(), s.lat_column_mat));
  EXPECT_TRUE(
      almost_equal(structure_data->atom_coordinate_cart, s.atom_info.coords));
  EXPECT_EQ(structure_data->atom_type, s.atom_info.names);
  EXPECT_EQ(structure_data->structure_factor_group.size(), 16);
  EXPECT_EQ(structure_data->structure_crystal_point_group.size(), 16);
}
