#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/mapping/SearchData.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;
using namespace CASM::mapping;

TEST(PrimSearchDataTest, Test1) {
  auto s =
      std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
  std::optional<std::vector<xtal::SymOp>> override_prim_factor_group =
      std::nullopt;
  bool enable_symmetry_breaking_atom_cost = true;
  auto prim_data = std::make_shared<PrimSearchData const>(
      s, override_prim_factor_group, enable_symmetry_breaking_atom_cost);

  EXPECT_TRUE(almost_equal(prim_data->prim_lattice.lat_column_mat(),
                           s->lattice().lat_column_mat()));
  Eigen::MatrixXd expected(3, 1);
  expected.col(0) = s->basis()[0].const_cart();
  EXPECT_TRUE(almost_equal(prim_data->prim_site_coordinate_cart, expected));
  EXPECT_EQ(prim_data->N_prim_site, 1);
  EXPECT_EQ(prim_data->prim_allowed_atom_types,
            std::vector<std::vector<std::string>>({{"A", "B"}}));
  EXPECT_EQ(prim_data->prim_factor_group.size(), 48);
  EXPECT_EQ(prim_data->prim_crystal_point_group.size(), 48);
  EXPECT_EQ(prim_data->prim_sym_invariant_displacement_modes.has_value(), true);
  EXPECT_EQ(prim_data->prim_sym_invariant_displacement_modes->size(), 0);
}

TEST(PrimSearchDataTest, Test2) {
  auto s = std::make_shared<xtal::BasicStructure const>(test::ZrO_prim());
  std::optional<std::vector<xtal::SymOp>> override_prim_factor_group =
      std::nullopt;
  bool enable_symmetry_breaking_atom_cost = true;
  auto prim_data = std::make_shared<PrimSearchData const>(
      s, override_prim_factor_group, enable_symmetry_breaking_atom_cost);

  EXPECT_TRUE(almost_equal(prim_data->prim_lattice.lat_column_mat(),
                           s->lattice().lat_column_mat()));

  {
    Eigen::MatrixXd expected(3, s->basis().size());
    Index i = 0;
    for (auto const &site : s->basis()) {
      expected.col(i++) = site.const_cart();
    }
    EXPECT_TRUE(almost_equal(prim_data->prim_site_coordinate_cart, expected));
  }

  EXPECT_EQ(prim_data->N_prim_site, 4);

  {
    std::vector<std::vector<std::string>> expected(
        {{"Zr"}, {"Zr"}, {"Va", "O"}, {"Va", "O"}});
    EXPECT_EQ(prim_data->prim_allowed_atom_types, expected);
  }

  EXPECT_EQ(prim_data->prim_factor_group.size(), 24);
  EXPECT_EQ(prim_data->prim_crystal_point_group.size(), 24);
  EXPECT_EQ(prim_data->prim_sym_invariant_displacement_modes.has_value(), true);
  EXPECT_EQ(prim_data->prim_sym_invariant_displacement_modes->size(), 0);
}

TEST(PrimSearchDataTest, Test3) {
  // test prim_sym_invariant_displacement_modes

  using namespace CASM::xtal;

  // lattice vectors as cols
  Eigen::Matrix3d lat;
  lat << 2.0, 0.0, 0.0,  //
      0.0, 2.0, 0.0,     //
      0.0, 0.0, 3.8;     //

  BasicStructure struc{Lattice{lat}};
  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");
  struc.push_back(
      Site(Coordinate(Eigen::Vector3d(0., 0., 0.), struc.lattice(), FRAC),
           std::vector<Molecule>{A}));
  struc.push_back(
      Site(Coordinate(Eigen::Vector3d(0.5, 0.5, 0.23), struc.lattice(), FRAC),
           std::vector<Molecule>{B}));
  struc.push_back(
      Site(Coordinate(Eigen::Vector3d(0.5, 0.5, 0.77), struc.lattice(), FRAC),
           std::vector<Molecule>{B}));

  auto s = std::make_shared<BasicStructure const>(struc);
  std::optional<std::vector<xtal::SymOp>> override_prim_factor_group =
      std::nullopt;
  bool enable_symmetry_breaking_atom_cost = true;
  auto prim_data = std::make_shared<PrimSearchData const>(
      s, override_prim_factor_group, enable_symmetry_breaking_atom_cost);

  EXPECT_TRUE(almost_equal(prim_data->prim_lattice.lat_column_mat(),
                           s->lattice().lat_column_mat()));

  {
    Eigen::MatrixXd expected(3, s->basis().size());
    Index i = 0;
    for (auto const &site : s->basis()) {
      expected.col(i++) = site.const_cart();
    }
    EXPECT_TRUE(almost_equal(prim_data->prim_site_coordinate_cart, expected));
  }

  EXPECT_EQ(prim_data->N_prim_site, 3);

  {
    std::vector<std::vector<std::string>> expected({{"A"}, {"B"}, {"B"}});
    EXPECT_EQ(prim_data->prim_allowed_atom_types, expected);
  }

  EXPECT_EQ(prim_data->prim_factor_group.size(), 16);
  EXPECT_EQ(prim_data->prim_crystal_point_group.size(), 16);
  EXPECT_EQ(prim_data->prim_sym_invariant_displacement_modes.has_value(), true);
  auto const &modes = *prim_data->prim_sym_invariant_displacement_modes;
  EXPECT_EQ(modes.size(), 1);
  {
    Eigen::MatrixXd expected(3, 3);
    expected << 0., 0., 0.,       //
        0., 0., 0.,               //
        0., 0.707107, -0.707107;  //
    EXPECT_TRUE(almost_equal(modes[0], expected));
  }
}
