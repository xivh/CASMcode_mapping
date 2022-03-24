#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/mapping/SearchData.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;
using namespace CASM::xtal;
using namespace CASM::mapping;

namespace CASM {
namespace mapping {
namespace mapping_impl {

// Defined in SearchData.cc
Eigen::MatrixXd make_supercell_site_coordinate_cart(
    xtal::UnitCellCoordIndexConverter const &unitcellcoord_index_converter,
    Eigen::MatrixXd const &prim_site_coordinate_cart,
    xtal::Lattice const &prim_lattice);

}  // namespace mapping_impl
}  // namespace mapping
}  // namespace CASM

namespace test {

/// \brief Modify r2 and atom_type (which correspond to F
/// (r1[i] + disp) = r2[i]) to be consistent with F (r1[i] + disp) = r2[perm[i]]
/// + trans
void permute_and_translate(Eigen::MatrixXd &r2,
                           std::vector<std::string> &atom_type,
                           std::vector<Index> const &perm,
                           Eigen::Vector3d const &trans) {
  // F (r1[i] + disp) = r2[perm[i]] + trans

  Eigen::MatrixXd r2_tmp = r2;
  for (Index i = 0; i < r2.cols(); ++i) {
    r2_tmp.col(perm[i]) = r2.col(i) - trans;
  }
  r2 = r2_tmp;

  std::vector<std::string> atom_type_tmp = atom_type;
  for (Index i = 0; i < atom_type.size(); ++i) {
    atom_type_tmp[perm[i]] = atom_type[i];
  }
  atom_type = atom_type_tmp;
}

// binary, BCC
std::shared_ptr<PrimSearchData const> make_prim_1() {
  // lattice vectors as cols
  Eigen::Matrix3d lat;
  lat << 4.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 4.0;

  BasicStructure prim{Lattice{lat}};

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  prim.push_back(
      Site(Coordinate(Eigen::Vector3d(0., 0., 0.), prim.lattice(), FRAC),
           std::vector<Molecule>{A, B}));
  prim.push_back(
      Site(Coordinate(Eigen::Vector3d(0.5, 0.5, 0.5), prim.lattice(), FRAC),
           std::vector<Molecule>{A, B}));

  return std::make_shared<PrimSearchData const>(
      std::make_shared<BasicStructure const>(prim));
}

struct SearchTestData {
  std::shared_ptr<PrimSearchData const> prim_data;

  // F * L1 * T * N = L2

  Eigen::Matrix3d L1;
  Eigen::Matrix3d F;
  Eigen::Matrix3d T;
  Eigen::Matrix3d N;
  Eigen::Matrix3d L2;
  LatticeMapping lattice_mapping;

  // F (r1_supercell[i] + disp) = r2[perm[i]] + trans

  Eigen::MatrixXd r1;
  Eigen::MatrixXd r1_supercell;
  Eigen::MatrixXd disp;
  std::vector<std::string> atom_type;
  Eigen::MatrixXd r2;
  std::vector<Index> perm;
  Eigen::Vector3d trans;

  std::shared_ptr<StructureSearchData const> structure_data;

  SearchTestData(std::shared_ptr<PrimSearchData const> _prim_data,
                 Eigen::Matrix3d const &_F, Eigen::Matrix3d const &_T,
                 Eigen::Matrix3d const &_N, Eigen::MatrixXd const &_disp,
                 std::vector<std::string> _atom_type, std::vector<Index> _perm,
                 Eigen::Vector3d const &_trans)
      : prim_data(std::move(_prim_data)),
        L1(prim_data->prim_lattice.lat_column_mat()),
        F(_F),
        T(_T),
        N(_N),
        L2(F * L1 * T * N),
        lattice_mapping(F, T, N),
        r1(prim_data->prim_site_coordinate_cart),
        r1_supercell(mapping_impl::make_supercell_site_coordinate_cart(
            UnitCellCoordIndexConverter(lround(T * N), r1.cols()), r1,
            prim_data->prim_lattice)),
        disp(_disp),
        atom_type(std::move(_atom_type)),
        r2(F * (r1_supercell + disp)),
        perm(std::move(_perm)),
        trans(_trans) {
    // permute and translate
    test::permute_and_translate(r2, atom_type, perm, trans);
    structure_data =
        std::make_shared<StructureSearchData const>(Lattice(L2), r2, atom_type);
  }
};

}  // namespace test

// Test LatticeMappingSearchData construction using a structure
// generated from the prim unit cell
TEST(LatticeMappingSearchDataTest, Test1) {
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

  test::SearchTestData d(test::make_prim_1(), F, T, N, disp, atom_type, perm,
                         trans);

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

// Test LatticeMappingSearchData constructiong using a structure
// generated from a supercell of the prim
// (T = [[1, 0, 0], [0, 1, 0], [0, 0, 2]])
TEST(LatticeMappingSearchDataTest, Test2) {
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

  test::SearchTestData d(test::make_prim_1(), F, T, N, disp, atom_type, perm,
                         trans);

  auto lattice_mapping_data = std::make_shared<LatticeMappingSearchData const>(
      d.prim_data, d.structure_data, d.lattice_mapping);

  EXPECT_TRUE(almost_equal(lattice_mapping_data->transformation_matrix_to_super,
                           lround(T * N)));
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
