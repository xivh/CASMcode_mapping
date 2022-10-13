#ifndef CASM_mapping_unittest_SearchTestData
#define CASM_mapping_unittest_SearchTestData

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/mapping/SearchData.hh"
#include "casm/misc/CASM_Eigen_math.hh"

using namespace CASM;

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
inline void permute_and_translate(Eigen::MatrixXd &r2,
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

struct SearchTestData {
  std::shared_ptr<mapping::PrimSearchData const> prim_data;

  // F * L1 * T * N = L2

  Eigen::Matrix3d L1;
  Eigen::Matrix3d F;
  Eigen::Matrix3d T;
  Eigen::Matrix3d N;
  Eigen::Matrix3d L2;
  mapping::LatticeMapping lattice_mapping;

  // F (r1_supercell[i] + disp) = r2[perm[i]] + trans
  // structure1_supercell_atom_type[i] = structure2_atom_type[perm[i]]

  Eigen::MatrixXd r1;
  Eigen::MatrixXd r1_supercell;
  Eigen::MatrixXd disp;
  std::vector<std::string> structure1_supercell_atom_type;
  std::vector<std::string> structure2_atom_type;
  Eigen::MatrixXd r2;
  std::vector<Index> perm;
  Eigen::Vector3d trans;

  std::shared_ptr<mapping::StructureSearchData const> structure_data;

  SearchTestData(std::shared_ptr<mapping::PrimSearchData const> _prim_data,
                 Eigen::Matrix3d const &_F, Eigen::Matrix3d const &_T,
                 Eigen::Matrix3d const &_N, Eigen::MatrixXd const &_disp,
                 std::vector<std::string> _structure1_supercell_atom_type,
                 std::vector<Index> _perm, Eigen::Vector3d const &_trans)
      : prim_data(std::move(_prim_data)),
        L1(prim_data->prim_lattice.lat_column_mat()),
        F(_F),
        T(_T),
        N(_N),
        L2(F * L1 * T * N),
        lattice_mapping(F, T, N),
        r1(prim_data->prim_site_coordinate_cart),
        r1_supercell(mapping::mapping_impl::make_supercell_site_coordinate_cart(
            xtal::UnitCellCoordIndexConverter(lround(T * N), r1.cols()), r1,
            prim_data->prim_lattice)),
        disp(_disp),
        structure1_supercell_atom_type(
            std::move(_structure1_supercell_atom_type)),
        structure2_atom_type(structure1_supercell_atom_type),
        r2(F * (r1_supercell + disp)),
        perm(std::move(_perm)),
        trans(_trans) {
    // permute and translate
    permute_and_translate(r2, structure2_atom_type, perm, trans);
    structure_data = std::make_shared<mapping::StructureSearchData const>(
        xtal::Lattice(L2), r2, structure2_atom_type);
  }
};

// binary, primitive BCC
inline std::shared_ptr<mapping::PrimSearchData const>
make_search_prim_binary_BCC(double a) {
  using namespace xtal;
  // lattice vectors as cols
  Eigen::Matrix3d lat;
  lat << -a / 2., a / 2., a / 2.,  //
      a / 2., -a / 2., a / 2.,     //
      a / 2., a / 2., -a / 2.;     //

  BasicStructure prim{Lattice{lat}};

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  prim.push_back(
      Site(Coordinate(Eigen::Vector3d(0., 0., 0.), prim.lattice(), FRAC),
           std::vector<Molecule>{A, B}));

  return std::make_shared<mapping::PrimSearchData const>(
      std::make_shared<BasicStructure const>(prim));
}

// binary + Va, primitive BCC
inline std::shared_ptr<mapping::PrimSearchData const>
make_search_prim_binary_vacancy_BCC(double a) {
  using namespace xtal;
  // lattice vectors as cols
  Eigen::Matrix3d lat;
  lat << -a / 2., a / 2., a / 2.,  //
      a / 2., -a / 2., a / 2.,     //
      a / 2., a / 2., -a / 2.;     //

  BasicStructure prim{Lattice{lat}};

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");
  Molecule Va = Molecule::make_vacancy();

  prim.push_back(
      Site(Coordinate(Eigen::Vector3d(0., 0., 0.), prim.lattice(), FRAC),
           std::vector<Molecule>{A, B, Va}));

  return std::make_shared<mapping::PrimSearchData const>(
      std::make_shared<BasicStructure const>(prim));
}

// binary, primitive FCC
inline std::shared_ptr<mapping::PrimSearchData const>
make_search_prim_binary_FCC(double a) {
  using namespace xtal;
  // lattice vectors as cols
  Eigen::Matrix3d lat;
  lat << 0., a / 2., a / 2.,  //
      a / 2., 0., a / 2.,     //
      a / 2., a / 2., 0.;     //

  BasicStructure prim{Lattice{lat}};

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  prim.push_back(
      Site(Coordinate(Eigen::Vector3d(0., 0., 0.), prim.lattice(), FRAC),
           std::vector<Molecule>{A, B}));

  return std::make_shared<mapping::PrimSearchData const>(
      std::make_shared<BasicStructure const>(prim));
}

// binary, primitive HCP
inline std::shared_ptr<mapping::PrimSearchData const>
make_search_prim_binary_HCP(double a, double c) {
  using namespace xtal;
  // lattice vectors as rows
  Eigen::Matrix3d lat;
  lat << a, 0., 0.,                    // a
      -a / 2., a * sqrt(3.) / 2., 0.,  // a
      0.0, 0.0, c;                     // c

  BasicStructure prim{Lattice{lat.transpose()}};

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  prim.push_back(
      Site(Coordinate(Eigen::Vector3d(0., 0., 0.), prim.lattice(), FRAC),
           std::vector<Molecule>{A, B}));
  prim.push_back(Site(Coordinate(Eigen::Vector3d(2. / 3., 1. / 3., 1. / 2.),
                                 prim.lattice(), FRAC),
                      std::vector<Molecule>{A, B}));

  return std::make_shared<mapping::PrimSearchData const>(
      std::make_shared<BasicStructure const>(prim));
}

// binary, conventional BCC
inline std::shared_ptr<mapping::PrimSearchData const>
make_search_prim_binary_conventional_BCC(double a) {
  using namespace xtal;
  // lattice vectors as cols
  Eigen::Matrix3d lat;
  lat << a, 0.0, 0.0, 0.0, a, 0.0, 0.0, 0.0, a;

  BasicStructure prim{Lattice{lat}};

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  prim.push_back(
      Site(Coordinate(Eigen::Vector3d(0., 0., 0.), prim.lattice(), FRAC),
           std::vector<Molecule>{A, B}));
  prim.push_back(
      Site(Coordinate(Eigen::Vector3d(0.5, 0.5, 0.5), prim.lattice(), FRAC),
           std::vector<Molecule>{A, B}));

  return std::make_shared<mapping::PrimSearchData const>(
      std::make_shared<BasicStructure const>(prim));
}

}  // namespace test

#endif
