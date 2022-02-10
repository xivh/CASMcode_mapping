#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

xtal::Lattice make_lattice(Eigen::Ref<Eigen::Matrix3d const> const &L) {
  return xtal::Lattice{L};
}

// {"axis_names", "basis"}
typedef std::pair<std::vector<std::string>, Eigen::MatrixXd> DoFSetInput;

/// \brief Construct DoFSetInput
///
/// \param axis_names DoFSet axis names. Size equals number of columns in basis.
/// \param basis Basis vectors, as columns. `x_standard = basis * x_prim`
///
DoFSetInput make_dofsetinput(std::vector<std::string> const &axis_names,
                             Eigen::MatrixXd const &basis) {
  return DoFSetInput(axis_names, basis);
}

std::map<std::string, xtal::SpeciesProperty> make_species_properties(
    std::map<std::string, Eigen::MatrixXd> species_properties) {
  std::map<std::string, xtal::SpeciesProperty> result;
  for (auto const &pair : species_properties) {
    result.emplace(pair.first, xtal::SpeciesProperty{AnisoValTraits(pair.first),
                                                     pair.second});
  }
  return result;
}

xtal::AtomPosition make_atom_position(
    std::string name, Eigen::Vector3d pos,
    std::map<std::string, xtal::SpeciesProperty> properties = {}) {
  xtal::AtomPosition atom(pos, name);
  atom.set_properties(properties);
  return atom;
}

std::map<std::string, Eigen::MatrixXd> get_atom_position_properties(
    xtal::AtomPosition const &atom) {
  std::map<std::string, Eigen::MatrixXd> result;
  for (auto const &pair : atom.properties()) {
    result.emplace(pair.first, pair.second.value());
  }
  return result;
}

xtal::Molecule make_molecule(
    std::string name, std::vector<xtal::AtomPosition> atoms = {},
    bool divisible = false,
    std::map<std::string, xtal::SpeciesProperty> properties = {}) {
  xtal::Molecule mol(name, atoms, divisible);
  mol.set_properties(properties);
  return mol;
}

std::map<std::string, Eigen::MatrixXd> get_molecule_properties(
    xtal::Molecule const &mol) {
  std::map<std::string, Eigen::MatrixXd> result;
  for (auto const &pair : mol.properties()) {
    result.emplace(pair.first, pair.second.value());
  }
  return result;
}

/// \brief Construct xtal::BasicStructure from JSON string
xtal::BasicStructure basicstructure_from_json(std::string const &prim_json_str,
                                              double xtal_tol) {
  jsonParser json{prim_json_str};
  ParsingDictionary<AnisoValTraits> const *aniso_val_dict = nullptr;
  return read_prim(json, xtal_tol, aniso_val_dict);
}

/// \brief Format xtal::BasicStructure as JSON string
std::string basicstructure_to_json(xtal::BasicStructure const &prim) {
  jsonParser json;
  write_prim(prim, json, FRAC);
  std::stringstream ss;
  ss << json;
  return ss.str();
}

xtal::BasicStructure make_basicstructure(
    Eigen::Ref<Eigen::Matrix3d const> const &L,
    Eigen::Ref<Eigen::MatrixXd const> const &B_frac,
    std::vector<std::vector<std::string>> const &occ_dof,
    std::vector<std::map<std::string, DoFSetInput>> const &local_dof,
    std::map<std::string, DoFSetInput> const &global_dof,
    std::map<std::string, xtal::Molecule> const &molecules) {
  // construct prim
  xtal::BasicStructure prim{L};

  // set basis sites
  for (Index b = 0; b < B_frac.cols(); ++b) {
    xtal::Coordinate coord{B_frac.col(b), prim.lattice(), FRAC};
    std::vector<xtal::Molecule> site_occ;
    for (std::string label : occ_dof[b]) {
      site_occ.push_back(molecules.at(label));
    }
    std::vector<xtal::SiteDoFSet> site_dofsets;
    for (auto const &pair : local_dof[b]) {
      std::string const &dofname = pair.first;
      DoFSetInput const &dofsetinfo = pair.second;
      site_dofsets.emplace_back(AnisoValTraits(dofname), dofsetinfo.first,
                                dofsetinfo.second,
                                std::unordered_set<std::string>{});
    }
    xtal::Site site{coord, site_occ, site_dofsets};
    prim.push_back(site, FRAC);
  }
  prim.set_unique_names(occ_dof);

  // set global dof
  std::vector<xtal::DoFSet> global_dofsets;
  for (auto const &pair : global_dof) {
    std::string const &dofname = pair.first;
    DoFSetInput const &dofsetinfo = pair.second;
    global_dofsets.emplace_back(AnisoValTraits(dofname), dofsetinfo.first,
                                dofsetinfo.second);
  }
  prim.set_global_dofs(global_dofsets);

  return prim;
}

// void get_basicstructure(
//     xtal::BasicStructure const &prim,
//     Eigen::Ref<Eigen::Matrix3d> &L,
//     Eigen::Ref<Eigen::MatrixXd> &B_frac,
//     std::vector<std::vector<std::string>> &occ_dof,
//     std::vector<std::map<std::string, DoFSetInput>> &local_dof,
//     std::map<std::string, DoFSetInput> &global_dof,
//     std::map<std::string, xtal::Molecule> &molecules) {
//
//   L = prim.lattice().lat_column_mat();
//
//   B_frac.resize(3, prim.basis().size());
//   Index b = 0;
//   for (auto const &site : prim.basis()) {
//     B_frac.col(b) = site.const_frac();
//     ++b;
//   }
//
//   occ_dof = xtal::allowed_molecule_names(prim);
//
//   local_dof.clear();
//   Index b = 0;
//   for (auto const &site : prim.basis()) {
//     std::map<std::string, DoFSetInput> site_dof;
//     for (auto const &pair : site.dofs()) {
//       std::string const &dofname = pair.first;
//       xtal::SiteDoFSet const &dofset = pair.second;
//       site_dof.emplace(
//           pair.first,
//           DoFSetInput(
//               dofset.component_names(),
//               dofset.basis()));
//     }
//     local_dof.push_back(site_dof);
//     ++b;
//   }
//
//   global_dof.clear();
//   for (auto const &pair : prim.global_dofs()) {
//     std::string const &dofname = pair.first;
//     xtal::DoFSet const &dofset = pair.second;
//     site_dof.emplace(
//         pair.first,
//         DoFSetInput(
//             dofset.component_names(),
//             dofset.basis()));
//   }
//
//   molecules.clear();
//   std::vector<std::vector<std::string>> mol_names = prim.unique_names();
//   if (mol_names.empty()) {
//     mol_names = xtal::allowed_molecule_unique_names(prim);
//   }
//   for (Index b = 0; b < mol_names.size(); ++b) {
//     for (Index i = 0; i < mol_names[b].size(); ++i) {
//       std::string const &name =  mol_names[b][i];
//       if (!molecules.count(name)) {
//         molecules.emplace(name, prim.basis()[b].occupant_dof()[i]);
//       }
//     }
//   }
//
// }

}  // namespace CASMpy

PYBIND11_MODULE(_xtal, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        CASMcode_crystallography bindings module
        ----------------------------------------
        .. currentmodule:: _xtal
        .. autosummary::
           :toctree: _generate
           MutablePrim
    )pbdoc";

  py::class_<xtal::Lattice>(m, "Lattice")
      .def(py::init<Eigen::Matrix3d const &, double, bool>(),
           py::arg("lat_column_mat"), py::arg("xtal_tol") = TOL,
           py::arg("force") = false)
      .def("lat_column_mat", &xtal::Lattice::lat_column_mat,
           py::return_value_policy::reference_internal);

  m.def("make_lattice", &make_lattice, "Construct a Lattice");

  py::class_<xtal::AtomPosition>(m, "AtomComponent")
      .def(py::init(&make_atom_position))
      .def("name", &xtal::AtomPosition::name)
      .def("coordinate", &xtal::AtomPosition::cart)
      .def("properties", &get_atom_position_properties);

  m.def("make_atom_component", &make_atom_position,
        "Construct an AtomComponent object");

  py::class_<xtal::Molecule>(m, "Molecule")
      .def(py::init(&make_molecule))
      .def("name", &xtal::Molecule::name)
      .def("is_vacancy", &xtal::Molecule::is_vacancy)
      .def("is_atomic", &xtal::Molecule::is_atomic)
      .def("is_divisible", &xtal::Molecule::is_divisible)
      .def("atoms", &xtal::Molecule::atoms)
      .def("atom", &xtal::Molecule::atom)
      .def("size", &xtal::Molecule::size)
      .def("properties", &get_molecule_properties);

  m.def("make_molecule", &make_molecule, "Construct a Molecule object");
  m.def("make_vacancy", &xtal::Molecule::make_vacancy,
        "Construct a vacancy (Molecule object)");
  m.def("make_atom", &xtal::Molecule::make_atom,
        "Construct a atom (Molecule object)");

  py::class_<DoFSetInput>(m, "DoFSetInput")
      .def(py::init(&make_dofsetinput))
      .def_readonly("axis_names", &DoFSetInput::first)
      .def_readonly("basis", &DoFSetInput::second);

  m.def("make_dofsetinput", &make_dofsetinput,
        "Construct a DoFSetInput object");

  py::class_<xtal::BasicStructure>(m, "Prim")
      .def(py::init(&make_basicstructure))
      .def("lattice", &xtal::BasicStructure::lattice);

  m.def("make_basicstructure", &make_basicstructure, "Construct a MutablePrim");
  m.def("basicstructure_from_json", &basicstructure_from_json,
        "Construct a MutablePrim from a JSON-formatted string");
  m.def("basicstructure_to_json", &basicstructure_to_json,
        "Construct a JSON-formatted string from a MutablePrim");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
