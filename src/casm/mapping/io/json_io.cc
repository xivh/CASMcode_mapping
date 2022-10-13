#include "casm/mapping/io/json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/mapping/AtomMapping.hh"
#include "casm/mapping/LatticeMapping.hh"
#include "casm/mapping/StructureMapping.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {

// LatticeMapping

jsonParser &to_json(mapping::LatticeMapping const &m, jsonParser &json) {
  json["deformation_gradient"] = m.deformation_gradient;
  json["transformation_matrix_to_supercell"] =
      lround(m.transformation_matrix_to_super).cast<int>();
  json["reorientation"] = lround(m.reorientation).cast<int>();
  json["right_stretch"] = m.right_stretch;
  json["isometry"] = m.isometry;
  json["left_stretch"] = m.left_stretch;

  return json;
}

void from_json(mapping::LatticeMapping &m, jsonParser const &json) {
  m = jsonConstructor<mapping::LatticeMapping>::from_json(json);
}

mapping::LatticeMapping jsonConstructor<mapping::LatticeMapping>::from_json(
    jsonParser const &json) {
  Eigen::Matrix3d deformation_gradient;
  Eigen::Matrix3l transformation_matrix_to_super;
  Eigen::Matrix3l reorientation;
  CASM::from_json(deformation_gradient, json["deformation_gradient"]);
  CASM::from_json(transformation_matrix_to_super,
                  json["transformation_matrix_to_supercell"]);
  CASM::from_json(reorientation, json["reorientation"]);

  return mapping::LatticeMapping(deformation_gradient,
                                 transformation_matrix_to_super.cast<double>(),
                                 reorientation.cast<double>());
}

// ScoredLatticeMapping

jsonParser &to_json(mapping::ScoredLatticeMapping const &value,
                    jsonParser &json) {
  json["lattice_cost"] = value.lattice_cost;
  to_json(static_cast<mapping::LatticeMapping const &>(value), json);
  return json;
}

void from_json(mapping::ScoredLatticeMapping &m, jsonParser const &json) {
  m = jsonConstructor<mapping::ScoredLatticeMapping>::from_json(json);
}

mapping::ScoredLatticeMapping jsonConstructor<
    mapping::ScoredLatticeMapping>::from_json(jsonParser const &json) {
  return mapping::ScoredLatticeMapping(
      json["lattice_cost"].get<double>(),
      jsonConstructor<mapping::LatticeMapping>::from_json(json));
}

// LatticeMappingResults

jsonParser &to_json(mapping::LatticeMappingResults const &results,
                    jsonParser &json) {
  to_json(results.data, json);
  return json;
}

void from_json(mapping::LatticeMappingResults &results,
               jsonParser const &json) {
  from_json(results.data, json);
}

// AtomMapping

jsonParser &to_json(mapping::AtomMapping const &m, jsonParser &json) {
  json["displacement"] = m.displacement.transpose();
  json["permutation"] = m.permutation;
  to_json(m.translation, json["translation"], jsonParser::as_array());

  return json;
}

void from_json(mapping::AtomMapping &m, jsonParser const &json) {
  m = jsonConstructor<mapping::AtomMapping>::from_json(json);
}

mapping::AtomMapping jsonConstructor<mapping::AtomMapping>::from_json(
    jsonParser const &json) {
  Eigen::MatrixXd displacement_t;
  std::vector<Index> permutation;
  Eigen::Vector3d translation;
  CASM::from_json(displacement_t, json["displacement"]);
  CASM::from_json(permutation, json["permutation"]);
  CASM::from_json(translation, json["translation"]);

  return mapping::AtomMapping(displacement_t.transpose(), permutation,
                              translation);
}

// ScoredAtomMapping

jsonParser &to_json(mapping::ScoredAtomMapping const &value, jsonParser &json) {
  json["atom_cost"] = value.atom_cost;
  to_json(static_cast<mapping::AtomMapping const &>(value), json);
  return json;
}

void from_json(mapping::ScoredAtomMapping &m, jsonParser const &json) {
  m = jsonConstructor<mapping::ScoredAtomMapping>::from_json(json);
}

mapping::ScoredAtomMapping
jsonConstructor<mapping::ScoredAtomMapping>::from_json(jsonParser const &json) {
  return mapping::ScoredAtomMapping(
      json["atom_cost"].get<double>(),
      jsonConstructor<mapping::AtomMapping>::from_json(json));
}

// AtomMappingResults

jsonParser &to_json(mapping::AtomMappingResults const &results,
                    jsonParser &json) {
  to_json(results.data, json);
  return json;
}

void from_json(mapping::AtomMappingResults &results, jsonParser const &json) {
  from_json(results.data, json);
}

// StructureMapping

jsonParser &to_json(mapping::StructureMapping const &m, jsonParser &json) {
  to_json(m.lattice_mapping, json);
  to_json(m.atom_mapping, json);

  return json;
}

void from_json(mapping::StructureMapping &m, jsonParser const &json,
               std::shared_ptr<xtal::BasicStructure const> const &prim) {
  m = jsonConstructor<mapping::StructureMapping>::from_json(json, prim);
}

mapping::StructureMapping jsonConstructor<mapping::StructureMapping>::from_json(
    jsonParser const &json,
    std::shared_ptr<xtal::BasicStructure const> const &prim) {
  return mapping::StructureMapping(
      prim, jsonConstructor<mapping::LatticeMapping>::from_json(json),
      jsonConstructor<mapping::AtomMapping>::from_json(json));
}

// StructureMappingCost

jsonParser &to_json(mapping::StructureMappingCost const &cost,
                    jsonParser &json) {
  to_json(cost.lattice_cost, json["lattice_cost"]);
  to_json(cost.atom_cost, json["atom_cost"]);
  to_json(cost.total_cost, json["total_cost"]);

  return json;
}

void from_json(mapping::StructureMappingCost &cost, jsonParser const &json) {
  cost = jsonConstructor<mapping::StructureMappingCost>::from_json(json);
}

mapping::StructureMappingCost jsonConstructor<
    mapping::StructureMappingCost>::from_json(jsonParser const &json) {
  return mapping::StructureMappingCost(json["lattice_cost"].get<double>(),
                                       json["atom_cost"].get<double>(),
                                       json["total_cost"].get<double>());
}

// ScoredStructureMapping

jsonParser &to_json(mapping::ScoredStructureMapping const &value,
                    jsonParser &json) {
  to_json(value.lattice_cost, json["lattice_cost"]);
  to_json(value.atom_cost, json["atom_cost"]);
  to_json(value.total_cost, json["total_cost"]);
  to_json(static_cast<mapping::StructureMapping const &>(value), json);

  return json;
}

void from_json(mapping::ScoredStructureMapping &m, jsonParser const &json,
               std::shared_ptr<xtal::BasicStructure const> const &prim) {
  m = jsonConstructor<mapping::ScoredStructureMapping>::from_json(json, prim);
}

mapping::ScoredStructureMapping
jsonConstructor<mapping::ScoredStructureMapping>::from_json(
    jsonParser const &json,
    std::shared_ptr<xtal::BasicStructure const> const &prim) {
  return mapping::ScoredStructureMapping(
      json["lattice_cost"].get<double>(), json["atom_cost"].get<double>(),
      json["total_cost"].get<double>(),
      jsonConstructor<mapping::StructureMapping>::from_json(json, prim));
}

// StructureMappingResults

jsonParser &to_json(mapping::StructureMappingResults const &results,
                    jsonParser &json) {
  to_json(results.data, json);
  return json;
}

void from_json(mapping::StructureMappingResults &results,
               jsonParser const &json,
               std::shared_ptr<xtal::BasicStructure const> const &prim) {
  from_json(results.data, json, prim);
}

}  // namespace CASM
