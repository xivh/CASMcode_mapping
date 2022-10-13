#ifndef CASM_mapping_json_io
#define CASM_mapping_json_io

#include <memory>

namespace CASM {

class jsonParser;
template <typename T>
struct jsonConstructor;

namespace xtal {

class BasicStructure;
}

namespace mapping {

struct LatticeMapping;
struct ScoredLatticeMapping;
struct LatticeMappingResults;
struct AtomMapping;
struct ScoredAtomMapping;
struct AtomMappingResults;
struct StructureMapping;
struct StructureMappingCost;
struct ScoredStructureMapping;
struct StructureMappingResults;
}  // namespace mapping

// LatticeMapping

jsonParser &to_json(mapping::LatticeMapping const &m, jsonParser &json);

void from_json(mapping::LatticeMapping &m, jsonParser const &json);

template <>
struct jsonConstructor<mapping::LatticeMapping> {
  static mapping::LatticeMapping from_json(jsonParser const &json);
};

// ScoredLatticeMapping

jsonParser &to_json(mapping::ScoredLatticeMapping const &value,
                    jsonParser &json);

void from_json(mapping::ScoredLatticeMapping &m, jsonParser const &json);

template <>
struct jsonConstructor<mapping::ScoredLatticeMapping> {
  static mapping::ScoredLatticeMapping from_json(jsonParser const &json);
};

// LatticeMappingResults

jsonParser &to_json(mapping::LatticeMappingResults const &results,
                    jsonParser &json);

void from_json(mapping::LatticeMappingResults &results, jsonParser const &json);

// AtomMapping

jsonParser &to_json(mapping::AtomMapping const &m, jsonParser &json);

void from_json(mapping::AtomMapping &m, jsonParser const &json);

template <>
struct jsonConstructor<mapping::AtomMapping> {
  static mapping::AtomMapping from_json(jsonParser const &json);
};

// ScoredAtomMapping

jsonParser &to_json(mapping::ScoredAtomMapping const &value, jsonParser &json);

void from_json(mapping::ScoredAtomMapping &m, jsonParser const &json);

template <>
struct jsonConstructor<mapping::ScoredAtomMapping> {
  static mapping::ScoredAtomMapping from_json(jsonParser const &json);
};

// AtomMappingResults

jsonParser &to_json(mapping::AtomMappingResults const &results,
                    jsonParser &json);

void from_json(mapping::AtomMappingResults &results, jsonParser const &json);

// StructureMapping

jsonParser &to_json(mapping::StructureMapping const &m, jsonParser &json);

void from_json(mapping::StructureMapping &m, jsonParser const &json,
               std::shared_ptr<xtal::BasicStructure const> const &prim);

template <>
struct jsonConstructor<mapping::StructureMapping> {
  static mapping::StructureMapping from_json(
      jsonParser const &json,
      std::shared_ptr<xtal::BasicStructure const> const &prim);
};

// StructureMappingCost

jsonParser &to_json(mapping::StructureMappingCost const &cost,
                    jsonParser &json);

void from_json(mapping::StructureMappingCost &cost, jsonParser const &json);

template <>
struct jsonConstructor<mapping::StructureMappingCost> {
  static mapping::StructureMappingCost from_json(jsonParser const &json);
};

// ScoredStructureMapping

jsonParser &to_json(mapping::ScoredStructureMapping const &value,
                    jsonParser &json);

void from_json(mapping::ScoredStructureMapping &m, jsonParser const &json,
               std::shared_ptr<xtal::BasicStructure const> const &prim);

template <>
struct jsonConstructor<mapping::ScoredStructureMapping> {
  static mapping::ScoredStructureMapping from_json(
      jsonParser const &json,
      std::shared_ptr<xtal::BasicStructure const> const &prim);
};

// StructureMappingResults

jsonParser &to_json(mapping::StructureMappingResults const &results,
                    jsonParser &json);

void from_json(mapping::StructureMappingResults &results,
               jsonParser const &json,
               std::shared_ptr<xtal::BasicStructure const> const &prim);

}  // namespace CASM

#endif
