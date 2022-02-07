#ifndef CASM_mapping_StrucMapping_json_io
#define CASM_mapping_StrucMapping_json_io

namespace CASM {

namespace mapping_v1 {
struct MappingNode;
}

class jsonParser;

jsonParser &to_json(mapping_v1::MappingNode const &mapping_node,
                    jsonParser &json);

}  // namespace CASM

#endif
