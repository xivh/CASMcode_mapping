#ifndef CASM_mapping_StrucMapping_json_io
#define CASM_mapping_StrucMapping_json_io

namespace CASM {

namespace mapping_impl {
struct MappingNode;
}

class jsonParser;

jsonParser &to_json(mapping_impl::MappingNode const &mapping_node,
                    jsonParser &json);

}  // namespace CASM

#endif
