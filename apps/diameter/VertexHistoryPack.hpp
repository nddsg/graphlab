#include <cstdint>
#include <vector>
#include <iostream>

#define VHPSETWIDTH (32)
#define VHPWIDTHSHIFT (5)
#define VHPSETBITMASK (31)
#define VHPSOLIDMASK (0xFFFFFFFF)

typedef uint32_t VHPPACKTYPE;

struct VertexHistoryPack {
  std::vector<uint32_t> v;
//  VertexHistoryPack() : v() {}
  VertexHistoryPack() : v() {
  }
  // inline VertexHistoryPack& operator+=(VertexHistoryPack& rvalue) {
    // v += rvalue.v;
    // return *this;
  // }
  bool get (uint32_t pos) const{
    uint32_t actual = pos >> VHPWIDTHSHIFT;
    // if (v.size() <= actual) {
      // return 0;
    // }
    VHPPACKTYPE packed = v[actual];
    VHPPACKTYPE internal = 1 << (pos & VHPSETBITMASK);
    return (packed & internal) != 0;
  }
  void set(uint32_t pos, bool value) {
    uint32_t actual = pos >> VHPWIDTHSHIFT;
    // if (v.size() <= actual) {
      // v.resize(actual + 1, 0);
    // }
    VHPPACKTYPE internal = 1 << (pos & VHPSETBITMASK);
    VHPPACKTYPE packed = v[actual];
    if (value) {
      v[actual] = packed | internal;
    }
    else {
      v[actual] = packed & ~internal; 
    }
  }
  void clear() {
    for (unsigned int i = 0; i < v.size(); i++) {
      v[i] = 0;
    }
  }
  void resize(uint32_t size) {
    v.resize(size, 0);
  }
  void dump() const{
    std::cout << "RAW:\n";
    for (uint32_t i = 0; i < v.size(); i++) {
      std::cout << "\t[" << i << "] => " << v[i] << "\n";
    }
    std::cout << "EXPANDED:\n";
    for (uint32_t i = 0; i < v.size() * VHPSETWIDTH; i++) {
      if (get(i)) {
        std::cout << "\t[" << i << "] => 1\n";
      }
    }
  }
#ifdef GRAPHLAB_SERIALIZE_HPP
  void save(graphlab::oarchive& oarc) const {
     oarc << v;
  }
  void load(graphlab::iarchive& iarc) {
     iarc >> v;
  }
#endif
};
