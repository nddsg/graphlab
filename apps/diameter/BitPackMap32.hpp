#include <cstdint>
#include <map>
#include <iostream>

#define SETWIDTH (32)
#define WIDTHSHIFT (5)
#define SETBITMASK (31)
#define SOLIDMASK (0xFFFFFFFF)

typedef std::map< uint32_t, uint32_t > mapint;

class BitPackMap {
  bool defaultValue;
  protected:
  mapint bitmap;
  public:
  class iterator {
    public:
    mutable uint32_t first;
    mutable bool second;

    private:
    mutable mapint::iterator internal_iterator;
    mutable mapint * internal_map;
    mutable int current_bit;
    void find_next() const {
      mapint::iterator internal_end = internal_map->end();
      if (internal_map->size() == 0) {
        internal_iterator = internal_end;
        current_bit = -1;
      }
      if (current_bit != -1) {
        current_bit++;
      }
      if (current_bit >= SETWIDTH) {
        current_bit = -1;
        internal_iterator++;
      }
      while (((current_bit == -1) && (internal_iterator != internal_end))
        || ((internal_iterator != internal_end) && (((1 << current_bit) & internal_iterator->second) == 0))) {
        if (current_bit == -1) {
          while ((internal_iterator->second == 0) && (internal_iterator != internal_end)) {
            // Skip blocks which have no set members
            internal_iterator++;
          }
        }
        current_bit++;
        if (current_bit >= SETWIDTH) {
          current_bit = -1;
          internal_iterator++;
        }
      }
      first = (internal_iterator->first << WIDTHSHIFT) + current_bit;
      second = (internal_iterator->second & (1 << current_bit)) != 0;
    }

    public:
    iterator& operator++() {
      find_next();
      return *this;
    }
    iterator operator++(int total) {
      find_next();
      return *this;
    }
    iterator* operator->() {
      return this;
    }
    iterator& operator*() {
      return *this;
    }
    bool operator==(iterator const & rhs) const {
      return (internal_iterator == rhs.internal_iterator) && (current_bit == rhs.current_bit);
    }
    bool operator!=(iterator const & rhs) const {
      return !(*this == rhs);
    }
    iterator& operator=(const iterator & i) {
      internal_iterator = i.internal_iterator;
      current_bit = i.current_bit;
      first = i.first;
      second = i.second;
      internal_map = i.internal_map;
      return *this;
    }
    const iterator* begin() const {
      internal_iterator = internal_map->begin();
      current_bit = -1;
      find_next();
      return this;
    }
    const iterator* end() const {
      internal_iterator = internal_map->end();
      current_bit = -1;
      return this;
    }
    iterator(const mapint & m) {
      internal_map = const_cast<mapint *>(&m);
      begin();
    }
  };

  BitPackMap() : defaultValue(false), bitmap() {}
  BitPackMap(bool defaultvalue) : defaultValue(defaultvalue), bitmap() {}
  ~BitPackMap(){}
  bool operator[] (uint32_t pos) const {
    uint32_t actual = pos >> WIDTHSHIFT;
    auto i = bitmap.find(actual);
    if (i == bitmap.end()) {
      return defaultValue;
    }
    uint32_t packed = i->second;
    uint32_t internal = 1 << (pos & SETBITMASK);
    return (packed & internal) != 0;
  }
  void set(uint32_t pos, bool value) {
    uint32_t actual = pos >> WIDTHSHIFT;
    uint32_t internal = 1 << (pos & SETBITMASK);
    uint32_t packed;
    auto i = bitmap.find(actual);
    if (i == bitmap.end()) {
      packed = defaultValue ? SOLIDMASK : 0;
    }
    else {
      packed = i->second;
    }
    if (value) {
      bitmap[actual] = packed | internal;
    }
    else {
      bitmap[actual] = packed & ~internal; 
    }
  }
  inline BitPackMap& operator+=(const BitPackMap& rvalue) {
    for (auto & iterator : rvalue.bitmap) {
      uint32_t temp;
      auto i = bitmap.find(iterator.first);
      if (i == bitmap.end()) {
        temp = defaultValue ? SOLIDMASK : 0;
      }
      else {
        temp = i->second;
      }
      this->bitmap[iterator.first] = temp | iterator.second;
    }
    return *this;
  }
  void clear() {
    this->bitmap.clear();
  }
  bool empty() const {
    return this->bitmap.empty();
  }
  void insert(uint32_t pos) {
    set(pos, true);
  }
  iterator begin() const {
    return *iterator(this->bitmap);
  }
  iterator end() const {
    iterator temp(this->bitmap);
    return *temp.end();
  }
  uint32_t size() const {
    uint32_t size = 0;
    for (iterator iterator = this->begin(); iterator != this->end(); iterator++) {
      size++;
    }
    return size;
  }

  void dump() const{
    for (iterator iterator = this->begin(); iterator != this->end(); iterator++) {
      std::cout << "\t[" << (iterator.first) << "] => " << iterator.second << std::endl;
    }
  }
#ifdef GRAPHLAB_SERIALIZE_HPP
  void save(graphlab::oarchive& oarc) const {
    //oarc << ;
  }
  void load(graphlab::iarchive& iarc) {
    //iarc >> ;
  }
#endif
};
