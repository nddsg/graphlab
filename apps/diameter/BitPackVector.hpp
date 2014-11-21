#include <cstdint>
#include <map>
#include <iostream>

#include <graphlab/serialization/oarchive.hpp>

#define SETWIDTH (32)
#define WIDTHSHIFT (5)
#define SETBITMASK (31)
#define SOLIDMASK (0xFFFFFFFF)

typedef std::vector< uint32_t > mapvector;

class BitPackVector {
  bool defaultValue;
  protected:
  mapvector bitvector;
  public:
  class iterator {
    public:
    mutable uint32_t first;
    mutable bool second;

    private:
    mutable uint32_t internal_iterator;
    mutable mapvector * internal_vector;
    mutable int current_bit;
    void find_next() const {
      uint32_t internal_end = internal_vector->size();
      if (internal_vector->size() == 0) {
        internal_iterator = internal_end;
        current_bit = -1;
      }
      if (current_bit != -1) {
        current_bit++;
      }
      while (((current_bit == -1) && (internal_iterator != internal_end))
        || ((((1 << current_bit) & (*internal_vector)[internal_iterator]) == 0) 
          && (internal_iterator != internal_end))) {
        if (current_bit == -1) {
          while (((*internal_vector)[internal_iterator] == 0) && (internal_iterator != internal_end)) {
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
      first = (internal_iterator << WIDTHSHIFT) + current_bit;
      second = ((*internal_vector)[internal_iterator] & (1 << current_bit)) != 0;
    }

    public:
    iterator& operator++() {
      find_next();
      return *this;
    }
    iterator operator++(int total) {
      iterator tmp = *this;
      find_next();
      return tmp;
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
      *internal_vector = *i.internal_vector;
      return *this;
    }
    const iterator* begin() const {
      internal_iterator = 0;
      current_bit = -1;
      find_next();
      return this;
    }
    const iterator* end() const {
      internal_iterator = internal_vector->size();
      current_bit = -1;
      return this;
    }
    iterator(const mapvector & m) {
      *internal_vector = *const_cast<mapvector *>(&m);
      begin();
    }
  };

  BitPackVector() : defaultValue(false), bitvector() {}
  BitPackVector(bool defaultvalue) : defaultValue(defaultvalue), bitvector() {}
  ~BitPackVector(){}
  bool operator[] (uint32_t pos) const {
    uint32_t actual = pos >> WIDTHSHIFT;
    uint32_t packed = bitvector[actual];
    uint32_t internal = 1 << ((pos & SETBITMASK) - 1);
    return (packed & internal) != 0;
  }
  void set(uint32_t pos, bool value) {
    uint32_t actual = pos >> WIDTHSHIFT;
    uint32_t internal = 1 << (pos & SETBITMASK);
    uint32_t packed = bitvector[actual];
    if (value) {
      bitvector[actual] = packed | internal;
    }
    else {
      bitvector[actual] = packed & ~internal; 
    }
  }
  inline BitPackVector& operator+=(const BitPackVector& rvalue) {
    if (!rvalue.bitvector.empty()) {
      if (bitvector.empty()) bitvector = rvalue.bitvector;
      else {
        for(size_t t = 0; t < rvalue.bitvector.size(); ++t) bitvector[t] = bitvector[t] | rvalue.bitvector[t];
      }
    }
    return *this;
  }
  void reserve(uint32_t size) {
    this->bitvector.reserve(size);
  }
  void clear() {
    for (uint32_t i = 0; i < this->bitvector.size(); i++) {
      this->bitvector[i] = defaultValue;
    }
  }
  bool empty() const {
    return this->bitvector.empty();
  }
  void insert(uint32_t pos) {
    set(pos, true);
  }
  iterator begin() const {
    return *iterator(this->bitvector);
  }
  iterator end() const {
    iterator temp(this->bitvector);
    return *temp.end();
  }

  void dump() {
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
