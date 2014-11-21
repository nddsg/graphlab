#include "VertexHistoryPack.hpp"


int main() {
  VertexHistoryPack a;
  a.resize(5);
  for (int i = 3; i < 8; i++) {
    a.set(i, true);
  }
  a.dump();
  return 0;
}
