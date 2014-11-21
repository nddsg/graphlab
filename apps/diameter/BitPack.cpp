#include "BitPackMap32.hpp"

int main() {
  BitPackMap a, b;
  // a.set(32, true);
  // b.set(30, true);
  // b.set(64, true);
  // a += b;
  for (int i = 30; i < 40; i++) {
    a.set(i, true);
  }
  for (int i = 20; i < 50; i++) {
    if (a[i]) {
      std::cout << "\t\t[" << i << "] => " << 1 << "\n";
    }
  }
  // a.set(32, false);
  a.dump();
  return 0;
}
