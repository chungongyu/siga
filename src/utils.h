#ifndef utils_h_
#define utils_h_

#include <cassert>
#include <sstream>
#include <string>

#ifndef SIZEOF_BITS
#define SIZEOF_BITS(x) (8 * sizeof(x))
#endif

#ifndef SIZEOF_ARRAY
#define SIZEOF_ARRAY(x)  (sizeof(x) / sizeof(x[0]))
#endif

#ifndef IS_POWER_OF_2
#define IS_POWER_OF_2(x) ((x) & ((x) - 1)) == 0
#endif

// return the x % y given that y is a power of 2
#ifndef MOD_POWER_2
#define MOD_POWER_2(x, y) (x) & ((y) - 1)
#endif

#ifndef SAFE_DELETE
#define SAFE_DELETE(x) if (x != NULL) {delete x; x = NULL;}
#endif

#ifndef SAFE_DELETE_ARRAY
#define SAFE_DELETE_ARRAY(x) if (x != NULL) {delete[] x; x = NULL;}
#endif

namespace Utils {
  // random
  void srand();
  int rand();
  inline int rand(int start, int stop) {
    assert(stop > start);
    return rand() % (stop - start) + start;
  }
  inline int rand(int stop) {
    return rand(0, stop);
  }

  // streams
  std::istream* ifstream(const std::string& filename);
  std::ostream* ofstream(const std::string& filename);
  std::string stem(const std::string& filename);
};  // namespace Utils

#endif  // utils_h_
