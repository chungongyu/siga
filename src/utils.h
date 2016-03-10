#ifndef utils_h_
#define utils_h_

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

namespace PairEnd {
    std::string basename(const std::string& name);
    std::string id(const std::string& name);
    std::string index(const std::string& name);
};

#endif // utils_h_
