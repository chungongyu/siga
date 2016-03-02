#ifndef utils_h_
#define utils_h_

#include <string>

#ifndef SIZEOF_BITS
#define SIZEOF_BITS(x) (8 * sizeof(x))
#endif

#ifndef SIZEOF_ARRAY
#define SIZEOF_ARRAY(x)  (sizeof(x) / sizeof(x[0]))
#endif

namespace PairEnd {
    std::string basename(const std::string& name);
    std::string id(const std::string& name);
    std::string index(const std::string& name);
};

#endif // utils_h_
