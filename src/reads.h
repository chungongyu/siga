#ifndef reads_h_
#define reads_h_

#include <iostream>
#include <vector>

namespace PairEnd {
    std::string basename(const std::string& name);
    std::string id(const std::string& name);
    std::string index(const std::string& name);
};

struct ReadInfo {
    ReadInfo() : length(0) {
    }
    ReadInfo(const std::string& name, size_t length) : name(name), length(length) {
    }
    std::string name;
    size_t length;
};

typedef std::vector< ReadInfo > ReadInfoList;

std::istream& operator>>(std::istream& stream, ReadInfoList& infos);

#endif // reads_h_
