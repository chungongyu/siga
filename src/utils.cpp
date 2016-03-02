#include "utils.h"

#include <cassert>

namespace PairEnd {
    std::string basename(const std::string& name) {
        assert(!name.empty());
        size_t pos = name.find_last_of('/');
        if (pos == std::string::npos) {
            return name;
        }
        return name.substr(0, pos);
    };

    std::string id(const std::string& name) {
        assert(!name.empty());
        std::string pid(name);

        size_t li = name.length() - 1;
        char last = name[li];

        if(last == 'A')
            pid[li] = 'B';
        else if(last == 'B')
            pid[li] = 'A';
        else if(last == '1')
            pid[li] = '2';
        else if(last == '2')
            pid[li] = '1';
        else if(last == 'f')
            pid[li] = 'r';
        else if(last == 'r')
            pid[li] = 'f';
        else
            pid = "";
        return pid;
    };

    std::string index(const std::string& name) {
        return name;
    };
};
