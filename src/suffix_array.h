#ifndef suffix_array_h_
#define suffix_array_h_

#include "kseq.h"

#include <iostream>
#include <vector>

class SuffixArray {
public:
    struct Elem {
        Elem() : i(-1), j(-1) {
        }
        Elem(size_t i, size_t j) : i(i), j(j) {
        }
        bool empty() const {
            return i == -1 || j == -1;
        }
        operator bool() const {
            return not empty();
        }
        size_t i;
        size_t j;
    };
    typedef std::vector< Elem > ElemList;

    SuffixArray() : _strings(0) {
    }
private:
    friend std::ostream& operator<<(std::ostream& stream, const SuffixArray& sa);
    friend std::istream& operator>>(std::istream& stream, SuffixArray& sa);
    friend class SAReader;
    friend class SAWriter;

    ElemList _elems;
    size_t _strings;
};

#endif // suffix_array_h_
