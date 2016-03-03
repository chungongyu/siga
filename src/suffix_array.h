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
        bool full() const {
            return j == 0;
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
    SuffixArray(size_t strings, size_t suffixes) : _strings(strings), _elems(suffixes) {
    }

    size_t size() const {
        return _elems.size();
    }
    const Elem& operator[](size_t i) const {
        return _elems[i];
    }
    Elem& operator[](size_t i) {
        return _elems[i];
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
