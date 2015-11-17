#include "bwt.h"
#include "suffix_array.h"

BWT::BWT(const SuffixArray& sa, const DNASeqList& sequences) : _data(sa.size(), ' ') {
    for (size_t i = 0; i < sa.size(); ++i) {
        const SuffixArray::Elem& elem = sa[i];
        const DNASeq& read = sequences[elem.i];
        char c = (elem.j == 0) ? '$' : read.seq[elem.j - 1];
        _data[i] = c;
    }
}

std::ostream& operator<<(std::ostream& stream, const BWT& bwt) {
    stream << bwt._data;
    return stream;
}

std::istream& operator>>(std::istream& stream, BWT& bwt) {
    stream >> bwt._data;
    return stream;
}
