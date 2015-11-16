#include "bwt.h"
#include "suffix_array.h"

BWT::BWT(const SuffixArray& as, const DNASeqList& sequences) {
}

class BWTWriter {
};

std::ostream& operator<<(std::ostream& stream, const BWT& bwt) {
    return stream;
}

std::istream& operator>>(std::istream& stream, const BWT& bwt) {
    return stream;
}
