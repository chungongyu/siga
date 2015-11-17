#include "fmindex.h"
#include "bwt.h"

std::ostream& operator<<(std::ostream& stream, const FMIndex& index) {
    stream << index._bwt;
    return stream;
}

std::istream& operator>>(std::istream& stream, FMIndex& index) {
    stream >> index._bwt;
    return stream;
}
