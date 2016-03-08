#include "coord.h"

// Interval
std::ostream& operator<<(std::ostream& stream, const Interval& r) {
    stream << r.start << " " << r.end;
    return stream;
}

std::istream& operator>>(std::istream& stream, Interval& r) {
    stream >> r.start >> r.end;
    return stream;
}

// SeqCoord
std::ostream& operator<<(std::ostream& stream, const SeqCoord& c) {
    stream << c.interval << " " << c.seqlen;
    return stream;
}

std::istream& operator>>(std::istream& stream, SeqCoord& c) {
    stream >> c.interval >> c.seqlen;
    return stream;
}

// Match
std::ostream& operator<<(std::ostream& stream, const Match& m) {
    stream << m.coords[0] << " " << m.coords[1] << " " << m.isReverse << " " << m.numDiff;
    return stream;
}

std::istream& operator>>(std::istream& stream, Match& m) {
    stream >> m.coords[0] >> m.coords[1] >> m.isReverse >> m.numDiff;
    return stream;
}

// Overlap
std::ostream& operator<<(std::ostream& stream, const Overlap& o) {
    stream << o.id[0] << " " << o.id[1] << " " << o.match;
    return stream;
}

std::istream& operator>>(std::istream& stream, Overlap& o) {
    stream >> o.id[0] >> o.id[1] >> o.match;
    return stream;
}
