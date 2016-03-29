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
SeqCoord SeqCoord::complement() const {
    size_t s = 0, e = -1;
    if (isFull()) {
        e = seqlen - 1;
    } else if (isEmpty()) {
    } else if (isLeftExtreme()) {
        s = interval.end + 1;
        e = seqlen - 1;
    } else {
        assert(isRightExtreme());
        e = interval.start - 1;
    }
    return SeqCoord(s, e, seqlen);
}

std::ostream& operator<<(std::ostream& stream, const SeqCoord& c) {
    stream << c.interval << " " << c.seqlen;
    return stream;
}

std::istream& operator>>(std::istream& stream, SeqCoord& c) {
    stream >> c.interval >> c.seqlen;
    return stream;
}

//
// Match
//

// Calculation the translation offset to shift
// a coord[1] position to a coord[0]. This must be calculated
// using canonical coordinates
size_t Match::translate10() const {
    assert(coords[0].length() == coords[1].length());
    if (isRC) {
        SeqCoord c = coords[0].flip();
        return c.interval.start - coords[1].interval.start;
    }
    return coords[0].interval.start - coords[1].interval.start;
}

SeqCoord Match::translate10(const SeqCoord& c) const {
    assert(coords[0].length() == coords[1].length());
    size_t t = translate10();
    SeqCoord r(c.interval.start + t, c.interval.end + t, coords[0].seqlen);
    if (isRC) {
        r.flip();
    }
    return r;
}

std::ostream& operator<<(std::ostream& stream, const Match& m) {
    stream << m.coords[0] << " " << m.coords[1] << " " << m.isRC << " " << m.numDiff;
    return stream;
}

std::istream& operator>>(std::istream& stream, Match& m) {
    stream >> m.coords[0] >> m.coords[1] >> m.isRC >> m.numDiff;
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
