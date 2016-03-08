#ifndef coord_h_
#define coord_h_

#include <iostream>

//
// Interval - A pair of integers denoting a closed interval
//
class Interval {
public:
    Interval() : start(0), end(0) {
    }

    Interval(size_t s, size_t e) : start(s), end(e) {
    }

    size_t length() const {
        return end + 1 - start;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Interval& i);
    friend std::istream& operator>>(std::istream& stream, Interval& i);    

    // data members
    size_t start;
    size_t end;
};

//
// SeqCoord - A data structure holding the coordinate of a substring of a sequence
//            which consists of an interval and the length of the string. Used to 
//            build matches and overlaps
//
class SeqCoord {
public:
    SeqCoord() : seqlen(0) {
    }
    SeqCoord(size_t s, size_t e, size_t l) : interval(s, e), seqlen(l) {
    }

    bool isLeftExtreme() const {
        return interval.start == 0;
    }
    bool isRightExtreme() const {
        return interval.end + 1 == seqlen;
    }
    bool isExtreme() const {
        return isLeftExtreme() || isRightExtreme();
    }
    bool isContained() const {
        return isLeftExtreme() && isRightExtreme();
    }
    bool isFull() const {
        interval.length() == seqlen;
    }

    friend std::ostream& operator<<(std::ostream& stream, const SeqCoord& c);
    friend std::istream& operator>>(std::istream& stream, SeqCoord& c);

    Interval interval;
    size_t seqlen;
};

//
// Match - A pair of coordinates representing the overlapping regions of two sequences
//
class Match {
public:
    Match() {
    }
    Match(const SeqCoord& c1, const SeqCoord& c2, bool isRC, size_t nd) : isReverse(isRC), numDiff(nd) {
        coords[0] = c1;
        coords[1] = c2;
    }
    Match(size_t s1, size_t e1, size_t l1, size_t s2, size_t e2, size_t l2, bool isRC, size_t nd) : isReverse(isRC), numDiff(nd) {
        coords[0] = SeqCoord(s1, e1, l1);
        coords[1] = SeqCoord(s2, e2, l2);
    }

    friend std::ostream& operator<<(std::ostream& stream, const Match& m);
    friend std::istream& operator>>(std::istream& stream, Match& m);

    SeqCoord coords[2];
    bool isReverse;
    size_t numDiff;
};

// 
// Overlap
// 
class Overlap {
public:
    Overlap() {
    }
    Overlap(const std::string& i1, const std::string& i2, const Match& m) : match(m) {
        id[0] = i1;
        id[1] = i2;
    }
    Overlap(const std::string& i1, const SeqCoord& c1, const std::string& i2, const SeqCoord& c2, bool isRC, size_t nd) : match(c1, c2, isRC, nd) {
        id[0] = i1;
        id[1] = i2;
    }
    Overlap(const std::string& i1, size_t s1, size_t e1, size_t l1, const std::string& i2, size_t s2, size_t e2, size_t l2, bool isRC, size_t nd) : match(s1, e1, l1, s2, e2, l2, isRC, nd) {
        id[0] = i1;
        id[1] = i2;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Overlap& o);
    friend std::istream& operator>>(std::istream& stream, Overlap& o);

    std::string id[2];
    Match match;
};

#endif // coord_h_
