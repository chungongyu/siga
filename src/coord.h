#ifndef coord_h_
#define coord_h_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

//
// Interval - A pair of integers denoting a closed interval
//
class Interval {
 public:
  Interval() : start(0), end(-1) {
  }

  Interval(size_t s, size_t e) : start(s), end(e) {
  }

  size_t length() const {
    return end + 1 - start;
  }

  void offset(size_t delta) {
    start += delta;
    end += delta;
  }

  // Flip a single position p to the reverse strand for a sequence of length l
  static size_t flip(size_t p, size_t l) {
    assert(l > p);
    return l - p - 1;
  }
  void flip(size_t l) {
    size_t t = start;
    start = flip(end, l);
    end = flip(t, l);
  }

  // Precondition: s1 >= e1 and s2 >= e2
  // Return true if the coordinates intersect
  static bool isIntersecting(size_t s1, size_t e1, size_t s2, size_t e2) {
    assert(s1 <= e1 && s2 <= e2);
    return !(s1 > e2 || s2 > e1);
  }

  friend std::ostream& operator<<(std::ostream& stream, const Interval& i);
  friend std::istream& operator>>(std::istream& stream, Interval& i);

  // data members
  size_t start;
  size_t end;
};

//
// SeqCoord - A data structure holding the coordinate of a substring of a sequence
//      which consists of an interval and the length of the string. Used to
//      build matches and overlaps
//
class SeqCoord {
 public:
  SeqCoord() : seqlen(0) {
  }
  SeqCoord(const Interval& i, size_t l) : interval(i), seqlen(l) {
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
    return interval.length() == seqlen;
  }
  bool isEmpty() const {
    return interval.start == 0 && interval.end == -1;
  }
  size_t length() const {
    return interval.length();
  }
  void extend(size_t len) {
    if (isLeftExtreme()) {
      interval.end += len;
    } else {
      assert(isRightExtreme() && interval.start >= len);
      interval.start -= len;
    }
  }
  void stretch(size_t len) {
    seqlen += len;
    interval.end += len;
  }
  // Flip mirrors the coordinates so they are on the other strand
  void flip() {
    interval.flip(seqlen);
  }
  SeqCoord flip() const {
    SeqCoord c(interval, seqlen);
    c.flip();
    return c;
  }

  SeqCoord complement() const;

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
  Match() : isRC(false), numDiff(0) {
  }
  Match(const SeqCoord& c1, const SeqCoord& c2, bool isRC, size_t nd) : isRC(isRC), numDiff(nd) {
    coords[0] = c1;
    coords[1] = c2;
  }
  Match(size_t s1, size_t e1, size_t l1, size_t s2, size_t e2, size_t l2, bool isRC, size_t nd)
      : isRC(isRC), numDiff(nd) {
    coords[0] = SeqCoord(s1, e1, l1);
    coords[1] = SeqCoord(s2, e2, l2);
  }

  size_t length() const {
    assert(coords[0].length() == coords[1].length());
    return coords[0].length();
  }

  bool isContainment() const {
    return coords[0].isContained() || coords[1].isContained();
  }

  size_t translate10() const;
  SeqCoord translate10(const SeqCoord& coord) const;

  friend std::ostream& operator<<(std::ostream& stream, const Match& m);
  friend std::istream& operator>>(std::istream& stream, Match& m);

  SeqCoord coords[2];
  bool isRC;
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
  Overlap(const std::string& i1, const SeqCoord& c1, const std::string& i2,
          const SeqCoord& c2, bool isRC, size_t nd)
      : match(c1, c2, isRC, nd) {
    id[0] = i1;
    id[1] = i2;
  }
  Overlap(const std::string& i1, size_t s1, size_t e1, size_t l1,
         const std::string& i2, size_t s2, size_t e2, size_t l2, bool isRC, size_t nd)
      : match(s1, e1, l1, s2, e2, l2, isRC, nd) {
    id[0] = i1;
    id[1] = i2;
  }

  bool isContainment() const {
    return match.isContainment();
  }
  size_t containedIdx() const {
    // The verts are mutually contained, return the lexographically lower id
    if (match.coords[0].isContained() && match.coords[1].isContained()) {
      return id[0] < id[1] ? 1 : 0;
    } else if (match.coords[0].isContained()) {
      return 0;
    }
    assert(match.coords[1].isContained());
    return 1;
  }

  friend std::ostream& operator<<(std::ostream& stream, const Overlap& o);
  friend std::istream& operator>>(std::istream& stream, Overlap& o);

  std::string id[2];
  Match match;
};

typedef std::vector<Overlap> OverlapList;

#endif  // coord_h_
