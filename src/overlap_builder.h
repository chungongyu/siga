#ifndef overlap_builder_h_
#define overlap_builder_h_

#include "fmindex.h"
#include "kseq.h"
#include "utils.h"

#include <vector>

//
// Flags indicating how a given read was aligned to the FM-index
//
struct AlignFlags {
public:
    AlignFlags() {
    }
    AlignFlags(bool qr, bool tr, bool qc) {
    }
private:
    static const size_t QUERYREV_BIT  = 0;
    static const size_t TARGETREV_BIT = 1;
    static const size_t QUERYCOMP_BIT = 2;
    uint8_t _data;
};

struct OverlapResult {
    OverlapResult() : substring(false), aborted(false) {
    }
    bool substring;
    bool aborted;
};

//
// A pair of intervals used for bidirectional searching a FM-index/reverse FM-index
//
class IntervalPair {
public:
    IntervalPair() {
    }
    bool valid() const {
        for (size_t i = 0; i < SIZEOF_ARRAY(_intervals); ++i) {
            if (!_intervals[i].valid()) {
                return false;
            }
        }
        return true;
    }
    FMIndex::Interval& operator[](size_t i) {
        return _intervals[i];
    }
    const FMIndex::Interval& operator[](size_t i) const {
        return _intervals[i];
    }

    void init(char c, const FMIndex* index, const FMIndex* rindex) {
        _intervals[0].init(c,  index);
        _intervals[1].init(c, rindex);
    }
    void updateL(char c, const FMIndex* index) {
        // Update the left index using the difference between the AlphaCounts in the reverse table
        DNAAlphabet::AlphaCount64 l = index->getOcc(_intervals[1].lower - 1);
        DNAAlphabet::AlphaCount64 u = index->getOcc(_intervals[1].upper);
        updateL(c, index, l, u);
    } 
    void updateR(char c, const FMIndex* index) {
    }
private:
    void updateL(char c, const FMIndex* index, const DNAAlphabet::AlphaCount64& l, const DNAAlphabet::AlphaCount64& u) {
        DNAAlphabet::AlphaCount64 diff = u - l;
        // Update the left index using the difference between the AlphaCounts in the reverse table
        _intervals[1].lower = _intervals[1].lower + std::accumulate(&diff[0], &diff[0] + DNAAlphabet::torank(c), 0);
        _intervals[1].upper = _intervals[1].lower + diff[DNAAlphabet::torank(c)] - 1;

        // Update the left index directly
        size_t pb = index->getPC(c);
        _intervals[0].lower = pb + l[DNAAlphabet::torank(c)];
        _intervals[0].upper = pb + u[DNAAlphabet::torank(c)] - 1;
    }

    FMIndex::Interval _intervals[2];
};

struct OverlapBlock {
    OverlapBlock(const IntervalPair& probe, const IntervalPair& ranges, size_t length) : probe(probe), ranges(ranges), length(length) {
    }

    IntervalPair ranges;
    IntervalPair probe;
    size_t length;
};


typedef std::vector< OverlapBlock > OverlapBlockList;

//
// OverlapBuilder - Implements all the logic for finding
//    and outputting overlaps for sequence reads
//
class OverlapBuilder {
public:
    OverlapBuilder(const FMIndex* fmi, const FMIndex* rfmi) : _fmi(fmi), _rfmi(rfmi) {
    }

    OverlapResult overlap(const DNASeq& read, size_t minOverlap, OverlapBlockList* blocks) const;
private:
    OverlapResult overlap(const std::string& seq, const FMIndex* fmi, const FMIndex* rfmi, size_t minOverlap, OverlapBlockList* overlaps, OverlapBlockList* contains) const;

    const FMIndex* _fmi;
    const FMIndex* _rfmi;
};

#endif // overlap_builder_h_
