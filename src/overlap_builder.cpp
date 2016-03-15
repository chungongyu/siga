#include "overlap_builder.h"
#include "utils.h"

//
// Interval holds a pair of integers which delineate an alignment of some string
//       to a BWT/FM-index
//
class Interval {
public:
    Interval() : lower(0), upper(0) {
    }
    bool valid() const {
        return upper > lower;
    }
    void init(char c, const FMIndex* index) {
        lower = index->getPC(c);
        upper = lower + index->getOcc(c, index->length() - 1) - 1;
    }

    size_t lower;
    size_t upper;
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
    Interval& operator[](size_t i) {
        return _intervals[i];
    }
    const Interval& operator[](size_t i) const {
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

    Interval _intervals[2];
};

struct OverlapBlock {
    OverlapBlock(const IntervalPair& probe, const IntervalPair& ranges, size_t length) : probe(probe), ranges(ranges), length(length) {
    }

    IntervalPair ranges;
    IntervalPair probe;
    size_t length;
};

void OverlapBuilder::overlap(const DNASeq& read, size_t minOverlap, OverlapBlockList* blocks) const {
    // The complete set of overlap blocks are collected in obWorkingList
    // The filtered set (containing only irreducible overlaps) are placed into pOBOut
    // by calculateIrreducibleHits
    const std::string& seq = read.seq;

    overlap(seq, _fmi, _rfmi, minOverlap, NULL, NULL);
    //overlap(make_complement_dna(seq), _fmi, _rfmi, minOverlap, NULL, NULL);
}

// Calculate the ranges in pBWT that contain a prefix of at least minOverlap basepairs that
// overlaps with a suffix of w. The ranges are added to the pOBList
void OverlapBuilder::overlap(const std::string& seq, const FMIndex* fmi, const FMIndex* rfmi, size_t minOverlap, OverlapBlockList* overlaps, OverlapBlockList* contains) const {
    // The algorithm is as follows:
    // We perform a backwards search using the FM-index for the string w.
    // As we perform the search we collect the intervals 
    // of the significant prefixes (len >= minOverlap) that overlap w.
    IntervalPair ranges;
    size_t l = seq.length();
    ranges.init(seq[l - 1], fmi, rfmi);

    // Collect the OverlapBlocks
    for (size_t i = l - 1; i > 1; --i) {
        // Compare the range of the suffix seq[i, l]
        ranges.updateL(seq[i - 1], fmi);

        if (l - i >= minOverlap) {
            // Calculate which of the prefixes that match w[i, l] are terminal
            // These are the proper prefixes (they are the start of a read)
            IntervalPair probe = ranges;
            probe.updateL('$', fmi);

            // The probe interval contains the range of proper prefixes
            if (probe[1].valid()) {
                overlaps->push_back(OverlapBlock(probe, ranges, l - i));
            }
        }
    }

    // Determine if this sequence is contained and should not be processed further
    ranges.updateL(seq[0], fmi);

    // Ranges now holds the interval for the full-length read
    // To handle containments, we output the overlapBlock to the final overlap block list
    // and it will be processed later
    // Two possible containment cases:
    // 1) This read is a substring of some other read
    // 2) This read is identical to some other read

    // Case 1 is indicated by the existance of a non-$ left or right hand extension
    // In this case we return no alignments for the string
    DNAAlphabet::AlphaCount64 lext = fmi->getOcc(ranges[0].upper) - fmi->getOcc(ranges[0].lower - 1);
    DNAAlphabet::AlphaCount64 rext = fmi->getOcc(ranges[1].upper) - fmi->getOcc(ranges[1].lower - 1);
    if (lext.hasDNA() || rext.hasDNA()) {
    } else {
        IntervalPair probe = ranges;
        probe.updateL('$', fmi);
        if (probe.valid()) {
            // terminate the contained block and add it to the contained list
            probe.updateR('$', rfmi);
            assert(probe.valid());
            contains->push_back(OverlapBlock(probe, ranges, l));
        }
    }
}
