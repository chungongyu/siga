#ifndef fmindex_h_
#define fmindex_h_

#include "alphabet.h"
#include "bwt.h"
#include "utils.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>

//
// FMMarkers - Marker classes used in the FM-index implementation
//

//
// LargeMarker - To allow random access to the 
// BWT symbols and implement the occurrence array
// we keep a vector of symbol counts every D1 symbols.
// These counts are the absolute number of times each
// symbol has been seen up to that point.
// 
class LargeMarker {
public:
    LargeMarker() : unitIndex(0) {
        memset(counts, 0, SIZEOF_ARRAY(counts));
    }

    size_t position() const {
        return std::accumulate(counts, counts + DNAAlphabet::ALL_SIZE, 0);
    }

    bool operator==(const LargeMarker& m) {
        for (size_t i = 0; i < DNAAlphabet::ALL_SIZE; ++i) {
            if (counts[i] != m.counts[i]) {
                return false;
            }
        }
        return true;
    }

    size_t counts[DNAAlphabet::ALL_SIZE];
    // The index in the RLString of the run that starts after
    // this marker. That is, if C = getActualPosition(), then
    // the run containing the B[C] is at unitIndex. This is not necessary
    // a valid index if there is a marker after the last symbol in the BWT
    size_t unitIndex;
};

typedef std::vector< LargeMarker > LargeMarkerList;

class FMIndex {
public:
    FMIndex() {
        initialize();
    }
    FMIndex(const BWT& bwt) : _bwt(bwt) {
        initialize();
    }
    FMIndex(const SuffixArray& sa, const DNASeqList& sequences) : _bwt(sa, sequences) {
        initialize();
    }
private:
    void initialize();

    friend std::ostream& operator<<(std::ostream& stream, const FMIndex& index);
    friend std::istream& operator>>(std::istream& stream, FMIndex& index);

    BWT _bwt;
    size_t _pred[DNAAlphabet::ALL_SIZE];
};

#endif // fmindex_h_
