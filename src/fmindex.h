#ifndef fmindex_h_
#define fmindex_h_

#include "alphabet.h"
#include "bwt.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>

//
// FMMarkers - Marker classes used in the FM-index implementation
//

template< class Storage >
class Marker {
public:
    Marker() : unitIndex(0) {
    }

    size_t total() const {
        return std::accumulate(&counts[0], &counts[0] + counts.size(), 0);
    }

    bool operator==(const Marker& m) {
        return counts == m.counts && unitIndex == m.unitIndex;
    }

    DNAAlphabet::AlphaCount< Storage > counts;
    Storage unitIndex;
};

//
// LargeMarker - To allow random access to the 
// BWT symbols and implement the occurrence array
// we keep a vector of symbol counts every D1 symbols.
// These counts are the absolute number of times each
// symbol has been seen up to that point.
// 
// unit - The index in the RLString of the run that starts after
// this marker. That is, if C = getActualPosition(), then
// the run containing the B[C] is at unitIndex. This is not necessary
// a valid index if there is a marker after the last symbol in the BWT
//
typedef Marker< uint64_t > LargeMarker;
typedef Marker< uint16_t > SmallMarker;

typedef std::vector< LargeMarker > LargeMarkerList;
typedef std::vector< SmallMarker > SmallMarkerList;

const size_t DEFAULT_SAMPLE_RATE_SMALL = 128;
const size_t DEFAULT_SAMPLE_RATE_LARGE = 8192;

class FMIndex {
public:
    FMIndex(size_t sampleRate = DEFAULT_SAMPLE_RATE_SMALL) : _sampleRate(sampleRate) {
        initialize();
    }
    FMIndex(const BWT& bwt, size_t sampleRate = DEFAULT_SAMPLE_RATE_SMALL) : _bwt(bwt), _sampleRate(sampleRate) {
        initialize();
    }
    FMIndex(const SuffixArray& sa, const DNASeqList& sequences, size_t sampleRate = DEFAULT_SAMPLE_RATE_SMALL) : _bwt(sa, sequences), _sampleRate(sampleRate) {
        initialize();
    }

    size_t getPC(char c) const {
        return _pred[DNAAlphabet::torank(c)];
    }
    size_t getOcc(char c, size_t i) const;
    DNAAlphabet::AlphaCount64 getOcc(size_t i) const;

    void info() const;
private:
    void initialize();

    friend std::ostream& operator<<(std::ostream& stream, const FMIndex& index);
    friend std::istream& operator>>(std::istream& stream, FMIndex& index);

    BWT _bwt;
    DNAAlphabet::AlphaCount64 _pred;
    LargeMarkerList _lmarkers;
    SmallMarkerList _smarkers;
    size_t _sampleRate;
};

#endif // fmindex_h_
