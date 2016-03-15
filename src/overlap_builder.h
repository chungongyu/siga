#ifndef overlap_builder_h_
#define overlap_builder_h_

#include "fmindex.h"
#include "kseq.h"

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


struct OverlapBlock;

typedef std::vector< OverlapBlock > OverlapBlockList;

//
// OverlapBuilder - Implements all the logic for finding
//    and outputting overlaps for sequence reads
//
class OverlapBuilder {
public:
    OverlapBuilder(const FMIndex* fmi, const FMIndex* rfmi) : _fmi(fmi), _rfmi(rfmi) {
    }

    void overlap(const DNASeq& read, size_t minOverlap, OverlapBlockList* blocks) const;
private:
    void overlap(const std::string& seq, const FMIndex* fmi, const FMIndex* rfmi, size_t minOverlap, OverlapBlockList* overlaps, OverlapBlockList* contains) const;

    const FMIndex* _fmi;
    const FMIndex* _rfmi;
};

#endif // overlap_builder_h_
