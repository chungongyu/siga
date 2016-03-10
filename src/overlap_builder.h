#ifndef overlap_builder_h_
#define overlap_builder_h_

#include "kseq.h"
#include "fmindex.h"

#include <vector>

struct OverlapBlock {
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

    void overlap(const DNASeq& read, size_t minOverlap, OverlapBlockList* blocks) const;
private:
    void overlap(const std::string& seq, const FMIndex* fmi, const FMIndex* rfmi, size_t minOverlap, OverlapBlockList* overlaps, OverlapBlockList* contains) const;

    const FMIndex* _fmi;
    const FMIndex* _rfmi;
};

#endif // overlap_builder_h_
