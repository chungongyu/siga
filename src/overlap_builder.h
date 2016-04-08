#ifndef overlap_builder_h_
#define overlap_builder_h_

#include "fmindex.h"
#include "kseq.h"
#include "utils.h"

#include <bitset>
#include <iostream>
#include <vector>

//
// Flags indicating how a given read was aligned to the FM-index
//
struct AlignFlags {
public:
    AlignFlags() {
    }
    AlignFlags(bool qr, bool tr, bool qc) {
        _data.set(QUERYREV_BIT, qr);
        _data.set(TARGETREV_BIT, tr);
        _data.set(QUERYCOMP_BIT, qc);
    }
private:
    static const size_t QUERYREV_BIT  = 0;
    static const size_t TARGETREV_BIT = 1;
    static const size_t QUERYCOMP_BIT = 2;
    std::bitset< 3 > _data;
};

struct OverlapResult;
struct OverlapBlock;
typedef std::vector< OverlapBlock > OverlapBlockList;

//
// OverlapBuilder - Implements all the logic for finding
//    and outputting overlaps for sequence reads
//
class OverlapBuilder {
public:
    OverlapBuilder(const FMIndex* fmi, const FMIndex* rfmi, const std::string& prefix="default") : _fmi(fmi), _rfmi(rfmi), _prefix(prefix) {
    }

    bool build(DNASeqReader& reader, size_t minOverlap, std::ostream& output, size_t threads=1, size_t* processed=NULL) const;
    bool build(const std::string& input, size_t minOverlap, const std::string& output, size_t threads=1, size_t* processed=NULL) const;
    
    OverlapResult overlap(const DNASeq& read, size_t minOverlap, OverlapBlockList* blocks) const;
private:
    OverlapResult overlap(const std::string& seq, const FMIndex* fmi, const FMIndex* rfmi, size_t minOverlap, OverlapBlockList* overlaps, OverlapBlockList* contains) const;
    bool hits2asqg(std::istream& input, std::ostream& output) const;

    const FMIndex* _fmi;
    const FMIndex* _rfmi;
    std::string _prefix;
};

#endif // overlap_builder_h_
