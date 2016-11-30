#ifndef overlap_builder_h_
#define overlap_builder_h_

#include "fmindex.h"
#include "kseq.h"

#include <iostream>
#include <list>

struct OverlapResult;
struct OverlapBlock;
typedef std::list< OverlapBlock > OverlapBlockList;

//
// OverlapBuilder - Implements all the logic for finding
//    and outputting overlaps for sequence reads
//
class OverlapBuilder {
public:
    OverlapBuilder(const FMIndex* fmi, const FMIndex* rfmi, const std::string& prefix="default", bool irreducible=true, bool rc=true) : _fmi(fmi), _rfmi(rfmi), _prefix(prefix), _irreducible(irreducible), _rc(rc) {
    }

    bool build(DNASeqReader& reader, size_t minOverlap, std::ostream& output, size_t threads=1, size_t* processed=NULL) const;
    bool build(const std::string& input, size_t minOverlap, const std::string& output, size_t threads=1, size_t* processed=NULL) const;

    bool rmdup(DNASeqReader& reader, std::ostream& output, std::ostream& duplicates, size_t threads=1, size_t* processed=NULL) const;
    bool rmdup(const std::string& input, const std::string& output, const std::string& duplicates, size_t threads=1, size_t* processed=NULL) const;
    
    OverlapResult overlap(const DNASeq& read, size_t minOverlap, OverlapBlockList* blocks) const;
    OverlapResult duplicate(const DNASeq& read, OverlapBlockList* blocks) const;

private:
    const FMIndex* _fmi;
    const FMIndex* _rfmi;
    std::string _prefix;
    bool _irreducible;
    bool _rc; // reverse complement
};

#endif // overlap_builder_h_
