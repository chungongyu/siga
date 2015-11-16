#ifndef suffix_array_builder_h_
#define suffix_array_builder_h_

#include "kseq.h"

class SuffixArray;

class SuffixArrayBuilder {
public:
    virtual bool build(const DNASeqList& sequences, SuffixArray& sa) = 0;
};

#endif // suffix_array_builder_h_
