#ifndef suffix_array_builder_h_
#define suffix_array_builder_h_

#include "kseq.h"

#include <string>

class SuffixArray;

class SuffixArrayBuilder {
public:
    static SuffixArrayBuilder* create(const std::string& algorithm);

    virtual bool build(const DNASeqList& sequences, SuffixArray* sa) = 0;
};

#endif // suffix_array_builder_h_
