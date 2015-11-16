#ifndef bwt_h_
#define bwt_h_

#include "kseq.h"

#include <iostream>

class SuffixArray;

//
// Run-length encoded Burrows Wheeler transform
//
class BWT {
public:
    BWT(const SuffixArray& sa, const DNASeqList& sequences);

private:
    friend std::ostream& operator<<(std::ostream& stream, const BWT& bwt);
    friend std::istream& operator>>(std::istream& stream, BWT& bwt);
};

#endif // bwt_h_
