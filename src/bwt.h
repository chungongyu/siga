#ifndef bwt_h_
#define bwt_h_

#include <iostream>

#include "kseq.h"
#include "rlstring.h"

class SuffixArray;

//
// Run-length encoded Burrows Wheeler transform
//
class BWT {
 public:
  BWT() : _strings(0), _suffixes(0) {
  }
  BWT(const SuffixArray& sa, const DNASeqList& sequences);

  const RLString& str() const {
    return _runs;
  }
  size_t length() const {
    return (size_t)_suffixes;
  }

 private:
  friend std::ostream& operator<<(std::ostream& stream, const BWT& bwt);
  friend std::istream& operator>>(std::istream& stream, BWT& bwt);
  friend class BWTReader;
  friend class BWTWriter;

  RLString _runs;   // The run-length encoded string
  uint64_t _strings;  // The number of strings in the collection
  uint64_t _suffixes;   // The total length of the bw string
};

#endif  // bwt_h_
