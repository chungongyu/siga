#ifndef bwt_h_
#define bwt_h_

#include "kseq.h"
#include "rlstring.h"

#include <iostream>

class SuffixArray;

//
// Run-length encoded Burrows Wheeler transform
//
class BWT {
public:
    BWT() {
    }
    BWT(const SuffixArray& sa, const DNASeqList& sequences);

private:
    friend std::ostream& operator<<(std::ostream& stream, const BWT& bwt);
    friend std::istream& operator>>(std::istream& stream, BWT& bwt);

    std::string _data;
};

enum BWFlag {       
    BWF_NOFMI = 0,
    BWF_HASFMI
};

class BWTWriter {
public:
    BWTWriter(std::ostream& stream) : _stream(stream), _numRuns(0), _posRun(0) {
    }
    
    bool write(const SuffixArray& sa, const DNASeqList& sequences);

private:
    bool writeHeader(size_t num_strings, size_t num_suffixes, BWFlag flag);
    bool writeChar(char c);
    bool writeRun(const RLUnit& run);
    bool finalize();

    RLUnit _currRun;
    size_t _numRuns;
    std::streampos _posRun;

    std::ostream& _stream;
};

#endif // bwt_h_
