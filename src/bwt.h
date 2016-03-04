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
    BWT() : _strings(0), _suffixes(0) {
    }
    BWT(const SuffixArray& sa, const DNASeqList& sequences);

private:
    friend std::ostream& operator<<(std::ostream& stream, const BWT& bwt);
    friend std::istream& operator>>(std::istream& stream, BWT& bwt);
    friend class BWTReader;
    friend class BWTWriter;

    RLList _runs;       // The run-length encoded string
    size_t _strings;    // The number of strings in the collection
    size_t _suffixes;   // The total length of the bw string
};

enum BWFlag {       
    BWF_NOFMI = 0,
    BWF_HASFMI
};

//
// Read a run length encoded binary BWT file from disk
//
class BWTReader {
public:
    BWTReader(std::istream& stream) : _stream(stream) {
    }

    bool read(BWT& bwt);
private:
    bool readHeader(size_t& num_strings, size_t& num_suffixes, BWFlag& flag);
    bool readRuns(RLList& runs, size_t numRuns);

    std::istream& _stream;
    size_t _numRuns;
};

//
// Write a run-length encoded BWT to a binary file
//
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
