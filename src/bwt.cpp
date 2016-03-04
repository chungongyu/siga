#include "bwt.h"
#include "rlstring.h"
#include "suffix_array.h"

enum BWFlag {       
    BWF_NOFMI = 0,
    BWF_HASFMI
};

const uint16_t BWT_FILE_MAGIC = 0xCACA;

BWT::BWT(const SuffixArray& sa, const DNASeqList& sequences) : _data(sa.size(), ' ') {
    for (size_t i = 0; i < sa.size(); ++i) {
        const SuffixArray::Elem& elem = sa[i];
        const DNASeq& read = sequences[elem.i];
        char c = (elem.j == 0) ? '$' : read.seq[elem.j - 1];
        _data[i] = c;
    }
}

class BWTWriter {
public:
    BWTWriter(std::ostream& stream) : _stream(stream), _numRuns(0), _posRun(0) {
    }
    
    bool write(const SuffixArray& sa, const DNASeqList& sequences) {
        size_t num_strings = sa.strings(), num_suffixes = sa.size();

        if (!writeHeader(num_strings, num_suffixes, BWF_NOFMI)) {
            return false;
        }
        for (size_t i = 0; i < num_suffixes; ++i) {
            const SuffixArray::Elem& elem = sa[i];
            const DNASeq& read = sequences[elem.i];
            char c = (elem.j == 0 ? '$' : read.seq[elem.j - 1]);
            if (!writeChar(c)) {
                return false;
            }
        }
        if (!finalize()) {
            return false;
        }
        return true;
    }

private:
    bool writeHeader(size_t num_strings, size_t num_suffixes, BWFlag flag) {
        if (!_stream.write((const char *)&BWT_FILE_MAGIC, sizeof(BWT_FILE_MAGIC))) {
            return false;
        }
        if (!_stream.write((const char *)&num_strings, sizeof(num_strings))) {
            return false;
        }
        if (!_stream.write((const char *)&num_suffixes, sizeof(num_suffixes))) {
            return false;
        }

        // Here we do not know the number of runs that are going to be written to the file
        // so we save the offset in the file and write a dummy value. After the bwt string
        // has been written, we return here and fill in the correct value
        _posRun = _stream.tellp();
        _numRuns = 0;
        _stream.write((const char *)&_numRuns, sizeof(_numRuns));

        return true;
    }
    bool writeChar(char c) {
        if (_currRun.initialized()) {
            if (_currRun == c && !_currRun.full()) {
                ++_currRun;
            } else {
                // Write out the old run and start a new one
                if (!writeRun(_currRun)) {
                    return false;
                }
                _currRun = RLUnit(c);
            }
        } else {
            // Start a new run 
            _currRun = RLUnit(c);
        }
        return true;
    }
    bool finalize() {
        if (_currRun.initialized()) {
            if (!writeRun(_currRun)) {
                return false;
            }
        }
        _stream.seekp(_posRun);
        _stream.write((const char *)&_numRuns, sizeof(_numRuns));
        _stream.seekp(std::ios_base::end);
        return _stream;
    }

    bool writeRun(const RLUnit& run) {
        if (!_stream.write((const char *)&run.data, sizeof(run.data))) {
            return false;
        }
        ++_numRuns;
        return true;
    }

    RLUnit _currRun;
    size_t _numRuns;
    size_t _posRun;

    std::ostream& _stream;
};

std::ostream& operator<<(std::ostream& stream, const BWT& bwt) {
    stream << bwt._data;
    return stream;
}

std::istream& operator>>(std::istream& stream, BWT& bwt) {
    stream >> bwt._data;
    return stream;
}
