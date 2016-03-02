#include "suffix_array_builder.h"

#include <cstring>

//
// Implementation of induced copying algorithm by Nong, Zhang, Chan
// Follows implementation given as an appendix to their 2008 paper
// '\0' is the sentinenl in this algorithm
//
class SAISBuilder : public SuffixArrayBuilder {
public:
    bool build(const DNASeqList& reads, SuffixArray* sa) {
        size_t num_strings = reads.size();

        // In the multiple strings case, we need a 2D bit array
        // to hold the L/S types for the suffixes
        char** type_array = new char*[num_strings];
        for (size_t i = 0; i < num_strings; ++i) {
            const DNASeq& read = reads[i];
            size_t num_bytes = read.seq.length() / 8 + 1;
            type_array[i] = new char[num_bytes];
            memset(type_array[i], 0, num_bytes);
        }

        // Classify each suffix as being L or S type
        for (size_t i = 0; i < num_strings; ++i) {
            const DNASeq& read = reads[i];
            size_t len = read.seq.length() + 1;

            // The empty suffix ($) for each string is defined to be S type
            // and hence the next suffix must be L type
            setBit(type_array, i, len - 1, 1);
            setBit(type_array, i, len - 2, 0);
            for (size_t j = len - 2; j > 0; --j) {
                char curr = read.seq[j - 1], next = read.seq[j];
                bool type = (curr < next || (curr == next && getBit(type_array, i, j + 1) == 1));
                setBit(type_array, i, j, type);
            }
        }

        return true;
    }

private:
    // set the element to b
    void setBit(char** bit_array, size_t i, size_t j, bool b) {
        char* ba = bit_array[i];
        size_t block = j / 8, offset = j % 8;
        ba[block] = (b ? (_MASK[offset] | ba[block]) : (~_MASK[offset] & ba[block]));
    }
    bool getBit(char** bit_array, size_t i, size_t j) {
        return bit_array[i][j / 8] & _MASK[j % 8] ? true : false;
    }

    static unsigned char _MASK[8];
};

unsigned char SAISBuilder::_MASK[8] = {0x80,0x40,0x20,0x10,0x08,0x04,0x02,0x01};

SuffixArrayBuilder* SuffixArrayBuilder::create(const std::string& algorithm) {
    if (algorithm == "SAIS") {
        return new SAISBuilder();
    }
    return NULL;
}

