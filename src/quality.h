#ifndef quality_h_
#define quality_h_

#include <cassert>
#include <cmath>

//
// Quality - functions for manipulating quality values
//
namespace Quality {
    namespace Phred {
        static const int DEFAULT_SCORE = 15;
        static const int PHRED64_DIFF = 31;

        // Convert the quality character from phred64 to phred33 encoding
        inline char _64to33(char c) {
            return (int)c - PHRED64_DIFF;
        }

        // Returns true if the character c is a valid Phred33 encoding of
        // a quality value in the range [0, 60]
        inline bool isValid(char c) {
            int p = (int)c - 33;
            return 0 <= p and p <= 60;
        }

        // Phred score transformations
        inline int fromchar(char b) {
            uint8_t v = b;
            assert(v >= 33);
            return v - 33;
        }
        inline char tochar(int p) {
            uint8_t v = (p <= 93) ? p : 93;
            char c = v + 33;
            return c;
        }
        inline int fromprob(double p) {
            return (int)(std::round(-10.0f * std::log10(p)));
        }
    };
};

#endif // quality_h_
