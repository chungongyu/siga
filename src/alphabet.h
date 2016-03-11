#ifndef alphabet_h_
#define alphabet_h_

#include <cstdint>
#include <cstring>

namespace DNAAlphabet {
    const size_t ALL_SIZE = 5;
    const char DNA_ALL[ALL_SIZE] = {'$', 'A', 'C', 'G', 'T'};

    inline int torank(char c) {
        static int RANK_ALL[256] = {
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,
            0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
        };
        return RANK_ALL[(size_t)c];
    }

    inline char tochar(int rank) {
        return DNA_ALL[rank];
    }

    template< class Storage >
    class AlphaCount {
    public:
        AlphaCount() {
            memset(_data, 0, sizeof(_data));
        }

        Storage& operator[](char c) {
            return _data[torank(c)];
        }
        const Storage& operator[](char c) const {
            return _data[torank(c)];
        }
        size_t size() const {
            return ALL_SIZE;
        }
        const Storage* array() const {
            return _data;
        }
    private:
        Storage _data[ALL_SIZE];
    };

    typedef AlphaCount< uint64_t > AlphaCount64;
    typedef AlphaCount< uint32_t > AlphaCount32;
    typedef AlphaCount< uint16_t > AlphaCount16;
};

#endif // alphabet_h_

