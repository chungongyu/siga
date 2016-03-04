#ifndef rlstring_h_
#define rlstring_h_

#include "alphabet.h"

#include <cassert>
#include <string>
#include <vector>

#define RL_COUNT_MASK   0x1F //00011111
#define RL_SYMBOL_MASK  0xE0 //11100000
#define RL_FULL_COUNT   31
#define RL_SYMBOL_SHIFT 5

//
// RLUnit - A run-length encoded unit of the FM-index
//
class RLUnit {
public:
    RLUnit() : data(0) {
    }
    RLUnit(char c) : data(1) {
        // Clear the current symbol
         data &= RL_COUNT_MASK;

         size_t code = DNAAlphabet::torank(c);
         code <<= RL_SYMBOL_SHIFT;
         data |= code;
    }

    // Returns true if the count cannot be incremented
    bool full() const {
         return count() == RL_FULL_COUNT;
    }
    bool empty() const {
        return count() == 0;
    }
    size_t count() const {
        return data & RL_COUNT_MASK;
    }
    bool initialized() const {
        return data > 0;
    }

    operator char() const {
        size_t code = data & RL_SYMBOL_MASK;
        code >>= RL_SYMBOL_SHIFT;
        return DNAAlphabet::tochar(code);
    }
    
    RLUnit& operator++() {
        assert(!full());
        ++data;
        return *this;
    }
    RLUnit& operator--() {
        assert(!empty());
        --data;
        return *this;
    }

    uint8_t data;
};

typedef std::vector< RLUnit > RLList;

class RLString {
public:
};

#endif // rlstring_h_
