#include "fmindex.h"
#include "bwt.h"

void FMIndex::initialize() {
    size_t counts[DNAAlphabet::ALL_SIZE];
    memset(counts, 0, SIZEOF_ARRAY(counts));
    size_t total = 0;

    const RLString& runs = _bwt.strings();
    for (size_t i = 0; i < runs.size(); ++i) {
        const RLUnit& run = runs[i];
        char c = (char)run;
        size_t len = run.count();

        // Update the count and advance the running total
        counts[DNAAlphabet::torank(c)] += len;
        total += len;

        // Check whether to place a new large marker

    }

    // Initialize C(a)
    memset(_pred, 0, SIZEOF_ARRAY(_pred));
    _pred[DNAAlphabet::torank('$')] = 0;
    for (size_t i = 1; i < DNAAlphabet::ALL_SIZE; ++i) {
        _pred[i] = _pred[i - 1] + counts[i - 1];
    }
}

std::ostream& operator<<(std::ostream& stream, const FMIndex& index) {
    stream << index._bwt;
    return stream;
}

std::istream& operator>>(std::istream& stream, FMIndex& index) {
    stream >> index._bwt;
    index.initialize();
    return stream;
}
