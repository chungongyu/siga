#define BOOST_TEST_MODULE "siga.index"
#include <boost/test/included/unit_test.hpp>

#include "alphabet.h"
#include "rlstring.h"
#include "suffix_array.h"

BOOST_AUTO_TEST_SUITE(Indexer);

BOOST_AUTO_TEST_CASE(Alphabet_torank) {
    BOOST_CHECK_EQUAL(DNAAlphabet::torank('$'), 0);
    BOOST_CHECK_EQUAL(DNAAlphabet::torank('A'), 1);
    BOOST_CHECK_EQUAL(DNAAlphabet::torank('C'), 2);
    BOOST_CHECK_EQUAL(DNAAlphabet::torank('G'), 3);
    BOOST_CHECK_EQUAL(DNAAlphabet::torank('T'), 4);
}
BOOST_AUTO_TEST_CASE(Alphabet_tochar) {
    BOOST_CHECK_EQUAL(DNAAlphabet::tochar(0), '$');
    BOOST_CHECK_EQUAL(DNAAlphabet::tochar(1), 'A');
    BOOST_CHECK_EQUAL(DNAAlphabet::tochar(2), 'C');
    BOOST_CHECK_EQUAL(DNAAlphabet::tochar(3), 'G');
    BOOST_CHECK_EQUAL(DNAAlphabet::tochar(4), 'T');
}
BOOST_AUTO_TEST_CASE(Alphabet_func) {
    for (size_t i = 0; i < DNAAlphabet::ALL_SIZE; ++i) {
        BOOST_CHECK_EQUAL(i, DNAAlphabet::torank(DNAAlphabet::tochar(i)));
    }
}

BOOST_AUTO_TEST_CASE(RLUnit_test) {
    {
        RLUnit unit;
        BOOST_CHECK(!unit.initialized());
        BOOST_CHECK(unit.empty());
    }
    {
        RLUnit unit('A');
        BOOST_CHECK(!unit.empty());
        BOOST_CHECK_EQUAL(unit.count(), 1);
        BOOST_CHECK_EQUAL((char)unit, 'A');
        ++unit;
        BOOST_CHECK_EQUAL(unit.count(), 2);
        --unit;
        BOOST_CHECK_EQUAL(unit.count(), 1);
    }
}

BOOST_AUTO_TEST_CASE(RLString_test) {
    DNASeq read("test", "AAACGGGTA");

    RLString runs;

    RLUnit run;
    for (const auto& c : read.seq) {
        if (run.initialized()) {
            if (run == c && !run.full()) {
                ++run;
            } else {
                runs.push_back(run);
                run = RLUnit(c);
            }
        } else {
            run = RLUnit(c);
        }
    }
    if (run.initialized()) {
        runs.push_back(run);
    }

    BOOST_CHECK_EQUAL(runs.size(), 5);
    BOOST_CHECK_EQUAL(runs[0], 'A');
    BOOST_CHECK_EQUAL(runs[0].count(), 3);
    BOOST_CHECK_EQUAL(runs[1], 'C');
    BOOST_CHECK_EQUAL(runs[1].count(), 1);
    BOOST_CHECK_EQUAL(runs[2], 'G');
    BOOST_CHECK_EQUAL(runs[2].count(), 3);
    BOOST_CHECK_EQUAL(runs[3], 'T');
    BOOST_CHECK_EQUAL(runs[3].count(), 1);
    BOOST_CHECK_EQUAL(runs[4], 'A');
    BOOST_CHECK_EQUAL(runs[4].count(), 1);
}

BOOST_AUTO_TEST_SUITE_END();
