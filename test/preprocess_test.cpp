#define BOOST_TEST_MODULE "siga.preprocess"
#include <boost/test/included/unit_test.hpp>

#include "primer_screen.h"

BOOST_AUTO_TEST_SUITE(preprocess);

BOOST_AUTO_TEST_CASE(PrimerScreen_contains) {
    {
        std::string pcr_free_a = "AATGATACGGCGACCACCGAGATCTACA";
        BOOST_CHECK(PrimerScreen::containsPrimer(pcr_free_a));
    }
    {
        std::string pcr_free_b = "GATCGGAAGAGCGGTTCAGCAGGAATGC";
        BOOST_CHECK(PrimerScreen::containsPrimer(pcr_free_b));
    }
    {
        std::string pcr_free_b = "AGATCGGAAGAGCGGTTCAGCAGGAATGC";
        BOOST_CHECK(!PrimerScreen::containsPrimer(pcr_free_b));
    }
}

BOOST_AUTO_TEST_SUITE_END();
