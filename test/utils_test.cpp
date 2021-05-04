#define BOOST_TEST_MODULE "siga.utils"
#include <boost/test/included/unit_test.hpp>

#include "utils.h"

BOOST_AUTO_TEST_SUITE(utils);

BOOST_AUTO_TEST_CASE(Utils_rand) {
  Utils::srand();
  {
    int v = Utils::rand();
    BOOST_CHECK_GE(v, 0);
    BOOST_CHECK_LE(v, RAND_MAX);
  }
  {
    int v = Utils::rand(10);
    BOOST_CHECK_GE(v, 0);
    BOOST_CHECK_LT(v, 10);
  }
  {
    int v = Utils::rand(1, 10);
    BOOST_CHECK_GE(v, 1);
    BOOST_CHECK_LT(v, 10);
  }
  {
    int v = Utils::rand(9, 10);
    BOOST_CHECK_GE(v, 9);
    BOOST_CHECK_LT(v, 10);
  }
}

BOOST_AUTO_TEST_CASE(Utils_stem) {
  BOOST_CHECK_EQUAL(Utils::stem("a.txt"), "a");
  BOOST_CHECK_EQUAL(Utils::stem("a.txt.gz"), "a");
  BOOST_CHECK_EQUAL(Utils::stem("a.txt.bz2"), "a");
}

BOOST_AUTO_TEST_SUITE_END();
