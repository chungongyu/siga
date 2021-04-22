#define BOOST_TEST_MODULE "siga.overlap"
#include <boost/test/included/unit_test.hpp>

#include "asqg.h"
#include "overlap_builder.h"

BOOST_AUTO_TEST_SUITE(overlap);

BOOST_AUTO_TEST_CASE(ASQG_fmt) {
  {
    ASQG::IntTagValue tag;
    BOOST_CHECK(tag.fromstring("test:i:100"));
    BOOST_CHECK((bool)tag); // initialized
    BOOST_CHECK_EQUAL(100, (int)tag);
  }
  {
    ASQG::FloatTagValue tag;
    BOOST_CHECK(tag.fromstring("test:f:100.0"));
    BOOST_CHECK((bool)tag); // initialized
    BOOST_CHECK_CLOSE(100.0, (float)tag, 1e-5);
  }
  {
    ASQG::StringTagValue tag;
    BOOST_CHECK(tag.fromstring("test:Z:100.0"));
    BOOST_CHECK((bool)tag); // initialized
    BOOST_CHECK_EQUAL("100.0", (std::string)tag);
  }
}

BOOST_AUTO_TEST_CASE(OverlapBuilder) {
  int i = 1;
  BOOST_CHECK(i);
}

BOOST_AUTO_TEST_SUITE_END();
