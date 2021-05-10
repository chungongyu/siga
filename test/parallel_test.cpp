#define BOOST_TEST_MODULE "siga.parallel"

#include <vector>
#include <map>

#include <boost/test/included/unit_test.hpp>
#include <boost/format.hpp>

#include "parallel_framework.h"

BOOST_AUTO_TEST_SUITE(parallel);

BOOST_AUTO_TEST_CASE(parallel_foreach) {
  {
    std::vector<int> x = {0, 1, 2, 4, 3, 5, 7, 6};
    int sum = 0;

    auto proc = [&](int tid, std::vector<int>::iterator i) {
        std::cout << boost::format("func: tid=%d, idx=%lu, val=%d") % tid % 0 % *i << std::endl;
        return *i;
      };
    auto postproc = [&](std::vector<int>::iterator, int x) {
        sum += x;
      };
    size_t threads = 0;
    parallel::foreach<int>(x.begin(), x.end(), proc, postproc, threads);
    std::cout << boost::format("result: sum=%d") % sum << std::endl;
    BOOST_CHECK(true);
  }
}

BOOST_AUTO_TEST_CASE(parallel_threads) {
  {
#ifdef _OPENMP
    std::cout << boost::format("threads: %lu") % parallel::threads(0) << std::endl;
#else
    std::cout << boost::format("threads: %lu") % parallel::threads(0) << std::endl;
#endif  // _OPENMP
    BOOST_CHECK(true);
  }
}

BOOST_AUTO_TEST_SUITE_END();
