#define BOOST_TEST_MODULE "siga.parallel"

#include <vector>
#include <map>

#include <boost/test/included/unit_test.hpp>
#include <boost/format.hpp>

#include "parallel_framework.h"

BOOST_AUTO_TEST_SUITE(parallel);

BOOST_AUTO_TEST_CASE(parallel_foreach) {
  {
    std::vector<int> x = {0, 1, 2, 4, 3};
    
    std::map<int, std::vector<int> > r;
    auto func = [&](int tid, std::vector<int>::iterator i) {
        std::cout << boost::format("func: tid=%d, idx=%lu, val=%d") % tid % 0 % *i << std::endl;
        auto& k = r[tid];
        k.push_back(*i);
      };
    size_t threads = 3;
    parallel::foreach(x.begin(), x.end(), [&](int tid, std::vector<int>::iterator i) {
          std::cout << boost::format("func: tid=%d, val=%d") % tid % *i << std::endl;
        }, threads);
    for (auto& i : r) {
      std::cout << boost::format("result: tid=%d, val_list: %s") % i.first % "" << std::endl;
    }
    BOOST_CHECK(true);
  }
}

BOOST_AUTO_TEST_SUITE_END();
