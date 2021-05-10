#ifndef paralell_framework_h_
#define paralell_framework_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_OMP_H
#include <omp.h>
#endif

#include <cassert>
#include <vector>

namespace parallel {

template <typename Iterator, typename Func>
void foreach(Iterator begin, Iterator end, Func func, size_t threads, size_t chunk = 1) {
#ifdef _OPENMP
  if (threads != 1) {
    if (threads > 1) {
      omp_set_num_threads(threads);
    }
    threads = omp_get_num_threads();

    bool done = false;
    std::vector<Iterator> inputs;
    while (!done) {
      if (begin != end) {
        inputs.push_back(begin++);
      } else {
        done = true;
      }

      // Once all buffers are full or the input is finished, dispatch the work to the threads
      if (inputs.size() == threads*chunk || done) {
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < inputs.size(); ++i) {
          int tid = omp_get_thread_num();
          func(tid, inputs[i]);
        }
        inputs.clear();
      }
    }
  } else {
#endif  // _OPENMP
    while (begin != end) {  // Serial version or threads == 1
      func(0, begin++);
    }
#ifdef _OPENMP
  }
#endif  // _OPENMP
}

template <typename Iterator, typename Func>
void foreach(Iterator begin, Iterator end, Func func) {
  foreach<Iterator, Func>(begin, end, func, 0);
}

};  // namespace parallel

#endif  // paralell_framework_h_
