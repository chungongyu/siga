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

template <typename Input, typename Output, typename Generator, typename Processor, typename PostProcessor>
void foreach(Generator& generator, Processor proc, PostProcessor postproc,
    size_t threads = 0, size_t batch = 1000) {
  Input input;
#ifdef _OPENMP
  if (threads != 1) {
    if (threads > 1) {
      omp_set_num_threads(threads);
    }
    bool done = false;
    std::vector<Input> inputs;
    while (!done) {
      if (generator.generate(input)) {
        inputs.push_back(input);
      } else {
        done = true;
      }

      // Once all buffers are full or the input is finished, dispatch the work to the threads
      if (!inputs.empty() && (inputs.size() >= threads*batch || done)) {
        std::vector<Output> outputs(inputs.size());

        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < inputs.size(); ++i) {
          outputs[i] = proc(omp_get_thread_num(), inputs[i]);
        }

        for (size_t i = 0; i < inputs.size(); ++i) {
          postproc(inputs[i], outputs[i]);
        }

        inputs.clear();
      }
    }
  } else {
#endif  // _OPENMP
    while (generator.generate(input)) {  // Serial version or threads == 1
      Output output = proc(0, input);
      postproc(input, output);
    }
#ifdef _OPENMP
  }
#endif  // _OPENMP
}

template <typename Iterator>
struct IteratorGenerator {
 IteratorGenerator(Iterator begin, Iterator end) : _begin(begin), _end(end) {
 }
 bool generate(Iterator& it) {
   if (_begin != _end) {
     it = _begin++;
     return true;
   }
   return false;
 }
 private:
  Iterator _begin;
  Iterator _end;
};

template <typename Output, typename Iterator, typename Processor, typename PostProcessor>
void foreach(Iterator begin, Iterator end, Processor proc, PostProcessor postproc,
    size_t threads = 0, size_t batch = 1000) {
  IteratorGenerator<Iterator> generator(begin, end);
  foreach<Iterator, Output, IteratorGenerator<Iterator>, Processor, PostProcessor>(
      generator, proc, postproc, threads, batch);
}

inline size_t threads(size_t threads) {
  if (threads != 1) {
#ifdef _OPENMP
    if (threads > 1) {
      return threads;
    }
    return omp_get_max_threads();
#endif  // _OPENMP
  }
  return 1;
}

};  // namespace parallel

#endif  // paralell_framework_h_
