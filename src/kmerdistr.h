#ifndef kmerdistr_h__
#define kmerdistr_h__

#include <map>

#include "fmindex.h"

class KmerDistribution {
 public:
  KmerDistribution() {
  }
  virtual ~KmerDistribution() {
  }

  static size_t sample(const FMIndex* index, size_t k, size_t n, KmerDistribution* distr);

  void add(int count) {
    auto i = _data.find(count);
    if (i != _data.end()) {
      ++i->second;
    } else {
      _data[count] = 1;
    }
  }
 private:
  // int -> size_t map of the number of times a kmer with multiplicty N has been seen
  std::map<int, size_t> _data;
};

#endif  // kmerdistr_h__
