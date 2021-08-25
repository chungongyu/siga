#ifndef correct_processor_h_
#define correct_processor_h_

#include <iostream>
#include <memory>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include "fmindex.h"
#include "kseq.h"
#include "suffix_array.h"

//
// Default params for corret command
//
#define kCorrectThreads         1
#define kCorrectBatchSize       1000
#define kCorrectAlgorithm       "kmer"
#define kCorrectKmerSize        31
#define kCorrectKmerThreshold   3
#define kCorrectKmerRounds      10
#define kCorrectKmerCountOffset 1

//
// CorrectProcessor
//
//
class CorrectProcessor {
 public:
  typedef boost::property_tree::ptree Options;
  struct Index {
    bool load(const std::string& prefix, const std::string& input, const Options& options);

    FMIndex fmi;
    std::shared_ptr<SuffixArray> sa;  // optional
    DNASeqList reads;  // optional
  };

  CorrectProcessor(const Options& options) : _options(options) {
  }

  bool process(const Index& index,
      DNASeqReader& reader, std::ostream& output) const;
  bool process(const Index& index,
      const std::string& input, const std::string& output) const;

 private:
  const Options& _options;
};

#endif  // correct_processor_h_
