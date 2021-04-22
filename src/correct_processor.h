#ifndef correct_processor_h_
#define correct_processor_h_

#include <iostream>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include "fmindex.h"
#include "kseq.h"

//
// Default params for corret command
//
#define kCorrectThreads         1
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

  CorrectProcessor(const Options& options) : _options(options) {
  }

  bool process(const FMIndex& index,
      DNASeqReader& reader, std::ostream& output,
      size_t threads = 1, size_t* processed = nullptr) const;
  bool process(const FMIndex& index,
      const std::string& input, const std::string& output,
      size_t threads = 1, size_t* processed = nullptr) const;

 private:
  const Options& _options;
};

#endif  // correct_processor_h_
