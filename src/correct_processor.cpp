#include "correct_processor.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <memory>
#include <vector>
#include <unordered_map>

#include <boost/format.hpp>

#include <log4cxx/logger.h>

#include "parallel_framework.h"
#include "utils.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.CorrectProcessor"));

//
// CorrectThreshold
//
class CorrectThreshold {
 public:
  static CorrectThreshold* get() {
    static CorrectThreshold inst;
    return &inst;
  }
  void minSupport(size_t minSupport) {
    _low = minSupport;
    _hight = minSupport + 1;
  }
  size_t requiredSupport(int phred) const {
    if (phred >= _cutoff) {
      return _hight;
    }
    return _low;
  }
 private:
  CorrectThreshold() : _cutoff(20) {
    minSupport(kCorrectKmerThreshold);
  }
  int _low;
  int _hight;
  int _cutoff;
};

//
// CorrectResult
//
struct CorrectResult {
  std::string seq;
  bool validQC;
};

class AbstractCorrector {
 public:
  virtual ~AbstractCorrector() {
  }
  virtual CorrectResult process(const DNASeqWorkItem& iterm) const = 0;

  static AbstractCorrector* create(const FMIndex& index, const CorrectProcessor::Options& options);
 protected:
  AbstractCorrector(const FMIndex& index, const CorrectProcessor::Options& options)
      : _index(index), _options(options) {
  }

  const CorrectProcessor::Options& _options;
  const FMIndex& _index;
};

class KmerCorrector : public AbstractCorrector {
 public:
  KmerCorrector(const FMIndex& index, const CorrectProcessor::Options& options) : AbstractCorrector(index, options) {
    _kmerSize = options.get<size_t>("kmer-size", kCorrectKmerSize);
    _maxAttempts = options.get<size_t>("kmer-rounds", kCorrectKmerRounds);
    _countOffset = options.get<size_t>("kmer-count-offset", kCorrectKmerCountOffset);
    CorrectThreshold::get()->minSupport(options.get<int>("kmer-threshold", kCorrectKmerThreshold));
  }

  CorrectResult process(const DNASeqWorkItem& item) const {
    CorrectResult r;

    // check read length
    if (item.read.seq.length() < _kmerSize) {
      r.seq = item.read.seq;
      r.validQC = false;
      return r;
    }

    std::string seq = item.read.seq;
    size_t k = _kmerSize;
    size_t n = seq.length();
    // For each kmer, calculate the minimum phred score seen in the bases
    // of the kmer
    std::vector<int> minPhredVector(n - k + 1);
    for (size_t i = k; i <= n; ++i) {
      int ps = std::numeric_limits<int>::max();
      for (size_t j = i - k; j < i; ++j) {
        ps = std::min(ps, item.read.score(j));
      }
      minPhredVector[i - k] = ps;
    }

    std::unordered_map<std::string, size_t> kmerCache;
    bool allSolid = false;
    size_t rounds = 0;
    bool done = false;
    while (!done) {
      // Compute the kmer counts across the read
      // and determine the positions in the read that are not covered by any solid kmers
      // These are the candidate incorrect bases
      std::vector<int> countVector(n - k + 1, 0);
      std::vector<int> solidVector(n, 0);

      for (size_t i = k; i <= n; ++i) {
        std::string kmer = seq.substr(i - k, k);

        // First check if this kmer is in the cache
        // If its not, find its count from the fm-index and cache it
        size_t count = 0;
        auto iter = kmerCache.find(kmer);
        if (iter != kmerCache.end()) {
          count = iter->second;
        } else {
          count = FMIndex::Interval::occurrences(kmer, &_index);
          kmerCache[kmer] = count;
        }

        // Get the phred score for the last base of the kmer
        int phred = minPhredVector[i - k];
        // Determine whether the base is solid or not based on phred scores
        if (count >= CorrectThreshold::get()->requiredSupport(phred)) {
          for (size_t j = 0; j < k; ++j) {
            solidVector[i - k + j] = 1;
          }
        }
      }

      allSolid = true;
      for (size_t i = 0; i < n; ++i) {
        LOG4CXX_DEBUG(logger, boost::format("Position[%d] = %d") % i % solidVector[i]);
        if (!solidVector[i]) {
          allSolid = false;
        }
      }
      LOG4CXX_DEBUG(logger, boost::format("Read %s solid = %d") % item.read.name % allSolid);
      // Stop if all kmers are well represented or we have exceeded the number of correction rounds
      if (allSolid || ++rounds > _maxAttempts) {
        break;
      }

      // Attempt to correct the leftmost potentially incorrect base
      bool corrected = false;
      for (size_t i = 0; i < n; ++i) {
        if (!solidVector[i]) {
          int phred = item.read.score(i);
          size_t threshold = CorrectThreshold::get()->requiredSupport(phred);

          // Attempt to correct the base using the leftmost covering kmer
          size_t leftIdx = (i + 1 >= _kmerSize ? i + 1 - _kmerSize : 0);
          if ((corrected = try2Correct(i, leftIdx, std::max(countVector[leftIdx] + _countOffset, threshold), seq))) {
            break;
          }
          // base was not corrected, try using the rightmost covering kmer
          size_t rightIdx = std::min(i, n - _kmerSize);
          if ((corrected = try2Correct(i, rightIdx, std::max(countVector[rightIdx] + _countOffset, threshold), seq))) {
            break;
          }
        }
      }

      // If no base in the read was corrected, stop the correction process
      if (!corrected) {
        assert(!allSolid);
        done = true;
      }
    }
    if (allSolid) {
      r.seq = seq;
      r.validQC = true;
    } else {
      r.seq = item.read.seq;
      r.validQC = false;
    }
    return r;
  }

 private:
  // Attempt to correct the base at position baseIdx in readSequence. Returns true if a correction was made
  // The correction is made only if the count of the corrected kmer is at least minCount
  bool try2Correct(size_t baseIdx, size_t kmerIdx, size_t minCount, std::string& read) const {
    assert(kmerIdx <= baseIdx && baseIdx < kmerIdx + _kmerSize);
    size_t deltaIdx = baseIdx - kmerIdx;
    char currBase = read[baseIdx];
    std::string kmer = read.substr(kmerIdx, _kmerSize);
    size_t bestCount = 0;
    char bestBase = '$';

    LOG4CXX_DEBUG(logger, boost::format("baseIdx: %d kmerIdx: %d %s %s")
        % baseIdx % kmerIdx % kmer % make_dna_reverse_complement_copy(kmer));

    for (size_t i = 0; i < DNAAlphabet::size; ++i) {
      char c = DNAAlphabet::DNA[i];
      if (c != currBase) {
        kmer[deltaIdx] = c;
        size_t count = FMIndex::Interval::occurrences(kmer, &_index);
        LOG4CXX_DEBUG(logger, boost::format("%c %lu") % c % count);
        if (count >= minCount) {
          if (bestBase != '$') {
            return false;
          }
          bestBase = c;
          bestCount = count;
        }
      }
    }
    if (bestCount >= minCount) {
      assert(bestBase != '$');
      read[baseIdx] = bestBase;
      return true;
    }
    return false;
  }

  size_t _kmerSize;
  size_t _maxAttempts;
  size_t _countOffset;
};

class OverlapCorrector : public AbstractCorrector {
 public:
  OverlapCorrector(const FMIndex& index, const CorrectProcessor::Options& options)
      : AbstractCorrector(index, options) {
  }
  CorrectResult process(const DNASeqWorkItem& item) const {
    CorrectResult r;
    return r;
  }
};

AbstractCorrector* AbstractCorrector::create(const FMIndex& index, const CorrectProcessor::Options& options) {
  std::string algorithm = options.get<std::string>("algorithm", kCorrectAlgorithm);
  if (algorithm == "kmer") {
    return new KmerCorrector(index, options);
  } else if (algorithm == "overlap") {
    return new OverlapCorrector(index, options);
  }
  return NULL;
}

class PostCorrector {
 public:
  PostCorrector(std::ostream& stream, std::ostream* discard = nullptr)
      : _stream(stream), _discard(discard) {
  }
  void operator()(const DNASeqWorkItem& workItem, const CorrectResult& result) {
    if (result.validQC) {
      DNASeq read = workItem.read;
      read.seq = result.seq;
      _stream << read;
    }
  }

 private:
  std::ostream& _stream;
  std::ostream* _discard;
};

bool CorrectProcessor::process(const FMIndex& index, DNASeqReader& reader, std::ostream& output,
      size_t threads, size_t* processed) const {
  {
    threads = parallel::threads(threads);
    DNASeqWorkItemGenerator<DNASeqWorkItem> generator(reader);
    std::vector<AbstractCorrector *> proclist(threads);
    for (size_t i = 0; i < threads; ++i) {
      proclist[i] = AbstractCorrector::create(index, _options);
      if (proclist[i] == nullptr) {
        return false;
      }
    }

    PostCorrector postproc(output);
    parallel::foreach<DNASeqWorkItem, CorrectResult>(generator,
        [&](int tid, const DNASeqWorkItem& it) {
          assert(tid < threads);
          return proclist[tid]->process(it);
        }, postproc, threads);

    for (size_t i = 0; i < threads; ++i) {
      delete proclist[i];
    }
  }
  return true;
}

bool CorrectProcessor::process(const FMIndex& index, const std::string& input, const std::string& output,
      size_t threads, size_t* processed) const {
  // DNASeqReader
  std::shared_ptr<std::istream> reads(Utils::ifstream(input));
  if (!reads) {
    LOG4CXX_ERROR(logger, boost::format("Failed to read file %s") % input);
    return false;
  }
  std::shared_ptr<DNASeqReader> reader(DNASeqReaderFactory::create(*reads));
  if (!reader) {
    LOG4CXX_ERROR(logger, boost::format("Failed to create DNASeqReader %s") % input);
    return false;
  }

  std::shared_ptr<std::ostream> out(Utils::ofstream(output));
  if (!out) {
    LOG4CXX_ERROR(logger, boost::format("Failed to create DNASeqWriter %s") % output);
    return false;
  }
  return process(index, *reader, *out, threads, processed);
}
