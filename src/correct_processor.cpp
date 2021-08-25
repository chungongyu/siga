#include "correct_processor.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <set>
#include <map>
#include <memory>
#include <vector>
#include <unordered_map>

#include <boost/format.hpp>

#include <log4cxx/logger.h>
#include <edlib.h>

#include "constant.h"
#include "kseq.h"
#include "suffix_array.h"
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
  std::string cigar;
};

class AbstractCorrector {
 public:
  virtual ~AbstractCorrector() {
  }
  virtual CorrectResult process(const DNASeqWorkItem& iterm) const = 0;

  static AbstractCorrector* create(const CorrectProcessor::Index& index, const CorrectProcessor::Options& options);
 protected:
  AbstractCorrector(const CorrectProcessor::Index& index, const CorrectProcessor::Options& options)
      : _index(index), _options(options) {
  }

  const CorrectProcessor::Options& _options;
  const CorrectProcessor::Index& _index;
};

class KmerCorrector : public AbstractCorrector {
 public:
  KmerCorrector(const CorrectProcessor::Index& index, const CorrectProcessor::Options& options) : AbstractCorrector(index, options) {
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
      std::vector<size_t> countVector(n - k + 1, 0);
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
          count = FMIndex::Interval::occurrences(kmer, &_index.fmi);
          kmerCache[kmer] = count;
        }

        countVector[i - k] = count;

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
        size_t count = FMIndex::Interval::occurrences(kmer, &_index.fmi);
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
  OverlapCorrector(const CorrectProcessor::Index& index, const CorrectProcessor::Options& options)
      : AbstractCorrector(index, options) {
    double errorRate = _options.get<double>("error-rate", 0.001);
    // _kmerSize = std::max((size_t)kCorrectKmerSize, (size_t)(log(0.75)/log(1.0-errorRate)));
    _kmerSize = options.get<size_t>("kmer-size", kCorrectKmerSize);
    assert(_kmerSize > 1);
    _kmerThreshold = options.get<size_t>("kmer-threshold", kCorrectKmerThreshold);
    LOG4CXX_INFO(logger, boost::format("[OverlapCorrector] error-rate: %lf kmer-size: %lu kmer-threshold: %lu") % errorRate % _kmerSize % _kmerThreshold);
    _alignConfig = edlibNewAlignConfig(3*_kmerSize*errorRate*(1-errorRate), EDLIB_MODE_HW, EDLIB_TASK_PATH, nullptr, 0);
  }

  CorrectResult process(const DNASeqWorkItem& item) const {
    CorrectResult r = {item.read.seq, true};

    if (item.read.seq.length() >= _kmerSize) {
      for (size_t loop = 0; loop < 10; ++loop) {
        bool validQC = true;

      size_t i = 0;
      while (i < r.seq.length()) {
        int delta = 0; // number of insertion(+) and deletion(-) corrections;
        size_t k = std::min(r.seq.length() - _kmerSize, i);
        std::string kmer = r.seq.substr(k, _kmerSize);
        size_t count = FMIndex::Interval::occurrences(kmer, &_index.fmi);
        if (count <= _kmerThreshold) {
          LOG4CXX_INFO(logger, boost::format("[OverlapCorrector] idx=%lu name=%s length=%lu kmer=%lu count=%lu") % item.idx % item.read.name % item.read.seq.length() % k % count);
          std::vector<KmerAlign> alignments;
          std::function<void(const std::string&)> inserter = [&](const std::string& fragment) {
                auto r = edlibAlign(kmer.c_str(), kmer.length(), fragment.c_str(), fragment.length(), &_alignConfig);
                if (r.status == EDLIB_STATUS_OK && 0 < r.editDistance && r.editDistance <= _alignConfig.k) {
                  alignments.emplace_back(fragment, r);

                  char* cigar = edlibAlignmentToCigar(r.alignment, r.alignmentLength, EDLIB_CIGAR_EXTENDED);
                  LOG4CXX_INFO(logger, boost::format("[OverlapCorrector] idx=%lu kmer=%lu count=%lu align: ed=%d locations=%d len=%d cigar=%s") % item.idx % k % count % r.editDistance % r.numLocations % r.alignmentLength % (cigar ? cigar : ""));
                  if (cigar) {
                    free(cigar);
                  }
                }
              };
          KmerMatch::get(this, r.seq, item.idx, k, inserter);
          if (alignments.size() == 1) {
            delta = alignments[0].correct(this, &r.seq, k);
            LOG4CXX_INFO(logger, boost::format("[OverlapCorrector] kmer=%s, target=%s, correct=%s") % kmer % alignments[0].target % r.seq.substr(k, _kmerSize));
          } else {
            validQC = false;
          }
          LOG4CXX_INFO(logger, boost::format("[OverlapCorrector] idx=%lu kmer=%lu count=%lu match: %lu") % item.idx % k % count % alignments.size());
          for (auto& r : alignments) {
            edlibFreeAlignResult(&r.align);
          }
        }
        i += _kmerSize - 1 + delta;
      }

        if (validQC) {
          break;
        }
      }
    }

    return r;
  }

 private:
  struct KmerAlign {
    KmerAlign() {}
    KmerAlign(const std::string& t, const EdlibAlignResult& r) : target(t), align(r) {
    }

    int correct(const OverlapCorrector* ec, std::string* query, size_t k) {
      int delta = 0;

      assert(align.alignment);
      assert(align.numLocations >= 1);
      size_t m = align.startLocations[0];
      size_t qidx = k, tidx = m;
      int istart = 0, iend = align.alignmentLength - 1;

      for (; istart < iend && align.alignment[istart] == EDLIB_EDOP_INSERT; ++istart) {++qidx;}
      // for (; istart < iend && align.alignment[istart] == EDLIB_EDOP_DELETE; ++istart) {++tidx;}
      // for (; istart < iend && align.alignment[iend] == EDLIB_EDOP_INSERT; --iend) {}
      for (; istart < iend && align.alignment[iend] == EDLIB_EDOP_DELETE; --iend) {}

      // Query #0 (8 residues): score = 2
      // T: AA--GCTA (0 - 5)
      //    ||  ||||
      // Q: AAAAGCTA (0 - 7)
	  //
	  // FIX:
      // T: --AAGCTA (0 - 5)
      //      ||||||
      // Q: AAAAGCTA (0 - 7)
      {
        int state = 0;
        size_t ins = 0;
        for (int i = istart; i <= iend; ++i) {
          if (state == 0) {
            if (align.alignment[i] == EDLIB_EDOP_INSERT) {
              state = 1;
            } else if (align.alignment[i] != EDLIB_EDOP_MATCH) {
              break;
            }
          }
          if (state == 1) {
            if (align.alignment[i] != EDLIB_EDOP_INSERT) {
              if (query->substr(qidx + ins, i - istart - ins) == target.substr(tidx, i - istart - ins)) {
                for (size_t j = 0; j < ins; ++j) {
                  std::swap(align.alignment[istart + j], align.alignment[i - ins + j]);
                }
              }
              break;
            }
            ++ins;
          }
        }
      }
      for (; istart < iend && align.alignment[istart] == EDLIB_EDOP_INSERT; ++istart) {++qidx;}

      for (int i = istart; i <= iend; ++i) {
        if (align.alignment[i] == EDLIB_EDOP_INSERT) {
          query->erase(qidx, 1);
          --delta;
        } else if (align.alignment[i] == EDLIB_EDOP_DELETE) {
          query->insert(qidx++, 1, target[tidx++]);
          ++delta;
        } else if (align.alignment[i] == EDLIB_EDOP_MISMATCH) {
          (*query)[qidx++] = target[tidx++];
        } else {
          assert(align.alignment[i] == EDLIB_EDOP_MATCH);
          ++qidx; ++tidx;
        }
      }

      return delta;
    }

    std::string target; // the location in the query sequence of this kmer
    EdlibAlignResult align;
  };

  // Struct to hold a partial match in the FM-index
  struct KmerMatch {
    size_t kmer; // the location in the query sequence of this kmer
    size_t idx;  // an index into the BWT.
    size_t query;  // the location in the query sequence of the solid kmer
    size_t target;  // the location in the target sequence of the solid kmer
    bool isRC;  // flag indicates the strand of the partial match

    bool operator<(const KmerMatch& o) const {
      if (idx != o.idx) {
        return idx < o.idx;
      }
      if (kmer != o.kmer) {
        return kmer < o.kmer;
      }
      return isRC < o.isRC;
    }
    bool operator==(const KmerMatch& o) const {
      return idx == o.idx && kmer == o.kmer &&
          isRC == o.isRC;
    }

    std::string getString(const OverlapCorrector* ec, const std::string& read) const {
      std::string seq;
      if (ec->_index.sa) {
        const std::string& r = ec->_index.reads[idx].seq;
        if (query < kmer) {
          size_t lr = kmer - query;
          size_t ll = std::min(std::min(query, target), ec->_kmerSize - lr);
          if (read.substr(query - ll, ll) == r.substr(target - ll, ll)) {
            seq = r.substr(target + lr, ec->_kmerSize);
          }
        } else if (query > kmer) {
          size_t ll = query - kmer;
          size_t lr = std::min(
              std::min(read.length() - query, r.length() - target) - ec->_kmerSize,
              ec->_kmerSize - ll);
          if (read.substr(query + ec->_kmerSize, lr) == r.substr(target + ec->_kmerSize, lr)) {
            seq = r.substr(target - std::min(target, ll), std::min(target, ll) + lr);
          }
        }
      }
      return seq;
    }

    static void get(const OverlapCorrector* ec, const std::string& read, size_t idx, size_t i,
        std::function<void(const std::string&)>& callback) {
      std::vector<std::pair<std::string, size_t> > fragments;
      auto is_overlap = [](const std::string& x, const std::string& y, size_t distance = 5) {
            for (size_t d = 0; d <= distance; ++d) {
              size_t i = d, j = 0;
              while (i < x.length() && j < y.length() && x[i++] == y[j++]) {}
              if (i == x.length() || j == y.length()) {
                return true;
              }
            }
            return false;
          };
      std::function<void(const std::string&, size_t)> inserter = [&](const std::string& fragment, size_t coverage) {
            auto it = std::find_if(fragments.begin(), fragments.end(), [&](const std::pair<std::string, size_t>& o) {
                  return is_overlap(o.first, fragment) || is_overlap(fragment, o.first);
                });
            if (it != fragments.end()) {
              if (it->second < coverage || it->first.length() < fragment.length()) {
                it->first = fragment;
              }
              it->second += coverage;
            } else {
              fragments.emplace_back(fragment, coverage);
              callback(fragment);
            }
          };
      KmerMatch::get(ec, read, idx, i, inserter);
    }

    static void get(const OverlapCorrector* ec, const std::string& read, size_t idx, size_t i,
        std::function<void(const std::string&, size_t)>& callback) {
      std::set<KmerMatch> matches;
      std::function<void(const KmerMatch&)> inserter = [&](const KmerMatch& match) {
          if (match.idx != idx) {
            matches.insert(match);
          }};

      KmerMatch::match(ec, read, i, inserter);

      // Refine the matches by computing proper overlaps between the sequences
      // Use the overlaps that meet the thresholds to build a multiple alignment
      std::unordered_map<std::string, size_t> fragments;
      for (const auto& match : matches) {
        std::string seq = match.getString(ec, read);
        if (seq.length() == ec->_kmerSize) {
          ++fragments[seq];
          //LOG4CXX_DEBUG(logger, boost::format("[OverlapCorrector] idx=%lu kmer=%lu match: query=%lu, target=%lu, idx=%lu, isRC=%d, seq=%s, seq=%s") % idx % k % match.query % match.target % match.idx % match.isRC % kmer % seq);
        }
      }
      for (const auto& it : fragments) {
        if (it.second > ec->_kmerThreshold) {
          //LOG4CXX_DEBUG(logger, boost::format("[OverlapCorrector] idx=%lu kmer=%lu match: seq=%s, seq=%s coverage=%lu") % idx % i % read.substr(i, ec->_kmerSize) % it.first % it.second);
          callback(it.first, it.second);
        }
      }
    }

    static void match(const OverlapCorrector* ec, const std::string& seq, size_t i,
        std::function<void(const KmerMatch&)>& callback) {
      size_t start = i - std::min(i, ec->_kmerSize), end = std::min(i + 2*ec->_kmerSize, seq.length());
      assert(ec->_kmerSize <= end - start);

      std::map<KmerMatch, bool> matches;
      auto matcher = [&](const std::string& kmer, size_t query, bool isRC) {
          FMIndex::Interval interval = FMIndex::Interval::get(kmer, &ec->_index.fmi);
          if (interval.valid() && interval.length() > ec->_kmerThreshold) {
            for (size_t idx = interval.lower; idx <= interval.upper; ++idx) {
              KmerMatch match = {i, idx, query, 0, isRC};
              matches.emplace(match, false);
            }
          }
        };

      for (size_t query = start; query <= end - ec->_kmerSize; ++query) {
        std::string kmer = seq.substr(query, ec->_kmerSize);
        matcher(kmer, query, false);  // forward
        // make_dna_reverse_complement(kmer);
        // matcher(kmer, query, true);  // opposite
      }

      // Backtrack through the kmer indices to turn them into read indices.
      for (auto it = matches.begin(); it != matches.end(); ++it) {
        if (it->second) {
          continue;  // This index has been visited
        }
        it->second = true;  // Mark this as visited

        KmerMatch match = it->first;
        // Backtrack the index until we hit the starting symbol
        while (true) {
          char c = ec->_index.fmi.getChar(match.idx);
          match.idx = ec->_index.fmi.getPC(c) + ec->_index.fmi.getOcc(c, match.idx - 1);
          // Check if the hash indicates we have visited this index. If so, stop the backtrack
          auto k = matches.find(match);
          if (k != matches.end()) {
            if (k->second) {  // We have processed this index already
              break;
            }
            k->second = true;
          }
          if (c == '$') {
            // We've found the lexicographic index for this read. Turn it into a proper ID
            if (ec->_index.sa) {
              assert(match.idx < ec->_index.sa->size());
              match.idx = (*ec->_index.sa)[match.idx].i;
              callback(match);
            }
            break;
          }
          ++match.target;
        }
      }
    }
  };

  size_t _kmerSize;
  size_t _kmerThreshold;
  EdlibAlignConfig _alignConfig;
};

bool CorrectProcessor::Index::load(const std::string& prefix, const std::string& input, const CorrectProcessor::Options& options) {
  if (!FMIndex::load(prefix + BWT_EXT, fmi)) {
    LOG4CXX_ERROR(logger, boost::format("failed to load FMIndex from %s") % prefix);
    return false;
  }
  std::string algorithm = options.get<std::string>("algorithm", kCorrectAlgorithm);
  if (algorithm == "overlap") {
    if (!(sa = std::shared_ptr<SuffixArray>(SuffixArray::load(prefix + SAI_EXT)))) {
      LOG4CXX_ERROR(logger, boost::format("failed to load suffix array from %s") % prefix);
      return false;
    }
    if (!ReadDNASequences(input, reads, (~kSeqWithComment)&(~kSeqWithComment))) {
      LOG4CXX_ERROR(logger, boost::format("failed to read sequences from %s") % input);
      return false;
    }
  }
  return true;
}

AbstractCorrector* AbstractCorrector::create(const CorrectProcessor::Index& index, const CorrectProcessor::Options& options) {
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
      if (!result.cigar.empty()) {
        std::string sep;
        if (!read.comment.empty()) {
          sep = " ";
        }
        read.comment += boost::str(boost::format("%sCG:Z:%s") % sep % result.cigar);
      }
      read.seq = result.seq;
      read.quality.clear();
      _stream << read;
    }
  }

 private:
  std::ostream& _stream;
  std::ostream* _discard;
};

bool CorrectProcessor::process(const Index& index, DNASeqReader& reader, std::ostream& output) const {
  {
    size_t threads = parallel::threads(_options.get<size_t>("threads", kCorrectThreads));
    size_t batch = _options.get<size_t>("batch-size", kCorrectBatchSize);

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
        }, postproc, threads, batch);

    for (size_t i = 0; i < threads; ++i) {
      delete proclist[i];
    }
  }
  return true;
}

bool CorrectProcessor::process(const Index& index, const std::string& input, const std::string& output) const {
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
  return process(index, *reader, *out);
}
