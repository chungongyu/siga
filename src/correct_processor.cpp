#include "correct_processor.h"
#include "sequence_process_framework.h"
#include "utils.h"

#include <fstream>
#include <memory>
#include <unordered_map>

#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.CorrectProcessor"));

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
    virtual CorrectResult process(const SequenceProcessFramework::SequenceWorkItem& iterm) const = 0;

    static AbstractCorrector* create(const CorrectProcessor::Options& options);
protected:
    AbstractCorrector(const CorrectProcessor::Options& options) : _options(options) {
    }

    const CorrectProcessor::Options& _options;
};

class KmerCorrector : public AbstractCorrector {
public:
    KmerCorrector(const CorrectProcessor::Options& options) : AbstractCorrector(options) {
        _kmerSize = options.get< size_t >("kmer-size", 31);
        _maxAttempts = options.get< size_t >("kmer-rounds", 10);
    }

    CorrectResult process(const SequenceProcessFramework::SequenceWorkItem& item) const {
        CorrectResult r;

        // check read length
        if (item.read.seq.length() < _options.get< size_t >("xx", -1)) {
            r.seq = item.read.seq;
            r.validQC = false;
            return r;
        }

        size_t k = _options.get< size_t >("xx", 17);
        size_t n = item.read.seq.length();
        if (k <= n) {
            // For each kmer, calculate the minimum phred score seen in the bases
            // of the kmer
            std::vector< int > minPhredVector(n - k + 1);
            for (size_t i = k; i <= n; ++i) {
                int ps = std::numeric_limits< int >::max();
                for (size_t j = i - k; j < i; ++j) {
                    ps = std::min(ps, item.read.score(j));
                }
                minPhredVector[i - k] = ps;
            }

            std::unordered_map< std::string, size_t > kmerCache;
            size_t rounds = 0;
            bool done = false;
            while (!done) {
                // Compute the kmer counts across the read
                // and determine the positions in the read that are not covered by any solid kmers
                // These are the candidate incorrect bases
                std::vector<int> countVector(n - k + 1, 0);
                std::vector<int> solidVector(n, 0);

                for (size_t i = k; i <= n; ++i) {
                    std::string kmer = item.read.seq.substr(i - k, k);

                    // First check if this kmer is in the cache
                    // If its not, find its count from the fm-index and cache it
                    size_t count = 0;
                    std::unordered_map< std::string, size_t >::iterator iter = kmerCache.find(kmer);
                    if (iter != kmerCache.end()) {
                        count = iter->second;
                    } else {
                        kmerCache[kmer] = count;
                    }

                    // Get the phred score for the last base of the kmer
                    int phred = minPhredVector[i - k];
                    // Determine whether the base is solid or not based on phred scores
                    if (count >= _kmerThreshold) {
                        for (size_t j = 0; j < k; ++j) {
                            solidVector[i - k + j] = 1;
                        }
                    }
                }

                bool allSolid = true;
                for (size_t i = 0; i < n; ++i) {
                    LOG4CXX_DEBUG(logger, boost::format("Position[%d] = %d") % i % solidVector[i]);
                    if (!solidVector[i]) {
                        allSolid = false;
                    }
                }
                LOG4CXX_DEBUG(logger, boost::format("Read %s is solid = %d") % item.read.name % allSolid);
                // Stop if all kmers are well represented or we have exceeded the number of correction rounds
                if(allSolid || ++rounds > _maxAttempts) {
                    break;
                }

                // Attempt to correct the leftmost potentially incorrect base
                bool corrected = false;
                for (size_t i = 0; i < n; ++i) {
                    if (!solidVector[i]) {
                        // Attempt to correct the base using the leftmost covering kmer
                        int phred = item.read.score(i);
                    }
                }
            }
        }
        return r;
    }

private:
    // Attempt to correct the base at position idx in readSequence. Returns true if a correction was made
    // The correction is made only if the count of the corrected kmer is at least minCount
    bool try2Correct(size_t baseIdx, size_t kmerIdx, size_t minCount, std::string& read) {
        assert(kmerIdx <= baseIdx && baseIdx < kmerIdx + _kmerSize);
        return false;
    }

    size_t _kmerSize;
    int _kmerThreshold;
    size_t _maxAttempts;
};

AbstractCorrector* AbstractCorrector::create(const CorrectProcessor::Options& options) {
    std::string algorithm = options.get< std::string >("algorithm", "kmer");
    if (algorithm == "kmer") {
        return new KmerCorrector(options);
    }
    return NULL;
}

class PostCorrector {
public:
    PostCorrector(std::ostream& stream, std::ostream* discard=NULL) : _stream(stream), _discard(discard) {
    }
    void process(const SequenceProcessFramework::SequenceWorkItem& workItem, const CorrectResult& result) {
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

bool CorrectProcessor::process(DNASeqReader& reader, std::ostream& output, size_t threads, size_t* processed) const {
    if (threads <= 1) {
        std::shared_ptr< AbstractCorrector > proc(AbstractCorrector::create(_options));
        if (!proc) {
            return false;
        }

        PostCorrector postproc(output);
        SequenceProcessFramework::SerialWorker<
            SequenceProcessFramework::SequenceWorkItem, 
            CorrectResult, 
            SequenceProcessFramework::SequenceWorkItemGenerator< SequenceProcessFramework::SequenceWorkItem >, 
            AbstractCorrector, 
            PostCorrector
            > worker;
        size_t num = worker.run(reader, proc.get(), &postproc);
        if (processed != NULL) {
            *processed = num;
        }
        return true;
    } else {
#ifdef _OPENMP
        std::vector< AbstractCorrector* > proclist(threads);
        for (size_t i = 0; i < threads; ++i) {
            proclist[i] = AbstractCorrector::create(_options);
            if (proclist[i] == NULL) {
                return false;
            }
        }

        PostCorrector postproc(output);
        SequenceProcessFramework::ParallelWorker<
            SequenceProcessFramework::SequenceWorkItem, 
            CorrectResult, 
            SequenceProcessFramework::SequenceWorkItemGenerator< SequenceProcessFramework::SequenceWorkItem >, 
            AbstractCorrector, 
            PostCorrector
            > worker;
        size_t num = worker.run(reader, &proclist, &postproc);
        if (processed != NULL) {
            *processed = num;
        }

        for (size_t i = 0; i < threads; ++i) {
            delete proclist[i];
        }
        return true;
#else
        LOG4CXX_ERROR(logger, "failed to load OpenMP");
        assert(false);
#endif // _OPENMP
    }
    return false;
}

bool CorrectProcessor::process(const std::string& input, const std::string& output, size_t threads, size_t* processed) const {

    // DNASeqReader
    std::ifstream reads(input);
    std::shared_ptr< DNASeqReader > reader(DNASeqReaderFactory::create(reads, &input));
    if (!reader) {
        LOG4CXX_ERROR(logger, boost::format("Failed to create DNASeqReader %s") % input);
        return false;
    }

    std::ofstream out(output);
    return process(*reader, out, threads, processed);
}
