#include "correct_processor.h"
#include "sequence_process_framework.h"
#include "utils.h"

#include <fstream>
#include <memory>

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
    }

    CorrectResult process(const SequenceProcessFramework::SequenceWorkItem& item) const {
        CorrectResult r;
        return r;
    }
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
    if (threads < 1) {
        std::shared_ptr< AbstractCorrector > proc(AbstractCorrector::create(_options));
        if (proc) {
            PostCorrector postproc(output);
            SequenceProcessFramework::SerialWorker<
                SequenceProcessFramework::SequenceWorkItem, 
                CorrectResult, 
                SequenceProcessFramework::WorkItemGenerator< SequenceProcessFramework::SequenceWorkItem >, 
                AbstractCorrector, 
                PostCorrector
                > worker;
            size_t num = worker.run(reader, proc.get(), &postproc);
            if (processed != NULL) {
                *processed = num;
            }
        }
    } else {
        assert(false);
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
