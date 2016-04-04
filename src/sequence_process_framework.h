#ifndef sequence_process_framework_h_
#define sequence_process_framework_h_

#include "kseq.h"

#include <fstream>
#include <iostream>
#include <memory>

#include <boost/format.hpp>

#include <log4cxx/logger.h>

namespace SequenceProcessFramework {
    //
    // SequenceWorkItem - Definition of the data structure used in the generic
    // functions to process a sequence
    //
    struct SequenceWorkItem {
        SequenceWorkItem() : idx(0) {
        }
        SequenceWorkItem(size_t idx, const DNASeq& read) : idx(idx), read(read) {
        }
        size_t idx;
        DNASeq read;
    };

    // Genereic class to generate work items using a seq reader
    template< class Input >
    class WorkItemGenerator {
    public:
        WorkItemGenerator(DNASeqReader& reader) : _reader(reader), _consumed(0) {
        }

        bool generate(SequenceWorkItem& item) {
            if (_reader.read(item.read)) {
                item.idx = _consumed;
                ++_consumed;
                return true;
            }
            return false;
        }
    
    private:
        DNASeqReader& _reader;
        size_t _consumed;
    };

    // Generic function to process n work items from a file. 
    // With the default value of -1, n becomes the largest value representable for
    // a size_t and all values will be read
    template< class Input, class Output, class Generator, class Processor, class PostProcessor >
    class SerialWorker {
    public:
        size_t run(Generator& generator, Processor* proc, PostProcessor* postproc, size_t n = -1) {
            Input workItem;
            while (generator.consumed() < n && generator.generate(workItem)) {
                Output output = proc->process(workItem);
                postproc->process(workItem, output);
            }
            
            LOG4CXX_INFO(logger, boost::format("processed %d sequences") % generator.consumed());

            return generator.consumed();
        }

        size_t run(DNASeqReader& reader, Processor* proc, PostProcessor* postproc, size_t n = -1) {
            WorkItemGenerator< Input > generator(reader);
            return run(generator, proc, postproc, n);
        }

        size_t run(std::istream& stream, Processor* proc, PostProcessor* postproc, size_t n = -1) {
            std::shared_ptr< DNASeqReader > reader(DNASeqReaderFactory::create(stream));
            if (reader) {
                return run(*reader, proc, postproc, n);
            }
            return 0;
        }
        size_t run(const std::string& filename, Processor* proc, PostProcessor* postproc, size_t n = -1) {
            std::ofstream stream(filename.c_str());
            return run(stream, proc, postproc, n);
        }

    private:
        static log4cxx::LoggerPtr logger;
    };

    template< class Input, class Output, class Generator, class Processor, class PostProcessor >
    log4cxx::LoggerPtr SerialWorker< Input, Output, Generator, Processor, PostProcessor >::logger(log4cxx::Logger::getLogger("arcs.SequenceProcessFramework"));
};

#endif // sequence_process_framework_h_
