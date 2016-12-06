#ifndef sequence_process_framework_h_
#define sequence_process_framework_h_

#include "kseq.h"

#include <omp.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

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
    
        size_t consumed() const {
            return _consumed;
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
                if (postproc != NULL) {
                    postproc->process(workItem, output);
                }
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

#if _OPENMP
    template< class Input, class Output, class Generator, class Processor, class PostProcessor >
    class ParallelWorker {
    public:
        size_t run(Generator& generator, std::vector< Processor* >* proclist, PostProcessor* postproc, size_t n = -1) {
            size_t threads = proclist->size();
            omp_set_num_threads(threads);

            bool done = false;
            std::vector< Input > inputs;
            while (!done) {
                Input workItem;
                if (generator.consumed() < n && generator.generate(workItem)) {
                    inputs.push_back(workItem);
                } else {
                    done = true;
                }

                // Once all buffers are full or the input is finished, dispatch the work to the threads
                if (inputs.size() == threads * 1000 || done) { 
                    std::vector< Output > outputs(inputs.size());

                    #pragma omp parallel for schedule(dynamic, 256)
                    for (size_t i = 0; i < inputs.size(); ++i) {
                        size_t tid = omp_get_thread_num();
                        outputs[i] = (*proclist)[tid]->process(inputs[i]);
                    }

                    for (size_t i = 0; i < inputs.size(); ++i) {
                        if (postproc != NULL) {
                            postproc->process(inputs[i], outputs[i]);
                        }
                    }
                    inputs.clear();
                }
            }
            
            LOG4CXX_INFO(logger, boost::format("processed %d sequences") % generator.consumed());

            return generator.consumed();
        }

        size_t run(DNASeqReader& reader, std::vector< Processor* >* proclist, PostProcessor* postproc, size_t n = -1) {
            WorkItemGenerator< Input > generator(reader);
            return run(generator, proclist, postproc, n);
        }

        size_t run(std::istream& stream, std::vector< Processor* >* proclist, PostProcessor* postproc, size_t n = -1) {
            std::shared_ptr< DNASeqReader > reader(DNASeqReaderFactory::create(stream));
            if (reader) {
                return run(*reader, proclist, postproc, n);
            }
            return 0;
        }
        size_t run(const std::string& filename, std::vector< Processor* >* proclist, PostProcessor* postproc, size_t n = -1) {
            std::ofstream stream(filename.c_str());
            return run(stream, proclist, postproc, n);
        }
    private:
        static log4cxx::LoggerPtr logger;
    };

    template< class Input, class Output, class Generator, class Processor, class PostProcessor >
    log4cxx::LoggerPtr ParallelWorker< Input, Output, Generator, Processor, PostProcessor >::logger(log4cxx::Logger::getLogger("arcs.SequenceProcessFramework"));
#endif // _OPENMP
};

#endif // sequence_process_framework_h_
