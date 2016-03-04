#include "bwt.h"
#include "config.h"
#include "constant.h"
#include "kseq.h"
#include "runner.h"
#include "suffix_array.h"
#include "suffix_array_builder.h"

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Indexer"));

#define SAI_EXT  ".sai"
#define RSAI_EXT ".rsai"
#define BWT_EXT  ".bwt"
#define RBWT_EXT ".rbwt"

class Indexer : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        std::string input = arguments[0];
        LOG4CXX_INFO(logger, boost::format("input: %s") % input);

        std::string output = boost::filesystem::path(input).stem().string();
        if (options.find("prefix") != options.not_found()) {
            output = options.get< std::string >("prefix");
        }
        LOG4CXX_INFO(logger, boost::format("output: %s") % output);

        std::string algorithm = options.get< std::string >("algorithm", "sais");
        LOG4CXX_INFO(logger, boost::format("algorithm: %s") % algorithm);

        DNASeqList reads;
        if (ReadDNASequences(input, reads)) {
            std::shared_ptr< SuffixArrayBuilder > builder(SuffixArrayBuilder::create(algorithm));
            if (builder) {
                // forward
                {
                    std::shared_ptr< SuffixArray > sa(builder->build(reads));
                    if (sa) {
                        // suffix array
                        {
                            boost::filesystem::ofstream out(output + SAI_EXT);
                            out << *sa;
                        }
                        // bwt
                        {
                            boost::filesystem::ofstream out(output + BWT_EXT);
                            BWTWriter w(out);
                            w.write(*sa, reads);
                        }
                    }
                }

                BOOST_FOREACH(DNASeq& read, reads) {
                    read.make_reverse();
                }

                // reverse
                {
                    std::shared_ptr< SuffixArray > sa(builder->build(reads));
                    if (sa) {
                        // suffix array
                        {
                            boost::filesystem::ofstream out(output + RSAI_EXT);
                            out << *sa;
                        }
                        // bwt
                        {
                            boost::filesystem::ofstream out(output + RBWT_EXT);
                            BWTWriter w(out);
                            w.write(*sa, reads);
                        }
                    }
                }
            } else {
                LOG4CXX_ERROR(logger, boost::format("Failed to create suffix array builder algorithm %s") % algorithm);
            }
        } else {
            LOG4CXX_ERROR(logger, boost::format("Failed to open input stream %s") % input);
            r = -1;
        }

        return r;
    }

private:
    Indexer() : Runner("c:s:a:t:p:g:h", boost::assign::map_list_of('a', "algorithm")('t', "threads")('p', "prefix")('g', "gap-array")) {
        RUNNER_INSTALL("index", this, "build the BWT and FM-index for a set of reads");
    }
    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("h") != options.not_found() || arguments.size() != 1) {
            return printHelps();
        }
        return 0;
    }
    int printHelps() const {
        std::cout << boost::format(
                "%s index [OPTION] ... READSFILE\n"
                "Index the reads in READSFILE using a suffixarray/bwt\n"
                "\n"
                "      -h, --help                       display this help and exit\n"
                "\n"
                "      -a, --algorithm=STR              BWT construction algorithm. STR can be:\n"
                "                                       sais - induced sort algorithm, slower but works for very long sequences (default)\n"
                "                                       ropebwt - very fast and memory efficient. use this for short (<200bp) reads\n"
                "      -t, --threads=NUM                use NUM threads to construct the index (default: 1)\n"
                "      -c, --check                      validate that the suffix array/bwt is correct\n"
                "      -p, --prefix=PREFIX              write index to file using PREFIX instead of prefix of READSFILE\n"
                "          --no-reverse                 suppress construction of the reverse BWT. Use this option when building the index\n"
                "                                       for reads that will be error corrected using the k-mer corrector, which only needs the forward index\n"
                "          --no-forward                 suppress construction of the forward BWT. Use this option when building the forward and reverse index separately\n"
                "          --no-sai                     suppress construction of the SAI file. This option only applies to -a ropebwt\n"
                "      -g, --gap-array=N                use N bits of storage for each element of the gap array. Acceptable values are 4,8,16 or 32. Lower\n"
                "                                       values can substantially reduce the amount of memory required at the cost of less predictable memory usage.\n"
                "                                       When this value is set to 32, the memory requirement is essentially deterministic and requires ~5N bytes where\n"
                "                                       N is the size of the FM-index of READS2.\n"
                "                                       The default value is 8.\n"
                "\n"
                ) % PACKAGE_NAME << std::endl;
        return 256;
    }

    static Indexer _runner;
};

Indexer Indexer::_runner;

