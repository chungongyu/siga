#include "config.h"
#include "constant.h"
#include "runner.h"

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>

#include <divsufsort.h>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Indexer"));

class Indexer : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        LOG4CXX_DEBUG(logger, "indexer reads begin");

        std::vector< std::string > filelist;
        std::copy(arguments.begin(), arguments.end(), std::back_inserter(filelist));
        LOG4CXX_DEBUG(logger, boost::format("input: %s") % boost::algorithm::join(filelist, ":"));

        char* text = "abracadabra";
        size_t n = strlen(text);
        int* sa = new int[n];

        divsufsort((const sauchar_t *)text, sa, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = sa[i]; j < n; ++j) {
                std::cout << text[j];
            }
            std::cout << std::endl;
        }

        delete[] sa;

        LOG4CXX_DEBUG(logger, "indexer reads end");
        return r;
    }

private:
    Indexer() : Runner("c:s:n:d:i:o:t:ESe:h", boost::assign::map_list_of('e', "READ_LENGTH_CUTOFF")('t', "THRESHOLD")) {
        RUNNER_INSTALL("index", this, "build the BWT and FM-index for a set of reads");
    }
    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("h") != options.not_found()) {
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
                "      -d, --disk=NUM                   use disk-based BWT construction algorithm. The suffix array/BWT will be constructed\n"
                "                                       for batchs of NUM reads at a time. To construct the suffix array of 200 megabases of sequence\n"
                "                                       requires ~2GB of memory, set this parameter accordingly.\n"
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

