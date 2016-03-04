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

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Correct"));

class Correct : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        return r;
    }

private:
    Correct() : Runner("c:s:a:t:p:g:h", boost::assign::map_list_of('a', "algorithm")('t', "threads")('p', "prefix")('g', "gap-array")) {
        RUNNER_INSTALL("correct", this, "build the BWT and FM-index for a set of reads");
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

    static Correct _runner;
};

Correct Correct::_runner;

