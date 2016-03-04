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

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Overlap"));

class Overlap : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        return r;
    }

private:
    Overlap() : Runner("c:s:a:t:p:g:h", boost::assign::map_list_of('a', "algorithm")('t', "threads")('p', "prefix")('m', "min-overlap")) {
        RUNNER_INSTALL("overlap", this, "compute overlaps between reads");
    }
    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("h") != options.not_found() || arguments.size() != 1) {
            return printHelps();
        }
        return 0;
    }
    int printHelps() const {
        std::cout << boost::format(
                "%s overlap [OPTION] ... READSFILE\n"
                "Compute pairwise overlap between all the sequences in READS\n"
                "\n"
                "      -h, --help                       display this help and exit\n"
                "\n"
                "      -t, --threads=NUM                use NUM threads to construct the index (default: 1)\n"
                "      -m, --min-overlap=LEN            minimum overlap required between two reads (default: 45)\n"
                "      -p, --prefix=PREFIX              write index to file using PREFIX instead of prefix of READSFILE\n"
                "\n"
                ) % PACKAGE_NAME << std::endl;
        return 256;
    }

    static Overlap _runner;
};

Overlap Overlap::_runner;

