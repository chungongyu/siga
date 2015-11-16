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
        std::cout << "arcs index -o [output] -e [read_cutoff] -n [buckets] [<inputs>]" << std::endl;
        return 256;
    }

    static Indexer _runner;
};

Indexer Indexer::_runner;

