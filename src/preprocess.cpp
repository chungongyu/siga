#include "constant.h"
#include "runner.h"

#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Preprocess"));

class Preprocess : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        LOG4CXX_DEBUG(logger, "preprocess reads begin");

        std::vector< std::string > filelist;
        std::copy(arguments.begin(), arguments.end(), std::back_inserter(filelist));
        LOG4CXX_DEBUG(logger, boost::format("input: %s") % boost::algorithm::join(filelist, ":"));

        LOG4CXX_DEBUG(logger, "preprocess reads end");
        return r;
    }

private:
    Preprocess() : Runner("c:s:n:d:i:o:t:ESe:h", boost::assign::map_list_of('e', "READ_LENGTH_CUTOFF")('t', "THRESHOLD")) {
        RUNNER_INSTALL("preprocess", this, "filter and quality-trim reads");
    }

    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("h") != options.not_found()) {
            return printHelps();
        }
        return 0;
    }
    int printHelps() const {
        std::cout << "arcs preprocess -o [output] -e [read_cutoff] -n [buckets] [<inputs>]" << std::endl;
        return 256;
    }

    static Preprocess _runner;
};

Preprocess Preprocess::_runner;
