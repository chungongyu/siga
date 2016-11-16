#include "constant.h"
#include "runner.h"

#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.GapFiller"));

class GapFiller : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        LOG4CXX_DEBUG(logger, "gapfill reads begin");

        std::vector< std::string > filelist;
        std::copy(arguments.begin(), arguments.end(), std::back_inserter(filelist));
        LOG4CXX_DEBUG(logger, boost::format("input: %s") % boost::algorithm::join(filelist, ":"));

        LOG4CXX_DEBUG(logger, "gapfill reads end");
        return r;
    }

private:
    GapFiller(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts) : Runner(shortopts, longopts) {
        RUNNER_INSTALL(name, this, description, kGapFill);
    }
    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("help") != options.not_found()) {
            return printHelps();
        }
        return 0;
    }
    int printHelps() const {
        std::cout << "arcs gapfill -o [output] -e [read_cutoff] -n [buckets] [<inputs>]" << std::endl;
        return 256;
    }

    static GapFiller _runner;
};

static const std::string shortopts = "c:s:n:d:i:o:t:ESe:h";
enum { OPT_HELP = 1 };
static const option longopts[] = {
    {"help",                no_argument,        NULL, 'h'}, 
    {NULL, 0, NULL, 0}, 
};
GapFiller GapFiller::_runner(
        "gapfill", 
        "filter and quality-trim reads", 
        shortopts,  
        longopts 
        );
