#include "config.h"
#include "constant.h"
#include "fmindex.h"
#include "overlap_builder.h"
#include "runner.h"

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.DuplicateRemove"));

class DuplicateRemove : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        std::string input = arguments[0];
        LOG4CXX_INFO(logger, boost::format("input: %s") % input);
        //input = boost::filesystem::path(input).stem().string();

        std::string output = boost::filesystem::path(input).stem().string();
        if (options.find("prefix") != options.not_found()) {
            output = options.get< std::string >("prefix");
        }
        LOG4CXX_INFO(logger, boost::format("output: %s") % output);

        FMIndex fmi, rfmi;
        if (FMIndex::load(output + BWT_EXT, fmi) && FMIndex::load(output + RBWT_EXT, rfmi)) {
            OverlapBuilder builder(&fmi, &rfmi, output);
            if (!builder.rmdup(input, output + RMDUP_EXT + ".fa", output + RMDUP_EXT + ".dups.fa")) {
                LOG4CXX_ERROR(logger, boost::format("Failed to remove duplicates from reads %s") % input);
                r = -1;
            }
        } else {
            LOG4CXX_ERROR(logger, boost::format("Failed to load FMIndex from %s") % input);
            r = -1;
        }

        return r;
    }

private:
    DuplicateRemove(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts) : Runner(shortopts, longopts) {
        RUNNER_INSTALL(name, this, description, kRmDup);
    }
    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("h") != options.not_found() || arguments.size() != 1) {
            return printHelps();
        }
        return 0;
    }
    int printHelps() const {
        std::cout << boost::format(
                "%s rmdup [OPTION] ... READSFILE\n"
                "Remove duplicated reads from the data set\n"
                "\n"
                "      -h, --help                       display this help and exit\n"
                "      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
                "      -t, --threads=N                  use N threads (default: 1)\n"
                "      -d, --sample-rate=N              sample the symbol counts every N symbols in the FM-index. Higher values use significantly\n"
                "                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
                "\n"
                ) % PACKAGE_NAME << std::endl;
        return 256;
    }

    static DuplicateRemove _runner;
};

static const std::string shortopts = "c:s:t:p:d:h";
enum { OPT_HELP = 1 };
static const option longopts[] = {
    {"prefix",              required_argument,  NULL, 'o'}, 
    {"threads",             required_argument,  NULL, 't'}, 
    {"sample-rate",         required_argument,  NULL, 'd'}, 
    {"help",                no_argument,        NULL, 'h'}, 
    {NULL, 0, NULL, 0}, 
};
DuplicateRemove DuplicateRemove::_runner(
        "rmdup", 
        "duplicate reads removal", 
        shortopts,  
        longopts 
        );

