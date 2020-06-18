#include "asqg.h"
#include "config.h"
#include "constant.h"
#include "fmindex.h"
#include "overlap_builder.h"
#include "runner.h"

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Overlap"));

class Overlapping : public Runner {
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
            output = options.get<std::string>("prefix");
        }
        LOG4CXX_INFO(logger, boost::format("output: %s%s%s") % output % ASQG_EXT % GZIP_EXT);

        FMIndex fmi, rfmi;
        if (FMIndex::load(output + BWT_EXT, fmi) && FMIndex::load(output + RBWT_EXT, rfmi)) {
            OverlapBuilder builder(&fmi, &rfmi, output, options.find("exhaustive") == options.not_found(), options.find("no-opposite-strand") == options.not_found());
            if (!builder.build(input, options.get<size_t>("min-overlap", 10), output + ASQG_EXT + GZIP_EXT, options.get<size_t>("threads", 1), options.get<size_t>("batch-size", 1000))) {
                LOG4CXX_ERROR(logger, boost::format("Failed to build overlaps from reads %s") % input);
                r = -1;
            }
        } else {
            LOG4CXX_ERROR(logger, boost::format("Failed to load FMIndex from %s") % input);
            r = -1;
        }

        return r;
    }

private:
    Overlapping(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts) : Runner(shortopts, longopts) {
        RUNNER_INSTALL(name, this, description, kOverlap);
    }
    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("help") != options.not_found() || arguments.size() != 1) {
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
                "          --batch-size=NUM             use NUM batches for each thread (default: 1000)\n"
                "      -m, --min-overlap=LEN            minimum overlap required between two reads (default: 45)\n"
                "      -p, --prefix=PREFIX              write index to file using PREFIX instead of prefix of READSFILE\n"
                "      -x, --exhaustive                 output all overlaps, including transitive edges\n"
                "          --no-opposite-strand         treat all reads as forward strand\n"
                "\n"
                ) % PACKAGE_NAME << std::endl;
        return 256;
    }

    static Overlapping _runner;
};

static const std::string shortopts = "c:s:t:p:m:xh";
enum { OPT_HELP = 1, OPT_BATCH_SIZE, OPT_NO_RC };
static const option longopts[] = {
    {"log4cxx",             required_argument,  NULL, 'c'}, 
    {"ini",                 required_argument,  NULL, 's'}, 
    {"prefix",              required_argument,  NULL, 'p'}, 
    {"threads",             required_argument,  NULL, 't'}, 
    {"batch-size",          required_argument,  NULL, OPT_BATCH_SIZE}, 
    {"min-overlap",         required_argument,  NULL, 'm'}, 
    {"exhaustive",          no_argument,        NULL, 'x'}, 
    {"no-opposite-strand",  no_argument,        NULL, OPT_NO_RC}, 
    {"help",                no_argument,        NULL, 'h'}, 
    {NULL, 0, NULL, 0}, 
};
Overlapping Overlapping::_runner(
        "overlap", 
        "compute overlaps between reads", 
        shortopts, 
        longopts
        );

