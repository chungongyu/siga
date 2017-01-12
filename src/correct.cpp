#include "config.h"
#include "constant.h"
#include "correct_processor.h"
#include "fmindex.h"
#include "runner.h"

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Correct"));

class Correct : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        std::string input = arguments[0];
        LOG4CXX_INFO(logger, boost::format("input: %s") % input);
        //input = boost::filesystem::path(input).stem().string();

        std::string prefix = boost::filesystem::path(input).stem().string();
        if (options.find("prefix") != options.not_found()) {
            prefix = options.get< std::string >("prefix");
        }

        FMIndex fmi;
        if (FMIndex::load(prefix + BWT_EXT, fmi)) {
            // Prepare parameters
            CorrectProcessor::Options parms;

            CorrectProcessor processor(parms);
            if (!processor.process(input, prefix, options.get< size_t >("threads", 1))) {
                LOG4CXX_ERROR(logger, boost::format("Failed to do error correction for reads %s") % input);
                r = -1;
            }
        } else {
            LOG4CXX_ERROR(logger, boost::format("Failed to load FMIndex from %s") % prefix);
            r = -1;
        }

        return r;
    }

private:
    Correct(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts) : Runner(shortopts, longopts) {
        RUNNER_INSTALL(name, this, description, kCorrect);
    }
    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("help") != options.not_found() || arguments.size() != 1) {
            return printHelps();
        }
        return 0;
    }
    int printHelps() const {
        std::cout << boost::format(
                "%s correct [OPTION] ... READSFILE\n"
                "Correct sequencing errors in all the reads in READSFILE\n"
                "\n"
                "      -h, --help                       display this help and exit\n"
                "\n"
                "      -p, --prefix=PREFIX              use PREFIX instead of prefix of READSFILE for the names of the index files\n"
                "      -t, --threads=NUM                use NUM threads for the computation (default: 1)\n"
                "\n"
                "      -k, --kmer-size=N                the length of the kmer to user (default: 31)\n"
                "      -x, --kmer-threshold=N           attempt to correct kmers that are seen less than N times (default: 3)\n"
                "      -i, --kmer-rounds=N              perform up to N rounds of kmer correction (default: 10)\n"
                "\n"
                ) % PACKAGE_NAME << std::endl;
        return 256;
    }

    static Correct _runner;
};

static const std::string shortopts = "c:s:p:t:k:x:i:h";
enum { OPT_HELP = 1 };
static const option longopts[] = {
    {"prefix",              required_argument,  NULL, 'o'}, 
    {"prefip",              required_argument,  NULL, 'p'}, 
    {"threads",             required_argument,  NULL, 't'}, 
    {"kmer-size",           required_argument,  NULL, 'k'}, 
    {"kmer-threshold",      required_argument,  NULL, 'x'}, 
    {"kmer-rounds",         required_argument,  NULL, 'i'}, 
    {"help",                no_argument,        NULL, 'h'}, 
    {NULL, 0, NULL, 0}, 
};
Correct Correct::_runner(
        "correct", 
        "correct sequencing errors in a set of reads", 
        shortopts, 
        longopts
        );

