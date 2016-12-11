#include "config.h"
#include "constant.h"
#include "runner.h"

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>
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
                "\n"
                ) % PACKAGE_NAME << std::endl;
        return 256;
    }

    static Correct _runner;
};

static const std::string shortopts = "c:s:a:t:p:g:h";
enum { OPT_HELP = 1 };
static const option longopts[] = {
    {"prefix",              required_argument,  NULL, 'o'}, 
    {"threads",             required_argument,  NULL, 't'}, 
    {"algorithm",           required_argument,  NULL, 'a'}, 
    {"help",                no_argument,        NULL, 'h'}, 
    {NULL, 0, NULL, 0}, 
};
Correct Correct::_runner(
        "correct", 
        "correct sequencing errors in a set of reads", 
        shortopts, 
        longopts
        );

