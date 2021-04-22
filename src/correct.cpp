#include <config.h>

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

#include "constant.h"
#include "correct_processor.h"
#include "fmindex.h"
#include "runner.h"
#include "utils.h"

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
    // input = Utils::stem(input);
    std::string stem = Utils::stem(input);

    std::string outfile = stem + EC_EXT + FA_EXT;
    if (options.find("outfile") != options.not_found()) {
      outfile = options.get<std::string>("outfile");
    }

    std::string prefix = stem;
    if (options.find("prefix") != options.not_found()) {
      prefix = options.get<std::string>("prefix");
    }

    FMIndex fmi;
    if (FMIndex::load(prefix + BWT_EXT, fmi)) {
      // Prepare parameters
      CorrectProcessor::Options parms(options);

      CorrectProcessor processor(parms);
      if (!processor.process(fmi, input, outfile, options.get<size_t>("threads", kCorrectThreads))) {
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
  Correct(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts)
      : Runner(shortopts, longopts) {
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
        "      -o, --outfile=FILE               write the corrected reads to FILE (default READFILE%s%s)\n"
        "      -t, --threads=NUM                use NUM threads for the computation (default: %d)\n"
        "      -a, --algorithm=STR              specify the correction algorithm to use. STR must be one of kmer,overlap. (default: %s)\n"
        "\n"
        "      -k, --kmer-size=N                the length of the kmer to user (default: %d)\n"
        "      -x, --kmer-threshold=N           attempt to correct kmers that are seen less than N times (default: %d)\n"
        "      -i, --kmer-rounds=N              perform up to N rounds of kmer correction (default: %d)\n"
        "      -O, --kmer-count-offset=N        when correcting a kmer, require the count of the new kmer is at least +N higher than the count of the old kmer. (default: %d)\n"
        "\n"
        ) % PACKAGE_NAME % EC_EXT % FA_EXT % kCorrectThreads % kCorrectAlgorithm % kCorrectKmerSize % kCorrectKmerThreshold % kCorrectKmerRounds % kCorrectKmerCountOffset << std::endl;
    return 256;
  }

  static Correct _runner;
};

static const std::string shortopts = "c:s:p:o:t:a:k:x:i:O:h";
enum { OPT_HELP = 1};
static const option longopts[] = {
    {"log4cxx",             required_argument,  nullptr, 'c'},
    {"ini",                 required_argument,  nullptr, 's'},
    {"prefix",              required_argument,  nullptr, 'p'},
    {"outfile",             required_argument,  nullptr, 'o'},
    {"threads",             required_argument,  nullptr, 't'},
    {"algorithm",           required_argument,  nullptr, 'a'},
    {"kmer-size",           required_argument,  nullptr, 'k'},
    {"kmer-threshold",      required_argument,  nullptr, 'x'},
    {"kmer-rounds",         required_argument,  nullptr, 'i'},
    {"kmer-count-offset",   required_argument,  nullptr, 'O'},
    {"help",                no_argument,        nullptr, 'h'},
    {nullptr, 0, nullptr, 0},
  };
Correct Correct::_runner(
    "correct",
    "correct sequencing errors in a set of reads",
    shortopts,
    longopts);
