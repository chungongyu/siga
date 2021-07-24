#include <config.h>

#include <fstream>
#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

#include "constant.h"
#include "fmindex.h"
#include "kseq.h"
#include "runner.h"
#include "utils.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Match"));

class Match : public Runner {
 public:
  int run(const Properties& options, const Arguments& arguments) {
    int r = 0;

    if ((r = checkOptions(options, arguments)) != 0) {
      return r;
    }

    LOG4CXX_INFO(logger, boost::format("input: %s") % boost::algorithm::join(arguments, ","));
    // input = Utils::stem(input);

    std::string prefix = Utils::stem(arguments[0]);
    if (options.find("prefix") != options.not_found()) {
      prefix = options.get<std::string>("prefix");
    }

    FMIndex fmi;
    if (FMIndex::load(prefix + BWT_EXT, fmi)) {
      for (const auto& input : arguments) {
        std::shared_ptr<std::istream> stream(Utils::ifstream(input));
        if (!stream) {
          LOG4CXX_ERROR(logger, boost::format("Failed to read file %s") % input);
          r = -1;
          break;
        }
        std::shared_ptr<DNASeqReader> reader(DNASeqReaderFactory::create(*stream));
        if (!reader) {
          LOG4CXX_ERROR(logger, boost::format("Failed to create DNASeqReader %s") % input);
          r = -1;
          break;
        }
        DNASeq read;
        while (reader->read(read)) {
          if (read.seq.length() > options.get<size_t>("max-length", -1)) {
            std::string start = read.seq.substr(0, options.get<size_t>("max-length", -1));
            std::string end = read.seq.substr(read.seq.length() - options.get<size_t>("max-length", -1));
            std::cout << boost::format("VT\t0\t%s\t%s\t%d\n") % read.name % read.seq % (FMIndex::Interval::occurrences(start, &fmi) + (options.find("no-opposite-strand") == options.not_found() ? FMIndex::Interval::occurrences(make_dna_reverse_complement_copy(start), &fmi) : 0));
            std::cout << boost::format("VT\t1\t%s\t%s\t%d\n") % read.name % read.seq % (FMIndex::Interval::occurrences(end, &fmi) + (options.find("no-opposite-strand") == options.not_found() ? FMIndex::Interval::occurrences(make_dna_reverse_complement_copy(end), &fmi) : 0));
          } else {
            std::cout << boost::format("VT\t0\t%s\t%s\t%d\n") % read.name % read.seq % (FMIndex::Interval::occurrences(read.seq, &fmi) + (options.find("no-opposite-strand") == options.not_found() ? FMIndex::Interval::occurrences(make_dna_reverse_complement_copy(read.seq), &fmi) : 0));
          }
        }
      }
    } else {
      LOG4CXX_ERROR(logger, boost::format("Failed to load FMIndex from %s") % prefix);
      r = -1;
    }

    return r;
  }

 private:
  Match(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts) : Runner(shortopts, longopts) {
    RUNNER_INSTALL(name, this, description, kMatch);
  }
  int checkOptions(const Properties& options, const Arguments& arguments) const {
    if (options.find("help") != options.not_found() || arguments.empty()) {
      return printHelps();
    }
    return 0;
  }
  int printHelps() const {
    std::cout << boost::format(
        "%s match [OPTION] ... READSFILE\n"
        "Match reads in READSFILE with ref\n"
        "\n"
        "      -h, --help                       display this help and exit\n"
        "\n"
        "      -p, --prefix=PREFIX              use PREFIX instead of prefix of READSFILE for the names of the index files\n"
        "          --no-opposite-strand         treat all reads as forward strand\n"
        "\n"
        ) % PACKAGE_NAME << std::endl;
    return 256;
  }

  static Match _runner;
};

static const std::string shortopts = "c:s:p:t:l:h";
enum { OPT_HELP = 1, OPT_NO_RC };
static const option longopts[] = {
    {"log4cxx",             required_argument,  nullptr, 'c'},
    {"ini",                 required_argument,  nullptr, 's'},
    {"prefix",              required_argument,  nullptr, 'p'},
    {"threads",             required_argument,  nullptr, 't'},
    {"max-length",          required_argument,  nullptr, 'l'},
    {"no-opposite-strand",  no_argument,        nullptr, OPT_NO_RC},
    {"help",                no_argument,        nullptr, 'h'},
    {nullptr, 0, nullptr, 0},
  };
Match Match::_runner(
    "match",
    "match a set of reads with ref",
    shortopts,
    longopts);
