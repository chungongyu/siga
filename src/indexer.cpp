#include <config.h>

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

#include "bwt.h"
#include "constant.h"
#include "kseq.h"
#include "runner.h"
#include "suffix_array.h"
#include "suffix_array_builder.h"
#include "utils.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Indexer"));

class Indexer : public Runner {
 public:
  int run(const Properties& options, const Arguments& arguments) {
    int r = 0;

    if ((r = checkOptions(options, arguments)) != 0) {
      return r;
    }

    std::string input = arguments[0];
    LOG4CXX_INFO(logger, boost::format("input: %s") % input);

    std::string output = Utils::stem(input);
    if (options.find("prefix") != options.not_found()) {
      output = options.get<std::string>("prefix");
    }
    LOG4CXX_INFO(logger, boost::format("output: %s.(%s|%s|%s|%s)") % output % SAI_EXT % BWT_EXT % RSAI_EXT % RBWT_EXT);

    std::string algorithm = options.get<std::string>("algorithm", "sais2");
    LOG4CXX_INFO(logger, boost::format("algorithm: %s") % algorithm);

    LOG4CXX_INFO(logger, boost::format("loading reads: %s") % input);
    DNASeqList reads;
    if (ReadDNASequences(input, reads, (~kSeqWithQuality)&(~kSeqWithComment))) {  // Strip Qualit & Comment
      std::shared_ptr<SuffixArrayBuilder> builder(SuffixArrayBuilder::create(algorithm));
      if (builder) {
        size_t threads = options.get<size_t>("threads", 1);

        // forward
        if (options.find("no-forward") == options.not_found()) {
          LOG4CXX_INFO(logger, boost::format("building forward index using %d threads") % threads);
          build(builder.get(), reads, threads, output + SAI_EXT, output + BWT_EXT);
        }

        // reverse
        if (options.find("no-reverse") == options.not_found()) {
          LOG4CXX_INFO(logger, boost::format("building reverse index using %d threads") % threads);
          for (auto& read : reads) {
            read.make_reverse();
          }

          build(builder.get(), reads, threads, output + RSAI_EXT, output + RBWT_EXT);
        }
        LOG4CXX_INFO(logger, boost::format("cpu: %.3f sec, max rss: %.3f GB") % Utils::cputime() % Utils::maxrss());
      } else {
        LOG4CXX_ERROR(logger, boost::format("Failed to create suffix array builder algorithm %s") % algorithm);
        r = -1;
      }
    } else {
      LOG4CXX_ERROR(logger, boost::format("Failed to open input file %s") % input);
      r = -1;
    }

    return r;
  }

 private:
  bool build(SuffixArrayBuilder* builder, const DNASeqList& reads, size_t threads, const std::string& safile, const std::string& bwtfile) {
    std::shared_ptr<SuffixArray> sa(builder->build(reads, threads));
    if (!sa) {
      return false;
    }

    // suffix array
    {
      boost::filesystem::ofstream out(safile);
      out << *sa;
      if (!out) {
        return false;
      }
    }
    // bwt
    {
      BWT bwt(*sa, reads);
      boost::filesystem::ofstream out(bwtfile);
      out << bwt;
      if (!out) {
        return false;
      }
    }
    return true;
  }

  Indexer(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts)
      : Runner(shortopts, longopts) {
    RUNNER_INSTALL(name, this, description, kIndex);
  }
  int checkOptions(const Properties& options, const Arguments& arguments) const {
    if (options.find("help") != options.not_found() || arguments.size() != 1) {
      return printHelps();
    }
    return 0;
  }
  int printHelps() const {
    std::cout << boost::format(
        "%s index [OPTION] ... READSFILE\n"
        "Index the reads in READSFILE using a suffixarray/bwt\n"
        "\n"
        "      -h, --help                       display this help and exit\n"
        "\n"
        "      -a, --algorithm=STR              BWT construction algorithm. STR can be:\n"
        "                                       sais - induced sort algorithm, slower but works for very long sequences\n"
        "                                       sais2 - very fast and works for very long sequences (default)\n"
        "                                       ropebwt - very fast and memory efficient. use this for short (<200bp) reads\n"
        "      -t, --threads=NUM                use NUM threads to construct the index (default: 1)\n"
        "      -p, --prefix=PREFIX              write index to file using PREFIX instead of prefix of READSFILE\n"
        "          --no-forward                 suppress construction of the forward BWT.\n"
        "          --no-reverse                 suppress construction of the reverse BWT.\n"
        "\n"
        ) % PACKAGE_NAME << std::endl;
    return 256;
  }

  static Indexer _runner;
};

static const std::string shortopts = "c:s:a:t:p:g:h";
enum { OPT_HELP = 1, OPT_NO_REVERSE, OPT_NO_FORWARD };
static const option longopts[] = {
    {"log4cxx",             required_argument,  nullptr, 'c'},
    {"ini",                 required_argument,  nullptr, 's'},
    {"prefix",              required_argument,  nullptr, 'p'},
    {"threads",             required_argument,  nullptr, 't'},
    {"algorithm",           required_argument,  nullptr, 'a'},
    {"no-reverse",          no_argument,        nullptr, OPT_NO_REVERSE},
    {"no-forward",          no_argument,        nullptr, OPT_NO_FORWARD},
    {"help",                no_argument,        nullptr, 'h'},
    {nullptr, 0, nullptr, 0},
  };
Indexer Indexer::_runner(
    "index",
    "build the BWT and FM-index for a set of reads",
    shortopts,
    longopts);
