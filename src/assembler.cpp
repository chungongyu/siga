#include <config.h>

#include <iostream>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

#include "bigraph.h"
#include "bigraph_visitors.h"
#include "constant.h"
#include "runner.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Assembler"));

class Assembler : public Runner {
 public:
  int run(const Properties& options, const Arguments& arguments) {
    int r = 0;

    if ((r = checkOptions(options, arguments)) != 0) {
      return r;
    }

    std::string input = arguments[0];
    LOG4CXX_INFO(logger, boost::format("input: %s") % input);

    std::string output = options.get<std::string>("prefix", "default");
    LOG4CXX_INFO(logger, boost::format("output: %s") % output);

    size_t minOverlap = options.get<size_t>("min-overlap", 40);
    int peMode = options.get<int>("pe-mode", 0);

    LOG4CXX_INFO(logger, "parameters:");
    LOG4CXX_INFO(logger, boost::format("min-overlap: %d") % minOverlap);
    if (options.find("max-distance") != options.not_found()) {
      LOG4CXX_INFO(logger, boost::format("max-distance: %d") % options.get<size_t>("max-distance"));
    }
    if (options.find("insert-size") != options.not_found()) {
      LOG4CXX_INFO(logger, boost::format("insert-size: %d") % options.get<size_t>("insert-size"));
    }
    if (options.find("insert-size-delta") != options.not_found()) {
      LOG4CXX_INFO(logger, boost::format("insert-size-delta: %d") % options.get<size_t>("insert-size-delta"));
    }

    Bigraph g;
    if (Bigraph::load(input, minOverlap, true, options.get<size_t>("max-edges", -1), &g)) {
      g.validate();
      LOG4CXX_INFO(logger, "load ok");

      // Visitors
      ChimericVisitor chVisit(options.get<size_t>("min-chimeric-length", 0),
                              options.get<size_t>("max-chimeric-coverage", -1),
                              options.get<size_t>("max-chimeric-delta", -1),
                              options.get<size_t>("num-reads", 0),
                              options.get<size_t>("genome-size", 0),
                              options.get<double>("uniq-threshold", 0.0));
      ContainRemoveVisitor containVisit;
      LoopRemoveVisitor loopVisit;
      MaximumOverlapVisitor moVisit(options.get<size_t>("max-overlap-delta", 0),
                                    options.find("max-overlap-carefully") != options.not_found(),
                                    options.get<size_t>("num-reads", 0),
                                    options.get<size_t>("genome-size", 0),
                                    options.get<double>("uniq-threshold", 13.0));
      StatisticsVisitor statsVisit;
      TrimVisitor trimVisit(options.get<size_t>("min-branch-length", 150), options.get<size_t>("min-branch-coverage", -1));

      // Pre-assembly graph stats
      LOG4CXX_INFO(logger, "[Stats] Input graph:");
      g.visit(&statsVisit);

      if (peMode == 1) {
        size_t average = options.get("insert-size", 500);
        size_t delta = options.get("insert-size-delta", 100);
        if (options.find("insert-size") == options.not_found() || options.find("insert-size-delta") == options.not_found()) {
          size_t mu = 500, sigma = 100;
          InsertSizeEstimateVisitor iseVisit(mu, sigma);
          g.visit(&iseVisit);
          if (options.find("insert-size") == options.not_found()) {
            average = mu;
          }
          if (options.find("insert-size-delta") == options.not_found()) {
            delta = sigma;
          }
        }
        PairedReadVisitor prVisit(options.get<size_t>("max-distance", 100), average, delta, options.get<size_t>("max-search-nodes", 100), options.get<size_t>("threads", 1), options.get<size_t>("batch-size", 1000));
        g.visit(&prVisit);
      } else {
        // LOG4CXX_INFO(logger, "Removing contained vertices from graph");
        // while (g.containment()) {
        //     g.visit(&containVisit);
        // }
      }

      // Compact together unbranched chains of vertices
      g.simplify();

#ifdef HAVE_MLPACK
      // std::unique_ptr<AIModel<mlpack::tree::DecisionTree<> > > model;
      std::unique_ptr<BaggingModel<mlpack::tree::DecisionTree<> > > model;
      // ai
      if (options.find("ai-model") != options.not_found()) {
        model = std::unique_ptr<BaggingModel<mlpack::tree::DecisionTree<> > >(new BaggingModel<mlpack::tree::DecisionTree<> >());
        if (model->load(options.get<std::string>("ai-model"))) {
          LOG4CXX_INFO(logger, "load ai model ok");
          // AIVisitor aiVisit(&model, 
          //                 options.get<size_t>("num-reads", 0),
          //                 options.get<size_t>("genome-size", 0));
          // g.visit(&aiVisit);

          // // Compact together unbranched chains of vertices
          // g.simplify();
        } else {
          model.reset();
          LOG4CXX_INFO(logger, "load ai model failed");
        }
      }
      std::unique_ptr<AIVisitor> aiVisit;
      if (model) {
        aiVisit = std::unique_ptr<AIVisitor>(new AIVisitor(model.get(), 
                            options.get<size_t>("num-reads", 0),
                            options.get<size_t>("genome-size", 0)));
      }
#endif  // HAVE_MLPACK

      // 10x!
      if (options.find("with-index") != options.not_found()) {
        LinkedReadVisitor lrVisit;
        g.visit(&lrVisit);

        // Compact together unbranched chains of vertices
        g.simplify();
      }

      if (peMode == 1) {
        // Trimming
        size_t trimRound = 0, numTrimRounds = options.get<size_t>("cut-terminal", 10);
        while (trimRound < numTrimRounds) {
          LOG4CXX_INFO(logger, boost::format("[Trim] Trim round: %d") % (trimRound + 1));
          bool modified = false;

          LOG4CXX_INFO(logger, "Trimming tips");
          if (g.visit(&trimVisit)) {
            modified = true;
            g.simplify();
          }

          if (!modified) {
            break;
          }

          LOG4CXX_INFO(logger, "[Stats] After trimming:");
          g.visit(&statsVisit);

          ++trimRound;
        }
      } else {
        // Pre-assembly graph stats
        LOG4CXX_INFO(logger, "[Stats] After removing contained vertices:");
        g.visit(&statsVisit);

        LinkedReadVisitor lrVisit(options.get<size_t>("min-linkedread-length", -1), options.get<size_t>("min-linkedread-coverage", -1));
        // Trimming
        size_t trimRound = 0, numTrimRounds = options.get<size_t>("cut-terminal", 10);
        while (trimRound++ < numTrimRounds) {
          LOG4CXX_INFO(logger, boost::format("[Trim] Trim round: %d") % (trimRound));
          bool modified = false;

#ifdef HAVE_MLPACK
          if (aiVisit) {
            LOG4CXX_INFO(logger, "Removing edges by ai model");
            if (g.visit(aiVisit.get())) {
              modified = true;
              g.simplify();
            }
          }
#endif  // HAVE_MLPACK

          LOG4CXX_INFO(logger, "Removing loops");
          if (g.visit(&loopVisit)) {
            modified = true;
            g.simplify();
          }
          LOG4CXX_INFO(logger, "Removing non-maximum overlap edges from graph");
          if (g.visit(&moVisit)) {
            modified = true;
            g.simplify();
          }
          LOG4CXX_INFO(logger, "Trimming tips");
          if (g.visit(&trimVisit)) {
            modified = true;
            g.simplify();
          }

          // 10x!
          if (options.find("with-index") != options.not_found()) {
            LOG4CXX_INFO(logger, "Linked Reads");
            if (g.visit(&lrVisit)) {
              modified = true;
              g.simplify();
            }
          }

          if (options.get<size_t>("min-chimeric-length", 0) > 0) {
            LOG4CXX_INFO(logger, "Removing chimerics:");
            if (g.visit(&chVisit)) {
              modified = true;
              g.simplify();
            }
          }

          if (!modified) {
            break;
          }

          LOG4CXX_INFO(logger, "[Stats] After trimming:");
          g.visit(&statsVisit);
        }
      }

      LOG4CXX_INFO(logger, "[Stats] Final graph:");
      g.visit(&statsVisit);

      // Write the results
      {
        boost::filesystem::ofstream stream(output + "-contigs.fa");
        if (stream) {
          FastaVisitor faVist(stream);
          g.visit(&faVist);
        } else {
          LOG4CXX_ERROR(logger, boost::format("failed to create file %s-contigs.fa") % output);
          r = 1;
        }
      }
      if (!Bigraph::save(output + "-graph" + ASQG_EXT + GZIP_EXT, &g)) {
        LOG4CXX_ERROR(logger, boost::format("failed to create file %s-graph.asqg.gz") % output);
        r = 1;
      }
    } else {
      LOG4CXX_ERROR(logger, boost::format("failed to open file %s") % input);
      r = 1;
    }

    return r;
  }

private:
  Assembler(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts) : Runner(shortopts, longopts) {
    RUNNER_INSTALL(name, this, description, kAssemble);
  }

  int checkOptions(const Properties& options, const Arguments& arguments) const {
    if (options.find("help") != options.not_found() || arguments.size() != 1) {
      return printHelps();
    }
    return 0;
  }
  int printHelps() const {
    std::cout << boost::format(
        "%s assemble [OPTION] ... ASQGFILE\n"
        "Create contigs from the assembly graph ASQGFILE.\n"
        "\n"
        "      -h, --help                     display this help and exit\n"
        "\n"
        "      -p, --prefix=NAME              use NAME as the prefix of the output files (output files will be NAME-contigs.fa, etc)\n"
        "      -m, --min-overlap=LEN          only use overlaps of at least LEN. This can be used to filter\n"
        "                                     the overlap set so that the overlap step only needs to be run once\n"
        "          --max-edges=N              limit each vertex to a maximum of N edges. For highly repetitive regions\n"
        "                                     this helps save memory by culling excessive edges around unresolvable repeats (default: inf)\n"
        "      -t, --threads=NUM              use NUM threads to construct the paired graph (default: 1)\n"
        "          --batch-size=NUM           use NUM batches for each thread (default: 1000)\n"
        "\n"
        "Paired reads parameters:\n"
        "          --pe-mode=INT              0 - do not treat reads as paired (default)\n"
        "                                     1 - treat reads as paired\n"
        "          --with-index               treat as 10x linked read data\n"
        "          --max-distance=INT         treat reads as connected whose distance is less than INT (default: 100)\n"
        "          --insert-size=INT          treat reads as paired with insert size INT (default: learned from paired reads)\n"
        "          --insert-size-delta=INT    treat reads as paired with insert size delta INT (default: learned from paired reads)\n"
        "\n"
        "Trimming parameters:\n"
        "      -x, --cut-terminal=N           cut off terminal branches in N rounds (default: 10)\n"
        "      -n, --min-branch-length=LEN    remove terminal branches only if they are less than LEN bases in length (default: 150)\n"
        "\n"
        "Maximal overlap parameters:\n"
        "      -d, --max-overlap-delta=LEN    remove branches only if they are less than LEN bases in length (default: 0)\n"
        "          --max-overlap-carefully    remove non-maximum overlap edges carefully\n"
        "\n"
        ) % PACKAGE_NAME << std::endl;
    return 256;
  }

  static Assembler _runner;
};

static const std::string shortopts = "c:s:p:t:m:M:x:n:C:l:A:a:d:N:G:T:h";
enum { OPT_HELP = 1, OPT_BATCH_SIZE, OPT_PEMODE, OPT_WITH_IDX, OPT_MAXDIST, OPT_INSERTSIZE, OPT_INSERTSIZE_DELTA, OPT_MAXEDGES, OPT_MAX_OVERLAP_CAREFULLY };
static const option longopts[] = {
    {"log4cxx",              required_argument,  NULL, 'c'},
    {"ini",                  required_argument,  NULL, 's'},
    {"prefix",               required_argument,  NULL, 'p'},
    {"min-overlap",          required_argument,  NULL, 'm'},
    {"max-edges",            required_argument,  NULL, OPT_MAXEDGES},
    {"threads",              required_argument,  NULL, 't'},
    {"batch-size",           required_argument,  NULL, OPT_BATCH_SIZE},
    {"pe-mode",              required_argument,  NULL, OPT_PEMODE},
    {"with-index",           no_argument,        NULL, OPT_WITH_IDX},
    {"max-distance",         required_argument,  NULL, OPT_MAXDIST},
    {"insert-size",          required_argument,  NULL, OPT_INSERTSIZE},
    {"insert-size-delta",    required_argument,  NULL, OPT_INSERTSIZE_DELTA},
    {"cut-terminal",         required_argument,  NULL, 'x'},
    {"min-branch-length",    required_argument,  NULL, 'n'},
    {"min-branch-coverage",  required_argument,  NULL, 'C'},
    {"max-overlap-delta",    required_argument,  NULL, 'd'},
    {"max-overlap-carefully",required_argument,  NULL, OPT_MAX_OVERLAP_CAREFULLY},
    {"min-chimeric-length",  required_argument,  NULL, 'l'},
    {"min-chimeric-coverage",required_argument,  NULL, 'A'},
    {"max-chimeric-delta",   required_argument,  NULL, 'a'},
    {"num-reads",            required_argument,  NULL, 'N'},
    {"genome-size",          required_argument,  NULL, 'G'},
    {"uniq-threshold",       required_argument,  NULL, 'T'},
#ifdef HAVE_MLPACK
    {"ai-model",             required_argument,  NULL, 'M'},
#endif  // HAVE_MLPACK
    {"help",                 no_argument,        NULL, 'h'},
    {NULL, 0, NULL, 0},
};
Assembler Assembler::_runner(
    "assemble", 
    "generate contigs from an assembly graph", 
    shortopts,  
    longopts);
