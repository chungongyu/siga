#include "bigraph.h"
#include "bigraph_visitors.h"
#include "constant.h"
#include "runner.h"

#include <iostream>
#include <memory>

#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

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

        std::string output = options.get< std::string >("prefix", "default");
        LOG4CXX_INFO(logger, boost::format("output: %s") % output);

        Bigraph g;
        if (loadASQG(input, options.get< size_t >("min-overlap", 0), false, options.get< size_t >("max-edges", 128), &g)) {
            g.validate();
            LOG4CXX_INFO(logger, "load ok");

            // Visitors
            ChimericVisitor chVisit(options.get< size_t >("min-chimeric-length", 0), options.get< size_t >("max-chimeric-delta", -1));
            ContainRemoveVisitor containVisit;
            MaximalOverlapVisitor moVisit(options.get< size_t >("max-overlap-delta", 0));
            SmoothingVisitor smoothVisit;
            StatisticsVisitor statsVisit;
            TrimVisitor trimVisit(options.get< size_t >("min-branch-length", 150));

            // Pre-assembly graph stats
            LOG4CXX_INFO(logger, "[Stats] Input graph:");
            g.visit(&statsVisit);

            LOG4CXX_INFO(logger, "Removing contained vertices from graph");
            while (g.containment()) {
                g.visit(&containVisit);
            }

            // Pre-assembly graph stats
            LOG4CXX_INFO(logger, "[Stats] After removing contained vertices:");
            g.visit(&statsVisit);

            // Compact together unbranched chains of vertices
            g.simplify();

            // Trimming
            size_t numTrimRounds = options.get< size_t >("cut-terminal", 10);
            for (size_t numTrim = 0; numTrim < numTrimRounds; ++numTrim) {
                LOG4CXX_INFO(logger, "Trimming tips");
                if (g.visit(&trimVisit)) {
                    g.simplify();
                }

                LOG4CXX_INFO(logger, "Performing variation smoothing");
                if (g.visit(&smoothVisit)) {
                    g.simplify();
                }

                if (options.get< size_t >("min-chimeric-length", 0) > 0) {
                    LOG4CXX_INFO(logger, "removing chimerics:");
                    if (g.visit(&chVisit)) {
                        g.simplify();
                    }
                }

                LOG4CXX_INFO(logger, "Removing non-maximal overlap edges from graph");
                if (g.visit(&moVisit)) {
                    g.simplify();
                }
            }

            if (numTrimRounds > 0) {
                LOG4CXX_INFO(logger, "[Stats] Graph after trimming:");
                g.visit(&statsVisit);
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
                    LOG4CXX_ERROR(logger, boost::format("failed to open stream %s-contigs.fa") % output);
                    r = 1;
                }
            }
            if (!saveASQG(output + "-graph.asqg.gz", &g)) {
                LOG4CXX_ERROR(logger, boost::format("failed to open stream %s-graph.asqg.gz") % output);
                r = 1;
            }
        } else {
            LOG4CXX_ERROR(logger, boost::format("failed to open stream %s") % input);
            r = 1;
        }

        return r;
    }

private:
    Assembler() : Runner("c:s:o:t:m:x:n:l:a:N:b:h", boost::assign::map_list_of('o', "prefix")('t', "threads")('m', "min-overlap")('x', "cut-terminal")('n', "min-branch-length")('l', "min-chimeric-length")('a', "max-chimeric-delta")('N', "max-edges")) {
        RUNNER_INSTALL("assemble", this, "generate contigs from an assembly graph");
    }

    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("h") != options.not_found() || arguments.size() != 1) {
            return printHelps();
        }
        return 0;
    }
    int printHelps() const {
        std::cout << boost::format(
                "%s assemble [OPTION] ... ASQGFILE\n"
                "Create contigs from the assembly graph ASQGFILE.\n"
                "\n"
                "      -h, --help                       display this help and exit\n"
                "\n"
                "      -o, --prefix=NAME            use NAME as the prefix of the output files (output files will be NAME-contigs.fa, etc)\n"
                "      -t, --threads=NUM                use NUM threads to construct the index (default: 1)\n"
                "      -m, --min-overlap=LEN            only use overlaps of at least LEN. This can be used to filter\n"
                "      -N, --max-edges=N                limit each vertex to a maximum of N edges. For highly repetitive regions\n"
                "                                       this helps save memory by culling excessive edges around unresolvable repeats (default: 128)\n"
                "\n"
                "Bubble/Variation removal parameters:\n"
                "      -b, --bubble=N                   perform N bubble removal steps (default: 3)\n"
                "\n"
                "Trimming parameters:\n"
                "      -x, --cut-terminal=N             cut off terminal branches in N rounds (default: 10)\n"
                "      -n, --min-branch-length=LEN      remove terminal branches only if they are less than LEN bases in length (default: 150)\n"
                "\n"
                "Chimeric parameters:\n"
                "      -l, --min-chimeric-length=LEN    remove chimerics only if they are less than LEN bases in length (default: 0)\n"
                "      -a, --max-chimeric-delta=LEN     remove chimerics only if they are less than LEN bases in length (default: 0)\n"
                "\n"
                ) % PACKAGE_NAME << std::endl;
        return 256;
    }

    static Assembler _runner;
};

Assembler Assembler::_runner;
