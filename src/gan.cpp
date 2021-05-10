#include "bigraph.h"
#include "bigraph_visitors.h"
#include "constant.h"
#include "fmindex.h"
#include "runner.h"

#include <iostream>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.GAN"));

class GANet : public Runner {
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

        size_t minOverlap = options.get<size_t>("min-overlap", 0);
        Bigraph g(options.get<size_t>("init-vertex-capacity", 0));
        if (Bigraph::load(input, minOverlap, true, options.get<size_t>("max-edges", 128), &g)) {
            g.validate();
            LOG4CXX_INFO(logger, "load ok");

            std::unique_ptr<FMIndex> ref;
            if (options.find("ref") != options.not_found()) {
                ref = std::unique_ptr<FMIndex>(new FMIndex());
                if (!FMIndex::load(options.get<std::string>("ref") + BWT_EXT, *ref)) {
                    ref = std::unique_ptr<FMIndex>();
                }
                if (ref) {
                    LOG4CXX_INFO(logger, "load ref ok");
                }
            }

            // Pre-assembly graph stats
            LOG4CXX_INFO(logger, "[Stats] Input graph:");
            StatisticsVisitor statsVisit;
            g.visit(&statsVisit);

            // GANVisitor
            GANVisitor ganVisit(std::cout, ref.get());

            // Trimming
            size_t trimRound = 0, numTrimRounds = options.get<size_t>("cut-terminal", 10);
            while (trimRound < numTrimRounds) {
                LOG4CXX_INFO(logger, boost::format("[Trim] Trim round: %d") % (trimRound + 1));
                bool modified = false;

                if (g.visit(&ganVisit)) {
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

            LOG4CXX_INFO(logger, "[Stats] Final graph:");
            g.visit(&statsVisit);

            // Write the results
            {
                boost::filesystem::ofstream stream(output + "-gan.fa");
                if (stream) {
                    FastaVisitor faVist(stream);
                    g.visit(&faVist);
                } else {
                    LOG4CXX_ERROR(logger, boost::format("failed to open stream %s-contigs.fa") % output);
                    r = 1;
                }
            }
            if (!Bigraph::save(output + "-gan" + ASQG_EXT + GZIP_EXT, &g)) {
                LOG4CXX_ERROR(logger, boost::format("failed to open stream %s-gan.asqg.gz") % output);
                r = 1;
            }
        } else {
            LOG4CXX_ERROR(logger, boost::format("failed to open stream %s") % input);
            r = 1;
        }

        return r;
    }

private:
    GANet(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts) : Runner(shortopts, longopts) {
        RUNNER_INSTALL(name, this, description, kUnknown);
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
                "      -h, --help                       display this help and exit\n"
                "\n"
                "      -p, --prefix=NAME                use NAME as the prefix of the output files (output files will be NAME-contigs.fa, etc)\n"
                "      -m, --min-overlap=LEN            only use overlaps of at least LEN. This can be used to filter\n"
                "                                       the overlap set so that the overlap step only needs to be run once\n"
                "          --max-edges=N                limit each vertex to a maximum of N edges. For highly repetitive regions\n"
                "                                       this helps save memory by culling excessive edges around unresolvable repeats (default: 128)\n"
                "          --init-vertex-capacity=INT   the initail capacity for veritices in bigraph INT (default 0)\n"
                "      -t, --threads=NUM                use NUM threads to construct the paired graph (default: 1)\n"
                "          --batch-size=NUM             use NUM batches for each thread (default: 1000)\n"
                "\n"
                "Paired reads parameters:\n"
                "          --pe-mode=INT                0 - do not treat reads as paired (default)\n"
                "                                       1 - treat reads as paired\n"
                "          --with-index                 treat as 10x linked read data\n"
                "          --max-distance=INT           treat reads as connected whose distance is less than INT (default: 100)\n"
                "          --insert-size=INT            treat reads as paired with insert size INT (default: learned from paired reads)\n"
                "          --insert-size-delta=INT      treat reads as paired with insert size delta INT (default: learned from paired reads)\n"
                "\n"
                "Trimming parameters:\n"
                "      -x, --cut-terminal=N             cut off terminal branches in N rounds (default: 10)\n"
                "      -n, --min-branch-length=LEN      remove terminal branches only if they are less than LEN bases in length (default: 150)\n"
                "\n"
                "Maximal overlap parameters:\n"
                "      -d, --max-overlap-delta=LEN      remove branches only if they are less than LEN bases in length (default: 0)\n"
                "\n"
                ) % PACKAGE_NAME << std::endl;
        return 256;
    }

    static GANet _runner;
};

static const std::string shortopts = "c:s:p:t:m:x:n:C:l:A:a:b:d:N:G:T:L:X:R:h";
enum { OPT_HELP = 1, OPT_BATCH_SIZE, OPT_PEMODE, OPT_WITH_IDX, OPT_MAXDIST, OPT_INSERTSIZE, OPT_INSERTSIZE_DELTA, OPT_MAXEDGES, OPT_INIT_VERTEX_CAPACITY };
static const option longopts[] = {
    {"log4cxx",             required_argument,  NULL, 'c'}, 
    {"ini",                 required_argument,  NULL, 's'}, 
    {"prefix",              required_argument,  NULL, 'p'}, 
    {"ref",                 required_argument,  NULL, 'R'}, 
    {"min-overlap",         required_argument,  NULL, 'm'}, 
    {"max-edges",           required_argument,  NULL, OPT_MAXEDGES}, 
    {"init-vertex-capacity",required_argument,  NULL, OPT_INIT_VERTEX_CAPACITY}, 
    {"threads",             required_argument,  NULL, 't'}, 
    {"batch-size",          required_argument,  NULL, OPT_BATCH_SIZE}, 
    {"pe-mode",             required_argument,  NULL, OPT_PEMODE}, 
    {"with-index",          no_argument,        NULL, OPT_WITH_IDX}, 
    {"max-distance",        required_argument,  NULL, OPT_MAXDIST}, 
    {"insert-size",         required_argument,  NULL, OPT_INSERTSIZE}, 
    {"insert-size-delta",   required_argument,  NULL, OPT_INSERTSIZE_DELTA}, 
    {"bubble",              required_argument,  NULL, 'b'}, 
    {"min-branch-length",   required_argument,  NULL, 'n'}, 
    {"min-branch-coverage", required_argument,  NULL, 'C'}, 
    {"max-overlap-delta",   required_argument,  NULL, 'd'}, 
    {"cut-terminal",        required_argument,  NULL, 'x'}, 
    {"min-chimeric-length", required_argument,  NULL, 'l'}, 
    {"min-chimeric-coverage",required_argument, NULL, 'A'}, 
    {"max-chimeric-delta",  required_argument,  NULL, 'a'}, 
    {"num-reads",           required_argument,  NULL, 'N'}, 
    {"genome-size",         required_argument,  NULL, 'G'}, 
    {"uniq-threshold",      required_argument,  NULL, 'T'}, 
    {"min-linkedread-length",required_argument, NULL, 'L'}, 
    {"min-linkedread-coverage",required_argument, NULL, 'X'}, 
    {"help",                no_argument,        NULL, 'h'}, 
    {NULL, 0, NULL, 0}, 
};
GANet GANet::_runner(
        "gan", 
        "generate an assembly network", 
        shortopts,  
        longopts 
        );
