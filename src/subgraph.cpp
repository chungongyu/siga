#include "bigraph.h"
#include "config.h"
#include "constant.h"
#include "runner.h"

#include <iostream>
#include <memory>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Subgraph"));

class SubgraphExtractor {
public:
    SubgraphExtractor(const Bigraph* g) : _graph(g) {
    }

    void extract(const Vertex* root, size_t span, Bigraph* sub) {
        assert(root != NULL);

        // copy graph parameters from the main graph
        sub->containment(_graph->containment());

        // add root to the subgraph
        addVertex(root, sub);

        if (span > 0) {
            EdgeCreator creator(sub, true, -1);
            // These are the edges in the main graph
            EdgePtrList edges = root->edges();
            BOOST_FOREACH(Edge* edge, edges) {
                if (edge->color() != GC_BLACK) {
                    Vertex* child = edge->end();
                    addVertex(child, sub);

                    Overlap overlap(root->id(), child->id(), edge->match());
                    creator.create(overlap);

                    edge->color(GC_BLACK);
                    edge->twin()->color(GC_BLACK);

                    extract(child, span - 1, sub);
                }
            }
        }
    }
private:
    void addVertex(const Vertex* vertex, Bigraph* sub) {
        // Make sure the vertex hasn't been added yet
        if (sub->getVertex(vertex->id()) == NULL) {
            sub->addVertex(new Vertex(vertex->id(), vertex->seq(), vertex->contained()));
        }
    }

    const Bigraph* _graph;
};

class Subgraph : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        const std::string rootId = arguments[0];
        const std::string input = arguments[1];
        const std::string output = options.get< std::string >("out", "subgraph.asqg.gz");
        LOG4CXX_INFO(logger, boost::format("input: %s") % input);

        Bigraph g;
        if (Bigraph::load(input, options.get< size_t >("min-overlap", 0), true, options.get< size_t >("max-edges", 128), &g)) {
            Bigraph sub;
            const Vertex* root = g.getVertex(rootId);
            if (root != NULL) {
                // add root to the subgraph
                SubgraphExtractor extractor(&g);
                extractor.extract(root, options.get< size_t >("size"), &sub);
                if (!Bigraph::save(output, &sub)) {
                    LOG4CXX_ERROR(logger, boost::format("failed to write stream ") % output);
                }
            } else {
                LOG4CXX_ERROR(logger, boost::format("Vertex %s not found in the graph.") % rootId);
            }
        } else {
            LOG4CXX_ERROR(logger, boost::format("failed to open stream %s") % input);
            r = -1;
        }

        return r;
    }

private:
    Subgraph(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts) : Runner(shortopts, longopts) {
        RUNNER_INSTALL(name, this, description, kSubgraph);
    }
    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("help") != options.not_found() || arguments.size() != 2) {
            return printHelps();
        }
        return 0;
    }
    int printHelps() const {
        std::cout << boost::format(
                "%s subgraph [OPTION] ... ID ASQGFILE\n"
                "Extract the subgraph around sequence with ID from an asqg file\n"
                "\n"
                "      -h, --help                       display this help and exit\n"
                "\n"
                "      -o, --out=FILE                   write the subgraph to FILE(default: subgraph.asqg.gz)\n"
                "      -m, --min-overlap=LEN            only use overlaps of at least LEN. This can be used to filter\n"
                "\n"
                "          --size=N                     the size of the subgraph to extract, all vertices that are at most N hops\n"
                "                                       away from the root will be included (default: 5)\n"
                "\n"
                ) % PACKAGE_NAME << std::endl;
        return 256;
    }

    static Subgraph _runner;
};

static const std::string shortopts = "o:m:h";
enum { OPT_HELP = 1, OPT_SIZE };
static const option longopts[] = {
    {"out",                 required_argument,  NULL, 'o'}, 
    {"min-overlap",         required_argument,  NULL, 'm'}, 
    {"size",                required_argument,  NULL, OPT_SIZE}, 
    {"help",                no_argument,        NULL, 'h'}, 
    {NULL, 0, NULL, 0}, 
};
Subgraph Subgraph::_runner(
        "subgraph", 
        "extract a subgraph from a graph", 
        shortopts, 
        longopts
        );

