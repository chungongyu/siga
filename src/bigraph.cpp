#include "bigraph.h"
#include "asqg.h"

#include <fstream>

#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Bigraph"));

//
// Vertex
//
void Vertex::addEdge(Edge* edge) {
}

//
// Bigraph
//
Bigraph::~Bigraph() {
    for (VertexTable::iterator i = _vertices.begin(); i != _vertices.end(); ++i) {
        delete i->second;
    }
}

bool Bigraph::addVertex(Vertex* vertex) {
    VertexTable::iterator i = _vertices.find(vertex->id);
    if (i != _vertices.end()) {
        return false;
    }
    _vertices[vertex->id] = vertex;
    return true;
}

Vertex* Bigraph::getVertex(const Vertex::Id& id) const {
    VertexTable::const_iterator i = _vertices.find(id);
    if (i != _vertices.end()) {
        return i->second;
    }
    return NULL;
}

void Bigraph::addEdge(Vertex* vertex, Edge* edge) {
    vertex->addEdge(edge);
}

class EdgeCreator {
public:
    EdgeCreator(Bigraph* g, bool allowContainments, size_t maxEdges) : _graph(g), _allowContainments(allowContainments), _maxEdges(maxEdges) {
    }

    bool create(const Overlap& overlap) {
        // Initialize data and perform checks
        bool isContainment = overlap.match.isContainment();

        Vertex* verts[2];
        for (size_t i = 0; i < 2; ++i) {
            verts[i] = _graph->getVertex(overlap.id[i]);

            // If one of the vertices is not in the graph, skip this edge
            // This can occur if one of the verts is a strict substring of some other vertex so it will
            // never be added to the graph
            if (verts[i] == NULL) {
                return false;
            }
        }

        // Check if this is a substring containment, if so mark the contained read
        // but do not create edges
        for (size_t i = 0; i < 2; ++i) {
            if (!overlap.match.coords[i].isExtreme()) {
                return false;
            }
        }

        // If either vertex has the maximum number of edges,
        // do not add any more. This is to protect against ultra-dense
        // regions of the graph inflating memory usage. The nodes that reach
        // this limit, and nodes connected to them are marked as super repeats.
        // After loading the graph, all edges to super repeats are cut to prevent
        // misassembly.
        {
        }

        if (!isContainment) {
            Edge* edges[2];
            for (size_t i = 0; i < 2; ++i) {
            }

            _graph->addEdge(verts[0], edges[0]);
            _graph->addEdge(verts[1], edges[1]);
        } else {
        }

        return true;
    }
private:
    Bigraph* _graph;
    bool _allowContainments;
    size_t _maxEdges;
};

bool loadASQG(std::istream& stream, size_t minOverlap, bool allowContainments, size_t maxEdges, Bigraph* g) {
    enum {
        STAGE_HEAD, 
        STAGE_VERTEX, 
        STAGE_EDGE
    };
    
    int stage = STAGE_HEAD;
    std::string line;
    while (std::getline(stream, line)) {
        ASQG::RecordType rt = ASQG::recordtype(line);
        switch(rt) {
            case ASQG::RT_HEADER: {
                if (stage != STAGE_HEAD) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected header record found at line %s") % line);
                    return false;
                }
                ASQG::HeaderRecord record;
                if (!ASQG::HeaderRecord::parse(line, record)) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected header record found at line %s") % line);
                    return false;
                }
                break;
            }
            case ASQG::RT_VERTEX: {
                if (stage == STAGE_HEAD) {
                    stage = STAGE_VERTEX;
                }
                if (stage != STAGE_VERTEX) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected vertex record found at line %s") % line);
                    return false;
                }
                ASQG::VertexRecord record;
                if (!ASQG::VertexRecord::parse(line, record)) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected vertex record found at line %s") % line);
                    return false;
                }
                Vertex* vertex = new Vertex(record.id, record.seq);
                if (!g->addVertex(vertex)) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Attempted to insert vertex into graph with a duplicate id: %s") % vertex->id);
                    LOG4CXX_ERROR(logger, "All reads must have a unique identifier");
                    delete vertex;
                    return false;
                }
                break;
            }
            case ASQG::RT_EDGE: {
                if (stage == STAGE_VERTEX) {
                    stage = STAGE_EDGE;
                }
                if (stage != STAGE_EDGE) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected edge record found at line %s") % line);
                    return false;
                }
                ASQG::EdgeRecord record;
                if (!ASQG::EdgeRecord::parse(line, record)) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected edge record found at line %s") % line);
                    return false;
                }
                const Overlap& ovr = record.overlap();
                // Add the edge to the graph
                if (ovr.match.length() >= minOverlap) {
                    EdgeCreator creator(g, allowContainments, maxEdges);
                    if (!creator.create(ovr)) {
                        return false;
                    }
                }

                break;
            }
            default: {
                return false;
            }
        }
    }

    return true;
}

bool loadASQG(const std::string& filename, size_t minOverlap, bool allowContainments, size_t maxEdges, Bigraph* g) {
    std::ifstream stream(filename.c_str());
    return loadASQG(stream, minOverlap, allowContainments, maxEdges, g);
}

