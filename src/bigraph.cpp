#include "bigraph.h"
#include "asqg.h"
#include "kseq.h"
#include "utils.h"

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Bigraph"));

#define GZIP_EXT ".gz"

//
// Edge
//

Edge::Dir Edge::EDGE_DIRECTIONS[2] = { Edge::ED_SENSE, Edge::ED_ANTISENSE };

std::string Edge::label() const {
    // get the unmatched coordinates in end vertex
    const SeqCoord& coord = _twin->coord();
    SeqCoord unmatched = coord.complement();
    const std::string& seq = _end->seq();
    std::string label = seq.substr(unmatched.interval.start, unmatched.length());
    if (comp() == EC_REVERSE) {
        make_reverse_complement_dna(label);
    }
    return label;
}

void Edge::join(Edge* edge) {
    // Update the match coordinate
    Match match = edge->match();
    _coord = match.translate10(_coord);

    if (edge->_comp == Edge::EC_REVERSE) {
        _comp = (Edge::Comp)(Edge::EC_COUNT - _comp - 1);
        _dir = (Edge::Dir)(Edge::ED_COUNT - _dir - 1);
    }

    _twin->extend(edge->_twin);
}

void Edge::extend(Edge* edge) {
    if (edge->_comp == Edge::EC_REVERSE) {
        _comp = (Edge::Comp)(Edge::EC_COUNT - _comp - 1);
    }
    _end = edge->end();
}

bool Edge::operator==(const Edge& edge) const {
    return _end->id() == edge._end->id() && _dir == edge._dir && _comp == edge._comp;
}

void Edge::validate() const {
    std::string v1 = start()->seq(), v2 = end()->seq();
    std::string m1 = v1.substr(_coord.interval.start, _coord.length());
    std::string m2 = v2.substr(_twin->_coord.interval.start, _twin->_coord.length());
    if (_comp == EC_REVERSE) {
        make_reverse_complement_dna(m2);
    }
    if (m1 != m2) {
        LOG4CXX_ERROR(logger, "Error, matching strings are not the same length");
        LOG4CXX_ERROR(logger, boost::format("V1M: %s") % m1);
        LOG4CXX_ERROR(logger, boost::format("V2M: %s") % m2);
        LOG4CXX_ERROR(logger, boost::format("V1: %s") % v1)
        assert(false);
    }
}

//
// Vertex
//
Vertex::~Vertex() {
    for (EdgePtrList::iterator i = _edges.begin(); i != _edges.end(); ++i) {
        delete *i;
    }
}

void Vertex::merge(Edge* edge) {
    // Merging two string vertices has two parts
    // First, the sequence of the vertex is extended
    // by the the content of the edge label
    // Then, all the edges that are pointing to this node
    // must be updated to contain the extension of the vertex
    Edge* twin = edge->twin();

    // Merge the sequence
    std::string label = edge->label();
    bool prepend = false;

    if (edge->dir() == Edge::ED_SENSE) {
        _seq += label;
    } else {
        _seq = label + _seq;
        prepend = true;
    }

    // Update the coverage value of the vertex
    _coverage += edge->end()->coverage();

    ////////////////////////////////////////////////////
    // Extend match
    // 
    // V1: |________________++++++|
    // V2:                  |++++++________________|
    //
    // V1->merge(edge)
    //
    // V1: |________________++++++++++++++++++++++|
    // V2:                 |++++++++++++++++++++++|
    //
    // V2->merge(edge)
    // prepend
    //
    // V1: |++++++++++++++++++++++|
    // V2: |++++++++++++++++++++++________________|
    //
    ////////////////////////////////////////////////////
    edge->coord().stretch(label.length());
    twin->coord().extend(label.length());

    // All the SeqCoords for the edges must have their seqlen field updated
    // Also, if we prepended sequence to this edge, all the matches in the 
    // SENSE direction must have their coordinates offset
    for (EdgePtrList::iterator i = _edges.begin(); i != _edges.end(); ++i) {
        SeqCoord& coord = (*i)->coord();
        coord.seqlen = _seq.length();
        if (prepend && (*i)->dir() == Edge::ED_SENSE && edge != *i) {
            coord.interval.offset(label.length());
        }
    }
}

void Vertex::addEdge(Edge* edge) {
    assert(edge->start() == this);
    _edges.push_back(edge);
}

void Vertex::removeEdge(Edge* edge) {
    for (EdgePtrList::iterator i = _edges.begin(); i != _edges.end(); ++i) {
        if (*i == edge) {
            _edges.erase(i);
            return;
        }
    }
    assert(false);
}

bool Vertex::hasEdge(Edge* edge) const {
    for (EdgePtrList::const_iterator i = _edges.begin(); i != _edges.end(); ++i) {
        if (**i == *edge) {
            return true;
        }
    }
    return false;
}

size_t Vertex::sweepEdges(GraphColor c) {
    size_t num = 0;

    EdgePtrList::iterator i = _edges.begin();
    while (i != _edges.end()) {
        Edge* edge = *i;
        if (edge->color() == c) {
            // Remove the edges pointing to this Vertex
            i = _edges.erase(i);
            SAFE_DELETE(edge);
            ++num;
        } else {
            ++i;
        }
    }

    return num;
}

void Vertex::deleteEdges() {
    for (EdgePtrList::iterator i = _edges.begin(); i != _edges.end(); ++i) {
        Edge* edge = *i;
        Edge* twin = edge->twin();

        Vertex* partner = edge->end();
        partner->removeEdge(twin);

        SAFE_DELETE(twin);
        SAFE_DELETE(edge);
    }
    _edges.clear();
}

void Vertex::validate() const {
    for (EdgePtrList::const_iterator i = _edges.begin(); i != _edges.end(); ++i) {
        (*i)->validate();
    }
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
    VertexTable::iterator i = _vertices.find(vertex->id());
    if (i != _vertices.end()) {
        return false;
    }
    _vertices[vertex->id()] = vertex;
    return true;
}

Vertex* Bigraph::getVertex(const Vertex::Id& id) const {
    VertexTable::const_iterator i = _vertices.find(id);
    if (i != _vertices.end()) {
        return i->second;
    }
    return NULL;
}

void Bigraph::removeVertex(Vertex* vertex) {
    const Vertex::Id& vid = vertex->id();
    _vertices.erase(vid);
}

size_t Bigraph::sweepVertices(GraphColor c) {
    size_t num = 0;

    VertexTable::iterator i = _vertices.begin();
    while (i != _vertices.end()) {
        VertexTable::iterator next = i;
        ++next;
        if (i->second->color() == c) {
            // Remove the edges pointing to this Vertex
            i->second->deleteEdges();

            // Remove the vertex from the collection
            removeVertex(i->second);
            SAFE_DELETE(i->second);
            ++num;
        }
        i = next;
    }

    return num;
}

void Bigraph::addEdge(Vertex* vertex, Edge* edge) {
    vertex->addEdge(edge);
}

size_t Bigraph::sweepEdges(GraphColor c) {
    size_t num = 0;

    for (VertexTable::iterator i = _vertices.begin(); i != _vertices.end(); ++i) {
        num += i->second->sweepEdges(c);
    }

    return num;
}

void Bigraph::simplify() {
    simplify(Edge::ED_SENSE);
    simplify(Edge::ED_ANTISENSE);
}

void Bigraph::simplify(Edge::Dir dir) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (VertexTable::iterator i = _vertices.begin(); i != _vertices.end(); ++i) {
            // Get the edges for this direction
            EdgePtrList edges = i->second->edges(dir);

            // If there is a single edge in this direction, merge the vertices
            // Don't merge singular self edges though
            if (edges.size() == 1 && !edges[0]->isSelf()) {
                // Check that the edge back is singular as well
                Edge* single = edges[0];
                Edge* twin = single->twin();
                Vertex* end = single->end();
                if (end->degrees(twin->dir()) == 1) {
                    merge(i->second, single);
                    changed = true;
                }
            }
        }
    }
}

void Bigraph::merge(Vertex* v1, Edge* edge) {
    Vertex* v2 = edge->end();

    // Merge the data
    v1->merge(edge);

    // Get the twin edge (the edge in V2 that points to V1)
    Edge* twin = edge->twin();

    // Ensure V2 has the twin edge
    assert(v2->hasEdge(twin));

    // Get the edge set opposite of the twin edge (which will be the new edges in this direction for V1)
    EdgePtrList transEdges = v2->edges((Edge::Dir)(Edge::ED_COUNT - twin->dir() - 1));

    // Move the edges from V2 to V1
    for (EdgePtrList::iterator i = transEdges.begin(); i != transEdges.end(); ++i) {
        Edge* transEdge = *i;

        // Remove the edge from V2, this does not destroy the edge
        v2->removeEdge(transEdge);

        // Join edge to the start of transEdge
        // This updates the starting point of transEdge to be V1
        // This calls Edge::extend on the twin edge
        transEdge->join(edge);
        assert(transEdge->dir() == edge->dir());
        v1->addEdge(transEdge);
    }

    // Remove the edge from V1 to V2
    v1->removeEdge(edge);
    SAFE_DELETE(edge);

    // Remove the edge from V2 to V1
    v2->removeEdge(twin);
    SAFE_DELETE(twin);

    // Remove V2
    // It is guarenteed to not be connected
    removeVertex(v2);
    SAFE_DELETE(v2);
}

void Bigraph::validate() const {
    for (VertexTable::const_iterator i = _vertices.begin(); i != _vertices.end(); ++i) {
        i->second->validate();
    }
}

bool Bigraph::visit(BigraphVisitor* visitor) {
    bool modified = false;

    visitor->previsit(this);
    for (VertexTable::const_iterator i = _vertices.begin(); i != _vertices.end(); ++i) {
        modified |= visitor->visit(this, i->second);
    }
    visitor->postvisit(this);

    return modified;
}

void Bigraph::color(GraphColor c) {
    for (VertexTable::const_iterator i = _vertices.begin(); i != _vertices.end(); ++i) {
        i->second->color(c);
    }
}

class EdgeCreator {
public:
    EdgeCreator(Bigraph* g, bool allowContainments, size_t maxEdges) : _graph(g), _allowContainments(allowContainments), _maxEdges(maxEdges) {
    }

    bool create(const Overlap& overlap) {
        // Initialize data and perform checks
        bool isContainment = overlap.match.isContainment();
        Edge::Comp comp = (overlap.match.isRC) ? Edge::EC_REVERSE : Edge::EC_SAME;

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
            size_t degrees0 = verts[0]->degrees(), degrees1 = verts[1]->degrees();
            if (degrees0 > _maxEdges || degrees1 > _maxEdges) {
                return false;
            }
        }

        if (!isContainment) {
            Edge* edges[2];
            for (size_t i = 0; i < 2; ++i) {
                const SeqCoord& coord = overlap.match.coords[i];
                Edge::Dir dir = coord.isLeftExtreme() ? Edge::ED_ANTISENSE : Edge::ED_SENSE;
                edges[i] = new Edge(verts[1 - i], dir, comp, coord);
            }

            edges[0]->twin(edges[1]);
            edges[1]->twin(edges[0]);

            _graph->addEdge(verts[0], edges[0]);
            _graph->addEdge(verts[1], edges[1]);
        } else {
            // Contained edges don't have a direction, they can be travelled from
            // one vertex to the other in either direction. Hence, we 
            // add two edges per vertex. Later during the contain removal
            // algorithm this is important to determine transitivity
            Edge* edges[4];
            for (size_t i = 0; i < 2; ++i) {
                const SeqCoord& coord = overlap.match.coords[i];
                edges[i    ] = new Edge(verts[1 - i], Edge::ED_SENSE,     comp, coord);
                edges[i + 2] = new Edge(verts[1 - i], Edge::ED_ANTISENSE, comp, coord);
            }

            // Twin the edges and add them to the graph
            edges[0]->twin(edges[1]);
            edges[1]->twin(edges[0]);

            edges[2]->twin(edges[3]);
            edges[3]->twin(edges[2]);

            _graph->addEdge(verts[0], edges[0]);
            _graph->addEdge(verts[1], edges[1]);

            _graph->addEdge(verts[0], edges[2]);
            _graph->addEdge(verts[1], edges[3]);
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
                const ASQG::IntTagValue& containment = record.containment();
                if (containment) {
                    g->containment((int)containment);
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
                Vertex* vertex = new Vertex(record.id, record.seq, record.substring ? (int)record.substring : false);
                if (!g->addVertex(vertex)) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Attempted to insert vertex into graph with a duplicate id: %s") % vertex->id());
                    LOG4CXX_ERROR(logger, "All reads must have a unique identifier");
                    SAFE_DELETE(vertex);
                    return false;
                }
                if (vertex->contained()) {
                    g->containment(true);
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
                    creator.create(ovr);
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
    try {
        boost::iostreams::filtering_istreambuf buf;
        if (boost::algorithm::ends_with(filename, GZIP_EXT)) {
            buf.push(boost::iostreams::gzip_decompressor());
        }
        buf.push(boost::iostreams::file_descriptor_source(filename));

        std::istream stream(&buf);
        return loadASQG(stream, minOverlap, allowContainments, maxEdges, g);
    } catch (...) {
        return false;
    }
}

bool saveASQG(std::ostream& stream, const Bigraph* g) {
    // Header
    {
        ASQG::HeaderRecord record;
        stream << record << '\n';
        if (!stream) {
            return false;
        }
    }
    // Vertices
    for (VertexTable::const_iterator i = g->_vertices.begin(); i != g->_vertices.end(); ++i) {
        ASQG::VertexRecord record(i->second->id(), i->second->seq());
        stream << record << '\n';
        if (!stream) {
            return false;
        }
    }
    // Edges
    for (VertexTable::const_iterator i = g->_vertices.begin(); i != g->_vertices.end(); ++i) {
        EdgePtrList edges = i->second->edges();
        for (EdgePtrList::const_iterator j = edges.begin(); j != edges.end(); ++j) {
            // We write one record for every bidirectional edge so only write edges
            // that are in canonical form (where id1 < id2)
            Edge* edge = (*j);
            Edge* twin = edge->twin();

            Overlap overlap(edge->start()->id(), edge->coord(), edge->end()->id(), twin->coord(), edge->comp() == Edge::EC_REVERSE, 0);
            if (overlap.id[0] <= overlap.id[1]) {
                // Containment edges are in both directions so only output one
                // record if it is a containment
                if (!overlap.match.isContainment() || edge->dir() == Edge::ED_SENSE) {
                    ASQG::EdgeRecord record(overlap);
                    stream << record << '\n';
                    if (!stream) {
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

bool saveASQG(const std::string& filename, const Bigraph* g) {
    try {
        boost::iostreams::filtering_ostreambuf buf;
        if (boost::algorithm::ends_with(filename, GZIP_EXT)) {
            buf.push(boost::iostreams::gzip_compressor());
        }
        buf.push(boost::iostreams::file_descriptor_sink(filename));

        std::ostream stream(&buf);
        return saveASQG(stream, g);
    } catch (...) {
        return false;
    }
}
