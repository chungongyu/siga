#ifndef bigraph_h_
#define bigraph_h_

#include "coord.h"

#include <cassert>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

enum GraphColor {
    GC_WHITE = 0, 
    GC_GRAY, 
    GC_BLACK, 
    GC_BLUE, 
    GC_RED
};

class Vertex;

//
// Edge
//
class Edge {
public:
    // The directions an edge can take.
    // Do not change the values
    enum Dir {
        ED_SENSE = 0,
        ED_ANTISENSE = 1,
        ED_COUNT = 2
    };

    // Flag indicating whether the sequences linked by an edge
    // are from the same strand or not.
    // Do not change the values
    enum Comp {
        EC_SAME = 0,
        EC_REVERSE = 1
    };

    Edge(Vertex* end, Dir dir, Comp comp, const SeqCoord& coord) : _end(end), _dir(dir), _comp(comp), _coord(coord) {
    }

    Vertex* start() const {
        assert(_twin != NULL);
        return _twin->end();
    }
    Vertex* end() const {
        return _end;
    }
    Dir dir() const {
        return _dir;
    }
    Comp comp() const {
        return _comp;
    }
    void twin(Edge* twin) {
        _twin = twin;
    }
    Edge* twin() const {
        return _twin;
    }
    void color(GraphColor color) {
        _color = color;
    }
    GraphColor color() const {
        return _color;
    }
private:
    Edge* _twin;
    Vertex* _end;
    Dir _dir;
    Comp _comp;
    SeqCoord _coord;
    GraphColor _color;
};

typedef std::vector< Edge* > EdgePtrList;

//
// Vertex
//
class Vertex {
public:
    typedef std::string Id;

    Vertex(const Id& id, const std::string& seq) : id(id), seq(seq), color(GC_WHITE), coverage(1) {
    }
    ~Vertex();

    // Merge another vertex into this vertex, as specified by pEdge
    void merge(Edge* edge);

    // Edge list operations
    void addEdge(Edge* edge);
    void removeEdge(Edge* edge);

    EdgePtrList edges() const {
        return _edges;
    }
    EdgePtrList edges(Edge::Dir dir) const {
        EdgePtrList ev;
        for (EdgePtrList::const_iterator i = _edges.begin(); i != _edges.end(); ++i) {
            if ((*i)->dir() == dir) {
                ev.push_back(*i);
            }
        }
        return ev;
    }

    size_t degrees() const {
        return _edges.size();
    }
    size_t degrees(Edge::Dir dir) const {
        EdgePtrList ev = edges(dir);
        return ev.size();
    }

    Id id;
    GraphColor color;
    std::string seq;
    size_t coverage; // Number of vertices that have been merged into this one

private:
    EdgePtrList _edges;
};

typedef std::unordered_map< Vertex::Id, Vertex* > VertexTable;

//
// Bigraph
//
class Bigraph {
public:
    Bigraph() {
    }
    ~Bigraph();

    bool addVertex(Vertex* vertex);
    Vertex* getVertex(const Vertex::Id& id) const;

    void addEdge(Vertex* vertex, Edge* edge);
private:
    VertexTable _vertices;
};

bool loadASQG(const std::istream& stream, size_t minOverlap, bool allowContainments, size_t maxEdges, Bigraph* g);
bool loadASQG(const std::string& filename, size_t minOverlap, bool allowContainments, size_t maxEdges, Bigraph* g);

#endif // bigraph_h_
