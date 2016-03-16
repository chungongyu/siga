#ifndef bigraph_h_
#define bigraph_h_

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

class Edge {
};

typedef std::vector< Edge* > EdgePtrList;

class Vertex {
public:
    enum {
        WHITE = 0, 
        GRAY, 
        BLACK, 
        BLUE, 
        RED
    };
    typedef uint8_t Color;
    typedef std::string Id;

    Vertex(const Id& id, const std::string& seq) : id(id), seq(seq), color(WHITE), coverage(1) {
    }

    void merge(Edge* edge);
    void addEdge(Edge* edge);

    Id id;
    Color color;
    std::string seq;
    size_t coverage; // Number of vertices that have been merged into this one

private:
    EdgePtrList _edges;
};

typedef std::unordered_map< Vertex::Id, Vertex* > VertexTable;

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
