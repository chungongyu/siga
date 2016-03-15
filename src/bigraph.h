#ifndef bigraph_h_
#define bigraph_h_

#include <string>
#include <unordered_map>

class Edge {
};

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

    Id id;
    Color color;
    std::string seq;
    size_t coverage; // Number of vertices that have been merged into this one
private:
};

typedef std::unordered_map< Vertex::Id, Vertex* > VertaexTable;

class Bigraph {
public:
    Bigraph() {
    }
    ~Bigraph();
private:
    VertaexTable _vertices;
};

#endif // bigraph_h_
