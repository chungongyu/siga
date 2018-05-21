#ifndef bigraph_visitors_h_
#define bigraph_visitors_h_

#include "bigraph.h"

#include <iostream>
#include <unordered_map>

class BigraphVisitor {
public:
    virtual void previsit(Bigraph* graph) {
    }
    virtual bool visit(Bigraph* graph, Vertex* vertex) = 0;
    virtual void postvisit(Bigraph* graph) {
    }
};

// Detects and removes small "chimeric" vertices from the graph
// when they are less than minLength in size
class ChimericVisitor : public BigraphVisitor {
public:
    ChimericVisitor(size_t minLength, size_t delta) : _minLength(minLength), _delta(delta) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _minLength;

    size_t _delta;
    size_t _chimeric;
};

// Remove contained vertices from the graph
class ContainRemoveVisitor : public BigraphVisitor {
public:
    ContainRemoveVisitor() : _contained(0) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _contained;
};

// Visit each node, writing it to a file as a fasta record
class FastaVisitor : public BigraphVisitor {
public:
    FastaVisitor(std::ostream& stream) : _stream(stream) {
    }
    bool visit(Bigraph* graph, Vertex* vertex);
private:
    std::ostream& _stream;
};

// Run the YU LIN's maximal overlap algorithm on each node
class MaximalOverlapVisitor : public BigraphVisitor {
public:
    MaximalOverlapVisitor(size_t delta) : _delta(delta), _dummys(0) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    static bool isSenseEdge(const Edge* edge);
    static bool isAntiSenseEdge(const Edge* edge);

    size_t _delta;
    size_t _dummys;
};

// Visit each paired node via zigzag, sweep false positive edges.
class PairedReadVisitor : public BigraphVisitor {
public:
    PairedReadVisitor(size_t minOverlap, size_t maxDistance, size_t maxNodes, size_t threads=1, size_t batch=1000) : _minOverlap(minOverlap), _maxDistance(maxDistance), _maxNodes(maxNodes), _threads(threads), _batch(batch) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    typedef std::vector< std::pair< const Vertex*, int > > DistanceList;
    typedef std::unordered_map< Vertex::Id, int > DistanceMap;
    typedef std::unordered_map< Vertex::Id, DistanceMap > LinkList;

    class VertexGenerator {
    public:
        VertexGenerator(const std::vector< const Vertex* >& vertices) : _vertices(vertices), _consumed(0) {
        }

        bool generate(const Vertex*& item) {
            if (_consumed < _vertices.size()) {
                item = _vertices[_consumed];
                ++_consumed;
                return true;
            }
            return false;
        }
        size_t consumed() const {
            return _consumed;
        }
    private:
        const std::vector< const Vertex* >& _vertices;
        size_t _consumed;
    };
    class VertexProcess {
    public:
        VertexProcess(Bigraph* graph, PairedReadVisitor* vistor) : _graph(graph), _visitor(vistor) {
        }
        DistanceList process(const Vertex* vertex);
    private:
        Bigraph* _graph;
        PairedReadVisitor* _visitor;
    };
    class VertexPostProcess {
    public:
        VertexPostProcess(LinkList* links) : _links(links) {
        }
        void process(const Vertex* vertex, DistanceList& links);
    private:
        void addLink(const Vertex::Id& v1, const Vertex::Id& v2, int distance, LinkList* links);

        LinkList* _links;
    };

    size_t _minOverlap;
    size_t _maxDistance;
    size_t _maxNodes;
    size_t _threads;
    size_t _batch;

    void addEdge(const Vertex::Id& v1, const Vertex::Id& v2, int distance, Bigraph* graph);

    std::vector< const Vertex* > _vertices;
    friend class VertexProcess;
};


// Smooth out variation in the graph
class SmoothingVisitor : public BigraphVisitor {
public:
    SmoothingVisitor() {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _simple;
    size_t _complex;
};

// Compile summary statistics for the graph 
class StatisticsVisitor : public BigraphVisitor {
public:
    StatisticsVisitor() {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _terminal;
    size_t _island;
    size_t _monobranch;
    size_t _dibranch;
    size_t _simple;
    size_t _edges;
    size_t _vertics;
};

// Detects and removes small "tip" vertices from the graph
// when they are less than minLength in size
class TrimVisitor : public BigraphVisitor {
public:
    TrimVisitor(size_t minLength) : _minLength(minLength) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _minLength;

    size_t _island;
    size_t _terminal;
};

#endif // bigraph_visitors_h_
