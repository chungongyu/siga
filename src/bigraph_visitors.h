#ifndef bigraph_visitors_h_
#define bigraph_visitors_h_

#include <iostream>

class Vertex;
class Bigraph;

class BigraphVisitor {
public:
    virtual void previsit(Bigraph* graph) {
    }
    virtual bool visit(Bigraph* graph, Vertex* vertex) = 0;
    virtual void postvisit(Bigraph* graph) {
    }
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

#endif // bigraph_visitors_h_
