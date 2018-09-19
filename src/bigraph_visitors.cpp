#include "bigraph_visitors.h"
#include "kseq.h"
#include "reads.h"
#include "sequence_process_framework.h"
#include "utils.h"

#include <functional>
#include <unordered_set>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.BigraphVisitor"));

//
// EdgeColorVisitor
//
class EdgeColorVisitor : public BigraphVisitor {
public:
    EdgeColorVisitor(GraphColor c) : _color(c) {
    }

    void previsit(Bigraph* graph) {
    }
    bool visit(Bigraph* graph, Vertex* vertex) {
        const EdgePtrList& edges = vertex->edges();
        BOOST_FOREACH(Edge* edge, edges) {
            edge->color(_color);
        }
        return true;
    }
    void postvisit(Bigraph* graph) {
    }
private:
    GraphColor _color;
};

//
// ChimericVisitor
//
void ChimericVisitor::previsit(Bigraph* graph) {
    _chimeric = 0;
    graph->color(GC_WHITE);
}

bool ChimericVisitor::visit(Bigraph* graph, Vertex* vertex) {
    // Check if this node is chimeric
    const std::string& seq = vertex->seq();
    if (vertex->degrees(Edge::ED_SENSE) == 1 && vertex->degrees(Edge::ED_ANTISENSE) == 1 && seq.length() < _minLength) {
        Edge* prevEdge = vertex->edges(Edge::ED_ANTISENSE)[0];
        Edge* nextEdge = vertex->edges(Edge::ED_SENSE)[0];
        Vertex* prevVert = prevEdge->end();
        Vertex* nextVert = nextEdge->end();

        bool chimeric = true;
        if (chimeric) {
            chimeric &= (prevVert->degrees(Edge::ED_SENSE) >= 2);
        }
        if (chimeric) {
            chimeric &= (nextVert->degrees(Edge::ED_ANTISENSE) >= 2);
        }
        if (chimeric) {
            // smallest?
            bool smallest = false;
            {
                EdgePtrList edges = prevVert->edges(Edge::ED_SENSE);
                for (size_t k = 0; k < edges.size(); ++k) {
                    if (edges[k]->coord().length() > prevEdge->coord().length() && edges[k]->coord().length() - prevEdge->coord().length() >= _delta) {
                        smallest = true;
                    }
                }
            }
            {
                EdgePtrList edges = nextVert->edges(Edge::ED_ANTISENSE);
                for (size_t k = 0; k < edges.size(); ++k) {
                    if (edges[k]->coord().length() > nextEdge->coord().length() && edges[k]->coord().length() - nextEdge->coord().length() >= _delta) {
                        smallest = true;
                    }
                }
            }
            chimeric &= smallest;
        }
        if (chimeric) {
            vertex->color(GC_BLACK);
            ++_chimeric;
            return true;
        }
    }
    return false;
}

void ChimericVisitor::postvisit(Bigraph* graph) {
    graph->sweepVertices(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[ChimericVisitor]: Removed %d chimeric") % _chimeric);
}
//
// ContainRemoveVisitor
//
void ContainRemoveVisitor::previsit(Bigraph* graph) {
    graph->color(GC_WHITE);

    // Clear the containment flag, if any containments are added
    // during this algorithm the flag will be reset and another
    // round must be re-run
    graph->containment(false);

    _contained = 0;
}

bool ContainRemoveVisitor::visit(Bigraph* graph, Vertex* vertex) {
    if (vertex->contained()) {
        // Add any new irreducible edges that exist when pToRemove is deleted
        // from the graph
        const EdgePtrList& edges = vertex->edges();

        // Delete the edges from the graph
        for (EdgePtrList::const_iterator i = edges.begin(); i != edges.end(); ++i) {
           Edge* edge = *i;
           Edge* twin = edge->twin();
           Vertex* end = edge->end();

           end->removeEdge(twin);
           vertex->removeEdge(edge);

           SAFE_DELETE(twin);
           SAFE_DELETE(edge);
        }

        vertex->color(GC_BLACK);
        ++_contained;

        return true;
    }
    return false;
}

void ContainRemoveVisitor::postvisit(Bigraph* graph) {
    graph->sweepVertices(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[ContainRemoveVisitor] Removed %d containment vertices") % _contained);
}

//
// FastaVisitor
//
bool FastaVisitor::visit(Bigraph* graph, Vertex* vertex) {
    DNASeq seq(vertex->id(), vertex->seq());
    _stream << seq;
    return false;
}

//
// LoopRemoveVisitor
//
void LoopRemoveVisitor::previsit(Bigraph* graph) {
    _loops.clear();
}

bool LoopRemoveVisitor::visit(Bigraph* graph, Vertex* vertex) {
    bool modified = false;
    if (vertex->degrees(Edge::ED_SENSE) == 1 && vertex->degrees(Edge::ED_ANTISENSE) == 1) {
        Edge* prevEdge = vertex->edges(Edge::ED_ANTISENSE)[0];
        Edge* nextEdge = vertex->edges(Edge::ED_SENSE)[0];
        Vertex* prevVert = prevEdge->end();
        Vertex* nextVert = nextEdge->end();
        if (prevVert == nextVert) {
            //vertex->color(GC_BLACK);
            _loops.push_back(vertex);
            modified = true;
        }
    }
    return modified;
}

void LoopRemoveVisitor::postvisit(Bigraph* graph) {
    LOG4CXX_INFO(logger, boost::format("[LoopRemoveVisitor] Removed %lu loop vertices") % _loops.size());
    BOOST_FOREACH(Vertex* vertex, _loops) {
        assert(vertex->degrees(Edge::ED_SENSE) == 1 && vertex->degrees(Edge::ED_ANTISENSE) == 1);
        Edge* prevEdge = vertex->edges(Edge::ED_ANTISENSE)[0];
        Edge* nextEdge = vertex->edges(Edge::ED_SENSE)[0];
        Vertex* prevVert = prevEdge->end();
        Vertex* nextVert = nextEdge->end();
        assert(prevVert == nextVert);

        LOG4CXX_INFO(logger, boost::format("[LoopRemoveVisitor] Remove loop vertex: id=%s, coverage=%d, prev vertex: id=%s, coverage=%d") % vertex->id() % vertex->coverage() % prevVert->id() % prevVert->coverage());


        Edge* nextTwin = nextEdge->twin();
        vertex->merge(nextEdge);
        vertex->removeEdge(nextEdge);
        SAFE_DELETE(nextEdge);
        nextVert->removeEdge(nextTwin);
        SAFE_DELETE(nextTwin);

        Edge* prevTwin = prevEdge->twin();
        prevVert->merge(prevTwin);
        prevVert->removeEdge(prevTwin);
        SAFE_DELETE(prevTwin);
        vertex->removeEdge(prevEdge);
        SAFE_DELETE(prevEdge);

        graph->removeVertex(vertex);
        SAFE_DELETE(vertex);
    }
}

//
// MaximalOverlapVisitor
//
void MaximalOverlapVisitor::previsit(Bigraph* graph) {
    // The graph must not have containments
    assert(!graph->containment());

    // Set all the edges in the graph to "vacant"
    EdgeColorVisitor ecVisit(GC_WHITE);
    graph->visit(&ecVisit);

    _dummys = 0;
}

class OverlapCmp {
public:
    bool operator()(const Edge* x, const Edge* y) const {
        return x->coord().length() > y->coord().length();
    }
};

class EdgeDirCmp {
public:
    EdgeDirCmp(const Edge* edge) : _edge(edge) {
    }
    bool operator()(const Edge* edge) const {
        return _edge->twin()->dir() == (edge->dir() + 1) % Edge::ED_COUNT;
    }
private:
    const Edge* _edge;
};

bool MaximalOverlapVisitor::visit(Bigraph* graph, Vertex* vertex) {
    bool modified = false;

    for (size_t i = 0; i < Edge::ED_COUNT; ++i) {
        Edge::Dir dir = Edge::EDGE_DIRECTIONS[i];
        EdgePtrList fwdlist = vertex->edges(dir);

        std::sort(fwdlist.begin(), fwdlist.end(), OverlapCmp());

        for (size_t j = 1; j < fwdlist.size(); ++j) {
            if (fwdlist[j]->color() == GC_BLACK) {
                continue;
            }

            if (fwdlist[0]->coord().length() - fwdlist[j]->coord().length() < _delta) {
                continue;
            }

            EdgePtrList revlist = fwdlist[j]->end()->edges();
            EdgePtrList::iterator last = std::remove_if(revlist.begin(), revlist.end(), EdgeDirCmp(fwdlist[j]));
            if (last != revlist.end()) {
                revlist.resize(std::distance(revlist.begin(), last));
            }
            assert(!revlist.empty());

            std::sort(revlist.begin(), revlist.end(), OverlapCmp());

            bool largest = revlist[0]->end() == vertex;
            for (size_t k = 1; k < revlist.size() && revlist[k]->coord().length() - revlist[0]->coord().length() < _delta && !largest; ++k) {
                largest = revlist[k]->end() == vertex; 
            }
            if (largest) {
                continue;
            }

            if (dir == Edge::ED_SENSE) {
                LOG4CXX_INFO(logger, boost::format("[MaximalOverlapVisitor] remove edge %s->%s (%d)") % fwdlist[j]->start()->id() % fwdlist[j]->end()->id() % fwdlist[j]->coord().length());
            } else {
                LOG4CXX_INFO(logger, boost::format("[MaximalOverlapVisitor] remove edge %s->%s (%d)") % fwdlist[j]->end()->id() % fwdlist[j]->start()->id() % fwdlist[j]->coord().length());
            }
            fwdlist[j]->color(GC_BLACK);
            fwdlist[j]->twin()->color(GC_BLACK);
            ++_dummys;
            modified = true;
        }
    }

    return modified;
}

void MaximalOverlapVisitor::postvisit(Bigraph* graph) {
    graph->sweepEdges(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[MaximalOverlapVisitor] Removed %d dummy edges") % _dummys);
}

// 
// InsertSizeEstimateVisitor
//
void InsertSizeEstimateVisitor::previsit(Bigraph* graph) {
    graph->color(GC_GREEN);

    _samples.clear();
}

bool InsertSizeEstimateVisitor::visit(Bigraph* graph, Vertex* vertex) {
    if (vertex->color() == GC_GREEN) {

        std::unordered_map< Vertex::Id, int > distancelist = boost::assign::map_list_of(vertex->id(), 0);
        vertex->color(GC_RED);

        ////////////////////////////////////////
        // search forward & backward
        //
        // cases:
        //              c
        //      f    fc | rc   f
        //    R1-->R2-->R3-->R4-->R5
        //    ^    |^   |^   |^   |
        //     \__/ \__/ \__/ \__/
        //      r    fc   rc   r

        //              c    c
        //      f    fc | r  | rc   f
        //    R1-->R2-->R3-->R4-->R5-->R6
        //    ^    |^   |^   |^   |^   |
        //     \__/ \__/ \__/ \__/ \__/
        //      r    fc   f    rc   r
        //
        // the state machine:
        //                     ---------    
        //                    |         |
        //                    V         | f
        //    -----  f/rc  ---------    |
        //   |start| ---> | forward | --
        //    -----        ---------   
        //      |          ^     |fc   
        //      |          |rc   V       
        //       \  fc/r  ----------   
        //         ----> | backword | --
        //                ----------    |
        //                    ^         | r
        //                    |         |
        //                     ---------
        //
        ////////////////////////////////////////
        for (size_t idx = 0; idx < Edge::ED_COUNT; idx++) {
            Edge::Dir searchDir = Edge::EDGE_DIRECTIONS[idx];
            int distance = 0, flag = (searchDir == Edge::ED_SENSE ? 1 : -1);
            Vertex* p = vertex;

            while (true) {
                EdgePtrList straight = this->edges(p, searchDir); // forward & backword edges

                // Don't walk thru singular self edges
                if (straight.size() != 1 || straight[0]->isSelf() || straight[0]->end()->color() == GC_RED) {
                    break;
                }
                Edge* single = straight[0];
                Edge* twin = single->twin();
                Vertex* end = single->end();
                EdgePtrList opposite = this->edges(end, twin->dir());
                if (opposite.size() != 1) {
                    break;
                }

                p = end;
                if (searchDir == Edge::ED_SENSE) {      // forward
                    distance += flag * single->coord().complement().length();
                } else {                                // backward
                    distance += flag * twin->coord().complement().length();
                }
                if (single->comp() == Edge::EC_REVERSE) {
                    searchDir = Edge::EDGE_DIRECTIONS[Edge::EC_COUNT - searchDir - 1];
                }
                distancelist[p->id()] = distance;
                LOG4CXX_DEBUG(logger, boost::format("InsertSizeEstimateVisitor::visit addVertex(%s,%s,%d)") % vertex->id() % p->id() % distance);

                p->color(GC_RED);
            }
        }

        for (std::unordered_map< Vertex::Id, int >::const_iterator i = distancelist.begin(); i != distancelist.end(); ++i) {
            Vertex::Id  pairId = PairEnd::id(i->first);
            if (i->first < pairId) {
                std::unordered_map< Vertex::Id, int >::const_iterator j = distancelist.find(pairId);
                if (j != distancelist.end()) {
                    size_t distance = abs(j->second - i->second);
                    _samples.push_back(distance);
                    LOG4CXX_DEBUG(logger, boost::format("InsertSizeEstimateVisitor::visit addLink(%s,%s,%d)") % i->first % j->first % distance);
                }
            }
        }
    }
    return false;
}

void InsertSizeEstimateVisitor::postvisit(Bigraph* graph) {
    LOG4CXX_INFO(logger, boost::format("InsertSizeEstimateVisitor::samples=%d") % _samples.size());

    typedef boost::accumulators::accumulator_set< double, boost::accumulators::stats< boost::accumulators::tag::count, boost::accumulators::tag::mean, boost::accumulators::tag::moment< 2 > > > Accumulator;
    Accumulator acc;
    std::for_each(_samples.begin(), _samples.end(), std::ref(acc));
    if (boost::accumulators::count(acc) > 0) {
        _average = (size_t)boost::accumulators::mean(acc);
        _delta = std::sqrt(
                boost::accumulators::moment< 2 >(acc) - std::pow(boost::accumulators::mean(acc), 2)
                );

        LOG4CXX_INFO(logger, boost::format("InsertSizeEstimateVisitor::average=%d, delta=%d") % _average % _delta);
    }

    graph->color(GC_GREEN);
}

EdgePtrList InsertSizeEstimateVisitor::edges(const Vertex* vertex, Edge::Dir dir) { // reduced edges
    EdgePtrList edges = vertex->edges(dir);
    std::sort(edges.begin(), edges.end(), OverlapCmp());
    {
        size_t k = 0;
        for (size_t i = 0; i < edges.size(); ++i) {
            if (edges[i]->coord().length() != edges[i]->coord().seqlen) {
                edges[k++] = edges[i];
            }
        }
        edges.resize(k);
    }
    {
        size_t k = 0;
        for (size_t i = 1; i < edges.size(); ++i) {
            if (edges[i]->coord().length() != edges[k]->coord().length() || edges[i]->label() != edges[k]->label()) {
                edges[++k] = edges[i];
            }
        }
        if (!edges.empty()) {
            edges.resize(k + 1);
        }
    }
    return edges;
}

// 
// BigraphSearchTree
//
class BigraphSearchTree {
public:
    class Node {
    public:
        Node(const Vertex* vertex, int distance) : vertex(vertex), distance(distance) {
        }
        ~Node() {
        }

        const Vertex* vertex;
        int distance;
    };
    typedef std::shared_ptr< Node > NodePtr;
    typedef std::vector< NodePtr > NodePtrList;

    struct NodePtrCmp {
        bool operator()(const NodePtr& x, const NodePtr& y) const {
            return x->vertex->id() == y->vertex->id() && x->distance == y->distance;
        }
    };

    template< class T >
    struct NodePtrHash {
        size_t operator()(const T& node) const {
            std::hash<std::string> hasher;
            return hasher(node->vertex->id());
        }
    };

    static void buildNew(Bigraph* graph, const Vertex* start, size_t maxDistance, size_t maxNodes, NodePtrList* leaves) {
        const Vertex* paired_v1 = graph->getVertex(PairEnd::id(start->id()));
        assert(paired_v1 != NULL);

        typedef std::deque< NodePtr > NodePtrQueue;

        std::unordered_set< NodePtr, NodePtrHash< NodePtr >, NodePtrCmp > visited;

        NodePtrQueue Q = boost::assign::list_of(NodePtr(new Node(start, 0)));
        while (!Q.empty() && leaves->size() < maxNodes && Q.size() < 1000*maxNodes) {
            NodePtr curr = Q.front();
            Q.pop_front();

            if (visited.find(curr) != visited.end()) {
                continue;
            }
            visited.insert(curr);

            if (abs(curr->distance) < maxDistance) {
                if (curr->distance != 0 && curr->vertex != start) {
                    const Vertex* paired_v2 = graph->getVertex(PairEnd::id(curr->vertex->id()));
                    assert(paired_v2 != NULL);
                    LOG4CXX_DEBUG(logger, boost::format("vertex1: %s<->%s, vertex2: %s<->%s\n") % start->id() % paired_v1->id() % curr->vertex->id() % paired_v2->id());

                    leaves->push_back(curr);
                }

                const EdgePtrList& edges = curr->vertex->edges(Edge::ED_SENSE);
                for (EdgePtrList::const_iterator i = edges.begin(); i != edges.end(); ++i) {
                    Edge* edge = *i;
                    int distance = curr->distance + edge->coord().interval.start;
                    NodePtr child(new Node(edge->end(), distance));
                    Q.push_back(child);
                }
            }
        }
    }

    static void build(const Vertex* start, const Vertex* end, Edge::Dir searchDir, size_t maxDistance, size_t maxNodes, NodePtrList* leaves) {
        typedef std::deque< NodePtr > NodePtrQueue;

        std::unordered_set< NodePtr, NodePtrHash< NodePtr >, NodePtrCmp > visited;

        NodePtrQueue Q = boost::assign::list_of(NodePtr(new Node(start, 0)));
        while (!Q.empty() && leaves->size() < maxNodes && Q.size() < 5*maxDistance) {
            NodePtr curr = Q.front();
            Q.pop_front();

            if (visited.find(curr) != visited.end()) {
                continue;
            }
            visited.insert(curr);

            if (abs(curr->distance) < maxDistance) {
                if (end == NULL) {
                    if (curr->distance != 0 && curr->vertex != start) {
                        leaves->push_back(curr);
                    }
                } else if (end->id() == curr->vertex->id()) {
                    leaves->push_back(curr);
                    return;
                }

                const EdgePtrList& edges = curr->vertex->edges();
                for (EdgePtrList::const_iterator i = edges.begin(); i != edges.end(); ++i) {
                    Edge* edge = *i;
                    if (edge->dir() == searchDir) {
                        int distance = curr->distance;
                        if (searchDir == Edge::ED_SENSE) {
                            distance += edge->coord().interval.start;
                        } else {
                            distance -= edge->twin()->coord().interval.start;
                        }
                        NodePtr child(new Node(edge->end(),  distance));
                        Q.push_back(child);
                    }
                }
            }
        }
    }
};

//
// PairedReadVisitor
//

template< class T >
struct DistanceListCmp {
    bool operator()(const std::pair< T, int >& x, const std::pair< T, int >& y) const {
        return x.second < y.second;
    }
};

PairedReadVisitor::DistanceList PairedReadVisitor::VertexProcess::process(const Vertex* vertex1) {
    PairedReadVisitor::DistanceList linklist;

    const Vertex* paired_v1 = _graph->getVertex(PairEnd::id(vertex1->id()));
    assert(paired_v1 != NULL);

    BigraphSearchTree::NodePtrList adjacents;
    BigraphSearchTree::build(vertex1, NULL, Edge::ED_SENSE, vertex1->seq().length() - _visitor->_minOverlap, _visitor->_maxNodes, &adjacents);
    //BigraphSearchTree::build(vertex1, NULL, Edge::ED_SENSE, _visitor->_maxDistance, _visitor->_maxNodes, &adjacents);

    size_t numNodes = 0;
    BOOST_FOREACH(const BigraphSearchTree::NodePtr& node1, adjacents) {
        if (numNodes >= 3) break;
        const Vertex* paired_v2 = _graph->getVertex(PairEnd::id(node1->vertex->id()));
        assert(paired_v2 != NULL);
        LOG4CXX_DEBUG(logger, boost::format("vertex1: %s<->%s, vertex2: %s<->%s") % vertex1->id() % paired_v1->id() % node1->vertex->id() % paired_v2->id());

        BigraphSearchTree::NodePtrList faraways;
        BigraphSearchTree::build(paired_v1, paired_v2, Edge::ED_SENSE, _visitor->_maxDistance, 1, &faraways);
        BigraphSearchTree::build(paired_v1, paired_v2, Edge::ED_ANTISENSE, _visitor->_maxDistance, 2, &faraways);
        BOOST_FOREACH(const BigraphSearchTree::NodePtr& node2, faraways) {
            if (numNodes >= 3) break;
            linklist.push_back(std::make_pair(node1->vertex, node1->distance));
            LOG4CXX_DEBUG(logger, boost::format("paired_read_all\t%s\t%s\t%d\t%s\t%s\t%d") % vertex1->id() % node1->vertex->id() % node1->distance % paired_v1->id() % node2->vertex->id() % node2->distance);
            ++numNodes;
        }
    }

    return linklist;
}

void PairedReadVisitor::VertexPostProcess::process(const Vertex* vertex1, DistanceList& linklist) {
    std::sort(linklist.begin(), linklist.end(), DistanceListCmp< const Vertex* >());
    for (size_t i = 0; i < linklist.size(); ++i) {
        const std::pair< const Vertex*, int >& xi = linklist[i];
        addLink(vertex1->id(), xi.first->id(), xi.second, _links);
//#if 0
        for (size_t j = i + 1; j < linklist.size(); ++j) {
            const std::pair< const Vertex*, int >& xj = linklist[j];
            const std::string& seqi = xi.first->seq();
            const std::string& seqj = xj.first->seq();

            if (xi.second != xj.second && boost::algorithm::equals(seqi.substr(xj.second - xi.second), seqj.substr(0, seqi.length() - (xj.second - xi.second)))) {
                addLink(xi.first->id(), xj.first->id(), xj.second - xi.second, _links);
            } else {
                LOG4CXX_DEBUG(logger, boost::format("paired_read_branch\t%s\t%s\t%d\t%s\t%d") % vertex1->id() % xi.first->id() % xi.second % xj.first->id() % xj.second);
            }
        }
//#endif
    }
}

void PairedReadVisitor::VertexPostProcess::addLink(const Vertex::Id& v1, const Vertex::Id& v2, int distance, LinkList* links) {
    LinkList::iterator i = links->find(v1);
    if (i == links->end()) {
        DistanceMap tbl = boost::assign::map_list_of(v2, distance); 
        (*links)[v1] = tbl;
    } else {
        DistanceMap::const_iterator k = i->second.find(v2);
        if (k != i->second.end() && k->second != distance) {
            LOG4CXX_WARN(logger, boost::format("PairedReadVisitor::addLink(%s, %s, %d) != %d") % v1 % v2 % distance % k->second);
        }
        if (k == i->second.end() || k->second > distance) {
            i->second[v2] = distance;
        }
    }
}

void PairedReadVisitor::previsit(Bigraph* graph) {
    _vertices.clear();
}

bool PairedReadVisitor::visit(Bigraph* graph, Vertex* vertex1) {
    const Vertex* paired_v1 = graph->getVertex(PairEnd::id(vertex1->id()));
    assert(paired_v1 != NULL);

    //if (PairEnd::id(paired_v1->id()) < PairEnd::id(vertex1->id())) {
        _vertices.push_back(vertex1);
    //}
    return false;
}

void PairedReadVisitor::postvisit(Bigraph* graph) {
    LinkList links;

    VertexGenerator generator(_vertices);
    std::vector< VertexProcess* > proclist(_threads);
    VertexPostProcess postproc(&links);
    for (size_t i = 0; i < proclist.size(); ++i) {
        proclist[i] = new VertexProcess(graph, this);
    }
#ifdef _OPENMP
    if (_threads > 1) {
        SequenceProcessFramework::ParallelWorker<
            const Vertex*, 
            DistanceList, 
            VertexGenerator, 
            VertexProcess, 
            VertexPostProcess
            > worker;
        size_t num = worker.run(generator, &proclist, &postproc, _batch);
    } else { // single thread
#endif // _OPENMP
        SequenceProcessFramework::SerialWorker<
            const Vertex*, 
            DistanceList, 
            VertexGenerator, 
            VertexProcess, 
            VertexPostProcess
            > worker;
        size_t num = worker.run(generator, proclist[0], &postproc);
#ifdef _OPENMP
    }
#endif // _OPENMP
    for (size_t i = 0; i < proclist.size(); ++i) {
        delete proclist[i];
    }

    EdgeColorVisitor ecVist(GC_BLACK);
    graph->visit(&ecVist);
    //graph->color(GC_BLACK);

    // simplify
    for (LinkList::const_iterator i = links.begin(); i != links.end(); ++i) {
        std::vector< std::pair< Vertex::Id, int > > nodelist;

        // sorted by distance
        std::copy(i->second.begin(), i->second.end(), std::back_inserter(nodelist));
        std::sort(nodelist.begin(), nodelist.end(), DistanceListCmp< Vertex::Id >());

        const std::pair< Vertex::Id, int > x0 = nodelist[0];
        LOG4CXX_DEBUG(logger, boost::format("paired_read_simplify\t%s\t%s\t%d") % i->first % x0.first % x0.second);
        addEdge(i->first, x0.first, x0.second, graph);

        for (size_t j = 1; j < nodelist.size(); ++j) {
            const std::pair< Vertex::Id, int >& xj = nodelist[j];
            LinkList::const_iterator v = links.find(x0.first);
            DistanceMap::const_iterator w;
            if (v != links.end() && (w = v->second.find(xj.first)) != v->second.end() && x0.second + w->second == xj.second) {
                continue;
            }
            Vertex* v0 = graph->getVertex(x0.first);
            Vertex* vj = graph->getVertex(xj.first);
            const std::string& seq0 = v0->seq();
            const std::string& seqj = vj->seq();

            if (boost::algorithm::equals(seq0.substr(xj.second - x0.second), seqj.substr(0, seq0.length() - (xj.second - x0.second)))) {
                continue;
            }

            LOG4CXX_DEBUG(logger, boost::format("paired_read_simplify\t%s\t%s\t%d") % i->first % xj.first % xj.second);
            addEdge(i->first, xj.first, xj.second, graph);
        }
    }

    //LOG4CXX_INFO(logger, boost::format("[PairedReadVisitor] Removed %d dummy edges") % _dummys);

    graph->sweepEdges(GC_BLACK);
}

void PairedReadVisitor::addEdge(const Vertex::Id& v1, const Vertex::Id& v2, int distance, Bigraph* graph) {
    Vertex* vertex[2] = {
        graph->getVertex(v1), 
        graph->getVertex(v2), 
    };

    Vertex* vertex1 = graph->getVertex(v1);
    Vertex* vertex2 = graph->getVertex(v2);

    bool found = false;
    const EdgePtrList& edges = vertex1->edges();
    BOOST_FOREACH(Edge* edge, edges) {
        if (edge->dir() == Edge::ED_SENSE && edge->end() == vertex2) {
            if (edge->coord().interval.start == distance) {
                edge->color(GC_WHITE);
                edge->twin()->color(GC_WHITE);
                found = true;
                break;
            }
            LOG4CXX_WARN(logger, boost::format("PairedReadVisitor::addEdge(%s, %s, %d) != %d") % v1 % v2 % distance % edge->coord().interval.start);
        }
    }
    if (!found) {
        Edge* edges[2] = {0};
        const std::string& seq1 = vertex[0]->seq();
        const std::string& seq2 = vertex[1]->seq();
        SeqCoord coord[2] = {
            SeqCoord(distance, seq1.length() - 1, seq1.length()),
            SeqCoord(0, seq1.length() - distance - 1, seq2.length())
        };

        Edge::Dir dir[2] = {
            Edge::ED_SENSE, 
            Edge::ED_ANTISENSE
        };
        for (size_t i = 0; i < 2; ++i) {
            edges[i] = new Edge(vertex[1 - i], dir[i], Edge::EC_SAME, coord[i]);
            edges[i]->color(GC_WHITE);
        }
        edges[0]->twin(edges[1]);
        edges[1]->twin(edges[0]);
        for (size_t i = 0; i < 2; ++i) {
            graph->addEdge(vertex[i], edges[i]);
        }
    }
}

//
// SmoothingVisitor
//
void SmoothingVisitor::previsit(Bigraph* graph) {
    graph->color(GC_WHITE);
    _simple = 0;
    _complex = 0;
}

bool SmoothingVisitor::visit(Bigraph* graph, Vertex* vertex) {
    return false;
}

void SmoothingVisitor::postvisit(Bigraph* graph) {
    graph->sweepVertices(GC_RED);
    LOG4CXX_INFO(logger, boost::format("[SmoothingVisitor] Removed %s simple and %d complex bubbles") % _simple % _complex);
}

// 
// StatisticsVisitor
//
void StatisticsVisitor::previsit(Bigraph*) {
    _terminal = 0;
    _island = 0;
    _monobranch = 0;
    _dibranch = 0;
    _simple = 0;
    _edges = 0;
    _vertics = 0;
}

bool StatisticsVisitor::visit(Bigraph* graph, Vertex* vertex) {
    size_t fdegrees = vertex->degrees(Edge::ED_SENSE);
    size_t rdegrees = vertex->degrees(Edge::ED_ANTISENSE);

    if (fdegrees == 0 && rdegrees == 0) {
        ++_island;
    } else if (fdegrees == 0 || rdegrees == 0) {
        ++_terminal;
    }

    if (fdegrees > 1 && rdegrees > 1) {
        ++_dibranch;
    } else if (fdegrees > 1 || rdegrees > 1) {
        ++_monobranch;
    }

    if (fdegrees == 1 || rdegrees == 1) {
        ++_simple;
    }

    _edges += fdegrees + rdegrees;
    ++_vertics;

    return false;
}

void StatisticsVisitor::postvisit(Bigraph*) {
    LOG4CXX_INFO(logger, boost::format("[StatisticsVisitor] Vertices: %d Edges: %d Islands: %d Tips: %d Monobranch: %d Dibranch: %d Simple: %d") % _vertics % _edges % _island % _terminal % _monobranch % _dibranch % _simple);
}

//
// TrimVisitor
//
void TrimVisitor::previsit(Bigraph* graph) {
    _island = 0;
    _terminal = 0;
    graph->color(GC_WHITE);
}

bool TrimVisitor::visit(Bigraph* graph, Vertex* vertex) {
    bool modified = false;

    const std::string& seq = vertex->seq();
    if (vertex->degrees() == 0) {
        // Is an island, remove if the sequence length is less than the threshold
        if (seq.length() < _minLength) {
            LOG4CXX_TRACE(logger, boost::format("[TrimVisitor] island %s %d") % vertex->id() % seq.length());
            vertex->color(GC_BLACK);
            ++_island;
            modified = true;
        }
    } else {
        // Check if this node is a dead-end
        for (size_t idx = 0; idx < Edge::ED_COUNT; idx++) {
            Edge::Dir dir = Edge::EDGE_DIRECTIONS[idx];
            if (vertex->degrees(dir) == 0 && seq.length() < _minLength) {
                LOG4CXX_TRACE(logger, boost::format("[TrimVisitor] terminal %s %d") % vertex->id() % seq.length());
                vertex->color(GC_BLACK);
                ++_terminal;
                modified = true;
                break;
            }
        }
    }

    return modified;
}

void TrimVisitor::postvisit(Bigraph* graph) {
    graph->sweepVertices(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[TrimVisitor] Removed %d island and %d dead-end short vertices") % _island % _terminal);
}
