#include "bigraph_visitors.h"

#include <functional>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include <log4cxx/logger.h>

#include <ssw_cpp.h>

#include "asqg.h"
#include "bigraph_search.h"
#include "fmindex.h"
#include "kseq.h"
#include "parallel_framework.h"
#include "reads.h"
#include "utils.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.BigraphVisitor"));

//
// Helpers
//
namespace Repeatness {
  double calc(const Vertex* vertex, size_t n, size_t G) {
    double delta = (double)vertex->seq().length();
    double k = (double)vertex->coverage();
    return delta*n/G - k*std::log(2.0);
  }
};  // namespace Repeatness

namespace Point {
  double avg(size_t c, size_t l) {
    return (double)(std::max(c, (size_t)1) - 1)/std::max(l, (size_t)1);
  }
  double avg(const Vertex* vertex) {
    return avg(vertex->coverage(), vertex->seq().length());
  }
};  // namespace Point

//
// EdgeColorVisitor
//
class EdgeColorVisitor : public BigraphVisitor {
public:
    EdgeColorVisitor(GraphColor c, bool twin=false) : _color(c), _twin(twin) {
    }
    EdgeColorVisitor(GraphColor c, std::function<bool(const Vertex*, const Edge* edge)> filter, bool twin=false) : _color(c), _filter(filter), _twin(twin) {
    }

    void previsit(Bigraph* graph) {
    }
    bool visit(Bigraph* graph, Vertex* vertex) {
        bool modified = false;
        const EdgePtrList& edges = vertex->edges();
        for (auto& edge : edges) {
            if (!_filter || _filter(vertex, edge)) {
                edge->color(_color);
                if (_twin) {
                    edge->twin()->color(_color);
                }
                modified = true;
            }
        }
        return modified;
    }
    void postvisit(Bigraph* graph) {
    }
private:
    GraphColor _color;
    bool _twin;
    std::function<bool(const Vertex*, const Edge*)> _filter;
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
    if (vertex->degrees(Edge::ED_SENSE) == 1 && vertex->degrees(Edge::ED_ANTISENSE) == 1 && seq.length() <= _minLength && Point::avg(vertex) <= Point::avg(_minCoverage, _minLength)) {
        Edge* prevEdge = vertex->edges(Edge::ED_ANTISENSE)[0];
        Edge* nextEdge = vertex->edges(Edge::ED_SENSE)[0];
        Vertex* prevVert = prevEdge->end();
        Vertex* nextVert = nextEdge->end();

        size_t n = _N > 0 ? _N : 1751447;
        size_t G = _G > 0 ? _G : 59128983;

        bool chimeric = true;
        if (chimeric) {
            chimeric &= (prevVert->degrees(Edge::ED_SENSE) >= 2);
        }
        if (chimeric) {
            size_t k = prevVert->coverage();
            size_t delta = prevVert->seq().length();
            double score = (n-k)*(log(G-delta)-log(G>2*delta ? G-2*delta : 0.001)) - k*log(2.0);
            //LOG4CXX_INFO(logger, boost::format("ChimericVisitor::vertex=%s,prevVert=%s,delta=%ld,score=%lf") % vertex->id() % prevVert->id() % delta % score);
            //chimeric &= (score >= _T);
        }
        if (chimeric) {
            chimeric &= (nextVert->degrees(Edge::ED_ANTISENSE) >= 2);
        }
        if (chimeric) {
            size_t k = nextVert->coverage();
            size_t delta = nextVert->seq().length();
            double score = (n-k)*(log(G-delta)-log(G>2*delta ? G-2*delta : 0.001)) - k*log(2.0);
            //LOG4CXX_INFO(logger, boost::format("ChimericVisitor::vertex=%s,nextVert=%s,delta=%ld,score=%lf") % vertex->id() % nextVert->id() % delta % score);
            //chimeric &= (score > _T);
        }
        if (chimeric) {
            auto smallest = [&](const EdgePtrList& edges, const Edge* edge)->bool {
                for (size_t k = 0; k < edges.size(); ++k) {
                    if (edges[k]->coord().length() < edge->coord().length() || edges[k]->coord().length() - edge->coord().length() < _delta) {
                        return false;
                    }
                }
               return true;
            };
            auto smallest_length = [&](const EdgePtrList& edges, const Edge* edge)->bool {
                for (size_t k = 0; k < edges.size(); ++k) {
                    if (edges[k]->end()->id() == vertex->id()) {
                        continue;
                    }
                    if (edges[k]->end()->seq().length() <= seq.length() + _delta) {
                        //LOG4CXX_INFO(logger, boost::format("ChimericVisitor::vertex=%s,xx=%s,seq.length=%d,xx.length=%d,delta=%ld") % vertex->id() % edges[k]->end()->id() % seq.length() % edges[k]->end()->seq().length() % _delta);
                        return false;
                    }
                }
               return true;
            };
            auto smallest_coverage = [&](const EdgePtrList& edges, const Edge* edge)->bool {
                for (size_t k = 0; k < edges.size(); ++k) {
                    if (edges[k]->end()->id() == vertex->id()) {
                        continue;
                    }
                    if (edges[k]->end()->coverage() <= vertex->coverage() + 3) {
                        //LOG4CXX_INFO(logger, boost::format("ChimericVisitor::smallest_coverage::vertex=%s,xx=%s,seq.length=%d,xx.length=%d,delta=%ld,xx.coverage=%ld,seq.coverage=%ld") % vertex->id() % edges[k]->end()->id() % seq.length() % edges[k]->end()->seq().length() % _delta % edges[k]->end()->coverage() % vertex->coverage());
                        return false;
                    }
                }
               return true;
            };
            auto smallest_new = [&](const EdgePtrList& edges, const Edge* edge)->bool {
                Vertex* linkVert = edge->end();
                size_t k = linkVert->coverage();
                size_t delta = linkVert->seq().length();
                double score = (n-k)*(log(G-delta)-log(G>2*delta ? G-2*delta : 0.001)) - k*log(2.0);
                //LOG4CXX_INFO(logger, boost::format("ChimericVisitor::vertex=%s,linkVert=%s,delta=%ld,score=%lf") % vertex->id() % linkVert->id() % delta % score);
                if (score < _T) {
                    return false;
                }
               return smallest_length(edges, edge) || smallest_coverage(edges, edge);
            };
            // smallest?
            // bool smallest = false;
            // {
            //     EdgePtrList edges = prevVert->edges(Edge::ED_SENSE);
            //     for (size_t k = 0; k < edges.size(); ++k) {
            //         if (edges[k]->coord().length() > prevEdge->coord().length() && edges[k]->coord().length() - prevEdge->coord().length() >= _delta) {
            //             smallest = true;
            //         }
            //     }
            // }
            // {
            //     EdgePtrList edges = nextVert->edges(Edge::ED_ANTISENSE);
            //     for (size_t k = 0; k < edges.size(); ++k) {
            //         if (edges[k]->coord().length() > nextEdge->coord().length() && edges[k]->coord().length() - nextEdge->coord().length() >= _delta) {
            //             smallest = true;
            //         }
            //     }
            // }
            chimeric &= (smallest_new(prevVert->edges(Edge::ED_SENSE), prevEdge) || smallest_new(nextVert->edges(Edge::ED_ANTISENSE), nextEdge));
            //LOG4CXX_INFO(logger, boost::format("ChimericVisitor::vertex=%s,smallest=%d,smallest=%d,chimeric=%d") % vertex->id() % smallest_new(prevVert->edges(Edge::ED_SENSE), prevEdge) % smallest_new(nextVert->edges(Edge::ED_ANTISENSE), nextEdge) % chimeric);
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
        for (auto edge : edges) {
           edge->color(GC_NONE);
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
    std::string index = vertex->index();
    if (vertex->coverage() > 1) {
        ASQG::IntTagValue cr(vertex->coverage());
        if (!seq.comment.empty()) {
            seq.comment += ' ';
        }
        seq.comment += cr.tostring(ASQG::COVERAGE_TAG);
    }
    if (!index.empty()) {
        ASQG::StringTagValue bx(index);
        if (!seq.comment.empty()) {
            seq.comment += ' ';
        }
        seq.comment += bx.tostring(ASQG::BARCODE_TAG);
    }
    std::string ext = vertex->extension();
    if (!ext.empty()) {
        ASQG::StringTagValue extension(ext);
        if (!seq.comment.empty()) {
            seq.comment += ' ';
        }
        seq.comment += extension.tostring(ASQG::EXTENSION_TAG);
    }
    _stream << seq;
    return false;
}

//
// IdenticalRemoveVisitor
//
void IdenticalRemoveVisitor::previsit(Bigraph* graph) {
    graph->color(GC_WHITE);
    _count = 0;
}

bool IdenticalRemoveVisitor::visit(Bigraph* graph, Vertex* vertex) {
    bool modified = false;
    if (vertex->contained()) {
        // Check if this vertex is identical to any other vertex
        const EdgePtrList& edges = vertex->edges();
        for (const auto& edge : edges) {
            Vertex* other = edge->end();
            if (vertex->seq().length() != other->seq().length()) {
                continue;
            }
            Overlap ovr(edge->start()->id(), edge->end()->id(), edge->match());
            if (!ovr.isContainment() || ovr.containedIdx() != 0) {
                continue;
            }
            if (vertex->seq() == other->seq()) {
                vertex->color(GC_BLACK);
                LOG4CXX_DEBUG(logger, boost::format("[IdenticalRemoveVisitor] Remove identical vertex: %s/%s") % vertex->id() % other->id());
                ++_count;
                break;
            }
        }
    }
    return modified;
}

void IdenticalRemoveVisitor::postvisit(Bigraph* graph) {
    graph->sweepVertices(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[IdenticalRemoveVisitorRemoveVisitor] Removed %lu identical vertices") % _count);
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
        if (!prevEdge->isSelf() && !nextEdge->isSelf() && prevVert == nextVert) {
            //vertex->color(GC_BLACK);
            _loops.push_back(vertex);
            modified = true;
        }
    }
    return modified;
}

void LoopRemoveVisitor::postvisit(Bigraph* graph) {
    /////////////////////////////////////////
    // R4 is a Loop node
    //
    //   R1-->R2-->R3
    //      /  \
    //      \  /
    //       R4
    //
    // We change it to:
    //   R1-->R2-->R4-->R2-->R3
    //
    /////////////////////////////////////////
    for (auto& vertex : _loops) {
        assert(vertex->degrees(Edge::ED_SENSE) == 1 && vertex->degrees(Edge::ED_ANTISENSE) == 1);
        Edge* prevEdge = vertex->edges(Edge::ED_ANTISENSE)[0];
        Edge* nextEdge = vertex->edges(Edge::ED_SENSE)[0];
        Vertex* prevVert = prevEdge->end();
        Vertex* nextVert = nextEdge->end();
        assert(prevVert == nextVert);

        LOG4CXX_INFO(logger, boost::format("[LoopRemoveVisitor] Remove loop vertex: id=%s, coverage=%d, seq=%d, prev vertex: id=%s, coverage=%d, seq=%d") % vertex->id() % vertex->coverage() % vertex->seq().length() % prevVert->id() % prevVert->coverage() % prevVert->seq().length());

        /////////////////////////////////////////
        //
        // R1-->R2-->R3
        //       \
        //        R4-->R2
        //
        ////////////////////////////////////////
        Edge* nextTwin = nextEdge->twin();
        vertex->merge(nextEdge);
        vertex->removeEdge(nextEdge);
        SAFE_DELETE(nextEdge);
        nextVert->removeEdge(nextTwin);
        SAFE_DELETE(nextTwin);

        /////////////////////////////////////////
        //
        //   R1-->R2-->R4-->R2-->R3
        //
        /////////////////////////////////////////
        Edge* prevTwin = prevEdge->twin();
        std::string label = prevTwin->label();
        bool prepend = false;
        if (prevTwin->dir() == Edge::ED_ANTISENSE) {
            prepend = true;
        }
        prevVert->merge(prevTwin);
        {
            EdgePtrList transEdges = prevVert->edges(Edge::EDGE_DIRECTIONS[Edge::ED_COUNT - prevEdge->dir() - 1]);
            for (auto& transEdge : transEdges) {
                if (transEdge != prevTwin && !prepend) {
                    SeqCoord& coord = transEdge->coord();
                    coord.interval.offset(label.length());
                }
            }
        }
        prevVert->removeEdge(prevTwin);
        SAFE_DELETE(prevTwin);
        vertex->removeEdge(prevEdge);
        SAFE_DELETE(prevEdge);

        graph->removeVertex(vertex);
        SAFE_DELETE(vertex);
    }
    LOG4CXX_INFO(logger, boost::format("[LoopRemoveVisitor] Removed %lu loop vertices") % _loops.size());
}

//
// MaximumOverlapVisitor
//
void MaximumOverlapVisitor::previsit(Bigraph* graph) {
  // The graph must not have containments
  assert(!graph->containment());

  // Set all the edges in the graph to "vacant"
  EdgeColorVisitor ecVisit(GC_WHITE, true);
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

bool MaximumOverlapVisitor::visit(Bigraph* graph, Vertex* vertex) {
  {
    size_t n = _N > 0 ? _N : 1751447;
    size_t G = _G > 0 ? _G : 59128983;
    size_t k = vertex->coverage();
    size_t delta = vertex->seq().length();
  
    double score = (n-k)*(log(G-delta)-log(G>2*delta ? G-2*delta : 0.001)) - k*log(2.0);
    if (score < _T) {
        return false;
    }
  }
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

      if (_carefully) {
        if (!fwdlist[j]->isSelf()) {  // Not a self
          EdgePtrList revlist = fwdlist[j]->end()->edges();
          auto last = std::remove_if(revlist.begin(), revlist.end(), EdgeDirCmp(fwdlist[j]));
          if (last != revlist.end()) {
            revlist.resize(std::distance(revlist.begin(), last));
          }
          assert(!revlist.empty());

          std::sort(revlist.begin(), revlist.end(), OverlapCmp());

          bool largest = revlist[0]->end() == vertex;
          for (size_t k = 1; k < revlist.size() && !largest; ++k) {
            if (revlist[0]->coord().length() - revlist[k]->coord().length() < _delta) {
              largest = revlist[k]->end() == vertex; 
            }
          }
          if (largest) {
            continue;
          }
        } else if (fwdlist[0]->isSelf()) {
          continue;
        }
      }

      if (dir == Edge::ED_SENSE) {
        LOG4CXX_INFO(logger, boost::format("[MaximumOverlapVisitor] remove edge %s->%s (%d)") % fwdlist[j]->start()->id() % fwdlist[j]->end()->id() % fwdlist[j]->coord().length());
      } else {
        LOG4CXX_INFO(logger, boost::format("[MaximumOverlapVisitor] remove edge %s->%s (%d)") % fwdlist[j]->end()->id() % fwdlist[j]->start()->id() % fwdlist[j]->coord().length());
      }

      fwdlist[j]->color(GC_BLACK);
      fwdlist[j]->twin()->color(GC_BLACK);

      ++_dummys;
      modified = true;
    }
  }

  return modified;
}

void MaximumOverlapVisitor::postvisit(Bigraph* graph) {
  graph->sweepEdges(GC_BLACK);
  LOG4CXX_INFO(logger, boost::format("[MaximumOverlapVisitor] Removed %d dummy edges") % _dummys);
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

        std::unordered_map<Vertex::Id, int> distancelist = {{vertex->id(), 0}};
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
                const EdgePtrList& straight = this->edges(p, searchDir); // forward & backword edges

                // Don't walk thru singular self edges
                if (straight.size() != 1 || straight[0]->isSelf() || straight[0]->end()->color() == GC_RED) {
                    break;
                }
                Edge* single = straight[0];
                Edge* twin = single->twin();
                Vertex* end = single->end();
                const EdgePtrList& opposite = this->edges(end, twin->dir());
                if (opposite.size() != 1) {
                    break;
                }

                p = end;
                if (searchDir == Edge::ED_SENSE) {      // forward
                    const SeqCoord& coord = single->coord();
                    distance += flag * (coord.seqlen - coord.length());
                    //distance += flag * single->coord().complement().length();
                } else {                                // backward
                    const SeqCoord& coord = twin->coord();
                    distance += flag * (coord.seqlen - coord.length());
                    //distance += flag * twin->coord().complement().length();
                }
                if (single->comp() == Edge::EC_REVERSE) {
                    searchDir = Edge::EDGE_DIRECTIONS[Edge::EC_COUNT - searchDir - 1];
                }
                distancelist[p->id()] = distance;
                LOG4CXX_DEBUG(logger, boost::format("InsertSizeEstimateVisitor::visit addVertex(%s,%s,%d)") % vertex->id() % p->id() % distance);

                p->color(GC_RED);
            }
        }

        for (auto i = distancelist.begin(); i != distancelist.end(); ++i) {
            Vertex::Id  pairId = PairEnd::id(i->first);
            if (i->first < pairId) {
                auto j = distancelist.find(pairId);
                if (j != distancelist.end()) {
                    size_t distance = std::abs(j->second - i->second);
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

    typedef boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::count, boost::accumulators::tag::mean, boost::accumulators::tag::moment<2> > > Accumulator;
    Accumulator acc;
    std::for_each(_samples.begin(), _samples.end(), std::ref(acc));
    if (boost::accumulators::count(acc) > 0) {
        _average = (size_t)boost::accumulators::mean(acc);
        _delta = std::sqrt(
                boost::accumulators::moment<2>(acc) - std::pow(boost::accumulators::mean(acc), 2)
                );

        LOG4CXX_INFO(logger, boost::format("InsertSizeEstimateVisitor::average=%d, delta=%d") % _average % _delta);
    }

    graph->color(GC_GREEN);
}

EdgePtrList InsertSizeEstimateVisitor::edges(const Vertex* vertex, Edge::Dir dir) { // reduced edges
    EdgePtrList edges = vertex->edges(dir);
    std::sort(edges.begin(), edges.end(), OverlapCmp());
    // remove if the edge links to the same seq with `vertex` or it's subseq
    {
        size_t k = 0;
        for (size_t i = 0; i < edges.size(); ++i) {
            if (!edges[i]->coord().isContained() || !edges[i]->coord().isExtreme()) { // containment
                edges[k++] = edges[i];
            }
        }
        edges.resize(k);
    }
    // remove the other edges if they links to the same seq
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
// PairedReadVisitor
//

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

typedef std::unordered_map<Vertex::Id, BigraphWalk::DistanceAttr> PairedDistanceMap;
typedef std::unordered_map<Vertex::Id, PairedDistanceMap> PairedLinkList;

class PairedVertexProcess {
public:
    PairedVertexProcess(Bigraph* graph, PairedReadVisitor* vistor) : _graph(graph), _visitor(vistor) {
    }
    BigraphWalk::NodePtrList operator()(int tid, const Vertex* vertex1) {
        return process(vertex1);
    }
    BigraphWalk::NodePtrList process(const Vertex* vertex1) {
        BigraphWalk::NodePtrList linklist;

        const Vertex* paired_v1 = _graph->getVertex(PairEnd::id(vertex1->id()));
        assert(paired_v1 != NULL);

        if (vertex1->id() == "ST-K00126:7:H5W53BBXX:553782:575696:10120:10533/1") {
            assert(true);
        }
        // BFS the adjacents of vertex1
        BigraphWalk::NodePtrList adjacents;
        if (vertex1->seq().length() > _visitor->_maxDistance) {
            std::function<bool(const Edge* edge)> filter = [](const Edge* edge) -> bool {
                    if (edge->dir() == Edge::ED_SENSE || edge->comp() == Edge::EC_REVERSE) {
                        const Edge* e = NULL;
                        if (edge->dir() == Edge::ED_SENSE) {
                            e = edge;
                        } else {
                            e = edge->twin();
                        }
                        const SeqCoord& coord = e->coord();
                        return coord.seqlen > coord.length();
                    }
                    return false;
                };
            BigraphWalk::build(vertex1, filter, NULL, 0, _visitor->_maxDistance, _visitor->_maxNodes, &adjacents);
        }
        std::sort(adjacents.begin(), adjacents.end(), [](const BigraphWalk::NodePtr& x, const BigraphWalk::NodePtr& y) -> bool {
                    return std::abs(x->attr.distance) < std::abs(y->attr.distance);
                }); // sort by distance

        // Match each virtual read with paired vertex1
        for (const auto& node1 : adjacents) {
            const Vertex* paired_v2 = _graph->getVertex(PairEnd::id(node1->vertex->id()));
            assert(paired_v2 != NULL);
            LOG4CXX_DEBUG(logger, boost::format("vertex1: %s<->%s, vertex2: %s<->%s, attr.distance=%d") % vertex1->id() % paired_v1->id() % node1->vertex->id() % paired_v2->id() % node1->attr.distance);

            BigraphWalk::NodePtrList faraways;
            for (size_t i = 0; i < Edge::ED_COUNT && faraways.empty(); ++i) {
                Edge::Dir dir = Edge::EDGE_DIRECTIONS[i];
                BigraphWalk::build(paired_v1, [dir](const Edge* edge)->bool {
                            return edge->dir() == dir;
                        }, paired_v2, 0, std::abs(node1->attr.distance) + _visitor->_insertDelta*4, 1, &faraways);
            }
            for (const auto& node2 : faraways) {
                linklist.push_back(node1);
                LOG4CXX_DEBUG(logger, boost::format("paired_read_all\t%s\t%s\t%d\t%s\t%s\t%d") % vertex1->id() % node1->vertex->id() % node1->attr.distance % paired_v1->id() % node2->vertex->id() % node2->attr.distance);
                //break;
            }
        }

        return linklist;
    }
private:
    Bigraph* _graph;
    PairedReadVisitor* _visitor;
};

class PairedVertexPostProcess {
public:
    PairedVertexPostProcess(PairedLinkList* links) : _links(links) {
    }
    void operator()(std::vector<const Vertex*>::iterator it, BigraphWalk::NodePtrList& linklist) {
        process(*it, linklist);
    }
    void process(const Vertex* vertex1, BigraphWalk::NodePtrList& linklist) {
        std::sort(linklist.begin(), linklist.end(), [](const BigraphWalk::NodePtr& x, const BigraphWalk::NodePtr& y) -> bool {
                    return std::abs(x->attr.distance) < std::abs(y->attr.distance);
                });
        for (size_t i = 0; i < linklist.size(); ++i) {
            const BigraphWalk::NodePtr& xi = linklist[i];

            BigraphWalk::DistanceAttr e = BigraphWalk::attrLink(xi->attr);
            addLink(vertex1->id(), xi->vertex->id(), e, _links);

            for (size_t j = i + 1; j < linklist.size(); ++j) {
                const BigraphWalk::NodePtr& xj = linklist[j];
                if (BigraphWalk::diffDir(xi->attr, xj->attr) || xi->attr.distance == xj->attr.distance) { // different search dir or equal
                    continue;
                }

                assert(std::abs(xi->attr.distance) < std::abs(xj->attr.distance));
                BigraphWalk::DistanceAttr e = BigraphWalk::attrLink(xi->attr, xj->attr);
                if (BigraphWalk::hasLink(xi->vertex, xj->vertex, e)) {
                    addLink(xi->vertex->id(), xj->vertex->id(), e, _links);
                } else {
                    LOG4CXX_DEBUG(logger, boost::format("paired_read_branch\t%s\t%s\t%d\t%s\t%d") % vertex1->id() % xi->vertex->id() % xi->attr.distance % xj->vertex->id() % xj->attr.distance);
                }
            }
        }
    }
private:
    void addLink(const Vertex::Id& v1, const Vertex::Id& v2, const BigraphWalk::DistanceAttr& e, PairedLinkList* links) {
        //LOG4CXX_TRACE(logger, boost::format("PairedLinkList::addLink\t%s\t%s\t%d\t%d\t%d") % v1 % v2 % e.distance % e.dir % e.comp);
        if (e.distance < 0) {
            BigraphWalk::DistanceAttr t = e.twin();
            t.distance = -t.distance;
            addLink(v2, v1, t, links);
        } else {
            auto i = links->find(v1);
            if (i == links->end()) {
                PairedDistanceMap tbl = {{v2, e}}; 
                (*links)[v1] = tbl;
            } else {
                auto k = i->second.find(v2);
                if (k != i->second.end() && k->second.distance != e.distance) {
                    LOG4CXX_WARN(logger, boost::format("PairedReadVisitor::addLink(%s, %s, %d) != %d") % v1 % v2 % e.distance % k->second.distance);
                }
                if (k == i->second.end() || k->second.distance > e.distance) {
                    i->second[v2] = e;
                }
            }
        }
    }

    PairedLinkList* _links;
};

class PairedEdgeCreator {
public:
    PairedEdgeCreator(Bigraph* graph) : _graph(graph) {
    }
    void create(const Vertex::Id& v1, const Vertex::Id& v2, const BigraphWalk::DistanceAttr& attr, GraphColor color) {
        assert(attr.distance > 0);

        Vertex* vertex[2] = {
            _graph->getVertex(v1), 
            _graph->getVertex(v2), 
        };

        bool found = false;
        const EdgePtrList& edges = vertex[0]->edges();
        for (auto& edge : edges) {
            //if (edge->dir() == attr.dir && edge->comp() == attr.comp && edge->end() == vertex[1]) {
            if (edge->dir() == attr.dir && edge->end() == vertex[1]) {
                if (edge->comp() == attr.comp && edge->coord().complement().length() == attr.distance) {
                    edge->color(color);
                    edge->twin()->color(color);
                    found = true;
                    break;
                }
                LOG4CXX_WARN(logger, boost::format("PairedReadVisitor::addEdge(%s, %s, %d) != %d") % v1 % v2 % attr.distance % edge->coord().interval.start);
            }
        }
        if (!found) {
            // sequnces
            const std::string& seq1 = vertex[0]->seq();
            const std::string& seq2 = vertex[1]->seq();

            // coordinates
            SeqCoord coord[2] = {
                SeqCoord(attr.distance, seq1.length() - 1, seq1.length()),
                SeqCoord(0, seq1.length() - attr.distance - 1, seq2.length())
            };
            if (attr.dir == Edge::ED_ANTISENSE) {
                coord[0] = SeqCoord(0, seq2.length() - attr.distance - 1, seq1.length());
                coord[1] = SeqCoord(attr.distance, seq2.length() - 1, seq2.length());
            }
            if (attr.comp == Edge::EC_REVERSE) {
                coord[1].flip();
            }

            EdgeCreator creator(_graph, true, -1);
            Overlap ovr(v1, coord[0], v2, coord[1], attr.comp == Edge::EC_REVERSE, 0);
            creator.create(ovr, color);
        }
    }
private:
    Bigraph* _graph;
};

struct PairedEdgeFilter {
    PairedEdgeFilter(GraphColor c) : _color(c), _vertex(NULL) {
        for (size_t i = 0; i < Edge::ED_COUNT; ++i) {
            _hasColor[i] = false;
        }
    }
    bool operator()(const Vertex* vertex, const Edge* edge) {
        //return edge->color() != _color;
        assert(vertex != NULL);
        if (vertex != _vertex) {
            _vertex = vertex;
            for (size_t i = 0; i < Edge::ED_COUNT; ++i) {
                _hasColor[i] = false;
            }
            const EdgePtrList& edges = vertex->edges();
            for (const auto& e : edges) {
                if (e->color() == _color) {
                    _hasColor[e->dir()] = true;
                }
            }
        }

        return (_hasColor[edge->dir()] && edge->color() != _color) || edge->coord().isFull();
    }
private:
    GraphColor _color;
    const Vertex* _vertex;
    bool _hasColor[Edge::ED_COUNT];
};

void PairedReadVisitor::postvisit(Bigraph* graph) {
    PairedLinkList links;

    {
        size_t threads = parallel::threads(_threads);
        std::vector<PairedVertexProcess *> proclist(_threads);
        for (size_t i = 0; i < proclist.size(); ++i) {
            proclist[i] = new PairedVertexProcess(graph, this);
        }
        PairedVertexPostProcess postproc(&links);

        parallel::foreach<BigraphWalk::NodePtrList>(_vertices.begin(), _vertices.end(),
            [&](int tid, std::vector<const Vertex*>::iterator it) {
                assert(tid < threads);
                return proclist[tid]->process(*it);
            }, postproc, threads);

        for (size_t i = 0; i < proclist.size(); ++i) {
            delete proclist[i];
        }
    }

    ///////////////////////////////////////////
    // simplify
    //
    // NOTE that all the distance>=0
    //
    // GC_GRAY  => unkown
    // GC_WHITE => true positives
    // GC_BLACK => false positives
    //
    //////////////////////////////////////////
    {
        EdgeColorVisitor ecVist(GC_GRAY);
        graph->visit(&ecVist);
    }

    // create true positive edges.
    PairedEdgeCreator creator(graph);
    for (auto i = links.begin(); i != links.end(); ++i) {
        std::vector<std::pair<Vertex::Id, BigraphWalk::DistanceAttr> > nodelist;

        // sorted by distance
        std::copy(i->second.begin(), i->second.end(), std::back_inserter(nodelist));
        std::sort(nodelist.begin(), nodelist.end(), [](const std::pair<Vertex::Id, BigraphWalk::DistanceAttr>& x, const std::pair<Vertex::Id, BigraphWalk::DistanceAttr>& y) -> bool {
                    return x.second.distance < y.second.distance;
                });

        for (size_t j = 0; j < nodelist.size(); ++j) {
            const std::pair<Vertex::Id, BigraphWalk::DistanceAttr>& xj = nodelist[j];

            bool hasLink = false;
            for (size_t k = 0; k < j; ++k) {
                const std::pair<Vertex::Id, BigraphWalk::DistanceAttr>& xk = nodelist[k];
                //fixme: check correctness
                if (xk.second.dir == xj.second.dir && BigraphWalk::hasLink(graph->getVertex(xk.first), xk.second, graph->getVertex(xj.first), xj.second)) {
                    hasLink = true;
                    break;
                }
            }
            if (!hasLink) {
                creator.create(i->first, xj.first, xj.second, GC_WHITE);
                LOG4CXX_INFO(logger, boost::format("paired_read_simplify\t%s\t%s\t%d") % i->first % xj.first % xj.second.distance);
            }
        }
    }

    // cleanup false positive edges.
    {
        PairedEdgeFilter peFilter(GC_WHITE);
        EdgeColorVisitor ecVisit(GC_BLACK, peFilter, true);
        graph->visit(&ecVisit);
    }

    graph->sweepEdges(GC_BLACK);
}

//
// LinkedReadVisitor
//

void LinkedReadVisitor::previsit(Bigraph* graph) {
    EdgeColorVisitor ecVisit(GC_WHITE);
    graph->visit(&ecVisit);

    _dummys = 0;
}

bool LinkedReadVisitor::visit(Bigraph* graph, Vertex* vertex) {
    if (vertex->seq().length() < _minLength || vertex->coverage() < _minCoverage) {
        return false;
    }
    bool modified = false;
    const auto& indexTbl1 = vertex->indexTbl();
    auto check = [&](Edge::Dir dir) {
        for (auto& edge : vertex->edges(dir)) {
            for (const auto& index : edge->end()->indexTbl()) {
                auto it = indexTbl1.equal_range(index.first);
            }
        }
    };
    const EdgePtrList& edges = vertex->edges();
    for (auto& edge : edges) {
        const Vertex::IndexTable& indexTbl2 = edge->end()->indexTbl();

        size_t fragment = 0;
        for (const auto& index : indexTbl2) {
            if (indexTbl1.find(index.first) != indexTbl1.end()) {
                ++fragment;
            }
        }
        if (fragment <= 1) {
            edge->color(GC_BLACK);
            edge->twin()->color(GC_BLACK);
            LOG4CXX_INFO(logger, boost::format("LinkedReadVisitor::visit(%s, %s)") % vertex->id() % edge->end()->id());
            ++_dummys;
        }
    }
    return modified;
}

void LinkedReadVisitor::postvisit(Bigraph* graph) {
    graph->sweepEdges(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[LinkedReadVisitor] Removed %d dummy edges") % _dummys);
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
// TransitiveReductionVisitor
//
void TransitiveReductionVisitor::previsit(Bigraph* graph) {
}

bool TransitiveReductionVisitor::visit(Bigraph* graph, Vertex* vertex) {
    bool modified = false;
    return modified;
}

void TransitiveReductionVisitor::postvisit(Bigraph* graph) {
}

namespace SeqAlign {
  double align(const std::string& x, const std::string& y, double match = 0.0, double mismatch = 1.0, double gap = 2.0) {
    std::vector<std::vector<double> > scores(x.length() + 1);
    for (size_t i = 0; i <= x.length(); ++i) {
      scores[i].resize(y.length() + 1);
      scores[i][0] = i*gap;
    }
    for (size_t j = 0; j <= y.length(); ++j) {
      scores[0][j] = j*gap;
    }
    for (size_t i = 0; i < x.length(); ++i) {
      for (size_t j = 0; j < y.length(); ++j) {
        scores[i+1][j+1] = scores[i][j] + (x[i] == y[j] ? match : mismatch);
        scores[i+1][j+1] = std::min(scores[i+1][j+1], scores[i+1][j] + gap);
        scores[i+1][j+1] = std::min(scores[i+1][j+1], scores[i][j+1] + gap);
      }
    }
    return scores[x.length()][y.length()];
  }
};  // SeqAlign

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

    auto similarity = [&](const std::string& x, const std::string& y) {
        return 0.0;
      };

    const std::string& seq = vertex->seq();
    if (vertex->degrees() == 0) {
        // Is an island, remove if the sequence length is less than the threshold
        if (seq.length() <= _minLength && Point::avg(vertex) <= Point::avg(_minCoverage, _minLength)) {
            LOG4CXX_TRACE(logger, boost::format("[TrimVisitor] island %s length=%d coverage=%d") % vertex->id() % seq.length() % vertex->coverage());
            vertex->color(GC_BLACK);
            ++_island;
            modified = true;
        }
    } else {
        // Check if this node is a dead-end
        for (size_t idx = 0; idx < Edge::ED_COUNT; idx++) {
            Edge::Dir dir = Edge::EDGE_DIRECTIONS[idx];
            if (vertex->degrees(dir) == 0 && seq.length() <= _minLength && Point::avg(vertex) <= Point::avg(_minCoverage, _minLength)) {
                LOG4CXX_TRACE(logger, boost::format("[TrimVisitor] terminal %s length=%d coverage=%d") % vertex->id() % seq.length() % vertex->coverage());
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

namespace AIFeat {

class HiFiParser {
 public:
  HiFiParser(const Vertex* start, const Vertex* end) : _start(start), _end(end) {
  }

  bool operator()(const std::map<std::string, double>* feats) const {
    size_t cnt[kFields] = {0};

    const auto& x = _start->indexTbl();
    const auto& y = _end->indexTbl();
    auto i = x.begin(), j = y.begin();
    while (i != x.end() && j != y.end()) {
      if (i->first < j->first) {
        ++i;
        ++cnt[kLinkAllX];
      } else if (i->first > j->first) {
        ++j;
        ++cnt[kLinkAllY];
      } else {
        Vertex::IndexTable::key_type key = i->first;
        Vertex::IndexTable::mapped_type vmin = i->second, vmax = i->second;
        while (i != x.end() && i->first == key) {
          ++i;
          ++cnt[kLinkAllX];
          ++cnt[kLinkCommonX];
          vmin = std::min(vmin, i->second);
          vmax = std::max(vmax, i->second);
        }
        while (j != y.end() && j->first == key) {
          ++j;
          ++cnt[kLinkAllY];
          ++cnt[kLinkCommonY];
          if (j->second < vmin) {
            ++cnt[kLinkMin];
          } else if (vmax < j->second) {
            ++cnt[kLinkMax];
          } else {
            ++cnt[kLinkMid];
          }
        }
      }
    }
    while (i != x.end()) {
      ++i;
      ++cnt[kLinkAllX];
    }
    while (j != y.end()) {
      ++j;
      ++cnt[kLinkAllY];
    }

    // (*feats)[] = cnt[kLinkAllX];
    return true;
  }
 private:
  enum {
    kLinkAllX = 0, 
    kLinkAllY,
    kLinkCommonX, 
    kLinkCommonY,
    kLinkMin, 
    kLinkMax,
    kLinkMid,
    kFields
  };
  const Vertex* _start;
  const Vertex* _end;
};

bool extract(const Edge* edge) {
  return true;
};

};  // namespace AIFeat

namespace HiFiParser {
 
enum {
  kLinkAllX = 0, 
  kLinkAllY,
  kLinkCommonX, 
  kLinkCommonY,
  kLinkMin, 
  kLinkMax,
  kLinkMid,
  kFields
};

void parse(const Vertex* start, const Vertex* end, size_t* cnt) {
  const auto& x = start->indexTbl();
  const auto& y = end->indexTbl();
  auto i = x.begin(), j = y.begin();
  while (i != x.end() && j != y.end()) {
    if (i->first < j->first) {
      ++i;
      ++cnt[kLinkAllX];
    } else if (i->first > j->first) {
      ++j;
      ++cnt[kLinkAllY];
    } else {
      Vertex::IndexTable::key_type key = i->first;
      Vertex::IndexTable::mapped_type vmin = i->second, vmax = i->second;
      while (i != x.end() && i->first == key) {
        ++i;
        ++cnt[kLinkAllX];
        ++cnt[kLinkCommonX];
        vmin = std::min(vmin, i->second);
        vmax = std::max(vmax, i->second);
      }
      while (j != y.end() && j->first == key) {
        ++j;
        ++cnt[kLinkAllY];
        ++cnt[kLinkCommonY];
        if (j->second < vmin) {
          ++cnt[kLinkMin];
        } else if (vmax < j->second) {
          ++cnt[kLinkMax];
        } else {
          ++cnt[kLinkMid];
        }
      }
    }
  }
  while (i != x.end()) {
    ++i;
    ++cnt[kLinkAllX];
  }
  while (j != y.end()) {
    ++j;
    ++cnt[kLinkAllY];
  }
}

void parse(const Edge* edge, size_t* cnt) {
  parse(edge->start(), edge->end(), cnt);
}

void parse(const Vertex* vertex, Edge::Dir dir, size_t* cnt) {
  auto sequenced = [&](const size_t* cnt) {
      return cnt[kLinkMin] + cnt[kLinkMid] + cnt[kLinkMax];
    };
  for (const auto& edge : vertex->edges(dir)) {
    size_t tmp[kFields] = {0};
    parse(edge, tmp);
    if (sequenced(tmp) >= sequenced(cnt)) {
      std::copy(tmp, tmp + kFields, cnt);
    }
  }
}

double linkr(size_t c, size_t x, size_t y) {
  if (x == 0 || y == 0) {
    return 0.0;
  }
  return (double)c/std::min(x, y);
}

};  // namespace HiFiParser

#ifdef HAVE_MLPACK
//
// AIVisitor
//
void AIVisitor::previsit(Bigraph* graph) {
  EdgeColorVisitor ecVist(GC_GRAY, true);
  graph->visit(&ecVist);

  _blacks = 0;
  _whites = 0;
  _grays = 0;
}

bool AIVisitor::visit(Bigraph* graph, Vertex* vertex) {
  bool modified = false;

  auto m = (BaggingModel<mlpack::tree::DecisionTree<> > *)_model;
  auto edges = vertex->edges(Edge::ED_SENSE);
  std::sort(edges.begin(), edges.end(), OverlapCmp());
  for (size_t i = 0; i < edges.size(); ++i) {
    Edge* edge = edges[i];
    const Vertex* end = edge->end();

    size_t j = 0, k = 0;
    {
      for (auto e : end->edges(Edge::ED_ANTISENSE)) {
        if (e != edge->twin()) {
          if (e->coord().length() >= edge->coord().length()) {
            ++j;
          }
          if (e->start()->seq().length() >= vertex->seq().length()) {
            ++k;
          }
        }
      }
    }

    size_t v2xcnt[HiFiParser::kFields] = {0}, x2ycnt[HiFiParser::kFields] = {0}, y2wcnt[HiFiParser::kFields] = {0};
    HiFiParser::parse(vertex, Edge::ED_ANTISENSE, v2xcnt);
    HiFiParser::parse(vertex, end, x2ycnt);
    HiFiParser::parse(end, Edge::ED_SENSE, y2wcnt);

    std::vector<double> vec = {
        (double)vertex->seq().length(),  // lenx
        (double)vertex->coverage(),      // coveragex
        Repeatness::calc(vertex, _n, _G),// repeatx
        (double)vertex->degrees(Edge::ED_ANTISENSE),  // indegreex
        (double)vertex->degrees(Edge::ED_SENSE),  // outdegreex
        (double)i,                       // orankx
        (double)end->seq().length(),     // leny
        (double)end->coverage(),         // coveragey
        Repeatness::calc(end, _n, _G),   // repeaty
        (double)end->degrees(Edge::ED_ANTISENSE),  // indegreey
        (double)end->degrees(Edge::ED_SENSE),  // outdegreey
        (double)j,                       // oranky
        (double)v2xcnt[HiFiParser::kLinkAllX], // alinkx
        (double)v2xcnt[HiFiParser::kLinkAllY], // alinky
        (double)v2xcnt[HiFiParser::kLinkCommonX], // clinkx
        (double)v2xcnt[HiFiParser::kLinkCommonY], // clinky
        (double)v2xcnt[HiFiParser::kLinkMin], // minlink
        (double)v2xcnt[HiFiParser::kLinkMax], // maxlink
        (double)v2xcnt[HiFiParser::kLinkMid], // midlink
        HiFiParser::linkr(v2xcnt[HiFiParser::kLinkMin], v2xcnt[HiFiParser::kLinkCommonX], v2xcnt[HiFiParser::kLinkCommonY]), 
        HiFiParser::linkr(v2xcnt[HiFiParser::kLinkMax], v2xcnt[HiFiParser::kLinkCommonX], v2xcnt[HiFiParser::kLinkCommonY]), 
        HiFiParser::linkr(v2xcnt[HiFiParser::kLinkMid], v2xcnt[HiFiParser::kLinkCommonX], v2xcnt[HiFiParser::kLinkCommonY]), 
        (double)x2ycnt[HiFiParser::kLinkAllX], // alinkx
        (double)x2ycnt[HiFiParser::kLinkAllY], // alinky
        (double)x2ycnt[HiFiParser::kLinkCommonX], // clinkx
        (double)x2ycnt[HiFiParser::kLinkCommonY], // clinky
        (double)x2ycnt[HiFiParser::kLinkMin], // minlink
        (double)x2ycnt[HiFiParser::kLinkMax], // maxlink
        (double)x2ycnt[HiFiParser::kLinkMid], // midlink
        HiFiParser::linkr(x2ycnt[HiFiParser::kLinkMin], x2ycnt[HiFiParser::kLinkCommonX], x2ycnt[HiFiParser::kLinkCommonY]), 
        HiFiParser::linkr(x2ycnt[HiFiParser::kLinkMax], x2ycnt[HiFiParser::kLinkCommonX], x2ycnt[HiFiParser::kLinkCommonY]), 
        HiFiParser::linkr(x2ycnt[HiFiParser::kLinkMid], x2ycnt[HiFiParser::kLinkCommonX], x2ycnt[HiFiParser::kLinkCommonY]), 
        (double)y2wcnt[HiFiParser::kLinkAllX], // alinkx
        (double)y2wcnt[HiFiParser::kLinkAllY], // alinky
        (double)y2wcnt[HiFiParser::kLinkCommonX], // clinkx
        (double)y2wcnt[HiFiParser::kLinkCommonY], // clinky
        (double)y2wcnt[HiFiParser::kLinkMin], // minlink
        (double)y2wcnt[HiFiParser::kLinkMax], // maxlink
        (double)y2wcnt[HiFiParser::kLinkMid], // midlink
        HiFiParser::linkr(y2wcnt[HiFiParser::kLinkMin], y2wcnt[HiFiParser::kLinkCommonX], y2wcnt[HiFiParser::kLinkCommonY]), 
        HiFiParser::linkr(y2wcnt[HiFiParser::kLinkMax], y2wcnt[HiFiParser::kLinkCommonX], y2wcnt[HiFiParser::kLinkCommonY]), 
        HiFiParser::linkr(y2wcnt[HiFiParser::kLinkMid], y2wcnt[HiFiParser::kLinkCommonX], y2wcnt[HiFiParser::kLinkCommonY]), 
        (double)k,                       // sranky
        (double)edge->coord().length(),  // overlap between x and y
        Point::avg(vertex),              // pointx
        Point::avg(end)                  // pointy
      };
    size_t l = m->Classify(vec);
    // LOG4CXX_INFO(logger, boost::format("[AIVisitor] start=%s end=%s label=%lu") % vertex->id() % end->id() % l);
    if (l < 1) {
      ++_blacks;
      edge->color(GC_BLACK);
      edge->twin()->color(GC_BLACK);
      modified = true;
    } else if (l > m->Size()/2) {
      ++_whites;
      edge->color(GC_WHITE);
      edge->twin()->color(GC_WHITE);
      // modified = true;
    } else {
      ++_grays;
    }
  }

  return modified;
}

void AIVisitor::postvisit(Bigraph* graph) {
  graph->sweepEdges(GC_BLACK);
  LOG4CXX_INFO(logger, boost::format("[AIVisitor]: Removed %d edges: whites=%lu, gray=%lu") % _blacks % _whites % _grays);
}
#endif  // HAVE_MLPACK

//
// UnitigVisitor
//
void UnitigVisitor::previsit(Bigraph* graph) {
  _unitigs = 0;
}

bool UnitigVisitor::visit(Bigraph* graph, Vertex* vertex) {
  bool modified = false;
  if (Repeatness::calc(vertex, _n, _G) >= _T) {
    for (size_t i = 0; i < Edge::ED_COUNT; ++i) {
      auto edges = vertex->edges(Edge::EDGE_DIRECTIONS[i]);
      if (edges.size() == 1 and Repeatness::calc(edges[0]->end(), _n, _G) < Repeatness::calc(vertex, _n, _G)) {
        Vertex* end = edges[0]->end();
        if (end->degrees(Edge::ED_SENSE) <= 1 and end->degrees(Edge::ED_ANTISENSE) <= 1) {
          continue;
        }

        size_t x2ycnt[HiFiParser::kFields] = {0};
        HiFiParser::parse(vertex, end, x2ycnt);

        Vertex tmp(end->id() + "_copy", end->seq(), end->contained(), end->index(), end->coverage(), end->extension());
        for (auto edge : end->edges(Edge::EDGE_DIRECTIONS[i])) {
          Vertex* verts[2] = {&tmp, edge->end()};
          Edge* arr[2] = {edge, edge->twin()};
          Edge* vec[2];
          for (size_t j = 0; j < 2; ++j) {
            vec[j] = new Edge(verts[1 - j], arr[j]->dir(), arr[j]->comp(), arr[j]->coord());
            vec[j]->color(arr[j]->color());
          }
          vec[0]->twin(vec[1]);
          vec[1]->twin(vec[0]);

          graph->addEdge(verts[0], vec[0]);
          graph->addEdge(verts[1], vec[1]);
        }

        Vertex* verts[2] = {vertex, &tmp};
        Edge* arr[2] = {edges[0], edges[0]->twin()};
        Edge* vec[2];
        for (size_t j = 0; j < 2; ++j) {
          vec[j] = new Edge(verts[1 - j], arr[j]->dir(), arr[j]->comp(), arr[j]->coord());
          vec[j]->color(arr[j]->color());
        }
        vec[0]->twin(vec[1]);
        vec[1]->twin(vec[0]);

        graph->addEdge(verts[0], vec[0]);
        graph->addEdge(verts[1], vec[1]);

        vertex->removeEdge(edges[0]);
        Edge* twin = edges[0]->twin();
        end->removeEdge(twin);
        SAFE_DELETE(edges[0]);
        SAFE_DELETE(twin);

        LOG4CXX_INFO(logger, boost::format("[UnitigVisitor] vertex: %s, length: %lu, end: %s, length: %lu, dir: %d, degree: %lu alinkx:%lu alinky:%lu clinkx:%lu clinky:%lu minlink:%lu maxlink:%lu midlink:%lu") % vertex->id() % vertex->seq().length() % end->id() % end->seq().length() % i % end->degrees(Edge::EDGE_DIRECTIONS[i])
                % x2ycnt[HiFiParser::kLinkAllX] % x2ycnt[HiFiParser::kLinkAllY] % x2ycnt[HiFiParser::kLinkCommonX] % x2ycnt[HiFiParser::kLinkCommonY] % x2ycnt[HiFiParser::kLinkMin] % x2ycnt[HiFiParser::kLinkMax] % x2ycnt[HiFiParser::kLinkMid]);

        assert(vertex->degrees(Edge::EDGE_DIRECTIONS[i]) == 1);
        graph->merge(vertex, vec[0]);

        ++_unitigs;
        modified = true;
      }
    }
  }
  return modified;
}

void UnitigVisitor::postvisit(Bigraph* graph) {
  LOG4CXX_INFO(logger, boost::format("[UnitigVisitor]: Walked %lu unitigs") % _unitigs);
}

//
// GANVisitor
//
void GANVisitor::previsit(Bigraph* graph) {
  EdgeColorVisitor ecVist(GC_GRAY, true);
  graph->visit(&ecVist);

  _blacks = 0;
  _whites = 0;
  _grays = 0;
}

bool GANVisitor::visit(Bigraph* graph, Vertex* vertex) {
  bool modified = false;

  static boost::regex cigar("^(\\d*)=$");
  typedef std::tuple<std::string, std::string, size_t, size_t> Alignment;
  auto AlignmentParser = [&](const std::string& text, Alignment* align) {
      std::vector<std::string> vec;
      boost::algorithm::split(vec, text, boost::is_any_of("|"));
      if (vec.size() < 4) {
        return false;
      }
      *align = std::make_tuple(vec[0], vec[1], boost::lexical_cast<size_t>(vec[2]), boost::lexical_cast<size_t>(vec[3]));
      return true;
    };
  auto AlignmentSorter = [&](const Alignment& x, const Alignment& y, size_t offsetx = 0, size_t offsety = 0) {
      if (std::get<1>(x) == std::get<1>(y)) {
        return std::get<2>(x) + offsetx < std::get<2>(y) + offsety;
      }
      return std::get<1>(x) < std::get<1>(y);
    };
  auto AlignmentMatcher = [&](const Edge* edge, const std::vector<Alignment>& curr, const std::vector<Alignment>& next, 
        std::vector<std::pair<Alignment, Alignment> >* mathes) {
      size_t i = 0, j = 0, x = vertex->seq().length(), y = edge->coord().length();
      while (i < curr.size() && j < next.size()) {
        if (AlignmentSorter(curr[i], next[j], x, y)) {
          ++i;
        } else if (AlignmentSorter(next[j], curr[i], y, x)) {
          ++j;
        } else {
          mathes->push_back(std::make_pair(curr[i], next[j]));
          ++i;
          ++j;
        }
      }
    };
  auto ExtParser = [&](const std::string& text, std::vector<Alignment>* alignments) {
      std::vector<std::string> vec;
      boost::algorithm::split(vec, text, boost::is_any_of(","));
      for (const auto& item : vec) {
        if (!item.empty()) {
          Alignment alignment;
          if (!AlignmentParser(item, &alignment)) {
            return false;
          }
          alignments->push_back(alignment);
        }
      }
      return true;
    };
  auto LabelParser = [&](std::pair<Alignment, Alignment>& m) {
      return boost::regex_match(std::get<0>(m.first), cigar) && boost::regex_match(std::get<0>(m.second), cigar);
    };

  auto RefMatcher = [&](const Edge* edge, size_t maxLength = -1) {
      if (_ref) {
        const std::string& start = edge->start()->seq();
        const std::string& end = edge->end()->seq();
        size_t overlap = edge->coord().length();
        assert(overlap <= start.length() && overlap <= end.length());
        std::string seq;
        if (start.length() > std::max(maxLength, overlap)) {
          seq += start.substr(start.length() - std::max(maxLength, overlap));
        } else {
          seq += start;
        }
        if (end.length() > std::max(maxLength, overlap)) {
          seq += end.substr(overlap, std::max(maxLength, overlap) - overlap);
        } else {
          seq += end.substr(overlap);
        }
        return FMIndex::Interval::occurrences(seq, _ref) > 0;
      }
      return false;
    };

  std::vector<Alignment> curralignments;
  if (!ExtParser(vertex->extension(), &curralignments)) {
    LOG4CXX_ERROR(logger, boost::format("[GANVisitor] ExtParser this id=%s") % vertex->id());
    return false;
  }
  // assert(!curralignments.empty());
  std::sort(curralignments.begin(), curralignments.end(), AlignmentSorter);

  auto edges = vertex->edges(Edge::ED_SENSE);
  std::sort(edges.begin(), edges.end(), OverlapCmp());
  size_t i = 0;
  for (size_t p = 0; p < edges.size(); ++p) {
    Edge* edge = edges[p];
    const Vertex* end = edge->end();

    if (p > 0 && edge->coord().length() < edges[p - 1]->coord().length()) {
      ++i;
    }

    size_t j = 0, k = 0;
    {
      for (auto e : end->edges(Edge::ED_ANTISENSE)) {
        if (e != edge->twin()) {
          if (e->coord().length() > edge->coord().length()) {
            ++j;
          }
          if (e->start()->seq().length() > vertex->seq().length()) {
            ++k;
          }
        }
      }
    }

    std::vector<Alignment> nextalignments;
    if (!ExtParser(end->extension(), &nextalignments)) {
      LOG4CXX_ERROR(logger, boost::format("[GANVisitor] ExtParser next id=%s") % end->id());
      continue;
    }
    // assert(!nextalignments.empty());
    std::sort(nextalignments.begin(), nextalignments.end(), AlignmentSorter);

    std::vector<std::pair<Alignment, Alignment> > pairs;
    AlignmentMatcher(edge, curralignments, nextalignments, &pairs);
    std::string seq = vertex->seq() + end->seq().substr(edge->coord().length());

    int label = 0;
    if (std::any_of(pairs.begin(), pairs.end(), LabelParser)) {
      label = 1;
    } else if (_ref && FMIndex::Interval::occurrences(seq, _ref) > 0) {
      label = 1;
    } else if (RefMatcher(edge, 1000)) {
      //label = 1;
    }
    if (label == 0) {
      if (Utils::rand(1000) < 1001) {
        ++_blacks;
        edge->color(GC_BLACK);
        edge->twin()->color(GC_BLACK);
      } else {
        ++_grays;
      }
      modified = true;
    } else if (label == 1) {
      ++_whites;
      edge->color(GC_WHITE);
      edge->twin()->color(GC_WHITE);
    }
    size_t v2xcnt[HiFiParser::kFields] = {0}, x2ycnt[HiFiParser::kFields] = {0}, y2wcnt[HiFiParser::kFields] = {0};

    HiFiParser::parse(vertex, Edge::ED_ANTISENSE, v2xcnt);
    HiFiParser::parse(vertex, end, x2ycnt);
    HiFiParser::parse(end, Edge::ED_SENSE, y2wcnt);

    LOG4CXX_INFO(logger, boost::format("[GANVisitor] %d %s %s lenx:%lu coveragex:%lu indegreex:%lu outdegreex:%lu orankx:%lu cigarx:%s chrx:%s leny:%lu coveragey:%lu indegreey:%lu outdegreey:%lu oranky:%lu sranky:%lu cigary:%s chry:%s is_self:%d overlap:%lu seq:%s alinkvx:%lu alinvy:%lu clinkvx:%lu clinkvy:%lu minlinkv:%lu maxlinkv:%lu midlinkv:%lu alinkx:%lu alinky:%lu clinkx:%lu clinky:%lu minlink:%lu maxlink:%lu midlink:%lu alinkwx:%lu alinkwy:%lu clinkwx:%lu clinkwy:%lu minlinkw:%lu maxlinkw:%lu midlinkw:%lu %lu %lu %lu %lu %lu %s %s %s %s") % label % vertex->id() % end->id()
        % vertex->seq().length() % vertex->coverage() % vertex->degrees(Edge::ED_ANTISENSE) % vertex->degrees(Edge::ED_SENSE) % i % (!curralignments.empty() ? std::get<0>(curralignments[0]) : "-") % (!curralignments.empty() ? std::get<1>(curralignments[0]) : "0")
        % end->seq().length() % end->coverage() % end->degrees(Edge::ED_ANTISENSE) % end->degrees(Edge::ED_SENSE) % j % k % (!nextalignments.empty() ? std::get<0>(nextalignments[0]) : "-") % (!nextalignments.empty() ? std::get<1>(nextalignments[0]) : "0")
        % edge->isSelf() % edge->coord().length()% seq
        % v2xcnt[HiFiParser::kLinkAllX] % v2xcnt[HiFiParser::kLinkAllY] % v2xcnt[HiFiParser::kLinkCommonX] % v2xcnt[HiFiParser::kLinkCommonY] % v2xcnt[HiFiParser::kLinkMin] % v2xcnt[HiFiParser::kLinkMax] % v2xcnt[HiFiParser::kLinkMid]
        % x2ycnt[HiFiParser::kLinkAllX] % x2ycnt[HiFiParser::kLinkAllY] % x2ycnt[HiFiParser::kLinkCommonX] % x2ycnt[HiFiParser::kLinkCommonY] % x2ycnt[HiFiParser::kLinkMin] % x2ycnt[HiFiParser::kLinkMax] % x2ycnt[HiFiParser::kLinkMid]
        % y2wcnt[HiFiParser::kLinkAllX] % y2wcnt[HiFiParser::kLinkAllY] % y2wcnt[HiFiParser::kLinkCommonX] % y2wcnt[HiFiParser::kLinkCommonY] % y2wcnt[HiFiParser::kLinkMin] % y2wcnt[HiFiParser::kLinkMax] % y2wcnt[HiFiParser::kLinkMid]
        % (!curralignments.empty() ? std::get<3>(curralignments[0]) : 0)
        % (!nextalignments.empty() ? std::get<3>(nextalignments[0]) : 0)
        % curralignments.size()
        % nextalignments.size()
        % pairs.size()
        % vertex->extension()
        % end->extension()
        % vertex->seq()
        % end->seq()
      );
  }
  return modified;
}

void GANVisitor::postvisit(Bigraph* graph) {
  graph->sweepEdges(GC_BLACK);
  LOG4CXX_INFO(logger, boost::format("[GANVisitor2]: Removed %d edges: whites=%lu, gray=%lu") % _blacks % _whites % _grays);
}
