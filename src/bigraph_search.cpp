#include "bigraph_search.h"
#include "kseq.h"

#include <unordered_set>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

size_t BigraphWalk::build(NodePtrQueue& Q, const Vertex* end, size_t minDistance, size_t maxDistance, size_t maxNodes, NodePtrList* leaves) {
    size_t num = 0;
    std::unordered_set< NodePtr, NodePtrHash< NodePtr >, NodePtrCmp > visited;
    while (!Q.empty() && num < maxNodes && Q.size() < 5*maxDistance) {
        std::pair<NodePtr, int> curr = Q.front();
        Q.pop_front();

        if (visited.find(curr.first) != visited.end()) {
            continue;
        }
        visited.insert(curr.first);

        if (std::abs(curr.first->attr.distance) < maxDistance) {
            if (std::abs(curr.first->attr.distance) >= minDistance) {
                if (end == NULL) {
                    if (curr.first->attr.distance != 0) {
                        ++num;
                        if (leaves != NULL) {
                            leaves->push_back(curr.first);
                        }
                    }
                } else if (end->id() == curr.first->vertex->id()) {
                    ++num;
                    if (leaves != NULL) {
                        leaves->push_back(curr.first);
                    }
                    break;
                }
            }

            Edge::Dir dir = curr.first->attr.dir;
            if (curr.first->attr.comp == Edge::EC_REVERSE) {
                dir = Edge::EDGE_DIRECTIONS[Edge::EC_COUNT - dir - 1];
            }
            const EdgePtrList& edges = curr.first->vertex->edges();
            BOOST_FOREACH(const Edge* edge, edges) {
                if (edge->dir() == dir) {
                    //int distance = curr.first->attr.distance;
                    int distance = 0;
                    if (dir == Edge::ED_SENSE) {
                        const SeqCoord& coord = edge->coord();
                        distance = coord.seqlen - coord.length();
                        //distance += curr.second * edge->coord().complement().length();
                    } else {
                        const SeqCoord& coord = edge->twin()->coord();
                        distance = coord.seqlen - coord.length();
                        //distance += curr.second * edge->twin()->coord().complement().length();
                    }
                    distance = curr.first->attr.distance + curr.second * distance;
                    NodePtr child(new Node(edge->end(), distance, dir, edge->comp()));
                    Q.push_back(std::make_pair(child, curr.second));
                }
            }
        }
    }

    return num;
}

size_t BigraphWalk::build(const Vertex* start, std::function< bool(const Edge* edge) > filter, const Vertex* end, size_t minDistance, size_t maxDistance, size_t maxNodes, NodePtrList* leaves) {
    const EdgePtrList& edges = start->edges();

    BigraphWalk::NodePtrQueue Q;
    BOOST_FOREACH(const Edge* edge, edges) {
        int flag = 1, distance = 0;
        if (edge->dir() == Edge::ED_SENSE) {
            const SeqCoord& coord = edge->coord();
            distance = coord.seqlen - coord.length();
            //distance += curr.second * edge->coord().complement().length();
        } else {
            const SeqCoord& coord = edge->twin()->coord();
            distance = coord.seqlen - coord.length();
            //distance += curr.second * edge->twin()->coord().complement().length();
            flag = -1;
        }
        if (!filter || filter(edge)) {
            NodePtr child(new Node(edge->end(), flag * distance, edge->dir(), edge->comp()));
            Q.push_back(std::make_pair(child, flag));
        }
    }

    return build(Q, end, minDistance, maxDistance, maxNodes, leaves);
}

bool BigraphWalk::hasLink(const Vertex* v1, const Vertex* v2, int distance, Edge::Dir dir, Edge::Comp comp) {
    if (distance < 0) {
        if (comp == Edge::EC_SAME) {
            return hasLink(v2, v1, -distance, Edge::EDGE_DIRECTIONS[Edge::ED_COUNT - dir - 1], comp);
        } else {
            return hasLink(v2, v1, -distance, dir, comp);
        }
    }
    assert(distance >= 0);
    std::string seq1 = v1->seq(), seq2 = v2->seq();
    if (comp == Edge::EC_REVERSE) {
        make_reverse_complement_dna(seq2);
    }
    return (
                dir == Edge::ED_SENSE && distance < seq1.length() && boost::algorithm::starts_with(seq2, seq1.substr(distance))
            ) || (
                dir == Edge::ED_ANTISENSE && distance < seq2.length()  && boost::algorithm::starts_with(seq1, seq2.substr(distance))
            );
}
