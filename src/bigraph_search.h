#ifndef bigraph_search_h_
#define bigraph_search_h_

#include "bigraph.h"

#include <deque>
#include <functional>
#include <memory>
#include <string>

// A walk on the string graph is given by the starting vertex
// then a vector of edges used in the walk
class BigraphWalk {
public:
    class DistanceAttr {
    public:
        DistanceAttr() {
        }
        DistanceAttr(int distance, Edge::Dir dir, Edge::Comp comp) : distance(distance), dir(dir), comp(comp) {
        }
        DistanceAttr twin() const {
            DistanceAttr o(distance, dir, comp);
            if (comp == Edge::EC_SAME) {
                o.dir = Edge::EDGE_DIRECTIONS[Edge::ED_COUNT - dir - 1];
            }
            return o;
        }
        Edge::Dir dir;
        Edge::Comp comp;
        int distance;
    };
    class Node {
    public:
        Node(const Vertex* vertex, int distance, Edge::Dir dir, Edge::Comp comp) : vertex(vertex), attr(distance, dir, comp) {
        }
        ~Node() {
        }

        const Vertex* vertex;
        DistanceAttr attr;
    };
    typedef std::shared_ptr<Node> NodePtr;
    typedef std::vector<NodePtr> NodePtrList;

    struct NodePtrCmp {
        bool operator()(const NodePtr& x, const NodePtr& y) const {
            return x->vertex->id() == y->vertex->id() && x->attr.distance == y->attr.distance;
        }
    };

    template <class T>
    struct NodePtrHash {
        size_t operator()(const T& node) const {
            std::hash<std::string> hasher;
            return hasher(node->vertex->id());
        }
    };

    typedef std::deque<std::pair<NodePtr, int> > NodePtrQueue;

    static size_t build(NodePtrQueue& Q, const Vertex* end, size_t minDistance, size_t maxDistance, size_t maxNodes, NodePtrList* leaves);
    static size_t build(const Vertex* start, std::function<bool(const Edge* edge)> filter, const Vertex* end, size_t minDistance, size_t maxDistance, size_t maxNodes, NodePtrList* leaves);

    static bool hasLink(const Vertex* v1, const Vertex* v2, int distance, Edge::Dir dir, Edge::Comp comp);
    static bool hasLink(const Vertex* v1, const Vertex* v2, int distance) {
        assert(distance >= 0);
        if (distance > 0) {
            return hasLink(v1, v2, distance, Edge::ED_SENSE, Edge::EC_SAME) 
                || hasLink(v1, v2, distance, Edge::ED_SENSE, Edge::EC_REVERSE) 
                || hasLink(v1, v2, distance, Edge::ED_ANTISENSE, Edge::EC_SAME)
                || hasLink(v1, v2, distance, Edge::ED_ANTISENSE, Edge::EC_REVERSE)
                ;
        }
        return false;
    }

    static void attrLink(const DistanceAttr& e1, const DistanceAttr& e2, DistanceAttr* e) {
        // ----------------------
        // (f , f )=>f
        // (f , fc)=>fc
        // (f , r )=>fc
        // (f , rc)=>f
        // (fc, f )=>rc
        // (fc, fc)=>r
        // (fc, r )=>r
        // (fc, rc)=>rc
        // (r , f )=>rc
        // (r , fc)=>r
        // (r , r )=>r
        // (r , rc)=>rc
        // (rc, f )=>f
        // (rc, fc)=>fc
        // (rc, r )=>fc
        // (rc, rc)=>f
        // ----------------------
        e->distance = e2.distance - e1.distance;
        if (e1.comp == Edge::EC_SAME) {
            e->dir = e1.dir;
        } else {
            e->dir = Edge::EDGE_DIRECTIONS[Edge::ED_COUNT - e1.dir - 1];
        }
        DistanceAttr t1 = e1.twin(), t2 = e2.twin();
        if (t1.dir == t2.dir) {
            e->comp = Edge::EC_SAME;
        } else {
            e->comp = Edge::EC_REVERSE;
        }
    }

    static DistanceAttr attrLink(const DistanceAttr& e1, const DistanceAttr& e2) {
        DistanceAttr e; attrLink(e1, e2, &e); return e;
    }

    static void attrLink(const DistanceAttr& e1, DistanceAttr* e) {
        DistanceAttr e0(0, e1.distance < 0 ? Edge::ED_ANTISENSE : Edge::ED_SENSE, Edge::EC_SAME);
        return attrLink(e0, e1, e);
    }

    static DistanceAttr attrLink(const DistanceAttr& e1) {
        DistanceAttr e; attrLink(e1, &e); return e;
    }

    static bool diffDir(const DistanceAttr& e1, const DistanceAttr& e2) {
        return (e1.distance < 0 || e2.distance < 0) && (e1.distance >= 0 || e2.distance >= 0);
    }

    static bool hasLink(const Vertex* v1, const Vertex* v2, const DistanceAttr& e) {
        return hasLink(v1, v2, e.distance, e.dir, e.comp);
    }

    static bool hasLink(const Vertex* v1, const DistanceAttr& e1, const Vertex* v2, const DistanceAttr& e2) {
        assert(!diffDir(e1, e2));
        if (std::abs(e1.distance) > std::abs(e2.distance)) {
            return hasLink(v2, e2, v1, e1);
        }

        DistanceAttr e = attrLink(e1, e2);
        return hasLink(v1, v2, e);
    }

    static bool hasLink(const NodePtr& v1, const NodePtr& v2) {
        return hasLink(v1->vertex, v1->attr, v2->vertex, v2->attr);
    }
};

#endif // bigraph_search_h_

