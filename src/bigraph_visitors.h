#ifndef bigraph_visitors_h_
#define bigraph_visitors_h_

#include <iostream>
#include <unordered_map>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H
#ifdef HAVE_MLPACK
#include "mlpack.h"
#endif  // HAVE_MLPACK
#include "bigraph.h"

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
  ChimericVisitor(size_t minLength, size_t minCoverage, size_t delta, size_t N, size_t G, double T = 13.0) : _minLength(minLength), _minCoverage(minCoverage), _delta(delta), _N(N), _G(G), _T(T) {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);
 private:
  size_t _minLength;
  size_t _minCoverage;
  size_t _N; // num of reads
  size_t _G; // genome size
  double _T; // threshold

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

// Remove identical vertices from the graph
class IdenticalRemoveVisitor : public BigraphVisitor {
 public:
  IdenticalRemoveVisitor() : _count(0) {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);
 private:
  size_t _count;
};

// Remove loops from the graph
class LoopRemoveVisitor : public BigraphVisitor {
 public:
  LoopRemoveVisitor() {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);
 private:
  std::vector<Vertex *> _loops;
};

// Run the YU LIN's maximum overlap algorithm on each node
class MaximumOverlapVisitor : public BigraphVisitor {
 public:
  MaximumOverlapVisitor(size_t delta, bool carefully, size_t N, size_t G, double T = 13.0) : _delta(delta), _carefully(carefully), _N(N), _G(G), _T(T), _dummys(0) {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);
 private:
  size_t _delta;
  size_t _dummys;
  bool _carefully;

  size_t _N; // num of reads
  size_t _G; // genome size
  double _T; // threshold
};

// Visit each paired node to estimate insert size variant under normal distribution assumption.
class InsertSizeEstimateVisitor : public BigraphVisitor {
 public:
  InsertSizeEstimateVisitor(size_t& average, size_t& delta) : _average(average), _delta(delta) {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);
 private:
  EdgePtrList edges(const Vertex* vertex, Edge::Dir dir);
  std::vector<size_t> _samples;

  size_t& _average;
  size_t& _delta;
};

// Visit each paired node via zigzag, sweep false positive edges.
class PairedReadVisitor : public BigraphVisitor {
 public:
  PairedReadVisitor(size_t maxDistance, size_t insertSize, size_t insertDelta, size_t maxNodes, size_t threads=1, size_t batch=1000) : _maxDistance(maxDistance), _insertSize(insertSize), _insertDelta(insertDelta), _maxNodes(maxNodes), _threads(threads), _batch(batch) {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);
 private:
  size_t _maxDistance;
  size_t _insertSize;
  size_t _insertDelta;
  size_t _maxNodes;
  size_t _threads;
  size_t _batch;

  std::vector<const Vertex *> _vertices;
  friend class PairedVertexProcess;
  friend class PairedVertexPostProcess;
};

// Visit each node via 10x linked read
class LinkedReadVisitor : public BigraphVisitor {
 public:
  LinkedReadVisitor(size_t minLength = -1, size_t minCoverage = -1) : _minLength(minLength), _minCoverage(minCoverage), _dummys(0) {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);
 private:
  size_t _minLength;
  size_t _minCoverage;
  size_t _dummys;
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

// Run the Myers transitive reduction algorithm on each node
class TransitiveReductionVisitor : public BigraphVisitor {
 public:
  TransitiveReductionVisitor() {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);
 private:
};

// Detects and removes small "tip" vertices from the graph
// when they are less than minLength in size
class TrimVisitor : public BigraphVisitor {
 public:
  TrimVisitor(size_t minLength, size_t minCoverage) : _minLength(minLength), _minCoverage(minCoverage) {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);
 private:
  size_t _minLength;
  size_t _minCoverage;

  size_t _island;
  size_t _terminal;
};

#ifdef HAVE_MLPACK
// Detects and removes all false positive edges from the graph
class AIVisitor : public BigraphVisitor {
 public:
  AIVisitor(void* model, size_t n, size_t G) : _model(model), _n(n), _G(G), _blacks(0), _whites(0), _grays(0) {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);

 private:
  template <typename T>
  T* model() {
    if (_model != nullptr) {
      return &((AIModel<T> *)_model)->model;
    }
    return nullptr;
  }
  void* _model;
  size_t _n;
  size_t _G;

  size_t _blacks;
  size_t _whites;
  size_t _grays;
};
#endif  // HAVE_MLPACK

class UnitigVisitor : public BigraphVisitor {
 public:
  UnitigVisitor(size_t n, size_t G, size_t T) : _n(n), _G(G), _T(T), _unitigs(0) {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);

 private:
  size_t _n;
  size_t _G;
  double _T; // threshold

  size_t _unitigs;
};

class FMIndex;

class GANVisitor : public BigraphVisitor {
 public:
  GANVisitor(std::ostream& stream, FMIndex* ref = nullptr) : _stream(stream), _ref(ref), _blacks(0), _whites(0), _grays(0) {
  }
  void previsit(Bigraph* graph);
  bool visit(Bigraph* graph, Vertex* vertex);
  void postvisit(Bigraph* graph);

 private:
  std::ostream& _stream;
  FMIndex* _ref;
  size_t _blacks;
  size_t _whites;
  size_t _grays;
};

#endif  // bigraph_visitors_h_
