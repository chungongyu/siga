#include "suffix_array_builder.h"

#include <cstring>
#include <numeric>
#include <iterator>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

#include <bcr.h>
#include <sais.hxx>

#include "alphabet.h"
#include "mkqs.h"
#include "suffix_array.h"
#include "utils.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.SuffixArrayBuilder"));

//
// Implementation of induced copying algorithm by Nong, Zhang, Chan
// Follows implementation given as an appendix to their 2008 paper
// '\0' is the sentinenl in this algorithm
//
class SAISBuilder : public SuffixArrayBuilder {
 public:
  SuffixArray* build(const DNASeqList& reads, size_t threads = 1) {
    assert(!reads.empty());

    size_t num_strings = reads.size();

    // In the multiple strings case, we need a 2D bit array
    // to hold the L/S types for the suffixes
    char** type_array = new char*[num_strings];
    for (size_t i = 0; i < num_strings; ++i) {
      const DNASeq& read = reads[i];
      size_t num_bytes = (read.seq.length() + 1) / 8 + 1;
      type_array[i] = new char[num_bytes];
      memset(type_array[i], 0, num_bytes);
    }

    // Classify each suffix as being L or S type
    for (size_t i = 0; i < num_strings; ++i) {
      const DNASeq& read = reads[i];
      size_t len = read.seq.length() + 1;

      // The empty suffix ($) for each string is defined to be S type
      // and hence the next suffix must be L type
      setBit(type_array, i, len - 1, 1);
      if (!read.seq.empty()) {
        setBit(type_array, i, len - 2, 0);
        for (size_t j = len - 2; j > 0; --j) {
          char curr = read.seq[j - 1], next = read.seq[j];
          bool type = (curr < next || (curr == next && getBit(type_array, i, j) == 1));
          setBit(type_array, i, j - 1, type);
        }
      }
    }

    // setup buckets
    size_t bucket_counts[DNAAlphabet::ALL_SIZE];
    size_t buckets[DNAAlphabet::ALL_SIZE];

    // find the ends of the buckets
    countBuckets(reads, bucket_counts, DNAAlphabet::ALL_SIZE);
    // getBuckets(bucket_counts, buckets, DNAAlphabet::ALL_SIZE, true);

    // Initialize the suffix array
    size_t num_suffixes = std::accumulate(&bucket_counts[0], &bucket_counts[0] + DNAAlphabet::ALL_SIZE, (size_t)0);
    LOG4CXX_DEBUG(logger, boost::format("initialize SA, strings: %d, suffixes: %d") % num_strings % num_suffixes);

    SuffixArray* sa = new SuffixArray(num_strings, num_suffixes);

    // Copy all the LMS substrings into the first n1 places in the SA
    size_t n1 = 0;
    for (size_t i = 0; i < num_strings; ++i) {
      const DNASeq& read = reads[i];
      for (size_t j = 0; j < read.seq.length() + 1; ++j) {
        if (isLMS(type_array, i, j)) {
          SuffixArray::Elem& ele = (*sa)[n1++];
          ele.i = i;
          ele.j = j;
        }
      }
    }

    // Call MKQS, first on the sequence and then on the index in the read table
    LOG4CXX_DEBUG(logger, boost::format("calling mkqs on %d of %d suffixes(%f), using %d threads")
        % n1 % num_suffixes % ((double)n1 / num_suffixes) % threads);
    {
      SuffixRadixCmp radixcmp(reads);
      SuffixIndexCmp indexcmp;
      if (threads <= 1) {
        mkqs2(&(*sa)[0], n1, 0, radixcmp, indexcmp);
      } else {
        mkqs_parallel(&(*sa)[0], n1, threads, radixcmp, indexcmp);
      }
    }
    LOG4CXX_DEBUG(logger, "mkqs finished");

    // Induction sort the remaining suffixes
    for (size_t i = n1; i < num_suffixes; ++i) {
      (*sa)[i] = SuffixArray::Elem();
    }

    // Find the ends of the buckets
    getBuckets(bucket_counts, buckets, DNAAlphabet::ALL_SIZE, true);

    for (size_t i = n1; i > 0; --i) {
      SuffixArray::Elem elem = (*sa)[i - 1];
      (*sa)[i - 1] = SuffixArray::Elem();  // empty
      const DNASeq& read = reads[elem.i];
      char c = read.seq[elem.j];
      (*sa)[--buckets[DNAAlphabet::torank(c)]] = elem;
    }

    induceSAl(reads, sa, type_array, bucket_counts, buckets, num_suffixes, DNAAlphabet::ALL_SIZE, false);
    induceSAs(reads, sa, type_array, bucket_counts, buckets, num_suffixes, DNAAlphabet::ALL_SIZE, true);

    // deallocate t array
    for (size_t i = 0; i < num_strings; ++i) {
      SAFE_DELETE_ARRAY(type_array[i]);
    }
    SAFE_DELETE_ARRAY(type_array);
    return sa;
  }

 private:
  class SuffixRadixCmp {
   public:
    SuffixRadixCmp(const DNASeqList& reads) : _reads(reads) {
    }

    // Get the character at position d for the SAElem
    char getChar(const SuffixArray::Elem& x, int d) const {
      const char* suffix = getChrPtr(x);
      return *(suffix + d);
    }
    // Get the suffix character string corresponding to this element
    const char* getChrPtr(const SuffixArray::Elem& x) const {
      const DNASeq& read = _reads[x.i];
      return read.seq.c_str() + x.j;
    }
   private:
    const DNASeqList& _reads;
  };
  // Compare two suffixes by their index in the read table
  // This is used for the final pass, after suffixes has been compared by sequence
  class SuffixIndexCmp {
   public:
    bool operator()(const SuffixArray::Elem& x, const SuffixArray::Elem& y) const {
      return x.i < y.i;
    }
  };

  void induceSAl(const DNASeqList& reads, SuffixArray* sa, char** type_array, size_t* counts, size_t* buckets, size_t n, size_t K, bool end) {
    getBuckets(counts, buckets, K, end);
    for (size_t i = 0; i < n; ++i) {
      const SuffixArray::Elem& ielem = (*sa)[i];
      if (!ielem.empty() && ielem.j > 0) {
        LOG4CXX_TRACE(logger, boost::format("Curr: %d %d") % ielem.i % ielem.j);

        SuffixArray::Elem jelem(ielem.i, ielem.j - 1);
        if (!getBit(type_array, jelem.i, jelem.j)) {
          const DNASeq& read = reads[jelem.i];
          char c = read.seq[jelem.j];
          LOG4CXX_TRACE(logger,  boost::format("<iSA1>Placing %d %d at position %d")
              % jelem.i % jelem.j % buckets[DNAAlphabet::torank(c)]);
          (*sa)[buckets[DNAAlphabet::torank(c)]++] = jelem;
        }
      }
    }
  }

  void induceSAs(const DNASeqList& reads, SuffixArray* sa, char** type_array, size_t* counts, size_t* buckets, size_t n, size_t K, bool end) {
    getBuckets(counts, buckets, K, end);
    for (size_t i = n; i > 0; --i) {
      const SuffixArray::Elem& ielem = (*sa)[i - 1];
      if (!ielem.empty() && ielem.j > 0) {
        LOG4CXX_TRACE(logger, boost::format("Curr: %d %d") % ielem.i % ielem.j);

        SuffixArray::Elem jelem(ielem.i, ielem.j - 1);
        if (getBit(type_array, jelem.i, jelem.j)) {
          const DNASeq& read = reads[jelem.i];
          char c = read.seq[jelem.j];
          LOG4CXX_TRACE(logger,  boost::format("<iSA1>Placing %d %d at position %d")
              % jelem.i % jelem.j % (buckets[DNAAlphabet::torank(c)] - 1));
          (*sa)[--buckets[DNAAlphabet::torank(c)]] = jelem;
        }
      }
    }
  }

  // Calculate the number of items that should be in each bucket
  void countBuckets(const DNASeqList& reads, size_t* counts, size_t K) {
    for (size_t i = 0; i < K; ++i) {
      counts[i] = 0;
    }
    for (size_t i = 0; i < reads.size(); ++i) {
      const DNASeq& read = reads[i];
      size_t len = read.seq.length();
      for (size_t j = 0; j < len; ++j) {
        ++counts[DNAAlphabet::torank(read.seq[j])];
      }
      ++counts[DNAAlphabet::torank('\0')];
    }
  }
  // If end is true, calculate the end of the buckets, otherwise
  // calculate the starts
  void getBuckets(size_t* counts, size_t* buckets, size_t K, bool end) {
    for (size_t i = 0; i < K; ++i) {
      buckets[i] = 0;
    }

    size_t sum = 0;
    for (size_t i = 0; i < K; ++i) {
      sum += counts[i];
      buckets[i] = end ? sum : sum - counts[i];
    }
  }

  // set the element to b
  void setBit(char** bit_array, size_t i, size_t j, bool b) {
    char* ba = bit_array[i];
    size_t block = j / 8, offset = j % 8;
    ba[block] = (b ? (_MASK[offset] | ba[block]) : (~_MASK[offset] & ba[block]));
  }
  bool getBit(char** bit_array, size_t i, size_t j) {
    return bit_array[i][j / 8] & _MASK[j % 8] ? 1 : 0;
  }

  bool isLMS(char** bit_array, size_t i, size_t j) {
    return j > 0 && getBit(bit_array, i, j) && !getBit(bit_array, i, j - 1);
  }

  static unsigned char _MASK[8];
};

unsigned char SAISBuilder::_MASK[8] = {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};

class NewSAISBuilder : public SuffixArrayBuilder {
 public:
  SuffixArray* build(const DNASeqList& reads, size_t /*threads = 1*/) {
    assert(!reads.empty());

    Helper h(reads);
    return h.build();
  }

 private:
  struct Helper {
    typedef int64_t MyInt;

    Helper(const DNASeqList& reads) : _reads(reads) {
      size_t suffixes = 0;

      for (size_t i = 0; i < _reads.size(); ++i) {
        suffixes += _reads[i].seq.length() + 1;
      }

      LOG4CXX_INFO(logger, boost::format("initialize index, strings: %lu, suffixes: %lu") % _reads.size() % suffixes);
      size_t k = 0;
      _indices.resize(suffixes);
      for (size_t i = 0; i < _reads.size(); ++i) {
        const auto& read = _reads[i];
        for (size_t j = 0; j < read.seq.length() + 1; ++j) {
          auto& elem = _indices[k++];
          elem.i = i;
          elem.j = j;
        }
      }
      assert(k == suffixes);
    }

    SuffixArray* build() {
      size_t n = _indices.size();
      Iterator T(this);
      std::vector<MyInt> S(n);
      int r = saisxx(T, &S[0], (MyInt)n, (MyInt)DNAAlphabet::ALL_SIZE);
      LOG4CXX_INFO(logger, boost::format("saisxx: %d") % r);
      if (r != 0) {
        return nullptr;
      }
      SuffixArray* sa = new SuffixArray(_reads.size(), _indices.size());
      for (size_t i = 0; i < n; ++i) {
        (*sa)[i] = _indices[S[i]];
      }
      return sa;
    }

   private:
    struct DNA {
      DNA(const Helper* h = nullptr, MyInt idx = 0) : _h(h), _e(nullptr) {
        if (_h != nullptr) {
          _e = &_h->_indices[idx];
        }
      }
      bool operator==(const DNA& o) const {
        if (_e == o._e) {
          return true;
        }
        char x = tochar(), y = o.tochar();
        return x == y && x != '\0';
      }
      bool operator!=(const DNA& o) const {
        return !(*this == o);
      }
      bool operator<(const DNA& o) const {
        int x = DNAAlphabet::torank(tochar()), y = DNAAlphabet::torank(o.tochar());
        return x < y || (x == 0 && _e->i < o._e->i);
      }
      bool operator>(const DNA& o) const {
        int x = DNAAlphabet::torank(tochar()), y = DNAAlphabet::torank(o.tochar());
        return x > y || (y == 0 && _e->i > o._e->i);
      }
      bool operator<=(const DNA& o) const {
        return !(*this > o);
      }
      bool operator>=(const DNA& o) const {
        return !(*this < o);
      }
      operator size_t() const {
        assert((_h && _e) || (!_h && !_e));
        if (_h != nullptr) {
          return DNAAlphabet::torank(tochar());
        }
        return 0;
      }

     private:
      char tochar() const {
        if (_h != nullptr) {
          assert(_e != nullptr);
          assert(_e->i < _h->_reads.size());
          const auto& read = _h->_reads[_e->i];
          if (_e->j < read.seq.length()) {
            return read.seq[_e->j];
          }
          assert(_e->j == read.seq.length());
        }
        return '\0';
      }

      const Helper* _h;
      const SuffixArray::Elem* _e;
    };

    class Iterator : public std::iterator<std::random_access_iterator_tag, DNA, MyInt, MyInt*, DNA> {
     public:
      Iterator(const Helper* h = nullptr, MyInt idx = 0) : _h(h), _idx(idx) {
      }

      reference operator*() const {
        return DNA(_h, _idx);
      }
      bool operator==(const Iterator& o) const {
        return _idx == o._idx;
      }
      bool operator!=(const Iterator& o) const {
        return !(*this == o);
      }

      // Forward iterator requirements
      Iterator& operator++() {
        ++_idx;
        return *this;
      }
      Iterator operator++(int) {
        Iterator o = *this;
        ++(*this);
        return o;
      }
      // Bidirectional iterator requirements
      Iterator& operator--() {
        --_idx;
        return *this;
      }
      Iterator operator--(int) {
        Iterator o = *this;
        --(*this);
        return o;
      }
      // Random access iterator requirements
      reference operator[](const difference_type& n) const {
        return DNA(_h, _idx + n);
      }
      Iterator operator+=(const difference_type& n) {
        _idx += n;
        return *this;
      }
      Iterator operator+(const difference_type& n) const {
        return Iterator(_h, _idx + n);
      }
      Iterator operator-=(const difference_type& n) {
        _idx -= n;
        return *this;
      }
      Iterator operator-(const difference_type& n) const {
        return Iterator(_h, _idx - n);
      }

     private:
      const Helper* _h;
      MyInt _idx;
    };

    const DNASeqList& _reads;
    std::vector<SuffixArray::Elem> _indices;
  };
};

class RopeBuilder : public SuffixArrayBuilder {
 public:
  SuffixArray* build(const DNASeqList& reads, size_t threads = 1) {
    assert(!reads.empty());
    return nullptr;
  }
};

SuffixArrayBuilder* SuffixArrayBuilder::create(const std::string& algorithm) {
  if (boost::algorithm::iequals(algorithm, "sais")) {
    return new SAISBuilder();
  } else if (boost::algorithm::iequals(algorithm, "sais2")) {
    return new NewSAISBuilder();
  } else if (boost::algorithm::iequals(algorithm, "rope")) {
    return new RopeBuilder();
  }
  return nullptr;
}
