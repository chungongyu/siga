#ifndef alphabet_h_
#define alphabet_h_

#include <cstdint>
#include <cstring>

#include <algorithm>
#include <iterator>
#include <utility>

namespace DNAAlphabet {

const size_t ALL_SIZE = 5;
const char DNA_ALL[ALL_SIZE] = {'$', 'A', 'C', 'G', 'T'};

static const size_t size = 4;
static const char DNA[size] = {'A', 'C', 'G', 'T'};

inline int torank(char c) {
  static int RANK_ALL[256] = {
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,
      0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
  return RANK_ALL[(size_t)c];
}

inline char tochar(size_t rank) {
  return DNA_ALL[rank];
}

template <class Storage>
class AlphaCount {
 public:
  AlphaCount() {
    memset(_data, 0, sizeof(_data));
  }

  Storage& operator[](size_t i) {
    return _data[i];
  }
  const Storage& operator[](size_t i) const {
    return _data[i];
  }
  size_t size() const {
    return ALL_SIZE;
  }
  bool hasDNA() const  {
    for (size_t i = 1; i < size(); ++i) {
      if (_data[i] > 0) {
        return true;
      }
    }
    return false;
  }
  void complement() {
    std::swap(_data[torank('A')], _data[torank('T')]);
    std::swap(_data[torank('C')], _data[torank('G')]);
  }
  AlphaCount<Storage> operator-(const AlphaCount<Storage>& c) const {
    AlphaCount<Storage> r;
    for (size_t i = 0; i < size(); ++i) {
      r._data[i] = _data[i] - c._data[i];
    }
    return r;
  }
  AlphaCount<Storage> operator+(const AlphaCount<Storage>& c) const {
    AlphaCount<Storage> r;
    for (size_t i = 0; i < size(); ++i) {
      r._data[i] = _data[i] + c._data[i];
    }
    return r;
  }
  AlphaCount<Storage>& operator-=(const AlphaCount<Storage>& c) {
    for (size_t i = 0; i < size(); ++i) {
      _data[i] -= c._data[i];
    }
    return *this;
  }
  AlphaCount<Storage>& operator+=(const AlphaCount<Storage>& c) {
    for (size_t i = 0; i < size(); ++i) {
      _data[i] += c._data[i];
    }
    return *this;
  }

 private:
  template <class T>
  friend std::ostream& operator<<(std::ostream& stream, const AlphaCount<T>& c);
  template <class T>
  friend std::istream& operator>>(std::istream& stream, AlphaCount<T>& c);

  Storage _data[ALL_SIZE];
};

typedef AlphaCount<uint64_t> AlphaCount64;
typedef AlphaCount<uint32_t> AlphaCount32;
typedef AlphaCount<uint16_t> AlphaCount16;

template <class T>
std::ostream& operator<<(std::ostream& stream, const AlphaCount<T>& c) {
  std::copy(c._data, c._data + c.size(), std::ostream_iterator<T>(stream, " "));
  return stream;
}

};  // namespace DNAAlphabet

#endif  // alphabet_h_
