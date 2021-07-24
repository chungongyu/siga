#include "suffix_array.h"

#include <fstream>

#include "utils.h"

static const uint16_t FILE_MAGIC = 0xCACA;

//
// Write a suffix array file to disk
//
class SAWriter {
 public:
  SAWriter(std::ostream& stream) : _stream(stream) {
  }
  bool write(const SuffixArray& sa) {
    if (!writeHeader(sa.strings(), sa.strings())) {
      return false;
    }
    for (size_t i = 0; i < sa.size(); ++i) {
      if (!writeElem(sa[i])) {
        return false;
      }
    }
    return true;
  }

 private:
  bool writeHeader(size_t strings, size_t elems) {
    if (_stream) {
      _stream << FILE_MAGIC << "\n";
      _stream << strings << "\n" << elems << "\n";
    }
    return (bool)_stream;
  }
  bool writeElem(const SuffixArray::Elem& elem) {
    if (_stream) {
      if (elem.full()) {
        _stream << elem.i << ' ' << elem.j << "\n";
      }
    }
    return (bool)_stream;
  }
  std::ostream& _stream;
};

std::ostream& operator<<(std::ostream& stream, const SuffixArray& sa) {
  SAWriter writer(stream);
  writer.write(sa);
  return stream;
}

//
// Read a suffix array file from disk
//
class SAReader {
 public:
  SAReader(std::istream& stream) : _stream(stream) {
  }
  bool read(SuffixArray& sa) {
    size_t elems = 0;
    if (!readHeader(sa._strings, elems)) {
      return false;
    }
    sa._elems.resize(elems);
    for (auto& elem : sa._elems) {
      if (!readElem(elem)) {
        return false;
      }
    }
    return true;
  }

 private:
  bool readHeader(size_t& strings, size_t& elems) {
    if (_stream) {
      uint16_t magic = 0;
      _stream >> magic;
      if (magic != FILE_MAGIC) {
        return false;
      }
      _stream >> strings >> elems;
    }
    return (bool)_stream;
  }
  bool readElem(SuffixArray::Elem& elem) {
    if (_stream) {
      _stream >> elem.i;
      _stream >> elem.j;
    }
    return (bool)_stream;
  }

  std::istream& _stream;
};

std::istream& operator>>(std::istream& stream, SuffixArray& sa) {
  SAReader reader(stream);
  reader.read(sa);
  return stream;
}

SuffixArray* SuffixArray::load(const std::string& filename) {
  std::ifstream stream(filename.c_str());
  return load(stream);
}

SuffixArray* SuffixArray::load(std::istream& stream) {
  SuffixArray* sa = new SuffixArray();

  SAReader reader(stream);
  if (reader.read(*sa)) {
    return sa;
  }

  SAFE_DELETE(sa);
  return nullptr;
}
