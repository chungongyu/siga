#include "suffix_array.h"

#include <boost/foreach.hpp>

template< const DNASeqList& reads >
struct DNAIndex {
    DNAIndex(size_t i) : _reads(reads) {
    }
private:
    SuffixArray::Elem _index;
    const DNASeqList& _reads;
};

class ReadList {
public:
    ReadList(const DNASeqList& reads) : _reads(reads) {
    }
private:
    const DNASeqList& _reads;
};

static const uint16_t FILE_MAGIC = 0xCACA;

class SAWriter {
public:
    SAWriter(std::ostream& stream) : _stream(stream) {
    }
    bool write(const SuffixArray& sa) {
        if (!writeHeader(sa._strings, sa._elems.size())) {
            return false;
        }
        BOOST_FOREACH(const SuffixArray::Elem& elem, sa._elems) {
            if (!writeElem(elem)) {
                return false;
            }
        }
        return true;
    }

private:
    bool writeHeader(size_t strings, size_t elems) {
        _stream << FILE_MAGIC;
        _stream << strings << elems;
        return true;
    }
    bool writeElem(const SuffixArray::Elem& elem) {
        _stream << elem.i << elem.j;
        return true;
    }
    std::ostream& _stream;
};

std::ostream& operator<<(std::ostream& stream, const SuffixArray& sa) {
    SAWriter writer(stream);
    writer.write(sa);
    return stream;
}

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
        BOOST_FOREACH(SuffixArray::Elem& elem, sa._elems) {
            if (!readElem(elem)) {
                return false;
            }
        }
        return true;
    }

private:
    bool readHeader(size_t& strings, size_t& elems) {
        uint16_t magic = 0;
        _stream >> magic;
        if (magic != FILE_MAGIC) {
            return false;
        }
        _stream >> strings >> elems;
        return true;
    }
    bool readElem(SuffixArray::Elem& elem) {
        _stream >> elem.i >> elem.j;
        return true;
    }

    std::istream& _stream;
};

std::istream& operator>>(std::istream& stream, SuffixArray& sa) {
    SAReader reader(stream);
    reader.read(sa);
    return stream;
}

