#include "bwt.h"

#include <memory>

#include "suffix_array.h"

BWT::BWT(const SuffixArray& sa, const DNASeqList& sequences)
    : _strings(sa.strings()), _suffixes(0) {
  RLUnit run;

  char c = '\0';
  std::unique_ptr<SuffixArray::BWTTraveller> travel(sa.travel());
  assert(travel);
  while (travel->next(c)) {
    ++_suffixes;
    if (run.initialized()) {
      if (run == c && !run.full()) {
        ++run;
      } else {
        // Write out the old run and start a new one
        _runs.push_back(run);
        run = RLUnit(c);
      }
    } else {
      // Start a new run
      run = RLUnit(c);
    }
  }
  if (run.initialized()) {
    _runs.push_back(run);
  }
}

const uint16_t BWT_FILE_MAGIC = 0xCACA;

enum BWFlag {
  BWF_NOFMI = 0,
  BWF_HASFMI
};

//
// Read a run length encoded binary BWT file from disk
//
class BWTReader {
 public:
  BWTReader(std::istream& stream) : _stream(stream) {
  }

  bool read(BWT& bwt);

 private:
  bool readHeader(uint64_t& num_strings, uint64_t& num_suffixes, BWFlag& flag);
  bool readRuns(RLString& runs, uint64_t numRuns);

  std::istream& _stream;
  uint64_t _numRuns;
};

bool BWTReader::read(BWT& bwt) {
  BWFlag flag;
  if (!readHeader(bwt._strings, bwt._suffixes, flag)) {
    return false;
  }
  if (!readRuns(bwt._runs, _numRuns)) {
    return false;
  }
  return true;
}

bool BWTReader::readHeader(uint64_t& num_strings, uint64_t& num_suffixes, BWFlag& flag) {
  uint16_t magic;
  if (!_stream.read((char *)&magic, sizeof(magic)) || magic != BWT_FILE_MAGIC) {
    return false;
  }
  if (!_stream.read((char *)&num_strings, sizeof(num_strings))) {
    return false;
  }
  if (!_stream.read((char *)&num_suffixes, sizeof(num_suffixes))) {
    return false;
  }
  if (!_stream.read((char *)&_numRuns, sizeof(_numRuns))) {
    return false;
  }
  if (!_stream.read((char *)&flag, sizeof(flag))) {
    return false;
  }
  return true;
}

bool BWTReader::readRuns(RLString& runs, uint64_t numRuns) {
  runs.resize(numRuns);
  if (!runs.empty()) {
    if (!_stream.read((char *)&runs[0], numRuns * sizeof(runs[0]))) {
      return false;
    }
  }
  return true;
}

//
// Write a run-length encoded BWT to a binary file
//
class BWTWriter {
 public:
  BWTWriter(std::ostream& stream) : _stream(stream), _numRuns(0), _posRun(0) {
  }

  bool write(const BWT& bwt);

 private:
  bool writeHeader(uint64_t num_strings, uint64_t num_suffixes, BWFlag flag);
  bool writeRun(const RLUnit& run);
  bool finalize();

  uint64_t _numRuns;
  uint64_t _posRun;

  std::ostream& _stream;
};

bool BWTWriter::write(const BWT& bwt) {
  uint64_t num_strings = bwt._strings, num_suffixes = bwt._suffixes;
  if (!writeHeader(num_strings, num_suffixes, BWF_NOFMI)) {
    return false;
  }
  for (const auto& run : bwt._runs) {
    if (!writeRun(run)) {
      return false;
    }
  }
  if (!finalize()) {
    return false;
  }
  return true;
}

bool BWTWriter::writeHeader(uint64_t num_strings, uint64_t num_suffixes, BWFlag flag) {
  if (!_stream.write((const char *)&BWT_FILE_MAGIC, sizeof(BWT_FILE_MAGIC))) {
    return false;
  }
  if (!_stream.write((const char *)&num_strings, sizeof(num_strings))) {
    return false;
  }
  if (!_stream.write((const char *)&num_suffixes, sizeof(num_suffixes))) {
    return false;
  }

  // Here we do not know the number of runs that are going to be written to the file
  // so we save the offset in the file and write a dummy value. After the bwt string
  // has been written, we return here and fill in the correct value
  _posRun = _stream.tellp();
  _numRuns = 0;
  if (!_stream.write((const char *)&_numRuns, sizeof(_numRuns))) {
    return false;
  }

  assert(flag == BWF_NOFMI);
  if (!_stream.write((const char *)&flag, sizeof(flag))) {
    return false;
  }

  return true;
}

bool BWTWriter::finalize() {
  _stream.seekp(_posRun);
  _stream.write((const char *)&_numRuns, sizeof(_numRuns));
  _stream.seekp(std::ios_base::end);
  return (bool)_stream;
}

bool BWTWriter::writeRun(const RLUnit& run) {
  if (!_stream.write((const char *)&run.data, sizeof(run.data))) {
    return false;
  }
  ++_numRuns;
  return true;
}

std::ostream& operator<<(std::ostream& stream, const BWT& bwt) {
  BWTWriter w(stream);
  w.write(bwt);
  return stream;
}

std::istream& operator>>(std::istream& stream, BWT& bwt) {
  BWTReader r(stream);
  r.read(bwt);
  return stream;
}
