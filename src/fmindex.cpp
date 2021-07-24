#include "fmindex.h"

#include <fstream>

#include <boost/format.hpp>

#include <log4cxx/logger.h>

#include "bwt.h"
#include "utils.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.FMIndex"));

//
// MarkerFill
//
template <class MarkerList>
class MarkerFill {
 public:
  MarkerFill(MarkerList& markers, size_t symbols, size_t sampleRate)
      : _markers(markers), _sampleRate(sampleRate) {
    initialize(symbols);
  }
  ~MarkerFill() {
    assert(_currIdx == _markers.size());
  }

  virtual void fill(const DNAAlphabet::AlphaCount64& counts, uint64_t total, size_t unitIndex, bool lastOne) = 0;

 protected:
  void initialize(size_t n) {
    // we place a marker at the beginning (with no accumulated counts), every sampleRate
    // bases and one at the very end (with the total counts)
    size_t required_markers = (n % _sampleRate == 0) ? (n / _sampleRate) + 1 : (n / _sampleRate) + 2;
    _markers.resize(required_markers);

    // Place a blank markers at the start of the data
    if (!_markers.empty()) {
      _markers[0].unitIndex = 0;
    }
    _currIdx = 1;
    _nextPos = _sampleRate;
  }

  MarkerList& _markers;
  size_t _sampleRate;

  size_t _currIdx;
  size_t _nextPos;
};

class LargeMarkerFill : public MarkerFill<LargeMarkerList> {
 public:
  LargeMarkerFill(LargeMarkerList& markers, size_t symbols, size_t sampleRate)
      : MarkerFill<LargeMarkerList>(markers, symbols, sampleRate) {
  }

  void fill(const DNAAlphabet::AlphaCount64& counts, uint64_t total, size_t unitIndex, bool lastOne) {
    // Check whether to place a new marker
    while (total >= _nextPos || lastOne) {
      // Sanity checks
      // The marker position should always be less than the running total unless
      // the number of symbols is smaller than the sample rate
      size_t expectedPos = _currIdx * _sampleRate;
      assert(expectedPos <= total || lastOne);
      assert(_currIdx < _markers.size());

      LargeMarker& marker = _markers[_currIdx++];
      marker.unitIndex = unitIndex;
      marker.counts = counts;

      _nextPos += _sampleRate;
      lastOne = lastOne && _currIdx < _markers.size();
    }
  }
};

class SmallMarkerFill : public MarkerFill<SmallMarkerList> {
 public:
  SmallMarkerFill(const LargeMarkerList& lmarkers, SmallMarkerList& smarkers, size_t symbols, size_t sampleRate)
      : MarkerFill<SmallMarkerList>(smarkers, symbols, sampleRate), _lmarkers(lmarkers) {
  }

  void fill(const DNAAlphabet::AlphaCount64& counts, uint64_t total, size_t unitIndex, bool lastOne) {
    while (total >= _nextPos || lastOne) {
      // Sanity checks
      // The marker position should always be less than the running total unless
      // the number of symbols is smaller than the sample rate
      size_t expectedPos = _currIdx * _sampleRate;
      assert(expectedPos <= total || lastOne);
      assert(_currIdx < _markers.size());

      // Calculate the large marker to set the relative count from
      // This is generally the most previously placed large block except it might
      // be the second-previous in the case that we placed the last large marker.
      size_t largeIndex = expectedPos / DEFAULT_SAMPLE_RATE_LARGE;
      const LargeMarker& lmarker = _lmarkers[largeIndex];
      SmallMarker& smarker = _markers[_currIdx++];

      for (size_t i = 0; i < lmarker.counts.size(); ++i) {
        // assert(counts[i] - lmarker.counts);
        smarker.counts[i] = counts[i] - lmarker.counts[i];
      }
      smarker.unitIndex = unitIndex - lmarker.unitIndex;

      _nextPos += _sampleRate;
      lastOne = lastOne && _currIdx < _markers.size();
    }
  }

  const LargeMarkerList& _lmarkers;
};

std::ostream& operator<<(std::ostream& stream, const FMIndex::Interval& interval) {
  stream << interval.lower << ' ' << interval.upper;
  return stream;
}

std::istream& operator>>(std::istream& stream, FMIndex::Interval& interval) {
  stream >> interval.lower >> interval.upper;
  return stream;
}

void FMIndex::initialize() {
  assert(IS_POWER_OF_2(_sampleRate));

  DNAAlphabet::AlphaCount64 counts;
  uint64_t total = 0;

  // Fill in the marker values
  // We wish to place markers every sampleRate symbols however since a run may
  // not end exactly on sampleRate boundaries, we place the markers AFTER
  // the run crossing the boundary ends

  LargeMarkerFill f1(_lmarkers, _bwt.length(), DEFAULT_SAMPLE_RATE_LARGE);
  SmallMarkerFill f2(_lmarkers, _smarkers, _bwt.length(), _sampleRate);

  const RLString& runs = _bwt.str();
  for (size_t i = 0; i < runs.size(); ++i) {
    const RLUnit& run = runs[i];
    char c = (char)run;
    size_t len = run.count();

    // Update the count and advance the running total
    counts[DNAAlphabet::torank(c)] += len;
    total += len;

    // Check whether to place a new large marker
    f1.fill(counts, total, i + 1, i == runs.size() - 1);

    // Check whether to place a new small marker
    f2.fill(counts, total, i + 1, i == runs.size() - 1);
  }

  // Initialize C(a)
  _pred[DNAAlphabet::torank(DNAAlphabet::DNA_ALL[0])] = 0;
  for (size_t i = 1; i < DNAAlphabet::ALL_SIZE; ++i) {
    size_t curr = DNAAlphabet::torank(DNAAlphabet::DNA_ALL[i]), prev = DNAAlphabet::torank(DNAAlphabet::DNA_ALL[i - 1]);
    _pred[curr] = _pred[prev] + counts[prev];
  }
}

void FMIndex::info() const {
  const RLString& runs = _bwt.str();
  LOG4CXX_INFO(logger, "fm-index info:");
  LOG4CXX_INFO(logger, boost::format("large sample rate: %d") % DEFAULT_SAMPLE_RATE_LARGE);
  LOG4CXX_INFO(logger, boost::format("small sample rate: %d") % _sampleRate);
  LOG4CXX_INFO(logger, boost::format("contains %lu symbols in %lu runs (%1.4lf symbols per run)")
      % _bwt.length() % runs.size()
      % (runs.empty() ? 0 : ((double)_bwt.length() / runs.size())));
  LOG4CXX_INFO(logger, boost::format("marker memory -- small markers: %lu, large markers: %lu")
      % _smarkers.size() % _lmarkers.size());
}

//
// MarkerFind
//
class MarkerFind {
 public:
  MarkerFind(const RLString& runs, const LargeMarkerList& lmarkers, const SmallMarkerList& smarkers, size_t sampleRate)
      : _runs(runs), _lmarkers(lmarkers), _smarkers(smarkers), _sampleRate(sampleRate) {
  }
  size_t find(char c, size_t i) const {
    DNAAlphabet::AlphaCount64 counts = find(i);
    return counts[DNAAlphabet::torank(c)];
  }

  DNAAlphabet::AlphaCount64 find(size_t i) const {
    // The counts in the marker are not inclusive (unlike the Occurrence class)
    // so we increment the index by 1.
    ++i;

    LargeMarker lmarker = nearest(i);
    size_t position = lmarker.total();

    DNAAlphabet::AlphaCount64 counts = lmarker.counts;
    size_t currIdx = lmarker.unitIndex;

    // Search forwards (towards 0) until idx is found
    while (position < i) {
      size_t delta = i - position;

      assert(currIdx < _runs.size());

      const RLUnit& run = _runs[currIdx++];
      size_t n = run.count();
      if (n > delta) {
        n = delta;
      }
      counts[DNAAlphabet::torank((char)run)] += n;
      position += n;
    }
    // Search backwards (towards 0) until idx is found
    while (position > i) {
      size_t delta = position - i;

      assert(currIdx <= _runs.size());

      const RLUnit& run = _runs[--currIdx];
      size_t n = run.count();
      if (n > delta) {
        n = delta;
      }
      counts[DNAAlphabet::torank((char)run)] -= n;
      position -= n;
    }

    assert(position == i);

    return counts;
  }

  char getChar(size_t i) const {
    LargeMarker amarker = upper(i);
    size_t k = amarker.total();
    assert(k >= i);
    size_t unitIndex = amarker.unitIndex;
    while (k > i) {
      assert(unitIndex != 0);
      const RLUnit& run = _runs[--unitIndex];
      k -= run.count();
    }
    const RLUnit& run = _runs[unitIndex];
    assert(k <= i && k + (size_t)run >= i);
    return (char)run;
  }

 private:
  LargeMarker nearest(size_t i) const {
    size_t baseIdx = i / _sampleRate;
    size_t offset = MOD_POWER_2(i, _sampleRate);  // equivalent to position % sampleRate
    if (offset >= _sampleRate>>1) {
      ++baseIdx;
    }
    return interpolated(baseIdx);
  }

  // Get the lowest interpolated marker whose position is strictly greater than position
  LargeMarker upper(size_t i) const {
    size_t baseIdx = (i / _sampleRate) + 1;
    return interpolated(baseIdx);
  }

  // Return a LargeMarker with values that are interpolated by adding
  // the relative count nearest to the requested position to the last
  // LargeMarker
  LargeMarker interpolated(size_t smallIdx) const {
    // Calculate the position of the LargeMarker closest to the target SmallMarker
    size_t largeIdx = smallIdx * _sampleRate / DEFAULT_SAMPLE_RATE_LARGE;

    LargeMarker absolute = _lmarkers[largeIdx];
    const SmallMarker& relative = _smarkers[smallIdx];
    for (size_t i = 0; i < absolute.counts.size(); ++i) {
      absolute.counts[i] += relative.counts[i];
    }
    absolute.unitIndex += relative.unitIndex;

    return absolute;
  }

  const RLString& _runs;
  const LargeMarkerList& _lmarkers;
  const SmallMarkerList& _smarkers;
  size_t _sampleRate;
};

char FMIndex::getChar(size_t i) const {
  MarkerFind finder(_bwt.str(), _lmarkers, _smarkers, _sampleRate);
  return finder.getChar(i);
}

std::string FMIndex::getString(size_t i) const {
  assert(i < length());

  // The range [0,n) in the BWT contains all the terminal
  // symbols for the reads. Search backwards from one of them
  // until the '$' is found gives a full string.
  std::string out;

  Interval interval(i, i);
  while (true) {
    assert(interval.valid());
    char c = getChar(interval.lower);
    if (c == '$') {
      break;
    }
    out += c;
    interval.update(c, this);
  }

  std::reverse(out.begin(), out.end());
  return out;
}

size_t FMIndex::getOcc(char c, size_t i) const {
  MarkerFind finder(_bwt.str(), _lmarkers, _smarkers, _sampleRate);
  return finder.find(c, i);
}

DNAAlphabet::AlphaCount64 FMIndex::getOcc(size_t i) const {
  MarkerFind finder(_bwt.str(), _lmarkers, _smarkers, _sampleRate);
  return finder.find(i);
}

std::ostream& operator<<(std::ostream& stream, const FMIndex& index) {
  stream << index._bwt;
  /*
  stream << "lmarkers" << std::endl;
  for (size_t i = 0; i < index._lmarkers.size(); ++i) {
    const LargeMarker& marker = index._lmarkers[i];
    stream << "--------------" << std::endl;
    stream << marker.unitIndex << std::endl;
    stream << marker.counts << std::endl;
  }
  stream << "smarkers" << std::endl;
  for (size_t i = 0; i < index._smarkers.size(); ++i) {
    const SmallMarker& marker = index._smarkers[i];
    stream << "--------------" << std::endl;
    stream << marker.unitIndex << std::endl;
    stream << marker.counts << std::endl;
  }
  stream << index._pred << std::endl;
  */
  return stream;
}

std::istream& operator>>(std::istream& stream, FMIndex& index) {
  stream >> index._bwt;
  index.initialize();
  return stream;
}

bool FMIndex::load(std::istream& stream, FMIndex& fmi) {
  try {
    stream >> fmi;
    fmi.info();
  } catch (...) {
    return false;
  }
  return true;
}

bool FMIndex::load(const std::string& filename, FMIndex& fmi) {
  std::ifstream stream(filename.c_str());
  return load(stream, fmi);
}
