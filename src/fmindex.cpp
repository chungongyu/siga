#include "fmindex.h"
#include "bwt.h"
#include "utils.h"

#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.FMIndex"));

static size_t _NumRequiredMarkers(size_t n, size_t d) {
    // we place a marker at the beginning (with no accumulated counts), every sampleRate
    // bases and one at the very end (with the total counts)
    return (n % d == 0) ? (n / d) + 1 : (n / d) + 2;
}

template< class MarkerList >
class MarkerFill {
public:
    MarkerFill(MarkerList& markers, size_t symbols, size_t sampleRate) : _markers(markers), _sampleRate(sampleRate) {
        initialize(symbols);
    }

    virtual void fill(size_t* counts, size_t total, size_t unitIndex, bool lastOne) = 0;
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
    }

    MarkerList& _markers;
    size_t _sampleRate;
};

class LargeMarkerFill : public MarkerFill< LargeMarkerList > {
public:
    LargeMarkerFill(LargeMarkerList& markers, size_t symbols, size_t sampleRate) : MarkerFill< LargeMarkerList >(markers, symbols, sampleRate) {
        _currIdx = 1;
        _nextPos = _sampleRate;
    }

    void fill(size_t* counts, size_t total, size_t unitIndex, bool lastOne) {
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
            memcpy(marker.counts, counts, SIZEOF_ARRAY(marker.counts));

            _nextPos += _sampleRate;
            lastOne = lastOne && total >= _nextPos;
        }
    }

    size_t _currIdx;
    size_t _nextPos;
};

class SmallMarkerFill : public MarkerFill< SmallMarkerList > {
public:
    SmallMarkerFill(const LargeMarkerList& lmarkers, SmallMarkerList& smarkers, size_t symbols, size_t sampleRate) : MarkerFill< SmallMarkerList >(smarkers, symbols, sampleRate), _lmarkers(lmarkers) {
    }

    void fill(size_t* counts, size_t total, size_t unitIndex, bool lastOne) {
    };

    const LargeMarkerList& _lmarkers;
};

void FMIndex::initialize() {
    assert(IS_POWER_OF_2(_sampleRate));

    size_t counts[DNAAlphabet::ALL_SIZE];
    memset(counts, 0, sizeof(counts));
    size_t total = 0;

    // Fill in the marker values
    // We wish to place markers every sampleRate symbols however since a run may
    // not end exactly on sampleRate boundaries, we place the markers AFTER
    // the run crossing the boundary ends

    LargeMarkerFill f1(_lmarkers, _bwt.length(), DEFAULT_SAMPLE_RATE_LARGE);
    //SmallMarkerFill f2(_lmarkers, _smarkers, _bwt.length(), _sampleRate);

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
        //f2.fill(counts, total, i + 1, i == runs.size() - 1);
    }

    // Initialize C(a)
    memset(_pred, 0, sizeof(_pred));
    _pred[DNAAlphabet::torank(DNAAlphabet::DNA_ALL[0])] = 0;
    for (size_t i = 1; i < DNAAlphabet::ALL_SIZE; ++i) {
        size_t curr = DNAAlphabet::torank(DNAAlphabet::DNA_ALL[i]), prev = DNAAlphabet::torank(DNAAlphabet::DNA_ALL[i - 1]);
        _pred[curr] = _pred[prev] + counts[prev];
    }
}

void FMIndex::info() const {
    const RLString& runs = _bwt.str();
    LOG4CXX_INFO(logger, "FMIndex info:");
    LOG4CXX_INFO(logger, boost::format("Large Sample rate: %d") % DEFAULT_SAMPLE_RATE_LARGE);
    LOG4CXX_INFO(logger, boost::format("Small Sample rate: %d") % _sampleRate);
    LOG4CXX_INFO(logger, boost::format("Contains %d symbols in %d runs (%1.4lf symbols per run)") % _bwt.length() % runs.size() % (runs.empty() ? 0 : ((double)_bwt.length() / runs.size())));
    LOG4CXX_INFO(logger, boost::format("Marker Memory -- Small Markers: %d (%.1lf MB) Large Markers: %d (%.1lf MB)") % _smarkers.size() % 0 % _lmarkers.size() % 0);
}

std::ostream& operator<<(std::ostream& stream, const FMIndex& index) {
    stream << index._bwt;
    return stream;
}

std::istream& operator>>(std::istream& stream, FMIndex& index) {
    stream >> index._bwt;
    index.initialize();
    return stream;
}
