#include "overlap_builder.h"
#include "asqg.h"
#include "constant.h"
#include "sequence_process_framework.h"

#include <iostream>
#include <memory>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.OverlapBuilder"));

//
// A pair of intervals used for bidirectional searching a FM-index/reverse FM-index
//
class IntervalPair {
public:
    IntervalPair() {
    }
    bool valid() const {
        for (size_t i = 0; i < SIZEOF_ARRAY(_intervals); ++i) {
            if (!_intervals[i].valid()) {
                return false;
            }
        }
        return true;
    }
    FMIndex::Interval& operator[](size_t i) {
        return _intervals[i];
    }
    const FMIndex::Interval& operator[](size_t i) const {
        return _intervals[i];
    }

    void init(char c, const FMIndex* index, const FMIndex* rindex) {
        _intervals[0].init(c,  index);
        _intervals[1].init(c, rindex);
    }
    void updateL(char c, const FMIndex* index) {
        // Update the left index using the difference between the AlphaCounts in the reverse table
        DNAAlphabet::AlphaCount64 l = index->getOcc(_intervals[1].lower - 1);
        DNAAlphabet::AlphaCount64 u = index->getOcc(_intervals[1].upper);
        updateL(c, index, l, u);
    } 
    void updateR(char c, const FMIndex* index) {
    }
private:
    friend std::ostream& operator<<(std::ostream& stream, const IntervalPair& pair);
    friend std::istream& operator>>(std::istream& stream, IntervalPair& pair);

    void updateL(char c, const FMIndex* index, const DNAAlphabet::AlphaCount64& l, const DNAAlphabet::AlphaCount64& u) {
        DNAAlphabet::AlphaCount64 diff = u - l;
        // Update the left index using the difference between the AlphaCounts in the reverse table
        _intervals[1].lower = _intervals[1].lower + std::accumulate(&diff[0], &diff[0] + DNAAlphabet::torank(c), 0);
        _intervals[1].upper = _intervals[1].lower + diff[DNAAlphabet::torank(c)] - 1;

        // Update the left index directly
        size_t pb = index->getPC(c);
        _intervals[0].lower = pb + l[DNAAlphabet::torank(c)];
        _intervals[0].upper = pb + u[DNAAlphabet::torank(c)] - 1;
    }

    FMIndex::Interval _intervals[2];
};

std::ostream& operator<<(std::ostream& stream, const IntervalPair& pair) {
    stream << pair._intervals[0] << ' ' << pair._intervals[1];
    return stream;
}

std::istream& operator>>(std::istream& stream, IntervalPair& pair) {
    stream >> pair._intervals[0] >> pair._intervals[1];
    return stream;
}

//
// OverlapBlock
//
struct OverlapBlock {
    OverlapBlock(const IntervalPair& probe, const IntervalPair& ranges, size_t length) : probe(probe), ranges(ranges), length(length) {
    }
    friend std::ostream& operator<<(std::ostream& stream, const OverlapBlock& block);
    friend std::istream& operator>>(std::istream& stream, OverlapBlock& block);

    IntervalPair ranges;
    IntervalPair probe;
    size_t length;
};

std::ostream& operator<<(std::ostream& stream, const OverlapBlock& block) {
    stream << block.ranges << ' ' << block.probe << ' ' << block.length;
    return stream;
}

std::istream& operator>>(std::istream& stream, OverlapBlock& block) {
    stream >> block.ranges >> block.probe >> block.length;
    return stream;
}

//
// OverlapResult
//
struct OverlapResult {
    OverlapResult() : substring(false), aborted(false) {
    }
    bool substring;
    bool aborted;
};

//
// OverlapProcess
//
class OverlapProcess {
public:
    OverlapProcess(const OverlapBuilder* builder, size_t minOverlap, std::ostream& stream) : _builder(builder), _minOverlap(minOverlap), _stream(stream) {
    }

    OverlapResult process(const SequenceProcessFramework::SequenceWorkItem& workItem) {
        OverlapBlockList blocks;
        OverlapResult result = _builder->overlap(workItem.read, _minOverlap, &blocks);

        //
        // Write overlap blocks out to a file
        //
        // Write the header info
        _stream << workItem.idx << ' ' << result.substring << ' ' << blocks.size() << ' ';
        BOOST_FOREACH(const OverlapBlock& block, blocks) {
            _stream << block << ' ';
        }
        _stream << '\n';

        return result;
    }

private:
    const OverlapBuilder* _builder;
    size_t _minOverlap;
    std::ostream& _stream;
};

//
// OverlapPostProcess
//
class OverlapPostProcess {
public:
    OverlapPostProcess(std::ostream& stream) : _stream(stream) {
    }

    void process(const SequenceProcessFramework::SequenceWorkItem& workItem, const OverlapResult& result) {
        ASQG::VertexRecord record(workItem.read.name, workItem.read.seq);
        record.substring = result.substring ? 1 : 0;
        _stream << record << '\n';
    }

private:
    std::ostream& _stream;
};

bool OverlapBuilder::build(DNASeqReader& reader, size_t minOverlap, std::ostream& output, size_t threads, size_t* processed) const {
    std::vector< std::string > hits;

    if (threads <= 1) { // single thread
        std::string hit = _prefix + HITS_EXT + GZIP_EXT;
        std::shared_ptr< std::streambuf > buf(ASQG::ofstreambuf(hit));
        if (!buf) {
            LOG4CXX_ERROR(logger, boost::format("failed to create hits %s") % hit);
            return false;
        }
        std::ostream stream(buf.get());
        OverlapProcess proc(this, minOverlap, stream);
        OverlapPostProcess postproc(output);

        SequenceProcessFramework::SerialWorker<
            SequenceProcessFramework::SequenceWorkItem, 
            OverlapResult, 
            SequenceProcessFramework::WorkItemGenerator< SequenceProcessFramework::SequenceWorkItem >, 
            OverlapProcess, 
            OverlapPostProcess
            > worker;
        size_t num = worker.run(reader, &proc, &postproc);
        if (processed != NULL) {
            *processed = num;
        }

        hits.push_back(hit);
    } else { // multi thread
        assert(false);
    }

    // Convert hits to ASQG
    BOOST_FOREACH(const std::string& filename, hits) {
        std::shared_ptr< std::streambuf > buf(ASQG::ifstreambuf(filename));
        if (!buf) {
            LOG4CXX_ERROR(logger, boost::format("failed to read hits %s") % filename);
            return false;
        }
        std::istream stream(buf.get());
        if (!hits2asqg(stream, output)) {
            LOG4CXX_ERROR(logger, boost::format("failed to convert hits to asqg %s") % filename);
            return false;
        }
    }

    return true;
}

bool OverlapBuilder::build(const std::string& input, size_t minOverlap, const std::string& output, size_t threads, size_t* processed) const {
    // DNASeqReader
    std::ifstream reads(input);
    std::shared_ptr< DNASeqReader > reader(DNASeqReaderFactory::create(reads));
    if (!reader) {
        LOG4CXX_ERROR(logger, boost::format("Failed to create DNASeqReader %s") % input);
        return false;
    }

    // ASQG
    std::shared_ptr< std::streambuf > buf(ASQG::ofstreambuf(output));
    if (!buf) {
        LOG4CXX_ERROR(logger, boost::format("Failed to create ASQG %s") % output);
        return false;
    }
    std::ostream asqg(buf.get());

    // Build
    return build(*reader, minOverlap, asqg, threads, processed);
}

OverlapResult OverlapBuilder::overlap(const DNASeq& read, size_t minOverlap, OverlapBlockList* blocks) const {
    // The complete set of overlap blocks are collected in obWorkingList
    // The filtered set (containing only irreducible overlaps) are placed into pOBOut
    // by calculateIrreducibleHits
    const std::string& seq = read.seq;

    return overlap(seq, _fmi, _rfmi, minOverlap, NULL, NULL);
    //overlap(make_complement_dna(seq), _fmi, _rfmi, minOverlap, NULL, NULL);
}

// Calculate the ranges in pBWT that contain a prefix of at least minOverlap basepairs that
// overlaps with a suffix of w. The ranges are added to the pOBList
OverlapResult OverlapBuilder::overlap(const std::string& seq, const FMIndex* fmi, const FMIndex* rfmi, size_t minOverlap, OverlapBlockList* overlaps, OverlapBlockList* contains) const {
    OverlapResult result;
    // The algorithm is as follows:
    // We perform a backwards search using the FM-index for the string w.
    // As we perform the search we collect the intervals 
    // of the significant prefixes (len >= minOverlap) that overlap w.
    IntervalPair ranges;
    size_t l = seq.length();
    ranges.init(seq[l - 1], fmi, rfmi);

    // Collect the OverlapBlocks
    for (size_t i = l - 1; i > 1; --i) {
        // Compare the range of the suffix seq[i, l]
        ranges.updateL(seq[i - 1], fmi);

        if (l - i >= minOverlap) {
            // Calculate which of the prefixes that match w[i, l] are terminal
            // These are the proper prefixes (they are the start of a read)
            IntervalPair probe = ranges;
            probe.updateL('$', fmi);

            // The probe interval contains the range of proper prefixes
            if (probe[1].valid()) {
                overlaps->push_back(OverlapBlock(probe, ranges, l - i));
            }
        }
    }

    // Determine if this sequence is contained and should not be processed further
    ranges.updateL(seq[0], fmi);

    // Ranges now holds the interval for the full-length read
    // To handle containments, we output the overlapBlock to the final overlap block list
    // and it will be processed later
    // Two possible containment cases:
    // 1) This read is a substring of some other read
    // 2) This read is identical to some other read

    // Case 1 is indicated by the existance of a non-$ left or right hand extension
    // In this case we return no alignments for the string
    DNAAlphabet::AlphaCount64 lext = fmi->getOcc(ranges[0].upper) - fmi->getOcc(ranges[0].lower - 1);
    DNAAlphabet::AlphaCount64 rext = fmi->getOcc(ranges[1].upper) - fmi->getOcc(ranges[1].lower - 1);
    if (lext.hasDNA() || rext.hasDNA()) {
    } else {
        IntervalPair probe = ranges;
        probe.updateL('$', fmi);
        if (probe.valid()) {
            // terminate the contained block and add it to the contained list
            probe.updateR('$', rfmi);
            assert(probe.valid());
            contains->push_back(OverlapBlock(probe, ranges, l));
        }
    }

    return result;
}

bool OverlapBuilder::hits2asqg(std::istream& input, std::ostream& output) const {
    return true;
}

