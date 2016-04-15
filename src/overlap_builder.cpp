#include "overlap_builder.h"
#include "asqg.h"
#include "constant.h"
#include "reads.h"
#include "sequence_process_framework.h"
#include "suffix_array.h"
#include "utils.h"

#include <iostream>
#include <memory>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.OverlapBuilder"));

//
// Flags indicating how a given read was aligned to the FM-index
//
struct AlignFlags {
public:
    AlignFlags() {
    }
    AlignFlags(bool qr, bool tr, bool qc) {
        _data.set(QUERYREV_BIT, qr);
        _data.set(TARGETREV_BIT, tr);
        _data.set(QUERYCOMP_BIT, qc);
    }
private:
    friend std::ostream& operator<<(std::ostream& stream, const AlignFlags& af);
    friend std::istream& operator>>(std::istream& stream, AlignFlags& af);

    static const size_t QUERYREV_BIT  = 0;
    static const size_t TARGETREV_BIT = 1;
    static const size_t QUERYCOMP_BIT = 2;
    std::bitset< 3 > _data;
};

static const AlignFlags kSuffixPrefixAF(false, false, false);
static const AlignFlags kPrefixPrefixAF(false, true, true);
static const AlignFlags kSuffixSuffixAF(true, false, true);
static const AlignFlags kPrefixSuffixAF(true, true, false);

std::ostream& operator<<(std::ostream& stream, const AlignFlags& af) {
    stream << af._data;
    return stream;
}

std::istream& operator>>(std::istream& stream, AlignFlags& af) {
    stream >> af._data;
    return stream;
}

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
        assert(i < SIZEOF_ARRAY(_intervals));
        return _intervals[i];
    }
    const FMIndex::Interval& operator[](size_t i) const {
        assert(i < SIZEOF_ARRAY(_intervals));
        return _intervals[i];
    }

    void init(char c, const FMIndex* index, const FMIndex* rindex) {
        _intervals[0].init(c,  index);
        _intervals[1].init(c, rindex);
    }
    void updateL(char c, const FMIndex* index) {
        // Update the left index using the difference between the AlphaCounts in the reverse table
        DNAAlphabet::AlphaCount64 l = index->getOcc(_intervals[0].lower - 1);
        DNAAlphabet::AlphaCount64 u = index->getOcc(_intervals[0].upper);
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
    OverlapBlock() : length(0) {
    }
    OverlapBlock(const IntervalPair& probe, const IntervalPair& ranges, size_t length, const AlignFlags& af) : probe(probe), ranges(ranges), length(length), af(af) {
    }

    Overlap overlap(const ReadInfo& query, const ReadInfo& target) const {
        return Overlap(
                query.name, 
                query.length - length,
                query.length - 1, 
                query.length, 
                target.name, 
                0, 
                length - 1, 
                target.length, 
                false, 
                0
                );
    }

    friend std::ostream& operator<<(std::ostream& stream, const OverlapBlock& block);
    friend std::istream& operator>>(std::istream& stream, OverlapBlock& block);

    IntervalPair probe;
    IntervalPair ranges;
    size_t length;
    AlignFlags af;
};

std::ostream& operator<<(std::ostream& stream, const OverlapBlock& block) {
    stream << block.probe << ' ' << block.ranges << ' ' << block.length << ' ' << block.af;
    return stream;
}

std::istream& operator>>(std::istream& stream, OverlapBlock& block) {
    stream >> block.probe >> block.ranges >> block.length >> block.af;
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

class Hits2ASQGConverter {
public:
    Hits2ASQGConverter(const SuffixArray& sa, const SuffixArray& rsa, DNASeqReader& reader) : _sa(sa), _rsa(rsa) {
        reader.reset();

        DNASeq read;
        while (reader.read(read)) {
            _readinfo.push_back(ReadInfo(read.name, read.seq.length()));
        }
    }

    bool convert(const std::string& hits, std::ostream& asqg) const {
        std::shared_ptr< std::streambuf > buf(ASQG::ifstreambuf(hits));
        if (!buf) {
            LOG4CXX_ERROR(logger, boost::format("failed to read hits %s") % hits);
            return false;
        }
        std::istream stream(buf.get());
        return convert(stream, asqg);
    }

    bool convert(std::istream& hits, std::ostream& asqg) const {
        while (hits) {
            // Read the overlap block for a read
            try {
                size_t count = 0, idx = 0;
                bool substring = false;
                hits >> idx >> substring >> count;

                for (size_t i = 0; i < count; ++i) {
                    OverlapBlock block;
                    hits >> block;
                    
                    // Iterate thru the range and write the overlaps
                    for (size_t j = block.probe[0].lower; j <= block.probe[0].upper; ++j) {
                        const ReadInfo& query = _readinfo[idx];

                        const ReadInfo& target = _readinfo[_sa[j].i];
                        if (query.name != target.name) {
                            Overlap o = block.overlap(query, target);

                            // The alignment logic above has potential to produce duplicate alignments
                            // To avoid this, we skip overlaps where the id of the first coord is lexo. lower than
                            // the second or the match is containment and the query is reversed (containments can be
                            // output up to 4 times total).
                            if (o.id[0] < o.id[1]) {
                                ASQG::EdgeRecord recod(o);
                                asqg << recod << '\n';
                            }
                        }
                    }
                }
            } catch (...) {
                return false;
            }
        }
        return true;
    }

private:
    const SuffixArray& _sa;
    const SuffixArray& _rsa;
    ReadInfoList _readinfo;
};

bool OverlapBuilder::build(DNASeqReader& reader, size_t minOverlap, std::ostream& output, size_t threads, size_t* processed) const {
    std::vector< std::string > hits;

    // Build and write the ASQG header
    {
        ASQG::HeaderRecord record;
        record.overlap(minOverlap);
        record.containment(1);
        output << record << '\n';
    }

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
    {
        std::shared_ptr< SuffixArray > sa(SuffixArray::load(_prefix + SAI_EXT)), rsa(SuffixArray::load(_prefix + RSAI_EXT));
        if (!sa || !rsa) {
            LOG4CXX_ERROR(logger, boost::format("failed to load suffix array index %s") % _prefix);
            return false;
        }

        Hits2ASQGConverter converter(*sa, *rsa, reader);
        BOOST_FOREACH(const std::string& filename, hits) {
            LOG4CXX_INFO(logger, boost::format("parsing file %s") % filename);
            if (!converter.convert(filename, output)) {
                LOG4CXX_ERROR(logger, boost::format("failed to convert hits to asqg %s") % filename);
                return false;
            }
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

class OverlapBlockFinder {
public:
    OverlapBlockFinder(const FMIndex* fmi, const FMIndex* rfmi, size_t minOverlap) : _fmi(fmi), _rfmi(rfmi), _minOverlap(minOverlap) {
    }

    // Calculate the ranges in FMI that contain a prefix of at least minOverlap basepairs that
    // overlaps with a suffix of w. The ranges are added to the pOBList
    void find(const std::string& seq, const AlignFlags& af, OverlapBlockList* overlaps, OverlapBlockList* contains, OverlapResult* result) const {
        assert(!seq.empty());
        // The algorithm is as follows:
        // We perform a backwards search using the FM-index for the string w.
        // As we perform the search we collect the intervals 
        // of the significant prefixes (len >= minOverlap) that overlap w.
        IntervalPair ranges;
        size_t l = seq.length();
        ranges.init(seq[l - 1], _fmi, _rfmi);

        // Collect the OverlapBlocks
        for (size_t i = l - 1; i > 1; --i) {
            // Compare the range of the suffix seq[i, l]
            ranges.updateL(seq[i - 1], _fmi);

            if (l - i + 1 >= _minOverlap) {
                // Calculate which of the prefixes that match w[i, l] are terminal
                // These are the proper prefixes (they are the start of a read)
                IntervalPair probe = ranges;
                probe.updateL('$', _fmi);

                // The probe interval contains the range of proper prefixes
                if (probe[1].valid()) {
                    //std::cout << "overlaps\t" << OverlapBlock(probe, ranges, l - i + 1, af) << std::endl;
                    overlaps->push_back(OverlapBlock(probe, ranges, l - i + 1, af));
                }
            }
        }

        // Determine if this sequence is contained and should not be processed further
        ranges.updateL(seq[0], _fmi);

        // Ranges now holds the interval for the full-length read
        // To handle containments, we output the overlapBlock to the final overlap block list
        // and it will be processed later
        // Two possible containment cases:
        // 1) This read is a substring of some other read
        // 2) This read is identical to some other read

        // Case 1 is indicated by the existance of a non-$ left or right hand extension
        // In this case we return no alignments for the string
        DNAAlphabet::AlphaCount64 lext =  _fmi->getOcc(ranges[0].upper) -  _fmi->getOcc(ranges[0].lower - 1);
        DNAAlphabet::AlphaCount64 rext = _rfmi->getOcc(ranges[1].upper) - _rfmi->getOcc(ranges[1].lower - 1);
        if (lext.hasDNA() || rext.hasDNA()) {
            result->substring = true;
        } else {
            IntervalPair probe = ranges;
            probe.updateL('$', _fmi);
            if (probe.valid()) {
                // terminate the contained block and add it to the contained list
                probe.updateR('$', _rfmi);
                assert(probe.valid());
                //std::cout << "contains\t" << OverlapBlock(probe, ranges, l, af) << std::endl;
                contains->push_back(OverlapBlock(probe, ranges, l, af));
            }
        }
    }

private:
    const FMIndex* _fmi;
    const FMIndex* _rfmi;
    size_t _minOverlap;
};

OverlapResult OverlapBuilder::overlap(const DNASeq& read, size_t minOverlap, OverlapBlockList* blocks) const {
    // The complete set of overlap blocks are collected in workinglist
    // The filtered set (containing only irreducible overlaps) are placed into blocks
    // by calculateIrreducibleHits
    const std::string& seq = read.seq;
    OverlapResult result;

    OverlapBlockFinder finder(_fmi, _rfmi, minOverlap), rfinder(_rfmi, _fmi, minOverlap);

    // Match the suffix of seq to prefixes
    finder.find(seq, kSuffixPrefixAF, blocks, blocks, &result);
    //rfinder.find(make_complement_dna(seq), kPrefixPrefixAF, blocks, blocks, &result);

    // Match the prefix of seq to suffixes
    //finder.find(make_reverse_complement_dna(seq), kSuffixSuffixAF, blocks, blocks, &result);
    //rfinder.find(make_reverse_dna(seq), kPrefixPrefixAF, blocks, blocks, &result);

    return result;
}
