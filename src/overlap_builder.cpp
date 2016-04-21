#include "overlap_builder.h"
#include "asqg.h"
#include "constant.h"
#include "reads.h"
#include "sequence_process_framework.h"
#include "suffix_array.h"
#include "utils.h"

#include <bitset>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
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
    bool test(size_t pos) const {
        return _data.test(pos);
    }

    static const size_t QUERYREV_BIT  = 0;
    static const size_t TARGETREV_BIT = 1;
    static const size_t QUERYCOMP_BIT = 2;
private:
    friend std::ostream& operator<<(std::ostream& stream, const AlignFlags& af);
    friend std::istream& operator>>(std::istream& stream, AlignFlags& af);

    std::bitset< 3 > _data;
};

static const AlignFlags kSuffixPrefixAF(false, false, false);
static const AlignFlags kSuffixSuffixAF(false, true, true);
static const AlignFlags kPrefixPrefixAF(true, false, true);
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
        // Update the left index using the difference between the AlphaCounts in the reverse table
        DNAAlphabet::AlphaCount64 l = index->getOcc(_intervals[1].lower - 1);
        DNAAlphabet::AlphaCount64 u = index->getOcc(_intervals[1].upper);
        updateR(c, index, l, u);
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
    void updateR(char c, const FMIndex* index, const DNAAlphabet::AlphaCount64& l, const DNAAlphabet::AlphaCount64& u) {
        DNAAlphabet::AlphaCount64 diff = u - l;

        _intervals[0].lower = _intervals[0].lower + std::accumulate(&diff[0], &diff[0] + DNAAlphabet::torank(c), 0);
        _intervals[0].upper = _intervals[0].lower + diff[DNAAlphabet::torank(c)] - 1;

        // Update the right index directly
        size_t pb = index->getPC(c);
        _intervals[1].lower = pb + l[DNAAlphabet::torank(c)];
        _intervals[1].upper = pb + u[DNAAlphabet::torank(c)] - 1;
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

    const FMIndex* index(const FMIndex* index, const FMIndex* rindex) const {
        return af.test(AlignFlags::TARGETREV_BIT) ? rindex : index;
    }

    DNAAlphabet::AlphaCount64 ext(const FMIndex* fmi, const FMIndex* rfmi) const {
        DNAAlphabet::AlphaCount64 count = probe[1].ext(index(fmi, rfmi));
        if (af.test(AlignFlags::QUERYCOMP_BIT)) {
            count.complement();
        }
        return count;
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
// Hit
//
struct Hit {
    Hit(size_t idx=-1, bool substring=false) : idx(idx), substring(substring) {
    }
    Hit(size_t idx, bool substring, const OverlapBlockList& blocks) : idx(idx), substring(substring), blocks(blocks) {
    }

    size_t idx;
    bool substring;
    OverlapBlockList blocks;
};

std::ostream& operator<<(std::ostream& stream, const Hit& hit) {
    // Write the header info
    stream << hit.idx << ' ' << hit.substring << ' ' << hit.blocks.size() << ' ';
    BOOST_FOREACH(const OverlapBlock& block, hit.blocks) {
        stream << block << ' ';
    }
    return stream;
}

std::istream& operator>>(std::istream& stream, Hit& hit) {
    size_t count = 0;

    stream >> hit.idx >> hit.substring >> count;
    for (size_t i = 0; i < count; ++i) {
        OverlapBlock block;
        stream >> block;
        hit.blocks.push_back(block);
    }

    return stream;
}

//
// OverlapProcess
//
class OverlapProcess {
public:
    OverlapProcess(const OverlapBuilder* builder, size_t minOverlap, std::ostream& stream) : _builder(builder), _minOverlap(minOverlap), _stream(stream) {
    }

    OverlapResult process(const SequenceProcessFramework::SequenceWorkItem& workItem) {
        Hit hit(workItem.idx);
        OverlapResult result = _builder->overlap(workItem.read, _minOverlap, &hit.blocks);
        hit.substring = result.substring;

        //
        // Write overlap blocks out to a file
        //
        _stream << hit << '\n';

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
        std::string line;
        while (std::getline(hits, line)) {
            boost::algorithm::trim(line);
            if (!line.empty()) {
                // Read the overlap block for a read
                Hit hit;
                std::stringstream ss(line);
                ss >> hit;

                BOOST_FOREACH(const OverlapBlock& block, hit.blocks) {
                    // Iterate thru the range and write the overlaps
                    for (size_t j = block.probe[0].lower; j <= block.probe[0].upper; ++j) {
                        const ReadInfo& query = _readinfo[hit.idx];

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

//
// DuplicateRemoveProcess
//
class DuplicateRemoveProcess {
public:
    DuplicateRemoveProcess(const OverlapBuilder* builder, std::ostream& stream) : _builder(builder), _stream(stream) {
    }

    OverlapResult process(const SequenceProcessFramework::SequenceWorkItem& workItem) {
        Hit hit(workItem.idx);
        OverlapResult result = _builder->duplicate(workItem.read, &hit.blocks);
        hit.substring = result.substring;
        _stream << workItem.read.name << '\t' << workItem.read.seq << '\t' << hit << '\n';

        return result;
    }

private:
    const OverlapBuilder* _builder;
    std::ostream& _stream;
};

//
// DuplicateRemovePostProcess
//
class DuplicateRemovePostProcess {
public:
    DuplicateRemovePostProcess() {
    }

    void process(const SequenceProcessFramework::SequenceWorkItem& workItem, const OverlapResult& result) {
    }
};

class Hits2FastaConverter {
public:
    Hits2FastaConverter(const SuffixArray& sa, const SuffixArray& rsa, DNASeqReader& reader) : _sa(sa), _rsa(rsa) {
        reader.reset();

        DNASeq read;
        while (reader.read(read)) {
            _readinfo.push_back(ReadInfo(read.name, read.seq.length()));
        }
    }

    bool convert(const std::string& hits, std::ostream& fasta) const {
        std::shared_ptr< std::streambuf > buf(ASQG::ifstreambuf(hits));
        if (!buf) {
            LOG4CXX_ERROR(logger, boost::format("failed to read hits %s") % hits);
            return false;
        }
        std::istream stream(buf.get());
        return convert(stream, fasta);
    }

    bool convert(std::istream& hits, std::ostream& fasta) const {
        std::string line;
        while (std::getline(hits, line)) {
            boost::algorithm::trim(line);
            if (!line.empty()) {
                // Read the overlap block for a read
                DNASeq item;
                Hit hit;
                std::stringstream ss(line);
                ss >> item.name >> item.seq >> hit;

                size_t numCopies = 0;
                BOOST_FOREACH(const OverlapBlock& block, hit.blocks) {
                    // Iterate thru the range and write the overlaps
                    assert(block.probe[0].lower <= block.probe[0].upper);
                    numCopies += block.probe[0].upper - block.probe[0].lower + 1;
                }

                bool isContained = false;
                if (hit.substring) {
                    isContained = true;
                }

                std::string meta = boost::str(boost::format("%s NumDuplicates=%d") % item.name % numCopies);
                if (isContained) {
                    // The read's index in the sequence data base
                    // is needed when removing it from the FM-index.
                    // In the output fasta, we set the reads ID to be the index
                    // and record its old id in the fasta header.
                } else {
                    item.name = boost::str(boost::format("%s %s") % item.name % meta);
                    fasta << item;
                }
            }
        }
        return true;
    }

private:
    const SuffixArray& _sa;
    const SuffixArray& _rsa;
    ReadInfoList _readinfo;
};

bool OverlapBuilder::rmdup(DNASeqReader& reader, std::ostream& output, size_t threads, size_t* processed) const {
    std::vector< std::string > hits;

    if (threads <= 1) { // single thread
        std::string hit = _prefix + RMDUP_EXT + HITS_EXT + GZIP_EXT;
        std::shared_ptr< std::streambuf > buf(ASQG::ofstreambuf(hit));
        if (!buf) {
            LOG4CXX_ERROR(logger, boost::format("failed to create hits %s") % hit);
            return false;
        }
        std::ostream stream(buf.get());
        DuplicateRemoveProcess proc(this, stream);
        DuplicateRemovePostProcess postproc;

        SequenceProcessFramework::SerialWorker<
            SequenceProcessFramework::SequenceWorkItem, 
            OverlapResult, 
            SequenceProcessFramework::WorkItemGenerator< SequenceProcessFramework::SequenceWorkItem >, 
            DuplicateRemoveProcess, 
            DuplicateRemovePostProcess
            > worker;
        size_t num = worker.run(reader, &proc, &postproc);
        if (processed != NULL) {
            *processed = num;
        }

        hits.push_back(hit);
    } else { // multi thread
        assert(false);
    }

    // Convert hits to fasta
    {
        std::shared_ptr< SuffixArray > sa(SuffixArray::load(_prefix + SAI_EXT)), rsa(SuffixArray::load(_prefix + RSAI_EXT));
        if (!sa || !rsa) {
            LOG4CXX_ERROR(logger, boost::format("failed to load suffix array index %s") % _prefix);
            return false;
        }

        Hits2FastaConverter converter(*sa, *rsa, reader);
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

bool OverlapBuilder::rmdup(const std::string& input, const std::string& output, size_t threads, size_t* processed) const {
    // DNASeqReader
    std::ifstream reads(input);
    std::shared_ptr< DNASeqReader > reader(DNASeqReaderFactory::create(reads));
    if (!reader) {
        LOG4CXX_ERROR(logger, boost::format("Failed to create DNASeqReader %s") % input);
        return false;
    }

    // FASTA
    std::ofstream fasta(output.c_str());
    if (!fasta) {
        LOG4CXX_ERROR(logger, boost::format("Failed to create FASTA %s") % output);
        return false;
    }

    // Build
    return rmdup(*reader, fasta, threads, processed);
}

class IrreducibleBlockListExtractor {
public:
    IrreducibleBlockListExtractor(const FMIndex* fmi, const FMIndex* rfmi) : _fmi(fmi), _rfmi(rfmi) {
    }

    bool extract(OverlapBlockList* blocks) {
        assert(blocks != NULL);

        // We store the overlap blocks in groups of blocks that have the same right-extension.
        // When a branch is found, the groups are split based on the extension
        typedef std::list< OverlapBlockList > BlockGroups;

        BlockGroups groups = boost::assign::list_of(*blocks);
        while (!groups.empty()) {
            // Perform one extenion round for each group.
            // If the top-level block has ended, push the result
            // to the final list and remove the group from processing
            BlockGroups incomings; // Branched blocks are placed here
            for (BlockGroups::iterator i = groups.begin(); i != groups.end(); ++i) {
                OverlapBlockList& blocklist = *i;
                bool eraseGroup = true;

                // Count the extensions in the top level (longest) blocks first
                DNAAlphabet::AlphaCount64 exts;
                size_t topLength = blocklist.front().length;
                for (OverlapBlockList::const_iterator j = blocklist.begin(); j != blocklist.end() && j->length == topLength; ++j) {
                    exts += j->ext(_fmi, _rfmi);
                }

                // Three cases:
                // 1) The top level block has ended as it contains the extension $. Output TLB and end.
                // 2) There is a singular unique extension base for all the blocks. Update the blocks and continue.
                // 3) There are multiple extension bases, split the block group and continue.
                // If some block other than the TLB ended, it must be contained within the TLB and it is not output
                // or considered further. 
                // Likewise if multiple distinct strings in the TLB ended, we only output the top one. The rest
                // must have the same sequence as the top one and are hence considered to be contained with the top element.
                if (exts[DNAAlphabet::torank('$')] > 0) {
                    // An irreducible overlap has been found. It is possible that there are two top level blocks
                    // (one in the forward and reverse direction). Since we can't decide which one
                    // contains the other at this point, we output hits to both. Under a fixed 
                    // length string assumption one will be contained within the other and removed later.
                    for (OverlapBlockList::const_iterator j = blocklist.begin(); j != blocklist.end() && j->length == topLength; ++j) {
                        DNAAlphabet::AlphaCount64 test = j->ext(_fmi, _rfmi);
                        if (test[DNAAlphabet::torank('$')] == 0) {
                            LOG4CXX_ERROR(logger, "substring read found during overlap computation.");
                            LOG4CXX_ERROR(logger, "Please run rmdup before  overlap.");
                            return false;
                        }

                        // Perform the final right-update to make the block terminal
                        OverlapBlock branched = *j;
                        branched.probe.updateR('$', branched.index(_fmi, _rfmi));
                        blocks->push_back(branched);

                        LOG4CXX_DEBUG(logger, boost::format("TLB of length %d has ended") % branched.length);
                    }
                } else {
                    // Count the extension for the rest of the blocks
                    for (OverlapBlockList::const_iterator j = blocklist.begin(); j != blocklist.end() && j->length != topLength; ++j) {
                        exts += j->ext(_fmi, _rfmi);
                    }

                    // If only one of the DNA characters has a non-zero count
                    if (std::count_if(&exts[0], &exts[0] + exts.size(), std::bind2nd(std::greater< uint64_t >(), 0)) == 1) {
                        // Update all the blocks using the unique extension character
                        // This character is in the canonical representation wrt to the query
                        char c = DNAAlphabet::tochar(
                                std::find_if(&exts[0], &exts[0] + exts.size(), std::bind2nd(std::greater< uint64_t >(), 0)) - &exts[0]
                                );
                        updateR(c, &blocklist);

                        // Set the flag to erase this group, it is finished
                        eraseGroup = false;
                    } else {
                        for (size_t j = 0; j < exts.size(); ++j) {
                            if (exts[j] > 0) {
                                OverlapBlockList branched = blocklist;
                                updateR(DNAAlphabet::tochar(j), &branched);
                                incomings.push_back(branched);
                            }
                        }
                    }
                }

                if (eraseGroup) {
                    i = groups.erase(i);
                } else {
                    ++i;
                }
            }

            // Splice in the newly branched blocks, if any
            std::copy(incomings.begin(), incomings.end(), std::back_inserter(groups));
        }

        return true;
    }
private:
    void updateR(char c, OverlapBlockList* blocks) {
        assert(blocks != NULL);
        OverlapBlockList::iterator i = blocks->begin();
        while (i != blocks->end()) {
            char b = i->af.test(AlignFlags::QUERYCOMP_BIT) ? make_complement_dna(c) : c;
            i->probe.updateR(b, i->index(_fmi, _rfmi));

            // remove the block from the list if its no longer valid
            if (!i->probe.valid()) {
                i = blocks->erase(i);
            } else {
                ++i;
            }
        }
    }

    const FMIndex* _fmi;
    const FMIndex* _rfmi;
};

class OverlapBlockFinder {
public:
    OverlapBlockFinder(const FMIndex* fmi, const FMIndex* rfmi, size_t minOverlap) : _fmi(fmi), _rfmi(rfmi), _minOverlap(minOverlap) {
    }

    // Calculate the ranges in FMI that contain a prefix of at least minOverlap basepairs that
    // overlaps with a suffix of w. The ranges are added to the pOBList
    void find(const std::string& seq, const AlignFlags& af, OverlapBlockList* overlaps, OverlapBlockList* contains, OverlapResult* result) const {
        assert(!seq.empty());
        // The algorithm is as follows:
        // We perform a backwards search using the FM-index for the string seq.
        // As we perform the search we collect the intervals 
        // of the significant prefixes (len >= minOverlap) that overlap seq.
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
                    if (overlaps != NULL) {
                        overlaps->push_back(OverlapBlock(probe, ranges, l - i + 1, af));
                    }
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
        // DNAAlphabet::AlphaCount64 lext =  _fmi->getOcc(ranges[0].upper) -  _fmi->getOcc(ranges[0].lower - 1);
        DNAAlphabet::AlphaCount64 lext =  ranges[0].ext(_fmi);
        // DNAAlphabet::AlphaCount64 rext = _rfmi->getOcc(ranges[1].upper) - _rfmi->getOcc(ranges[1].lower - 1);
        DNAAlphabet::AlphaCount64 rext =  ranges[1].ext(_rfmi);
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
                if (contains != NULL) {
                    contains->push_back(OverlapBlock(probe, ranges, l, af));
                }
            }
        }
    }

private:
    const FMIndex* _fmi;
    const FMIndex* _rfmi;
    size_t _minOverlap;
};

class SubMaximalBlockFilter {
public:
    SubMaximalBlockFilter(const FMIndex* fmi, const FMIndex* rfmi) : _fmi(fmi), _rfmi(rfmi) {
    }

    void filter(OverlapBlockList* blocks) {
        assert(blocks != NULL);
        // This algorithm removes any sub-maximal OverlapBlocks from pList
        // The list is sorted by the left coordinate and iterated through
        // if two adjacent blocks overlap they are split into maximal contiguous regions
        // with resolveOverlap. The resulting list is merged back into pList. This process
        // is repeated until each block in pList is a unique range
        // The bookkeeping in the intersecting case could be more efficient 
        // but the vast vast majority of the cases will not have overlapping 
        // blocks.
        if (!blocks->empty()) {
            IntervalLeftSorter sorter;
            blocks->sort(sorter);
            OverlapBlockList::iterator prev = blocks->begin();
            OverlapBlockList::iterator curr = std::next(prev);
            while (curr != blocks->end()) {
                // Check if prev and curr overlaps
                if (Interval::isIntersecting(prev->probe[0].lower, prev->probe[1].upper, curr->probe[0].lower, curr->probe[0].upper)) {

                    // Merge the new elements in and start back from the beginning of the list
                    OverlapBlockList resolved;
                    resolve(*prev, *curr, &resolved);

                    blocks->erase(curr);
                    blocks->erase(prev);
                    blocks->merge(resolved, sorter);

                    prev = blocks->begin();
                } else {
                    ++prev;
                }
                curr = std::next(prev);
            }
        }
    }
private:
    void resolve(const OverlapBlock& x, const OverlapBlock& y, OverlapBlockList* resolve) {
    }

    class IntervalLeftSorter {
    public:
        bool operator()(const OverlapBlock& x, const OverlapBlock& y) const {
            return x.probe[0].lower < y.probe[0].lower;
        }
    };
    const FMIndex* _fmi;
    const FMIndex* _rfmi;
};

OverlapResult OverlapBuilder::overlap(const DNASeq& read, size_t minOverlap, OverlapBlockList* blocks) const {
    // The complete set of overlap blocks are collected in workinglist
    // The filtered set (containing only irreducible overlaps) are placed into blocks
    // by calculateIrreducibleHits
    OverlapResult result;

    const std::string& seq = read.seq;
    OverlapBlockFinder finder(_fmi, _rfmi, minOverlap), rfinder(_rfmi, _fmi, minOverlap);

    // Match the suffix of seq to prefixes
    finder.find(seq, kSuffixPrefixAF, blocks, blocks, &result);
    rfinder.find(make_complement_dna(seq), kSuffixSuffixAF, blocks, blocks, &result);

    // Match the prefix of seq to suffixes
    rfinder.find(make_reverse_dna(seq), kPrefixSuffixAF, blocks, blocks, &result);
    finder.find(make_reverse_complement_dna(seq), kPrefixPrefixAF, blocks, blocks, &result);

    SubMaximalBlockFilter filter(_fmi, _rfmi);
    filter.filter(blocks);
    
    IrreducibleBlockListExtractor extractor(_fmi, _rfmi);
    extractor.extract(blocks);

    return result;
}

OverlapResult OverlapBuilder::duplicate(const DNASeq& read, OverlapBlockList* blocks) const {
    OverlapResult result;

    const std::string& seq = read.seq;
    size_t minOverlap = seq.length();
    OverlapBlockFinder finder(_fmi, _rfmi, minOverlap), rfinder(_rfmi, _fmi, minOverlap);

    finder.find(seq, kSuffixPrefixAF, NULL, blocks, &result);
    rfinder.find(make_complement_dna(seq), kSuffixSuffixAF, NULL, blocks, &result);

    return result;
}

