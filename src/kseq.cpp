#include "kseq.h"

#include <fstream>
#include <map>
#include <memory>
#include <numeric>
#include <unordered_map>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.DNASeq"));

char make_complement_dna(char c) {
    static std::map< char, char > mapping = {
        {'A', 'T'},
        {'C', 'G'},
        {'G', 'C'},
        {'T', 'A'},
        {'N', 'N'},
        };
    return mapping[c];
}

void make_complement_dna(std::string& sequence) {
    static std::map< char, char > mapping = {
        {'A', 'T'},
        {'C', 'G'},
        {'G', 'C'},
        {'T', 'A'},
        {'N', 'N'},
        };

    size_t N = sequence.length();
    for (size_t i = 0; i < N; ++i) {
        sequence[i] = mapping[sequence[i]];
    }
}

std::string make_complement_dna_copy(const std::string& sequence) {
    std::string complement = sequence;
    make_complement_dna(complement);
    return complement;
}

void make_reverse_dna(std::string& sequence) {
    std::reverse(sequence.begin(), sequence.end());
}

std::string make_reverse_dna_copy(const std::string& sequence) {
    std::string reverse = sequence;
    make_reverse_dna(reverse);
    return reverse;
}

void make_reverse_complement_dna(std::string& sequence) {
    make_complement_dna(sequence);
    make_reverse_dna(sequence);
}

std::string make_reverse_complement_dna_copy(const std::string& sequence) {
    std::string complement = make_complement_dna_copy(sequence);
    make_reverse_dna(complement);
    return complement;
}

void DNASeq::make_complement() {
    make_complement_dna(seq);
}

void DNASeq::make_reverse() {
    std::reverse(seq.begin(), seq.end());
    if (!quality.empty()) {
        std::reverse(quality.begin(), quality.end());
    }
}

void DNASeq::make_reverse_complement() {
    make_complement();
    make_reverse();
}

std::ostream& operator << (std::ostream& os, const DNASeq& seq) {
    if (seq.quality.empty()) {
        os << '>' << seq.name << '\n';
        os << seq.seq << '\n';
    } else {
        os << '@' << seq.name << '\n';
        os << seq.seq << '\n';
        os << '+' << '\n';
        os << seq.quality << '\n';
    }
    return os;
}

DNASeqReader* DNASeqReaderFactory::create(std::istream& stream, const void* extra) {
    if (stream) {
        int c = stream.peek();
        if (c == '@') {
            return new FASTQReader(stream, extra);
        } else if (c == '>') {
            return new FASTAReader(stream, extra);
        }
    }
    return NULL;
}

bool FASTQReader::read(DNASeq& sequence) {
    enum {
        kName, 
        kSequence, 
        kName2, 
        kQuality, 
    };

    if (_stream) {
        int state = kName;
        std::string buf;

        while (std::getline(_stream, buf)) {
            boost::algorithm::trim(buf);
            if (buf.empty()) continue;
            if (state == kName) {
                if (boost::algorithm::starts_with(buf, "@")) {
                    sequence.name = buf.substr(1);
                    state = kSequence;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fastq=>invalid line for sequence name: %s") % buf);
                    return false;
                }
            } else if (state == kSequence) {
                sequence.seq = buf;
                state = kName2;
            } else if (state == kName2) {
                if (boost::algorithm::starts_with(buf, "+") && (buf.length() == 1 || boost::algorithm::ends_with(buf, sequence.name))) {
                    state = kQuality;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fastq=>names aren't equal: %s") % buf);
                    return false;
                }
            } else if (state == kQuality) {
                if (buf.length() == sequence.seq.length()) {
                    sequence.quality = buf;
                    // name
                    size_t i = sequence.name.find_first_of(" \t");
                    if (i != std::string::npos) {
                        sequence.name.resize(i);
                    }
                    return true;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fastq=>length of sequence and quality are not equal: %s") % buf);
                    return false;
                }
            }
        }
    }

    return false;
}

bool FASTAReader::read(DNASeq& sequence) {
    if (_stream) {
        std::string seq;

        std::string line;
        while (std::getline(_stream, line)) {
            boost::algorithm::trim(line);
            if (line.empty()) continue;
            if (boost::algorithm::starts_with(line, ">")) {
                if (!seq.empty() && !_name.empty()) {
                    // name
                    size_t i = _name.find_first_of(" \t");
                    if (i != std::string::npos) {
                        _name.resize(i);
                    }
                    sequence.name = _name;
                    sequence.seq = seq;
                    _name = line.substr(1);
                    return true;
                } else if (!_name.empty()) {
                    LOG4CXX_WARN(logger, boost::format("fastq=>invalid line for sequence name: %s") % line);
                    return false;
                }
                _name = line.substr(1);
            } else {
                seq += line;
            }
        }

        // the last one
        if (!seq.empty() && !_name.empty()) {
            // name
            size_t i = _name.find_first_of(" \t");
            if (i != std::string::npos) {
                _name.resize(i);
            }
            sequence.name = _name;
            sequence.seq = seq;
            return true;
        }
    }

    return false;
}

bool ReadDNASequences(std::istream& stream, DNASeqList& sequences) {
    std::shared_ptr< DNASeqReader > reader(DNASeqReaderFactory::create(stream));
    if (!reader) {
        return false;
    }
    DNASeq seq;
    while (reader->read(seq)) {
        sequences.push_back(seq);
    }
    return true;
}

bool ReadDNASequences(const std::string& file, DNASeqList& sequences) {
    std::ifstream stream(file.c_str());
    return ReadDNASequences(stream, sequences);
}

bool ReadDNASequences(const std::vector< std::string >& filelist, DNASeqList& sequences) {
    for (const auto& file : filelist) {
        if (!ReadDNASequences(file, sequences)) {
            return false;
        }
    }
    return true;
}
