#include "kseq.h"

#include <fstream>
#include <map>
#include <memory>
#include <numeric>
#include <unordered_map>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

#include "utils.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.DNASeq"));

char make_dna_complement(char c) {
  static std::map<char, char> mapping = {
      {'A', 'T'},
      {'C', 'G'},
      {'G', 'C'},
      {'T', 'A'},
      {'N', 'N'},
    };
  return mapping[c];
}

void make_dna_complement(std::string& sequence) {
  static std::map<char, char> mapping = {
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

std::string make_dna_complement_copy(const std::string& sequence) {
  std::string complement = sequence;
  make_dna_complement(complement);
  return complement;
}

void make_dna_reverse(std::string& sequence) {
  std::reverse(sequence.begin(), sequence.end());
}

std::string make_dna_reverse_copy(const std::string& sequence) {
  std::string reverse = sequence;
  make_dna_reverse(reverse);
  return reverse;
}

void make_dna_reverse_complement(std::string& sequence) {
  make_dna_complement(sequence);
  make_dna_reverse(sequence);
}

std::string make_dna_reverse_complement_copy(const std::string& sequence) {
  std::string complement = make_dna_complement_copy(sequence);
  make_dna_reverse(complement);
  return complement;
}

void make_seq_name(std::string& name, std::string& comment) {
  size_t i = name.find_first_of(" \t");
  if (i != std::string::npos) {
    comment = name.substr(i + 1);
    name.resize(i);
  } else {
    comment.clear();
  }
}

DNASeq::DNASeq(const std::string& name, const std::string& seq) : name(name), seq(seq) {
  make_seq_name(this->name, this->comment);
}

DNASeq::DNASeq(const std::string& name, const std::string& seq, const std::string& quality)
    : name(name), seq(seq), quality(quality) {
  make_seq_name(this->name, this->comment);
}

void DNASeq::make_complement() {
  make_dna_complement(seq);
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
    os << '>' << seq.name;
    if (!seq.comment.empty()) {
      os << ' ' << seq.comment;
    }
    os << '\n';
    os << seq.seq << '\n';
  } else {
    os << '@' << seq.name;
    if (!seq.comment.empty()) {
      os << ' ' << seq.comment;
    }
    os << '\n';
    os << seq.seq << '\n';
    os << '+' << '\n';
    os << seq.quality << '\n';
  }
  return os;
}

DNASeqReader* DNASeqReaderFactory::create(std::istream& stream) {
  if (stream) {
    int c = stream.peek();
    if (c == '@') {
      return new FASTQReader(stream);
    } else if (c == '>') {
      return new FASTAReader(stream);
    }
    LOG4CXX_ERROR(logger, boost::format("failed to create dna reader: %c") % c);
  }
  return nullptr;
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
        if (boost::algorithm::starts_with(buf, "+") && (buf.length() == 1 ||
            boost::algorithm::ends_with(buf, sequence.name))) {
          state = kQuality;
        } else {
          LOG4CXX_WARN(logger, boost::format("fastq=>names aren't equal: %s") % buf);
          return false;
        }
      } else if (state == kQuality) {
        if (buf.length() == sequence.seq.length()) {
          sequence.quality = buf;
          // name
          make_seq_name(sequence.name, sequence.comment);
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
          sequence.name = _name;
          make_seq_name(sequence.name, sequence.comment);
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
      sequence.name = _name;
      make_seq_name(sequence.name, sequence.comment);
      sequence.seq = seq;
      return true;
    }
  }

  return false;
}

bool ReadDNASequences(std::istream& stream, DNASeqList& sequences, uint32_t flags) {
  std::shared_ptr<DNASeqReader> reader(DNASeqReaderFactory::create(stream));
  if (!reader) {
    return false;
  }
  DNASeq seq;
  while (reader->read(seq)) {
    if (!(flags & kSeqWithQuality)) {
      seq.quality.clear();
    }
    if (!(flags & kSeqWithComment)) {
      seq.comment.clear();
    }
    sequences.push_back(seq);
  }
  if (flags & kSeqShrink) {
    sequences.shrink_to_fit();
  }
  return true;
}

bool ReadDNASequences(const std::string& file, DNASeqList& sequences, uint32_t flags) {
  std::shared_ptr<std::istream> stream(Utils::ifstream(file));
  if (stream) {
    return ReadDNASequences(*stream, sequences, flags);
  }
  return false;
}

bool ReadDNASequences(const std::vector<std::string>& filelist, DNASeqList& sequences, uint32_t flags) {
  for (const auto& file : filelist) {
    if (!ReadDNASequences(file, sequences, flags)) {
      return false;
    }
  }
  return true;
}
