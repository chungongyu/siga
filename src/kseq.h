#ifndef kseq_h_
#define kseq_h_

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "quality.h"

char make_complement_dna(char c);
void make_complement_dna(std::string& dna);
std::string make_complement_dna_copy(const std::string& dna);
void make_reverse_dna(std::string& dna);
std::string make_reverse_dna_copy(const std::string& dna);
void make_reverse_complement_dna(std::string& dna);
std::string make_reverse_complement_dna_copy(const std::string& dna);

//
// DNASeq represents a DNA sequence.
//
class DNASeq {
public:
    DNASeq() {}
    DNASeq(const std::string& name, const std::string& seq) : name(name), seq(seq) {}
    DNASeq(const std::string& name, const std::string& seq, const std::string& quality) : name(name), seq(seq), quality(quality) {}
    virtual ~DNASeq() {}

    std::string name;
    std::string seq;
    std::string quality;

    int score(size_t i) const {
        if (!quality.empty()) {
            assert(i < quality.length());
            return Quality::Phred::fromchar(quality[i]);
        }
        return Quality::Phred::DEFAULT_SCORE;
    }

    void make_complement();
    void make_reverse();
    void make_reverse_complement();
private:
    friend class FASTQReader;
    friend std::ostream& operator << (std::ostream& os, const DNASeq& seq);
};

typedef std::vector<DNASeq> DNASeqList;
bool ReadDNASequences(std::istream& stream, DNASeqList& sequences);
bool ReadDNASequences(const std::string& file, DNASeqList& sequences);
bool ReadDNASequences(const std::vector<std::string>& filelist, DNASeqList& sequences);

class DNASeqReader {
public:
    DNASeqReader(std::istream& stream) : _stream(stream) {
        _pos = 0;//_stream.tellg();
    }
    virtual ~DNASeqReader() {
    }

    virtual void reset() {
        _stream.clear();
        _stream.seekg(_pos);
    }
    virtual bool read(DNASeq& sequence) = 0;

    bool hasAttr(const std::string& key) const {
        return _attrs.find(key) != _attrs.end();
    }
    std::string getAttr(const std::string& key, const std::string& defval="") const {
        auto it = _attrs.find(key);
        if (it != _attrs.end()) {
            return it->second;
        }
        return defval;
    }
    void setAttr(const std::string& key, const std::string& val) {
        _attrs[key] = val;
    }
protected:
    std::istream& _stream;
    size_t _pos;
    std::map<std::string, std::string> _attrs;
};

class DNASeqReaderFactory {
public:
    static DNASeqReader* create(std::istream& stream);
};

//
// FASTQ Format Specification
// Syntax
//    <fastq>	:=	<block>+
//    <block>	:=	@<seqname>\n<seq>\n+[<seqname>]\n<qual>\n
//    <seqname>	:=	[A-Za-z0-9_.:-]+
//    <seq>	:=	[A-Za-z\n\.~]+
//    <qual>	:=	[!-~\n]+
// Requirements
//    1. The <seqname> following '+' is optional, but if it appears right after '+', it should be 
//    identical to the <seqname> following '@'.
//    2. The length of <seq> is identical the length of <qual>. Each character in <qual> represents 
//    the phred quality of the corresponding nucleotide in <seq>.
//    3. If the Phred quality is $Q, which is a non-negative integer, the corresponding quality 
//    character can be calculated with the following Perl code:
//        $q = chr(($Q<=93? $Q : 93) + 33);
//    where chr() is the Perl function to convert an integer to a character based on the ASCII  
//    table.
//    4. Conversely, given a character $q, the corresponding Phred quality can be calculated with:
//        $Q = ord($q) - 33;
//    where ord() gives the ASCII code of a character.
// 
class FASTQReader : public DNASeqReader {
public:
    FASTQReader(std::istream& stream) : DNASeqReader(stream) {
    }
    
    bool read(DNASeq& sequence);
};

//
// FASTA Format Specification
// Syntax
//    <fasta>	:=	<block>+
//    <block>	:=	><seqname>\n<seq>\n
//    <seqname>	:=	[A-Za-z0-9_.:-]+
//    <seq>	:=	[A-Za-z\n\.~]+
// 
class FASTAReader : public DNASeqReader {
public:
    FASTAReader(std::istream& stream) : DNASeqReader(stream) {
    }

    void reset() {
        _name.clear();
        DNASeqReader::reset();
    }
    bool read(DNASeq& sequence);
private:
    std::string _name;
};

#endif // kseq_h_
