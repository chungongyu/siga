#ifndef asqg_h_
#define asqg_h_

#include <sstream>
#include <string>
#include <vector>

#include "coord.h"

//
// The ASQG format describes an assembly graph
//
namespace ASQG {
    const char FIELD_SEP = '\t';
    const char TAG_SEP = ':';

    enum RecordType {
        RT_NONE = -1, 
        RT_HEADER = 0,
        RT_VERTEX,
        RT_EDGE
    };

    void tokenize(std::vector<std::string>& container, const std::string& text, char sep);

    template <class T>
    class TagValue {
    public:
         TagValue() : _initialized(false) {
         }
         TagValue(const T& value) : _value(value), _initialized(true) {
         }
         std::string tostring(const std::string& key) const {
             std::stringstream ss;
             ss << key << TAG_SEP << typecode(_value) << TAG_SEP << _value;
             return ss.str();
         }
         bool fromstring(const std::string& text) {
             std::vector<std::string> tokens;
             tokenize(tokens, text, TAG_SEP);
             if (tokens.size() != 3) {
                 return false;
             }
             if (tokens[1].length() != 1 || tokens[1][0] != typecode(_value)) {
                 return false;
             }
             std::stringstream ss(tokens[2]);
             ss >> _value;
             _initialized = true;
             return true;
         }
         operator bool() const {
             return _initialized;
         }
         operator T() const {
             assert(_initialized);
             return _value;
         }
    private:
         char typecode(int) const {
             return 'i';
         }
         char typecode(char) const {
             return 'A';
         }
         char typecode(const std::string&) const {
             return 'Z';
         }
         char typecode(float) const {
             return 'f';
         }

         T _value;
         bool _initialized;
    };

    typedef TagValue<int>         IntTagValue;
    typedef TagValue<float>       FloatTagValue;
    typedef TagValue<std::string> StringTagValue;

    // A header record is just a tag:value pairs
    class HeaderRecord {
    public:
        HeaderRecord();

        static bool parse(const std::string& text, HeaderRecord& record);

        void overlap(int overlap) {
            _overlap = overlap;
        }
        int overlap() const {
            return _overlap;
        }
        void infile(const std::string& filename) {
            _infile = filename;
        }
        void errorRate(float errorRate) {
            _errorRate = errorRate;
        }
        const FloatTagValue& errorRate() const {
            return _errorRate;
        }
        void containment(int v) {
            _containment = v;
        }
        const IntTagValue& containment() const {
            return _containment;
        }
        void transitive(int v) {
            _transitive = v;
        }
        const IntTagValue& transitive() const {
            return _transitive;
        }
    private:
        friend std::ostream& operator<<(std::ostream& stream, const HeaderRecord& record);
        friend std::istream& operator>>(std::istream& stream, HeaderRecord& record);

        IntTagValue _version;
        FloatTagValue _errorRate;
        StringTagValue _infile;
        IntTagValue _overlap;
        IntTagValue _containment;
        IntTagValue _transitive;
    };

    // A vertex record is an id, sequence and an array of
    // tag:value 
    class VertexRecord {
    public:
        VertexRecord() {
        }
        VertexRecord(const std::string& id, const std::string& seq) : id(id), seq(seq) {
        }

        static bool parse(const std::string& text, VertexRecord& record);

        std::string id;
        std::string seq;
        IntTagValue substring;
    private:
        friend std::ostream& operator<<(std::ostream& stream, const VertexRecord& record);
        friend std::istream& operator>>(std::istream& stream, VertexRecord& record);
    };

    // An edge record is just an overlap object and tag:values
    class EdgeRecord {
    public:
        EdgeRecord() {
        }
        EdgeRecord(const Overlap& overlap) : _overlap(overlap) {
        }

        static bool parse(const std::string& text, EdgeRecord& record);

        void cigar(const std::string& v) {
            _cigar = v;
        }
        void identity(float pi) {
            _identity = pi;
        }
        const Overlap& overlap() const {
            return _overlap;
        }
    private:
        friend std::ostream& operator<<(std::ostream& stream, const EdgeRecord& record);
        friend std::istream& operator>>(std::istream& stream, EdgeRecord& record);

        Overlap _overlap;
        StringTagValue _cigar;
        FloatTagValue _identity;
    };

    RecordType recordtype(const std::string& record);
};

#endif // asqg_h_
