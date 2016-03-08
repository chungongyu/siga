#ifndef asqg_h_
#define asqg_h_

#include <map>
#include <string>

#include <boost/any.hpp>

#include "coord.h"

//
// The ASQG format describes an assembly graph
//
namespace ASQG {
    const char FIELD_SEP = '\t';
    const char TAG_SEP = ':';

    enum RecordType {
        RT_HEADER = 0,
        RT_VERTEX,
        RT_EDGE
    };

    // A header record is just a tag:value pairs
    class HeaderRecord {
    public:
        HeaderRecord();
        HeaderRecord(const std::string& record);

    private:
        std::map< std::string, boost::any > _values;
    };

    // A vertex record is an id, sequence and an array of
    // tag:value 
    class VertexRecord {
    public:
    private:
        std::string _id;
        std::string _seq;
    };

    // An edge record is just an overlap object and tag:values
    class EdgeRecord {
    public:
    private:
        Overlap _overlap;
    };
};

#endif // asqg_h_
