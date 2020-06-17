#include "asqg.h"
#include "constant.h"
#include "utils.h"

#include <cassert>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.ASQG"));

namespace ASQG {
    static const int HEADER_VERSION = 1;
    // Record ID tags
    static const std::string HEAD_TAG("HT");
    static const std::string VERTEX_TAG("VT");
    static const std::string EDGE_TAG("ED");

    // Header tags
    static const std::string VERSION_TAG("VN");
    static const std::string OVERLAP_TAG("OL");
    static const std::string INFILE_TAG("IN");
    static const std::string ERRRATE_TAG("ER");
    static const std::string CONTAINMENT_TAG("CN");
    static const std::string TRANSITIVE_TAG("TE");

    // Vertex tags
    static const std::string SUBSTRING_TAG("SS");

    // Edge tags
    static const std::string CIGAR_TAG("CG");
    static const std::string PERCENT_IDENTITY_TAG("PI");

    void tokenize(std::vector< std::string >& container, const std::string& text, char sep) {
        boost::algorithm::split(container, text, boost::is_from_range(sep, sep));
    }

    //
    // HeaderRecord
    //
    HeaderRecord::HeaderRecord() : _version(HEADER_VERSION) {
    }

    bool HeaderRecord::parse(const std::string& text, HeaderRecord& record) {
        std::vector< std::string > fields;
        tokenize(fields, text, FIELD_SEP);

        if (fields.empty()) {
            LOG4CXX_ERROR(logger, "Error: Header record is incomplete.");
            LOG4CXX_ERROR(logger, boost::format("Record: %s") % text);
            return false;
        }
        if (!boost::algorithm::equals(fields[0], HEAD_TAG)) {
            LOG4CXX_ERROR(logger, "Error: Header does not have a header tag");
            LOG4CXX_ERROR(logger, boost::format("Record: %s") % text);
            return false;
        }
        for (size_t i = 1; i < fields.size(); ++i) {
            if (boost::algorithm::starts_with(fields[i], VERSION_TAG)) {
                if (!record._version.fromstring(fields[i])) {
                    return false;
                }
            } else if (boost::algorithm::starts_with(fields[i], OVERLAP_TAG)) {
                if (!record._overlap.fromstring(fields[i])) {
                    return false;
                }
            } else if (boost::algorithm::starts_with(fields[i], INFILE_TAG)) {
                if (!record._infile.fromstring(fields[i])) {
                    return false;
                }
            } else if (boost::algorithm::starts_with(fields[i], ERRRATE_TAG)) {
                if (!record._errorRate.fromstring(fields[i])) {
                    return false;
                }
            } else if (boost::algorithm::starts_with(fields[i], CONTAINMENT_TAG)) {
                if (!record._containment.fromstring(fields[i])) {
                    return false;
                }
            } else if (boost::algorithm::starts_with(fields[i], TRANSITIVE_TAG)) {
                if (!record._transitive.fromstring(fields[i])) {
                    return false;
                }
            }
        }
        
        return true;
    }

    std::ostream& operator<<(std::ostream& stream, const HeaderRecord& record) { // Version
        std::vector< std::string > fields;
        assert(record._version);
        fields.push_back(record._version.tostring(VERSION_TAG));

        if (record._errorRate) {
            fields.push_back(record._errorRate.tostring(ERRRATE_TAG));
        }
        if (record._overlap) {
            fields.push_back(record._overlap.tostring(OVERLAP_TAG));
        }
        if (record._infile) {
            fields.push_back(record._infile.tostring(INFILE_TAG));
        }
        if (record._containment) {
            fields.push_back(record._containment.tostring(CONTAINMENT_TAG));
        }
        if (record._transitive) {
            fields.push_back(record._transitive.tostring(TRANSITIVE_TAG));
        }

        stream << HEAD_TAG;
        for (const auto& item : fields) {
            stream << FIELD_SEP;
            stream << item;
        }
        return stream;
    }

    std::istream& operator>>(std::istream& stream, HeaderRecord& record) {
        std::string line;
        std::getline(stream, line);
        bool r = HeaderRecord::parse(line, record);
        assert(r);
        return stream;
    }

    //
    // VertexRecord
    //
    bool VertexRecord::parse(const std::string& text, VertexRecord& record) {
        std::vector< std::string > fields;
        tokenize(fields, text, FIELD_SEP);
        if (fields.size() < 3) {
            LOG4CXX_ERROR(logger, "Error: Vertex record is incomplete.");
            LOG4CXX_ERROR(logger, boost::format("Record: %s") % text);
            return false;
        }
        if (!boost::algorithm::equals(fields[0], VERTEX_TAG)) {
            LOG4CXX_ERROR(logger, "Error: Record does not have a vertex tag");
            LOG4CXX_ERROR(logger, boost::format("Record: %s") % text);
            return false;
        }
        record.id = fields[1];
        record.seq = fields[2];

        for (size_t i = 3; i < fields.size(); ++i) {
            if (boost::algorithm::starts_with(fields[i], SUBSTRING_TAG)) {
                if (!record.substring.fromstring(fields[i])) {
                    return false;
                }
            }
        }
        return true;
    }

    std::ostream& operator<<(std::ostream& stream, const VertexRecord& record) {
        stream << VERTEX_TAG << FIELD_SEP << record.id << FIELD_SEP << record.seq;
        if (record.substring) {
            stream << FIELD_SEP << record.substring.tostring(SUBSTRING_TAG);
        }
        return stream;
    }

    std::istream& operator>>(std::istream& stream, VertexRecord& record) {
        std::string line;
        std::getline(stream, line);
        bool r = VertexRecord::parse(line, record);
        assert(r);
        return stream;
    }

    //
    // EdgeRecord
    //
    bool EdgeRecord::parse(const std::string& text, EdgeRecord& record) {
        std::vector< std::string > fields;
        tokenize(fields, text, FIELD_SEP);
        if (fields.size() < 2) {
            LOG4CXX_ERROR(logger, "Error: Edge record is incomplete.");
            LOG4CXX_ERROR(logger, boost::format("Record: %s") % text);
            return false;
        }
        if (!boost::algorithm::equals(fields[0], EDGE_TAG)) {
            LOG4CXX_ERROR(logger, "Error: Edge does not have a edge tag");
            LOG4CXX_ERROR(logger, boost::format("Record: %s") % text);
            return false;
        }
        std::stringstream ss(fields[1]);
        ss >> record._overlap;
        for (size_t i = 2; i < fields.size(); ++i) {
            if (boost::algorithm::starts_with(fields[i], CIGAR_TAG)) {
                if (!record._cigar.fromstring(fields[i])) {
                    return false;
                }
            } else if (boost::algorithm::starts_with(fields[i], PERCENT_IDENTITY_TAG)) {
                if (!record._identity.fromstring(fields[i])) {
                    return false;
                }
            }
        }
        return true;
    }

    std::ostream& operator<<(std::ostream& stream, const EdgeRecord& record) {
        stream << EDGE_TAG << FIELD_SEP << record._overlap;
        if (record._cigar) {
            stream << FIELD_SEP << record._cigar.tostring(CIGAR_TAG);
        }
        if (record._identity) {
            stream << FIELD_SEP << record._identity.tostring(PERCENT_IDENTITY_TAG);
        }
        return stream;
    }

    std::istream& operator>>(std::istream& stream, EdgeRecord& record) {
        std::string line;
        std::getline(stream, line);
        bool r = EdgeRecord::parse(line, record);
        assert(r);
        return stream;
    }

    RecordType recordtype(const std::string& record) {
        if (boost::algorithm::starts_with(record, HEAD_TAG)) {
            return RT_HEADER;
        } else if (boost::algorithm::starts_with(record, VERTEX_TAG)) {
            return RT_VERTEX;
        } else if (boost::algorithm::starts_with(record, EDGE_TAG)) {
            return RT_EDGE;
        }
        return RT_NONE;
    }

    std::streambuf* ifstreambuf(const std::string& filename) {
        boost::iostreams::filtering_istreambuf* buf = new boost::iostreams::filtering_istreambuf();
        try {
            if (boost::algorithm::ends_with(filename, GZIP_EXT)) {
                buf->push(boost::iostreams::gzip_decompressor());
            }
            buf->push(boost::iostreams::file_descriptor_source(filename));
        } catch (...) {
            SAFE_DELETE(buf);
        }
        return buf;
    }

    std::streambuf* ofstreambuf(const std::string& filename) {
        boost::iostreams::filtering_ostreambuf* buf = new boost::iostreams::filtering_ostreambuf();
        try {
            if (boost::algorithm::ends_with(filename, GZIP_EXT)) {
                buf->push(boost::iostreams::gzip_compressor());
            }
            buf->push(boost::iostreams::file_descriptor_sink(filename));
        } catch (...) {
            SAFE_DELETE(buf);
        }
        return buf;
    }
};
