#include "utils.h"
#include "constant.h"

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Utils"));

namespace Utils {
    std::istream* ifstream(const std::string& filename) {
        boost::iostreams::filtering_istream* stream = new boost::iostreams::filtering_istream();
        try {
            if (boost::algorithm::ends_with(filename, GZIP_EXT)) {
                stream->push(boost::iostreams::gzip_decompressor());
            }
            stream->push(boost::iostreams::file_descriptor_source(filename));
        } catch (...) {
            SAFE_DELETE(stream);
        }
        return stream;
    }
    std::ostream* ofstream(const std::string& filename) {
        boost::iostreams::filtering_ostream* stream = new boost::iostreams::filtering_ostream();
        try {
            if (boost::algorithm::ends_with(filename, GZIP_EXT)) {
                stream->push(boost::iostreams::gzip_compressor());
            }
            stream->push(boost::iostreams::file_descriptor_sink(filename));
        } catch (...) {
            SAFE_DELETE(stream);
        }
        return stream;
    }
};
