#include "utils.h"

#include <cstdlib>
#include <ctime>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <log4cxx/logger.h>

#include "constant.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Utils"));

namespace Utils {
  void srand() {
    ::srand(time(NULL));
  }
  int rand() {
    return ::rand();
  }

  std::istream* ifstream(const std::string& filename) {
    boost::iostreams::filtering_istream* stream = new boost::iostreams::filtering_istream();
    try {
      if (boost::algorithm::ends_with(filename, GZIP_EXT)) {
        stream->push(boost::iostreams::gzip_decompressor());
      } else if (boost::algorithm::ends_with(filename, BZIP_EXT)) {
        stream->push(boost::iostreams::bzip2_decompressor());
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
      } else if (boost::algorithm::ends_with(filename, BZIP_EXT)) {
        stream->push(boost::iostreams::bzip2_compressor());
      }
      stream->push(boost::iostreams::file_descriptor_sink(filename));
    } catch (...) {
      SAFE_DELETE(stream);
    }
    return stream;
  }

  std::string stem(const std::string& filename) {
    if (boost::algorithm::ends_with(filename, GZIP_EXT)) {
      return stem(filename.substr(0, filename.length() - strlen(GZIP_EXT)));
    } else if (boost::algorithm::ends_with(filename, BZIP_EXT)) {
      return stem(filename.substr(0, filename.length() - strlen(BZIP_EXT)));
    }
    return boost::filesystem::path(filename).stem().string();
  }
};  // namespace Utils
