#include "utils.h"

#include <cassert>
#include <cstdlib>
#include <ctime>
#include <memory>

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

template <typename Alloc = std::allocator<char> >
struct basic_gzip_decompressor_delegator : boost::iostreams::basic_gzip_decompressor<Alloc> {
 public:
  typedef char char_type;
  struct category : boost::iostreams::basic_gzip_decompressor<Alloc>::category, boost::iostreams::input_seekable {};
  basic_gzip_decompressor_delegator() : impl_(new boost::iostreams::basic_gzip_decompressor<Alloc>()) {
  }
  template <typename Sink>
  std::streamsize write(Sink& snk, const char_type* s, std::streamsize n) {
    return impl_->write(snk, s, n);
  }
  template <typename Source>
  std::streamsize read(Source& src, char_type* s, std::streamsize n) {
    return impl_->read(src, s, n);
  }
  template <typename Source>
  void close(Source& src, BOOST_IOS::openmode m) {
    impl_->close(src, m);
  }
  template <typename Source>
  void close(Source& src) {
    this->close(src, BOOST_IOS::in);
  }
  template <typename Source>
  std::streampos seek(Source& src, boost::iostreams::stream_offset off, BOOST_IOS::seekdir way,
      BOOST_IOS::openmode which = BOOST_IOS::in | BOOST_IOS::out) {
    assert(off == 0);
    src.pubseekoff(off, way);
    impl_.reset(new boost::iostreams::basic_gzip_decompressor<Alloc>());
    return 0;
  }

  std::string file_name() const { return impl_->file_name(); }
  std::string comment() const { return impl_->comment(); }
  bool text() const { return impl_->text(); }
  int os() const { return impl_->os(); }
  std::time_t mtime() const { return impl_->mtime(); }

 private:
  std::shared_ptr<boost::iostreams::basic_gzip_decompressor<Alloc> > impl_;
};

BOOST_IOSTREAMS_PIPABLE(basic_gzip_decompressor_delegator, 1)
typedef basic_gzip_decompressor_delegator<> gzip_decompressor_delegator;

template <typename Alloc = std::allocator<char> >
struct basic_bzip2_decompressor_delegator : boost::iostreams::basic_bzip2_decompressor<Alloc> {
 public:
  typedef char char_type;
  struct category : boost::iostreams::basic_bzip2_decompressor<Alloc>::category, boost::iostreams::input_seekable {};
  basic_bzip2_decompressor_delegator() : impl_(new boost::iostreams::basic_bzip2_decompressor<Alloc>()) {
  }
  template <typename Sink>
  std::streamsize write(Sink& snk, const char_type* s, std::streamsize n) {
    return impl_->write(snk, s, n);
  }
  template <typename Source>
  std::streamsize read(Source& src, char_type* s, std::streamsize n) {
    return impl_->read(src, s, n);
  }
  template <typename Source>
  void close(Source& src, BOOST_IOS::openmode m) {
    impl_->close(src, m);
  }
  template <typename Source>
  void close(Source& src) {
    this->close(src, BOOST_IOS::in);
  }
  template <typename Source>
  std::streampos seek(Source& src, boost::iostreams::stream_offset off, BOOST_IOS::seekdir way,
      BOOST_IOS::openmode which = BOOST_IOS::in | BOOST_IOS::out) {
    assert(off == 0);
    src.pubseekoff(off, way);
    impl_.reset(new boost::iostreams::basic_bzip2_decompressor<Alloc>());
    return 0;
  }

 private:
  std::shared_ptr<boost::iostreams::basic_bzip2_decompressor<Alloc> > impl_;
};

BOOST_IOSTREAMS_PIPABLE(basic_bzip2_decompressor_delegator, 1)
typedef basic_bzip2_decompressor_delegator<> bzip2_decompressor_delegator;

std::istream* ifstream(const std::string& filename) {
  boost::iostreams::filtering_istream* stream = new boost::iostreams::filtering_istream();
  try {
    if (boost::algorithm::ends_with(filename, GZIP_EXT)) {
      stream->push(gzip_decompressor_delegator());
      // stream->push(boost::iostreams::gzip_decompressor());
    } else if (boost::algorithm::ends_with(filename, BZIP_EXT)) {
      stream->push(bzip2_decompressor_delegator());
      // stream->push(boost::iostreams::bzip2_decompressor());
    }
    stream->push(boost::iostreams::file_descriptor_source(filename, std::ios_base::in|std::ios_base::binary));
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
