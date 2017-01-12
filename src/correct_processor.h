#ifndef correct_processor_h_
#define correct_processor_h_

#include "kseq.h"

#include <iostream>
#include <string>

#include <boost/property_tree/ptree.hpp>

//
// CorrectProcessor
//
//
class CorrectProcessor {
public:
    typedef boost::property_tree::ptree Options;

    CorrectProcessor(const Options& options) : _options(options) {
    }

    bool process(DNASeqReader& reader, std::ostream& output, size_t threads=1, size_t* processed=NULL) const;
    bool process(const std::string& input, const std::string& output, size_t threads=1, size_t* processed=NULL) const;

private:
    const Options& _options;
};

#endif // correct_processor_h_
