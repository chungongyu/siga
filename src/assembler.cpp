#include "assembler.h"
#include "constant.h"

#include <iostream>

#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Assembler"));

int Assembler::run(const Properties& options, const Arguments& arguments) {
    int r = 0;

    if ((r = checkOptions(options)) != 0) {
        return r;
    }

    LOG4CXX_DEBUG(logger, "assemble begin");

    LOG4CXX_DEBUG(logger, "assemble end");
    return r;
}

Assembler Assembler::_runner;

Assembler::Assembler() : Runner("c:s:o:h") {
    RUNNER_INSTALL("assemble", this, "generate contigs from an assembly graph");
}

int Assembler::printHelps() const {
    std::cout << "arcs assemble -o [input] -d [workdir]" << std::endl;
    return 256;
}

int Assembler::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
