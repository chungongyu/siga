#include "constant.h"
#include "runner.h"

#include <iostream>

#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Assembler"));

class Assembler : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        LOG4CXX_DEBUG(logger, "assemble begin");

        LOG4CXX_DEBUG(logger, "assemble end");
        return r;
    }

private:
    Assembler() : Runner("c:s:o:h") {
        RUNNER_INSTALL("assemble", this, "generate contigs from an assembly graph");
    }

    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("h") != options.not_found()) {
            return printHelps();
        }
        return 0;
    }
    int printHelps() const {
        std::cout << "arcs assemble -o [input] -d [workdir]" << std::endl;
        return 256;
    }

    static Assembler _runner;
};

Assembler Assembler::_runner;
