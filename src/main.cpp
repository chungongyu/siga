#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

#include "constant.h"
#include "runner.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.main"));

int main(int argc, char* argv[]) {
    // Synchronizing iostreams with printf-style I/O can be costly. 
    // std::cin and std::cout are by default synchronized with printf.
    std::ios::sync_with_stdio(false);

    if (argc < 2) {
        return RunnerManager::get()->help(argc, argv);
    }

    RunnerPtr runner = RunnerManager::get()->create(argv[1]);
    if (!runner) {
        return RunnerManager::get()->help(argc, argv);
    }

    Properties options;
    {
        // command line options
        Properties cmd;
        const option* longopts = NULL;
        const std::string& shortopts = runner->options(&longopts);
        int opt = -1;
        while ((opt = getopt_long(argc - 1, argv + 1, shortopts.c_str(), longopts, NULL)) != -1) {
            int flag = 0;
            std::string key = runner->transform((char)opt, &flag);
            if (flag != required_argument) {
                cmd.put(key, NULL);
            } else {
                std::string val = optarg;
                // multiple key=value
                if (cmd.find(key) != cmd.not_found()) {
                    val = boost::str(boost::format("%s:%s") % cmd.get<std::string>(key) % val);
                }
                cmd.put(key, val);
            }
        }

        // config log4cxx.
        const std::string log_config = cmd.get<std::string>("log4cxx", "log4cxx.properties");
        if (boost::filesystem::exists(log_config)) {
            log4cxx::PropertyConfigurator::configure(log_config);
        } else {
            log4cxx::BasicConfigurator::configure();
        }
        
        // load ini options
        if (cmd.find("ini") != cmd.not_found()) {
            const std::string file_config = cmd.get<std::string>("ini");
            try {
                boost::property_tree::read_ini(file_config, options);
            } catch (const boost::property_tree::ini_parser_error& e) {
                LOG4CXX_ERROR(logger, boost::format("load %s failed(%s).") % file_config % e.what());
                return 1;
            }
        }

        // merge options
        for (auto it = cmd.begin(); it != cmd.end(); it++){
            options.put(it->first,it->second.data());
        }
    }

    Arguments arguments;
    std::copy(argv + optind + 1, argv + argc, std::back_inserter(arguments));

    return runner->run(options, arguments);
}
