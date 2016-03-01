#include "config.h"
#include "constant.h"
#include "kseq.h"
#include "runner.h"

#include <fstream>
#include <memory>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Preprocess"));

class Preprocess : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        // parameters
        LOG4CXX_INFO(logger, "Parameters:");
        LOG4CXX_INFO(logger, boost::format("Min length: %d") % options.get< size_t >("min-length", 40));

        // input
        std::vector< std::string > filelist;
        std::copy(arguments.begin(), arguments.end(), std::back_inserter(filelist));
        LOG4CXX_DEBUG(logger, boost::format("input: %s") % boost::algorithm::join(filelist, ":"));

        // output
        std::ostream* out = &std::cout;
        if (options.find("out") != options.not_found()) {
            std::string file = options.get< std::string >("out");
            out = new std::ofstream(file.c_str());
        }

        // process
        if (*out) {
            BOOST_FOREACH(const std::string& file, filelist) {
                LOG4CXX_INFO(logger, boost::format("Processing %s") % file);

                std::ifstream stream(file.c_str());
                std::shared_ptr< DNASeqReader > reader(DNASeqReaderFactory::create(stream));
                if (reader) {
                    DNASeq seq;
                    while (reader->read(seq)) {
                        if (processRead(options, seq)) {
                            *out << seq;
                        }
                    }
                } else {
                    LOG4CXX_ERROR(logger, boost::format("Failed to open %s") % file);
                }
            }
        } else {
            LOG4CXX_ERROR(logger, "Failed to open output stream");
        }

        if (out != &std::cout) {
            delete out;
        }

        return r;
    }

private:
    bool processRead(const Properties& options, DNASeq& record) const {
        // Ensure sequence is entirely ACGT
        if (record.seq.find_first_not_of("ACGT") != std::string::npos) {
            return false;
        }
        if (record.seq.length() < options.get< size_t >("min-length", 40)) {
            return false;
        }
        return true;
    }

    Preprocess() : Runner("c:s:o:p:q:m:h", boost::assign::map_list_of('o', "out")('p', "pe-mode")('m', "min-length")) {
        RUNNER_INSTALL("preprocess", this, "filter and quality-trim reads");
    }

    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("h") != options.not_found()) {
            return printHelps();
        }
        return 0;
    }
    int printHelps() const {
        std::cout << boost::format(
                "%s preprocess [OPTION] READS1 READS2 ...\n"
                "Prepare READS1, READS2, ... data files for assembly\n"
                "If pe-mode is turned on (pe-mode=1) then if a read is discarded its pair will be discarded as well.\n"
                "\n"
                "      -h, --help                       display this help and exit\n"
                "          --seed                       set random seed\n"
                "\n"
                "Input/Output options:\n"
                "      -o, --out=FILE                   write the reads to FILE (default: stdout)\n"
                "      -p, --pe-mode=INT                0 - do not treat reads as paired (default)\n"
                "                                       1 - reads are paired with the first read in READS1 and the second\n"
                "                                       read in READS2. The paired reads will be interleaved in the output file\n"
                "                                       2 - reads are paired and the records are interleaved within a single file.\n"
                "          --pe-orphans=FILE            if one half of a read pair fails filtering, write the passed half to FILE\n"
                "\n"
                "Conversions/Filtering:\n"
                "          --phred64                    convert quality values from phred-64 to phred-33.\n"
                "          --discard-quality            do not output quality scores\n"
                "      -q, --quality-trim=INT           perform Heng Li's BWA quality trim algorithm. \n"
                "                                       Reads are trimmed according to the formula:\n"
                "                                       argmax_x{\\sum_{i=x+1}^l(INT-q_i)} if q_l<INT\n"
                "                                       where l is the original read length.\n"
                "      -f, --quality-filter=INT         discard the read if it contains more than INT low-quality bases.\n"
                "                                       Bases with phred score <= 3 are considered low quality. Default: no filtering.\n"
                "                                       The filtering is applied after trimming so bases removed are not counted.\n"
                "                                       Do not use this option if you are planning to use the BCR algorithm for indexing.\n"
                "      -m, --min-length=INT             discard sequences that are shorter than INT\n"
                "                                       this is most useful when used in conjunction with --quality-trim. Default: 40\n"
                "      -h, --hard-clip=INT              clip all reads to be length INT. In most cases it is better to use\n"
                "                                       the soft clip (quality-trim) option.\n"
                "      --permute-ambiguous              Randomly change ambiguous base calls to one of possible bases.\n"
                "                                       If this option is not specified, the entire read will be discarded.\n"
                "      -s, --sample=FLOAT               Randomly sample reads or pairs with acceptance probability FLOAT.\n"
                "      --dust                           Perform dust-style filtering of low complexity reads.\n"
                "      --dust-threshold=FLOAT           filter out reads that have a dust score higher than FLOAT (default: 4.0).\n"
                "      --suffix=SUFFIX                  append SUFFIX to each read ID\n"
                "\n"
                "Adapter/Primer checks:\n"
                "          --no-primer-check            disable the default check for primer sequences\n"
                "      -r, --remove-adapter-fwd=STRING\n"
                "      -c, --remove-adapter-rev=STRING  Remove the adapter STRING from input reads.\n"
                "\n"
                ) % PACKAGE_NAME << std::endl;

        return 256;
    }

    static Preprocess _runner;
};

Preprocess Preprocess::_runner;
