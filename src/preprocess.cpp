#include "config.h"
#include "constant.h"
#include "kseq.h"
#include "primer_screen.h"
#include "quality.h"
#include "reads.h"
#include "runner.h"
#include "utils.h"

#include <fstream>
#include <memory>
#include <iostream>
#include <set>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Preprocess"));

static int LOW_QUALITY_PHRED_SCORE = 3;

class Preprocess : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        // parameters
        LOG4CXX_INFO(logger, "Parameters:");
        LOG4CXX_INFO(logger, boost::format("PE Mode: %d") % options.get<int>("pe-mode", 0));
        if (options.find("pe-orientation") != options.not_found()) {
            LOG4CXX_INFO(logger, boost::format("PE Orientation: %d") % options.get<std::string>("pe-orientation", "fr"));
        }
        if (options.find("min-length") != options.not_found()) {
            LOG4CXX_INFO(logger, boost::format("Min length: %d") % options.get<int>("min-length"));
        }
        if (options.find("hard-clip") != options.not_found()) {
            LOG4CXX_INFO(logger, boost::format("hard clip: %d") % options.get<int>("hard-clip"));
        }
        if (options.find("sample-rate") != options.not_found()) {
            LOG4CXX_INFO(logger, boost::format("sample rate: %f") % options.get<float>("sample-rate"));
        }
        if (options.find("quality-trim") != options.not_found()) {
            LOG4CXX_INFO(logger, boost::format("Quality Trim: %d") % options.get<int>("quality-trim"));
        }
        if (options.find("quality-filter") != options.not_found()) {
            LOG4CXX_INFO(logger, boost::format("Quality Filter: %d") % options.get<int>("quality-filter"));
        }

        srand(time(NULL));

        // input
        std::vector<std::string> filelist;
        std::copy(arguments.begin(), arguments.end(), std::back_inserter(filelist));
        LOG4CXX_INFO(logger, boost::format("input: %s") % boost::algorithm::join(filelist, ":"));

        // output
        std::ostream* out = &std::cout;
        if (options.find("out") != options.not_found()) {
            std::string file = options.get<std::string>("out");
            out = Utils::ofstream(file);
            LOG4CXX_INFO(logger, boost::format("output: %s") % file);
        }

        // process
        if (*out) {
            Statistics stats;
            if ((r = processReads(options, filelist, *out, stats)) == 0) {
                // report statistics
                LOG4CXX_INFO(logger, "Preprocess stats:");
                LOG4CXX_INFO(logger, boost::format("Reads parsed:\t%d") % stats.numReadsRead);
                if (stats.numReadsRead > 0) {
                    LOG4CXX_INFO(logger, boost::format("Reads kept:\t%d(%f)") % stats.numReadsKept % ((double)stats.numReadsKept / stats.numReadsRead));
                    LOG4CXX_INFO(logger, boost::format("Reads failed primer screen:\t%d(%e)") % stats.numReadsPrimer % ((double)stats.numReadsPrimer / stats.numReadsRead));
                }
                LOG4CXX_INFO(logger, boost::format("Bases parsed:\t%d") % stats.numBasesRead);
                if (stats.numBasesRead) {
                    LOG4CXX_INFO(logger, boost::format("Bases kept:\t%d(%f)") % stats.numBasesKept % ((double)stats.numBasesKept / stats.numBasesRead));
                }
                LOG4CXX_INFO(logger, boost::format("Number of incorrectly paired reads that were discarded: %d") % stats.numInvalidPE);
            }
        } else {
            LOG4CXX_ERROR(logger, "Failed to open output stream");
        }

        if (out != &std::cout) {
            SAFE_DELETE(out);
        }

        return r;
    }

private:
    struct Statistics {
        Statistics() : numReadsRead(0), numReadsKept(0), numBasesRead(0), numBasesKept(0), numReadsPrimer(0), numInvalidPE(0) {
        }
        size_t numReadsRead;
        size_t numReadsKept;
        size_t numBasesRead;
        size_t numBasesKept;
        size_t numReadsPrimer;
        size_t numInvalidPE;
    };

    bool samplePass(const Properties& options) const {
        if (options.find("sample-rate") == options.not_found()) {
            return true;
        }
        float sampleRate = options.get<float>("sample-rate");
        return rand() / (RAND_MAX + 1.0f) < sampleRate;
    }

    int processReads(const Properties& options, const std::vector<std::string>& inputs, std::ostream& output, Statistics& stats) {
        int peMode = options.get<int>("pe-mode", 0);

        if (peMode == 0) {
            return processSingleEnds(options, inputs, output, stats);
        } else if (peMode == 1) {
            return processPairEnds1(options, inputs, output, stats);
        } else if (peMode == 2) {
            return processPairEnds2(options, inputs, output, stats);
        }

        LOG4CXX_ERROR(logger, boost::format("Invalid pe mode parameter: %d") % peMode);
        return -1;
    }

    int processSingleEnds(const Properties& options, const std::vector<std::string>& inputs, std::ostream& output, Statistics& stats) {
        for (const auto& file : inputs) {
            LOG4CXX_INFO(logger, boost::format("Processing %s") % file);

            int r = -1;
            if (file == "-") {
                r = processSingleEnds(options, std::cin, output, stats);
            } else {
                std::shared_ptr<std::istream> stream(Utils::ifstream(file));
                if (stream) {
                    r = processSingleEnds(options, *stream, output, stats);
                }
            }
            if (r != 0) {
                LOG4CXX_ERROR(logger, boost::format("Failed to open input stream %s") % file);
                return r;
            }
        }
        return 0;
    }

    int processSingleEnds(const Properties& options, std::istream& input, std::ostream& output, Statistics& stats) {
        std::shared_ptr<DNASeqReader> reader(DNASeqReaderFactory::create(input));
        return processSingleEnds(options, reader.get(), output, stats);
    }

    int processSingleEnds(const Properties& options, DNASeqReader* reader, std::ostream& output, Statistics& stats) {
        if (reader) {
            DNASeq read;
            while (reader->read(read)) {
                if (processRead(options, read, stats) && samplePass(options)) {
                    output << read;
                    // Statistics
                    ++stats.numReadsKept;
                    stats.numBasesKept += read.seq.length();
                }
            }
            return 0;
        }
        return -1;
    }

    int processPairEnds1(const Properties& options, const std::vector<std::string>& inputs, std::ostream& output, Statistics& stats) {
        if (inputs.size() % 2 != 0) {
            LOG4CXX_ERROR(logger, boost::format("An even number of files must be given for pe-mode 1"));
            return -1;
        }

        size_t i = 0;
        while (i < inputs.size()) {
            std::string file1 = inputs[i++];
            std::string file2 = inputs[i++];
            LOG4CXX_INFO(logger, boost::format("Processing %s,%s") % file1 % file2);

            std::shared_ptr<std::istream> stream1(Utils::ifstream(file1)), stream2(Utils::ifstream(file2));
            if (!stream1 || !stream2) {
                LOG4CXX_ERROR(logger, boost::format("Failed to read pair ends: %s,%s") % file1 % file2);
                return -1;
            }
            std::shared_ptr<DNASeqReader> reader1(DNASeqReaderFactory::create(*stream1));
            std::shared_ptr<DNASeqReader> reader2(DNASeqReaderFactory::create(*stream2));
            int r = processPairEnds(options, reader1.get(), reader2.get(), output, stats);
            if (r != 0) {
                LOG4CXX_ERROR(logger, boost::format("Failed to process pair ends: %s,%s") % file1 % file2);
                return r;
            }
        }

        return 0;
    }

    int processPairEnds2(const Properties& options, const std::vector<std::string>& inputs, std::ostream& output, Statistics& stats) {
        for (const auto& file : inputs) {
            LOG4CXX_INFO(logger, boost::format("Processing %s") % file);

            int r = -1;
            if (file == "-") {
                r = processPairEnds(options, std::cin, std::cin, output, stats);
            } else {
                std::shared_ptr<std::istream> stream(Utils::ifstream(file));
                if (stream) {
                    r = processPairEnds(options, *stream, *stream, output, stats);
                }
            }
            if (r != 0) {
                LOG4CXX_ERROR(logger, boost::format("Failed to process pair ends: %s") % file);
                return r;
            }
        }
        return 0;
    }

    int processPairEnds(const Properties& options, std::istream& input1, std::istream& input2, std::ostream& output, Statistics& stats) {
        std::shared_ptr<DNASeqReader> reader1(DNASeqReaderFactory::create(input1));
        std::shared_ptr<DNASeqReader> reader2(DNASeqReaderFactory::create(input2));
        return processPairEnds(options, reader1.get(), reader2.get(), output, stats);
    }

    int processPairEnds(const Properties& options, DNASeqReader* reader1, DNASeqReader* reader2, std::ostream& output, Statistics& stats) {
        if (reader1 != NULL && reader2 != NULL) {
            std::string orientation = options.get("pe-orientation", "fr");

            DNASeq read1, read2;
            while (reader1->read(read1) && reader2->read(read2)) {
                // If the names of the records are the same, append a /1 and /2 to them
                if (read1.name == read2.name) {
                    read1.name += "/1";
                    read2.name += "/2";
                }

                // Ensure the read names are sensible
                std::string expectedID2 = PairEnd::id(read1.name);
                std::string expectedID1 = PairEnd::id(read2.name);

                if (expectedID1 != read1.name || expectedID2 != read2.name) {
                    LOG4CXX_WARN(logger, "Pair names do not match (expected format /1,/2 or /A,/B)");
                    LOG4CXX_WARN(logger, boost::format("Read1 name: %s") % read1.name);
                    LOG4CXX_WARN(logger, boost::format("Read2 name: %s") % read2.name);
                    // Statistics
                    stats.numInvalidPE += 2;
                }

                bool passed1 = processRead(options, read1, stats);
                bool passed2 = processRead(options, read2, stats);
                if (passed1 && passed2 && samplePass(options)) {
                    if (boost::algorithm::iequals(orientation, "fr")) {
                        read2.make_reverse_complement();
                    } else if (boost::algorithm::iequals(orientation, "rf")) {
                        read1.make_reverse_complement();
                    }
                    output << read1 << read2;
                    // Statistics
                    stats.numReadsKept += 2;
                    stats.numBasesKept += (read1.seq.length() + read2.seq.length());
                }
            }
            return 0;
        }
        return -1;
    }

    bool processRead(const Properties& options, DNASeq& record, Statistics& stats) const {
        // Statistics
        ++stats.numReadsRead;
        stats.numBasesRead += record.seq.length();

        // Ensure sequence is entirely ACGT
        if (record.seq.find_first_not_of("ACGT") != std::string::npos) {
            LOG4CXX_DEBUG(logger, boost::format("%s is not entirely ACGT") % record.name);
            return false;
        }

        // Validate the quality string (if present) and
        // perform any necessary transformations
        if (!record.quality.empty()) {
            bool phred64 = options.find("phred64") != options.not_found();

            // Calculate the range of phred scores for validation
            bool allValid = true;
            for (size_t i = 0; i < record.quality.length(); ++i) {
                if (phred64) {
                    record.quality[i] = Quality::Phred::_64to33(record.quality[i]);
                }
                allValid = Quality::Phred::isValid(record.quality[i]) && allValid;
            }

            if (!allValid) {
                LOG4CXX_ERROR(logger, boost::format("Error: read %s has out of range quality values.") % record.name);
                LOG4CXX_ERROR(logger, boost::format("Expected phred%d.") % 33);
                LOG4CXX_ERROR(logger, boost::format("Quality string: %s") % record.quality);
                LOG4CXX_ERROR(logger, boost::format("Check your data and re-run preprocess with the correct quality scaling flag."));
            }
        }

        // Hard clip
        {
            int maxLength = options.get<int>("hard-clip", 0);
            if (maxLength > 0) {
                hardClip(maxLength, record);
            }
        }

        // Quality trim
        {
            int qualityTrim = options.get<int>("quality-trim", 0);
            if (qualityTrim > 0 && !record.quality.empty()) {
                softClip(qualityTrim, record);
            }
        }

        // Quality filter
        {
            int qualityFilter = options.get<int>("quality-filter", -1);
            if (qualityFilter >= 0 && !record.quality.empty()) {
                size_t numLowQuality = countLowQuality(record);
                if (numLowQuality >= qualityFilter) {
                    return false;
                }
            }
        }

        // Primer screen
        if (options.find("no-primer-check") == options.not_found()) {
            if (PrimerScreen::containsPrimer(record.seq)) {
                // Statistics
                ++stats.numReadsPrimer;
                return false;
            }
        }

        // Min length
        if (record.seq.length() < options.get<size_t>("min-length", 40)) {
            return false;
        }
        return true;
    }

    // Count the number of low quality bases in the read
    size_t countLowQuality(const DNASeq& record) const {
        assert(record.seq.length() == record.quality.length());

        size_t n = 0;
        for (auto q : record.quality) {
            int ps = Quality::Phred::fromchar(q);
            if (ps <= LOW_QUALITY_PHRED_SCORE) {
                ++n;
            }
        }
        return n;
    }

    // Perform a soft-clipping of the sequence by removing low quality bases from the
    // 3' end using Heng Li's algorithm from bwa
    void softClip(int qualityTrim, DNASeq& record) const {
        assert(record.seq.length() == record.quality.length());

        size_t i = record.seq.length();
        int terminalScore = Quality::Phred::fromchar(record.quality[i - 1]);
        // Only perform soft-clipping if the last base has qual less than qualTrim
        if (terminalScore < qualityTrim) {
            size_t endpoint = 0; // not inclusive
            int max = 0;

            int subSum = 0;
            while (i > 0) {
                int ps = Quality::Phred::fromchar(record.quality[i - 1]);
                int score = qualityTrim - ps;
                subSum += score;
                if (subSum > max) {
                    max = subSum;
                    endpoint = i;
                }
                --i;
            }

            // Clip the read
            hardClip(endpoint, record);
        }
    }

    // Perform a hard-clipping
    void hardClip(size_t endpoint, DNASeq& record) const {
        if (record.seq.length() > endpoint) {
            record.seq.resize(endpoint);
        }
        if (record.quality.length() > endpoint) {
            record.quality.resize(endpoint);
        }
    }

    Preprocess(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts) : Runner(shortopts, longopts) {
        RUNNER_INSTALL(name, this, description, kPreprocess);
    }

    int checkOptions(const Properties& options, const Arguments& arguments) const {
        std::set<std::string> orientations = {"fr", "ff", "rf"};
        if (options.find("help") != options.not_found() 
                || orientations.find(options.get<std::string>("pe-orientation", "fr")) == orientations.end() 
                || arguments.empty()) {
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
                "\n"
                "Input/Output options:\n"
                "      -o, --out=FILE                   write the reads to FILE (default: stdout)\n"
                "          --pe-mode=INT                0 - do not treat reads as paired (default)\n"
                "                                       1 - reads are paired with the first read in READS1 and the second\n"
                "                                       read in READS2. The paired reads will be interleaved in the output file\n"
                "                                       2 - reads are paired and the records are interleaved within a single file.\n"
                "          --pe-orientation=STR         orientation of reads for paired-end reads. orientation = fr(default), rf, ff.\n"
                "\n"
                "Conversions/Filtering:\n"
                "          --phred64                    convert quality values from phred-64 to phred-33.\n"
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
                "          --hard-clip=INT              clip all reads to be length INT. In most cases it is better to use\n"
                "                                       the soft clip (quality-trim) option.\n"
                "          --sample-rate=FLOAT          randomly sample reads or pairs with acceptance probability FLOAT.\n"
                "\n"
                "Adapter/Primer checks:\n"
                "          --no-primer-check            disable the default check for primer sequences\n"
                "\n"
                ) % PACKAGE_NAME << std::endl;

        return 256;
    }

    static Preprocess _runner;
};

static const std::string shortopts = "c:s:o:p:q:f:m:h";
enum { OPT_HELP = 1, OPT_PE_MODE, OPT_PE_ORIENTATION, OPT_PHRED64, OPT_HARD_CLIP, OPT_SAMPLE_RATE, OPT_NO_PRIMER_CHECK };
static const option longopts[] = {
    {"log4cxx",             required_argument,  NULL, 'c'}, 
    {"ini",                 required_argument,  NULL, 's'}, 
    {"out",                 required_argument,  NULL, 'o'}, 
    {"pe-mode",             required_argument,  NULL, OPT_PE_MODE}, 
    {"pe-orientation",      required_argument,  NULL, OPT_PE_ORIENTATION}, 
    {"phred64",             no_argument,        NULL, OPT_PHRED64}, 
    {"quality-trim",        required_argument,  NULL, 'q'}, 
    {"quality-filter",      required_argument,  NULL, 'f'}, 
    {"min-length",          required_argument,  NULL, 'm'}, 
    {"hard-clip",           required_argument,  NULL, OPT_HARD_CLIP}, 
    {"sample-rate",         required_argument,  NULL, OPT_SAMPLE_RATE}, 
    {"no-primer-check",     no_argument,        NULL, OPT_NO_PRIMER_CHECK}, 
    {"help",                no_argument,        NULL, 'h'}, 
    {NULL, 0, NULL, 0}, 
};
Preprocess Preprocess::_runner(
        "preprocess", 
        "filter and quality-trim reads", 
        shortopts,  
        longopts 
        );
