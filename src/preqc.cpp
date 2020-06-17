#include "config.h"
#include "constant.h"
#include "fmindex.h"
#include "kmerdistr.h"
#include "kseq.h"
#include "quality.h"
#include "reads.h"
#include "runner.h"
#include "utils.h"

#include <fstream>
#include <memory>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.PreQC"));

#define kPreQCThreads   1

class PreQC : public Runner {
public:
    typedef rapidjson::StringBuffer JSONBuffer;
    typedef rapidjson::PrettyWriter< JSONBuffer > JSONWriter;

    int run(const Properties& options, const Arguments& arguments) {
        int r = 0;

        if ((r = checkOptions(options, arguments)) != 0) {
            return r;
        }

        // parameters
        LOG4CXX_INFO(logger, "Parameters:");

        // input
        std::string input = arguments[0];
        LOG4CXX_INFO(logger, boost::format("input: %s") % input);

        JSONBuffer qc;
        {
            JSONWriter writer(qc);

            // Top-level document
            writer.StartObject();

            // In simple mode we only compute metrics that do not need the FM-index
            QualityStatistics statistics(options.get< double >("sample-rate", 0.05));
            r = statistics.stats(input, &writer);

            if (r == 0 && options.find("simple") == options.not_found()) {
                std::string prefix = boost::filesystem::path(input).stem().string();
                // Load the FM-index and compute the rest of the metrics if requested
                LOG4CXX_INFO(logger, boost::format("Loading FM-index of %s") % prefix);

                GenomeEstimator estimator(options.get< size_t >("kmer", 31), options.find("diploid-reference-mode") != options.not_found());
                r = estimator.estimate(prefix, &writer);
            }

            // End document
            writer.EndObject();
        }

        std::cout << qc.GetString() << std::endl;

        return r;
    }

private:
    class QualityStatistics {
    public:
        struct Statistics {
            Statistics() : qualitysum(0), q30num(0), count(0) {
            }
            size_t qualitysum;
            size_t q30num;
            size_t count;
        };

        QualityStatistics(double sampleRate) : _sampleRate(sampleRate) {
        }

        int stats(const std::string& file, JSONWriter* writer) const {
            int r = 0;

            LOG4CXX_INFO(logger, boost::format("Processing %s") % file);
            if (file == "-") {
                r = stats(std::cin, writer);
            } else {
                std::ifstream stream(file.c_str());
                r = stats(stream, writer);
            }
            if (r != 0) {
                LOG4CXX_ERROR(logger, boost::format("Failed to open input stream %s") % file);
            }
            return r;
        }

        int stats(std::istream& stream, JSONWriter* writer) const {
            std::shared_ptr< DNASeqReader > reader(DNASeqReaderFactory::create(stream));
            return stats(reader.get(), writer);
        }

        int stats(DNASeqReader* reader, JSONWriter* writer) const {
            if (reader) {
                std::vector< Statistics > bases;

                DNASeq read;
                while (reader->read(read)) {
                    if ((double)rand() / RAND_MAX < _sampleRate && read.seq.length() == read.quality.length()) {
                        size_t l = read.seq.length();
                        if (l > bases.size()) {
                            bases.resize(l);
                        }
                        for (size_t i = 0; i < l; ++i) {
                            Statistics& base = bases[i];
                            int q = Quality::Phred::fromchar(read.quality[i]);
                            ++base.count;
                            base.qualitysum += q;
                            base.q30num += (q >= 30) ? 1 : 0;
                        }
                    }
                }

                writer->String("QualityScores");
                writer->StartObject();
                {
                    writer->String("mean_quality");
                    writer->StartArray();
                    for (const auto& base : bases) {
                        writer->Double((double)base.qualitysum / base.count);
                    }
                    writer->EndArray();

                    writer->String("fraction_q30");
                    writer->StartArray();
                    for (const auto& base : bases) {
                        writer->Double((double)base.q30num / base.count);
                    }
                    writer->EndArray();
                }
                writer->EndObject();

                return 0;
            }
            return -1;
        }

    private:
        double _sampleRate;
    };

    class GenomeEstimator {
    public:
        GenomeEstimator(size_t kmer = 31, bool diploid = false) : _kmer(kmer), _diploid(diploid), _samples(50000) {
        }

        int estimate(const std::string& prefix, JSONWriter* writer) const {
            int r = 0;
            if ((r = estimateSize(NULL, writer)) != 0) {
                return r;
            }
            return r;
        }
    private:
        int estimateSize(const FMIndex* index, JSONWriter* writer) const {
            KmerDistribution kmerdistr;
            size_t nl = KmerDistribution::sample(index, _kmer, _samples, &kmerdistr);

            // calculate the k-mer count model parameters from the distribution
            // this gives us the estimated proportion of kmers that contain errors

            writer->String("GenomeSize");
            writer->StartObject();
            {
                writer->String("k");
                writer->Int64(_kmer);
                writer->String("size");
                writer->Int64(0);
            }
            writer->EndObject();
            return 0;
        }
        size_t _samples;
        size_t _kmer;
        bool _diploid;
    };

    PreQC(const std::string& name, const std::string& description, const std::string& shortopts, const option* longopts) : Runner(shortopts, longopts) {
        RUNNER_INSTALL(name, this, description, kPreQC);
    }
    int checkOptions(const Properties& options, const Arguments& arguments) const {
        if (options.find("help") != options.not_found() || arguments.empty()) {
            return printHelps();
        }
        return 0;
    }
    int printHelps() const {
        std::cout << boost::format(
                "%s preqc [OPTION] READSFILE\n"
                "Preform pre-assembly quality checks\n"
                "\n"
                "      -h, --help                       display this help and exit\n"
                "\n"
                "      -t, --threads=NUM                use NUM threads (default: %d)\n"
                "          --simple                     only compute the metrics that do not need the FM-index\n"
                "\n"
                ) % PACKAGE_NAME % kPreQCThreads << std::endl;

        return 256;
    }

    static PreQC _runner;
};

static const std::string shortopts = "c:s:o:t:h";
enum { OPT_HELP = 1, OPT_SIMPLE };
static const option longopts[] = {
    {"prefix",              required_argument,  NULL, 'o'}, 
    {"threads",             required_argument,  NULL, 't'}, 
    {"simple",              required_argument,  NULL, OPT_SIMPLE}, 
    {"help",                no_argument,        NULL, 'h'}, 
    {NULL, 0, NULL, 0}, 
};
PreQC PreQC::_runner(
        "preqc", 
        "preform pre-assembly quality checks", 
        shortopts,  
        longopts 
        );
