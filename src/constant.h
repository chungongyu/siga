#ifndef constant_h_
#define constant_h_

#include <cstdio>

#define SAI_EXT   ".sai"
#define RSAI_EXT  ".rsai"
#define BWT_EXT   ".bwt"
#define RBWT_EXT  ".rbwt"
#define ASQG_EXT  ".asqg"
#define HITS_EXT  ".hits"
#define GZIP_EXT  ".gz"
#define RMDUP_EXT ".rmdup"
#define EC_EXT    ".ec"
#define FA_EXT    ".fa"

// common
extern const char* kLogConfig;
extern const char* kWorkDir;
extern const size_t kKmerSize;

// preprocess
extern const size_t kBuckets;

// gap_filling
extern const int MAX_MERGE_GAP;
extern const int MIN_EDGE_COUNT_FOR_TRAINING;
extern const int MAX_CHOICE;
extern const int CANDI_THRESHOLD;

// comand sorting
enum {
    kPreprocess, 
    kIndex, 
    kCorrect, 
    kOverlap, 
    kAssemble, 
    kRmDup, 
    kScaffold, 
    kGapFill,
    kPreQC, 
    kUnkown = 1000
};

#endif // constant_h_
