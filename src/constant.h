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
#define BZIP_EXT  ".bz2"
#define RMDUP_EXT ".rmdup"
#define EC_EXT  ".ec"
#define FA_EXT  ".fa"

// command sorting
enum {
  kPreprocess,
  kIndex,
  kCorrect,
  kOverlap,
  kAssemble,
  kSubgraph,
  kRmDup,
  kPreQC,
  kMatch,
  kUnkown = 1000
};

#endif  // constant_h_
