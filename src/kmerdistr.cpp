#include "kmerdistr.h"

#include <string>

size_t KmerDistribution::sample(const FMIndex* index, size_t k, size_t n, KmerDistribution* distr) {
    size_t L = 0;
    size_t N = index->length();

    // Learn k-mer occurrence distribution for this value of k
    for (size_t i = 0; i < n; ++i) {
        size_t idx = rand() % N;
        std::string s = index->getString(idx);
        if (s.length() < k) {
            continue;
        }

        for (size_t j = k; j < s.length(); ++j) {
            std::string kmer = s.substr(j - k, k);
            size_t count = 0;
            if (distr != NULL) {
                distr->add(count);
            }
        }
        
        L += s.length();
    }

    return L;
}
