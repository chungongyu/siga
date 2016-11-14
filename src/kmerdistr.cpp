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
            std::string w = s.substr(j - k, k);
            std::string v = make_reverse_dna_copy(w);

            size_t count = 0;
            count += FMIndex::Interval::occurrences(w, index);
            count += FMIndex::Interval::occurrences(v, index);

            if (distr != NULL) {
                distr->add(count);
            }
        }
        
        L += s.length();
    }

    return L;
}
