#define BOOST_TEST_MODULE "siga.index"
#include <boost/test/included/unit_test.hpp>
#include <boost/format.hpp>

#include <ssw_cpp.h>

#include "alphabet.h"
#include "rlstring.h"
#include "suffix_array.h"

BOOST_AUTO_TEST_SUITE(Indexer);

BOOST_AUTO_TEST_CASE(Alphabet_torank) {
  BOOST_CHECK_EQUAL(DNAAlphabet::torank('$'), 0);
  BOOST_CHECK_EQUAL(DNAAlphabet::torank('A'), 1);
  BOOST_CHECK_EQUAL(DNAAlphabet::torank('C'), 2);
  BOOST_CHECK_EQUAL(DNAAlphabet::torank('G'), 3);
  BOOST_CHECK_EQUAL(DNAAlphabet::torank('T'), 4);
}
BOOST_AUTO_TEST_CASE(Alphabet_tochar) {
  BOOST_CHECK_EQUAL(DNAAlphabet::tochar(0), '$');
  BOOST_CHECK_EQUAL(DNAAlphabet::tochar(1), 'A');
  BOOST_CHECK_EQUAL(DNAAlphabet::tochar(2), 'C');
  BOOST_CHECK_EQUAL(DNAAlphabet::tochar(3), 'G');
  BOOST_CHECK_EQUAL(DNAAlphabet::tochar(4), 'T');
}
BOOST_AUTO_TEST_CASE(Alphabet_func) {
  for (size_t i = 0; i < DNAAlphabet::ALL_SIZE; ++i) {
    BOOST_CHECK_EQUAL(i, DNAAlphabet::torank(DNAAlphabet::tochar(i)));
  }
}

BOOST_AUTO_TEST_CASE(RLUnit_test) {
  {
    RLUnit unit;
    BOOST_CHECK(!unit.initialized());
    BOOST_CHECK(unit.empty());
  }
  {
    RLUnit unit('A');
    BOOST_CHECK(!unit.empty());
    BOOST_CHECK_EQUAL(unit.count(), 1);
    BOOST_CHECK_EQUAL((char)unit, 'A');
    ++unit;
    BOOST_CHECK_EQUAL(unit.count(), 2);
    --unit;
    BOOST_CHECK_EQUAL(unit.count(), 1);
  }
}

BOOST_AUTO_TEST_CASE(RLString_test) {
  DNASeq read("test", "AAACGGGTA");

  RLString runs;

  RLUnit run;
  for (const auto& c : read.seq) {
    if (run.initialized()) {
      if (run == c && !run.full()) {
        ++run;
      } else {
        runs.push_back(run);
        run = RLUnit(c);
      }
    } else {
      run = RLUnit(c);
    }
  }
  if (run.initialized()) {
    runs.push_back(run);
  }

  BOOST_CHECK_EQUAL(runs.size(), 5);
  BOOST_CHECK_EQUAL(runs[0], 'A');
  BOOST_CHECK_EQUAL(runs[0].count(), 3);
  BOOST_CHECK_EQUAL(runs[1], 'C');
  BOOST_CHECK_EQUAL(runs[1].count(), 1);
  BOOST_CHECK_EQUAL(runs[2], 'G');
  BOOST_CHECK_EQUAL(runs[2].count(), 3);
  BOOST_CHECK_EQUAL(runs[3], 'T');
  BOOST_CHECK_EQUAL(runs[3].count(), 1);
  BOOST_CHECK_EQUAL(runs[4], 'A');
  BOOST_CHECK_EQUAL(runs[4].count(), 1);
}

SuffixArray::Elem toelem(const std::vector<std::string>& reads, const std::vector<size_t>& anchors, size_t idx) {
  auto it = std::lower_bound(anchors.begin(), anchors.end(), idx);
  assert(it != anchors.end());
  size_t i = std::distance(anchors.begin(), it);
  size_t delta = *it - idx;
  size_t j = reads[i].length() - delta;
  return SuffixArray::Elem(i, j);
}

BOOST_AUTO_TEST_CASE(SSW_test) {
  const std::string ref = "CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAA";
  const std::string query = "CTGAGCCGGTAAATC";

  // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner(query);
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref
  aligner.Align(ref, &alignment);

  BOOST_CHECK_EQUAL(alignment.ref_begin, 8);
  BOOST_CHECK_EQUAL(alignment.ref_end, 21);
  BOOST_CHECK_EQUAL(alignment.query_begin, 0);
  BOOST_CHECK_EQUAL(alignment.query_end, 14);
  BOOST_CHECK_EQUAL(alignment.cigar_string.c_str(), "4=1X4=1I5=");

  std::cout << "===== SSW result =====" << std::endl;
  std::cout << "Best Smith-Waterman score:\t" << alignment.sw_score << std::endl
            << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << std::endl
            << "Reference start:\t" << alignment.ref_begin << std::endl
            << "Reference end:\t" << alignment.ref_end << std::endl
            << "Query start:\t" << alignment.query_begin << std::endl
            << "Query end:\t" << alignment.query_end << std::endl
            << "Next-best reference end:\t" << alignment.ref_end_next_best << std::endl
            << "Number of mismatches:\t" << alignment.mismatches << std::endl
            << "Cigar: " << alignment.cigar_string << std::endl;
  std::cout << "======================" << std::endl;
}

BOOST_AUTO_TEST_CASE(SuffixArray_toelem) {
  std::vector<std::string> reads = {"ACGT", "ACTGN", "CTGANA"};

  size_t suffixes = 0;
  std::vector<size_t> anchors(reads.size());
  std::vector<SuffixArray::Elem> indices;
  for (size_t i = 0; i < reads.size(); ++i) {
    suffixes += reads[i].length() + 1;
    anchors[i] = suffixes - 1;
    for (size_t j = 0; j < reads[i].length() + 1; ++j) {
      indices.push_back(SuffixArray::Elem(i, j));
    }
  }
  BOOST_CHECK_EQUAL(anchors.back() + 1, suffixes);

  std::cout << "reads:" << std::endl;
  std::copy(reads.begin(), reads.end(), std::ostream_iterator<std::string>(std::cout, "\n"));
  std::cout << std::endl;

  std::cout << boost::format("suffixes: %lu") % suffixes << std::endl;

  std::cout << "anchors:" << std::endl;
  std::copy(anchors.begin(), anchors.end(), std::ostream_iterator<size_t>(std::cout, ","));
  std::cout << std::endl;

  for (size_t i = 0; i < suffixes; ++i) {
    SuffixArray::Elem elem = toelem(reads, anchors, i);
    BOOST_CHECK_EQUAL(elem, indices[i]);
    std::cout << boost::format("%lu\t(%lu,%lu)\t(%lu,%lu)") % i % elem.i % elem.j % indices[i].i % indices[i].j << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END();
