// ssw_cpp.h
// Created by Wan-Ping Lee
// Last revision by Mengyao Zhao on 2017-05-30

#ifndef COMPLETE_STRIPED_SMITH_WATERMAN_CPP_H_
#define COMPLETE_STRIPED_SMITH_WATERMAN_CPP_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace StripedSmithWaterman {

struct Alignment {
  uint16_t sw_score;           // The best alignment score
  uint16_t sw_score_next_best; // The next best alignment score
  size_t  ref_begin;          // Reference begin position of the best alignment
  size_t  ref_end;            // Reference end position of the best alignment
  size_t  query_begin;        // Query begin position of the best alignment
  size_t  query_end;          // Query end position of the best alignment
  size_t  ref_end_next_best;  // Reference end position of the next best alignment
  size_t  mismatches;         // Number of mismatches of the alignment
  std::string cigar_string;    // Cigar string of the best alignment
  std::vector<uint32_t> cigar; // Cigar stored in the BAM format
                               //   high 28 bits: length
			       //   low 4 bits: M/I/D/S/X (0/1/2/4/8);
  void Clear() {
    sw_score           = 0;
    sw_score_next_best = 0;
    ref_begin          = 0;
    ref_end            = 0;
    query_begin        = 0;
    query_end          = 0;
    ref_end_next_best  = 0;
    mismatches         = 0;
    cigar_string.clear();
    cigar.clear();
  };
};

struct Filter {
  // NOTE: No matter the filter, those five fields of Alignment will be given anyway.
  //       sw_score; sw_score_next_best; ref_end; query_end; ref_end_next_best.
  // NOTE: Only need score of alignments, please set 'report_begin_position'
  //       and 'report_cigar' false.

  bool report_begin_position;    // Give ref_begin and query_begin.
                                 //   If it is not set, ref_begin and query_begin are -1.
  bool report_cigar;             // Give cigar_string and cigar.
                                 //   report_begin_position is automatically TRUE.

  // When *report_cigar* is true and alignment passes these two filters,
  //   cigar_string and cigar will be given.
  uint16_t score_filter;         // score >= score_filter
  uint16_t distance_filter;      // ((ref_end - ref_begin) < distance_filter) &&
                                 // ((query_end - read_begin) < distance_filter)

  Filter()
    : report_begin_position(true)
    , report_cigar(true)
    , score_filter(0)
    , distance_filter(32767)
  {}

  Filter(const bool& pos, const bool& cigar, const uint16_t& score, const uint16_t& dis)
    : report_begin_position(pos)
    , report_cigar(cigar)
    , score_filter(score)
    , distance_filter(dis)
    {}
};

class Aligner {
 public:
  // =========
  // @function Construct an Aligner by assigning ref and scores.
  //             The function will build the {A.C,G,T,N} aligner.
  //             If you target for other character aligners, then please
  //             use the other constructor and pass the corresponding matrix in.
  // =========
  Aligner(const std::string& query, uint8_t match_score = 2, uint8_t mismatch_penalty = 2,
	  uint8_t gap_opening_penalty = 3, uint8_t gap_extending_penalty = 1);

  virtual ~Aligner();

  // =========
  // @function Align target to query.
  //             SetReferenceSequence.
  // @param    target    The target sequence.
  // @param    filter    The filter for the alignment.
  // @param    alignment The container contains the result.
  // @param    maskLen   The distance between the optimal and suboptimal alignment ending position will >= maskLen. We suggest to 
  //                     use readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function 
  //                     will NOT return the suboptimal alignment information.
  // @return   True: succeed; false: fail.
  // =========
  bool Align(const std::string& target, const Filter& filter, int32_t maskLen, Alignment* alignment) const;
  bool Align(const std::string& target, int32_t maskLen, Alignment* alignment) const;
  bool Align(const std::string& target, Alignment* alignment) const;

 private:
  void TranslateBase(const std::string& seq, std::vector<int8_t>* translated) const;
  std::vector<int8_t>* TranslateBase(const std::string& seq) const;

  std::shared_ptr<std::vector<int8_t> > translated_query_;
  std::shared_ptr<std::vector<int8_t> > score_matrix_;
  uint8_t gap_opening_penalty_;   // default: 3
  uint8_t gap_extending_penalty_; // default: 1
  std::shared_ptr<void> handler_;
}; // class Aligner

} // namespace StripedSmithWaterman

#endif // COMPLETE_STRIPED_SMITH_WATERMAN_CPP_H_
