// ssw_cpp.cpp
// Created by Wan-Ping Lee
// Last revision by Mengyao Zhao on 2017-05-30

#include "ssw_cpp.h"
#include "ssw.h"

#include <cassert>
#include <algorithm>
#include <sstream>

namespace {

static const int8_t kBaseTranslation[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    // A     C            G
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    //           T
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    // a     c            g
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    //           t
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

std::vector<int8_t>* BuildDNAScoreMatrix(uint8_t match_score, uint8_t mismatch_penalty) {
  // The score matrix looks like
  //                 // A,  C,  G,  T,  N
  //  score_matrix_ = { 2, -2, -2, -2, -2, // A
  //                   -2,  2, -2, -2, -2, // C
  //                   -2, -2,  2, -2, -2, // G
  //                   -2, -2, -2,  2, -2, // T
  //                   -2, -2, -2, -2, -2};// N

  std::vector<int8_t>* matrix = new std::vector<int8_t>(5*5);
  size_t k = 0;
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      (*matrix)[k++] = ((i == j) ? match_score : static_cast<int8_t>(-mismatch_penalty));
    }
    (*matrix)[k++] = static_cast<int8_t>(-mismatch_penalty); // For N
  }

  for (size_t i = 0; i < 5; ++i) {
    (*matrix)[k++] = static_cast<int8_t>(-mismatch_penalty); // For N
  }
  return matrix;
}

void ConvertAlignment(const s_align& s_al,
                      const int& query_len,
                      StripedSmithWaterman::Alignment* al) {
  al->sw_score           = s_al.score1;
  al->sw_score_next_best = s_al.score2;
  al->ref_begin          = s_al.ref_begin1;
  al->ref_end            = s_al.ref_end1;
  al->query_begin        = s_al.read_begin1;
  al->query_end          = s_al.read_end1;
  al->ref_end_next_best  = s_al.ref_end2;

  al->cigar.clear();
  al->cigar_string.clear();

  if (s_al.cigarLen > 0) {
    std::ostringstream cigar_string;
    if (al->query_begin > 0) {
      uint32_t cigar = to_cigar_int(al->query_begin, 'S');
      al->cigar.push_back(cigar);
      cigar_string << al->query_begin << 'S';
    }

    for (int i = 0; i < s_al.cigarLen; ++i) {
      al->cigar.push_back(s_al.cigar[i]);
      cigar_string << cigar_int_to_len(s_al.cigar[i]) << cigar_int_to_op(s_al.cigar[i]);
    }

    int end = query_len - al->query_end - 1;
    if (end > 0) {
      uint32_t cigar = to_cigar_int(end, 'S');
      al->cigar.push_back(cigar);
      cigar_string << end << 'S';
    }

    al->cigar_string = cigar_string.str();
  } // end if
}

// @Function:
//     Calculate the length of the previous cigar operator
//     and store it in new_cigar and new_cigar_string.
//     Clean up in_M (false), in_X (false), length_M (0), and length_X(0).
void CleanPreviousMOperator(
    bool* in_M,
    bool* in_X,
    uint32_t* length_M,
    uint32_t* length_X,
    std::vector<uint32_t>* new_cigar,
    std::ostringstream* new_cigar_string) {
  if (*in_M) {
    uint32_t match = to_cigar_int(*length_M, '=');
    new_cigar->push_back(match);
    (*new_cigar_string) << *length_M << '=';
  } else if (*in_X){ //in_X
    uint32_t match = to_cigar_int(*length_X, 'X');
    new_cigar->push_back(match);
    (*new_cigar_string) << *length_X << 'X';
  }

  // Clean up
  *in_M = false;
  *in_X = false;
  *length_M = 0;
  *length_X = 0;
}

// @Function:
//     1. Calculate the number of mismatches.
//     2. Modify the cigar string:
//         differentiate matches (M) and mismatches(X).
// @Return:
//     The number of mismatches.
int CalculateNumberMismatch(
    StripedSmithWaterman::Alignment* al,
    int8_t const *ref,
    int8_t const *query,
    const int& query_len) {

  ref   += al->ref_begin;
  query += al->query_begin;
  int mismatch_length = 0;

  std::vector<uint32_t> new_cigar;
  std::ostringstream new_cigar_string;

  if (al->query_begin > 0) {
    uint32_t cigar = to_cigar_int(al->query_begin, 'S');
    new_cigar.push_back(cigar);
    new_cigar_string << al->query_begin << 'S';
  }

  bool in_M = false; // the previous is match
  bool in_X = false; // the previous is mismatch
  uint32_t length_M = 0;
  uint32_t length_X = 0;

  for (unsigned int i = 0; i < al->cigar.size(); ++i) {
    char op = cigar_int_to_op(al->cigar[i]);
    uint32_t length = cigar_int_to_len(al->cigar[i]);
    if (op == 'M') {
      for (uint32_t j = 0; j < length; ++j) {
  if (*ref != *query) {
    ++mismatch_length;
          if (in_M) { // the previous is match; however the current one is mismatche
      uint32_t match = to_cigar_int(length_M, '=');
      new_cigar.push_back(match);
      new_cigar_string << length_M << '=';
    }
    length_M = 0;
    ++length_X;
    in_M = false;
    in_X = true;
  } else { // *ref == *query
    if (in_X) { // the previous is mismatch; however the current one is matche
      uint32_t match = to_cigar_int(length_X, 'X');
      new_cigar.push_back(match);
      new_cigar_string << length_X << 'X';
    }
    ++length_M;
    length_X = 0;
    in_M = true;
    in_X = false;
  } // end of if (*ref != *query)
  ++ref;
  ++query;
      }
    } else if (op == 'I') {
      query += length;
      mismatch_length += length;
      CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);
      new_cigar.push_back(al->cigar[i]);
      new_cigar_string << length << 'I';
    } else if (op == 'D') {
      ref += length;
      mismatch_length += length;
      CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);
      new_cigar.push_back(al->cigar[i]);
      new_cigar_string << length << 'D';
    }
  }

  CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);

  int end = query_len - al->query_end - 1;
  if (end > 0) {
    uint32_t cigar = to_cigar_int(end, 'S');
    new_cigar.push_back(cigar);
    new_cigar_string << end << 'S';
  }

  al->cigar_string.clear();
  al->cigar.clear();
  al->cigar_string = new_cigar_string.str();
  al->cigar = new_cigar;

  return mismatch_length;
}

void SetFlag(const StripedSmithWaterman::Filter& filter, uint8_t* flag) {
  if (filter.report_begin_position) *flag |= 0x08;
  if (filter.report_cigar) *flag |= 0x0f;
}

} // namespace


namespace StripedSmithWaterman {

Aligner::Aligner(const std::string& query, uint8_t match_score, uint8_t mismatch_penalty,
        uint8_t gap_opening_penalty, uint8_t gap_extending_penalty)
    : translated_query_(TranslateBase(query))
    , score_matrix_(BuildDNAScoreMatrix(match_score, mismatch_penalty))
    , gap_opening_penalty_(gap_opening_penalty)
    , gap_extending_penalty_(gap_extending_penalty) {
  assert(!translated_query_->empty());
  const int8_t score_size = 2;
  handler_ = std::shared_ptr<void>(ssw_init(&(*translated_query_)[0], translated_query_->size(),
      &(*score_matrix_)[0], 5, score_size), [&](void* h) {
        init_destroy((s_profile *)h);
      });
}

Aligner::~Aligner(){
}

bool Aligner::Align(const std::string& target, const Filter& filter, int32_t maskLen,
                    Alignment* alignment) const {
  if (!handler_ || target.empty()) {
    return false;
  }

  std::vector<int8_t> translated_target;
  TranslateBase(target, &translated_target);

  uint8_t flag = 0;
  SetFlag(filter, &flag);
  s_align* s_al = ssw_align((s_profile *)handler_.get(),
      &translated_target[0], translated_target.size(),
      static_cast<int>(gap_opening_penalty_),
      static_cast<int>(gap_extending_penalty_),
      flag, filter.score_filter, filter.distance_filter, maskLen);

  alignment->Clear();
  ConvertAlignment(*s_al, translated_query_->size(), alignment);
  alignment->mismatches = CalculateNumberMismatch(&*alignment, &translated_target[0], &(*translated_query_)[0], translated_query_->size());

  // Free memory
  align_destroy(s_al);

  return true;
}

bool Aligner::Align(const std::string& target, int32_t maskLen, Alignment* alignment) const {
  static Filter filter;
  return Align(target, filter, maskLen, alignment);
}

bool Aligner::Align(const std::string& target, Alignment* alignment) const {
  int32_t maskLen = std::max(15, (int32_t)translated_query_->size() / 2);
  return Align(target, maskLen, alignment);
}

void Aligner::TranslateBase(const std::string& seq, std::vector<int8_t>* translated) const {
  std::for_each(seq.begin(), seq.end(), [&](char c) {
        translated->push_back(kBaseTranslation[c]);
      });
}

std::vector<int8_t>* Aligner::TranslateBase(const std::string& seq) const {
  std::vector<int8_t>* translated = new std::vector<int8_t>();
  TranslateBase(seq, translated);
  return translated;
}

}  // namespace StripedSmithWaterman
