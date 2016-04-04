#ifndef overlap_process_h_
#define overlap_process_h_

#include "overlap_builder.h"
#include "sequence_process_framework.h"

#include <iostream>

//
// OverlapProcess - Compute the overlap blocks for reads
//
class OverlapProcess {
public:
    OverlapProcess(OverlapBuilder* builder, size_t minOverlap, std::ostream& stream) : _builder(builder), _minOverlap(minOverlap), _stream(stream) {
    }

    OverlapResult process(const SequenceProcessFramework::SequenceWorkItem& workItem);

private:
    OverlapBuilder* _builder;
    size_t _minOverlap;
    std::ostream& _stream;
};

#endif // overlap_process_h_
