#include "overlap_process.h"

//
// OverlapProcess
//
OverlapResult OverlapProcess::process(const SequenceProcessFramework::SequenceWorkItem& workItem) {
    OverlapBlockList blocks;
    OverlapResult result = _builder->overlap(workItem.read, _minOverlap, &blocks);
    return result;
}
