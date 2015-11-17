#ifndef fmindex_h_
#define fmindex_h_

#include "bwt.h"

#include <iostream>

class FMIndex {
public:
    FMIndex() {
    }
private:
    friend std::ostream& operator<<(std::ostream& stream, const FMIndex& index);
    friend std::istream& operator>>(std::istream& stream, FMIndex& index);

    BWT _bwt;
};

#endif // fmindex_h_
