#ifndef assemble_h_
#define assemble_h_

#include "runner.h"

class Assembler : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments);

private:
    Assembler();
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static Assembler _runner;
};

#endif // assemble_h_
