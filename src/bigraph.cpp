#include "bigraph.h"

Bigraph::~Bigraph() {
    for (VertaexTable::iterator i = _vertices.begin(); i != _vertices.end(); ++i) {
        delete i->second;
    }
}
