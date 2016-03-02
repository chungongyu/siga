//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// PrimerScreen - Singleton class to filter sequences 
// that match a database of primer sequences
//
#ifndef primer_screen_h_
#define primer_screen_h_

#include <string>
#include <vector>

class PrimerScreen {
public:
    // Return true if the sequence fails the primer check
    static bool containsPrimer(const std::string& seq);

private:
    PrimerScreen(); 
    std::vector< std::string > _db;
};

#endif // primer_screen_h_
