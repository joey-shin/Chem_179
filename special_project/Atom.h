#ifndef SPECIAL_PROJECT_ATOM_H
#define SPECIAL_PROJECT_ATOM_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

struct Atom{
    int atomic_number;
    vec coord;
};

void construct_Atoms(vector<Atom>& Atoms, const string& file);
void print_Atoms(const vector<Atom>& Atoms);


#endif //SPECIAL_PROJECT_ATOM_H
