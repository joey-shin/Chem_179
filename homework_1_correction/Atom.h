#ifndef HOMEWORK_1_CORRECTION_ATOM_H
#define HOMEWORK_1_CORRECTION_ATOM_H

#include <iostream>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

struct Atom{
    int atomic_number;
    vec coord(3);
};

void construct_Atoms(const vector<Atom>& Atoms, const string& file);

#endif //HOMEWORK_1_CORRECTION_ATOM_H
