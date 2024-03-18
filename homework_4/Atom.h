#ifndef HOMEWORK_4_ATOM_H
#define HOMEWORK_4_ATOM_H

#include <iostream>
#include <vector>
#include <armadillo>

#include "GTO.h"
#include "STO3G.h"

struct Atom{
    vector<STO3G> AOs;
    vec coord;
    int atomic_number;
};

struct Basis{
    vector<STO3G> basis;
    int n;
    int p;
    int q;
};


void construct_Atoms(vector<Atom>& atoms, const string& mol_file);

void construct_basis(Basis& basis, const vector<Atom>& Atoms);

#endif //HOMEWORK_4_ATOM_H
