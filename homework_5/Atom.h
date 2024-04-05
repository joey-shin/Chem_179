#ifndef HOMEWORK_5_ATOM_H
#define HOMEWORK_5_ATOM_H

#include <iostream>
#include <vector>
#include <armadillo>

#include "GTO.h"
#include "STO3G.h"

//main structure, contains Atomic orbitals
struct Atom{
    vector<STO3G> AOs;
};

//expansion of every STO3G in every Atom into vector for convenience
struct Basis{
    vector<STO3G> basis;
    vector<int> atom_index;
    int n;
    int p;
    int q;
};


void construct_Atoms(vector<Atom>& atoms, const string& mol_file);

void construct_basis(Basis& basis, const vector<Atom>& Atoms);

#endif //HOMEWORK_5_ATOM_H
