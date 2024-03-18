#ifndef HOMEWORK_3_BASIS_H
#define HOMEWORK_3_BASIS_H

#include <iostream>
#include <vector>
#include <armadillo>

#include "GTO.h"

struct STO3G{
    vector<GTO> primitive;
    vector<int> L; //angular quantum number vector
    int total_L;
    vector<double> coord;
    int atomic_number;
};

void construct_basis(vector<STO3G>& basis, const string& H_STO3G_file, const string& C_STO3G_file, const string& mol_file, int& n);

#endif //HOMEWORK_3_BASIS_H