#ifndef HOMEWORK_5_STO3G_H
#define HOMEWORK_5_STO3G_H

#include <iostream>
#include <vector>
#include <armadillo>

#include "GTO.h"

struct STO3G{
    vector<GTO> primitive;
    vector<int> L_vec; //angular quantum number vector
    vec coord;
    int atomic_number;
    int total_L;
    double IA;
    double beta;
    int Z;
};

void construct_AO(STO3G& AO, int& atomic_number, int& total_L, const vec& coord);

#endif //HOMEWORK_5_STO3G_H
