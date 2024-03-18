#ifndef HOMEWORK_4_STO3G_1_H
#define HOMEWORK_4_STO3G_1_H

#include <iostream>
#include <vector>
#include <armadillo>

#include "GTO.h"

struct STO3G{
    vector<GTO> primitive;
    vector<int> L_vec; //angular quantum number vector
    vec coord;
    int total_L;
    double IA;
    double beta;
    int Z;
};

void construct_AO(STO3G& AO, int& atomic_number, int& total_L, const vec& coord);

#endif //HOMEWORK_4_STO3G_1_H
