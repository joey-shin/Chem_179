#ifndef HOMEWORK_4_GTO_1_H
#define HOMEWORK_4_GTO_1_H

#include <iostream>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

struct GTO{
    double exp;
    double d; //contraction coefficient corresponding to GTO
    double N; //Normalization constant corresponding to GTO
};

void read_STO3G(vector<GTO>& primitive, const int& total_L, string& basis_file_name);

void normalize_GTO(vector<GTO>& primitive, vector<int>& L, const vec& coord);

#endif //HOMEWORK_4_GTO_1_H
