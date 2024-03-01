#ifndef HOMEWORK_3_GTO_H
#define HOMEWORK_3_GTO_H

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

void read_STO3G(vector<GTO>& primitive, const int& atomic_number, const int& total_L,
                const string& H_STO3G_file, const string& C_STO3G_file);

void normalize_GTO(vector<GTO>& primitive, vector<int>& L, vector<double>& coord);


#endif //HOMEWORK_3_GTO_H
