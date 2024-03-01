#ifndef HOMEWORK_3_EVAL_H
#define HOMEWORK_3_EVAL_H

#include <iostream>
#include <vector>
#include <armadillo>
#include <math.h>

#include "STO3G.h"

using namespace std;
using namespace arma;

double overlap_3D_ana_int(double exp_a, double exp_b, vector<int> l_a, vector<int> l_b, vector<double> coord_a, vector<double> coord_b);

void overlap_matrix(mat& S, vector<STO3G>& basis);

void hamiltonian_matrix(mat& H, const mat& S, vector<STO3G>& basis);

void orthogonalization_matrix(mat& X, mat& S);

void coefficient_matrix_energy(mat& H, mat& X, mat& C, double& E, int& n);



#endif //HOMEWORK_3_EVAL_H
