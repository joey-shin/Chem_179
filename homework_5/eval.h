#ifndef HOMEWORK_5_EVAL_H
#define HOMEWORK_5_EVAL_H

#include <iostream>
#include <vector>
#include <armadillo>
#include <math.h>

#include "STO3G.h"
#include "Atom.h"

using namespace std;
using namespace arma;


double overlap_3D_ana_int(double exp_a, double exp_b, vector<int> l_a, vector<int> l_b, vec coord_a,
                          vec coord_b);

void gamma_matrix(mat& gamma, vector<Atom>& molecule);

void overlap_matrix(mat& S, Basis& basis);

void core_hamiltonian_matrix(mat& H, const mat& S, const Basis& basis, const vector<Atom>& molecule);

//mat coeff_occ(const mat& C, const int pq);

void density_matrix(mat &P, const mat& C, const int pq);

void fock_matrix(mat& F, const mat& P_alpha, const mat& P_beta, const mat& S, const Basis& basis,
                 const vector<Atom>& molecule);

void SCF(mat& F_A, mat& F_B, mat& P_A, mat& P_B, const mat& P_A_i, const mat& P_B_i, const mat& S, const Basis& basis,
         const vector<Atom>& molecule);

double electron_E(const mat& H, const mat& F_A, const mat& F_B, const mat& P_A, const mat& P_B, Basis& basis);

double nuc_repulsion_E(vector<Atom>& molecule);

void V_nuc_R_matrix(mat& V_nuc_R, const vector<Atom>& molecule);

void X_matrix(mat& X, const Basis& basis, const mat& P_A, const mat& P_B);

void gamma_R_matrix(mat& gamma_R, const vector<Atom>& molecule);

void S_R_matrix(mat& S_R, const Basis& basis);

void Y_matrix(mat& Y, const vector<Atom>& molecule, const Basis& basis, const mat& P_A, const mat& P_B);


#endif //HOMEWORK_5_EVAL_H
