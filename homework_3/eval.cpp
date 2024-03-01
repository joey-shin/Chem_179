#include <iostream>
#include <vector>
#include <armadillo>
#include <math.h>

#include "eval.h"
#include "STO3G.h"

using namespace std;
using namespace arma;

int factorial(int l){
    int total = 1;

    for(int i = l; i > 0; i--){
        total *= i;
    }
    return total;
}

int double_factorial(int l){
    int total = 1;

    if((l % 2) == 0){
        for(int i = l; i > 0; i--){
            if((i % 2) == 0){
                total *= i;
            }
        }
    } else{
        for(int i = l; i > 0; i--){
            if((i % 2) == 1){
                total *= i;
            }
        }
    }
    return total;
}

int binomial(int m, int n){
    return factorial(m) / (factorial(n) * factorial(m - n));
}


//returns overlap integration with analytical integration in 1D
double overlap_1D_ana_int(double exp_a, double exp_b, int l_a, int l_b, double coord_a, double coord_b){
    double sum = 0.0;
    double coord_p = ((exp_a * coord_a) + (exp_b * coord_b)) / (exp_a + exp_b);
    double prefactor = exp((-1 * exp_a * exp_b * pow((coord_a - coord_b),2)) / (exp_a + exp_b));

    for(int i = 0; i <= l_a; i++){
        for(int j = 0; j <= l_b; j++){
            if(((i + j) % 2) == 0){
                sum += binomial(l_a, i) * binomial(l_b, j) * ((double_factorial(i + j - 1) *
                                                               pow((coord_p - coord_a),(l_a - i)) * pow((coord_p - coord_b),(l_b - j))) /
                                                              pow((2 * (exp_a + exp_b)),((i + j) * 0.5)));
            }
        }
    }
    return prefactor * sqrt((M_PI) / (exp_a + exp_b)) * sum;
}

//returns overlap integration with analytical integration in 3D
double overlap_3D_ana_int(double exp_a, double exp_b, vector<int> l_a, vector<int> l_b, vector<double> coord_a, vector<double> coord_b){
    double total_overlap = 1;

    for(int i = 0; i < 3; i++){
        total_overlap *= overlap_1D_ana_int(exp_a, exp_b, l_a[i], l_b[i], coord_a[i], coord_b[i]);
    }
    return total_overlap;
}

double STO3G_overlap(STO3G A, STO3G B){
    double overlap = 0.0;

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            overlap += A.primitive[i].d * B.primitive[j].d * A.primitive[i].N * B.primitive[j].N *
                    overlap_3D_ana_int(A.primitive[i].exp, B.primitive[j].exp, A.L, B.L,
                                       A.coord, B.coord);
        }
    }
    return overlap;
}


//construct overlap matrix
void overlap_matrix(mat& S, vector<STO3G>& basis){
    for(int i = 0; i < basis.size(); i++){
        for(int j = 0; j < basis.size(); j++){
            S(i, j) = STO3G_overlap(basis[i], basis[j]);
        }
    }
}


//returns diagonals of the hamiltonian matrix
double return_hamiltonian(STO3G basis){
    double h_H1s = -13.6;
    double h_C2s = -21.4;
    double h_C2p = -11.4;

    if(basis.atomic_number == 1){
        return h_H1s;
    } else if(basis.atomic_number == 6){
        if(basis.total_L == 0){
            return h_C2s;
        } else if(basis.total_L == 1){
            return h_C2p;
        }
    }
    throw runtime_error("molecule is not a hydrocarbon. ");
}


//construct hamiltonian matrix
void hamiltonian_matrix(mat& H, const mat& S, vector<STO3G>& basis){
    double K = 1.75;

    for(int i = 0; i < S.n_rows; i++){
        for(int j = 0; j < S.n_cols; j++){
            //setting diagonals of matrix
            if(i == j){
                H(i, j) = return_hamiltonian(basis[i]);
            } else{ //setting off diagonals
                H(i, j) = 0.5 * K * (return_hamiltonian(basis[i]) + return_hamiltonian(basis[j])) * S(i, j);
            }
        }
    }
}


//construct orthogonalizing matrix
void orthogonalization_matrix(mat& X, mat& S){
    mat S_evec, S_eval_inv_sqrt;
    vec S_eval;

    eig_sym(S_eval, S_evec, S);
    S_eval_inv_sqrt = inv(diagmat(sqrt(S_eval)));
    X = S_evec * S_eval_inv_sqrt * S_evec.t();
}


//construct coefficient matrix and total molecule energy
void coefficient_matrix_energy(mat& H, mat& X, mat& C, double& E, int& n){
    //form Hamltonian in orthogonal basis
    mat H_ortho = X.t() * H * X;

    //
    mat U;
    vec E_eval;
    eig_sym(E_eval, U, H_ortho);

    C = X * U;

    for(int i = 0; i < n; i++){
        E += 2 * E_eval(i);
    }
}