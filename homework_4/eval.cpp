#include <iostream>
#include <vector>
#include <armadillo>
#include <math.h>

#include "eval.h"
#include "STO3G.h"
#include "Atom.h"

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
                                                               pow((coord_p - coord_a),(l_a - i)) *
                                                               pow((coord_p - coord_b),(l_b - j))) /
                                                              pow((2 * (exp_a + exp_b)),((i + j) * 0.5)));
            }
        }
    }
    return prefactor * sqrt((M_PI) / (exp_a + exp_b)) * sum;
}

//returns overlap integration with analytical integration in 3D
double overlap_3D_ana_int(double exp_a, double exp_b, vector<int> l_a, vector<int> l_b, vec coord_a,
                          vec coord_b){
    double total_overlap = 1;

    for(int i = 0; i < 3; i++){
        total_overlap *= overlap_1D_ana_int(exp_a, exp_b, l_a[i], l_b[i],
                                            coord_a(i), coord_b(i));
    }
    return total_overlap;
}


//calculates Boys function
double boysFunction(double alpha_k, double alpha_k_p, double beta_k, double beta_k_p, double distance){
    double sigma_A = pow((alpha_k + alpha_k_p), -1.0);
    double sigma_B = pow((beta_k + beta_k_p), -1.0);
    double U_A = pow((M_PI * sigma_A), 1.5);
    double U_B = pow((M_PI * sigma_B), 1.5);
    double V_square = pow((sigma_A + sigma_B), -1.0);
    double T = V_square * distance;

    if(T == 0){
        return 27.211 * U_A * U_B * sqrt(2 * V_square) * sqrt(2.0 / M_PI);
    }
    return 27.211 * U_A * U_B * sqrt(pow(distance, -1.0)) * erf(sqrt(T));
}

//calculate gamma element
double gamma_element(Atom A, Atom B){
    double total = 0;

    for(int i = 0; i < 3; i++){ //S orbital always in the [0] index
        double d_p_k = A.AOs[0].primitive[i].d * A.AOs[0].primitive[i].N;
        for(int j = 0; j < 3; j++){
            double d_p_k_p = A.AOs[0].primitive[j].d * A.AOs[0].primitive[j].N;
            for(int k = 0; k < 3; k++){
                double d_p_l = B.AOs[0].primitive[k].d * B.AOs[0].primitive[k].N;
                for(int l = 0; l < 3; l++){
                    double d_p_l_p = B.AOs[0].primitive[l].d * B.AOs[0].primitive[l].N;

                    total += d_p_k * d_p_k_p * d_p_l * d_p_l_p * boysFunction(
                            A.AOs[0].primitive[i].exp, A.AOs[0].primitive[j].exp,
                            B.AOs[0].primitive[k].exp, B.AOs[0].primitive[l].exp,
                            pow(norm(A.AOs[0].coord - B.AOs[0].coord), 2));
                }
            }
        }
    }

    return total;
}

//construct gamma matrix
void gamma_matrix(mat& gamma, vector<Atom>& molecule){
    for(int i = 0; i < molecule.size(); i++){
        for(int j = 0; j < molecule.size(); j++){
            gamma(i, j) = gamma_element(molecule[i], molecule[j]);
        }
    }
}


//calculates two STO3G overlap
double STO3G_overlap(STO3G A, STO3G B){
    double overlap = 0.0;

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            overlap += A.primitive[i].d * B.primitive[j].d * A.primitive[i].N * B.primitive[j].N *
                       overlap_3D_ana_int(A.primitive[i].exp, B.primitive[j].exp, A.L_vec, B.L_vec,
                                          A.coord, B.coord);
        }
    }
    return overlap;
}

//construct overlap matrix with basis set
void overlap_matrix(mat& S, Basis& basis){
    for(int i = 0; i < basis.basis.size(); i++){
        for(int j = 0; j < basis.basis.size(); j++){
            S(i, j) = STO3G_overlap(basis.basis[i], basis.basis[j]);
        }
    }
}

//construct core hamiltonian matrix
void core_hamiltonian_matrix(mat& H, const mat& S, const Basis& basis, const vector<Atom>& molecule){
    double term_1, term_2, term_3;

    for(int i = 0; i < basis.basis.size(); i++){
        term_1 = 0.0, term_2 = 0.0, term_3 = 0.0;

        Atom A = molecule[basis.atom_index[i]];
        for(int j = 0; j < basis.basis.size(); j++){
            if(i == j){
                //terms 1 and 2 calculation
                term_1 = basis.basis[i].IA; // 0.5 not in the first term???
                term_2 = (A.AOs[0].Z - 0.5) * gamma_element(A, A);

                //term 3 calculation
                for(int k = 0; k < molecule.size(); k++){
                    Atom B = molecule[k];
                    term_3 += B.AOs[0].Z * gamma_element(A, B);
                }
                term_3 -= A.AOs[0].Z * gamma_element(A, A);

                H(i, i) = -term_1 - term_2 - term_3;
            } else{
                H(i, j) = 0.5 * (basis.basis[i].beta + basis.basis[j].beta) * S(i, j);
            }
        }
    }
}

//Returns eigenvectors of a symmetric matrix
mat diagonalize_coeff(const mat& A){
    vec lambda;
    mat X;
    eig_sym(lambda, X, A);

    return X;
}

//Returns eigenvalues of a symmetric matrix
vec diagonalize_eigen(const mat& A){
    vec lambda;
    mat X;
    eig_sym(lambda, X, A);

    return lambda;
}

//construct density matrix from occupied MO coefficient matrix
void density_matrix(mat &P, const mat& C, const int pq){
    mat C_occ(C.n_rows, pq);

    for(int i = 0; i < C.n_rows; i++){
        for(int j = 0; j < pq; j++){
            C_occ(i, j) = C(i, j);
        }
    }

    P = C_occ * C_occ.t();
}

//construct vector of densities of atom overlap
void density_atom_overlap(vec& P_atom_overlap, const mat& P_total, const vector<Atom>& molecule,
                          const Basis& basis){
    double total = 0.0;

    for(int i = 0; i < molecule.size(); i++){
        total = 0.0;
        for(int j = 0; j < basis.basis.size(); j++){
            if(basis.atom_index[j] == i){
                total += P_total(j, j);
            }
        }
        P_atom_overlap[i] = total;
    }
}

//construct fock matrix
void fock_matrix(mat& F, const mat& P_alpha, const mat& P_beta, const mat& S, const Basis& basis,
                 const vector<Atom>& molecule){
    double term_1, term_2, term_3;

    mat P_total = P_alpha + P_beta;
    vec P_atom_overlap(molecule.size());
    density_atom_overlap(P_atom_overlap, P_total, molecule, basis);

    for(int i = 0; i < basis.basis.size(); i++){
        term_1 = 0.0, term_2 = 0.0, term_3 = 0.0;

        Atom A = molecule[basis.atom_index[i]];
        for(int j = 0; j < basis.basis.size(); j++){
            if(i == j){
                //terms 1 and 2 calculation
                term_1 = basis.basis[i].IA; // 0.5 not in the first term???
                term_2 = ((P_atom_overlap(basis.atom_index[i]) - A.AOs[0].Z) -
                        (P_alpha(i, i) - 0.5)) * gamma_element(A, A);

                //term 3 calculation
                for(int k = 0; k < molecule.size(); k++){
                    Atom B = molecule[k];
                    term_3 += (P_atom_overlap(k) - B.AOs[0].Z) * gamma_element(A, B);
                }
                term_3 -= (P_atom_overlap(basis.atom_index[i]) - A.AOs[0].Z) * gamma_element(A, A);

                F(i, i) = -term_1 + term_2 + term_3;
            } else{
                Atom B = molecule[basis.atom_index[j]];
                term_1 = 0.5 * (basis.basis[i].beta + basis.basis[j].beta) * S(i, j);
                term_2 = P_alpha(i, j) * gamma_element(A, B);
                F(i, j) = term_1 - term_2;

            }
        }
    }
}

//SCF convergence algorithm, return Fock matrix and density matrix
void SCF(mat& F_A, mat& F_B, mat& P_A, mat& P_B, const mat& P_A_i, const mat& P_B_i, const mat& S, const Basis& basis,
         const vector<Atom>& molecule){

    //initial condition setup
    mat P_A_new(basis.basis.size(), basis.basis.size()), P_B_new(basis.basis.size(), basis.basis.size()),
    P_A_old(basis.basis.size(), basis.basis.size()), P_B_old(basis.basis.size(), basis.basis.size()),
    P_total(basis.basis.size(), basis.basis.size()),
    F_A_new(basis.basis.size(), basis.basis.size()), F_B_new(basis.basis.size(), basis.basis.size()),
    F_A_old(basis.basis.size(), basis.basis.size()), F_B_old(basis.basis.size(), basis.basis.size());

    vec P_atom_overlap(molecule.size());

    mat C_A, C_B, E_A(basis.basis.size(), basis.basis.size()), E_B(basis.basis.size(), basis.basis.size());

    P_A_new = P_A_i;
    P_B_new = P_B_i;

    bool converge = false;
    int i = 0;

    while(!converge){
        cout << '\n' << '\n' << "iteration: " << i << endl;

        fock_matrix(F_A_new, P_A_new, P_B_new, S, basis, molecule);
        fock_matrix(F_B_new, P_B_new, P_A_new, S, basis, molecule);

        cout << "F_alpha" << '\n' << F_A_new << endl;
        cout << "F_beta" << '\n' << F_B_new << endl;

        C_A = diagonalize_coeff(F_A_new);
        C_B = diagonalize_coeff(F_B_new);

        cout << "C_alpha" << '\n' << C_A << endl;
        cout << "C_beta" << '\n' << C_B << endl;

        P_A_old = P_A_new;
        P_B_old = P_B_new;

        density_matrix(P_A_new, C_A, basis.p);
        density_matrix(P_B_new, C_B, basis.q);
        P_total = P_A_new + P_B_new;
        density_atom_overlap(P_atom_overlap, P_total, molecule, basis);

        cout << "P_alpha" << '\n' << P_A_new << endl;
        cout << "P_beta" << '\n' << P_B_new << endl;
        cout << "P_total" << '\n' << P_total << endl;
        cout << "P_atom_overlap" << '\n' << P_atom_overlap << endl;


        F_A_old = F_A_new;
        F_B_old = F_B_new;

        if(approx_equal(P_A_new,P_A_old,"absdiff",1e-6)){
            converge = true;
        }

        i++;
    }
    E_A = diagonalize_eigen(F_A_new);
    E_B = diagonalize_eigen(F_B_new);

    cout << '\n' << "Final Convergence: " << endl;
    cout << "E_alpha" << '\n' << E_A << endl;
    cout << "E_beta" << '\n' << E_B << endl;
    cout << "C_alpha" << '\n' << C_A << endl;
    cout << "C_beta" << '\n' << C_B << endl;

    F_A = F_A_new;
    F_B = F_B_new;
    P_A = P_A_new;
    P_B = P_B_new;
}

//returns electron energy term of energy
double electron_E(const mat& H, const mat& F_A, const mat& F_B, const mat& P_A, const mat& P_B, Basis& basis){
    double alpha_E = 0.0, beta_E = 0.0;

    for(int i = 0; i < basis.basis.size(); i++){
        for(int j = 0; j < basis.basis.size(); j++){
            alpha_E += P_A(i, j) * (H(i, j) + F_A(i, j));
            beta_E += P_B(i, j) * (H(i, j) + F_B(i, j));
        }
    }
    return 0.5 * (alpha_E + beta_E); //27.211 eV/A.U. conversion not needed?
}

//return nuclear repulsion term of energy
double nuc_repulsion_E(vector<Atom>& molecule){
    double E = 0.0;

    for(int i = 0; i < molecule.size(); i++){
        Atom A = molecule[i];
        for(int j = 0; j < molecule.size(); j++){
            Atom B = molecule[j];
            if(i != j){
                E += (A.AOs[0].Z * B.AOs[0].Z) / (norm(A.AOs[0].coord - B.AOs[0].coord));
            }
        }
    }
    return E * 27.211 * 0.5;
}