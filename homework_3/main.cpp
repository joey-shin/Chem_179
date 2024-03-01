#include <iostream>
#include <armadillo>

#include "GTO.h"
#include "STO3G.h"
#include "eval.h"

//console input structure: ./[executable file] basis/H_STO3G.txt basis/C_STO3G.txt sample_input/[molecule].txt
int main(int argc, char* argv[]) {

    vector<STO3G> basis;
    int n = 0; //number of electron pairs

    //constructs basis
    construct_basis(basis, argv[1], argv[2], argv[3], n);

    cout << "number of basis: " << basis.size() << endl;
    cout << "number of electrons pairs: " << n << '\n' << endl;

    for(int i = 0; i < basis.size(); i++){
        cout << '\n' << "atomic number: " << basis[i].atomic_number << '\n';
        for(int j = 0; j < 3; j++){
            cout << "basis [" << i << "] GTO [" << j << "] exp: " << basis[i].primitive[j].exp << endl;
            cout << "basis [" << i << "] GTO [" << j << "] contract coeff: " << basis[i].primitive[j].d << endl;
            cout << "basis [" << i << "] GTO [" << j << "] normalization constant: " << basis[i].primitive[j].N << endl;
        }
        cout << "L: ";
        for(int j = 0; j < 3; j++){
            cout << basis[i].L[j] << ' ';
        }
        cout << endl;

        cout << "Coord: ";
        for(int j = 0; j < 3; j++){
            cout << basis[i].coord[j] << ' ';
        }
        cout << '\n' << '\n' << endl;
    }


    mat S(basis.size(), basis.size());
    overlap_matrix(S, basis);

    cout << "Overlap Matrix: " << endl;
    cout << S << endl;


    mat H(basis.size(), basis.size());
    hamiltonian_matrix(H, S, basis);

    cout << "Hamiltonian Matrix: " << endl;
    cout << H << endl;


    mat X(basis.size(), basis.size());
    orthogonalization_matrix(X, S);

    cout << "Orthogonalization Matrix: " << endl;
    cout << X << endl;

    mat C(basis.size(), basis.size());
    double E = 0.0;
    coefficient_matrix_energy(H, X, C, E, n);

    cout << "Coefficient Matrix: " << endl;
    cout << C << endl;
    cout << "Total Energy: " << E << endl;

    return 0;
}
