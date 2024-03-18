#include <iostream>
#include <armadillo>

#include "GTO.h"
#include "STO3G.h"
#include "Atom.h"
#include "eval.h"

int main(int argc, char* argv[]) {

    //construct vector of atoms
    vector<Atom> atoms;
    Basis basis;
    construct_Atoms(atoms, argv[1]);
    construct_basis(basis, atoms);\

    //print all values of Atom structure and sub structures
    cout << "number of Atoms: " << atoms.size() << endl;
    for(int i = 0; i < atoms.size(); i++){
        cout << '\n' << "Atom " << i << '\n';
        cout << '\n' << "atomic number: " << atoms[i].atomic_number << '\n';
        cout << "Coord: ";
        for(int j = 0; j < 3; j++){
            cout << atoms[i].coord[j] << ' ';
        }

        for(int j = 0; j < atoms[i].AOs.size(); j++){
            cout << '\n' << '\n' << "STO3G " << j << '\n';
            cout << "total_L: " << atoms[i].AOs[j].total_L << '\n';
            cout << "IA: " << atoms[i].AOs[j].IA << '\n';
            cout << "beta: " << atoms[i].AOs[j].beta << '\n';
            cout << "Z: " << atoms[i].AOs[j].Z << '\n';
            cout << "L: ";
            for(int k = 0; k < 3; k++){
                cout << atoms[i].AOs[j].L_vec[k] << ' ';
            }
            cout << '\n' << endl;

            for(int k = 0; k < atoms[i].AOs[j].primitive.size(); k++){
                cout << "primitive " << k << '\n';
                cout << "exp " << atoms[i].AOs[j].primitive[k].exp << '\n';
                cout << "d " << atoms[i].AOs[j].primitive[k].d << '\n';
                cout << "N " << atoms[i].AOs[j].primitive[k].N << '\n' << endl;
            }
        }
    }

    //print all values of Basis structure
    cout << "number of basis: " << basis.basis.size() << endl;
    for(int i = 0; i < basis.basis.size(); i++){
        cout << '\n' << "basis " << i << '\n';
        cout << "L: ";
        for(int k = 0; k < 3; k++){
            cout << basis.basis[i].L_vec[k] << ' ';
        }
        cout << '\n' << endl;
    }

    cout  << "n: " << basis.n << endl;
    cout  << "p: " << basis.p << endl;
    cout  << "q: " << basis.q << endl;


    mat gamma(atoms.size(), atoms.size());
    gamma_matrix(gamma, atoms);
    cout << "gamma matrix: " << endl;
    cout << gamma << endl;

    /*
    vector<STO3G> basis;
    int n = 0; //number of electron pairs

    //constructs basis
    //format for console inputs
    //argv[1]: sample_input/[filename]
    construct_basis(basis, argv[1], n);

    //basis information output
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

    //construct overlap matrix
    mat S(basis.size(), basis.size());
    overlap_matrix(S, basis);

    cout << "Overlap Matrix: " << endl;
    cout << S << endl;

    //construct hamiltonian matrix
    mat H(basis.size(), basis.size());
    hamiltonian_matrix(H, S, basis);

    cout << "Hamiltonian Matrix: " << endl;
    cout << H << endl;

    //construct orthogonalizing matrix
    mat X(basis.size(), basis.size());
    orthogonalization_matrix(X, S);

    cout << "Orthogonalization Matrix: " << endl;
    cout << X << endl;

    //construct coefficient matrix and calculate energy
    mat C(basis.size(), basis.size());
    double E = 0.0;
    coefficient_matrix_energy(H, X, C, E, n);

    cout << "Coefficient Matrix: " << endl;
    cout << C << endl;
    cout << "Total Energy: " << E << endl;
    */

    return 0;
}
