#include <iostream>
#include <armadillo>

#include "GTO.h"
#include "STO3G.h"
#include "Atom.h"
#include "eval.h"

int main(int argc, char* argv[]) {

    //construct vector of atoms
    vector<Atom> molecule;
    Basis basis;
    construct_Atoms(molecule, argv[1]);
    construct_basis(basis, molecule);

    /*
     *
     *
    //print all values of Atom structure and sub structures
    //for debugging purpose
    cout << "numb er of Atoms: " << molecule.size() << endl;
    for(int i = 0; i < molecule.size(); i++){
        cout << '\n' << "Atom " << i << '\n';

        for(int j = 0; j < molecule[i].AOs.size(); j++){
            cout << '\n' << "STO3G " << j << '\n';
            cout << "total_L: " << molecule[i].AOs[j].total_L << '\n';
            cout << "atomic number: " << molecule[i].AOs[j].atomic_number << '\n';
            cout << "IA: " << molecule[i].AOs[j].IA << '\n';
            cout << "beta: " << molecule[i].AOs[j].beta << '\n';
            cout << "Z: " << molecule[i].AOs[j].Z << '\n';
            cout << "L: ";
            for(int k = 0; k < 3; k++){
                cout << molecule[i].AOs[j].L_vec[k] << ' ';
            }
            cout << endl;
            cout << "Coord: ";
            for(int k = 0; k < 3; k++){
                cout << molecule[i].AOs[j].coord[k] << ' ';
            }
            cout << '\n' << endl;

            for(int k = 0; k < molecule[i].AOs[j].primitive.size(); k++){
                cout << "primitive " << k << '\n';
                cout << "exp " << molecule[i].AOs[j].primitive[k].exp << '\n';
                cout << "d " << molecule[i].AOs[j].primitive[k].d << '\n';
                cout << "N " << molecule[i].AOs[j].primitive[k].N << '\n' << endl;
            }
        }
    }

    //print all values of Basis structure
    cout << '\n' << "number of basis: " << basis.basis.size() << endl;
    for(int i = 0; i < basis.basis.size(); i++){
        cout << '\n' << "basis " << i << '\n';
        cout << "L: ";
        for(int k = 0; k < 3; k++){
            cout << basis.basis[i].L_vec[k] << ' ';
        }
        cout << endl;
        cout << "L: ";
        for(int k = 0; k < 3; k++){
            cout << basis.basis[i].coord[k] << ' ';
        }
        cout << '\n' << endl;
    }
    cout << "atom index: " << endl;
    for(int i = 0; i < basis.atom_index.size(); i++){
        cout << basis.atom_index[i] << " ";
    }
    cout << '\n' << endl;
    cout  << "n: " << basis.n << endl;
    cout  << "p: " << basis.p << endl;
    cout  << "q: " << basis.q << '\n' << endl;
     */

    mat gamma(molecule.size(), molecule.size());
    mat S(basis.basis.size(), basis.basis.size());
    mat H(basis.basis.size(), basis.basis.size());

    mat P_A_i(basis.basis.size(), basis.basis.size(), fill::zeros);
    mat P_B_i(basis.basis.size(), basis.basis.size(), fill::zeros);

    mat F_A(basis.basis.size(), basis.basis.size());
    mat F_B(basis.basis.size(), basis.basis.size());
    mat P_A(basis.basis.size(), basis.basis.size());
    mat P_B(basis.basis.size(), basis.basis.size());

    //gamma, S, H matrix construction
    gamma_matrix(gamma, molecule);
    overlap_matrix(S, basis);
    core_hamiltonian_matrix(H, S, basis, molecule);

    //SCF convergence algorithm
    SCF(F_A, F_B, P_A, P_B, P_A_i, P_B_i, S, basis, molecule);

    //energy calculation
    double electron_Energy = electron_E(H, F_A, F_B, P_A, P_B, basis);
    double nuclear_repulsion_E = nuc_repulsion_E(molecule);
    double total_energy = electron_Energy + nuclear_repulsion_E;



    cout << "gamma matrix: " << endl;
    cout << gamma << endl;

    cout << "S: " << endl;
    cout << S << endl;

    cout << "core H: " << endl;
    cout << H << endl;

    //print energy values
    cout << "electron energy: " << electron_Energy << '\n';
    cout << "nuclear repulsion energy: " << nuclear_repulsion_E << '\n';
    cout << "total energy: " << total_energy << '\n' << endl;

    mat V_nuc_R;
    V_nuc_R_matrix(V_nuc_R, molecule);

    cout << "V_nuc_R: " << endl;
    cout << V_nuc_R << endl;


    mat X(basis.basis.size(), basis.basis.size());
    X_matrix(X, basis, P_A, P_B);

    cout << "X: " << endl;
    cout << X << endl;


    mat gamma_R;
    gamma_R_matrix(gamma_R, molecule);

    cout << "gamma_R: " << endl;
    cout << gamma_R << endl;

    return 0;
}
