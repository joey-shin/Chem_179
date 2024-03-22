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
    cout << "number of Atoms: " << molecule.size() << endl;
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

    //gamma matrix construction
    mat gamma(molecule.size(), molecule.size());
    gamma_matrix(gamma, molecule);
    cout << "gamma matrix: " << endl;
    cout << gamma << endl;

    //overlap matrix construction
    mat S(basis.basis.size(), basis.basis.size());
    overlap_matrix(S, basis);
    cout << "S: " << endl;
    cout << S << endl;

    //hailtonian matrix construction
    mat H(basis.basis.size(), basis.basis.size());
    core_hamiltonian_matrix(H, S, basis, molecule);
    cout << "core H: " << endl;
    cout << H << endl;

    //initial density matrix construction (P_alpha_initial = P_beta_initial = 0)
    mat P_A_i(basis.basis.size(), basis.basis.size(), fill::zeros);
    mat P_B_i(basis.basis.size(), basis.basis.size(), fill::zeros);

    mat F_A(basis.basis.size(), basis.basis.size());
    mat F_B(basis.basis.size(), basis.basis.size());
    mat P_A(basis.basis.size(), basis.basis.size());
    mat P_B(basis.basis.size(), basis.basis.size());

    //SCF convergence algorithm
    SCF(F_A, F_B, P_A, P_B, P_A_i, P_B_i, S, basis, molecule);

    //energy calculation
    double electron_Energy = electron_E(H, F_A, F_B, P_A, P_B, basis);
    cout << "electron energy: " << electron_Energy << '\n';

    double nuclear_repulsion_E = nuc_repulsion_E(molecule);
    cout << "nuclear repulsion energy: " << nuclear_repulsion_E << '\n';

    double total_energy = electron_Energy + nuclear_repulsion_E;
    cout << "total energy: " << total_energy << '\n';

    return 0;
}
