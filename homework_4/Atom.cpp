#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <armadillo>
#include <string>
#include <math.h>


#include "STO3G.h"
#include "GTO.h"
#include "Atom.h"

//set angular momentum matrix will all possible orientations given total quantum number
void cartesian_angular_momentum(mat &L, const int &total_L){
    int total_sum;

    for(int i = 0; i <= total_L; i++){
        for(int j = 0; j <= total_L; j++){
            for(int k = 0; k <= total_L; k++){
                total_sum = i + j + k;
                if (total_sum == total_L){
                    L.resize(L.n_rows + 1, 3);
                    L(L.n_rows - 1, 0) = k;
                    L(L.n_rows - 1, 1) = j;
                    L(L.n_rows - 1, 2) = i;
                }
            }
        }
    }
}


//reads file and constructs basis
void construct_Atoms(vector<Atom>& atoms, const string& mol_file){
    ifstream in_mol;

    in_mol.open(mol_file, ios::in); //opening file + catch errors
    if(!in_mol.is_open()){
        throw runtime_error(mol_file + " didn't open!");
    }

    string line;
    int number_of_atoms, charge, atomic_number;

    //read number of atoms and charge
    getline(in_mol, line);
    istringstream isstring(line);
    isstring >> number_of_atoms >> charge;

    //set up primitive GTO vector and read atomic number and coordinates
    for(int i = 0; i < number_of_atoms; i++){
        Atom instance_Atom;

        getline(in_mol, line);
        istringstream isstring(line);

        //get atomic number
        isstring >> atomic_number;
        if(atomic_number <= 0){
            throw runtime_error("atomic number invalid");
        }

        instance_Atom.atomic_number = atomic_number;

        mat L_matrix;

        //set up L matrix
        if(atomic_number == 1){
            cartesian_angular_momentum(L_matrix, 0);
        } else if(atomic_number >= 3 && atomic_number <= 9){
            cartesian_angular_momentum(L_matrix, 0);
            cartesian_angular_momentum(L_matrix, 1);
        }

        //create instance coord vector
        double coord_temp;

        //vector<double> current_coord;
        instance_Atom.coord.set_size(3);

        for(int j = 0; j < 3; j++){
            isstring >> coord_temp;
            instance_Atom.coord(j) = coord_temp;
        }

        //create new STO3G for each set of quantum numbers
        for(int j = 0; j < L_matrix.n_rows; j++){
            //create new instance STO3G basis
            STO3G instance_AO;

            int total_L = 0;
            vector<int> L_Vec;

            //create L vectors in instance_AO
            for(int k = 0; k < 3; k++){
                total_L += L_matrix(j, k);
                instance_AO.L_vec.push_back(L_matrix(j, k));
                instance_AO.coord(k) = instance_Atom.coord(k);
            }
            instance_AO.total_L = total_L;

            //set up atomic number, coord vector, L vector, total_L in current STO3G basis
            construct_AO(instance_AO, atomic_number, total_L, instance_Atom.coord);
            instance_Atom.AOs.push_back(instance_AO);
        }
        atoms.push_back(instance_Atom);
    }

    in_mol.close();
}

void construct_basis(Basis& basis, const vector<Atom>& atoms){
    basis.n = 0;

    for(int i = 0; i < atoms.size(); i++){
        for(int j = 0; j < atoms[i].AOs.size(); j++){
            basis.basis.push_back(atoms[i].AOs[j]);
        }
        basis.n += atoms[i].AOs[0].Z; //Z in all AO contained in atom identical, using index 0
    }

    basis.p = ceil(basis.n / 2.0);
    basis.q = floor(basis.n / 2.0);
}