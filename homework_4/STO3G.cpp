#include <iostream>
#include <vector>
#include <armadillo>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <math.h>

using namespace std;
using namespace arma;

#include "STO3G.h"
#include "GTO.h"


/*
//reads file and constructs basis
void construct_basis(vector<STO3G>& basis, const string& mol_file, int& n){
    ifstream in_mol;

    in_mol.open(mol_file, ios::in); //opening file + catch errors
    if(!in_mol.is_open()){
        throw runtime_error(mol_file + " didn't open!");
    }

    string line;
    int number_of_atoms, charge, atomic_number, total_L;
    int a = 0, b = 0;

    //read number of atoms and charge
    getline(in_mol, line);
    istringstream isstring(line);
    isstring >> number_of_atoms >> charge;

    //set up primitive GTO vector and read atomic number and coordinates
    for(int i = 0; i < number_of_atoms; i++){
        getline(in_mol, line);
        istringstream isstring(line);
        mat L_matrix;

        //get atomic number
        isstring >> atomic_number;
        if(atomic_number <= 0){
            throw runtime_error("atomic number invalid");
        }

        //set up L
        if(atomic_number == 1){
            cartesian_angular_momentum(L_matrix, 0);
        } else if(atomic_number >= 3 && atomic_number <= 9){
            cartesian_angular_momentum(L_matrix, 0);
            cartesian_angular_momentum(L_matrix, 1);
        }

        //calculating number of MO
        if(atomic_number == 1){
            b += 1;
        } else if(atomic_number == 6){
            a += 1;
        }


        //create instance coord vector
        double coord_temp;
        vector<double> current_coord;
        for(int j = 0; j < 3; j++){
            isstring >> coord_temp;
            current_coord.push_back(coord_temp);
        }

        //create new STO3G for each set of quantum numbers
        for(int j = 0; j < L_matrix.n_rows; j++){
            //create new instance STO3G basis
            STO3G instance_basis;

            //set up atomic number, coord vector, L vector, total_L in current STO3G basis
            instance_basis.atomic_number = atomic_number;

            for(int k = 0; k < 3; k++){
                instance_basis.coord.push_back(current_coord[k]);
            }

            total_L = 0;
            for(int k = 0; k < 3; k++){
                total_L += L_matrix(j, k);
                instance_basis.L.push_back(L_matrix(j, k));
            }
            instance_basis.total_L = total_L;

            //set up primitive GTO vector in current STO3G basis
            for(int k = 0; k < 3; k++){
                GTO gto;
                instance_basis.primitive.push_back(gto);
            }

            //set up GTOs in GTO vector
            read_STO3G(instance_basis.primitive, instance_basis.atomic_number, instance_basis.total_L);
            normalize_GTO(instance_basis.primitive, instance_basis.L, instance_basis.coord);

            basis.push_back(instance_basis);
        }
    }

    //calculate number of electron pairs
    n = ((4 * a) + b) / 2;

    //catch if n is not an integer
    if(floor(n) != n){
        throw runtime_error("number of electrons not a whole number");
    }

    in_mol.close();
}
 */


//reads file and constructs a single basis
void construct_AO(STO3G& AO, int& atomic_number, int& total_L, const vec& coord){
    AO.total_L = total_L;
    AO.coord = coord;
    AO.atomic_number = atomic_number;

    for(int i = 0; i < 3; i++){
        GTO gto;
        AO.primitive.push_back(gto);
    }

    //assign basis file name
    string basis_file_name;
    switch (atomic_number) {
        case 1:
            basis_file_name = "basis/H_STO3G.txt";
            AO.Z = 1;
            AO.IA = 7.176;
            AO.beta = -9;
            break;
        case 6:
            basis_file_name = "basis/C_STO3G.txt";
            AO.Z = 4;
            if(total_L == 0){
                AO.IA = 14.051;
            } else if(total_L == 1){
                AO.IA = 5.572;
            }
            AO.beta = -21;
            break;
        case 7:
            basis_file_name = "basis/N_STO3G.txt";
            AO.Z = 5;
            if(total_L == 0){
                AO.IA = 19.361;
            } else if(total_L == 1){
                AO.IA = 7.275;
            }
            AO.beta = -25;
            break;
        case 8:
            basis_file_name = "basis/O_STO3G.txt";
            AO.Z = 6;
            if(total_L == 0){
                AO.IA = 25.390;
            } else if(total_L == 1){
                AO.IA = 9.111;
            }
            AO.beta = -31;
            break;
        case 9:
            basis_file_name = "basis/F_STO3G.txt";
            AO.Z = 7;
            if(total_L == 0){
                AO.IA = 32.272;
            } else if(total_L == 1){
                AO.IA = 11.080;
            }
            AO.beta = -39;
            break;
        default:
            throw runtime_error("Atomic number " + to_string(atomic_number) + " isn't valid nor defined");
    }

    read_STO3G(AO.primitive, AO.total_L, basis_file_name);
    normalize_GTO(AO.primitive, AO.L_vec, coord);
}