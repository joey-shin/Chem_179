#include <iostream>
#include <vector>
#include <math.h>
#include <armadillo>

#include "Atom.h"

using namespace std;
using namespace arma;

const double e_Au = 5.29; //binding energy kcal/mol
const double sigma_Au = 2.951; //eq bond dist

//Lennard-Jones potential Atom input
//OBSOLETE
/*
double E_lj(const vector<Atom> & Atoms){ //not needed
    double total_E = 0.0;
    double R_ij;

    for(int i = 0; i < Atoms.size(); i++){
        for(int j = i + 1; j < Atoms.size(); j++){
            R_ij = norm(Atoms[i].coord - Atoms[j].coord);
            total_E += e_Au * (pow((sigma_Au/R_ij), 12.0) - (2 * pow((sigma_Au/R_ij), 6.0)));
        }
    }

    return total_E;
}
 */

//Lennard-Jones potential coord matrix input
double E_lj(const mat& coord){
    double total_E = 0.0;
    double R_ij;

    for(int i = 0; i < coord.n_cols; i++){
        for(int j = i + 1; j < coord.n_cols; j++){
            R_ij = norm(coord.col(i) - coord.col(j));
            total_E += e_Au * (pow((sigma_Au/R_ij), 12.0) - (2 * pow((sigma_Au/R_ij), 6.0)));
        }
    }

    return total_E;
}


//take analytical gradient element, spec. atom and DOF
double F_lj_element(double R_ij, double x_ij){
    double deriv_elem = e_Au * ((12 * pow(sigma_Au, 12.0)/pow(R_ij, 13.0)) - (12 * pow(sigma_Au, 6.0)/pow(R_ij, 7.0)));
    double chain_elem = x_ij / R_ij;
    return deriv_elem * chain_elem;
}


//construct analytical gradient, Atoms input
//OBSOLETE
/*
void F_lj(mat& F, const vector<Atom> & Atoms){ //not needed
    double R_ij;
    double x_ij;

    for(int i = 0; i < Atoms.size(); i++){ //run through all pairs of atom interactions
        for(int j = 0; j < Atoms.size(); j++){
            if(i != j){
                R_ij = norm(Atoms[i].coord - Atoms[j].coord);
                for(int k = 0; k < 3; k++){ //run through each DOF
                    x_ij = Atoms[i].coord(k) - Atoms[j].coord(k);
                    F(k,i) += F_lj_element(R_ij, x_ij);
                }
            }
        }
    }
}

 */

//construct analytical gradient coord matrix input
void F_lj(mat& F, const mat& coord){
    double R_ij;
    double x_ij;
    F.zeros(3, coord.n_cols);

    for(int i = 0; i < coord.n_cols; i++){ //run through all pairs of atom interactions
        for(int j = 0; j < coord.n_cols; j++){
            if(i != j){
                R_ij = norm(coord.col(i) - coord.col(j));
                for(int k = 0; k < 3; k++){ //run through each DOF
                    x_ij = coord.col(i)(k) - coord.col(j)(k);
                    F(k,i) += F_lj_element(R_ij, x_ij);
                }
            }
        }
    }
}