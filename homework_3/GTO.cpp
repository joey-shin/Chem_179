#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <armadillo>
#include <string>
#include <math.h>

#include "GTO.h"
#include "eval.h"


using namespace std;
using namespace arma;



//reads STO3G file and constructs vector of GTO
void read_STO3G(vector<GTO>& primitive, const int& atomic_number, const int& total_L,
                const string& H_STO3G_file, const string& C_STO3G_file){
    ifstream in_H_STO3G;
    ifstream in_C_STO3G;

    in_H_STO3G.open(H_STO3G_file, ios::in); //opening file + catch errors
    if(!in_H_STO3G.is_open()){
        throw runtime_error(H_STO3G_file + " didn't open!");
    }

    in_C_STO3G.open(C_STO3G_file, ios::in); //opening file + catch errors
    if(!in_C_STO3G.is_open()){
        throw runtime_error(C_STO3G_file + " didn't open!");
    }

    string line, temp;

    //sets contraction coefficients of each GTO
    if(atomic_number == 1){
        for(int i = 0; i < 3; i++){
            getline(in_H_STO3G, line);
            istringstream isstring(line);
            isstring >> primitive[i].exp >> primitive[i].d;
        }
    } else if(atomic_number == 6){
        if(total_L == 0){
            for(int i = 0; i < 3; i++){
                getline(in_C_STO3G, line);
                istringstream isstring(line);
                isstring >> primitive[i].exp >> primitive[i].d;
            }
        } else if(total_L == 1){
            for(int i = 0; i < 3; i++){
                getline(in_C_STO3G, line);
                istringstream isstring(line);
                isstring >> primitive[i].exp >> temp >> primitive[i].d;
            }
        }
    }
}


//sets normalization constant N of each GTO
void normalize_GTO(vector<GTO>& primitive, vector<int>& L, vector<double>& coord){
    double overlap;
    for(int i = 0; i < 3; i ++){
        overlap = overlap_3D_ana_int(primitive[i].exp, primitive[i].exp, L, L, coord,
                                     coord);
        primitive[i].N = pow(sqrt(overlap), -1.0);
    }
}