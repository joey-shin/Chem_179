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
void read_STO3G(vector<GTO>& primitive, const int& atomic_number, const int& total_L, string& basis_file_name){
    string line, temp;
    ifstream in;

    in.open(basis_file_name, ios::in); //opening file + catch errors
    if(!in.is_open()){
        throw runtime_error(basis_file_name + " didn't open!");
    }

    if(total_L == 0){
        for(int i = 0; i < 3; i++){
            getline(in, line);
            istringstream isstring(line);
            isstring >> primitive[i].exp >> primitive[i].d;
        }
    } else if(total_L == 1){
        for(int i = 0; i < 3; i++){
            getline(in, line);
            istringstream isstring(line);
            isstring >> primitive[i].exp >> temp >> primitive[i].d;
        }
    }
}


//sets normalization constant N of each GTO
void normalize_GTO(vector<GTO>& primitive, vector<int>& L, const vec& coord){
    double overlap;
    for(int i = 0; i < 3; i ++){
        overlap = overlap_3D_ana_int(primitive[i].exp, primitive[i].exp, L, L, coord,
                                     coord);
        primitive[i].N = pow(sqrt(overlap), -1.0);
    }
}