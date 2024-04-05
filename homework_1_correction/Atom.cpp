#include <iostream>
#include <vector>
#include <armadillo>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "Atom.h"

void construct_Atoms(const vector<Atom>& Atoms, const string& file){
    ifstream in;

    in.open(file, ios::in);
    if(!in.is_open()){
        throw runtime_error(file + " didn't open!");
    }

    string line;
    int number_of_atoms;

    //read number of atoms and charge
    getline(in, line);
    istringstream isstring(line);
    isstring >> number_of_atoms;

    for(int i = 0; i < number_of_atoms; i++){
        Atom ref_atom;

        getline(in, line);
        istringstream isstring(line);
        isstring >> atom.atomic_number >> atom.coord(0) >> atom.coord(1) >> atom.coord(2);

        Atoms.push_back(ref_atom);
    }
}

