#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <armadillo>

#include "Atom.h"

using namespace std;
using namespace arma;

//built vector of Atom (assuming that all atoms are Gold)
void construct_Atoms(vector<Atom>& Atoms, const string& file){
    ifstream in;

    in.open(file, ios::in); //opening file + catch errors
    if(!in.is_open()){
        throw runtime_error(file + " didn't open!");
    }

    string line;
    int number_of_atoms, atomic_number;
    vec coord(3);

    getline(in, line);
    number_of_atoms = stoi(line);

    for(int i = 0; i < number_of_atoms; i++){
        Atom instance_Atom;

        getline(in, line);
        istringstream isstring(line);

        //get atomic number
        isstring >> atomic_number >> coord(0) >> coord(1) >> coord(2);
        instance_Atom.atomic_number = atomic_number;
        instance_Atom.coord = coord;

        Atoms.push_back(instance_Atom);
    }

    in.close();
}

//print data structures
void print_Atoms(const vector<Atom>& Atoms){
    cout << "number of atoms: " << Atoms.size() << endl;
    for(int i = 0; i < Atoms.size(); i++){
        cout << "Atom " << i << ": " << endl;
        cout << "atomic_number: " << Atoms[i].atomic_number << endl;
        cout << "coord: " << '\n' << Atoms[i].coord << endl;
    }
}