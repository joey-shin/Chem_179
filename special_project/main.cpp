#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <armadillo>
#include <chrono>
#include <thread>

#include "Atom.h"
#include "eval.h"
#include "algo_one.h"
#include "algo_two.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {

    string file = argv[1];
    vector<Atom> Atoms;

    construct_Atoms(Atoms, file);
    print_Atoms(Atoms);

    mat coord;
    coord_matrix(coord, Atoms);

    double stepsize = 0.3;
    double epsilon = 1e-2;
    double epsilon_f = 1e-4;
    double epsilon_e = 1e-2;

    mat inital_Hessian = eye(3 * coord.n_cols, 3 * coord.n_cols);


    //output for gradient descent
    cout << "gradient descent" << endl;
    auto start = chrono::steady_clock::now();
    gradient_descent(Atoms, stepsize, epsilon);
    auto end = chrono::steady_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Execution time: " << duration.count() << " ms" << '\n' << endl;


    //output for golden_ratio_line_search
    cout << "golden_ratio_line_search" << endl;
    start = chrono::steady_clock::now();
    golden_ratio_line_search(Atoms, stepsize, epsilon_f, epsilon_e);
    end = chrono::steady_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Execution time: " << duration.count() << " ms" << '\n' << endl;


    //output for Newton Step
    cout << "Newton_step" << endl;
    start = chrono::steady_clock::now();
    newton_method(Atoms, epsilon_f, epsilon_e);
    end = chrono::steady_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Execution time: " << duration.count() << " ms" << '\n' << endl;


    //output for golden_ratio_line_search + newton step



    //output for BFGS
    cout << "BFGS" << endl;
    start = chrono::steady_clock::now();
    BFGS(Atoms, inital_Hessian, epsilon_f, stepsize);
    end = chrono::steady_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Execution time: " << duration.count() << " ms" << '\n' << endl;


}
