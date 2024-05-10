#include <iostream>
#include <vector>
#include <math.h>
#include <stdexcept>
#include <armadillo>

#include "Atom.h"
#include "eval.h"
#include "algo_one.h"

using namespace std;
using namespace arma;

const double e_Au = 5.29; //binding energy kcal/mol
const double sigma_Au = 2.951; //eq bond dist 2.951


mat vec_to_coord(vec& x){
    int n_cols = x.size() / 3;
    mat coord(3, n_cols);
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < n_cols; j++){
            coord(i, j) = x((3 * j) + i);
        }
    }

    return coord;
}

//Hessian matrix with central differences
void H_central_diff(mat& H, const mat& coord, const double& stepsize){
    vec x = vectorise(coord);
    H.zeros(x.size(), x.size());

    for(int i = 0; i < x.size(); i++){ //run through each DOF in each atom pair interactions
        for(int j = 0; j < x.size(); j++){
            if(i == j){ //central diff diagonal elements
                //construct forward and backwards steps
                vec x_for = x;
                x_for(i) += stepsize;
                vec x_back = x;
                x_back(i) -= stepsize;

                //calculate energies
                double E_for = E_lj(vec_to_coord(x_for));
                double E_back = E_lj(vec_to_coord(x_back));
                double E_x = E_lj(vec_to_coord(x));

                //final formula
                H(i, j) = (E_for + E_back - (2 * E_x)) / (pow(stepsize, 2.0));
            } else{ //central diff off diagonal elements
                //construct forward and backwards steps
                vec xi_for_xj_for = x;
                xi_for_xj_for(i) += stepsize;
                xi_for_xj_for(j) += stepsize;

                vec xi_for_xj_back = x;
                xi_for_xj_back(i) += stepsize;
                xi_for_xj_back(j) -= stepsize;

                vec xi_back_xj_for = x;
                xi_back_xj_for(i) -= stepsize;
                xi_back_xj_for(j) += stepsize;

                vec xi_back_xj_back = x;
                xi_back_xj_back(i) -= stepsize;
                xi_back_xj_back(j) -= stepsize;

                //calculate energies
                double E_xi_for_xj_for = E_lj(vec_to_coord(xi_for_xj_for));
                double E_xi_for_xj_back = E_lj(vec_to_coord(xi_for_xj_back));
                double E_xi_back_xj_for = E_lj(vec_to_coord(xi_back_xj_for));
                double E_xi_back_xj_back = E_lj(vec_to_coord(xi_back_xj_back));

                //final formula
                H(i, j) = (E_xi_for_xj_for - E_xi_for_xj_back - E_xi_back_xj_for + E_xi_back_xj_back) /
                        (4 * pow(stepsize, 2.0));
            }
        }
    }
}

//Hessian matrix with analytical derivative
void H_analytical(mat& H, const mat& coord){
    H.zeros(coord.n_cols * 3, coord.n_cols * 3);
    double R_ij, x_k, x_l;
    double term_1, term_2;

    for(int i = 0; i < coord.n_cols; i++){ //run through each atom pair interaction
        for(int j = 0; j < coord.n_cols; j++){
            if(i != j){ //excluding i != j (avoiding R = 0)
                R_ij = norm(coord.col(i) - coord.col(j));
                for(int k = 0; k < 3; k++){ //run through each chain rule difference components
                    for(int l = 0; l < 3; l++){
                        x_k = coord.col(i)(k) - coord.col(j)(k);
                        x_l = coord.col(i)(l) - coord.col(j)(l);
                        if(k == l){ //diagonal elements
                            term_1 = pow(sigma_Au, 12) * ((168 * pow(x_k, 2) / pow(R_ij, 16)) - (12 / pow(R_ij, 14)));
                            term_2 = 2 * pow(sigma_Au, 6) * ((48 * pow(x_k, 2) / pow(R_ij, 10)) - (6 / pow(R_ij, 8)));
                        } else { //off-diagonal elements
                            term_1 = pow(sigma_Au, 12) * 168 * x_k * x_l / pow(R_ij, 16);
                            term_2 = pow(sigma_Au, 6) * 96 * x_k * x_l / pow(R_ij, 10);
                        }
                        //off diagonal block matrix (i != j)
                        H(3*i + k, 3*j + l) += -1 * e_Au * (term_1 - term_2);
                        //diagonal back matrix (i == j), sum of contributions form off diagonal block matricies
                        H(3*i + k, 3*i + l) += e_Au * (term_1 - term_2);
                    }
                }
            }
        }
    }
}

//newton_method algorithm
void newton_method(vector<Atom> Atoms, const double& epsilon_f, const double& epsilon_e){
    vec delta, g, x;
    mat coord_old, coord_new, F, H, H_i;
    double E_lj_old, E_lj_new;

    coord_matrix(coord_old, Atoms);
    F_lj(F, coord_old);
    g = vectorise(F);
    E_lj_old = E_lj(coord_old);

    int counter = 0;

    while (norm(g) > epsilon_f || abs(E_lj_new - E_lj_old) > epsilon_e){

        //construct vectorized coordinate matrix
        x = vectorise(coord_old);

        //construct vectorized gradient matrix
        g = vectorise(F);

        //construct Hessian matrix
        H_analytical(H, coord_old);

        //construct inverse Hessian matrix
        try{
            H_i = inv(H);
        } catch (const runtime_error& error){
            cerr << "Caught exception: " << error.what() << endl;
            cout << error.what() << ", Algorithm Terminated" << '\n' << endl;
            break;
        }

        //calculate step
        delta = H_i * g;

        //update coordinates
        x += delta;
        coord_new = vec_to_coord(x);
        E_lj_new = E_lj(coord_new);

        coord_old = coord_new;
        E_lj_old = E_lj_new;
        F_lj(F, coord_new);
        g = vectorise(F);

        cout << "Optimized Energy: " << E_lj_new << endl;

        counter++;
    }

    cout << "Iterations: " << counter << endl;
    cout << "Optimized Energy: " << E_lj_new << endl;
    cout << "Optimized Coordinates: " << endl;
    cout << coord_new << endl;
}


/*
//backtrack_line_Search algorithm with Armijo Condition and Curvature Condition
//DOES NOT WORK
double backtrack_line_search(vec& x, vec& g, vec& direction){
    vec x_new, g_new, dot_product_1, dot_product_2, dot_product_3;
    mat F_new;
    double dot_1, dot_2, dot_3;

    double a = 1;
    double c1 = 1e-4;
    double c2 = 0.9;

    double E_lj_old = E_lj(vec_to_coord(x));
    x_new = x + a * direction;
    double E_lj_new = E_lj(vec_to_coord(x_new));
    F_lj(F_new, vec_to_coord(x_new));
    g_new = vectorise(F_new);

    dot_product_1 = g.t() * direction;
    dot_1 = dot_product_1(0);
    dot_product_2 = g_new.t() * direction;
    dot_2 = dot_product_2(0);
    dot_product_3 = g.t() * direction;
    dot_3 = dot_product_3(0);

    cout << "E_lj_new: " << E_lj_new << endl;
    cout << "c1 * a * dot_1: " << c1 * a * dot_1 << endl;
    cout << "dot_2: " << dot_2 << endl;
    cout << "c2 * dot_3: " << c2 * dot_3 << '\n' << endl;

    int counter = 0;

    //Armijo Condition and Curvature Condition
    while((E_lj_new >= (E_lj_old + (c1 * a * dot_1))) && ((dot_2) <= (c2 * dot_3))){
        a *= 0.5;
        x_new = x + a * direction;
        F_lj(F_new, vec_to_coord(x_new));
        g_new = vectorise(F_new);

        dot_product_2 = g_new.t() * direction;
        dot_2 = dot_product_2(0);

        counter++;

        cout << "b_line_search iteration: " << counter << endl;
        cout << "a: " << a << endl;

        cout << "E_lj_new: " << E_lj_new << endl;
        cout << "c1 * a * dot_1: " << c1 * a * dot_1 << endl;
        cout << "dot_2: " << dot_2 << endl;
        cout << "c2 * dot_3: " << c2 * dot_3 << '\n' << endl;
    }

    cout << "a: " << a << endl;
    return a;
}
 */


//BFGS algorithm
//not working
void BFGS(vector<Atom> Atoms, const mat& inital_Hessian, const double& epsilon_f, double stepsize){
    vec delta, x_old, x_new, g_old, g_new, y, direction, left, right, dot_product;
    mat coord, F, H, I, H_i, left_mat, right_mat;
    double a, r;

    coord_matrix(coord, Atoms);
    x_old = vectorise(coord);

    F_lj(F, coord);
    g_old = vectorise(F);

    H = inital_Hessian;
    I = eye(H.n_cols, H.n_cols);

    a = stepsize;
    int counter = 0;

    while(norm(g_old) > epsilon_f){
        direction = -1 * H * g_old;
        //a = backtrack_line_search(x_old, g_old, direction);
        delta = a * direction;
        x_new = x_old + delta;

        F_lj(F, vec_to_coord(x_new));
        g_new = vectorise(F);
        y = g_new - g_old;

        dot_product = y.t() * delta;
        r = pow(dot_product(0), -1);

        left_mat = (I - (r * (delta * y.t())));
        right_mat = (I - (r * (y * delta.t())));

        H_i = left_mat * H * right_mat;
        H = H_i + (r * (delta * delta.t()));

        x_old = x_new;
        g_old = g_new;

        counter++;

        /*
        cout << "BFGS iteration: " << counter << endl;
        cout << "inverse Hessian: " << endl;
        cout << H << endl;
         */
    }

    cout << "Iterations: " << counter << endl;
    cout << "Optimized Energy: " << E_lj(vec_to_coord(x_new)) << endl;
    cout << "Optimized Coordinates: " << endl;
    cout << vec_to_coord(x_new) << endl;
}

