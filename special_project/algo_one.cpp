#include <iostream>
#include <vector>
#include <math.h>
#include <armadillo>

#include "Atom.h"
#include "eval.h"

using namespace std;
using namespace arma;


//extract atom coord into matrix form
//***coord matrix dim doesn't need to be initialized
void coord_matrix(mat& coord, const vector<Atom> & Atoms){
    for(int i = 0; i < Atoms.size(); i++){
        coord = join_rows(coord, Atoms[i].coord);
    }
}


//takes in atoms and updates to optimized structure
void gradient_descent(vector<Atom> Atoms, double stepsize, const double& epsilon){
    mat coord_old, coord_new, F;
    double E_lj_old, E_lj_new;

    coord_matrix(coord_old, Atoms);
    E_lj_old = E_lj(coord_old);
    F_lj(F, coord_old);

    int counter = 0;

    while(norm(F, "fro") > epsilon){

        coord_new = coord_old + stepsize * F/norm(F, "fro");
        E_lj_new = E_lj(coord_new);

        if(E_lj_new < E_lj_old){
            coord_old = coord_new;
            E_lj_old = E_lj_new;
            F_lj(F, coord_new);
            stepsize *= 1.2;
        } else{
            stepsize *= 0.5;
        }

        counter++;
    }

    cout << "Iterations: " << counter << endl;
    cout << "Optimized Energy: " << E_lj_old << endl; // E_lj_old == E_lj_new
    cout << "Optimized Coordinates: " << endl;
    cout << coord_new << endl;
}


//imported form reference code for HW1
void golden_ratio_line_search(vector<Atom> Atoms, double stepsize, const double& epsilon_f, const double& epsilon_e){
    double golden_ratio = 0.38197;
    mat A, B, C, D, opt_coord, F, unit_F;

    double AB, BD, AD;
    double golden_tol = 1e-7;
    double init_stepsize = stepsize;

    coord_matrix(A, Atoms);
    double E_lj_A = E_lj(A);
    F_lj(F, A);

    double old_e = 1e308;
    double cur_e = E_lj_A;

    int counter = 0;
    //main loop
    while(norm(F,"fro") > epsilon_f || abs(cur_e - old_e) > epsilon_e){
        unit_F = F/norm(F,"fro");
        counter++;

        //search for point B
        int bracket_count = 0;
        bool FindB = true;
        B = A + stepsize * unit_F;
        double E_lj_B = E_lj(B);
        while (E_lj_B > E_lj_A) {
            stepsize /= 2; //Halve the step size if B is not found in this iteration
            B = A + stepsize * unit_F;
            E_lj_B = E_lj(B);
            bracket_count++;
            //break the loop if we cannot find B in finite iterations
            if (bracket_count > 100) {
                FindB = false;
                break;
            }
        }

        //search for point D if point B is found
        bracket_count = 0;
        bool FindD = true;
        if (FindB) {
            D = B + stepsize * unit_F;
            double E_lj_D = E_lj(D);
            while (E_lj_D < E_lj_B) {
                stepsize *= 1.2; //Increase the stepsize if D is not found in this iteration
                D = B + stepsize * unit_F;
                E_lj_D = E_lj(D);
                bracket_count++;
                //break the loop if we cannot find D in finite iterations
                if (bracket_count > 100) {
                    FindD = false;
                    break;
                }
            }
        }

        // if we fail to find B or D, it is possible there's no local minima along the
        // negative gradient direction, in this case we move one step towards the negative
        // gradient to decrease energy, then use line search in the future steps
        if (!FindB || !FindD) {
            stepsize = init_stepsize;
            mat new_point = A + stepsize * unit_F;
            while (E_lj(new_point) > E_lj_A){
                stepsize /= 2;
                new_point = A - stepsize * unit_F;
            }
            A = new_point;
            opt_coord = A;
            F_lj(F, A);
            old_e = cur_e;
            cur_e = E_lj(A);
            continue;
        }
        else { //If we find both B and D, start the Golden section search
            AB = norm(B - A, "fro");
            BD = norm(D - B, "fro");
            AD = norm(D - A, "fro");
            if (AB < BD) {
                C = D + golden_ratio * (A - D);
            }
            else {
                C = B;
                B = A + golden_ratio * (D - A);
            }
            while (AD > golden_tol) {
                if (E_lj(B) > E_lj(C)) {
                    A = B;
                    B = C;
                } else {
                    D = C;
                }
                AB = norm(B - A, "fro");
                BD = norm(D - B, "fro");
                AD = norm(D - A, "fro");
                if (AB < BD) {
                    C = D + golden_ratio*(A - D);
                } else {
                    C = B;
                    B = A + golden_ratio*(D - A);
                }
            }
            //renew the point according to the relative energy between B and C
            if (E_lj(B) > E_lj(C)) {
                opt_coord = C;
            } else {
                opt_coord = B;
            }
            F_lj(F, opt_coord);
            old_e = cur_e;
            cur_e = E_lj(opt_coord);
        }
    }

    cout << "Iterations: " << counter << endl;
    cout << "Optimized Energy: " << E_lj(opt_coord) << endl;
    cout << "Optimized Coordinates: " << endl;
    cout << opt_coord << endl;
}


//timing function runtime
//auto start = std::chromo::high_resolution_clock::now();