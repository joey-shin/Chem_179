#ifndef SPECIAL_PROJECT_ALGO_ONE_H
#define SPECIAL_PROJECT_ALGO_ONE_H

#include <iostream>
#include <vector>
#include <armadillo>

#include "Atom.h"
#include "eval.h"

void coord_matrix(mat& coord, const vector<Atom> & Atoms);
void gradient_descent(vector<Atom> Atoms, double stepsize, const double& epsilon);
void golden_ratio_line_search(vector<Atom> Atoms, double stepsize, const double& epsilon_f, const double& epsilon_e);


#endif //SPECIAL_PROJECT_ALGO_ONE_H
