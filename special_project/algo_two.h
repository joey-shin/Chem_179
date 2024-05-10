#ifndef SPECIAL_PROJECT_ALGO_TWO_H
#define SPECIAL_PROJECT_ALGO_TWO_H

#include <iostream>
#include <vector>
#include <armadillo>

#include "Atom.h"
#include "eval.h"

void H_central_diff(mat& H, const mat& coord, const double& stepsize);
void H_analytical(mat& H, const mat& coord);
void newton_method(vector<Atom> Atoms, const double& epsilon_f, const double& epsilon_e);
void BFGS(vector<Atom> Atoms, const mat& inital_Hessian, const double& epsilon_f, double stepsize);


#endif //SPECIAL_PROJECT_ALGO_TWO_H
