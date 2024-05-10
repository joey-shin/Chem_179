#ifndef SPECIAL_PROJECT_EVAL_H
#define SPECIAL_PROJECT_EVAL_H

#include <iostream>
#include <vector>
#include <armadillo>

#include "Atom.h"

//double E_lj(const vector<Atom> & Atoms);
double E_lj(const mat& coord);
//void F_lj(mat& F, const vector<Atom> & Atoms);
void F_lj(mat& F, const mat& coord);

#endif //SPECIAL_PROJECT_EVAL_H
