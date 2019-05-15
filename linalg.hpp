//
//  linalg.hpp
//  project
//
//  Created by Tony Cao on 12/6/16.
//  Copyright © 2016 Tony Cao. All rights reserved.
//

#ifndef linalg_hpp
#define linalg_hpp

//#include <stdio.h>
#include <valarray>
#include <functional>

#define matrix std::valarray<std::valarray<double>>
#define vec std::valarray<double>

void conjugate_gradient(double * A, vec d);

matrix toMatrix(double * A, int N);

vec matrixMultiply(matrix & A, vec & d);

double dot(vec & a, vec & b);

vec applyPreconditioner(vec r);

#endif /* linalg_hpp */
