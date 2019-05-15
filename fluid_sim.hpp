//
//  fluid_sim.hpp
//  project
//
//  Created by Tony Cao on 12/4/16.
//  Copyright © 2016 Tony Cao. All rights reserved.
//

#ifndef fluid_sim_hpp
#define fluid_sim_hpp

#include <stdio.h>

void set_boundary(int N, int btype, double *d);

void sim_step(int N, double dt, double viscosity, double diffusion, double sigma, double vf, double vd, double vg,
	double *u, double *u_prev, double *v, double *v_prev, double *p, double *p_prev, double *p_blur,
	double *p_target, double *p_target_blur, double *tmp1, double *tmp2);

#endif /* fluid_sim_hpp */
