//
//  target_driven.hpp
//
//  Created by Thanapat Katapornsiri on 1 May 2019.
//  Copyright © 2019 Thanapat Katapornsiri. All rights reserved.
//

#ifndef target_driven_hpp
#define target_driven_hpp

#include <stdio.h>

void guassian_blur(int N, double *p_blur, double *p, double *p_blur_h, double sigma);
void add_driving_force(int N, double *u, double *v, double *u_prev, double *v_prev, double *p_blur, double *p_target_blur, double vf, double dt);
void attenuate(int N, double *u, double *v, double *u_prev, double *v_prev, double vd, double dt);
void gather(int N, double *p, double *p_prev, double *p_target, double *p_target_blur, 
	double *error_flux_x, double *error_flux_y, double vg, double dt);
double get_total(int N, double *d);

#endif /* target_driven_hpp */
