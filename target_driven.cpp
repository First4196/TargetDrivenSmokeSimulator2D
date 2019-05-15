//
//  target_driven.cpp
//
//  Created by Thanapat Katapornsiri on 1 May 2019.
//  Copyright © 2019 Thanapat Katapornsiri. All rights reserved.
//

#include "target_driven.hpp"
#include "fluid_sim.hpp"
#include <iostream> 
#include <stdio.h>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>

#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x,y) {double *_=x; x=y; y=_;}

void guassian_blur(int N, double *p_blur, double *p, double *p_blur_h, double sigma){
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			p_blur_h[IX(i, j)] = 0;
			p_blur[IX(i, j)] = 0;
		}
	}
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			for (int k = std::max(1, i - 10); k <= std::min(N, i + 10); k++) {
				p_blur_h[IX(k, j)] += p[IX(i, j)] * exp(-(k - i)*(k - i)/(N * N * sigma * sigma));
			}
		}
	}
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			for (int k = std::max(1, j - 10); k <= std::min(N, j + 10); k++) {
				p_blur[IX(i, k)] += p_blur_h[IX(i, j)] * exp(-(k - j)*(k - j) / (N * N * sigma * sigma));
			}
		}
	}
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			p_blur[IX(i, j)] /= 2 * M_PI * sigma * sigma;
		}
	}
	set_boundary(N, 0, p_blur);
}

void add_driving_force(int N, double *u, double *v, double *u_prev, double *v_prev, double *p_blur, double *p_target_blur, double vf, double dt) {
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			double norm_grad_x = 0.5f * N * (p_target_blur[IX(i + 1, j)] - p_target_blur[IX(i - 1, j)]) / p_target_blur[IX(i, j)];
			double norm_grad_y = 0.5f * N * (p_target_blur[IX(i, j + 1)] - p_target_blur[IX(i, j - 1)]) / p_target_blur[IX(i, j)];
			u[IX(i, j)] = u_prev[IX(i, j)] + dt * vf * norm_grad_x * p_blur[IX(i, j)];
			v[IX(i, j)] = v_prev[IX(i, j)] + dt * vf * norm_grad_y * p_blur[IX(i, j)];
		}
	}
	set_boundary(N, 1, u);
	set_boundary(N, 2, v);
}

void attenuate(int N, double *u, double *v, double *u_prev, double *v_prev, double vd, double dt) {
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			u[IX(i, j)] = u_prev[IX(i, j)] - dt * vd * u_prev[IX(i, j)];
			v[IX(i, j)] = v_prev[IX(i, j)] - dt * vd * v_prev[IX(i, j)];
		}
	}
	set_boundary(N, 1, u);
	set_boundary(N, 2, v);
}

void gather(int N, double *p, double *p_prev, double *p_target, double *p_target_blur, 
	double *error_flux_x, double *error_flux_y, double vg, double dt) {
	
	// error_flux
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			error_flux_x[IX(i, j)] = 0.5f * N * 
				(p_prev[IX(i+1, j)] - p_target[IX(i+1, j)]) - (p_prev[IX(i-1, j)] - p_target[IX(i-1, j)]);
			error_flux_y[IX(i, j)] = 0.5f * N *
				(p_prev[IX(i, j+1)] - p_target[IX(i, j+1)]) - (p_prev[IX(i, j-1)] - p_target[IX(i, j-1)]);
		}
	}
	set_boundary(N, 0, error_flux_x);
	set_boundary(N, 0, error_flux_y);

	// gather smoke
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			double x0 = p_prev[IX(i-1, j)] * p_target_blur[IX(i-1, j)] * error_flux_x[IX(i-1, j)];
			double x1 = p_prev[IX(i+1, j)] * p_target_blur[IX(i+1, j)] * error_flux_x[IX(i+1, j)];
			double y0 = p_prev[IX(i, j-1)] * p_target_blur[IX(i, j-1)] * error_flux_y[IX(i, j-1)];
			double y1 = p_prev[IX(i, j+1)] * p_target_blur[IX(i, j+1)] * error_flux_y[IX(i, j+1)];
			double divergence = 0.5 * N * (x1 - x0 + y1 - y0);
			p[IX(i, j)] = p_prev[IX(i, j)] + dt * vg * divergence;
			if (p[IX(i, j)] < 0) {
				p[IX(i, j)] = 0;
			}
			if (p[IX(i, j)] > 1) {
				p[IX(i, j)] = 1;
			}
		}
	}
	set_boundary(N, 0, p);
}

double get_total(int N, double *d) {
	double total = 0;
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			total += d[IX(i, j)];
		}
	}
	return total;
}