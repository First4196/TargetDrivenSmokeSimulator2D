//
//  fluid_sim.cpp
//  project
//
//  Created by Tony Cao on 12/4/16.
//  Copyright © 2016 Tony Cao. All rights reserved.
//

#include "fluid_sim.hpp"
#include "target_driven.hpp"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <algorithm>

#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x,y) {double *_=x; x=y; y=_;}

void add_source ( int N, double *d, double *s, double dt ){
    int size=(N+2)*(N+2);
	for (int i = 0; i < size; i++) {
		d[i] += dt * s[i];
	}
}

void set_boundary ( int N, int btype, double *d ){
	// p : btype 0
	// u : btype 1
	// v : btype 2
    for (int i=1 ; i<=N ; i++ ) {
        d[IX(0  ,i)] = btype == 1 ? -d[IX(1,i)] : d[IX(1,i)];
        d[IX(N+1,i)] = btype == 1 ? -d[IX(N,i)] : d[IX(N,i)];
        d[IX(i,0  )] = btype == 2 ? -d[IX(i,1)] : d[IX(i,1)];
        d[IX(i,N+1)] = btype == 2 ? -d[IX(i,N)] : d[IX(i,N)];
    }
    d[IX(0  ,0  )] = 0.5f * (d[IX(1,0  )]+d[IX(0  ,1)]);
    d[IX(0  ,N+1)] = 0.5f * (d[IX(1,N+1)]+d[IX(0  ,N)]);
    d[IX(N+1,0  )] = 0.5f * (d[IX(N,0  )]+d[IX(N+1,1)]);
    d[IX(N+1,N+1)] = 0.5f * (d[IX(N,N+1)]+d[IX(N+1,N)]);
}

// iterative solver 
// laplace(d) = d0
void lin_solve(int N, int btype, double *d, double *d0, double a, double c){
    for (int k=0; k<20; k++) {
		for (int i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) {
				d[IX(i, j)] = (d0[IX(i, j)] + a * (d[IX(i - 1, j)] + d[IX(i + 1, j)] + d[IX(i, j - 1)] + d[IX(i, j + 1)])) / c;
			}
		}
		set_boundary(N, btype, d);
    }
}

void diffuse ( int N, int btype, double *d, double *d0, double diff, double dt ){
    double a = dt * diff * N * N;
    lin_solve ( N, btype, d, d0, a, 1+4*a);
}

void advect ( int N, int btype, double *d, double *d_prev, double *u, double *v, double dt ){
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			double x = i - dt * N * u[IX(i, j)];
			double y = j - dt * N * v[IX(i, j)];
			if (x < 0.5f) {
				x = 0.5f;
			}
			if (x > N + 0.5f) {
				x = N + 0.5f;
			}
			int i0 = (int)x; 
			int i1 = i0 + 1;
			if (y < 0.5f) {
				y = 0.5f;
			}
			if (y > N + 0.5f) {
				y = N + 0.5f;
			}
			int j0 = (int)y; 
			int j1 = j0 + 1;
			double s1 = x - i0; 
			double s0 = 1.0f - s1;
			double t1 = y - j0;
			double t0 = 1.0f - t1;
			d[IX(i, j)] = s0 * (t0 * d_prev[IX(i0, j0)] + t1 * d_prev[IX(i0, j1)]) +
				s1 * (t0 * d_prev[IX(i1, j0)] + t1 * d_prev[IX(i1, j1)]);
		}
	}
	set_boundary(N, btype, d);
}

void project (int N, double *u, double *v, double *p, double *divergence){
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			divergence[IX(i, j)] = -0.5f*(u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
			p[IX(i, j)] = 0;
		}
	}
	set_boundary(N, 0, divergence);
	set_boundary(N, 0, p);
    
    lin_solve(N, 0, p, divergence, 1, 4);
    
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
			v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
		}
	}
	set_boundary(N, 1, u);
	set_boundary(N, 2, v);
}

void sim_step(int N, double dt, double viscosity, double diffusion, double sigma, double vf, double vd, double vg,
	double *u, double *u_prev, double *v, double *v_prev, double *p, double *p_prev, double *p_blur, 
	double *p_target, double *p_target_blur, double *tmp1, double *tmp2){
	
	// calculate p_blur
	guassian_blur(N, p_blur, p, tmp1, sigma);

	// vel step
	SWAP(u_prev, u);
	SWAP(v_prev, v);
	diffuse(N, 1, u, u_prev, viscosity, dt);
	diffuse(N, 2, v, v_prev, viscosity, dt);
	project(N, u, v, tmp1, tmp2);
	SWAP(u_prev, u); 
	SWAP(v_prev, v);
	add_driving_force(N, u, v, u_prev, v_prev, p_blur, p_target_blur, vf, dt);
	SWAP(u_prev, u);
	SWAP(v_prev, v);
	attenuate(N, u, v, u_prev, v_prev, vd, dt);
	SWAP(u_prev, u);
	SWAP(v_prev, v);
	advect(N, 1, u, u_prev, u_prev, v_prev, dt); 
	advect(N, 2, v, v_prev, u_prev, v_prev, dt);
	project(N, u, v, tmp1, tmp2);

	// p step
	SWAP(p_prev, p);
	diffuse(N, 0, p, p_prev, diffusion, dt);
	SWAP(p_prev, p); 
	advect(N, 0, p, p_prev, u, v, dt);
	SWAP(p_prev, p);
	gather(N, p, p_prev, p_target, p_target_blur, tmp1, tmp2, vg, dt);
}