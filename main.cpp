//
//  main.cpp
//  project
//
//  Created by Tony Cao on 12/4/16.
//  Copyright © 2016 Tony Cao. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <windows.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "linalg.hpp"
#include "vecmath.h"
#include "GL/GL.H"
#include "GL/GLAux.h"
#include "GL/GLU.H"
#include "GL/glui.h"
#include "GL/glut.h"

#include "fluid_sim.hpp"
#include "target_driven.hpp"

#define IX(i,j) ((i)+(N+2)*(j))
#define NUM_DISPLAY_MODES 2
#define DISPLAY_MODE_SMOKE 0
#define DISPLAY_MODE_VELOCITY 1

// params
static int N; // grid width and height
static double dt; // timestep
static double viscosity;
static double diffusion;

// target driven params
static double sigma;
static double vf;
static double vd;
static double vg;
static int num_targets;
static std::vector<double> targets_time;

// target driven varaibles
static int current_target;
static double current_time;

// grid storage
static double *u, *u_prev; // velocity in x-direction
static double *v, *v_prev; // velocity in y-direction
static double *p, *p_prev, *p_blur; // density
static double *tmp1, *tmp2; // temps
static std::vector<double*> p_targets, p_targets_blur; // target densities

// display
static int win_id;
static int win_x, win_y;
static int display_mode;

// mouse input
static int mouse_down[3];
static int mouse_x, prev_mouse_x;
static int mouse_y, prev_mouse_y;

// const
static double MAX_COLOR = 255.0f;

// free data
static void free_data ( void ){
	if (u) free(u);
	if (u_prev) free(u_prev);
	if (v) free(v);
	if (v_prev) free(v_prev);
	if (p) free(p);
	if (p_prev) free(p_prev);
	if (p_blur) free(p_blur);
	if (tmp1) free(tmp1);
	if (tmp2) free(tmp2);
	for (int k = 0; k < p_targets.size(); k++) {
		if (p_targets[k]) free(p_targets[k]);
	}
	p_targets.clear();
	for (int k = 0; k < p_targets_blur.size(); k++) {
		if (p_targets_blur[k]) free(p_targets_blur[k]);
	}
	p_targets_blur.clear();
}

// clear data
static void clear_data ( void ){
	int size = (N + 2)*(N + 2);
    for(int i=0; i<size; i++) {
        u[i] = u_prev[i] = v[i] = v_prev[i] = p[i] = p_prev[i] = p_blur[i] = tmp1[i] = tmp2[i] = 0.0f;
    }
	for (int k = 0; k < p_targets.size(); k++) {
		for (int i = 0; i < size; i++) {
			p_targets[k][i] = 0.0f;
		}
	}
	for (int k = 0; k < p_targets_blur.size(); k++) {
		for (int i = 0; i < size; i++) {
			p_targets_blur[k][i] = 0.0f;
		}
	}
}

// allocate data
static int allocate_data ( void ){
	int size = (N + 2)*(N + 2);
  
	u = (double *)malloc(size * sizeof(double));
	u_prev = (double *)malloc(size * sizeof(double));
	v = (double *)malloc(size * sizeof(double));
	v_prev = (double *)malloc(size * sizeof(double));
	p = (double *)malloc(size * sizeof(double));
	p_prev = (double *)malloc(size * sizeof(double));
	p_blur = (double *)malloc(size * sizeof(double));
	tmp1 = (double *)malloc(size * sizeof(double));
	tmp2 = (double *)malloc(size * sizeof(double));
	p_targets.resize(num_targets + 1);
	for (int k = 0; k < p_targets.size(); k++) {
		p_targets[k] = (double *)malloc(size * sizeof(double));
	}
	p_targets_blur.resize(num_targets + 1);
	for (int k = 0; k < p_targets_blur.size(); k++) {
		p_targets_blur[k] = (double *)malloc(size * sizeof(double));
	}
	return 1;
}

// display
static void draw_modes() {
	glColor3f(0.0f, 1.0f, 0.0f);
	glRasterPos2f(.005, .98);
	std::string display_mode_str;
	if (display_mode == DISPLAY_MODE_SMOKE) {
		display_mode_str = "smoke";
	}
	else {
		display_mode_str = "velocity";

	}
	for (char c : "Display mode : " + display_mode_str) {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, c);
	}
}
static void draw_velocity(void){
	double h = 1.0f / N;

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);

	for (int i = 1; i <= N; i++) {
		double x = (i - 0.5f)*h;
		for (int j = 1; j <= N; j++) {
			double y = (j - 0.5f)*h;

			glVertex2f(x, y);
			glVertex2f(x + u[IX(i, j)], y + v[IX(i, j)]);
		}
	}

	glEnd();
}
static void draw_smoke(void){
    glBegin ( GL_QUADS );
    double h = 1.0f/N;

    for (int i=0; i<=N; i++) {
        double x = (i-0.5f)*h;

        for (int j=0; j<=N; j++) {
            double y = (j-0.5f)*h;
            
            double p00 = p[IX(i,j)];
            double p01 = p[IX(i,j+1)];
            double p10 = p[IX(i+1,j)];
            double p11 = p[IX(i+1,j+1)];
            
            glColor3f(p00, p00, p00 ); glVertex2f(x, y);
            glColor3f(p10, p10, p10 ); glVertex2f(x+h, y);
            glColor3f(p11, p11, p11 ); glVertex2f(x+h, y+h);
            glColor3f(p01, p01, p01 ); glVertex2f(x, y+h);
        }
    }
    
    glEnd ();
}
static void pre_display(void) {
	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}
static void post_display(void) {
	glutSwapBuffers();
}
static void display_func(void) {
	pre_display();
	if (display_mode == DISPLAY_MODE_SMOKE) {
		draw_smoke();
	}
	else {
		draw_velocity();
	}
	draw_modes();
	post_display();
}
static void reshape_func(int width, int height) {
	glutSetWindow(win_id);
	glutReshapeWindow(width, height);
	win_x = width;
	win_y = height;
}

// keyboard inputs
static void special_key_func ( int key, int x, int y ){
    switch(key) {
        case GLUT_KEY_LEFT:
            break;
        case GLUT_KEY_RIGHT:
            break;
        case GLUT_KEY_UP:
            break;
        case GLUT_KEY_DOWN:
            break;
    }
}
static void key_func ( unsigned char key, int x, int y ){
	switch(key){
		case 'c':
		case 'C':
			clear_data();
			break;
        
		case 'q':
		case 'Q':
			free_data();
			exit(0);
			break;

		case 'v':
		case 'V':
			display_mode = (display_mode + 1) % NUM_DISPLAY_MODES;
			break;
	}
}

// mouse input
static void mouse_func ( int button, int state, int x, int y ){
    mouse_x = x;
    mouse_y = y;
    mouse_down[button] = state == GLUT_DOWN;
}
static void motion_func ( int x, int y ){
    mouse_x = x;
    mouse_y = y;
}
// relate mouse input forces
static void get_from_UI(double * d, double * u, double * v, double * temp, int curr_color_mode) {
	/*int i, j, size = (N+2)*(N+2);

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = d[i] = temp[i] = 0.0f;
	}
	if (!mouse_down[0] && !mouse_down[2]) return;

	i = (int)((       mx /(double)win_x)*N+1);
	j = (int)(((win_y-my)/(double)win_y)*N+1);
	if ( i<1 || i>N || j<1 || j>N ) return;

	if (mouse_down[0]){

		//i = (int)((       mx /(double)win_x)*N+1);
		//j = (int)(((win_y-my)/(double)win_y)*N+1);

		//if ( i<1 || i>N || j<1 || j>N ) return;

		if ( add_mode == VELOCITY_MODE ) {
			if (curr_color_mode == 3){
				u[IX(i,j)] = force * (mx-(omx));
				v[IX(i,j)] = force * ((omy)-my);
				std::cout << mx << omx << std::endl;
				std::cout << my << omy << std::endl;
				omx = mx;
				omy = my;

			}
		}

		else if ( add_mode == SMOKE_MODE ) {
			if (color_mode == curr_color_mode){
				d[IX(i,j)] = smoke_source;
			}
		}

		else if ( add_mode == TEMPERATURE_MODE ) {
			temp[IX(i, j)] = temp_source;
		}

		//omx = mx;
		//omy = my;
	}
	if (mouse_down[2]) {
		for (int start = 0; start < block_width; start++) {
			for (int end = 0; end < block_height; end++) {
				boundary[IX(i + start, j + end)] = 1.0f;
			}
		}
	}*/
	return;
}

// main loop function
static void idle_func ( void ){    
    //get_from_UI(p_prev, u_prev, v_prev);
    
    sim_step(N, dt, viscosity, diffusion, sigma, vf, vd, vg, u, u_prev, v, v_prev, p, p_prev, p_blur,
		p_targets[current_target], p_targets_blur[current_target], tmp1, tmp2);
	
	current_time += dt;
	std::cout << "Current_time : " << current_time << std::endl;
	std::cout << "Current_target : " << current_target << std::endl;
	if (current_target < num_targets && current_time > targets_time[current_target]) {
		current_target++;
	}
    
    glutSetWindow ( win_id );
    glutPostRedisplay ();
}

// open glut window
static void open_glut_window ( void ){
    glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );
    
    glutInitWindowPosition ( 0, 0 );
    glutInitWindowSize ( win_x, win_y );
    win_id = glutCreateWindow ( "Alias | wavefront" );
    
    glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
    glClear ( GL_COLOR_BUFFER_BIT );
    glutSwapBuffers ();
    glClear ( GL_COLOR_BUFFER_BIT );
    glutSwapBuffers ();
    
    pre_display ();
    
    glutKeyboardFunc ( key_func );
    glutSpecialFunc ( special_key_func );
    glutMouseFunc ( mouse_func );
    glutMotionFunc ( motion_func );
    glutReshapeFunc ( reshape_func );
    glutIdleFunc ( idle_func );
    glutDisplayFunc ( display_func );
}

// main
int main ( int argc, char ** argv ){
    glutInit ( &argc, argv );

	{
		std::ifstream cfg("C:\\Users\\first\\Desktop\\Workspace\\realtime_physics_simulation\\TargetDrivenSmokeSimulator2D\\data\\cfg.txt");

		cfg >> N >> dt >> viscosity >> diffusion;
		cfg >> sigma >> vf >> vd >> vg;
		cfg >> num_targets;
		for (int i = 0; i < num_targets; i++) {
			double tmp;
			cfg >> tmp;
			targets_time.push_back(tmp);
		}

		std::cout << "N : " << N << std::endl;
		std::cout << "dt : " << dt << std::endl;
		std::cout << "viscosity : " << viscosity << std::endl;
		std::cout << "diffusion : " << diffusion << std::endl;
		std::cout << "sigma : " << sigma << std::endl;
		std::cout << "vf : " << vf << std::endl;
		std::cout << "vd : " << vd << std::endl;
		std::cout << "vg : " << vg << std::endl;
		std::cout << "num_targets : " << num_targets << std::endl;
		std::cout << "targets_time : ";
		for (int i = 0; i < num_targets; i++) {
			std::cout << targets_time[i] << " ";
		}
		std::cout << std::endl;

		current_target = 0;
		current_time = 0.0;
	}

    if (!allocate_data()) exit(1);
    clear_data ();

	// init density from p_init.txt
	/*{
		std::ifstream p_init("C:\\Users\\first\\Desktop\\Workspace\\realtime_physics_simulation\\TargetDrivenSmokeSimulator2D\\data\\p_init.txt");
		for (int i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) {
				double tmp;
				p_init >> tmp;
				p[IX(i, j)] = tmp;
			}
		}
	}*/

	// init density as rect
	{
		int xmin = N * 1 / 8;
		int xmax = N * 2 / 8;
		int ymin = N * 1 / 8;
		int ymax = N * 7 / 8;
		double _p = 0.5;
		for (int i = xmin; i <= xmax; i++) {
			for (int j = ymin; j <= ymax; j++) {
				p[IX(i, j)] = _p;
				p_targets[0][IX(i, j)] = _p;
			}
		}
	}

	// load target density from p_target_{i}.txt
	/*{
		for (int k = 0; k < num_targets; k++) {
			std::ifstream p_target("C:\\Users\\first\\Desktop\\Workspace\\realtime_physics_simulation\\TargetDrivenSmokeSimulator2D\\data\\p_target" + std::to_string(k) + ".txt");
			for (int i = 1; i <= N; i++) {
				for (int j = 1; j <= N; j++) {
					double tmp;
					p_target >> tmp;
					p_targets[k][IX(i, j)] = tmp;
				}
			}
		}
	}*/

	// set target density as donut
	{
		int cenx = N / 2;
		int ceny = N / 2;
		double r = 0; // N / 8
		double R = N / 4;
		double _p = 0.5;
		for (int i = cenx - R; i <= cenx + R; i++) {
			for (int j = ceny - R; j <= ceny + R; j++) {
				double dx = i - cenx;
				double dy = j - ceny;
				double dr = sqrt(dx * dx + dy * dy);
				if (r <= dr && dr <= R) {
					p_targets[1][IX(i, j)] = _p;
				}
			}
		}
	}
	{
		int xmin = N * 6 / 8;
		int xmax = N * 7 / 8;
		int ymin = N * 1 / 8;
		int ymax = N * 7 / 8;
		double _p = 0.5;
		for (int i = xmin; i <= xmax; i++) {
			for (int j = ymin; j <= ymax; j++) {
				p_targets[2][IX(i, j)] = _p;
			}
		}
	}

	// pre compute p_targets_blur
	for (int k = 0; k <= num_targets; k++) {
		guassian_blur(N, p_targets_blur[k], p_targets[k], tmp1, sigma);
	}

    win_x = 512;
    win_y = 512;
    open_glut_window ();
    
    glutMainLoop ();
    
    exit ( 0 );
}
