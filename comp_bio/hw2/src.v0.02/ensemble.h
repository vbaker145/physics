
/*
 * Code to accompany the PHYS 462/562 Computational Biophysics course.
 * 
 * Copyright (c) 2009-2010 Luis R Cruz Cruz, Philadelphia, PA
 *
 * Please address questions to:
 *  Luis R Cruz Cruz
 *  Physics Department
 *  Drexel University
 *  3141 Chestnut St.
 *  Philadelphia, PA 19104
 */


#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <iostream>
#include <fstream>

#include "atom.h"

#include "defs.h"
#include "macros.h"

using namespace std;

class ensemble
{
  friend ostream &operator<<(ostream &os, const ensemble &z);

  public:
  	ensemble(void);
  	~ensemble(void);
  private:
  	ensemble(const ensemble &z);
  	ensemble &operator=(const ensemble &z);
  public:
  	void readFile(char*,ofstream&);

	void init_positions(int a);	// ASSIGN POSITIONS TO PARTICLES
	void init_velocities(int b);// ASSIGN VELOCITIES TO PARTICLES

	void output_properties(int counter, double real_time, ofstream&);
	void output_positions(void);
	void output_positions(char*);
	void output_position(int a);
	void output_velocities(void);
	void output_velocities(char*);
	void output_speeds(char* outfile);
	void output_total_linear_momentum(void);
	void output_total_linear_momentum(ofstream&);

	void clear_forces(void);
	void clear_forces(int);

	void init_lattice(void);	// PUTS PARTICLES IN A LATTICE
	void init_random_pos_check_energy(void); // PUTS PARTICLES AT RANDOM POSITIONS
										// WHERE CHANGES IN ENERGY ARE CHECKED AGAINST NEW POSITIONS.
	void init_random_vel(void); // GIVES RANDOM VELS TO PARTICLES

	double calc_temperature(void); // UPDATES INSTANTANEOUS TEMPERATURE
	void calc_pressure(void); // UPDATES INSTANTANEOUS PRESSURE

	long N;					// TOTAL NUMBER OF ATOMS
	atom *atoms;			// ARRAY OF ATOMS

	double total_time;		// TOTAL TIME IN REAL UNITS FOR THE SIMULATION
	long total_num_steps;	// TOTAL NUMBER OF STEPS, CORRESPONDING TO THE TOTAL TIME 
							// DIVIDED BY THE TIME STEP
	double delta_t;			// TIME STEP DT
	double density;			// TOTAL DENSITY OF THE SYSTEM
	double Lx, Ly, Lz;		// SIZE OF THE SYSTEM
	double inv_Lx, inv_Ly, inv_Lz;		// INVERSE OF SIZE OF THE SYSTEM
	double temperature;		// TEMPERATURE FOR THE SYSTEM
	double requested_temperature;		// REQUESTED TEMPERATURE FOR THE SYSTEM
	double pot_energy;			// POTENTIAL ENERGY OF THE SYSTEM
	double kin_energy;			// KINETIC ENERGY OF THE SYSTEM
	double energy;			// ENERGY OF THE SYSTEM
	double virSum;			// SUM OF THE VIRIAL, FOR THE PRESSURE CALC
	double pressure;			// PRESSURE FOR THE SYSTEM
};

#endif // ENSEMBLE_H
