
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

#include "infobjarray.h"

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
  	void readFile(char*);
  	void writeFile(char*);
	void apply_inits(ofstream&);		// USE PARAMETERS READ FROM FILE TO INITIALIZE

	void init_positions(int a);	// ASSIGN POSITIONS TO PARTICLES
	void init_velocities(int b);// ASSIGN VELOCITIES TO PARTICLES

	void output_properties(int counter, double real_time, ofstream&);
	void output_positions(void);
	void output_positions(int, char*);
	void  input_positions(char*);
	void output_position(int a);
	vec  get_image_position(int a);
	void output_velocities(void);
	void output_velocities(char*);
	void output_total_linear_momentum(void);
	void output_total_linear_momentum(ofstream&);
	double output_pressure(void);
	double output_temperature(void);

	void clear_forces(void);
	void clear_forces(int);

	void init_lattice(void);	// PUTS PARTICLES IN A LATTICE
	void init_random_pos_check_energy(void); // PUTS PARTICLES AT RANDOM POSITIONS
										// WHERE CHANGES IN ENERGY ARE CHECKED AGAINST NEW POSITIONS.
	void init_random_pos_two_particles(void); // PUTS TWO PARTICLES AT RANDOM
	void init_xtal_pos_two_particles(void); // PUTS TWO PARTICLES ON THE X-AXIS
	void init_xtal_polymer(void); // PUTS POLYMERS IN A XTAL STRUCTURE
	void init_from_file(void); // READS ATOM POSITIONS FROM FILE
	void init_zero_vel(void);	// GIVES ZERO TO ALL VELOCITIES
	void init_random_vel(void); // GIVES RANDOM VELS TO PARTICLES
	void init_two_part_vel(void); // VELOCITIES FOR JUST A TWO PARTICLE SYSTEM
	void init_polymers_vel(void); // VELOCITIES FOR ATOMS IN A POLYMER
	void scale_vel(void);	// SCALE VELOCITIES BY SQUARE ROOT OF INV TEMPERATURE
	void init_random_angles(void); // GIVES RANDOM ANGLES TO PARTICLES
	void adjust_box_sizes(double lx, double ly, double lz);
	void scale_positions(double sx, double sy, double sz);

	double calc_temperature(void); // UPDATES INSTANTANEOUS TEMPERATURE
	void calc_pressure(void); // UPDATES INSTANTANEOUS PRESSURE
	void calc_density(void); // UPDATES INSTANTANEOUS DENSITY FOR NPT
	double order_param_LJXY(void); // ORIENTATIONAL ORDER PARAMETER IN THE LJXY CASE
	double order_param_LJLL(void); // ORIENTATIONAL ORDER PARAMETER IN THE LJLL CASE

	void fix_atoms_beyond_radius(double radius); // FIXES ATOMS THAT ARE OUTSIDE A CIRCLE OF RADIUS radius
	long num_not_fixed;

	long N;					// TOTAL NUMBER OF ATOMS
	atom *atoms;			// ARRAY OF ATOMS

	int do_polymers;		// CONSIDER ATOMS BONDED INTO POLYMERS
	long num_chains;		// TOTAL NUMBER OF POLYMERS
	long atoms_per_chain;	// NUM ATOMS PER POLYMER
	double length_bonds;	// LENGTH OF BONDS IN THE POLYMER

    infobjarray<vec,1000> tag_particle_pos;

	void reset_tag_pos(void);
	void add_tag_pos(vec &r);
	void record_tag_particles_position(void);
	void write_tag_pos(void);
	long tag_particle;
	short do_tag;

	double total_time;		// TOTAL TIME IN REAL UNITS FOR THE SIMULATION
	double simulation_time;		// 	SIMULATION TIME IN REAL UNITS FOR THE SIMULATION
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
	double pressure;			// FIXED PRESSURE FOR THE SYSTEM, USED IN NPT
	double run_pressure;			// INSTANTANEOUS PRESSURE FOR THE SYSTEM
	double rc;				// CUTOFF DISTANCE
	double rl;				// LIST CUTOFF DISTANCE

	short init_pos;			// FLAG TO INDICATE INITIAL POSITIONS
					// 0 - LATTICE
					// 1 - RANDOM
					// 2 - TWO PARTICLES
	short init_vel;			// FLAG TO INDICATE INITIAL VELOCITY DISTRIBUTION
					// 0 - ZERO
					// 1 - RANDOM
	short dynamics;			// FLAG TO INDICATE DYNAMICS
							// 0 - REGULAR MD, N*(N-1)/2 CALCULATION
							// 1 - LJ POTENTIAL CUTOFF AT Rc
							// 2 - LJ POTENTIAL CUTOFF AT Rc WITH VERLET LIST
							// 3 - HARD SPHERES DMD
							// 4 - LANGEVIN DYNAMICS
							// 5 - LANGEVIN DYNAMICS, 2
							// 6 - MC of LJ
							// 7 - MC of LJ+XY
							// 8 - MC of LJ+LL
	double damping_constant;
	double time_meas;		// TIME OVER WHICH MEAN SQUARE DISPLACEMENT AND PAIR CORRELATION FUNCTION 
							// ARE CALCULATED
	int nhis_pcf;			// NUMBER OF HISTOGRAM BINS FOR THE PAIR CORRELATION FUNCTION
	short do_vacf;			// FLAG TO DO OR NOT VACF
	short do_pcf;			// FLAG TO DO OR NOT PCF
	short do_R;				// FLAG TO DO OR NOT END-TO-END AND RG
	double threshold_energy; 	// ENERGY THRESHOLD TO USE WHEN PUTTING RANDOM INIT POSITIONS
	double mc_dr;			// DISTANCE STEP FOR MONTE CARLO
	short do_mc_adjust;			// ADJUST SIZE OF MC STEP BY ACCEPTANCE
	double XY_J;			// STRENGTH OF XY MODEL INTERACTION
	double theta_o;			// ANGULAR STEP SIZE FOR MC
	double two_part_vel;	// FOR A TWO-PARTICLE SYSTEM, THIS IS THE UP/DOWN
							// SPEED TO TRY TO MAKE AN ORBIT

	short do_minimize;		// MINIMIZE POSITIONS OR NOT
};

#endif // ENSEMBLE_H
