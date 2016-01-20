
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

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "atom.h"
#include "ensemble.h"
#include "integrator.h"
#include "minima.h"
#include "forces.h"

/*
 * PROGRAM TO RUN SIMULATIONS OF N-BODY PROBLEMS
 */

using namespace std;

int main(int argc, char* argv[])
{
  // =============================================
  // CHECK FOR CORRECT NUMBER OF INPUT PARAMETERS,
  // EXIT IF MISSING PARAMETERS
  if (argc<2)
  {
	cout << "Usage:\n";
	cout << "\trun_me "
	  << " <InputFile> "
	  << "\n\n";
	exit(1);
  }

  int i,j;
  // =============================================
  // CREATE AND INITIALIZE:
  // - THE SYSTEM
  ensemble sys;
  // - ACTION OBJECT
  forces physics;
  // - FINITE DIFF. INTEGRATOR
  integrator integrate;
  // - MINIMIZATION OBJECT
  minima minimize;
  // - OUTPUT FILES
  ofstream warn("warnings.txt");
  ofstream out_f("thermo.dat");

  // =============================================
  //  READ INITIAL PARAMETERS
  sys.readFile(argv[1],warn);

  // =============================================
  // SET SEED FOR RANDOM GENERATOR
  long seed;
  long time_sub = 1234234670;
  // USE TIME AS SEED
  seed = - long( time(NULL) - time_sub );
  // SEED THE GENERATOR
  srand(seed);

  // =============================================
  // INITIALIZE POSITIONS : 1 - LATTICE, 2 - RANDOM WITH ENERGY CHECK
  sys.init_positions(1);
  
  // =============================================
  // OPTIONAL OUTPUTS TO VERIFY INITIAL SETUP 
  //
  //char dummy[40] = "out_init_r.xyz";
  //sys.output_positions(dummy);
  //sys.output_position(5);
  //
  //char dummy2[40] = "out_init_v.dat";
  //sys.output_velocities(dummy2);
  //
  //sys.output_total_linear_momentum();

  // =============================================
  // OPTIONAL MINIMIZATION STEP
  //int    particle  = 5;
  //double step_size = 0.1;
  //double criterium = 1.0e-4;
  // - MOVE ONE PARTICLE
  //minimize.steepest_desdents_arbitrary_step_one_particle(particle, step_size, criterium, physics, sys);
  // - MOVE ALL PARTICLES
  //for (i=0; i<20; i++)
  //minimize.steepest_desdents_arbitrary_step(step_size, criterium, physics, sys);
  // - PRINT POSITIONS AFTER MINIMIZING
  //char dummy3[40] = "out_minim_r.xyz";
  //sys.output_positions(dummy3);

  // =============================================
  // INITIALIZE RANDOM VELOCITIES
  sys.init_velocities(1);

  // =============================================
  // PRINT INITIAL SPEEDS
  char dummy6[40] = "out_initial_v_3D.dat";
  sys.output_speeds(dummy6);
  //sys.output_velocities(dummy6);

  // =============================================
  // MAIN TIME LOOP 
  long n       = 0;
  int turn_off = 0;
  while (n < sys.total_num_steps && !turn_off)
  {
	// **************************
	// CALCULATE FORCES
	physics.calc_LJ_forces(sys);

	// **************************
	// INTEGRATE EQ. OF MOTION
	// 	- VERLET ALGORITHM
	integrate.verlet(sys);
	// 	- VELOCITY VERLET ALGORITHM
	//integrate.velocity_verlet(physics, sys);

	// **************************
	// UPDATE PRESSURE
	sys.calc_pressure();

	// **************************
	// PERIODIC OUTPUT OF THERMODYNAMIC PROPERTIES
	sys.output_properties(n, n*sys.delta_t, out_f);

	// **************************
	// PRINT TIME
	if (n % 10 ==0) PRINT5("time = ",n*sys.delta_t,": (",n*sys.delta_t/sys.total_time*100.,"%)");

	// **************************
	// INCREMENT LOOP COUNTER
	n++;
  }

  // =============================================
  // LAST PRINTOUTS BEFORE EXITING
  //
  char dummy4[40] = "out_final_v_3D.dat";
  sys.output_speeds(dummy4);
  //sys.output_velocities(dummy4);
  //
  char dummy5[40] = "out_final_r.xyz";
  sys.output_positions(dummy5);

  warn << "-- Finished! --\n";

  warn.close();
  out_f.close();
  return 1;
}

