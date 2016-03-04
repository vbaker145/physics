
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
#include "diffusion.h"
#include "paircorrfunc.h"
#include "R.h"
#include "verletlist.h"
#include "stepwise.h"
#include "langevin.h"
#include "montecarlo.h"

/*
 * PROGRAM TO RUN SIMULATIONS OF N-BODY PROBLEMS
 */

using namespace std;

int mainFunc(int argc, char* argv[], double temp)
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
  // - LANGEVIN INTEGRATOR
  langevin brownian;
  // - MINIMIZATION OBJECT
  minima minimize;
  // - STEPWISE INTERACTION ALGORITHM
  stepwise hardsphere;
  // - MC ALGORITHM
  montecarlo mc;
  // - OUTPUT FILES
  OUTPUT_FILES

  // =============================================
  //  READ INITIAL PARAMETERS
  sys.readFile(argv[1]);

  //VJB: modify gamma for multi-parameter runs
  //sys.damping_constant = gamma;
  sys.temperature = temp;
  out_f.close();

  char therm_fname[80];
  sprintf(therm_fname, "results/thermo_%d.dat", (int) sys.temperature);
  out_f.open((const char *) therm_fname);

  // =============================================
  //  WRITE INITIAL PARAMETERS
  char dummy8[40] = "results/input.txt";
  sys.writeFile(dummy8);
  // =============================================
  //  APPLY PARAMETERS
  sys.apply_inits(warn);
  params << sys;
  params.close();

  // =============================================
  // SET SEED FOR RANDOM GENERATOR
  SET_SEED 

  // =============================================
  // INITIALIZE POSITIONS : 0 - LATTICE, 1 - RANDOM WITH ENERGY CHECK
  sys.init_positions(sys.init_pos);
  
  // =============================================
  // VERLET NEIGHBOR LIST
  verletlist neighborsverlet(sys.N, sys.rc, sys.rl);

  // =============================================
  // OPTIONAL OUTPUTS TO VERIFY INITIAL SETUP 
  //
  char dummy[40] = "results/out_init_r.xyz";
  sys.output_positions(1,dummy);
  //sys.output_position(5);
  //
  //char dummy2[40] = "results/out_init_v.dat";
  //sys.output_velocities(dummy2);
  //
  //sys.output_total_linear_momentum();
  //
  //char dummy7[40] = "results/enrg_lndscp.dat";
  //physics.energy_landscape(sys,dummy7);

  // =============================================
  // DEFINE CONFINED SYSTEM
#if CONFINED == 1
  double min_L = MIN(sys.Lx,sys.Ly);
  min_L = MIN(min_L,sys.Lz);
  double rad;
  //rad = min_L + 5.;
  //rad = min_L / 2.;
  //rad = min_L / 2.5;
  //rad = min_L / 4.;
  //rad = min_L / 6.;
  //rad = min_L / 8.;
  rad = 3.75;
  sys.fix_atoms_beyond_radius(rad);
  PRINT2("Confined to radius = ",rad);
#endif

  // =============================================
  // OPTIONAL MINIMIZATION STEP
  int    particle  = 5;
  double step_size = 0.1;
  //double step_size = 0.5;
  double criterium = 1.0e-4;
  //double criterium = 1.0e-6;
  // - MOVE ONE PARTICLE
  //minimize.steepest_descents_arbitrary_step_one_particle(particle, step_size, criterium, physics, sys);
  // - MOVE ALL PARTICLES
  //int steps_min = 1; 
  int steps_min = 10; 
  if (sys.do_minimize)
  {
	for (i=0; i<steps_min; i++)
	  minimize.steepest_descents_arbitrary_step(step_size, criterium, physics, sys);
  }
  /*
  */
  // - PRINT POSITIONS AFTER MINIMIZING
  char dummy3[40] = "results/out_minim_r.xyz";
  sys.output_positions(1,dummy3);

  // =============================================
  // INITIALIZE RANDOM VELOCITIES
  sys.init_velocities(sys.init_vel);
  //char dummy6[40] = "results/out_init_v.dat";
  //sys.output_velocities(dummy6);

  // =============================================
  // CREATE DIFFUSION OBJECT
  diffusion diff(sys);
  // CREATE PAIR CORRELATION FUNCTION OBJECT
  paircorrfunc pcf(sys);
  // CREATE END-TO-END AND RG FUNCTION OBJECT
  R_meas R(sys);

  // =============================================
  // UPDATE VERLET LISTS FOR THE FIRST STEP
  neighborsverlet.update_verletlist(sys);

  // =============================================
  // INITIALIZE TAG PARTICLE ARRAY
  sys.reset_tag_pos();

  // =============================================
  // MAIN TIME LOOP 
  long n       = 0;
  sys.simulation_time = 0.;
  int turn_off = 0;

#if BENCH_IT == 1
  // **************************
  //	FOR MEASURING EXECUTION TIME
  long bench_t = time(NULL);
  int  bench_h = 31;
#endif

  //while (n < sys.total_num_steps && !turn_off)
  while (sys.simulation_time <= sys.total_time && !turn_off)
  {
	// **************************
	// CHOOSE STEPS ACCORDING TO DYNAMICS
	switch(sys.dynamics)
	{
	  case 0:
		{
		  // ##################################
		  // REGULAR LJ MD, N*(N-1)/2 CALCULATION
		  // ##################################

		  // **************************
		  // CALCULATE FORCES
		  double virSum;
		  physics.calc_LJ_forces(virSum, sys);

		  // **************************
		  // INTEGRATE EQ. OF MOTION
		  // 	- VERLET ALGORITHM
		  integrate.verlet(sys);
		  // 	- VELOCITY VERLET ALGORITHM
		  //integrate.velocity_verlet(physics, sys);

		  break;
		}
	  case 1:
		{
		  // ##################################
		  // LJ POTENTIAL CUTOFF AT Rc
		  // ##################################

		  // **************************
		  // CALCULATE FORCES
		  //physics.calc_LJ_forces_pot_cutoff(sys);
		  // **************************
		  // CALCULATE FORCES HERE ONLY THE FIRST TIME
		  if (n==0)
		  {
			physics.calc_LJ_forces_pot_cutoff(sys);
		  }

		  // **************************
		  // INTEGRATE EQ. OF MOTION PT 1	
		  // 	- VERLET ALGORITHM
		  //integrate.verlet(sys);
		  // 	- VELOCITY VERLET ALGORITHM
		  integrate.velocity_verlet_pt1(physics, sys);

		  // **************************
		  //  CALCULATE FORCES IN BETWEEN
		  physics.calc_LJ_forces_pot_cutoff(sys);

		  // **************************
		  // INTEGRATE EQ. OF MOTION PT 2	
		  // 	- VELOCITY VERLET ALGORITHM
		  integrate.velocity_verlet_pt2(physics, sys);

		  break;
		}
	  case 2:
		{
		  // ##################################
		  // LJ POTENTIAL CUTOFF AT Rc
		  // WITH VERLET LISTS
		  // ##################################

		  // **************************
		  // CALCULATE FORCES
		  physics.calc_LJ_forces_verletlists(neighborsverlet, sys);

		  // **************************
		  // INTEGRATE EQ. OF MOTION
		  // 	- VERLET ALGORITHM
		  integrate.verlet(sys);
		  // 	- VELOCITY VERLET ALGORITHM
		  //integrate.velocity_verlet(physics, sys);

		  // **************************
		  // CHECK IF NEED TO UPDATE VERLET LISTS
		  //if (1)
		  if (neighborsverlet.check_verletlist(sys))
		  //if ( neighborsverlet.check_verletlist(sys) || (n%3==0) )
		  {
			//PRINT2("in update",n);
			neighborsverlet.update_verletlist(sys);
		  }
		  /*
		   */

		  break;
		}
	  case 3:
		{
		  // ##################################
		  // HARD SPHERES DMD
		  // ##################################

		  // **************************
		  // USE HARD SPHERE STEPWISE ALGORITHM
		  hardsphere.calc_hardsphere_stepwise(sys);

		  break;
		}
	  case 4:
		{
		  // ##################################
		  // LANGEVIN DYNAMICS
		  // ##################################

		  // **************************
		  // CALCULATE FORCES
		  physics.calc_LJ_forces_pot_cutoff(sys);
		  //physics.calc_LJ_forces_verletlists(neighborsverlet, sys);

		  // **************************
		  // INTEGRATE EQ. OF MOTION
		  //  - LANGEVIN DYNAMICS
		  brownian.integrate_ermak(sys);

		  // **************************
		  // CHECK IF NEED TO UPDATE VERLET LISTS
		  /*
		  if (neighborsverlet.check_verletlist(sys))
		  {
			//PRINT2("in update",n);
			neighborsverlet.update_verletlist(sys);
		  }
		   */

		  break;
		}
	  case 5:
		{
		  // ##################################
		  // LANGEVIN DYNAMICS, INTEGRATION METHOD
		  // WITH FORCE NOT CONSTANT
		  // ##################################

		  // **************************
		  // CALCULATE FORCES HERE ONLY THE FIRST TIME
		  if (n==0)
		  {
			physics.calc_LJ_forces_pot_cutoff(sys);
			//physics.calc_LJ_forces_verletlists(neighborsverlet, sys);
		  }

		  // **************************
		  // INTEGRATE EQ. OF MOTION
		  //  - LANGEVIN DYNAMICS STEP 1
		  //    UPDATE POSITIONS
		  brownian.integrate_other_pt1(sys);

		  // **************************
		  //  CALCULATE FORCES IN BETWEEN
		  physics.calc_LJ_forces_pot_cutoff(sys);

		  // **************************
		  // INTEGRATE EQ. OF MOTION
		  //  - LANGEVIN DYNAMICS STEP 2
		  //    UPDATE VELOCITIES
		  brownian.integrate_other_pt2(sys);

		  // **************************
		  // CHECK IF NEED TO UPDATE VERLET LISTS
		  /*
		  if (neighborsverlet.check_verletlist(sys))
		  {
			//PRINT2("in update",n);
			neighborsverlet.update_verletlist(sys);
		  }
		   */

		  break;
		}
	  case 6:
		{
		  // ##################################
		  // MC OF LJ PARTICLES, NO ROTATION
		  // ##################################

		  // **************************
		  // INITIALIZE THE POTENTIAL ENERGY BY 
		  // CALCULATING FORCES BUT ONLY ON THE FIRST STEP
		  if (n==0)
		  {
			double virSum;
			physics.calc_LJ_forces(virSum, sys);
			// USING SHIFTED POTENTIAL INTRODUCES ALREADY
			// THE SHIFT AS AN ADDITIONAL TERM IN THE INITIAL
			// ENERGY
			//physics.calc_LJ_forces_pot_cutoff(sys);
		  }

		  // **************************
		  // DO MC STEP
		  double acc;
		  acc = mc.mc_step(physics, sys);
		  cout <<  "accepted = " << acc 
		       <<  " mc_step = " << sys.mc_dr 
			<< "\n";

		  // **************************
		  // ADJUST MC STEP IF REQUESTED
		  if (sys.do_mc_adjust)
			 sys.mc_dr = mc.adjust_mc_step_size(acc, sys.mc_dr);

		  break;
		}
	  case 7:
		{
		  // ##################################
		  // MC OF LJ + XY MODEL PARTICLES
		  // ##################################

		  // **************************
		  // INITIALIZE THE POTENTIAL ENERGY BY 
		  // CALCULATING FORCES BUT ONLY ON THE FIRST STEP
		  if (n==0)
		  {
			// INITIALIZE ANGLES
			sys.init_random_angles();

			// THIS IS IGNORING THE XY MODEL TERM!
			double virSum;
			physics.calc_LJ_forces(virSum, sys);
			// USING SHIFTED POTENTIAL INTRODUCES ALREADY
			// THE SHIFT AS AN ADDITIONAL TERM IN THE INITIAL
			// ENERGY
			//physics.calc_LJ_forces_pot_cutoff(sys);
		  }

		  // **************************
		  // DO MC STEP
		  double acc;
		  acc = mc.mc_step_2D_angle(0, physics, sys);
		  cout <<  "accepted = " << acc 
		       <<  " mc_step = " << sys.mc_dr 
		       <<  " mc_theta = " << DEGREE(sys.theta_o)
			<< "\n";

		  // **************************
		  // ADJUST MC STEP IF REQUESTED
		  if (sys.do_mc_adjust)
			 sys.mc_dr = mc.adjust_mc_step_size_theta(acc, sys.mc_dr, sys.theta_o);

		  break;
		}
	  case 8:
		{
		  // ##################################
		  // MC OF LJ + LL MODEL PARTICLES
		  // ##################################

		  // **************************
		  // INITIALIZE THE POTENTIAL ENERGY BY 
		  // CALCULATING FORCES BUT ONLY ON THE FIRST STEP
		  if (n==0)
		  {
			// INITIALIZE ANGLES
			sys.init_random_angles();

			// THIS IS IGNORING THE LL MODEL TERM!
			double virSum;
			physics.calc_LJ_forces(virSum, sys);
			// USING SHIFTED POTENTIAL INTRODUCES ALREADY
			// THE SHIFT AS AN ADDITIONAL TERM IN THE INITIAL
			// ENERGY
			//physics.calc_LJ_forces_pot_cutoff(sys);
		  }

		  // **************************
		  // DO MC STEP
		  double acc;
		  acc = mc.mc_step_2D_angle(1, physics, sys);
		  cout <<  "accepted = " << acc 
		       <<  " mc_step = " << sys.mc_dr 
		       <<  " mc_theta = " << DEGREE(sys.theta_o)
			<< "\n";

		  // **************************
		  // ADJUST MC STEP IF REQUESTED
		  if (sys.do_mc_adjust)
			 sys.mc_dr = mc.adjust_mc_step_size_theta(acc, sys.mc_dr, sys.theta_o);

		  break;
		}
	  case 9:
		{
		  // ##################################
		  // MC OF LJ PARTICLES, NO ROTATION
		  // NPT ENSEMBLE
		  // ##################################

		  // **************************
		  // INITIALIZE THE POTENTIAL ENERGY BY 
		  // CALCULATING FORCES BUT ONLY ON THE FIRST STEP
		  if (n==0)
		  {
			double virSum;
			physics.calc_LJ_forces(virSum, sys);
			// USING SHIFTED POTENTIAL INTRODUCES ALREADY
			// THE SHIFT AS AN ADDITIONAL TERM IN THE INITIAL
			// ENERGY
			//physics.calc_LJ_forces_pot_cutoff(sys);
		  }

		  // **************************
		  // DO MC STEP
		  double acc;
		  acc = mc.mc_step_NPT(physics, sys);
		  cout <<  "accepted = " << acc 
		       <<  " mc_step = " << sys.mc_dr 
			<< "\n";

		  // **************************
		  // ADJUST MC STEP IF REQUESTED
		  if (sys.do_mc_adjust)
			 sys.mc_dr = mc.adjust_mc_step_size(acc, sys.mc_dr);

		  break;
		}
	  case 10:
		{
		  // ##################################
		  // POLYMER DYNAMICS WITH
		  // REGULAR LJ MD, N*(N-1)/2 CALCULATION
		  // ##################################

		  // **************************
		  // CALCULATE FORCES
		  double virSum;
		  physics.calc_LJ_forces(virSum, sys);

		  // **************************
		  // INTEGRATE EQ. OF MOTION
		  // 	- VERLET ALGORITHM
		  //integrate.verlet(sys);
		  // 	- VELOCITY VERLET ALGORITHM
		  //integrate.velocity_verlet(physics, sys);
		  // 	- VERLET ALGORITHM USING SHAKE CONSTRAINTS
		  integrate.verlet_polymers(virSum, sys);

		  break;
		}
	  case 11:
		{
		  // ##################################
		  // POLYMER DYNAMICS WITH
		  // LJ POTENTIAL CUTOFF AT Rc
		  // WITH VERLET LISTS
		  // ##################################

		  // **************************
		  // CALCULATE FORCES
		  physics.calc_LJ_forces_verletlists(neighborsverlet, sys);

		  // **************************
		  // INTEGRATE EQ. OF MOTION
		  // 	- VERLET ALGORITHM
		  //integrate.verlet(sys);
		  // 	- VELOCITY VERLET ALGORITHM
		  //integrate.velocity_verlet(physics, sys);
		  // 	- VERLET ALGORITHM USING SHAKE CONSTRAINTS
		  double virSum;
		  integrate.verlet_polymers(virSum, sys);

		  // **************************
		  // CHECK IF NEED TO UPDATE VERLET LISTS
		  if (neighborsverlet.check_verletlist(sys))
		  //if (neighborsverlet.check_verletlist(sys) || (n%3==0))
		  {
			//PRINT2("in update",n);
			neighborsverlet.update_verletlist(sys);
		  }
		  /*
		   */

		  break;
		}
	  default:
		PRINT2("ERROR: NO CHOSEN DYNAMICS",sys.dynamics);
		exit(0);
	}


	// **************************
	// UPDATE PRESSURE
	sys.calc_pressure();

	// **************************
	// PERIODIC OUTPUT OF THERMODYNAMIC PROPERTIES
	sys.output_properties(n, sys.simulation_time, out_f);
	//sys.output_properties(n, n*sys.delta_t, out_f);
	//sys.output_position(5);
	sys.output_total_linear_momentum(out_p);
	
	// **************************
	// CALCULATE MEAN SQUARE DISP. AND VACF
	if (sys.do_vacf)
		diff.accum(n,sys);
	// **************************
	// CALCULATE PAIR CORRELATION FUNCTION
	if (sys.do_pcf)
		pcf.accum(n,sys);
	// **************************
	// CALCULATE END-TO-END AND RG FUNCTION
	if (sys.do_R)
		R.calc(sys,out_rg,out_rg2,out_eoe,out_eoe_vec,out_eoe2);

	// **************************
	// RECORD TAG PARTICLE'S LOCATION
	if (sys.do_tag)
	  sys.record_tag_particles_position();

	// **************************
	// PRINT TIME
	//if (n % 10 ==0) PRINT5("time = ",n*sys.delta_t,": (",n*sys.delta_t/sys.total_time*100.,"%)");
	int print_time = 10;
	//int print_time = 100;
	//int print_time = 1000;
	if (sys.dynamics != 3)
		if (n % print_time ==0) PRINT5("time = ",sys.simulation_time,": (",sys.simulation_time/sys.total_time*100.,"%)");

	// **************************
	// INCREMENT LOOP COUNTER
	n++;
	// INCREMENT SIMULATION TIME
	// - EXCLUDE FOR DMD THAT CALCULATES ITS OWN
	// - SIMULATION TIME
	if(sys.dynamics!=3)
		sys.simulation_time = n*sys.delta_t;

#if BENCH_IT == 1
	// **************************
	//	OUTPUT EXECUTION TIME 
	if (n % bench_h == 0)

	{
	  cout << "  RUNNING AT " 
		<< float(time(NULL) - bench_t) / float(bench_h*sys.delta_t) 
		<< " seconds per unit time\n";
	  bench_t = time(NULL);
	}
#endif
  }

  // =============================================
  // LAST PRINTOUTS BEFORE EXITING
  //
  //char dummy4[40] = "results/out_final_v.dat";
  //sys.output_velocities(dummy4);
  //
  //char dummy5[40] = "results/out_final_r.xy";
  //char dummy5[40] = "results/out_final_r.xyz";
  //sys.output_positions(dummy5);
  /*
  cout 
	<< sys.output_temperature()  << " "
	<< sys.output_pressure()  << " "
	<< "\n";
	*/
  sys.write_tag_pos();

  warn << "-- Finished! --\n";

  warn.close();
  out_f.close();
  out_p.close();
  return 1;
}

//Refactored code to automate multi-parameter runs
int main(int argc, char* argv[] ) {

	for( int i=1; i<=10; i++ ) {
		mainFunc(argc, argv, (double) i);
	}

}

