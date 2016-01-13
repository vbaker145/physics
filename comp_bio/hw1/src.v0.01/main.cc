
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
#include <stdlib.h>

#include "atom.h"
#include "ensemble.h"
#include "integrator.h"
#include "forces.h"

/*
 * PROGRAM TO RUN SIMULATIONS OF N-BODY PROBLEMS
 */

using namespace std;

int main(int argc, char* argv[])
{
  // =============================================
  // CHECK FOR CORRECT NUMBER OF INPUT PARAMETERS
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
  // CREATE THE SYSTEM
  ensemble sys;
  // CREATE ACTION OBJECT
  forces physics;
  // CREATE FINITE DIFF. INTEGRATOR
  integrator integrate;

  // =============================================
  //  READ SIMULATION INITIAL PARAMETERS INTO THE DATA OBJECT	
  sys.readFile(argv[1]);

  // =============================================
  // SET SEED FOR RANDOM GENERATOR
  long seed;
  long time_sub = 1234234670;
  // USE TIME AS SEED
  seed = - long( time(NULL) - time_sub );
  // SEED THE GENERATOR
  srand(seed);

  // =============================================
  // INITIALIZE POSITIONS AND VELOCITIES
  //   sys.init(i,j):
  //   	i = 1 : PUT PARTICLES IN A LATTICE
  //   	i = 2 : PUT PARTICLES IN RANDOM POSITIONS, BUT CHECK WITH ENERGY
  //   	===================
  //   	j = 1 : ASSIGN RANDOM VELOCITIES
  //   	j = 2 :
  sys.init(1,1);
  //sys.init(2,1);
  
  // OPTIONAL OUTPUTS TO VERIFY INITIAL SETUP 
  //cout << sys;
  //sys.output_positions();
  //sys.output_total_linear_momentum();


  // =============================================
  // MAIN TIME LOOP 
  long n = 0;
  while (n < sys.total_num_steps)
  {
	physics.calc_LJ_forces(sys);
	integrate.verlet(sys);

	// CHOOSE YOUR OUTPUT BY UNCOMMENTING
	// A LINE
	sys.output_properties(n, n*sys.delta_t);

	// .. OF ALL ATOMS
	//sys.output_positions();
	// .. OF ONLY ONE ATOM
	//sys.output_positions(2);
	
    	//sys.output_total_linear_momentum();

	n++;
  }

  return 1;
}

