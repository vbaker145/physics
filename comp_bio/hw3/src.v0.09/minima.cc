
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

#include <fstream>
#include <math.h>

#include "defs.h"
#include "macros.h"
#include "wait.h"

#include "minima.h"

#include "ensemble.h"
#include "forces.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const minima &z)
{
  // FRIEND FUNCTION TO WRITE ELEMENTS OF DATA
//os << z.N  << " "
//<< "\n"
//;

  return os;
}


/*
 * ----------------------------------------------------------------------
 *                 INTEGRATOR CLASS
 */

minima::minima(void)
{
}


minima::~minima(void) { }

void minima::steepest_descents_arbitrary_step_one_particle(int part, double step_size, double crit, forces &physics, ensemble& sys)
{
  // STEEPEST DESCENT WHERE WE MOVE A GIVEN PARTICLE IN SEARCH OF A LOCAL MINIMA

  int i,j;
  int n = 0;
  int keep_going = 1;
  vec direction;
  vec new_position;
  double curr_energy;
  double prev_energy;
  double virSum;

  while (keep_going)
  {
		curr_energy = physics.calc_LJ_forces_single_particle(part, direction, virSum, sys);

		// INITIALIZE THE PREVIOUS ENERGY VARIABLE
		// IF THIS IS THE FIRST TIME AROUND
		if (n==0) prev_energy = curr_energy + 1.;
		n = 1;

		//cout << "direction " << direction << "\n";
	 	//PRINT3("energy",curr_energy,prev_energy);
	 	//PRINT2("local energy",curr_energy);
		
		// INCREMENT POSITION IN THE DIRECTION OF direction
		// WHICH SHOULD BE MINIMIZING THE ENERGY
		new_position           = sys.atoms[part].r + step_size * direction;
		sys.atoms[part].r_old  = sys.atoms[part].r;
		sys.atoms[part].r      = new_position;

		// UPDATE STEP SIZE IF START OSCILLATING AROUND A MINIMA
		if (curr_energy > prev_energy) step_size *= 0.9;

		// STOP IF ACHIEVE DESIRED CONVERGENCE
		if (ABS(curr_energy-prev_energy)<crit) keep_going = 0;

		// SAVE VALUES OF THE ENERGY FOR NEXT ITERATION
		prev_energy = curr_energy;
  }

}

void minima::steepest_descents_arbitrary_step(double step_size, double crit, forces &physics, ensemble& sys)
{
  // STEEPEST DESCENT WHERE WE MOVE ALL PARTICLES IN SEARCH OF A MINIMA
  
  int i,j;
  int n = 0;
  int keep_going = 1;
  int part;
  vec direction;
  vec new_position;
  double curr_energy;
  double prev_energy;
  double virSum;

  wait(1);
  // LOOP OVER ALL PARTICLES AND MOVE SUCH THAT
  // WE GET LOWER ENERGY
  PRINT1("MINIMIZING:");
  for (part=0; part<sys.N; part++)
  {
	// EXCLUDE FIXED PARTICLES
	if (sys.atoms[part].is_not_fixed())
	{
	  // INITIALIZE
	  n          = 0;
	  keep_going = 1;

	  // PRINT PROGRESS
	  if ( (part % 10 == 0) && part>0 )
	  {
	  	PRINT2(100.*double(part+1)/double(sys.N),"%");
	  }
	  else
	  {
		cout << ".";
	  }

	  // MOVE GIVEN PARTICLE
	  while (keep_going)
	  {
		curr_energy = physics.calc_LJ_forces_single_particle(part, direction, virSum, sys);
		//PRINT2("energy",curr_energy);
		//PRINT4("particle",part,"energy",curr_energy);

		// INITIALIZE THE PREVIOUS ENERGY VARIABLE
		// IF THIS IS THE FIRST TIME AROUND
		if (n==0) prev_energy = curr_energy + 1.;
		n = 1;

		// INCREMENT POSITION IN THE DIRECTION OF direction
		// WHICH SHOULD BE MINIMIZING THE ENERGY
		new_position           = sys.atoms[part].r + step_size * direction;
		sys.atoms[part].r_old  = sys.atoms[part].r;
		sys.atoms[part].r      = new_position;

		// UPDATE STEP SIZE IF START OSCILLATING AROUND A MINIMA
		if (curr_energy > prev_energy) step_size *= 0.9;

		// STOP IF ACHIEVE DESIRED CONVERGENCE
		if (ABS(curr_energy-prev_energy)<crit) keep_going = 0;

		// SAVE VALUES OF THE ENERGY FOR NEXT ITERATION
		prev_energy = curr_energy;
	  }
	}

	double total_energy;
	// PRINT OUT TOTAL ENERGY (FOR THE WHOLE SYSTEM)
    total_energy = physics.calc_LJ_forces(virSum, sys);
	//PRINT2("TOTAL ENERGY",total_energy);
  }
  PRINT1("done minimizing!");
}


/*
 * ----------------------------------------------------------------------
 */
