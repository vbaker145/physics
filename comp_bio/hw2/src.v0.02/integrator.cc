
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

#include "integrator.h"

#include "ensemble.h"
#include "forces.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const integrator &z)
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

integrator::integrator(void)
{
}


integrator::~integrator(void) { }

void integrator::verlet(ensemble& sys)
{
  vec sum_vel;
  double sum_vel2 = 0.;

  vec new_r;

  for (int i=0; i<sys.N; i++)
  {
	// GET NEW POSITIONS
	new_r = 2.*sys.atoms[i].r - sys.atoms[i].r_old + sys.delta_t*sys.delta_t * sys.atoms[i].f / sys.atoms[i].m ;

	// GET VELOCITY
	sys.atoms[i].v = (new_r - sys.atoms[i].r_old) / (2.*sys.delta_t);

	// VELOCITY OF CENTER OF MASS
	sum_vel   = sum_vel + sys.atoms[i].v;
	// TOTAL KINETIC
	sum_vel2 += sys.atoms[i].v.norm2();

	// UPDATE
	sys.atoms[i].r_old = sys.atoms[i].r;
	sys.atoms[i].r     = new_r;
  }

  // UPDATE INSTANTANEOUS TEMPERATURE
  sys.kin_energy  = sum_vel2 / 2.;
#if NDIM == 2
  sys.temperature = sum_vel2 / (2.*sys.N-2.);
#else
  sys.temperature = sum_vel2 / (3.*sys.N-3.);
#endif
  // UPDATE TOTAL ENERGY
  //sys.energy = (sys.pot_energy + 0.5 * sum_vel2) / sys.N;
  //sys.energy = (sys.pot_energy + sys.kin_energy) / sys.N;
  sys.energy = sys.pot_energy + sys.kin_energy;
}

void integrator::velocity_verlet(forces &physics, ensemble& sys)
{
  // NOTE THAT THIS velocity verlet ALGORITHM IS HIGHLY INEFICIENT!
  // IT CALCULATES THE FORCE INSIDE HERE, BUT THEN THE TIME LOOP
  // RECALCULATES THIS FORCE WHILE NOT REUSING THESE NUMBERS.

  int i;
  vec sum_vel;
  double sum_vel2 = 0.;

  vec new_r;
  vec middle_v;

  for (i=0; i<sys.N; i++)
  {
	// GET NEW POSITIONS
	// NOTE THAT r_old ARE NOT USED HERE
	sys.atoms[i].r = sys.atoms[i].r + sys.delta_t*sys.atoms[i].v + sys.delta_t*sys.delta_t*sys.atoms[i].f/(2.0*sys.atoms[i].m) ;

	// GET HALF-STEP VELOCITY
	sys.atoms[i].v_new = sys.atoms[i].v + 0.5 * sys.delta_t *sys.atoms[i].f/sys.atoms[i].m;
  }

  // CALCULATE FORCES (AGAIN), ALGORITHM REQUIRES IT.
  // THIS IS INEFICIENT AS THESE SAME FORCES ARE AGAIN CALCULATED
  // OUTSIDE OF THIS LOOP - IN THE MAIN PROGRAM.
  physics.calc_LJ_forces(sys);	

  for (i=0; i<sys.N; i++)
  {
	// GET WHOLE-STEP VELOCITY USING THE NEW FORCES
	sys.atoms[i].v = sys.atoms[i].v_new + 0.5 * sys.delta_t *sys.atoms[i].f/sys.atoms[i].m;

	// VELOCITY OF CENTER OF MASS
	sum_vel   = sum_vel + sys.atoms[i].v;
	// TOTAL KINETIC
	sum_vel2 += sys.atoms[i].v.norm2();
  }

  // UPDATE INSTANTANEOUS TEMPERATURE
  sys.kin_energy  = sum_vel2 / 2.;
#if NDIM == 2
  sys.temperature = sum_vel2 / (2.*sys.N-2.);
#else
  sys.temperature = sum_vel2 / (3.*sys.N-3.);
#endif
  // UPDATE TOTAL ENERGY
  //sys.energy = (sys.pot_energy + 0.5 * sum_vel2) / sys.N;
  //sys.energy = (sys.pot_energy + sys.kin_energy) / sys.N;
  sys.energy = sys.pot_energy + sys.kin_energy;
}

void integrator::rescale_T_simple(ensemble& sys)
{
  double rescale_factor = sqrt( sys.requested_temperature / sys.temperature );
  PRINT2("rescaling by ",rescale_factor);
  for (int i=0; i<sys.N; i++)
  {
	// RESCALE VELOCITIES
	sys.atoms[i].v = sys.atoms[i].v * rescale_factor;
  }
}

void integrator::zero_velocities(ensemble& sys)
{
  for (int i=0; i<sys.N; i++)
  {
	// ZERO OUT ALL VELOCITIES
	sys.atoms[i].v.reset();
  }
}

/*
 * ----------------------------------------------------------------------
 */
