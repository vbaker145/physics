
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

#include "forces.h"

#include "ensemble.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const forces &z)
{
  // FRIEND FUNCTION TO WRITE ELEMENTS OF DATA
//os << z.N  << " "
//<< "\n"
//;

  return os;
}


/*
 * ----------------------------------------------------------------------
 *                 FORCES CLASS
 */

forces::forces(void) { }

forces::~forces(void) { }

void forces::calc_LJ_forces(ensemble& sys)
{
  // CALCULATE THE LENNARD JONES FORCES ON ALL THE PARTICLES

  // CLEAR ENERGY
  sys.pot_energy = 0.;
  // CLEAR FORCES
  sys.clear_forces();
  // CLEAR VIRIAL
  sys.virSum = 0.;

  vec dr;
  double r2, inv_r2, inv_r6;
  double force_val;

  // PAIRWISE INTERACTIONS
  for (int i=0; i<sys.N-1; i++)
	for (int j=i+1; j<sys.N; j++)
	{
	  dr = sys.atoms[i].r - sys.atoms[j].r;
	  // USE IMAGE DISTANCES
	  dr.x = IMAGE(dr.x, sys.inv_Lx, sys.Lx);
	  dr.y = IMAGE(dr.y, sys.inv_Ly, sys.Ly);
	  dr.z = IMAGE(dr.z, sys.inv_Lz, sys.Lz);

	  r2 = dr.norm2();
	  inv_r2 = 1.0 / r2;
	  inv_r6 = inv_r2 * inv_r2 * inv_r2;
	  // FORCE IN REDUCED UNITS
	  force_val = 48. * inv_r2*inv_r6 * ( inv_r6 - 0.5 );

	  sys.atoms[i].f = sys.atoms[i].f + force_val * dr;
	  sys.atoms[j].f = sys.atoms[j].f - force_val * dr;

	  // VIRIAL
	  sys.virSum += force_val * r2;
	
	  sys.pot_energy += 4.*inv_r6*(inv_r6 - 1.);
	  /*
	  if (sqrt(r2)<1.)
	    PRINT5("POT IJ ", i, j, sqrt(r2), 4.*inv_r6*(inv_r6 - 1.));
		*/
	}
}

double forces::calc_LJ_forces_single_particle(int ii, vec &direction, ensemble& sys)
{
  // CALCULATE THE LENNARD JONES FORCES ON A SINGLE GIVEN PARTICLE
  // AND RETURN THE ENERGY CONTRINBUTION DUE TO ONLY THAT PARTICLE.
  // RETURN DIRECTION OF THE FORCE IN direction

  // CLEAR ENERGY
  double pot_energy = 0.;
  // CLEAR FORCES
  sys.clear_forces(ii);

  vec dr;
  double r2, inv_r2, inv_r6;
  double force_val;

  // PAIRWISE INTERACTIONS BETWEEN PARTICLE ii
  // AND THE REST
  for (int i=0; i<sys.N; i++)
  {
	if (i!=ii)
	{
	  dr = sys.atoms[i].r - sys.atoms[ii].r;
	  // USE IMAGE DISTANCES
	  dr.x = IMAGE(dr.x, sys.inv_Lx, sys.Lx);
	  dr.y = IMAGE(dr.y, sys.inv_Ly, sys.Ly);
	  dr.z = IMAGE(dr.z, sys.inv_Lz, sys.Lz);

	  r2 = dr.norm2();
	  inv_r2 = 1.0 / r2;
	  inv_r6 = inv_r2 * inv_r2 * inv_r2;
	  // FORCE IN REDUCED UNITS
	  force_val = 48. * inv_r2*inv_r6 * ( inv_r6 - 0.5 );

	  sys.atoms[ii].f = sys.atoms[ii].f - force_val * dr;

	  pot_energy += 4.*inv_r6*(inv_r6 - 1.);
	}
  }
  // SAVE THE DIRECTION FOR RETURN 
  direction = sys.atoms[ii].f / sys.atoms[ii].f.norm(); 

  return pot_energy;
}


/*
 * ----------------------------------------------------------------------
 */
