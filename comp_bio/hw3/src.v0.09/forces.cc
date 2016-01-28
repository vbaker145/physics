
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
#include "verletlist.h"

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

double forces::calc_LJ_forces(double &virSum, ensemble& sys)
{
  // CALCULATE THE LENNARD JONES FORCES ON ALL THE PARTICLES
  // N*(N-1)/2 CALCULATION

  // CLEAR ENERGY
  sys.pot_energy = 0.;
  // CLEAR FORCES
  sys.clear_forces();
  // CLEAR VIRIAL
  sys.virSum = 0.;
  virSum     = 0.;

  vec dr;
  double r2, inv_r2, inv_r6;
  double force_val;

  // PAIRWISE INTERACTIONS
  for (int i=0; i<sys.N-1; i++)
	for (int j=i+1; j<sys.N; j++)
	{
	  // EXCLUDE FORCE CALCULATION IF BOTH PARTICLES
	  // ARE FIXED IN SPACE
	  if ( sys.atoms[i].is_not_fixed() || sys.atoms[j].is_not_fixed() )
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
	  }
	}
  virSum = sys.virSum;

  return sys.pot_energy;
}

void forces::calc_LJ_forces_pot_cutoff(ensemble& sys)
{
  // CALCULATE THE LENNARD JONES FORCES ON ALL THE PARTICLES.
  // USE A CUTOFF FOR THE INTERACTIONS AND COMPENSATE BY
  // ADDING A CONSTANT TO THE ENERGY.
  // THIS IS THE 'SHIFTED POTENTIAL' TRICK.
  // N*Nc/2 CALCULATION.
  int i,j;

  // SAVE FORCES AS OLD
  for (i=0; i<sys.N; i++)
	sys.atoms[i].f_old = sys.atoms[i].f;

  // CLEAR ENERGY
  sys.pot_energy = 0.;
  // CLEAR FORCES
  sys.clear_forces();
  // CLEAR VIRIAL
  sys.virSum = 0.;

  vec dr;
  double r2, inv_r2, inv_r6;
  double force_val;
  double rc2  = sys.rc * sys.rc;
  double ecut = 4.*( 1./pow(rc2,6) - 1./pow(rc2,3) );

  // PAIRWISE INTERACTIONS
  for (i=0; i<sys.N-1; i++)
	for (j=i+1; j<sys.N; j++)
	{
	  // EXCLUDE FORCE CALCULATION IF BOTH PARTICLES
	  // ARE FIXED IN SPACE
	  if ( sys.atoms[i].is_not_fixed() || sys.atoms[j].is_not_fixed() )
	  {
		dr = sys.atoms[i].r - sys.atoms[j].r;
		// USE IMAGE DISTANCES
		dr.x = IMAGE(dr.x, sys.inv_Lx, sys.Lx);
		dr.y = IMAGE(dr.y, sys.inv_Ly, sys.Ly);
		dr.z = IMAGE(dr.z, sys.inv_Lz, sys.Lz);

		r2 = dr.norm2();

		if (r2 < rc2)
		{
		  inv_r2 = 1.0 / r2;
		  inv_r6 = inv_r2 * inv_r2 * inv_r2;
		  // FORCE IN REDUCED UNITS
		  force_val = 48. * inv_r2*inv_r6 * ( inv_r6 - 0.5 );

		  sys.atoms[i].f = sys.atoms[i].f + force_val * dr;
		  sys.atoms[j].f = sys.atoms[j].f - force_val * dr;

		  // VIRIAL
		  sys.virSum += force_val * r2;

		  // CORRECTED POTENTIAL FOR THE AMOUNT AT V(rc)
		  sys.pot_energy += 4.*inv_r6*(inv_r6 - 1.) - ecut;
		}
	  }

	}
}

void forces::calc_LJ_forces_verletlists(verletlist &vl, ensemble& sys)
{
  // CALCULATE THE LENNARD JONES FORCES ON ALL THE PARTICLES
  // IN THE VERLET NEIGHBORHOOD
  // N*Nv CALCULATION, WHERE Nv NUMBER OF PARTICLES INSIDE
  // THE SKIN SPHERE CUTOFF

  // CLEAR ENERGY
  sys.pot_energy = 0.;
  // CLEAR FORCES
  sys.clear_forces();
  // CLEAR VIRIAL
  sys.virSum = 0.;

  int i,j;
  vec dr;
  double r2, inv_r2, inv_r6;
  double force_val;
  //double rv2  = vl.r_skin*vl.r_skin;
  //double ecut = 4.*( 1./pow(rv2,6) - 1./pow(rv2,3) );
  double rc2  = sys.rc * sys.rc;
  double ecut = 4.*( 1./pow(rc2,6) - 1./pow(rc2,3) );

  // PAIRWISE INTERACTIONS.
  // LOOP OVER VERLET NEIGHBORS
  for (i=0; i<sys.N; i++)
	for (j=0; j<vl.number_list[i]; j++)
	{
	  // EXCLUDE FORCE CALCULATION IF BOTH PARTICLES
	  // ARE FIXED IN SPACE
	  if ( sys.atoms[i].is_not_fixed() || sys.atoms[j].is_not_fixed() )
	  {
		dr = sys.atoms[i].r - sys.atoms[vl.list[i][j]].r;

		// USE IMAGE DISTANCES
		dr.x = IMAGE(dr.x, sys.inv_Lx, sys.Lx);
		dr.y = IMAGE(dr.y, sys.inv_Ly, sys.Ly);
		dr.z = IMAGE(dr.z, sys.inv_Lz, sys.Lz);

		r2 = dr.norm2();

		// ONLY INTERACT WITH PARTICLES CLOSER THAN rc
		// EVEN WHEN THEY ARE INSIDE THE LIST.
		if (r2 < rc2)
		{
		  inv_r2 = 1.0 / r2;
		  inv_r6 = inv_r2 * inv_r2 * inv_r2;
		  // FORCE IN REDUCED UNITS
		  force_val = 48. * inv_r2*inv_r6 * ( inv_r6 - 0.5 );

		  sys.atoms[i].f = sys.atoms[i].f + force_val * dr;
		  // BECAUSE EACH PARTICLE HAS A LIST, HAVE TO
		  // PREVENT DOUBLE COUNTING, SO NO NEWTON'S 3RD
		  // ADVANTAGE
		  //sys.atoms[j].f = sys.atoms[j].f - force_val * dr;

		  // VIRIAL
		  //sys.virSum += force_val * r2;
		  // THE 1/2 TO ACCOUNT FOR DOUBLE COUNTING,
		  // SAME REASON AS ABOVE
		  sys.virSum += force_val * r2 * 0.5;

		  // THE 1/2 TO ACCOUNT FOR DOUBLE COUNTING
		  //sys.pot_energy += 4.*0.5*inv_r6*(inv_r6 - 1.);
		  sys.pot_energy += ( 4.*inv_r6*(inv_r6 - 1.) - ecut ) / 2.;
		  //sys.pot_energy += 4.*0.5*inv_r6*(inv_r6 - 1.) - ecut;
		  //sys.pot_energy += 4.*inv_r6*(inv_r6 - 1.) - ecut;
		}
	  }
	}
}

double forces::calc_LJ_forces_single_particle(int ii, vec &direction, double &virSum, ensemble& sys)
{
  // CALCULATE THE LENNARD JONES FORCES ON A SINGLE GIVEN PARTICLE
  // AND RETURN THE ENERGY CONTRINBUTION DUE TO ONLY THAT PARTICLE.
  // RETURN DIRECTION OF THE FORCE IN direction
  // RETURN CONTRIBUTION TO VIRIAL IN virSum

  // CLEAR ENERGY
  double pot_energy = 0.;
  // CLEAR FORCES
  sys.clear_forces(ii);
  // CLEAR VIRIAL
  virSum = 0.;

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

	  // VIRIAL
	  virSum += force_val * r2;

	  pot_energy += 4.*inv_r6*(inv_r6 - 1.);
	}
  }
  // SAVE THE DIRECTION FOR RETURN 
  direction = sys.atoms[ii].f / sys.atoms[ii].f.norm(); 

  return pot_energy;
}

double forces::calc_LJXY_forces_single_particle(int ii, vec &direction, double &virSum, ensemble& sys)
{
  // CALCULATE THE LENNARD JONES FORCES AND XY MODEL ANGLE ENERGY TERM
  // ON A SINGLE GIVEN PARTICLE
  // AND RETURN THE ENERGY CONTRIBUTION DUE TO ONLY THAT PARTICLE.
  // RETURN DIRECTION OF THE FORCE IN direction
  // RETURN CONTRIBUTION TO VIRIAL IN virSum
  //
  // NOTE: VIRIAL MIGHT BE OFF SINCE WE ARE NOT TAKING INTO ACCOUNT THE TORQUE

  // CLEAR ENERGY
  double pot_energy = 0.;
  // CLEAR FORCES
  sys.clear_forces(ii);
  // CLEAR VIRIAL
  virSum = 0.;

  vec dr;
  double r2, inv_r2, inv_r6;
  double rc2  = sys.rc*sys.rc;
  double force_val;
  double XY_term;

  // VALUE OF INTERACTION ENERGY
  //double XY_J = 1.;
  //double XY_J = 10.;

  // PAIRWISE INTERACTIONS BETWEEN PARTICLE ii
  // AND THE REST
  XY_term = 0.;
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

	  // VIRIAL
	  virSum += force_val * r2;

	  // XY MODEL TERM ONLY ON PARTICLES WITHIN CUTOFF
	  if (r2 < rc2)
	  {
		XY_term += -sys.XY_J * cos(sys.atoms[i].theta - sys.atoms[ii].theta);
	  }


	  pot_energy += 4.*inv_r6*(inv_r6 - 1.);
	}
  }
  pot_energy += XY_term;

  // SAVE THE DIRECTION FOR RETURN 
  direction = sys.atoms[ii].f / sys.atoms[ii].f.norm(); 

  return pot_energy;
}

double forces::calc_LJLL_forces_single_particle(int ii, vec &direction, double &virSum, ensemble& sys)
{
  // CALCULATE THE LENNARD JONES FORCES AND LL MODEL (COSINE SQUARE) ANGLE ENERGY TERM
  // ON A SINGLE GIVEN PARTICLE
  // AND RETURN THE ENERGY CONTRIBUTION DUE TO ONLY THAT PARTICLE.
  // RETURN DIRECTION OF THE FORCE IN direction
  // RETURN CONTRIBUTION TO VIRIAL IN virSum
  //
  // NOTE: VIRIAL MIGHT BE OFF SINCE WE ARE NOT TAKING INTO ACCOUNT THE TORQUE

  // CLEAR ENERGY
  double pot_energy = 0.;
  // CLEAR FORCES
  sys.clear_forces(ii);
  // CLEAR VIRIAL
  virSum = 0.;

  vec dr;
  double r2, inv_r2, inv_r6;
  double rc2  = sys.rc*sys.rc;
  double force_val;
  double LL_term;
  double cos2;

  // VALUE OF INTERACTION ENERGY
  //double XY_J = 1.;
  //double XY_J = 10.;

  // PAIRWISE INTERACTIONS BETWEEN PARTICLE ii
  // AND THE REST
  LL_term = 0.;
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

	  // VIRIAL
	  virSum += force_val * r2;

	  // LL MODEL TERM ONLY ON PARTICLES WITHIN CUTOFF
	  if (r2 < rc2)
	  {
		cos2 = cos(sys.atoms[i].theta - sys.atoms[ii].theta);
		cos2 = CUAD(cos2);
		LL_term += -sys.XY_J * cos2;
	  }


	  pot_energy += 4.*inv_r6*(inv_r6 - 1.);
	}
  }
  pot_energy += LL_term;

  // SAVE THE DIRECTION FOR RETURN 
  direction = sys.atoms[ii].f / sys.atoms[ii].f.norm(); 

  return pot_energy;
}

// ---------- ENERGY LANDSCAPES ------------

void forces::energy_landscape(ensemble& sys, char* outfile)
{
  // CALCULATES ENERGY LANDSCAPE
  int i,j,k;
  //double step = .1;
  double step = .05;
  int num_x   = int( sys.Lx / step );
  int num_y   = int( sys.Ly / step );
  int num_z   = int( sys.Lz / step );
  double xx,yy,zz;
  vec pos;
  double curr_energy;

  ofstream out_f(outfile);

#if NDIM == 3
  PRINT1("No energy landscape in three dimensions!");
#else
  // THE ORDERING OF THE WRITING OF THE ENERGY MATRIX
  // CORRESONDS TO A CORRECT VISUALIZATION BY THE vis
  // PROGRAM.
  //for (j=1; j<num_y; j++)
  for (j=num_y; j>=0; j--)
  {
	for (i=0; i<num_x; i++)
	{
	  xx = i*step - sys.Lx / 2.;
	  yy = j*step - sys.Ly / 2.;
	  pos = vec( xx , yy , 0. );
	  curr_energy = calc_LJ_forces_energy_landscape(pos, sys);

	  //out_f << xx << " " << yy << " " << curr_energy << "\n";
	  out_f << curr_energy << " ";
	}
	out_f << "\n";
  }
#endif

  out_f.close();
}

double forces::calc_LJ_forces_energy_landscape(vec &position, ensemble& sys)
{
  // CALCULATE THE LENNARD JONES TOTAL POTENTIAL ENERGY AT 
  // POSITION position FROM ALL THE PARTICLES.
  // N*(N-1)/2 CALCULATION.

  // CLEAR ENERGY
  sys.pot_energy = 0.;
  // CLEAR FORCES
  sys.clear_forces();
  // CLEAR VIRIAL
  sys.virSum = 0.;

  vec dr;
  double r2, inv_r2, inv_r6;
  double force_val;
  int consider = 1;
  double eps;

  // PAIRWISE INTERACTIONS
  for (int i=0; i<sys.N; i++)
  {
	dr = sys.atoms[i].r - position;

	// USE IMAGE DISTANCES
	dr.x = IMAGE(dr.x, sys.inv_Lx, sys.Lx);
	dr.y = IMAGE(dr.y, sys.inv_Ly, sys.Ly);
	dr.z = IMAGE(dr.z, sys.inv_Lz, sys.Lz);

	// CHECK FOR POSITIONS TOO CLOSE TO THE 
	// HARD CORE
	//eps = 0.0 * sys.atoms[i].d;
	//eps = 0.1 * sys.atoms[i].d;
	eps = 0.5 * sys.atoms[i].d;
	if ( (dr.norm2() >= CUAD(sys.atoms[i].d - eps)) && consider )
	{
	  r2 = dr.norm2();
	  inv_r2 = 1.0 / r2;
	  inv_r6 = inv_r2 * inv_r2 * inv_r2;

	  sys.pot_energy += 4.*inv_r6*(inv_r6 - 1.);
	}
	else
	{
	  // OVERLAP. TIME TO EXIT WITH 0 VALUE.
	  sys.pot_energy = 0.;
	  consider       = 0;
	}
  }

  return sys.pot_energy;
}

/*
 * ----------------------------------------------------------------------
 */
