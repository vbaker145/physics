
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
	// CHECK TO SEE IF ATOM IS FIXED
	if (sys.atoms[i].is_not_fixed())
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
  }

  // UPDATE INSTANTANEOUS TEMPERATURE
  sys.kin_energy  = sum_vel2 / 2.;
#if NDIM == 2
  //sys.temperature = sum_vel2 / (2.*sys.N-2.);
  sys.temperature = sum_vel2 / (2.*sys.num_not_fixed-2.);
#else
  //sys.temperature = sum_vel2 / (3.*sys.N-3.);
  sys.temperature = sum_vel2 / (3.*sys.num_not_fixed-3.);
#endif
  // UPDATE TOTAL ENERGY
  //sys.energy = (sys.pot_energy + 0.5 * sum_vel2) / sys.N;
  //sys.energy = (sys.pot_energy + sys.kin_energy) / sys.N;
  sys.energy = sys.pot_energy + sys.kin_energy;
}

void integrator::verlet_polymers(double &virSum, ensemble& sys)
{
  // SHAKE ALGORITHM FOR CONSTRAINT DYNAMICS OF A CHAIN MOLECULE
  // FROM Allen and Tildesley, F.08 FORTRAN CODE
  int i,j,k;
  int atomi;
  int it, isdone, a, b;
  vec sum_vel;
  int moving[sys.atoms_per_chain];
  int moved[sys.atoms_per_chain];
  double sum_vel2 = 0.;
  double length_bonds2 = sys.length_bonds * sys.length_bonds;
  double PABSQ, RABSQ, DIFFSQ, RPAB;
  double RMA, RMB, GAB;

  // OPEN CHAINS (i.e. NOT CYCLIC)
  int num_bonds = sys.atoms_per_chain - 1;

  // ADJUSTABLE PARAMETERS
  int maxit = 1000;
  double tolerance = 1.e-7;
  double RPTOL     = 1.e-7;

  vec *new_r = new vec[sys.atoms_per_chain];
  vec new_pr;
  vec new_rr;
  vec DD;

  // LOOP OVER CHAINS
  for (j=0; j<sys.num_chains; j++)
  {
	atomi = j * sys.atoms_per_chain;

	// LOOP OVER ATOMS IN THE CHAIN
	for (i=0; i<sys.atoms_per_chain; i++)
	{
		// GET NEW UNBOUNDED POSITIONS
		new_r[i]  = 2.*sys.atoms[i+atomi].r - sys.atoms[i+atomi].r_old + sys.delta_t*sys.delta_t * sys.atoms[i+atomi].f / sys.atoms[i+atomi].m ;
		//cout << "pos " << sys.atoms[i+atomi].r << " +--+ " << sys.atoms[i+atomi].r_old << "\n";
		moving[i] = 0;
		moved[i]  = 1;
	}

	it     = 0;
	isdone = 0;

	// BEGIN ITERATIVE LOOP FOR THE SHAKE ALGORITHM
	while ( (!isdone) && (it <= maxit) )
	{
	  isdone = 1;

	  for (a=0; a<num_bonds; a++)
	  {
		b = a + 1;
		if (b>sys.atoms_per_chain) b = 1;

		if ( moved[a] || moved[b] )
		{
		  new_pr   = new_r[a] - new_r[b];
		  new_pr.x = IMAGE(new_pr.x, sys.inv_Lx, sys.Lx);
		  new_pr.y = IMAGE(new_pr.y, sys.inv_Ly, sys.Ly);
		  new_pr.z = IMAGE(new_pr.z, sys.inv_Lz, sys.Lz);

		  PABSQ = new_pr.norm2();
		  RABSQ = length_bonds2;
		  DIFFSQ = RABSQ - PABSQ;
		  //cout << "   DIFFSQ between " << a << " and " << b << " = " << DIFFSQ << "\n";
		  //cout << "      PABSQ " << PABSQ << "\n";
		  //cout << "      RABSQ " << RABSQ << "\n";
		  //cout << "      crit " << RABSQ * 2 * tolerance << "\n";

		  if (ABS(DIFFSQ) > (RABSQ * 2 * tolerance) )
		  {
			new_rr = sys.atoms[a+atomi].r - sys.atoms[b+atomi].r;
			new_rr.x = IMAGE(new_rr.x, sys.inv_Lx, sys.Lx);
			new_rr.y = IMAGE(new_rr.y, sys.inv_Ly, sys.Ly);
			new_rr.z = IMAGE(new_rr.z, sys.inv_Lz, sys.Lz);

			RPAB = new_rr.x*new_pr.x + new_rr.y*new_pr.y + new_rr.z*new_pr.z;

			if ( RPAB < (RABSQ * RPTOL) )
			{
			  cout << "CONSTRAINT FAILURE\n";
			  cout << "RPAB = " << RPAB << " RABSQ * RPTOL = " << RABSQ * RPTOL << "\n";
			  exit(1);
			}

			RMA = 1. / sys.atoms[a+atomi].m;
			RMB = 1. / sys.atoms[b+atomi].m;
			GAB = DIFFSQ / ( 2.0 * ( RMA + RMB ) * RPAB );
			sys.virSum += GAB * RABSQ;
			DD = new_rr * GAB;
			//cout << " ->DD " << DD << "\n";

			new_r[a] += RMA * DD;
			new_r[b] -= RMB * DD;

			moving[a] = 1;
			moving[b] = 1;
			isdone = 0;

			//cout << "mol  " << j << "\n";
			//cout << "atoms " << a << " " << b << "\n";
			//cout << "atomi " << atomi << "\n";
			//cout << "posa " << sys.atoms[a+atomi].r << "\n";
			//cout << "posb " << sys.atoms[b+atomi].r << "\n";
		  }
		}
	  }

	  for (a=0; a<num_bonds; a++)
	  {
		moved[a]  = moving[a];
		moving[a] = 0;
	  }
	  ++it;
	  if (it > maxit) cout << "MAX NUM ITERATIONS EXCEEDED IN SHAKE\n";
	  //cout << " 00 NUM ITERATIONS " << it << "\n";
	}
	//cout << " 01 TOTAL NUM ITERATIONS MOL " << j << " = " << it << "\n";
	// END ITERATIVE LOOP

	if (!isdone)
	{
	  cout << "TOO MANY CONSTRAINT ITERATIONS\n";
	  cout << "MOLECULE " << j << "\n";
	  exit(1);
	}

	//
	for (i=0; i<sys.atoms_per_chain; i++)
	{
		// GET VELOCITY
		sys.atoms[i+atomi].v = (new_r[i] - sys.atoms[i+atomi].r_old) / (2.*sys.delta_t);

		// VELOCITY OF CENTER OF MASS
		sum_vel   = sum_vel + sys.atoms[i+atomi].v;
		// TOTAL KINETIC
		sum_vel2 += sys.atoms[i+atomi].v.norm2();

		// UPDATE
		//cout << " ============== \n";
		//cout << " UPDATING POSITIONS \n";
		//cout << "   OLD = " << sys.atoms[i+atomi].r_old << " : " << sys.atoms[i+atomi].r << "\n";
		sys.atoms[i+atomi].r_old = sys.atoms[i+atomi].r;
		//cout << "   NEW = " << new_r[i] << "\n";
		sys.atoms[i+atomi].r     = new_r[i];
		//cout << " ============== \n";
	}
  }

  // UPDATE INSTANTANEOUS TEMPERATURE
  sys.kin_energy  = sum_vel2 / 2.;
#if NDIM == 2
  //sys.temperature = sum_vel2 / (2.*sys.N-2.);
  sys.temperature = sum_vel2 / (2.*sys.num_not_fixed-2.);
#else
  //sys.temperature = sum_vel2 / (3.*sys.N-3.);
  sys.temperature = sum_vel2 / (3.*sys.num_not_fixed-3.);
#endif
  // UPDATE TOTAL ENERGY
  //sys.energy = (sys.pot_energy + 0.5 * sum_vel2) / sys.N;
  //sys.energy = (sys.pot_energy + sys.kin_energy) / sys.N;
  sys.energy = sys.pot_energy + sys.kin_energy;

  // VIRIAL
  sys.virSum = sys.virSum / (sys.delta_t*sys.delta_t) / 3.0;

  delete [] new_r;
}

void integrator::velocity_verlet_pt1(forces &physics, ensemble& sys)
{
  // VELOCITY VERLET PART 1
  // HAVE TO CALCULATE FORCES BEFORE GOING TO PART 2
  //  ======================= 1 ======================

  int i;
  vec sum_vel;
  double sum_vel2 = 0.;

  vec new_r;
  vec middle_v;

  for (i=0; i<sys.N; i++)
  {
	// CHECK TO SEE IF ATOM IS FIXED
	if (sys.atoms[i].is_not_fixed())
	{
	  // GET NEW POSITIONS
	  // NOTE THAT r_old ARE NOT USED HERE
	  sys.atoms[i].r = sys.atoms[i].r + sys.delta_t*sys.atoms[i].v + sys.delta_t*sys.delta_t*sys.atoms[i].f/(2.0*sys.atoms[i].m) ;

	  // GET HALF-STEP VELOCITY
	  sys.atoms[i].v_new = sys.atoms[i].v + 0.5 * sys.delta_t *sys.atoms[i].f/sys.atoms[i].m;
	}
  }
}

void integrator::velocity_verlet_pt2(forces &physics, ensemble& sys)
{
  // VELOCITY VERLET PART 2
  //  ======================= 2 ======================

  int i;
  vec sum_vel;
  double sum_vel2 = 0.;

  vec new_r;
  vec middle_v;

  for (i=0; i<sys.N; i++)
  {
	// CHECK TO SEE IF ATOM IS FIXED
	if (sys.atoms[i].is_not_fixed())
	{
	  // GET WHOLE-STEP VELOCITY USING THE NEW FORCES
	  sys.atoms[i].v = sys.atoms[i].v_new + 0.5 * sys.delta_t *sys.atoms[i].f/sys.atoms[i].m;

	  // VELOCITY OF CENTER OF MASS
	  sum_vel   = sum_vel + sys.atoms[i].v;
	  // TOTAL KINETIC
	  sum_vel2 += sys.atoms[i].v.norm2();
	}
  }

  // UPDATE INSTANTANEOUS TEMPERATURE
  sys.kin_energy  = sum_vel2 / 2.;
#if NDIM == 2
  //sys.temperature = sum_vel2 / (2.*sys.N-2.);
  sys.temperature = sum_vel2 / (2.*sys.num_not_fixed-2.);
#else
  //sys.temperature = sum_vel2 / (3.*sys.N-3.);
  sys.temperature = sum_vel2 / (3.*sys.num_not_fixed-3.);
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
	// RESCALE VELOCITIES.
	// FIXED ATOMS SHOULD NOT HAVE A VELOCITY AT THIS POINT
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
