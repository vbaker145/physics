
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

#include "stepwise.h"

#include "ensemble.h"
#include "verletlist.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const stepwise &z)
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

stepwise::stepwise(void) { }

stepwise::~stepwise(void) { }

void stepwise::calc_hardsphere_stepwise(ensemble& sys)
{
  // CALCULATE COLLISION TIMES AND MOVE
  // FOR HARD SPHERE SYSTEM.
  // DOES NOT USE ANY KIND OF LIST.
  // N*(N-1)/2 CALCULATION PER COLLISION

  // CLEAR ENERGY
  sys.pot_energy = 0.;
  // CLEAR VIRIAL
  sys.virSum = 0.;

  int i_col, j_col;
  vec dr;
  vec rij;
  vec dv;
  double r2, v2, bij, bij2, c2;
  double deter, tij1, tij2;
  double t_test_col, r_min;
  double fac_delta_v;

  // ========== 1 ===============================
  // LOOP OVER PAIRS TO FIND THE NEXT
  // COLLISION TIME
  double t_col = sys.total_time;
  for (int i=0; i<sys.N-1; i++)
	for (int j=i+1; j<sys.N; j++)
	{
	  // EXCLUDE COLLITION TIME CALCULATION IF BOTH PARTICLES
	  // ARE FIXED IN SPACE
	  if ( sys.atoms[i].is_not_fixed() || sys.atoms[j].is_not_fixed() )
	  {
		dr = sys.atoms[i].r - sys.atoms[j].r;
		// USE IMAGE DISTANCES
		dr.x = IMAGE(dr.x, sys.inv_Lx, sys.Lx);
		dr.y = IMAGE(dr.y, sys.inv_Ly, sys.Ly);
		dr.z = IMAGE(dr.z, sys.inv_Lz, sys.Lz);
		r2 = dr.norm2();

		// VELOCITY DIFFERENCE
		dv = sys.atoms[i].v - sys.atoms[j].v;
		v2 = dv.norm2();

		// COLLISION DETERMINANT
		bij  = ScalarProd(dr,dv);
		bij2 = bij * bij;
		c2   = v2*(r2-1.);

		// TEST FOR POSITIVE ROOTS
		if ( (bij<0.) && (bij2>=c2) )
		{
		  deter = sqrt( bij2 - c2 );
		  tij1  = (-bij - deter)/v2;
		  tij2  = (-bij + deter)/v2;
		  // COLLISION TIME IS THE SMALLEST ROOT
		  if (tij1>0.&&tij2<0.)
			t_test_col = tij1;
		  else if (tij2>0.&&tij1<0.)
			t_test_col = tij2;
		  else
			t_test_col = MIN(tij1, tij2);
		}
		else
		{
		  t_test_col = sys.total_time;
		}

		// TEST WHETHER THIS COLLISION TIME IS
		// SMALLER THAN THE ONE IN RECORD
		if (t_test_col<t_col)
		{
		  i_col = i;
		  j_col = j;
		  t_col = t_test_col;
		}
	  }
	}

  // THIS IS NOT RIGHT!
  //sys.virSum += force_val * r2;
  sys.virSum     += 0.;
  sys.pot_energy += 0.;

  // ========== 2 ===============================
  // CALCULATE VELOCITY FACTOR AND RECALCULATE
  // THE B FACTOR FOR THE VELOCITY UPDATES
  vec new_ri = sys.atoms[i_col].r + sys.atoms[i_col].v * t_col;
  vec new_rj = sys.atoms[j_col].r + sys.atoms[j_col].v * t_col;
  rij = new_ri - new_rj;
  rij.x = IMAGE(rij.x, sys.inv_Lx, sys.Lx);
  rij.y = IMAGE(rij.y, sys.inv_Ly, sys.Ly);
  rij.z = IMAGE(rij.z, sys.inv_Lz, sys.Lz);
  //
  vec new_v  = sys.atoms[i_col].v - sys.atoms[j_col].v;
  //
  bij  = ScalarProd(rij,new_v);
  fac_delta_v = -bij/CUAD(1.);
  vec delta_v = fac_delta_v * rij;
  //cout << "FAC DELTA V " << delta_v << " * " << fac_delta_v << " * " << rij << "\n";

  // ========== 3 ===============================
  // ADVANCE EVERY PARTICLE BALLISTICALLY.
  // PERFORM COLLISION ON THE COLLIDING PAIR.
  double sum_vel2 = 0.;
  vec new_r;
  for (int i=0; i<sys.N; i++)
  {
	if ( sys.atoms[i].is_not_fixed() )
	{
	  // GET NEW POSITION
	  new_r = sys.atoms[i].r + sys.atoms[i].v * t_col;
	  //cout << i << " position " << sys.atoms[i].r << "\n";
	  //cout << "     velocity " << sys.atoms[i].v << "\n";

	  // VELOCITY REMAINS THE SAME EXCEPT
	  // FOR COLLIDING PARTICLES
	  if (i==i_col) 
		sys.atoms[i].v = sys.atoms[i].v + delta_v;
	  else if (i==j_col)
		sys.atoms[i].v = sys.atoms[i].v - delta_v;

	  // TOTAL KINETIC
	  sum_vel2 += sys.atoms[i].v.norm2();

	  // UPDATE
	  sys.atoms[i].r_old = sys.atoms[i].r;
	  sys.atoms[i].r     = new_r;
	}
  }

  // ========== 4 ===============================
  // UPDATE RUNNING CLOCK
  sys.simulation_time += t_col;
  cout << sys.simulation_time << " " << t_col << "\n";

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




/*
 * ----------------------------------------------------------------------
 */
