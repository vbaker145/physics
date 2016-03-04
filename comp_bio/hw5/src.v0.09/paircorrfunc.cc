
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

#include "paircorrfunc.h"

#include "ensemble.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const paircorrfunc &z)
{
  // FRIEND FUNCTION TO WRITE ELEMENTS OF DATA
//os << z.N  << " "
//<< "\n"
//;

  return os;
}


/*
 * ----------------------------------------------------------------------
 *                 PAIR CORRELATION FUNCTION CLASS
 */

paircorrfunc::paircorrfunc(void) { }

paircorrfunc::~paircorrfunc(void) 
{ 
  delete [] num_bins;
  delete [] gr;
}

paircorrfunc::paircorrfunc(ensemble& sys) 
{ 
  // NUMBER OF HISTOGRAM BINS
  nhis = sys.nhis_pcf;

  // TOTAL TIME FOR EACH MEASUREMENT
  duration = sys.time_meas;
  // TOTAL STEPS FOR EACH MEASUREMENT
  i_duration = int(duration / sys.delta_t);
  // CURRENT COUNTER
  curr_time = 0;

  // HOW MANY SAMPLES
  howmany = int(sys.total_num_steps / (1.0*i_duration));

  // - GET SMALLEST DIMENSION OF BOX
  minL = MIN(sys.Lx,sys.Ly);
#if NDIM == 3
#else
  minL = MIN(minL,sys.Lz);
#endif
  // - SIZE OF THE BIN
  delta_gr = minL / (2.0*nhis);

  // SPACE FOR THE ARRAYS 
  num_bins  = new int[nhis];
  gr        = new double[nhis];
}

void paircorrfunc::reset(ensemble& sys)
{
  int i,j;
  
  // INITIALIZE THE LOCAL COUNTER
  curr_time = 0;

  // INITIALIZE ARRAYS
  for (i=0; i<nhis; i++)
  {
  	num_bins[i]  = 0;
  	gr[i]  = 0.;
  }
}

void paircorrfunc::accum(int n, ensemble& sys)
{
  // MAIN FUNCTION THAT EITHER RESETS OR ACCUMULATES
  // MEASUREMENTS
  int i,j, i_gr;
  double r, vol_bin, num_ideal;
  vec dr;

  // CHECK WITHIN MEASUREMENT CYCLE.
  // ACCUMULATION LOOP.
  if (1.0*n/i_duration <= howmany)
  {
	// **********************
	// RESET IF NEEDED.
	if ( n % i_duration == 0 )
	{
	  // IF THE SIMULATION COUNTER COINCIDES WITH THE 
	  // DURATION OF THE MEASUREMENT, WE ARE STARTING A
	  // NEW MEAS. CYCLE 
	  if (n!=0)
	  {
		// APPEND RESULTS FROM THE PREVIOUS CYCLE.
		ofstream out_gr("results/gr.dat",ios_base::app);

		for (j=1; j<nhis; j++)
		{
		  r = delta_gr * (j+.5);
#if NDIM == 3
		  vol_bin = ( pow(j+1,3)-pow(j,3) ) * delta_gr*delta_gr*delta_gr;
		  num_ideal = (4./3.) * PI * vol_bin * sys.density;
#else
		  vol_bin = ( pow(j+1,2)-pow(j,2) ) * delta_gr*delta_gr;
		  num_ideal = PI * vol_bin * sys.density;
#endif
		  //out_gr  << j*delta_gr << " " << gr[j]/( i_duration * sys.N * num_ideal )  << "\n";
		  // NORMALIZE BY THE NUMBER OF BINS PER DISTANCE. USED WHEN
		  // PARTICLES ARE FIXED IN SPACE
		  // -- NOT QUITE WORKING FOR FIXED PARTICLES!
		  //out_gr  << j*delta_gr << " " << gr[j]/( i_duration * num_bins[j] * num_ideal )  << "\n";
		  //out_gr  << j*delta_gr << " " << gr[j]/( i_duration *  num_ideal )  << "\n";
		  out_gr  << j*delta_gr << " " << gr[j]/( i_duration * sys.num_not_fixed * num_ideal )  << "\n";
		}
		out_gr  << "\n";

		out_gr.close();
	  }

	  // RESET VALUES
	  this->reset(sys);	
	}
	// END RESET
	// **********************

	// **********************
	// START ACCUMULATE RESULTS.
	for (int i=0; i<sys.N-1; i++)
	  for (int j=i+1; j<sys.N; j++)
	  {
		// ONLY MEASURE IF BOTH PARTICLES ARE NOT FIXED.
		// -- NOT SURE THIS IS THE WAY TO GO!
		if ( sys.atoms[i].is_not_fixed() && sys.atoms[j].is_not_fixed() )
		{
		  dr = sys.atoms[i].r - sys.atoms[j].r;
		  // USE IMAGE DISTANCES
		  dr.x = IMAGE(dr.x, sys.inv_Lx, sys.Lx);
		  dr.y = IMAGE(dr.y, sys.inv_Ly, sys.Ly);
		  dr.z = IMAGE(dr.z, sys.inv_Lz, sys.Lz);

		  // DISTANCE
		  r = dr.norm();

		  // ACCUMULATE IF LESS THAN HALF OF THE
		  // BOX SIZE
		  if ( r < (minL/2.0) )
		  {
			i_gr     = int( r / delta_gr );
			gr[i_gr] = gr[i_gr] + 2.;
			// THIS IS FOR THE FIXED PARTICLES, BUT NOT
			// SURE THIS IS WORKING...
			num_bins[i_gr]++;
		  }
		}

	  }
	// END ACCUMULATE RESULTS.
	// **********************

	curr_time++;
  }
}

/*
 * ----------------------------------------------------------------------
 */
