
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

#include "diffusion.h"

#include "ensemble.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const diffusion &z)
{
  // FRIEND FUNCTION TO WRITE ELEMENTS OF DATA
//os << z.N  << " "
//<< "\n"
//;

  return os;
}


/*
 * ----------------------------------------------------------------------
 *                 DIFFUSION CLASS
 */

diffusion::diffusion(void) { }

diffusion::~diffusion(void) 
{ 
  delete [] r2t;
  delete [] vacf;
  delete [] rt0;
  delete [] vt0;
}

diffusion::diffusion(ensemble& sys) 
{ 
  // TOTAL TIME FOR EACH MEASUREMENT
  duration = sys.time_meas;
  // TOTAL STEPS FOR EACH MEASUREMENT
  i_duration = int(duration / sys.delta_t);
  // CURRENT COUNTER
  curr_time = 0;

  // HOW MANY SAMPLES
  howmany = int(sys.total_num_steps / (1.0*i_duration));

  // SPACE FOR THE ARRAYS 
  r2t  = new double[i_duration];
  vacf = new double[i_duration];

  rt0 = new vec[sys.N];
  vt0 = new vec[sys.N];

  //PRINT3(duration,i_duration,howmany);
}

void diffusion::reset(ensemble& sys)
{
  int i,j;
  
  // INITIALIZE THE LOCAL COUNTER
  curr_time = 0;

  // INITIALIZE ARRAYS
  for (i=0; i<i_duration; i++)
  {
  	r2t[i]  = 0.;
  	vacf[i] = 0.;
  }

  // SAVE INITIAL VALUES OF POSITION AND VELOCITY
  for (i=0; i<sys.N; i++)
  {
	rt0[i] = sys.atoms[i].r;
	vt0[i] = sys.atoms[i].v;
  }
}

void diffusion::accum(int n, ensemble& sys)
{
  // MAIN FUNCTION THAT EITHER RESETS OR ACCUMULATES
  // MEASUREMENTS
  int i,j;
  vec dist;

  // CHECK WITHIN MEASUREMENT CYCLE.
  // ACCUMULATION LOOP.
  if (1.0*n/i_duration <= howmany)
  {
	// RESET IF NEEDED.
	if ( n % i_duration == 0 )
	{
	  // IF THE SIMULATION COUNTER COINCIDES WITH THE 
	  // DURATION OF THE MEASUREMENT, WE ARE STARTING A
	  // NEW MEAS. CYCLE 

	  if (n!=0)
	  {
		// WRITE RESULTS FROM THE PREVIOUS CYCLE.
		//PRINT1("*** PRINT RESULTS FOR R2T AND VACF");
		ofstream out_r2t("results/r2t.dat",ios_base::app);
		ofstream out_vacf("results/vacf.dat",ios_base::app);

		for (j=0; j<i_duration; j++)
		{
		  // TAKE CARE OF FIXED ATOMS
		  //out_r2t  << j*sys.delta_t << " " << r2t[j]/sys.N  << "\n";
		  out_r2t  << j*sys.delta_t << " " << r2t[j]/sys.num_not_fixed  << "\n";
		  //out_vacf << j*sys.delta_t << " " << vacf[j]/sys.N << "\n";
		  out_vacf << j*sys.delta_t << " " << vacf[j]/sys.num_not_fixed << "\n";
		}
	  	out_r2t  << "\n";
	  	out_vacf << "\n";

		out_r2t.close();
		out_vacf.close();
	  }

	  // RESET VALUES
	  this->reset(sys);	
	}

	// ACCUMULATE RESULTS.
	// AVERAGE OVER ALL PARTICLES
	// FOR THIS PARTICULAR TIME.
	for (i=0; i<sys.N; i++)
	{
	  // ONLY COUNT MOBILE ATOMS
	  if (sys.atoms[i].is_not_fixed())
	  {
		dist = rt0[i] - sys.atoms[i].r;
		r2t[curr_time]  = r2t[curr_time]  + ScalarProd(dist,dist);
		vacf[curr_time] = vacf[curr_time] + ScalarProd(vt0[i],sys.atoms[i].v) ;
	  }
	}

	curr_time++;
  }
}

/*
 * ----------------------------------------------------------------------
 */
