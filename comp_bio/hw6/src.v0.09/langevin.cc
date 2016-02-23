
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

#include "langevin.h"

#include "ensemble.h"
#include "forces.h"
#include "ran01.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const langevin &z)
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

langevin::langevin(void)
{
}


langevin::~langevin(void) { }

void langevin::integrate_ermak(ensemble& sys)
{
  // INTEGRATION ALGORITHM FROM:
  // Ermak, D.L. and H. Buckholtz, J. Comp. Phys. 35 (1980) 169.
  // REFERED FROM THE HINCHLIFFE BOOK, P. 374 

  vec sum_vel;
  double sum_vel2 = 0.;
  double gamma = sys.damping_constant;

  vec new_r;
  vec new_v;

  // c PARAMETERS
  double gt = gamma * sys.delta_t; 
  double c0 = exp( -gt );
  double c1 = (1. - c0) / gt;
  double c2 = (1. - c1) / gt;

  // STDEVS FOR THE GAUSSIANS
  double rt = sys.requested_temperature;
  double dt = sys.delta_t;
  double sigma2_r = CUAD(dt) * rt * (2. - (3.-4*exp(-gt)+exp(-2*gt))/gt )/gt; 
  double sigma2_v = rt * (1. - exp(-gt) );
  double sigma_r = sqrt(sigma2_r);
  double sigma_v = sqrt(sigma2_v);

  for (int i=0; i<sys.N; i++)
  {
	// CHECK TO SEE IF ATOM IS FIXED
	if (sys.atoms[i].is_not_fixed())
	{
	  // GET NEW POSITIONS
	  new_r = sys.atoms[i].r 
		+ (c1 * dt) * sys.atoms[i].v 
		+ c2 * sys.atoms[i].f * CUAD(dt) / sys.atoms[i].m 
#if NDIM == 2
		+ vec( ran_gaussian_02_mu_sigma(0.,sigma_r), ran_gaussian_02_mu_sigma(0.,sigma_r), 0. );
#else
		+ vec( ran_gaussian_02_mu_sigma(0.,sigma_r), ran_gaussian_02_mu_sigma(0.,sigma_r), ran_gaussian_02_mu_sigma(0.,sigma_r) );
#endif

	  // GET VELOCITY
	  new_v = c0 * sys.atoms[i].v 
		+ (c1 * dt) * sys.atoms[i].f  / sys.atoms[i].m 
#if NDIM == 2
		+ vec( ran_gaussian_02_mu_sigma(0.,sigma_v), ran_gaussian_02_mu_sigma(0.,sigma_v), 0. );
#else
		+ vec( ran_gaussian_02_mu_sigma(0.,sigma_v), ran_gaussian_02_mu_sigma(0.,sigma_v), ran_gaussian_02_mu_sigma(0.,sigma_v) );
#endif

	  // UPDATE V
	  sys.atoms[i].v = new_v;

	  // VELOCITY OF CENTER OF MASS
	  sum_vel   = sum_vel + sys.atoms[i].v;
	  // TOTAL KINETIC
	  sum_vel2 += sys.atoms[i].v.norm2();

	  // UPDATE R
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

void langevin::integrate_other_pt1(ensemble& sys)
{
  // INTEGRATION ALGORITHM USING TWO STEPS, FOR DISTANCE
  // UPDATES AND FOR VELOCITY UPDATES
  //     * PART 1 *

  double gamma = sys.damping_constant;
  vec new_r;

  // c PARAMETERS
  double gt = gamma * sys.delta_t; 
  double c0 = exp( -gt );
  double c1 = (1. - c0) / gt;
  double c2 = (1. - c1) / gt;

  // STDEVS FOR THE GAUSSIANS
  double rt = sys.requested_temperature;
  double dt = sys.delta_t;
  double sigma2_r = CUAD(dt) * rt * (2. - (3.-4*exp(-gt)+exp(-2*gt))/gt )/gt; 
  double sigma_r = sqrt(sigma2_r);

  for (int i=0; i<sys.N; i++)
  {
	// CHECK TO SEE IF ATOM IS FIXED
	if (sys.atoms[i].is_not_fixed())
	{
	  // GET NEW POSITIONS
	  new_r = sys.atoms[i].r 
		+ (c1 * dt) * sys.atoms[i].v 
		+ c2 * sys.atoms[i].f * CUAD(dt) / sys.atoms[i].m 
#if NDIM == 2
		+ vec( ran_gaussian_02_mu_sigma(0.,sigma_r), ran_gaussian_02_mu_sigma(0.,sigma_r), 0. );
#else
		+ vec( ran_gaussian_02_mu_sigma(0.,sigma_r), ran_gaussian_02_mu_sigma(0.,sigma_r), ran_gaussian_02_mu_sigma(0.,sigma_r) );
#endif

	  // UPDATE R
	  sys.atoms[i].r_old = sys.atoms[i].r;
	  sys.atoms[i].r     = new_r;
	}
  }
}

void langevin::integrate_other_pt2(ensemble& sys)
{
  // INTEGRATION ALGORITHM USING TWO STEPS, FOR DISTANCE
  // UPDATES AND FOR VELOCITY UPDATES
  //     * PART 2 *

  double gamma = sys.damping_constant;
  vec sum_vel;
  double sum_vel2 = 0.;
  vec new_v;

  // c PARAMETERS
  double gt = gamma * sys.delta_t; 
  double c0 = exp( -gt );
  double c1 = (1. - c0) / gt;
  double c2 = (1. - c1) / gt;

  // STDEVS FOR THE GAUSSIANS
  double rt = sys.requested_temperature;
  double dt = sys.delta_t;
  double sigma2_v = rt * (1. - exp(-gt) );
  double sigma_v = sqrt(sigma2_v);

  for (int i=0; i<sys.N; i++)
  {
	// CHECK TO SEE IF ATOM IS FIXED
	if (sys.atoms[i].is_not_fixed())
	{
	  // GET VELOCITY
	  new_v = c0 * sys.atoms[i].v 
		+ ( (c1-c2) * dt) * sys.atoms[i].f_old / sys.atoms[i].m 
		+ (  c2     * dt) * sys.atoms[i].f     / sys.atoms[i].m 
#if NDIM == 2
		+ vec( ran_gaussian_02_mu_sigma(0.,sigma_v), ran_gaussian_02_mu_sigma(0.,sigma_v), 0. );
#else
		+ vec( ran_gaussian_02_mu_sigma(0.,sigma_v), ran_gaussian_02_mu_sigma(0.,sigma_v), ran_gaussian_02_mu_sigma(0.,sigma_v) );
#endif

	  // UPDATE V
	  sys.atoms[i].v = new_v;

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



/*
 * ----------------------------------------------------------------------
 */
