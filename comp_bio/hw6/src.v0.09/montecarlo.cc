
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

#include "montecarlo.h"

#include "ensemble.h"
#include "forces.h"
#include "ran01.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const montecarlo &z)
{
  // FRIEND FUNCTION TO WRITE ELEMENTS OF DATA
//os << z.N  << " "
//<< "\n"
//;

  return os;
}


/*
 * ----------------------------------------------------------------------
 *                 MONTE CARLO CLASS
 */

montecarlo::montecarlo(void)
{
}


montecarlo::~montecarlo(void) { }

double montecarlo::mc_step(forces &physics, ensemble& sys)
{
  // ========================================
  // FUNCTION TO PERFORM AN MC STEP

  vec new_r;
  vec old_r;
  vec rand_r;
  int accepted_steps = 0;
  double eo, en;
  double virSum, virSumTemp;
  vec direction;
  double beta = 1./sys.temperature;
  double delta_e, delta_be;

  // MC STEP
  for (int i=0; i<sys.N; i++)
  {
	// CHECK TO SEE IF ATOM IS FIXED
	if (sys.atoms[i].is_not_fixed())
	{
	  // ========================================
	  // CALCULATE CURRENT ENERGY OF PARTICLE i
	  eo = physics.calc_LJ_forces_single_particle(i, direction, virSumTemp, sys);
	  // SAVE VIRIAL AND OLD POSITION
	  virSum = virSumTemp;
	  old_r  = sys.atoms[i].r;

	  // ========================================
	  // GENERATE RANDOM VECTOR
#if NDIM == 2
	  rand_r = vec( (rand01()-0.5) * sys.mc_dr, (rand01()-0.5) * sys.mc_dr, 0. );
#else
	  rand_r = vec( (rand01()-0.5) * sys.mc_dr, (rand01()-0.5) * sys.mc_dr, (rand01()-0.5) * sys.mc_dr );
#endif
	  // GENERATE NEW TRIAL POSITION
	  sys.atoms[i].r += rand_r;

	  // ========================================
	  // CALCULATE TRIAL ENERGY OF PARTICLE i
	  en = physics.calc_LJ_forces_single_particle(i, direction, virSumTemp, sys);

	  // ========================================
	  // METROPOLIS ACCEPTANCE CRITERIA
	  delta_e  = en - eo;
	  delta_be = beta * delta_e;
	  if (delta_be < 75.)
	  {
		if (delta_e < 0.)
		{
		  // ACCEPTED.
		  // CHANGE ENERGY AND VIRIAL BY THEIR DIFFERENCES
		  sys.pot_energy += delta_e;
		  sys.virSum      = sys.virSum + virSumTemp - virSum;
		  accepted_steps++;
		}
		else if ( exp(-delta_be) > rand01() )
		{
		  // ACCEPTED.
		  // CHANGE ENERGY AND VIRIAL BY THEIR DIFFERENCES
		  sys.pot_energy += delta_e;
		  sys.virSum      = sys.virSum + virSumTemp - virSum;
		  accepted_steps++;
		}
		else
		{
		  // REJECTED
		  sys.atoms[i].r = old_r;
		}
	  }
	  else
	  {
		// REJECTED
		sys.atoms[i].r = old_r;
	  }

	}
  }

  // MAKE SURE THERE IS NOT KINETIC ENERGY
  sys.kin_energy  = 0.;

  // UPDATE TOTAL ENERGY
  //sys.energy = sys.pot_energy + sys.kin_energy;
  sys.energy = sys.pot_energy;

  return (accepted_steps * 1.0 / sys.N);
}

double montecarlo::mc_step_2D_angle(short which_pot, forces &physics, ensemble& sys)
{
  int accepted_steps = 0;
#if NDIM==3
  PRINT1("NO MC WITH ANGLE IN 3D YET!");
#else
  // ========================================
  // FUNCTION TO PERFORM AN MC STEP,
  // PARTICLES HAVE ANGLES
  // 'which_pot' = 0 : XY MODEL
  // 'which_pot' = 1 : LL MODEL

  vec new_r;
  vec old_r;
  vec rand_r;
  double old_angle;
  double rand_o;
  double eo, en;
  double virSum, virSumTemp;
  vec direction;
  double beta = 1./sys.temperature;
  double delta_e, delta_be;


  // MC STEP
  for (int i=0; i<sys.N; i++)
  {
	// CHECK TO SEE IF ATOM IS FIXED
	if (sys.atoms[i].is_not_fixed())
	{
	  // ========================================
	  // CALCULATE CURRENT ENERGY OF PARTICLE i
	  if (which_pot == 0)
	  	eo = physics.calc_LJXY_forces_single_particle(i, direction, virSumTemp, sys);
	  else if (which_pot == 1)
	  	eo = physics.calc_LJLL_forces_single_particle(i, direction, virSumTemp, sys);
	  // SAVE VIRIAL AND OLD POSITION
	  virSum    = virSumTemp;
	  old_r     = sys.atoms[i].r;
	  old_angle = sys.atoms[i].theta;

	  // ========================================
	  // GENERATE RANDOM POSITION VECTOR AND ORIENTATION
	  rand_r = vec( (rand01()-0.5) * sys.mc_dr, (rand01()-0.5) * sys.mc_dr, 0. );
	  rand_o = sys.theta_o * (rand01()-0.5);

	  // GENERATE NEW TRIAL POSITION
	  sys.atoms[i].r     += rand_r;
	  sys.atoms[i].theta += rand_o;

	  // ========================================
	  // CALCULATE TRIAL ENERGY OF PARTICLE i
	  if (which_pot == 0)
	  	en = physics.calc_LJXY_forces_single_particle(i, direction, virSumTemp, sys);
	  else if (which_pot == 1)
	  	en = physics.calc_LJLL_forces_single_particle(i, direction, virSumTemp, sys);

	  // ========================================
	  // METROPOLIS ACCEPTANCE CRITERIA
	  delta_e  = en - eo;
	  delta_be = beta * delta_e;
	  if (delta_be < 75.)
	  {
		if (delta_e < 0.)
		{
		  // ACCEPTED.
		  // CHANGE ENERGY AND VIRIAL BY THEIR DIFFERENCES
		  sys.pot_energy += delta_e;
		  sys.virSum      = sys.virSum + virSumTemp - virSum;
		  accepted_steps++;
		}
		else if ( exp(-delta_be) > rand01() )
		{
		  // ACCEPTED.
		  // CHANGE ENERGY AND VIRIAL BY THEIR DIFFERENCES
		  sys.pot_energy += delta_e;
		  sys.virSum      = sys.virSum + virSumTemp - virSum;
		  accepted_steps++;
		}
		else
		{
		  // REJECTED
		  sys.atoms[i].r     = old_r;
		  sys.atoms[i].theta = old_angle;
		}
	  }
	  else
	  {
		// REJECTED
		sys.atoms[i].r     = old_r;
		sys.atoms[i].theta = old_angle;
	  }

	}
  }

  // MAKE SURE THERE IS NOT KINETIC ENERGY
  sys.kin_energy  = 0.;

  // UPDATE TOTAL ENERGY
  //sys.energy = sys.pot_energy + sys.kin_energy;
  sys.energy = sys.pot_energy;
#endif
  return (accepted_steps * 1.0 / sys.N);
}

double montecarlo::mc_step_NPT(forces &physics, ensemble& sys)
{
  // ========================================
  // FUNCTION TO PERFORM AN MC STEP
  // IN THE NPT ENSEMBLE

  vec new_r;
  vec old_r;
  vec rand_r;
  int accepted_steps = 0;
  double eo, en;
  double vo, vn;
  double virSum, virSumTemp;
  vec direction;
  double beta = 1./sys.temperature;
  double delta_e, delta_be;
  int i,ii;
  double lnvn;
  double rat_xy, rat_xz;
  double lx, ly, lz;
  double lxo, lyo, lzo;
  // MAX SIZE OF VOLUME CHANGE
  double delta_v = 0.1;

  // ========================================
  // MC STEP,
  // DO N+1 LOOPS TO ALLOW FOR VOLUME CHANGE
  for (int ii=0; ii<=sys.N; ii++)
  {
	// ========================================
	// DECIDE TO DO PARTICLE OR VOLUME TRIAL
	//
	// N+1 POSSIBILITIES COUNTING THE VOLUME CHANGE,
	// INT(RAND()*(N+1)) GIVES VALUES BETWEEN [0-N] EXCEPT FOR 
	// THE EXCEPTION AT RAND()=1,
	//
	if ( int(rand01()*(sys.N+1)) <= sys.N+1-2 )
	{
	  // ========================================
	  // DO PARTICLE MOVE

	  // CHOOSE PARTICLE AT RANDOM,
	  // BE CAREFUL WITH POSSIBILITY OF i=N THAT DOES NOT EXIST
	  i = MIN( int(rand01()*sys.N), sys.N-1 );
	  //PRINT2("doing particle move",i);
	  // CHECK TO SEE IF ATOM IS FIXED
	  if (sys.atoms[i].is_not_fixed())
	  {
		// ========================================
		// CALCULATE CURRENT ENERGY OF PARTICLE i
		eo = physics.calc_LJ_forces_single_particle(i, direction, virSumTemp, sys);
		// SAVE VIRIAL AND OLD POSITION
		virSum = virSumTemp;
		old_r  = sys.atoms[i].r;

		// ========================================
		// GENERATE RANDOM VECTOR
#if NDIM == 2
		rand_r = vec( (rand01()-0.5) * sys.mc_dr, (rand01()-0.5) * sys.mc_dr, 0. );
#else
		rand_r = vec( (rand01()-0.5) * sys.mc_dr, (rand01()-0.5) * sys.mc_dr, (rand01()-0.5) * sys.mc_dr );
#endif
		// GENERATE NEW TRIAL POSITION
		sys.atoms[i].r += rand_r;

		// ========================================
		// CALCULATE TRIAL ENERGY OF PARTICLE i
		en = physics.calc_LJ_forces_single_particle(i, direction, virSumTemp, sys);

		// ========================================
		// METROPOLIS ACCEPTANCE CRITERIA
		delta_e  = en - eo;
		delta_be = beta * delta_e;
		if (delta_be < 75.)
		{
		  if (delta_e < 0.)
		  {
			// ACCEPTED.
			// CHANGE ENERGY AND VIRIAL BY THEIR DIFFERENCES
			sys.pot_energy += delta_e;
			sys.virSum      = sys.virSum + virSumTemp - virSum;
			accepted_steps++;
		  }
		  else if ( exp(-delta_be) > rand01() )
		  {
			// ACCEPTED.
			// CHANGE ENERGY AND VIRIAL BY THEIR DIFFERENCES
			sys.pot_energy += delta_e;
			sys.virSum      = sys.virSum + virSumTemp - virSum;
			accepted_steps++;
		  }
		  else
		  {
			// REJECTED
			sys.atoms[i].r = old_r;
		  }
		}
		else
		{
		  // REJECTED
		  sys.atoms[i].r = old_r;
		}
	  }
	}
	else
	{
	  //PRINT1("doing volume move");
	  // ========================================
	  // DO VOLUME CHANGE

	  // ========================================
	  // CALCULATE CURRENT ENERGY AT CURRENT VOLUME
	  eo = physics.calc_LJ_forces(virSumTemp, sys);
	  // SAVE VIRIAL AND OLD POSITION
	  virSum = virSumTemp;
	  // SAVE OLD VOLUME AND RATIOS
#if NDIM == 2
	  vo     = sys.Lx * sys.Ly;
#else
	  vo     = sys.Lx * sys.Ly * sys.Lz;
#endif
	  rat_xy = sys.Lx / sys.Ly;
	  rat_xz = sys.Lx / sys.Lz;
	  lxo = sys.Lx;
	  lyo = sys.Ly;
	  lzo = sys.Lz;
	  //PRINT3("  current energy virial", eo, virSum);
	  //PRINT4("      distance A:",lxo,lyo,lzo);

	  // ========================================
	  // GENERATE RANDOM VOLUME CHANGE
	  lnvn = log(vo) + (rand01()-0.5) * delta_v;
	  // GENERATE NEW TRIAL LENGTHS
	  vn = exp(lnvn);
#if NDIM == 2
	  lx = pow(vn * rat_xy,0.5);
	  lz = lzo;
#else
	  lx = pow(vn * rat_xy * rat_xz,(1./3.));
	  lz = lx / rat_xz;
#endif
	  ly = lx / rat_xy;
	  //PRINT4("      distance B:",lx,ly,lz);
	  // READJUST BOX SIZES IN SIMULATION
	  sys.adjust_box_sizes(lx, ly, lz);
	  sys.scale_positions(lx/lxo, ly/lyo, lz/lzo);

	  // ========================================
	  // CALCULATE TRIAL ENERGY OF PARTICLE i
	  en = physics.calc_LJ_forces(virSumTemp, sys);

	  // ========================================
	  // METROPOLIS ACCEPTANCE CRITERIA,
	  // GIBBS FREE ENERGY
	  delta_e  = (en-eo) + sys.pressure*(vn-vo) - (sys.N+1)*log(vn/vo)/beta ;
	  delta_be = beta * delta_e;
	  //PRINT2("   delta_e",delta_e);
	  //PRINT7("distances",lxo,lx,lyo,ly,lzo,lz);
	  if (delta_be < 75.)
	  {
		if (delta_e < 0.)
		{
		  // ACCEPTED.
		  // ENERGY AND VIRIAL HAVE ALREADY BEEN RECALCULATED.
		  // UPDATE DENSITY
		  sys.calc_density();
		  accepted_steps++;
		}
		else if ( exp(-delta_be) > rand01() )
		{
		  // ACCEPTED.
		  // ENERGY AND VIRIAL HAVE ALREADY BEEN RECALCULATED.
		  // UPDATE DENSITY
		  sys.calc_density();
		  accepted_steps++;
		}
		else
		{
		  // REJECTED.
		  // GO BACK TO OLD ENERGY AND VIRIAL
		  sys.pot_energy = eo;
		  sys.virSum     = virSum;
		  // GO BACK TO OLD SIZES AND POSITIONS
		  sys.adjust_box_sizes(lxo, lyo, lzo);
		  sys.scale_positions(lxo/lx, lyo/ly, lzo/lz);
		}
	  }
	  else
	  {
		// REJECTED.
		// GO BACK TO OLD ENERGY AND VIRIAL
		sys.pot_energy = eo;
		sys.virSum     = virSum;
		// GO BACK TO OLD SIZES AND POSITIONS
		sys.adjust_box_sizes(lxo, lyo, lzo);
		sys.scale_positions(lxo/lx, lyo/ly, lzo/lz);
	  }
	}
  }

  // MAKE SURE THERE IS NOT KINETIC ENERGY
  sys.kin_energy  = 0.;

  // UPDATE TOTAL ENERGY
  //sys.energy = sys.pot_energy + sys.kin_energy;
  sys.energy = sys.pot_energy;

  return (accepted_steps * 1.0 / (sys.N+1));
}

double montecarlo::adjust_mc_step_size(double acceptance, double mc_dr)
{
  // ADJUST SIZE OF MC STEP BY LOOKING AT
  // ACCEPTANCE RATIO.
  double crit_change = 0.5;
  double per_decr = 0.80;
  double per_incr = 1.20;

  if (acceptance > crit_change)
  {
	// DOING TOO GOOD, INCREASE STEP
	mc_dr *= per_incr;
  }
  else if (acceptance < crit_change)
  {
	// NEED ENCOURAGEMENT, DECRASE STEP
	mc_dr *= per_decr;
  }

  // NOTE THAT IF == crit_change, LEAVE ALONE
  return mc_dr;
}

double montecarlo::adjust_mc_step_size_theta(double acceptance, double mc_dr, double& theta)
{
  // ADJUST SIZE OF MC STEP BY LOOKING AT
  // ACCEPTANCE RATIO.
  double crit_change = 0.5;
  double per_decr = 0.80;
  double per_incr = 1.20;

  if (acceptance > crit_change)
  {
	// DOING TOO GOOD, INCREASE STEP
	mc_dr *= per_incr;
	theta *= per_incr;
  }
  else if (acceptance < crit_change)
  {
	// NEED ENCOURAGEMENT, DECRASE STEP
	mc_dr *= per_decr;
	theta *= per_decr;
  }

  // ADJUST THETA TO LIMITS
  // IN RADIANS
  if (theta<=0.)
	theta = 0.009;
  else if (theta>.6)
	theta = 0.5236;

  // NOTE THAT IF == crit_change, LEAVE ALONE
  return mc_dr;
}


/*
 * ----------------------------------------------------------------------
 */
