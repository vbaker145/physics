
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

#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <iostream>

using namespace std;

class ensemble;
class forces;

class montecarlo
{
  friend ostream &operator<<(ostream &os, const montecarlo &z);

  public:
  	montecarlo(void);
  	~montecarlo(void);
  private:
  	montecarlo(const montecarlo &z);
  	montecarlo &operator=(const montecarlo &z);
  public:
	// DO ONE MC STEP
 	double mc_step(forces&, ensemble& sys);
	// DO ONE MC STEP IN NPT
	double mc_step_NPT(forces &physics, ensemble& sys);
	// DO ONE MC STEP OF LJ+XY MODEL
 	double mc_step_2D_angle(short, forces&, ensemble& sys);
	// ADJUST SIZE OF MC STEP BASED ON ACCEPTANCE
	double adjust_mc_step_size(double acceptance, double mc_dr);
	double adjust_mc_step_size_theta(double acceptance, double mc_dr, double& theta);

};

#endif // MONTECARLO_H
