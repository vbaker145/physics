
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

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>

using namespace std;

class ensemble;

class integrator
{
  friend ostream &operator<<(ostream &os, const integrator &z);

  public:
  	integrator(void);
  	~integrator(void);
  private:
  	integrator(const integrator &z);
  	integrator &operator=(const integrator &z);
  public:
	// VERLET ALGORITHM
	void verlet(ensemble& sys);

	// VELOCITY RESCALING FUNCTIONS
	void rescale_T_simple(ensemble& sys);
	void zero_velocities(ensemble& sys);
};

#endif // INTEGRATOR_H
