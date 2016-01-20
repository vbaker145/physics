
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

#ifndef MINIMA_H
#define MINIMA_H

#include <iostream>

using namespace std;

class ensemble;
class forces;

class minima
{
  friend ostream &operator<<(ostream &os, const minima &z);

  public:
  	minima(void);
  	~minima(void);
  private:
  	minima(const minima &z);
  	minima &operator=(const minima &z);
  public:

	// FUNCTION THAT MOVES INDIVIDUAL PARTICLES ALONG THE DIRECTION
	// OF THE TOTAL FORCE (STEEPEST DESCENT) UNTIL A LOCAL MINIMA
	// IN THE ENERGY IS FOUND.
	// WORKS BY USING AN ARBITRARY STEP.
	// 	FOR A GIVEN PARTICLE
	void steepest_desdents_arbitrary_step_one_particle(int particle, double step_size, double crit, forces &physics, ensemble& sys);
	// 	FOR ALL PARTICLES
	void steepest_desdents_arbitrary_step(double step_size, double crit, forces &physics, ensemble& sys);
};

#endif // MINIMA_H
