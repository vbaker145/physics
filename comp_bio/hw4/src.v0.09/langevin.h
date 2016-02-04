
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

#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <iostream>

using namespace std;

class ensemble;
class forces;

class langevin
{
  friend ostream &operator<<(ostream &os, const langevin &z);

  public:
  	langevin(void);
  	~langevin(void);
  private:
  	langevin(const langevin &z);
  	langevin &operator=(const langevin &z);
  public:
	// INTEGRATION ALGORITHM FROM HINCHLIFFE
 	void integrate_ermak(ensemble& sys);

	// INTEGRATION ALGORITHM WITH FORCE NOT CONSTANT
	// BETWEEN TIME STEPS
	void integrate_other_pt1(ensemble& sys);
	void integrate_other_pt2(ensemble& sys);
};

#endif // LANGEVIN_H
