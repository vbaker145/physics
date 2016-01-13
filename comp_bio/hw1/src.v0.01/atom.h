
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

#ifndef ATOM_H
#define ATOM_H

#include <iostream>

#include "vec.h"

using namespace std;

class atom
{
  friend ostream &operator<<(ostream &os, const atom &z);

  public:
  	atom(void);
  	~atom(void);
  private:
  	atom(const atom &z);
  	atom &operator=(const atom &z);
  public:
	void output_position(void);
	void output_velocity(void);

	double m;					// ATOM MASS
	double d;                   // DIAMETER
	vec r;						// POSITION
	vec r_old;					// OLD POSITION
	vec v; 						// VELOCITIES
	vec f;						// FORCE
};

#endif // ATOM_H
