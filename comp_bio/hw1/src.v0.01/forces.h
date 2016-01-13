
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


#ifndef FORCES_H
#define FORCES_H

#include <iostream>

using namespace std;

class ensemble;

class forces
{
  friend ostream &operator<<(ostream &os, const forces &z);

  public:
  	forces(void);
  	~forces(void);
  private:
  	forces(const forces &z);
  	forces &operator=(const forces &z);
  public:
	void calc_LJ_forces(ensemble&);
};

#endif // FORCES_H
