
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
class verletlist;
class vec;

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
	double calc_LJ_forces(double&, ensemble&);
	double calc_LJ_forces_single_particle(int ii, vec &direction, double&, ensemble& sys);
	void calc_LJ_forces_verletlists(verletlist &vl, ensemble& sys);

	void calc_LJ_forces_pot_cutoff(ensemble&);

	double calc_LJXY_forces_single_particle(int ii, vec &direction, double &virSum, ensemble& sys);
	double calc_LJLL_forces_single_particle(int ii, vec &direction, double &virSum, ensemble& sys);

	void energy_landscape(ensemble& sys,char*);
	double calc_LJ_forces_energy_landscape(vec &position, ensemble& sys);
};

#endif // FORCES_H
