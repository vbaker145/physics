
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

#ifndef VERLETLIST_H
#define VERLETLIST_H

#include <iostream>

using namespace std;

class ensemble;
class vec;
class forces;

class verletlist
{
  friend ostream &operator<<(ostream &os, const verletlist &z);

  public:
  	verletlist(void);
  	~verletlist(void);
  private:
  	verletlist(const verletlist &z);
  	verletlist &operator=(const verletlist &z);
  public:
  	verletlist(long Np, double , double );
	// CHECKS TO SEE IF LIST NEEDS UPDATING
	int check_verletlist(ensemble&);
	// UPDATES THE LIST WITH NEW NEIGHBOR ID'S
	void update_verletlist(ensemble &sys);
	// RETURNS 1 IF PARTICLE j IS IN THE LIST
	// OF PARTICLE i
	int is_in_list(int i, int j);

	long N;
	long *number_list;			// ARRAY THAT HOLDS THE NUMBER OF NEIGHBOR
									// PARTICLES FOR A GIVEN PARTICLE
	vec *r_orig;					// ARRAY OF ORIGINAL POSITIONS AT THE MOMENT OF
									// THE SORTING OF THE LISTS
	double r_skin;					// RADIUS OF THE OUTER SKIN OF LIST
	double r_c;						// POTENTIAL CUTOFF DISTANCE
	long **list;					// VERLET NEIGHBOR LIST. 1ST INDEX IS PARTICLE ID,
									// 2ND INDEX RUNS OVER INTERNAL NEIGHBOR LABEL,
									// VALUE OF THE LIST IS THE ID OF THE NEIGHBOR.
};

#endif // VERLETLIST_H
