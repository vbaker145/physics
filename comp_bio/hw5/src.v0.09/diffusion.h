
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


#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <iostream>

using namespace std;

class ensemble;
class vec;

class diffusion
{
  friend ostream &operator<<(ostream &os, const diffusion &z);

  private:
  	diffusion(void);
  public:
  	~diffusion(void);
  private:
  	diffusion(const diffusion &z);
  	diffusion &operator=(const diffusion &z);
  public:
  	diffusion(ensemble&);	// CONSTRUCTOR

	void accum(int n, ensemble&);	// PERIODICALLY CALLED TO UPDATE
									// ARRAYS
	void reset(ensemble&);					// ZEROS THE ARRAYS

	double duration;	// DURATION TIME OF EACH MEASUREMENT
	int i_duration;		// DURATION STEPS OF EACH MEASUREMENT
	int howmany;		// HOW MANY TIMES WE WILL SAMPLE 
						// GIVEN A duration
	int curr_time;		// COUNTER FOR THE CURRENT LOCAL STEP NO.

	double *r2t;		// ARRAY OF THE MEAN SQUARE DISPLACEMENTS
	double *vacf;		// ARRAY OF THE VELOCITY AUTOCORRELATION FUNCTION
	vec *rt0;			// HOLDS THE INITIAL VALUE OF POSITION
	vec *vt0;			// HOLDS THE INITIAL VALUE OF VELOCITY
};

#endif // DIFFUSION_H
