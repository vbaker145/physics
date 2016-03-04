
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


#ifndef PAIRCORRFUNC_H
#define PAIRCORRFUNC_H

#include <iostream>

using namespace std;

class ensemble;
class vec;

class paircorrfunc
{
  friend ostream &operator<<(ostream &os, const paircorrfunc &z);

  private:
  	paircorrfunc(void);
  public:
  	~paircorrfunc(void);
  private:
  	paircorrfunc(const paircorrfunc &z);
  	paircorrfunc &operator=(const paircorrfunc &z);
  public:
  	paircorrfunc(ensemble&);	// CONSTRUCTOR

	void accum(int n, ensemble&);	// PERIODICALLY CALLED TO UPDATE
									// ARRAYS
	void reset(ensemble&);					// ZEROS THE ARRAYS

	double duration;	// DURATION TIME OF EACH MEASUREMENT
	int i_duration;		// DURATION STEPS OF EACH MEASUREMENT
	int howmany;		// HOW MANY TIMES WE WILL SAMPLE 
						// GIVEN A duration
	int curr_time;		// COUNTER FOR THE CURRENT LOCAL STEP NO.
	int nhis;			// NUMBER OF HISTOGRAM BINS
	double minL;		// MINIMUM SIZE OF THE CUBE
	double delta_gr;	// SIZE OF THE BIN

	int *num_bins;		// ARRAY OF THE NUMBER OF COUNTED BINS PER DISTANCE
	double *gr;		// ARRAY OF THE PAIR CORRELATION FUNCTION
};

#endif // PAIRCORRFUNC_H
