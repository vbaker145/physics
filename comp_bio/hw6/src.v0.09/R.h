
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


#ifndef R_MEASURE_H
#define R_MEASURE_H

#include <iostream>

using namespace std;

class ensemble;
class vec;

class R_meas
{
  friend ostream &operator<<(ostream &os, const R_meas &z);

  private:
  	R_meas(void);
  public:
  	~R_meas(void);
  private:
  	R_meas(const R_meas &z);
  	R_meas &operator=(const R_meas &z);
  public:
  	R_meas(ensemble&);	// CONSTRUCTOR

	void calc(ensemble&, ofstream&, ofstream&, ofstream&, ofstream&, ofstream&);	
					// PERIODICALLY CALLED TO CALCULATE AND
					// OUTPUT THE END-TO-END DISTANCE AND
					// RADIUS OF GYRATION
	double *Rg;		// RADIUS OF GYRATION
	double *Rg2;	// SQUARE OF RADIUS OF GYRATION
	double *Eoe; 	// END-TO-END DISTANCE
	vec    *Eoe_vec;// END-TO-END VECTOR
};

#endif // R_MEASURE_H
