
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


#ifndef STEPWISE_H
#define STEPWISE_H

#include <iostream>

using namespace std;

class ensemble;
class verletlist;
class vec;

class stepwise
{
  friend ostream &operator<<(ostream &os, const stepwise &z);

  public:
  	stepwise(void);
  	~stepwise(void);
  private:
  	stepwise(const stepwise &z);
  	stepwise &operator=(const stepwise &z);
  public:
	void calc_hardsphere_stepwise(ensemble&);
};

#endif // STEPWISE_H
